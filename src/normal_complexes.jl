using Polymake: polytope
using Combinatorics: powerset, permutations, levicivita
using LinearAlgebra: I, det
using Symbolics

# Normal complexes: the combinatorial hypercubes cut out of simplicial cones by
# a vector of heights h, their volumes, and their lattice point counts.
#
# Conventions: a cone is given by the matrix R whose rows are its rays, and h
# assigns a height h_i to the i-th ray. The region of the cone cut out by
# ⟨r_i, w⟩ <= h_i is a combinatorial hypercube whose corners are indexed by
# the subsets J of rays, encoded as the binary digits of an integer k.


####################
# COMBINATORIAL HYPERCUBES
####################

# Corner k of the combinatorial hypercube cut out by h inside the cone with rays
# the rows of R. Corner k corresponds to the set J of rays given by the binary
# digits of k, and is the point w in span(R_J) with R_J * w = h_J.
function hypercube_corner(R::AbstractMatrix{Rational{T}}, h, k) where T
    J = [(k >> i) % 2 == 1 for i in 0:size(R, 1)-1]
    # Corner 0 is the origin
    if iszero(J)
        return zero(R[1, :])
    end
    # Restrict h and R to the rays in J
    h_ = h[J]
    R_ = R[J, :]
    # Find the point that projects to h_ in the span of R_
    return R_' * inv(R_ * R_') * h_
end

# All 2^d corners of the combinatorial hypercube, in the order of hypercube_corner
function hypercube_corners(R::AbstractMatrix{Rational{T}}, h) where T <: Integer
    return map(k -> hypercube_corner(R, h, k), 0:2^size(R, 1) - 1)
end

# Signed volume of the simplex spanned by the origin and the rows of R
function simplex_signed_volume(R)
    return det(R) / factorial(size(R, 1))
end

# Volume of the combinatorial hypercube with vertices verts (indexed as in hypercube_corners).
# Uses the staircase triangulation of the cube: for each permutation σ, the
# simplex with vertices 0, e_σ1, e_σ1 + e_σ2, ..., weighted by the sign of σ so
# that all simplices are counted with the same orientation.
function hypercube_volume(verts)
    n = Int(log2(length(verts)))
    vol = 0
    for flag in permutations(1:n)
        edges = map(1:n) do i
            vert_index = sum(1 .<< (flag[1:i] .- 1)) + 1
            verts[vert_index] - verts[1]
        end
        R = reduce(hcat, edges)'
        vol += levicivita(flag) * simplex_signed_volume(R)
    end
    return vol
end


####################
# VOLUME POLYNOMIAL
####################

# Volume of the combinatorial hypercube in the cone with rays the rows of R,
# as a polynomial in vars, measured in the coordinates given by pairing with the rays
function cone_volume_polynomial(R, vars)
    verts = hypercube_corners(R, vars) .|> w -> Vector(R * w)
    sgn = sign(det(R * R'))

    return sgn * hypercube_volume(verts)
end

# Volume polynomial of the matroid M in the variables h_F: the sum of the cone
# volume polynomials over the maximal cones of the Bergman fan
function volume_polynomial(M)
    dim = maximum(length.(chains(M)))
    cones = filter(c -> length(c) == dim, chains(M))

    vol = 0
    for cone in cones
        R = bergman_cone_rays(M, cone)
        vars = Symbolics.variables(:h, cone)

        vol += cone_volume_polynomial(R, vars)
    end
    return Symbolics.expand(vol)
end

# Snapper polynomial of M computed from its volume polynomial: write the volume
# polynomial in the f basis, replace each monomial Π_F f_F^(a_F) by the product
# of rising factorials Π_F f_F (f_F + 1) ... (f_F + a_F - 1), and translate by
# f_E ↦ f_E + 1, with E the ground set. Agrees with snapper_polynomial(M).
function snapper_polynomial_from_volume(M)
    f = Symbolics.variables(:f, nonempty_flats(M))
    f_E = Symbolics.variable(:f, flat_index(M, Set(ground_set(M))))

    vol_f = Symbolics.expand(Symbolics.substitute(volume_polynomial(M), h_to_f(M)))
    rising = to_rising_factorials(vol_f, f)
    return Symbolics.expand(Symbolics.substitute(rising, Dict(f_E => f_E + 1)))
end


####################
# EHRHART POLYNOMIAL OF A CONE
####################

# Matrix whose rows are the vectors vs
function to_matrix(vs)
    return reduce(vcat, [v' for v in vs])
end

# Primitive integer vector on the ray through r
function primitive_vector(r)
    return r // gcd(r...)
end

# Tangent cones (rays as rows) at the vertices of a combinatorial hypercube with
# vertices verts (indexed as in hypercube_corners). The rays at a vertex are the primitive
# edge directions towards its d neighbours, which differ from it in one binary digit.
function hypercube_tangent_cones(verts::AbstractVector{<:AbstractVector{Rational{T}}}) where T <: Integer
    dim = Int(log2(length(verts)))
    cones = Matrix{Rational{T}}[]
    for i in 1:length(verts)
        base = verts[i]
        rays = []
        for j in 0:dim-1
            v = verts[xor(i - 1, 1 << j) + 1]
            push!(rays, primitive_vector(v - base))
        end
        push!(cones, to_matrix(rays))
    end
    return cones
end

# Polymake cone of the heights h for which the region cut out by h inside the
# cone with rays the rows of R is a genuine (pseudo-cubical) hypercube:
# h >= 0, and each corner w_J lies on the correct side of the remaining facets.
function pseudo_cubical_cone(R)
    n = size(R, 1)
    inequalities = []
    push!(inequalities, Matrix{Rational}(I, n, n))
    # We need each facet of the polytope to be a hypercube
    # so the orthogonal projection of the vertex on each cone must lie in its
    # relative interior
    for J in powerset(1:n, 1, n-1)
        # projection matrix onto subspace generated by facet J
        # For any w in the span, with Rw = h, we have w = Ah
        A = R[J, :]' * inv(R[J, :] * R[J, :]')
        for j in setdiff(1:n, J)
            # These inequalities come from Proposition 3.1, which essentially state that w_J belongs to the polytope cut out by the half-planes
            # The original inequality is r_j * w_J <= h_j for all j not in J
            # We have w_J = A_J * h_J and here we also take the transpose since we work with rows
            ineq = zeros(Rational, n)
            ineq[j] = 1
            ineq[J] = -A' * R[j, :]
            push!(inequalities, ineq)
        end
    end
    return polytope.Cone(INEQUALITIES=reduce(hcat, inequalities)')
end

# Number of lattice points of the combinatorial hypercube cut out by the
# heights vars inside the cone with rays the rows of R, as a polynomial in vars.
# Points are measured in the coordinates B' * w. The tangent cones at the
# vertices are taken from the hypercube at the pseudo-cubical heights h_pc
# (see pseudo_cubical_cone), which has the same combinatorial type.
function cone_ehrhart_polynomial(B::AbstractMatrix{Rational{T}}, R::AbstractMatrix{Rational{T}}, vars, h_pc) where T <: Integer
    verts_pc = hypercube_corners(R, h_pc) .|> w -> Vector{Rational{T}}(B' * w)
    cones = hypercube_tangent_cones(verts_pc)
    verts = hypercube_corners(R, vars) .|> w -> Vector(B' * w)

    return lattice_count(cones, verts)
end


####################
# EHRHART POLYNOMIAL OF A MATROID
####################

# Pseudo-cubical heights for each cone (chain of flats) of the Bergman fan of M:
# a relative interior point of pseudo_cubical_cone for the rays of the cone
function pseudo_cubical_heights(M; B=I)
    return Dict(cone => Vector{Rational}(pseudo_cubical_cone(bergman_cone_rays(M, cone; B=B)).REL_INT_POINT)
                for cone in chains(M))
end

# Lattice point counts of the faces of the hypercubes of the Bergman fan of M.
# The entry (cone, face) is the count for the hypercube of the face (a subchain
# of cone), measured in the coordinates given by pairing with the rays of cone.
function face_ehrhart_polynomials(M; B=I, h_pc=pseudo_cubical_heights(M; B=B))
    polynomials = Dict()
    for cone in chains(M)
        r = bergman_cone_rays(M, cone; B=B)
        vars = Symbolics.variables(:h, cone)
        for J in powerset(1:length(cone), 1)
            face = cone[J]
            polynomials[(cone, face)] = cone_ehrhart_polynomial(r', r[J, :], vars[J], h_pc[face])
        end
    end
    return polynomials
end

# Contribution of a cone c of dimension 1 or 2 to the Ehrhart polynomial, where
# P = face_ehrhart_polynomials(M):
#   dimension 1:  P(c, c) - 1
#   dimension 2:  P(c, c) - 1/2 Σ_{rays ρ of c} (P(c, ρ) + P(ρ, ρ)) + 1
function ehrhart_contribution(cone, polynomials)
    if length(cone) == 1
        return polynomials[(cone, cone)] - 1
    elseif length(cone) == 2
        result = polynomials[(cone, cone)] + 1
        for ray in (cone[1:1], cone[2:2])
            result -= 1//2 * (polynomials[(cone, ray)] + polynomials[(ray, ray)])
        end
        return result
    else
        error("ehrhart_contribution is only implemented for cones of dimension at most 2")
    end
end

# Ehrhart polynomial of the normal complex of M in the variables h_F: the sum of
# the contributions of all cones of the Bergman fan, plus 1 for the origin.
# Only implemented when all cones have dimension at most 2 (M of rank at most 3).
function ehrhart_polynomial(M; B=I, polynomials=nothing)
    cones = chains(M)
    if maximum(length.(cones)) > 2
        error("ehrhart_polynomial is only implemented for Bergman fans with cones of dimension at most 2")
    end
    if polynomials === nothing
        polynomials = face_ehrhart_polynomials(M; B=B)
    end

    ehr = 1 + sum(ehrhart_contribution(cone, polynomials) for cone in cones)
    return Symbolics.expand(ehr)
end
