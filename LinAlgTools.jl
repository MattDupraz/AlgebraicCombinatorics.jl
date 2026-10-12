using LinearAlgebra, DataStructures
using Base.Iterators

# BinomialSeries and Series (used in brion_coefficient), unless already included
isdefined(@__MODULE__, :BinomialSeries) || include("Series.jl")

# Exact linear algebra over ℚ and ℤ, lattice reduction, and Barvinok-style
# lattice point counting for rational cones and polygons.
#
# Conventions: vectors are stored as the ROWS of a matrix (a cone is given by
# the matrix whose rows are its rays). Integer input is promoted to rationals
# wherever division is needed, so all results are exact.

# Element type in which division is exact (integers are promoted to rationals)
field_type(T) = T <: Integer ? Rational{T} : T


####################
# GRAM-SCHMIDT AND LLL
####################

# Gram-Schmidt orthogonalisation of the rows b_i of B, without normalisation.
# Returns (B_star, mu) where the rows b*_i of B_star are pairwise orthogonal and
#   b_i = b*_i + Σ_{j<i} mu[i, j] b*_j.
# Assumes the rows of B are linearly independent.

function gram_schmidt(B::AbstractMatrix{T}) where T
    m, n = size(B)
    F = field_type(T)
    B_star = Matrix{F}(B)
    mu = zeros(F, (m, m))
    norms = zeros(F, m)  # |b*_j|^2

    for i in 1:m
        for j in 1:i-1
            mu[i, j] = dot(@view(B[i, :]), @view(B_star[j, :])) / norms[j]
            @views B_star[i, :] .-= mu[i, j] .* B_star[j, :]
        end
        norms[i] = dot(@view(B_star[i, :]), @view(B_star[i, :]))
    end

    return B_star, mu
end

# LLL reduction of the lattice spanned by the rows of B.
# Returns a reduced basis B' and a unimodular matrix U such that U * B = B'.
# delta ∈ (1/4, 1) is the Lovász parameter.

function lll_reduction(B::AbstractMatrix{T}, delta=3//4) where T
    B = copy(B)
    m = size(B, 1)
    delta = Rational{BigInt}(delta)
    # Gram-Schmidt coefficients grow quickly, so compute them in BigInt to avoid overflow.
    # They are computed once and then updated incrementally (Cohen, Algorithm 2.6.3).
    B_star, mu = gram_schmidt(big.(B))
    norms = [dot(r, r) for r in eachrow(B_star)]  # |b*_i|^2
    U = Matrix{T}(I, m, m)  # Unimodular transformation matrix
    k = 2

    while k <= m
        # Size reduction: b_k -= q b_j, making |mu[k, j]| <= 1/2
        for j in (k-1):-1:1
            if abs(mu[k, j]) > 1//2
                q = round(mu[k, j])
                @views B[k, :] .-= T(q) .* B[j, :]
                @views U[k, :] .-= T(q) .* U[j, :]
                mu[k, j] -= q
                for i in 1:j-1
                    mu[k, i] -= q * mu[j, i]
                end
            end
        end

        # Lovász condition
        if norms[k] >= (delta - mu[k, k-1]^2) * norms[k-1]
            k += 1
        else
            # Swap b_k and b_{k-1}, then update mu and norms in place
            B[[k, k-1], :] = B[[k-1, k], :]
            U[[k, k-1], :] = U[[k-1, k], :]
            for j in 1:k-2
                mu[k, j], mu[k-1, j] = mu[k-1, j], mu[k, j]
            end
            nu = mu[k, k-1]
            new_norm = norms[k] + nu^2 * norms[k-1]
            mu[k, k-1] = nu * norms[k-1] / new_norm
            norms[k] = norms[k-1] * norms[k] / new_norm
            norms[k-1] = new_norm
            for i in k+1:m
                t = mu[i, k]
                mu[i, k] = mu[i, k-1] - nu * t
                mu[i, k-1] = t + mu[k, k-1] * mu[i, k]
            end
            k = max(k-1, 2)
        end
    end

    return B, U
end

# Approximate shortest vector of the lattice spanned by the rows of B: the
# vector of an LLL-reduced basis with the smallest ∞-norm.
# Returns (v, c) with v = c' * B and c an integer vector.

function shortest_vector_approximation(B::AbstractMatrix{T}) where T
    reduced_basis, U = lll_reduction(B)

    infinity_norms = [maximum(abs, reduced_basis[i, :]) for i in 1:size(reduced_basis, 1)]
    shortest_index = argmin(infinity_norms)
    shortest_vector = reduced_basis[shortest_index, :]
    coordinates = U[shortest_index, :]

    return shortest_vector, coordinates
end


####################
# ROW REDUCTION AND SUBSPACES
####################

# Reduced row echelon form, computed exactly over field_type(T)

function reduced_row_echelon_form(A::AbstractMatrix{T}) where T
    A = Matrix{field_type(T)}(A)  # Work on a copy to avoid modifying the original matrix
    rows, cols = size(A)
    lead = 1  # Column of the current pivot

    for r in 1:rows
        if lead > cols
            break
        end

        # Find the next pivot: the first row i >= r with A[i, lead] != 0,
        # moving to the next column if this one is zero below row r
        p = r
        while A[p, lead] == 0
            if p == rows
                p = r
                lead += 1
                if lead > cols
                    return A
                end
            else
                p += 1
            end
        end

        # Swap rows
        A[[r, p], :] = A[[p, r], :]

        # Normalize the pivot row
        lv = A[r, lead]
        @views A[r, :] ./= lv

        # Clear the rest of the pivot column
        for i in 1:rows
            if i != r && !iszero(A[i, lead])
                f = A[i, lead]
                @views A[i, :] .-= f .* A[r, :]
            end
        end

        lead += 1
    end

    return A
end

# Returns a maximal linearly independent subset of the rows of vs.
# A row v is kept if it is not fixed by the projection P = B'(BB')⁻¹B onto
# the span of the rows B kept so far.

function extract_basis(vs::AbstractMatrix{T}) where T
    B = vs[1:0, :]
    for row in 1:size(vs, 1)
        v = vs[row:row, :]
        iszero(v) && continue
        if isempty(B) || v * B' * inv(Matrix{field_type(T)}(B * B')) * B != v
            B = [B; v]
        end
    end
    return B
end

# Basis of the row space of A: the non-zero rows of its RREF

function row_space_basis(A::AbstractMatrix{T}) where T
    ref_A = reduced_row_echelon_form(A)
    nonzero_rows = [i for i in 1:size(ref_A, 1) if !iszero(@view ref_A[i, :])]
    return ref_A[nonzero_rows, :]
end

# Basis (as rows) of the orthogonal complement of the row space of A.
# Projects each standard basis vector onto the complement and row-reduces the result.

function orthogonal_complement(A::AbstractMatrix{T}) where T
    m, n = size(A)
    A, _ = gram_schmidt(A)
    Q = Matrix{eltype(A)}(I, n, n)  # Start with the identity matrix as the initial basis
    norms = [dot(r, r) for r in eachrow(A)]

    # Get orthogonal part of each vector
    for i in 1:n
        for j in 1:m
            c = dot(@view(A[j, :]), @view(Q[i, :])) / norms[j]
            @views Q[i, :] .-= c .* A[j, :]
        end
    end

    return row_space_basis(Q)
end


####################
# HERMITE NORMAL FORM AND LATTICES
####################

# Column-style Hermite normal form of an integer matrix A.
# Returns (H, U) with U unimodular and A * U = H, where H is in column echelon
# form: each pivot is positive, entries right of a pivot are zero, and entries
# left of a pivot lie in [0, pivot).

function hermite_normal_form(A::AbstractMatrix{T}) where T <: Integer
    H = copy(A)
    m, n = size(H)
    U = Matrix{T}(I, n, n)  # Unimodular transformation matrix

    pivot = 1 # Row of the current pivot
    j = 1     # Column of the current pivot

    while j <= n && pivot <= m
        # Euclidean algorithm on the columns: repeat until H[pivot, j+1:end] is zero
        while j < n && !iszero(@view H[pivot, j+1:end])
            # Find the column with the smallest non-zero entry in the current row
            min_col = j
            for i in j+1:n
                if H[pivot, i] != 0 &&
                        (H[pivot, min_col] == 0 || abs(H[pivot, i]) < abs(H[pivot, min_col]))
                    min_col = i
                end
            end

            # Swap columns to bring the smallest entry to the pivot
            if min_col != j
                H[:, [j, min_col]] = H[:, [min_col, j]]
                U[:, [j, min_col]] = U[:, [min_col, j]]
            end

            # Make the pivot positive
            if H[pivot, j] < 0
                @views H[:, j] .*= -1
                @views U[:, j] .*= -1
            end

            # Reduce entries to the right of the pivot
            for i in j+1:n
                if H[pivot, i] != 0
                    q = div(H[pivot, i], H[pivot, j])
                    @views H[:, i] .-= q .* H[:, j]
                    @views U[:, i] .-= q .* U[:, j]
                end
            end
        end

        if H[pivot, j] != 0
            # Make the pivot positive
            if H[pivot, j] < 0
                @views H[:, j] .*= -1
                @views U[:, j] .*= -1
            end
            # Reduce entries to the left of the pivot into [0, pivot)
            for i in 1:j-1
                if H[pivot, i] != 0
                    q = fld(H[pivot, i], H[pivot, j])
                    @views H[:, i] .-= q .* H[:, j]
                    @views U[:, i] .-= q .* U[:, j]
                end
            end
            j += 1
            pivot += 1
        else
            # No pivot in this row; try the next row in the same column
            pivot += 1
        end
    end

    return H, U
end

# Basis (as rows) of the lattice ℤⁿ ∩ span(rows of B).
# This lattice is the integer kernel of the normal vectors B_perp, which is
# spanned by the columns of U corresponding to the zero columns of the HNF of B_perp.

function lattice_basis_inside_span(B::AbstractMatrix{Rational{T}}) where T <: Integer
    B = row_space_basis(B)
    r, n = size(B)  # r = dimension of the span

    if r == n
        return Matrix{T}(I, r, r)
    end

    # Integer normal vectors (dividing by the gcd clears denominators)
    B_perp = orthogonal_complement(B)
    B_perp = Matrix{T}(B_perp//gcd(B_perp))

    # B_perp has rank n - r, so the last r columns of H are zero
    H, U = hermite_normal_form(B_perp)
    return U[:, n-r+1:n]'
end


####################
# CONES AND BARVINOK DECOMPOSITION
####################

# A cone (rays as rows) with a sign, for signed decompositions
struct SignedCone{T}
    matrix::Matrix{Rational{T}}
    sign::Int
end
SignedCone(matrix::AbstractMatrix{Rational{T}}, sign) where T = SignedCone{T}(matrix, sign)

# A signed cone together with the determinant of its ray matrix (its index,
# up to sign, for a full-dimensional simplicial cone)
struct ConeData{T}
    matrix::Matrix{Rational{T}}
    sign::Int
    det::Rational{T}
end
ConeData(matrix::AbstractMatrix{Rational{T}}, sign, det) where T = ConeData{T}(matrix, sign, det)

# Dual of the full-dimensional simplicial cone with rays the rows of A.
# The i-th dual ray is the primitive integer normal of the facet spanned by all
# rays except the i-th, oriented so that it pairs positively with the i-th ray.

function simplicial_cone_dual(A::AbstractMatrix{Rational{T}}) where T <: Integer
    if size(A, 1) != size(A, 2)
        error("matrix defining cone must be square")
    end

    A_det = det(A)

    rays = Matrix{Rational{T}}[]
    for i in 1:size(A, 1)
        facet = A[setdiff(1:end, i), :]
        normal = orthogonal_complement(facet)
        ray_gen = normal//gcd(normal)

        # Orient the normal to point to the same side as the i-th ray
        if sign(det([A[1:i-1, :]; ray_gen; A[i+1:end, :]])) != sign(A_det)
            ray_gen = -ray_gen
        end

        push!(rays, ray_gen)
    end

    return reduce(vcat, rays)
end

# Barvinok's signed decomposition of the simplicial cone C (rays as rows) into
# unimodular cones, modulo lower-dimensional cones.
# The cone is first written in a lattice basis of its span, making it
# full-dimensional. The decomposition is done on the dual cone, where the
# lower-dimensional cones that get dropped do not affect lattice point counts
# of the primal; the resulting unimodular cones are then dualised back.
# Assumes that C is simplicial.

function barvinok_decomposition(C::AbstractMatrix{Rational{T}}) where T <: Integer
    B = Matrix{Rational{T}}(lattice_basis_inside_span(C))
    # Basis change into basis B
    P = B' * inv(B * B')
    # Cone in coordinates of lattice basis (full-dimensional now)
    A = C * P

    if size(A, 1) != size(A, 2)
        error("cone is not simplicial")
    end

    if abs(det(A)) == 1
        return [SignedCone(C, 1)]
    end

    A_dual = simplicial_cone_dual(A)

    # Subdivide until every cone is unimodular
    cone_queue = Queue{ConeData{T}}()
    enqueue!(cone_queue, ConeData{T}(A_dual, 1, det(A_dual)))

    unimodular_cones = ConeData{T}[]

    while !isempty(cone_queue)
        cone = dequeue!(cone_queue)
        subcones = subdivide_cone_once(cone)
        for subcone in subcones
            if abs(subcone.det) == 1
                push!(unimodular_cones, subcone)
            else
                enqueue!(cone_queue, subcone)
            end
        end
    end

    # Dualise back and return to the original coordinates
    primal_cones = unimodular_cones .|>
        dual_cone -> SignedCone(simplicial_cone_dual(dual_cone.matrix) * B,
                               dual_cone.sign)

    return primal_cones
end

# One step of Barvinok's subdivision of a full-dimensional simplicial cone with
# rays a_i (rows of A). Finds a short lattice vector z = Σ λ_i a_i (small λ_i,
# via LLL on A⁻¹), and replaces each ray a_i by z in turn. Each resulting cone
# has |det| = |λ_i| |det A|, which is smaller than |det A|.
# Cones with zero determinant are lower-dimensional and are dropped.

function subdivide_cone_once(cone::ConeData{T}) where T <: Integer
    if iszero(cone.det)
        error("cone is not full-dimensional")
    end

    A = cone.matrix
    A_inv = inv(A)
    lambda, z = shortest_vector_approximation(A_inv)
    # Prefer z on the same side as the cone
    if all(lambda .<= 0)
        lambda = -lambda
        z = -z
    end

    cones = ConeData{T}[]
    for i in 1:size(A, 1)
        C = copy(A)
        C[i, :] = z
        det_C = det(C)
        if det_C != 0
            # Sign is + if C has the same orientation as A, - otherwise
            sign_C = cone.sign * Int(sign(cone.det) * sign(det_C))
            push!(cones, ConeData(C, sign_C, det_C))
        end
    end
    return cones
end


####################
# LATTICE POINT COUNTING
####################

# Contribution of the unimodular cone `cone` shifted to the integer point
# `vertex` to the lattice point count, via Brion's theorem.
# The generating function x^v / Π (1 - x^{r_i}) is specialised along lambda by
# x^m ↦ (1 + t)^⟨m, lambda⟩, and the constant term of the resulting Laurent
# series in t is returned.

function brion_coefficient(cone::SignedCone{T}, vertex, lambda) where T <: Integer
    A = cone.matrix
    num = BinomialSeries(dot(vertex, lambda))
    denom = Series(1)
    for i in 1:size(A, 1)
        ray = A[i, :]
        denom *= 1 - BinomialSeries(dot(ray, lambda))
    end
    # denom = t^offset * (invertible series), so num / denom has a pole of order offset
    offset = first_nonzero_index(denom)
    quot = num / (denom << offset)
    return cone.sign * quot[offset]
end

# Enumerates ℕ^dim: the binary digits of k are dealt out to the coordinates in
# turn, so every non-negative integer vector is reached for some k

function get_an_integer_vector(k, dim)
    result = zeros(Int, dim)
    e = 0
    while k != 0
        for i in 1:dim
            result[i] += 2^e * (k % 2)
            k >>= 1
        end
        e += 1
    end
    return result
end

# Finds an integer vector lambda with ⟨r, lambda⟩ != 0 for every row r of A,
# so that no denominator in brion_coefficient vanishes

function find_good_lambda(A)
    k = 1
    while true
        lambda = get_an_integer_vector(k, size(A, 2))
        if all(A * lambda .!= 0)
            return lambda
        end
        k += 1
    end
end

# Number of lattice points of the polytope with the given integer vertices,
# where cones[i] is the tangent cone at vertices[i] (rays as rows)

function compute_lattice_count(cones::AbstractVector{<:AbstractMatrix{Rational{T}}},
        vertices::AbstractVector{<:AbstractVector}) where T <: Integer
    decompositions = cones .|> barvinok_decomposition
    lambda = find_good_lambda(reduce(vcat, [cone.matrix for cone in flatten(decompositions)]))
    result = 0
    for (decomposition, vertex) in zip(decompositions, vertices)
        for cone in decomposition
            result += brion_coefficient(cone, vertex, lambda)
        end
    end
    return result
end

# Number of lattice points of a convex lattice polygon with vertices given in
# cyclic order. The tangent cone at each vertex is spanned by the primitive
# edge directions towards its two neighbours.

function polytope_lattice_count(vertices)
    n = length(vertices)
    cones = map(1:n) do i
        r1 = vertices[(n+i-2)%n + 1] - vertices[i]  # Previous vertex
        r2 = vertices[i%n + 1] - vertices[i]        # Next vertex
        r1 = r1//gcd(r1)
        r2 = r2//gcd(r2)
        [r1' ; r2']
    end
    return compute_lattice_count(cones, vertices)
end
