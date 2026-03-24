module LinAlgTools

using ..Series
using LinearAlgebra, DataStructures
using Base.Iterators

export gram_schmidt, lll_reduction, shortest_vector_approximation, extract_basis,
    row_space_basis, orthogonal_complement, hermite_normal_form, lattice_basis_inside_span,
    barvinok_decomposition, simplicial_cone_dual, compute_lattice_count
export SignedCone

function gram_schmidt(B::AbstractMatrix{T}) where T
    m, n = size(B)
    B_star = copy(B)
    mu = zeros(T, (m, m))

    for i in 1:m
        for j in 1:i-1
            mu[i, j] = dot(B[i, :], B_star[j, :]) / dot(B_star[j, :], B_star[j, :])
            B_star[i, :] -= mu[i, j] * B_star[j, :]
        end
    end

    return B_star, mu
end

# Returns a reduced basis of B (say B') and a unimodular matrix U such that U * B = B'
function lll_reduction(B::AbstractMatrix{T}, delta=0.75) where T
    B = copy(B)
    n = size(B, 2)
    B_star, mu = gram_schmidt(B)
    U = Matrix{T}(I, n, n)  # Unimodular transformation matrix
    k = 2

    while k <= n
        for j in (k-1):-1:1
            if abs(mu[k, j]) > 0.5
                B[k, :] -= round(mu[k, j]) * B[j, :]
                U[k, :] -= round(mu[k, j]) * U[j, :]
                B_star, mu = gram_schmidt(B)
            end
        end

        lhs = dot(B_star[k, :], B_star[k, :])
        rhs = (delta - mu[k, k-1]^2) * dot(B_star[k-1, :], B_star[k-1, :])
        if lhs >= rhs
            k += 1
        else
            B[[k, k-1], :] = B[[k-1, k], :]
            U[[k, k-1], :] = U[[k-1, k], :]
            B_star, mu = gram_schmidt(B)
            k = max(k-1, 2)
        end
    end

    return B, U
end

# Returns shortest vector and its coordinates in basis B
function shortest_vector_approximation(B::AbstractMatrix{T}) where T
    reduced_basis, U = lll_reduction(B)

    infinity_norms = [maximum(abs, reduced_basis[i, :]) for i in 1:size(reduced_basis, 2)]
    shortest_index = argmin(infinity_norms)
    shortest_vector = reduced_basis[shortest_index, :]
    coordinates = U[shortest_index, :]

    return shortest_vector, coordinates
end

function reduced_row_echelon_form(A::AbstractMatrix{T}) where T
    A = copy(A)  # Work on a copy to avoid modifying the original matrix
    rows, cols = size(A)
    lead = 1
    row_count = 1

    for r in 1:rows
        if lead > cols
            break
        end

        # Find next pivot
        i = r
        while A[i, lead] == 0
            if i == rows
                i = r
                lead += 1
                if lead > cols
                    return A
                end
            else
                i += 1
            end
        end

        # Swap rows
        A[[r, i], :] = A[[i, r], :]

        # Normalize the pivot row
        lv = A[r, lead]
        A[r, :] = A[r, :] / lv

        for i in 1:rows
            if i != r
                A[i, :] -= A[i, lead] * A[r, :]
            end
        end

        lead += 1
    end

    return A
end

function extract_basis(A::AbstractMatrix{T}) where T
    B = vs[:, 1:1]
    for row in 2:size(vs, 1)
		v = vs[row, :]'
        if v * B' * inv(B * B') * B != v
            B = [B; v]
        end
    end
    return B
end

function row_space_basis(A::Matrix{T}) where T
    ref_A = reduced_row_echelon_form(A)
    basis = []

    for row in eachrow(ref_A)
        if any(row .!= 0)  # Check if the row is non-zero
            push!(basis, row)
        end
    end

    return hcat(basis...)'
end

function orthogonal_complement(A::AbstractMatrix{T}) where T
    m, n = size(A)
    A, _ = gram_schmidt(A)
    Q = Matrix{T}(I, n, n)  # Start with the identity matrix as the initial basis

    # Get orthogonal part of each vector
    for i in 1:n
        for j in 1:m
            Q[i, :] -= (dot(A[j, :], Q[i, :]) // dot(A[j, :], A[j, :])) * A[j, :]
        end
    end

    return row_space_basis(Q)
end

function hermite_normal_form(A::AbstractMatrix{T}) where T <: Integer
    H = copy(A)
    m, n = size(H)
    U = Matrix{T}(I, n, n)  # Unimodular transformation matrix

    pivot = 1 # row of first non-zero coefficient in column
    j = 1

    while j <= n && pivot <= m
        while j<n && !iszero(H[pivot, j+1:end])
            # Find the column with the smallest non-zero entry in the current row
            min_col = j
            for i in j+1:n
                if H[pivot, i] != 0 &&
                        (H[pivot, min_col] == 0 || H[pivot, i] < H[pivot, min_col])
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
                H[:, j] = -H[:, j]
                U[:, j] = -U[:, j]
            end

            # Reduce entries to the right of the pivot
            for i in j+1:n
                if H[pivot, i] != 0
                    q = div(H[pivot, i], H[pivot, j])
                    H[:, i] -= q * H[:, j]
                    U[:, i] -= q * U[:, j]
                end
            end
        end

        # Reduce entries to the left of the pivot
        if H[pivot, j] != 0
            for i in 1:j-1
                if H[pivot, i] != 0
                    q = div(H[pivot, i], H[pivot, j])
                    H[:, i] -= q * H[:, j]
                    U[:, i] -= q * U[:, j]
                end
            end
            j += 1
            pivot += 1
        else 
            pivot += 1
        end
    end

    return H, U
end

function lattice_basis_inside_span(B::AbstractMatrix{Rational{T}}) where T <: Integer
    B = row_space_basis(B)
    r = size(B, 1)
        
    if size(B, 1) == size(B, 2)
        return Matrix{T}(I, r, r)
    end
    B_perp= orthogonal_complement(B)
    B_perp = Matrix{T}(B_perp//gcd(B_perp))
    H, U = hermite_normal_form(B_perp)

    r = size(B, 1)
    n = size(U, 2)
    basis = U[:, n-r+1:n]'
    return basis
end

struct SignedCone{T}
    matrix::AbstractMatrix{Rational{T}}
    sign::Int
end

struct ConeData{T}
    matrix::AbstractMatrix{Rational{T}}
    sign::Int
    det::Rational{T}
end

function simplicial_cone_dual(A::AbstractMatrix{Rational{T}}) where T <: Integer
    if size(A, 1) != size(A, 2)
        error("matrix defining cone must be square")
    end

    A_det = det(A)

    rays = []
    for i in 1:size(A, 1)
        facet = A[setdiff(1:end, i), :]
        normal = orthogonal_complement(facet)
        ray_gen = normal//gcd(normal)
        
        if sign(det([A[1:i-1, :]; ray_gen; A[i+1:end, :]])) != sign(A_det)
            ray_gen = -ray_gen
        end

        push!(rays, ray_gen)
    end

    return reduce(vcat, rays)
end

# Assumes that C is simplicial
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

    cone_queue = Queue{ConeData{T}}()
    enqueue!(cone_queue, ConeData{T}(A_dual, 1, det(A_dual)))

    unimodular_cones = []

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
    
    primal_cones = unimodular_cones .|>
        dual_cone -> SignedCone(simplicial_cone_dual(dual_cone.matrix) * B,
                               dual_cone.sign)

    return primal_cones
end

# Assumes A is full dimensional
function subdivide_cone_once(cone::ConeData{T}) where T <:Integer
    if iszero(cone.det)
        error("cone is not full-dimensional")
    end

    A = cone.matrix
    A_inv = inv(A) 
    lambda, z = shortest_vector_approximation(A_inv)
    if all(lambda .< 0)
        lambda = -lambda
        z = -z
    end

    cones = []
    for i in 1:size(A, 1)
        C = copy(A)
        C[i, :] = z
        det_C = det(C)
        if det_C != 0
            sign_C = cone.sign * Int(sign(cone.det) * sign(det_C))
            push!(cones, ConeData(C, sign_C, det_C))
        end
    end
    return cones
end

function brion_coefficient(cone::SignedCone{T}, vertex, lambda) where T <: Integer
    A = cone.matrix
    num = BinomialSeries(dot(vertex, lambda))
    denom = FiniteSeries([1])
    for i in 1:size(A, 1)
        ray = A[i, :]
        denom *= 1 - BinomialSeries(dot(ray, lambda))
    end
    offset = first_nonzero_index(denom)
    quot = num / (denom << offset)
    return cone.sign * get_coeff(quot, offset)
end

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

function polytope_lattice_count(vertices)
    cones = []
    n = length(vertices)
    for i in 1:n
        r1 = vertices[(n+i-2)%n + 1] - vertices[i]
        r2 = vertices[i%n + 1] - vertices[i]
        r1 = r1//gcd(r1)
        r2 = r2//gcd(r2)
        push!(cones, [r1' ; r2'])
    end
    return compute_lattice_count(cones, vertices)
end

end

