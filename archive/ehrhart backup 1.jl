ENV["CRYPTOGRAPHY_OPENSSL_NO_LEGACY"] = 1
using Symbolics, LinearAlgebra, Random, Oscar, PyCall
_arango = pyimport("arango")

function binomial(x::Number, k::Integer)
    k < 0 && throw(DomainError("k=$k must be non-negative"))
    k == 0 && return oneunit(x)/one(k)
    return prod(i -> (x-i+1)/i, 1:k)
end

function nonempty_flats(M)
    return filter(!isempty, flats(M))
end

function nonempty_proper_flats(M)
    return filter(F -> F ≠ matroid_groundset(M), nonempty_flats(M))
end

function chains(sets; tail=[])
	cs = [tail]
    for S in sets
        if isempty(tail) || S ⊊ first(tail)
            append!(cs, chains(sets; tail=[S, tail...]))
		end
	end
	return cs
end

function nonempty_chains(sets)
    return filter(!isempty, chains(sets))
end

function vectors_less_than(vec)
    return Iterators.product([0:n for n in vec]...)
end

function join(M, Fs)
    return closure(M, union(Fs...))
end

function snapper_polynomial(M)
    ne_flats = nonempty_flats(M)

    f = Symbolics.variables(:f, 1:length(ne_flats))

    ranks_vec = [rank(M, F) for F in ne_flats]
    ranks = Dict(nonempty_flats(M) .=> ranks_vec)
    ms = vectors_less_than(max.(0, ranks_vec .- 1))
	snap = 0

    closures = Dict()

	for m in ms
		dragon_hall_rado = true
		for m_ in vectors_less_than(m)
            if all(iszero, m_)
                continue
            end
            S = union([F for (i, F) in enumerate(ne_flats)
                           if m_[i] > 0]...)
            
            F = 0
            if haskey(closures, S)
                F = closures[S]
            else
                F = closure(M, S)
                closures[S] = F
            end

            if ranks[F] <= sum(m_)
				dragon_hall_rado = false
				break
			end
		end
		if dragon_hall_rado
            m = collect(m)
            snap += prod(binomial.(m + f .- 1, m))
		end
	end

	return Symbolics.expand(snap)
end

struct Line
	origin
	ray_gen
	Line(origin::Vector, ray_gen::Vector) = new(origin, ray_gen)
	Line(ray_gen::Vector) = new(zeros(size(ray_gen)), ray_gen)
end

function Base.:*(A::Matrix, l::Line)::Line
	return Line(A*l.origin, A*l.ray_gen)
end

function Base.:*(l::Line, A::Matrix)::Line
	return Line(l.origin*A, l.ray_gen*A)
end

function get_ray_gen(M, F)
    v = [e in F for e in matroid_groundset(M)] * 1
	v_proj = v[1:length(v) - 1] .- last(v)
	return v_proj
end

function f_to_h(M)
    n = length(nonempty_proper_flats(M))
    f = Symbolics.variables(:f, 1:n + 1)
    h = Symbolics.variables(:h, 1:n)

	A = Rational.(reduce(hcat, [-1 * [F ⊆ G for G in nonempty_proper_flats(M)]
                                for F in nonempty_proper_flats(M)
	])) # M is a change of basis from (L_F - H) to D_F
    A * (f[1:n] .- last(f))
	return Dict(f .=> Rational.([Matrix(I, n, n); -ones(n)']) * inv(A) * h)
end

function h_to_f(M)
    n = length(nonempty_proper_flats(M))
    f = Symbolics.variables(:f, 1:n + 1)
    h = Symbolics.variables(:h, 1:n)
	
    # Hyperplane class
    H = [matroid_groundset(M)[1] ∈ F for F in nonempty_proper_flats(M)] 
	A = Rational.(reduce(hcat, [H - [F ⊆ G for G in nonempty_proper_flats(M)]
                                for F in nonempty_flats(M)
	])) # M is a change of basis from L_F to D_F
	return Dict(h .=> A * f)
end

function parallel(r1, r2)
	return det([r1 r2]) == 0
end

function intersect(l1::Line, l2::Line)
	if parallel(l1.ray_gen, l2.ray_gen)
		return nothing
	end
	A = [-l1.ray_gen l2.ray_gen]
	result = l1.origin + l1.ray_gen * (inv(A) * (l1.origin - l2.origin))[1]
	return Symbolics.expand.(result)
end

function ray_gen(r)
	return r // gcd(r...)
end

function divide(v, r)
	if r[1] != 0
		return v[1]//r[1]
	elseif r[2] != 0
		return v[2]//r[2]
	else
		throw(DomainError("r is zero"))
	end
end

function parallelogram_ehr(r1, r2, d1, d2)
	return Symbolics.expand(det([d1*r1 d2*r2]) + d1 + d2 + 1)
end

function quadrilateral_ehr(r1, r2, r3, r4, d1, d2, d3, d4)
	return Symbolics.expand(
		(parallelogram_ehr(r1, r2, d1, d2)
		+ parallelogram_ehr(r3, r4, d3, d4))//2
	)
end

function quadrilateral_vol(P1, P2, P3, P4)
    vol = (det([P2 - P1 P3 - P1]) + det([P3 - P4 P2 - P4]))//2;
	return Symbolics.expand(vol)
end

function cone2d_vol(r1, r2, h1, h2)
	cone_basis = Rational[r1 r2]

	r1_perp = r2 - dot(r1, r2)//dot(r1, r1) * r1
	r2_perp = -(r1 - dot(r1, r2)//dot(r2, r2) * r2)

	s1 = cone_basis' * r1
	s2 = cone_basis' * r2
	s1_perp = cone_basis' * r1_perp
	s2_perp = cone_basis' * r2_perp
	
	P1 = h1 * s1 // dot(r1, r1)
	P2 = h2 * s2 // dot(r2, r2)

	L1 = Line(P1, s1_perp)
	L2 = Line(P2, s2_perp)
	
	Q = intersect(L1, L2)

    vol = quadrilateral_vol(zeros(Rational, size(Q)), P1, P2, Q)
    return Symbolics.expand(vol)
end

function cone2d_ehr(r1, r2, h1, h2; exclude_r1=false, exclude_r2=false)
	cone_basis = Rational[r1 r2]

	r1_perp = r2 - dot(r1, r2)//dot(r1, r1) * r1
	r2_perp = -(r1 - dot(r1, r2)//dot(r2, r2) * r2)

	s1 = cone_basis' * r1
	s2 = cone_basis' * r2
	s1_perp = cone_basis' * r1_perp
	s2_perp = cone_basis' * r2_perp

	g1 = ray_gen(s1)
	g2 = ray_gen(s2)
	g1_perp = ray_gen(s1_perp)
	g2_perp = ray_gen(s2_perp)
	
	P1 = h1 * s1 // dot(r1, r1)
	P2 = h2 * s2 // dot(r2, r2)

	L1 = Line(P1, g1_perp)
	L2 = Line(P2, g2_perp)
	
	Q = intersect(L1, L2)

	d1 = divide(P1, g1)
	d2 = divide(P2, g2)
	
	d1_perp = divide(P1 - Q, -g1_perp)
	d2_perp = divide(P2 - Q, g2_perp)

	correction = 0
	if exclude_r1 && exclude_r2
		correction = d1 + d2 + 1
	elseif exclude_r1
		correction = d1 + 1
	elseif exclude_r2
		correction = d2 + 1
	end
	
	ehr = quadrilateral_ehr(g1, g2, g2_perp, -g1_perp, d1, d2, d2_perp, d1_perp) - correction

    #vol = (h1 * h2 //(dot(r1, r1) * dot(r2, r2))) *
    #        det(cone_basis' * cone_basis)//2 + 
    #        (h1 - h2 * dot(r1, r2)//dot(r2, r2)) * 
    #        (h2 - h1 * dot(r1, r2)//dot(r1, r1))//2

    #println(Symbolics.expand(vol - quadrilateral_vol(zeros(Rational, size(Q)), P1, P2, Q)))

    #ehr2 = 1//2*((h1 - h2 * dot(r1, r2)//dot(r2, r2)) +
    #        (h2 - h1 * dot(r1, r2)//dot(r1, r1)) +
    #        h1//dot(r1, r1)*gcd(dot(r1, r1), dot(r1, r2)) +
    #        h2//dot(r2, r2)*gcd(dot(r2, r2), dot(r1, r2))) + vol + 1

    #println(Symbolics.expand(ehr + correction - ehr2))
    
    return Symbolics.expand(ehr)
end

# Returns the number of lattice points in the polytope delimited by the cone
# over the vectors in r and cut out by the normal hyperplanes
# B is a change of basis (for the lattice)
function cone_lattice_count(B, r, h)
    function get_w(k)
        J = [i+1 for i in 0:length(r)-1 if (k >> i) % 2 == 1]
        if isempty(J)
            return r[1] * 0
        end
        h_ = [h[j] for j in J]
        R_ = reduce(hcat, [r[j] for j in J])
        return R_ * inv(R_' * R_) * h_
    end
    w = map(get_w, 0:2^length(r) - 1)
    y = [B * u for u in w]
    return length(lattice_points(convex_hull(y)))
end

function conv2rat(A::Array{QQFieldElem})
    return map(x -> x.num//x.den, A)
end

function extract_basis(vs)
    B = vs[1]
    for v in vs[2:length(vs)]
        if B * inv(B' * B) * B' * v != v
            B = hcat(B, v)
        end
    end
    return B
end

function cone_ehr(B, r, h)
	#monomials = vec(Iterators.product(
	#	[0:deg for var in h]...
	#) |> collect |> filter(m -> sum(m) <= deg))
    
    B = reduce(hcat, r)'
    R = reduce(hcat, [B * u for u in r])
    
    cone = positive_hull(R)
    dual_cone = polarize(cone)
    pseudo_cubical = Oscar.intersect(cone, dual_cone)
    pc_rays = map(conv2rat ∘ collect, rays(pseudo_cubical))
    pc_basis = extract_basis(pc_rays)

    denoms = [lcm(map(denominator, B * pc_basis)[i, :]...)
              for i in 1:size(pc_basis,1)]


    return pc_rays, pc_basis, B*(pc_basis * denoms)
end

function interpolate_polynomial(f, vars, deg, xs_basis)
	monomials = vec(Iterators.product(
		[0:deg for var in vars]...
	) |> collect |> filter(m -> sum(m) <= deg))
	xs = [sum(xs_basis .* m) for m in monomials] 
	ys = f.(xs)
	M = Symbolics.value.(
		reduce(hcat, [[prod(x.^m) for m in monomials] for x in xs])
	)'
	return dot(inv(M) * ys, [prod(vars .^ m) for m in monomials])
end

function ray_ehr(r1, r2, h1, h2)
	cone_basis = Rational[r1 r2]

	s1 = cone_basis' * r1
	g1 = ray_gen(s1)
	P1 = h1 * s1 // dot(r1, r1)
	d1 = divide(P1, g1)

	return d1 + 1
end

function random_GL(n)
	A = Matrix{Rational}(I, n, n)
	k = rand(1:100)
	for i in 1:k
		p = randperm(n)
		E = Matrix{Rational}(I, n, n)
		if rand(Bool)
			Base.permutecols!!(E, p)
		else
			E[p[1], p[2]] += -1 * rand(Bool)
		end
		A = A*E
	end
	return A
end

function ehrhart_polynomial(M)
    n = length(nonempty_proper_flats(M))
    h = Dict(nonempty_proper_flats(M) .=> Symbolics.variables(:h, 1:n))

    A = random_GL(length(matroid_groundset(M)) - 1)
    rays = Dict(F => A * get_ray_gen(M, F) for F in nonempty_proper_flats(M))

    cones_2d = filter(c -> length(c) == 2,
                      nonempty_chains(nonempty_proper_flats(M)))
    ehr = 0
    for cone in cones_2d
        r1 = rays[cone[1]]
        r2 = rays[cone[2]]
        h1 = h[cone[1]]
        h2 = h[cone[2]]
        ehr += cone2d_ehr(r1, r2, h1, h2; exclude_r1=true, exclude_r2=true)
        ehr += 1//2 * (ray_ehr(r1, r2, h1, h2) - (1 + h1))
        ehr += 1//2 * (ray_ehr(r2, r1, h2, h1) - (1 + h2))
    end

    for F in nonempty_proper_flats(M)
        ehr += h[F]
    end

    ehr += 1

    ehr = Symbolics.expand(ehr)
    return ehr
end

function volume_polynomial(M)
    n = length(nonempty_proper_flats(M))
    h = Dict(nonempty_proper_flats(M) .=> Symbolics.variables(:h, 1:n))

    A = random_GL(length(matroid_groundset(M)) - 1)
    rays = Dict(F => A * get_ray_gen(M, F) for F in nonempty_proper_flats(M))

    cones_2d = filter(c -> length(c) == 2,
                      nonempty_chains(nonempty_proper_flats(M)))
    vol = 0
    for cone in cones_2d
        r1 = rays[cone[1]]
        r2 = rays[cone[2]]
        h1 = h[cone[1]]
        h2 = h[cone[2]]
        vol += 2*cone2d_vol(r1, r2, h1, h2)
    end

    vol = Symbolics.expand(vol)
    return vol
end

function vol_to_ehr(M)
    n = length(nonempty_proper_flats(M))
    h = Symbolics.variables(:h, 1:n)

    vol = volume_polynomial(M)//2
    result = 1 + vol
    for var in h
        result += Symbolics.derivative(vol, var)//2
    end
    return Symbolics.expand(result)
end

function compare_polynomials(M)
    snap = snapper_polynomial(M)
    ehr = ehrhart_polynomial(M)
    return Symbolics.expand(substitute(ehr, h_to_f(M)) - snap)
end

function to_h_dict(M, line_bundle)
    n = length(nonempty_flats(M))

    f = Symbolics.variables(:f, 1:n)
    L = Symbolics.variables(:L, 1:n)
    h = Symbolics.variables(:h, 1:n - 1)
    D = Symbolics.variables(:D, 1:n - 1)

    vals = [Symbolics.coeff(line_bundle, D[i]) for i in 1:n-1]

    dict = h_to_f(M)
    for j in 1:n
        f_to_coeff = Dict(f .=> Num(0))
        f_to_coeff[f[j]] = Symbolics.coeff(line_bundle, L[j])
        vals += [substitute(dict[h[i]], f_to_coeff)
                 for i in 1:n-1]
    end

    return Dict(h .=> vals)
end

function eval_h(M, p, line_bundle)
    return Symbolics.expand(substitute(p, to_h_dict(M, line_bundle)))
end

function get_matroids_db()
    client = _arango.ArangoClient(hosts="https://matroid.mathematik.uni-siegen.de:443")
    db = client.db("MatroidsDB", username="matroid", password="matroid")
    matroids = db.collection("matroids_split_public")
    return db
end

function query_db(db; representable=false, min_atoms=1, max_atoms=11, limit=100)
    db.aql.execute("FOR doc IN matroids_split_public FILTER doc.NumberOfAtoms >= $min_atoms && doc.NumberOfAtoms <= $max_atoms && doc.IsRepresentable == $representable LIMIT $limit RETURN doc")
end

function matroid_from_doc(doc)
    return matroid_from_hyperplanes(doc["ListOfCoatoms"], doc["NumberOfAtoms"])
end

