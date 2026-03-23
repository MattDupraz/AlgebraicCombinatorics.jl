### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a3e5efb2-f8ef-45e7-be41-b1b05a6253bb
using Polymake

# ╔═╡ 5b2c43a4-dd97-11ef-1c79-1b5007c54b8d
using Symbolics, LinearAlgebra, Random

# ╔═╡ 719a18c6-4dd9-4633-9d5d-771bbd858ab4
using XML

# ╔═╡ 189b646c-1812-4b35-a433-718e2580346a
using Serialization

# ╔═╡ 76e48c1c-913a-4d39-ab5b-25654c3daa87
macro remember(str, expr)
	filename = "cache/$str.dat"
	try
		return deserialize(filename)
	catch
		try
			mkdir("cache")
		catch
		end
		return :(
			result = $expr;
			serialize($filename, result);
			result
		)
	end
end

# ╔═╡ bded0dea-feb6-4156-9c5e-69e95325089a
function get_latte_string(P)
	result = ""
	facets = Matrix(P.FACETS)
	linearity = Matrix(P.AFFINE_HULL)
	A = [facets; linearity]
	denoms = reduce(lcm, denominator.(A);
					init=BigInt(1),
					dims=2)
	A = A .* denoms
	d = size(A)
	result *= "$(d[1]) $(d[2])\n"
	for i in 1:size(A, 1)
		for j in 1:size(A, 2)
			result *= "$(A[i, j]) "
		end
		result *= "\n"
	end
	n_facets = size(facets, 1)
	n_lin = size(linearity, 1)
	if n_lin >= 1
		result *= "linearity $n_lin "
		for i in 1:n_lin
			result *= "$(n_facets + i) "
		end
	end
	return result
end

# ╔═╡ 3c9f89cf-c45e-4943-a375-99cb35c630ad
function count_lattice_points(P)
	str = get_latte_string(P)
	uid = rand(UInt32)

	dir = "temp-$uid"
	mkdir(dir)
	
	write("$dir/polytope", str)

	run(pipeline(Cmd(`latte-count polytope`, dir=dir); stdout=devnull, stderr=devnull))
	result = parse(BigInt, read("$dir/numOfLatticePoints", String))

	rm(dir; recursive=true, force=true)
	return result
end

# ╔═╡ 6ae0a0b0-0134-4cc3-a52d-fda551ace3c3
struct SVGVisual
	obj::Polymake.LibPolymake.BigObjectAllocated
end

# ╔═╡ 56fea39d-0a9d-4f70-8892-934a88264123
function Base.show(io::IO, ::MIME"image/svg+xml", x::Polymake.LibPolymake.BigObjectAllocated)
	xml_str = Polymake._get_visual_string_svg(Polymake.visual(Polymake.Visual, x))
	doc = parse(Node, xml_str)
	doc[3]["style"] = "background: #ffffff"
	doc[3]["width"] = "100%"
	doc[3]["height"] = ""
	print(io, XML.write(doc))
end

# ╔═╡ db24b03a-02dc-4d9b-834d-b4387b46c875
function Base.show(io::IO, ::MIME"image/svg+xml", svg::SVGVisual)
	show(io, MIME("image/svg+xml"), svg.obj)
end

# ╔═╡ 65b8f904-b394-464f-a8c7-cd6ac88ca6e0
function binomial(x::Number, k::Integer)
    k < 0 && throw(DomainError("k=$k must be non-negative"))
    k == 0 && return oneunit(x)//one(k)
    return prod(i -> (x-i+1)//i, 1:k)
end

# ╔═╡ cda53267-e690-4891-9a53-581e1ff683dd
struct Matroid
	n_elements::Int
	rank::Int
	flats::Vector{Set{Int}}
	ranks::Vector{Int}
	pm_matroid::Polymake.LibPolymake.BigObjectAllocated
	function Matroid(M::Polymake.LibPolymake.BigObjectAllocated)
		arr = Polymake.Array(M.LATTICE_OF_FLATS.DECORATION)
		flats = [Polymake.decoration_face(dec) for dec in arr]
		ranks = [Polymake.decoration_rank(dec) for dec in arr]
		new(M.N_ELEMENTS, M.RANK, flats, ranks, M)
	end
end

# ╔═╡ 002a96fe-1395-4200-9074-dcc591b39fc1
function svg_visual(M::Matroid)
	SVGVisual(M.pm_matroid.LATTICE_OF_FLATS)
end

# ╔═╡ 5174ba72-74e8-4ae2-b504-fa4a36f84b01
function groundset(M::Matroid)
	return 0:M.n_elements - 1
end

# ╔═╡ dfdeb32d-6c6e-4704-8006-ad7ab1233ed7
function get_flat_index(M, S)
	findfirst(map(F -> issetequal(S, F), M.flats))
end

# ╔═╡ eb4541cd-29fc-4f41-b668-97701ba9efbe
function meet(M, Fs)
	S = intersect([M.flats[i] for i in Fs]...)
	get_flat_index(M, S)
end

# ╔═╡ 2f0ed098-5d7e-4f3e-8ef0-97cb09d9dad4
function closure(M, S)
	get_flat_index(M, intersect(filter(F -> S ⊆ F, M.flats)...))
end

# ╔═╡ 4fb7b3e5-86d5-4b6a-94f8-b155ab1b6e2d
function join(M, Fs)
	S = union([M.flats[i] for i in Fs]...)
	closure(M, S)
end

# ╔═╡ 8176fa85-c220-4c41-a221-9f48eb090e74
function flats_of_rank_range(M, r, s)
	return [i for i in 1:length(M.flats)
		if r <= M.ranks[i] && M.ranks[i] <= s ]
end

# ╔═╡ ff29b293-0bb7-491c-85da-1112d61a5f32
function nonempty_flats(M)
    return flats_of_rank_range(M, 1, M.rank)
end

# ╔═╡ ceab87f9-4a40-4de4-9b16-49141064393a
function nonempty_proper_flats(M)
    return flats_of_rank_range(M, 1, M.rank - 1)
end

# ╔═╡ 4b85ad51-8cc6-4382-966e-dc6ce3172690
# Return the chains of flats from Fs
function chains(M, Fs; tail=[])
	cs = [tail]
    for F in Fs
        if isempty(tail) || M.flats[F] ⊊ M.flats[first(tail)]
            append!(cs, chains(M, Fs; tail=[F; tail]))
		end
	end
	return cs
end

# ╔═╡ 00681b2c-6ae0-48b5-8cfd-29d90384a7dd
# Returns the non-empty chains of proper nonempty flats of M
function chains(M)
	return filter(!isempty, chains(M, nonempty_proper_flats(M)))
end

# ╔═╡ 960fc95e-b0f4-4fce-8481-8ffec865bf33
function f_to_h(M)
	nep_flats = nonempty_proper_flats(M)
    n = length(nep_flats)
    f = Symbolics.variables(:f, nonempty_flats(M))
    h = Symbolics.variables(:h, nep_flats)

	A = Rational.(reduce(hcat, [-1 * [M.flats[F] ⊆ M.flats[G] for G in nep_flats]
                                for F in nep_flats
	])) # M is a change of basis from (L_F - H) to D_F
	return Dict(f .=> Rational.([Matrix(I, n, n); -ones(n)']) * inv(A) * h)
end

# ╔═╡ e27d30dc-5f52-468d-b482-2edf6aa67c88
function h_to_f(M)
	nep_flats = nonempty_proper_flats(M)
    f = Symbolics.variables(:f, nonempty_flats(M))
    h = Symbolics.variables(:h, nep_flats)
	
    # Hyperplane class
    H = [groundset(M)[1] ∈ M.flats[F] for F in nep_flats] 
	A = Rational.(reduce(hcat, [H - [M.flats[F] ⊆ M.flats[G] for G in nep_flats]
                                for F in nonempty_flats(M)
	])) # M is a change of basis from L_F to D_F
	return Dict(h .=> A * f)
end

# ╔═╡ af89df1d-7afd-47ba-b93d-035d33e1c75c
# returns all monomials in n_var variables of degree at most deg, represented as vectors of size n_var
function m_vectors_iter(max_m, max_sum)
	if length(max_m) == 1
		return 0:min(max_m[1], max_sum)
	end
	return Iterators.flatmap(
		n -> Iterators.product(n, m_vectors_iter(max_m[2:length(max_m)], max_sum - n))
				.|> t -> [t[1]; t[2]],
		0:min(max_m[1], max_sum))
end

# ╔═╡ 2bd8e1ac-471c-435c-af62-717c1d50d8e0
#=╠═╡
function snapper_polynomial(M)
    ne_flats = nonempty_flats(M)

    f = Symbolics.variables(:f, ne_flats)

    ms = m_vectors_iter(max.(0, M.ranks[ne_flats] .- 1), M.rank) |> collect
	snap = 0

    closures = Dict()

	for (i, m) in enumerate(ms)
		dragon_hall_rado = true
		write("progress", "$i/$(length(ms)) $m")
		for m_ in collect(vectors_less_than(m))
            if all(iszero, m_)
                continue
            end
            S = join(M, [F for (i, F) in enumerate(ne_flats)
                           if m_[i] > 0])

            if M.ranks[S] <= sum(m_)
				dragon_hall_rado = false
				break
			end
		end
		if dragon_hall_rado
            m_rat = Vector{Rational}(m)
            snap += prod(binomial.(m_rat + f .- 1//1, m))
		end
	end

	return Symbolics.expand(snap)
end
  ╠═╡ =#

# ╔═╡ 4a021ce4-d60a-4289-a45f-e7437384918f
function get_ray_gen(M, F)
    v = Rational[e in M.flats[F] for e in groundset(M)] * BigInt(1)
	v_proj = v[1:length(v) - 1] .- last(v)
	return v_proj
end

# ╔═╡ 8b252e90-fb1b-4dc6-85e0-81b2fedc61aa
function ray_gen(r)
	return r // gcd(r...)
end

# ╔═╡ 4f722a34-aaf7-4ef3-8d9e-46102be79e09
function homogeneize(A)
	return [ones(Rational, size(A, 1)) A]
end

# ╔═╡ dcbe98d9-25a6-4666-a7ff-57b28a6d1eb2
# the ws are the corners of the combinatorial hypercube
function get_w(r, h, k)
	J = [(k >> i) % 2 == 1 for i in 0:size(r, 1)-1]
	if iszero(J)
		return zero(r[1, :])
	end
	h_ = h[J]
	r_ = r[J, :]
	return r_' * inv(r_ * r_') * h_
end

# ╔═╡ d0adbc23-f2e2-46de-af36-55c2e5a4f276
function get_ws(r, h)
    return hcat(map(k -> get_w(r, h, k), 0:2^size(r, 1) - 1)...)'
end

# ╔═╡ eb00d49e-7107-4a43-a4fb-1c359138d252
function cone_lattice_count(B, r, h)
	w = get_ws(r, h)
    y = w * B
	P = polytope.Polytope(POINTS=homogeneize(y))
    return count_lattice_points(P)
end

# ╔═╡ 9fcad0b6-ac34-4286-b54f-f26b6a70291e
function get_polytope(B, r, h)
	w = get_ws(r, h)
    y = w * B
	P = polytope.Polytope(POINTS=homogeneize(y))
	return P
end

# ╔═╡ ee4558c0-af44-4aa5-a863-426ca84c0a8d
function compute_lattice_counts_parallel(B, r, hs)
	polytopes = hs .|> h -> get_polytope(B, r, h)
	latte_strs = polytopes .|> P -> get_latte_string(P)
	results = zeros(BigInt, size(polytopes))

	Threads.@threads for i in 1:length(polytopes)
		uid = rand(UInt32)
		dir = "temp-$uid"
		mkdir(dir)
		
		write("$dir/polytope", latte_strs[i])
		
		run(pipeline(Cmd(`latte-count polytope`, dir=dir); stdout=devnull, stderr=devnull))
		result = parse(BigInt, read("$dir/numOfLatticePoints", String))
		
		rm(dir; recursive=true, force=true)
		results[i] = result
	end
	return results
end

# ╔═╡ dfefa444-3fbc-42a4-96aa-c2af75f00425
function extract_basis(vs)
    B = vs[1:1, :]
    for row in 2:size(vs, 1)
		v = vs[row, :]'
        if v * B' * inv(B * B') * B != v
            B = [B; v]
        end
    end
    return B
end

# ╔═╡ 211351c0-cdc2-4635-b75c-cdfab56dec1b
# returns all monomials in n_var variables of degree at most deg, represented as vectors of size n_var
function monomials_iter(deg, n_vars)
	if n_vars == 1
		return 0:deg
	end
	return Iterators.flatmap(
		n -> Iterators.product(n, monomials_iter(deg - n, n_vars - 1))
				.|> t -> [t[1]; t[2]],
		0:deg)
end

# ╔═╡ 0863c1fb-eac8-4839-808f-1c0c69937e19
function interpolate_polynomial(xs, ys, vars, deg)
	monomials = monomials_iter(deg, length(vars))
	# M records the values of the monomials at the sample points xs
	M = (
		reduce(hcat, [[prod(x.^m) for m in monomials] for x in xs])
	)'
	# M should be invertible and we should have M * [a_1, ..., a_n] = ys, where the a_i are the coefficients of each monomial in the polynomial
	return dot(inv(M) * ys, [prod(vars .^ m) for m in monomials])
end

# ╔═╡ 36721cd7-608f-4ba6-9e4c-096f913cb882
function cone_ehr(B, r, vars)    
    cone = polytope.Cone(INPUT_RAYS=r)
    dual_cone = polytope.polarize(cone)
	# We want the outer corner to project positively on each ray,
	# so that the polyhedron is well-defined.
    pseudo_cubical = polytope.intersection(cone, dual_cone)
	pc_rays = Matrix{Rational}(pseudo_cubical.RAYS)
	# Rescale each ray generator by gcd to get primitive ray generators
	pc_rays = pc_rays .// reduce(gcd, pc_rays; init=0//1, dims=2)
	# Extract a basis, which will help us sample independent values
    pc_basis = extract_basis(pc_rays)

	# Dimension of the cone (should be the same as the number of rays in r)
	dim = size(pc_basis, 1)
	# Samples of outer corner in the interior of the cone
	samples = monomials_iter(dim, length(vars)) .|> m ->  pc_basis' * (m .+ 1)
	# Corresponding h values giving rise to these points
	h_values = samples .|> p -> r * p
	# Calculate all the corners for each h
	ws = h_values .|> h -> get_ws(r, h) * B
	# Want to assure all of these points are lattice points, so scale everything by common denominator
	denoms = ws .|> w -> reduce(lcm, denominator.(w); init =BigInt(1), dims=1)
	denom = lcm((Iterators.flatten(ws) .|> denominator)...)
	h_values_scaled = map(h -> denom * h, h_values)
	# Calculate lattice counts for selected samples
	#lattice_counts = h_values .|> h -> cone_lattice_count(B, r, denom*h)
	lattice_counts = compute_lattice_counts_parallel(B, r, h_values_scaled)
	# Interpolate the samples with a polynomial of expected degree
	f = interpolate_polynomial(h_values_scaled, lattice_counts, vars, dim)
	
    return f#, ws, denoms, denom, h_values_scaled, lattice_counts, (h_values_scaled .|> h -> get_ws(r, h) * B)
end

# ╔═╡ 845a0f6b-6628-4ee6-a89b-0b6ec31b6187
function compute_polynomials(M, B=I)	
	cones = chains(M)
	polynomials = Dict()
	
    for (i, cone) in enumerate(cones)
		r = reduce(hcat, [B * get_ray_gen(M, F) for F in cone])'
		vars = Symbolics.variables(:h, cone)
		for k in 1:2^length(cone)-1
			J = [(k >> i) % 2 == 1 for i in 0:length(cone)-1]
			subcone = cone[J]
			subcone_r = r[J, :]
			subcone_vars = vars[J]

			f = cone_ehr(r', subcone_r, subcone_vars)
			polynomials[(cone, subcone)] = f
		end
		
		write("progress", "$i/$(length(cones))")
    end
	
    return polynomials
end

# ╔═╡ 52873363-960a-4cd1-9e70-ac7055f928d7
x = Symbolics.variables(:x, 1:10)

# ╔═╡ 097e42ca-56ed-4a03-8ba8-0d8faafa849d
function ehrhart_polynomial(M, polynomials=nothing)
	if polynomials == nothing
		polynomials = compute_polynomials(M)
	end
	
    n = length(nonempty_proper_flats(M))

	cones = chains(M)
	dim = maximum(length.(cones))
	cones_by_dim = [filter(c -> length(c) == d, cones) for d in 1:dim]	
    ehr = 0
	if (dim >= 3)
	    for (i, cone) in enumerate(cones_by_dim[3])
	        r = reduce(hcat, [get_ray_gen(M, F) for F in cone])'
	        vars = Symbolics.variables(:h, cone)
	
			result = 0
			result += polynomials[(cone, cone)]
			for j in 1:3
				J = [k != j for k in 1:3]
				subcone = cone[J]
				result -= 1//2 * (polynomials[(cone, subcone)] + polynomials[(subcone, subcone)])
	
				for k in (1:3)[J]
					ray = cone[k:k]
					result += 1//6 * (polynomials[(cone, ray)] + polynomials[(subcone, ray)] + polynomials[(ray, ray)])
				end
	
				#result += x[4] * polynomials[(cone, cone[j:j])]
			end
			result -= 1
	
			ehr += result
	    end
	end
	
    for (i, cone) in enumerate(cones_by_dim[2])
		result = 0
		result += polynomials[(cone, cone)]
		result -= 1//2 * (polynomials[(cone, cone[1:1])] +
						polynomials[(cone[1:1], cone[1:1])])
		result -= 1//2 * (polynomials[(cone, cone[2:2])] +
						polynomials[(cone[2:2], cone[2:2])])
		result += 1

		ehr += result
    end

    for cone in cones_by_dim[1]
        ehr += polynomials[(cone, cone)] - 1
    end

    ehr += 1

    ehr = Symbolics.expand(ehr)
    return ehr
end

# ╔═╡ 5393fae8-ee75-4647-9205-56faeb53b176
# ╠═╡ disabled = true
#=╠═╡
fano = Matroid(matroid.fano_matroid())
  ╠═╡ =#

# ╔═╡ 078c01ce-886d-4e14-89fa-d201affd8ff0
#=╠═╡
chains(fano)
  ╠═╡ =#

# ╔═╡ 7f55c842-2e3f-4764-be9a-9c510626d362
#=╠═╡
fano_snap = snapper_polynomial(fano)
  ╠═╡ =#

# ╔═╡ a5206183-58a5-430e-89b9-1e7e87d623c5
#=╠═╡
fano_polynomials = compute_polynomials(fano)
  ╠═╡ =#

# ╔═╡ bb25f4cc-ff1f-4f67-990e-530e0a369b2f
# ╠═╡ disabled = true
#=╠═╡
Symbolics.expand(Symbolics.substitute(fano_snap, f_to_h(fano)) - ehrhart_polynomial(fano, fano_polynomials))
  ╠═╡ =#

# ╔═╡ be10a598-918c-4bec-af4d-017e9522e948
@time cone_ehr(r', r, Symbolics.variables(:x, 1:2))

# ╔═╡ 1c071329-005a-496e-ad78-f535286e7d6c
#=╠═╡
M.flats
  ╠═╡ =#

# ╔═╡ 39164e66-13d0-483c-b0f4-9cedb5d030d4
#=╠═╡
closure(M, Set([0, 4]))
  ╠═╡ =#

# ╔═╡ 686dd79f-4815-455c-8123-f60e7656357c
#=╠═╡
chains(M)
  ╠═╡ =#

# ╔═╡ 6fdd1ad0-3fa3-41eb-8e6c-a11f63db2fdc
#=╠═╡
svg_visual(M)
  ╠═╡ =#

# ╔═╡ 5fbbb173-fe2c-4538-aa4e-379396331c05
M_polynomials = @remember("M_polynomials", compute_polynomials(M));

# ╔═╡ d444f5c6-8686-4ab6-ac31-7e9f90ac9d35
#=╠═╡
M_ehrhart = ehrhart_polynomial(M, M_polynomials);
  ╠═╡ =#

# ╔═╡ f6254849-8c30-49e4-a736-0b10577b46f0
M_snap = @remember("M_snap", snapper_polynomial(M));

# ╔═╡ b403367b-93a8-4a73-bd13-48b372d23a52
M_ehrhart_true = @remember("M_ehrhart_true", 
	Symbolics.expand(Symbolics.substitute(M_snap, f_to_h(M))));

# ╔═╡ 5ac3cc96-7e85-45a9-972d-b7fb1bc4a3bf
#=╠═╡
M_ehrhart - M_ehrhart_true
  ╠═╡ =#

# ╔═╡ 82b0d9ec-ffe8-4491-b707-5379fee808c4
#=╠═╡
chains(M)
  ╠═╡ =#

# ╔═╡ 300cdd94-3f84-4240-bd1c-6fa21d252664
U = Matroid(matroid.uniform_matroid(4,5));

# ╔═╡ 475bd32e-2eba-4a44-9d21-c7ae214a3438
svg_visual(U)

# ╔═╡ ecc9e495-347f-479c-94a6-3578c63f0eae
U.flats

# ╔═╡ de72a522-9c20-4f27-a4b6-ea86cdef68cf
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

# ╔═╡ 1b20a227-d07a-4a3d-9a0a-703b3ced4177
U_polynomials = @remember("U_polynomials", compute_polynomials(U));

# ╔═╡ 2074cd20-8e28-481c-aaa0-24bdbb1bd41d
A = BigInt[1 0 0 0; 2 1 0 0; 1 0 1 0; 0 3 2 1]

# ╔═╡ 4b986231-8c14-40fe-9934-18089e11bc76
U_polynomials2 = @remember("U_polynomials2", compute_polynomials(U, A))

# ╔═╡ a242f773-e719-4401-afed-90669b4442f8
filter(p->Symbolics.coeff(p.second, 1) != 1, U_polynomials2)

# ╔═╡ bca104b9-852b-46ce-8ce8-dadd6d45d418
U_ehrhart = ehrhart_polynomial(U, U_polynomials2);

# ╔═╡ 18273213-efb2-4e64-9651-5983780e0e09
U_snap = @remember("U_snap", snapper_polynomial(U));

# ╔═╡ 6e4a6519-aa95-42dc-941e-41ef445c6c6a
U_ehrhart_true = @remember("U_ehrhart_true",
	Symbolics.expand(Symbolics.substitute(U_snap, f_to_h(U))));

# ╔═╡ c9c6d515-5e92-48fe-a990-15c6e85d4e28
f = U_ehrhart_true - U_ehrhart

# ╔═╡ 254166fa-2fa6-4fad-b3fb-57d8d9100eee
#=╠═╡
hs = Symbolics.variables(:h, nonempty_proper_flats(M));
  ╠═╡ =#

# ╔═╡ 3b58ab04-48a6-4da8-87f7-6df41b0c3e6c
#=╠═╡
monomials = Set([hs [h1 * h2 for (h1, h2) in Iterators.product(hs, hs)]]);
  ╠═╡ =#

# ╔═╡ 80967c37-7227-45e6-ae49-48cc5053c8a5
#=╠═╡
gs = [Symbolics.substitute(Symbolics.coeff(f, m), Dict(hs .=> 0)) for m in monomials]
  ╠═╡ =#

# ╔═╡ 1596e4d3-fd60-4ba4-852e-24486ba3df0a
#=╠═╡
gs[5] - gs[6]
  ╠═╡ =#

# ╔═╡ 86cc36e2-1766-4a63-9b8b-92ec2b589ebb
#=╠═╡
gs[15] - gs[16]
  ╠═╡ =#

# ╔═╡ 3bcdc4cc-5b49-4617-a362-adfeda9ae1c5
U34 = Matroid(matroid.uniform_matroid(3, 4))

# ╔═╡ b607abfe-f106-4849-8d67-4cdb01f864b7
#=╠═╡
cone_ehr(r', r, x[1:3])
  ╠═╡ =#

# ╔═╡ baa2348f-4be2-4acf-861b-325601dd62bf
U34_polyns = compute_polynomials(U34)

# ╔═╡ 48b4e6ce-a283-4f0b-92fe-4c0ded6571dd
ehrhart_polynomial(U34, U34_polyns)

# ╔═╡ 0781d38c-6eb2-48ca-afdf-bc520034eea6
#=╠═╡
U34_snapper = snapper_polynomial(U34)
  ╠═╡ =#

# ╔═╡ 93f2a35d-01da-411e-9563-53e40d83084d
#=╠═╡
U34_ehr_true = Symbolics.expand(Symbolics.substitute(U34_snapper, f_to_h(U34)))
  ╠═╡ =#

# ╔═╡ 26172237-30a0-46e8-a825-b9ab3ebd112f
#=╠═╡
M = Matroid(matroid.matroid_from_graph(graph.graph_from_edges([[0, 1], [1, 2], [0, 2], [1, 3], [2, 3], [3, 4], [1, 4]])))
  ╠═╡ =#

# ╔═╡ f25b25e5-42ee-4f03-9ba8-2f756a74569d
# ╠═╡ disabled = true
#=╠═╡
r = [get_ray_gen(U, 4) get_ray_gen(U, 8) get_ray_gen(U, 17)]' * A'
  ╠═╡ =#

# ╔═╡ e7806bb9-3fa0-4395-99ea-28f0b9f2b98c
# ╠═╡ disabled = true
#=╠═╡
function vectors_less_than(vec)
	return Iterators.product([0:n for n in vec]...)
end
  ╠═╡ =#

# ╔═╡ a2ddce7c-22ab-4b22-a46a-188c123bed82
# ╠═╡ disabled = true
#=╠═╡
r = Rational[1 0 1 1; 0 1 0 1]
  ╠═╡ =#

# ╔═╡ 4d2d10ad-55e0-4f69-a52f-fac7ad0c3b5d
# ╠═╡ disabled = true
#=╠═╡
M = Matroid(matroid.fano_matroid())
  ╠═╡ =#

# ╔═╡ 73a9d5f0-cade-432f-aeeb-6f530532342d
#=╠═╡
function vectors_less_than(vec)
	if length(vec) == 1
		return 0:vec[1]
	end
	return Iterators.flatmap(
		n -> Iterators.product(n, vectors_less_than(vec[2:length(vec)]))
				.|> t -> [t[1]; t[2]],
		0:vec[1])
end
  ╠═╡ =#

# ╔═╡ 3952d683-6446-445d-8b02-3922a9125d35
#=╠═╡
r = reduce(hcat, [get_ray_gen(U34, F) for F in nonempty_proper_flats(U34)])'
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Polymake = "d720cf60-89b5-51f5-aff5-213f193123e7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Serialization = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
XML = "72c71f33-b9b6-44de-8c94-c961784809e2"

[compat]
Polymake = "~0.11.24"
Symbolics = "~6.24.0"
XML = "~0.3.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "cabb40dc88d030228b294ba20945663b93e80e17"

[[deps.ADTypes]]
git-tree-sha1 = "ab3815fc06f9c171ee81277f00ecee19b47e0b91"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.12.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.ASL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6252039f98492252f9e47c312c8ffda0e3b9e78d"
uuid = "ae81ac8f-d209-56e5-92de-9978fef736f9"
version = "0.1.3+0"

[[deps.AbstractAlgebra]]
deps = ["LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "4d02060414d46e141e79efbd1484901c0fb4ce24"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.43.12"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "0ba8f4c1f06707985ffb4804fdad1bf97b233897"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.41"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "50c3c56a52972d78e8be9fd135bfb91c9574c140"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.1.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "017fcb757f8e921fb44ee063a7aafe5f89b86dd1"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.18.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "d8b0439d2be438a5f2cd68ec158fe08a7b2595b7"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.9"

[[deps.BinaryWrappers]]
deps = ["JLLWrappers", "Scratch"]
git-tree-sha1 = "7fea8f658689fa5062b23f4400eda888b7ae2aaa"
uuid = "f01c122e-0ea1-4f85-ad8f-907073ad7a9f"
version = "0.1.3"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8873e196c2eb87962a2048b3b8e08946535864a1"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+4"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.CxxWrap]]
deps = ["Libdl", "MacroTools", "libcxxwrap_julia_jll"]
git-tree-sha1 = "9208ab6d13777dace5981b451171891faf59e6a2"
uuid = "1f15a43c-97ca-5a2a-ae31-89f07a497df4"
version = "0.16.0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DecFP]]
deps = ["DecFP_jll", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "4a10cec664e26d9d63597daf9e62147e79d636e3"
uuid = "55939f99-70c6-5e9b-8bb0-5071ed7d61fd"
version = "1.3.2"

[[deps.DecFP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e9a8da19f847bbfed4076071f6fef8665a30d9e5"
uuid = "47200ebd-12ce-5be5-abb7-8e082af23329"
version = "2.0.3+1"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "03aa5d44647eaec98e1920635cdfed5d5560a8b9"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.117"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "a7e9f13f33652c533d49868a534bfb2050d1365f"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.15"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "9a3ae38b460449cc9e7dd0cfb059c76028724627"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.1"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Expronicon]]
deps = ["MLStyle", "Pkg", "TOML"]
git-tree-sha1 = "fc3951d4d398b5515f91d7fe5d45fc31dccb3c9b"
uuid = "6b7a57c9-7cc1-4fdf-b7f5-e857abae3636"
version = "0.8.5"

[[deps.FLINT_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "MPFR_jll", "OpenBLAS32_jll"]
git-tree-sha1 = "9327ebef3e04034b2aa1b9c59869f36717843ce5"
uuid = "e134572f-a0d5-539d-bddf-3cad8db41a82"
version = "300.100.301+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLPK_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6aa6294ba949ccfc380463bf50ff988b46de5bc7"
uuid = "e8aa6df9-e6ca-548a-97ff-1f85fc5b8b98"
version = "5.0.1+1"

[[deps.GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"
version = "6.3.0+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "50aedf345a709ab75872f80a2779568dc0bb461b"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.11.2+3"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "2bd56245074fab4015b9174f24ceba8293209053"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.27"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.Ipopt_jll]]
deps = ["ASL_jll", "Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "MUMPS_seq_jll", "SPRAL_jll", "libblastrampoline_jll"]
git-tree-sha1 = "4f55ad688c698a4f77d892a1cb673f7e8a30f178"
uuid = "9cc047cb-c261-5740-88fc-0cf96f7bdcc7"
version = "300.1400.1700+0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "191686b1ac1ea9c89fc52e996ad15d1d241d1e33"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.10.1+0"

[[deps.METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2eefa8baa858871ae7770c98c3c2a7e46daba5b4"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.3+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MPFR_jll]]
deps = ["Artifacts", "GMP_jll", "Libdl"]
uuid = "3a97d323-0669-5f0c-9066-3539efd106a3"
version = "4.2.1+0"

[[deps.MUMPS_seq_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "libblastrampoline_jll"]
git-tree-sha1 = "0eab12f94948ca67908aec14b9f2ebefd17463fe"
uuid = "d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"
version = "500.700.301+0"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MongoC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "snappy_jll"]
git-tree-sha1 = "0ab1e13479581def03ea0b19733b8b830177889a"
uuid = "90100e71-7732-535a-9be7-2e9affd1cfc1"
version = "1.28.1+0"

[[deps.Mongoc]]
deps = ["Dates", "DecFP", "MongoC_jll", "Serialization"]
git-tree-sha1 = "9b0ee6f5843a8df838d76794a5d50e77b8a3b08f"
uuid = "4fe8b98c-fc19-5c23-8ec2-168ff83495f2"
version = "0.9.2"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "8d39779e29f80aa6c071e7ac17101c6e31f075d7"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "43122df26d27424b23577d59e2d8020f28386516"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.2"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "fe891aea7ccd23897520db7f16931212454e277e"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.1"

[[deps.Ncurses_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3690e6c58c16ba676bcc9b5654762fe8a05db1c7"
uuid = "68e3532b-a499-55ff-9963-d1c0c0748b3a"
version = "6.5.0+1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ninja_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d9c6d1dca818768f1671244dcdf8e8e08b2b8937"
uuid = "76642167-d241-5cee-8c94-7a494e8cb7b7"
version = "1.12.1+0"

[[deps.OffsetArrays]]
git-tree-sha1 = "5e1897147d1ff8d98883cda2be2187dcf57d8f0c"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.15.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ece4587683695fe4c5f20e990da0ed7e83c351e7"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+3"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "12f1439c4f986bb868acda6ea33ebc78e19b95ad"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.7.0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "966b85253e959ea89c53a9abebbf2e964fbf593b"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.32"

[[deps.PPL_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "dc69ba2b374b5bd6e876fa3a8441563c0eeafb9e"
uuid = "80dd9cbb-8b87-5171-a280-372cc418f402"
version = "1.2.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Perl_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Readline_jll"]
git-tree-sha1 = "7ab65a258bcf6da373cab49af462aead452d3960"
uuid = "83958c19-0796-5285-893e-a1267f8ec499"
version = "5.34.1+0"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.Polymake]]
deps = ["AbstractAlgebra", "BinaryWrappers", "CxxWrap", "Downloads", "JSON", "Libdl", "Mongoc", "NetworkOptions", "Ninja_jll", "Perl_jll", "Pidfile", "Pkg", "REPL", "Scratch", "SparseArrays", "TOPCOM_jll", "lib4ti2_jll", "libpolymake_julia_jll", "polymake_jll", "polymake_oscarnumber_jll"]
git-tree-sha1 = "0a606d37843442f11b972046610eb94216b6089c"
uuid = "d720cf60-89b5-51f5-aff5-213f193123e7"
version = "0.11.24"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "b8a399e95663485820000f26b6a43c794e166a49"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.4"

[[deps.Readline_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ncurses_jll"]
git-tree-sha1 = "69684dc9c2c69f7c515097841991362cca0739ea"
uuid = "05236dd9-4125-5232-aa7c-9ec0c9b2c25a"
version = "8.2.1+1"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "32f824db4e5bab64e25a12b22483a30a6b813d08"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.27.4"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SCIP_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "GMP_jll", "Ipopt_jll", "JLLWrappers", "Libdl", "Pkg", "Readline_jll", "Zlib_jll", "bliss_jll", "boost_jll"]
git-tree-sha1 = "cf69186eb29307fbb2319b90e6133797bad983ce"
uuid = "e5ac4fe4-a920-5659-9bf8-f9f73e9e79ce"
version = "800.0.301+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SPRAL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "Libdl", "METIS_jll", "libblastrampoline_jll"]
git-tree-sha1 = "11f3da4b25efacd1cec8e263421f2a9003a5e8e0"
uuid = "319450e9-13b8-58e8-aa9f-8fd1420848ab"
version = "2024.5.8+0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "Expronicon", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "80396dc3496a1718fc356797fa54034b6c533a18"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.72.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "6149620767866d4b0f0f7028639b6e661b6a1e44"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.12"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "0444a37a25fab98adbd90baa806ee492a3af133a"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.6.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "47091a0340a675c738b1304b58161f3b0839d454"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.10"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "fd2d4f0499f6bb4a0d9f5030f5c7d61eed385e03"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.37"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fabf4650afe966a2ba646cabd924c3fd43577fc3"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TermInterface", "TimerOutputs", "Unityper", "WeakValueDicts"]
git-tree-sha1 = "f631a4de315d3486ec768f9e9a6a3a85660427e8"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.11.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "IfElse", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "OffsetArrays", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "1211b5df046c13c955ba90c210c6394b99b2d1eb"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.24.0"

    [deps.Symbolics.extensions]
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TOPCOM_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg", "cddlib_jll"]
git-tree-sha1 = "adbe178144e762ae7057fcb8e26564de7ee2e36c"
uuid = "36f60fef-b880-50dc-9289-4aaecee93cc3"
version = "0.17.8+0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "d7298ebdfa1654583468a487e8e83fae9d72dac3"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.26"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.WeakValueDicts]]
git-tree-sha1 = "98528c2610a5479f091d470967a25becfd83edd0"
uuid = "897b6980-f191-5a31-bcb0-bf3c4585e0c1"
version = "0.1.0"

[[deps.XML]]
deps = ["Mmap", "OrderedCollections"]
git-tree-sha1 = "360c178fdee6816bcc9d436c6f324d582e84f1d1"
uuid = "72c71f33-b9b6-44de-8c94-c961784809e2"
version = "0.3.4"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "622cf78670d067c738667aaa96c553430b65e269"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+0"

[[deps.bliss_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f8b75e896a326a162a4f6e998990521d8302c810"
uuid = "508c9074-7a14-5c94-9582-3d4bc1871065"
version = "0.77.0+1"

[[deps.boost_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7a89efe0137720ca82f99e8daa526d23120d0d37"
uuid = "28df3c45-c428-5900-9ff8-a3135698ca75"
version = "1.76.0+1"

[[deps.cddlib_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "08f5df03703af917b9bfec47b9767eb943220d08"
uuid = "f07e07eb-5685-515a-97c8-3014f6152feb"
version = "0.94.14+0"

[[deps.lib4ti2_jll]]
deps = ["Artifacts", "GLPK_jll", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79116185def638e73ea3d88ca6c10e210a1dc183"
uuid = "1493ae25-0f90-5c0e-a06c-8c5077d6d66f"
version = "1.6.10+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libcxxwrap_julia_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b386562c39065094e2d1e1939ea289dbd448fe3e"
uuid = "3eaa8342-bff7-56a5-9981-c04077f7cee7"
version = "0.13.3+0"

[[deps.libpolymake_julia_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "FLINT_jll", "JLLWrappers", "Libdl", "TOPCOM_jll", "lib4ti2_jll", "libcxxwrap_julia_jll", "polymake_jll"]
git-tree-sha1 = "4885c6d17996ea3fa5492cfe8c68f0f52be82820"
uuid = "4d8266f6-2b3b-57e3-ad7a-d431eaaac945"
version = "0.13.0+0"

[[deps.lrslib_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9169c4d5823195a7a2173fc6b479468478857438"
uuid = "3873f7d0-7b7c-52c3-bdf4-8ab39b8f337a"
version = "0.3.3+0"

[[deps.nauty_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "fc7dc197df0648cd5f965801bfe086abd9325add"
uuid = "55c6dc9b-343a-50ca-8ff2-b71adb3733d5"
version = "2.6.13+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.normaliz_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "FLINT_jll", "GMP_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl", "MPFR_jll", "nauty_jll"]
git-tree-sha1 = "b84641c3add3c9e435546b1f2867105ce209791f"
uuid = "6690c6e9-4e12-53b8-b8fd-4bffaef8839f"
version = "300.1000.200+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.polymake_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "FLINT_jll", "GMP_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl", "MPFR_jll", "MongoC_jll", "PPL_jll", "Perl_jll", "SCIP_jll", "bliss_jll", "boost_jll", "cddlib_jll", "lrslib_jll", "normaliz_jll"]
git-tree-sha1 = "bfd75c077ece1a66325ec61a6cb6b0b23420f4b2"
uuid = "7c209550-9012-526c-9264-55ba7a78ba2c"
version = "400.1300.2+0"

[[deps.polymake_oscarnumber_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl", "libcxxwrap_julia_jll", "libpolymake_julia_jll", "polymake_jll"]
git-tree-sha1 = "0f1f590da4a2323defeecbf50a1d21448e7108fc"
uuid = "10f31823-b687-53e6-9f29-edb9d4da9f9f"
version = "0.3.3+1"

[[deps.snappy_jll]]
deps = ["Artifacts", "JLLWrappers", "LZO_jll", "Libdl", "Lz4_jll", "Zlib_jll"]
git-tree-sha1 = "a6c92533edec9755b3c1083b028da19a41f92e25"
uuid = "fe1e1685-f7be-5f59-ac9f-4ca204017dfd"
version = "1.2.1+3"
"""

# ╔═╡ Cell order:
# ╠═a3e5efb2-f8ef-45e7-be41-b1b05a6253bb
# ╠═5b2c43a4-dd97-11ef-1c79-1b5007c54b8d
# ╠═719a18c6-4dd9-4633-9d5d-771bbd858ab4
# ╠═76e48c1c-913a-4d39-ab5b-25654c3daa87
# ╠═bded0dea-feb6-4156-9c5e-69e95325089a
# ╠═3c9f89cf-c45e-4943-a375-99cb35c630ad
# ╠═6ae0a0b0-0134-4cc3-a52d-fda551ace3c3
# ╠═56fea39d-0a9d-4f70-8892-934a88264123
# ╠═db24b03a-02dc-4d9b-834d-b4387b46c875
# ╠═65b8f904-b394-464f-a8c7-cd6ac88ca6e0
# ╠═cda53267-e690-4891-9a53-581e1ff683dd
# ╠═002a96fe-1395-4200-9074-dcc591b39fc1
# ╠═5174ba72-74e8-4ae2-b504-fa4a36f84b01
# ╠═4d2d10ad-55e0-4f69-a52f-fac7ad0c3b5d
# ╠═dfdeb32d-6c6e-4704-8006-ad7ab1233ed7
# ╠═eb4541cd-29fc-4f41-b668-97701ba9efbe
# ╠═2f0ed098-5d7e-4f3e-8ef0-97cb09d9dad4
# ╠═4fb7b3e5-86d5-4b6a-94f8-b155ab1b6e2d
# ╠═8176fa85-c220-4c41-a221-9f48eb090e74
# ╠═ff29b293-0bb7-491c-85da-1112d61a5f32
# ╠═ceab87f9-4a40-4de4-9b16-49141064393a
# ╠═1c071329-005a-496e-ad78-f535286e7d6c
# ╠═39164e66-13d0-483c-b0f4-9cedb5d030d4
# ╠═4b85ad51-8cc6-4382-966e-dc6ce3172690
# ╠═00681b2c-6ae0-48b5-8cfd-29d90384a7dd
# ╠═686dd79f-4815-455c-8123-f60e7656357c
# ╠═960fc95e-b0f4-4fce-8481-8ffec865bf33
# ╠═e27d30dc-5f52-468d-b482-2edf6aa67c88
# ╠═af89df1d-7afd-47ba-b93d-035d33e1c75c
# ╠═73a9d5f0-cade-432f-aeeb-6f530532342d
# ╠═e7806bb9-3fa0-4395-99ea-28f0b9f2b98c
# ╠═2bd8e1ac-471c-435c-af62-717c1d50d8e0
# ╠═4a021ce4-d60a-4289-a45f-e7437384918f
# ╠═8b252e90-fb1b-4dc6-85e0-81b2fedc61aa
# ╠═4f722a34-aaf7-4ef3-8d9e-46102be79e09
# ╠═dcbe98d9-25a6-4666-a7ff-57b28a6d1eb2
# ╠═d0adbc23-f2e2-46de-af36-55c2e5a4f276
# ╠═eb00d49e-7107-4a43-a4fb-1c359138d252
# ╠═9fcad0b6-ac34-4286-b54f-f26b6a70291e
# ╠═ee4558c0-af44-4aa5-a863-426ca84c0a8d
# ╠═dfefa444-3fbc-42a4-96aa-c2af75f00425
# ╠═211351c0-cdc2-4635-b75c-cdfab56dec1b
# ╠═0863c1fb-eac8-4839-808f-1c0c69937e19
# ╠═36721cd7-608f-4ba6-9e4c-096f913cb882
# ╠═845a0f6b-6628-4ee6-a89b-0b6ec31b6187
# ╠═52873363-960a-4cd1-9e70-ac7055f928d7
# ╠═097e42ca-56ed-4a03-8ba8-0d8faafa849d
# ╠═5393fae8-ee75-4647-9205-56faeb53b176
# ╠═078c01ce-886d-4e14-89fa-d201affd8ff0
# ╠═7f55c842-2e3f-4764-be9a-9c510626d362
# ╠═a5206183-58a5-430e-89b9-1e7e87d623c5
# ╠═bb25f4cc-ff1f-4f67-990e-530e0a369b2f
# ╠═a2ddce7c-22ab-4b22-a46a-188c123bed82
# ╠═be10a598-918c-4bec-af4d-017e9522e948
# ╠═26172237-30a0-46e8-a825-b9ab3ebd112f
# ╠═6fdd1ad0-3fa3-41eb-8e6c-a11f63db2fdc
# ╠═189b646c-1812-4b35-a433-718e2580346a
# ╠═5fbbb173-fe2c-4538-aa4e-379396331c05
# ╠═d444f5c6-8686-4ab6-ac31-7e9f90ac9d35
# ╠═f6254849-8c30-49e4-a736-0b10577b46f0
# ╠═b403367b-93a8-4a73-bd13-48b372d23a52
# ╠═5ac3cc96-7e85-45a9-972d-b7fb1bc4a3bf
# ╠═82b0d9ec-ffe8-4491-b707-5379fee808c4
# ╠═300cdd94-3f84-4240-bd1c-6fa21d252664
# ╠═475bd32e-2eba-4a44-9d21-c7ae214a3438
# ╠═ecc9e495-347f-479c-94a6-3578c63f0eae
# ╠═de72a522-9c20-4f27-a4b6-ea86cdef68cf
# ╠═1b20a227-d07a-4a3d-9a0a-703b3ced4177
# ╠═2074cd20-8e28-481c-aaa0-24bdbb1bd41d
# ╠═f25b25e5-42ee-4f03-9ba8-2f756a74569d
# ╠═b607abfe-f106-4849-8d67-4cdb01f864b7
# ╠═4b986231-8c14-40fe-9934-18089e11bc76
# ╠═a242f773-e719-4401-afed-90669b4442f8
# ╠═bca104b9-852b-46ce-8ce8-dadd6d45d418
# ╠═18273213-efb2-4e64-9651-5983780e0e09
# ╠═6e4a6519-aa95-42dc-941e-41ef445c6c6a
# ╠═c9c6d515-5e92-48fe-a990-15c6e85d4e28
# ╠═254166fa-2fa6-4fad-b3fb-57d8d9100eee
# ╠═3b58ab04-48a6-4da8-87f7-6df41b0c3e6c
# ╠═80967c37-7227-45e6-ae49-48cc5053c8a5
# ╠═1596e4d3-fd60-4ba4-852e-24486ba3df0a
# ╠═86cc36e2-1766-4a63-9b8b-92ec2b589ebb
# ╠═3bcdc4cc-5b49-4617-a362-adfeda9ae1c5
# ╠═3952d683-6446-445d-8b02-3922a9125d35
# ╠═baa2348f-4be2-4acf-861b-325601dd62bf
# ╠═48b4e6ce-a283-4f0b-92fe-4c0ded6571dd
# ╠═0781d38c-6eb2-48ca-afdf-bc520034eea6
# ╠═93f2a35d-01da-411e-9563-53e40d83084d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
