### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ 312dfd9d-e4d0-4fac-a2eb-0191020ff66a
begin
	import Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	using AlgebraicCombinatorics
	using Polymake: matroid, graph
	using Combinatorics: powerset
	using LinearAlgebra: I
	import Symbolics as Sym
end

# ╔═╡ 845a0f6b-6628-4ee6-a89b-0b6ec31b6187
function compute_polynomials(M; B=I)	
	cones = chains(M)
	polynomials = Dict()

	pseudocubical_values = Dict()

	for cone in cones
		R = reduce(hcat, [B * bergman_ray(M, F) for F in cone])'
		pseudo_cubical = pseudo_cubical_cone(R)
		h_pc = Vector{Rational}(pseudo_cubical.REL_INT_POINT)
		pseudocubical_values[cone] = h_pc
	end

	for i in 1:length(cones)
		cone = cones[i]
		r = reduce(hcat, [B * bergman_ray(M, F) for F in cone])'
		vars = Sym.variables(:h, cone)
		for J in powerset(1:length(cone), 1)
			subcone = cone[J]
			subcone_r = r[J, :]
			subcone_vars = vars[J]

			f = cone_ehrhart_polynomial(r', subcone_r, subcone_vars, pseudocubical_values[subcone])
			polynomials[(cone, subcone)] = f
			for J_ in powerset(J, 1)
				subsubcone = cone[J_]
				subsubcone_r = r[J_, :]
				subsubcone_vars = vars[J_]

				lat_basis = lattice_basis_inside_span(subcone_r * r')
	
				f = cone_ehrhart_polynomial(r' * lat_basis', subsubcone_r, subsubcone_vars, pseudocubical_values[subsubcone])
				polynomials[(cone, subcone, subsubcone)] = f
			end
		end
    end
	
    return polynomials
end

# ╔═╡ 634a0421-25b5-4bb6-963f-916f310ab4ed
x = Sym.variables(:x, 1:10)

# ╔═╡ 097e42ca-56ed-4a03-8ba8-0d8faafa849d
function ehrhart_polynomial(M; B=I, polynomials=nothing)
	if polynomials == nothing
		polynomials = compute_polynomials(M; B=B)
	end
	
    n = length(nonempty_proper_flats(M))

	cones = chains(M)
	dim = maximum(length.(cones))
	cones_by_dim = [filter(c -> length(c) == d, cones) for d in 1:dim]	
    ehr = 0
	if (dim >= 3)
	    for (i, cone) in enumerate(cones_by_dim[3])
	        r = reduce(hcat, [bergman_ray(M, F) for F in cone])'
	        vars = Sym.variables(:h, cone)
	
			result = 0
			result += polynomials[(cone, cone)]
			for j in 1:3
				J = [k != j for k in 1:3]
				subcone = cone[J]
				result -= 1//2 * (polynomials[(cone, subcone)] + polynomials[(subcone, subcone)])
	
				for k in (1:3)[J]
					ray = cone[k:k]
					#result += 1//6 * (polynomials[(cone, ray)] + polynomials[(subcone, ray)] + polynomials[(ray, ray)])
					#result += x[1] * (polynomials[(cone, subcone, ray)])
					result += 1//6 * polynomials[(subcone, ray)]
				end
	
				result += 1//3 *(polynomials[(cone, cone[j:j])] + polynomials[(cone[j:j], cone[j:j])])
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

    ehr = Sym.expand(ehr)
    return ehr
end

# ╔═╡ f7f4ca84-bf8f-4962-a898-b92412db4bdd
# ╠═╡ disabled = true
#=╠═╡
M = Matroid(matroid.matroid_from_graph(graph.graph_from_edges([[0, 1], [1, 2], [0, 2], [1, 3], [2, 3], [3, 4], [1, 4]])))
  ╠═╡ =#

# ╔═╡ 5c91579c-34a0-448a-b5b4-08d95693ac10
#=╠═╡
svg_visual(M)
  ╠═╡ =#

# ╔═╡ 0bddbbd1-efa4-4d4e-9cb9-33ef6f83aae5
#=╠═╡
M.flats
  ╠═╡ =#

# ╔═╡ d65b475d-79f7-4227-a33c-045f83b3e5d0
#=╠═╡
[(F, M.flats[F], bergman_ray(M, F)) for F in nonempty_flats(M)]
  ╠═╡ =#

# ╔═╡ a06ec159-edaf-441e-92ae-ebe60330bbb4
#=╠═╡
polynomials = compute_polynomials(M)
  ╠═╡ =#

# ╔═╡ 8436dbb4-cbde-428b-a1e1-ca0c19270ab5
#=╠═╡
snap = snapper_polynomial(M);
  ╠═╡ =#

# ╔═╡ d8b1d015-016c-4a43-a85f-43c113f1c2f9
#=╠═╡
snap_h = Sym.expand(Sym.substitute(snap, f_to_h(M)));
  ╠═╡ =#

# ╔═╡ 37b7624b-f0d1-434f-83c6-1dbff3265525
#=╠═╡
p = snap_h - ehrhart_polynomial(M; polynomials=polynomials)
  ╠═╡ =#

# ╔═╡ 89e3aa92-021f-4a9e-8c5a-e17c01d51e0c
#=╠═╡
a = [Sym.coeff(p, Sym.variable(:h, i)) for i in 2:26]
  ╠═╡ =#

# ╔═╡ 1b3f732b-e2e1-409d-95b6-927a70936a78
#=╠═╡
collect(zip(a, [bergman_ray(M, F) for F in 2:26], M.flats[2:26]))[a .!= 0]
  ╠═╡ =#

# ╔═╡ 0f398106-1db4-491c-8ed2-24978cf9cc21
hypercube_volume(Vector{Rational}[[0, 0, 0], [1, 0,0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]])

# ╔═╡ 63910e9a-8c8c-4144-b6ca-ca15dc59bbde
hypercube_volume(Vector{Rational}[[0, 0], [1, 0], [0, 1], [2, 1]])

# ╔═╡ 464c9484-8dc0-4d36-beb2-54abdcdf329e
F = Matroid(matroid.fano_matroid())

# ╔═╡ fb6c5e16-11fd-414e-8347-64dc12de97b7


# ╔═╡ e9b461a5-231f-4776-ad44-9344a8bb182c
P = Matroid(matroid.uniform_matroid(4,4))

# ╔═╡ 2e6ca780-e780-46cc-8d15-f8504ea63dc9
g = Sym.variables(:g, 1:15)

# ╔═╡ d05c8cc6-95c1-4707-a483-91f412dcb368
sym_dict = Dict(Sym.variables(:h, 2:15) .=> [g[1], g[1], g[1], g[1], g[2], g[2], g[2], g[2], g[2], g[2], g[3], g[3], g[3], g[3]])

# ╔═╡ 0f846dbf-194e-4052-a148-ad93fe93a7ce
volume_polynomial(P)

# ╔═╡ 6ffc1041-21db-4030-88fc-82d35ceb5584
f = Sym.expand(Sym.substitute(volume_polynomial(P), sym_dict))

# ╔═╡ 9463d81b-f700-4560-a12f-488cbc323867
Sym.polynomial_coeffs(f, Sym.variables(:h, 2:15))

# ╔═╡ 0665a407-0de6-40b0-b312-cdb09a76fc3a
f_poch = Sym.expand(to_falling_factorials(f, Sym.variables(:g, 1:3)))

# ╔═╡ 6eb16052-ee49-44fd-921b-90a0b47db1d7
h = f_poch - f

# ╔═╡ 1b5990cf-17ab-4676-b63e-0317602738f6
df = [Sym.derivative(f, Sym.variable(:g, i)) for i in 1:3]

# ╔═╡ e4dbe57f-062b-4ea5-a8c7-443c3e54491b
Sym.expand(sum(df .* [-1, -0, -1]))

# ╔═╡ 94195c83-7a5d-4cec-ad1d-4a2cae1e4006
res = Sym.expand(h - sum(df .* [1//2, 0, 1//2]))

# ╔═╡ 748443e6-d85f-46cb-b7d7-83a942e8e841
sqs2 = [Sym.coeff(Sym.expand(g - sum(df .* Sym.variables(:a, 2:15))), Sym.variable(:h, i)*Sym.variable(:h, 2)) for i in 2:15]

# ╔═╡ a36264ce-ca20-4008-afd1-e849b680596b
Sym.symbolic_solve(sqs2 .~ 0, Sym.variables(:a, [2, 6,7,8, 12, 13, 14]))

# ╔═╡ 253b9439-e316-4185-b9d2-3223ec733872
sqs = vcat([Sym.coeff(Sym.expand(h - sum(df .* Sym.variables(:a, 1:3))), Sym.variable(:g, i)* Sym.variable(:g, j)) for i in 1:3 for j in 1:3]
)

# ╔═╡ e5c106d6-0002-46d2-9ff6-d21bc0adfc7d
eqs = [exp ~ 0 for exp in sqs if !Sym._iszero(exp)]

# ╔═╡ 2eac58c8-5c22-4426-af72-06f81c11ba5c
Sym.symbolic_solve(eqs[1:3], Sym.variables(:a, 1:3))

# ╔═╡ 82407c48-4256-4ea2-ac1f-accbc99c7b42
A = hcat([[Sym.coeff(a, Sym.variable(:a, i)) for i in 2:15] for a in sqs]...)

# ╔═╡ b908088e-324c-42a0-b019-150eec489d5b
A * A'

# ╔═╡ eb0939cc-dcc0-49e0-b964-134bbfac24cd
#=╠═╡
vol_M = volume_polynomial(M)
  ╠═╡ =#

# ╔═╡ 53bdbe99-d1b0-460f-b116-d8d480e6f938
#=╠═╡
diff = vol_M - snap_h
  ╠═╡ =#

# ╔═╡ f59dfc1a-71bc-48f0-982c-c7a2a984a8e0
#=╠═╡
fs = [Sym.coeff(diff, Sym.variable(:h, i)) for i in 1:34]
  ╠═╡ =#

# ╔═╡ a0a43e95-9d74-458f-a8af-49f0a4138245
#=╠═╡
b = reduce(hcat, [[Sym.coeff(f, Sym.variable(:x, j)) for j in 1:3] for f in fs])
  ╠═╡ =#

# ╔═╡ 2217761d-f38b-4e79-8874-c3c679b622ed
#=╠═╡
y = [Sym.coeff(f, 1) for f in fs]
  ╠═╡ =#

# ╔═╡ 97c74a06-c53e-444c-bdd8-4335e2cb5826
#=╠═╡
y = [Sym.coeff(a, 1) for a in sqs]
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═312dfd9d-e4d0-4fac-a2eb-0191020ff66a
# ╠═845a0f6b-6628-4ee6-a89b-0b6ec31b6187
# ╠═634a0421-25b5-4bb6-963f-916f310ab4ed
# ╠═097e42ca-56ed-4a03-8ba8-0d8faafa849d
# ╠═f7f4ca84-bf8f-4962-a898-b92412db4bdd
# ╠═5c91579c-34a0-448a-b5b4-08d95693ac10
# ╠═0bddbbd1-efa4-4d4e-9cb9-33ef6f83aae5
# ╠═d65b475d-79f7-4227-a33c-045f83b3e5d0
# ╠═a06ec159-edaf-441e-92ae-ebe60330bbb4
# ╠═8436dbb4-cbde-428b-a1e1-ca0c19270ab5
# ╠═d8b1d015-016c-4a43-a85f-43c113f1c2f9
# ╠═37b7624b-f0d1-434f-83c6-1dbff3265525
# ╠═89e3aa92-021f-4a9e-8c5a-e17c01d51e0c
# ╠═1b3f732b-e2e1-409d-95b6-927a70936a78
# ╠═0f398106-1db4-491c-8ed2-24978cf9cc21
# ╠═63910e9a-8c8c-4144-b6ca-ca15dc59bbde
# ╠═9463d81b-f700-4560-a12f-488cbc323867
# ╠═464c9484-8dc0-4d36-beb2-54abdcdf329e
# ╠═fb6c5e16-11fd-414e-8347-64dc12de97b7
# ╠═e9b461a5-231f-4776-ad44-9344a8bb182c
# ╠═2e6ca780-e780-46cc-8d15-f8504ea63dc9
# ╠═d05c8cc6-95c1-4707-a483-91f412dcb368
# ╠═0f846dbf-194e-4052-a148-ad93fe93a7ce
# ╠═6ffc1041-21db-4030-88fc-82d35ceb5584
# ╠═0665a407-0de6-40b0-b312-cdb09a76fc3a
# ╠═6eb16052-ee49-44fd-921b-90a0b47db1d7
# ╠═1b5990cf-17ab-4676-b63e-0317602738f6
# ╠═e4dbe57f-062b-4ea5-a8c7-443c3e54491b
# ╠═94195c83-7a5d-4cec-ad1d-4a2cae1e4006
# ╠═748443e6-d85f-46cb-b7d7-83a942e8e841
# ╠═a36264ce-ca20-4008-afd1-e849b680596b
# ╠═253b9439-e316-4185-b9d2-3223ec733872
# ╠═e5c106d6-0002-46d2-9ff6-d21bc0adfc7d
# ╠═2eac58c8-5c22-4426-af72-06f81c11ba5c
# ╠═82407c48-4256-4ea2-ac1f-accbc99c7b42
# ╠═97c74a06-c53e-444c-bdd8-4335e2cb5826
# ╠═b908088e-324c-42a0-b019-150eec489d5b
# ╠═eb0939cc-dcc0-49e0-b964-134bbfac24cd
# ╠═53bdbe99-d1b0-460f-b116-d8d480e6f938
# ╠═f59dfc1a-71bc-48f0-982c-c7a2a984a8e0
# ╠═a0a43e95-9d74-458f-a8af-49f0a4138245
# ╠═2217761d-f38b-4e79-8874-c3c679b622ed
