### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ a3e5efb2-f8ef-45e7-be41-b1b05a6253bb
using Polymake: Polymake, polytope, matroid, graph, LibPolymake.BigObjectAllocated as PMBigObject, @pm, common

# ╔═╡ 1c9ecfc8-a36a-477b-9670-822fe20bd873
using DataStructures

# ╔═╡ fae50ce3-b791-4012-a202-57c6cdf8db6f
using Combinatorics:powerset, permutations

# ╔═╡ afd7d500-7975-48e9-9467-2c6c59ae824c
using SimplePlutoInclude

# ╔═╡ 747f9307-4f18-4d8d-b85d-691da228aeb1
using LinearAlgebra: I, det

# ╔═╡ 1300785b-a1b4-469b-af15-d74269b8f790
using SHA

# ╔═╡ df693074-1d7b-4280-a7d2-7aa1ed25b691
using Serialization

# ╔═╡ f054e974-f9cc-4065-97ae-3b0aa1ab4522
using XML

# ╔═╡ 5b2c43a4-dd97-11ef-1c79-1b5007c54b8d
import Symbolics as Sym

# ╔═╡ 312dfd9d-e4d0-4fac-a2eb-0191020ff66a
@plutoinclude "NormalComplexes.jl" all=true

# ╔═╡ 81fdb0ae-8327-4ffa-8b40-6af55468bd01
@plutoinclude "Pochhammer.jl" all=true

# ╔═╡ fe958ff4-3294-4516-a97c-d93f48a765d1
@plutoinclude "Cache.jl"

# ╔═╡ 845a0f6b-6628-4ee6-a89b-0b6ec31b6187
function compute_polynomials(M; B=I)	
	cones = chains(M)
	polynomials = Dict()

	pseudocubical_values = Dict()

	for cone in cones
		R = reduce(hcat, [B * get_ray_gen(M, F) for F in cone])'
		pseudo_cubical = get_pseudo_cubical_cone(R)
		h_pc = Vector{Rational}(pseudo_cubical.REL_INT_POINT)
		pseudocubical_values[cone] = h_pc
	end

	for i in 1:length(cones)
		cone = cones[i]
		r = reduce(hcat, [B * get_ray_gen(M, F) for F in cone])'
		vars = Sym.variables(:h, cone)
		for J in powerset(1:length(cone), 1)
			subcone = cone[J]
			subcone_r = r[J, :]
			subcone_vars = vars[J]

			f = cone_ehr(r', subcone_r, subcone_vars, pseudocubical_values[subcone])
			polynomials[(cone, subcone)] = f
			for J_ in powerset(J, 1)
				subsubcone = cone[J_]
				subsubcone_r = r[J_, :]
				subsubcone_vars = vars[J_]

				lat_basis = lattice_basis_inside_span(subcone_r * r')
	
				f = cone_ehr(r' * lat_basis', subsubcone_r, subsubcone_vars, pseudocubical_values[subsubcone])
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
	        r = reduce(hcat, [get_ray_gen(M, F) for F in cone])'
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
[(F, M.flats[F], get_ray_gen(M, F)) for F in nonempty_flats(M)]
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
collect(zip(a, [get_ray_gen(M, F) for F in 2:26], M.flats[2:26]))[a .!= 0]
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
f_poch = Sym.expand(pochhammer(f, Sym.variables(:g, 1:3)))

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Polymake = "d720cf60-89b5-51f5-aff5-213f193123e7"
SHA = "ea8e919c-243c-51af-8825-aaa63cd721ce"
Serialization = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
SimplePlutoInclude = "6f00a2c5-ea4a-46bf-9183-91b7b57a087f"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
XML = "72c71f33-b9b6-44de-8c94-c961784809e2"

[compat]
Combinatorics = "~1.0.2"
DataStructures = "~0.18.20"
Latexify = "~0.16.6"
Polymake = "~0.11.25"
SimplePlutoInclude = "~0.2.0"
Symbolics = "~6.29.1"
XML = "~0.3.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.1"
manifest_format = "2.1"
project_hash = "00c4b1d5a8146e85ee81d8f4f81c41745147dfe9"

[[deps.ADTypes]]
git-tree-sha1 = "fb97701c117c8162e84dfcf80215caa904aef44f"
registries = "General"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.13.0"

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
registries = "General"
uuid = "ae81ac8f-d209-56e5-92de-9978fef736f9"
version = "0.1.3+0"

[[deps.AbstractAlgebra]]
deps = ["LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays"]
git-tree-sha1 = "1c765891be37479d2ab6f4fe1faff30fea27fc2f"
registries = "General"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.44.8"
weakdeps = ["Requires", "Test"]

    [deps.AbstractAlgebra.extensions]
    TestExt = "Test"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
registries = "General"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
registries = "General"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

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
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cd8b948862abee8f3d3e9b73a102a9ca924debb0"
registries = "General"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.2.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
registries = "General"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "017fcb757f8e921fb44ee063a7aafe5f89b86dd1"
registries = "General"
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
registries = "General"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.9"

[[deps.BinaryWrappers]]
deps = ["JLLWrappers", "Scratch"]
git-tree-sha1 = "7fea8f658689fa5062b23f4400eda888b7ae2aaa"
registries = "General"
uuid = "f01c122e-0ea1-4f85-ad8f-907073ad7a9f"
version = "0.1.3"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
registries = "General"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
registries = "General"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
registries = "General"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
registries = "General"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
registries = "General"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
registries = "General"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
registries = "General"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
registries = "General"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
registries = "General"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.CxxWrap]]
deps = ["Libdl", "MacroTools", "libcxxwrap_julia_jll"]
git-tree-sha1 = "12349409599af22cfe0a7675857d969fe71b7767"
registries = "General"
uuid = "1f15a43c-97ca-5a2a-ae31-89f07a497df4"
version = "0.16.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
registries = "General"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
registries = "General"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
registries = "General"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DecFP]]
deps = ["DecFP_jll", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "88e521a871a1b11488a1f48b8c085b4be8f71be5"
registries = "General"
uuid = "55939f99-70c6-5e9b-8bb0-5071ed7d61fd"
version = "1.4.1"

[[deps.DecFP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e9a8da19f847bbfed4076071f6fef8665a30d9e5"
registries = "General"
uuid = "47200ebd-12ce-5be5-abb7-8e082af23329"
version = "2.0.3+1"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
registries = "General"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "03aa5d44647eaec98e1920635cdfed5d5560a8b9"
registries = "General"
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
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "a7e9f13f33652c533d49868a534bfb2050d1365f"
registries = "General"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.15"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "9a3ae38b460449cc9e7dd0cfb059c76028724627"
registries = "General"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.1"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
registries = "General"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
registries = "General"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
registries = "General"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FLINT_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "MPFR_jll", "OpenBLAS32_jll"]
git-tree-sha1 = "9327ebef3e04034b2aa1b9c59869f36717843ce5"
registries = "General"
uuid = "e134572f-a0d5-539d-bddf-3cad8db41a82"
version = "300.100.301+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
registries = "General"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
registries = "General"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
registries = "General"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
registries = "General"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLPK_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6aa6294ba949ccfc380463bf50ff988b46de5bc7"
registries = "General"
uuid = "e8aa6df9-e6ca-548a-97ff-1f85fc5b8b98"
version = "5.0.1+1"

[[deps.GMP_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"
version = "6.3.0+2"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
registries = "General"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f93a9ce66cd89c9ba7a4695a47fd93b4c6bc59fa"
registries = "General"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.12.0+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "2bd56245074fab4015b9174f24ceba8293209053"
registries = "General"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.27"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
registries = "General"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
registries = "General"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
registries = "General"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
registries = "General"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.Ipopt_jll]]
deps = ["ASL_jll", "Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "MUMPS_seq_jll", "SPRAL_jll", "libblastrampoline_jll"]
git-tree-sha1 = "4f55ad688c698a4f77d892a1cb673f7e8a30f178"
registries = "General"
uuid = "9cc047cb-c261-5740-88fc-0cf96f7bdcc7"
version = "300.1400.1700+0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
registries = "General"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
registries = "General"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
registries = "General"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
registries = "General"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cd714447457c660382fe634710fb56eb255ee42e"
registries = "General"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.6"

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
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.18.0+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.1+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
registries = "General"
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

[[deps.METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2eefa8baa858871ae7770c98c3c2a7e46daba5b4"
registries = "General"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.3+0"

[[deps.MPFR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "GMP_jll", "Libdl"]
uuid = "3a97d323-0669-5f0c-9066-3539efd106a3"
version = "4.2.2+0"

[[deps.MUMPS_seq_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "libblastrampoline_jll"]
git-tree-sha1 = "0eab12f94948ca67908aec14b9f2ebefd17463fe"
registries = "General"
uuid = "d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"
version = "500.700.301+0"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
registries = "General"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
registries = "General"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MongoC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "snappy_jll"]
git-tree-sha1 = "0ab1e13479581def03ea0b19733b8b830177889a"
registries = "General"
uuid = "90100e71-7732-535a-9be7-2e9affd1cfc1"
version = "1.28.1+0"

[[deps.Mongoc]]
deps = ["Dates", "DecFP", "MongoC_jll", "Serialization"]
git-tree-sha1 = "9b0ee6f5843a8df838d76794a5d50e77b8a3b08f"
registries = "General"
uuid = "4fe8b98c-fc19-5c23-8ec2-168ff83495f2"
version = "0.9.2"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "453de0fc2be3d11b9b93ca4d0fddd91196dcf1ed"
registries = "General"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.5"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "8d39779e29f80aa6c071e7ac17101c6e31f075d7"
registries = "General"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "491bdcdc943fcbc4c005900d7463c9f216aabf4c"
registries = "General"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "cc0a5deefdb12ab3a096f00a6d42133af4560d71"
registries = "General"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.2"

[[deps.Ncurses_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b5e7e7ad16adfe5f68530f9f641955b5b0f12bbb"
registries = "General"
uuid = "68e3532b-a499-55ff-9963-d1c0c0748b3a"
version = "6.5.1+0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Ninja_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d9c6d1dca818768f1671244dcdf8e8e08b2b8937"
registries = "General"
uuid = "76642167-d241-5cee-8c94-7a494e8cb7b7"
version = "1.12.1+0"

[[deps.OffsetArrays]]
git-tree-sha1 = "5e1897147d1ff8d98883cda2be2187dcf57d8f0c"
registries = "General"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.15.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ece4587683695fe4c5f20e990da0ed7e83c351e7"
registries = "General"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
registries = "General"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "966b85253e959ea89c53a9abebbf2e964fbf593b"
registries = "General"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.32"

[[deps.PPL_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "dc69ba2b374b5bd6e876fa3a8441563c0eeafb9e"
registries = "General"
uuid = "80dd9cbb-8b87-5171-a280-372cc418f402"
version = "1.2.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Perl_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Readline_jll"]
git-tree-sha1 = "7ab65a258bcf6da373cab49af462aead452d3960"
registries = "General"
uuid = "83958c19-0796-5285-893e-a1267f8ec499"
version = "5.34.1+0"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
registries = "General"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.Polymake]]
deps = ["AbstractAlgebra", "BinaryWrappers", "CxxWrap", "Downloads", "JSON", "Libdl", "Mongoc", "NetworkOptions", "Ninja_jll", "Perl_jll", "Pidfile", "Pkg", "REPL", "Scratch", "SparseArrays", "TOPCOM_jll", "lib4ti2_jll", "libpolymake_julia_jll", "polymake_jll", "polymake_oscarnumber_jll"]
git-tree-sha1 = "5addbf9d589c4d052ec0aef4bd5cc7e886d12c00"
registries = "General"
uuid = "d720cf60-89b5-51f5-aff5-213f193123e7"
version = "0.11.25"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
registries = "General"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
registries = "General"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
registries = "General"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "b8a399e95663485820000f26b6a43c794e166a49"
registries = "General"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.4"

[[deps.Readline_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ncurses_jll"]
git-tree-sha1 = "6044f482a91c7aa2b82ab614aedd726be633ad05"
registries = "General"
uuid = "05236dd9-4125-5232-aa7c-9ec0c9b2c25a"
version = "8.2.13+0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
registries = "General"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "e96b644f7bfbf1015f8e42a7c7abfae2a48fafbf"
registries = "General"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.31.0"

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
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
registries = "General"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
registries = "General"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
registries = "General"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SCIP_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "GMP_jll", "Ipopt_jll", "JLLWrappers", "Libdl", "Pkg", "Readline_jll", "Zlib_jll", "bliss_jll", "boost_jll"]
git-tree-sha1 = "cf69186eb29307fbb2319b90e6133797bad983ce"
registries = "General"
uuid = "e5ac4fe4-a920-5659-9bf8-f9f73e9e79ce"
version = "800.0.301+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.SPRAL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "Libdl", "METIS_jll", "libblastrampoline_jll"]
git-tree-sha1 = "11f3da4b25efacd1cec8e263421f2a9003a5e8e0"
registries = "General"
uuid = "319450e9-13b8-58e8-aa9f-8fd1420848ab"
version = "2024.5.8+0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "ee305515b0946db5f56af699e8b5804fee04146c"
registries = "General"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.75.1"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "6149620767866d4b0f0f7028639b6e661b6a1e44"
registries = "General"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.12"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "0444a37a25fab98adbd90baa806ee492a3af133a"
registries = "General"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.6.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
registries = "General"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
registries = "General"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SimplePlutoInclude]]
deps = ["Dates"]
git-tree-sha1 = "a98dee6bdc63f648d92b28348dcb9b59c2881fbf"
registries = "General"
uuid = "6f00a2c5-ea4a-46bf-9183-91b7b57a087f"
version = "0.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
registries = "General"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.13.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
registries = "General"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "e3be13f448a43610f978d29b7adf78c76022467a"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.12"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
registries = "General"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
registries = "General"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
registries = "General"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
registries = "General"
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
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.10.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "d6c04e26aa1c8f7d144e1a8c47f1c73d3013e289"
registries = "General"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.38"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fabf4650afe966a2ba646cabd924c3fd43577fc3"
registries = "General"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "ExproniconLite", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "TimerOutputs", "Unityper", "WeakValueDicts"]
git-tree-sha1 = "e893939e258de0de23dc233ec6a618a4b53554ff"
registries = "General"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.27.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "IfElse", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "OffsetArrays", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "e021ac0b4dc29ce0a5b38fd91c27bd36d6d7df13"
registries = "General"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.29.1"

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
registries = "General"
uuid = "36f60fef-b880-50dc-9289-4aaecee93cc3"
version = "0.17.8+0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
registries = "General"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
registries = "General"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "d155450e6dff2a8bc2fcb81dcb194bd98b0aeb46"
registries = "General"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.2"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
registries = "General"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3832505b94c1868baea47764127e6d36b5c9f29e"
registries = "General"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.27"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

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
registries = "General"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.WeakValueDicts]]
git-tree-sha1 = "98528c2610a5479f091d470967a25becfd83edd0"
registries = "General"
uuid = "897b6980-f191-5a31-bcb0-bf3c4585e0c1"
version = "0.1.0"

[[deps.XML]]
deps = ["Mmap", "OrderedCollections"]
git-tree-sha1 = "360c178fdee6816bcc9d436c6f324d582e84f1d1"
registries = "General"
uuid = "72c71f33-b9b6-44de-8c94-c961784809e2"
version = "0.3.4"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.bliss_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f8b75e896a326a162a4f6e998990521d8302c810"
registries = "General"
uuid = "508c9074-7a14-5c94-9582-3d4bc1871065"
version = "0.77.0+1"

[[deps.boost_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7a89efe0137720ca82f99e8daa526d23120d0d37"
registries = "General"
uuid = "28df3c45-c428-5900-9ff8-a3135698ca75"
version = "1.76.0+1"

[[deps.cddlib_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "08f5df03703af917b9bfec47b9767eb943220d08"
registries = "General"
uuid = "f07e07eb-5685-515a-97c8-3014f6152feb"
version = "0.94.14+0"

[[deps.lib4ti2_jll]]
deps = ["Artifacts", "GLPK_jll", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79116185def638e73ea3d88ca6c10e210a1dc183"
registries = "General"
uuid = "1493ae25-0f90-5c0e-a06c-8c5077d6d66f"
version = "1.6.10+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libcxxwrap_julia_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70091a11c756455c15bc6380ea000c16f38f2d3e"
registries = "General"
uuid = "3eaa8342-bff7-56a5-9981-c04077f7cee7"
version = "0.13.4+0"

[[deps.libpolymake_julia_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "FLINT_jll", "JLLWrappers", "Libdl", "TOPCOM_jll", "lib4ti2_jll", "libcxxwrap_julia_jll", "polymake_jll"]
git-tree-sha1 = "d05957bb0eacb2747e6d042d56034b005c220561"
registries = "General"
uuid = "4d8266f6-2b3b-57e3-ad7a-d431eaaac945"
version = "0.13.1+0"

[[deps.lrslib_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9169c4d5823195a7a2173fc6b479468478857438"
registries = "General"
uuid = "3873f7d0-7b7c-52c3-bdf4-8ab39b8f337a"
version = "0.3.3+0"

[[deps.nauty_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "fc7dc197df0648cd5f965801bfe086abd9325add"
registries = "General"
uuid = "55c6dc9b-343a-50ca-8ff2-b71adb3733d5"
version = "2.6.13+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.67.1+0"

[[deps.normaliz_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "FLINT_jll", "GMP_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl", "MPFR_jll", "nauty_jll"]
git-tree-sha1 = "b84641c3add3c9e435546b1f2867105ce209791f"
registries = "General"
uuid = "6690c6e9-4e12-53b8-b8fd-4bffaef8839f"
version = "300.1000.200+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.2+0"

[[deps.polymake_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "FLINT_jll", "GMP_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl", "MPFR_jll", "MongoC_jll", "PPL_jll", "Perl_jll", "SCIP_jll", "bliss_jll", "boost_jll", "cddlib_jll", "lrslib_jll", "normaliz_jll"]
git-tree-sha1 = "bfd75c077ece1a66325ec61a6cb6b0b23420f4b2"
registries = "General"
uuid = "7c209550-9012-526c-9264-55ba7a78ba2c"
version = "400.1300.2+0"

[[deps.polymake_oscarnumber_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl", "libcxxwrap_julia_jll", "libpolymake_julia_jll", "polymake_jll"]
git-tree-sha1 = "14b654205ab1599612521c587f5b11c662ba0ed1"
registries = "General"
uuid = "10f31823-b687-53e6-9f29-edb9d4da9f9f"
version = "0.3.4+0"

[[deps.snappy_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8bb12ecfc3ffde2ef72a957cf2fab27cc193d251"
registries = "General"
uuid = "fe1e1685-f7be-5f59-ac9f-4ca204017dfd"
version = "1.2.2+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╠═a3e5efb2-f8ef-45e7-be41-b1b05a6253bb
# ╠═1c9ecfc8-a36a-477b-9670-822fe20bd873
# ╠═5b2c43a4-dd97-11ef-1c79-1b5007c54b8d
# ╠═fae50ce3-b791-4012-a202-57c6cdf8db6f
# ╠═afd7d500-7975-48e9-9467-2c6c59ae824c
# ╠═747f9307-4f18-4d8d-b85d-691da228aeb1
# ╠═1300785b-a1b4-469b-af15-d74269b8f790
# ╠═df693074-1d7b-4280-a7d2-7aa1ed25b691
# ╠═f054e974-f9cc-4065-97ae-3b0aa1ab4522
# ╠═312dfd9d-e4d0-4fac-a2eb-0191020ff66a
# ╠═81fdb0ae-8327-4ffa-8b40-6af55468bd01
# ╠═fe958ff4-3294-4516-a97c-d93f48a765d1
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
