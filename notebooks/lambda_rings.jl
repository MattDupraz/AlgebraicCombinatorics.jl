### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 1529505c-cc19-11ef-3358-6bc94b79ee9e
begin
	import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.instantiate()
	using NormalComplexes
	using Symbolics
end

# ╔═╡ d8cdba7d-df75-4231-b516-59a5e0cd6a96
function lambda(p)
	result = Series(1)
	for (monom, coeff) in polynomial_coeffs(p, Symbolics.get_variables(p))[1]
		result *= BinomialSeries(coeff) ⊙ GeometricSeries(monom)
	end
	return result
end

# ╔═╡ 279499c0-d576-492a-831e-e1af7837f30d
function apply_operator(op, polyn)
	op = expand(op)
	result = 0
	vars = Symbolics.get_variables(op)
	for (monom, coeff) in polynomial_coeffs(op, vars)[1]
		result += coeff * Symbolics.substitute(polyn, Dict(
			s => s + Symbolics.degree(monom, s) for s in vars
		))
	end
	return expand(result)
end

# ╔═╡ 63338baa-fa61-46b9-a4d6-7746d8f3f436
first_coeffs(LogSeries(), 5)

# ╔═╡ cb0bc7d6-c563-40e4-ac8a-414967da2acc
@variables x y z

# ╔═╡ 001b2d51-59e9-4ed7-9ea0-597e4093a4ff
first_coeffs(lambda(1//2), 5)

# ╔═╡ d4d1003d-728f-4962-ba54-77aab9ffeab9
apply_operator((x - 1)^2 - 2 * (y - 1), x^2 + y)

# ╔═╡ 642d37f9-5c10-4bae-bae1-c62bb13f4fc1
typeof(1 + x)

# ╔═╡ 201a9b6a-ff47-47ac-a914-a25277c50b8f
lambda((1 - x)*(1 - y))

# ╔═╡ 52a6ec22-c3f6-4b94-8b96-904c7a9da837
D(x) = 1 - x

# ╔═╡ c6e3138a-9f6e-4858-81eb-4c9a0324dbeb
expand(substitute(lambda(D(x) + D(x)^2 * (1//2) + D(x)^3 * (1//3) + D(x)^4 * (1//4) + D(x)^5 * (1//5) + D(x)^6 * (1//6) + D(x)^7 * (1//7) + D(x)^8 * (1//8))[1], Dict(x => 1 -x, y => 1 - y)))

# ╔═╡ b722d699-0ec1-4e0f-8c4b-64a9b2c90d0f
apply_operator(lambda(-2*D(x) - D(x)^2)[1], x^2 + x + x*y)

# ╔═╡ 695f16c9-cd75-41a2-9d3d-9948e5858a19
expand(substitute(lambda(D(x) - D(y))[2], Dict(x => 1 - x, y => 1 - y)))

# ╔═╡ 9a0ce425-00c5-4821-afb4-d51eca70db1f
apply_operator(lambda(D(x) + D(y))[1], x-y)

# ╔═╡ 08a461bb-3feb-4bdf-9e3b-3f07546e1ebf
expand(substitute(lambda(D(x)^2 + 2*D(y))[2], Dict(x => 1 - x, y => 1 - y, z=>1-z)))

# ╔═╡ 8ad78a55-9e1a-4949-9790-534d69c47531
apply_operator(lambda(D(x)^2 + 2*D(y))[1], x^2 + y)

# ╔═╡ 6c7c678a-e691-4705-bb4f-f20129ee09a0
expand(substitute(lambda(D(x)^2 - 2*D(y))[2], Dict(x => 1 - x, y => 1 - y, z=>1-z)))

# ╔═╡ b6eec44d-f5ae-4582-af1a-0b2600a1b5d4
apply_operator(lambda(D(x) + D(y) + D(z))[3], x + y - 2 * z)

# ╔═╡ ce2abd81-2c0d-485e-a25a-79f78b8807d1
lambda(-6D(x)-3D(x)^2-1)

# ╔═╡ cb729f7d-fe0a-4dee-a029-3a63fa2aa528
expand(substitute(lambda(-24D(x)-12D(x)^2-4D(x)^3-D(x)^4)[4], Dict(x => 1 - x, y => 1 - y, z=>1-z)))

# ╔═╡ Cell order:
# ╠═1529505c-cc19-11ef-3358-6bc94b79ee9e
# ╠═d8cdba7d-df75-4231-b516-59a5e0cd6a96
# ╠═279499c0-d576-492a-831e-e1af7837f30d
# ╠═63338baa-fa61-46b9-a4d6-7746d8f3f436
# ╠═cb0bc7d6-c563-40e4-ac8a-414967da2acc
# ╠═001b2d51-59e9-4ed7-9ea0-597e4093a4ff
# ╠═d4d1003d-728f-4962-ba54-77aab9ffeab9
# ╠═642d37f9-5c10-4bae-bae1-c62bb13f4fc1
# ╠═201a9b6a-ff47-47ac-a914-a25277c50b8f
# ╠═c6e3138a-9f6e-4858-81eb-4c9a0324dbeb
# ╠═52a6ec22-c3f6-4b94-8b96-904c7a9da837
# ╠═b722d699-0ec1-4e0f-8c4b-64a9b2c90d0f
# ╠═695f16c9-cd75-41a2-9d3d-9948e5858a19
# ╠═9a0ce425-00c5-4821-afb4-d51eca70db1f
# ╠═08a461bb-3feb-4bdf-9e3b-3f07546e1ebf
# ╠═8ad78a55-9e1a-4949-9790-534d69c47531
# ╠═6c7c678a-e691-4705-bb4f-f20129ee09a0
# ╠═b6eec44d-f5ae-4582-af1a-0b2600a1b5d4
# ╠═ce2abd81-2c0d-485e-a25a-79f78b8807d1
# ╠═cb729f7d-fe0a-4dee-a029-3a63fa2aa528
