# The Snapper polynomial of a matroid can be computed from its volume polynomial:
# write the volume polynomial in the f basis, replace each monomial
# Π_F f_F^(a_F) by the product of rising factorials Π_F f_F^(a_F rising), and
# translate by 1 in the direction f_E (E the ground set): f_E ↦ f_E + 1.
# This is snapper_polynomial_from_volume; the steps are spelled out below.
#
# Matroid on a, b, c, d with maximal proper flats ad, bd, cd, abc: a, b, c are
# collinear and d lies off their line.
#
# Run from the repository root with: julia --project=. examples/snapper_volume.jl

using AlgebraicCombinatorics
using Polymake: matroid
using Symbolics

# Realisation: a, b, c in the plane z = 0 (with c = a + b), d outside it
M = Matroid(matroid.Matroid(VECTORS=[1 0 0; 0 1 0; 1 1 0; 0 0 1]))

# Name each flat by its elements, e.g. flat {0, 3} is "ad"
flat_name(F) = join(sort(collect(M.flats[F])) .|> e -> "abcd"[e + 1])

# Rename the variables f_F after their flats, for printing
readable = Dict(Symbolics.variable(:f, F) => Symbolics.variable(Symbol("f_", flat_name(F))) for F in nonempty_flats(M))
show_readable(p) = Symbolics.expand(Symbolics.substitute(p, readable))

f = Symbolics.variables(:f, nonempty_flats(M))
f_E = Symbolics.variable(:f, flat_index(M, Set(ground_set(M))))

vol_f = Symbolics.expand(Symbolics.substitute(volume_polynomial(M), h_to_f(M)))
rising = to_rising_factorials(vol_f, f)
translated = Symbolics.expand(Symbolics.substitute(rising, Dict(f_E => f_E + 1)))
snap = snapper_polynomial(M)

println("Volume polynomial in the f basis:\n  ", show_readable(vol_f))
println("\nIn rising factorials:\n  ", show_readable(rising))
println("\nTranslated by f_abcd ↦ f_abcd + 1:\n  ", show_readable(translated))
println("\nSnapper polynomial:\n  ", show_readable(snap))
println("\nDifference: ", Symbolics.expand(translated - snap))

# The same holds in rank 4, where ehrhart_polynomial is not implemented
U = Matroid(matroid.uniform_matroid(4, 5))
println("\nU(4, 5), difference: ", Symbolics.expand(snapper_polynomial_from_volume(U) - snapper_polynomial(U)))
