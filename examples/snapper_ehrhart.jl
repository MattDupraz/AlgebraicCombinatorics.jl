# The Ehrhart polynomial of the normal complex of a matroid agrees with its Snapper
# polynomial, after the change of variables f_to_h.
#
# Matroid on a, b, c, d with maximal proper flats ad, bd, cd, abc: a, b, c are
# collinear and d lies off their line.
#
# Run from the repository root with: julia --project=. examples/snapper_ehrhart.jl

using NormalComplexes
using Polymake: matroid
using Symbolics

# Realisation: a, b, c in the plane z = 0 (with c = a + b), d outside it
M = Matroid(matroid.Matroid(VECTORS=[1 0 0; 0 1 0; 1 1 0; 0 0 1]))

# Name each flat by its elements, e.g. flat {0, 3} is "ad"
flat_name(F) = join(sort(collect(M.flats[F])) .|> e -> "abcd"[e + 1])
println("Maximal proper flats: ", [flat_name(F) for F in flats_of_rank_range(M, M.rank - 1, M.rank - 1)])

# Rename the variables f_F and h_F after their flats, for printing
readable = merge(
    Dict(Symbolics.variable(:f, F) => Symbolics.variable(Symbol("f_", flat_name(F))) for F in nonempty_flats(M)),
    Dict(Symbolics.variable(:h, F) => Symbolics.variable(Symbol("h_", flat_name(F))) for F in nonempty_proper_flats(M)))
show_readable(p) = Symbolics.expand(Symbolics.substitute(p, readable))

snap = snapper_polynomial(M)
ehr = ehrhart_polynomial(M)

println("\nSnapper polynomial:\n  ", show_readable(snap))
println("\nEhrhart polynomial:\n  ", show_readable(ehr))

println("\nChange of variables f_to_h:")
substitution = f_to_h(M)
for F in nonempty_flats(M)
    f = Symbolics.variable(:f, F)
    println("  ", show_readable(f), " = ", show_readable(substitution[f]))
end

# The Snapper polynomial in the h variables is the Ehrhart polynomial, and the
# Ehrhart polynomial in the f variables is the Snapper polynomial
snap_h = Symbolics.expand(Symbolics.substitute(snap, f_to_h(M)))
ehr_f = Symbolics.expand(Symbolics.substitute(ehr, h_to_f(M)))
println("\nSnapper(f_to_h) - Ehrhart = ", Symbolics.expand(snap_h - ehr))
println("Ehrhart(h_to_f) - Snapper = ", Symbolics.expand(ehr_f - snap))
