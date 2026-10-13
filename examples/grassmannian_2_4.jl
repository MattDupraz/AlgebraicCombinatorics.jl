# Todd class of the Grassmannian Gr(2, 4), in the Chern classes of the tautological bundle.
# Run from the repository root with: julia --project=. examples/grassmannian_2_4.jl

using NormalComplexes
using Symbolics

k = 2
n = 4
T = todd_class_grassmannian(k, n)

# Since [S ⊗ Q] = [O], we know that the product of the chern classes is 1 in the Chow ring
S_inv = inv(Series([1; Symbolics.variables(:s, 1:2)]))
q_to_s = Dict([Symbolics.variable(:q, i) => S_inv[i] for i in 1:2])

T = substitute(T, q_to_s)

println(first_coeffs(T, k*(n - k)))
