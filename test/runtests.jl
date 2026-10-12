using AlgebraicCombinatorics
using Test
using LinearAlgebra: det, I
using Polymake: matroid
using Symbolics

# Symbolic equality after expansion
≂(a, b) = iszero(Num(Symbolics.expand(a - b)))

@testset "AlgebraicCombinatorics" begin

@testset "Series" begin
    @variables x t
    @test first_coeffs(BinomialSeries(5), 6) == [binomial(5, i) for i in 0:6]
    @test BinomialSeries(x)[2] ≂ x * (x - 1) / 2
    @test first_coeffs(ExpSeries() * ExpSeries(-1), 5) == [1; zeros(Int, 5)]
    # log(1 + (exp(t) - 1)) = t
    @test first_coeffs(LogSeries() ∘ (ExpSeries() - 1), 5) == [0, 1, 0, 0, 0, 0]
    @test first_coeffs(GeometricSeries(2) ⊙ GeometricSeries(3), 3) == [1, 6, 36, 216]
    @test first_nonzero_index(Series(0, 0, 3)) == 2
end

@testset "SymmetricPolynomials" begin
    xs = Symbolics.variables(:x, 1:3)
    e = Symbolics.variables(:e, 1:3)
    p = sum(xs .^ 2)
    @test to_elementary(p, xs) ≂ e[1]^2 - 2e[2]
    @test substitute_elementary(to_elementary(p, xs), xs) ≂ p
end

@testset "Pochhammer" begin
    @variables a b
    @test falling_factorial(a, 3) ≂ a * (a - 1) * (a - 2)
    @test rising_factorial(a, 3) ≂ a * (a + 1) * (a + 2)
    @test to_falling_factorials(a^2 * b + a, [a, b]) ≂ a * (a - 1) * b + a
    @test to_rising_factorials(a^2, [a]) ≂ a * (a + 1)
    xs = Symbolics.variables(:x, 1:3)
    @test vandermonde(xs) ≂ (xs[1] - xs[2]) * (xs[1] - xs[3]) * (xs[2] - xs[3])
end

@testset "Grassmannians" begin
    # Todd class of a line bundle: 1 + c/2 + c^2/12 - c^4/720
    c = Symbolics.variable(:c, 1)
    T = todd_class(1)
    @test all(first_coeffs(T, 4) .≂ [1, c / 2, c^2 / 12, 0, -c^4 / 720])
end

@testset "Cache" begin
    cd(mktempdir()) do
        @test (@cache 1 + 1) == 2
        @test (@cache 1 + 1) == 2  # read back from disk
    end
end

@testset "LinAlgTools" begin
    B = [1 1 1; -1 0 2; 3 5 6]
    R, U = lll_reduction(B)
    @test U * B == R
    @test abs(det(Rational{BigInt}.(U))) == 1
    @test shortest_vector_approximation(B)[1] == R[argmin([maximum(abs, r) for r in eachrow(R)]), :]

    A = [1 2 3; 2 4 6; 1 0 1]
    @test reduced_row_echelon_form(A) == [1 0 1; 0 1 1; 0 0 0]
    @test row_space_basis(A) == [1 0 1; 0 1 1]
    @test extract_basis(A) == [1 2 3; 1 0 1]
    Q = orthogonal_complement([1//1 2 3])
    @test size(Q, 1) == 2 && iszero([1 2 3] * Q')

    A = [2 4 6; 3 5 7]
    H, V = hermite_normal_form(A)
    @test A * V == H
    @test abs(det(Rational{BigInt}.(V))) == 1
    @test H == [2 0 0; 0 1 0]

    L = lattice_basis_inside_span([1//1 1 0; 0 1 1])
    @test size(L, 1) == 2 && iszero(L * [1, -1, 1])

    # Lattice points of polygons, checked against Pick's theorem
    @test polygon_lattice_count([[0, 0], [3, 0], [3, 2], [0, 2]]) == 12
    @test polygon_lattice_count([[0, 0], [4, 0], [0, 3]]) == 11
    @test polygon_lattice_count([[0, 0], [5, 1], [2, 7]]) == 20
end

@testset "Matroids" begin
    M = Matroid(matroid.uniform_matroid(2, 3))
    @test M.rank == 2
    @test length(nonempty_proper_flats(M)) == 3
    @test join(M, [2, 3]) == 5
    @test meet(M, [2, 3]) == 1
    @test closure(M, Set([0])) == 2
    @test join(["a", "b"], ",") == "a,b"
    f5 = Symbolics.variable(:f, 5)
    @test snapper_polynomial(M) ≂ 1 + f5

    # f_to_h and h_to_f are inverse changes of variables
    M = Matroid(matroid.uniform_matroid(3, 4))
    vol = volume_polynomial(M)
    @test Symbolics.substitute(Symbolics.substitute(vol, h_to_f(M)), f_to_h(M)) ≂ vol
end

@testset "NormalComplexes" begin
    unit_cube = Vector{Rational{Int}}[[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
                                      [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
    @test hypercube_volume(unit_cube) == 1
    # Quadrilateral (0,0), (1,0), (2,1), (0,1), corners listed in hypercube_corners order: area 3/2
    @test hypercube_volume(Vector{Rational{Int}}[[0, 0], [1, 0], [0, 1], [2, 1]]) == 3//2

    h = Symbolics.variables(:h, 2:4)
    M = Matroid(matroid.uniform_matroid(2, 3))
    @test volume_polynomial(M) ≂ sum(h)

    # Lattice points of the box [0, h1] × [0, h2] in the positive quadrant
    R = Matrix{Rational{Int}}(I, 2, 2)
    @variables h1 h2
    p = cone_ehrhart_polynomial(R, R, [h1, h2], [1//1, 1//1])
    @test p ≂ (h1 + 1) * (h2 + 1)
end

end
