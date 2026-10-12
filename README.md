# AlgebraicCombinatorics

Julia tools for symbolic power series, symmetric polynomials, exact lattice
computations, matroids and their normal complexes (volume and Ehrhart polynomials).

| File | Contents |
|---|---|
| `src/Series.jl` | Formal power series, with named series (`BinomialSeries`, `ExpSeries`, …) |
| `src/SymmetricPolynomials.jl` | Rewriting symmetric polynomials in elementary symmetric polynomials |
| `src/Pochhammer.jl` | Falling and rising factorials |
| `src/Grassmannians.jl` | Todd classes of bundles and Grassmannians |
| `src/LinAlgTools.jl` | Exact linear algebra, LLL, Hermite normal form, Barvinok decomposition, lattice point counts |
| `src/Matroids.jl` | Matroids via polymake, flats, Bergman fan rays, Snapper polynomial |
| `src/NormalComplexes.jl` | Volume and Ehrhart polynomials of normal complexes |
| `src/Cache.jl` | `@cache` macro storing results in `.cache/` |

## Installation

```julia
pkg> add git@git.imp.fu-berlin.de:matthed69/AlgebraicCombinatorics.jl.git
```

To edit the code, use `dev` instead of `add`, or clone the repository and run
`pkg> dev path/to/AlgebraicCombinatorics.jl`.

## Usage

```julia
using AlgebraicCombinatorics
using Polymake: matroid

M = Matroid(matroid.uniform_matroid(3, 4))
volume_polynomial(M)
```

The Pluto notebooks are in `notebooks/`. `ehrhart.jl` and `lambda_rings.jl`
activate this project in their first cell, so they work from a clone of the
repository. `notebook_1.jl` to `notebook_6.jl` are archived early experiments
and self-contained. `examples/` contains standalone scripts, run with
`julia --project=. examples/<script>.jl`.

## Tests

```julia
pkg> activate .
pkg> test
```
