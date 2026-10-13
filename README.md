# NormalComplexes

Julia code accompanying my research on matroids and their normal complexes,
in particular their volume, Ehrhart and Snapper polynomials.

This is research code rather than a polished library: the interface may change
without notice.

Most of this code was originally a collection of scripts used in my research, which
was later reorganized into a package with the help of an LLM. 

| File | Contents |
|---|---|
| `src/series.jl` | Formal power series, with named series (`BinomialSeries`, `ExpSeries`, …) |
| `src/symmetric_polynomials.jl` | Rewriting symmetric polynomials in elementary symmetric polynomials |
| `src/pochhammer.jl` | Falling and rising factorials |
| `src/grassmannians.jl` | Todd classes of bundles and Grassmannians |
| `src/linalg_tools.jl` | Exact linear algebra, LLL, Hermite normal form, Barvinok decomposition, lattice point counts |
| `src/matroids.jl` | Matroids via polymake, flats, Bergman fan rays, Snapper polynomial |
| `src/normal_complexes.jl` | Volume and Ehrhart polynomials of normal complexes |
| `src/cache.jl` | `@cache` macro storing results in `.cache/` |

## Installation

```julia
pkg> add https://github.com/MattDupraz/NormalComplexes.jl
```

To edit the code, use `dev` instead of `add`, or clone the repository and run
`pkg> dev path/to/NormalComplexes.jl`.

## Usage

```julia
using NormalComplexes
using Polymake: matroid

M = Matroid(matroid.uniform_matroid(3, 4))
volume_polynomial(M)
```

The Pluto notebooks are in `notebooks/`. They activate this project in their
first cell, so they work from a clone of the repository. `examples/` contains
standalone scripts, run with `julia --project=. examples/<script>.jl`.

## Tests

```julia
pkg> activate .
pkg> test
```

## License

MIT, see [`LICENSE`](LICENSE).
