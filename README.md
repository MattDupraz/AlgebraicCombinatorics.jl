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

## References

The volume computation of normal complexes follows [NR], the Snapper
polynomial of a matroid comes from [LLPP], and the formula expressing the
Snapper polynomial in terms of the volume polynomial is from [DGM]. The lattice
point counting is based on the Barvinok decomposition as described in [DHTY].

- [NR] A. Nathanson, D. Ross, *Tropical fans and normal complexes: putting the
  "volume" back in "volume polynomials"*, Adv. Math. **420** (2023), 108981.
  [doi:10.1016/j.aim.2023.108981](https://doi.org/10.1016/j.aim.2023.108981),
  [arXiv:2110.08647](https://arxiv.org/abs/2110.08647),
  [Zbl 1536.14059](https://zbmath.org/?q=an:1536.14059)
- [LLPP] M. Larson, S. Li, S. Payne, N. Proudfoot, *K-rings of wonderful
  varieties and matroids*, Adv. Math. **441** (2024), 109554.
  [doi:10.1016/j.aim.2024.109554](https://doi.org/10.1016/j.aim.2024.109554),
  [arXiv:2210.03169](https://arxiv.org/abs/2210.03169),
  [Zbl 1535.05266](https://zbmath.org/?q=an:1535.05266)
- [DGM] M. Dupraz, A. Gross, L. Monin, *A Riemann-Roch theorem for Frobenius
  quotients*, preprint (2026).
  [arXiv:2607.26176](https://arxiv.org/abs/2607.26176)
- [DHTY] J. A. De Loera, R. Hemmecke, J. Tauzer, R. Yoshida, *Effective lattice
  point counting in rational convex polytopes*, J. Symbolic Comput. **38**
  (2004), no. 4, 1273–1302.
  [doi:10.1016/j.jsc.2003.04.003](https://doi.org/10.1016/j.jsc.2003.04.003),
  [Zbl 1137.52303](https://zbmath.org/?q=an:1137.52303)

## License

MIT, see [`LICENSE`](LICENSE).
