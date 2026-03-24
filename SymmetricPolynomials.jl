"""
    SymmetricPolynomials

Express a symmetric polynomial as a polynomial in the elementary symmetric
polynomials e₁, e₂, …, eₙ using the classical leading-term cancellation
algorithm (fundamental theorem of symmetric polynomials).

Dependencies:
    julia> using Pkg; Pkg.add(["Symbolics", "Combinatorics"])
"""
module SymmetricPolynomials

using Symbolics
using Combinatorics

export elem_sym_polys, to_elementary, substitute_elementary


# ── Internal helpers ───────────────────────────────────────────────────────────

"""
    _D(f, x) -> Num

Differentiate Symbolics expression `f` with respect to `x` and expand.
"""
_D(f, x) = expand(Symbolics.derivative(f, x))


"""
    _poly_degree(f, x) -> Int

Return the polynomial degree of `f` in variable `x` via repeated symbolic
differentiation. Returns -1 if `f` is identically zero.
"""
function _poly_degree(f, x)
    g = expand(f)
    iszero(g) && return -1
    d = 0
    while true
        dg = _D(g, x)
        iszero(dg) && return d
        g  = dg
        d += 1
    end
end


"""
    _monomial_coeff(f, vars, exps) -> Num

Extract the coefficient of the monomial ∏ᵢ vars[i]^exps[i] from polynomial `f`.

Uses the formula:
    coeff = (1 / ∏ᵢ exps[i]!) · ∂^(∑exps[i]) f / ∏ᵢ ∂vars[i]^exps[i]  |_{vars=0}
"""
function _monomial_coeff(f, vars, exps)
    g = expand(f)
    for (x, a) in zip(vars, exps)
        a == 0 && continue
        for _ in 1:a
            g = _D(g, x)
        end
        g = expand(g / factorial(a))
        g = expand(substitute(g, Dict(x => 0)))
    end
    return g     # a pure numeric constant at this point
end


"""
    _leading_exponents(f, vars) -> Vector{Int}

Return the exponent vector [a₁, a₂, …, aₙ] of the lexicographically leading
monomial of the symmetric polynomial `f`, with ordering vars[1] > … > vars[n].

For a symmetric polynomial the vector satisfies a₁ ≥ a₂ ≥ … ≥ aₙ ≥ 0.

**Algorithm:**
- a₁ = deg(f, vars[1])
- a₂ = deg(coeff of vars[1]^a₁ in f,  vars[2])
- aᵢ = deg(coeff of vars[1]^a₁…vars[i-1]^a_{i-1},  vars[i])

Each successive coefficient is extracted using the _monomial_coeff formula,
truncated to the running prefix.
"""
function _leading_exponents(f, vars)
    n    = length(vars)
    exps = zeros(Int, n)
    g    = expand(f)

    for i in 1:n
        iszero(g) && break

        aᵢ = _poly_degree(g, vars[i])
        aᵢ <= 0 && break           # degree 0 ⟹ constant (symmetric ⟹ all vars[j] absent)

        exps[i] = aᵢ

        # Reduce g to the coefficient of vars[i]^aᵢ:
        #   g ← (1/aᵢ!) · ∂^aᵢ g / ∂vars[i]^aᵢ
        # After this, g is a polynomial in vars[i+1], …, vars[n] only.
        for _ in 1:aᵢ
            g = _D(g, vars[i])
        end
        g = expand(g / factorial(aᵢ))
    end

    return exps
end


# ── Public API ─────────────────────────────────────────────────────────────────

"""
    elem_sym_polys(vars) -> Vector{Num}

Return the elementary symmetric polynomials [e₁, e₂, …, eₙ] in `vars`:

    eₖ = Σ_{|S|=k} ∏_{i ∈ S} vars[i]

# Example
```julia
@variables x y z
elem_sym_polys([x, y, z])
# [x+y+z, x*y + x*z + y*z, x*y*z]
```
"""
function elem_sym_polys(vars::AbstractVector)
    n = length(vars)
    return [sum(prod(vars[collect(s)]) for s in combinations(1:n, k)) for k in 1:n]
end


"""
    to_elementary(f, vars; max_iter=500, sym=:e) -> Num

Express the symmetric polynomial `f` in variables `vars` as a polynomial in the
elementary symmetric polynomials e₁, …, eₙ.

Returns a Symbolics expression in symbolic proxies `sym[1], …, sym[n]`,
where `sym[k]` represents eₖ = Σ_{|S|=k} ∏_{i∈S} vars[i].

# Arguments
- `f`        – A symmetric polynomial (Symbolics.Num)
- `vars`     – Vector of symbolic variables w.r.t. which `f` is symmetric
- `max_iter` – Safety cap on reduction steps (default 500)
- `sym` – Symbol used for the output proxies (default `:e`)

# Returns
`Num` expression in `sym[1], …, sym[n]`.  To recover `f` in the
original variables, call `substitute_elementary(result, vars; sym)`.

# Algorithm  (leading-term cancellation)
Given the lex-leading monomial c·x^a with a₁ ≥ a₂ ≥ … ≥ aₙ, subtract

    c · e₁^(a₁-a₂) · e₂^(a₂-a₃) · ⋯ · eₙ^(aₙ)

which has exactly the same leading monomial.  Each step strictly reduces
the leading exponent in lex order, guaranteeing termination.

# Examples
```julia
@variables x y z

# Newton's identity: p₂ = e₁² - 2e₂  (default sym)
to_elementary(x^2 + y^2 + z^2, [x, y, z])
# ⟹  e[1]^2 - 2*e[2]

# Custom symbol σ
to_elementary(x^2 + y^2 + z^2, [x, y, z], sym = :σ)
# ⟹  σ[1]^2 - 2*σ[2]

# e₃ directly
to_elementary(x*y*z, [x, y, z], sym = :p)
# ⟹  p[3]

# Schur polynomial s_{(2,1)}
to_elementary(expand(x^2*y + x^2*z + y^2*x + y^2*z + z^2*x + z^2*y),
                         [x, y, z])
# ⟹  e[1]*e[2] - 3*e[3]
```
"""
function to_elementary(f, vars::AbstractVector;
                                  max_iter::Int    = 500,
                                  sym::Symbol = :e)
    n = length(vars)

    e_polys = elem_sym_polys(vars)                        # actual eₖ as polynomials in vars
    E = Symbolics.variables(sym, 1:n)          # symbolic stand-ins for the output

    result = zero(Num)
    g      = expand(f)

    for _ in 1:max_iter
        g = expand(g)
        iszero(g) && break

        a = _leading_exponents(g, vars)

        # ── Constant remainder ─────────────────────────────────────────────
        # all(iszero, a)  ⟺  g does not depend on any variable  ⟺  g is a constant
        if all(iszero, a)
            result = expand(result + g)
            g      = zero(Num)
            break
        end

        # ── Leading coefficient ────────────────────────────────────────────
        c = _monomial_coeff(g, vars, a)

        if iszero(c)
            error("Encountered zero leading coefficient for non-zero exponent $a. " *
                  "Is `f` genuinely symmetric in all of `vars`?")
        end

        # ── Build  e₁^(a₁-a₂) · e₂^(a₂-a₃) · ⋯ · eₙ^(aₙ) ───────────────
        e_term = one(Num)   # in original variables
        E_term = one(Num)   # in symbolic e[k] proxies

        for i in 1:n
            dᵢ = a[i] - (i < n ? a[i+1] : 0)   # successive differences
            if dᵢ > 0
                e_term = expand(e_term * e_polys[i]^dᵢ)
                E_term = E_term * E[i]^dᵢ
            end
        end

        result = expand(result + c * E_term)
        g      = expand(g - c * e_term)
    end

    iszero(expand(g)) ||
        error("to_elementary: remainder $g ≠ 0 after $max_iter steps.\n" *
              "Verify that `f` is symmetric in every variable in `vars`.")

    return result
end


"""
    substitute_elementary(expr, vars; sym=:e) -> Num

Replace the symbolic proxies `sym[k]` in `expr` (output of
`to_elementary`) with the actual elementary symmetric polynomials
of `vars`, recovering the original polynomial.

`sym` must match the one passed to `to_elementary`.

Useful for round-trip verification.
"""
function substitute_elementary(expr, vars::AbstractVector; sym::Symbol = :e)
    n       = length(vars)
    e_polys = elem_sym_polys(vars)
    E       = Symbolics.variables(sym, 1:n)
    sub     = Dict(E[k] => e_polys[k] for k in 1:n)
    return expand(substitute(expr, sub))
end


# ── Demo ───────────────────────────────────────────────────────────────────────

function demo()
    @variables x y z

    tests = [
        ("e₁ = x+y+z",
         x + y + z),

        ("e₂ = xy+xz+yz",
         x*y + x*z + y*z),

        ("p₂ = x²+y²+z²  (Newton's identity: e₁²−2e₂)",
         x^2 + y^2 + z^2),

        ("p₃ = x³+y³+z³  (Newton: e₁³−3e₁e₂+3e₃)",
         x^3 + y^3 + z^3),

        ("e₃ = xyz",
         x*y*z),

        ("Schur s_{(2,1)}: Σ x²y  (= e₁e₂−3e₃)",
         expand(x^2*y + x^2*z + y^2*x + y^2*z + z^2*x + z^2*y)),

        ("(x+y+z)⁴ expanded  (= e₁⁴)",
         expand((x + y + z)^4)),
    ]

    println("═"^65)
    println("  Symmetric polynomial → elementary symmetric polynomials")
    println("  Variables: x, y, z")
    println("═"^65)

    for (label, f) in tests
        expr = to_elementary(f, [x, y, z])
        # Verify round-trip
        recovered = substitute_elementary(expr, [x, y, z])
        ok = isequal(expand(recovered - f), 0) ? "✓" : "✗ MISMATCH"
        println("\n  f = $label")
        println("    = $expr  [$ok]")
    end

    println("\n" * "═"^65)
    println("  Two-variable example")
    println("═"^65)
    @variables a b
    f2 = a^3 + b^3
    expr2 = to_elementary(f2, [a, b])
    println("\n  f = a³+b³  =  $expr2")
    println("  (= e[1]³ - 3·e[1]·e[2]  ✓)")

    println("\n" * "═"^65)
    println("  Custom sym examples")
    println("═"^65)

    f3 = x^2 + y^2 + z^2

    for sym in [:e, :σ, :p, :h]
        expr = to_elementary(f3, [x, y, z]; sym = sym)
        recovered = substitute_elementary(expr, [x, y, z]; sym = sym)
        ok = isequal(expand(recovered - f3), 0) ? "✓" : "✗ MISMATCH"
        println("\n  sym = :$sym  →  $expr  [$ok]")
    end
end

end
