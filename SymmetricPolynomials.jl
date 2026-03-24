# SymmetricPolynomials
#
# Express a symmetric polynomial as a polynomial in the elementary symmetric
# polynomials e₁, e₂, …, eₙ using the leading-term cancellation algorithm.

using Symbolics
using Combinatorics

# ── Internal helpers ───────────────────────────────────────────────────────────

_D(f::Num, x::Num) = expand(Symbolics.derivative(f, x))

# Degree of f in variable x via repeated differentiation
function _poly_degree(f::Num, x::Num)
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

# Coefficient of ∏ᵢ vars[i]^exps[i] in f
function _monomial_coeff(f::Num, vars::AbstractVector{Num}, exps::Vector{Int})
    g = expand(f)
    for (x, a) in zip(vars, exps)
        a == 0 && continue
        for _ in 1:a
            g = _D(g, x)
        end
        g = expand(g / factorial(a))
        g = expand(substitute(g, Dict(x => 0)))
    end
    return g
end

# Exponent vector of the lex-leading monomial of symmetric polynomial f
function _leading_exponents(f::Num, vars::AbstractVector{Num})
    n    = length(vars)
    exps = zeros(Int, n)
    g    = expand(f)
    for i in 1:n
        iszero(g) && break
        aᵢ = _poly_degree(g, vars[i])
        aᵢ <= 0 && break
        exps[i] = aᵢ
        for _ in 1:aᵢ
            g = _D(g, vars[i])
        end
        g = expand(g / factorial(aᵢ))
    end
    return exps
end

# ── Public API ─────────────────────────────────────────────────────────────────

# Elementary symmetric polynomials e₁, …, eₙ in vars
function elem_sym_polys(vars::AbstractVector{Num})
    n = length(vars)
    return [sum(prod(vars[collect(s)]) for s in combinations(1:n, k)) for k in 1:n]
end

# Express symmetric polynomial f in vars as a polynomial in e₁, …, eₙ.
# Output uses symbolic variables sym_1, …, sym_n (default sym=:e).
function to_elementary(f::Num, vars::AbstractVector{Num}; max_iter::Int=500, sym::Symbol=:e)
    n       = length(vars)
    e_polys = elem_sym_polys(vars)
    E       = Symbolics.variables(sym, 1:n)
    result  = zero(Num)
    g       = expand(f)

    for _ in 1:max_iter
        g = expand(g)
        iszero(g) && break

        a = _leading_exponents(g, vars)

        if all(iszero, a)
            result = expand(result + g)
            g      = zero(Num)
            break
        end

        c = _monomial_coeff(g, vars, a)
        iszero(c) && error("zero leading coefficient for exponent $a — is f symmetric?")

        # Subtract c · e₁^(a₁-a₂) · e₂^(a₂-a₃) · ⋯ · eₙ^(aₙ)
        e_term = one(Num)
        E_term = one(Num)
        for i in 1:n
            dᵢ = a[i] - (i < n ? a[i+1] : 0)
            if dᵢ > 0
                e_term = expand(e_term * e_polys[i]^dᵢ)
                E_term = E_term * E[i]^dᵢ
            end
        end

        result = expand(result + c * E_term)
        g      = expand(g - c * e_term)
    end

    iszero(expand(g)) ||
        error("remainder $g ≠ 0 after $max_iter steps — is f symmetric in all vars?")

    return result
end

# Substitute e₁, …, eₙ back into expr to recover the original polynomial
function substitute_elementary(expr::Num, vars::AbstractVector{Num}; sym::Symbol=:e)
    n       = length(vars)
    e_polys = elem_sym_polys(vars)
    E       = Symbolics.variables(sym, 1:n)
    return expand(substitute(expr, Dict(E[k] => e_polys[k] for k in 1:n)))
end
