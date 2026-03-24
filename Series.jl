module Series

using Symbolics
using ..SymmetricPolynomials
import ..SymmetricPolynomials: to_elementary
import Symbolics: substitute

export Series, cached
export GeometricSeries, BinomialSeries, ExpSeries, LogSeries
export substitute, to_elementary
export first_coeffs, first_nonzero_index, ⊙

import Base: +, /, *, ^, -, <<, >>, show, getindex, inv, map

#####################
# SERIES TYPE
#####################

# A series is a function ℕ → value

struct Series
    f::Function
    Series(f::Function)            = new(f)
    Series(coeffs::AbstractVector) = new(i -> (i < 0 || i + 1 > length(coeffs)) ? 0 : coeffs[i + 1])
    Series(coeffs::Number...)      = new(i -> (i < 0 || i + 1 > length(coeffs)) ? 0 : coeffs[i + 1])
end

getindex(s::Series, i) = s.f(i)

function first_coeffs(s::Series, n)
    return [s[i] for i in 0:n]
end

function first_nonzero_index(s::Series; limit=1000)
    for i in 0:limit
        s[i] != 0 && return i
    end
    return nothing
end

# Returns a memoized version of a series or function

function cached(f::Function)
    cache = []
    Series(function(i)
        while length(cache) < i + 1
            push!(cache, f(length(cache)))
        end
        return cache[i + 1]
    end)
end

cached(s::Series) = cached(i -> s[i])

######################
# NAMED SERIES
######################

# Geometric series: elem^i

GeometricSeries(elem) = Series(i -> elem^i)

# Binomial series: coefficients of (1 + t)^x

BinomialSeries(x) = Series(i -> i == 0 ? oneunit(x)//1 : prod(j -> (x-j+1)//j, 1:i))

# Exponential series: exp(t) around t=0

ExpSeries(t=1) = Series(i -> t^i // factorial(i))

# Log series: log(1 + t) around t=0

LogSeries(t=1) = Series(function(deg)
    if deg == 0
        return 0
    elseif deg % 2 == 1
        return (t^deg)//deg
    else
        return -(t^deg)//deg
    end
end)

####################
# OPERATIONS
####################

# Applies a function to each coefficient

map(f, s::Series) = Series(i -> f(s[i]))

# Applies a substitution dictionary to each coefficient

substitute(s::Series, dict) = map(x -> expand(substitute(x, dict)), s)

# Expresses each coefficient as a polynomial in the elementary symmetric
# polynomials of vars, using sym for the output symbols

to_elementary(s::Series, vars; sym=:e) =
    map(x -> to_elementary(x, vars; sym=sym), s)

# Shift left

<<(s::Series, n::Int) = Series(i -> s[i + n])

# Shift right

>>(s::Series, n::Int) = Series(i -> s[i - n])

# Additive inverse

-(s::Series) = Series(i -> -s[i])

# Addition

+(a::Series, b::Series) = Series(i -> a[i] + b[i])
-(a::Series, b::Series) = Series(i -> a[i] - b[i])

# Hadamard product (element-wise multiplication)

⊙(a::Series, b::Series) = Series(i -> a[i] * b[i])

# Multiplication by a constant

*(a::Number, b::Series) = Series(i -> a * b[i])
*(a::Series, b::Number) = Series(i -> a[i] * b)

# Cauchy product of two series

*(a::Series, b::Series) = Series(d -> expand(sum(a[k] * b[d-k] for k in 0:d)))

# Division of series (requires invertible constant term)

function /(a::Series, b::Series)
    iszero(b[0]) && error("constant coefficient of divisor is not invertible")
    ref = Ref{Series}()
    ref[] = cached(function(d)
        d == 0 && return a[0] / b[0]
        s = sum(ref[][k] * b[d-k] for k in 0:d-1; init=zero(Num))
        return expand((a[d] - s) / b[0])
    end)
    return ref[]
end

inv(a::Series) = 1 / a

# Exponentiation

function ^(a::Series, e::Integer)
    if e == 0
        return Series(1)
    elseif e > 0
        return a * (a ^ (e - 1))
    else
        return 1 / a^(-e)
    end
end

# Promotion rules

/(a::Number, b::Series) = Series(a) / b
/(a::Series, b::Number) = a * inv(b)
+(a::Number, b::Series) = Series(a) + b
+(a::Series, b::Number) = a + Series(b)
-(a::Number, b::Series) = Series(a) - b
-(a::Series, b::Number) = a - Series(b)

end
