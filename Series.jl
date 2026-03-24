module Series

using Symbolics
using ..SymmetricPolynomials

export FiniteSeries, FuncSeries
export GeometricSeries, BinomialSeries, ExpSeries, LogSeries
export SubstituteSeries, SymmetricSeries
export get_coeff, first_coeffs, first_nonzero_index, ⊙

import Base: +, /, *, ^, -, <<, >>, show, getindex, inv, map

#####################
# ABSTRACT TYPE
#####################

abstract type AbstractSeries end
# expects get_coeff(series, i) to be defined

function getindex(series::AbstractSeries, n)
    return get_coeff(series, n)
end

function first_coeffs(series::AbstractSeries, n)
    return [series[i] for i in 0:n]
end

function first_nonzero_index(series::AbstractSeries; limit=1000)
    for i in 0:limit
        if series[i] != 0
            return i
        end
    end
    return nothing
end

######################
# CONCRETE TYPES
######################

# A series defined by a function ℕ → value

struct FuncSeries <: AbstractSeries
    f
end

function get_coeff(s::FuncSeries, index)
    return s.f(index)
end

# A series defined by a function ℕ → value, with memoization

struct CachedSeries <: AbstractSeries
    cache::Vector
    f
    CachedSeries(f) = new([], f)
end

function get_coeff(s::CachedSeries, index)
    while length(s.cache) < index + 1
        push!(s.cache, s.f(length(s.cache)))
    end
    return s.cache[index + 1]
end

# A series defined by a finite list of coefficients

struct FiniteSeries <: AbstractSeries
    coeffs::Vector
end

function get_coeff(series::FiniteSeries, index)
    if index < 0 || index + 1 > length(series.coeffs)
        return 0
    else
        return series.coeffs[index + 1]
    end
end

######################
# BASIC SERIES
######################

# Geometric series: elem^i

GeometricSeries(elem) = FuncSeries(i -> elem^i)

# Binomial series: coefficients of (1 + t)^x

BinomialSeries(x) = FuncSeries(i -> i == 0 ? oneunit(x)//1 : prod(j -> (x-j+1)//j, 1:i))

# Exponential series: exp(t) around t=0

ExpSeries(t=1) = FuncSeries(i -> t^i // factorial(i))

# Log series: log(1 + t) around t=0

LogSeries(t=1) = FuncSeries(function(deg)
    if deg % 2 == 0
        return (t^deg)//deg
    else
        return -(t^deg)//deg
    end
end)

# Applies a function to each coefficient of a series

map(f, s::AbstractSeries) = FuncSeries(i -> f(s[i]))

# Applies a substitution dictionary to each coefficient

SubstituteSeries(series, dict) = map(x -> expand(substitute(x, dict)), series)

# Expresses each coefficient as a polynomial in the elementary symmetric
# polynomials of vars, using sym_name for the output symbols

SymmetricSeries(vars, series, sym_name=:e) =
    map(x -> symmetric_to_elementary(x, vars; sym_name=sym_name), series)

####################
# OPERATIONS
####################

# Shift left

<<(s::AbstractSeries, n::Int) = FuncSeries(i -> s[i + n])

# Shift right

>>(s::AbstractSeries, n::Int) = FuncSeries(i -> s[i - n])

# Additive inverse

-(s::AbstractSeries) = FuncSeries(i -> -s[i])

# Addition

+(a::AbstractSeries, b::AbstractSeries) = FuncSeries(i -> a[i] + b[i])
-(a::AbstractSeries, b::AbstractSeries) = FuncSeries(i -> a[i] - b[i])

# Hadamard product (element-wise multiplication)

⊙(a::AbstractSeries, b::AbstractSeries) = FuncSeries(i -> a[i] * b[i])

# Multiplication by a constant

*(a::T, b::AbstractSeries) where T <: Number = FuncSeries(i -> a * b[i])
*(a::AbstractSeries, b::T) where T <: Number = FuncSeries(i -> a[i] * b)

# Cauchy product of two series (requires caching)

function *(a::AbstractSeries, b::AbstractSeries)
    CachedSeries(d -> expand(sum(a[k] * b[d-k] for k in 0:d)))
end

# Division of series (requires invertible constant term and caching)

function /(a::AbstractSeries, b::AbstractSeries)
    iszero(b[0]) && error("constant coefficient of divisor is not invertible")
    ref = Ref{CachedSeries}()
    ref[] = CachedSeries(function(d)
        d == 0 && return a[0] / b[0]
        s = sum(ref[][k] * b[d-k] for k in 0:d-1; init=zero(Num))
        return expand((a[d] - s) / b[0])
    end)
    return ref[]
end

function inv(a::AbstractSeries)
    return 1 / a
end

# Exponentiation of series

function ^(a::AbstractSeries, e::Integer)
    if e == 0
        return FiniteSeries([1])
    elseif e > 0
        return a * (a ^ (e - 1))
    else
        return 1 / a^(-e)
    end
end

# Promotion rules

/(a::T, b::AbstractSeries) where T <: Number = FiniteSeries([a]) / b
/(a::AbstractSeries, b::T) where T <: Number = a * inv(b)
+(a::T, b::AbstractSeries) where T <: Number = FiniteSeries([a]) + b
+(a::AbstractSeries, b::T) where T <: Number = a + FiniteSeries([b])
-(a::T, b::AbstractSeries) where T <: Number = FiniteSeries([a]) - b
-(a::AbstractSeries, b::T) where T <: Number = a - FiniteSeries([b])

end
