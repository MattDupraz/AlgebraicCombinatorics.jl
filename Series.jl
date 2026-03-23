module Series

using Symbolics

export FiniteSeries, GeometricSeries, BinomialSeries, ExpSeries, LogSeries
export get_coeff, first_coeffs, first_nonzero_index, ⊙

import Base: +, /, *, ^, -, <<, >>, show, getindex, inv

#####################
# ABSTRACT TYPES
#####################

abstract type AbstractSeries end
# expects get_coeff(series, i) to be defined

# Overloads [] operator
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

# Cached series

abstract type CacheSeries <: AbstractSeries end
# expects compute_coeff(series, i) to be defined

function get_coeff(series::CacheSeries, index)
	while length(series.cache) < index + 1
		push!(series.cache, compute_coeff(series, length(series.cache)))
	end
	return series.cache[index + 1]
end

######################
# BASIC SERIES
######################

# Finite series

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

# Geometric series

struct GeometricSeries <: AbstractSeries
	elem
end

function get_coeff(series::GeometricSeries, deg)
	return series.elem^deg
end

# Binomial series
# series (1 + t)^power

struct BinomialSeries <: AbstractSeries
	power
end

function binomial(x::Number, k::Integer)
    k < 0 && throw(DomainError("k=$k must be non-negative"))
    k == 0 && return oneunit(x)//one(k)
    return prod(i -> (x-i+1)//i, 1:k)
end

function get_coeff(series::BinomialSeries, deg)
	return binomial(series.power, deg)
end

# Exponential series
# series exp(t) around t=0

struct ExpSeries <: AbstractSeries
end

function get_coeff(series::ExpSeries, deg)
	return 1//factorial(deg)
end

# Log series
# series of log(1 + t) around t = 0

struct LogSeries <: AbstractSeries
end

function get_coeff(series::LogSeries, deg)
	if deg % 2 == 0
		return 1//deg
	else
		return -1//deg
	end
end

####################
# OPERATIONS
####################


# Shift left

struct ShiftLeft <: AbstractSeries
	series::AbstractSeries
	shift::Int
end

function get_coeff(s::ShiftLeft, i)
    return s.series[i + s.shift]
end

function <<(s::AbstractSeries, i::Int)
	return ShiftLeft(s, i)
end

# Shift right

struct ShiftRight <: AbstractSeries
    series::AbstractSeries
    shift::Int
end

function get_coeff(s::ShiftRight, i)
    return s.series[i - s.shift]
end

function >>(s::AbstractSeries, i::Int)
	return ShiftRight(s, i)
end

# Additive inverse series

struct Negative <: AbstractSeries
	orig::AbstractSeries
end

function get_coeff(s::Negative, index)
    return -s.orig[index]
end

function -(s::AbstractSeries)
	return Negative(s)
end

# Addition of two series

struct Add <: AbstractSeries
	a::AbstractSeries
	b::AbstractSeries
end

function get_coeff(s::Add, index)
    return s.a[index] + s.b[index]
end

function +(a::AbstractSeries, b::AbstractSeries)
	return Add(a, b)
end

function -(a::AbstractSeries, b::AbstractSeries)
	return a + (-b)
end

# Multiplication series

struct Mult <: CacheSeries
	cache::Vector

	a::AbstractSeries
	b::AbstractSeries
    Mult(a, b) = new([], a, b)
end

# Open question: does it actually make sense to cache computations here?

function compute_coeff(series::Mult, d)
    result = 0
    for k in 0:d
        result += series.a[k] * series.b[d-k]
    end
    return expand(result)
end

function *(a::AbstractSeries, b::AbstractSeries)
	return Mult(a, b)
end

# Multiplication by a constant

struct MultConst <: AbstractSeries
    a::AbstractSeries
    b
end

function get_coeff(series::MultConst, index)
    return series.a[index] * series.b
end

*(a::T, b::AbstractSeries) where T <: Number = MultConst(b, a)
*(a::AbstractSeries, b::T) where T <: Number = MultConst(a, b)


# Hadamard product (element-wise multiplication)

struct ElemMult <: AbstractSeries
	a::AbstractSeries
	b::AbstractSeries
end

function get_coeff(series::ElemMult, index)
    return series.a[index] * series.b[index]
end

function ⊙(a::AbstractSeries, b::AbstractSeries)
	return ElemMult(a, b)
end

# Division of series (requires invertible constant coefficient)

struct Div <: CacheSeries
    cache::Vector

    a::AbstractSeries
    b::AbstractSeries
    Div(a, b) = new([], a, b)
end

function compute_coeff(series::Div, d)
    if iszero(series.b[0])
        error("constant coefficient of divisor is not invertible")
    end

    if d < 0
        return 0
    elseif d == 0
        return series.a[0] / series.b[0]
    else
        s = 0
        for k in 0:d-1
            s += series[k] * series.b[d-k]
        end
        return expand((series.a[d] - s) / series.b[0])
    end
end

function /(a::AbstractSeries, b::AbstractSeries)
	return Div(a, b)
end

function inv(a::AbstractSeries)
    return 1 / a
end

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
