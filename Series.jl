module Series

using Symbolics
using Latexify

export GenericSeries, FiniteSeries, PowerSeries, PolynomialSeries, BinomialSeries
export get_coeff

import Base: +, /, *, ^, ∘, -, <<, show

abstract type AbstractSeries end
# expects get_coeff(series, i) to be defined

function first_coeffs(series::AbstractSeries, n)
	return [get_coeff(series, i) for i in 0:n]
end

function polyn_head(series::AbstractSeries, d)
	return sum(first_coeffs(series, d) .* [Symbolics.variable(:t)^i for i in 0:d]; init=0)
end

function print_series(io::IO, s::AbstractSeries, deg::Int)
	print(io, "\$\$" * latexify(polyn_head(s, deg)) * "+ o(t^{$deg}) \$\$")
end

function show(io::IO, m::MIME"text/latex", s::AbstractSeries)
	print_series(io, s, 5)
end

# Generic series

struct GenericSeries <: AbstractSeries
	symbol::Symbol
end

function get_coeff(series::GenericSeries, deg)
	if (deg == 0) return 1 end
	return Symbolics.variable(series.symbol, deg)
end

# Finite series

struct FiniteSeries <: AbstractSeries
	coeffs::Vector
end

function get_coeff(series::FiniteSeries, index)
	if index + 1 > length(series.coeffs)
		return 0
	else
		return series.coeffs[index + 1]
	end
end

# Polynomial series

struct PolynomialSeries <: AbstractSeries
	coeffs
	var
	PolynomialSeries(polyn, var) = new(Symbolics.polynomial_coeffs(polyn, [var])[1], var)
end

function get_coeff(s::PolynomialSeries, i)
	return get(s.coeffs, s.var^i, 0)
end

# Shift left

struct ShiftLeftSeries <: AbstractSeries
	series::AbstractSeries
	shift::Int
end

function get_coeff(s::ShiftLeftSeries, i)
	return get_coeff(s.series, i + s.shift)
end

function <<(s::AbstractSeries, i::Int)
	return ShiftLeftSeries(s, i)
end


# Additive inverse series

struct MinusSeries <: AbstractSeries
	orig::AbstractSeries
end

function get_coeff(s::MinusSeries, index)
	return -get_coeff(s.orig, index)
end

function -(s::AbstractSeries)
	return MinusSeries(s)
end

# Addition of two series

struct AddSeries <: AbstractSeries
	a::AbstractSeries
	b::AbstractSeries
end

function get_coeff(s::AddSeries, index)
	return get_coeff(s.a, index) + get_coeff(s.b, index)
end

function +(a::AbstractSeries, b::AbstractSeries)
	return AddSeries(a, b)
end

function -(a::AbstractSeries, b::AbstractSeries)
	return AddSeries(a, MinusSeries(b))
end

# Geometric series

struct PowerSeries <: AbstractSeries
	a
end

function get_coeff(series::PowerSeries, deg)
	return series.a^deg
end

# Element-wise multiplication

struct ElemMultSeries <: AbstractSeries
	a::AbstractSeries
	b::AbstractSeries
end

function get_coeff(series::ElemMultSeries, deg)
	return get_coeff(series.a, deg) * get_coeff(series.b, deg)
end

function ∘(a::AbstractSeries, b::AbstractSeries)
	return ElemMultSeries(a, b)
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
	return binomial(series.power, deg)# * series.a^(series.power - deg) * series.b^deg 
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

# Inverse series

struct InverseSeries <: CacheSeries
	orig::AbstractSeries
	cache::Vector
	InverseSeries(orig) = new(orig, [])
end

function compute_coeff(series::InverseSeries, index)
	if index == 0
		return 1//get_coeff(series.orig, 0)
	end
	d = index - 1
	p = polyn_head(series.orig, d + 1)
	q = polyn_head(series, d)
	r = expand(p * q)
	return -expand(Symbolics.coeff(r, Symbolics.variable(:t)^index)//get_coeff(series.orig, 0))
end

# Multiplication series

struct MultSeries <: CacheSeries
	a::AbstractSeries
	b::AbstractSeries
	cache::Vector
	MultSeries(a, b) = new(a, b, [])
end

function compute_coeff(series::MultSeries, d)
	p = polyn_head(series.a, d)
	q = polyn_head(series.b, d)
	r = expand(p*q)
	return expand(Symbolics.coeff(r, Symbolics.variable(:t)^d))
end

function *(a::AbstractSeries, b::AbstractSeries)
	return MultSeries(a, b)
end

function /(a::AbstractSeries, b::AbstractSeries)
	return MultSeries(a, InverseSeries(b))
end

function ^(a::AbstractSeries, e::Integer)
	if e == 0
		return FiniteSeries([1])
	elseif e > 0
		return a * (a ^ (e - 1))
	else
		inv = InverseSeries(a)
		return inv ^ (-e)
	end
end

end
