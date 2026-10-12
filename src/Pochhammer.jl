using Symbolics

function falling_factorial(x, n::Int)
    if n < 0
        error("Falling factorial exponent must be non-negative.")
    elseif n == 0
        return 1
    elseif n == 1
        return x
    else
        term = x
        for i in 1:(n-1)
            term *= (x - i)
        end
        return term
    end
end

function rising_factorial(x, n::Int)
    if n < 0
        error("Rising factorial exponent must be non-negative.")
    elseif n == 0
        return 1
    elseif n == 1
        return x
    else
        term = x
        for i in 1:(n-1)
            term *= (x + i)
        end
        return term
    end
end

function to_falling_factorials(expr, vars)
	coeffs = Symbolics.polynomial_coeffs(expr, vars)[1]
	result = 0
	for (monom, coeff) in coeffs
		poch_monom = 1
		for var in vars
			poch_monom *= falling_factorial(var, Symbolics.degree(monom, var))
        end
        #println("coeff: ", coeff)
        #println("monom: ", monom)
        #println("poch_monom: ", poch_monom)

		result += coeff * poch_monom
	end
	return Symbolics.expand(Num(result))
end

function to_rising_factorials(expr, vars)
	coeffs = Symbolics.polynomial_coeffs(expr, vars)[1]
	result = 0
	for (monom, coeff) in coeffs
		poch_monom = 1
		for var in vars
			poch_monom *= rising_factorial(var, Symbolics.degree(monom, var))
        end
        #println("coeff: ", coeff)
        #println("monom: ", monom)
        #println("poch_monom: ", poch_monom)

		result += coeff * poch_monom
	end
	return Symbolics.expand(Num(result))
end

function vandermonde(vars)
	result = 1
	for i in 1:length(vars)
		for j in i+1:length(vars)
			result *= vars[i] - vars[j]
		end
	end
	
	return result
end
