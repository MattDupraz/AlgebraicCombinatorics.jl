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

function pochhammer(expr, vars)
	coeffs = Sym.polynomial_coeffs(expr, vars)[1]
	result = 0
	for (coeff, monom) in coeffs
		poch_monom = 1
		for var in vars
			poch_monom *= falling_factorial(var, Sym.degree(monom, var))
		end
		result += coeff * poch_monom
	end
	return Sym.simplify(result; expand=true)
end

function vandermonde(vars)
	result = 1
	for i in 1:length(x)
		for j in i+1:length(x)
			result *= x[i] - x[j]
		end
	end
	
	return result
end
