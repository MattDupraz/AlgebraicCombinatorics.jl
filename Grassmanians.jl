using Symbolics

include("SymmetricPolynomials.jl")

include("Series.jl")
using .Series


ToddSeries() = inv((1 - ExpSeries(-1)) << 1)
ToddSeries(w) = ToddSeries() ⊙  GeometricSeries(w)


function todd_class(k, n)
    A = Symbolics.variables(:a, 1:k) # Chern roots of the tautological line bundle
    B = Symbolics.variables(:b, 1:(n - k)) # Chern roots of the quotient bundle
    
    T = SymmetricSeries(B,
                        SymmetricSeries(A, prod([ToddSeries(-a + b) for a in A, b in B]), :s),
                        :q)
    return T
end
