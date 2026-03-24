using Symbolics

include("SymmetricPolynomials.jl")
include("Series.jl")


ToddSeries() = inv((1 - ExpSeries(-1)) << 1)
ToddSeries(w) = ToddSeries() ⊙  GeometricSeries(w)

function todd_class(d)
    A = Symbolics.variables(:a, 1:d)

    T = to_elementary(prod([ToddSeries(a) for a in A]), A; sym=:c)
end

function todd_class_grassmanian(k, n)
    A = Symbolics.variables(:a, 1:k) # Chern roots of the tautological line bundle
    B = Symbolics.variables(:b, 1:(n - k)) # Chern roots of the quotient bundle
    
    T = to_elementary(
            to_elementary(prod([ToddSeries(-a + b) for a in A, b in B]), A; sym=:s),
            B; sym=:q)
    return T
end

# Compute coeffs for Grassmanian(2, 4)
k = 2
n = 4
T = todd_class_grassmanian(k, n)

# Since [S ⊗ Q] = [O], we know that the product of the chern classes is 1 in the Chow ring
S_inv = inv(Series([1; Symbolics.variables(:s, 1:2)]))
q_to_s = Dict([Symbolics.variable(:q, i) => S_inv[i] for i in 1:2])

T = substitute(T, q_to_s)

first_coeffs(T, k*(n - k))
