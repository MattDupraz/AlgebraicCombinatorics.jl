using Symbolics

# Todd classes, expressed in Chern classes via symmetric polynomials

# Todd series t / (1 - exp(-t)), and its version in the variable w*t
ToddSeries() = inv((1 - ExpSeries(-1)) << 1)
ToddSeries(w) = ToddSeries() ⊙  GeometricSeries(w)

# Todd class of a rank d bundle, in its Chern classes c_1, ..., c_d
function todd_class(d)
    A = Symbolics.variables(:a, 1:d)

    T = to_elementary(prod([ToddSeries(a) for a in A]), A; sym=:c)
end

# Todd class of the Grassmannian Gr(k, n), in the Chern classes s_i of the
# tautological bundle and q_i of the quotient bundle
function todd_class_grassmannian(k, n)
    A = Symbolics.variables(:a, 1:k) # Chern roots of the tautological line bundle
    B = Symbolics.variables(:b, 1:(n - k)) # Chern roots of the quotient bundle

    T = to_elementary(
            to_elementary(prod([ToddSeries(-a + b) for a in A, b in B]), A; sym=:s),
            B; sym=:q)
    return T
end
