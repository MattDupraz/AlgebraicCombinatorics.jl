using Polymake: Polymake, polytope, matroid, graph, LibPolymake.BigObjectAllocated as PMBigObject, @pm, common
using XML
using Symbolics
using LinearAlgebra: I

# Matroids via their lattice of flats (computed by polymake), together with the
# rays of the Bergman fan and the Snapper polynomial. The volume polynomial is
# in normal_complexes.jl.
#
# Conventions: a flat is referred to by its index into M.flats, and the ground
# set is 0:n-1 (polymake's numbering). Polynomials use Symbolics variables
# f_F (one per non-empty flat F) and h_F (one per non-empty proper flat F).


####################
# MATROID TYPE
####################

# Wrapper so that a polymake object is displayed as an SVG (e.g. in Pluto)
struct SVGVisual
    obj::PMBigObject
end

function Base.show(io::IO, ::MIME"image/svg+xml", x::PMBigObject)
    xml_str = Polymake._get_visual_string_svg(Polymake.visual(Polymake.Visual, x))
    doc = parse(Node, xml_str)
    doc[3]["style"] = "background: #ffffff"
    doc[3]["width"] = "100%"
    doc[3]["height"] = ""
    print(io, XML.write(doc))
end

function Base.show(io::IO, ::MIME"image/svg+xml", svg::SVGVisual)
    show(io, MIME("image/svg+xml"), svg.obj)
end

# A matroid with its flats and their ranks, read from a polymake matroid
struct Matroid
    n_elements::Int
    rank::Int
    flats::Vector{Set{Int}}
    ranks::Vector{Int}       # ranks[i] is the rank of flats[i]
    pm_matroid::PMBigObject
    function Matroid(M::PMBigObject)
        arr = Polymake.Array(M.LATTICE_OF_FLATS.DECORATION)
        flats = [Polymake.decoration_face(dec) for dec in arr]
        ranks = [Polymake.decoration_rank(dec) for dec in arr]
        new(M.N_ELEMENTS, M.RANK, flats, ranks, M)
    end
end

# Picture of the lattice of flats
function svg_visual(M::Matroid)
    SVGVisual(M.pm_matroid.LATTICE_OF_FLATS)
end

function ground_set(M::Matroid)
    return 0:M.n_elements - 1
end


####################
# LATTICE OF FLATS
####################

# Index of the flat equal to the set S (nothing if S is not a flat)
function flat_index(M, S)
    findfirst(map(F -> issetequal(S, F), M.flats))
end

# Meet (intersection) of the flats with indices Fs
function meet(M, Fs)
    S = intersect([M.flats[i] for i in Fs]...)
    flat_index(M, S)
end

# Closure of the set S: the smallest flat containing it
function closure(M, S)
    flat_index(M, intersect(filter(F -> S ⊆ F, M.flats)...))
end

# Join of the flats with indices Fs: the closure of their union.
# A method of Base.join, so that join on strings keeps working.
function Base.join(M::Matroid, Fs)
    S = union([M.flats[i] for i in Fs]...)
    closure(M, S)
end

# Indices of the flats with rank in [r, s]
function flats_of_rank_range(M, r, s)
    return [i for i in 1:length(M.flats)
            if r <= M.ranks[i] && M.ranks[i] <= s]
end

function nonempty_flats(M)
    return flats_of_rank_range(M, 1, M.rank)
end

function nonempty_proper_flats(M)
    return flats_of_rank_range(M, 1, M.rank - 1)
end

# Returns all chains of flats from Fs that extend `tail` downwards, including
# `tail` itself. Chains are listed from smallest to largest flat.
function chains(M, Fs; tail=[])
    cs = [tail]
    for F in Fs
        if isempty(tail) || M.flats[F] ⊊ M.flats[first(tail)]
            append!(cs, chains(M, Fs; tail=[F; tail]))
        end
    end
    return cs
end

# Returns the non-empty chains of proper nonempty flats of M
# (the cones of the Bergman fan)
function chains(M)
    return filter(!isempty, chains(M, nonempty_proper_flats(M)))
end


####################
# SNAPPER POLYNOMIAL
####################

# Iterates over the integer vectors m with 0 <= m[i] <= max_m[i] and sum(m) <= max_sum
# (as plain integers rather than vectors if length(max_m) == 1)
function bounded_vectors(max_m, max_sum)
    if length(max_m) == 1
        return 0:min(max_m[1], max_sum)
    end
    return Iterators.flatmap(
        n -> Iterators.product(n, bounded_vectors(max_m[2:length(max_m)], max_sum - n))
                .|> t -> [t[1]; t[2]],
        0:min(max_m[1], max_sum))
end

# Iterates over the integer vectors (as tuples) between 0 and vec, componentwise
function vectors_below(vec)
    return Iterators.product([0:n for n in vec]...)
end

# Snapper polynomial of M in the variables f_F, F a non-empty flat:
#   Σ_m Π_F binomial(f_F + m_F - 1, m_F)
# over the vectors m with m_F <= rk(F) - 1 and |m| <= rk(M) satisfying the
# Dragon-Hall-Rado condition: rk(∨{F : m'_F > 0}) > |m'| for all 0 != m' <= m.
# Writes the current progress to the file "progress".
function snapper_polynomial(M)
    ne_flats = nonempty_flats(M)

    f = Symbolics.variables(:f, ne_flats)

    ms = bounded_vectors(max.(0, M.ranks[ne_flats] .- 1), M.rank) |> collect
    snap = 0

    for (i, m) in enumerate(ms)
        #write("progress", "$i/$(length(ms)) $m")

        # Check the Dragon-Hall-Rado condition for m
        dragon_hall_rado = true
        for m_ in collect(vectors_below(m))
            if all(iszero, m_)
                continue
            end
            S = join(M, [F for (j, F) in enumerate(ne_flats) if m_[j] > 0])

            if M.ranks[S] <= sum(m_)
                dragon_hall_rado = false
                break
            end
        end

        if dragon_hall_rado
            # binomial(f_F + m_F - 1, m_F) = Π_{i=0}^{m_F-1} (f_F + i) / (i + 1),
            # expanded by hand since Base.binomial stays unevaluated on symbolic input
            snap += prod(prod((f[j] + i) * (1//(i+1)) for i in 0:m[j]-1; init=1)
                         for j in eachindex(m))
        end
    end

    return Symbolics.expand(snap)
end


####################
# BERGMAN FAN
####################

# Ray of the Bergman fan for the flat F: the indicator vector of F modulo the
# all-ones vector, represented by subtracting the last coordinate and dropping it
function bergman_ray(M, F)
    v = Rational[e in M.flats[F] for e in ground_set(M)] * BigInt(1)
    v_proj = v[1:length(v) - 1] .- last(v)
    return v_proj
end

# Rays (as rows) of the cone of the Bergman fan given by a chain of flats,
# optionally in the coordinates given by the matrix B
function bergman_cone_rays(M, cone; B=I)
    return reduce(hcat, [B * bergman_ray(M, F) for F in cone])'
end


####################
# CHANGE OF VARIABLES
####################

# Substitution expressing each f_F (F non-empty) in terms of the h_G (G non-empty proper)
function f_to_h(M)
    nep_flats = nonempty_proper_flats(M)
    n = length(nep_flats)
    f = Symbolics.variables(:f, nonempty_flats(M))
    h = Symbolics.variables(:h, nep_flats)

    # A is a change of basis from (L_F - H) to D_F
    A = Rational.(reduce(hcat, [-1 * [M.flats[F] ⊆ M.flats[G] for G in nep_flats]
                                for F in nep_flats]))
    return Dict(f .=> Rational.([Matrix(I, n, n); -ones(n)']) * inv(A) * h)
end

# Substitution expressing each h_G (G non-empty proper) in terms of the f_F (F non-empty)
function h_to_f(M)
    nep_flats = nonempty_proper_flats(M)
    f = Symbolics.variables(:f, nonempty_flats(M))
    h = Symbolics.variables(:h, nep_flats)

    # Hyperplane class
    H = [ground_set(M)[1] ∈ M.flats[F] for F in nep_flats]
    # A is a change of basis from L_F to D_F
    A = Rational.(reduce(hcat, [H - [M.flats[F] ⊆ M.flats[G] for G in nep_flats]
                                for F in nonempty_flats(M)]))
    return Dict(h .=> A * f)
end
