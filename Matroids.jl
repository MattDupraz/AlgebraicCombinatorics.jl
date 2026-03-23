module Matroids

using Polymake: Polymake, polytope, matroid, graph, LibPolymake.BigObjectAllocated as PMBigObject, @pm, common
using XML

export Matroid
export svg_visual, groundset, get_flat_index, meet, closure, join, flats_of_rank_range, nonempty_flats, nonempty_proper_flats, chains

# ╔═╡ 6ae0a0b0-0134-4cc3-a52d-fda551ace3c3
struct SVGVisual
	obj::PMBigObject
end

# ╔═╡ 56fea39d-0a9d-4f70-8892-934a88264123
function Base.show(io::IO, ::MIME"image/svg+xml", x::PMBigObject)
	xml_str = Polymake._get_visual_string_svg(Polymake.visual(Polymake.Visual, x))
	doc = parse(Node, xml_str)
	doc[3]["style"] = "background: #ffffff"
	doc[3]["width"] = "100%"
	doc[3]["height"] = ""
	print(io, XML.write(doc))
end

# ╔═╡ db24b03a-02dc-4d9b-834d-b4387b46c875
function Base.show(io::IO, ::MIME"image/svg+xml", svg::SVGVisual)
	show(io, MIME("image/svg+xml"), svg.obj)
end

# ╔═╡ cda53267-e690-4891-9a53-581e1ff683dd
struct Matroid
	n_elements::Int
	rank::Int
	flats::Vector{Set{Int}}
	ranks::Vector{Int}
	pm_matroid::PMBigObject
	function Matroid(M::PMBigObject)
		arr = Polymake.Array(M.LATTICE_OF_FLATS.DECORATION)
		flats = [Polymake.decoration_face(dec) for dec in arr]
		ranks = [Polymake.decoration_rank(dec) for dec in arr]
		new(M.N_ELEMENTS, M.RANK, flats, ranks, M)
	end
end

# ╔═╡ 002a96fe-1395-4200-9074-dcc591b39fc1
function svg_visual(M::Matroid)
	SVGVisual(M.pm_matroid.LATTICE_OF_FLATS)
end

# ╔═╡ 5174ba72-74e8-4ae2-b504-fa4a36f84b01
function groundset(M::Matroid)
	return 0:M.n_elements - 1
end

# ╔═╡ dfdeb32d-6c6e-4704-8006-ad7ab1233ed7
function get_flat_index(M, S)
	findfirst(map(F -> issetequal(S, F), M.flats))
end

# ╔═╡ eb4541cd-29fc-4f41-b668-97701ba9efbe
function meet(M, Fs)
	S = intersect([M.flats[i] for i in Fs]...)
	get_flat_index(M, S)
end

# ╔═╡ 2f0ed098-5d7e-4f3e-8ef0-97cb09d9dad4
function closure(M, S)
	get_flat_index(M, intersect(filter(F -> S ⊆ F, M.flats)...))
end

# ╔═╡ 4fb7b3e5-86d5-4b6a-94f8-b155ab1b6e2d
function join(M, Fs)
	S = union([M.flats[i] for i in Fs]...)
	closure(M, S)
end

# ╔═╡ 8176fa85-c220-4c41-a221-9f48eb090e74
function flats_of_rank_range(M, r, s)
	return [i for i in 1:length(M.flats)
		if r <= M.ranks[i] && M.ranks[i] <= s ]
end

# ╔═╡ ff29b293-0bb7-491c-85da-1112d61a5f32
function nonempty_flats(M)
    return flats_of_rank_range(M, 1, M.rank)
end

# ╔═╡ ceab87f9-4a40-4de4-9b16-49141064393a
function nonempty_proper_flats(M)
    return flats_of_rank_range(M, 1, M.rank - 1)
end

# ╔═╡ 4b85ad51-8cc6-4382-966e-dc6ce3172690
# Return the chains of flats from Fs
function chains(M, Fs; tail=[])
	cs = [tail]
    for F in Fs
        if isempty(tail) || M.flats[F] ⊊ M.flats[first(tail)]
            append!(cs, chains(M, Fs; tail=[F; tail]))
		end
	end
	return cs
end

# ╔═╡ 00681b2c-6ae0-48b5-8cfd-29d90384a7dd
# Returns the non-empty chains of proper nonempty flats of M
function chains(M)
	return filter(!isempty, chains(M, nonempty_proper_flats(M)))
end

end
