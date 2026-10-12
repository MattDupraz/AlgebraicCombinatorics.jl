module AlgebraicCombinatorics

# Symbolic power series and symmetric polynomials
include("Series.jl")
include("SymmetricPolynomials.jl")
include("Pochhammer.jl")
include("Grassmannians.jl")

# Caching of expensive results to disk
include("Cache.jl")

# Exact linear algebra, lattices and lattice point counting
include("LinAlgTools.jl")

# Matroids and their normal complexes
include("Matroids.jl")
include("NormalComplexes.jl")

# Series.jl
export Series, first_coeffs, first_nonzero_index, cached,
    GeometricSeries, BinomialSeries, ExpSeries, LogSeries, ⊙

# SymmetricPolynomials.jl
export elem_sym_polys, to_elementary, substitute_elementary

# Pochhammer.jl
export falling_factorial, rising_factorial, to_falling_factorials, to_rising_factorials,
    vandermonde

# Grassmannians.jl
export ToddSeries, todd_class, todd_class_grassmannian

# Cache.jl
export @cache

# LinAlgTools.jl
export gram_schmidt, lll_reduction, shortest_vector_approximation,
    reduced_row_echelon_form, extract_basis, row_space_basis, orthogonal_complement,
    hermite_normal_form, lattice_basis_inside_span,
    SignedCone, simplicial_cone_dual, barvinok_decomposition, brion_coefficient,
    lattice_count, polygon_lattice_count

# Matroids.jl
export Matroid, svg_visual, ground_set, flat_index, meet, closure,
    flats_of_rank_range, nonempty_flats, nonempty_proper_flats, chains,
    snapper_polynomial, bergman_ray, bergman_cone_rays, f_to_h, h_to_f

# NormalComplexes.jl
export hypercube_corner, hypercube_corners, simplex_signed_volume, hypercube_volume,
    cone_volume_polynomial, volume_polynomial,
    to_matrix, primitive_vector, hypercube_tangent_cones, pseudo_cubical_cone,
    cone_ehrhart_polynomial,
    pseudo_cubical_heights, face_ehrhart_polynomials, ehrhart_contribution, ehrhart_polynomial

end
