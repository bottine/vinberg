using AbstractAlgebra
using LinearAlgebra
using JSON
using StaticArrays
using ToggleableAsserts
using Hecke




using Markdown
using InteractiveUtils




Qx, x = QQ["x"]
f = x^2 - 5
K, a = NumberField(f, "a"); @assert istotally_real(K)
OK = ring_of_integers(K)
M = matrix(K, 2, 2, [-a, 0, 0 , 1]) # The bilinear form
V = quadratic_space(K, M)
L = lattice(V)


function get_possible_root_norms(
	number_field,# The number field in which we play
	algebraic_integers, # Its ring of algebraic integers
	quad_space, # A quadratic space 
	lattice, # A corresponding lattice
)
	# Probably all arguments can be recovered from the lattice only -> TODO
	K = number_field
	OK = algebraic_integers
	
	U, mU = unit_group(OK)
	Q, mQ = quo(U, 2) # quo(U,2) yields Q = U/U² and the quotient morphism mQ: U -> Q
	representatives_of_units_up_to_squares = [ mU(preimage(mQ, q)) for q in Q]
	
	twice_discriminant = ideal(OK,2*discriminant(V))
	twice_discriminant_factors = factor(twice_discriminant)
	@assert all(Hecke.isprincipal_fac_elem(idl)[1] for idl in keys(twice_discriminant_factors))
	
	unit_factors_of_root_lengths = Dict([u=>1 for u in representatives_of_units_up_to_squares])
	prime_factors_of_root_lengths = Dict([evaluate(Hecke.isprincipal_fac_elem(idl)[2]) => mul for (idl,mul) in twice_discriminant_factors])
	all_factors_of_root_lengths = merge(prime_factors_of_root_lengths, unit_factors_of_root_lengths)
	
	all_root_norms = 
	filter(
		l -> istotally_positive(K(l)),
		products(all_factors_of_root_lengths)
	)
	return unique(all_root_norms)
	
end

