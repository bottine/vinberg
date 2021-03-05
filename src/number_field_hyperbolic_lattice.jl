using LinearAlgebra
using JSON
using StaticArrays
using ToggleableAsserts
using Hecke



function ideal_gcd(ring,elems)
    idls = [ideal(ring,ring(e)) for e in elems]
    gcd_idls = reduce(gcd,idls)
    return gcd_idls
end

function elem_gcd(ring,elems)
    is_principal,gcd_elems = Hecke.isprincipal(ideal_gcd(ring,elems))
    @assert is_principal 
    return gcd_elems
    
end

function products(d)
    ks = collect(keys(d))
    vs = collect(values(d))
    partial_products = Base.product([collect(0:v) for v in vs]...)
    return [prod([k^p for (k,p) in zip(ks,pp)]) for pp in partial_products]
end

function get_possible_root_norms(
    quad_space::Hecke.QuadSpace
)

    number_field = quad_space.K
    algebraic_integers = ring_of_integers(number_field)
    
    return get_possible_root_norms(number_field,algebraic_integers,quad_space)
end

function get_possible_root_norms(
	quad_space, # A quadratic space 
	algebraic_integers, # Its ring of algebraic integers
)
    ring = algebraic_integers
    space = quad_space
    field = space.K

    units, morphism_units = unit_group(ring)
    twice_units, morphism_twice_units = quo(units, 2) # quo(U,2) yields Q = U/U² and the quotient morphism mQ: U -> Q
    representatives_of_units_up_to_squares = [ morphism_units(preimage(morphism_twice_units, q)) for q in twice_units]

    # is it true in general that the lengths divide twice the last invariant factor?
    # PROOF??? TODO

    # last invariant factor is det(G) / gcd(all elements of the cofactor matrix)
    # We use ideals to compute the gcd
    gram = space.gram 
    cofactors = det(gram) * inv(gram)
    gcd_cofactors = ideal_gcd(ring,collect(cofactors))


    twice_last_invariant_factor = ideal(ring,2)*colon(ideal(ring,ring(det(gram))),gcd_cofactors)
    twice_last_invariant_factor_factors = Hecke.factor(twice_last_invariant_factor)
    @assert all(Hecke.isprincipal(idl)[1] for idl in keys(twice_last_invariant_factor_factors))

    unit_factors_of_root_lengths = Dict([u=>1 for u in representatives_of_units_up_to_squares])
    prime_factors_of_root_lengths = Dict([Hecke.isprincipal(idl)[2] => mul for (idl,mul) in twice_last_invariant_factor_factors])
    all_factors_of_root_lengths = merge(prime_factors_of_root_lengths, unit_factors_of_root_lengths)

    all_root_norms = 
    filter(
        l -> istotally_positive(field(l)),
        products(all_factors_of_root_lengths)
    )
    return unique(all_root_norms)

end

function is_integral(space,ring,vector)
    field = space.K
    
    for c in vector
        if !in(field(c),ring)
            return false
        end
    end
    return true
end

function has_positive_norm(space,ring,vector)
    field = space.K
    
    l = inner_product(space,vector,vector)

    # positive norm
    P = infinite_places(field)
    if !ispositive(l, P[1]) # Taking P[1] means considering our canonical embedding in ℝ
        return false
    end
    return true
end

function is_primitive(space,ring,vector)
    return isunit(ideal_gcd(ring,vector))
end

function crystallographic_condition(space,ring,vector)
    l = inner_product(space,vector,vector)

    for b in eachcol(I(length(vector)))
        if !in(2*inner_product(space,collect(b),vector),ideal(ring,l))
            return false
        end
    end
    return true
end

function is_root(space,ring,vector)

    !has_positive_norm(space,ring,vector) && return false

    !is_integral(space,ring,vector) && return false

    !is_primitive(space,ring,vector) && return false

    !crystallographic_condition(space,ring,vector) && return false 

    @assert inner_product(space,vector,vector) ∈ get_possible_root_norms(space,ring)

    return true

end

#=
function alg_integers_with_bounded_T2_norm(OK,max_norm)
    G = trace_matrix(OK)
    candidates = Hecke.__enumerate_gram(G, max_norm)
    res = [dot(g[1],basis(OK)) for g in candidates]
    @assert all(t2(r) ≤ max_norm for r in res)
    return res
end

function next_k₀(OK,V,l,remaining,max)
    
    m = 

    while isempty(remaining)
        max = max+1
       
        # max is the max value of k that we look for
        # Then |σ(k²/l*⟨v₀,v₀⟩)| ≤ 0 means |σ(k)|² ≤ |σ(l*⟨v₀,v₀⟩)|² ≤ ⟨v₀,v₀⟩*t₂(l)

        max_T2_norm = max^2 + ?? 

        candidates = alg_integers_with_bounded_T2_norm(OK)

    end

    k = pop!(remaining)
    return k
end

function enumerate_ratios(G)
    
    #Need to enumerate k₀>0 increasingly

    # Let l denote the possible norm²: need to have σ(k₀²/l⟨v₀,v₀⟩) ≤ 1 for all non-identity Galois conjugates σ (see talk by Mark, timestamp 15:21)
    # This means we can bound the T₂ norm of k₀

end
=#


##########################

Qx, x = QQ["x"]
f = x^2 - 5
K, a = NumberField(f, "a"); @assert istotally_real(K)
OK = ring_of_integers(K)
M = matrix(K, 2, 2, [-a, 0, 0 , 1]) # The bilinear form
V = quadratic_space(K,M)
L = lattice(V)


########################
#
QK,Qa = Hecke.rationals_as_number_field()
OQK = ring_of_integers(QK)
QM = matrix(QK, 4, 4, [-1,0,0,0,  0,2,0,0,  0,0,6,0,  0,0,0,6])
QV = quadratic_space(QK,QM)

