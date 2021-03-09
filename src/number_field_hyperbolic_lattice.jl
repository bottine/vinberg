using LinearAlgebra
using JSON
using StaticArrays
using ToggleableAsserts
using Hecke
using AbstractAlgebra



function Base.isless(x::nf_elem,y::nf_elem)
    @assert parent(x) == parent(y)
    field = parent(x)
    @assert istotally_real(field)
    places = infinite_places(field)
    
    if x-y == 0
       return false
    end
    return isnegative(x-y,places[1])
end

Base.isless(x::nf_elem,y) =  Base.isless(x,parent(x)(y))
Base.isless(x,y::nf_elem) = Base.isless(parent(y)(x),y)

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem},y) =  Base.isless(x,parent(x).nf(y))
Base.isless(x,y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(parent(y).nf(x),y)

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem}, y::nf_elem) = Base.isless(parent(x).nf(x),y) 
Base.isless(x::nf_elem, y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(x,parent(y).nf(y))

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem}, y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(parent(x).nf(x),parent(y).nf(y)) 

Maybe{T} = Union{T,Nothing}

# Everything goes there, eventually
mutable struct VinbergData
    field#::AnticNumberField
    ring#::Maybe{NfAbsOrd{AnticNumberField,nf_elem}}
    gram_matrix#::AbstractAlgebra.Generic.MatSpaceElem{NfAbsOrdElem{AnticNumberField,nf_elem}}
    quad_space#::Maybe{Hecke.QuadSpace{AnticNumberField,AbstractAlgebra.Generic.MatSpaceElem{nf_elem}}}
    lattice#::Maybe{QuadLat{AnticNumberField,AbstractAlgebra.Generic.MatSpaceElem{nf_elem},Hecke.PMat{nf_elem,Hecke.NfAbsOrdFracIdl{AnticNumberField,nf_elem}}}}

    least_k_by_root_length::Maybe{Dict}
    
end

"""
    is_feasible(field,matrix)

    Check that the quadratic form defined by the matrix `matrix` has signature ``(n,1)`` and that all its (non-trivial) Galois conjugates have signature ``(n+1,0)`` (Equivalently, all are definite, the original is hyperbolic and the (non-trivial) conjugates are positive definite).
"""
function is_feasible(V::Hecke.QuadSpace)
    
    P = infinite_places(V.K)
    n = dim(V)

    return signature(V,P[1]) == (dim(V)-1,1,0) && all(signature(V,p) == (dim(V),0,0) for p in P[2:end])
end


function signature(V::Hecke.QuadSpace,P::InfPlc)
       diag = [inner_product(V,v,v) for v in collect(Vector{nf_elem},eachrow(Matrix(orthogonal_basis(V))))]
       filter!(≠(0),diag)
       sig =  (count([ispositive(d,P) for d in diag]),count([isnegative(d,P) for d in diag]),dim(V)-length(diag))

       @assert sum(sig) == dim(V)
       return sig
end

function diag_signature(field,diagonal,place)
    non_zeros = filter(≠(0),diagonal)
    (count([ispositive(field(d),place) for d in non_zeros]),count([isnegative(field(d),place) for d in non_zeros]),length(diagonal)-length(non_zeros))
end

function is_diago_and_feasible(field,matrix)
    
    (m,n) = size(matrix)
    
    m ≠ n && return false
    !isdiagonal(matrix) && return false

    diagonal = [matrix[i,i] for i in 1:n]
   
    P = infinite_places(field)

    diag_signature(field,diagonal,P[1]) ≠ (n-1,1,0) && return false
    any(diag_signature(field,diagonal,p) ≠ (n,0,0) for p in P[2:end]) && return false

    # the diagonal is ordered increasingly (this is not really necessary)
    for i in 1:n-1
        diagonal[i] > diagonal[i+1] && return false 
    end

    return true

end

function VinbergData(field,matrix)

    println("hello")
    ring = maximal_order(field)
    quad_space = quadratic_space(field, matrix)
    quad_lattice = Hecke.lattice(quad_space)
    println("hello")

    @assert all(field(c) ∈ ring for c in matrix) "The Gram matrix must have coefficients in the ring of integers."

    println("OK")
    @assert is_diago_and_feasible(field,matrix)

    println("world")

    least_k_by_root_length = Dict([ring(l) => (ring(0),NfAbsOrdElem{AnticNumberField,nf_elem}[]) for l in possible_root_norms_up_to_squared_units(ring,field,quad_space)])
    
    vd = VinbergData(field,ring,matrix,quad_space,quad_lattice,least_k_by_root_length)

    return vd

end

# useless but kept for reference
function Base.getproperty(vd::VinbergData,::Val{:ring})
    if isnothing(vd.ring)
        vd.ring = maximal_order(vd.field)
    end
    return vd.ring
end

"""

Goal: enumerate all ``k=k₀`` such that ``r = (k₀,…)`` can form a root of length ``l`` (i.e. such that there exist coeffs ``k₁,…,k_n`` forming a root).
`α` is the negative eigenvalue of the diag form, with corresponding vector `v₀`.
`ls` contain the  inner products ``r⊙r`` that we need to consider. 
We need to have, for all non identity Galois automorphism ``σ`` that ``|σ(k²/lα)|  ≤ 1``, which means
`` σ(k)² ≤ |σ(lα)|.

We can therefore take the max ``M = \\max_{σ≠id}|σ(lα)|``
and we will know that if ``0≤k≤n``, then ``t₂(k) ≤ k² + GM``, where ``G`` is one minus the number of conjugates.

"""
function enumerate_k(vd::VinbergData,l,k_min,k_max)
    
    println("hello enumerate k for $l")

    ring = vd.ring
    field = vd.field

    α = vd.gram_matrix[1,1]
    @assert α < 0 

    M = maximum(abs.(convert.(Float64,conjugates_real(field(l*α),32)[2:end])))
    P = infinite_places(field)
    num_non_trivial_conjugates = length(P)-1
    k_min_squared_approx = Float64(conjugates_real(field(k_min^2),32)[1])
    k_max_squared_approx = Float64(conjugates_real(field(k_max^2),32)[1])


    trace_mat = trace_matrix(ring)
    ring_basis = basis(ring)
    candidates = [ring(dot(ring_basis,v)) for (v,magnitude) in Hecke.short_vectors(Zlattice(gram = trace_mat), k_min_squared_approx , k_max_squared_approx  + num_non_trivial_conjugates*M + 1)]
    # the -1 and +1 are just here to have a security margin
    println("Got ", length(candidates), " candidates.")
    println([conjugates_real(field(c),8)[1] for c in candidates])
    println([c for c in candidates])
    # Only k≥0
    filter!(
        k-> field(k) ≥ 0,
        candidates,    
    )
    
    println("Got ", length(candidates), " candidates.")

    #  |σ(k²/lα)| ≤ 1 ⇔ σ(k²) ≤ |σ(lα)| ⇔ σ(k⁴) ≤ σ(l²α²) ⇔ σ(l²α² - k⁴) ≥ 0  
    filter!(
        k -> all(
            field((l*α)^2 - k^4) ≥ 0 
            for p in P[2:end]
        ),
        candidates,
    )
    println("Got ", length(candidates), " candidates.")
    return candidates
end

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


function possible_root_norms_up_to_squared_units(
    ring,
    field,
    space,
)

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

    # this is the ideal ⟨2*last_invariant_factor⟩ = ⟨2*det(gram)/gcd(cofactors(gram))⟩
    # constructed by taking the ideal I :=⟨gcd(cofactors(gram))⟩
    # and the ideal                   J := ⟨det(gram)⟩
    # then taking the ideal (I:J) = \{x: xI ⊆ J\}
    # then multiplying by the ideal ⟨2⟩
    # Probably one can do something cleaner
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
    return [ring(l) for l in unique(all_root_norms)]

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

    #Not gonna work since it is only up to squared units: neeed to use the morphisms constructed in possible_root_norms_up_to_squared_units
    #@assert inner_product(space,vector,vector) ∈ possible_root_norms_up_to_squared_units(space,ring)

    return true

end

function is_root(vd::VinbergData,vector)
    return is_root(vd.quad_space,vd.ring,vector)
end

function next_min_k_for_l!(vd,l)

    (current_k,remaining) = vd.least_k_by_root_length[l] 

    k_min,k_max = current_k,current_k+1

    while isempty(remaining)
        println("remaining is empty!")
        remaining = enumerate_k(vd,l,k_min,k_max)
        filter!(>(current_k),remaining)
        sort!(remaining)
        k_min,k_max = k_max,k_max+1
    end
    new_k = popfirst!(remaining)
    vd.least_k_by_root_length[l] = (new_k,remaining)
end

function next_min_ratio!(vd::VinbergData)
    field = vd.field
    min_pair = nothing
    for (l,(k,remaining)) in vd.least_k_by_root_length
        if isnothing(min_pair) || field(k)//field(l) < field(min_pair[1])//field(min_pair[2])
            min_pair = (k,l)
        end
    end
    next_min_k_for_l!(vd,min_pair[2])
    return min_pair
end

function extend_root_stem(vd::VinbergData,stem)
    
    
    
end


Qx, x = Hecke.QQ["x"]
f = x^2 - 5
K, a = Hecke.NumberField(f, "a"); @assert istotally_real(K)
M = matrix(K, 2, 2, [a, 0, 0 , 1]) # The bilinear form



#
#
#########################
##
#QK,Qa = Hecke.rationals_as_number_field()
#OQK = ring_of_integers(QK)
#QM = matrix(QK, 4, 4, [-1,0,0,0,  0,2,0,0,  0,0,6,0,  0,0,0,6])
#QV = quadratic_space(QK,QM)
#
