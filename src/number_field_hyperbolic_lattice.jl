using LinearAlgebra
using JSON
using StaticArrays
using ToggleableAsserts
using Hecke
using AbstractAlgebra
using Convex, Cbc, COSMO
import MathOptInterface

#include("util.jl")
include("diagrams.jl")

using Main.Diagrams: build_diagram_and_subs, extend!, is_finite_volume







###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#
#    Misc
#
###

function products(d)
    ks = collect(keys(d))
    vs = collect(values(d))
    partial_products = Base.product([collect(0:v) for v in vs]...)
    return [prod([k^p for (k,p) in zip(ks,pp)]) for pp in partial_products]
end

"""
    Coxeter_coeff(r₁,r₂)

Compute the label corresponding to the edge between the vertices corresponding to the hyperplanes defined by `r₁` and `r₂` respectively.

# Remarks

* This is incomplete, and only works for acute angles I think.
* The code is copied from B&P, and I haven't yet made in more general.
"""
function Coxeter_coeff(space, ring, r₁, r₂)
    @toggled_assert is_root(space, ring, r₁) && is_root(space, ring, r₂) "The elements must be roots"

    angle = Hecke.inner_product(space,r₁,r₂)
    cos² = approx(angle^2//(Hecke.inner_product(space,r₁,r₁)*Hecke.inner_product(space,r₂,r₂)))
    if cos² == 0
        return 2
    elseif cos² == 1
        return 0
    elseif cos² > 1
        return 1
    else

        #   cos(π/m)² = r₁⋅r₂ / r₁² r₂² 
        # ⇒ cos(π/m) = √(r₁⋅r₂/r₁r₂)
        # ⇒ m = π/acos(√(r₁⋅r₂/r₁r₂))

        # TODO Is this rounding dangerous?
        #      Can the angle be different from a submultiple of π? if yes, how to deal with it?

        m = round(Int, π/acos(√cos²))

        return m
    end

end

get_Coxeter_matrix(space, ring, roots) = reduce(hcat,[[Coxeter_coeff(space, ring, r₁,r₂) for r₁ in roots] for r₂ in roots])

# TODO: make exact
function is_necessary_halfspace(cone_roots,root) 
   


    float_cone_roots = Vector{Vector{Float64}}([approx.(cone_root) for cone_root in cone_roots])    
    float_root = Vector{Float64}(approx.(root))
    
    n = length(root) 

    # x' * (A * r) ≤ 0 ∀ r
    # (A * r)' * x ≤ 0 ∀ r

    #x = Variable(n, IntVar)
    x = Variable(n)
    p = satisfy()       # satisfiability question 
    for cone_root in float_cone_roots
        p.constraints += x' * cone_root ≤ 0 # hyperplanes defining the cone
    end
    p.constraints += x' * float_root ≥ 1 # other side of the half space defined by root
    # it should only be strictly bigger than zero, but Convex.jl does not do "strictly", so we change it to ≥ 1 (and since we have a cone, it should be the same result)

    
    Convex.solve!(p,Cbc.Optimizer(verbose=0,loglevel=0), verbose=false, warmstart=false)
    #solve!(p,COSMO.Optimizer(verbose=false), verbose=false)
   

    if p.status == MathOptInterface.INFEASIBLE 
        return false
    elseif p.status == MathOptInterface.OPTIMAL
        #println(p.optval)
        return true
    else
        println("can't tell! ($(p.status))")
        println("             $(p))")
    end

end


function drop_redundant_halfspaces(
    roots
) 
    
    
    for i in length(roots):-1:1
       
        rr = copy(roots)
        r = popat!(rr,i)

        if ! is_necessary_halfspace(rr,r)
            return drop_redundant_halfspaces(rr) 
        end
    end
    
    return roots
    
end










###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#
#    Algebra stuff
#
###


function Base.:>(x::nf_elem,y::nf_elem,p::InfPlc)
    x == y && return false
    return ispositive(x-y,p)
end
function Base.:≥(x::nf_elem,y::nf_elem,p::InfPlc)
    x == y && return true
    return ispositive(x-y,p)
end

Base.:>(x::nf_elem,y,p::InfPlc) = >(x,parent(x)(y),p)
Base.:>(x,y::nf_elem,p::InfPlc) = >(parent(y)(x),y,p)
Base.:≥(x::nf_elem,y,p::InfPlc) = ≥(x,parent(x)(y),p)
Base.:≥(x,y::nf_elem,p::InfPlc) = ≥(parent(y)(x),y,p)

Base.:<(x,y,p::InfPlc) = >(y,x,p)
Base.:≤(x,y,p::InfPlc) = ≥(y,x,p)

function Base.isless(x::nf_elem,y::nf_elem)
    @assert parent(x) == parent(y)
    field = parent(x)
    @assert istotally_real(field)
    places = infinite_places(field)
    
    return <(x,y,places[1])
end

Base.isless(x::nf_elem,y) =  Base.isless(x,parent(x)(y))
Base.isless(x,y::nf_elem) = Base.isless(parent(y)(x),y)

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem},y) =  Base.isless(x,parent(x).nf(y))
Base.isless(x,y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(parent(y).nf(x),y)

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem}, y::nf_elem) = Base.isless(parent(x).nf(x),y) 
Base.isless(x::nf_elem, y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(x,parent(y).nf(y))

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem}, y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(parent(x).nf(x),parent(y).nf(y)) 


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

divides(a,b,ring) = b ∈ ideal(ring,a)

function t2_exact(x::S) where S <: NumFieldElem
    @assert istotally_real(parent(x))
    return trace(x^2)
end
function t2_exact(x::NfAbsOrdElem)
  return t2_exact(x.elem_in_nf)
end

function short_t2_elems(O::NfAbsOrd, lb, ub)
    @assert istotally_real(nf(O))

    trace = Hecke.trace_matrix(O)
    basis = Hecke.basis(O)

    lat = Hecke.short_vectors(Zlattice(gram = trace), lb, ub)
    candidates = [O(Hecke.dot(basis,v)) for (v,t) in lat]



    @assert all(lb-1 ≤ t2(c) && t2(c) ≤ ub+1 for c in candidates)
    return candidates
end

function non_neg_short_t2_elems(O::NfAbsOrd, lb, ub)
    candidates = short_t2_elems(O,lb,ub)
    if lb == 0
        push!(candidates,O(0))
    end
    map!(c -> (c<0 ? -c : c), candidates, candidates)
    return candidates
end

approx(x) = Float64(conjugates_real(x)[1])
approx(x::NfAbsOrdElem{AnticNumberField,nf_elem}) = conjugates_real(x.elem_in_nf)[1]

function diagm(K::AnticNumberField,diag)
    n = length(diag)
    M = fill(K(0),n,n)
    for i in 1:n
        M[i,i] = K(diag[i])
    end
    return matrix(K,M)
end

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#
#    Quad Space & Lattice stuff
#
###


"""
    is_feasible(field,matrix)

    Check that the quadratic form defined by the matrix `matrix` has signature ``(n,1)`` and that all its (non-trivial) Galois conjugates have signature ``(n+1,0)`` (Equivalently, all are definite, the original is hyperbolic and the (non-trivial) conjugates are positive definite).
"""
function is_feasible(V::Hecke.QuadSpace)
    
    P = infinite_places(V.K)
    n = dim(V)

    return signature(V,P[1]) == (dim(V)-1,1,0) && all(signature(V,p) == (dim(V),0,0) for p in P[2:end])
end


function signature(
    V::Hecke.QuadSpace,
    P::InfPlc,
)
    diag = diagonal(V) 
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



###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#
#    Roots stuff
#
###


function is_integral(space,ring,vector)
    field = space.K
    
    for c in vector
        if !Hecke.in(field(c),ring)
            return false
        end
    end
    return true
end

function has_positive_norm(space,ring,vector)
    field = space.K
    
    l = Hecke.inner_product(space,vector,vector)

    if l ≤ 0 
        return false
    end
    return true
end

function is_primitive(space,ring,vector)
    return isunit(ideal_gcd(ring,vector))
end

function crystallographic_condition(space,ring,vector)
    l = Hecke.inner_product(space,vector,vector)

    for b in eachcol(I(length(vector)))
        if  !divides(l,2*Hecke.inner_product(space,collect(b),vector),ring)
            return false
        end
    end
    return true
end


function is_root(space,ring,vector)

    @debug "is_root($vector)"

    !has_positive_norm(space,ring,vector) && return false
    
    @debug "✓ positive length"

    !is_integral(space,ring,vector) && return false
    
    @debug "✓ integral"

    !is_primitive(space,ring,vector) && return false
    
    @debug "✓ primitive"

    !crystallographic_condition(space,ring,vector) && return false 
    
    @debug "✓ crystallographic"

    #Not gonna work since it is only up to squared units: neeed to use the morphisms constructed in possible_root_norms_up_to_squared_units
    #@assert Hecke.inner_product(space,vector,vector) ∈ possible_root_norms_up_to_squared_units(space,ring)

    return true

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













###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#
#    Vinberg Algo stuff
#
###


@enum VinbergDataStage Initialized ConeRoots MoreRoots Done


Maybe{T} = Union{T,Nothing}

# Everything goes there, eventually
mutable struct VinbergData

    stage::VinbergDataStage 

    dim::Int
    field::AnticNumberField
    ring::NfAbsOrd{AnticNumberField,nf_elem}
    gram_matrix::AbstractAlgebra.Generic.MatSpaceElem{nf_elem}
    quad_space::Hecke.QuadSpace{AnticNumberField,AbstractAlgebra.Generic.MatSpaceElem{nf_elem}}
    #lattice::QuadLat{AnticNumberField,AbstractAlgebra.Generic.MatSpaceElem{nf_elem},Hecke.PMat{nf_elem,Hecke.NfAbsOrdFracIdl{AnticNumberField,nf_elem}}}


    least_k_by_root_length::Dict{NfAbsOrdElem{AnticNumberField,nf_elem},Tuple{NfAbsOrdElem{AnticNumberField,nf_elem},Array{NfAbsOrdElem{
    AnticNumberField,nf_elem},1}}}
    candidate_roots::Vector{Vector{nf_elem}}
    accepted_roots::Vector{Vector{nf_elem}}
    
end



function VinbergData(field,matrix)

    (n,m) = size(matrix)
    @assert n == m
    dim = n

    ring = maximal_order(field)
    quad_space = quadratic_space(field, matrix)
    #quad_lattice = Hecke.lattice(quad_space)

    println("VinbergData($field,$matrix)")
    println("matrix of type $(typeof(matrix))")

    @assert all(field(c) ∈ ring for c in matrix) "The Gram matrix must have coefficients in the ring of integers."

    @assert is_diago_and_feasible(field,matrix)


    least_k_by_root_length = Dict([ring(l) => (ring(0),NfAbsOrdElem{AnticNumberField,nf_elem}[]) for l in possible_root_norms_up_to_squared_units(ring,field,quad_space)])
    
    vd = VinbergData(Initialized,dim,field,ring,matrix,quad_space,#=quad_lattice,=#least_k_by_root_length,[],[])

    return vd

end

diag_coeff(vd::VinbergData,i::Int) = vd.gram_matrix[i,i] 

basepoint(vd::VinbergData) = vd.field.(vcat([1],zeros(Int,vd.dim-1)))
vector_length(vd,root) = Hecke.inner_product(vd.quad_space,root,root)
fake_dist_to_basepoint(vd,root) = (root[1]^2//vector_length(vd,root))
⊙(v,vd::VinbergData,w) = Hecke.inner_product(vd.quad_space,v,w)

# useless but kept for reference
function Base.getproperty(vd::VinbergData,::Val{:ring})
    if isnothing(vd.ring)
        vd.ring = maximal_order(vd.field)
    end
    return vd.ring
end


# TODO make that better
sum_at_places(val,first_place_idx) = sum(convert.(Float64,conjugates_real(val,32)[first_place_idx:end]))

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
    
    ring = vd.ring
    field = vd.field

    α = vd.gram_matrix[1,1]
    @assert α < 0 

    M = sum_at_places(field(l//α),2)
    P = infinite_places(field)
    k_min_squared_approx = Float64(conjugates_real(field(k_min^2),32)[1])
    k_max_squared_approx = Float64(conjugates_real(field(k_max^2),32)[1])


    candidates = non_neg_short_t2_elems(ring, k_min_squared_approx-1 , k_max_squared_approx  + M + 1)
    # the -1 and +1 are just here to have a security margin
    
    candidates = vcat(candidates, .- candidates)

    #Only k≥0
    filter!(
        k -> field(k) ≥ 0,
        candidates,    
    )
    
   
    # crystallographic_condition
    filter!(
        k -> divides(l,2*k*α,ring),
        candidates,    
    )
    
    # conjugates are pos-def and bounded length
    filter!(
        k -> all(≤(α*k^2,l,p) for p in P[2:end]),
        candidates,
    )
      
    return candidates
end


function next_min_k_for_l(vd,current_k,remaining_k,l)

    k_min,k_max = current_k,current_k+10

    while isempty(remaining_k)
        remaining_k = enumerate_k(vd,l,k_min,k_max)
        filter!(>(current_k),remaining_k)
        sort!(remaining_k)
        k_min,k_max = k_max,k_max+10
    end
    new_k = popfirst!(remaining_k)
    return (new_k,remaining_k)
end

function next_min_k_for_l!(vd,l)

    (current_k,remaining_k) = vd.least_k_by_root_length[l] 
    (new_k,new_remaining_k) = next_min_k_for_l(vd,current_k,remaining_k,l)
    vd.least_k_by_root_length[l] = (new_k,new_remaining_k)

end

function next_min_ratio(vd)
    field = vd.field
    min_pair = nothing

    val(k,l) = field(k^2)//field(l)
    

    for (l,(k,remaining)) in vd.least_k_by_root_length
        if isnothing(min_pair) || val(k,l) < val(min_pair...)
            min_pair = (k,l)
        end
    end
    return min_pair
end

function next_min_ratio!(vd)
    
    min_pair = next_min_ratio(vd)
    next_min_k_for_l!(vd,min_pair[2])
    return min_pair
end

zero_from(v,i) = all(==(0),v[i:end])

function extend_root_stem(vd::VinbergData,stem,root_length,bounds=nothing)
   
    j = length(stem) + 1
    
    if !isnothing(bounds) && any( bound < 0 && zero_from(root,j) for (root, bound) in zip(vd.accepted_roots,bounds))
        println("out of bounds!")
        return Vector{nf_elem}()
    end

    @assert isdiagonal(vd.gram_matrix)
    field = vd.field
    ring = vd.ring
    space = vd.quad_space
    d = [vd.gram_matrix[i,i] for i in 1:vd.dim]
    P = infinite_places(field)


    # stem = [k₀,…,k_j]

    l = root_length
    tab = "  "^j
    #@info tab * "extend_root_stem($stem, $root_length)"

    if j == vd.dim + 1

        #@info tab * "stem is complete"

        if is_root(space,ring,stem) && Hecke.inner_product(vd.quad_space,stem,stem) == l
            #@info tab * "and it's a root of length $l"
            return [stem]
        else
            #@info tab * "and it's bad (length is $(Hecke.inner_product(vd.quad_space,stem,stem)))"
            return Vector{nf_elem}()
        end
    else
        
        #@info tab * "stem is not complete"
        

        α = d[j]
        S_j = l - sum([d[i]*stem[i]^2 for i in 1:length(stem)]) 
        candidates_k_j = non_neg_short_t2_elems(vd.ring, 0,sum_at_places(field(S_j//α),1)+1) # TODO the +1 here is because of inexact computations --> can we make everything exact? --> yes in this case since sum_at_places ranges over all places, so probably just a trace or an exact t2 computation

        candidates_k_j = vcat(candidates_k_j, .- candidates_k_j)
        
        #@info tab * "candidates for j=$j are $candidates_k_j"

        #filter!(≥(0),candidates_k_j)
        #@info tab * "only pos            are $candidates_k_j"
        
        filter!(
            k-> all(≤(field(α*k^2),field(S_j),p) for p in P),
            candidates_k_j
        )
        #@info tab * "OK norm             are $candidates_k_j"
        filter!(
            k -> divides(l,2*k*α,ring),
            candidates_k_j,    
        )
        #@info tab * "OK crystal          are $candidates_k_j"
        bounds_updated(k) = if bounds === nothing
            nothing
        else
            [b - diag_coeff(vd,j)*k*r[j] for (b,r) in zip(bounds,vd.accepted_roots)]
        end
        return vcat([extend_root_stem(vd,vcat(stem,[k]),root_length,bounds_updated(k)) for k in candidates_k_j]...)
    end
    
end

function next_roots!(vd::VinbergData)
    (k,l) = next_min_ratio!(vd)
    @info "next_roots for l=$l and k = $k"
    return extend_root_stem(vd,[k],l,[-k*diag_coeff(vd,1)*r[1] for r in vd.accepted_roots])
end



function next_root!(vd::VinbergData)
    while isempty(vd.candidate_roots)
        vd.candidate_roots = [[c.elem_in_nf for c in r] for r in next_roots!(vd)]
    end
    return pop!(vd.candidate_roots)
end


function cone_roots!(vd::VinbergData)

    @assert vd.stage === Initialized "Too late to compute cone roots!"

    roots_at_distance_zero = []
    
    while true
        root = next_root!(vd)
        
        if root[1] ≠ 0
            pushfirst!(vd.candidate_roots,root)
            break
        end
        
        push!(roots_at_distance_zero,root)

    end
    
    # We put first the roots with integer coordinates to maximize the chance of having them in the cone roots
    # It's not necessary but easier to analyze the output and compare with rgug then
    integer_roots = filter!(r -> all(isinteger,vd.field.(r)), roots_at_distance_zero) 
    sort!(integer_roots)
    prepend!(integer_roots,roots_at_distance_zero)


    cone_roots = []
    @debug "starting with $(length(roots_at_distance_zero)) at dist zero"


    for r in roots_at_distance_zero
        @debug "looking at $r"
        if  all((-1)*r ≠ cr for cr in cone_roots)
            @debug "so far so good"
            if is_necessary_halfspace(cone_roots,-vd.gram_matrix.entries*r)
                @debug "degeneration"
                push!(cone_roots,r)
            end
        
        end
        @debug "have $(length(cone_roots)) cone roots" 
    end

    cone_roots = drop_redundant_halfspaces(cone_roots)
    @info "have $(length(cone_roots)) cone roots" 
    return cone_roots
end

function enumerate_roots!(vd)

    vd.accepted_roots = cone_roots!(vd)
    
    Coxeter_matrix = get_Coxeter_matrix(vd.quad_space, vd.ring, vd.accepted_roots) 
    diagram = build_diagram_and_subs(Coxeter_matrix,vd.dim-1)
    

    while length(vd.accepted_roots) < 100
        root = next_root!(vd)
        
        @info "Candidate $root"

        if Hecke.inner_product(vd.quad_space,basepoint(vd),root) ≤ 0 && all(Hecke.inner_product(vd.quad_space,prev,root) ≤  0 for prev in vd.accepted_roots)

            vd.stage = MoreRoots

            println("New root : ", root)

            extend!(diagram,[Coxeter_coeff(vd.quad_space, vd.ring, r,root) for r in vd.accepted_roots])
            push!(vd.accepted_roots,root)
            
            println("Matrix is")
            display(get_Coxeter_matrix(vd.quad_space, vd.ring, vd.accepted_roots))

            if is_finite_volume(diagram)
                vd.stage = Done
                @info "And the diagram has finite volume."
                break
            end

        end
        
    end
    println("Matrix is")
    display(get_Coxeter_matrix(vd.quad_space, vd.ring, vd.accepted_roots))
    
    return vd.accepted_roots

end



Qx, x = Hecke.QQ["x"]
f = x^2 - 2
K, a = Hecke.NumberField(f, "a"); @assert istotally_real(K)
M = matrix(K, 3, 3, [-1+a,0,0, 0,1,0, 0,0,1]) # The bilinear form

#=
r@fedora ~/U/v/r/A/build (master)> ./alvin -k=Q[sqrt 2] -qf -1-T,1,1
Quadratic form (2,1): -1 - 1 * T(2), 1, 1
Field of definition: Q[ sqrt(2) ]

Vectors: 
	e1 = (0, -1, 1)
	e2 = (0, 0, -1)
	e3 = (1, 1 + 1 * T(2), 0)
Algorithm ended

Graph written in the file: 
	output/2-1+T2,1,1.coxiter

Computation time: 0.0736732s
=# 

#
#
#########################
##
QK,Qa = Hecke.rationals_as_number_field()
QM = matrix(QK, 4, 4, [-1,0,0,0,  0,2,0,0,  0,0,6,0,  0,0,0,6])

