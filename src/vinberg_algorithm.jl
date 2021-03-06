using IterTools

include("util.jl")
include("diagrams.jl")
include("hyperbolic_lattices.jl")
include("root_enumeration.jl")


"""
Given the collection of roots `rr` and a root `r`, (and the corresponding partial products `rr_pp` and `r_pp`) tests whether the halfspace ``ℋ_r⁻`` is redundant for the convex cone 
```math
    ⋂_{r₀∈rr} ℋ_{r₀}⁻
```
(with ``ℋ_v^± = \\{u: ±u⊙v≥0\\}``).

""" 
function is_necessary_halfspace(
    rr#=::Vector{HyperbolicLatticeElement{n}}=#,
    rr_pp#=::Vector{SVector{n,Int}}=#,
    r#=::HyperbolicLatticeElement{n}=#,
    r_pp#=::SVector{n,Int}=#
) #where {n}
    
    #TODO: put signature back in a way that the compiler likes

    L = r.L
    @toggled_assert all(ro.L == r.L for ro in rr)
    @toggled_assert is_root(r)
    @toggled_assert L.G*r.vec == r_pp
    @toggled_assert all(is_root(ro) for ro in rr)
    @toggled_assert all(L.G*ro.vec == ro_pp for (ro,ro_pp) in zip(rr,rr_pp))

    return is_necessary_halfspace(rr_pp,r_pp)

end


"""
    drop_redundant_halfspaces(roots,roots_pp)

Given the collection of roots `roots` defining a cone and the corresponding partial products `roots_pp`, drop all which are redundant in the halfspace decomposition of the cone. 
"""
function drop_redundant_halfspaces(
    roots#=::Vector{HyperbolicLatticeElement{n}}=#,
    roots_pp#=::Vector{SVector{n,Int}}=#
)# where {n}
    
    #TODO: put signature back in a way that the compiler likes
    
    for i in reverse(collect(1:length(roots)))
       
        rr = copy(roots)
        r = popat!(rr,i)
        rr_pp = copy(roots_pp)
        r_pp = popat!(rr_pp,i)

        if ! is_necessary_halfspace(rr_pp,r_pp)
            return drop_redundant_halfspaces(rr,rr_pp) 
        end
    end
    
    return (roots,roots_pp)
    
end


"""
    roots_of_fundamental_cone(VL,roots_at_distance_zero)

Given the collection `roots_at_distance_zero` of roots at distance zero to `VL.v₀` (technically, roots with induced hyperplane at distance zero from ``v₀``), find a minimal subset defining a fundamental cone containing ``v₀``.
"""
function roots_of_fundamental_cone(
    VL::VinbergLattice,
    roots_at_distance_zero::Vector{HyperbolicLatticeElement}
)
    
    @toggled_assert (length(roots_at_distance_zero) == length(Set(roots_at_distance_zero))) "No root should appear more than once (sanity check)."

    # Store the candidate roots, and their partial product with G, for efficiency's sake, hopefully ("pp" stands for "partial product")
    candidate_roots = roots_at_distance_zero
    candidate_roots_pp = [VL.L.G * r.vec for r in candidate_roots]
    @info "Starting with $(length(candidate_roots)) candidate roots (at distance zero) for fundamental cone."

    # Start with an empty set of roots (and corresponding partial products) in the cone.
    cone_roots::Vector{HyperbolicLatticeElement} = Vector()
    cone_roots_pp::Vector{SVector{rk(VL.L),Int}}= Vector()

    # Iterate over all roots (and the precomputed corresponding partial product)
    for (r,r_pp) in zip(candidate_roots,candidate_roots_pp)
        @debug "Considering candidate root $r."

        # If the opposite of the root already is in the cone roots, adding this one would make the cone degenerate.
        # 
        # TODO: I also first added a “acute angle” condition, but:
        #
        # * B&P don't use it;
        # * It actually yields wrong results.
        # 
        # Why is that?
        if  all((-1)*r ≠ cr for cr in cone_roots) # && all(cr.vec' *r_pp < 0 for cr in cone_roots)
            @debug "Its opposite is not in the previous cone_roots."
           
            # Test that adding the halfspace ``ℋ_r⁻`` defined by `r` doesn't make the resulting cone degenerate (checked by `is_necessary_halfspace(cone_roots, -r`).
            # Indeed, the intersection becomes degenerate iff the intersection with ``ℋ_{-r}⁻ = ℋ_r⁺`` is /not/ strictly smaller than the original cone, which is iff ``-r`` defines a necessary halfspace.
            #if is_necessary_halfspace(cone_roots::Vector{HyperbolicLatticeElement{rk(VL.L)}},cone_roots_pp::Vector{SVector{rk(VL.L),Int}},-r::HyperbolicLatticeElement{rk(VL.L)},-r_pp::SVector{rk(VL.L),Int})
            if is_necessary_halfspace(cone_roots,cone_roots_pp,-r,-r_pp)
                @debug "And it does not make the cone degenerate."
                
                push!(cone_roots,r)
                push!(cone_roots_pp,r_pp)
            end
        
        end
        
    end
    
    (cone_roots,cone_roots_pp) = drop_redundant_halfspaces(cone_roots, cone_roots_pp)
    @info "Returning $(length(cone_roots)) cone roots."
    return (cone_roots, cone_roots_pp)

end


"""
    Vinberg_Algorithm(G [,v₀vec,rounds])

Run the Vinberg for `rounds` rounds (or indefinitely if `isnothing(rounds)`) 

* on the lattice given by the quadratic form `G`,
* with basepoint given by `v₀vec` (or chosen by diagonalization if `isnothing(v₀vec)`).
"""
function Vinberg_Algorithm(
    G::Array{Int,2};
    v₀_vec::Union{Array{Int,1},Nothing}=nothing,
    rounds=nothing
)
    VL = VinbergLattice(G;v₀_vec=v₀_vec)
    return Vinberg_Algorithm(VL;rounds=rounds)
end


"""
    Vinberg_Algorithm(VL::VinbergLattice[;rounds=nothing])

Run the Vinberg algorithm on the lattice ``VL.L`` with basepoint ``VL.v₀`` until `rounds` roots have been considered (if `rounds=nothing`, then the algorithm runs indefinitely long).
"""
function Vinberg_Algorithm(
    VL::VinbergLattice;
    rounds=nothing
)

    
    v₀ = VL.v₀
    G = VL.L.G
    n = rk(VL.L)
   
    # We make our iterator peekable so that we can look at the next root without consuming it: 
    new_roots_iterator = roots_by_distance(VL)
    @info "Initialized root iterator."

    roots_at_distance_zero = Vector{HyperbolicLatticeElement}()

    # Getting next (actually the first!) root
    while fake_dist(VL,(local new_root=new_roots_iterator())) == 0
        push!(roots_at_distance_zero,new_root)
    end
    # All roots at distance zero have now been considered
    # `new_root` is now the first root not at distance zero
    
    # Sort the roots, so as to have a predictable choice of fundamental cone
    sort!(roots_at_distance_zero)
    @info "Got all roots at distance zero."
    
    # Get all roots defining a fundamental cone for `v₀`
    # Note that there is some degree of freedom here since `v₀` can lie in different cones.
    # But once a cone is fixed, all the other roots appearing afterwards are uniquely defined.
    (roots, roots_pp) = roots_of_fundamental_cone(VL,roots_at_distance_zero)
    @info "Got the roots of a fundamental cone."
    for r in roots
        @info "                         : $r."
    end
    # Construct the corresponding `DiagramAndSubs` object containing all subdiagrams of the Coxeter diagram defined by the roots
    Coxeter_matrix = reduce(hcat,[[Coxeter_coeff(r₁,r₂) for r₁ in roots] for r₂ in roots])
    diagram = build_diagram_and_subs(Coxeter_matrix,n-1)
    @info "Built the corresponding Coxeter diagram."

    round_num = 0 
    for new_root in Iterators.flatten([[new_root],new_roots_iterator]) 

        !isnothing(rounds) && round_num ≥ rounds && break 
        round_num += 1

        @debug "Trying new root $new_root."
        
        @toggled_assert iff(
                            (all(r_pp' * new_root.vec ≤ 0 for r_pp in roots_pp) && times_v₀(VL,new_root) < 0),
                            (all((new_root⊙r) ≤ 0 for r in roots) && VL.v₀ ⊙ new_root < 0)
                           ) "Precomputed products should be the same as non-precomputed ones."
        
        if any( (r_pp' * new_root.vec > 0)::Bool for r_pp in roots_pp)
            # Acute angle condition ?
            @debug "It doesn't satisfy the acute angle condition; discarding it."
        
        elseif times_v₀(VL,new_root) ≥ 0
            # Normally, the way the code is written now makes it so that this condition is always satisfied: we only get roots containing v₀ in their correct halfspace

            # v₀ on the correct side of the halfspace ?
            @debug "The corresponding halfspace doesn't contain v₀; discarding it."
        
        else
            
            # new_root now satisfies all the requirements to be added to our set of roots.
            new_root_pp = G*new_root.vec
            extend!(diagram,[Coxeter_coeff(r,new_root) for r in roots])
            push!(roots,new_root)
            push!(roots_pp,new_root_pp)


            @info "Found new satisfying root: $new_root."
            if is_finite_volume(diagram)
                @info "And the diagram has finite volume."
                break
            end

        end
       
    end
   
    println("Decision ($(rounds)) :", is_finite_volume(diagram))
    @toggled_assert length(roots) == length(drop_redundant_halfspaces(roots,roots_pp)[1])
    return sort([r.vec for r in roots])

end


