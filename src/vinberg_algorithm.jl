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
    rr::Vector{HyperbolicLatticeElement{n}},
    r::HyperbolicLatticeElement{n},
)::Bool where {n}
    
    #TODO: put signature back in a way that the compiler likes

    lat = r.lat
    @toggled_assert all(ro.lat == r.lat for ro in rr)
    @toggled_assert is_root(r)
    @toggled_assert all(is_root(ro) for ro in rr)
    
    return is_necessary_halfspace(lat.D,[ro.diag_coordinates for ro in rr],r.diag_coordinates)

end


"""
    drop_redundant_halfspaces(roots,roots_pp)

Given the collection of roots `roots` defining a cone and the corresponding partial products `roots_pp`, drop all which are redundant in the halfspace decomposition of the cone. 
"""
function drop_redundant_halfspaces(
    roots#=::Vector{HyperbolicLatticeElement{n}}=#,
)# where {n}
    
    #TODO: put signature back in a way that the compiler likes
    
    for i in length(roots):-1:1
       
        rr = copy(roots)
        r = popat!(rr,i)

        if ! is_necessary_halfspace(rr,r)
            return drop_redundant_halfspaces(rr) 
        end
    end
    
    return roots
    
end


"""
    roots_of_fundamental_cone(lat,roots_at_distance_zero)

Given the collection `roots_at_distance_zero` of roots at distance zero to `lat.v₀` (technically, roots with induced hyperplane at distance zero from ``v₀``), find a minimal subset defining a fundamental cone containing ``v₀``.
"""
function roots_of_fundamental_cone(
    lat::HyperbolicLattice{n},
    roots_at_distance_zero::Vector{HyperbolicLatticeElement{n}}
)::Vector{HyperbolicLatticeElement{n}} where {n}
    
    @toggled_assert (length(roots_at_distance_zero) == length(Set(roots_at_distance_zero))) "No root should appear more than once (sanity check)."

    # Store the candidate roots, and their partial product with G, for efficiency's sake, hopefully ("pp" stands for "partial product")
    candidate_roots = roots_at_distance_zero
    @info "Starting with $(length(candidate_roots)) candidate roots (at distance zero) for fundamental cone."

    # Start with an empty set of roots (and corresponding partial products) in the cone.
    cone_roots::Vector{HyperbolicLatticeElement{n}} = Vector()

    # Iterate over all roots (and the precomputed corresponding partial product)
    for r in candidate_roots
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
            #if is_necessary_halfspace(cone_roots::Vector{HyperbolicLatticeElement{rk(lat)}},cone_roots_pp::Vector{SVector{rk(lat),Int}},-r::HyperbolicLatticeElement{rk(lat)},-r_pp::SVector{rk(lat),Int})
            if is_necessary_halfspace(cone_roots,-r)
                @debug "And it does not make the cone degenerate."
                
                push!(cone_roots,r)
            end
        
        end
        
    end
    
    cone_roots = drop_redundant_halfspaces(cone_roots)
    @info "Returning $(length(cone_roots)) cone roots."
    return cone_roots

end

"""
    Vinberg_Algorithm(lat::HyperbolicLattice[;rounds=nothing])

Run the Vinberg algorithm on the lattice ``lat`` with basepoint ``lat.v₀`` until `rounds` roots have been considered (if `rounds=nothing`, then the algorithm runs indefinitely long).
"""
function Vinberg_Algorithm(
    lat::HyperbolicLattice{n};
    rounds=nothing
)::Vector{HyperbolicLatticeElement{n}} where {n}

    
    G = lat.G
   
    # We make our iterator peekable so that we can look at the next root without consuming it: 
    @info "Getting roots at distance zero"
    roots_at_distance_zero = Vector{HyperbolicLatticeElement{n}}([r for r in roots_by_distance(lat,x->true,==(0))])
    # Sort the roots, so as to have a predictable choice of fundamental cone
    sort!(roots_at_distance_zero)
    @info "Got all roots at distance zero."
    
    # Get all roots defining a fundamental cone for `v₀`
    # Note that there is some degree of freedom here since `v₀` can lie in different cones.
    # But once a cone is fixed, all the other roots appearing afterwards are uniquely defined.
    roots::Vector{HyperbolicLatticeElement{n}} = roots_of_fundamental_cone(lat,roots_at_distance_zero)
    @info "Got the roots of a fundamental cone."
    for r in roots
        @info "                         : $r."
    end
    # Construct the corresponding `DiagramAndSubs` object containing all subdiagrams of the Coxeter diagram defined by the roots
    Coxeter_matrix = get_Coxeter_matrix(roots) 
    diagram = build_diagram_and_subs(Coxeter_matrix,n-1)
    @info "Built the corresponding Coxeter diagram."

    round_num = 0

    roots_by_distance_iter = roots_by_distance(lat,>(0),x->true,roots)

    new_roots = nothing # those are fed to the iterator 
    while true
   
        new_root = roots_by_distance_iter(new_roots) 
        new_roots = nothing # those are fed to the iterator 

        !isnothing(rounds) && round_num ≥ rounds && break 
        round_num += 1


        @debug "Trying new root $new_root."
        
        if any( (r ⊙ new_root > 0)::Bool for r in roots)
            # Acute angle condition ?
            @debug "It doesn't satisfy the acute angle condition; discarding it."
        
        elseif times_v₀(lat,new_root) ≥ 0
            # Normally, the way the code is written now makes it so that this condition is always satisfied: we only get roots containing v₀ in their correct halfspace

            # v₀ on the correct side of the halfspace ?
            @debug "The corresponding halfspace doesn't contain v₀; discarding it."
        
        else
            
            # new_root now satisfies all the requirements to be added to our set of roots.
            extend!(diagram,[Coxeter_coeff(r,new_root) for r in roots])
            push!(roots,new_root)
            new_roots = [new_root]

            @info "Found new satisfying root: $new_root."
            if is_finite_volume(diagram)
                @info "And the diagram has finite volume."
                break
            end

        end
       
    end
   
    println("Decision ($(rounds)) :", is_finite_volume(diagram))
    @toggled_assert length(roots) == length(drop_redundant_halfspaces(roots))
    return roots

end

function Vinberg_Algorithm(
    G,
    rounds=nothing
) where {n}
    return Vinberg_Algorithm(HyperbolicLattice(G),rounds=rounds)
end


