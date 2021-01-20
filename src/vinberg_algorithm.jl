include("util.jl")
include("qsolve.jl")
include("diagrams.jl")
include("hyperbolic_lattices.jl")
include("root_enumeration.jl")

# Code adapted from N. V. Bogachev and A. Yu. Perepechko:
#
#   https://github.com/aperep/vinberg-algorithm


# TODO
#
# * See here for memoization: https://github.com/JuliaCollections/Memoize.jl
# * Extensive unit tests!
#

"""
    is_necessary_halfspace(rr,r[,rrpp,rpp])

Given the collection of roots `rr` and a root `r`, tests whether the halfspace ``ℋ_r⁻`` is redundant for thhe convex cone cone 
```math
    ⋂_{r₀∈rr} ℋ_{r₀}⁻
```
(with ``ℋ_v^± = \\{u: ±u⊙v≥0\\}``).

"""
function is_necessary_halfspace(rr::Vector{HyperbolicLatticeElement},r::HyperbolicLatticeElement)
    
    L = r.L
    @toggled_assert all(ro.L == r.L for ro in rr)
    @toggled_assert is_root(r)
    @toggled_assert all(is_root(ro) for ro in rr)

    return is_necessary_halfspace(Vector{SVector{rk(L),Int}}([ro.vec for ro in rr]),r.L.G,r.vec)

end

"""
Given the collection of roots `rr` and a root `r`, (and the corresponding partial products `rr_pp` and `r_pp`) tests whether the halfspace ``ℋ_r⁻`` is redundant for thhe convex cone cone 
```math
    ⋂_{r₀∈rr} ℋ_{r₀}⁻
```
(with ``ℋ_v^± = \\{u: ±u⊙v≥0\\}``).

""" 
function is_necessary_halfspace(
    rr#=::Vector{HyperbolicLatticeElement{n}}=#,
    rrpp#=::Vector{SVector{n,Int}}=#,
    r#=::HyperbolicLatticeElement{n}=#,
    rpp#=::SVector{n,Int}=#
) #where {n}
    
    #TODO: put signature back in a way that the compiler likes

    L = r.L
    @toggled_assert all(ro.L == r.L for ro in rr)
    @toggled_assert is_root(r)
    @toggled_assert L.G*r.vec == rpp
    @toggled_assert all(is_root(ro) for ro in rr)
    @toggled_assert all(L.G*ro.vec == ropp for (ro,ropp) in zip(rr,rrpp))

    return is_necessary_halfspace(rrpp,rpp)

end


"""
    drop_redundant_halfspaces(roots)

Given the collection of roots `roots` defining a cone, drop all which are redundant in the halfspace decomposition of the cone. 
"""
function drop_redundant_halfspaces(roots::Vector{HyperbolicLatticeElement})
    
    for i in reverse(collect(1:length(roots)))
       
        rr = copy(roots)
        r = popat!(rr,i)

        if ! is_necessary_halfspace(rr,r)
            return drop_redundant_halfspaces(rr) 
        end
    end
    
    return roots
    
end


"""
    drop_redundant_halfspaces(roots,roots_pp)

Given the collection of roots `roots` defining a cone and the corresponding partial products `roots_pp`, drop all which are redundant in the halfspace decomposition of the cone. 
"""
function drop_redundant_halfspaces(
    roots#=::Vector{HyperbolicLatticeElement{n}}=#,
    rootspp#=::Vector{SVector{n,Int}}=#
)# where {n}
    
    #TODO: put signature back in a way that the compiler likes
    
    for i in reverse(collect(1:length(roots)))
       
        rr = copy(roots)
        r = popat!(rr,i)
        rrpp = copy(rootspp)
        rpp = popat!(rrpp,i)

        if ! is_necessary_halfspace(rrpp,rpp)
            return drop_redundant_halfspaces(rr,rrpp) 
        end
    end
    
    return (roots,rootspp)
    
end

"""
    roots_of_fundamental_cone(VL,roots_at_distance_zero)

Given the collection `roots_at_distance_zero` of roots at distance zero to `VL.v₀` (technically, roots with induced hyperplane at distance zero from ``v₀``), find a minimal subset defining a fundamental cone containing ``v₀``.
"""
function roots_of_fundamental_cone(VL::VinbergLattice,roots_at_distance_zero::Vector{HyperbolicLatticeElement})
    
    @toggled_assert (length(roots_at_distance_zero) == length(Set(roots_at_distance_zero))) "No root should appear more than once."

    # Store the candidate roots, and their partial product with G, for efficiency's sake, hopefully ("pp" stands for "partial product")
    candidate_roots = roots_at_distance_zero
    candidate_roots_pp = [VL.L.G * r.vec for r in candidate_roots]
    @info "Found candidate roots for fundamental cone (numbering $(length(candidate_roots)))."

    # Start with an empty set of roots (and partial products) in the cone.
    cone_roots::Vector{HyperbolicLatticeElement} = Vector()
    cone_roots_pp::Vector{SVector{rk(VL.L),Int}}= Vector()

    # Iterate over all roots (and the precomputed corresponding partial product)
    for (r,r_pp) in zip(candidate_roots,candidate_roots_pp)
        @info "Considering candidate root $r."

        # If the opposite of the root already is in the cone roots, adding this one would make the cone degenerate.
        # 
        # TODO: I also first added a “acute angle” condition, but:
        #
        # * B&P don't use it;
        # * It actually yields wrong results.
        # 
        # Why is that?
        if  all((-1)*r ≠ cr for cr in cone_roots) #  && all(cr⊙r ≤ 0 for cr in cone_roots)
            @info "Its opposite is not in the previous cone_roots."
           
            # Test that adding the halfspace ``ℋ_r⁻`` defined by `r` doesn't make the resulting cone degenerate (checked by `is_necessary_halfspace(cone_roots, -r`).
            # Indeed, the intersection becomes degenerate iff the intersection with ``ℋ_{-r}⁻ = ℋ_r⁺`` is /not/ strictly smaller than the original cone, which is iff ``-r`` defines a necessary halfspace.
            #if is_necessary_halfspace(cone_roots::Vector{HyperbolicLatticeElement{rk(VL.L)}},cone_roots_pp::Vector{SVector{rk(VL.L),Int}},-r::HyperbolicLatticeElement{rk(VL.L)},-r_pp::SVector{rk(VL.L),Int})
            if is_necessary_halfspace(cone_roots,cone_roots_pp,-r,-r_pp)
                @info "And it does not make the cone degenerate."
                
                push!(cone_roots,r)
                push!(cone_roots_pp,r_pp)

                @info "Dropping now redundant halfspaces."
                (cone_roots,cone_roots_pp) = drop_redundant_halfspaces(cone_roots, cone_roots_pp)
            end
        
        end
        
    end
    
    #cone_roots = drop_redundant_halfspaces(cone_roots)
    @info "Dropping rendundant halfspaces; returning cone roots (amounting to: $(length(cone_roots)))."
    return (cone_roots, cone_roots_pp)

end

function is_finite_volume(roots::Array{HyperbolicLatticeElement,(1)},VL::VinbergLattice)
    Gram_matrix = Int.(reduce(hcat,[[r1⊙r2 for r1 in roots] for r2 in roots]))
    Coxeter_matrix = Gram_to_Coxeter(Gram_matrix)
    return isnothing(Coxeter_matrix) ? false : is_fin_vol(Coxeter_matrix,rk(VL.L)-1)
    
end


function Vinberg_Algorithm(G::Array{Int,2};v₀_vec::Union{Array{Int,1},Nothing}=nothing,rounds=nothing)
    VL = VinbergLattice(G;v₀_vec=v₀_vec)
    return Vinberg_Algorithm(VL;rounds=rounds)
end

"""
    Vinberg_Algorithm(VL::VinbergLattice[;rounds=nothing])

Run the Vinberg algorithm on the lattice `VL.L` with basepoint `VL.v₀` until `rounds` roots have been considered (if `rounds=nothing`, then the algorithm runs indefinitely long).
"""
function Vinberg_Algorithm(VL::VinbergLattice;rounds=nothing)

    
    # `r` counts the numer of remaining rounds the algorithm is allowed to use.
    # By convention here, nothing can be understood as ∞, so that:
    # 
    # * decrease(∞) = ∞;
    # * ∞ ≥ 0 is true.
    @inline decrease(r) = (r === nothing ? nothing : r-1)
    @inline geqzero(r) = (r=== nothing ? true : r ≥ 0)

    v₀ = VL.v₀
    G = VL.L.G
    
    # initiate our root enumeration object
    new_roots_iterator = RootsByDistance(VL)
    @info "Initialized root iterator."

    # get all roots at distance zero (corresponding to hyperplanes containing `v₀`)
    roots_at_distance_zero::Vector{HyperbolicLatticeElement} = next_at_distance_zero!(new_roots_iterator)
    # and sort them (this is not necessary for the algorithm, but ensures a level of predictability in the output)
    sort!(roots_at_distance_zero)
    @info "Got all roots at distance zero."
    
    # Get all roots defining a fundamental cone for `v₀`
    # Note that there is some degree of freedom here since `v₀` can lie in different cones.
    # But once a cone is fixed, all the other roots appearing afterwards are uniquely defined.
    (roots, roots_pp) = roots_of_fundamental_cone(VL,roots_at_distance_zero)
    @info "Got the roots of a fundamental cone."

    
    start = true

    new_root = next!(new_roots_iterator)
    @info "Trying with the root $new_root."
    
    while start || ( ! is_finite_volume(roots,VL) && geqzero(rounds))
        start = false
        
        

        @toggled_assert iff(
                            (all(r_pp' * new_root.vec ≤ 0 for r_pp in roots_pp) && times_v₀(VL,new_root) < 0),
                            (all((new_root⊙r) ≤ 0 for r in roots) && VL.v₀ ⊙ new_root < 0)
                           )
        
        while !(all(r_pp' * new_root.vec ≤ 0 for r_pp in roots_pp) && times_v₀(VL,new_root) < 0) && geqzero(rounds) 
            rounds = decrease(rounds)
            @info "The root either does not face away from v₀ or doesn't satisfy the accute angle condition."

            new_root = next!(new_roots_iterator)
            @info "Trying with the root $new_root."
        end

        new_root_pp = G*new_root.vec
        @info "Testing whether the root defines a necessary hyperplane."
        if is_necessary_halfspace(roots, roots_pp, new_root, new_root_pp)
            # seems like this is always satisfied?
            #
            push!(roots,new_root)
            push!(roots_pp,new_root_pp)
            @info "It does; adding it to our collection."
        end

        # TODO: Here instead of only checking whether new_root defines a necessary hyperplane, we should call 
        #roots = drop_redundant_halfspaces(roots)
        # but this implies that we should also correct partial_times
        # Hence: that's a TODO!

    end
   
    println("Decision ($(rounds)) :", is_finite_volume(roots,VL))
    println("can we drop hyperplanes? $(length(roots)) vs $(length(drop_redundant_halfspaces(roots)))")
    return [r.vec for r in roots]

end


