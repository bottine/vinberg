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
    is_necessary_halfspace(rr,r)

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
    drop_redundant_halfspaces(roots)

Given the collection of roots `roots` defining a cone, drop all which are redundant in the halfspace decomposition of the cone. 
"""
function drop_redundant_halfspaces(roots::Vector{HyperbolicLatticeElement})
    
    @toggled_assert all(is_root(r) for r in roots)

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
    roots_of_fundamental_cone(VL,roots_at_distance_zero)

Given the collection `roots_at_distance_zero` of roots at distance zero to `VL.v₀` (technically, roots with induced hyperplane at distance zero from ``v₀``), find a minimal subset defining a fundamental cone containing ``v₀``.
"""
function roots_of_fundamental_cone(VL::VinbergLattice,roots_at_distance_zero::Vector{HyperbolicLatticeElement})
    
    @toggled_assert (length(roots_at_distance_zero) == length(Set(roots_at_distance_zero))) "No root should appear more than once."

    candidate_roots = roots_at_distance_zero
    @info "Found candidate roots for fundamental cone (numbering $(length(candidate_roots)))."

    # Start with an empty set of roots in the cone.
    cone_roots::Vector{HyperbolicLatticeElement} = Vector()

    # Iterate over all roots (and the precomputed corresponding partial product)
    for r in candidate_roots
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
           
            # Test that the halfspace defined by `r` cuts the cone into two non-degenerate cones, i.e. 
            # * `r` is necessary (checked by `is_necessary_halfspace(cone_roots, r)`, and,
            # * adding the halfspace ``ℋ_r⁻`` doesn't make the resulting cone degenerate (checked by `is_necessary_halfspace(cone_roots, -r`) 
            if is_necessary_halfspace(cone_roots, r) && is_necessary_halfspace(cone_roots, -r)
                @info "And it cuts the cone into two non-degenerate parts; adding it."
                
                push!(cone_roots,r)

                @info "Dropping now redundant halfspaces."
                cone_roots = drop_redundant_halfspaces(cone_roots)
            end
        
        end
        
    end
    
    cone_roots = drop_redundant_halfspaces(cone_roots)
    @info "Dropping rendundant halfspaces; returning cone roots (amounting to: $(length(cone_roots)))."
    return cone_roots

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
    roots::Array{HyperbolicLatticeElement,(1)} = roots_of_fundamental_cone(VL,roots_at_distance_zero)
    # Store the result of r.vec' * G for all roots in the fundamental cone for performance gains (hopefully)
    partial_times = [r.vec' * G for r in roots]
    @info "Got the roots of a fundamental cone."

    
    start = true

    new_root = next!(new_roots_iterator)
    @info "Trying with the root $new_root."
    
    while start || ( ! is_finite_volume(roots,VL) && geqzero(rounds))
        start = false
        
        

        @toggled_assert iff(
                            (all(pt * new_root.vec ≤ 0 for pt in partial_times) && times_v₀(VL,new_root) < 0),
                            (all((new_root⊙r) ≤ 0 for r in roots) && VL.v₀ ⊙ new_root < 0)
                           )
        
        while !(all(pt * new_root.vec ≤ 0 for pt in partial_times) && times_v₀(VL,new_root) < 0) && geqzero(rounds) 
            rounds = decrease(rounds)
            @info "The root either does not face away from v₀ or doesn't satisfy the accute angle condition."

            new_root = next!(new_roots_iterator)
            @info "Trying with the root $new_root."
        end

        @info "Testing whether the root defines a necessary hyperplane."
        if is_necessary_halfspace(roots,new_root)
            # seems like this is always satisfied?
            #
            push!(roots,new_root)
            push!(partial_times,new_root.vec' * G)
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


