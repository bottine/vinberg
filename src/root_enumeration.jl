using DataStructures


include("util.jl")
include("qsolve.jl")
include("hyperbolic_lattices.jl")

"""
Given a `VinbergLattice` with underlying lattice ``L``, basepoint ``v₀``, orthogonal complement ``V₁=v₀^⟂`` and coset representatives ``W`` for the sublattice ``⟨v₀⟩⊕V₁``, any element ``u`` of ``L`` can be uniquely represented as 
```math
    u = a₀v₀ + v₁ + w
``` 
with:

* ``a₀∈ℤ``,
* ``v₁∈V₁``,
* ``w∈W``.


The fake distance between ``v₀`` and the hyperplane ``H_u`` is then ``(u⊙v₀)²/(u⊙u) = ((a₀v₀ + v₁ + w)⊙v₀)²/(u⊙u) = (a₀v₀⊙v₀ + w⊙v₀)²/(u⊙u)``.
We call "Decomposition Pattern" the data of ``a₀,w`` _and_ ``k=u⊙u``; it follows that the fake distance is therefore determined by the decomposition pattern (we don't need to know what exactly ``u`` is, but only its norm ``k`` and its components ``a₀v₀`` and ``k``).

The type `RootDecompositionPattern` is just a synonym for the tuple `(k,w,a₀)`. 
"""
RootDecompositionPattern = Tuple{Int,                       # = k
                                 HyperbolicLatticeElement,  # = w
                                 Int                        # = a
                                }

"""
`Pattern` is a shorthand for `RootDecompositionPattern`
"""
Pattern = RootDecompositionPattern

"""
Computes the fake distance between the hyperplane between a root matching the root decomposition pattern `p` and `VL.v₀`.
"""
function fake_dist(p::RootDecompositionPattern,VL::VinbergLattice)
    (k,w,a) = p
    (a*VL.v₀_norm + (times_v₀(VL,w)))^2//k
end

"""
Recall that if ``u = a₀v₀ + v₁ + w`` is a root with ``k=u⊙u``, the distance from ``v₀`` to its hyperplane ``H_u`` is entirely defined by ``(k,w,a)``.
Thus, to be able to iterate over all roots by increasing distance, we first iterate over all patterns by distance.
Note that:

* The set ``W`` of coset representatives for ``⟨v₀⟩⊕V₁`` is finite.
* The set of norms ``k=u⊙u`` is also finite since the norm of a root must divide the last invariant factor.
* ``a₀`` can take any integer value.

Furthermore, for a fixed pair ``(k,w)`` minimizing the distance amounts to minimizing the quadratic (in ``a₀``) function
```math
    (u⊙v₀)²/(u⊙u) = (a₀v₀⊙v₀ + v₀⊙w)²/k.
```
The minimum of this function is attained at ``a = -v₀⊙w/v₀⊙v₀``, and its values increases with ``a`` getting farther from this minimum in both direction.

It follows that for a given fixed pair ``(k,w)`` iterating over the possible ``a₀``s by increasing distance amounts to finding this minimum, and keeping track of both possible directions for ``a`` to go.
Since the number of pairs ``(k,w)`` is finite, we simply keep, for each such pair, the two "iterators" ``a₋`` (decreasing from ``a₀``) and ``a₊`` (increasing).
"""
mutable struct RootDecompositionPatternsByDistance
    # The VinbergLattice in which we're working
    VL::VinbergLattice
    # For each pair `(w,k)` the "iterators" a₋ and a₊
    next_least_as_for_w_and_k::Union{
                                     Nothing,  # In case we haven't started iterating
                                     SortedDict{
                                                Tuple{HyperbolicLatticeElement,Int} # Keys are (w,k)
                                                ,Tuple{Int,Int}}}                   # Values are (a₋,a₊)
    # We also keep a copy of the current fake_dist, for verification purposes
    current_fake_dist::Union{Nothing,Rational{Int}}
end

"""
Constructs an instance of `RootDecompositionPatternsByDistance` by iterating over pairs `(k,w)`, and for each such pair:

* Finding the ``a₀`` (non-necessarily integer) attaining the minimum distance for the pair ``(k,w)``.
* Taking its `ceil` and `floor` to define `a₊` and `a₋` respectively (the "iterators").
"""
function RootDecompositionPatternsByDistance(VL::VinbergLattice)
    
    v₀ = VL.v₀
    W = VL.W
    
    # This is the fake dist to be minimized
    # Only define it for verification purposes
    # We can't just use `fake_dist(::Pattern,::VinbergLattice)` because we use it on `a₀`, not necessarily an integer.
    f_dist(k,w,a) = (a*VL.v₀_norm + (times_v₀(VL,w)))^2//k

    next_least_as_for_w_and_k = SortedDict{Tuple{HyperbolicLatticeElement,Int},Tuple{Int,Int}}()
     
    for w in W, k in root_lengths(VL.L)
        # Find the minimal ``a₀``
        a₀::Rational{Int} = -times_v₀(VL,w)//VL.v₀_norm; 
        @toggled_assert -(w⊙v₀)//(v₀⊙v₀) == -times_v₀(VL,w)//VL.v₀_norm "The optimized computation should equal the full one."
        
        # its integer approximations in both directions 
        a₋::Int = floor(a₀)
        a₊::Int = ceil(a₀)

        # verify that they are actually not as good
        @toggled_assert f_dist(k,w,a₋) ≥ f_dist(k,w,a₀)
        @toggled_assert f_dist(k,w,a₊) ≥ f_dist(k,w,a₀)
        @debug "For (w=$w,k=$k), (a₋,a₀,a₊) is ($a₋,$a₀,$a₊) with values $((f_dist(k,w,a₋),f_dist(k,w,a₀),f_dist(k,w,a₊)))."
        
        # push the result
        push!(next_least_as_for_w_and_k, (w,k) => (a₋,a₊))
    end

    current_fake_dist::Union{Nothing,Rational{Int}} = nothing

    return RootDecompositionPatternsByDistance(VL,next_least_as_for_w_and_k,current_fake_dist)


end

"""
    next!(pats::RootDecompositionPattern)

Returns:

* A `Vector` of `Pattern`s, all at the same distance.
* The distance.

The function proceeds as follows:
* Iterate over all `(w,k)`, and for each get the "iterators" `a₊,a₋`.
* Find the minimum distance over all those.
* Collect the `Pattern`s `(w,k,a)` (with `a` either `a₊` or `a₋`) attaining this minimum
* Return all of those, and the distance.
"""
function next!(pats::RootDecompositionPatternsByDistance)::Tuple{Vector{RootDecompositionPattern},Rational{Int}}
    
    v₀ = pats.VL.v₀
    VL = pats.VL
    
    # The fake distance, redefined here so that we don't have to give `VL` explicitely at each call.
    f_dist(k,w,a) = (a*VL.v₀_norm + (times_v₀(VL,w)))^2//k
        
    # No min_val yet, no patterns either
    min_val = nothing
    min_patterns=Vector{Pattern}()
    for  ((w,k),(a₋,a₊)) in pats.next_least_as_for_w_and_k
        
        # Verify that everything we see now is farther away than the previous fake_dist.
        @toggled_assert pats.current_fake_dist === nothing || f_dist(k,w,a₊) ≥ pats.current_fake_dist 
        @toggled_assert pats.current_fake_dist === nothing || f_dist(k,w,a₋) ≥ pats.current_fake_dist 
      
        # We first treat `a₊`, then `a₋`.

        # If `min_val` not set yet, or set but we are lower than it, push our pattern to `min_patterns`.
        # Note that `min_val === nothing` iff `min_patterns == []`.
        if min_val === nothing || f_dist(k,w,a₊) < min_val
            min_patterns = Vector{Pattern}([(k,w,a₊)])
            min_val = f_dist(k,w,a₊)

        # If we're exactly at `min_val`, push our `Pattern`
        elseif f_dist(k,w,a₊) == min_val
            push!(min_patterns,(k,w,a₊))
        end

        # Same as for `a₊`
        if min_val === nothing || f_dist(k,w,a₋) < min_val
            min_patterns = Vector{Pattern}([(k,w,a₋)])
            min_val = f_dist(k,w,a₋)
        elseif f_dist(k,w,a₋) == min_val
            # Here we differ from `a₊`:
            # It may happen that a₋ == a₊ (only at the very beginning, i.e. distance zero from v₀).
            # In this case, if a₊ has already been added, we shouldn't not add it again as `a₋`.
            # So, we ensure that a₋ is not added if it was already added as a₊.
            if isempty(min_patterns) || min_patterns[end] ≠ (k,w,a₋)
                push!(min_patterns,(k,w,a₋))
            end
        end
    end
   
    # Store the new minimal distance.
    pats.current_fake_dist = min_val

    # For each `a` appearing in a pattern, "discard it" by advancing the corresponding iterator `a₊` or `a⁻`.
    for (k,w,a) in min_patterns
        (a₋,a₊) = pats.next_least_as_for_w_and_k[(w,k)]
        if a₊ == a
            a₊ += 1
        end
        if a₋ == a
            a₋ -= 1
        end
        pats.next_least_as_for_w_and_k[(w,k)] = (a₋,a₊)
    end
    
    # Return the patterns, and their distance
    return (min_patterns,min_val)

end

"""
    Allows iterating over roots by batches corresponding either to a common distance, or a common pattern.
"""
mutable struct BatchedRoots
    # The pattern iterator
    pats::RootDecompositionPatternsByDistance
    # The patterns we have "in store"
    current_patterns::Union{Nothing,Vector{Pattern}}
end

function BatchedRoots(VL::VinbergLattice)
    pats = RootDecompositionPatternsByDistance(VL)
    current_patterns = nothing
    return BatchedRoots(pats,current_patterns)
end

"""
    Returns the next batch of roots all, at the same (next least) distance.
    It may be that the next least distance actually contains no root, in which case we return an empty `Vector`.
"""
function next_dist_batch!(br::BatchedRoots)
    
    VL = br.pats.VL
    v₀ = VL.v₀
    
    # If there is no pattern in store:
    if br.current_patterns === nothing || isempty(br.current_patterns)
        # then, ask for the next ones from our `RootDecompositionPatternsByDistance`.
        (br.current_patterns,dist) = next!(br.pats)
    end
    
    # Find all roots corresponding to all patterns in store
    roots_batches = [roots_decomposed_into(VL,a*v₀ + w,k) for (k,w,a) in br.current_patterns]
    this_batch =  vcat(roots_batches...)

    # Empty our pattern store
    br.current_patterns = []
    
    # Ensure that the distances for our roots are what br.pats thinks they are
    @toggled_assert all(fake_dist(VL,r) == br.pats.current_fake_dist for r in this_batch)
   
    return this_batch

end

"""
    Returns the next batch of roots, all corresponding to the same (next least distance) pattern. 
    It may be that the pattern actually contains no root, in which case we return an empty `Vector`.
"""
function next_pattern_batch!(br::BatchedRoots)

    VL = br.pats.VL
    v₀ = VL.v₀

    # If there is no pattern in store:
    if br.current_patterns === nothing || isempty(br.current_patterns)
        # then, ask for the next ones from our `RootDecompositionPatternsByDistance`.
        (br.current_patterns, dist) = next!(br.pats)
    end
    
    # Pop the next pattern
    (k,w,a) = pop!(br.current_patterns)
    
    # Find the corresponding roots
    this_batch = roots_decomposed_into(VL,a*v₀ + w,k)
    
    # Ensure that the distances for our roots are what br.pats thinks they are
    @toggled_assert all(fake_dist(VL,r) == br.pats.current_fake_dist for r in this_batch)
    
    return this_batch 

end

"""
Allows iterating over roots (with increasing distance) with the following possibilities:

* Either get all roots at distance zero (can only be done on a fresh instance)
* Get one next root.
* Get some "natural" batch (i.e. what `BatchedRoots` spits out with `next_pattern_batch!`).
"""
mutable struct RootsByDistance
    # Batched roots iterator
    br::BatchedRoots
    # What roots we have asked `br` but haven spitted out ourselves yet (i.e. a buffer)
    remaining::Vector{HyperbolicLatticeElement}
    # Whether this is a fresh instance
    fresh::Bool
end

function RootsByDistance(VL::VinbergLattice)
    return RootsByDistance(BatchedRoots(VL),[],true)
end

"""
    All roots at distance zero (that is, such that ``v₀`` is contained in the corresponding hyperplane).
"""
function next_at_distance_zero!(rbd::RootsByDistance) 
    @toggled_assert rbd.fresh == true "Only allowed to get roots at distance zero on a fresh instance."
    rbd.fresh = false

    # Since this is a fresh instance, we haven't gotten anything yet from `rbd.br`, so that `next_dist_batch!` returns all roots at the next possible distance, that is those at distance zero.
    return next_dist_batch!(rbd.br)
end

"""
Some (non-empty) batch of roots.
"""
function next_batch!(rbd::RootsByDistance)
    rbd.fresh = false
    
    # If nothing in our buffer, fill it again (and try as long as `next_pattern_batch` gives us nothing)
    while isempty(rbd.remaining)
        rbd.remaining = next_pattern_batch!(rbd.br)
    end

    # Return the content of rbd.remaining, and empty it.
    batch = rbd.remaining
    rbd.remaining = []
    return batch
end

"""
Exactly one root.
"""
function next!(rbd::RootsByDistance)
    rbd.fresh = false
    
    # If nothing in our buffer, fill it again (and try as long as `next_pattern_batch` gives us nothing)
    while isempty(rbd.remaining)
        rbd.remaining = next_pattern_batch!(rbd.br)
    end
    
    # The buffer is not empty: pop one root out of it.
    return pop!(rbd.remaining)
end

"""
    roots_decomposed_into(VL,a,k)

Returns all roots ``r`` of ``L`` that can be written as 
```math
    r = a + v₁
```
and satisfying ``r⊙r = k`` with ``v₁∈V₁``.

# Strategy

We first look for arbitrary elements ``r`` satisfying these conditions, and then filter on roots (using `is_root`).
If ``r = a + v₁``, and ``r⊙r=k``, by unrolling, we get
```math
    v₁⊙v₁ + 2v₁⊙a + a⊙a = k,
```
that is
```math
    v₁'Gv₁ + 2v₁'Ga + a'Ga = k
```
with ``v₁∈V₁``.
Recall that ``V₁ = v₀^⟂``.
We then transform the equation by taking a basis for ``V₁``:
Let ``M₁`` be the matrix whose columns are the basis elements for ``V₁``, and ``b₁`` the representation of ``v₁`` under this basis, so that ``v₁ = M₁b₁``.
Plugging this in the equation above yields:
```math
    b₁'M₁'GM₁b₁ + 2b₁'M₁Ga + a'Ga = k,
```
and now ``M₁'GM₁`` is positive definite (since ``V₁`` is the orthogonal complement of ``v₀`` and ``G`` is of signature ``(n,1)``).
This results in a "positive definite" quadratic equation for ``b₁``, which we know how to solve. 
It suffices then to transform back
"""
function roots_decomposed_into(VL::VinbergLattice, a::HyperbolicLatticeElement, k::Int)
    
    r = rk(VL.L)
    M₁ = VL.V₁_basis_matrix
    
    # Translate our problem into a positive definite quadratic equation, solvable by `qsolve()`
    A::SMatrix{r-1,r-1,Int} = M₁' * VL.L.G * M₁
    b::SVector{r-1,Int} = 2 *( M₁' * VL.L.G * a.vec)
    γ::Int = a⊙a - k
    
    # Finds solutions
    solutions = qsolve(A, b, γ)
    
    # Translate them back
    #solutions_in_L::Array{HyperbolicLatticeElement,1} = (b -> VL.L(M₁ * b + a.vec)).(solutions) # This is prettier, but performance-wise horrible. TODO: investigate
    solutions_in_L::Array{HyperbolicLatticeElement,1} = (b -> HyperbolicLatticeElement(VL.L,M₁ * b + a.vec)).(solutions)
    @toggled_assert all(norm(u) == k for u in solutions_in_L) "``u⊙u`` must be equal to k"
    @toggled_assert all((u-a)⊙VL.v₀ == 0 for u in solutions_in_L) "``u-a`` must lie in `V₁`"


    # Filter to get roots only (we know that `k` is the norm of `u`, so feed it to `is_root` already)
    return filter(u->is_root(u,k),solutions_in_L)
end


