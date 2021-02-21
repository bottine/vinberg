using DataStructures
using ResumableFunctions

include("util.jl")
include("qsolve_diag_constrained.jl")
include("hyperbolic_lattices.jl")

"""
Given a `HyperbolicLattice` with underlying lattice ``L``, basepoint ``v₀``, orthogonal complement ``V₁=v₀^⟂`` and coset representatives ``W`` for the sublattice ``⟨v₀⟩⊕V₁``, any element ``u`` of ``L`` can be uniquely represented as 
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
Computes the fake distance between the hyperplane between a root matching the root decomposition pattern `p` and `v₀(lat)`.
"""
function fake_dist(p::RootDecompositionPattern,lat::HyperbolicLattice)
    (k,w,a) = p
    (a*v₀_norm(lat) + (times_v₀(lat,w)))^2//k
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
The minimum is exactly achieved when ``v₀`` is in the halfspace defined by ``u``, and with ``a₀`` getting smaller, ``v₀`` gets in the negative halfspace, while with ``a₀`` getting bigger, ``v₀`` gets in the positive halfspace.

It follows that for a given fixed pair ``(k,w)`` iterating over the possible ``a₀``s by increasing distance amounts to finding this minimum, and keeping track of both possible directions for ``a`` to go.
Since the number of pairs ``(k,w)`` is finite, we can simply keep, for each such pair, the two "iterators" ``a₋`` (decreasing from ``a₀``) and ``a₊`` (increasing).
Note that if ``a₀`` is not an integer, we'll need to take the closest one, obviously.

**SINCE** we only care about roots containing ``v₀`` in their **positive** halfspace, we can simply discard the ``a₋`` part, and only iterate over ``a₊``.
That is, we iterate by distance to `v₀`, but only one "one side"
"""
mutable struct RootDecompositionPatternsByDistance
    # The HyperbolicLattice in which we're working
    lat::HyperbolicLattice
    # For each pair `(w,k)` the "iterator" a₊
    next_least_a_for_w_and_k::Union{
                                     Nothing,  # In case we haven't started iterating
                                     SortedDict{
                                                Tuple{HyperbolicLatticeElement,Int} # Keys are (w,k)
                                                ,Int}}                   # Value is a₊
    # We also keep a copy of the current fake_dist, for verification purposes
    current_fake_dist::Union{Nothing,Rational{Int}}
end

"""
Constructs an instance of `RootDecompositionPatternsByDistance` by iterating over pairs `(k,w)`, and for each such pair:

* Finding the ``a₀`` (non-necessarily integer) attaining the minimum distance for the pair ``(k,w)``.
* Taking its `ceil`to define `a₊`(the "iterator").
"""
function RootDecompositionPatternsByDistance(lat::HyperbolicLattice)
    
    W = lat.W
    
    # This is the fake dist to be minimized
    # Only define it for verification purposes
    # We can't just use `fake_dist(::Pattern,::HyperbolicLattice)` because we use it on `a₀`, not necessarily an integer.
    f_dist(k,w,a) = (a*v₀_norm(lat) + (times_v₀(lat,w)))^2//k

    next_least_a_for_w_and_k = SortedDict{Tuple{HyperbolicLatticeElement,Int},Int}()
     
    for w_ in W, k in root_lengths(lat)
        w = lat(w_) 
        # Find the minimal ``a₀``
        a₀::Rational{Int} = -times_v₀(lat,w)//v₀_norm(lat); 
        @toggled_assert -(w⊙v₀(lat))//(v₀(lat)⊙v₀(lat)) == -times_v₀(lat,w)//v₀_norm(lat) "The optimized computation should equal the full one."
        
        # Its integer approximation in the positive direction 
        a₊::Int = ceil(a₀)

        # verify that they are actually not as good
        @toggled_assert f_dist(k,w,a₊) ≥ f_dist(k,w,a₀)
        @debug "For (w=$w,k=$k), (a₀,a₊) is ($a₀,$a₊) with values $((f_dist(k,w,a₀),f_dist(k,w,a₊)))."
        
        # push the result
        push!(next_least_a_for_w_and_k, (w,k) => a₊)
    end

    current_fake_dist::Union{Nothing,Rational{Int}} = nothing

    return RootDecompositionPatternsByDistance(lat,next_least_a_for_w_and_k,current_fake_dist)


end

"""
    next!(pats::RootDecompositionPattern)

Returns:

* A `Vector` of `Pattern`s, all at the same distance.
* The distance.

The function proceeds as follows:
* Iterate over all `(w,k)`, and for each get the "iterator" `a₊`.
* Find the minimum distance over all those.
* Collect the `Pattern`s `(w,k,a)` (with `a` equal to `a₊`) attaining this minimum
* Return all of those, and the distance.
"""
function next!(pats::RootDecompositionPatternsByDistance)::Tuple{Vector{RootDecompositionPattern},Rational{Int}}
    
    lat = pats.lat
    
    # The fake distance, redefined here so that we don't have to give `lat` explicitely at each call.
    f_dist(k,w,a) = (a*v₀_norm(lat) + (times_v₀(lat,w)))^2//k
        
    # No min_val yet, no patterns either
    min_val = nothing
    min_patterns=Vector{Pattern}()
    for  ((w,k),a₊) in pats.next_least_a_for_w_and_k
        
        # Verify that everything we see now is farther away than the previous fake_dist.
        @toggled_assert pats.current_fake_dist === nothing || f_dist(k,w,a₊) ≥ pats.current_fake_dist 
      
        # We first treat `a₊`.

        # If `min_val` not set yet, or set but we are lower than it, push our pattern to `min_patterns`.
        # Note that `min_val === nothing` iff `min_patterns == []`.
        if min_val === nothing || f_dist(k,w,a₊) < min_val
            min_patterns = Vector{Pattern}([(k,w,a₊)])
            min_val = f_dist(k,w,a₊)

        # If we're exactly at `min_val`, push our `Pattern`
        elseif f_dist(k,w,a₊) == min_val
            push!(min_patterns,(k,w,a₊))
        end
    end
   
    # Store the new minimal distance.
    pats.current_fake_dist = min_val

    # For each `a` appearing in a pattern, "discard it" by advancing the corresponding iterator `a₊` `.
    for (k,w,a) in min_patterns
        a₊ = pats.next_least_a_for_w_and_k[(w,k)]
        if a₊ == a
            a₊ += 1
        end
        pats.next_least_a_for_w_and_k[(w,k)] = a₊
    end
    # Return the patterns, and their distance
    return (min_patterns,min_val)

end


"""
    roots_decomposed_into(lat,a,k)

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
function roots_decomposed_into(lat::HyperbolicLattice, a::HyperbolicLatticeElement, k::Int) end

@resumable function roots_decomposed_into(
    lat::HyperbolicLattice{r},
    a₀::Int,
    w::HyperbolicLatticeElement{r},
    k::Int,
    constraints=[],
)::HyperbolicLatticeElement{r} where {r}

    a = a₀*v₀(lat) + w
    
    w_coordinate = lat.W_coordinates[findall(x->x==vec(w),lat.W)[1]]

    M₁ = lat.P[:,2:end]
    
    println("All the constraints $constraints")

    # Translate our problem into a positive definite quadratic equation, solvable by `qsolve()`
    A = M₁' * lat.G * M₁
    b = 2 *( M₁' * lat.G * vec(a))
    γ = a⊙a - k
  
    sols = Vector{HyperbolicLatticeElement{r}}()
    # Finds solutions, translate them back to the lattice
    for u in qsolve_diag_con(diag(A),b,γ,constraints)
        uu = lat(M₁ * u + vec(a))
        @toggled_assert norm(uu) == k  "``u⊙u`` must be equal to k"
        @toggled_assert (uu-a)⊙v₀(lat) == 0 "``u-a`` must lie in `V₁`"
        if is_root(uu,k)
            @yield uu 
        end
    end
end




@resumable function roots_by_distance(
    lat::HyperbolicLattice{r},
    interval_lower_bound=(x->true), # must match a non-empty interval of positive numbers
    interval_upper_bound=(x->true),
    known_roots=[],
) :: HyperbolicLatticeElement{r} where {r} 
    
    patterns = RootDecompositionPatternsByDistance(lat)


    while true
        
        (current_patterns,dist) = next!(patterns)

        if !interval_lower_bound(dist)
            continue
        end
        if !interval_upper_bound(dist)
            break
        end

        for (k,w,a) in current_patterns
            constraints = [(root[2:end] .* D[2:end],-root⊙(a*w + v₀(lat))) for root in known_roots]
            for root in roots_decomposed_into(lat,a,w,k,constraints)
                @toggled_assert fake_dist(lat,root) == patterns.current_fake_dist
                @yield root
            end
        end
    end
end
