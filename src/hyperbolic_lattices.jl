using AbstractAlgebra
using LinearAlgebra
using JSON
using StaticArrays
using ToggleableAsserts

import Base: vec, convert

include("util.jl")

ScaledInt = Int

Rat = Rational{Int}

"""
A hyperbolic lattice, defined by a symmetric matrix with integer coefficients, of signature `(n,1)`.

Stored in the structure are:

* `G`, the symmetric matrix defining the scalar product.
  It has to be symmetric of signature `(n,1)`.

For convenience/performance, also stored are:

* The rank `r` of the lattice.
  This is just `r == size(G)[1] == size(G)[2]`.
* The two components `P,D` of a diagonalization of `G`, i.e. a matrix `P` and a vector `D` satisfying `P'*G*P == diagm(D)`.
* The last invariant factor of the lattice, which can be computed as `abs(Int(det(rG)//gcd(cofactors_G)))`.
* The possible root lengths of the lattice, i.e. a finite list of integers that a root can have as a norm (roots and length to be defined below)
"""
struct HyperbolicLattice{rank}
    # Rank
    r::Int
    
    # Matrix defining the bilinear form: We store it as a `SMatrix` for performance reasons
    G::SMatrix{rank,rank,Int}
    
    # Diagonalisation pair.
    P::SMatrix{rank,rank,Int}
    Pinv::SMatrix{rank,rank,Rat}
    D::SVector{rank,Int}
    W::Vector{SVector{rank,Int}}
    common_denominator::Int
    W_coordinates::Vector{SVector{rank,ScaledInt}}

    # latast invariant factor and possible root lengths.
    last_invariant_factor::Int
    root_lengths::Vector{Int}
end


"""
    Asserts that G is a symmetric matrix of signature `(n,1)`.
"""
function assert_sig_n_1_matrix(G)
   
    rank = size(G)[1]
    @assert size(G) == (rank,rank)          "G must be a square matrix."
    @assert issymmetric(G)                  "G must be symmetric."
    @assert signature(G) == (rank-1,1)      "G must have signature (n,1)." 
end

"""
    Constructor for hyperbolic lattices.
    Only takes a matrix as its argument, and computes the other values.
"""
function HyperbolicLattice(G)

    assert_sig_n_1_matrix(G)
    
    rank = size(G)[1]

    # G can be essentially in any "matrix-like" form, so we convert it to a SMatrix (the small 's' stands for "static"
    sG = SMatrix{rank,rank,Int}(G)

    # Diagonalize G
    D,P = diagonalize(sG)
    @assert P'*G*P == D ""
    @assert isdiag(D)
    @assert D[1] < 0
    # Get the diagonal vector of `D`
    D = diag(D) 
   
  
    absdetP = Int(abs(det(Rational{Int}.(P))))

    W_,W_coordinates_ = get_integer_points_with_coordinates(P,absdetP)
    W = Vector{SVector{rank,Int}}(W_) # coset representatives under "natural" coordinates

    common_denominator = lcm([lcm(denominator.(w)) for w in W_coordinates_])
    W_coordinates = Vector{SVector{rank,Int}}([Int.(common_denominator*w) for w in W_coordinates_]) # coset represnetative under "diagonal" coordinates

    rG = Rational{Int}.(sG)
    cofactorsG = det(rG) * inv(rG) # 
    last_invariant_factor = abs(Int(det(rG)//gcd(cofactorsG)))
    
    # The possible lengths of roots are the divisor of `2*last_invariant_factor`.
    twice_LIF = 2*last_invariant_factor
    root_lengths =  [k for k in 1:twice_LIF if twice_LIF%k == 0]
    
    return HyperbolicLattice(
        rank,
        sG,
        P,
        inv(Rat.(P)),
        D,
        W,
        common_denominator,
        W_coordinates,
        last_invariant_factor,
        root_lengths
    )
end


struct FixedDenom{denominator}
    numerator::Int
end

function Base.:(+)(a::FixedDenom{d},b::FixedDenom{d}) where {d}
    return FixedDenom{d}(a.numerator+b.numerator)
end

struct HyperbolicLatticeElement{rk}
    lat::HyperbolicLattice{rk}
    diag_coordinates::SVector{rk,ScaledInt} 
end

function Base.:(==)(lat1::HyperbolicLattice,lat2::HyperbolicLattice)
    lat1.G == lat2.G
end

function root_lengths(lat::HyperbolicLattice)
    return lat.root_lengths
end

"""
The rank of the lattice.
"""
function rk(lat::HyperbolicLattice)
    return lat.r
end


function v₀(lat::HyperbolicLattice{rk}) where {rk}
    return HyperbolicLatticeElement(lat,SVector{rk,ScaledInt}(vcat([lat.common_denominator],[0 for i in 1:rk-1])))
end

function times_v₀(lat::HyperbolicLattice,v::HyperbolicLatticeElement)::Int
    return lat.D[1]*v.diag_coordinates[1]/lat.common_denominator
end

"""
    Create a hyperbolic lattice element out of a lattice `L` and a vector `vec` by first converting `vec` to a `SVector`.
    Allows using the syntax `L(vec)`.
"""
function (lat::HyperbolicLattice)(vec)
    svec = SVector{rk(lat),Int}(vec)
    HyperbolicLatticeElement(lat,Int.(lat.common_denominator*(lat.Pinv*svec)))
end

function Base.show(io::IO,v::HyperbolicLatticeElement) 
    show(io,vec(v))
end

function Base.isequal(v1::HyperbolicLatticeElement,v2::HyperbolicLatticeElement)
    v1.lat == v2.lat && v1.diag_coordinates == v2.diag_coordinates
end

function Base.:(==)(v1::HyperbolicLatticeElement,v2::HyperbolicLatticeElement)
    v1.lat == v2.lat && v1.diag_coordinates == v2.diag_coordinates
end

function diag_rep(v::HyperbolicLatticeElement)
    return (v.diag_coordinates,v.lat.common_denominator)
end

function diag_rep_rational(v::HyperbolicLatticeElement)
    return (1//lat.common_denominator) * diag_coordinates
end

function std_rep(v::HyperbolicLatticeElement{r})::SVector{r,Int} where {r}
    return (v.lat.Pinv * v.diag_coordinates)//v.lat.common_denominator
end

"""
    inner_product(v,w)

Compute the inner_product of its two arguments in their common lattice.
Obtained as `v' * G * w`.
"""
function inner_product(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement)::Int
    @toggled_assert v.lat == w.lat "The elements must belong to the same lattice"
    D = v.lat.D
    vv = v.diag_coordinates
    ww = w.diag_coordinates
    return Int(sum(vv .* D .*  ww)/v.lat.common_denominator^2)
end

⊙(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement) = inner_product(v,w)

function norm(v::HyperbolicLatticeElement)
    return inner_product(v,v)
end

"""
    The standard basis of the lattice, that is, given by the vectors (1,0,…,0),…,(0,…,0,1,0…,0),…,(0,…,0,1).
"""
function standard_basis(lat::HyperbolicLattice)
    return lat.(SVector{rk(lat),Int}.(collect(eachcol(I(rk(lat))))))
end

function Base.:-(v::HyperbolicLatticeElement) 
    return HyperbolicLatticeElement(v.lat,-v.diag_coordinates)
end 

function Base.:+(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement)
    @toggled_assert v.lat == w.lat "The elements must belong to the same lattice"
    return HyperbolicLatticeElement(v.lat,v.diag_coordinates+w.diag_coordinates)
end 

function Base.:-(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement)
    @toggled_assert v.lat == w.lat "The elements must belong to the same lattice"
    return v + (-w) 
end 

function Base.:*(k::Int, v::HyperbolicLatticeElement)
    return HyperbolicLatticeElement(v.lat,k*v.diag_coordinates)
end 

function v₀_norm(lat::HyperbolicLattice)
    lat.D[1]
end

"""
    Coxeter_coeff(r₁,r₂)

Compute the label corresponding to the edge between the vertices corresponding to the hyperplanes defined by `r₁` and `r₂` respectively.

# Remarks

* This is incomplete, and only works for acute angles I think.
* The code is copied from B&P, and I haven't yet made in more general.
* TODO: probably deserves using `r₁_pp` and `r₂_pp` to skip matrix multiplications where possible.
"""
function Coxeter_coeff(r₁::HyperbolicLatticeElement,r₂::HyperbolicLatticeElement)
    @toggled_assert r₁.lat == r₂.lat "The elements must belong to the same lattice"
    @toggled_assert is_root(r₁) && is_root(r₂) "The elements must be roots"

    angle = r₁⊙r₂
    cos² = angle^2//(norm(r₁)*norm(r₂))
    if cos² == 0
        return 2
    elseif cos² == 1//4
        return 3
    elseif cos² == 1//2
        return 4
    elseif cos² == 3//4
        return 6
    elseif cos² == 1
        return 0
    elseif cos² > 1
        return 1
    else
        println("ERROR")
        println("cos² is $cos²")
        println("r₁ is $r₁ and r₂ is $r₂")
        return nothing
    end

end

"""
Ordering on `HyperboliclatticeElement` simply defined as the ordering on their representing vectors.
"""
function Base.isless(v::HyperbolicLatticeElement, w::HyperbolicLatticeElement)
    @toggled_assert v.lat == w.lat "The elements must belong to the same lattice"
    return Base.isless(vec(v),vec(w))
end

"""
    crystallographic_condition(v[, norm_v])

Test whether ``v`` satisfies the crystallographic condition in its lattice ``L``, that is, whether ``v⊙w`` is in ``L`` for any ``w`` in ``L``.
If `norm_v` is not `nothing`, use it as the norm of `v`, i.e. assume that `norm(v) == norm_v`.
"""
function crystallographic_condition(v::HyperbolicLatticeElement,norm_v::Int)
    
    @toggled_assert norm_v == norm(v)

    @toggled_assert iff(
                        all(2*(e⊙v) % (v⊙v) == 0 for e in standard_basis(v.lat)), 
                        all(2*(row⋅vec(v)) % norm_v == 0 for row in eachrow(v.lat.G))
                       ) "Optimized check for the crystallographic condition should be equivalent to the complete one."
    return all( (((2*x) % norm_v) == 0)::Bool for x in v.lat.G * vec(v) )

end

function crystallographic_condition(v::HyperbolicLatticeElement)
    crystallographic_condition(v,norm(v))
end
    
function is_primitive(v::HyperbolicLatticeElement)
    return abs(gcd(vec(v))) == 1
end

"""
    is_root(v[, norm_v])

Test whether  `v` is a root of its lattice.
If `norm_v` is not nothing, assume it is equal to the norm of `v`, i.e. `norm_v == norm(v)`

We follow B&P in our definition of a root; `v` is a root if:

* `v` has positive norm.
* `v` is primitive, in that the gcd of its coefficients must be 1.
* `v` satisfies the crystallographic condition, meaning that reflecting along `v` preserves the lattice.
  This can be checked by verifying that the reflection of elements of the standard basis are still in the lattice.
"""
function is_root(v::HyperbolicLatticeElement,norm_v::Int)::Bool
    
    @toggled_assert norm(v) == norm_v
    
    # A root has  positive norm.
    # TODO: check that it's indeed "positive" and not "non-negative"
    norm_v ≤ 0 && return false

    # A root is a primitive vector, i.e. the gcd of its entries is 1
    !is_primitive(v) && return false

    # A root has length dividing twice the last invariant factor (see B&P) but 
    # this condition is actually not by definition, so we could skip it and get the same result afaik
    norm_v ∉ v.lat.root_lengths && return false

    # A root respects the crystallographic condition
    !crystallographic_condition(v,norm_v) && return false

    return true

end

function is_root(v::HyperbolicLatticeElement)::Bool
    return is_root(v,norm(v))
end

function vec(v::HyperbolicLatticeElement)
    return Int.(v.lat.P*v.diag_coordinates/v.lat.common_denominator) 
end

"""
    fake_dist(VL,e) == (e⊙v₀)²//e⊙e

where `Vlat::VinbergLattice` and `e::HyperbolicLatticeElement` is assumed to be a root.
The actual distance between the vector `v₀` and the hyperplane `e^⟂` in the hyperbolic space `ℍ^n` is monotonous with this "fake distance".
We rather use the "fake distance" since it does not involve approximate computations.
"""
function fake_dist(lat::HyperbolicLattice,e::HyperbolicLatticeElement)

    @toggled_assert is_root(e)
    @toggled_assert e.lat == lat
    @toggled_assert times_v₀(lat,e)^2//(e⊙e) == (e⊙v₀(lat))^2//(e⊙e)
    
    return times_v₀(lat,e)^2//(e⊙e)

end

function zero_elem(lat::HyperbolicLattice)
    HyperbolicLatticeElement(lat,SVector{rk(lat),Int}(zeros(rk(lat))))
end

"""
    Computes the reflection of `v` along `r`.
"""
function reflection(r::HyperbolicLatticeElement,v::HyperbolicLatticeElement)
   
    @toggled_assert v.lat == r.lat
    @toggled_assert is_root(r)

    return v - (2*(r⊙v)/(r⊙r))*r

end

"""
    Computes `sinh(the distance between v and the hyperplane defined by e)`.
    `e` must be a root.
"""
function sinh_distance_to_hyperplane(v::HyperbolicLatticeElement, e::HyperbolicLatticeElement)
    
    @toggled_assert v.lat == e.lat
    @toggled_assert is_root(e)
    
    return sqrt( - (e⊙v)^2 / ( norm(e)*norm(v) ) ) 
end 

distance_to_hyperplane(
    v::HyperbolicLatticeElement,
    e::HyperbolicLatticeElement
) =  asinh( sinh_distance_to_hyperplane(v,e) )


