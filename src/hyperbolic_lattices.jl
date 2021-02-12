using AbstractAlgebra
using LinearAlgebra
using JSON
using StaticArrays
using ToggleableAsserts

import Base: vec, convert

include("util.jl")
include("qsolve.jl")




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
    D::SVector{rank,Int}
    
    # Last invariant factor and possible root lengths.
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
    # Get the diagonal vector of `D`
    D = diag(D) 
    
    # Compute the last invariant factor:
    # * taking `rG` should force the determinant and inverse computations to be exact.
    # * we use [this idea](https://stackoverflow.com/questions/58149645/julia-how-can-we-compute-the-adjoint-or-classical-adjoint-linear-algebra) to compute the cofactor matrix.i
    #
    # TODO: do we **need** the `abs(Int(…))` part?
    rG = Rational{Int}.(sG)
    cofactorsG = det(rG) * inv(rG) # 
    last_invariant_factor = abs(Int(det(rG)//gcd(cofactorsG)))
    
    # The possible lengths of roots are the divisor of `2*last_invariant_factor`.
    twice_LIF = 2*last_invariant_factor
    root_lengths =  [k for k in 1:twice_LIF if twice_LIF%k == 0]
    
    return HyperbolicLattice(rank,sG,P,D,last_invariant_factor,root_lengths)
end

function Base.:(==)(L1::HyperbolicLattice,L2::HyperbolicLattice)
    L1.G == L2.G
end

function root_lengths(L::HyperbolicLattice)
    return L.root_lengths
end

"""
The rank of the lattice.
"""
function rk(L::HyperbolicLattice)
    return L.r
end


"""
An element of a hyperbolic lattice, represented by:
* `L::HyperbolicLattice`, the lattice it is a member of.
* `vec`, the vector of its coordinate in the standard basis.
    
"""
struct HyperbolicLatticeElement{rank}
    L::HyperbolicLattice{rank}
    vec::SVector{rank,Int}
end


"""
    Create a hyperbolic lattice element out of a lattice `L` and a vector `vec` by first converting `vec` to a `SVector`.
    Allows using the syntax `L(vec)`.
"""
function (L::HyperbolicLattice)(vec)
    svec = SVector{rk(L),Int}(vec)
    HyperbolicLatticeElement(L,svec)
end

function Base.show(io::IO,v::HyperbolicLatticeElement) 
    show(io,v.vec)
end

function Base.isequal(v1::HyperbolicLatticeElement,v2::HyperbolicLatticeElement)
    v1.L == v2.L && v1.vec == v2.vec
end

function Base.:(==)(v1::HyperbolicLatticeElement,v2::HyperbolicLatticeElement)
    v1.L == v2.L && v1.vec == v2.vec
end

"""
    inner_product(v,w)

Compute the inner_product of its two arguments in their common lattice.
Obtained as `v' * G * w`.
"""
function inner_product(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement)
    @toggled_assert v.L == w.L "The elements must belong to the same lattice"
    G = v.L.G
    vv = v.vec
    ww = w.vec
    return vv⋅(G*ww)
end

⊙(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement) = inner_product(v,w)

function norm(v::HyperbolicLatticeElement)
    return inner_product(v,v)
end

"""
    The standard basis of the lattice, that is, given by the vectors (1,0,…,0),…,(0,…,0,1,0…,0),…,(0,…,0,1).
"""
function standard_basis(L::HyperbolicLattice)
    return L.(SVector{rk(L),Int}.(collect(eachcol(I(rk(L))))))
end

function Base.:-(v::HyperbolicLatticeElement) 
    return HyperbolicLatticeElement(v.L,-v.vec)
end 

function Base.:+(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement)
    @toggled_assert v.L == w.L "The elements must belong to the same lattice"
    return HyperbolicLatticeElement(v.L,v.vec+w.vec)
end 

function Base.:-(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement)
    @toggled_assert v.L == w.L "The elements must belong to the same lattice"
    return v + (-w) 
end 

function Base.:*(k::Int, v::HyperbolicLatticeElement)
    return HyperbolicLatticeElement(v.L,k*v.vec)
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
    @toggled_assert r₁.L == r₂.L "The elements must belong to the same lattice"
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
        return nothing
    end

end

"""
Ordering on `HyperboliclatticeElement` simply defined as the ordering on their representing vectors.
"""
function Base.isless(v::HyperbolicLatticeElement, w::HyperbolicLatticeElement)
    @toggled_assert v.L == w.L "The elements must belong to the same lattice"
    return Base.isless(v.vec,w.vec)
end

"""
    crystallographic_condition(v[, norm_v])

Test whether ``v`` satisfies the crystallographic condition in its lattice ``L``, that is, whether ``v⊙w`` is in ``L`` for any ``w`` in ``L``.
If `norm_v` is not `nothing`, use it as the norm of `v`, i.e. assume that `norm(v) == norm_v`.
"""
function crystallographic_condition(v::HyperbolicLatticeElement,norm_v=nothing::Union{Nothing,Int})
     
    nv::Int = norm_v === nothing ? norm(v) : norm_v
    @toggled_assert iff(
                        all(2*(e⊙v) % (v⊙v) == 0 for e in standard_basis(v.L)), 
                        all(2*(row⋅v.vec) % nv == 0 for row in eachrow(v.L.G))
                       ) "Optimized check for the crystallographic condition should be equivalent to the complete one."
    return all( (((2*x) % nv) == 0)::Bool for x in v.L.G * v.vec )

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
function is_root(v::HyperbolicLatticeElement,norm_v=nothing::Union{Nothing,Int})
    
    nv::Int = norm_v === nothing ? norm(v) : norm_v
    
    # A root has  positive norm.
    # TODO: check that it's indeed "positive" and not "non-negative"
    if nv ≤ 0
        return false
    end

    # A root is a primitive vector, i.e. the gcd of its entries is 1
    vv = v.vec
    if abs(gcd(vv)) ≠ 1
        return false
    end

    # A root has length dividing twice the last invariant factor (see B&P) but 
    # this condition is actually not by definition, so we could skip it and get the same result afaik
    if nv ∉ v.L.root_lengths
        return false
    end

    # A root respects the crystallographic condition
    if !crystallographic_condition(v,nv) 
        return false
    end

    return true

end


"""
A `VinbergLattice` is simply a `HyperbolicLattice` `L` plus a choice of vector `v₀` of the lattice, of negative norm.
Given `L` and `v₀`, some more objects can be computed:
Let `V₁` be the subgroup of `L` consisting of elements orthogonal to v₀, i.e. `V₁ = v₀^⟂`.
We can compute a basis of `V₁`, using methods of the package `AbstractAlgebra`.
The sublattice `⟨v₀⟩⊕V₁` is of finite index in L, and we can therefore find a finite set of (coset) representatives.

We therefore store:

* The lattice `L`,
* The negative norm vector `v₀`

And the related data for convenience/performance:

* `V₁_basis_matrix` is rank×rank_minus_1 matrix whose columns are a basis for `V₁`
* `W` a set of coset representatives for the sublattice `⟨v₀⟩⊕V₁`.
* `v₀_vec_times_G == v₀.vec' * G` corresponds to the computation v₀.vec ⋅ G, which makes multiplication of v₀ with other elements faster:
  `v₀⊙e = v₀.vec' * L.G * e.vec = v₀_vec_times_G * e.vec` so that we can skip one matrix multiplication
* `v₀_norm == v₀⊙v₀` is the norm of v₀, which we precompute sincet it is used often.

Note that the type VinbergLattice{rank,rank_minus_1} is parametrized by the two parameters `rank,rank_minus_1`, but they must satisfy `rank-1 == rank_minus_1`.
This useless second parameter is due to a [limitation of julia](https://discourse.julialang.org/t/addition-to-parameter-of-parametric-type/20059/5).
"""
struct VinbergLattice{rank,rank_minus_1}
    
    # The lattice
    L::HyperbolicLattice{rank}
   
    # The negative norm vector "basepoint"
    v₀::HyperbolicLatticeElement{rank}

    # A basis of V₁, in the form of a matrix
    V₁_basis_matrix::SMatrix{rank,rank_minus_1,Int}

    # The representatives for `⟨v₀⟩⊕V₁` 
    W::Vector{HyperbolicLatticeElement{rank}}

    # Precomputations
    v₀_vec_times_G::SVector{rank,Int} # used to hopefully make computations of the form v₀⊙something faster
    v₀_norm::Int # the norm of v₀ computed once and for all
end

"""
The result of multiplying `e` with `v₀`, i.e. `e⊙v₀`.
"""
function times_v₀(VL::VinbergLattice,e::HyperbolicLatticeElement)
    @toggled_assert e.L == VL.L

    @toggled_assert VL.v₀_vec_times_G ⋅ e.vec == VL.v₀⊙e "Optimized computation should not differ from official one"
    return VL.v₀_vec_times_G ⋅ e.vec
end

function Base.:(==)(L1::VinbergLattice,L2::VinbergLattice)
    @assert false "Let's not compare vinberg lattices yet"
end

"""
    VinbergLattice(G[;v₀_vec])

Construct a `VinbergLattice`, with an optional negative norm vector `v₀`.
If none is provided, one is found by diagonalizing the form.
"""
function VinbergLattice(G::Array{Int,2};v₀_vec::Union{Array{Int,1},Nothing}=nothing)
   
    assert_sig_n_1_matrix(G)
    
    rk = size(G)[1]

    L = HyperbolicLattice(G)
    v₀ = (v₀_vec == nothing ? negative_vector(L) : HyperbolicLatticeElement(L,SVector{rk,Int}(v₀_vec)))
    
    @toggled_assert norm(v₀) < 0 "v₀ needs to have negative norm"
    
    # Compute a matrix whose columns form a basis of `V₁`
    V₁_basis_matrix = SMatrix{rk,rk-1,Int}(matrix_basis_of_orthogonal_complement(v₀))

    # Computes a matrix `M` whose columns are the basis of `V₁` computed before, plus `v₀`
    M = hcat(v₀.vec,V₁_basis_matrix)
    sM = SMatrix{rk,rk}(M)
   
    # Compute the representatives of `⟨v₀⟩⊕V₁`
    W = get_integer_points(M)
    W = (x -> HyperbolicLatticeElement(L,SVector{rk,Int}(x))).(W)
  
    v₀_norm = norm(v₀)

    @toggled_assert length(W) == abs(det(Rational{Int}.(M))) "The number of representative, i.e. the index of the sublattice must be the same as the determinant of the matrix whose columns are a basis for the sublattice.\n Here we have $(length(W)) ≠ $(abs(det(Rational{Int}.(M))))"
    @toggled_assert v₀_norm % length(W) == 0 "According to B&P, the norm of v₀ must divide the index of the sublattice"
   
    v₀_vec_times_G = SVector{rk,Int}(v₀.vec' *  G)

    @info "Initialized Vinberg Lattice:"
    @info "v₀ is $v₀"
    @info "V₁_basis_matrix is $V₁_basis_matrix"
    @info "W is $W"

    return VinbergLattice(L,v₀,V₁_basis_matrix,W,v₀_vec_times_G,v₀_norm)   
end

function Base.convert(v::HyperbolicLatticeElement) 
    return v.vec::Array{Int,1} 
end

"""
    fake_dist(VL,e) == (e⊙v₀)²//e⊙e

where `VL::VinbergLattice` and `e::HyperbolicLatticeElement` is assumed to be a root.
The actual distance between the vector `v₀` and the hyperplane `e^⟂` in the hyperbolic space `ℍ^n` is monotonous with this "fake distance".
We rather use the "fake distance" since it does not involve approximate computations.
"""
function fake_dist(VL::VinbergLattice,e::HyperbolicLatticeElement)

    @toggled_assert is_root(e)
    @toggled_assert e.L == VL.L
    @toggled_assert (VL.v₀_vec_times_G ⋅ e.vec)^2//(e⊙e) == (e⊙VL.v₀)^2//(e⊙e)
    @toggled_assert times_v₀(VL,e)^2//(e⊙e) == (e⊙VL.v₀)^2//(e⊙e)
    
    return times_v₀(VL,e)^2//(e⊙e)

end

function zero_elem(L::HyperbolicLattice)
    HyperbolicLatticeElement(L,SVector{rk(L),Int}(zeros(rk(L))))
end

function zero_elem(VL::VinbergLattice)
    zero_elem(VL.L) 
end

"""
    Computes the reflection of `v` along `r`.
"""
function reflection(r::HyperbolicLatticeElement,v::HyperbolicLatticeElement)
   
    @toggled_assert v.L == r.L
    @toggled_assert is_root(r)

    return v - (2*(r⊙v)/(r⊙r))*r

end

"""
Computes a basis (over ℤ) for the orthogonal complement `v^⟂` of `v`.
An element `e` lies in `v^⟂` iff `e.vec'*G*v.vec == e⊙v == 0`, which is iff `e.vec` is in the kernel of `G*v.vec`

"""
function matrix_basis_of_orthogonal_complement(v::HyperbolicLatticeElement)
   
    rank = rk(v.L)

    # We make use of `AbstractAlgebra`, so first need to convert to its types.
    # We create a space of mutrices over ℤ and of dimensions `1×rk(L)`
    S = MatrixSpace(ZZ, 1, rank)

    # Politely ask AbstractAlgebra for a basis
    Mv = v.L.G * v.vec
    MMv = S(reduce(vcat,Mv))
    @toggled_assert MMv == S(Mv')
    right_ker_rank, right_ker = right_kernel(MMv)
    
    @toggled_assert right_ker_rank == rank - 1 "We need a basis for V₁; its dimension must be rk(L)-1"

    # convert back to a type we understand
    return SMatrix{rank,rank-1}(Matrix(right_ker))

end

"""
    Computes `sinh(the distance between v and the hyperplane defined by e)`.
    `e` must be a root.
"""
function sinh_distance_to_hyperplane(v::HyperbolicLatticeElement, e::HyperbolicLatticeElement)
    
    @toggled_assert v.L == e.L
    @toggled_assert is_root(e)
    
    return sqrt( - (e⊙v)^2 / ( norm(e)*norm(v) ) ) 
end 

distance_to_hyperplane(v::HyperbolicLatticeElement, e::HyperbolicLatticeElement) =  asinh( sinh_distance_to_hyperplane(v,e) )

"""
Find a vector of negative norm in the lattice, by diagonalization.
"""
function negative_vector(L::HyperbolicLattice)
    
    rank = rk(L)
    
    # Find `P` and `D` diagonal such that `P'GP = D` is diagonal
    G,D,P = L.G,L.D,L.P 
    @toggled_assert P'*G*P == Diagonal(D)

    # D necessarily has one negative entry since the signature of `G` is `(rk-1,1)`
    # Find the element corresponding to this negative entry
    v₀_vec = [P'[i,:] for i in eachindex(D) if D[i]<0][1]
    # Normalize it 
    v₀_vec = (1//gcd(v₀_vec))*v₀_vec
    
    v₀ = HyperbolicLatticeElement(L,SVector{rank,Int}(v₀_vec))
    

    # verify that it indeed has negative norm
    @toggled_assert norm(v₀) < 0 "v₀ must have negative norm"

    return v₀

end



