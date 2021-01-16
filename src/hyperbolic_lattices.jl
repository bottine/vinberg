using AbstractAlgebra
using LinearAlgebra
using JSON
using StaticArrays

import Base: vec, convert

include("util.jl")
include("qsolve.jl")

struct HyperbolicLattice{rank}
    r::Int
    G::SMatrix{rank,rank,Int}
    P::SMatrix{rank,rank,Int}
    D::SVector{rank,Int}
    last_invariant_factor::Int
    root_lengths::Vector{Int}
end

# Presumably any root must have length dividing twice the last invariant factor
function root_lengths(L::HyperbolicLattice)
    return L.root_lengths
end



function HyperbolicLattice(G)

    assert_sig_n_1_matrix(G)
    
    rank = size(G)[1]

    sG = SMatrix{rank,rank,Int}(G)

    D,P = diagonalize(sG)
    @assert P'*G*P == D
    @assert isdiag(D)
    D = diag(D) # Get the diagonal vector of `D`

    rG = Rational{Int}.(sG)
    cofactorsG = det(rG) * inv(rG) # https://stackoverflow.com/questions/58149645/julia-how-can-we-compute-the-adjoint-or-classical-adjoint-linear-algebra
    last_invariant_factor = abs(Int(det(rG)//gcd(cofactorsG)))

    twice_LIF = 2*last_invariant_factor
    root_lengths =  [k for k in 1:twice_LIF if twice_LIF%k == 0]
    
    return HyperbolicLattice(rank,sG,P,D,last_invariant_factor,root_lengths)
end


function Base.:(==)(L1::HyperbolicLattice,L2::HyperbolicLattice)
    L1.G == L2.G
end



function rank(L::HyperbolicLattice)
    return L.r
end

struct HyperbolicLatticeElement{rank}
    L::HyperbolicLattice{rank}
    vec::SVector{rank,Int}
    #vec::Array{Int,1}
end
function Base.isequal(v1::HyperbolicLatticeElement,v2::HyperbolicLatticeElement)
    v1.L == v2.L && v1.vec == v2.vec
end

function Base.:(==)(v1::HyperbolicLatticeElement,v2::HyperbolicLatticeElement)
    v1.L == v2.L && v1.vec == v2.vec
end
struct VinbergLattice{rank,rank_minus_1}
    # https://discourse.julialang.org/t/addition-to-parameter-of-parametric-type/20059/5
    
    # By which we mean, a quadratic lattice along with
    # * An element of negative norm v₀
    # * A basis V1_basis of the sublattice v₀^\perp = V₁
    # * A set of representatives W_reps for the quotient of L by V₁⊕ℤv₀

    L::HyperbolicLattice{rank}
    v0::HyperbolicLatticeElement{rank}
    V1Mat::SMatrix{rank,rank_minus_1,Int} # ideally a svector of quadlatticeelements since we know it has rank-1 elements
    W_reps::Vector{HyperbolicLatticeElement{rank}}
    v0vec_times_G::SVector{rank,Int} # used to hopefully make computations of the form v₀⊙something faster
    v0norm::Int # the norm of v₀ computed once and for all
end

function times_v0(VL::VinbergLattice,e::HyperbolicLatticeElement)
    #@assert e.L == VL.L

    return VL.v0vec_times_G ⋅ e.vec
end

function Base.:(==)(L1::VinbergLattice,L2::VinbergLattice)
    @assert false "Let's not compare vinberg lattices yet"
end


function VinbergLattice(G::Array{Int,2};v0vec::Union{Array{Int,1},Nothing}=nothing)
   
    assert_sig_n_1_matrix(G)
    
    rank = size(G)[1]

    L = HyperbolicLattice(G)
    v0 = (v0vec == nothing ? negative_vector(L) : HyperbolicLatticeElement(L,SVector{rank,Int}(v0vec)))
    
    @assert norm(v0) < 0

    V1_basis = basis_of_orhogonal_complement(L,v0)

    V1Mat = SMatrix{rank,rank-1,Int}(vcat(V1_basis...))

    M = hcat(v0.vec,V1_basis...)
    sM = SMatrix{rank,rank}(M)
    W = get_integer_points(M)
    W_reps = (x -> HyperbolicLatticeElement(L,SVector{rank,Int}(x))).(W)
   
    @assert length(W_reps) == abs(det(M)) "only $(length(W_reps)) but need $(det(M))"
    @assert norm(v0) % length(W_reps) == 0
   
    v0vec_times_G = SVector{rank,Int}(v0.vec' *  G)
    v0norm = v0⊙v0


    return VinbergLattice(L,v0,V1Mat,W_reps,v0vec_times_G,v0norm)   
end

function Base.convert(v::HyperbolicLatticeElement) 
    return v.vec::Array{Int,1} 
end

same_quad_lattice(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement) = v.L == w.L

function inner_product(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement)
    #@assert same_quad_lattice(v,w) "The elements must belong to the same lattice"
    G = v.L.G
    vv = v.vec
    ww = w.vec
    return vv⋅(G*ww)
end

⊙(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement) = inner_product(v,w)

function norm(v::HyperbolicLatticeElement)
    return inner_product(v,v)
end

function standard_basis(L::HyperbolicLattice)
    I = LinearAlgebra.Diagonal(ones((rank(L))))
    basis = []
    for i in 1:rank(L)
        push!(basis,HyperbolicLatticeElement(L,SVector{rank(L),Int}(I[:,i]))) 
    end
    return basis
end

function Base.:-(v::HyperbolicLatticeElement)         ## the `Base.:` syntax is necessary to make julia undenstand that we want to extend the function `-` to our type
    return HyperbolicLatticeElement(-v.vec,v.L)
end 

function Base.:+(v::HyperbolicLatticeElement,w::HyperbolicLatticeElement)
    #@assert same_quad_lattice(v,w) "The elements must belong to the same lattice"
    return HyperbolicLatticeElement(v.L,v.vec+w.vec)
end 

function Base.:*(k::Int, v::HyperbolicLatticeElement)
    return HyperbolicLatticeElement(v.L,k*v.vec)
end 

function is_root(v::HyperbolicLatticeElement)

    nv = norm(v)
    # A root has (non strictly) positive norm
    if nv < 0
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
    #if ! all(2*(e⊙v) % (v⊙v) == 0 for e in standard_basis(v.L))
    if ! all(2*(col⋅v.vec) % nv == 0 for col in eachcol(v.L.G))
        return false
    end

    return true

end

function fake_dist(VL::VinbergLattice,e::HyperbolicLatticeElement)
    # computes the value - (e⊙v₀)² / e⊙e ,
    # which is monotonous with the distance between v₀ and H_e (the hyperplane defined by e)

    @assert e.L == VL.L
    # @assert is_root(e)
    
    return -(VL.v0vec_times_G ⋅ e.vec)^2//(e⊙e)

end


function zero_elem(L::HyperbolicLattice)
    HyperbolicLatticeElement(L,SVector{rank(L),Int}(zeros(rank(L))))
end

function zero_elem(VL::VinbergLattice)
    zero_elem(VL.L) 
end

function reflection(r::HyperbolicLatticeElement,v::HyperbolicLatticeElement)
   
    #@assert same_quad_lattice(v,r) "The elements must belong to the same lattice"
    #@assert is_root(r) "r needs to be a root"

    return v - (2*(r⊙v)/(r⊙r))*r

end


function vec(v::HyperbolicLatticeElement)
    return v.vec
end

# v₀ a vector of the lattice L
# Let V₁ the orthogonal subspace
# x ∈ V₁ iff x'*G*v₀  == 0, i.e. x⊙v₀ == 0
#
# WARNING: I'm not sure the right_kernel really gives us a basis…
#
function basis_of_orhogonal_complement(L::HyperbolicLattice, v0::HyperbolicLatticeElement)
    @assert v0.L == L

    S = MatrixSpace(ZZ, 1, rank(L))
    Mv0 = L.G * v0.vec
    MMv0 = S(reduce(vcat,Mv0))
    right_ker_rank, right_ker = right_kernel(MMv0)
    @assert right_ker_rank == rank(L) - 1 "We need a basis!"


    span = [u for u in eachcol(Matrix(right_ker))]
    return span

end


# The distance between the hyperplane H_{e} and the point v0
function sinh_distance_to_hyperplane(v0::HyperbolicLatticeElement, e::HyperbolicLatticeElement)
    
    #@assert same_quad_lattice(v0,e) "The elements must belong to the same lattice"
    #@assert is_root(e)
    
    return sqrt( - (e⊙v0)^2 / ( norm(e)*norm(v0) ) ) 
end 
# The distance between the hyperplane H_{e} and the point v0
function distance_to_hyperplane(v0::HyperbolicLatticeElement, e::HyperbolicLatticeElement)
    
    return asinh( sinh_distance_to_hyperplane(v0,e) )

end 


function signature(G)
    
    D,P = diagonalize(G)
    pos = filter(>(0),diag(D))
    neg = filter(<(0),diag(D))
    
    return (length(pos),length(neg))
end


function assert_sig_n_1_matrix(G)
    
    @assert length(size(G)) == 2            "G must be a matrix"
    @assert size(G)[1] == size(G)[2]        "G must be square"
    @assert issymmetric(G)                  "G must be symmetric"
    np = size(G)[1]
    @assert true "TODO"  (np-1,1) == signature(G)        "G must have signature (n,1)" 
end



function negative_vector(L::HyperbolicLattice)

    assert_sig_n_1_matrix(L.G)  
    
    # np = n+1
    np = rank(L)
    # Find P and D diagonal such that P'GP = D is diagonal
    G,D,P = L.G,L.D,L.P 
    @assert P'*G*P == Diagonal(L.D)

    # D Necessarily has one negative entry since of signature (n,1)
    # Find the element corresponding to this negative entry
    v0 = [P'[i,:] for i in eachindex(D) if D[i]<0][1]
    # normalize it 
    v0 = (1//gcd(v0))*v0
    
    v0 = HyperbolicLatticeElement(L,SVector{rank(L),Int}(v0))
    

    # verify that it indeed has negative norm
    @assert norm(v0) < 0 "v₀ must have negative norm"

    return v0

end



