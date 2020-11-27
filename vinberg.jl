#using Nemo
#using Hecke
#using AbstractAlgebra
using LinearAlgebra
using SymPy

# Code adapted from N. V. Bogachev and A. Yu. Perepechko:
#
#   https://github.com/aperep/vinberg-algorithm


# TODO
#
# * See here for memoization: https://github.com/JuliaCollections/Memoize.jl
# * Extensive unit tests!
#


include("util.jl")
include("qsolve.jl")

struct QuadLattice
    n_pos::Integer
    n_neg::Integer
    G::LinearAlgebra.Symmetric{Integer,AbstractMatrix{Integer}}
    P::AbstractMatrix{Integer}
    D::Array{Integer,1} 
end

function QuadLattice(G)

    assert_sig_n_1_matrix(G)

    G = Symmetric(G) 

    np = size(G)[1]
    n_pos = np-1
    n_neg = 1
    D,P = diagonalize(G)
    @assert P'*G*P == D
    @assert Diagonal(D) == D
    D = diag(D) # Get the diagonal vector of `D` 

    return QuadLattice(n_pos,n_neg,G,P,D)
end

function rank(L::QuadLattice)
    return L.n_pos + L.n_neg
end

struct QuadLatticeElement
    L::QuadLattice
    vec::Array{BigInt,1}
end

same_quad_lattice(v::QuadLatticeElement,w::QuadLatticeElement) = v.L == w.L

function inner_product(v::QuadLatticeElement,w::QuadLatticeElement)
    @assert same_quad_lattice(v,w) "The elements must belong to the same lattice"
    G = v.L.G
    vv = v.vec
    ww = w.vec
    return vv⋅(G*ww)
end

⊙(v::QuadLatticeElement,w::QuadLatticeElement) = inner_product(v,w)

function norm(v::QuadLatticeElement)
    return inner_product(v,v)
end

function standard_basis(L::QuadLattice)
    I = LinearAlgebra.diagonal(ones((rank(L))))
    basis = []
    for i in 1:rank(L)
        push!(basis,QuadLatticeElement(I[:,i],L)) 
    end
    return basis
end

function Base.:-(v::QuadLatticeElement)         ## the `Base.:` syntax is necessary to make julia undenstand that we want to extend the function `-` to our type
    @assert same_quad_lattice(v,w) "The elements must belong to the same lattice"
    return QuadLatticeElement(-v.vec,v.L)
end 

function Base.:+(v::QuadLatticeElement,w::QuadLatticeElement)
    @assert same_quad_lattice(v,w) "The elements must belong to the same lattice"
    return QuadLatticeElement(v.vec+w.vec,v.L)
end 

function Base.:*(k::Integer, v::QuadLatticeElement)
    return QuadLatticeElement(k*v.vec,v.L)
end 

function is_root(v::QuadLatticeElement)

    if norm(v) < 0
        return false
    end

    return all(2*(e⊙v) % v⊙v == 0 for e in standard_basis(v.L))

end

function reflection(r::QuadLatticeElement,v::QuadLatticeElement)
   
    @assert same_quad_lattice(v,r) "The elements must belong to the same lattice"
    @assert is_root(r) "r needs to be a root"

    return v - (2*(r⊙v)/(r⊙r))*r

end

# I **think** this is the last invariant factor (Prop 1 in Bogachev&Perepechko)
# TODO: check that? what does it mean?
function last_invariant_factor(L::QuadLattice)
    return abs(det(L.G)/gcd(adjoint(L.G)))
end

function vec(v::QuadLatticeElement)
    return v.vec
end

# The vectors [v0 v1 … v_n] = P are orthogonal, hence linearly independent, so that they span a finite index subgroup of L.
# This function should spit out a set of representatives 
function diagonalization_representatives(L::QuadLattice)
    
    P = L.P
    G = L.G
    
    # The columns of P define a basis for the subgroup, so any element in the subgroup must lie in the parallelipiped spanned by these vectors
    # get_integer_points returns exactly those elements

    return (x -> QuadLatticeElement(x,L)).(get_integer_points(P,rank(L)))

end 


# Presumably any root must have length dividing twice the last invariant factor
function root_lengths(L::QuadLattice)
    return [k for k in 1:last_invariant_factor(L) if last_invariant_factor(L)%k == 0]
end

# The distance between the hyperplane H_{e} and the point v0
function distance_to_hyperplane(v0::QuadLatticeElement, e::QuadLatticeElement)
    
    @assert same_quad_lattice(v0,e) "The elements must belong to the same lattice"
    @assert is_root(e)
    
    return asinh( sqrt( - (e⊙v0)^2 / ( norm(e)*norm(v0) ) ) )

end 

#function diagonalize(A::LinearAlgebra.Symmetric{Core.Integer,Array{Core.Integer,2}})


function signature(G)
    
    D,T = LinearAlgebra.eigen(G)
    
    pos = 0
    neg = 0
    for d in D
        if d < 0
           neg = neg + 1
        elseif d > 0
           pos = pos + 1
        end
    end
    
    return (pos,neg)

end


function assert_sig_n_1_matrix(G)
    
    @assert length(size(G)) == 2            "G must be a matrix"
    @assert size(G)[1] == size(G)[2]        "G must be square"
    @assert LinearAlgebra.issymmetric(G)    "G must be symmetric"
    np = size(G)[1]
    @assert (np-1,1) == signature(G)        "G must have signature (n,1)" 
    
end



function negative_vector(L::QuadLattice)

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
    
    v0 = QuadLatticeElement(L,v0)
    

    # verify that it indeed has negative norm
    @assert norm(v0) < 0 "v₀ must have negative norm"

    return v0

end

function roots_of_fundamental_cone(v0::QuadLatticeElement)

    return []

end



struct RootsByDistance
   v0::QuadLatticeElement
   norm::Integer
   current_batch::Array{QuadLattice,(1)}
end

function next!(r::RootsByDistance)
    return r.v0 # that's not a root 
end


function is_finite_volume(roots::Array{QuadLatticeElement,(1)})
    false
end

function Vinberg_Algorithm(G)


    L = QuadLattice(G)

    # We first start with just one point of the fundamental domain.
    v0 = negative_vector(L)
   
    roots::Array{QuadLatticeElement,(1)} = []

    append!(roots, roots_of_fundamental_cone(v0))


#    while ! is_finite_volume(roots)
#        new_root = next!()
#       push!(roots,new_root)
#        println("hello")
#    end


end


G0 = [-3 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 2;
     ]

G1= [-10 0  0 0; 
      0  2  -1 0; 
      0 -1  2 0; 
      0  0   0 1]
G2 = [-7 0   0 0; 
      0 2 -1 0; 
      0 -1 2 -1; 
      0 0 -1 2]
G = [1 2 3; 
     2 3 4; 
     3 4 5]



#test_diagonalize()
Vinberg_Algorithm(G1)
Vinberg_Algorithm(G2)
#Vinberg_Algorithm(G)
