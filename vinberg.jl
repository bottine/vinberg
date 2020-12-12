#using Nemo
#using Hecke
using AbstractAlgebra
using LinearAlgebra
# using SymPy

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
    @assert isdiag(D)
    D = diag(D) # Get the diagonal vector of `D` 

    return QuadLattice(n_pos,n_neg,G,P,D)
end

function rank(L::QuadLattice)
    return L.n_pos + L.n_neg
end

struct QuadLatticeElement
    L::QuadLattice
    vec::Array{BigInt,1}
    #vec::Array{Int,1}
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

function zero_elem(L)
    QuadLatticeElement(L,zeros(rank(L)))
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

# v₀ a vector of the lattice L
# Let V₁ the orthogonal subspace
# x ∈ V₁ iff x'*G*v₀  == 0, i.e. x⊙v₀ == 0
#
# WARNING: I'm not sure the right_kernel really gives us a basis…
#
function basis_of_orhogonal_complement(L::QuadLattice, v0::QuadLatticeElement)
    @assert v0.L == L

    S = MatrixSpace(ZZ, 1, rank(L))
    println("S is $S")
    Mv0 = L.G * v0.vec
    MMv0 = S(reduce(vcat,Mv0))
    right_ker_rank, right_ker = right_kernel(MMv0)


    span = [QuadLatticeElement(L,u) for u in eachcol(Matrix(right_ker))]
    return span

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

function roots_of_fundamental_cone(L::QuadLattice,v0::QuadLatticeElement, V1_basis::Array{QuadLatticeElement,1})
    
    @assert v0.L == L

    possible_roots = [v for k in root_lengths(L) for v in roots_decomposed_into(L,v0,V1_basis,zero_elem(L),k) if is_root(v)]
    println("possible roots are $possible_roots")
    filter!(is_root,possible_roots)
    println("after filtering $possible_roots")

    roots = []
    for r in possible_roots
        if is_necessary_hyperplane(roots,L.G, r)
            push!(roots,r)
        end
    end
    
    println("cone roots roots are now $roots")

    return roots

end



struct RootsByDistance
   v0::QuadLatticeElement
   norm::Integer
   current_batch::Array{QuadLattice,(1)}
end

function next!(r::RootsByDistance)
    return r.v0 # that's not a root 
end


function roots_decomposed_into(
    L::QuadLattice,
    v0::QuadLatticeElement,
    V1::Array{QuadLatticeElement,1},
    a::QuadLatticeElement,
    k::Integer)
    # cf eponimous function in B&P's code
    # with added assumption that a lies in span({v₀}), which simplifies the 
    # computations but only works for finding the cone roots

    # We are looking for a root ``v = a + v₁``
    # satisfying 
    # * ``v⊙v = k``
    # * ``v₁∈V₁``
    #
    # (v₁+a)⊙(v₁+a) = v₁⊙v₁ + 2 v₁⊙a + a⊙a
    # so
    # (v₁+a)⊙(v₁+a) = k iff  v₁⊙v₁ + 2 v₁⊙a = k-a⊙a,
    # and
    
    V1Mat = reduce(hcat,[v.vec for v in V1])

    solutions = qsolve(V1Mat' * L.G * V1Mat, V1Mat' * L.G * a.vec, a⊙a - k)

    solutions_in_L = (x -> QuadLatticeElement(L,V1Mat * x + a)).(solutions)
    
    return solutions_in_L

end


function is_finite_volume(roots::Array{QuadLatticeElement,(1)})
    false
end

function Vinberg_Algorithm(G)
    
    println("(n,1) form matrix: ")
    display(G)
    println()
    println("-------------------")

    L = QuadLattice(G)

    # We first start with just one point of the fundamental domain.
    v0 = negative_vector(L)

    println("negative norm vector v₀")
    display(v0.vec)
    println()
    println("-------------------")

    V1 = basis_of_orhogonal_complement(L,v0)


    roots::Array{QuadLatticeElement,(1)} = []

    append!(roots, roots_of_fundamental_cone(L,v0,V1))


#    while ! is_finite_volume(roots)
#        new_root = next!()
#       push!(roots,new_root)
#        println("hello")
#    end


end


G0 = [-3 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 2;
     ]

G1= [-10 0  0 0; 
      0  2  -1 0; 
      0 -1  2 0; 
      0  0   0 1]
G2 = [-7 0   0 0; 
      0 2 -1 0; 
      0 -1 2 -1; 
      0 0 -1 2]



#test_diagonalize()
Vinberg_Algorithm(G1)
Vinberg_Algorithm(G2)
#Vinberg_Algorithm(G)
