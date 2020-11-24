#using Nemo
#using Hecke
#using AbstractAlgebra
using LinearAlgebra
using SymPy

# Code adapted from N. V. Bogachev and A. Yu. Perepechko:
#
#   https://github.com/aperep/vinberg-algorithm

struct QuadLattice
    nplus::Integer
    nminus::Integer
    G::Symmetric{Integer,AbstractMatrix{Integer}}
    P::AbstractMatrix{Integer}
    D::Array{Integer,1} 
end

function rank(L::QuadLattice)
    return nplus + nminus
end

struct QuadLatticeElement
    L::QuadLattice
    vec::Array{Integer,1}
end

function inner_product(v::QuadLatticeElement,w::QuadLatticeElement)
    @assert v.L == w.L "The elements must belong to the same lattice"
    G = v.L.G
    vv = v.vec
    ww = w.vec
    return QuadLatticeElement(vv⋅(G*ww),v.L)
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

function -(v::QuadLatticeElement)
    @assert v.L == w.L "The elements must belong to the same lattice"
    return QuadLatticeElement(-v.vec,v.L)
end 

function +(v::QuadLatticeElement,w:QuadLatticeElement)
    @assert v.L == w.L "The elements must belong to the same lattice"
    return QuadLatticeElement(v.vec+w.vec,v.L)
end 

function *(k:Integer, v::QuadLatticeElement)
    return QuadLatticeElement(k*v.vec,v.L)
end 

function is_root(v::QuadLatticeElement)

    return all(2*(e⊙v) % v⊙v == 0 for e in standard_basis(v.L))

end

function reflection(r::QuadLatticeElement,v::QuadLatticeElement)
   
    @assert is_root(r) "r needs to be a root"

    return v - (2*(r⊙v)/(r⊙r))*r

end


#function diagonalize(A::LinearAlgebra.Symmetric{Core.Integer,Array{Core.Integer,2}})
function diagonalize(A)
    # returns T and D with G = TDT'
    # algorithm copied from there https://math.stackexchange.com/questions/1388421/reference-for-linear-algebra-books-that-teach-reverse-hermite-method-for-symmetr
    # plus a gcd step to reduce the growth of values 
    
    A = BigInt.(A)

    n = size(A)[1]
    i0 = 1
    M = [A I]
    while i0 ≤ n

        if M[i0,i0] ≠ 0
            
            for i in range(i0+1,stop=n)
                g =  gcd(M[i0,i0],M[i0,i])
                mizi = M[i0,i]//g
                miziz = M[i0,i0]//g
                M[i,:] = (-mizi*M[i0,:] + miziz*M[i,:])
                M[:,i] = (-mizi*M[:,i0] + miziz*M[:,i])
            end
            i0 = i0 + 1

        elseif any(M[k,k]≠0 for k in range(i0+1,stop=n)) 

            k = [k for k in range(i0+1,stop=n) if M[k,k]≠0][1]
            M[i0,:], M[k,:] = M[k,:], M[i0,:]
            M[:,i0], M[:,k] = M[:,k], M[:,i0]
        
        elseif any(M[i,j] ≠ 0 for i in range(i0,stop=n), j in range(i0,stop=n))
        
            (i,j) = [(i,j) for  i in range(i0,stop=n), j in range(i0,stop=n) if M[i,j]≠0][1]
            M[i,:] = M[i,:] + M[j,:]
            M[:,i] = M[:,i] + M[:,j]
        
        end

    end
   

    D = M[1:n,1:n]
    Q = M[1:n,n+1:2*n]
    P = Q'
    
    println(det(P))

    @assert LinearAlgebra.isdiag(D) "D is diagonal", D
    @assert P'*A*P == D "We have a diagonalization"

    
    return (D,P)

end

function test_diagonalize()
    
    for i in range(1,stop=10)
        println("-")
        for n in range(4,stop=10)
            print(n,", ")

            #M = rand(range(-20,stop=20),n,n)
            M = rand(range(-10,stop=10),n,n)
            M = M + M' # make it symmetric
            @assert LinearAlgebra.issymmetric(M)
            (D,P) = diagonalize(M)
            println(maximum(D),maximum(P))
            @assert P'*M*P == D

        end
    end

end 

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
    G,D,P = G,LinearAlgebra.diagonalize(L.D),L.P 
    @assert P'*G*P == D

    # D Necessarily has one negative entry since of signature (n,1)
    # Find the element corresponding to this negative entry
    v0 = [T[i,:] for i in eachindex(D) if D[i]<0][1]
    # normalize it 
    v0 = (1/gcd(v0))*v0
    
    v0 = QuadLatticeElement(L,v0)

    # verify that it indeed has negative norm
    @assert norm(v0) < 0 "v₀ must have negative norm"

    return v0

end

function first_roots(P,D)
    # Let's follow p.112 of Guglielmetti
    
#    seen = Set()
#    roots = Set()
#    for i in  eachindex(D) if D[i] ∉ seen
#        push!(seen,D[i])
#        
#        push!(roots, P[i,:])
#        last_j = i
#        for j in eachindex(D) if j > i && D[j] == D[i]
#                push!(roots, P[j,:] - P[last_j,:])
#                last_j = j
#        end end
#
#    end end
#    return roots
#
end

function get_integer_points(M,n)
    
    @assert (n,n) == size(M) "Matrix should be square"
    @assert (n,0) == signature(M)

    negative = [sum([min(v,0) for v in M[i,:]]) for i in range(1,stop=n)]
    positive = [sum([max(v,0) for v in M[i,:]]) for i in range(1,stop=n)]

    bounding_box = [ [] ]
    for i in range(1,stop=n)
        neg_i = negative[i]
        pos_i = positive[i]
        bounding_box = [hcat(vec, val) for vec in bounding_box for val in range(neg_i-1,stop=pos_i+1)] 
    end

    function parallelipiped_contains(v)
        Q = v/M
        return all(c < 1 && c >= 0 for c in Q)
    end
    
    return [v for v in bounding_box if parallelipiped_contains(v)]

end





function Vinberg_Algorithm(G)

    # Given a ``G`` is a ``(n+1)×(n+1)`` symmetric square matrix of signature ``(n,1)``, apply Vinberg's algorithm to (try to) find a finite volume fundamental polyhedron for the group ``O_r(L)`` of automorphisms (w.r.t the bilinear form) of ``ℤ^{n+1}`` generated by reflections along vectors of ``ℤ^{n+1}``.

    # We first verify that ``G`` really is a symmetric integral matrix of signature ``(n,1)``.
    assert_sig_n_1_matrix(G)
  
    # The bilinear form defined by ``G`` defines an inner product ``(a,b)↦a⋅(Gb))``.
    ⊙(a,b) = inner_product(a,b,G) 
    
    # We first start with just one point of the fundamental domain.
    v0 = negative_vector(G)
    
    M1 = G # should be gram matrix of v₀^⟂
    M2 = G # should be the matrix associated to V1, but I don't know what it is
    W = sort(get_integer_points([M2,v0]),by=w -> -(w⊙v0))
    
    En = abs(det(G)/gcd(adjoint(M)))
    root_lengths = [k for k in range(1,2*En+1) if (2*En)%k == 0]

    roots = first_roots(P,d)



end


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



test_diagonalize()
Vinberg_Algorithm(G1)
Vinberg_Algorithm(G2)
#Vinberg_Algorithm(G)
