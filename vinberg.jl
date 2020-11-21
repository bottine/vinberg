#using Nemo
#using Hecke
#using AbstractAlgebra
using LinearAlgebra
using SymPy

# Code adapted from N. V. Bogachev and A. Yu. Perepechko:
#
#   https://github.com/aperep/vinberg-algorithm
#
#


function diagonalize(A)
    # returns T and D with G = TDT'
    # algorithm copied from there https://math.stackexchange.com/questions/1388421/reference-for-linear-algebra-books-that-teach-reverse-hermite-method-for-symmetr
    #
    
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
                #M[i,:] = (-M[i,i0]*M[i0,:] + M[i0,i0]*M[i,:])//g
                #M[:,i] = (-M[i0,i]*M[:,i0] + M[i0,i0]*M[:,i])//g
            end

            i0 = i0 + 1

        elseif M[i0,i0] == 0 && any(M[k,k]≠0 for k in range(i0+1,stop=n)) 
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
    
    @assert LinearAlgebra.isdiag(D) "D is diagonal", D
    @assert P'*A*P == D "We have a diagonalization"

    
    return (D,P)

end

function test_diagonalize()
    
    for i in range(1,stop=100)
        println("-")
        for n in range(4,stop=10)
            print(n,", ")

            M = rand(range(-20,stop=20),n,n)
            M = M + M' # make it symmetric
            @assert LinearAlgebra.issymmetric(M)
            (D,P) = diagonalize(M)
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

function inner_product(a,b,G)
    
    assert_sig_n_1_matrix(G)
    return a⋅(G*b)

end

function negative_vector(G)

    assert_sig_n_1_matrix(G)  
    
    # np = n+1
    np = size(G)[1]
    

    # ⊙ = \odot
    ⊙(a,b) = inner_product(a,b,G)
   

    # need to find matrices:
    # * D diagonal
    # * T with T⁻¹DT = G
    v0 = 0
    M = 1 
    # M originated from the perceived need to multiply v0 (as returned by the following loop) by a suitably large
    # number for its integer approximation to be precise but maybe that's not needed…

    D,T = diagonalize(G)     
    

    for i in eachindex(D)
        if D[i] < 0
            v0 = T[i,:]
            break
        end
    end
    v0 = (1/gcd(v0))*v0
    
    @assert inner_product(v0,v0,G) < 0 "v₀ must have negative norm"

    return v0

end

function get_integer_points(M,n)
    
    @assert (n,n) == size(M) "Matrix should be square"
    @assert (n,0) == signature(M)

    negative = [sum([min(v,0) for v in M[i,:]]) for i in range(1,stop=n)]
    positive = [sum([max(v,0) for v in M[i,:]]) for i in range(1,stop=n)]

    bounding_box = []
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

function reflection(r,G,v)

    assert_sig_n_1_matrix(G)
    
    return v - 2*(inner_product(r,v,G)/inner_product(r,v,G))*r

end
    

function Vinberg_Algorithm(G)
    ### G is a (n+1)×(n+1) square matrix
    assert_sig_n_1_matrix(G)
    
    v0 = negative_vector(G)
    println(v0)


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
