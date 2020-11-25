## Copied from qsolve.py

using LinearAlgebra
using Base








function diagonalize(A)
    # returns T and D with D = T'GT
    # algorithm copied from there https://math.stackexchange.com/questions/1388421/reference-for-linear-algebra-books-that-teach-reverse-hermite-method-for-symmetr
    # plus a gcd step to reduce the growth of values 
    
    A = BigInt.(A)

    n = size(A)[1]
    i0 = 1
    M = [A I]
    while i0 ≤ n
        
        
        
        min_i = sort(i0:n,by=(idx->abs(M[idx,idx])))[1]
        k = min_i 
        if min_i ≠ i0
            #println("Case 0")
            M[i0,:], M[k,:] = M[k,:], M[i0,:]
            M[:,i0], M[:,k] = M[:,k], M[:,i0]

        elseif M[i0,i0] ≠ 0
            #println("case I") 
            for i in i0+1:n
                g =  gcd(M[i0,i0],M[i0,i])
                mizi = M[i0,i]//g
                miziz = M[i0,i0]//g
                M[i,:] = (-mizi*M[i0,:] + miziz*M[i,:])
                M[:,i] = (-mizi*M[:,i0] + miziz*M[:,i])
            end
            i0 = i0 + 1

        elseif any(M[k,k]≠0 for k in i0+1:n) 
            #println("case II") 

            k = [k for k in i0+1:n if M[k,k]≠0][1]
            M[i0,:], M[k,:] = M[k,:], M[i0,:]
            M[:,i0], M[:,k] = M[:,k], M[:,i0]
            
            # Now we're good for case I
        
        elseif any(M[i,j] ≠ 0 for i in i0:n, j in i0:n)
            #println("case III") 
        
            (i,j) = sort([(i,j) for  i in i0:n, j in i0:n if M[i,j]≠0], by=((i,j)-> abs(M[i,j])))[1]
            M[i,:] = M[i,:] + M[j,:]
            M[:,i] = M[:,i] + M[:,j]

            # Now we're good for case II
        
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
    
    for i in 1:10
        println("-")
        for n in 4:10
            print(n,", ")

            #M = rand(-20:20,n,n)
            M = rand(-10:10,n,n)
            M = M + M' # make it symmetric
            @assert LinearAlgebra.issymmetric(M)
            (D,P) = diagonalize(M)
            println(maximum(D),maximum(P))
            @assert P'*M*P == D

        end
    end

end 

function get_integer_points(M)
   
    n = size(M,1)
    @assert (n,n) == size(M) "Matrix should be square (and invertible)"

    minimum = [sum([min(v,0) for v in M[i,:]]) for i in 1:n]
    maximum = [sum([max(v,0) for v in M[i,:]]) for i in 1:n]

    bounding_box = [[]] 
    for i in 1:n
        bounding_box = [vcat(vec, [val]) for vec in bounding_box for val in minimum[i]:maximum[i]-1] 
    end
    

    function parallelipiped_contains(v)
        Q = inv(M)*v
        return all(c < 1 && c >= 0 for c in Q)
    end
    
    integer_points = [Vector(v) for v in bounding_box if parallelipiped_contains(Vector(v))]
    @assert length(integer_points) ≤ abs(det(M)) "discrete volume ≤ volume"
    # TODO is the discrete volume equal always?

    return integer_points

end


