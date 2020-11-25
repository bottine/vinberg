## Copied from qsolve.py

using LinearAlgebra
using Base

#include("vinberg.jl")

function qform_minimum(A,b,γ)
    minushalf_b = -0.5 * b
    x = A \ minushalf_b 
    # then Ax = minushalf_b 

    # returns a minimum vector x for the quadratic polynomial, and its value
    #
    return (x, x ⋅ (A*x) + b ⋅ x + γ)
end


function solve_quadratic_poly(a,b,c)

    Δ = b^2 - 4*a*c
    if Δ ≥ 0
        δ = sqrt(Δ)
        return [(-b-δ)/2,(-b+δ)/2]
    else
        return []
    end

end


function bounding_box_diago(A,b,γ)
    

    n = size(A,1)

    # A is diagonal

    value(x) = x' * Diagonal(A) * x + b' * x + γ

    x,minval = qform_minimum(Diagonal(A),b,γ)
    
    if minval ≥ 0
        return Set()
    end

    max = []
    min = []
    for i in 1:n 
        sols = solve_quadratic_poly(A[i],b[i],γ-minval) 
        @assert size(sols,1) == 2 "Should always have solutions because pos def I think"
        push!(min,sols[1]-1) 
        push!(max,sols[2]+1) 
    end

    bounding_box = [[]]
    for i in 1:n 
        bounding_box = [vcat(vec, val) for vec in bounding_box for val in min[i]:max[i]] 
    end

    return bounding_box

end


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
            
            for i in i0+1:n
                g =  gcd(M[i0,i0],M[i0,i])
                mizi = M[i0,i]//g
                miziz = M[i0,i0]//g
                M[i,:] = (-mizi*M[i0,:] + miziz*M[i,:])
                M[:,i] = (-mizi*M[:,i0] + miziz*M[:,i])
            end
            i0 = i0 + 1

        elseif any(M[k,k]≠0 for k in i0+1:n) 

            k = [k for k in i0+1:n if M[k,k]≠0][1]
            M[i0,:], M[k,:] = M[k,:], M[i0,:]
            M[:,i0], M[:,k] = M[:,k], M[:,i0]
        
        elseif any(M[i,j] ≠ 0 for i in i0:n, j in i0:n)
        
            (i,j) = [(i,j) for  i in i0:n, j in i0:n if M[i,j]≠0][1]
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

function get_integer_points(M)
    
    n = size(M,1)
    @assert (n,n) == size(M) "Matrix should be square"

    negative = [sum([min(v,0) for v in M[i,:]]) for i in 1:n]
    positive = [sum([max(v,0) for v in M[i,:]]) for i in 1:n]

    bounding_box = [ [] ]
    neg_1 = negative[1]
    pos_1 = positive[1]
    bounding_box = [vcat(vec, val) for vec in bounding_box for val in neg_1-1:pos_1+1] 
    for i in 2:n
        neg_i = negative[i]
        pos_i = positive[i]
        bounding_box = [vcat(vec, val) for vec in bounding_box for val in neg_i-1:pos_i+1] 
    end

    function parallelipiped_contains(v)
        Q = v/M
        return all(c < 1 && c >= 0 for c in Q)
    end
    
    return [v' for v in bounding_box if parallelipiped_contains(v')]

end

function qsolve_naive(A,b,γ)

    D,P = diagonalize(A)
    d = diag(D)

    @assert P'*A*P == D "We have a diagonalization"
    @assert inv(P')' == inv(P)

    Q = inv(P)
    @assert Q'*D*Q == A
    
    value(x) = x'*A*x + b'*x + γ

    # x'*A*x + b'*x + γ  ~~~> x'*(Q'*D*Q)*x + (b'*inv(Q)) Q*x + γ
    
    bounding_box_diag = bounding_box_diago(d,P'*b,γ)
    representatives    = get_integer_points(P)


    solutions = Set()
    for v0 in bounding_box_diag, r in representatives
        v = v0'/P + r
        if value(v') == 0
            push!(solutions,v)
        end
    end
    return solutions
end

#function qsolve_iterative(A,b,γ)
#    # A a matrix symmetric pos def
#    # B a vector
#    # γ a scalar
#    # All with coeffs in Z
#    #
#    # Looks for a solution to ⟨x,Ax⟩ + ⟨x,b⟩ + γ = 0
#    #                      ie  x^tAx +  xb^t + γ = 0
#
#    n = size(A)[1]
#    sols = Set()
#    
#	if n == 1
#	
#		a = A[1,1]
#		b = b[1]
#		c = γ
#	
#        Δ = b^2 - 4*a*c
#        if Δ ≥ 0
#            δ = sqrt(Δ)
#            x = floor((-b+δ)/2)
#            if a*x^2 + b*x + c == 0
#                push!(sols,x)
#            end
#            x = floor((-b-δ)/2)
#            if a*x^2 + b*x + c == 0
#                push!(sols,x)
#            end
#    	end
#		
#
#    else # n > 1
#
#        (x, val) = qform_minimum(A,b,γ)
#
#
#        if val > 0
#        	# If the minimum is > 0, the quadratic form has no zero
#			# do nothing
#        else
#
#        	a = floor(x[n])
#        	
#			
#        	v = a
#        	while true
#        		# Looking for a solution x for ⟨x,Ax⟩ + ⟨x,b⟩ + γ = 0 with x_n=v fixed yields 
#        		# a (n-1)-dimensional quadratic form:        	
#        		AA = A[1:n-1,1:n-1]
#        		bb = b[1:n-1]+ v*A[n,1:n-1] + v*A[1:n-1,n]
#        		γγ = A[n,n]*v*v + b[n] + γ
#        		
#        		sols_v = qsolve_iterative(AA,bb,γγ)
#        		if isempty(sols_v) 
#        			break
#        		else 
#        			for s in sols_v
#        				push!(sols, push!(s,v))
#        			end
#        		end
#        		
#        		v = v+1
#        	end
#        	# same idea, opposite direction
#        	v = a
#        	while true
#        		# Looking for a solution x for ⟨x,Ax⟩ + ⟨x,b⟩ + γ = 0 with x_n=v fixed yields 
#        		# a (n-1)-dimensional quadratic form:        	
#        		AA = A[1:n-1,:n-1]
#        		bb = b[1:n-1]+ v*A[n,1:n-1] + v*A[1:n-1,n]
#        		γγ = A[n,n]*v*v + b[n] + γ
#        		
#        		sols_v = qsolve_iterative(AA,bb,γγ)
#        		if isempty(sols_v) 
#        			break
#        		else 
#        			for s in sols_v
#        				push!(sols, push!(s,v))
#        			end
#        		end
#        		
#        		v = v-1
#        	end
#        	
#        
#        end
#
#    end 
#    return sols
#    
#end


A = [1 2 3; 
     2 1 4; 
     3 4 1]
b = [0; 0; 0]
γ = -2


println(eigen(A))
println(qsolve_naive(A,b,γ))

