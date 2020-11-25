## Copied from qsolve.py

using LinearAlgebra
using Base

include("util.jl")

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
    # since A is diago, this could be done by hand 
    
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


function qsolve_naive(A,b,γ)

    D,P = diagonalize(A)
    d = diag(D)

    @assert P'*A*P == D "We have a diagonalization"

    @assert inv(P')' == inv(P)
    Q = inv(Rational{Int}.(P)) # forces the inverse computation to be rational if I'm not mistaken
    @assert Q*P == I && P*Q == I
    @assert Q'*P' == I && P'*Q' == I 
    
    @assert Q'*D*Q == A
    
    value(x) = x'*A*x + b'*x + γ

    # x'*A*x + b'*x + γ  ~~~> x'*(Q'*D*Q)*x + (b'*inv(Q)) Q*x + γ
    
    bounding_box_diag = bounding_box_diago(d,P'*b,γ)
    representatives    = get_integer_points(P)

    println("# representatives: ", length(representatives))
    println("# bounding_box_diag: ", length(bounding_box_diag))

    solutions = Set()
    for v0 in bounding_box_diag, r in representatives
        v = inv(P)*v0 + r
        if value(v) == 0
            push!(solutions,v)
        end
    end
    return solutions
end



# println(qsolve_naive(A,b,γ))

