## Copied from qsolve.py

using Plots
using LinearAlgebra
using Base

include("util.jl")

function qform_minimum(A,b,γ)
    minushalf_b = -0.5 * b
    x = A \ minushalf_b 
    # then Ax = minushalf_b 

    # returns a minimum vector x for the quadratic polynomial, and its value
    #
    println("qform_minimum($A,$b,$γ)…")

    return (x, x ⋅ (A*x) + b ⋅ x + γ)
end


function solve_quadratic_poly(a,b,c)
    
    println("lookking for sols of $a x² + $b x + $c == 0")

    Δ = b^2 - 4*a*c
    if Δ ≥ 0
        δ = sqrt(Δ)
        return [(-b-δ)/2,(-b+δ)/2]
    else
        return []
    end

end


function bounding_box_diago(A,b,γ)
    
    @assert length(size(A)) == 1
    @assert isposdef(Diagonal(A))

    n = size(A,1)

    value(x) = x' * Diagonal(A) * x + b' * x + γ

    x,minval = qform_minimum(Diagonal(A),b,γ)
    println("returnning $x, $minval")
    # since A is diago, this could be done by hand 
    
    if minval ≥ 0
        return Set()
    end

    max = []
    min = []
    for i in 1:n 
        sols = solve_quadratic_poly(A[i],2*A[i]*x[i] + b[i],minval) 
        println("($i) solutions are $sols")
        @assert size(sols,1) == 2 "Should always have solutions because pos def I think"
        push!(min,round(x[i] + sols[1],digits=0)-1) 
        push!(max,round(x[i] + sols[2],digits=0)+1) 
    end

    bounding_box = [[]]
    for i in 1:n 
        bounding_box = [vcat(vec, val) for vec in bounding_box for val in min[i]:max[i]] 
    end
    
    println("Bounding box is:")
    display(bounding_box)
    println("\n ---------------")

    return bounding_box

end

qsolve(A,b,γ) = qsolve_naive(A,b,γ)

function qsolve_naive(A,b,γ)
    

    @assert issymmetric(A) "A must be symmetric"
    @assert isposdef(A) "A must be positive definite"

    n = size(A,1)
    
    D,P = if isdiag(A)
        A,Diagonal(ones(n))
    else
        diagonalize(A)
    end
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

    draw_2d(A,b,γ)

    return solutions
end


function draw_2d(A,b,γ,size=10)
    
    value(x) = x'*A*x + b'*x + γ
    
    #backend(:plotly)
    #VV = [value([i;j]) for i in -size:size, j in -size:size]
    #gui(heatmap(VV))

    dv(x) = if value(x) == 0 
        "0"
    elseif value(x) < 0
        "⋅"
    else
        "#"
    end

    VV = join([join([dv([i;j]) for i in -size:size])*"\n" for j in -size:size])
    println(VV)

end

function test_qsolve()
    D = [
    1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1]

P = [
    1 1 1 1;
    0 0 0 1;
    0 0 1 0;
    0 1 0 0]

    A = P'*D*P

    b = [1; 2; 3; 4]

    γ = 4



    #

end
