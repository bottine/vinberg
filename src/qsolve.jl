## Copied from qsolve.py

# using Plots
using LinearAlgebra
using AbstractAlgebra
using Base
using Polyhedra

include("util.jl")

function qform_minimum(A,b,γ)

    #@info "> qform_minimum(…)"
    minushalf_b = -0.5 * b
    x = A \ minushalf_b 
    # then Ax = minushalf_b 

    # returns a minimum vector x for the quadratic polynomial, and its value
    #
    #println("qform_minimum($A,$b,$γ)…")

    #@info "< qform_minimum(…)"
    return (x, x ⋅ (A*x) + b ⋅ x + γ)
end


function solve_quadratic_poly(a,b,c)
    
    #println("lookking for sols of $a x² + $b x + $c == 0")

    Δ = b^2 - 4*a*c
    if Δ ≥ 0
        δ = sqrt(Δ)
        return [(-b-δ)/2,(-b+δ)/2]
    else
        return []
    end

end

# seems like their (B&P) version is way better than diagonalization + my naive check
function qsolve_iterative(A,b,γ)

end


function bounding_box_diago(A,b,γ)
   
    #@info "> bounding_box_diago(…)"

    @assert length(size(A)) == 1
    @assert isposdef(Diagonal(A))

    n = size(A,1)

    value(x) = x' * Diagonal(A) * x + b' * x + γ

    function qform_minimum_diago(A,b,γ)
        x = []
        for i in eachindex(A)
            push!(x, -b[i]/(2*A[i]))
        end
        return (x,value(x))
    end


    x,minval = qform_minimum_diago(A,b,γ)
    # could also have used x,minval = qform_minimum(Diagonal(A),b,γ)
   

    if minval ≥ 0
        return Set()
    end

    max = []
    min = []

    #sol_plus = []
    #sol_minus = []
    
    for i in 1:n 
        sols = solve_quadratic_poly(A[i],2*A[i]*x[i] + b[i],minval) 
        #println("($i) solutions are $sols")
        @assert size(sols,1) == 2 "Should always have solutions because pos def I think"
        @assert sols[1] < sols[2]
        push!(min,BigInt(round(x[i] + sols[1],digits=0))-1) 
        push!(max,BigInt(round(x[i] + sols[2],digits=0))+1) 
        #push!(sol_minus,sols[1])
        #push!(sol_plus,sols[2])
        #println("$(min[end]) -- $(sol_minus[end]) -- $(sol_plus[end]) -- $(max[end])")
    end
    
#    quadrant(vv) = [vv[i] > x[i] ? (vv[i]-x[i],sol_plus[i]) : (-(vv[i]-x[i]),-sol_minus[i]) for i in 1:n]
#    
#    in_diamond(vv) = begin
#        println("$vv gives $(quadrant(vv))")
#        println(sum([ww[1]/ww[2] for ww in quadrant(vv)]))
#        return sum([(ww[1]+1)/ww[2] for ww in quadrant(vv)]) < 1
#    end

    bounding_box = [[]]
    for i in 1:n 
        bounding_box = [vcat(vec, val) for vec in bounding_box for val in min[i]:max[i]] 
    end


    #@info "< bounding_box_diago(…)"
    return [vv for vv in bounding_box]

end

function qsolve(A::Array{Int,2},b::Array{Int,1},γ::Int)
    return qsolve_naive(BigInt.(A),BigInt.(b),BigInt.(γ))
end

function qsolve(A,b,γ)
    return qsolve_naive(A,b,γ)
end

function qsolve_naive(A::Array{BigInt,2},b::Array{BigInt,1},γ::BigInt)
    
#    println("QSOLVE x'Ax + b'x + γ = 0 with A as:")
#    display(A)
#    println("")
#    println(" and b as")
#    display(b)
#    println("")
#    println(" and γ as $γ")
    
    #@info "> qsolve_naive(…)"

    @assert issymmetric(A) "A must be symmetric"
    @assert isposdef(A) "A must be positive definite"
    # TODO put this back!!!!! but make sure it actually works

    n = size(A,1)
    
    D::Array{BigInt,2},P::Array{BigInt,2} = if isdiag(A)
        A,Diagonal(ones(n))
    else
        diagonalize(A)
    end
    d = diag(D)

#    println("Diagonalization: P is ($(typeof(P)))")
#    display(P)
#    println("\n and D is($(typeof(D)))")
#    display(D)
#    println("\n")


    @assert P'*A*P == D "We have a diagonalization"

    @assert inv(P')' == inv(P)
    Q = inv(Rational{BigInt}.(P)) # forces the inverse computation to be rational if I'm not mistaken
    @assert Q*P == I && P*Q == I
    @assert Q'*P' == I && P'*Q' == I 
    
    @assert Q'*D*Q == A
    
    value(x) = x'*A*x + b'*x + γ

    # x'*A*x + b'*x + γ  ~~~> x'*(Q'*D*Q)*x + (b'*inv(Q)) Q*x + γ
    
    bounding_box_diag = bounding_box_diago(d,P'*b,γ)
    representatives    = get_integer_points(P)
#
#    println("# representatives: ", length(representatives))
#    println("# bounding_box_diag: ", length(bounding_box_diag))
#
    solutions = Set()
    for v0 in bounding_box_diag, r in representatives
        v = P*v0 + r
        #println("testing $v = P⁻¹$v0 + $r")
        if value(v) == 0
            push!(solutions,v)
        end
    end
    

    #@info "< qsolve_naive(…)"
    return [BigInt.(s) for s in solutions] 
end


function draw_2d(A,b,γ;size=10)
    
    value(x) = x'*A*x + b'*x + γ
    
    #backend(:plotly)
    #VV = [value([i;j]) for i in -size:size, j in -size:size]
    #gui(heatmap(VV))

    dv(x) = if x == [0;0]
        if value(x) == 0
            "∅"
        else
            "X"
        end
    elseif value(x) == 0
        println("value($x) = $(value(x)) !!")
        "0"
    elseif value(x) < 0
        "⋅"
    else
        " "
    end

    VV = join([join([dv([i;j]) for i in -size:size])*"\n" for j in -size:size])
    println(VV)

end

function test_qsolve()
    
    
    # This specific example should only return one value:
    D = [1 0;
         0 2]
    P = [2 2;
         1 3]
    b = [-1; -2]
    c = -220
    A = P'*D*P
    @assert qsolve(A,b,c) == Set(Any[BigFloat[8.0, -6.0]])

    D = [1 0;
         0 1]
    P = [1 0;
         0 1]
    b = [0; 0]
    c = -25
    A = P'*D*P
    @assert qsolve(A,b,c) == Set(Any[[-5.0, 0.0], [3.0, 4.0], [4.0, 3.0], [0.0, 5.0], [-4.0, -3.0], [3.0, -4.0], [0.0, -5.0], [4.0, -3.0], [-3.0, -4.0], [-3.0, 4.0], [5.0, 0.0], [-4.0, 3.0]])



    
    D = [1 0 0;
         0 1 0;
         0 0 1]
    
    P = [1 0 0;
         0 1 0;
         0 0 1]
    
    A = P'*D*P
    
    b = [0; 0; 0]
    c = -26
    
    println(qsolve(A,b,c))


end

   D = [4 0 ;
        0 1]
    
    P = [1 0;
         0 1]
    
    A = P'*D*P
    
    b = [20; 20]
    c = -26

