## Copied from (B&P's) qsolve.py

using LinearAlgebra
using AbstractAlgebra
using Base
using Memoize
using ToggleableAsserts

include("util.jl")

"""
    qsolve_minimum(A,b,γ)
    
Return the minimum (point) of the form ``x'Ax + b'x + γ`` and its value, where ``A`` is positive definite.
"""
function qform_minimum(
    A::SMatrix{rank,rank,Int},
    b::SVector{rank,Int},
    γ::Int
) where {rank}

    minushalf_b = -b//2
    
    x = A \ minushalf_b
    #=
    @toggled_assert floor.(inv(Rational{Int}.(A)) * minushalf_b)  == floor.(x)  """
    Approximate and rational computations should agree 
    ($(Int.(floor.(inv(Rational{Int}.(A)) * minushalf_b))) vs $(Int.(floor.(x))))
    ($(inv(Rational{Int}.(A)) * minushalf_b) vs $(x))
    """
    =# 
    # This check is problematic since on 0.9999999… the rational one gives the right floor (1) while the float one gives (0)
    # I think this is no problem but may need to ensure that, hence TODO

    
    return (x, x ⋅ (A*x) + b ⋅ x + γ)
end

"""
    solve_quadratic_poly(a,b,c)

Return integers solutions of the quadratic poly, if real solutions exist, or nothing otherwise
"""
@inline function int_solve_quadratic_poly(a::Int,b::Int,c::Int)
    
    Δ = b^2 - 4*a*c
    res = []
    if Δ ≥ 0
        δ = Int(round(sqrt(Δ),digits=0)) # actually faster than isqrt( ) it seems
        if δ^2 == Δ
            if (-b-δ) % (2*a) == 0
                push!(res,(-b-δ)//(2*a))
            end
            if δ≠0 && (-b+δ) % (2*a) == 0
                push!(res,(-b+δ)//(2*a))
            end
        end
    else
        return nothing
    end
    return res
end

"""
    qsolve_iterative(A::SMatrix{rank,rank,Int},b::SVector{rank, Int},γ::Int)

Solve the positive definite problem ``x'Ax + b'x + γ = 0`` with integral solutions, using an iterative approach.
This is the case where rank=1.
This method is pretty much entirely copied from B&P.
**Important.** By convention, if no solutions (even real ones) exist for the equation, returns nothing, and not just an empty list. 
"""
function qsolve_iterative(A::SMatrix{1,1,Int},b::SVector{1, Int},γ::Int) where {rank}
    
    a = A[1,1]
    b = b[1]
    return int_solve_quadratic_poly(a,b,γ)
   
end

"""
    qsolve_iterative(A::SMatrix{rank,rank,Int},b::SVector{rank, Int},γ::Int)

Solve the positive definite problem ``x'Ax + b'x + γ = 0`` with integral solutions, using an iterative approach.
This method is pretty much entirely copied from B&P.
**Important.** By convention, if no solutions (even real ones) exist for the equation, returns nothing, and not just an empty list. 
"""
@memoize Dict function qsolve_iterative(A::SMatrix{rank,rank,Int},b::SVector{rank, Int},γ::Int
) where {rank}


    @assert rank > 1 "Rank 1 case treated above"

    (min_point,min_val) = qform_minimum(A,b,γ)
    
    if min_val > 0.1 
        return nothing
    end
    x = min_point[end]
    sols::Vector{SVector{rank,Int}} = []

    x_floor = floor(x)
    x_ceil = x_floor+1

    A_(y)::SMatrix{rank-1,rank-1,Int} = A[1:end-1,1:end-1]
    b_(y)::SVector{rank-1,Int} = b[1:end-1] + y*A[end,1:end-1] + y*A[1:end-1,end]
    γ_(y)::Int = A[end,end]*y^2 +b[end]*y + γ

    y = x_floor
    sols_y = qsolve_iterative(A_(y),b_(y),γ_(y))
    while sols_y ≠ nothing
        append!(sols,[vcat(sol,SVector{1,Int}(y)) for sol in sols_y])
        y -= 1
        sols_y = qsolve_iterative(A_(y),b_(y),γ_(y))
    end
    
    y = x_ceil
    sols_y = qsolve_iterative(A_(y),b_(y),γ_(y))
    while sols_y ≠ nothing
        append!(sols,[vcat(sol,SVector{1,Int}(y)) for sol in sols_y])
        y += 1
        sols_y = qsolve_iterative(A_(y),b_(y),γ_(y))
    end
    
    @toggled_assert length(Set(sols)) == length(sols) "We should have not solution appearing twice"
    
    return sols 
end

"""
    qsolve(A,b,γ)

Use `qsolve_iterative` to find integral solutions to ``x'Ax + b'x + γ``, and return those.
The construction of `qsolve_iterative` is such that it can either return `nothing` or an empty list if there are no integer solutions: this function is just there to map theh case `nothing` to that of an empty list.
"""
function qsolve(
    A::SMatrix{rank,rank,Int},
    b::SVector{rank,Int},
    γ::Int
) where {rank}
    res = qsolve_iterative(A,b,γ)
    if res === nothing
        return []
    else 
        return res
    end
end

