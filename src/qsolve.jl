## Copied from (B&P's) qsolve.py

using LinearAlgebra
using AbstractAlgebra
using Base
using Memoize
using ToggleableAsserts
using LRUCache
using ResumableFunctions

include("util.jl")

"""
    qsolve_minimum(A,b,γ)
    
Return the minimum (point) of the form ``x'Ax + b'x + γ`` and its value, where ``A`` is positive definite.
"""
function qform_minimum(
    A::SMatrix{rk,rk,Int},
    b::SVector{rk,Int},
    γ::Int
   )::Tuple{SVector{rk,Float64},Float64} where {rk}

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
                push!(res,Int((-b-δ)//(2*a)))
            end
            if δ≠0 && (-b+δ) % (2*a) == 0
                push!(res,Int((-b+δ)//(2*a)))
            end
        end
    end
    return res
end

const max_cached_rk=2
const cache_sizes = [100000000,1000000000]
const qsolve_cache = [LRU{Tuple{SMatrix{i,i,Int},SVector{i,Int},Int,SVector{i,Float64}},Vector{SVector{i,Int}}}(maxsize=cache_sizes[i]) for i in 1:max_cached_rk]

"""
    qsolve_iterative(A::SMatrix{rk,rk,Int},b::SVector{rk, Int},γ::Int)

Solve the positive definite problem ``x'Ax + b'x + γ = 0`` with integral solutions, using an iterative approach.
This is the case where rk=1.
This method is pretty much entirely copied from B&P.
"""
function qsolve_iterative(
    A::SMatrix{1,1,Int},
    b::SVector{1, Int},
    γ::Int,
    x::SVector{1, Float64},
) where {rk}
    
    a = A[1,1]
    b = b[1]
    return [SVector{1,Int}(s) for s in int_solve_quadratic_poly(a,b,γ)]
   
end

"""
    qsolve_iterative(A::SMatrix{rk,rk,Int},b::SVector{rk, Int},γ::Int)

Solve the positive definite problem ``x'Ax + b'x + γ = 0`` with integral solutions, using an iterative approach.
This method is pretty much entirely copied from B&P.
"""
function qsolve_iterative(
    A::SMatrix{rk, rk,Int},
    b::SVector{rk, Int},
    γ::Int,
    x::SVector{rk, Float64},
) where {rk}


    @assert rk > 1 "Rank 1 case treated above"
    sols::Vector{SVector{rk,Int}} = []
    
    y = Int(floor(x[rk]))
    Aspecy,bspecy,γspecy = sub_A_b_γ(A,b,γ,y)
    yb = feasible(Aspecy,bspecy,γspecy)

    while yb ≠ nothing
        sols_y = cached_qsolve_iterative(Aspecy,bspecy,γspecy,yb)
        append!(sols,[vcat(sol,SVector{1,Int}(y)) for sol in sols_y])
        y -= 1
        Aspecy,bspecy,γspecy = sub_A_b_γ(A,b,γ,y)
        yb = feasible(Aspecy,bspecy,γspecy)
    end
    
    y = Int(floor(x[rk])+1)
    Aspecy,bspecy,γspecy = sub_A_b_γ(A,b,γ,y)
    yb = feasible(Aspecy,bspecy,γspecy)

    while yb ≠ nothing
        sols_y = cached_qsolve_iterative(Aspecy,bspecy,γspecy,yb)
        append!(sols,[vcat(sol,SVector{1,Int}(y)) for sol in sols_y])
        y += 1
        Aspecy,bspecy,γspecy = sub_A_b_γ(A,b,γ,y)
        yb = feasible(Aspecy,bspecy,γspecy)
    end   
    
    @toggled_assert length(Set(sols)) == length(sols) "No solution should be appearing twice"
    
    return sols 
end


function cached_qsolve_iterative(A::SMatrix{rk,rk,Int},b::SVector{rk, Int},γ::Int,x) where {rk}
    if rk > max_cached_rk
        return qsolve_iterative(A,b,γ,x)
    else
        return get!(qsolve_cache[rk],(A,b,γ,x)) do
            qsolve_iterative(A,b,γ,x)
        end
    end
end

"""
    qsolve(A,b,γ)

Use `qsolve_iterative` to find integral solutions to ``x'Ax + b'x + γ``, and return those.
The construction of `qsolve_iterative` is such that it can either return `nothing` or an empty list if there are no integer solutions: this function is just there to map theh case `nothing` to that of an empty list.
"""
function qsolve(
    A::SMatrix{rk,rk,Int},
    b::SVector{rk,Int},
    γ::Int
) where {rk}
    

    y = feasible(A,b,γ)
    if y ≠ nothing
        return cached_qsolve_iterative(A,b,γ,y)
    else
        return []
    end

end


function feasible(    
    A::SMatrix{rk,rk,Int},
    b::SVector{rk, Int},
    γ::Int,
)::Union{Nothing,SVector{rk,Float64}} where {rk}

    (min_point,min_val) = qform_minimum(A,b,γ)
    
    if min_val > 0.1 
        return nothing
    else
        return min_point::SVector{rk,Float64}
    end
end

@resumable function res_qsolve(
    A::SMatrix{rk,rk,Int},
    b::SVector{rk, Int},
    γ::Int,
)::SVector{rk,Int}  where {rk}
    x = feasible(A,b,γ)
    if ! isnothing(x) 
        for r in res_qsolve_iterative(A,b,γ,x)
            @yield r
        end
    end
end


function sub_A_b_γ(
    A::SMatrix{rk,rk,Int},
    b::SVector{rk, Int},
    γ::Int,
    z::Int,
) where {rk}
    return (
        SMatrix{rk-1,rk-1,Int}(A[1:rk-1,1:rk-1]),
        SVector{rk-1,Int}(b[1:rk-1] + z*A[rk,1:rk-1] + z*A[1:rk-1,rk]),
        (A[rk,rk]*z^2 +b[rk]*z + γ)::Int,
    )
end

@resumable function res_qsolve_iterative(
    A::SMatrix{rk,rk,Int},
    b::SVector{rk, Int},
    γ::Int,
    x::SVector{rk, Float64},
):: SVector{rk,Int}  where {rk}

    lower_qsolve = rk-1 ≤ max_cached_rk ? cached_qsolve_iterative : res_qsolve_iterative

    if rk == 1
        res = int_solve_quadratic_poly(A[1,1],b[1],γ)
        if res === nothing
            @yield nothing
        end
        for r in res
            @yield SVector{1,Int}(r)
        end
        return 
    end

    @assert rk > 1 "Rank 1 case treated above"
    y = Int(floor(x[rk]))
    
    Aspecy,bspecy,γspecy = sub_A_b_γ(A,b,γ,y)
    yb = feasible(Aspecy,bspecy,γspecy)

    while yb ≠ nothing
        
        for sol_y in lower_qsolve(Aspecy,bspecy,γspecy,yb)
            @yield vcat(sol_y,SVector{1,Int}(y))
        end
        y -= 1
        Aspecy,bspecy,γspecy = sub_A_b_γ(A,b,γ,y)
        yb = feasible(Aspecy,bspecy,γspecy)
    end
    
    y = Int(floor(x[rk])+1)

    Aspecy,bspecy,γspecy = sub_A_b_γ(A,b,γ,y)
    yb = feasible(Aspecy,bspecy,γspecy)
    
    while yb ≠ nothing

        for sol_y in lower_qsolve(Aspecy,bspecy,γspecy,yb)
            @yield vcat(sol_y,SVector{1,Int}(y))
        end
        y += 1
        Aspecy,bspecy,γspecy = sub_A_b_γ(A,b,γ,y)
        yb = feasible(Aspecy,bspecy,γspecy)
    end

end
