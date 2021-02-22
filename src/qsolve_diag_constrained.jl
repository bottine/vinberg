## Copied from (B&P's) qsolve.py

using LinearAlgebra
using AbstractAlgebra
using Base
using Memoize
using ToggleableAsserts
using LRUCache
using ResumableFunctions
using IterTools

include("util.jl")

Rat = Rational{Int}

val(D,b,γ,z) = dot(z,diagm(D),z) + dot(b,z) + γ

function min_quad(
    a::Int,
    b::Int,
    c::Int,
)::Rat
    
    @assert a>0
    min_point = -b//(2*a)
    return min_point
end

function min_quad(
    D::SVector{n,Int},
    b::SVector{n,Int},
    γ::Int,
)::SVector{n,Rat} where {n}

    @toggled_assert all(d > 0 for d in D) 
    min_point =  -(b .// (2*D))
    @toggled_assert min_point == - (inv(diagm(Rat.(D)))*b)//2
    return min_point
end

@inline function float_solve_quadratic_poly(a::Int,b::Int,c::Rat)
    
    Δ = b^2 - 4*a*c
    res = []
    if Δ ≥ 0
        δ = sqrt(Δ)
        push!(res,(-b-δ)/(2*a))
        if δ≠0 
            push!(res,(-b+δ)/(2*a))
        end
    end
    return res
end


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

function updateat(v::SVector{n,T},a::T,idx::Int) where {n,T}
    @toggled_assert 1 ≤ idx && idx ≤ n
    insert(deleteat(v,idx),idx,a)
end

function last_non_zero(v) 
    for i in length(v):-1:1
        if v[i] ≠ 0
            return i
        end
    end
    return -1
end

@resumable function qsolve_diag_con(
    D::SVector{n,Int},
    b::SVector{n,Int},
    γ::Int,
    con#::Vector{Tuple{SVector{n,Int},Int}},
)::SVector{n,Int} where {n}

    con = [filter(x -> last_non_zero(x[1]) == i, con) for i in 1:n]
    ccon = map(level -> Vector{AffineConstraint{n}}(map(c -> (SVector{length(c[1]),Int}(c[1]),c[2]), level)), con) 
    
    for s in qsolve_diag_constrained(SVector{n,Int}(D),SVector{n,Int}(b),γ,ccon)
        @yield s
    end

end

AffineConstraint{N} = Tuple{SVector{N,Int},Int}
# A pair vec, val representing the half-space {x: x⋅vec≤val}

ConstraintsByLastNonZero{N} = Vector{Vector{AffineConstraint{N}}}

satisfied(c::AffineConstraint) = c[2] ≥ 0

function sub_first_in_constraint(
    c::AffineConstraint{n},
    x::Int
)::AffineConstraint{n-1} where {n}
    vec,val = c
    return (popfirst(vec),val-vec[1]*x)
end

function sub_first_in_constraints(
    constraints::ConstraintsByLastNonZero{n},
    x::Int,
)::Tuple{Bool,ConstraintsByLastNonZero{n-1}} where {n}
    level0 = constraints[1]
    feasible = all([satisfied(sub_first_in_constraint(c,x)) for c in level0])
    !feasible && return (feasible,[])
    
    return (feasible, [map(c -> sub_first_in_constraint(c,x),level) for level in constraints[2:end]])

end

update_D_b_γ(D,b,γ,last_x) = (popfirst(D),popfirst(b),γ+D[1]*last_x^2 + last_x*b[1])
matching_D_b_γ(D,b,γ,min_but) = (D[1],b[1],γ+dot(min_but,diagm(D[2:end]),min_but) + dot(min_but,b[2:end]))

@resumable function qsolve_diag_constrained(
    D::SVector{n,Int},
    b::SVector{n,Int},
    γ::Int,
    constraints::ConstraintsByLastNonZero{n},
    depth=0,
)::SVector{n,Int} where {n}
    
    if any(length(level) > 0 for level in constraints)
        #println(" "^depth, "qsolve_diag_constrained:")   
        #println(" "^depth, "$D")   
        #println(" "^depth, "$b")   
        #println(" "^depth, "$γ")
        for level in constraints
            for c in level
                #println(" "^depth, c)  
            end
            #println(" "^depth, "---------------------------------------------------")
        end
    end
    @toggled_assert all(d > 0 for d in D) 
   

    if n == 1
    
        candidate_sols = SVector{1,Int}.(int_solve_quadratic_poly(D[1],b[1],γ))
        for s in candidate_sols
            feasible, remaining = sub_first_in_constraints(constraints,s[1])
            feasible && @yield s
        end

    else

        min_point = min_quad(D,b,γ) 

        min_val = val(D,b,γ,min_point)


        if min_val ≤ 0

            x_zeros = float_solve_quadratic_poly(matching_D_b_γ(D,b,γ,min_point[2:end])...)
            x_bottom = min(x_zeros...)
            x_top = max(x_zeros...)

            for x in Int(floor(x_bottom)):Int(ceil(x_top))
                (feasible,updated) = sub_first_in_constraints(constraints,x)
                if feasible
                    for sol_x in qsolve_diag_constrained(update_D_b_γ(D,b,γ,x)...,updated,depth+1)
                        
                        @yield pushfirst(sol_x,x)
                    end
                else
                    #println("---")
                end
            

            end
        end
    end

end
