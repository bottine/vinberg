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

function val( # zDz + zb + γ with all of those with indices ranging from n to N
    k::Int,
    D::SVector{N},
    b::SVector{N},
    γ::Int,
    z::SVector{N}
   ) where {N} 
    diag_product(k,z,D,z) + dot(k,b,z) + γ
end

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
    k::Int,
    D::SVector{N,Int},
    b::SVector{N,Int},
    γ::Int,
)::SVector{N-k+1,Rat} where {N}

    @toggled_assert all(d > 0 for d in D) 
    min_point =  -(b[k:N] .// D[k:N]) //2
    @toggled_assert min_point == - (inv(diagm(Rat.(D[k:N])))*b[k:N])//2
    return min_point
end

@inline function float_solve_quadratic_poly(a::Int,b::Int,c::Rat)
    
    Δ = b^2 - 4*a*c
    res = Float64[]
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
    res = Int[]
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

function unzip(pairs::Vector{Tuple{SVector{N,Int},Int}})::Tuple{Vector{SVector{N,Int}},Vector{Int}} where {N} 
    return ([p[1] for p in pairs],[p[2] for p in pairs])
end

@resumable function qsolve_diag_con(
    D::SVector{N,Int},
    b::SVector{N,Int},
    γ::Int,
    constraints::Vector{Tuple{SVector{N,Int},Int}},
    the_norm::Int,
    w_diag::SVector{N,Int},
    common_denom::Int,
)::SVector{N,Int} where {N}

    con = [filter(x -> last_non_zero(x[1]) == i, constraints) for i in 1:N]
    ccon = [unzip(level) for level in con] 


    for s in qsolve_diag_constrained(SVector{0,Int}(),D,b,γ,min_quad(1,D,b,γ),ccon,the_norm,w_diag,common_denom)
        @yield s
    end

end



function sub_kth_in_constraint(
    k::Int,
    x::Int,
    vec::SVector{N,Int},
    val::Int
)::Int where {N}
    return val-vec[k]*x
end

function sub_kth_in_constraints(
    k::Int,
    x::Int,
    vecs::Vector{SVector{N,Int}},
    vals::Vector{Int}
)::Vector{Int} where {N}
    return [sub_kth_in_constraint(k,x,vec,val) for (vec,val) in zip(vecs,vals)]
end

function sub_kth_in_constraints(
    k::Int,
    x::Int,
    constraints::Vector{Tuple{Vector{SVector{N,Int}},Vector{Int}}},
)::Tuple{Bool,Vector{Tuple{Vector{SVector{N,Int}},Vector{Int}}}} where {N}

    level0 = constraints[1]
    feasible = all([v ≥ 0 for v in sub_kth_in_constraints(k,x,level0[1],level0[2])])
    !feasible && return (feasible,[])
    
    return (feasible, [(level[1],sub_kth_in_constraints(k,x,level[1],level[2])) for level in constraints[2:end]])

end

update_γ(k,D,b,γ,last_x) = γ+D[k]*last_x^2 + last_x*b[k]
matching_D_b_γ(k,D,b,γ,min) = (D[k],b[k],γ+diag_product(k,min,D,min) + dot(k,min,b))

@resumable function qsolve_diag_constrained(
    prefix::SVector{k_minus_one,Int},
    D::SVector{N,Int},
    b::SVector{N,Int},
    γ::Int,
    min_point::SVector{N,Rat},
    constraints::Vector{Tuple{Vector{SVector{N,Int}},Vector{Int}}},
    the_norm::Int,
    w_diag::SVector{N,Int},
    common_denom::Int,
    depth=0,
)::SVector{N,Int} where {k_minus_one,N}
    
    k = k_minus_one + 1
    @toggled_assert all(d > 0 for d in D) 
  

    if k == N
    
        candidate_sols = SVector{1,Int}.(int_solve_quadratic_poly(D[k],b[k],γ))
        for s in candidate_sols
            feasible, remaining = sub_kth_in_constraints(k,s[1],constraints)
            feasible && @yield vcat(prefix,s)
        end

    else

        min_val = val(k,D,b,γ,min_point)

        if min_val ≤ 0

            x_zeros = float_solve_quadratic_poly(matching_D_b_γ(k,D,b,γ,min_point)...)
            x_bottom = min(x_zeros...)
            x_top = max(x_zeros...)

            for x in Int(floor(x_bottom)):Int(ceil(x_top))
                (feasible,updated_cons) = sub_kth_in_constraints(k,x,constraints)
               
                if feasible && Int(2*(x*D[k] + w_diag[k]*D[k]//common_denom)//common_denom) % the_norm == 0
                    for sol_x in qsolve_diag_constrained(push(prefix,x),D,b,update_γ(k,D,b,γ,x),min_point,updated_cons,the_norm,w_diag,common_denom,depth+1)
                        
                        @yield sol_x 
                    end
                else
                    # not interested, thanks 
                end
            

            end
        end
    end

end
