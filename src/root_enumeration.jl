
include("util.jl")
include("qsolve.jl")
include("hyperbolic_lattices.jl")



mutable struct RootsByDistance
    VL::VinbergLattice
    next_least_a0s_for_k_and_w::Union{Nothing,Vector{Tuple{Int,HyperbolicLatticeElement,Int,Int}}}
    current_a0_and_k_and_w::Union{Nothing,Tuple{Int,Int,HyperbolicLatticeElement}}
    roots_for_current_a0_and_k_and_w::Vector{HyperbolicLatticeElement}
    current_fake_dist::Union{Nothing,Rational{Int}}
end



function RootsByDistance(VL::VinbergLattice)
    
    v0 = VL.v0

    W_reps = VL.W_reps


    next_least_a0s_for_k_and_w = Vector{Tuple{Int,HyperbolicLatticeElement,Int,Int}}()
     
    for w in W_reps, k in root_lengths(VL.L)
        least::Rational{Int} = -times_v0(VL,w)//VL.v0norm # -(w⊙v0)/(v0⊙v0)
        #@info "least($(w.vec),$k) at $least"
        least_plus::Int = ceil(least)
        least_minus::Int = floor(least)
        push!(next_least_a0s_for_k_and_w, (k,w,least_plus,least_minus))
    end


    current_a0_and_k_and_w::Union{Nothing,Tuple{Int,Int,HyperbolicLatticeElement}} = nothing
    #current_a0_and_k_and_w = nothing;
    roots_for_current_a0_and_k_and_w::Vector{HyperbolicLatticeElement} = Vector{HyperbolicLatticeElement}()
    
    current_fake_dist::Union{Nothing,Rational{Int}} = nothing

    return RootsByDistance(VL,next_least_a0s_for_k_and_w,current_a0_and_k_and_w,roots_for_current_a0_and_k_and_w,current_fake_dist)

end


function next_with_decomposition!(r::RootsByDistance)
    #@info "> next!(roots_by_distance)"
    v0 = r.VL.v0
    VL = r.VL

    to_minimize(a0,k,w) = (a0*VL.v0norm + (times_v0(VL,w)))^2//k 
    

    while r.current_a0_and_k_and_w === nothing || isempty(r.roots_for_current_a0_and_k_and_w)
        
        min_val = nothing
        min_tuple = nothing
        min_idx = nothing
        for (idx, (k,w,a0plus,a0minus)) in enumerate(r.next_least_a0s_for_k_and_w)
            if min_tuple === nothing || to_minimize(a0plus,k,w) < min_val
                min_tuple = (a0plus,k,w)
                min_val = to_minimize(a0plus,k,w)
                min_idx = idx
            end
            if min_tuple === nothing || to_minimize(a0minus,k,w) < min_val
                min_tuple = (a0minus,k,w)
                min_val = to_minimize(a0minus,k,w)
                min_idx = idx
            end
        end
        
        # we got the triplet k,w,a0 minimizing to_minimizeance
        r.current_a0_and_k_and_w = min_tuple
         
        (a0,k,w) = r.current_a0_and_k_and_w
        
        # we update the dictionary
        (k,w,a0plus,a0minus) = r.next_least_a0s_for_k_and_w[min_idx]
        if a0plus == a0
            a0plus += 1
        end
        if a0minus == a0
            a0minus -= 1
        end
        
        r.next_least_a0s_for_k_and_w[min_idx] = (k,w,a0plus,a0minus)

        # we update fake_dist
        r.current_fake_dist = (a0*VL.v0norm + VL.v0vec_times_G ⋅ w.vec)^2 // k

        r.roots_for_current_a0_and_k_and_w = roots_decomposed_into(r.VL,a0*v0 + w,k)
    end
 
    #@info "< next!(roots_by_distance)"
    return (r.current_a0_and_k_and_w,pop!(r.roots_for_current_a0_and_k_and_w),r.current_fake_dist)

end


function next_with_dist!(r::RootsByDistance)
    
    res = next_with_decomposition!(r)
    return (res[2],res[3])

end

function next!(r::RootsByDistance)
    return next_with_decomposition!(r)[2]
end


function roots_decomposed_into(VL::VinbergLattice, a::HyperbolicLatticeElement, k::Int)
    # cf eponimous function in B&P's code

    # We are looking for a root ``v = a + v₁``
    # satisfying 
    # * ``v⊙v = k``
    # * ``v₁∈V₁``
    #
    # (v₁+a)⊙(v₁+a) = v₁⊙v₁ + 2 v₁⊙a + a⊙a
    # so
    # (v₁+a)⊙(v₁+a) = k iff  v₁⊙v₁ + 2 v₁⊙a = k-a⊙a,
    # and
   
    #@info "> roots_decomposed_into(VL, $(a.vec), $k)"
    #println("…  $((-(a⊙VL.v0))/(k^0.5))")
    
    r = rank(VL.L)
    V1Mat = VL.V1Mat

    #solutions = qsolve_naive(BigInt.(V1Mat' * VL.L.G * V1Mat), BigInt.(V1Mat' * VL.L.G * a.vec), BigInt(a⊙a - k))
    A::SMatrix{r-1,r-1,Int} = V1Mat' * VL.L.G * V1Mat
    b::SVector{r-1,Int} = 2 *( V1Mat' * VL.L.G * a.vec)
    γ::Int = a⊙a - k
    solutions = qsolve(A, b, γ)
    #println("qsolve ", V1Mat' * VL.L.G * V1Mat, " , ", 2 *(V1Mat' * VL.L.G * a.vec), " , " ,  a⊙a - k)
    #println(solutions)
    solutions_in_L::Array{HyperbolicLatticeElement,1} = (x -> HyperbolicLatticeElement(VL.L,V1Mat * x + a.vec)).(solutions)
     
    #@info "< roots_decomposed_into(VL, $(a.vec), $k)"
    return filter(is_root,solutions_in_L)
end

function enumerate_roots(VL;num=200,file=nothing)
    
    roots_by_distance = RootsByDistance(VL)
    
    roots = []
    last_fake_dist = 0
    open(file,"w") do io 

        while num > 0
            println("next root $num")
            r = next_with_decomposition!(roots_by_distance)
            #println("$(r.vec) : $((sinh_distance_to_hyperplane(VL.v0,r))^2)")
            @assert r[3] ≥ last_fake_dist
            last_fake_dist = r[3]
            push!(roots,r)
            println(io,"$(r[1][1]),$(r[1][2]),$(r[1][3].vec),$(r[2].vec),$(r[3])")
            num -= 1
        end
    end
    return roots
end

