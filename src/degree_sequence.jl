
Deg = Int

empty_deg = 0


#=
#   0  -> []
#   1  -> [3]
#   2  -> [3,3]
#   3  -> [3,3,3]
#   4  -> [3,3,3]
#   5  -> [4]
#   6  -> [4,3]
#   7  -> [4,3,3]
#   8  -> [4,4]
#   9  -> [5]
#   10 -> [5,3]
#   11 -> [6]
#   12 -> [6,3]
#   13 -> [∞]
#   n  -> [n-13]
#
#
#
=#


@inline push_three(n::Deg) = begin
    if n ≤ 3 || 5 ≤ n ≤ 6 || n == 9 || n == 11
        return n+1
    else
        return nothing
    end
end
@inline push_four(n::Deg) = begin
    if n ≤ 2 
        return n+5
    elseif n == 5
        return 8
    else
        return nothing
    end
end
@inline push_five(n::Deg) = begin
    if n ≤ 2
        return n+9
    else
        return nothing
    end
end
@inline push_six(n::Deg) = begin
    if n ≤ 2
        return n+11
    else
        return nothing
    end
end
@inline push_infty(n::Deg) = (n==0 ? 13 : nothing)
@inline push_big(n::Deg,l::Int) = (n==0 ? 13+l : nothing)

@inline function push_label(n::Deg,l::Int)::Union{Deg,Nothing} 
    if l == 3
        return push_three(n)
    elseif l == 4
        return push_four(n)
    elseif l == 5
        return push_five(n)
    elseif l == 6 
        return push_six(n)
    elseif l == 0
        return push_infty(n)
    elseif l ≥ 7
        return push_big(n,l)
    else
        return nothing
    end
end

function big_label(n::Deg)::Union{Int,Nothing}
    if n ≥ 20
        return n-13
    else
        return nothing
    end
end

# Degree Sequence: Each vertex has an associated "multi-degree" which is a multiset containing the labels of edges incident to the vertex, encoded as an int
# A DegSeq is the multiset of the mult-degrees associated to the vertices of a diagram
mutable struct DegSeq 
    content::Vector{Deg}
end

function short_vec_to_deg(vv::Vector{Int})::Union{Int,Nothing}
    @assert length(vv) ≤ 4
    res = empty_deg
    for v in vv
        res = push_label(res,v)
        if res === nothing
            return nothing
        end
    end
    return res 

end

function Base.push!(ds::DegSeq,v::Vector{Int})
    @assert length(v) ≤ 4
    push!(ds.content,short_vec_to_deg(v))
    sort!(ds.content)
    return ds
end

function Base.push!(ds::DegSeq,v::Deg)
    @assert length(v) ≤ 4
    push!(ds.content,v)
    sort!(ds.content)
    return ds
end

function Base.:+(ds1::DegSeq,ds2::DegSeq)
    DegSeq(sort(vcat(ds1.content,ds2.content)))
end

function Base.:*(k::Int,ds::DegSeq)
    @assert k≥0
    DegSeq(reduce(vcat,[[v for i in 1:k] for v in ds.content]))
end

function Base.length(ds::DegSeq)
    length(ds.content)
end

function Base.:(==)(ds1::DegSeq,ds2::DegSeq)
    (ds1.content == ds2.content)
end

# Given an array of arrays of ints, returns the associated degree sequence
function deg_seq(vv::Vector{Vector{Int}})
    @assert all(length(v) ≤ 4 for v in vv)
    sorted_vv = [short_vec_to_deg(v) for v in vv]
    return DegSeq(sort(sorted_vv))
end

# The degree sequences corresponding to each irreducible diagram types follow:

const deg_seq_a1 = deg_seq(Vector{Vector{Int}}([[]]))
deg_seq_a(n::Int) = begin
    @assert n≥2
    2*deg_seq([[3]]) + (n-2)*deg_seq([[3,3]])::DegSeq
end

const deg_seq_b2 = deg_seq([[4],[4]])
const deg_seq_b3 = deg_seq([[4],[4,3],[3]])
deg_seq_b(n)::DegSeq = begin
    @assert n≥3
    deg_seq([[4]]) + deg_seq([[4,3]]) + (n-3)*deg_seq([[3,3]]) + deg_seq([[3]])
end

deg_seq_d(n::Int)::DegSeq = begin
    @assert n≥4
    deg_seq([[3,3,3]]) + (n-4)*deg_seq([[3,3]]) + 3*deg_seq([[3]])
end

deg_seq_A(n::Int)::DegSeq = begin
    @assert n≥3
    n*deg_seq([[3,3]])
end

deg_seq_B4 = begin
    deg_seq([[3,3,4],]) + 2*deg_seq([[3],]) + deg_seq([[4],])
end

deg_seq_B(n::Int)::DegSeq = begin
    @assert n≥5
    deg_seq([[3,3,3]]) + 2*deg_seq([[3]]) + (n-5)*deg_seq([[3,3]]) + deg_seq([[3,4]])  + deg_seq([[4]])
end
deg_seq_C3 = begin
    deg_seq([[4,4]]) +  2*deg_seq([[4]])   
end
deg_seq_C(n::Int)::DegSeq = begin
    @assert n≥4
    2*deg_seq([[4,3]]) +  2*deg_seq([[4]])  + (n-4)*deg_seq([[3,3]])
end

const deg_seq_D5 = begin
    deg_seq([[3,3,3,3]]) + 4*deg_seq([[3]])
end
deg_seq_D(n::Int)::DegSeq = begin
    @assert n≥6
    2*deg_seq([[3,3,3]]) + 4*deg_seq([[3]]) + (n-6)*deg_seq([[3,3]])
end


const deg_seq_f4 = deg_seq([[3],[3],[3,4],[3,4]])
const deg_seq_F4 = deg_seq([[3],[3],[3,3],[3,4],[3,4]])
const deg_seq_h2 = deg_seq([[5],[5]])
const deg_seq_h3 = deg_seq([[5],[5,3],[3]])
const deg_seq_h4 = deg_seq([[5],[5,3],[3,3],[3]])
const deg_seq_g2 = deg_seq([[6],[6]])
const deg_seq_G2 = deg_seq([[6],[6,3],[3]])
const deg_seq_Iinfty = deg_seq([[0],[0]])
function deg_seq_i(n::Int)::DegSeq
    deg_seq([[n],[n]])
end

const deg_seq_e6 = deg_seq([[3,3,3],[3,3],[3,3],[3],[3],[3]])
const deg_seq_e7 = deg_seq([[3,3,3],[3,3],[3,3],[3,3],[3],[3],[3]])
const deg_seq_e8 = deg_seq([[3,3,3],[3,3],[3,3],[3,3],[3,3],[3],[3],[3]])
const deg_seq_E6 = deg_seq([[3,3,3],[3,3],[3,3],[3,3],[3],[3],[3]])
const deg_seq_E7 = deg_seq([[3,3,3],[3,3],[3,3],[3,3],[3,3],[3],[3],[3]])
const deg_seq_E8 = deg_seq([[3,3,3],[3,3],[3,3],[3,3],[3,3],[3,3],[3],[3],[3]])



