"""
    GenDeg 

In an integer-labelled graph, we def "generalized degree" of a vertex as the multiset of the labels of its incident edges.
In the case of **irreducible spherical/affine** Coxeter diagrams, this (generalized) degree can't take "many" forms, so we encode it as an integer, as follows:
```
    Int -> MultiSet
    0  -> []
    1  -> [3]
    2  -> [3,3]
    3  -> [3,3,3]
    4  -> [3,3,3]
    5  -> [4]
    6  -> [4,3]
    7  -> [4,3,3]
    8  -> [4,4]
    9  -> [5]
    10 -> [5,3]
    11 -> [6]
    12 -> [6,3]
    13 -> [∞]
    n  -> [n-13]
```    
"""
GenDeg = Int

"""
    empty_deg = 0

By convention, to `0` corresponds an isolated vertex.
"""
const empty_deg = 0




"""
    push_three(n::GenDeg)::Union{Nothing,GenDeg}

Add an edge with label three to a (generalized) degree, if it results in a legal degree for a spherical/affine irreducible diagram.
Otherwise, return `nothing`
"""
@inline push_three(n::GenDeg) = begin
    if n ≤ 3 || 5 ≤ n ≤ 6 || n == 9 || n == 11
        return n+1
    else
        return nothing
    end
end

"""
    push_four(n::GenDeg)::Union{Nothing,GenDeg}

Add an edge with label four to a degree, if it results in a legal degree for a spherical/affine irreducible diagram.
Otherwise, return `nothing`
"""
@inline push_four(n::GenDeg) = begin
    if n ≤ 2 
        return n+5
    elseif n == 5
        return 8
    else
        return nothing
    end
end

"""
    push_five(n::GenDeg)::Union{Nothing,GenDeg}

Add an edge with label five to a degree, if it results in a legal degree for a spherical/affine irreducible diagram.
Otherwise, return `nothing`
"""
@inline push_five(n::GenDeg) = begin
    if n ≤ 2
        return n+9
    else
        return nothing
    end
end

"""
    push_six(n::GenDeg)::Union{Nothing,GenDeg}

Add an edge with label six to a degree, if it results in a legal degree for a spherical/affine irreducible diagram.
Otherwise, return `nothing`
"""
@inline push_six(n::GenDeg) = begin
    if n ≤ 2
        return n+11
    else
        return nothing
    end
end

"""
    push_infty(n::GenDeg)::Union{Nothing,GenDeg}

Add an edge with label infty to a degree, if it results in a legal degree for a spherical/affine irreducible diagram.
Otherwise, return `nothing`
"""
@inline push_infty(n::GenDeg) = (n==0 ? 13 : nothing)

"""
    push_big(n::GenDeg,l::Int)::Union{Nothing,GenDeg}

Add an edge with label ``l>6`` to a degree, if it results in a legal degree for a spherical/affine irreducible diagram.
Otherwise, return `nothing`
"""
@inline push_big(n::GenDeg,l::Int) = (n==0 ? 13+l : nothing)

"""
    push_label(n::GenDeg,l::Int)::Union{Nothing,GenDeg}

Add an edge with label ``l`` to a degree, if it results in a legal degree for a spherical/affine irreducible diagram.
Otherwise, return `nothing`
"""
@inline function push_label(n::GenDeg,l::Int)::Union{GenDeg,Nothing} 
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

"""
    big_label(n::GenDeg)

If `n` encodes a (generalized) degree corresponding to just an edge with label ``∞>l>6``, return the corresponding ``l``, and `nothing` otherwise.
"""
function big_label(n::GenDeg)::Union{Int,Nothing}
    if n ≥ 20
        return n-13
    else
        return nothing
    end
end

"""
    GenDegSeq

The "degree sequence" of a combinatorial graph is the multiset (or ordered sequence) of the degree of its vertices.
For (irreducible spherical/affine) Coxeter diagram, `GenDegSeq` is similarly defined as the sorted vector of "generalized" degrees ([GenDeg](@ref)).
"""
mutable struct GenDegSeq 
    content::Vector{GenDeg}
end

"""
    short_vec_to_deg(vv::Vector{Int})

Encode the "generalized degree" given by `vv` into an object of type `GenDeg` if possible (if the resulting generalized degree can appear in a spherical/affine irreducible diagram).
"""
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

"""
    push!(ds,v)

Given a generalized degree sequence `ds` and a vector `v` defining a generalized degree, add `v` to the sequence `ds`.
"""
function Base.push!(ds::GenDegSeq,v::Vector{Int})
    @assert length(v) ≤ 4
    push!(ds.content,short_vec_to_deg(v))
    sort!(ds.content)
    return ds
end

"""
    push!(ds,v)

Given a generalized degree sequence `ds` and a `GenDeq` encoding a generalized degree, add `v` to the sequence `ds`.
"""
function Base.push!(ds::GenDegSeq,v::GenDeg)
    @assert length(v) ≤ 4
    push!(ds.content,v)
    sort!(ds.content)
    return ds
end

"""
    +(ds1,ds2)

Concatenate degree sequences `ds1` and `ds2`.
"""
function Base.:+(ds1::GenDegSeq,ds2::GenDegSeq)
    GenDegSeq(sort(vcat(ds1.content,ds2.content)))
end

"""
    *(k,ds)

Concatenate `k` copies of the degree sequence `ds` together.
"""
function Base.:*(k::Int,ds::GenDegSeq)
    @assert k≥0
    GenDegSeq(reduce(vcat,[[v for i in 1:k] for v in ds.content]))
end

"""
    length(ds)

The length of a degree sequence (corresponding to the number of vertices in the diagram it represents).
"""
function Base.length(ds::GenDegSeq)
    length(ds.content)
end

"""
    ==(ds1,ds2)

Whether `ds1` and `ds2` are equal degree sequence.
Note that since we keep them sorted (as for "usual" degree sequence), equality testing is just equality of vectors.
"""
function Base.:(==)(ds1::GenDegSeq,ds2::GenDegSeq)
    (ds1.content == ds2.content)
end

"""
    deg_seq(vv)

Given the `Vector` of `Vector` of `Int`s `vv`, return the associated degree sequence (each vector of integers represents a generalized degree).
"""
function deg_seq(vv::Vector{Vector{Int}})
    @assert all(length(v) ≤ 4 for v in vv)
    sorted_vv = [short_vec_to_deg(v) for v in vv]
    return GenDegSeq(sort(sorted_vv))
end




# We can now define the degree sequences corresponding to each possible irreducible spherical/affine diagram
# Note that for instance `deg_seq_A(n)` corresponds to the diagram `DT_A` with ``n`` vertices, which is usually called ``\tilde{A}_{n-1}``, hence a difference of 1 for all indices in affine diagrams.
# Let's document only `deg_seq_a(n::Int)`.

"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_a1 = deg_seq(Vector{Vector{Int}}([[]]))

"""
    deg_seq_a(n)

Compute the generalized degree sequence associated to the spherical diagram of type A with ``n`` vertices.
Beware that for all sporadic diagrams, the specified integer corresponds to the _rank_ of the diagram, while for affine non sporadic, it corresponds to the number of vertices.
"""
deg_seq_a(n::Int) = begin
    @assert n≥2
    2*deg_seq([[3]]) + (n-2)*deg_seq([[3,3]])::GenDegSeq
end

"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_b2 = deg_seq([[4],[4]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_b3 = deg_seq([[4],[4,3],[3]])
"""
Similar as [`deg_seq_a`](@ref)
"""
deg_seq_b(n)::GenDegSeq = begin
    @assert n≥3
    deg_seq([[4]]) + deg_seq([[4,3]]) + (n-3)*deg_seq([[3,3]]) + deg_seq([[3]])
end

"""
Similar as [`deg_seq_a`](@ref)
"""
deg_seq_d(n::Int)::GenDegSeq = begin
    @assert n≥4
    deg_seq([[3,3,3]]) + (n-4)*deg_seq([[3,3]]) + 3*deg_seq([[3]])
end

"""
Similar as [`deg_seq_a`](@ref)
"""
deg_seq_A(n::Int)::GenDegSeq = begin
    @assert n≥3
    n*deg_seq([[3,3]])
end

"""
Similar as [`deg_seq_a`](@ref)
"""
deg_seq_B4 = begin
    deg_seq([[3,3,4],]) + 2*deg_seq([[3],]) + deg_seq([[4],])
end

"""
Similar as [`deg_seq_a`](@ref)
"""
deg_seq_B(n::Int)::GenDegSeq = begin
    @assert n≥5
    deg_seq([[3,3,3]]) + 2*deg_seq([[3]]) + (n-5)*deg_seq([[3,3]]) + deg_seq([[3,4]])  + deg_seq([[4]])
end
"""
Similar as [`deg_seq_a`](@ref)
"""
deg_seq_C3 = begin
    deg_seq([[4,4]]) +  2*deg_seq([[4]])   
end
"""
Similar as [`deg_seq_a`](@ref)
"""
deg_seq_C(n::Int)::GenDegSeq = begin
    @assert n≥4
    2*deg_seq([[4,3]]) +  2*deg_seq([[4]])  + (n-4)*deg_seq([[3,3]])
end

"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_D5 = begin
    deg_seq([[3,3,3,3]]) + 4*deg_seq([[3]])
end
"""
Similar as [`deg_seq_a`](@ref)
"""
deg_seq_D(n::Int)::GenDegSeq = begin
    @assert n≥6
    2*deg_seq([[3,3,3]]) + 4*deg_seq([[3]]) + (n-6)*deg_seq([[3,3]])
end


"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_f4 = deg_seq([[3],[3],[3,4],[3,4]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_F4 = deg_seq([[3],[3],[3,3],[3,4],[3,4]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_h2 = deg_seq([[5],[5]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_h3 = deg_seq([[5],[5,3],[3]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_h4 = deg_seq([[5],[5,3],[3,3],[3]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_g2 = deg_seq([[6],[6]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_G2 = deg_seq([[6],[6,3],[3]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_I∞ = deg_seq([[0],[0]])
"""
Similar as [`deg_seq_a`](@ref)
"""
function deg_seq_i(n::Int)::GenDegSeq
    deg_seq([[n],[n]])
end

const deg_seq_e6 = deg_seq([[3,3,3],[3,3],[3,3],[3],[3],[3]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_e7 = deg_seq([[3,3,3],[3,3],[3,3],[3,3],[3],[3],[3]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_e8 = deg_seq([[3,3,3],[3,3],[3,3],[3,3],[3,3],[3],[3],[3]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_E6 = deg_seq([[3,3,3],[3,3],[3,3],[3,3],[3],[3],[3]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_E7 = deg_seq([[3,3,3],[3,3],[3,3],[3,3],[3,3],[3],[3],[3]])
"""
Similar as [`deg_seq_a`](@ref)
"""
const deg_seq_E8 = deg_seq([[3,3,3],[3,3],[3,3],[3,3],[3,3],[3,3],[3],[3],[3]])



