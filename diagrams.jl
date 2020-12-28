using MLStyle
using LightGraphs
using GraphPlot
using Colors
using Multisets
using Memoize
using StaticArrays

# TODO
#
# * maybe implement DegSeq more efficiently
# * memoize connected_diagram_type by hand
# * make computation of build_deg_seq_and_associated_data more efficient by using previously constructed deg_seqs of the previous components?
# * clean up structure
# * store subdiagrams in arrays rather than a dict since all we do is iterate
#
# * be efficient enough to run on the ./graphs/18-vinbxx examples

import Base.push!, Base.length

function gug_coxiter_path_to_matrix(path)
    s = open(path) do file
        read(file, String)
    end

    return gug_coxiter_to_matrix(s)
end

function gug_coxiter_to_matrix(descr)
   
    lines = split(descr,"\n")
    @assert length(lines) ≥ 2 "CoxIter graph description must have non trivial content!"
    
    num_vertices = nothing
    rank = nothing
    m = match(r"^(?<num_vertices>\d\d*)\s\s*(?<rank>\d\d*)", lines[1])
    if m === nothing
        println("can't match first line (we require the rank to be explicitely given):")
        println(lines[1])
        return nothing 
    end

    num_vertices = parse(Int,m[:num_vertices])
    rank = parse(Int,m[:rank])
     
    D = fill(2,num_vertices,num_vertices)
    for i in 1:num_vertices
        D[i,i] = 1
    end
   
    
    m = match(r"^vertices labels:(?<all_labels>(\s+\w+)*)\s*", lines[2])
    if m ≠ nothing # labelled vertice
        vertices = split(m[:all_labels])
        if length(vertices) ≠ num_vertices
            println("not enough vertex labels")
            return nothing
        end
        
        if length(lines) == 2
            println("not enough lines")
            return nothing
        end
        for line in lines[3:end]
            m = match(r"^(?<from>\w+)\s+(?<to>\w+)\s+(?<label>\d+)", line)
            if m === nothing
                continue
            end
            from = indexin([m[:from]],vertices)[1]
            to = indexin([m[:to]],vertices)[1]
            label = parse(Int,m[:label])
            println("$from and $to and $label") 
            if from ≠ nothing && to ≠ nothing && from ≠ to
                D[from,to] = label
                D[to,from] = label
            end
        end
  

    else # non labelled vertices
        for line in lines[2:end]
            m = match(r"^(?<from>\d\d*)\s\s*(?<to>\d\d*)\s\s*(?<label>\d\d*)", line)
            if m === nothing
                continue
            end
            from = parse(Int,m[:from])
            to = parse(Int,m[:to])
            label = parse(Int,m[:label])
          
            if from ≠ to && from ∈ 1:num_vertices && to ∈ 1:num_vertices
                D[from,to] = label
                D[to,from] = label
            end
        end
    end


    
    return D, rank

end




# Referring to https://en.wikipedia.org/wiki/Coxeter%E2%80%93Dynkin_diagram
# These are the different isomorphism types of irreducible (=connected) subdiagrams of spherical/affine Coxeter-Dynkin diagrams
@enum DiagramType begin
    # Lowercase == Spherical
    DT_a
    DT_b # == c 
    DT_d

    DT_e6
    DT_e7
    DT_e8
    DT_f4
    DT_g2
    DT_h2
    DT_h3
    DT_h4
    DT_in
    # Uppercase == Affine
    DT_A
    DT_B
    DT_C
    DT_D

    DT_E6
    DT_E7
    DT_E8
    DT_F4
    DT_G2
    DT_Iinfty
end


@memoize function is_spherical(dt::DiagramType)::Bool
    dt ∈ [DT_a,DT_b,DT_d,DT_e6,DT_e7,DT_e8,DT_f4,DT_g2,DT_h2,DT_h3,DT_h4,DT_in]
end

@memoize function is_affine(dt::DiagramType)::Bool
    dt ∈ [DT_A,DT_B,DT_C,DT_D,DT_E6,DT_E7,DT_E8,DT_F4,DT_G2,DT_Iinfty]
end

@memoize function is_sporadic(dt::DiagramType)::Bool
    dt ∈ [DT_e6,DT_e7,DT_e8,DT_f4,DT_g2,DT_h2,DT_h3,DT_h4,DT_in,DT_E6,DT_E7,DT_E8,DT_F4,DT_G2,DT_Iinfty]
end



#0  -> ()
#1  -> (3)
#2  -> (3,3)
#3  -> (3,3,3)
#4  -> (3,3,3,3)
#5  -> (4)
#6  -> (4,3)
#7  -> (4,3,3)
#8  -> (4,4)
#9  -> (5)
#10 -> (5,3)
#11 -> (6)
#12 -> (6,3)
#13 -> (∞)
#13+n -> (n)
Deg = Int

push_three(n::Deg) = begin
    if n ≤ 3 || 5 ≤ n ≤ 6 || n == 9 || n == 11
        return n+1
    else
        return nothing
    end
end
push_four(n::Deg) = begin
    if n ≤ 2 
        return n+5
    elseif n == 5
        return 8
    else
        return nothing
    end
end
push_five(n::Deg) = begin
    if n ≤ 2
        return n+9
    else
        return nothing
    end
end
push_six(n::Deg) = begin
    if n ≤ 2
        return n+11
    else
        return nothing
    end
end
push_infty(n::Deg) = (n==0 ? 13 : nothing)
push_big(n::Deg,l::Int) = (n==0 ? 13+l : nothing)


# Degree Sequence: Each vertex has an associated "multi-degree" which is a multiset containing the labels of edges incident to the vertex
# A DegSeq is the multiset of the mult-degrees associated to the vertices of a diagram
mutable struct DegSeq 
    content::Vector{SVector{4,Int}}
end

function short_vec_to_svec(v::Vector{Int})::SVector{4,Int}
    @assert length(v) ≤ 4
    filled = vcat(v, fill(2,(4-length(v),1)))
    return SVector{4,Int}(filled)

end

function Base.push!(ds::DegSeq,v::Vector{Int})
    @assert length(v) ≤ 4
    push!(ds.content,sort(short_vec_to_svec(v)))
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
    sorted_vv = [sort(short_vec_to_svec(v)) for v in vv]
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





# A connected induced subdiagram.
# connected = irreducible
# We store the vertices and the type
struct ConnectedInducedSubDiagram
    vertices::BitSet
    type::DiagramType
end

function CISD(vertices::BitSet,type)
    ConnectedInducedSubDiagram(vertices,type)
end
function CISD(vertices::Array{Int,1},type)
    ConnectedInducedSubDiagram(BitSet(vertices),type)
end


card(c::ConnectedInducedSubDiagram) = length(c.vertices)

is_empty(c::ConnectedInducedSubDiagram) = length(c.vertices) == 0

the_singleton(v::Int) = CISD(BitSet(v),DT_a)


# An arbitrary induced subdiagram
# Stored as a collection of its irreducible components, plus whether it is affine or spherical
struct InducedSubDiagram
    connected_components::Vector{ConnectedInducedSubDiagram}
    is_affine::Bool
    is_spherical::Bool
end

function InducedSubDiagram(connected_components::Vector{ConnectedInducedSubDiagram})
    this_is_affine = all(is_affine(c.type) for c in connected_components) 
    this_is_spherical = all(is_spherical(c.type) for c in connected_components) 
    
    @assert ! (length(connected_components) > 0 && this_is_affine && this_is_spherical) "A diagram can't both be spherical and affine."
    
    return InducedSubDiagram(connected_components, this_is_affine, this_is_spherical) 
end

function the_empty_isd() 
    return InducedSubDiagram(Vector{ConnectedInducedSubDiagram}())
end


struct DiagramAndSubs 
    D::Array{Int,(2)}
    subs::Vector{Vector{Tuple{BitSet,InducedSubDiagram}}}
end


function print_das(das::DiagramAndSubs)
    
    println("The matrix is given by:")
    display(das.D)
    println()
    println("and the subdiagrams:")
    for i in eachindex(das.subs)
        println("cardinality $(i-1):")
        for (sub_support, sub_diagram) in das.subs[i]
            println("$(collect(sub_support)) [affine=$(sub_diagram.is_affine),spherical=$(sub_diagram.is_spherical)]:")
            for component in sub_diagram.connected_components
                println("    $(component.type) : $(component.vertices)")
            end
        end
        println("")
    end
    println("***********************")

end

function is_compact(dag::DiagramAndSubs,n::Int)

    subs = dag.subs
   
    # just to ensure that `subs` is long enough
    # TODO clean this
    subs = vcat(subs, [[] for i in 1:n])

    if ! any(is_spherical(subdiagram) for (support,subdiagram) in subs[n+1])
        return false
    end

    for (support, subdiagram) in subs[n]

        if  is_spherical(subdiagram)
            num_extensions = 0
            for (sup,sub) in subs[n+1]
                if support < sup && is_spherical(sub)
                    num_extensions += 1
                end
            end
            if num_extensions ≠ 2
                @debug "the subdiagram of support $support has $(length(extensions)) affine/spherical extensions"
                return false
            end
        end
    
    end
    

    return true

end

function is_finite_volume(dag::DiagramAndSubs,n::Int)

    subs = dag.subs

    subs = vcat(subs, [[] for i in 1:n])
    
    if ! any(is_spherical(subdiagram) || is_affine(subdiagram) for (support,subdiagram) in subs[n+1])
        return false
    end

    for (support, subdiagram) in subs[n]

        if  is_spherical(subdiagram)
            num_extensions = 0
            for (sup,sub) in subs[n+1]
                if support < sup && (is_spherical(sub) || is_affine(sub))
                    num_extensions += 1
                end
            end
            if num_extensions ≠ 2
                @debug "the subdiagram of support $support has $(length(extensions)) affine/spherical extensions"
                return false
            end
        end
    
    end
    

    return true

end



Arg = Tuple{BitSet,Array{Int,2},Bool}
Base.hash(a::Arg, h::UInt) = hash(a[1], hash(:Arg, h))
Base.isequal(a::Arg, b::Arg) = Base.isequal(hash(a), hash(b))


@memoize Dict function connected_diagram_type(VS::BitSet,D::Array{Int,2}; only_sporadic::Bool=false)
    
    @assert true "The diagram here is assumed connected. maybe this deserves a check"
    

    deg_seq_and_assoc = build_deg_seq_and_associated_data(VS,D)
    
    joined = nothing
    if  length(VS) ≤ 9 
        joined = connected_sporadic_diagram_type(VS,D,deg_seq_and_assoc) 
    end
    if joined === nothing && !only_sporadic
        joined = connected_non_sporadic_diagram_type(VS,D,deg_seq_and_assoc)
    end
    
    return joined
end

function build_deg_seq_and_associated_data(VS::BitSet,D::Array{Int,2})

    @assert true "The diagram here is assumed connected. maybe this deserves a check"

    @debug "build_deg_seq_and_associated_data(…)"
    @debug "VS is $(collect(VS))"

    VSC = length(VS)
    
    deg1_vertices = Dict{Int,BitSet}()
    deg3_vertices = Dict{Int,BitSet}()
    
    deg_seqs::DegSeq = deg_seq(Vector{Vector{Int}}()) 
    for v in VS

        @debug "looking at $v"

        deg_seq_v = Vector{Int}() 
        simple_neighbors_v = BitSet()
        for u in VS if u ≠ v
            if D[u,v] == 3
                push!(simple_neighbors_v,u)
            end
            if D[u,v] ≠ 2 
                push!(deg_seq_v,D[u,v])
            end
        end end
        if deg_seq_v == [3,3,3]
            push!(deg3_vertices,v => simple_neighbors_v)
        elseif deg_seq_v == [3]
            push!(deg1_vertices,v => simple_neighbors_v)
        end
        push!(deg_seqs, deg_seq_v)
        
        @debug "$v has deg_seq = $deg_seq_v"
    end    
    ds = deg_seqs
    
    @debug "deg_seq is $ds versus"
    @debug "           $deg_seq_f4"
    @debug "end of build_deg_seq…"

    center::Union{Nothing,Int} = ( length(collect(keys(deg3_vertices))) > 0 ? collect(keys(deg3_vertices))[1] : nothing )
    center_neighbors::BitSet = ( center === nothing ? BitSet() : deg3_vertices[center] )
    extremities::BitSet = BitSet(keys(deg1_vertices))
    extremities_neighbors::BitSet = ( isempty(extremities) ? BitSet() : ∪(values(deg1_vertices)...) )
    
    return (ds, deg1_vertices, deg3_vertices, center, center_neighbors, extremities, extremities_neighbors)::Tuple{DegSeq,Dict{Int,BitSet},Dict{Int,BitSet},Union{Nothing,Int},BitSet,BitSet,BitSet}

end

function connected_non_sporadic_diagram_type(VS::BitSet,D::Array{Int,2},deg_seq_and_assoc::Tuple{DegSeq,Dict{Int,BitSet},Dict{Int,BitSet},Union{Nothing,Int},BitSet,BitSet,BitSet})
    
    @assert true "The diagram here is assumed connected. maybe this deserves a check"
    
    (ds, deg1_vertices, deg3_vertices, center, center_neighbors, extremities, extremities_neighbors) = deg_seq_and_assoc # build_deg_seq_and_associated_data(VS,D) 
    

    vertices = VS
    n = length(vertices)

    if false
        @assert false "+++"
    elseif ds == deg_seq_a1
        return CISD(vertices,DT_a)
    elseif n≥2 && ds == deg_seq_a(n)
        return CISD(vertices,DT_a)
    
    elseif ds == deg_seq_b2
        return CISD(vertices,DT_b)
    elseif ds == deg_seq_b3
        return CISD(vertices,DT_b)
    elseif n≥4 && ds == deg_seq_b(n)
        return CISD(vertices,DT_b)
    
    elseif  n≥4 && ds == deg_seq_d(n)  && length(extremities ∩ center_neighbors) ≥ 2
        return CISD(vertices,DT_d)

    elseif n≥3 && ds == deg_seq_A(n)    
        return CISD(vertices,DT_A)

    elseif ds == deg_seq_B4
        return CISD(vertices,DT_B)
    elseif n≥5 && ds == deg_seq_B(n)
        return CISD(vertices,DT_B)

    elseif ds == deg_seq_C3
        return CISD(vertices,DT_C)
    elseif n≥4 && ds == deg_seq_C(n)
        return CISD(vertices,DT_C)

    elseif ds == deg_seq_D5
        return CISD(vertices,DT_D)
    elseif n≥6 && ds == deg_seq_D(n)
        return CISD(vertices,DT_D)
    else
        return nothing

    end    
end

function connected_sporadic_diagram_type(VS::BitSet,D::Array{Int,2},deg_seq_and_assoc::Tuple{DegSeq,Dict{Int,BitSet},Dict{Int,BitSet},Union{Nothing,Int},BitSet,BitSet,BitSet})

    @assert true "The diagram here is assumed connected. maybe this deserves a check"


    (ds, deg1_vertices, deg3_vertices, center, center_neighbors, extremities, extremities_neighbors) = deg_seq_and_assoc # build_deg_seq_and_associated_data(VS,D) 


    if false
        @assert false "For alignment's sake"
      
    # les sporadiques **pas** de type DT_e ou DT_E
    elseif ds == deg_seq_f4
        return CISD(VS,DT_f4)
    elseif ds == deg_seq_F4
        return CISD(VS,DT_F4)
    elseif ds == deg_seq_h2
        return CISD(VS,DT_h2)
    elseif ds == deg_seq_h3
        return CISD(VS,DT_h3)
    elseif ds == deg_seq_h4
        return CISD(VS,DT_h4)
    elseif ds == deg_seq_g2
        return CISD(VS,DT_g2)
    elseif ds == deg_seq_G2
        return CISD(VS,DT_G2)
    elseif ds == deg_seq_Iinfty
        return CISD(VS,DT_Iinfty)
    elseif length(ds) ≥ 1 && length(ds.content[1]) ≥ 1 && ds.content[1][1] ≥ 6 && ds == deg_seq_i(ds.content[1][1])
        return CISD(VS,DT_in)


    elseif ds == deg_seq_e6 && length(center_neighbors∩extremities) == 1    
        return CISD(VS,DT_e6)
    elseif ds == deg_seq_e7 && length(center_neighbors∩extremities) == 1   
        return CISD(VS,DT_e7)
    elseif ds == deg_seq_e8 && length(center_neighbors∩extremities) == 1 && length(extremities_neighbors ∩ center_neighbors) == 1 
        return CISD(VS,DT_e8)
    elseif ds == deg_seq_E6 && length(center_neighbors∩extremities) == 0    
        return CISD(VS,DT_E6)
    elseif ds == deg_seq_E7 && length(center_neighbors∩extremities) == 1 && length(extremities_neighbors ∩ center_neighbors) == 0 
        return CISD(VS,DT_E7)
    elseif ds == deg_seq_E8 && length(center_neighbors∩extremities) == 1 && length(extremities_neighbors ∩ center_neighbors) == 1 
        return CISD(VS,DT_E8)
    end
    
    return nothing

end

function try_extend(VS::BitSet,S::InducedSubDiagram,D::Array{Int,2},v::Int)
    components = S.connected_components
    
    # joined should/will be of type ConnectedInducedSubDiagram
    joined = nothing # Here is the result
    
    # special case, no component in S, so we just return the singleton
    if length(components) == 0
        joined = the_singleton(v)::ConnectedInducedSubDiagram
        only_joined = Vector{ConnectedInducedSubDiagram}()::Vector{ConnectedInducedSubDiagram}
        push!(only_joined, joined)
        return InducedSubDiagram(only_joined)
    end
    
    freedom = 4

    neighboring_components::Vector{ConnectedInducedSubDiagram} = Vector{ConnectedInducedSubDiagram}()
    #neighboring_components_size::BitSet = BitSet()
    non_neighboring_components::Vector{ConnectedInducedSubDiagram} = Vector{ConnectedInducedSubDiagram}()

    total_size::Int = 1
    only_sporadic::Bool = false
    
    for c in components
        for u in c.vertices
            if D[u,v] == 1 # dotted edge => early out
                return nothing
            elseif D[u,v] ≠ 2 
                if is_affine(c.type) # can't extend affine types => early out
                    return nothing
                end
                if freedom ≤ 0  # can't have too many connections, neither too high degrees
                    return nothing
                end
                if D[u,v] == 0
                    freedom = 0
                else
                    freedom -= (min(D[u,v],5)-2)
                end
                if length(neighboring_components) == 0 || c≠neighboring_components[end]
                    push!(neighboring_components,c)
                end
                total_size += card(c)
                if is_sporadic(c.type)
                    only_sporadic = true
                end
                if D[u,v] > 4
                    only_sporadic = true
                end
            end
        end
        if length(neighboring_components) == 0 ||  c ≠ neighboring_components[end]
            push!(non_neighboring_components,c)
        end

    end

    
    vertices::BitSet = ∪(BitSet(v),[c.vertices for c in neighboring_components]...)

    joined = connected_diagram_type(vertices,D;only_sporadic=only_sporadic)

    if joined === nothing
        return nothing
    else
        new_components = non_neighboring_components 
        push!(new_components,joined)
        return InducedSubDiagram(new_components)
    end

end


is_affine(isd::InducedSubDiagram) = all(is_affine(c.type) for c in isd.connected_components)
is_spherical(isd::InducedSubDiagram) = all(is_spherical(c.type) for c in isd.connected_components)

function extend(das::DiagramAndSubs, v::Array{Int,1}; max_card::Union{Nothing,Int}=nothing)
  
    D = das.D
    old_subs = das.subs
    max_card = (max_card === nothing ?  length(v)+1 : max_card) # if set to nothing, give it high enough value that it won't restrict anything

    n = length(v) 
    @assert size(D) == (n,n)
    
    # Extend D with v
    D = [D v]
    D = [D; [v;1]']
   
    new_vertex = n+1
    
    new_subs::Vector{Vector{Tuple{BitSet,InducedSubDiagram}}} = [[] for i in eachindex(old_subs)]
    for i in eachindex(old_subs)
        if i ≤ max_card 
            for (support,subdiagram) in old_subs[i]
                extended_with_v = try_extend(support,subdiagram,D,new_vertex)
                if extended_with_v ≠ nothing 
                    push!(new_subs[i],(support∪BitSet(new_vertex),extended_with_v))
                end
            end
        end 
    end 
   
    new_aligned::Vector{Vector{Tuple{BitSet,InducedSubDiagram}}} = vcat([[]], new_subs) 
    old_aligned::Vector{Vector{Tuple{BitSet,InducedSubDiagram}}} = vcat(old_subs, [[]])
    return DiagramAndSubs(D,[vcat(old_aligned[i],new_aligned[i]) for i in eachindex(old_aligned)])
    
end


function build_diagram_and_subs(M::Array{Int,2};max_card::Union{Nothing,Int}=nothing)
   
    empty!(memoize_cache(connected_diagram_type))

    n = size(M,1)
    @assert size(M) == (n,n) "M must be square"
    @assert M == M' "M must be symmetric"
    @assert all(l ≥ 0 for l in M) "M must have non-negative entries"

    subs = Vector{Vector{Tuple{BitSet,InducedSubDiagram}}}([[]])
    push!(subs[1],(BitSet(), the_empty_isd()))

    das = DiagramAndSubs(reshape([],0,0),subs)
    for i in 1:n
        println("extending with vertex $i")
        das = extend(das,M[i,1:i-1];max_card=max_card)
        #print_das(das)
    end
    return das
end



# ##########################
#
# Misc user-facing functions
#
# ##########################
#

function is_compact_respectively_finvol(path::String)
    
    # Get file content in s as in https://en.wikibooks.org/wiki/Introducing_Julia/Working_with_text_files
    s = open(path) do file
        read(file, String)
    end

    ret = gug_coxiter_to_matrix(s)
    if ret === nothing
        println("Failed reading $path")
    else
        (D, rank) = ret
        if D === nothing || rank === nothing
            println("Error reading file probably")
        else
            das = build_diagram_and_subs(D;max_card=rank)
            print_das(das)
            return (is_compact(das, rank), is_finite_volume(das, rank))
        end
    end


end



function check_all_graphs(sub_directory="")
    

    for (root, dirs, files) in walkdir("./graphs/"*sub_directory)
        for path in joinpath.(root, files)
            if endswith(path,".coxiter")
                println("path: $path")
                println(is_compact_respectively_finvol(path))
            end
        end
    end

end


function check_some_graphs()

    @time begin
        for path in ["graphs/13-mcl11.coxiter"] 
                println("path: $path")
                println(is_compact_respectively_finvol(path))
        end
        for (root, dirs, files) in walkdir("./graphs/simplices")
            for path in joinpath.(root, files)
                if endswith(path,".coxiter")
                    println("path: $path")
                    println(is_compact_respectively_finvol(path))
                end
            end
        end
    end

end
