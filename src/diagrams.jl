using MLStyle
using LightGraphs
using GraphPlot
using Colors
using Multisets
using Memoize
using StaticArrays
using Debugger

include("sbitset.jl")

# TODO
#
# * maybe implement DegSeq more efficiently
# * memoize connected_diagram_type by hand
# * make computation of build_deg_seq_and_associated_data more efficient by using previously constructed deg_seqs of the previous components?
# * clean up structure
#
# * be efficient enough to run on the ./graphs/18-vinbxx examples

import Base.push!, Base.length, Base.copy


include("diagram_type.jl")
include("degree_sequence.jl")
include("coxiter_io.jl")




# A connected induced subdiagram.
# connected = irreducible
# We store the vertices and the type
struct ConnectedInducedSubDiagram
    vertices::SBitSet{1}
    type::DiagramType
end
function Base.:(==)(c1::ConnectedInducedSubDiagram,c2::ConnectedInducedSubDiagram)
    (c1.vertices == c2.vertices)
end

function CISD(vertices::SBitSet{1},type)
    ConnectedInducedSubDiagram(vertices,type)
end
function CISD(vertices::Array{Int,1},type)
    ConnectedInducedSubDiagram(SBitSet{1}(vertices),type)
end

function Base.copy(c::ConnectedInducedSubDiagram)
    return CISD(copy(c.vertices),c.type) 
end

card(c::ConnectedInducedSubDiagram) = length(c.vertices)

is_empty(c::ConnectedInducedSubDiagram) = length(c.vertices) == 0

the_singleton(v::Int) = CISD(SBitSet{1}(v),DT_a)


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


is_affine(isd::InducedSubDiagram) = all(is_affine(c.type) for c in isd.connected_components)
is_spherical(isd::InducedSubDiagram) = all(is_spherical(c.type) for c in isd.connected_components)

function Base.copy(c::InducedSubDiagram)
    return InducedSubDiagram(copy(c.connected_components),c.is_affine,c.is_spherical)
end


function the_empty_isd() 
    return InducedSubDiagram(Vector{ConnectedInducedSubDiagram}())
end


struct DiagramAndSubs 
    D::Array{Int,(2)}
    subs::Vector{Vector{Tuple{SBitSet{1},InducedSubDiagram}}}
end

function subdiagrams_of_card(das::DiagramAndSubs,n::Int)
    if n+1 > length(das.subs)
        return []
    else
        return das.subs[n+1]
    end
end

function dump_das(das::DiagramAndSubs;range=nothing)
   
    dense_bitset_str(b::SBitSet{1}) = *("[",[string(i)*"," for i in b]...,"]")

    println("### Subdiagrams for Coxeter matrix:")
    display(das.D)
    println()
    println("###")
    for i in eachindex(das.subs)
        if range === nothing || i ∈ range 
            println("Cardinality $(i-1):")
            for (sub_support, sub_diagram) in das.subs[i]
                print("    $(dense_bitset_str(sub_support)) $(sub_diagram.is_affine ? "A" : " ") $(sub_diagram.is_spherical ? "S" : " ") = ")
                for component in sub_diagram.connected_components
                    print("$(component.type)$(dense_bitset_str(component.vertices)) ∪ ")
                end
                println()
            end
            println()
        end
    end
    println("###")

end


function print_das(das::DiagramAndSubs;rank=nothing)
    
    println("The matrix is given by:")
    display(das.D)
    println()
    println("and the subdiagrams:")
    for i in eachindex(das.subs)
        if rank == nothing || i == rank || i == rank+1
            println("cardinality $(i-1):")
            for (sub_support, sub_diagram) in das.subs[i]
                println("$(collect(sub_support)) [affine=$(sub_diagram.is_affine),spherical=$(sub_diagram.is_spherical)]:")
                for component in sub_diagram.connected_components
                    println("    $(component.type) : $(component.vertices)")
                end
            end
            println("")
        end
    end
    println("***********************")

end

function is_compact(dag::DiagramAndSubs,dimension::Int)
    
    # The dimension is the rank
    #

   
    if ! any(is_spherical(subdiagram) for (support,subdiagram) in subdiagrams_of_card(dag,dimension))
        return false
    end

    for (support, subdiagram) in subdiagrams_of_card(dag,dimension-1)

        if  is_spherical(subdiagram)
            num_extensions = 0
            for (sup,sub) in subdiagrams_of_card(dag,dimension)
                if support ⊆ sup && is_spherical(sub)
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

function is_finite_volume(dag::DiagramAndSubs,dimension::Int)

    
    if ! any(is_spherical(subdiagram) || is_affine(subdiagram) for (support,subdiagram) in subdiagrams_of_card(dag,dimension))
        return false
    end

    for (support, subdiagram) in subdiagrams_of_card(dag,dimension-1)

        if  is_spherical(subdiagram)
            num_extensions = 0
            extensions = []
            for (sup,sub) in subdiagrams_of_card(dag,dimension)
                if support ⊆ sup && (is_spherical(sub) || is_affine(sub))
                    num_extensions += 1
                    push!(extensions,sub)
                end
            end
            if num_extensions ≠ 2
                println("the $(is_spherical(subdiagram) ? "spherical " : "") subdiagram $support has extensions")
                for e in extensions
                    println("$e   $(is_spherical(e) ? "(spherical)" : "")$(is_affine(e) ? "(spherical)" : "")" )
                end
                println("-----------------------------------------------------------")
                @debug "the subdiagram of support $support has $(length(extensions)) affine/spherical extensions"
                return false
            end
        end
    
    end
    

    return true

end



# dirty hack to memoize `connected_diagram_type` only on its first argument
Arg = Tuple{SBitSet{1},Array{Int,2},Bool}
Base.hash(a::Arg, h::UInt) = hash(a[1], hash(:Arg, h))
Base.isequal(a::Arg, b::Arg) = Base.isequal(hash(a), hash(b))


@memoize Dict function connected_diagram_type(VS::SBitSet{1},D::Array{Int,2}; only_sporadic::Bool=false)
    
    @assert true "The diagram here is assumed connected. maybe this deserves a check"


    deg_seq_and_assoc = build_deg_seq_and_associated_data(VS,D)
   
    if deg_seq_and_assoc === nothing
        return nothing
    end
    

    joined = nothing
    if  length(VS) ≤ 9
        joined = connected_sporadic_diagram_type(VS,D,deg_seq_and_assoc) 
    end
    if joined === nothing && !only_sporadic
        joined = connected_non_sporadic_diagram_type(VS,D,deg_seq_and_assoc)
    end
 


    return joined
end

@inline function build_deg_seq_and_associated_data(VS::SBitSet{1},D::Array{Int,2})

    
    @assert true "The diagram here is assumed connected. maybe this deserves a check"

    @debug "build_deg_seq_and_associated_data(…)"
    @debug "VS is $(collect(VS))"

    VSC = length(VS)
   
    deg1 = SBitSet{1}()
    deg1_neigh = SBitSet{1}()
    deg1_vertices = Dict{Int,SBitSet{1}}()
    deg3_vertices = Dict{Int,SBitSet{1}}()
   

    deg_seqs::DegSeq = DegSeq(Vector{Deg}()) 
    for v in VS

        @debug "looking at $v"

        deg_v = empty_deg 
        simple_neighbors_v = SBitSet{1}()
        for u in VS if u ≠ v
            if D[u,v] == 3
                simple_neighbors_v = simple_neighbors_v| SBitSet{1}(u)
            end
            if D[u,v] ≠ 2
                deg_v = push_label(deg_v,D[u,v])
                if deg_v === nothing
                    return nothing
                end
            end
        end end
        if deg_v == short_vec_to_deg([3,3,3])
            push!(deg3_vertices,v => simple_neighbors_v)
        elseif deg_v == short_vec_to_deg([3])
            deg1 = deg1 | SBitSet{1}(v)
            deg1_neigh = deg1_neigh | simple_neighbors_v
            push!(deg1_vertices,v => simple_neighbors_v)
        end
        # TODO add early exits here
        push!(deg_seqs, deg_v)
        
    end    
    ds = deg_seqs


    
    @debug "deg_seq is $ds versus"
    @debug "           $deg_seq_f4"
    @debug "end of build_deg_seq…"

    center_data::Union{Nothing,             #no data because we don't care
                       Tuple{Int,           #the center
                             SBitSet{1},    #its neighbors
                             SBitSet{1},    #the deg1 vertices
                             SBitSet{1}}    #their neighbors
                    } = nothing
    center::Union{Nothing,Int} = ( length(collect(keys(deg3_vertices))) > 0 ? collect(keys(deg3_vertices))[1] : nothing )
    center_neighbors::SBitSet{1} = ( center === nothing ? SBitSet{1}() : deg3_vertices[center] )
    extremities::SBitSet{1} = deg1 
    extremities_neighbors::SBitSet{1} = ( isempty(extremities) ? SBitSet{1}() : deg1_neigh ) 
    
    return (ds, deg1_vertices, deg3_vertices, center, center_neighbors, extremities, extremities_neighbors)::Tuple{DegSeq,Dict{Int,SBitSet{1}},Dict{Int,SBitSet{1}},Union{Nothing,Int},SBitSet{1},SBitSet{1},SBitSet{1}}

end

@inline function connected_non_sporadic_diagram_type(VS::SBitSet{1},D::Array{Int,2},deg_seq_and_assoc::Tuple{DegSeq,Dict{Int,SBitSet{1}},Dict{Int,SBitSet{1}},Union{Nothing,Int},SBitSet{1},SBitSet{1},SBitSet{1}})
    
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

@inline function connected_sporadic_diagram_type(VS::SBitSet{1},D::Array{Int,2},deg_seq_and_assoc::Tuple{DegSeq,Dict{Int,SBitSet{1}},Dict{Int,SBitSet{1}},Union{Nothing,Int},SBitSet{1},SBitSet{1},SBitSet{1}})

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
    elseif length(ds) == 2  &&
        big_label(ds.content[1]) ≠ nothing &&
        big_label(ds.content[1]) ≥ 7 &&
        ds == deg_seq_i(big_label(ds.content[1]))
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

function try_extend(VS::SBitSet{1},S::InducedSubDiagram,D::Array{Int,2},v::Int)
   
   
    if SBitSet{1}([1,2,3,4,5,6,8,9,10,11,12,13,14,15,16]) == VS && v == 17 
        println("does $VS + $v define a diagram?")
    end

    # special case, no component in S, so we just return the singleton
    if length(S.connected_components) == 0
        
        if SBitSet{1}([1,2,3,4,5,6,8,9,10,11,12,13,14,15,16]) == VS && v == 17 
            println("dead")
        end
        joined = the_singleton(v)::ConnectedInducedSubDiagram
        only_joined = Vector{ConnectedInducedSubDiagram}([joined])
        return InducedSubDiagram(only_joined)
    end
    
    # joined should/will be of type ConnectedInducedSubDiagram
    joined = nothing # Here is the result
    joined_vertices::SBitSet{1} = SBitSet{1}(v)

    components = copy(S.connected_components)
    
    freedom = 4


    total_size::Int = 1
    only_sporadic::Bool = false
    popped = false 
    idx = 1
    while idx ≤ length(components)
        c = components[idx]
            
        popped = false
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
                
                joined_vertices = joined_vertices | c.vertices
                if popped == false 
                    popat!(components,idx)
                    popped = true
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
        if popped == false
            idx+=1
        end

    end
    
    joined = connected_diagram_type(joined_vertices,D;only_sporadic=only_sporadic)
    if SBitSet{1}([1,2,3,4,5,6,8,9,10,11,12,13,14,15,16]) == VS && v == 17 
        println("joined is $joined")
    end
    
    
    if joined === nothing
        return nothing
    else
        push!(components,joined)
    if SBitSet{1}([1,2,3,4,5,6,8,9,10,11,12,13,14,15,16]) == VS && v == 17 
        println("conn_comps is $components")
    end
        return InducedSubDiagram(components)
    end

end



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
    
    function extend_one(sup_sub::Tuple{SBitSet{1},InducedSubDiagram})
        res = try_extend(sup_sub[1],sup_sub[2],D,new_vertex)
        if res === nothing
            return nothing
        else
            return (sup_sub[1]|SBitSet{1}(new_vertex),res)
        end
    end

    new_subs::Vector{Vector{Tuple{SBitSet{1},InducedSubDiagram}}} = map(x -> filter(y->!isnothing(y),extend_one.(x)),old_subs[1:min(max_card+1,end)])

    old_aligned::Vector{Vector{Tuple{SBitSet{1},InducedSubDiagram}}} = vcat(old_subs, [[]])
    new_aligned::Vector{Vector{Tuple{SBitSet{1},InducedSubDiagram}}} = vcat([[]],new_subs, [[] for i in 1:length(old_aligned) - length(new_subs)]) 
    

    old_and_new = [vcat(old_aligned[i],new_aligned[i]) for i in eachindex(old_aligned)]

    return DiagramAndSubs(D,old_and_new)
    
end


function build_diagram_and_subs(M::Array{Int,2};max_card::Union{Nothing,Int}=nothing)
   
    empty!(memoize_cache(connected_diagram_type))

    n = size(M,1)
    @assert size(M) == (n,n) "M must be square"
    @assert M == M' "M must be symmetric"
    @assert all(l ≥ 0 for l in M) "M must have non-negative entries"

    subs = Vector{Vector{Tuple{SBitSet{1},InducedSubDiagram}}}([[]])
    push!(subs[1],(SBitSet{1}(), the_empty_isd()))

    das = DiagramAndSubs(reshape([],0,0),subs)
    for i in 1:n
        das = extend(das,M[i,1:i-1];max_card=max_card)

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
            #dump_das(das;range=nothing)
            return (is_compact(das, rank), is_finite_volume(das, rank))
        end
    end


end



function check_all_graphs(sub_directory="")
    

    for (root, dirs, files) in walkdir("../graphs/"*sub_directory)
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
        for path in ["../graphs/13-mcl11.coxiter"] 
                println("path: $path")
                println(is_compact_respectively_finvol(path))
        end
        for (root, dirs, files) in walkdir("../graphs/simplices")
            for path in joinpath.(root, files)
                if endswith(path,".coxiter")
                    println("path: $path")
                    println(is_compact_respectively_finvol(path))
                end
            end
        end
    end

end
