using Memoize
using StaticArrays
using Debugger

include("sbitset.jl")

# TODO
#
# * memoize connected_diagram_type by hand
# * make computation of build_deg_seq_and_associated_data more efficient by using previously constructed deg_seqs of the previous components?
# * clean up structure
# * back from bitsets to lists?
#
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
    vertices::SBitSet{4}
    type::DiagramType
end
function Base.:(==)(c1::ConnectedInducedSubDiagram,c2::ConnectedInducedSubDiagram)
    (c1.vertices == c2.vertices)
end

function CISD(vertices::SBitSet{4},type)
    ConnectedInducedSubDiagram(vertices,type)
end
function CISD(vertices::Array{Int,1},type)
    ConnectedInducedSubDiagram(SBitSet{4}(vertices),type)
end

function Base.copy(c::ConnectedInducedSubDiagram)
    return CISD(copy(c.vertices),c.type) 
end

card(c::ConnectedInducedSubDiagram) = length(c.vertices)

is_empty(c::ConnectedInducedSubDiagram) = isempty(c.vertices)

the_singleton(v::Int) = CISD(SBitSet{4}(v),DT_a)


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


function cisd_rank(cisd::ConnectedInducedSubDiagram)
    if is_affine(cisd.type)
        return length(cisd.vertices) - 1
    elseif is_spherical(cisd.type)
        return length(cisd.vertices)
    else
        @assert false "unreachable"
    end
end

function isd_rank(isd::InducedSubDiagram)
    
    #@assert is_spherical(isd) || is_affine(isd) "We only care about the rank for spherical/affine diagrams"
    
    return sum(vcat([0], [cisd_rank(c) for c in isd.connected_components])) # should be equal cardinality - connected_compnents

end

mutable struct DiagramAndSubs 
    D::Array{Int,(2)}
    d::Int # = dimension
    subs::Vector{Tuple{SBitSet{4},InducedSubDiagram}}
    spherical_subs_rank_d_minus_1::Vector{Tuple{SBitSet{4},InducedSubDiagram}}
    spherical_subs_rank_d::Vector{Tuple{SBitSet{4},InducedSubDiagram}}
    affine_subs_rank_d_minus_1::Vector{Tuple{SBitSet{4},InducedSubDiagram}}
end




function dump_das(das::DiagramAndSubs;range=nothing)
   
    dense_bitset_str(b::SBitSet{4}) = *("[",[string(i)*"," for i in b]...,"]")

    println("### Subdiagrams for Coxeter matrix:")
    display(das.D)
    println()
    println("###")
    for i in eachindex(das.subs)
        if range === nothing || i ∈ range 
            println("Cardinality $(i-1):")
            for (sub_support, sub_diagram) in das.subs[i]
                print("    $(dense_bitset_str(sub_support)) $(sub_diagram.is_affine ? "A" : " ") $(sub_diagram.is_spherical ? "S" : " ") r=$(isd_rank(sub_diagram)) is ")
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




function is_compact(das::DiagramAndSubs)
    
    # The dimension is the rank
    #

    dense_bitset_str(b::SBitSet{4}) = *("[",[string(i)*"," for i in b]...,"]")
   
    if isempty(das.spherical_subs_rank_d)
        return false
    end

    for (support, subdiagram) in das.spherical_subs_rank_d_minus_1 

        extensions = []
        num_extensions = 0
        for (sup,sub) in das.spherical_subs_rank_d 
            if support ⊆ sup
                num_extensions += 1
                push!(extensions,(sup,sub))
            end
        end
        if num_extensions ≠ 2
            #println("the subdiagram of support $support has $(length(extensions)) affine/spherical extensions")
            #for (sup,sub) in extensions
            #    println("    $(dense_bitset_str(sup)) $(sub.is_affine ? "A" : " ") $(sub.is_spherical ? "S" : " ") r=$(isd_rank(sub)) = ")
            #end
            @debug "the subdiagram of support $support has $(length(extensions)) affine/spherical extensions"
            return false
        end
    
    end
    

    return true

end

function is_finite_volume(das::DiagramAndSubs)

    
    dense_bitset_str(b::SBitSet{4}) = *("[",[string(i)*"," for i in b]...,"]")
    
    if isempty(das.spherical_subs_rank_d) && isempty(das.affine_subs_rank_d_minus_1)
        return false
    end

    for (support, subdiagram) in das.spherical_subs_rank_d_minus_1

        num_extensions = 0
        extensions = []
        for (sup,sub) in das.spherical_subs_rank_d
            if support ⊆ sup
                num_extensions += 1
                push!(extensions,(sup,sub))
            end
        end
        for (sup,sub) in das.affine_subs_rank_d_minus_1
            if support ⊆ sup
                num_extensions += 1
                push!(extensions,(sup,sub))
            end
        end
        if num_extensions ≠ 2
            #println("the subdiagram of support $support has $(length(extensions)) affine/spherical extensions")
            #for (sup,sub) in extensions
            #    println("    $(dense_bitset_str(sup)) $(sub.is_affine ? "A" : " ") $(sub.is_spherical ? "S" : " ") r=$(isd_rank(sub))")
            #end
            return false
        end

    
    end
    

    return true

end

function is_compact_finite_volume(das::DiagramAndSubs)

    
    compact = true
    fin_vol = true
    dense_bitset_str(b::SBitSet{4}) = *("[",[string(i)*"," for i in b]...,"]")
    
    if isempty(das.spherical_subs_rank_d)
        compact = false
        if isempty(das.affine_subs_rank_d_minus_1)
            fin_vol = false
            return (compact, fin_vol)
        end
    end

    for (support, subdiagram) in das.spherical_subs_rank_d_minus_1

        num_sph_extensions = 0
        num_aff_extensions = 0
        for (sup,sub) in das.spherical_subs_rank_d
            if support ⊆ sup
                num_sph_extensions += 1
            end
        end
        for (sup,sub) in das.affine_subs_rank_d_minus_1
            if support ⊆ sup
                num_aff_extensions += 1
            end
        end
        if num_sph_extensions ≠ 2
            compact = false
        end
        if num_sph_extensions + num_aff_extensions ≠ 2
            fin_vol = false
        end
        if (compact,fin_vol) == (false,false)
            return (compact,fin_vol)
        end
    
    end

    return (compact,fin_vol)

end



# dirty hack to memoize `connected_diagram_type` only on its first argument
Arg = Tuple{SBitSet{4},Array{Int,2},Bool}
Base.hash(a::Arg, h::UInt) = hash(a[1], hash(:Arg, h))
Base.isequal(a::Arg, b::Arg) = Base.isequal(hash(a), hash(b))


@memoize Dict function connected_diagram_type(VS::SBitSet{4},D::Array{Int,2}; only_sporadic::Bool=false)
    
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

@inline function build_deg_seq_and_associated_data(VS::SBitSet{4},D::Array{Int,2})

    
    @assert true "The diagram here is assumed connected. maybe this deserves a check"

    @debug "build_deg_seq_and_associated_data(…)"
    @debug "VS is $(collect(VS))"

    VSC = length(VS)
   
    deg1 = SBitSet{4}()
    deg1_neigh = SBitSet{4}()
    deg3 = SBitSet{4}()
    deg3_neigh = SBitSet{4}()

    deg_seqs::DegSeq = DegSeq(Vector{Deg}()) 
    for v in VS

        @debug "looking at $v"

        deg_v = empty_deg 
        simple_neighbors_v = SBitSet{4}()
        for u in VS if u ≠ v
            if D[u,v] == 3
                simple_neighbors_v = simple_neighbors_v | SBitSet{4}(u)
            end
            if D[u,v] ≠ 2
                deg_v = push_label(deg_v,D[u,v])
                if deg_v === nothing
                    return nothing
                end
            end
        end end
        if deg_v == short_vec_to_deg([3,3,3])
            deg3 = deg3 | SBitSet{4}(v)
            deg3_neigh = deg3_neigh | simple_neighbors_v
        elseif deg_v == short_vec_to_deg([3])
            deg1 = deg1 | SBitSet{4}(v)
            deg1_neigh = deg1_neigh | simple_neighbors_v
        end
        # TODO add early exits here
        push!(deg_seqs, deg_v)
        
    end    
    ds = deg_seqs
    

    
    @debug "deg_seq is $ds versus"
    @debug "           $deg_seq_f4"
    @debug "end of build_deg_seq…"


    center::SBitSet{4} = deg3
    center_neighbors::SBitSet{4} = deg3_neigh 
    extremities::SBitSet{4} = deg1 
    extremities_neighbors::SBitSet{4} = deg1_neigh 
    
    return (ds, deg1, deg1_neigh, deg3, deg3_neigh)

end

@inline function connected_non_sporadic_diagram_type(VS::SBitSet{4},D::Array{Int,2},deg_seq_and_assoc)
    
    @assert true "The diagram here is assumed connected. maybe this deserves a check"
    
    (ds, extremities, extremities_neighbors, center, center_neighbors) = deg_seq_and_assoc # build_deg_seq_and_associated_data(VS,D) 
    
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
    elseif n≥5 && ds == deg_seq_B(n) && length(extremities ∩ center_neighbors) == 2
        return CISD(vertices,DT_B)

    elseif ds == deg_seq_C3
        return CISD(vertices,DT_C)
    elseif n≥4 && ds == deg_seq_C(n)
        return CISD(vertices,DT_C)

    elseif ds == deg_seq_D5 
        return CISD(vertices,DT_D)
    elseif n≥6 && ds == deg_seq_D(n) && length(extremities ∩ center_neighbors) ≥ 4
        return CISD(vertices,DT_D)
    else
        return nothing

    end    
end

@inline function connected_sporadic_diagram_type(VS::SBitSet{4},D::Array{Int,2},deg_seq_and_assoc)

    @assert true "The diagram here is assumed connected. maybe this deserves a check"


    (ds, extremities, extremities_neighbors, center, center_neighbors) = deg_seq_and_assoc # build_deg_seq_and_associated_data(VS,D) 
    


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
    elseif ds == deg_seq_E6 && isempty(center_neighbors∩extremities)    
        return CISD(VS,DT_E6)
    elseif ds == deg_seq_E7 && length(center_neighbors∩extremities) == 1 && isempty(extremities_neighbors ∩ center_neighbors) 
        return CISD(VS,DT_E7)
    elseif ds == deg_seq_E8 && length(center_neighbors∩extremities) == 1 && length(extremities_neighbors ∩ center_neighbors) == 1 
        return CISD(VS,DT_E8)
    end
    
    return nothing

end

function try_extend(VS::SBitSet{4},S::InducedSubDiagram,D::Array{Int,2},v::Int)::Union{Nothing,InducedSubDiagram}
   
   

    # special case, no component in S, so we just return the singleton
    if isempty(S.connected_components)
        joined = the_singleton(v)::ConnectedInducedSubDiagram
        only_joined = Vector{ConnectedInducedSubDiagram}([joined])
        return InducedSubDiagram(only_joined)
    end
    
    # joined should/will be of type ConnectedInducedSubDiagram
    joined = nothing # Here is the result
    joined_vertices::SBitSet{4} = SBitSet{4}(v)

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
    
    if joined === nothing
        return nothing
    else
        push!(components,joined)

        return InducedSubDiagram(components)::InducedSubDiagram
    end

end



function extend!(das::DiagramAndSubs, v::Array{Int,1})

    n = length(v) 
    @assert size(das.D) == (n,n)
    
    # Extend D with v
    das.D = [das.D v]
    das.D = [das.D; [v;1]']
   
    new_vertex = n+1
    singleton_v = SBitSet{4}(new_vertex)

    #= Can't be more efficient it seems
    @inline function extend_one(supsub::Tuple{SBitSet{4},InducedSubDiagram})
        sup,sub = supsub
        res = try_extend(sup,sub,das.D,new_vertex)
        if res ≠ nothing
            sup_v = sup|singleton_v
            sub_v = res
            if is_spherical(sub_v) && length(sup_v) == das.d-1
                push!(das.spherical_subs_rank_d_minus_1,(sup_v,sub_v))
            elseif is_spherical(sub_v) && length(sup_v) == das.d
                push!(das.spherical_subs_rank_d,(sup_v,sub_v))
            elseif is_affine(sub_v) && length(sup_v) - length(sub_v.connected_components) == das.d-1
                push!(das.affine_subs_rank_d_minus_1,(sup_v,sub_v))
            end
            return (sup_v,sub_v)
        else
            return nothing
        end
    end
    new_subs::Vector{Tuple{SBitSet{4},InducedSubDiagram}} =filter(x -> !isnothing(x),extend_one.(das.subs))
    append!(das.subs,new_subs)
    sizehint!(das.subs,length(das.subs)*2)
    =#
    
    #new_subs::Vector{Tuple{SBitSet{4},InducedSubDiagram}} = sizehint!([],length(das.subs)) 
    #sizehint!(das.subs,length(das.subs)*2)
    len = length(das.subs)
    i = len+1
    resize!(das.subs,2*len)
    for j in 1:len
        (sup,sub) = das.subs[j]
        res = try_extend(sup,sub,das.D,new_vertex)
        if res ≠ nothing
            sup_v = sup|singleton_v
            sub_v = res
            #push!(new_subs,(sup_v,sub_v))
            if is_spherical(sub_v) && length(sup_v) == das.d-1
                push!(das.spherical_subs_rank_d_minus_1,(sup_v,sub_v))
            elseif is_spherical(sub_v) && length(sup_v) == das.d
                push!(das.spherical_subs_rank_d,(sup_v,sub_v))
            elseif is_affine(sub_v) && length(sup_v) - length(sub_v.connected_components) == das.d-1
                push!(das.affine_subs_rank_d_minus_1,(sup_v,sub_v))
            end
            das.subs[i] = (sup_v,sub_v)
            i += 1
        end
    end

    resize!(das.subs,i-1)
    #append!(das.subs,new_subs)
   
end


function build_diagram_and_subs(M::Array{Int,2},dimension::Int)
   
    empty!(memoize_cache(connected_diagram_type))

    n = size(M,1)
    @assert size(M) == (n,n) "M must be square"
    @assert M == M' "M must be symmetric"
    @assert all(l ≥ 0 for l in M) "M must have non-negative entries"

    subs = Vector{Tuple{SBitSet{4},InducedSubDiagram}}([])
    push!(subs,(SBitSet{4}(), the_empty_isd()))

    das = DiagramAndSubs(reshape([],0,0),dimension,subs,[],[],[])
    for i in 1:n
        #println(i)
        extend!(das,M[i,1:i-1])
    end
    return das
end



# ##########################
#
# Misc user-facing functions
#
# ##########################
#

function is_fin_vol(D::Array{Int,2},dimension::Int)
    
    das = build_diagram_and_subs(D,dimension)
    return is_finite_volume(das)

end

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
            das = build_diagram_and_subs(D,rank)
            #dump_das(das;range=nothing)
            compact = is_compact(das)
            fin_vol = is_finite_volume(das)
            return is_compact_finite_volume(das) 
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
