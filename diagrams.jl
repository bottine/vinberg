using MLStyle
using LightGraphs
using GraphPlot
using Colors
using Multisets


include("diagram_matrices.jl")


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


function is_spherical(dt::DiagramType)::Bool
    dt ∈ [DT_a,DT_b,DT_d,DT_e6,DT_e7,DT_e8,DT_f4,DT_g2,DT_h2,DT_h3,DT_h4,DT_in]
end

function is_affine(dt::DiagramType)::Bool
    dt ∈ [DT_A,DT_B,DT_C,DT_D,DT_E6,DT_E7,DT_E8,DT_F4,DT_G2,DT_Iinfty]
end

function is_sporadic(dt::DiagramType)::Bool
    dt ∈ [DT_e6,DT_e7,DT_e8,DT_f4,DT_g2,DT_h2,DT_h3,DT_h4,DT_in,DT_E6,DT_E7,DT_E8,DT_F4,DT_G2,DT_Iinfty]
end





MS = Multiset

# Degree Sequence: Each vertex has an associated "multi-degree" which is a multiset containing the labels of edges incident to the vertex
# A DegSeq is the multiset of the mult-degrees associated to the vertices of a diagram
DegSeq = MS{MS{Int}}

# Given an array of arrays of ints, returns the associated degree sequence
function deg_seq(vv::Array{Array{Int,1},1})
    return MS(MS.(vv))::DegSeq
end

# The degree sequences corresponding to each irreducible diagram types follow:

const deg_seq_a1 = deg_seq(Array{Array{Int,1},1}([[]]))
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
    deg_seq([[3,3,4]]) + 2*deg_seq([[3]]) + deg_seq([[4]])
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
    vertices::Array{Int,1}
    type::DiagramType
end

function CISD(vertices,type)
    ConnectedInducedSubDiagram(vertices,type)
end


card(c::ConnectedInducedSubDiagram) = length(c.vertices)

is_empty(c::ConnectedInducedSubDiagram) = length(c.vertices) == 0

the_singleton(v::Int) = CISD([v],DT_a)


# An arbitrary induced subdiagram
# Stored as a collection of its irreducible components, plus whether it is affine or spherical
struct InducedSubDiagram
    connected_components::Set{ConnectedInducedSubDiagram}
    is_affine::Bool
    is_spherical::Bool
end

function InducedSubDiagram(connected_components::Set{ConnectedInducedSubDiagram})
    this_is_affine = all(is_affine(c.type) for c in connected_components) 
    this_is_spherical = all(is_spherical(c.type) for c in connected_components) 
    
    @assert ! (length(connected_components) > 0 && this_is_affine && this_is_spherical) "A diagram can't both be spherical and affine."
    
    return InducedSubDiagram(connected_components, this_is_affine, this_is_spherical) 
end

function the_empty_isd() 
    return InducedSubDiagram(Set{ConnectedInducedSubDiagram}())
end


struct DiagramAndSubs 
    D::Array{Int,(2)}
    subs::Dict{Set{Int},InducedSubDiagram}
end


function print_das(das::DiagramAndSubs)
    
    println("The matrix is given by:")
    display(das.D)
    println()
    println("and the subdiagrams are:")
    for (sub_support, sub_diagram) in das.subs
        println("$(collect(sub_support)) [affine=$(sub_diagram.is_affine),spherical=$(sub_diagram.is_spherical)]:")
        for component in sub_diagram.connected_components
            println("    $(component.type) : $(component.vertices)")
        end
    end
    println("***********************")

end



function is_finite_volume(dag::DiagramAndSubs,n::Int)
    # Has spherical/parabolic diagram of rank n-1
    has_spherical_sub_of_rank_n = false
    has_affine_sub_of_rank_n_minus_1 = false

    subs = dag.subs

    for (support, subdiagram) in subs
        
        if length(support) == n && is_spherical(subdiagram)
            has_spherical_sub_of_rank_n = true
        end
        if length(support) == n && is_affine(subdiagram)
            has_affine_sub_of_rank_n_minus_1 = true
        end

        if length(support) == n-1 && is_spherical(subdiagram)
            extensions = [
                ext_support for (ext_support, ext_subdiagram) in subs if 
                support < ext_support && 
                (
                    ( length(ext_support) == n && is_spherical(ext_subdiagram) ) || 
                    ( length(ext_support) == n && is_affine(ext_subdiagram) ) 

                ) 
            ]
            if length(extensions) ≠ 2
                @debug "the subdiagram of support $support has $(length(extensions)) affine/spherical extensions"
                return false
            end
        end
    
    end
    
    @debug "has_affine_sub_of_rank_n_minus_1 = $has_affine_sub_of_rank_n_minus_1"
    @debug "has_spherical_sub_of_rank_n = $has_spherical_sub_of_rank_n"

    return has_affine_sub_of_rank_n_minus_1 || has_spherical_sub_of_rank_n

end


function build_deg_seq_and_associated_data(VS::Set{Int},D::Array{Int,2})

    @assert true "The diagram here is assumed connected. maybe this deserves a check"

    @debug "build_deg_seq_and_associated_data(…)"
    @debug "VS is $(collect(VS))"

    VSC = length(VS)
    
    deg1_vertices = Dict{Int,Set{Int}}()
    deg3_vertices = Dict{Int,Set{Int}}()
    
    deg_seqs = MS{MS{Int}}()
    for v in VS

        @debug "looking at $v"

        deg_seq_v = MS{Int}()
        simple_neighbors_v = Set{Int}()
        for u in VS if u ≠ v
            if D[u,v] == 3
                push!(simple_neighbors_v,u)
            end
            if D[u,v] ≠ 2 
                push!(deg_seq_v,D[u,v])
            end
        end end
        if deg_seq_v == MS([3,3,3])
            push!(deg3_vertices,v => simple_neighbors_v)
        elseif deg_seq_v == MS([3])
            push!(deg1_vertices,v => simple_neighbors_v)
        end
        push!(deg_seqs, deg_seq_v)
        
        @debug "$v has deg_seq = $deg_seq_v"
    end    
    ds = deg_seqs

    center::Union{Nothing,Int} = ( length(collect(keys(deg3_vertices))) > 0 ? collect(keys(deg3_vertices))[1] : nothing )
    center_neighbors::Set{Int} = ( center === nothing ? Set() : deg3_vertices[center] )
    extremities::Set{Int} = Set(keys(deg1_vertices))
    extremities_neighbors::Set{Int} = ( isempty(extremities) ? Set() : ∪(values(deg1_vertices)...) )
    
    return (ds, deg1_vertices, deg3_vertices, center, center_neighbors, extremities, extremities_neighbors)::Tuple{DegSeq,Dict{Int,Set{Int}},Dict{Int,Set{Int}},Union{Nothing,Int},Set{Int},Set{Int},Set{Int}}

end

function small_connected_non_sporadic_diagram_type(VS::Set{Int},D::Array{Int,2},deg_seq_and_assoc::Tuple{DegSeq,Dict{Int,Set{Int}},Dict{Int,Set{Int}},Union{Nothing,Int},Set{Int},Set{Int},Set{Int}})
    
    @assert true "The diagram here is assumed connected. maybe this deserves a check"
    
    (ds, deg1_vertices, deg3_vertices, center, center_neighbors, extremities, extremities_neighbors) = deg_seq_and_assoc # build_deg_seq_and_associated_data(VS,D) 
    

    vertices = collect(VS)
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

function small_connected_sporadic_diagram_type(VS::Set{Int},D::Array{Int,2},deg_seq_and_assoc::Tuple{DegSeq,Dict{Int,Set{Int}},Dict{Int,Set{Int}},Union{Nothing,Int},Set{Int},Set{Int},Set{Int}})

    @assert true "The diagram here is assumed connected. maybe this deserves a check"


    (ds, deg1_vertices, deg3_vertices, center, center_neighbors, extremities, extremities_neighbors) = deg_seq_and_assoc # build_deg_seq_and_associated_data(VS,D) 


    VS = collect(VS)
    
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
    elseif length(ds) ≥ 1 && length(first(ds)) ≥ 1 && first(first(ds)) ≥ 6 && ds == deg_seq_i(first(first(ds)))
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

function try_extend(VS::Set{Int},S::InducedSubDiagram,D::Array{Int,2},v::Int)
    components = S.connected_components
    
    # joined should/will be of type ConnectedInducedSubDiagram
    joined = nothing # Here is the result
    
    # special case, no component in S, so we just return the singleton
    if length(components) == 0
        joined = the_singleton(v)
        return InducedSubDiagram(Set([joined]))
    end
    
    freedom = 4

    neighboring_components::Set{ConnectedInducedSubDiagram} = Set()
    neighboring_components_size::Set{Int} = Set()
    non_neighboring_components::Set{ConnectedInducedSubDiagram} = Set()

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
                freedom -= (min(D[u,v],5)-2)
                push!(neighboring_components,c)
                total_size += card(c)
                if is_sporadic(c.type)
                    only_sporadic = true
                end
                if D[u,v] > 4
                    only_sporadic = true
                end
            end
        end
        if c ∉ neighboring_components
            push!(non_neighboring_components,c)
        end

    end

    
    vertices::Set{Int} = ∪(Set(),[Set(c.vertices)::Set{Int} for c in neighboring_components]...)
    push!(vertices,v)

    deg_seq_and_assoc = build_deg_seq_and_associated_data(vertices,D)

    if joined === nothing && total_size ≤ 9 
        joined = small_connected_sporadic_diagram_type(vertices,D,deg_seq_and_assoc) 
    end
    if joined === nothing && !only_sporadic
        joined = small_connected_non_sporadic_diagram_type(vertices,D,deg_seq_and_assoc)
    end
    

    if joined === nothing
        return nothing
    else
        new_components = non_neighboring_components 
        push!(new_components,joined)
        return InducedSubDiagram(Set{ConnectedInducedSubDiagram}(new_components))
    end

end


is_affine(isd::InducedSubDiagram) = all(is_affine(c.type) for c in isd.connected_components)
is_spherical(isd::InducedSubDiagram) = all(is_spherical(c.type) for c in isd.connected_components)

function extend(das::DiagramAndSubs, v::Array{Int,1}; max_card::Union{Nothing,Int}=nothing)
  
    D = das.D
    subs = das.subs

    n = length(v) 
    @assert size(D) == (n,n)
    
    # Extend D with v
    D = [D v]
    D = [D; [v;1]']
   
    new_vertex = n+1

    new_subs::Dict{Set{Int},InducedSubDiagram} = Dict{Set{Int},InducedSubDiagram}()
    for (V,S) in subs
        if !isnothing(max_card) && length(V) ≤ max_card - 1
            S_and_v = try_extend(V,S,D,new_vertex)
            if S_and_v ≠ nothing 
                push!(new_subs,(V∪Set([new_vertex])) => S_and_v)
            end
        else
            # means the subdiagram is too big and we don't care (because it won't matter for finite volume/cocompactness verification) 
        end
    end
    return DiagramAndSubs(D,merge(subs,new_subs))
    
end


function build_diagram_and_subs(M::Array{Int,2};max_card::Union{Nothing,Int}=nothing)
    
    n = size(M,1)
    @assert size(M) == (n,n) "M must be square"
    @assert M == M' "M must be symmetric"
    @assert all(l ≥ 0 for l in M) "M must have non-negative entries"

    subs = Dict{Set{Int},InducedSubDiagram}()
    push!(subs,Set{Int}() => the_empty_isd())

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


function is_finite_volume(path::String)
    
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
            return is_finite_volume(build_diagram_and_subs(D;max_card=rank),rank)
        end
    end


end


function check_all_graphs()
    

    for (root, dirs, files) in walkdir("./graphs/")
        for path in joinpath.(root, files)
            if endswith(path,".coxiter")
                println("path: $path")
                println(is_finite_volume(path))
            end
        end
    end

end


function check_some_graphs()

    @time begin
        for path in ["graphs/13-mcl11.coxiter"] 
                println("path: $path")
                println(is_finite_volume(path))
        end
        for (root, dirs, files) in walkdir("./graphs/simplices")
            for path in joinpath.(root, files)
                if endswith(path,".coxiter")
                    println("path: $path")
                    println(is_finite_volume(path))
                end
            end
        end
    end

end
