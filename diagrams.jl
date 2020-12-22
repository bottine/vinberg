using MLStyle
using LightGraphs
using GraphPlot
using Colors
using Multisets

include("diagram_matrices.jl")

# Referring to https://en.wikipedia.org/wiki/Coxeter%E2%80%93Dynkin_diagram#Application_with_uniform_polytopes
# These are the different isomorphism types of subdiagrams of spherical/affine Coxeter-Dynkin diagrams
# 
# WARNING: for non sporadic diagrams, the number associated is the number of vertices, and not the rank
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

is_spherical(dt::DiagramType) = dt ∈ Set([DT_a,DT_b,DT_d,DT_e6,DT_e7,DT_e8,DT_f4,DT_g2,DT_h2,DT_h3,DT_h4,DT_in])
is_affine(dt::DiagramType) = dt ∈ Set([DT_A,DT_B,DT_C,DT_D,DT_F4,DT_G2,DT_Iinfty])

struct ConnectedInducedSubDiagram
    vertices::Array{Int,1}
    type::DiagramType
    ext_points::Set{Int}
    sec_points::Set{Int}
end

function CISD(vertices,type,ext_points=Set(),sec_points=Set())
    ConnectedInducedSubDiagram(vertices,type,ext_points,sec_points)
end


card(c::ConnectedInducedSubDiagram) = length(c.vertices)

struct InducedSubDiagram
    connected_components::Set{ConnectedInducedSubDiagram}
end 

function the_empty_isd() 
    return InducedSubDiagram(Set())
end


struct DiagramAndSubs 
    D::Array{Int,(2)}
    subs::Dict{BitSet,(InducedSubDiagram)}
end


function print_DAS(das::DiagramAndSubs)
    
    println("The matrix is given by:")
    display(das.D)
    println()
    println("and the subdiagrams are:")
    for (sub_support, sub_diagram) in das.subs
        println("$(collect(sub_support)):")
        for component in sub_diagram.connected_components
            println("    $(component.type) : $(component.vertices)")
        end
    end
    println("***********************")

end

function is_finite_volume(path::String)
    
    # Get file content in s as in https://en.wikibooks.org/wiki/Introducing_Julia/Working_with_text_files
    s = open(path) do file
        read(file, String)
    end

    (D, rank) = gug_coxiter_to_matrix(s)

    if D === nothing || rank === nothing
        println("Error reading file probably")
    else
        return is_finite_volume(build_diagram_and_subs(D),rank)
    end


end

function is_finite_volume(dag,dim)
    # Has spherical/parabolic diagram of rank n-1
    n = dim-1
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
                return false
            end
        end
    
    end

    return has_affine_sub_of_rank_n_minus_1 || has_spherical_sub_of_rank_n

end

is_empty(S1::ConnectedInducedSubDiagram) = length(S1.vertices) == 0
the_singleton(v::Integer) = CISD([v],DT_a)

function build_deg_seq_and_associated_data(VS::BitSet,D)

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
    println("ds is $ds")

    center = ( length(collect(keys(deg3_vertices))) > 0 ? collect(keys(deg3_vertices))[1] : nothing )
    center_neighbors = ( center === nothing ? Set() : deg3_vertices[center] )
    extremities = Set(keys(deg1_vertices))
    extremities_neighbors = ( isempty(extremities) ? Set() : ∪(values(deg1_vertices)...) )
    
    return (ds, deg1_vertices, deg3_vertices, center, center_neighbors, extremities, extremities_neighbors)

end

function small_connected_non_sporadic_diagram_type(VS::BitSet,D)
    
    @assert true "The diagram here is assumed connected. maybe this deserves a check"
    
    (ds, deg1_vertices, deg3_vertices, center, center_neighbors, extremities, extremities_neighbors) = build_deg_seq_and_associated_data(VS,D)
    

    vertices = collect(VS)
    n = length(vertices)

    if false
        @assert false "+++"
    elseif ds == deg_seq_a1
        return CISD(vertices,DT_a)
    elseif n≥2 && ds == deg_seq_a(n)
        return CISD(vertices,DT_a, Set(keys(deg1_vertices)))
    
    elseif ds == deg_seq_b2
        return CISD(vertices,DT_b)
    elseif ds == deg_seq_b3
        return CISD(vertices,DT_b)
    elseif n≥4 && ds == deg_seq_b(n)
        return CISD(vertices,DT_b, Set(keys(deg1_vertices)))
    
    elseif  n≥4 && ds == deg_seq_d(n)  && length(extremities ∩ center_neighbors) ≥ 2
        return CISD(vertices,DT_d,setdiff(extremities,center_neighbors))

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

function small_connected_sporadic_diagram_type(VS::BitSet,D)

    @assert true "The diagram here is assumed connected. maybe this deserves a check"


    (ds, deg1_vertices, deg3_vertices, center, center_neighbors, extremities, extremities_neighbors) = build_deg_seq_and_associated_data(VS,D)


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
    elseif ds == deg_seq_i(ds[1][1])
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

function try_extend(VS::BitSet,S::InducedSubDiagram,D,v::Int)
    components = S.connected_components
    
    # joined should/will be of type ConnectedInducedSubDiagram
    joined = nothing # Here is the result
    
    # special case, no component in S, so we just return the singleton
    if length(components) == 0
        joined = the_singleton(v)
        return InducedSubDiagram(Set([joined]))
    end


    # list of components, corresponding adjacencent vertices 
    neighbors_v = Set([u for u in eachindex(D[v,:]) if u≠v && u in VS && D[v,u]≠2])
    neighbors_v_labels = [D[u,v] for u in neighbors_v] 
    non_neighboring_components = []
    neighboring_components_data = [] 
    for c in components
        neighbors_in_c = [u for u in c.vertices if u in neighbors_v]
        if length(neighbors_in_c) > 0
            push!(neighboring_components_data,(c, [D[u,v] for u in neighbors_in_c]))
        else
            push!(non_neighboring_components,c)
        end
    end
   
    

    sort!(neighboring_components_data, by=(x->(x[1].type,card(x[1])))) # I think since the tuples have c.type as first entry and card(c) as second, the ordering is automatically the right one
    # TODO check that 
    ncd = neighboring_components_data
    @assert all(ncd[i][1].type ≤ ncd[i+1][1].type for i in 1:length(ncd)-1) "ordered according to type first"
    @assert all( ( ncd[i][1].type ≤ ncd[i+1][1].type ? card(ncd[i][1])≤card(ncd[i+1][1]) : true ) for i in 1:length(ncd)-1) "after type, ordered according to cardinality"
    

    neighboring_components_vertices = []
    neighboring_components_vertices_labels = []
    nv = neighbors_v 
    nc = [d[1] for d in ncd]            # neighboring_components 
    nct = [d[1].type for d in ncd]           # neighboring_components_types
    ncs  = [card(d[1]) for d in ncd]          # neighboring_components_sizes
    ncvl = [d[2] for d in ncd]          # neighboring_components_vertices_labels
    ncep = [d[1].ext_points for d in ncd]

    #=
    function swap(col::Set{T},from::T,to::T)
        @assert from ∈ col
        @assert to ≠ from
        return push(remove(col,from,to))
    end


    # This block covers only the non sporadic diagrams
    if false
        @assert false "For alignment's sake"

    elseif length(nv) == 0
        joined = the_singleton(v)
    
    elseif length(nv) == 1
         
        println("Hello length(nv) == 1")
        # v has only one neighbor, so extends exactly one (non sporadic) diagram with an edge
        # Since the affine diagrams can't be extended, the neighboring diagram of v is necessarily of type a,b,d.
        # The way those can be extended are:
        #
        #
        #     v ---- * ---- (the rest)          or         v ==== * ---- (the rest)           or            v ---- * ---- (the rest)
        #                                                                                                          |
        #                                                                                                          |
        #                                                                                                          *
        
        if false
            @assert false "For alignment's sake"
        
        # These correspond to "prolonging" a spherical (non-sporadic) diagram (i.e. v----*), thus resulting in a diagram of same type
        elseif nct[1] == DT_a && ncvl[1] == 3 && nv[1] ∈ ncep[1] 
            joined = CISD(push(nc.v,v),DT_a,ext_points=swap(ncep[1],nv[1]),v)
        elseif nct[1] == DT_b && ncvl[1] == 3 && nv[1] ∈ ncep[1]
            joined = CISD(push(nc.v,v),DT_b,ext_points=swap(ncep[1],nv[1],v))
        elseif nct[1] == DT_d && ncvl[1] == 3 && nv[1] ∈ ncep[1]
            joined = CISD(push(nc.v,v),DT_d,ext_points=swap(ncep[1],nv[1]),v)
        
        # These correspond to adding a  v ==== * edge to a spherical (non-sporadic) diagram, thus resulting in a b,B, or C diagram
        elseif nct[1] == DT_a && ncvl[1] == 4 && nv[1] ∈ ncep[1] 
            joined = CISD(push(nc.v,v),DT_b,ext_points=swap(ncep[1],nv[1],v)) 
        elseif nct[1] == DT_b && ncvl[1] == 4 && nv[1] ∈ ncep[1]
            joined = CISD(push(nc.v,v),DT_C)
        elseif nct[1] == DT_d && ncvl[1] == 4 && nv[1] ∈ ncep[1]
            joined = CISD(push(nc.v,v),DT_B)
        
        # These correspond to adding a  "third" edge to just before the end of a spherical (non-sporadic) diagram, thus resulting in a d,B or D diagram
        elseif nct[1] == DT_a && ncvl[1] == 3 && length(nc[1].left_end) == 2 && nv[1] == nc[1].left_end[2]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_d,[v,nc[1].left_end[1]],nc[1].right_end)
        elseif nct[1] == DT_a && ncvl[1] == 3 && length(nc[1].right_end) == 2 && nv[1] == nc[1].right_end[2]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_d,[v,nc[1].right_end[1]],nc[1].left_end) # By convention, the right end is the special one, hence the swap
        elseif nct[1] == DT_b && ncvl[1] == 3 && length(nc[1].left_end) == 2 && nv[1] == nc[1].left_end[2]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_B,[v,nc[1].left_end[1]],nc[1].right_end)
        elseif nct[1] == DT_d && ncvl[1] == 3 && length(nc[1].left_end) == 2 && nv[1] == nc[1].left_end[2]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_D,[v,nc[1].left_end[1]],nc[1].right_end)
        end

    elseif length(nv) == 2

       # v has 2 neighbors, either twice the same connected component, in which case it's of type a, and adding v makes it type A
       #                    or two connected components…
        
        println("Hello length(nv) == 2")
       
        if false
            @assert false "For alignment's sake"
        
        # type a to type A
        elseif length(nc) == 1 && nct[1] == DT_a && ncvl == [3,3] && nv == [nc[1].left_end[1],nc[1].right_end[1]]
            joined = ConnectedInducedSubDiagram(push(nc.vertices,v),DT_A,[])
        elseif length(nc) == 1 && nct[1] == DT_a && ncvl == [3,3] && nv == [nc[1].right_end[1],nc[1].left_end[1]]
            joined = ConnectedInducedSubDiagram(push(nc.vertices,v),DT_A,[])

        # v connects two components of type a each
        elseif length(nc) == 2 && nct == [DT_a,DT_a] && ncvl == [3,3] && nv == [nc[1].left_end[1],nc[2].right_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_a, nc[1].right_end, nc[2].left_end)
        elseif length(nc) == 2 && nct == [DT_a,DT_a] && ncvl == [3,3] && nv == [nc[1].right_end[1],nc[2].left_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_a, nc[1].left_end, nc[2].right_end)
        elseif length(nc) == 2 && nct == [DT_a,DT_a] && ncvl == [3,3] && nv == [nc[1].right_end[1],nc[2].right_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_a, nc[1].left_end, nc[2].left_end)
        elseif length(nc) == 2 && nct == [DT_a,DT_a] && ncvl == [3,3] && nv == [nc[1].left_end[1],nc[2].left_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_a, nc[1].right_end, nc[2].right_end)
           
        
        # v connects DT_a to DT_b
        elseif length(nc) == 2 && nct == [DT_a,DT_b] && ncvl == [3,3] && nv == [nc[1].left_end[1],nc[2].left_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_b, nc[1].right_end, nc[2].right_end)
        elseif length(nc) == 2 && nct == [DT_a,DT_b] && ncvl == [3,3] && nv == [nc[1].right_end[1],nc[2].left_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_b, nc[1].left_end, nc[2].right_end)
        
        # v connects DT_a to DT_d
        elseif length(nc) == 2 && nct == [DT_a,DT_d] && ncvl == [3,3] && nv == [nc[1].left_end[1],nc[2].left_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_d, nc[1].right_end, nc[2].right_end)
        elseif length(nc) == 2 && nct == [DT_a,DT_d] && ncvl == [3,3] && nv == [nc[1].right_end[1],nc[2].left_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_d, nc[1].left_end, nc[2].right_end)
           
        # v connects DT_b to DT_b resulting in DT_C
        elseif length(nc) == 2 && nct == [DT_b,DT_b] && ncvl == [3,3] && nv == [nc[1].left_end[1],nc[2].left_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_C, nc[1].right_end, nc[2].right_end)
        
        # v connects DT_b to DT_d resulting in DT_B
        elseif length(nc) == 2 && nct == [DT_b,DT_d] && ncvl == [3,3] && nv == [nc[1].left_end[1],nc[2].left_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_C, nc[1].right_end, nc[2].right_end)
        
            # v connects DT_d to DT_d resulting in DT_D
        elseif length(nc) == 2 && nct == [DT_d,DT_d] && ncvl == [3,3] && nv == [nc[1].left_end[1],nc[2].left_end[1]]
            joined = ConnectedInducedSubDiagram([nc[1].vertices nc[2].vertices [v]],DT_D, nc[1].right_end, nc[2].right_end)


        end


    elseif length(nv) == 3

        # v has 3 neighbors, hence 3 neighboring components too
        # Necessarily v yields a forky end, and thus corresponding to extensions of type a --> d, b --> B or d --> D
        
        println("Hello length(nv) == 3")

        if false
            @assert false "For alignment's sake"
        elseif nct == [DT_a,DT_a,DT_a] && ncvl == [3,3,3] && nv[1] == nc[1].left_end[1] && ncs == [ncs[1],1,1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_d,nv[2:3],nc[1].right_end)
        elseif nct == [DT_a,DT_a,DT_a] && ncvl == [3,3,3] && nv[1] == nc[1].right_end[1] && ncs == [ncs[1],1,1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_d,nc[1].left_end,nv[2:3])
        elseif nct == [DT_b,DT_a,DT_a] && ncvl == [3,3,3] && nv[1] == nc[1].left_end[1] && ncs == [ncs[1],1,1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_B,nv[2:3],nc[1].right_end)
        elseif nct == [DT_d,DT_a,DT_a] && ncvl == [3,3,3] && nv[1] == nc[1].left_end[1] && ncs == [ncs[1],1,1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_D,nv[2:3],nc[1].right_end)
        end

    end

    
    =#
    
    vertices = BitSet()
    for c in nc
        vertices = vertices ∪ BitSet(c.vertices)
    end
    push!(vertices,v)
    if joined === nothing && sum(ncs) + 1 ≤ 9 
        joined = small_connected_sporadic_diagram_type(vertices,D) 
    end
    if joined === nothing
        joined = small_connected_non_sporadic_diagram_type(vertices,D)
    end
    

    if joined === nothing
        return nothing
    else
        new_components = non_neighboring_components; 
        push!(new_components,joined)
        return InducedSubDiagram(Set(new_components))
    end

end


is_affine(isd::InducedSubDiagram) = all(is_affine(c.type) for c in isd.connected_components)
is_spherical(isd::InducedSubDiagram) = all(is_spherical(c.type) for c in isd.connected_components)

function extend(DAS::DiagramAndSubs, v::Array{Int,1})

    D = DAS.D
    subs = DAS.subs

    n = length(v) 
    @assert size(D) == (n,n)

    # Extend D with v
    D = [D v]
    D = [D; [v;1]']
   
    new_vertex = n+1

    new_subs::Dict{BitSet,InducedSubDiagram} = Dict()
    for (V,S) in subs
        S_and_v = try_extend(V,S,D,new_vertex)
        if S_and_v ≠ nothing 
            push!(new_subs,(V∪BitSet([new_vertex])) => S_and_v)
        end
    end
    return DiagramAndSubs(D,merge(subs,new_subs))
    
end


function build_diagram_and_subs(M)
    n = size(M,1)
    @assert size(M) == (n,n)
    # @assert issymetric(M)
    
    lol = Dict()
    push!(lol,BitSet() => the_empty_isd())

    DAS = DiagramAndSubs(reshape([],0,0),lol)
    for i in 1:n
        DAS = extend(DAS,M[i,1:i-1]) 
    end
    return DAS
end


is_finite_volume(build_diagram_and_subs([1 2 3;
                        2 1 3;
                        3 2 1]),2)
