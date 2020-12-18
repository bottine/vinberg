using MLStyle
using LightGraphs
using GraphPlot
using Colors

# Referring to https://en.wikipedia.org/wiki/Coxeter%E2%80%93Dynkin_diagram#Application_with_uniform_polytopes
# These are the different isomorphism types of subdiagrams of spherical/affine Coxeter-Dynkin diagrams
@enum DiagramType begin
    # Lowercase == Spherical
    DT_a
    DT_b # == C 
    DT_d
    DT_A
    DT_B
    DT_C
    DT_D
    DT_spherical_sporadic
    DT_affine_sporadic
end

is_spherical(dt::DiagramType) = dt ∈ Set([DT_a,DT_b,DT_d,DT_spherical_sporadic])
is_affine(dt::DiagramType) = dt ∈ Set([DT_A,DT_B,DT_C,DT_D,DT_affine_sporadic])


struct ConnectedInducedSubDiagram
    vertices::Array{Int,(1)}
    type::DiagramType
    left_end::Array{Int,1}
    right_end::Array{Int,1}
end


card(c::ConnectedInducedSubDiagram) = length(c.vertices)

function are_pairwise_disconnected(VV,D) 
    for vv in VV, ww in VV if vv ≠ ww
        for v in vv, w in ww
            if v == w || D[v,w] ≠ 2
                return false
            end
        end
    end end
    return true
end


function is_path(vv,D,first=(x -> x == 3),last=(x -> x == 3),middle=(val,idx -> val==3))
    n = length(vv)
    for i in 1:n, j in i+1:n
        if i == 1 && j == 2
            if !first(D[vv[i],vv[j]])
                return false
            end
        elseif i==n-1 && j == n 
            if !last(D[vv[i],vv[j]])
                return false
            end
        elseif j == i+1
            if !middle(D[vv[i],vv[j],i])##D[vv[i],vv[j]] ≠ 3
                return false
            end
        elseif j > i+1
            if D[vv[i],vv[j]] ≠ 2
                return false
            end
        else
            @assert false "unreachable"
        end
    end
    return true
end

function reindex(vv,along...)
    reduce(vcat,[vv[a] for a in along])
end

function drop_slice(vv,from,to) 
    return reindex(vv,1:from-1,to+1:length(vv))
end

function is_subdiagram(c::ConnectedInducedSubDiagram,D,n)
    @assert size(D) == (n,n)
    @assert D == D'

    # canonical "orientation" so that we can always assume the left end is the "lower" one (arbitrary order)
    if c.left_end > c.right_end  
        c = ConnectedInducedSubDiagram(reverse(c.vertices),c.right_end,c.left_end)
    end

    l = c.left_end
    r = c.right_end

    vv = c.vertices
    card = length(vv)

    return begin
        if false
            # Welcome to the poor person's pattern matching
            @assert false "unreachable"
        elseif (l,r) == (NoEnd,NoEnd) && is_path(vv,D) 
            (true, DT_a)
        elseif (l,r) == (NoEnd,ASingleFour) && card ≥ 2 && is_path(vv,D,last=(x -> x == 4)) 
            (true, DT_b)
        elseif (l,r) == (NoEnd, ASimpleBranch) && card ≥ 4 && is_path(drop_slice(vv,card,card),D) && is_path(drop_slice(vv,card-1,card-1),D) && D[vv[end-1],vv[end]] == 2  
            (true, DT_d)
        elseif (l,r) == (NoEnd, ASingleInfty) && card == 2
            (true, DT_I)
        elseif (l,r) == (NoEnd, ASingleNgeq5) && card == 2 && D[vv[1],vv[2]] == 5
            (true, DT_h2)
        elseif (l,r) == (NoEnd, ASingleNgeq5) && card == 2 && D[vv[1],vv[2]] == 6
            (true, DT_g2)
        elseif (l,r) == (NoEnd, ASingleNgeq5) && card == 2 && D[vv[1],vv[2]] ≥ 7
            (true, DT_i)
        elseif (l,r) == (NoEnd, ASingleNgeq5) && card == 3 && D[vv[2],vv[3]] == 5 && is_path(vv, D, last=(x -> x == 5))
            (true, DT_h3)
        elseif (l,r) == (NoEnd, ASingleNgeq5) && card == 4 && D[vv[3],vv[4]] == 5 && is_path(vv, D, last=(x -> x == 5))
            (true, DT_h4)
        elseif (l,r) == (NoEnd, ASingleNgeq5) && card == 3 && D[vv[2],vv[3]] == 6 && is_path(vv, D, last=(x -> x == 6))
            (true, DT_G2)
        elseif (l,r) == (NoEnd, AFourAndThen) && card == 4 && is_path(vv[1:end-1],last=(x-> x == 4)) && is_path(vv[2:end],first=(x-> x == 4)) && are_pairwise_disconnected([vv[1:2],vv[4:4]]) 
            (true, DT_f4)
        elseif (l,r) == (NoEnd, AFourAndThen) && card == 5 && is_path(vv[1:end-1],last=(x-> x == 4)) && is_path(vv[end-1:end],first=(x-> x == 4)) && are_pairwise_disconnected([vv[1:3],vv[5:5]]) # TODO needs more
            (true, DT_F4)
        elseif (l,r) == (ASimpleBranch, ASingleFour) && card ≥ 5 && is_path(drop_slice(vv,2,2),D,last=(x -> x == 4)) && is_path(drop_slice(vv,1,1),D,last=(x -> x == 4)) && D[vv[1],vv[2]] == 2
            (true, DT_B)
        elseif (l,r) == (ASimpleBranch, ASingleFour) && card ≥ 7 && is_path(vv[2:end-1],D) && are_pairwise_disconnected([vv[3,end-2],vv[1:1],vv[2:2],vv[end-1:end-1],vv[end:end]])
            (true, DT_D)
        elseif (l,r) == (ASingleFour, ASingleFour) && card ≥ 5 && is_path(vv,D,first=(x->x==4),last=(x->x==4))
            (true, DT_C)
        elseif (l,r) == (ALoop, ALoop) && card ≥ 3 && is_path(vv[1:end-1]) && is_path(vcat(vv[2:end], [vv[1]])) 
            (true, DT_F4)
        elseif (l,r) == (NoEnd, TwoOneBranch) && card == 6 && is_path(drop_slice(vv,3,3),D) && is_path(vv[3:end],D) && is_path(reindex(vv,1:2,4:-1:3))
            (true, DT_e6)
        elseif (l,r) == (NoEnd, TwoOneBranch) && card == 7 && is_path(drop_slice(vv,3,3),D) && is_path(vv[3:end],D) && is_path(reindex(vv,1:2,4:-1:3))
            (true, DT_e7)
        elseif (l,r) == (NoEnd, TwoOneBranch) && card == 8 && is_path(drop_slice(vv,3,3),D) && is_path(vv[3:end],D) && is_path(reindex(vv,1:2,4:-1:3))
            (true, DT_e8)
        elseif (l,r) == (NoEnd, TwoOneBranch) && card == 9 && is_path(drop_slice(vv,3,3),D) && is_path(vv[3:end],D) && is_path(reindex(vv,1:2,4:-1:3))
            (true, DT_E8)
        elseif (l,r) == (NoEnd, TwoTwoBranch) && card == 7 && is_path(drop_slice(vv,3,4),D) && is_path(vv[3:end],D) && is_path(reindex(vv,1:2,5:-1:3))
            (true, DT_E6)
        elseif (l,r) == (NoEnd, ThreeOneBranch) && card == 8 && is_path(drop_slice(vv,4,4),D) && is_path(vv[4:end],D) && is_path(reindex(vv,1:3,5:-1:4))
            (true, DT_E7)
        else 
            (false, "oh boy")
        end


    end
    
end


function is_simple_end(c::ConnectedInducedSubDiagram, v::Integer)
    if c.left_end == NoEnd && v == c.vertices[1]
        return true
    elseif c.right_end == NoEnd && v == c.vertices[-1]
        return true
    else
        return false
    end
end


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

function is_finite_volume(dag,dim)
    # Has spherical/parabolic diagram of rank n-1
    n = dim-1
    has_spherical_sub_of_rank_n = false
    has_affine_sub_of_rank_n_minus_1 = false

    subs = dag.subs

    for (support, subdiagram) in subs
        
        if length(support) == n && is_spherical(subdiagram.type)
            has_spherical_sub_of_rank_n = true
        end
        if length(support) == n && is_affine(subdiagram.type)
            has_affine_sub_of_rank_n_minus_1 = true
        end

        if length(support) == n-1 && is_spherical(subdiagram.type)
            extensions = [
                ext_support for (ext_support, ext_subdiagram) in subs if 
                support < ext_support && 
                (
                    ( length(ext_support) == n && is_spherical(ext_subdiagram.type) ) || 
                    ( length(ext_support) == n && is_affine(ext_subdiagram.type) ) 

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
the_singleton(v::Integer) = ConnectedInducedSubDiagram([v],NoEnd,NoEnd)

function small_diagram_type(VS::BitSet,D)
    @assert length(VS) ≤ 9 "only small diagrams"
    VS = Array(VS)
    
    if length(VS) == 2 && D[VS[1],VS[2]] == 3
        return ConnectedInducedSubDiagram(VS,NoEnd,NoEnd,DT_a)
    elseif length(VS) == 2 && D[VS[1],VS[2]] == 4
        return ConnectedInducedSubDiagram(VS,NoEnd,ASingleFour,DT_b)
    elseif length(VS) == 2 && D[VS[1],VS[2]] == 5
        return ConnectedInducedSubDiagram(VS,NoEnd,ASingleNgeq5,DT_h2)
    elseif length(VS) == 2 && D[VS[1],VS[2]] == 6
        return ConnectedInducedSubDiagram(VS,NoEnd,ASingleNgeq5,DT_g2)
    elseif length(VS) == 2 && D[VS[1],VS[2]] ≥ 7
        return ConnectedInducedSubDiagram(VS,NoEnd,ASingleNgeq5,DT_i)
    elseif length(VS) == 2 && D[VS[1],VS[2]] ≥ -1
        return ConnectedInducedSubDiagram(VS,NoEnd,ASingleNgeq5,DT_I)
    end
    
    if any(D[u,v] == 4 for u in VS for v in VS if u ≠ v)
        # f4 or F4 or bn or Bn or Cn
    end

    if any(D[u,v] == 5 for u in VS for v in VS if u ≠ v)
        # H3 or H4
    end
    
    if any(D[u,v] == 6 for u in VS for v in VS if u ≠ v)
        # G2
    end


end

function try_extend(VS::BitSet,S::InducedSubDiagram,D,v::Int)
    components = S.connected_components
    
    # joined should/will be of type ConnectedInducedSubDiagram
    joined = nothing # Here is the result
    
    # special case, no component in S, so we just return the singleton
    if length(components) == 0
        joined = singleton(v)
        return InducedSubDiagram(Set([joined]))
    end

    # list of components, corresponding adjacencent vertices 
    neighbors_v = Set([u for u in eachindex(D[v,:]) if u≠v && u in VS && D[v,u]≠2])
    neighbors_v_labels = [D[u,v] for u in neighbors_v] 
    non_neighboring_components = []
    neighboring_components = []
    neighboring_components_vertices = []
    neighboring_components_vertices_labels = []
    for c in components
        neighbors_in_c = [u for u in c.vertices if u in neighbors_v]
        if length(neighbors_in_c) > 0
            push!(neighboring_components,c)
            push!(neighboring_components_vertices,neighbors_in_c)
            push!(neighboring_components_vertices_labels,[D[u,v] for u in neighbors_in_c])
        else
            push!(non_neighboring_components,c)
        end
    end
   
    # TODO sort by types of connected components and then by cardinality 

    nv = neighbors_v 
    nc = neighboring_components
    nct = neighboring_components_types
    ncv = neighboring_components_vertices
    ncs = neighboring_components_sizes
    ncvl = neighboring_components_vertices_labels

    only_simple_edges =  all(l == 3 for l in neighbors_v_labels)



    if false
        @assert false "For alignment's sake"

    elseif sum(ncs) + 1 ≤ 9 # all sporadic & small subgraphs covered here 
        joined = small_diagram_type(reduce(∪,[c.vertices for c in nc],D) ∪ BitSet([v]))
    
    elseif length(nv) == 1
         
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
        elseif nct[1] == DT_a && ncvl[1] == 3 && nv == nc[1].left_end[1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_a,[v,nc[1].left_end[2]],nc[1].right_end)
        elseif nct[1] == DT_a && ncvl[1] == 3 && nv == nc[1].right_end[1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_a,nc[1].left_end,[v,nc[1].left_right[2]])
        elseif nct[1] == DT_b && ncvl[1] == 3 && nv == nc[1].left_end[1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_b,[v,nc[1].left_end[2]],nc[1].right_end)
        elseif nct[1] == DT_d && ncvl[1] == 3 && nv == nc[1].left_end[1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_d,[v,nc[1].left_end[2]],nc[1].right_end)
        
        # These correspond to adding a  v ==== * edge to a spherical (non-sporadic) diagram, thus resulting in a b,B, or C diagram
        elseif nct[1] == DT_a && ncvl[1] == 4 && nv == nc[1].left_end[1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_b,[v,nc[1].left_end[2]],nc[1].right_end) 
        elseif nct[1] == DT_a && ncvl[1] == 4 && nv == nc[1].right_end[1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_b,[v,nc[1].right_end[2]],nc[1].left_end) # By convention, the right end is the special one, hence the swap
        elseif nct[1] == DT_b && ncvl[1] == 4 && nv == nc[1].left_end[1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_C,[v,nc[1].left_end[2]],nc[1].right_end)
        elseif nct[1] == DT_d && ncvl[1] == 4 && nv == nc[1].left_end[1]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_B,[v,nc[1].left_end[2]],nc[1].right_end)
        
        # These correspond to adding a  "third" edge to just before the end of a spherical (non-sporadic) diagram, thus resulting in a d,B or D diagram
        elseif nct[1] == DT_a && ncvl[1] == 3 && length(nc[1].left_end) == 2 && nv == nc[1].left_end[2]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_d,[v,nc[1].left_end[1]],nc[1].right_end)
        elseif nct[1] == DT_a && ncvl[1] == 3 && length(nc[1].right_end) == 2 && nv == nc[1].right_end[2]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_d,[v,nc[1].right_end[1]],nc[1].left_end) # By convention, the right end is the special one, hence the swap
        elseif nct[1] == DT_b && ncvl[1] == 3 && length(nc[1].left_end) == 2 && nv == nc[1].left_end[2]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_B,[v,nc[1].left_end[1]],nc[1].right_end)
        elseif nct[1] == DT_d && ncvl[1] == 3 && length(nc[1].left_end) == 2 && nv == nc[1].left_end[2]
            joined = ConnectedInducedSubDiagram(push(nc.v,v),DT_D,[v,nc[1].left_end[1]],nc[1].right_end)
        end

    elseif length(nv) == 2

       # v has 2 neighbors, either twice the same connected component, in which case it's of type a, and adding v makes it type A
       #                    or two connected components…
        
       if false
            @assert false "For alignment's sake"
        
        # type a to type A
        elseif length(nc) == 1 && nct[1] == DT_a && ncvl == [3,3] && nv == [nc[1].left_end[1],nc[1].right_end[1]]
            joined = ConnectedInducedSubDiagram(push(nc.vertices,v),DT_A,[],[])
        elseif length(nc) == 1 && nct[1] == DT_a && ncvl == [3,3] && nv == [nc[1].right_end[1],nc[1].left_end[1]]
            joined = ConnectedInducedSubDiagram(push(nc.vertices,v),DT_A,[],[])

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


    else
        return nothing
    end

end

#function try_extend3(VS::BitSet,S::InducedSubDiagram,D,v::Integer)
#
#    components = S.connected_components
#    
#    # joined should/will be of type ConnectedInducedSubDiagram
#    joined = Nothing # Here is the result
#    
#    # special case, no component in S, so we just return the singleton
#    if length(components) == 0
#        joined = singleton(v)
#        return InducedSubDiagram(Set([joined]))
#    end
#    
#
#    # list of components, corresponding adjacencent vertices 
#    neighbors_v = Set([u for u in eachindex(D[v,:]) if u≠v && u in VS && D[v,u]≠2])
#    neighbors_v_labels = [D[u,v] for u in neighbors_v] 
#    non_neighboring_components = []
#    neighboring_components = []
#    neighboring_components_vertices = []
#    neighboring_components_vertices_labels = []
#    for c in components
#        neighbors_in_c = [u for u in c.vertices if u in neighbors_v]
#        if length(neighbors_in_c) > 0
#            push!(neighboring_components,c)
#            push!(neighboring_components_vertices,neighbors_in_c)
#            push!(neighboring_components_vertices_labels,[D[u,v] for u in neighbors_in_c])
#        else
#            push!(non_neighboring_components,c)
#        end
#    end
#   
#    # TODO sort by size of connected components, and then by edge labels
#
#    nv = neighbors_v 
#    nc = neighboring_components
#    nct = neighboring_components_types
#    ncv = neighboring_components_vertices
#    ncs = neighboring_components_sizes
#    ncvl = neighboring_components_vertices_labels
#
#    only_simple_edges =  all(l == 3 for l in neighbors_v_labels)
#
#    if false
#        @assert false "unreachable"
#    elseif sum(ncs) + 1 ≤ 9 # all sporadic subgraphs covered here 
#        joined = small_diagram_type(reduce(∪,[c.vertices for c in nc],D) ∪ BitSet([v]))
#    elseif ncs == [ncs[1],1,1] && nct == [DT_a,DT_a,DT_a] && ncvl == [3,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # d
#    elseif ncs == [ncs[1],1,1] && nct == [DT_b,DT_a,DT_a] && ncvl == [3,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # B
#    elseif ncs == [ncs[1],1,1] && nct == [DT_d,DT_a,DT_a] && ncvl == [3,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # D
#    elseif ncs == [ncs[1]] && nct == [DT_a] && ncvl == [3,3] && nct == [NoEnd,NoEnd]
#        # A
#        
#
#    
#    end
#
#
#end
#
#function try_extend2(VS::BitSet,S::InducedSubDiagram,D,v::Integer)
#    
#    components = S.connected_components
#    
#    # joined should/will be of type ConnectedInducedSubDiagram
#    joined = Nothing # Here is the result
#    
#    # special case, no component in S, so we just return the singleton
#    if length(components) == 0
#        joined = singleton(v)
#        return InducedSubDiagram(Set([joined]))
#    end
#    
#
#    # list of components, corresponding adjacencent vertices 
#    neighbors_v = Set([u for u in eachindex(D[v,:]) if u≠v && u in VS && D[v,u]≠2])
#    neighbors_v_labels = [D[u,v] for u in neighbors_v] 
#    non_neighboring_components = []
#    neighboring_components = []
#    neighboring_components_vertices = []
#    neighboring_components_vertices_labels = []
#    for c in components
#        neighbors_in_c = [u for u in c.vertices if u in neighbors_v]
#        if length(neighbors_in_c) > 0
#            push!(neighboring_components,c)
#            push!(neighboring_components_vertices,neighbors_in_c)
#            push!(neighboring_components_vertices_labels,[D[u,v] for u in neighbors_in_c])
#        else
#            push!(non_neighboring_components,c)
#        end
#    end
#   
#    # TODO sort by size of connected components, and then by edge labels
#
#    nv = neighbors_v 
#    nc = neighboring_components
#    nct = neighboring_components_types
#    ncv = neighboring_components_vertices
#    ncs = neighboring_components_sizes
#    ncvl = neighboring_components_vertices_labels
#
#    only_simple_edges =  all(l == 3 for l in neighbors_v_labels)
#
#    if false
#        @assert false "unreachable"
#    elseif ncs == [2,2,1] && nct == [DT_a,DT_a,DT_a] && ncvl == [3,3,3]
#        # e6
#    elseif ncs == [3,2,1] && nct == [DT_a,DT_a,DT_a] && ncvl == [3,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # e7
#    elseif ncs == [4,2,1] && nct == [DT_a,DT_a,DT_a] && ncvl == [3,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # e8
#    elseif ncs == [2,2,2] && nct == [DT_a,DT_a,DT_a] && ncvl == [3,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # E6
#    elseif ncs == [3,3,1] && nct == [DT_a,DT_a,DT_a] && ncvl == [3,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # E7
#    elseif ncs == [5,2,1] && nct == [DT_a,DT_a,DT_a] && ncvl == [3,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # E8
#    elseif ncs == [ncs[1],1,1] && nct == [DT_a,DT_a,DT_a] && ncvl == [3,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # d
#    elseif ncs == [1,1,1] && nct == [DT_a,DT_a,DT_a] && ncvl == [4,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # B4
#    elseif ncs == [ncs[1],1,1] && nct == [DT_b,DT_a,DT_a] && ncvl == [3,3,3] && nct == [NoEnd,NoEnd,NoEnd]
#        # B
#    elseif ncs == [ncs[1]] && nct == [DT_a] && ncvl == [3] && nct == [NoEnd]
#        # a
#    elseif ncs == [ncs[1]] && nct == [DT_a] && ncvl == [3,3] && nct == [NoEnd,NoEnd]
#        # A
#    elseif ncs == [ncs[1]] && nct == [DT_a] && ncvl == [3] && (nct == [2] || nct == [-2])
#        # e6, e7, e8 or E8
#    elseif ncs == [ncs[1]] && nct == [DT_a] && ncvl == [3] && (nct == [2] || nct == [-2])
#        # e6, e7, e8 or E8
#    end
#
#
#
#end
#
#
#function try_extend(VS::BitSet,S::InducedSubDiagram,D,v::Integer)
#    
#    println("____________ try_extend ________________")
#    println("$S with $v along")
#    display(D)
#    println()
#    println("________________________________________")
#
#    # We DO NOT DO the sporadic diagrams yet!
#    
#    components = S.connected_components
#    
#    # joined should/will be of type ConnectedInducedSubDiagram
#    joined = Nothing # Here is the result
#    
#    # special case, no component in S, so we just return the singleton
#    if length(components) == 0
#        joined = singleton(v)
#        return InducedSubDiagram(Set([joined]))
#    end
#
#    
#    neighbors_v = Set([u for u in eachindex(D[v,:]) if u≠v && u in VS && D[v,u]≠2])
#    non_neighboring_components = []
#    neighboring_components = []
#    for c in components
#        neighbors_in_c = [u for u in c.vertices if u in neighbors_v]
#        if length(neighbors_in_c) > 0
#            push!(neighboring_components,(c,neighbors_in_c))
#        else
#            push!(non_neighboring_components,c)
#        end
#    end
#
#     
#
#    if length(neighbors_v) > 3 # no subgraph has vertices of valency > 3
#        return Nothing
#    end
#
#    if length(neighbors_v) == 3 && any(D[u,v] ≠ 3 for u in neighbors_v) # if valency is == 3, labels must all be 3
#        return Nothing
#    end
#    
#    if length(neighbors_v) == 2 && all(D[u,v] ≠ 3 for u in neighbors_v) # if valency is == 2, at least one label must be 3
#        return Nothing
#    end
#   
#    if length(neighbors_v) == 3
#
#        for (c,neighbors_in_c) in neighboring_components 
#            if length(neighbors_in_c) ≥ 2 # If v is trivalent, it will make no loop/ in other words all components linked through v are linked via only one of their vertices
#                return Nothing
#            end
#            if D[neighbors_in_c[1],v] ≠ 3 # All labels must be equal 3
#                return Nothing
#            end
#        end
#        
#        @assert length(neighboring_components) == 3
#        
#        # TODO sporadic cases
#        #
#        sort!(neighboring_components,by=((c,u) -> card(c)))
#        
#        nc1,u1 = neighboring_components[1]
#        nc2,u2 = neighboring_components[2]
#        nc3,u3 = neighboring_components[3]
#
#        if card(nc2) > 1
#            return Nothing
#        else
#            if is_simple_end(nc3,u3)
#                if nc3.vertices[-1] == u3[1]
#                    joined = ConnectedInducedSubDiagram(nc3.vertices * [v] * nc1.vertices * nc2.vertices,nc3.right_end,ASimpleBranch)
#                elseif nc3.vertices[1] == u3[1]
#                    joined = ConnectedInducedSubDiagram( nc1.vertices * nc2.vertices * [v] * nc3.vertices, ASimpleBranch,nc3.right_end)
#                else
#                    @assert false "unreachable"
#                    return Nothing
#                end
#            else
#                return Nothing
#            end
#        end
#
#    
#
#    elseif length(neighbors_v) == 2 && length(neighboring_components) == 2
#
#        
#        C1,u1 = neighboring_components[1]
#        u1 = u1[1]
#        l1 = D[u1,v]
#        V1 = C1.vertices
#
#        C2,u2 = neighboring_components[2]
#        u2 = u2[1]
#        l2 = D[u2,v]
#        V2 = C2.vertices
#        
#        joined = Nothing
#
#        if l1 == 3 && l2 == 4 && card(C1) == 1
#            joined = ConnectedInducedSubDiagram([u1,v,u2],NoEnd,ASingleFour)
#        elseif l2 == 3 && l1 == 4 && card(C1) == 1
#            joined = ConnectedInducedSubDiagram([u1,v,u2],ASingleFour,NoEnd)
#        elseif l1 == 3 && l2 == 3 && u1 == V1[-1] && u2 == V2[1]
#            joined = ConnectedInducedSubDiagram([V1 ; [v]; V2],C1.left_end,C2.right_end)
#        elseif l1 == 3 && l2 == 3 && u1 == V1[1] && u2 == V2[-1]
#            joined = ConnectedInducedSubDiagram([V2 ; [v] ; V1],C2.left_end,C1.right_end)
#        elseif l1 == 3 && l2 == 3 && u1 == V1[1] && u2 == V2[1]
#            joined = ConnectedInducedSubDiagram([reverse(V1) ; [v] ; V1],C2.right_end,C2.right_end)
#        elseif l1 == 3 && l2 == 3 && u1 == V1[-1] && u2 == V2[-1]
#            joined = ConnectedInducedSubDiagram([V1 ; [v] ; reverse(V2)],C2.left_end,C2.left_end)
#        else
#            return Nothing
#        end
#
#        return 
#
#    elseif length(neighbors_v) == 2 && length(neighboring_components) == 1
#        println("neighbors of $v are $neighbors_v")
#        println("neighboring_components are $neighboring_components")
#        c,U = neighboring_components[1]
#        println("U is $U")
#        V = c.vertices
#        u1 = U[1]
#        u2 = U[2]
#        l1 = D[u1,v]
#        l2 = D[u2,v]
#        if u1 == V[1] && u2 == V[end] && l1 == 3 && l2 == 3
#            joined = ConnectedInducedSubDiagram([V ; [v]], ALoop, ALoop)
#        elseif u2 == V[1] && u1 == V[end] && l1 == 3 && l2 == 3
#            joined = ConnectedInducedSubDiagram([V ; [v]], ALoop, ALoop)
#        else
#            return Nothing
#        end
#    elseif length(neighbors_v) == 1
#        println("neighbors of $v are $neighbors_v")
#        println("neighboring_components are $neighboring_components")
#        C,U = neighboring_components[1]
#        V = C.vertices
#        u = U[1]
#        l = D[u,v]
#    
#        if is_simple_end(C,u) && l == 3
#            if C.left_end == NoEnd && V[1] == u
#                joined = ConnectedInducedSubDiagram([[v] ; V], NoEnd, C.right_end)
#            elseif C.right_end == NoEnd && V[-1] == u
#                joined = ConnectedInducedSubDiagram([V ; [v]], C.left_end, NoEnd)
#            else
#                @assert false "unreachable" 
#            end
#        elseif is_simple_end(C,u) && l == 4
#            if C.left_end == NoEnd && V[1] == u
#                joined = ConnectedInducedSubDiagram([[v] ; V], ASingleFour, C.right_end)
#            elseif C.right_end == NoEnd && V[-1] == u
#                joined = ConnectedInducedSubDiagram([V ; [v]], C.left_end, ASingleFour)
#            else
#                @assert false "unreachable" 
#            end
#        elseif card(C) ≥ 3 && V[2] == u && l == 3
#                joined = ConnectedInducedSubDiagram([[v] ; V], ASimpleBranch, C.right_end)
#        elseif card(C) ≥ 3 && V[-2] == u && l == 3
#                joined = ConnectedInducedSubDiagram([V ; [v]], C.left_end, ASimpleBranch)
#        else
#            return Nothing
#        end
#
#        # TODO 
#       
#    elseif length(neighbors_v) == 0
#        # Easy case:
#        joined = singleton(v)
#    else
#        @assert false "We shouldn't be here!"
#    end
#    
#    if joined == Nothing
#        return Nothing
#    end
#    
#    println("########## We have ##################")
#    display(D)
#    println()
#    println(joined)
#    println("########## ####### ##################")
#    @assert is_subdiagram(joined,D,size(D,1))
#    @assert is_subdiagram(joined,D,size(D,1))
#
#    return InducedSubDiagram(Set([joined])∪Set(non_neighboring_components))
#
#end

function extend(DAS::DiagramAndSubs, v::Array{Int,1})

    D = DAS.D
    subs = DAS.subs

    n = length(v) 
    @assert size(D) == (n,n)

    # Extend D with v
    D = [D v]
    D = [D; [v;0]']
   
    new_vertex = n+1

    new_subs::Dict{BitSet,InducedSubDiagram} = Dict()
    for (V,S) in subs
        S_and_v = try_extend(V,S,D,new_vertex)
        if S_and_v ≠ Nothing && S_and_v ≠ nothing 
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
    println(DAS)

end

function test_non_sporadic_connected_diagrams()
    
    # tripod with center the vertex 4
    tripod = [0 2 2 3;
              2 0 2 3;
              2 2 0 3;
              3 3 3 0]
    @assert(is_subdiagram(ConnectedInducedSubDiagram([1, 2, 4, 3],ASimpleBranch,NoEnd),tripod,4)[1])
    @assert(!is_subdiagram(ConnectedInducedSubDiagram([4, 1, 2, 3],ASimpleBranch,NoEnd),tripod,4)[1], "it's a tripod but the description is invalid")


end
