using MLStyle

# Referring to https://en.wikipedia.org/wiki/Coxeter%E2%80%93Dynkin_diagram#Application_with_uniform_polytopes
# These are the different isomorphism types of subdiagrams of spherical/affine Coxeter-Dynkin diagrams
@enum DiagramType begin
    # Lowercase == Spherical
    DT_a
    DT_b # == C 
    DT_d
    DT_e6
    DT_e7
    DT_e8
    DT_f4
    DT_g2
    DT_h2
    DT_i
    DT_h3
    DT_h4
    DT_# Uppercase == Affine
    DT_A
    DT_B
    DT_C
    DT_D
    DT_E6
    DT_E7
    DT_E8
    DT_F4
    DT_G2
    DT_I
end

@enum EndType begin
    NoEnd = 0
    ASimpleBranch = 1
    ASingleFour = 2
    ALoop = 3
    # Above are the only ones needed for non sporadic
    TwoOneBranch = 4    # Needed for E types
    TwoTwoBranch = 5   # Needed for E types
    ThreeOneBranch = 6 # Needed for E types
    ASingleNgeq5  = 7  
    ASingleInfty = 8
    AFourAndThen = 9   # Needed for F_4

end

struct ConnectedInducedSubDiagram
    vertices::Array{Int,(1)}
    left_end::EndType
    right_end::EndType
end

card(c::ConnectedInducedSubDiagram) = length(c.vertices)

macro needs(condition, comment)
    return quote 
        if ! $(esc(condition))
            if $(esc(comment)) ≠ Nothing
                println($(esc(comment)))
            end
            return false
        end
    end
end

macro iff(condition, comment)
    return quote 
        if $(esc(comment)) ≠ Nothing
            println($(esc(comment)))
        end
        if $(esc(condition))
            return true
        else
            return false
        end
    end
end

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


function is_path(vv,D,first=(x -> x == 3),last=(x -> x == 3))
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
            if D[vv[i],vv[j]] ≠ 3
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
    return vcat( vv[1:from-1],  vv[to+1:end])
end

function is_subdiagram2(c::ConnectedInducedSubDiagram,D,n)
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
        elseif (l,r) == (NoEnd,NoEnd) && card ≤ 3 && is_path(vv,D) 
            (true, Nothing)
        elseif (l,r) == (NoEnd,NoEnd) && card ≥ 4 && is_path(vv,D) 
            (true, DT_a)
        elseif (l,r) == (NoEnd,ASingleFour) && card ≤ 3 && is_path(vv,D,last=(x -> x == 4)) 
            (true, Nothing)
        elseif (l,r) == (NoEnd,ASingleFour) && card ≥ 4 && is_path(vv,D,last=(x -> x == 4)) 
            (true, DT_b)
        elseif (l,r) == (NoEnd, ASimpleBranch) && card ≤ 4 && is_path(drop_slice(vv,card,card),D) && is_path(drop_slice(vv,card-1,card-1),D) && D[vv[end-1],vv[end]] == 2  
            (true, Nothing)
        elseif (l,r) == (NoEnd, ASimpleBranch)  && card ≥ 5 && is_path(drop_slice(vv,card,card),D) && is_path(drop_slice(vv,card-1,card-1),D) && D[vv[end-1],vv[end]] == 2  
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
        elseif (l,r) == (NoEnd, AFourAndThen) && card == 3 && is_path(vv,first=(x -> x == 4))
            (true, Nothing)
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
        elseif (l,r) == (NoEnd, TwoOneBranch) && card == 5 && is_path(reindex(vv,1:2,4:5),D) && is_path(vv[3:5],D) && is_path(reindex(vv,1:2,4:-1:3))
            (true, Nothing) # a piece of an E diagram
        elseif (l,r) == (NoEnd, TwoOneBranch) && card == 6 && is_path(drop_slice(vv,3,3),D) && is_path(vv[3:end],D) && is_path(reindex(vv,1:2,4:-1:3))
            (true, DT_e6)
        elseif (l,r) == (NoEnd, TwoOneBranch) && card == 7 && is_path(drop_slice(vv,3,3),D) && is_path(vv[3:end],D) && is_path(reindex(vv,1:2,4:-1:3))
            (true, DT_e7)
        elseif (l,r) == (NoEnd, TwoOneBranch) && card == 8 && is_path(drop_slice(vv,3,3),D) && is_path(vv[3:end],D) && is_path(reindex(vv,1:2,4:-1:3))
            (true, DT_e8)
        elseif (l,r) == (NoEnd, TwoOneBranch) && card == 9 && is_path(drop_slice(vv,3,3),D) && is_path(vv[3:end],D) && is_path(reindex(vv,1:2,4:-1:3))
            (true, DT_E8)
        elseif (l,r) == (NoEnd, TwoTwoBranch) && card == 6 && is_path(drop_slice(vv,3,4),D) && is_path(vv[3:end],D) && is_path(reindex(vv,1:2,5:-1:3))
            (true, Nothing) # A piece of E6
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
    subs::Dict{BitSet,InducedSubDiagram}
end


empty(S1::ConnectedInducedSubDiagram) = length(S1.vertices) == 0
singleton(v::Integer) = ConnectedInducedSubDiagram([v],NoEnd,NoEnd)


function try_extend(VS::BitSet,S::InducedSubDiagram,D,v::Integer)
    
    println("____________ try_extend ________________")
    println("$S with $v along")
    display(D)
    println()
    println("________________________________________")

    # We DO NOT DO the sporadic diagrams yet!
    
    components = S.connected_components
    
    joined = Nothing # Here is the result

    if length(components) == 0
        joined = singleton(v)
        return InducedSubDiagram(Set([joined]))
    end

    
    neighbors_v = Set([u for u in eachindex(D[v,:]) if u≠v && u in VS && D[v,u]≠2])
    non_neighboring_components = []
    neighboring_components = []
    for c in components
        neighbors_in_c = [u for u in c.vertices if u in neighbors_v]
        if length(neighbors_in_c) > 0
            push!(neighboring_components,(c,neighbors_in_c))
        else
            push!(non_neighboring_components,c)
        end
    end

    if length(neighbors_v) > 3 # no subgraph has vertices of valency > 3
        return Nothing
    end

    if length(neighbors_v) == 3 && any(D[u,v] ≠ 3 for u in neighbors_v) # if valency is == 3, labels must all be 3
        return Nothing
    end
    
    if length(neighbors_v) == 2 && all(D[u,v] ≠ 3 for u in neighbors_v) # if valency is == 2, at least one label must be 3
        return Nothing
    end
   
    if length(neighbors_v) == 3

        for (c,neighbors_in_c) in neighboring_components 
            if length(neighbors_in_c) ≥ 2 # If v is trivalent, it will make no loop/ in other words all components linked through v are linked via only one of their vertices
                return Nothing
            end
            if D[neighbors_in_c[1],v] ≠ 3 # All labels must be equal 3
                return Nothing
            end
        end
        
        @assert length(neighboring_components) == 3
        
        # TODO sporadic cases
        #
        sort!(neighboring_components,by=((c,u) -> card(c)))
        
        nc1,u1 = neighboring_components[1]
        nc2,u2 = neighboring_components[2]
        nc3,u3 = neighboring_components[3]

        if card(nc2) > 1
            return Nothing
        else
            if is_simple_end(nc3,u3)
                if nc3.vertices[-1] == u3[1]
                    joined = ConnectedInducedSubDiagram(nc3.vertices * [v] * nc1.vertices * nc2.vertices,nc3.right_end,ASimpleBranch)
                elseif nc3.vertices[1] == u3[1]
                    joined = ConnectedInducedSubDiagram( nc1.vertices * nc2.vertices * [v] * nc3.vertices, ASimpleBranch,nc3.right_end)
                else
                    @assert false "unreachable"
                    return Nothing
                end
            else
                return Nothing
            end
        end

    

    elseif length(neighbors_v) == 2 && length(neighboring_components) == 2

        
        C1,u1 = neighboring_components[1]
        u1 = u1[1]
        l1 = D[u1,v]
        V1 = C1.vertices

        C2,u2 = neighboring_components[2]
        u2 = u2[1]
        l2 = D[u2,v]
        V2 = C2.vertices
        
        joined = Nothing

        if l1 == 3 && l2 == 4 && card(C1) == 1
            joined = ConnectedInducedSubDiagram([u1,v,u2],NoEnd,ASingleFour)
        elseif l2 == 3 && l1 == 4 && card(C1) == 1
            joined = ConnectedInducedSubDiagram([u1,v,u2],ASingleFour,NoEnd)
        elseif l1 == 3 && l2 == 3 && u1 == V1[-1] && u2 == V2[1]
            joined = ConnectedInducedSubDiagram([V1 ; [v]; V2],C1.left_end,C2.right_end)
        elseif l1 == 3 && l2 == 3 && u1 == V1[1] && u2 == V2[-1]
            joined = ConnectedInducedSubDiagram([V2 ; [v] ; V1],C2.left_end,C1.right_end)
        elseif l1 == 3 && l2 == 3 && u1 == V1[1] && u2 == V2[1]
            joined = ConnectedInducedSubDiagram([reverse(V1) ; [v] ; V1],C2.right_end,C2.right_end)
        elseif l1 == 3 && l2 == 3 && u1 == V1[-1] && u2 == V2[-1]
            joined = ConnectedInducedSubDiagram([V1 ; [v] ; reverse(V2)],C2.left_end,C2.left_end)
        else
            return Nothing
        end

        return 

    elseif length(neighbors_v) == 2 && length(neighboring_components) == 1
        println("neighbors of $v are $neighbors_v")
        println("neighboring_components are $neighboring_components")
        c,U = neighboring_components[1]
        println("U is $U")
        V = c.vertices
        u1 = U[1]
        u2 = U[2]
        l1 = D[u1,v]
        l2 = D[u2,v]
        if u1 == V[1] && u2 == V[end] && l1 == 3 && l2 == 3
            joined = ConnectedInducedSubDiagram([V ; [v]], ALoop, ALoop)
        elseif u2 == V[1] && u1 == V[end] && l1 == 3 && l2 == 3
            joined = ConnectedInducedSubDiagram([V ; [v]], ALoop, ALoop)
        else
            return Nothing
        end
    elseif length(neighbors_v) == 1
        println("neighbors of $v are $neighbors_v")
        println("neighboring_components are $neighboring_components")
        C,U = neighboring_components[1]
        V = C.vertices
        u = U[1]
        l = D[u,v]
    
        if is_simple_end(C,u) && l == 3
            if C.left_end == NoEnd && V[1] == u
                joined = ConnectedInducedSubDiagram([[v] ; V], NoEnd, C.right_end)
            elseif C.right_end == NoEnd && V[-1] == u
                joined = ConnectedInducedSubDiagram([V ; [v]], C.left_end, NoEnd)
            else
                @assert false "unreachable" 
            end
        elseif is_simple_end(C,u) && l == 4
            if C.left_end == NoEnd && V[1] == u
                joined = ConnectedInducedSubDiagram([[v] ; V], ASingleFour, C.right_end)
            elseif C.right_end == NoEnd && V[-1] == u
                joined = ConnectedInducedSubDiagram([V ; [v]], C.left_end, ASingleFour)
            else
                @assert false "unreachable" 
            end
        elseif card(C) ≥ 3 && V[2] == u && l == 3
                joined = ConnectedInducedSubDiagram([[v] ; V], ASimpleBranch, C.right_end)
        elseif card(C) ≥ 3 && V[-2] == u && l == 3
                joined = ConnectedInducedSubDiagram([V ; [v]], C.left_end, ASimpleBranch)
        else
            return Nothing
        end

        # TODO 
       
    elseif length(neighbors_v) == 0
        # Easy case:
        joined = singleton(v)
    else
        @assert false "We shouldn't be here!"
    end
    
    if joined == Nothing
        return Nothing
    end
    
    println("########## We have ##################")
    display(D)
    println()
    println(joined)
    println("########## ####### ##################")
    @assert is_subdiagram(joined,D,size(D,1))
    @assert is_subdiagram2(joined,D,size(D,1))

    return InducedSubDiagram(Set([joined])∪Set(non_neighboring_components))

end

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
    @assert(is_subdiagram2(ConnectedInducedSubDiagram([1, 2, 4, 3],ASimpleBranch,NoEnd),tripod,4)[1])
    @assert(!is_subdiagram2(ConnectedInducedSubDiagram([4, 1, 2, 3],ASimpleBranch,NoEnd),tripod,4)[1], "it's a tripod but the description is invalid")

    badtripod = [0 4 2 3;
              4 0 2 3;
              2 2 0 3;
              3 3 3 0]
    @assert(!is_subdiagram(ConnectedInducedSubDiagram([1, 2, 4, 3],ASimpleBranch,NoEnd),badtripod,4))
    @assert(!is_subdiagram(ConnectedInducedSubDiagram([4, 1, 2, 3],ASimpleBranch,NoEnd),badtripod,4))

end
