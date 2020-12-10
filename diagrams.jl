

# Referring to https://en.wikipedia.org/wiki/Coxeter%E2%80%93Dynkin_diagram#Application_with_uniform_polytopes
# These are the different isomorphism types of subdiagrams of spherical/affine Coxeter-Dynkin diagrams


@enum EndType begin
    ASimpleBranch
    ASingleFour
    ASingleInfty
    ALoop
    NoEnd
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
function is_non_sporadic_subdiagram(c::ConnectedInducedSubDiagram,D,n)
    # D is the Coxeter matrix
    # c the subdiagram
    #
    @assert size(D) == (n,n)
    @assert D == D'
    vv = c.vertices
    
    println(all(v ≥ 1 && v ≤ n for v in vv))
    @needs(all(v ≥ 1 && v ≤ n for v in vv), "Elements of c must all be vertices")

    @needs(length(vv) == length(Set(vv)), "All vertices must appear exactly once")
    
    if card(c) == 0
    
        return true
    
    elseif card(c) == 1
        
        return true
    
    elseif card(c) == 2
        
        u = vv[1]
        v = vv[2]

        @iff(D[u,v] == 3 || D[u,v] == 4, "If only two vertices, the edge must be a 4 or a 3")
    

    elseif card(c) == 3
        if c.left_end == NoEnd && c.right_end == NoEnd
            @iff((D[vv[1],vv[2]], D[vv[2],vv[3]],D[vv[1],vv[3]]) == (3,3,2), 
                 "must have *-3-*-3-*")
        elseif c.left_end == NoEnd && c.right_end == ASingleFour 
            @iff((D[vv[1],vv[2]], D[vv[2],vv[3]],D[vv[1],vv[3]]) == (3,4,2),
                 "must have *-3-*-4-*")
        elseif c.left_end == ASingleFour && c.right_end == NoEnd 
            @iff((D[vv[1],vv[2]], D[vv[2],vv[3]],D[vv[1],vv[3]]) == (4,3,2), 
                 "must have *-4-*-3-*")
        elseif c.left_end == ALoop && c.right_end == ALoop 
            @iff((D[vv[1],vv[2]], D[vv[2],vv[3]],D[vv[1],vv[3]]) == (3,3,3), 
                 "Loop of length three, i.e. triangle")
        else
            return false
        end

    else
        
        @assert card(c) ≥ 4

        #forbidden end configurations
        ends = Set([c.left_end,c.right_end])

        @needs(!(ends == Set([ASingleFour]) 
               || (length(ends) == 2 && ALoop ∈ ends) 
               || ends == Set([ASingleFour,ASimpleBranch])), "Not all end combinations are legal")

    
        if c.left_end == ALoop && c.right_end == ALoop
            @needs(D[1,end] == 3, 
                   "If loop, it must close")
            @needs(all(D[vv[i],vv[i+1]] == 3 for i in 1:card(c)-1), 
                   "If loop, needs path")
            @needs(all(D[vv[i],vv[j]] == 2 for i in 1:card(c)-1 for j in i+2:card(c)), 
                   "If loop, no edge anywhere else")

            return true

        end

        start = 1
        stop = card(c)

        if c.left_end == ASingleFour 
            @needs(D[vv[1],vv[2]] == 4,
                   "Left end is *-4-*")
            start = 2
            @needs(all(D[vv[1],vv[i]] == 2 for i in start+1:card(c)),
                  "No connection between leftmost vertex and all others")

        elseif c.left_end == ASimpleBranch
            @needs(D[vv[1],vv[2]] == 2 && D[vv[1],vv[3]] == 3 && D[vv[2],vv[3]] == 3,
                   "Left end is :>*, i.e. a branch")
            
            start = 3
            println([(D[vv[2],vv[i]],D[vv[1],vv[i]])  for i in start+1:card(c)])
            @needs(all((D[vv[2],vv[i]],D[vv[1],vv[i]]) == (2,2) for i in start+1:card(c)),
                  "No connection between leftmost vertices and all others")

        elseif c.left_end == NoEnd
            start = 1
        end 
        
        if c.right_end == ASingleFour
            @needs(D[vv[end],vv[end-1]] == 4,
                   "Right end is *-4-*")
            
            stop = card(c) - 1
            @needs(all(D[vv[end],vv[i]] == 2 for i in 1:stop-1),
                  "No connection between rightmost vertex and all others")
        elseif c.right_end == ASimpleBranch
            @needs(D[vv[end],vv[end-1]] == 2 && D[vv[end],vv[end-2]] == 3 && D[vv[end-1],vv[end-2]] == 3,
                   "Right end is :>*, i.e. a branch")
            
            stop = card(c) - 2
            @needs(all((D[vv[end-1],vv[i]],D[vv[end],vv[i]]) == (2,2) for i in 1:stop-1),
                  "No connection between rightmost vertices and all others")

        elseif c.right_end == NoEnd
            stop = card(c)
        end
       
        for i in start:stop, j in start:stop
            if j-i == 1 
                @needs(D[vv[i],vv[j]] == 3, "consecutive vertices linked through simple edge")
            end
            if abs(j-i) ≥ 2  
                @needs(D[vv[i],vv[j]] == 2, "non-consecutive vertices not adjacent")
            end
        end

        return true
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
    @assert is_non_sporadic_subdiagram(joined,D,size(D,1))

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
    @assert(is_non_sporadic_subdiagram(ConnectedInducedSubDiagram([1, 2, 4, 3],ASimpleBranch,NoEnd),tripod,4))
    @assert(!is_non_sporadic_subdiagram(ConnectedInducedSubDiagram([4, 1, 2, 3],ASimpleBranch,NoEnd),tripod,4), "it's a tripod but the description is invalid")

    badtripod = [0 4 2 3;
              4 0 2 3;
              2 2 0 3;
              3 3 3 0]
    @assert(!is_non_sporadic_subdiagram(ConnectedInducedSubDiagram([1, 2, 4, 3],ASimpleBranch,NoEnd),badtripod,4))
    @assert(!is_non_sporadic_subdiagram(ConnectedInducedSubDiagram([4, 1, 2, 3],ASimpleBranch,NoEnd),badtripod,4))

end
