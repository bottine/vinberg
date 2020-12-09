

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


function the_empty_cisb() 
    return ConnectedInducedSubDiagram([],NoEnd,NoEnd)
end

function <(S1::ConnectedInducedSubDiagram,S2::ConnectedInducedSubDiagram)
    S1.vertices ⊆ S2.vertices
end

struct DiagramAndSubs 
    D::Array{Int,(2)}
    subs::Dict{BitSet,ConnectedInducedSubDiagram}
end

function all_compatible_pairs(DAS::DiagramAndSubs)
    return [((v1,s1),(v2,s2)) for (v1,s1) in DAS.subs, (v2,s2) in DAS.subs if isdisjoint(v1,v2)]
end

empty(S1::ConnectedInducedSubDiagram) = length(S1.vertices) == 0



function try_join(S1::ConnectedInducedSubDiagram,S2::ConnectedInducedSubDiagram,D,v::Integer)
    
    adj(u,w) =  D[u,w] ≠ 2 ? true : false
    # TODO
    # We don't do the sporadics diagrams for now
    

    


    
    if empty(S1) && empty(S2)
        # both empty means we just return the vertex `v`
        return ConnectedInducedSubDiagram([v],NoEnd,NoEnd) 
   
    end
    
    # We assume S1 and S2 are at distance ≥2 to begin with, because otherwise their union is a ConnectedInducedSubDiagram
    # of its own, and should be treated separately
    if ! all(!adj(u,w)  for u in S1.vertices for w in S2.vertices)
        # "Induced subgraphs must be at distance ≥ 2"
        return Nothing
    end

    if empty(S1) || empty(S2)
    
        if empty(S1)
            S1,S2 = S2,S1
        end
        
        @assert empty(S2) && !empty(S1) # just to be sure

        V1 = S1.vertices
        # neighbors of v, with edge label
        neighbors = [(u,D[u,v]) for u in V1 if adj(u,v)]


        if length(neighbors) == 0 
            return Nothing
        elseif length(neighbors) == 1
            
            (u,l) = neighbors[1]

            if l == 3
                
                if length(S1.vertices) == 1
                    return ConnectedInducedSubDiagram([[v]; V1], NoEnd, NoEnd)
                elseif length(V1) == 2 
                    return ConnectedInducedSubDiagram([[v] ; V1], NoEnd, NoEnd)
                elseif length(V1) == 3
                    if u == V1[1]
                        return ConnectedInducedSubDiagram([[v] ; V1], NoEnd, NoEnd)
                    elseif u == V1[2]
                        return ConnectedInducedSubDiagram([[v] ; V1], ASimpleBranch, NoEnd)
                    elseif u == V1[-1]
                        return ConnectedInducedSubDiagram([V1 ; [v]], NoEnd, NoEnd)
                    else
                        return Nothing    
                    end
                else
                     if u == V1[1]
                        return ConnectedInducedSubDiagram([[v] ; V1], NoEnd, NoEnd)
                    elseif u == V1[2]
                        return ConnectedInducedSubDiagram([[v] ; V1], ASimpleBranch, NoEnd)
                    elseif u == V1[-1]
                        return ConnectedInducedSubDiagram([V1 ; [v]], NoEnd, NoEnd)
                    elseif u == V1[-2]
                        return ConnectedInducedSubDiagram([V1 ; [v]], NoEnd, ASimpleBranch)
                    else
                        return Nothing    
                    end                   
                end

            elseif l == 4
                if u == V1[1] && S1.left_end == NoEnd
                    return ConnectedInducedSubDiagram([[v] ; V1], ASingleFour,S1.right_end)
                if u == V1[-1] && S1.right_end == NoEnd
                    return ConnectedInducedSubDiagram([V1 ; [v]], S1.left_end, ASingleFour)
                else
                    return Nothing
                end
                
            else
                #sporadic cases but we don't do them yet
                return Nothing
            end


            

        elseif length(neighbors) == 2
            # the only case in which the resulting diagram is valid is if we're closing a loop, i.e. getting an \tilde{A}_n diagram
            (u1,l1),(u2,l2) = neighbors[1],neighbors[2]
            if u1 == V1[1] && u2 == V1[-1] && l1 == 3 && l2 == 3
                return ConnectedInducedSubDiagram([V1 ; [v]], ALoop, ALoop)
            elseif u2 == V1[1] && u1 == V1[-1] && l1 == 3 && l2 == 3
                return ConnectedInducedSubDiagram([V1 ; [v]], ALoop, ALoop)
            else
                return Nothing
            end
        else
            return Nothing
        end 
        end

    else

        @assert !empty(S1) && !empty(S2) 

        

        if all(!adj(u,v) for u in S1.vertices) || all(!adj(u,v) for u in S2.vertices)
            return Nothing # v can't do any joining
        end
        
        # Now is the tedious part: 
        # we have to decide when two subdiagrams can be joined into another
        # *legal* subdiagram

        neighbors1 = [(u,D[u,v]) for u in S1.vertices if adj(u,v)]
        neighbors2 = [(u,D[u,v]) for u in S2.vertices if adj(u,v)]
        @assert !isempty(neighbors1)
        @assert !isempty(neighbors2)

        if length(neighbors1) > 1 || length(neighbors2) > 1
            return Nothing
        end

        V1 = S1.vertices
        V2 = S2.vertices
        
        (u1,l1) = neighbors1[1]
        (u2,l2) = neighbors2[1]

        if l1 == 3 && l2 == 4 && length(V2) == 1
            return ConnectedInducedSubDiagram([u1,v,u2],NoEnd,ASingleFour)
        elseif l2 == 3 && l1 == 4 && length(V1) == 1
            return ConnectedInducedSubDiagram([u1,v,u2],ASingleFour,NoEnd)
        elseif l1 == 3 && l2 == 3 && u1 == V1[-1] && u2 == V2[1]
            return ConnectedInducedSubDiagram([V1 ; [v]; V2],S1.left_end,S2.right_end)
        elseif l1 == 3 && l2 == 3 && u1 == V1[1] && u2 == V2[-1]
            return ConnectedInducedSubDiagram([V2 ; [v] ; V1],S2.left_end,S1.right_end)
        elseif l1 == 3 && l2 == 3 && u1 == V1[1] && u2 == V2[1]
            return ConnectedInducedSubDiagram([reverse(V1) ; [v] ; V1],S2.right_end,S2.right_end)
        elseif l1 == 3 && l2 == 3 && u1 == V1[-1] && u2 == V2[-1]
            return ConnectedInducedSubDiagram([V1 ; [v] ; reverse(V2)],S2.left_end,S2.left_end)
            
        else
            return Nothing
        end


    end
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

    new_subs::Dict{BitSet,ConnectedInducedSubDiagram} = Dict()
    for ((V1,S1),(V2,S2)) in all_compatible_pairs(DAS)
        S1_v_S2 = try_join(S1,S2,D,new_vertex)
        if S1_v_S2 ≠ Nothing && S1_v_S2 ≠ nothing 
            push!(new_subs,(V1∪V2∪BitSet([new_vertex])) => S1_v_S2)
        end
    end
    return DiagramAndSubs(D,merge(subs,new_subs))
    
end


function build_diagram_and_subs(M)
    n = size(M,1)
    @assert size(M) == (n,n)
    # @assert issymetric(M)
    
    lol = Dict()
    push!(lol,BitSet() => the_empty_cisb())

    DAS = DiagramAndSubs(reshape([],0,0),lol)
    for i in 1:n
        DAS = extend(DAS,M[i,1:i-1]) 
    end
    println(DAS)

end
