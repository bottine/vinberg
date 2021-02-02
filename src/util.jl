## Copied from qsolve.py

using Convex, Cbc
using LinearAlgebra
using Base
using MathOptInterface
using BitIntegers
using StaticArrays
using Memoize
using ToggleableAsserts
# To do diagonalization on StaticArrays, we need fixed-size element types, so BigInt is out.
# But we still need to have big integers.
# Hopefully Int512 is enough


iff(a::Bool,b::Bool) = a&&b || (!a)&&(!b)

function signature(G)
    
    D,P = diagonalize(G)
    pos = filter(>(0),diag(D))
    neg = filter(<(0),diag(D))
    
    return (length(pos),length(neg))
end

function diagonalize(A::Array{Int,2}) 
    return diagonalize(SMatrix{size(A)[1],size(A)[2],Int}(A))
end

function diagonalize(A::SMatrix{rank,rank,Int}) where {rank}
    # returns T and D with D = T'GT
    # algorithm copied from there https://math.stackexchange.com/questions/1388421/reference-for-linear-algebra-books-that-teach-reverse-hermite-method-for-symmetr
    # plus a gcd step to reduce the growth of values
    # plus the "Case 0" step to further reduce, but it's not enough
    @assert issymmetric(A) "A must be symmetric."


    # NO BigInt is it dangerous?
    A = Int512.(A)
    #A = Int128.(A)
    #A = BigInt.(A)

    n = size(A)[1]
    i0 = 1
    I_ = SMatrix{rank,rank,Int128}(I)
    M::MMatrix{rank,2*rank} = [A I_]
    while i0 ≤ n
       
      
        # look at non zero diagonal entries
        non_zero_diag = [k for k in i0:n if M[k,k] ≠ 0]
        non_zero_diag = sort!(non_zero_diag,by=(k -> abs(M[k,k])))

#        println("====================")
#        println("i0 = $i0")
#        display(M)
#        println("")
#        println("non_zero_diag = $non_zero_diag")
#    
        if length(non_zero_diag) == 0
            non_zero_coordinates = [(i,j) for  i in i0:n, j in i0:n if M[i,j]≠0]
            if isempty(non_zero_coordinates)
                break
            end
            (i,j) = (sort(non_zero_coordinates, by=(x-> abs(M[x[1],x[2]]))))[1]
            M[i,:] = M[i,:] + M[j,:]
            M[:,i] = M[:,i] + M[:,j]
        else
            
            k = non_zero_diag[1]
            M[i0,:], M[k,:] = M[k,:], M[i0,:]
            M[:,i0], M[:,k] = M[:,k], M[:,i0]

            for i in i0+1:n
                g =  gcd(M[i0,i0],M[i0,i])
                mizi = M[i0,i]//g
                miziz = M[i0,i0]//g
                M[i,:] = (-mizi*M[i0,:] + miziz*M[i,:])
                M[:,i] = (-mizi*M[:,i0] + miziz*M[:,i])
            end
            i0 = i0 + 1
        end
    end
   

    D::SMatrix{rank,rank,Int} = M[1:n,1:n]
    Q::SMatrix{rank,rank,Int} = M[1:n,n+1:2*n]
    P::SMatrix{rank,rank,Int} = Q'
   


    @assert LinearAlgebra.isdiag(D) "D is diagonal", D
    @assert P'*A*P == D "We have a diagonalization"
    
    return (D,P)

end


function test_diagonalize()
   
   
    for n in 4:20
        println("\nn is $n")
        for i in 1:4*(20-n)
            print(".")
            #M = rand(-20:20,n,n)
            M = rand(-20:20,n,n)
            M = M + M' # make it symmetric
            @assert LinearAlgebra.issymmetric(M)
            (D,P) = diagonalize(M)
            @assert P'*M*P == D

        end
    end

end 

function get_integer_points(M)
  
    invM = inv(M) 
    # TODO: Maybe we want an exact inverse: invM = inv(Rational{Int}.(M)) 
    # But I think approximate is enough


    n = size(M,1)
    @assert (n,n) == size(M) "Matrix should be square (and invertible)"
    @assert det(M) ≠ 0 "Matrix should be square (and invertible)"

    minimum = [sum([min(v,0) for v in M[i,:]]) for i in 1:n]
    maximum = [sum([max(v,0) for v in M[i,:]]) for i in 1:n]

    bounding_box = [[]] 
    for i in 1:n
        bounding_box = [vcat(vec, [val]) for vec in bounding_box for val in minimum[i]-1:maximum[i]+1] 
    end
    
    function parallelipiped_contains(v)
        Q = invM*v
        return all(c < 1 && c >= 0 for c in Q)
    end
    
    integer_points = [Vector(v) for v in bounding_box if parallelipiped_contains(Vector(v))]
    @toggled_assert length(integer_points) == abs(det(Rational{Int}.((M)))) "index = |determinant| = volume (I think)\n but have $(length(integer_points)) ≠ $(abs(det(Rational{Int}.(M))))"
    # TODO is the discrete volume equal always?

    return integer_points

end


@memoize Dict function is_necessary_halfspace(cone_roots::Vector{SVector{rank,Int}},root::SVector{rank,Int}) where {rank}


    n = rank

    # x' * (A * r) ≤ 0 ∀ r
    # (A * r)' * x ≤ 0 ∀ r

    x = Variable(n, IntVar)
    p = satisfy()       # satisfiability question 
    for cone_root in cone_roots
        p.constraints += x' * cone_root ≤ 0 # hyperplanes defining the cone
    end
    p.constraints += x' * root ≥ 1 # other side of the half space defined by root
    # it should only be strictly bigger than zero, but Convex.jl does not do "strictly", so we change it to ≥ 1 (and since we have a cone, it should be the same result)

    
    solve!(p,Cbc.Optimizer(verbose=0,loglevel=0), verbose=false, warmstart=false)
   

    if p.status == MathOptInterface.INFEASIBLE 
        return false
    elseif p.status == MathOptInterface.OPTIMAL
        #println(p.optval)
        return true
    else
        println("can't tell!")
    end

end


function is_necessary_halfspace(cone_roots::Vector{SVector{rank,Int}},A::SMatrix{rank,rank,Int},root::SVector{rank,Int}) where {rank}
    
    return is_necessary_halfspace((r -> A*r).(cone_roots),A*root)

end


function Gram_to_Coxeter(G::Matrix{Int})
    # stolen from B&P
    function weight(M, i, j)
        cos2 = (M[i,j]*M[i,j])//(M[i,i]*M[j,j])
        if cos2 == 0
            return 2
        elseif cos2 == 1//4
            return 3
        elseif cos2 == 1//2
            return 4
        elseif cos2 == 3//4
            return 6
        elseif cos2 == 1
            return 0
        elseif cos2 > 1
            return 1
        else
            return nothing
        end
    end
    Coxeter_matrix = reduce(hcat,[[weight(G,i,j) for i in 1:size(G,1)] for j in 1:size(G,2)])
    if any(g == nothing for g in Coxeter_matrix)
        return nothing
    else
        return Coxeter_matrix
    end

end

function txtmat_path_to_matrix(path)
    s = open(path) do file
        read(file, String)
    end

    return txtmat_to_matrix(s)
end

function txtmat_to_matrix(descr)
   
    lines = split(descr,"\n")
    @assert length(lines) ≥ 2 "CoxIter graph description must have non trivial content!"
    
    M = []
    
    for line in lines
        m = match(r"^(?<entries>(\s*-?\d+)*)\s*(?<comment>#.*)?$", line)
        if m ≠ nothing # labelled vertice
            entries = split(m[:entries])
            if ! isempty(entries)
                int_entries = (x->parse(Int,x)).(entries)
                push!(M,int_entries)
            end
        end
    end
    
    Mat = hcat(M...)
    
    return Mat

end

function matrix_to_txtmat(Mat)
    return join([join(string.(line)," ") for line in eachrow(Mat)], "\n") 
end

function matrix_to_txtmat_path(Mat, path)

    s = open(path,"w") do file
        print(file, matrix_to_txtmat(Mat))
    end
end
