## Copied from qsolve.py

using Convex, Cbc, COSMO
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

Rat = Rational{Int}

iff(a::Bool,b::Bool) = a&&b || (!a)&&(!b)

function diag_product(x::SizedVector{N},y::SizedVector{N},z::SizedVector{N}) where {N}
    val = 0
    @inbounds for i in 1:N
       val += x[i]*y[i]*z[i]
    end
    return val
end
function diag_product(n::Int,x::SizedVector{N},y::SizedVector{N},z::SizedVector{N}) where {N}
    @assert 1 ≤ n && n ≤ N
    val = 0
    @inbounds for i in n:N
       val += x[i]*y[i]*z[i]
    end
    return val
end
function dot(n::Int,x::SizedVector{N},y::SizedVector{N}) where {N}
    @assert 1 ≤ n && n ≤ N
    val = 0
    @inbounds for i in n:N
       val += x[i]*y[i]
    end
    return val
end

function signature(G)
    
    D,P = diagonalize(G)
    pos = filter(>(0),diag(D))
    neg = filter(<(0),diag(D))
    
    return (length(pos),length(neg))
end

function diagonalize(A::Array{Int,2}) 
    return diagonalize(SizedMatrix{size(A)[1],size(A)[2],Int}(A))
end

function diagonalize(A::SizedMatrix{dim,dim,Int}) where {dim}
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
    I_ = SizedMatrix{dim,dim,Int128}(I)
    M::MMatrix{dim,2*dim} = [A I_]
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
   

    D::MMatrix{dim,dim,Int} = M[1:n,1:n]
    Q::MMatrix{dim,dim,Int} = M[1:n,n+1:2*n]
    P::MMatrix{dim,dim,Int} = Q'
     


    @assert LinearAlgebra.isdiag(D) "D is diagonal", D
    @assert P'*A*P == D "We have a diagonalization"
   

    # If the signature is (n,1), we put the negative eigenvalue first
    if dim > 1 && D[1,1] ≥ 0
        for i in 2:dim
            if D[i,i] < 0
                
                # Construct a transposition matrix
                T = MMatrix{dim,dim,Int}(I(dim))
                T[i,1] = 1
                T[1,i] = 1
                T[i,i] = 0
                T[1,1] = 0
                
                # transpose on the diagonal
                D[1,1],D[i,i] = D[i,i],D[1,1]
                # apply the transposition to P
                P = P*T
                
                break
            end
        end
    end
    
    @assert LinearAlgebra.isdiag(D) "D is diagonal", D
    @assert P'*A*P == D "We have a diagonalization"
    @assert A == inv(P')*D*inv(P) "We have a diagonalization"

    return (SizedMatrix{dim,dim,Int}(D),SizedMatrix{dim,dim,Int}(P))

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

function get_integer_points(
    M::SizedMatrix{rank,rank,Int}
)::Vector{SizedVector{rank,Int}} where {rank}
  
    invM = inv(M) 
    # TODO: Maybe we want an exact inverse: invM = inv(Rational{Int}.(M)) 
    # But I think approximate is enough


    n = size(M,1)
    @assert (n,n) == size(M) "Matrix should be square (and invertible)"
    @assert det(M) ≠ 0 "Matrix should be square (and invertible)"

    minimum = [sum([min(v,0) for v in M[i,:]]) for i in 1:n]
    maximum = [sum([max(v,0) for v in M[i,:]]) for i in 1:n]
    interval = [collect(min:max) for (min,max) in zip(minimum,maximum)]

    function parallelipiped_contains(v::SizedVector{rank,Int})
        Q = invM*v
        return all(c < 1 && c >= 0 for c in Q)
    end
    
    bb = SizedVector{rank,Int}.(Base.product(interval...))
    integer_points = [b for b in bb if parallelipiped_contains(b)] 
    
    dtr = Int(abs(det(Rational{Int}.((M)))))
    #integer_points = [v for v in bb if parallelipiped_contains(v)]
    @toggled_assert length(integer_points) == dtr "index = |determinant| = volume (I think)\n but have $(length(integer_points)) ≠ $(dtr)"
    # TODO is the discrete volume equal always?

    return integer_points
end

function get_integer_points_with_coordinates(
    M::SizedMatrix{rank,rank,Int},
    absdet::Int
)::Tuple{
    Vector{SizedVector{rank,Int}},
    Vector{SizedVector{rank,Rational{Int}}},
} where {rank}
  


    n = size(M,1)
    
    @toggled_assert (n,n) == size(M) "Matrix should be square (and invertible)"
    @toggled_assert abs(det(Rational{Int}.(M))) == absdet
    @toggled_assert absdet ≠ 0 "Matrix should be square (and invertible)"
    
    invM = inv(Rational{Int}.(M)) 

    minimum = [sum([min(v,0) for v in M[i,:]]) for i in 1:n]
    maximum = [sum([max(v,0) for v in M[i,:]]) for i in 1:n]
    interval = [collect(min:max) for (min,max) in zip(minimum,maximum)]

    function parallelipiped_contains(v::SizedVector{rank,Int})
        coordinates = invM*v
        if all(c < 1 && c >= 0 for c in coordinates) 
            return coordinates
        else
            return nothing
        end
    end
    
    bb = SizedVector{rank,Int}.(Base.product(interval...))
    integer_points = SizedVector{absdet,SizedVector{rank,Int}}(undef);
    integer_points_coordinates = SizedVector{absdet,SizedVector{rank,Rational{Int}}}(undef);
    idx = 1
    for b in bb
        c = parallelipiped_contains(b)
        if !isnothing(c)
            integer_points[idx] = b
            integer_points_coordinates[idx] = c
            idx+=1
        end
        if idx == absdet+1
            break
        end
    end
    
    #integer_points = [v for v in bb if parallelipiped_contains(v)]
    #@toggled_assert idx == absdet+1 "index = |determinant| = volume (I think)"

    return integer_points,integer_points_coordinates
end


function is_necessary_halfspace(
    D::SizedVector{rank,Int},
    cone_roots::Vector{SizedVector{rank,Int}},
    root::SizedVector{rank,Int},
   ) where {rank}

    return is_necessary_halfspace(Vector{SizedVector{rank,Int}}([D .* cr for cr in cone_roots]),D .* root)

end

function is_necessary_halfspace(cone_roots::Vector{SizedVector{rank,Int}},root::SizedVector{rank,Int}) where {rank}


    n = rank

    # x' * (A * r) ≤ 0 ∀ r
    # (A * r)' * x ≤ 0 ∀ r

    #x = Variable(n, IntVar)
    x = Variable(n)
    p = satisfy()       # satisfiability question 
    for cone_root in cone_roots
        p.constraints += x' * cone_root ≤ 0 # hyperplanes defining the cone
    end
    p.constraints += x' * root ≥ 1 # other side of the half space defined by root
    # it should only be strictly bigger than zero, but Convex.jl does not do "strictly", so we change it to ≥ 1 (and since we have a cone, it should be the same result)

    
    #solve!(p,Cbc.Optimizer(verbose=0,loglevel=0), verbose=false, warmstart=false)
    solve!(p,COSMO.Optimizer(verbose=false), verbose=false)
   

    if p.status == MathOptInterface.INFEASIBLE 
        return false
    elseif p.status == MathOptInterface.OPTIMAL
        #println(p.optval)
        return true
    else
        println("can't tell! ($(p.status))")
        println("             $(p))")
    end

end

#=
function is_necessary_halfspace(cone_roots::Vector{SizedVector{rank,Int}},A::SizedMatrix{rank,rank,Int},root::SizedVector{rank,Int}) where {rank}
    
    return is_necessary_halfspace((r -> A*r).(cone_roots),A*root)

end
=#

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
