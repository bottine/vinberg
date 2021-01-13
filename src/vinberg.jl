#using Nemo
#using Hecke
using AbstractAlgebra
using LinearAlgebra
using Polyhedra
using CDDLib


import Base: vec, convert


# Code adapted from N. V. Bogachev and A. Yu. Perepechko:
#
#   https://github.com/aperep/vinberg-algorithm


# TODO
#
# * See here for memoization: https://github.com/JuliaCollections/Memoize.jl
# * Extensive unit tests!
#


include("util.jl")
include("qsolve.jl")
include("diagrams.jl")
include("Some_Lattices.jl")

struct QuadLattice
    n_pos::Int
    n_neg::Int
    G::LinearAlgebra.Symmetric{Int,AbstractMatrix{Int}}
    P::AbstractMatrix{Int}
    D::Array{Int,1}
    last_invariant_factor::Int
    root_lengths::Vector{Int}
end

# Presumably any root must have length dividing twice the last invariant factor
function root_lengths(L::QuadLattice)
    return L.root_lengths
end



function QuadLattice(G)

    assert_sig_n_1_matrix(G)

    G = Symmetric(G) 

    np = size(G)[1]
    n_pos = np-1
    n_neg = 1
    D,P = diagonalize(G)
    @assert P'*G*P == D
    @assert isdiag(D)
    D = diag(D) # Get the diagonal vector of `D`

    rG = Rational{Int}.(G)
    cofactorsG = det(rG) * inv(rG) # https://stackoverflow.com/questions/58149645/julia-how-can-we-compute-the-adjoint-or-classical-adjoint-linear-algebra
    last_invariant_factor = abs(Int(det(rG)//gcd(cofactorsG)))

    twice_LIF = 2*last_invariant_factor
    root_lengths =  [k for k in 1:twice_LIF if twice_LIF%k == 0]
    
    return QuadLattice(n_pos,n_neg,G,P,D,last_invariant_factor,root_lengths)
end


function Base.:(==)(L1::QuadLattice,L2::QuadLattice)
    L1.G == L2.G
end



function rank(L::QuadLattice)
    return L.n_pos + L.n_neg
end

struct QuadLatticeElement
    L::QuadLattice
    vec::Array{Int,1}
    #vec::Array{Int,1}
end
function Base.isequal(v1::QuadLatticeElement,v2::QuadLatticeElement)
    v1.L == v2.L && v1.vec == v2.vec
end

function Base.:(==)(v1::QuadLatticeElement,v2::QuadLatticeElement)
    v1.L == v2.L && v1.vec == v2.vec
end
struct VinbergLattice
    # By which we mean, a quadratic lattice along with
    # * An element of negative norm v₀
    # * A basis V1_basis of the sublattice v₀^\perp = V₁
    # * A set of representatives W_reps for the quotient of L by V₁⊕ℤv₀

    L::QuadLattice
    v0::QuadLatticeElement
    V1_basis # ::Collection{QuadLatticeElement}
    W_reps # ::{QuadLatticeElement}
    v0vec_times_G::Adjoint{Int,Array{Int,1}} # used to hopefully make computations of the form v₀⊙something faster
    v0norm::Int # the norm of v₀ computed once and for all
end

function times_v0(VL::VinbergLattice,e::QuadLatticeElement)
    #@assert e.L == VL.L

    return VL.v0vec_times_G * e.vec
end

function Base.:(==)(L1::VinbergLattice,L2::VinbergLattice)
    @assert false "Let's not compare vinberg lattices yet"
end


function VinbergLattice(G::Array{Int,2};v0vec::Union{Array{Int,1},Nothing}=nothing)
   
    assert_sig_n_1_matrix(G)
    

    L = QuadLattice(G)
    v0 = (v0vec == nothing ? negative_vector(L) : QuadLatticeElement(L,v0vec))
    
    @assert norm(v0) < 0

    V1_basis = basis_of_orhogonal_complement(L,v0)

    @info "V1_basis is $([v.vec for v in V1_basis])"

    M = Matrix(reshape(v0.vec,(rank(L),1)))
    for v in V1_basis
        M = hcat(M,Matrix(reshape(v.vec,(rank(L),1))))
    end
    W = get_integer_points(M)
    W_reps = (x -> QuadLatticeElement(L,x)).(W)
   
    @assert length(W_reps) == abs(det(M)) "only $(length(W_reps)) but need $(det(M))"
    @assert norm(v0) % length(W_reps) == 0
   
    v0vec_times_G = v0.vec' *  G
    println("v0voc_times_G = $v0vec_times_G")
    v0norm = v0⊙v0


    return VinbergLattice(L,v0,V1_basis,W_reps,v0vec_times_G,v0norm)   
end

function Base.convert(v::QuadLatticeElement) 
    return v.vec::Array{Int,1} 
end

same_quad_lattice(v::QuadLatticeElement,w::QuadLatticeElement) = v.L == w.L

function inner_product(v::QuadLatticeElement,w::QuadLatticeElement)
    #@assert same_quad_lattice(v,w) "The elements must belong to the same lattice"
    G = v.L.G
    vv = v.vec
    ww = w.vec
    return vv⋅(G*ww)
end

⊙(v::QuadLatticeElement,w::QuadLatticeElement) = inner_product(v,w)

function norm(v::QuadLatticeElement)
    return inner_product(v,v)
end

function standard_basis(L::QuadLattice)
    I = LinearAlgebra.Diagonal(ones((rank(L))))
    basis = []
    for i in 1:rank(L)
        push!(basis,QuadLatticeElement(L,I[:,i])) 
    end
    return basis
end

function Base.:-(v::QuadLatticeElement)         ## the `Base.:` syntax is necessary to make julia undenstand that we want to extend the function `-` to our type
    return QuadLatticeElement(-v.vec,v.L)
end 

function Base.:+(v::QuadLatticeElement,w::QuadLatticeElement)
    #@assert same_quad_lattice(v,w) "The elements must belong to the same lattice"
    return QuadLatticeElement(v.L,v.vec+w.vec)
end 

function Base.:*(k::Int, v::QuadLatticeElement)
    return QuadLatticeElement(v.L,k*v.vec)
end 

function is_root(v::QuadLatticeElement)

    nv = norm(v)
    # A root has (non strictly) positive norm
    if nv < 0
        return false
    end

    # A root is a primitive vector, i.e. the gcd of its entries is 1
    vv = v.vec
    if abs(gcd(vv)) ≠ 1
        return false
    end

    # A root has length dividing twice the last invariant factor (see B&P) but 
    # this condition is actually not by definition, so we could skip it and get the same result afaik
    if nv ∉ v.L.root_lengths
        return false
    end

    # A root respects the crystallographic condition
    #if ! all(2*(e⊙v) % (v⊙v) == 0 for e in standard_basis(v.L))
    if ! all(2*(col⋅v.vec) % nv == 0 for col in eachcol(v.L.G))
        return false
    end

    return true

end

function zero_elem(L::QuadLattice)
    QuadLatticeElement(L,zeros(rank(L)))
end

function zero_elem(VL::VinbergLattice)
    zero_elem(VL.L) 
end

function reflection(r::QuadLatticeElement,v::QuadLatticeElement)
   
    #@assert same_quad_lattice(v,r) "The elements must belong to the same lattice"
    #@assert is_root(r) "r needs to be a root"

    return v - (2*(r⊙v)/(r⊙r))*r

end


function vec(v::QuadLatticeElement)
    return v.vec
end

# v₀ a vector of the lattice L
# Let V₁ the orthogonal subspace
# x ∈ V₁ iff x'*G*v₀  == 0, i.e. x⊙v₀ == 0
#
# WARNING: I'm not sure the right_kernel really gives us a basis…
#
function basis_of_orhogonal_complement(L::QuadLattice, v0::QuadLatticeElement)
    @assert v0.L == L

    S = MatrixSpace(ZZ, 1, rank(L))
    Mv0 = L.G * v0.vec
    MMv0 = S(reduce(vcat,Mv0))
    right_ker_rank, right_ker = right_kernel(MMv0)
    @assert right_ker_rank == rank(L) - 1 "We need a basis!"


    span = [QuadLatticeElement(L,u) for u in eachcol(Matrix(right_ker))]
    return span

end


# The distance between the hyperplane H_{e} and the point v0
function sinh_distance_to_hyperplane(v0::QuadLatticeElement, e::QuadLatticeElement)
    
    #@assert same_quad_lattice(v0,e) "The elements must belong to the same lattice"
    #@assert is_root(e)
    
    return sqrt( - (e⊙v0)^2 / ( norm(e)*norm(v0) ) ) 
end 
# The distance between the hyperplane H_{e} and the point v0
function distance_to_hyperplane(v0::QuadLatticeElement, e::QuadLatticeElement)
    
    return asinh( sinh_distance_to_hyperplane(v0,e) )

end 


function signature(G)
    
    D,T = LinearAlgebra.eigen(G)
    

    pos = 0
    neg = 0
    for d in D
        if d < 0
           neg = neg + 1
        elseif d > 0
           pos = pos + 1
        end
    end
    return (pos,neg)

end


function assert_sig_n_1_matrix(G)
    
    @assert length(size(G)) == 2            "G must be a matrix"
    @assert size(G)[1] == size(G)[2]        "G must be square"
    @assert LinearAlgebra.issymmetric(G)    "G must be symmetric"
    np = size(G)[1]
    @assert (np-1,1) == signature(G)        "G must have signature (n,1)" 
end



function negative_vector(L::QuadLattice)

    assert_sig_n_1_matrix(L.G)  
    
    # np = n+1
    np = rank(L)
    # Find P and D diagonal such that P'GP = D is diagonal
    G,D,P = L.G,L.D,L.P 
    @assert P'*G*P == Diagonal(L.D)

    # D Necessarily has one negative entry since of signature (n,1)
    # Find the element corresponding to this negative entry
    v0 = [P'[i,:] for i in eachindex(D) if D[i]<0][1]
    # normalize it 
    v0 = (1//gcd(v0))*v0
    
    v0 = QuadLatticeElement(L,v0)
    

    # verify that it indeed has negative norm
    @assert norm(v0) < 0 "v₀ must have negative norm"

    return v0

end

function roots_in_V1(VL::VinbergLattice)

    sort([v for k in root_lengths(VL.L) for v in roots_decomposed_into(VL,zero_elem(VL.L),k) if is_root(v)],by=(v->v.vec))

end

function is_necessary_hyperplane(rr::Vector{QuadLatticeElement},r::QuadLatticeElement)
    
    L = r.L
    #@assert all(ro.L == r.L for ro in rr)
    #@assert is_root(r)
    #@assert all(is_root(ro) for ro in rr)

    return is_necessary_hyperplane([ro.vec for ro in rr],Array{Int,2}(r.L.G),r.vec)

end

function drop_redundant_roots(roots::Vector{QuadLatticeElement})

    for i in 1:length(roots)
        r = roots[i]
        rr = vcat(roots[1:i-1], roots[i+1:end])
        

        if ! is_necessary_hyperplane(rr,r)
            return drop_redundant_roots(rr) 
        end
    end
    
    return roots
    
end

function roots_of_fundamental_cone(VL::VinbergLattice)




    possible_roots = roots_in_V1(VL) 
    @info "possible_roots are: $([r.vec for r in possible_roots])"

    cone_roots::Vector{QuadLatticeElement} = Vector()
    for r in possible_roots
        #println("looking at $(r.vec)")
        #println("cone_roots is $([r.vec for r in cone_roots])")
        #@assert ((-1)*r).vec == -r.vec
        #@assert QuadLatticeElement(r.L,-r.vec) == (-1)*r
        # @assert (r ∉ cone_roots) ⊻ all((-1)*r ≠ cr for cr in cone_roots) "yes?" TODO WTF
        if  all((-1)*r ≠ cr for cr in cone_roots) # &&  all(r⊙cr ≤ 0 for cr in cone_roots) # TODO seems like we can't test equality for our own type QuadLatticeElement
            #println("compatible with other roots")
            # TODO I added the condition r⊙cr≤0 but I'm not sure it's good to have it there
            if is_necessary_hyperplane(cone_roots, r) && is_necessary_hyperplane(cone_roots, (-1)*r)
                println("adding $(r.vec)")
                push!(cone_roots,r)
                cone_roots = drop_redundant_roots(cone_roots)
            end
        
        end
        
    end
    
    #return cone_roots
    return drop_redundant_roots(cone_roots)

end

mutable struct RootsByDistance
    VL::VinbergLattice
    next_least_a0_for_k_and_w::Union{Nothing,Dict{Tuple{Int,QuadLatticeElement},Tuple{Int,Int}}}
    current_a0_and_k_and_w::Union{Nothing,Tuple{Int,Int,QuadLatticeElement}}
    roots_for_current_a0_and_k_and_w::Vector{QuadLatticeElement}
    
end

function RootsByDistance(VL::VinbergLattice;no_distance_zero=false)
    
    @info "> RootsByDistance(…)"
    v0 = VL.v0

    W_reps = VL.W_reps
    # TODO Figure out how to throw away roots at distance zero from the get-go

    next_least_a0_for_k_and_w::Dict{Tuple{Int,QuadLatticeElement},Tuple{Int,Int}} = Dict()
        
    for w in W_reps, k in root_lengths(VL.L)
        least::Rational{Int} = times_v0(VL,w)//VL.v0norm # (w⊙v0)/(v0⊙v0)
        least_plus::Int = ceil(least)
        least_minus::Int = floor(least)
        push!(next_least_a0_for_k_and_w, (k,w) => (least_plus,least_minus))
    end

    @info "W_reps has cardinality $(length(VL.W_reps)) and the number of possible root lengths is $(length(root_lengths(VL.L)))"

    current_a0_and_k_and_w::Union{Nothing,Tuple{Int,Int,QuadLatticeElement}} = nothing;
    #current_a0_and_k_and_w = nothing;
    roots_for_current_a0_and_k_and_w::Vector{QuadLatticeElement} = Vector{QuadLatticeElement}()
    
    @info "< RootsByDistance(…)"
    
    return RootsByDistance(VL,next_least_a0_for_k_and_w,current_a0_and_k_and_w,roots_for_current_a0_and_k_and_w)

end



function next!(r::RootsByDistance)
    
    #@info "> next!(roots_by_distance)"
    v0 = r.VL.v0
    VL = r.VL

    #dist(a0,k,w) = abs(a0*(v0⊙v0) + (w⊙v0))/sqrt(-(k*(v0⊙v0))) # TODO: ensure that the minus sign is needed
    dist(a0,k,w) = abs(a0*VL.v0norm + (times_v0(VL,w)))/sqrt(-(k*VL.v0norm)) # TODO: ensure that the minus sign is needed

    while r.current_a0_and_k_and_w === nothing || isempty(r.roots_for_current_a0_and_k_and_w)
        
        min_val = nothing
        min_tuple = nothing
        for ((k,w),(a0plus,a0minus)) in r.next_least_a0_for_k_and_w 
            if min_tuple === nothing || dist(a0plus,k,w) < min_val
                min_tuple = (a0plus,k,w)
                min_val = dist(a0plus,k,w)
            elseif min_tuple === nothing || dist(a0minus,k,w) < min_val
                min_tuple = (a0minus,k,w)
                min_val = dist(a0minus,k,w)
            end
        end
        
        # we got the triplet k,w,a0 minimizing distance
        r.current_a0_and_k_and_w = min_tuple
         
        a0 = r.current_a0_and_k_and_w[1]
        k = r.current_a0_and_k_and_w[2]
        w = r.current_a0_and_k_and_w[3]

        # we update the dictionary
        (a0plus,a0minus) = r.next_least_a0_for_k_and_w[(k,w)]
        if a0plus == a0
            a0plus += 1
        end
        if a0minus == a0
            a0minus -= 1
        end
        r.next_least_a0_for_k_and_w[(k,w)] = (a0plus,a0minus)
        r.roots_for_current_a0_and_k_and_w = filter(is_root,roots_decomposed_into(r.VL,a0*v0 + w,k))
    end
 
    #@info "< next!(roots_by_distance)"
    return pop!(r.roots_for_current_a0_and_k_and_w)

end



function roots_decomposed_into(VL::VinbergLattice, a::QuadLatticeElement, k::Int)
    # cf eponimous function in B&P's code

    # We are looking for a root ``v = a + v₁``
    # satisfying 
    # * ``v⊙v = k``
    # * ``v₁∈V₁``
    #
    # (v₁+a)⊙(v₁+a) = v₁⊙v₁ + 2 v₁⊙a + a⊙a
    # so
    # (v₁+a)⊙(v₁+a) = k iff  v₁⊙v₁ + 2 v₁⊙a = k-a⊙a,
    # and
   
    #@info "> roots_decomposed_into(VL, $(a.vec), $k)"
    #println("…  $((-(a⊙VL.v0))/(k^0.5))")

    V1Mat::Array{Int,2} = reduce(hcat,[v.vec for v in VL.V1_basis])
   

    #solutions = qsolve_naive(BigInt.(V1Mat' * VL.L.G * V1Mat), BigInt.(V1Mat' * VL.L.G * a.vec), BigInt(a⊙a - k))
    solutions = qsolve_iterative(V1Mat' * VL.L.G * V1Mat, 2 *( V1Mat' * VL.L.G * a.vec), a⊙a - k)
    #println("qsolve ", V1Mat' * VL.L.G * V1Mat, " , ", 2 *(V1Mat' * VL.L.G * a.vec), " , " ,  a⊙a - k)
    #println(solutions)
    solutions_in_L::Array{QuadLatticeElement,1} = (x -> QuadLatticeElement(VL.L,V1Mat * x + a.vec)).(solutions)
     
    #@info "< roots_decomposed_into(VL, $(a.vec), $k)"
    return solutions_in_L
end

function enumerate_roots(VL;num=200,no_distance_zero=false)
    
    roots_by_distance = RootsByDistance(VL;no_distance_zero=no_distance_zero)

    while num > 0
        r = next!(roots_by_distance)
        println("$(r.vec) : $((sinh_distance_to_hyperplane(VL.v0,r))^2)")
        num -= 1
    end
end


function is_finite_volume(roots::Array{QuadLatticeElement,(1)},VL::VinbergLattice)
    Gram_matrix = Int.(reduce(hcat,[[r1⊙r2 for r1 in roots] for r2 in roots]))
    #println("gram matrix is ")
    #display(Gram_matrix)
    Coxeter_matrix = Gram_to_Coxeter(Gram_matrix)
    #println("Coxeter matrix is")
    #display(Coxeter_matrix)
    #println()
    #println("-----------------------------------")
    return isnothing(Coxeter_matrix) ? false : is_fin_vol(Coxeter_matrix,rank(VL.L)-1)
    
end

function Vinberg_Algorithm(G;v0vec=nothing,num_remaining_rounds=100)
        

    VL = VinbergLattice(G;v0vec=v0vec)
    
    v0 = VL.v0 
    

    println("v0 is $v0")

    roots::Array{QuadLatticeElement,(1)} = []
    partial_times = []
    #append!(roots, roots_of_fundamental_cone(VL))
    append!(roots, roots_of_fundamental_cone(VL))
    partial_times = [r.vec' * G for r in roots]

    println("cone roots are:")
    for r in roots
        println(r.vec)
    end
    println("----------------")
    
    

    new_roots_iterator = RootsByDistance(VL)

    

    start = true


    while start || ( ! is_finite_volume(roots,VL) && num_remaining_rounds > 0)
        start = false
        num_remaining_rounds -= 1

        new_root = next!(new_roots_iterator)
        while new_root⊙v0 == 0 # skip roots at distance zero since they have been covered in FundPoly
            new_root = next!(new_roots_iterator)
        end

        #println("($(length(roots)))trying $(new_root.vec)             [$(distance_to_hyperplane(v0,new_root))]")
        
        #while !(all((new_root⊙r) ≤ 0 for r in roots) && times_v0(VL,new_root) < 0) && num_remaining_rounds > 0
        while !(all(pt * new_root.vec ≤ 0 for pt in partial_times) && times_v0(VL,new_root) < 0) && num_remaining_rounds > 0
            num_remaining_rounds -= 1
            new_root = next!(new_roots_iterator)
            #println("($(length(roots)))trying $(new_root.vec)            [$(distance_to_hyperplane(v0,new_root))]")
        end

        #println("new root : $(new_root.vec)")
        push!(roots,new_root)
        push!(partial_times,new_root.vec' * G)
        #println("now have : $([r.vec for r in roots])")

    end
   
    println("decision?")
    println(is_finite_volume(roots,VL))
    println([r.vec for r in roots])

    println("can we drop hyperplanes? $(length(roots))")
    roots = drop_redundant_roots(roots)
    println(length(roots))
    println("remaining rounds: $num_remaining_rounds")
    return [r.vec for r in roots]

end


G0 = [-1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1;
     ]

G1= [-10 0  0 0; 
      0  2  -1 0; 
      0 -1  2 0; 
      0  0   0 1]
G2 = [-7 0   0 0; 
      0 2 -1 0; 
      0 -1 2 -1; 
      0 0 -1 2]


# Example 6.2.5 in Guglielmetti's thesis
Gug_625 = [-1 0 0 0;
            0 2 0 0;
            0 0 6 0;
            0 0 0 6]

include("Some_Lattices.jl")
function test_some_lattices()
    for (name,matrix,basepoint,roots,rounds) in Lattice_table
        println("Checking $name")
        @assert Vinberg_Algorithm(matrix;v0vec=basepoint,num_remaining_rounds=rounds) == roots 
    end
end



