include("util.jl")
include("qsolve.jl")
include("diagrams.jl")
include("hyperbolic_lattices.jl")
include("root_enumeration.jl")

# Code adapted from N. V. Bogachev and A. Yu. Perepechko:
#
#   https://github.com/aperep/vinberg-algorithm


# TODO
#
# * See here for memoization: https://github.com/JuliaCollections/Memoize.jl
# * Extensive unit tests!
#


function is_necessary_hyperplane(rr::Vector{HyperbolicLatticeElement},r::HyperbolicLatticeElement)
    
    L = r.L
    #@assert all(ro.L == r.L for ro in rr)
    #@assert is_root(r)
    #@assert all(is_root(ro) for ro in rr)

    return is_necessary_hyperplane(Vector{SVector{rank(L),Int}}([ro.vec for ro in rr]),r.L.G,r.vec)

end

function drop_redundant_roots(roots::Vector{HyperbolicLatticeElement})

    for i in reverse(collect(1:length(roots)))
       
        rr = copy(roots)
        r = popat!(rr,i)
        #r = roots[i]
        #rr = vcat(roots[1:i-1], roots[i+1:end])
        

        if ! is_necessary_hyperplane(rr,r)
            return drop_redundant_roots(rr) 
        end
    end
    
    return roots
    
end

function roots_of_fundamental_cone(VL::VinbergLattice)




    possible_roots = roots_in_V1(VL) 

    cone_roots::Vector{HyperbolicLatticeElement} = Vector()
    for r in possible_roots
        #println("looking at $(r.vec)")
        #println("cone_roots is $([r.vec for r in cone_roots])")
        #@assert ((-1)*r).vec == -r.vec
        #@assert HyperbolicLatticeElement(r.L,-r.vec) == (-1)*r
        # @assert (r ∉ cone_roots) ⊻ all((-1)*r ≠ cr for cr in cone_roots) "yes?" TODO WTF
        if  all((-1)*r ≠ cr for cr in cone_roots) # &&  all(r⊙cr ≤ 0 for cr in cone_roots) # TODO seems like we can't test equality for our own type HyperbolicLatticeElement
            #println("compatible with other roots")
            # TODO I added the condition r⊙cr≤0 but I'm not sure it's good to have it there
            if is_necessary_hyperplane(cone_roots, r) && is_necessary_hyperplane(cone_roots, (-1)*r)
                push!(cone_roots,r)
                cone_roots = drop_redundant_roots(cone_roots)
            end
        
        end
        
    end
    
    #return cone_roots
    return drop_redundant_roots(cone_roots)

end

function is_finite_volume(roots::Array{HyperbolicLatticeElement,(1)},VL::VinbergLattice)
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


function Vinberg_Algorithm(G::Array{Int,2};v0vec::Union{Array{Int,1},Nothing}=nothing,rounds=nothing)
    VL = VinbergLattice(G;v0vec=v0vec)
    return Vinberg_Algorithm(VL;rounds=rounds)
end

function Vinberg_Algorithm(VL::VinbergLattice;rounds=nothing)

   
    @inline decrease(r) = (r === nothing ? nothing : r-1)
    geqzero(r) = (r=== nothing ? true : r ≥ 0)

    v0 = VL.v0
    G = VL.L.G
    
    
    roots::Array{HyperbolicLatticeElement,(1)} = []
    partial_times = []
    #append!(roots, roots_of_fundamental_cone(VL))
    append!(roots, roots_of_fundamental_cone(VL))
    partial_times = [r.vec' * G for r in roots]


    new_roots_iterator = RootsByDistance(VL)

    start = true


    while start || ( ! is_finite_volume(roots,VL) && geqzero(rounds))
        start = false
        

        
        new_root = next!(new_roots_iterator)
        # TODO optimize to start directly at roots at distance > 0 rather than skip them
        while new_root⊙v0 == 0 # skip roots at distance zero since they have been covered in FundPoly
            new_root = next!(new_roots_iterator)
        end

        #println("($(length(roots)))trying $(new_root.vec)             [$(distance_to_hyperplane(v0,new_root))]")
        
        #while !(all((new_root⊙r) ≤ 0 for r in roots) && times_v0(VL,new_root) < 0) && num_remaining_rounds > 0
        while !(all(pt * new_root.vec ≤ 0 for pt in partial_times) && times_v0(VL,new_root) < 0) && geqzero(rounds) 
            
            rounds = decrease(rounds)

            new_root = next!(new_roots_iterator)
            #println("($(length(roots)))trying $(new_root.vec)            [$(distance_to_hyperplane(v0,new_root))]")
        end

        #println("new root : $(new_root.vec)")
        if true # is_necessary_hyperplane(roots,new_root)
            # seems like this is always satisfied?
            #
            push!(roots,new_root)
            push!(partial_times,new_root.vec' * G)
        end
        #println("now have : $([r.vec for r in roots])")

    end
   
    println("Decision ($(rounds)) :", is_finite_volume(roots,VL))
    println("can we drop hyperplanes? $(length(roots)) vs $(length(drop_redundant_roots(roots)))")
    return [r.vec for r in roots]

end


include("Some_Lattices.jl")
function test_some_lattices()
    for (name,matrix,basepoint,roots,rounds) in Lattice_table
        println("Checking $name")
        @assert Vinberg_Algorithm(matrix;v0vec=basepoint,rounds=rounds) == roots "$name failed!" 
    end
end

function write_all_lattices()
    for (name,matrix,basepoint,roots,rounds) in Lattice_table
        println("Checking $name")
        matrix_to_txtmat_path(matrix,"lattices/$name")
    end
end

function Vinberg_Algorithm_stats(folder,path;v0vec=nothing)
    G = txtmat_path_to_matrix(joinpath(folder,path))
    roots,time = @timed Vinberg_Algorithm(VinbergLattice(G,v0vec=v0vec))
    v0vec = VinbergLattice(G,v0vec=v0vec).v0.vec
    #return JSON.json(Dict("path"=>path,"matrix"=>G,"v0"=>v0vec,"roots"=>roots,"time"=>time))
    return Dict("path"=>path,"matrix"=>G,"v0"=>v0vec,"roots"=>roots,"time"=>time)
end


function all_json_output()
    outs = []
    for (name,matrix,basepoint,roots,rounds) in Lattice_table
        push!(outs,"\t" * Vinberg_Algorithm_JSON_output("lattices/"*name,v0vec=basepoint))
    end
    
    return "["*join(outs,", ")*"]"
end


function test_suite(label=nothing;cache_behaviour=:empty_batched)

    seen = []
    
    if cache_behaviour == :empty_batched
        println("caches are emptied between batched i.e. every time test_suite is run")
    end
    if cache_behaviour == :empty_singles
        println("caches are emptied between every single lattice")
    end

    if cache_behaviour == :empty_batched
            empty!(memoize_cache(qsolve_iterative))
            empty!(memoize_cache(is_necessary_hyperplane))
    end

    open("lattices/known_values.json", "r") do io
       
        # https://gist.github.com/silgon/0ba43e00e0749cdf4f8d244e67cd9d6a
        all_entries = read(io,String)  # file information to string
        known_values = JSON.parse(all_entries)
        new_values = []

        for entry in known_values
            path = String(entry["path"])
            matrix = Array{Int,2}(hcat(entry["matrix"]...))
            v0vec = Array{Int,1}(entry["v0"])
            roots = Array{Array{Int,1},1}(entry["roots"])
            time = Float64(entry["time"])
            time_history = "time_history" ∈ keys(entry) ? Dict{String,Float64}(entry["time_history"]) : Dict{String,Float64}()
            push!(seen,path)
        
            G = txtmat_path_to_matrix("lattices/"*path)
            @assert G == matrix
            println("Looking at $path:")
            
            if cache_behaviour == :empty_singles
                println("emptying caches")
                empty!(memoize_cache(qsolve_iterative))
                empty!(memoize_cache(is_necessary_hyperplane))
            end 

            my_roots,my_time = @timed Vinberg_Algorithm(VinbergLattice(G,v0vec=v0vec))
            v0vec = VinbergLattice(G,v0vec=v0vec).v0.vec
            
             
            if label≠nothing
                push!(time_history,label=>my_time)
            end

            if my_roots ≠ roots
                println("ERROR on roots")
                @assert false
            end

            println("Time change ratio (in %):                                   ", round(100*my_time/time,digits=1))
            if my_time > time
                println("Taking too long! ($my_time vs $time)")
            elseif label≠nothing
                entry["time"] = my_time
            end
            entry["time_history"] = time_history
            push!(new_values,entry)
            
            

        end
        
        println("looked at $seen")
        
        for (root, dirs, files) in walkdir("lattices/")
            for path in files
                if  endswith(path,".lat") && path ∉ seen
                    println("One more: $path")
                    push!(new_values, Vinberg_Algorithm_stats("lattices/",path))
                end
            end
        end
        
        println(JSON.json(new_values))

    end
end
