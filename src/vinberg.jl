using Logging

include("vinberg_algorithm.jl")
include("lattice_playground.jl")

# Code adapted from N. V. Bogachev and A. Yu. Perepechko:
#
#   https://github.com/aperep/vinberg-algorithm



function test_some_lattices()
    for (name, matrix, basepoint, roots, rounds) in Lattice_table
        println("Checking $name")
        @assert Vinberg_Algorithm(matrix;v₀_vec=basepoint,rounds=rounds) == roots "$name failed!" 
    end
end

function write_all_lattices()
    for (name, matrix, basepoint, roots, rounds) in Lattice_table
        println("Checking $name")
        matrix_to_txtmat_path(matrix, "lattices/$name")
    end
end

function Vinberg_Algorithm_stats(
    folder, 
    path;
    v₀_vec=nothing
)
    G = txtmat_path_to_matrix(joinpath(folder, path))
    roots, time = @timed Vinberg_Algorithm(VinbergLattice(G, v₀_vec=v₀_vec))
    v₀_vec = VinbergLattice(G, v₀_vec=v₀_vec).v₀.vec
    # return JSON.json(Dict("path"=>path,"matrix"=>G,"v₀"=>v₀_vec,"roots"=>roots,"time"=>time))
    return Dict("path" => path, "matrix" => G, "v₀" => v₀_vec, "roots" => roots, "time" => time)
end


function all_json_output()
    outs = []
    for (name, matrix, basepoint, roots, rounds) in Lattice_table
        push!(outs, "\t" * Vinberg_Algorithm_JSON_output("lattices/" * name, v₀_vec=basepoint))
    end
    
    return "[" * join(outs, ", ") * "]"
end


function test_suite(label=nothing;cache_behaviour=:empty_batched,log_location=nothing)

    if log_location ≠ nothing 
        println("Logging to $log_location")
        logger = SimpleLogger(open(log_location, "w+"))
        global_logger(logger)
    end

    seen = []
    
    if cache_behaviour == :empty_batched
        println("caches are emptied between batched i.e. every time test_suite is run")
    end
    if cache_behaviour == :empty_singles
        println("caches are emptied between every single lattice")
    end

    if cache_behaviour == :empty_batched
        empty!(memoize_cache(is_necessary_halfspace))
    end

    open("lattices/known_values.json", "r") do io
       
        # https://gist.github.com/silgon/0ba43e00e0749cdf4f8d244e67cd9d6a
        all_entries = read(io, String)  # file information to string
        known_values = JSON.parse(all_entries)
        new_values = []

        for entry in known_values
            path = String(entry["path"])
            matrix = Array{Int,2}(hcat(entry["matrix"]...))
            v₀_vec = Array{Int,1}(entry["v0"])
            roots = Array{Array{Int,1},1}(entry["roots"])
            time = Float64(entry["time"])
            time_history = "time_history" ∈ keys(entry) ? Dict{String,Float64}(entry["time_history"]) : Dict{String,Float64}()
            push!(seen, path)
        
            G = txtmat_path_to_matrix("lattices/" * path)
            @assert G == matrix
            println("Looking at $path:")
            
            if cache_behaviour == :empty_singles
                println("emptying caches")
                empty!(memoize_cache(is_necessary_halfspace))
            end 

            my_roots, my_time = @timed Vinberg_Algorithm(VinbergLattice(G, v₀_vec=v₀_vec))
            VL = VinbergLattice(G, v₀_vec=v₀_vec)
            v₀_vec = VL.v₀.vec
            
             
            if label ≠ nothing
                push!(time_history, label => my_time)
            end

            myr = [(fake_dist(VL, VL.L(r)),r) for r in my_roots]
            ofr = [(fake_dist(VL, VL.L(r)),r) for r in roots]

            sort!(myr)
            sort!(ofr)

            if myr ≠ ofr
                println("$(myr) vs $(ofr)")
                display(myr)
                println("vs")
                display(ofr)
                @assert false
            end

            println("Time change ratio (in %):                                   ", round(100 * my_time / time, digits=1))
            if my_time > time
                println("Taking too long! ($my_time vs $time)")
            elseif label ≠ nothing
                entry["time"] = my_time
            end
            entry["time_history"] = time_history
            push!(new_values, entry)
            
            

        end
        
        println("looked at $seen")
        
        for (root, dirs, files) in walkdir("lattices/")
            for path in files
                if endswith(path, ".lat") && path ∉ seen
                    println("One more: $path")
                    push!(new_values, Vinberg_Algorithm_stats("lattices/", path))
                end
            end
        end
        
        println(JSON.json(new_values))

    end
end
