"""
    gug_coxiter_path_to_matrix(path)

Try to load a "`CoxIter`-formatted" diagram description located at `path` and output the resulting matrix (if nothing goes wrong).

See also: [`gug_coxiter_to_matrix`](@ref)
"""
function gug_coxiter_path_to_matrix(path)
    s = open(path) do file
        read(file, String)
    end

    return gug_coxiter_to_matrix(s)
end

"""
    gug_coxiter_to_matrix(descr)

Try to decode a "`CoxIter`-formated" diagram description contained in `descr` and output the resulting matrix (if nothing goes wrong).
Not all options accepted by `CoxIter` are accepted here.
The following formats are accepted:

* First line: ":number_vertices :rank". Each remaining line ":vertex_number :vertex_number :edge_label".
* First line: ":number_vertices :rank". Second line: first the litteral "vertices labels:" then a sequence of identifier separated by spaces (must be as many as number_vertices). All reamining lines ":identifier :identifier :edge_label".
"""
function gug_coxiter_to_matrix(descr)
   
    lines = split(descr,"\n")
    @assert length(lines) ≥ 2 "CoxIter graph description must have non trivial content!"
    
    num_vertices = nothing
    rank = nothing
    m = match(r"^(?<num_vertices>\d\d*)\s\s*(?<rank>\d\d*)", lines[1])
    if m === nothing
        println("can't match first line (we require the rank to be explicitely given):")
        println(lines[1])
        return nothing 
    end

    num_vertices = parse(Int,m[:num_vertices])
    rank = parse(Int,m[:rank])
     
    D = fill(2,num_vertices,num_vertices)
    for i in 1:num_vertices
        D[i,i] = 1
    end
   
    
    m = match(r"^vertices labels:(?<all_labels>(\s+\w+)*)\s*", lines[2])
    if m ≠ nothing # labelled vertice
        vertices = split(m[:all_labels])
        if length(vertices) ≠ num_vertices
            println("not enough vertex labels")
            return nothing
        end
        
        if length(lines) == 2
            println("not enough lines")
            return nothing
        end
        for line in lines[3:end]
            m = match(r"^(?<from>\w+)\s+(?<to>\w+)\s+(?<label>\d+)", line)
            if m === nothing
                continue
            end
            from = indexin([m[:from]],vertices)[1]
            to = indexin([m[:to]],vertices)[1]
            label = parse(Int,m[:label])
            if from ≠ nothing && to ≠ nothing && from ≠ to
                D[from,to] = label
                D[to,from] = label
            end
        end
  

    else # non labelled vertices
        for line in lines[2:end]
            m = match(r"^(?<from>\d\d*)\s\s*(?<to>\d\d*)\s\s*(?<label>\d\d*)", line)
            if m === nothing
                continue
            end
            from = parse(Int,m[:from])
            to = parse(Int,m[:to])
            label = parse(Int,m[:label])
          
            if from ≠ to && from ∈ 1:num_vertices && to ∈ 1:num_vertices
                D[from,to] = label
                D[to,from] = label
            end
        end
    end


    
    return D, rank

end
