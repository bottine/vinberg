include("diagrams.jl")

# input file in first and only command line argument
args = string.(ARGS)

path = string(ARGS[1])

for path in args
println("Processing $path")
    (D,r) = gug_coxiter_path_to_matrix(path)
    das = build_diagram_and_subs(D,r)
    println("Cocompact     ? ", is_compact(das))
    println("Finite Volume ? ", is_finite_volume(das))
end
