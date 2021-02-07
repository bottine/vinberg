"""
    is_finite_volume(matrix,dimension)

* `matrix` is the Coxeter matrix 
* The dimension of ambient space
"""
def is_finite_volume(matrix,dimension):
    
    # Need to first install PyJulia following the readme here https://github.com/JuliaPy/pyjulia
    from julia import Main
    
    path_to_diagrams_dot_jl = "./diagrams.jl"
    Main.include(path_to_diagrams_dot_jl)
    
    Main.das = Main.build_diagram_and_subs(matrix,dimension)

    return Main.is_finite_volume(Main.das)


