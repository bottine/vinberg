module Diagrams

    import Base: zero, one, +, -, *, /, ^, sqrt, inv, iszero, ==, show, literal_pow

    export is_compact_respectively_finvol 

    include("sbitset.jl")
    include("degree_sequence.jl")
    include("diagrams.jl")

end    
