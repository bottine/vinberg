using Pkg
for pkg in ["ResumableFunctions", "AbstractAlgebra", "BitIntegers", "CSV", "Convex", "COSMO", "Cbc", "DataStructures", "LinearAlgebra", "Logging", "MathOptInterface", "Memoize", "StaticArrays", "ToggleableAsserts"]
    Pkg.add(pkg)
end
