using LinearAlgebra

include("util.jl")


"""
    save_lattice(lat,name)

Save the lattice `lat` in the file `lattices/name.lat`.
"""
function save_lattice(lat, name)
    
    assert_sig_n_1_matrix(lat)
    @assert match(r"^[A-Za-z0-9_-][A-Za-z0-9_-]*$",name) ≠ nothing "The name must consists only on alphanumerical characters, - and _."
    @assert name ∉ keys(L) "The name must not already exist."
    
    path = "lattices/" * name * ".lat" 
    s = open(path,"w") do file
        print(file, matrix_to_txtmat(lat))
    end
end

"""
    ⊕(M₁,M₂)

Compute the "block diagonal" matrix with blocks ``M₁`` and ``M₂``, assuming both are square symetric.
"""
function ⊕(M₁,M₂)
    
    n₁ = size(M₁)[1]
    n₂ = size(M₂)[1]
    @assert size(M₁) == (n₁,n₁)
    @assert size(M₂) == (n₂,n₂)
    @assert M₁' == M₁
    @assert M₂' == M₂

    P = zeros(Int,n₁+n₂,n₁+n₂)
    P[1:n₁,1:n₁] = M₁
    P[n₁+1:end,n₁+1:end] = M₂

    return P
end

# We now define some useful matrixes to build up lattices.
U_ = [0 1; 1 0]
A_(n) = Matrix(SymTridiagonal([2 for i in 1:n], [1 for i in 1:n-1]))
I_ = I(1)

# Using Vinberg's L(n) naming convention, we get, for instance:
Vin07_L16 = U_ ⊕ [2 1;1 4] 
Vin07_L19 = -8*I_ ⊕ A_(3)

# From B&P
BP18_12 = U_ ⊕ 36*I_ ⊕ 6*I_

# This one doesn't seem to be reflective, but it may well be
G2 = [-7 0   0 0; 
      0 2 -1 0; 
      0 -1 2 -1; 
      0 0 -1 2]

# Load all lattices in the `lattices/` folder
L = Dict()
for (root, dirs, files) in walkdir("lattices/")
    for path in files
        if  endswith(path,".lat")
            push!(L,path[1:end-4]=> txtmat_path_to_matrix("lattices/"*path))
        end
    end
end
   

  
