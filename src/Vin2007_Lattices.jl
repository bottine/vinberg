using LinearAlgebra

function ⊕(M1,M2)
    
    n1 = size(M1)[1]
    n2 = size(M2)[1]
    @assert size(M1) == (n1,n1)
    @assert size(M2) == (n2,n2)
    @assert M1' == M1
    @assert M2' == M2

    P = zeros(Int,n1+n2,n1+n2)
    P[1:n1,1:n1] = M1
    P[n1+1:end,n1+1:end] = M2

    return P
end

U = [0 1; 1 0]
A_(n) = Matrix(I(n))

# Using Vinberg's L(n) naming convention

Vin07_L1 = 4*U ⊕ 2*A_(1)  ⊕ 2*A_(1)
Vin07_L2 = 2*U ⊕ 2*A_(1)  ⊕ 2*A_(1)
Vin07_L3 = 1*U ⊕ 2*A_(1)  ⊕ 2*A_(1)
Vin07_L4 = 3*U ⊕ 2*A_(1)  ⊕ 2*A_(1)
Vin07_L5 = (-4)*A_(1) ⊕ A_(3)
Vin07_L6 = U ⊕ 2*A_(1)  ⊕ 8*A_(1)

Vin07_L15 = U ⊕ [2 1;1 10]
Vin07_L16 = U ⊕ [2 1;1 4] 
Vin07_L17 = U ⊕ [2 1;1 14] 
Vin07_L18 = U ⊕ [2 0;0 4] 
Vin07_L19 = 8*A_(1) ⊕ A_(3)
