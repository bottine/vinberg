""" 
    DiagramType

All possible types of **irreducible** spherical or affine [Coxeter diagram](https://en.wikipedia.org/wiki/Coxeter%E2%80%93Dynkin_diagram).
Lowercase represents spherical diagrams, and uppercase represents affine diagrams. 
"""
@enum DiagramType begin
    # Spherical non sporadic
    DT_a
    DT_b
    DT_d

    # Spherical sporadic
    DT_e6
    DT_e7
    DT_e8
    DT_f4
    DT_g2
    DT_h2
    DT_h3
    DT_h4
    DT_in

    # Affine non sporadic
    DT_A
    DT_B
    DT_C
    DT_D

    # Affine sporadic
    DT_E6
    DT_E7
    DT_E8
    DT_F4
    DT_G2
    DT_I∞
end

"""
    is_spherical(dt::DiagramType)

Whether the type is spherical.
"""
@inline function is_spherical(dt::DiagramType)::Bool
    #dt ∈ [DT_a,DT_b,DT_d,DT_e6,DT_e7,DT_e8,DT_f4,DT_g2,DT_h2,DT_h3,DT_h4,DT_in]
    dt == DT_a || dt == DT_b || dt == DT_d || dt == DT_e6 || dt == DT_e7 || dt == DT_e8 || dt == DT_f4 || dt == DT_g2 || dt == DT_h2 || dt == DT_h3 || dt == DT_h4 || dt == DT_in

end

"""
    is_affine(dt::DiagramType)

Whether the type is affine.
"""
@inline function is_affine(dt::DiagramType)::Bool
    #dt ∈ [DT_A,DT_B,DT_C,DT_D,DT_E6,DT_E7,DT_E8,DT_F4,DT_G2,DT_I∞]
    dt == DT_A || dt == DT_B || dt == DT_C || dt == DT_D || dt == DT_E6 || dt == DT_E7 || dt == DT_E8 || dt == DT_F4 || dt == DT_G2 || dt == DT_I∞
end

"""
    is_sporadic(dt::DiagramType)

Whether the type is sporadic.
"""
@inline function is_sporadic(dt::DiagramType)::Bool
    #dt ∈ [DT_e6,DT_e7,DT_e8,DT_f4,DT_g2,DT_h2,DT_h3,DT_h4,DT_in,DT_E6,DT_E7,DT_E8,DT_F4,DT_G2,DT_I∞]
    dt == DT_e6 || dt == DT_e7 || dt == DT_e8 || dt == DT_f4 || dt == DT_g2 || dt == DT_h2 || dt == DT_h3 || dt == DT_h4 || dt == DT_in || dt == DT_E6 || dt == DT_E7 || dt == DT_E8 || dt == DT_F4 || dt == DT_G2 ||dt == DT_I∞
end
