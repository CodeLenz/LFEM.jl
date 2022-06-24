#
# Sequencia de chamadas para a solução do problema harmônico
#
function Solve_harmonic(mesh::Mesh, w::Float64 ; x=Float64[], p=1.0)
  
    # Basic check
    w >= 0.0 || throw("Solve_modal:: angular frequency w must be positive")

    # Assembly K and M
    K = Global_K(mesh;x=x,p=p)
    M = Global_M(mesh;x=x)

    # Assembly F
    F = Point_load(mesh)

    # Free dofs
    free_dofs = mesh.free_dofs
    
    # Local views to the free dofs
    KV = @view  K[free_dofs, free_dofs]
    MV = @view  M[free_dofs, free_dofs]

    # Local damping
    CV = 1E-6*KV

    # Harmonic matrix
    KD = KV + w*im*CV - (w^2)*MV

    # Lu decomposition
    LU = lu(KD)
    
    # Harmonic displacement
    Ul = LU\F[free_dofs]

    # Expand 
    Ud = Expand_vector!(Ul,free_dofs)
    
    # Return
    return Ud, LU
    
 end
