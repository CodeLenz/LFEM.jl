#
# Sequencia de chamadas para a solução do problema de equilíbrio
#
function Solve_KU(mesh::Mesh; x=Float64[], p=1.0)
  
    # Assembly
    K = Global_K(mesh;x=x,p=p)
    F = Point_load(mesh)

    # Free dofs
    free_dofs = mesh.free_dofs
    
    # Solve just for free dofs
    Chol = cholesky(K[free_dofs,free_dofs])
    Ul = Chol\F[free_dofs]
    
    # Expand homogeneous ebc
    Us  = zeros(length(F))
    Us[free_dofs] .= Ul
    
    return Us, F, Chol
    
 end
