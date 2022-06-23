#
# Sequencia de chamadas para a solução do problema modal
#
# nev -> number of modes
# which=:SM -> Smaller in magnitude
function Solve_Modal(mesh::Mesh, nev=4, which=:SM; x=Float64[], p=1.0)
  
    # Assembly K and M
    K = Global_K(mesh;x=x,p=p)
    M = Global_M(mesh;x=x)

    # Free dofs
    free_dofs = mesh.free_dofs
    
    # Local views to the free dofs
    KV = @view Symmetric(K[free_dofs, free_dofs])
    MV = @view Symmetric(M[free_dofs, free_dofs])

    # Solve using Arpack
    λ, ϕ = eigs(KV,MV,nev=nev,which=which)

    # Expand the modes to the full mesh
    modes = zero(ϕ)
    for j=1:nev
      modes[free_dofs,j] .= ϕ[:,j]
    end

    # Return the eigenvalues and the eigenvectors
    return λ, modes
    
 end