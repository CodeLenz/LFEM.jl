#
# Sequencia de chamadas para a solução do problema modal
#
# nev -> number of modes
# which=:SM -> Smaller in magnitude
function Solve_modal(mesh::Mesh, nev=4, which=:SM; x=Float64[], p=1.0)
  
    # Assembly K and M
    K = Global_K(mesh;x=x,p=p)
    M = Global_M(mesh;x=x)

    # Free dofs
    free_dofs = mesh.free_dofs
    
    # Local views to the free dofs
    KV = @view  K[free_dofs, free_dofs]
    MV = @view  M[free_dofs, free_dofs]

    # Solve using Arpack
    λ, ϕ = eigs(Symmetric(KV),Symmetric(MV),nev=nev,which=which)

    # Expand the modes to the full mesh
    dim = 2
    if isa(mesh,Mesh3D)
        dim=3
    end
    modes = zeros(dim*mesh.bmesh.nn,nev)
    for j=1:nev
      modes[free_dofs,j] .= ϕ[:,j]
    end

    # Return the eigenvalues and the eigenvectors
    return λ, modes
    
 end
