"""
Solve the modal problem (M(x) - λK(x))ϕ(x) = 0

  Solve_modal(mesh::Mesh, x::Vector{Float64}, kparam::Function, mparam::Function, nev=4, which=:SM)

where 

    x is a ne x 1 vector of design varibles 
    kparam(xe): R->R is the material parametrization for K (SIMP like)
    mparam(xe): R->R is the material parametrization for M (SIMP like)
    nev is the number of eigenvalues and eigenvectors to compute
    which is the range (:SM is smaller in magnitude, for example)

Returns:

    λ = eigenvalues vector (nev x 1)
    modes = matrix dim*nn x nev with the eigenvectors
"""
function Solve_modal(mesh::Mesh, x::Vector{Float64}, kparam::Function, mparam::Function, nev=4, which=:SM)
  
    # Basic assertion
    length(x)==Get_ne(mesh) || throw("Solve_modal:: length of x must be ne")

    # Assembly K and M
    K = Global_K(mesh,x,kparam)
    M = Global_M(mesh,x,mparam)

    # Free dofs
    free_dofs = mesh.free_dofs
    
    # Local views to the free dofs
    KV = @view  K[free_dofs, free_dofs]
    MV = @view  M[free_dofs, free_dofs]

    # Solve using Arpack
    λ, ϕ = eigs(Symmetric(KV),Symmetric(MV),nev=nev,which=which)

    # Expand the modes to the full mesh
    dim = Get_dim(mesh)
    nn  = Get_nn(mesh)
    modes = zeros(dim*nn,nev)
    @inbounds for j=1:nev
      modes[free_dofs,j] .= ϕ[:,j]
    end

    # Return the eigenvalues and the eigenvectors
    return λ, modes
    
 end

 """
Solve the modal problem (M - λK)ϕ = 0

  Solve_modal(mesh::Mesh,nev=4, which=:SM)

where 
    nev is the number of eigenvalues and eigenvectors to compute
    which is the range (:SM is smaller in magnitude, for example)

Returns:

    λ = eigenvalues vector (nev x 1)
    modes = matrix dim*nn x nev with the eigenvectors
"""
function Solve_modal(mesh::Mesh, nev=4, which=:SM)

    # x->1.0 mapping
    dummy_f(x)=1.0

    # x is not used
    x = Vector{Float64}(undef,Get_ne(mesh))

    # Call Solve_modal
    Solve_modal(mesh, x, dummy_f, dummy_f, nev, which)
  
end
  