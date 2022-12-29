"""
Solve the modal problem (M(x) - λK(x))ϕ(x) = 0

  Solve_modal(mesh::Mesh, x::Vector{Float64}, kparam::Function; 
              mparam::Function, nev=4, which=:SM,  σ=1.0, loadcase::Int64=1)

where 

    x is a ne x 1 vector of design varibles 
    kparam(xe): R->R is the material parametrization for K (SIMP like)
    mparam(xe): R->R is the material parametrization for M (SIMP like)
    nev is the number of eigenvalues and eigenvectors to compute
    loadcase is the loadcase

Returns:

    λ = eigenvalues vector (nev x 1)
    modes = matrix dim*nn x nev with the eigenvectors
"""
function Solve_modal(mesh::Mesh, x::Vector{Float64}, kparam::Function, 
                     mparam::Function; nev=4, loadcase::Int64=1)
  
    # Basic assertions
    length(x)==Get_ne(mesh) || throw("Solve_modal:: length of x must be ne")
    0<=loadcase<=mesh.nload || throw("Solve_modal:: invalid loadcase")

    # Assembly K and M
    K = Global_K(mesh,x,kparam)
    M = Global_M(mesh,x,mparam)

    # Free dofs
    free_dofs = mesh.free_dofs[loadcase]
    
    # Local views to the free dofs
    K =  K[free_dofs, free_dofs]
    M =  M[free_dofs, free_dofs]

    # Solve using Arnoldi interface
    λ, ϕ = Solve_Eigen_(K,M,nev)

    # Total number of dofs
    dim = Get_dim(mesh)
    nn  = Get_nn(mesh)
    ngls = dim*nn
  
    # Alocate the matrix (full number of dofs)
    PHI = zeros(ngls,nev)

    # Expand the eigenvectores
    Usf  = zeros(ngls)
    for i=1:nev
 
        # Expande esse modo para os gls globais
        fill!(Usf,0.0)
        Expand_vector!(Usf,real.(ϕ[:,i]),free_dofs)
        PHI[:,i] .= Usf

    end
 
    # Return the eigenvalues and the eigenvectors
    return λ, PHI
    
 end

 """
Solve the modal problem (M - λK)ϕ = 0

  Solve_modal(mesh::Mesh ;nev=4, loadcase=1)

where 
    nev is the number of eigenvalues and eigenvectors to compute
    loadcase is the loadcase

Returns:

    λ = eigenvalues vector (nev x 1)
    modes = matrix dim*nn x nev with the eigenvectors
"""
function Solve_modal(mesh::Mesh; nev=4,  loadcase::Int64=1)

    # x->1.0 mapping
    dummy_f(x)=1.0

    # x is not used
    x = Vector{Float64}(undef,Get_ne(mesh))

    # Call Solve_modal
    Solve_modal(mesh, x, dummy_f, dummy_f, nev=nev, loadcase=loadcase)
  
end
  


#
#
#
# Organize the eigenvalues and eigenvector . To be used as a post-processor
# for the Modal Analysis.
#
#
#
function Organize_Eigen(lambda::Vector,phi::Matrix,ngls::Int64,free_dofs::Vector)

    # Convert to real numbers
    lamb_before = real.(lambda)

    # Number of eigenvalues
    nev_before = length(lamb_before)

    # Make sure to get only the positive eigenvalues
    # in crescent order
    n_effective = 0
    lamb_a  = Float64[]
    pos_lamb = Int64[]
    for i=1:nev_before
        if lamb_before[i] > 0.0
            push!(lamb_a,lamb_before[i])
            push!(pos_lamb,i) 
            n_effective += 1
        end
    end

    # Avoid the situation of no positive eigenvalues
    n_effective >=1 || error("Organize_Eigen:: there is no valid positive solution")

    # sort
    ordem = sortperm(lamb_a)
    lamb = lamb_a[ordem]

    # Convert the eigenvectors to real numbers
    phi_real = real.(phi[:,pos_lamb[ordem]])

    # Alocate the matrix (full number of dofs)
    PHI = zeros(ngls,n_effective)

    # Expand the eigenvectores
    Usf  = zeros(ngls)
    for i=1:n_effective

        # Expande esse modo para os gls globais
        fill!(Usf,0.0)
        Expand_vector!(Usf,phi_real[:,i],free_dofs)
        PHI[:,i] .= Usf

    end

    # Return the positive eigenvalues and their eigenvectors
    # but just the nev requested ones
    return lamb, PHI

end

