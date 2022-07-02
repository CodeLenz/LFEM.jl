"""
Solve the harmonic problem Kd(x,w)Ud(x,w) = F(w), 
where Kd(x,w)= K(x)-M(x)w^2 + im*w*C(x)

    Solve_harmonic(mesh::Mesh, w::Float64, x::Vector{Float64}, kparam::Function, mparam::Function)

where 

    w is the angular frequency
    x is a ne x 1 vector of design varibles 
    kparam(xe): R->R is the material parametrization for K (SIMP like)
    mparam(xe): R->R is the material parametrization for M (SIMP like)

Returns:

    Ud = displacement vector (ComplexF64) of size dim*nn x 1
    LU = LU factorization of Kd(x,w) (just free positions)
"""
function Solve_harmonic(mesh::Mesh, w::Float64, x::Vector{Float64}, kparam::Function, mparam::Function)
  
    # Basic checks
    w >= 0.0 || throw("Solve_harmonic:: angular frequency w must be >=0.0")
    length(x)==Get_ne(mesh) || throw("Solve_harmonic:: length of x must be ne")

    # Assembly K and M
    K = Global_K(mesh,x,kparam)
    M = Global_M(mesh,x,mparam)

    # total size
    nfull = size(K,1)
  
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
    Ud = Expand_vector(Ul,nfull,free_dofs)
    
    # Return
    return Ud, LU
    
 end

 """
Solve the harmonic problem Kd(w)Ud(w) = F(w), 
where Kd(w)= K-Mw^2 + im*w*C

    Solve_harmonic(mesh::Mesh, w::Float64)

where 

    w is the angular frequency
Returns:

    Ud = displacement vector (ComplexF64) of size dim*nn x 1
    LU = LU factorization of Kd(x,w) (just free positions)
"""
function Solve_harmonic(mesh::Mesh, w::Float64)

      # x->1.0 mapping
      dummy_f(x)=1.0

      # x is not used
      x = Vector{Float64}(undef,Get_ne(mesh))
  
     # Call Solve_harmonic
     Solve_harmonic(mesh, w, x, dummy_f, dummy_f)

end
