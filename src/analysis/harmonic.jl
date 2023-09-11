"""
Solve the harmonic problem Kd(x,w)Ud(x,w) = F(w), 
where Kd(x,w)= K(x)-M(x)w^2 + im*w*C(x)

    Solve_harmonic(mesh::Mesh, w::Float64, α_c::Float64, β_c::Float64, x::Vector{Float64}, 
                   kparam::Function, mparam::Function ; 
                   lumped=true,loadcase=1)

where 

    w is the angular frequency
    α_c and  β_c are the parameters for proportional damping
    x is a ne x 1 vector of design varibles 
    kparam(xe): R->R is the material parametrization for K (SIMP like)
    mparam(xe): R->R is the material parametrization for M (SIMP like)
    lumped is true for lumped mass matrices
    loadcase is the loadcase

Returns:

    Ud = displacement vector (ComplexF64) of size dim*nn x 1
    linsolve = LinearSolve object with factored linear problem

"""
function Solve_harmonic(mesh::Mesh, w::Float64, α_c::Float64, β_c::Float64,
                        x::Vector{Float64},
                        kparam::Function, mparam::Function; lumped=true, loadcase::Int64=1)
  
    # Basic checks
    w >= 0.0 || throw("Solve_harmonic:: angular frequency w must be >=0.0")
    length(x)==Get_ne(mesh) || throw("Solve_harmonic:: length of x must be ne")
    0<=loadcase<=mesh.nload || throw("Solve_harmonic:: invalid loadcase")

    # Assembly K, M and C
    K = Global_K(mesh,x,kparam)
    M = Global_M(mesh,x,mparam,lumped=lumped)
    C = Global_C(M,K,mesh,α_c,β_c)

    # total size
    nfull = size(K,1)
  
    # Assembly F
    F = Point_load(mesh,loadcase)

    # Free dofs
    free_dofs = mesh.free_dofs[loadcase]
    
    # Views
    K =  @view K[free_dofs, free_dofs]
    C =  @view C[free_dofs, free_dofs]
    M =  @view M[free_dofs, free_dofs]

    # Harmonic matrix 
    @inbounds KD = sparse(K) .+ (w*im).*sparse(C) .- (w^2).*sparse(M)

    # Create LinearSolve problem
    prob = LinearProblem(KD,complex.(F[free_dofs]))
    linsolve = init(prob)

    # Harmonic displacement
    Ul = solve(linsolve)

    # Expand 
    Ud = Expand_vector(Ul.u,nfull,free_dofs)
    
    # Return
    return Ud, linsolve
    
 end

 """
Solve the harmonic problem Kd(w)Ud(w) = F(w), 
where Kd(w)= K-Mw^2 + im*w*C

    Solve_harmonic(mesh::Mesh, w::Float64, α_c::Float64, β_c::Float64 ; loadcase=1)

where 

    w is the angular frequency  
    α_c and  β_c are the parameters for proportional damping
    lumped=true is for lumped mass matrix
    loadcase is the loadcase  

Returns:

    Ud = displacement vector (ComplexF64) of size dim*nn x 1  
    linsolve = LinearSolve object with factored linear problem
"""
function Solve_harmonic(mesh::Mesh, w::Float64, α_c::Float64, β_c::Float64 ;
                       lumped=true, loadcase::Int64=1)

      # x->1.0 mapping
      dummy_f(x)=1.0

      # x is not used
      x = Vector{Float64}(undef,Get_ne(mesh))
  
     # Call Solve_harmonic
     
     Solve_harmonic(mesh, w, α_c, β_c, x, dummy_f, dummy_f, lumped=lumped, loadcase=loadcase)

end
