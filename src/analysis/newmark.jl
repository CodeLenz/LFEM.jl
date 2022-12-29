"""
Solve the transient problem M(x)A(x,t) + C(x)V(x,t) + K(x,t)U(x,t) = F(t), 
using Newmark-beta method.

    Solve_newmark(mesh::Mesh, f!::Function, gls::Matrix{Int64}, 
                  ts::Tuple{Float64, Float64}, Δt::Float64,
                  x::Vector{Float64}, kparam::Function, mparam::Function,
                  verbose=false;
                  U0=Float64[], V0=Float64[],
                  β=1/4, γ=1/2, 
                  α_c=0.0, β_c=1E-6
                  loadcase=1)

where 

    ts is Tupple with initial and end time (Ti,Tf)  
    Δt is (fixed) time step  
    x is a ne x 1 vector of design varibles   
    kparam(xe): R->R is the material parametrization for K (SIMP like)  
    mparam(xe): R->R is the material parametrization for M (SIMP like)  
    verbose is false or true
    U0 and V0 are the initial conditions  
    β and γ are the parameters of the Newmark method
    α_c and  β_c are the coefficients for proportional damping C=α_cM + β_c*K
    loadcase is the loadcase
 
    f!(t,F,mesh,loadcase) must be a function of t, mesh and F where F is dim*nn x 1,
                Example: 
    
                function f!(t,F,mesh::Mesh,loadcase=1)
                            P = Point_load(mesh,loadcase)
                            F.= cos(2*t)*P 
                end   

    gls is a matrix with [node gl ;
                          node gl ...] to monitor

Return three arrays of size ng x nt, where ng is size(gls,1) and nt is the 
number of time steps (length of t0:Δt:tf)

    A_U displacements
    A_V velocities
    A_A accelerations
    
    A_t is a vector of size nt x 1 with discrete times
    dofs is a vector with the (global) monitored dofs
"""
function Solve_newmark(mesh::Mesh, f!::Function, gls::Matrix{Int64}, 
                      ts::Tuple{Float64, Float64}, Δt::Float64,
                      x::Vector{Float64}, kparam::Function, 
                      mparam::Function, verbose=false; 
                      U0=Float64[], V0=Float64[], β=1/4, γ=1/2,
                      α_c=0.0, β_c=1E-6,
                      loadcase::Int64=1)


    #
    #                              Initialization
    #             

    # Tspan    
    t0,tf = ts
    tspan = t0:Δt:tf-Δt

    # Number of time steps
    nt = length(tspan)

    # Alias
    ne = Get_ne(mesh)
    nn = Get_nn(mesh)

    # Check if loadcase is valid
    0<=loadcase<=mesh.nload || throw("Solve_newmark:: invalid loadcase")

    # free dofs
    free_dofs = mesh.free_dofs[loadcase]
    nfree = mesh.ngls[loadcase]

    # Dimension
    dim = Get_dim(mesh)

    # Dimension of full vectors
    nfull = dim*nn

    # Check x
    length(x)==ne || throw("Solve_newmark:: length of x must be ne")

    # Check initial conditons
    if isempty(U0)
        resize!(U0,nfull)
        fill!(U0,0.0)
    else
        length(U0)==dim*nn || throw("Newmark::U0 should be dim*nn")
    end

    if isempty(V0)
        resize!(V0,nfull)
        fill!(V0,0.0)
    else
        length(V0)==dim*nn || throw("Newmark::V0 should be dim*nn")
    end

    # Force vector
    F = zeros(nfull)

    # List with DOFs to monitor
    if isempty(gls)
        error("At least one valid pair of node/dof must be monitored (gls)")
    else

        # Lets make it easy to use it in this routine
        ndofs = size(gls,1)
        
        verbose && println("Newmark::monitoring $ndofs DOFs in $nt time steps")
        
        dofs = zeros(Int64,ndofs)
        for i=1:ndofs
            dofs[i] = dim*(gls[i,1]-1)+gls[i,2]
        end
    end

    #
    #                           Newmark L.H.S
    #

    # Global stiffness matrix
    K = Global_K(mesh, x, kparam)

    # Global mass matrix
    M = Global_M(mesh, x, mparam)

    # Global damping matrix
    C = Global_C(M,K,mesh,α_c,β_c)

    # Some views
    VK =  K[free_dofs,free_dofs]
    VM =  M[free_dofs,free_dofs]
    VC =  C[free_dofs,free_dofs]

    #
    # Initial acceleration
    #

    # Initial force
    f!(t0,F,mesh,loadcase)

    # Lets make a final consistency test
    @assert length(F)==nfull "Solve_newmark:: Function f!(t,F) must return a $nfull length vector F"

    rhs = F[free_dofs] .- VK*U0[free_dofs] .- VC*V0[free_dofs]
    A0f = VM\rhs

    # Expand A0f 
    A0 = Expand_vector(A0f,nfull,free_dofs)

    # Arrays to monitor the solution. The number of time points is 
    # tspan / dt + 1.
    A_t = Vector{Float64}(undef,nt+1)
    A_U = Array{Float64}(undef,nt+1,ndofs)
    A_V = Array{Float64}(undef,nt+1,ndofs)
    A_A = Array{Float64}(undef,nt+1,ndofs)

    # Store initial values (time zero)
    A_t[1]    = t0
    A_U[1,:] .= U0[dofs]
    A_V[1,:] .= V0[dofs]
    A_A[1,:] .= A0[dofs]

    # Define here to avoid allocations
    U = similar(F)
    V = similar(F)
    A = similar(F)
    b = similar(A0f)
    Af = similar(A0f)

    # Newmark operator
    MN = VM .+ β*VK*Δt^2 .+ γ*VC*Δt

    # Create a LinearPoblem
    prob = LinearProblem(Symmetric(MN),b)
    linsolve = init(prob)

    # Main Loop. At each t in the loop we are at t, evaluating for the next time steps
    # t + Δt.
    count = 2
    for t in tspan
        
        # Force in the next time step
        f!(t+Δt,F,mesh,loadcase)  

        # R.H.S in t+dt
        LinearSolve.set_b(linsolve,F[free_dofs] .- VK*U0[free_dofs] .-(VC .+Δt*VK)*V0[free_dofs] .- (VC*Δt*(1-γ) .+ VK*(1/2-β)*Δt^2)*A0[free_dofs])   

        # Solve for A in t+Δt
        sol = solve(linsolve)
        Af .= sol.u

        # Expand A0f 
        Expand_vector!(A,Af,free_dofs)

        # Velocity and displacement at t+Δt
        V .= V0 .+ Δt*( (1-γ)*A0 .+ γ*A )
        U .= U0 .+ Δt*V0 .+ ( (1/2-β)*A0 .+ β*A )*Δt^2

        # Store values at t+Δt
        A_t[count]    = t + Δt
        A_U[count,:] .= U[dofs]
        A_V[count,:] .= V[dofs]
        A_A[count,:] .= A[dofs]
        count += 1

        # translation
        A0 .= A
        V0 .= V
        U0 .= U


    end #t

    # Retorna os arrays de monitoramento
    return A_U, A_V, A_A, A_t, dofs
end 


"""
Solve the transient problem MA(t) + CV(t) + K(t)U(t) = F(t), 
using Newmark-beta method.

    Solve_newmark(mesh::Mesh, f!::Function, gls::Matrix{Int64}, 
                  ts::Tuple{Float64, Float64}, Δt::Float64,
                  verbose=false;
                  U0=Float64[], V0=Float64[], β=1/4, γ=1/2,loadcase=1)

where 

    ts is Tupple with initial and end time (Ti,Tf)
    Δt is (fixed) time steps
    verbose is false or true
    U0 and V0 are the initial conditions  
    β and γ are the parameters of the Newmar method
    loadcase is the loadcase

    f!(t,F,mesh,loadcase) must be a function of t, mesh and F where F is dim*nn x 1,
                Example: 
    
                function f!(t,F,mesh::Mesh,loadcase=1)
                            P = Point_load(mesh,loadcase)
                            F.= cos(2*t)*P 
                end   

    gls is a matrix with [node gl ;
                          node gl ...] to monitor

Return three arrays of size ng x nt, where ng is size(gls,1) and nt is the 
number of time steps (length of t0:Δt:tf)

    A_U displacements
    A_V velocities
    A_A accelerations
    
    A_t is a vector of size nt x 1 with discrete times
    dofs is a vector with the (global) monitored dofs
"""
function Solve_newmark(mesh::Mesh, f!::Function, gls::Matrix{Int64},
                       ts::Tuple{Float64, Float64}, Δt::Float64,
                       verbose=false;
                       U0=Float64[], V0=Float64[], β=1/4, γ=1/2,
                       loadcase::Int64=1)

      # x->1.0 mapping
      dummy_f(x)=1.0

      # x is not used
      x = Vector{Float64}(undef,Get_ne(mesh))

      # Call Solve_newmark
      Solve_newmark(mesh,f!,gls,ts,Δt,x,dummy_f,dummy_f,verbose,
                    U0=U0,V0=V0,β=β,γ=γ,loadcase=loadcase)
 
end   





"""
Solve the transient problem M(x)A(x,t) + C(x)V(x,t) + K(x,t)U(x,t) = F(t), 
using Newmark-beta method. Matrizes are given and no mesh information is 
used

    Solve_newmark(M::AbstractMatrix,C::AbstractMatrix,K::AbstractMatrix, f!::Function, gls::Matrix{Int64}, 
                  ts::Tuple{Float64, Float64}, Δt::Float64,
                  verbose=false;
                  U0=Float64[], V0=Float64[],
                  β=1/4, γ=1/2)

where 

    M, C and K are given matrices
    ts is Tupple with initial and end time (Ti,Tf)  
    Δt is (fixed) time step  
    verbose is false or true
    U0 and V0 are the initial conditions  
    β and γ are the parameters of the Newmark method
 
    f!(t,F) must be a function of t  and F where F is dim*nn x 1,
                Example: 
    
                function f!(t,F)
                            P = [0.0;10.0]
                            F.= cos(2*t)*P 
                end   

    gls is a vector gls to monitor

Return three arrays of size ng x nt, where ng is size(gls,1) and nt is the 
number of time steps (length of t0:Δt:tf)

    A_U displacements
    A_V velocities
    A_A accelerations
    
    A_t is a vector of size nt x 1 with discrete times
"""
function Solve_newmark(M::AbstractMatrix,C::AbstractMatrix,K::AbstractMatrix, f!::Function, gls::Vector{Int64}, 
                      ts::Tuple{Float64, Float64}, Δt::Float64,
                      verbose=false; 
                      U0=Float64[], V0=Float64[], β=1/4, γ=1/2)


    #
    #                              Initialization
    #             

    # Tspan    
    t0,tf = ts
    tspan = t0:Δt:tf-Δt

    # Number of time steps
    nt = length(tspan)

    # Problem's size
    nfull = size(M,1)

    # Newmark operator
    MN = M .+ β*K*Δt^2 .+ γ*C*Δt

    # Check initial conditons
    if isempty(U0)
        resize!(U0,nfull)
        fill!(U0,0.0)
    else
        length(U0)==nfull || throw("Newmark::U0 should be $nfull")
    end

    if isempty(V0)
        resize!(V0,nfull)
        fill!(V0,0.0)
    else
        length(V0)==nfull || throw("Newmark::V0 should be $nfull")
    end

    # Force vector
    F = zeros(nfull)

    # List with DOFs to monitor
    if isempty(gls)
        error("At least one valid dof  must be monitored (gls)")
    end


    # Initial acceleration
    f!(0.0,F)
    A0 = M\(F- K*U0 -C*V0)

    # Arrays to monitor the solution. The number of time points is 
    # tspan / dt + 1.
    ndofs = length(gls)
    A_t = Vector{Float64}(undef,nt+1)
    A_U = Array{Float64}(undef,nt+1,ndofs)
    A_V = Array{Float64}(undef,nt+1,ndofs)
    A_A = Array{Float64}(undef,nt+1,ndofs)

    # Store initial values (time zero)
    A_t[1]    = t0
    A_U[1,:] .= U0[gls]
    A_V[1,:] .= V0[gls]
    A_A[1,:] .= A0[gls]

    # Define here to avoid allocations inside the loop
    U = similar(F)
    V = similar(F)
    A = similar(F)
    b = similar(F)

    # Main loop
    count = 2
    for t in tspan

            f!(t+Δt,F)
            b .= F .- K*U0 .-(C .+Δt*K)*V0 .- (C*Δt*(1-γ) .+ K*(1/2-β)*Δt^2)*A0
            A .= MN\b
            V .= V0 .+ Δt*( (1-γ)*A0 .+ γ*A )
            U .= U0 .+ Δt*V0 .+ ( (1/2-β)*A0 .+ β*A )*Δt^2
            
           # Store values at t+Δt
            A_t[count]    = t + Δt
            A_U[count,:] .= U[gls]
            A_V[count,:] .= V[gls]
            A_A[count,:] .= A[gls]
            count += 1

            U0 .= U
            V0 .= V
            A0 .= A
    end

    # Return the values 
    return A_U, A_V, A_A, A_t

end


