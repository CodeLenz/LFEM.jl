
 
#
# f!(t,F,m) must be a function of t, mesh and F where F is dim*nn 
#
# function f!(t,F,m::Mesh)
#          P = Point_load(m)
#          F.= cos(2*t)*P 
# end   
#
#
# gls must be a matrix with [node gl ; node gl ...] to monitor
#
#
#
function Solve_newmark(mesh::Mesh, f!::Function, gls::Matrix{Int64}, ts::Tuple(Float64,Float64), Δt::Float64;
                 x=Float64[], U0=Float64[], V0=Float64[], β=1/4, γ=1/2, p=1.0)


    #
    #                              Initialization
    #             

    # Tspan
    
    t0,tf = ts
    tspan = t0:Δt:tf

    # Number of time steps
    nt = length(tspan)

    # Alias
    nn = mesh.bmesh.nn

    # free dofs
    free_dofs = mesh.free_dofs
    nfree = mesh.ngls

    # Dimension
    dim = 2
    if isa(mesh,Mesh3D)
        dim=3
    end

    # Dimension of full vectors
    nfull = dim*nn

    # Check x
    if isempty(x)
        x = ones(mesh.bmesh.ne)
    end

    # Check initial conditons
    if isempty(U0)
        resize!(U0,nfull)
        fill!(U0,0.0)
    else
        length(U0)==dim*nn || throw("Newmark::U0::wrong dimension")
    end

    if isempty(V0)
        resize!(V0,nfull)
        fill!(V0,0.0)
    else
        length(V0)==dim*nn || throw("Newmark::V0::wrong dimension")
    end

    # Force vectors
    F = zeros(nfull)

    # List with DOFs to monitor
    if isempty(gls)
        error("At least one valid pair of node/dof must be monitored (gls)")
    else

        # Lets make it easy to use it in this routine
        ndofs = size(gls,1)
        println("Newmark::monitoring $ndofs DOFs in $nt time steps")
        dofs = zeros(Int64,ndofs)
        for i=1:ndofs
            dofs[i] = dim*(gls[i,1]-1)+gls[i,2]
        end
    end

    #
    #                           Newmark L.H.S
    #

    # Global stiffness matrix
    K = Global_K(mesh, x=x, p=p)

    # Global mass matrix
    M = Global_M(mesh, x=x)

    # Just to play a little bit
    C = 1E-6*K

    # Newmark operator
    MN = M .+ β*K*Δt^2 .+ γ*C*Δt

    # Cholesky decomposition
    CMN = cholesky(Symmetric(MN[free_dofs,free_dofs]))

    #
    # Initial acceleration
    #

    # Initial force
    f!(t0,F,mesh)

    # Lets make a final consistency test
    @assert length(F)==nfull "Function f!(t,F) must return a $nfull length vector F"

    rhs = F .- K*U0 .- C*V0
    A0f = M[free_dofs,free_dofs]\rhs[free_dofs]

    # Expand A0f 
    A0 = Expand_vector(A0f,nfull,free_dofs)

    # Arrays to monitor the solution
    A_t = Vector{Float64}(undef,nt+1)
    A_U = Array{Float64}(undef,nt+1,ndofs)
    A_V = Array{Float64}(undef,nt+1,ndofs)
    A_A = Array{Float64}(undef,nt+1,ndofs)

    # Store initial values
    A_t[1]    = t0
    A_U[1,:] .= U0[dofs]
    A_V[1,:] .= V0[dofs]
    A_A[1,:] .= A0[dofs]

    # Main Loop
    count = 2
    for t in tspan
        
        # Force in the next time step
        f!(t+Δt,F,mesh)  

        # R.H.S in t+dt
        b = F .- K*U0 .-(C .+Δt*K)*V0 .- (C*Δt*(1-γ) .+ K*(1/2-β)*Δt^2)*A0

        # Solve for A
        Af = CMN\b[free_dofs]

        # Expand A0f 
        A = Expand_vector(Af,nfull,free_dofs)

        # Velocity and displacement at t+Δt
        V = V0 .+ Δt*( (1-γ)*A0 .+ γ*A )
        U = U0 .+ Δt*V0 .+ ( (1/2-β)*A0 .+ β*A )*Δt^2

        # translation
        A0 .= A
        V0 .= V
        U0 .= U

        A_t[count]    = t+Δt
        A_U[count,:] .= U0[dofs]
        A_V[count,:] .= V0[dofs]
        A_A[count,:] .= A0[dofs]
        count += 1

    end #t

    # Retorna os arrays de monitoramento
    return A_U, A_V, A_A, A_t

end 


