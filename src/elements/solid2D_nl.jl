"""
B Matricies for 2D elements
    Bs_solid2D_nl(r::T,s::T,x::Vector{T},y::Vector{T},ue::Vector{T}) where T

    ue is the 8×1 vector of displacements for this element

"""
function Bs_solid2D_nl(r::T,s::T,x::Vector{T},y::Vector{T},ue::Vector{T}) where T

    # Derivates of N (no bubble for now)
    dNrs = dN_solid2D(r,s)[1:4]

    # Jacobian Matrix
    J = Jacobian_solid2D(x,y,dNrs)

    # Determinant
    detJ = det(J)

    # Map derivatives w.r.t (r,s) to (x,y)
    # (1)   ... (4)
    # d/dx      d/dx
    # d/dy      d/dy
    #
    dNxy = J\dNrs

    # B0 matrix - Linear
    B0 = zeros(T,3,8)
    for j=1:4

        # columns
        c1 = 2*(j-1)+1
        c2 = c1+1

        # xx 
        B0[1,c1]=dNxy[1,j]
        # yy
        B0[2,c2]=dNxy[2,j]
        # xy
        B0[3,c1]=dNxy[2,j] 
        B0[3,c2]=dNxy[1,j]

    end

    # Recover u1 and u2
    u1 = ue[1:2:7]
    u2 = ue[2:2:8]

    # Evaluate u1,1 u1,2 u2,1 and u2,2
    # Derivative of component i wrt j
    # Bathe 565
    u11 = dot(dNxy[1,:],u1)
    u12 = dot(dNxy[2,:],u1)
    u21 = dot(dNxy[1,:],u2)
    u22 = dot(dNxy[2,:],u2)

    # B1 - Non linear terms
    B1 = zeros(T,3,8)
    for j=1:4

        # columns
        c1 = 2*(j-1)+1
        c2 = c1+1

        # xx 
        B1[1,c1]=u11*dNxy[1,j]
        B1[1,c2]=u21*dNxy[1,j]
        
        # yy
        B1[2,c1]=u12*dNxy[2,j]
        B1[2,c2]=u22*dNxy[2,j]
    
        # xy
        B1[3,c1]=u11*dNxy[2,j] + u12*dNxy[1,j]  
        B1[3,c2]=u21*dNxy[2,j] + u22*dNxy[1,j]

    end

    # G 
    G = zeros(T,4,8)
    for j=1:4

        # columns
        c1 = 2*(j-1)+1
        c2 = c1+1

        G[1,c1]=dNxy[1,j]
        G[2,c1]=dNxy[2,j]
        G[3,c2]=dNxy[1,j]
        G[4,c2]=dNxy[2,j]

    end

    # Return B's , G  and detJ
    return B0, B1, G, detJ

end

"""
Stiffness Matrix and internal force 2D element
    KF_solid2D_nl(m::Mesh2D,ele::Int64,U::Vector)

    U is the current Displacement vector for the structure
"""
function KF_solid2D_nl(m::Mesh2D,ele::Int64,U::Vector)

    # Coordinates
    x,y = Nodal_coordinates(m,ele)

    # Geometry
    geo = m.geo_ele[ele]

    # Thickness
    thick = m.geometries[geo].thickness

    # Constitutive matrix 
    C = Constitutive(m,ele)

    # Stresses
    S = zeros(3)

    # Matrix Shat
    Shat = zeros(4,4)

    # Gauss points
    G = Gauss_2D()

    # Global DOFs
    gls = DOFs(mesh,ele) 

    # Element displacements
    ug = U[gls]

    # Matrix
    K = zeros(8,8)

    # (internal) Force vector
    F = zeros(8)

    # Main loop
    for i=1:4

        # Gauss points
        r,s = G[:,i]

        # Stresses at rs (as a 3×1 vector)
        S .= Stress_solid2D_nl(r,s,mesh,ele,ug)

        # Matriz Shat
        Shat[1,1]=S[1]; Shat[1,2]=S[3]; Shat[2,1]=S[3]; Shat[2,2]=S[2]
        Shat[3,3]=S[1]; Shat[3,4]=S[3]; Shat[4,3]=S[3]; Shat[4,4]=S[2]
    
        # B matrices
        B0,B1,G,dJ = Bs_solid2D_nl(r,s,x,y,ug)

        # B
        B = B0 .+ B1

        # Add 
        K .= K .+ transpose(B)*C*(B)*dJ*thick .+ transpose(G)*Shat*G*dJ*thick
        F .= F .+ transpose(B)*S
    end

    return K,F

end

"""
Second PK Local stresses for solid 2D 
    Stress_solid2D_nl(r::Float64,s::Float64,mesh::Mesh2D,ele::Int64,U::Vector{T})

"""
function Stress_solid2D_nl(r::Float64,s::Float64,mesh::Mesh2D,ele::Int64,ug::Vector{T}) where T

    # Consitutive relation
    C = Constitutive(mesh,ele)

    # Coordinates
    x,y = Nodal_coordinates(mesh,ele)

    # B0,B1,... matrix
    B0,B1,_ = Bs_solid2D_nl(r,s,x,y,ug)
    
    # Second PK Stresses
    C*(B0.+0.5*B1)*ug
    
end
