"""
#D Gauss points for the trilinear element

    Gauss_3D()

return a 3x8 matrix with eigth Gauss points

"""
function Gauss_3D()

    pp = 1.0/sqrt(3.0)
    G = SMatrix{3,8,Float64}([-pp  pp  pp  -pp  -pp  pp  pp  -pp;  # r
                              -pp -pp  pp   pp  -pp -pp  pp   pp;  # s
                              -pp -pp -pp  -pp   pp  pp  pp   pp]) # t

end

"""
Derivative of N with respect to r,  s and t
    dN_solid3D(r::T,s::T,t::T) where T
"""
function dN_solid3D(r::T,s::T,t::T) where T

    dN1 = (1/8)*[-(1-s)*(1-t); -(1-r)*(1-t); -(1-r)*(1-s)]
    dN2 = (1/8)*[ (1-s)*(1-t); -(r+1)*(1-t); -(1+r)*(1-s)]
    dN3 = (1/8)*[ (s+1)*(1-t) ; (r+1)*(1-t); -(1+r)*(s+1)]
    dN4 = (1/8)*[-(s+1)*(1-t) ; (1-r)*(1-t); -(1-r)*(s+1)]

    dN5 = (1/8)*[-(1-s)*(1+t); -(1-r)*(1+t);  (1-r)*(1-s)]
    dN6 = (1/8)*[ (1-s)*(1+t); -(r+1)*(1+t);  (1+r)*(1-s)]
    dN7 = (1/8)*[ (s+1)*(1+t) ; (r+1)*(1+t);  (1+r)*(s+1)]
    dN8 = (1/8)*[-(s+1)*(1+t) ; (1-r)*(1+t);  (1-r)*(s+1)]
 
    # Internal nodes (for incompatible element)
    dN9  = [ -2*r ;  0.0  ; 0.0]
    dN10 = [  0.0 ;  -2*s ; 0.0]
    dN11 = [  0.0 ;  0.0  ; -2*t]
    
    return [dN1 dN2 dN3 dN4 dN5 dN6 dN7 dN8 dN9 dN10 dN11]

end

"""
Jacobian matrix for solid3D

    Jacobian_solid3D(x::Vector{T},y::Vector{T},z::Vector{T},
                     dN::Matrix{T}) where T
"""
function Jacobian_solid3D(x::Vector{T},y::Vector{T},z::Vector{T},
                          dN::Matrix{T}) where T

    # Jacobian 
    J11 = zero(T); J12 = zero(T); J13 = zero(T)
    J21 = zero(T); J22 = zero(T); J23 = zero(T)
    J31 = zero(T); J32 = zero(T); J33 = zero(T)
    for i=1:8
       J11 += x[i]*dN[1,i]; J12 += y[i]*dN[1,i]; J13 += z[i]*dN[1,i]
       J21 += x[i]*dN[2,i]; J22 += y[i]*dN[2,i]; J23 += z[i]*dN[2,i]
       J31 += x[i]*dN[3,i]; J32 += y[i]*dN[3,i]; J33 += z[i]*dN[3,i]
    end

    return SMatrix{3,3,T}([J11 J12 J13;
                           J21 J22 J23;
                           J31 J32 J33])

end

"""
B Matrix (with additional bublle functions) for 3D elements
    B_solid3D(r::T,s::T,t::T,x::Vector{T},y::Vector{T},z::Vector{T}) where T
"""
function B_solid3D(r::T,s::T,t::T,x::Vector{T},y::Vector{T},z::Vector{T}) where T

    # Derivates of N
    dNrst = dN_solid3D(r,s,t)

    # Jacobian Matrix
    J = Jacobian_solid3D(x,y,z,dNrst)

    # Determinant
    detJ = det(J)

    # Map derivatives w.r.t (r,s,t) to (x,y,z)
    dNxyz = J\dNrst

    # B matrix
    #
    # Lets use the traditional xx, yy, zz, yz, xz, xy
    #
    # xx xy xz
    #    yy yz
    #       zz

    B = @MMatrix zeros(T,6,33)
    for j=1:11

        # xx
        B[1,3*(j-1)+1] = dNxyz[1,j]
        # yy
        B[2,3*(j-1)+2] = dNxyz[2,j]
        # zz
        B[3,3*(j-1)+3] = dNxyz[3,j]
        # yz
        B[4,3*(j-1)+2] = dNxyz[3,j] 
        B[4,3*(j-1)+3] = dNxyz[2,j]
        # xz
        B[5,3*(j-1)+1] = dNxyz[3,j] 
        B[5,3*(j-1)+3] = dNxyz[1,j]
        # xy
        B[6,3*(j-1)+1] = dNxyz[2,j] 
        B[6,3*(j-1)+2] = dNxyz[1,j]

    end

    # Return B and detJ
    return B, detJ

end

"""
Stiffness Matrix for (incompatible) 3D element
    K_solid3D(m::Mesh3D,ele::Int64)
"""
function K_solid3D(m::Mesh3D,ele::Int64)

    # Coordinates
    x,y,z = Nodal_coordinates(m,ele)

    # Constitutive matrix
    C = Constitutive(m,ele)
    
    # Gauss points
    G = Gauss_3D()

    # Matrix
    K = @MMatrix zeros(33,33)

    # Main loop
    for i=1:8

        # Gauss points
        r,s,t = G[:,i]

        # B matrix
        B, dJ = B_solid3D(r,s,t,x,y,z)

        # Add 
        K .= K .+ transpose(B)*C*B*dJ
        
    end

    # Guyan reduction
    Kaa = @view K[1:24,1:24]
    Kab = @view K[1:24,25:33]
    Kbb = @view K[25:33,25:33]
    return MMatrix{24,24}(Kaa .- Kab*(Kbb\Kab'))

end

"""
Interpolation matrix for solid 3D (For M)
    N_solid3D(r::T,s::T,t::T) where T
"""
function N_solid3D(r::T,s::T,t::T) where T

    # Functions
    N1 = (1/8)*(1-r)*(1-s)*(1-t); 
    N2 = (1/8)*(1+r)*(1-s)*(1-t); 
    N3 = (1/8)*(1+r)*(1+s)*(1-t); 
    N4 = (1/8)*(1-r)*(1+s)*(1-t); 
    N5 = (1/8)*(1-r)*(1-s)*(1+t); 
    N6 = (1/8)*(1+r)*(1-s)*(1+t); 
    N7 = (1/8)*(1+r)*(1+s)*(1+t); 
    N8 = (1/8)*(1-r)*(1+s)*(1+t); 

    return SMatrix{3,24}([N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 0;
                           0 N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 ;
                           0 0 N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8])

end

"""
Consistent mass matrix for solid 3D
    M_solid3D(m::Mesh3D,ele::Int64,lumped=false)
"""
function M_solid3D(m::Mesh3D,ele::Int64,lumped=false)

    # We have to implement the lumped version at the end
    lumped || throw("M_solid3D:: lumped mass still not implemented")

    # Coordinates
    x,y,z = Nodal_coordinates(m,ele)

    # Material 
    mat = m.mat_ele[ele]

    # Density 
    dens = m.materials[mat].density
    
    # Gauss points
    G = Gauss_3D()

    # Matrix
    M = @MMatrix zeros(24,24)

    # Main loop
    for i=1:8

        # Gauss points
        r,s,t = G[:,i]

        # N matrix
        N = N_solid3D(r,s,t)

        # Derivates of N
        dNrst = dN_solid3D(r,s,t)

        # Jacobian matrix
        J = Jacobian_solid3D(x,y,z,dNrst)

        # Add 
        M .= M .+ transpose(N)*N*(det(J)*dens)
        
    end

    # Return M
    return M

end


"""
Local stress for solid 3D ((Not expanding bubble DOFs)
    Stress_solid3D(r::Float64,s::Float64,t::Float64,mesh::Mesh2D,ele::Int64,U::Vector{T})
"""
function Stress_solid3D(r::Float64,s::Float64,t::Float64,
                        mesh::Mesh3D,ele::Int64,U::Vector{T}) where T

    # Consitutive relation
    C = Constitutive(mesh,ele)

    # Coordinates
    x,y,z = Nodal_coordinates(mesh,ele)

    # B matrix
    B, _ = B_solid3D(r,s,t,x,y,z)

    # Global DOFs
    gls = DOFs(mesh,ele) 

    # Element displacements
    ug = SVector{24,T}(U[gls])
    
    # Stress
    MVector{6,T}(C*B[:,1:24]*ug)
    
end


#
# Volume 
#
"""
Return the volume of element ele

   Volume_solid3D(mesh::Mesh,ele::Int64)

"""
function Volume_solid3D(mesh::Mesh,ele::Int64)
    
    # Gauss points
    G = Gauss_3D()
    
    # Coordinates
    x,y,z = Nodal_coordinates(mesh,ele)

    volume = 0.0
    for i=1:8

        # Gauss points
        r,s,t = G[:,i]

        # N matrix
        N = N_solid3D(r,s,t)

        # Derivates of N
        dNrst = dN_solid3D(r,s,t)

        # Jacobian matrix
        J = Jacobian_solid3D(x,y,z,dNrst)
 
        # add
        volume += det(J)

    end
    return volume
    
end
