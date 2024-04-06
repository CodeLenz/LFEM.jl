
"""
2D Gauss points for the bilinear element

    Gauss_2D()

return a 2x4 matrix with four Gauss points

"""
function Gauss_2D()

    pp = 1.0/sqrt(3.0)
    G = SMatrix{2,4,Float64}([-pp  pp pp -pp ; 
                              -pp -pp pp  pp])
end

"""
Derivative of N with respect to r and s
    dN_solid2D(r::T,s::T) where T

"""
function dN_solid2D(r::T,s::T) where T

    dN1 = (1/4)*[ s-1 ;  r-1]
    dN2 = (1/4)*[ 1-s ; -r-1]
    dN3 = (1/4)*[ 1+s ;  r+1]
    dN4 = (1/4)*[-1-s ;  1-r]

    # Internal nodes (for incompatible element)
    dN5 = [ -2*r ;  0.0]
    dN6 = [  0.0 ;  -2*s]

    return [dN1 dN2 dN3 dN4 dN5 dN6]

end

"""
Jacobian matrix for solid2D

    Jacobian_solid2D(x::Vector{T},y::Vector{T},dN::Matrix{T}) where T

"""
function Jacobian_solid2D(x::Vector{T},y::Vector{T},dN::Matrix{T}) where T

    # Jacobian 
    J11 = zero(T)
    J12 = zero(T)
    J21 = zero(T)
    J22 = zero(T)
    @inbounds for i=1:4 
       J11 += x[i]*dN[1,i]
       J12 += y[i]*dN[1,i]
       J21 += x[i]*dN[2,i]
       J22 += y[i]*dN[2,i]
    end

    return SMatrix{2,2,T}([J11 J12 ; J21 J22])

end

"""
B Matrix (with additional bublle functions) for 2D elements
    B_solid2D(r::T,s::T,x::Vector{T},y::Vector{T}) where T

"""
function B_solid2D(r::T,s::T,x::Vector{T},y::Vector{T}) where T

    # Derivates of N
    dNrs = dN_solid2D(r,s)

    # Jacobian Matrix
    J = Jacobian_solid2D(x,y,dNrs)

    # Determinant
    detJ = det(J)

    # Map derivatives w.r.t (r,s) to (x,y)
    dNxy = J\dNrs

    # B matrix
    B = @MMatrix zeros(T,3,12)
    @inbounds for j=1:6

        # xx         
        B[1,2*(j-1)+1]=dNxy[1,j]
        # yy
        B[2,2*(j-1)+2]=dNxy[2,j]
        # xy
        B[3,2*(j-1)+1]=dNxy[2,j] 
        B[3,2*(j-1)+2]=dNxy[1,j]

    end

    # Return B and detJ
    return B, detJ

end

"""
Stiffness Matrix for (incompatible) 2D element
    K_solid2D(m::Mesh2D,ele::Int64)

"""
function K_solid2D(m::Mesh2D,ele::Int64)

    # Coordinates
    x,y = Nodal_coordinates(m,ele)

    # Geometry
    geo = m.geo_ele[ele]

    # Thickness
    thick = m.geometries[geo].thickness

    # Constitutive matrix 
    C = Constitutive(m,ele)
    
    # Gauss points
    G = Gauss_2D()

    # Matrix
    K = @MMatrix zeros(12,12)

    # Main loop
    @inbounds for i=1:4

        # Gauss points
        r,s = G[:,i]

        # B matrix
        B, dJ = B_solid2D(r,s,x,y)

        # Add 
        K .= K .+ transpose(B)*C*B*dJ*thick
        
    end

    # Guyan reduction
    Kaa = @view K[1:8,1:8]
    Kab = @view K[1:8,9:12]
    Kbb = @view K[9:12,9:12]

    return MMatrix{8,8}(Kaa .- Kab*(Kbb\Kab'))

end

"""
Interpolation matrix for solid 2D (For M)
    N_solid2D(r::T,s::T) where T

"""
function N_solid2D(r::T,s::T) where T

    # Functions
    N1 = (1/4)*(1-r)*(1-s)
    N2 = (1/4)*(1+r)*(1-s)
    N3 = (1/4)*(1+r)*(1+s)
    N4 = (1/4)*(1-r)*(1+s)

    return SMatrix{2,8}([N1 0 N2 0 N3 0 N4 0 ; 
                          0 N1 0 N2 0 N3 0 N4 ])

end


"""
Consistent mass matrix for solid 2D
    M_solid2D(m::Mesh2D,ele::Int64,lumped=false)

"""
function M_solid2D(m::Mesh2D,ele::Int64;lumped=false)

    # Coordinates
    x,y = Nodal_coordinates(m,ele)

    # Material 
    mat = m.mat_ele[ele]

    # Geometry
    geo = m.geo_ele[ele]

    # Density and thickness
    dens = m.materials[mat].density
    thick = m.geometries[geo].thickness

    # Gauss points
    G = Gauss_2D()

    # Matrix
    M = @MMatrix zeros(8,8)

    # Main loop
    mass = 0.0
    @inbounds for i=1:4

        # Gauss points
        r,s = G[:,i]

        # N matrix
        N = N_solid2D(r,s)

        # Derivates of N
        dNrs = dN_solid2D(r,s)

        # Jacobian matrix
        J = Jacobian_solid2D(x,y,dNrs)

        # Determinant 
        DJ = det(J)

        # Update mass
        mass += DJ*thick*dens

        # Add 
        M .= M .+ transpose(N)*N*(DJ*thick*dens)
          
    end

    # If lumped, we use the special technique in pg 445 Hugues
    if lumped

        # Dofs X
        glsx = [1;3;5;7]
        diagx = diag(M[glsx,glsx])
        alpx = mass / sum(diagx)
        diagx .*= alpx

        # Dofs Y
        glsy = [2;4;6;8]
        diagy = diag(M[glsy,glsy])
        alpy = mass / sum(diagy)
        diagy .*= alpy
      
        # Mix the X and Y and return the diagonal matrix
        dig = zeros(8)
        dig[glsx].= diagx
        dig[glsy].= diagy
    
        return diagm(dig)
    else
        return M
    end

end

"""
Local stress for solid 2D ((Not expanding bubble DOFs)
    Stress_solid2D(r::Float64,s::Float64,mesh::Mesh2D,ele::Int64,U::Vector{T})

"""
function Stress_solid2D(r::Float64,s::Float64,mesh::Mesh2D,ele::Int64,U::Vector{T}) where T

    # Consitutive relation
    C = Constitutive(mesh,ele)

    # Coordinates
    x,y = Nodal_coordinates(mesh,ele)

    # B matrix
    B, _ = B_solid2D(r,s,x,y)

    # Global DOFs
    gls = DOFs(mesh,ele) 

    # Element displacements
    ug = SVector{8,T}(U[gls])
    
    # Stress
    MVector{3,T}(C*B[:,1:8]*ug)
    
end

#
# Volume 
#
"""
Return the volume of element ele

   Volume_solid2D(mesh::Mesh,ele::Int64)

"""
function Volume_solid2D(mesh::Mesh2D,ele::Int64)
    
    # Gauss points
    G = Gauss_2D()

    # Coordinates
    x,y = Nodal_coordinates(mesh,ele)

    # Thickness
    thick = Get_geometry(mesh,ele).thickness

    volume = 0.0
    @inbounds for i=1:4

        # Gauss points
        r,s = G[:,i]

        # N matrix
        N = N_solid2D(r,s)

        # Derivates of N
        dNrs = dN_solid2D(r,s)

        # Jacobian matrix
        J = Jacobian_solid2D(x,y,dNrs)
 
        # add
        volume += det(J)*thick

    end
    return volume
    
end
