
#
# Derivative of N with respect to r and s
#
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

#
# Jacobian
#
function Jacobian_solid2D(x::Vector{T},y::Vector{T},dN::Matrix{T}) where T

    # Jacobian 
    J11 = zero(T)
    J12 = zero(T)
    J21 = zero(T)
    J22 = zero(T)
    for i=1:4 
       J11 += x[i]*dN[1,i]
       J12 += y[i]*dN[1,i]
       J21 += x[i]*dN[2,i]
       J22 += y[i]*dN[2,i]
    end

    return SMatrix{2,2,T}([J11 J12 ; J21 J22])

end

#
# B Matrix (with additional functions)
#
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
    B = zeros(T,3,12)
    for j=1:6

        # xx
        B[1,2*(j-1)+1]=dNxy[1,j]
        # yy
        B[2,2*(j-1)+2]=dNxy[2,j]
        # xy
        B[3,2*(j-1)+1]=dNxy[2,j] 
        B[3,2*(j-1)+2]=dNxy[1,j]

    end

    # Return B and detJ
    return SMatrix{3,12,T}(B), detJ

end

#
# Stiffness Matrix
#
function K_solid2D(m::Mesh2D,ele)

    # Alias
    bm = m.bmesh

    # Coordinates
    x,y = Nodal_coordinates(m,ele)

    # Material 
    mat = m.mat_ele[ele]

    # Geometry
    geo = m.geo_ele[ele]

    # Thickness
    thick = m.geometries[geo].thickness

    # Constitutive matrix (EPT)
    Ex = m.materials[mat].Ex
    νxy = m.materials[mat].νxy
    G = Ex/(2*(1+νxy))
    c = Ex/(1-νxy^2)
    C = SMatrix{3,3,Float64}([  c    νxy*c 0.0 ;
                               νxy*c  c    0.0 ;
                               0.0  0.0    G ]   )

    # Gauss points
    pp = 1.0/sqrt(3.0)
    G = [-pp  pp pp -pp ; 
         -pp -pp pp  pp]

    # Matrix
    K = @MMatrix zeros(12,12)

    # Main loop
    for i=1:4

        # Gauss points
        r,s = G[:,i]

        # B matrix
        B, dJ = B_solid2D(r,s,x,y)

        # Add 
        K .= K .+ dot(B,C,B)*dJ*thick
        
    end

    # Guyan reduction
    Kaa = @view K[1:8,1:8]
    Kab = @view K[1:8,9:12]
    Kbb = @view K[9:12,9:12]
    return Kaa .- Kab*(Kbb\Kab')

end

#
# Interpolation (For M)
#
function N_solid2D(r::T,s::T) where T

    # Functions
    N1 = (1/4)*(1-r)*(1-s)
    N2 = (1/4)*(1+r)*(1-s)
    N3 = (1/4)*(1+r)*(1+s)
    N4 = (1/4)*(1-r)*(1+s)

    return SMatrix{2,8}([N1 0 N2 0 N3 0 N4 0 ; 
                          0 N1 0 N2 0 N3 0 N4 ])

end


#
# Mass  Matrix (Consistent)
#
function M_solid2D(m::Mesh2D,ele,lumped=false)

    # Alias
    bm = m.bmesh

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
    pp = 1.0/sqrt(3.0)
    G = [-pp  pp pp -pp ; 
         -pp -pp pp  pp]

    # Matrix
    M = @MMatrix zeros(8,8)

    # Main loop
    for i=1:4

        # Gauss points
        r,s = G[:,i]

        # N matrix
        N = N_solid2D(r,s)

        # Derivates of N
        dNrs = dN_solid2D(r,s)

        # Jacobian matrix
        J = Jacobian_solid2D(x,y,dNrs)

        # Add 
        M .= M .+ dot(N,N)*(det(J)*thick*dens)
          
    end

    # Return M
    return M

end
