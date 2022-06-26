
#
# Derivative of N with respect to r, s and t
#
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

#
# Jacobian
#
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

    return [J11 J12 J13;
            J21 J22 J23;
            J31 J32 J33]

end

#
# B Matrix (with additional functions)
#
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

    B = zeros(T,6,33)
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

#
# Stiffness Matrix
#
function K_solid3D(m::Mesh3D,ele)

    # Alias
    bm = m.bmesh

    # Coordinates
    x,y,z = Nodal_coordinates(m,ele)

    # Material 
    mat = m.mat_ele[ele]

    # Constitutive matrix
    Ex = m.materials[mat].Ex
    νxy = m.materials[mat].νxy
    G = Ex/(2*(1+νxy))
    c0 = 1-2*νxy^2-νxy
    c2 = (Ex*νxy)/c0
    c1 = Ex/c0 - c2
    C = [  c1    c2  c2  0.0 0.0 0.0;
           c2    c1  c2  0.0 0.0 0.0;
           c2    c2  c1  0.0 0.0 0.0;
           0.0  0.0  0.0  G  0.0 0.0;
           0.0  0.0  0.0 0.0  G  0.0;
           0.0  0.0  0.0 0.0 0.0  G]   

    # Gauss points
    pp = 1.0/sqrt(3.0)
    G = [-pp  pp  pp  -pp  -pp  pp  pp  -pp; # r
         -pp -pp  pp   pp  -pp -pp  pp   pp; # s
         -pp -pp -pp  -pp   pp  pp  pp   pp] # t

    # Matrix
    K = zeros(33,33)

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
    return Kaa .- Kab*(Kbb\Kab')

end

#
# Interpolation (For M)
#
function N_solid3D(r::T,s::T,t::T) where T

    # Functions
    N1 = (1/8)*(1-r)*(1-s)*(1-t); 
    N2 = (1/8)*(1+r)*(1-s)*(1-t); 
    N3 = (1/8)*(1+r)*(1+s)*(1-t); 
    N4 = (1/8)*(1-r)*(1+s)*(1-t); 
    N5 = (1/8)*(1-r)*(1-s)*(1-t); 
    N6 = (1/8)*(1+r)*(1-s)*(1-t); 
    N7 = (1/8)*(1+r)*(1+s)*(1-t); 
    N8 = (1/8)*(1-r)*(1+s)*(1-t); 

    return [N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 0;
            0 N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 ;
            0 0 N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8]

end

#
# Mass  Matrix (Consistent)
#
function M_solid3D(m::Mesh3D,ele,lumped=false)

    # Alias
    bm = m.bmesh

    # Coordinates
    x,y,z = Nodal_coordinates(m,ele)

    # Material 
    mat = m.mat_ele[ele]

    # Density 
    dens = m.materials[mat].density
    
    # Gauss points
    pp = 1.0/sqrt(3.0)
    G = [-pp  pp  pp  -pp  -pp  pp  pp  -pp; # r
         -pp -pp  pp   pp  -pp -pp  pp   pp; # s
         -pp -pp -pp  -pp   pp  pp  pp   pp] # t

    # Matrix
    M = zeros(24,24)

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
