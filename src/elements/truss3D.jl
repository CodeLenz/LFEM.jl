function K_truss3D(E,A,L)

         SMatrix{6,6,Float64}((E*A/L)*[  1.0 0.0 0.0  -1.0 0.0 0.0;
                                         0.0 0.0 0.0   0.0 0.0 0.0; 
                                         0.0 0.0 0.0   0.0 0.0 0.0;
                                        -1.0 0.0 0.0   1.0 0.0 0.0;
                                         0.0 0.0 0.0   0.0 0.0 0.0;
                                         0.0 0.0 0.0   0.0 0.0 0.0] )

end

#
# Local mass matrix (Lumped)
#
function M_truss3D(dens,A,L,lumped=true)
    SMatrix{6,6,Float64}((dens*A*L/2)*[ 1.0 0.0 0.0  0.0 0.0 0.0;
                                        0.0 1.0 0.0  0.0 0.0 0.0; 
                                        0.0 0.0 1.0  0.0 0.0 0.0;
                                        0.0 0.0 0.0  1.0 0.0 0.0;
                                        0.0 0.0 0.0  0.0 1.0 0.0;
                                        0.0 0.0 0.0  0.0 0.0 1.0] )
end

#
# Local B matrix
#
function B_truss3D(L)
    SMatrix{1,6,Float64}([-1/L 0.0 0.0 1/L 0.0 0.0])
end

#
# Local stress
#
function Stress_truss3D(mesh::Mesh3D,ele::Int64,U::Vector{Float64};xe=1.0,p=1.0,q=0.0)

    # Descobre os dados do elemento
    mat = mesh.mat_ele[ele]
    Ee = mesh.materials[mat].Ex
    Le = BMesh.Length(mesh.bmesh,ele)

    # Matriz B do elemento
    B = B_truss3D(Le)

    # Matriz de rotação
    Te = T_matrix(mesh.bmesh,ele)

    # Global DOFs
    gls = DOFs(mesh.bmesh,ele) 

    # Element displacements in global reference
    ug = SVector{6,Float64}(U[gls])
         
    # Element displacements in local reference
    u = Te*ug

    # Stress
    (xe^(p-q))*Ee*B*u
    

end
