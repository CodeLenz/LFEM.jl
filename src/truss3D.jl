function K_truss3D(E,A,L)

         (E*A/L)*[  1.0 0.0 0.0  -1.0 0.0 0.0;
                    0.0 0.0 0.0   0.0 0.0 0.0; 
                    0.0 0.0 0.0   0.0 0.0 0.0;
                   -1.0 0.0 0.0   1.0 0.0 0.0;
                    0.0 0.0 0.0   0.0 0.0 0.0;
                    0.0 0.0 0.0   0.0 0.0 0.0]

end

#
# Local B matrix
#
function B_truss3D(L)
    [-1/L 0.0 0.0 1/L 0.0 0.0]
end

#
# Local stress
#
function Stress_truss3D(mesh::Mesh3D,ele::Int64,U::Vector{Float64};xe=1.0,p=1.0,q=0.0)

    # Descobre os dados do elemento
    Ee = mesh.materials[1].Ex
    Le = Length(mesh.bmesh,ele)

    # Matriz B do elemento
    B = B_truss3D(Le)

    # Matriz de rotação
    Te = T_matrix(mesh.bmesh,ele)

    # Global DOFs
    gls = DOFs(mesh.bmesh,ele) 

    # Local displacements
    u = Te*U[gls]

    # Stress
    (xe^(p-q))*Ee*B*u
    

end
