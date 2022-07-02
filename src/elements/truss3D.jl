"""
Local stiffness matrix for truss3D
       K_truss3D(mesh::Mesh3D,ele::Int64)
"""
function K_truss3D(mesh::Mesh3D,ele::Int64)
  
        # Element properties
        mat = mesh.mat_ele[ele]
        geo = mesh.geo_ele[ele]

        # Material and geometric properties
        Ee = mesh.materials[mat].Ex
        Ae = mesh.geometries[geo].A
        Le = LMesh.Length(mesh,ele)

        SMatrix{6,6,Float64}((Ee*Ae/Le)*[  1.0 0.0 0.0  -1.0 0.0 0.0;
                                           0.0 0.0 0.0   0.0 0.0 0.0; 
                                           0.0 0.0 0.0   0.0 0.0 0.0;
                                          -1.0 0.0 0.0   1.0 0.0 0.0;
                                           0.0 0.0 0.0   0.0 0.0 0.0;
                                           0.0 0.0 0.0   0.0 0.0 0.0] )

end

"""
Local mass matrix for truss3D (lumped)
       M_truss3D(mesh::Mesh3D,ele::Int64)
"""
function M_truss3D(mesh::Mesh,ele::Int64;lumped=true)

    # Element properties
    mat = mesh.mat_ele[ele]
    geo = mesh.geo_ele[ele]

    # Material and geometric properties
    De = mesh.materials[mat].density
    Ae = mesh.geometries[geo].A
    Le = LMesh.Length(mesh,ele)

    SMatrix{6,6,Float64}((De*Ae*Le/2)*[ 1.0 0.0 0.0  0.0 0.0 0.0;
                                        0.0 1.0 0.0  0.0 0.0 0.0; 
                                        0.0 0.0 1.0  0.0 0.0 0.0;
                                        0.0 0.0 0.0  1.0 0.0 0.0;
                                        0.0 0.0 0.0  0.0 1.0 0.0;
                                        0.0 0.0 0.0  0.0 0.0 1.0] )
end

"""
Local B matrix for truss3D
       B_truss3D(mesh::Mesh3D,ele::Int64)
"""
function B_truss3D(mesh::Mesh3D,ele::Int64)
    Le = LMesh.Length(mesh,ele)
    SMatrix{1,6,Float64}([-1/Le 0.0 0.0 1/Le 0.0 0.0])
end

"""
Local stress for truss3D
       Stress_truss3D(mesh::Mesh3D,ele::Int64,U::Vector{Float64})
       
It returns stress as [sxx] for compatibility with solid elements.
"""
function Stress_truss3D(mesh::Mesh3D,ele::Int64,U::Vector{Float64})

    # Descobre os dados do elemento
    mat = mesh.mat_ele[ele]
    Ee = mesh.materials[mat].Ex
    
    # Matriz B do elemento
    B = B_truss3D(mesh,ele)

    # Global DOFs
    gls = DOFs(mesh,ele) 

    # Element displacements in global reference
    ug = SVector{6,Float64}(U[gls])
         
    # Element displacements in local reference
    u = To_local(ug,mesh,ele)

    # Stress
    Ee*B*u
    
end
