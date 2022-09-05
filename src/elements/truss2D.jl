"""
Local stiffness matrix for truss2D
       K_truss2D(mesh::Mesh2D,ele::Int64)
"""
function K_truss2D(mesh::Mesh2D,ele::Int64)
  
       # Element properties
       mat = mesh.mat_ele[ele]
       geo = mesh.geo_ele[ele]

       # Material and geometric properties
       Ee = mesh.materials[mat].Ex
       Ae = mesh.geometries[geo].A
       Le = LMesh.Length(mesh,ele)
  
       SMatrix{4,4,Float64}( (Ee*Ae/Le)*[ 1.0 0.0 -1.0 0.0 ;
                                          0.0 0.0  0.0 0.0 ; 
                                         -1.0 0.0  1.0 0.0 ;
                                          0.0 0.0  0.0 0.0 ] )
end
  
"""
Geometric Local stiffness matrix for truss2D
       Ks_truss2D(mesh::Mesh2D,ele::Int64, s::Float64)
"""
function Ks_truss2D(mesh::Mesh2D,ele::Int64,s::Vector{Float64})
  
       # Element properties
       mat = mesh.mat_ele[ele]
       geo = mesh.geo_ele[ele]

       # Geometric properties
       Ae = mesh.geometries[geo].A
       Le = LMesh.Length(mesh,ele)

       SMatrix{4,4,Float64}( (s[1]*Ae/Le)*[ 0.0  0.0  0.0  0.0 ;
                                            0.0  1.0  0.0 -1.0 ; 
                                            0.0  0.0  0.0  0.0 ;
                                            0.0 -1.0  0.0  1.0 ] )
end


"""
Local mass matrix for truss2D (lumped)
       M_truss2D(mesh::Mesh2D,ele::Int64)
"""
function M_truss2D(mesh::Mesh2D,ele::Int64;lumped=true)

       # Element properties
       mat = mesh.mat_ele[ele]
       geo = mesh.geo_ele[ele]

       # Material and geometric properties
       De = mesh.materials[mat].density
       Ae = mesh.geometries[geo].A
       Le = LMesh.Length(mesh,ele)
  
        SMatrix{4,4,Float64}((De*Ae*Le/2)*[ 1.0 0.0  0.0 0.0 ;
                                            0.0 1.0  0.0 0.0 ; 
                                            0.0 0.0  1.0 0.0 ;
                                            0.0 0.0  0.0 1.0 ] )
end


"""
Local B matrix for truss2D
       B_truss2D(mesh::Mesh2D,ele::Int64)
"""
function B_truss2D(mesh::Mesh2D,ele::Int64)
       Le = LMesh.Length(mesh,ele)
       SMatrix{1,4,Float64}( [-1/Le 0.0 1/Le 0.0] )
end

"""
Local stress for truss2D
       Stress_truss2D(mesh::Mesh2D,ele::Int64,U::Vector{T})
       
It returns stress as [sxx] for compatibility with solid elements.
"""
function Stress_truss2D(mesh::Mesh2D,ele::Int64,U::Vector{T}) where T

       # Descobre os dados do elemento
       mat = mesh.mat_ele[ele]
       Ee = mesh.materials[mat].Ex

       # Matriz B do elemento
       B = B_truss2D(mesh,ele)

       # Global DOFs
       gls = DOFs(mesh,ele) 

       # Element displacements in global reference
       ug = SVector{4,T}(U[gls])
       
       # Element displacements in local reference
       u = To_local(ug,mesh,ele)

       # Stress
       Ee*B*u
       
end
