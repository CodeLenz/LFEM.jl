#
# Local stiffness matrix
#
function K_truss2D(mesh::Mesh,ele::Int64)
  
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
  
#
# Local mass matrix (Lumped)
#
function M_truss2D(mesh::Mesh,ele::Int64;lumped=true)

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


#
# Local B matrix
#
function B_truss2D(mesh::Mesh,ele)
       Le = LMesh.Length(mesh,ele)
       SMatrix{1,4,Float64}( [-1/Le 0.0 1/Le 0.0] )
end

#
# Local stress
#
function Stress_truss2D(mesh::Mesh2D,ele::Int64,U::Vector{Float64};xe=1.0,p=1.0,q=0.0)

       # Descobre os dados do elemento
       mat = mesh.mat_ele[ele]
       Ee = mesh.materials[mat].Ex

       # Matriz B do elemento
       B = B_truss2D(mesh,ele)

       # Global DOFs
       gls = DOFs(mesh,ele) 

       # Element displacements in global reference
       ug = SVector{4,Float64}(U[gls])
       
       # Element displacements in local reference
       u = To_local(ug,mesh,ele)

       # Stress
       (xe^(p-q))*Ee*B*u
       

end
