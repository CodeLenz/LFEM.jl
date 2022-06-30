#
# Overload some functions defined in BMesh
#
import BMesh:Conec
"""
  Return the connectivities of element ele

    Conect(mesh::Mesh,ele::Int64)

  as an Int64 vector.  
"""
function Conect(mesh::Mesh,ele::Int64)

    BMesh.Conect(mesh.bmesh,ele)

end


import BMesh:Coord
"""
  Return the coordinates of node 

    Coord(mesh::Mesh,node::Int64)

  as a vector [x;y;(z)]
"""
function Coord(mesh::Mesh,node::Int64)

   BMesh.Coord(mesh.bmesh,node)
    
end

import BMesh:Length
"""
  Return the distance between two nodes of the element

    Length(mesh::Mesh,ele::Int64;nodes=(1,2))

  where the default is the distance between (local) nodes 1 and 2
"""
function Length(mesh::Mesh,ele::Int64;nodes=(1,2))

    BMesh.Length(mesh.bmesh,ele;nodes=nodes)

end


import BMesh:DOFs
"""
Return a dim*n vector with the global DOFs of element ele

    DOFs(mesh::Mesh,ele::Int64)

where n is the number of nodes of ele and dim=2 or 3.
"""
function DOFs(mesh::Mesh,ele::Int64)

    BMesh.DOFs(mesh.bmesh,ele)

end

import BMesh:T_matrix
"""
Evaluate the rotation matrix T that maps from global to local reference systems. This version
evaluates Rotation internally.

    T_matrix(mesh::Mesh, ele::Int64, α=0.0)

"""
function T_matrix(mesh::Mesh, ele::Int64, α=0.0)
    
    BMesh.T_matrix(mesh.bmesh,ele,α)
    
end

