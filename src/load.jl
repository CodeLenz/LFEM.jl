#
# Point loads
#
function Point_load(mesh::Mesh)
   
    # Dimension
    dim = 2
    if isa(mesh,Mesh3D)
        dim = 3
    end

    # Create an empty vector
    F = zeros(dim*mesh.bmesh.nn)
    
    # Loop 
    for i=1:mesh.nnbc
        
        # Node, dof and value
        no = Int(mesh.nbc[i,1])
        gl = Int(mesh.nbc[i,2])
        val = mesh.nbc[i,3]
        
        # global position
        pos = dim*(no-1)+gl   
        
        # sum this contribution
        F[pos] = F[pos] + val
        
    end #i
    
    # Return the vector
    return F
    
end

