#
# Point loads
#
function Point_load(mesh::Mesh,loadcase::Int64=1)
   
    # Basic assertion
    0<=loadcase<=mesh.nload || throw("Point_load:: invalid loadcase")

    # Dimension
    dim = Get_dim(mesh)

    # Create an empty vector
    F = zeros(dim*mesh.bmesh.nn)
    
    # Loop 
    @inbounds for i=1:mesh.nnbc
        
        # Node, dof and value
        no   = Int(mesh.nbc[i,1])
        gl   = Int(mesh.nbc[i,2])
        load = Int(mesh.nbc[i,4])

        if load==loadcase

            val = mesh.nbc[i,3]
        
            # global position
            pos = dim*(no-1)+gl   
        
            # sum this contribution
            F[pos] = F[pos] + val

        end
        
    end #i
    
    # Return the vector
    return F
    
end

