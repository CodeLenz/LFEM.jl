
# Bilinear base
bibase(r,s) = @SMatrix [1.0 r s r*s]
 
# Trilinear base
tribase(r,s,t) = @SMatrix [1.0 r s t r*s*t]

#
# Given a list of natural coordinates in [-1,1], return the matrix A for the Least Square
#

#
# Coordinates must be np x dim array with the normalizes values in [-1,1]
#
function Matrix_A_patch(coordinates::Array)

 
   # Dimension
   dim = size(coordinates,2)

   # size
   s = ifelse(dim==2,4,5)

   # Initialize matriz A
   A = zeros(s,s)

   # Loop over the vectors
   for j=1:size(coordinates,1)

       # Evaluate the base at this point 
       if dim==2
          p = bibase(coordinates[j,1],coordinates[j,2])
       else
          p = tribase(coordinates[j,1],coordinates[j,2],coordinates[j,3])
       end

       # Add
       A .+= transpose(p)*p

   end

   # Return A
   return A
   
end


#
# components must be element-wise component, relative to the coordinates 
# centroids
#
function Vector_b_patch(coordinates::Array,component::Vector{T}) where T

  
  
    # Dimension
    dim = size(coordinates,2)
 
    # size
    s = ifelse(dim==2,4,5)
 
    # Initialize vector b
    b = zeros(s)
 
    # Loop over the vectors
    for j=1:size(coordinates,1)
 
        # Evaluate the base at this point 
        if dim==2
           p = bibase(coordinates[j,1],coordinates[j,2])
        else
           p = tribase(coordinates[j,1],coordinates[j,2],coordinates[j,3])
        end
 
        # Add
        b .+= transpose(p).*component[j]
 
    end
 
    # Return b
    return b
    
 end
 

 ############################# Two almost equal functions due to some problem
 ############################# when returning the smooth(rs)=dot(a,base(rs[1],rs[2],..))

#
# Given a Patch and a stress component, return the coefficients 
# of the patch
#
function Patch_smooth2D(patch::Patch,component::Vector{T}) where T

    # The dimension of the problem can be asserted from the dimension of lower_coord
    dim = 2

    # Number of neighbours in this patch
    nneig = length(patch.neighbors)

    # Array of coordinates
    coordinates = zeros(nneig,dim)

    # Components IN THE CENTROIDS of the PATCH
    vector_component = zeros(nneig)

    # First, we must loop along the elements of the patch
    # finding the coordinates of its centroid
    cont = 1
    for ele in patch.neighbors

        # Centroid of this element
        centroid = Centroid(mesh,ele)

        # Map the centroid to the natural coordinates of the patch
        rst = Map_coordinates_to_patch(patch,centroid)

        # Store the component
        vector_component[cont] = component[ele]

        # Add to the array of coordinates
        coordinates[cont,:] .= vec(rst)
        cont += 1  

    end # eles

    # Build A
    A = Matrix_A_patch(coordinates)

    # Vector b for this patch, component
    b = Vector_b_patch(coordinates,vector_component)

    # Coefficientes are
    a = A\b

    return smooth(rs)=dot(a,bibase(rs[1],rs[2]))
    
end


#
# Given a Patch and a stress component, return the coefficients 
# of the patch
#
function Patch_smooth3D(patch::Patch,component::Vector{T}) where T

    # The dimension of the problem can be asserted from the dimension of lower_coord
    dim = 3

    # Number of neighbours in this patch
    nneig = length(patch.neighbors)

    # Array of coordinates
    coordinates = zeros(nneig,dim)

    # Components IN THE CENTROIDS of the PATCH
    vector_component = zeros(nneig)

    # First, we must loop along the elements of the patch
    # finding the coordinates of its centroid
    cont = 1
    for ele in patch.neighbors

        # Centroid of this element
        centroid = Centroid(mesh,ele)

        # Map the centroid to the natural coordinates of the patch
        rst = Map_coordinates_to_patch(patch,centroid)

        # Store the component
        vector_component[cont] = component[ele]

        # Add to the array of coordinates
        coordinates[cont,:] .= vec(rst)
        cont += 1  

    end # eles

    # Build A
    A = Matrix_A_patch(coordinates)

    # Vector b for this patch, component
    b = Vector_b_patch(coordinates,vector_component)

    # Coefficientes are
    a = A\b

    return smooth(rst)=dot(a,tribase(rst[1],rst[2],rst[3]))
    
end

#
# Loop over all patches and project the interpolated component stresses 
# over the nodes
#
function Patches(mesh::Mesh, stresses::Array{T}) where T

    # Number of nodes
    nn = Get_nn(mesh)

    # Number of stress components
    nsc = size(stresses,2)

    # Nodal stresses
    nodal = zeros(nn,nsc)

    # Number of accesses
    acesses = zeros(Int64,nn)

    # Dimensionality
    dim = ifelse(nsc==3,2,3)

    # Loop over central elements
    for ele in mesh

        # Build the patch
        p = Patch(mesh,ele)

        # For each component, find the LS and 
        # project onto the nodes of the patch
        for j=1:nsc

            # Function 
            func = if dim==2
                      Patch_smooth2D(p,stresses[:,j]) 
                   else
                      Patch_smooth3D(p,stresses[:,j]) 
                   end

            # Nodes for this patch
            nodes = p.nodes

            # Accesses
            acesses[nodes] .+= 1

            # Loop over the nodes of the patch
            for node in p.nodes

                # Coordinates of this node
                coord = Coord(mesh,node)

                # Map to the patch
                rst = Map_coordinates_to_patch(p,coord)

                # Evaluate the component at the position
                comp = func(rst)

                # Add value to the nodal array
                nodal[node,j] += comp
                
            end # node

        end # j

    end # ele

    # Return the nodal average 
    return nodal ./ acesses

end