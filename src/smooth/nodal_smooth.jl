#
# Stress smoothing by simple nodal average
#
function Nodal_stress_smooth(mesh,stresses::Array{TF}) where TF

    # As we are considering stresses just at the center of each element
    # it is a good idea to check for the number of collumns in 
    # stresses
    nc = size(stresses,2)

    # Assertion
    (nc==3 || nc==6) || error("Nodal_stress_smooth:: implemented for centroidal stresses only")

    # Allocate the output array
    smooth = zeros(TF,mesh.bmesh.nn,nc)

    # First, we must find the elements sharing each node. BMesh has 
    # a function for it
    neighbours = BMesh.Elements_sharing_nodes(mesh.bmesh)
    
    # We are assuming that stresses were evaluated at the center
    # of each element. Thus, a simple average is performed at each node
    for node=1:mesh.bmesh.nn 

        # Neighbours 
        elements = neighbours[node]

        # "volumes" for each element
        vols =  [Volume_element(mesh,e) for e in elements]

        # Stresses for each element
        sigmas = stresses[elements,:]
        
        # Scale each stress for its volume
        scale_sigma = sigmas .* vols

        # For each element and each component, evaluate the mean
        smooth[node,:] .= vec(sum(scale_sigma,dims=1)) ./ sum(vols)

    end

    # Return smoothed stresses
    return smooth

end