#
# Evaluate the relative error for a given element. For simplicity,
# we assume that the element is not distorted 
#
# sigma: element (centroidal) stresses (non smoothed)
# nodal: nodal stresses obtained by some smoothing procedure
#
function Element_error_stress(mesh::Mesh2D,ele::Int,sigma::Matrix{T},nodal::Matrix{T}) where T

    # Since sigma is at the centroid and nodal is at the nodes,
    # we must interpolate the nodal stresses to the same point (0,0)
    N = N_smooth_solid2D(0.0,0.0)

    # Nodes of this element
    nodes = Connect(mesh,ele)

    # Lets interpolate the three components at the same time
    interpo = [N;N;N]*(nodal[nodes,:])

    # The "error" at the center is
    error = sigma[ele,:] .- interpo

    # Area for the element
    V = Volume_element(mesh,ele)

    # The "integral" is 
    e_s = sqrt(dot(error,error)*V)

end

function Element_error_stress(mesh::Mesh3D,ele::Int,sigma::Matrix{T},nodal::Matrix{T}) where T

    # Since sigma is at the centroid and nodal is at the nodes,
    # we must interpolate the nodal stresses to the same point (0,0)
    N = N_smooth_solid3D(0.0,0.0,0.0)

    # Nodes of this element
    nodes = Connect(mesh,ele)

    # Lets interpolate the three components at the same time
    interpo = [N;N;N;N;N;N]*(nodal[nodes,:])

    # The "error" at the center is
    error = sigma[ele,:] .- interpo

    # Area for the element
    V = Volume_element(mesh,ele)

    # The "integral" is 
    e_s = sqrt(dot(error,error)*V)

end