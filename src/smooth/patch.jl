
#
# Data type containing the basic data of a given Patch
#
struct Patch

    # Central element
    central_ele::Int64

    # List of neighbours (elements)
    neighbors::Vector{Int64}

    # List of nodes (all the nodes in the patch)
    nodes::Vector{Int64}

    # Extreme coordinates
    lower_coord::Vector{Float64}
    upper_coord::Vector{Float64}

    # Constructor
    function Patch(mesh::Mesh,central_ele::Int64)

        # Evaluate the Patch data
        neighbors, lower_coord, upper_coord, nodes = Patch_element(mesh,central_ele)

        # Create the patch
        new(central_ele,neighbors,nodes,lower_coord,upper_coord)

    end
end


#
# Given a central element, return the first order 
# neighborhood as a list of unique elements and also
# the extreme coordinates of the patch (bounding box)
#
function Patch_element(mesh::Mesh,central_ele::Int64)

    # Nodes of this element
    nodes = Connect(mesh,central_ele)
    
    # Allocate an empty vector for the first order
    # neighbors (elements)
    neighbors = Int64[]

    # Allocate an empty vector for all the nodes in the patch
    vector_nodes = Int64[]

    # Alias for connect
    connect = mesh.bmesh.connect

    # Extreme coordinates. Lets start with the extremes for the 
    # entire mesh. Yes...we start it "swapped"
    lower_coord = vec(maximum(mesh.bmesh.coord,dims=1))
    upper_coord = vec(minimum(mesh.bmesh.coord,dims=1))

    # Dimension
    nc = length(lower_coord)

    # Find all elements sharing nodes with the central element
    for node in nodes

        # Loop over the elements
        for ele in mesh 

            # Nodes for this element
            nodes_ele = connect[ele,:] 
 
            # If this element contains the node, we 
            # store it in the list
            if node in nodes_ele

                # Store the element 
                push!(neighbors,ele)

                # Store the nodes in vector_nodes
                vector_nodes = [vector_nodes;nodes_ele]

                # Coordinates of the nodes of this element
                # x, y, (z)
                coord = Nodal_coordinates(mesh,ele)

                # Loop over the coordinates
                for j=1:nc

                    lower_coord[j] = min(lower_coord[j],minimum(coord[j]))
                    upper_coord[j] = max(upper_coord[j],maximum(coord[j]))

                end #j

            end #if
        end #ele
    end #node

    # As there will be repeated entries, we return just que unique entries
    return sort(unique(neighbors)), lower_coord, upper_coord, sort(unique(vector_nodes))

end


#
# Map coordinates to the natural coordinates of a Patch
#
function Map_coordinates_to_patch(patch::Patch,coord)

    # Dimension
    dim = length(coord)

    # Output array
    output = zeros(dim)

    # Extreme coordinates 
    lower = patch.lower_coord
    upper = patch.upper_coord

    # We must assert that the coordinates are INSIDE the patch
    all(lower .<= coord .<= upper) || error("Map_coordinates_to_patch::Invalid coordinates")

    # general mapping for a coordinate
    gm(x,xl,xu) = -1 + 2*(x-xl)/(xu-xl)

    # For each coordinate, do the mapping
    for i=1:dim
        output[i] = gm(coord[i],lower[i],upper[i])
    end

    return output

end
