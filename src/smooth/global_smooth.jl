#
# Global smooth 
#

# Build the global "mass" matrix 
function Global_smooth_M(mesh::Mesh)

    # Element type
    etype = Get_etype(mesh)

    # Functions and number of nodes depend on the 
    # type of element
    if etype===:solid2D
        nn = 4
        f(m,e) = M_smooth_solid2D(m,e)
    elseif etype===:solid3D
        nn = 8
        f(m,e) = M_smooth_solid3D(m,e)
    else
        error("Global_smooth_M:: just for elasticity")
    end

    # Allocate auxiliar vectors
    VI = Int64[]; 
    VJ = Int64[]; 
    VV = Float64[]; 

    # Loop pelos elementos, calculando a matriz local Me de cada um
    # e posicionando na M
    for ele in mesh

        # Local "mass" matrix
        Me = f(mesh,ele)
        
        # Find element nodes
        nodes = Connect(mesh,ele)

        # Adiciona a matriz do elemento (rotacionada) à matriz Global
        @inbounds for i=1:nn
            gi = nodes[i]
            @inbounds for j=1:nn
                gj = nodes[j]
                push!(VI,gi)
                push!(VJ,gj)
                push!(VV, Me[i,j])
            end #j
        end #i

    end #ele

    # Generate the sparse matrix
    M = sparse(VI,VJ,VV)
    dropzeros!(M)

    # Return the global matrix
    return M
    
end


# Build the global "force" vector for a given stress component
function Global_smooth_F(mesh::Mesh,component::Vector{T}) where T

    # Element type
    etype = Get_etype(mesh)

    # Functions  depending on the type of element
    if etype===:solid2D
        fF(m,e,t) = F_smooth_solid2D(m,e,t)
    elseif etype===:solid3D
        fF(m,e,t) = F_smooth_solid3D(m,e,t)
    else
        error("Global_smooth_F:: just for elasticity")
    end

    # Number of nodes in the mesh
    nn = Get_nn(mesh)

    # Allocate the output vector
    F = zeros(T,nn)
    
    # Loop over elements
    for ele in mesh

        # Stress component for this element
        stress = component[ele]

        # Local "force" vector
        Fe = fF(mesh,ele,stress)
        
        # Find element nodes
        nodes = Connect(mesh,ele)

        # Add to the global vector
        F[nodes] .+= Fe

    end #ele

    # Return the global vector
    return F
    
end

#
# Global stress smoothing
#
function Global_stress_smooth(mesh,stresses::Array{TF}) where TF

    # As we are considering stresses just at the center of each element
    # it is a good idea to check for the number of collumns in 
    # stresses
    nc = size(stresses,2)

    # Assertion
    (nc==3 || nc==6) || error("Global_stress_smooth:: implemented for centroidal stresses only")

    # Allocate the output array
    smooth = zeros(TF,mesh.bmesh.nn,nc)

    # Build the global "mass" matrix (only once)
    M = Global_smooth_M(mesh)

    # Cholesky decomposition
    cM = cholesky(M)

    # For each component, build the "force" vector and solves the 
    # linear system
    for j=1:nc

        # "Force" vector
        F = Global_smooth_F(mesh,stresses[:,j])

        # Store the smoothed solution
        smooth[:,j] .= cM\F

    end

    # return the smoothed stresses
    return smooth

end