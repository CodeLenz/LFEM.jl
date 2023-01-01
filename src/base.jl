#
# Copy v to V[free_dofs] inplace
#
function Expand_vector!(V::Vector{T},v::Vector{T},free_dofs::Vector{Int64}) where {T}

    # Basic assertions
    length(v)==length(free_dofs)||throw("Expand_vector!:: v and free_dofs must have the same size")
    length(V)>=length(free_dofs)||throw("Expand_vector!:: length of V must be larger than the length of v")

    # Copy
    V[free_dofs].=v
end


function Expand_vector(v::Vector{T},dim::Int64,free_dofs::Vector{Int64}) where {T}

    V = zeros(T,dim)
    Expand_vector!(V,v,free_dofs)
    return V

end

#
# Function to consider the loadcase
#
function Expand_vector(v::Vector{T},mesh::Mesh,loadcase::Int64) where {T}

  # Output vector
  dimension = Get_nn(mesh)*Get_dim(mesh)
  V = zeros(T,dimension)

  # Free dofs for this loadcase
  0<=loadcase<=mesh.nload || throw("Expand_vector:: invalid loadcase")
  free_dofs = mesh.free_dofs[loadcase]

  Expand_vector!(V,v,free_dofs)
  return V

end


 #
 # Change of reference Local to Global
 #
 function To_global(M::AbstractMatrix,mesh::Mesh,ele::Int64)

    # Bail out if solid
    Get_eclass(mesh)===:solid && return M
    
    # Evaluate the rotation matrix for this element
    Te = T_matrix(mesh.bmesh,ele)
  
    # Rotaciona a matriz local para o sistema global 
    transpose(Te)*M*Te
   
  end
 
  #
  # Change of reference  Global to Local
  #
  function To_local(M::AbstractMatrix,mesh::Mesh,ele::Int64)
 
    # Bail out if solid
    Get_eclass(mesh)===:solid && return M
    
    # Evaluate the rotation matrix for this element
    Te = T_matrix(mesh.bmesh,ele)
  
    # Rotaciona a matriz local para o sistema global 
    Te*M*transpose(Te)
   
  end
 
 
  #
  # Change of reference Local to Global
  #
  function To_global(V::AbstractVector,mesh::Mesh,ele::Int64)
 
    # Bail out if solid
    Get_eclass(mesh)===:solid && return V
    
    # Evaluate the rotation matrix for this element
    Te = T_matrix(mesh.bmesh,ele)
  
    # Rotaciona o vetor local para o sistema global 
    transpose(Te)*V
   
  end
 
  #
  # Change of reference Global to Local
  #
  function To_local(V::AbstractVector,mesh::Mesh,ele::Int64)
 
    # Bail out if solid
    Get_eclass(mesh)===:solid && return V
    
    # Evaluate the rotation matrix for this element
    Te = T_matrix(mesh.bmesh,ele)
  
    # Rotaciona o vetor global para o local 
    Te*V
   
  end
 
