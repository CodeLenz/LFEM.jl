#
# Driver for stress (2D)
#
function Stress(r::Float64,s::Float64,mesh::Mesh2D,ele::Int64,U::Vector{T}) where T
  
   etype = Get_etype(mesh)
  
   if etype==:truss2D
         return  Stress_truss2D(mesh,ele,U)  
   elseif etype==:solid2D
       return  Stress_solid2D(r,s,mesh,ele,U)  
   else
       error("Stress:: element $etype:: not implemented")
   end
  
end

# Default  (r=s=0.0)
function Stress(mesh::Mesh2D,ele::Int64,U::Vector{T}) where T
   Stress(0.0,0.0,mesh,ele,U)
end


#
# Driver for stress (3D)
#
function Stress(r::Float64,s::Float64,t::Float64,
                mesh::Mesh3D,ele::Int64,U::Vector{T}) where T
  
   etype = Get_etype(mesh)
  
   if etype==:truss3D
       return  Stress_truss3D(mesh,ele,U)  
   elseif etype==:solid3D
       return  Stress_solid3D(r,s,t,mesh,ele,U)  
   else
       error("Stress:: element $etype:: not implemented")
   end
  
end


# Default  (r=s=t=0.0)
function Stress(mesh::Mesh3D,ele::Int64,U::Vector{T}) where T
   Stress(0.0,0.0,0.0,mesh,ele,U)
end


# von-Mises Equivalent stress
function Equivalent_stress(sigma, mesh::Mesh, eps=1E-5)

   # Voigt matrix
   V = Voigt_equivalent(mesh,1)

   sqrt(real(dot(sigma,V,sigma))+eps^2)

end

# von-Mises Equivalent stress
# V is the Voigt matrix
function Equivalent_stress(sigma,V,eps=1E-5)
   
   sqrt(real(dot(sigma,V,sigma))+eps^2)
   
end


#
# Driver for Stiffness (2D)
#
function Local_K(mesh::Mesh2D,ele::Int64)

   # Element type
   etype = Get_etype(mesh)
  
   # If solid
   if etype==:solid2D 

      Ke = K_solid2D(mesh,ele)

   elseif etype==:truss2D

      Ke = K_truss2D(mesh,ele)

   else
      
      throw("Local_K:: etype $etype is invalid")

   end
   
   return Ke

end

#
# Driver for Stiffness (3D)
#
function Local_K(mesh::Mesh3D,ele::Int64)

   # Element type
   etype = Get_etype(mesh)
  
   # If solid
   if etype==:solid3D 

      Ke = K_solid3D(mesh,ele)

   elseif etype==:truss3D  
      
      Ke = K_truss3D(mesh,ele)

   else
      
      throw("Local_K:: etype $etype is invalid")

   end
   
   return Ke

end


#
# Driver for Geometric Stiffness (2D)
#
function Local_Ks(mesh::Mesh2D,ele::Int64,stress::AbstractVector)

   # Element type
   etype = Get_etype(mesh)
  
   # If solid
   if etype==:solid2D 

      error("Local_Ks:: não implementado para elasticidade 2D")

   elseif etype==:truss2D

      Ke = Ks_truss2D(mesh,ele,stress)

   else
      
      throw("Local_Ks:: etype $etype is invalid")

   end
   
   return Ke

end

#
# Driver for Geometric Stiffness (3D)
#
function Local_Ks(mesh::Mesh3D,ele::Int64, stress::AbstractVector)

   # Element type
   etype = Get_etype(mesh)
  
   # If solid
   if etype==:solid3D 

      
      error("Local_Ks:: não implementado para elasticidade 3D")


   elseif etype==:truss3D  
      
       
      Ke = Ks_truss3D(mesh,ele,stress)


   else
      
      throw("Local_Ks:: etype $etype is invalid")

   end
   
   return Ke

end


#
# Driver for Mass (2D)
#
function Local_M(mesh::Mesh2D,ele::Int64)

   # Element type
   etype = Get_etype(mesh)
  
   # If solid
   if etype==:solid2D 

      Me = M_solid2D(mesh,ele)

   elseif etype==:truss2D

      Me = M_truss2D(mesh,ele)

   else
      
      throw("Local_M:: etype $etype is invalid")

   end
   
   return Me

end


#
# Driver for Mass (3D)
#
function Local_M(mesh::Mesh3D,ele::Int64)

   # Element type
   etype = Get_etype(mesh)
  
   # If solid
   if etype==:solid3D 

      Me = M_solid3D(mesh,ele)

   elseif etype==:truss3D  
      
      Me = M_truss3D(mesh,ele)

   else
      
      throw("Local_M:: etype $etype is invalid")

   end
   
   return Me

end



#
# Volume for truss2D and truss2D
#
"""
Return the volume of element ele

   Volume_truss(mesh::Mesh,ele::Int64)

Valid for :truss2D and :truss3D
"""
function Volume_truss(mesh::Mesh,ele::Int64)
       
   # Area
   geometry = Get_geometry(mesh,ele) 
   A = geometry.A

   # Length
   L = LMesh.Length(mesh,ele) 

   return L*A
end

#
# Return the volume
#
"""
Return the volume of element ele

   Volume_element(mesh::Mesh,ele::Int64)

Valid for all elements.
"""
function Volume_element(mesh::Mesh,ele::Int64)

   # Element type
   etype = Get_etype(mesh)
  
   volume = 0.0
   if etype==:truss2D || etype==:truss3D
      volume =  Volume_truss(mesh,ele)
   elseif etype==:solid2D
      volume =  Volume_solid2D(mesh,ele)
   else 
      volume = Volume_solid3D(mesh,ele)
   end
   return volume
end


"""
Return the B matrix of element ele - in the center
AND WITH NO BUBBLE

   B_element(mesh::Mesh2D,ele::Int64)

Valid for all elements.
"""
function B_element(mesh::Mesh2D,ele::Int64)

   # Element type
   etype = Get_etype(mesh)
  
   if etype===:truss2D 
      B =  B_truss2D(mesh,ele)
   elseif etype===:solid2D
      x,y = Nodal_coordinates(mesh,ele)
      B_,_ = B_solid2D(0.0,0.0,x,y)
      B = @view B_[:,1:8]
   else
     error("B_element::element type not defined")
   end
  
  return B
  
end


"""
Return the B matrix of element ele - in the center
AND WITH NO BUBBLE

   B_element(mesh::Mesh3D,ele::Int64)

Valid for all elements.
"""
function B_element(mesh::Mesh3D,ele::Int64)

   # Element type
   etype = Get_etype(mesh)
  
   if etype===:truss3D
      B =  B_truss3D(mesh,ele)
   elseif etype===:solid3D
      x,y,z = Nodal_coordinates(mesh,ele)
      B_,_= B_solid3D(0.0,0.0,0.0,x,y,z)
      B = @view B_[:,1:24]
   else
     error("B_element::element type not defined")
   end
  
  return B
  
end
