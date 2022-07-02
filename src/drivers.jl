#
# Driver for stress (2D)
#
function Stress(r::Float64,s::Float64,mesh::Mesh2D,ele::Int64,U::Vector{Float64})
  
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
function Stress(mesh::Mesh2D,ele::Int64,U::Vector{Float64})
   Stress(0.0,0.0,mesh,ele,U)
end


#
# Driver for stress (3D)
#
function Stress(r::Float64,s::Float64,t::Float64,
                mesh::Mesh3D,ele::Int64,U::Vector{Float64})
  
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
function Stress(mesh::Mesh3D,ele::Int64,U::Vector{Float64})
   Stress(0.0,0.0,0.0,mesh,ele,U)
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



