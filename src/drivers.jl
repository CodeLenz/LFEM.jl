#
# Driver for stress
#
function Stress(mesh::Mesh,ele::Int64,U::Vector{Float64};xe=1.0,p=1.0,q=0.0)
  
   if mesh.bmesh.etype==:truss2D
         return  Stress_truss2D(mesh,ele,U;xe=xe,p=p,q=q)  
   elseif mesh.bmesh.etype==:truss3D
       return  Stress_truss3D(mesh,ele,U;xe=xe,p=p,q=q)  
   else
       error("Stress:: element $mesh.bmesh.etype:: not implemented")
  end
  
end
    
#
# Driver for Stiffness
#
function Local_K(mesh::Mesh,ele::Int64)

   # Element type
   etype = mesh.bmesh.etype
  
   # If solid
   if etype==:solid2D 

      Ke = K_solid2D(mesh,ele)

   elseif etype==:solid3D

      Ke = K_solid3D(mesh,ele)   

   elseif etype==:truss2D

      Ke = K_truss2D(mesh,ele)

   elseif etype==:truss3D  
      
      Ke = K_truss3D(mesh,ele)

   else
      
      throw("Local_K:: etype $etype is invalid")

   end
   
   return Ke

end


#
# Driver for Mass
#
function Local_M(mesh::Mesh,ele::Int64)

   # Element type
   etype = mesh.bmesh.etype
  
   # If solid
   if etype==:solid2D 

      Me = M_solid2D(mesh,ele)

   elseif etype==:solid3D

      Me = M_solid3D(mesh,ele)   

   elseif etype==:truss2D

      Me = M_truss2D(mesh,ele)

   elseif etype==:truss3D  
      
      Me = M_truss3D(mesh,ele)

   else
      
      throw("Local_M:: etype $etype is invalid")

   end
   
   return Me

end


