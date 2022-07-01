#
# Driver for stress (2D)
#
function Stress(r::Float64,s::Float64,mesh::Mesh2D,ele::Int64,U::Vector{Float64};xe=1.0,p=1.0,q=0.0)
  
   etype = Get_etype(mesh)
  
   if etype==:truss2D
         return  Stress_truss2D(mesh,ele,U;xe=xe,p=p,q=q)  
   elseif etype==:solid2D
       return  Stress_solid2D(r,s,mesh,ele,U;xe=xe,p=p,q=q)  
   else
       error("Stress:: element $etype:: not implemented")
   end
  
end

# Default  (r=s=0.0)
function Stress(mesh::Mesh2D,ele::Int64,U::Vector{Float64};xe=1.0,p=1.0,q=0.0)
   Stress(0.0,0.0,mesh,ele,U;xe=xe,p=p,q)
end


#
# Driver for stress (3D)
#
function Stress(r::Float64,s::Float64,t::Float64,
                mesh::Mesh3D,ele::Int64,U::Vector{Float64};xe=1.0,p=1.0,q=0.0)
  
   etype = Get_etype(mesh)
  
   if etype==:truss3D
       return  Stress_truss3D(mesh,ele,U;xe=xe,p=p,q=q)  
   elseif etype==:solid3D
       return  Stress_solid3D(r,s,t,mesh,ele,U;xe=xe,p=p,q=q)  
   else
       error("Stress:: element $etype:: not implemented")
   end
  
end


# Default  (r=s=t=0.0)
function Stress(mesh::Mesh3D,ele::Int64,U::Vector{Float64};xe=1.0,p=1.0,q=0.0)
   Stress(0.0,0.0,0.0,mesh,ele,U;xe=xe,p=p,q)
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



#
# Evaluate stresses for the entire mesh (central points only)
#
"""
Return stresses for the entire mesh. This version evaluates only in the central point.

  Stresses(mesh::Mesh,U::Vector{T};x=Float64[],p=1.0,q=0.0)

The output is a matrix ne x ncol, where ncol is 1 for :truss2D and 3D, 3 for :solid2D and 6 for :solid3D


"""
function Stresses(mesh::Mesh,U::Vector{T};x=Float64[],p=1.0,q=0.0)  where T

   # Number of elements
   ne = Get_ne(mesh)

   # Check for empty x
   if isempty(x)
      x1 = ones(ne)
   else
      x1=x
   end   

   # Alocate the output array. It depends on the number of stresses 
   # for each element type
   etype = Get_etype(mesh)

   ncol = 1
   if etype==:solid2D
      ncol = 3
   elseif etype==:solid3D
      ncol = 6
   end

   # Alocate
   stresses = zeros(ne,ncol)

   # Loop 
   for ele in mesh
      s = Stress(mesh,ele,U,xe=x1[ele],p=p,q=q)
      stresses[ele,:].=s[:]
   end

   # Return stresses
   return stresses

end

