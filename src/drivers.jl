#
# Driver for stress
#
function Stress(mesh::Mesh,ele::Int64,U::Vector{Float64};xe=1.0,p=1.0,q=0.0)
  
   if mesh.bmesh.etype==:truss2D
         return  Stress_truss3D(mesh,ele,U;xe=xe,p=p,q=q)  
   elseif mesh.bmesh.etype==:truss3D
       return  Stress_truss3D(mesh,ele,U;xe=xe,p=p,q=q)  
   else
       error("Stress:: element $mesh.bmesh.etype:: not implemented")
  end
  
end
    
