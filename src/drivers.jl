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
# Driver for Keg
#
function Keg_truss(mesh::Mesh,ele::Int64)
  
  # Descobre os dados do elemento
  mat = mesh.mat_ele[ele]
  Ee = mesh.materials[mat].Ex
  
  geo = mesh.geo_ele[ele]
  Ae = mesh.geometries[geo].A
  Le = Length(mesh.bmesh,ele)
     
  # Element type
  etype = mesh.bmesh.etype
  
  # Monta a matriz local (4 × 4)
  if etype==:truss2D
     Ke = K_truss2D(Ee,Ae,Le)
  elseif etype==:truss3D
     Ke = K_truss3D(Ee,Ae,Le)
  else
     error("Keg_truss::elemento $etype ainda não implementado")
  end

  # Evaluate the rotation matrix for this element
  Te = T_matrix(mesh.bmesh,ele)

  # Rotaciona a matriz local para o sistema global 
  Ke .= transpose(Te)*Ke*Te
 
  # Retorna a matriz local no sistema global
  return Ke
  
end


#
# Driver for Meg
#
function Meg_truss(mesh::Mesh,ele::Int64)
  
   # Descobre os dados do elemento
   mat = mesh.mat_ele[ele]
   Ee = mesh.materials[mat].Ex
   dense = mesh.materials[mat].density
   
   geo = mesh.geo_ele[ele]
   Ae = mesh.geometries[geo].A
 
   Le = Length(mesh.bmesh,ele)
      
   # Element type
   etype = mesh.bmesh.etype
   
   # Monta a matriz local (4 × 4)
   if etype==:truss2D
      Me = M_truss2D(dense,Ae,Le)
   elseif etype==:truss3D
      Me = M_truss3D(dense,Ae,Le)
   else
      error("Meg_truss::elemento $etype ainda não implementado")
   end
 
   # Evaluate the rotation matrix for this element
   Te = T_matrix(mesh.bmesh,ele)
 
   # Rotaciona a matriz local para o sistema global 
   Me .= transpose(Te)*Me*Te
  
   # Retorna a matriz local no sistema global
   return Me
   
 end
 
