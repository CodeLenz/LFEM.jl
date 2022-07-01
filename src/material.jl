"""
Return the correct constitutive matrix for element ele

    Constitutive(mesh::Mesh,ele::Int64)
    
"""
function Constitutive(mesh::Mesh,ele::Int64)

   # Material data for this element
   mat = Get_material(mesh,ele)

   # Element type
   etype = Get_etype(mesh,ele)

   # Basic data for linear elastic
   Ex = mat.Ex
   νxy = mat.νxy
   model = mat.model

   # For custom material
   if model==:Custom
   
      # Size of input matrix    
      S1,S2 = size(mat.custom)

      # Convert to StaticMatrix      
      C = SMatrix{S1,S2}(mat.custom)

   elseif etype==:solid3D

        G = Ex/(2*(1+νxy))
        c0 = 1-2*νxy^2-νxy
        c2 = (Ex*νxy)/c0
        c1 = Ex/c0 - c2
        C = SMatrix{6,6,Float64}([  c1    c2  c2  0.0 0.0 0.0;
                                    c2    c1  c2  0.0 0.0 0.0;
                                    c2    c2  c1  0.0 0.0 0.0;
                                    0.0  0.0  0.0  G  0.0 0.0;
                                    0.0  0.0  0.0 0.0  G  0.0;
                                    0.0  0.0  0.0 0.0 0.0  G])


   elseif model==:EPT && etype==:solid2D

        G = Ex/(2*(1+νxy))
        c = Ex/(1-νxy^2)
        C = SMatrix{3,3,Float64}([  c    νxy*c 0.0 ;
                                νxy*c  c    0.0 ;
                                0.0  0.0    G ]   )

   elseif  model==:EPD && etype==:solid2D

       throw("Constitutive::EPD not implemented")

   else   

       throw("Constitutive::should not happen")
   end  
      
   # Return constitutive matrix
   return C


end