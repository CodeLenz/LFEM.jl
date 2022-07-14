"""
Return the correct constitutive matrix for element ele

    Constitutive(mesh::Mesh,ele::Int64)
    
"""
function Constitutive(mesh::Mesh,ele::Int64)

   # Material data for this element
   mat = Get_material(mesh,ele)

   # Element type
   etype = Get_etype(mesh)

   # Basic data for linear elastic
   Ex = mat.Ex
   νxy = mat.νxy
   model = mat.model

   # Special (simpler) case
   if etype===:truss2D || etype===:truss3D
      return SMatrix{1,1}([Ex])
   end

   # For custom material
   if model===:Custom
   
      # Size of input matrix    
      S1,S2 = size(mat.custom)

      # Convert to StaticMatrix      
      C = SMatrix{S1,S2}(mat.custom)

   elseif etype===:solid3D

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


   elseif model===:EPT && etype===:solid2D

        G = Ex/(2*(1+νxy))
        c = Ex/(1-νxy^2)
        C = SMatrix{3,3,Float64}([  c    νxy*c 0.0 ;
                                νxy*c  c    0.0 ;
                                0.0  0.0    G ]   )

   elseif  model===:EPD && etype===:solid2D

       throw("Constitutive::EPD not implemented")

   else   

       throw("Constitutive::should not happen")
   end  
      
   # Return constitutive matrix
   return C


end


###################################################################################################
#                        von Mises stress matrix
###################################################################################################
function Voigt_equivalent(mesh::Mesh)

    # Get element classe (truss/solid)
    eclass = Get_eclass(mesh)

    if eclass===:truss

        return SMatrix{1,1}([1.0])

    else    

        # Get element type
        etype = Get_etype(mesh)

        if etype===:solid3D

            error("Voigt_equivalent::3D - implementar")

        else 
            
            if model===:EPT
            
                return SMatrix{3,3}(  [1.0 -0.5 0.0;
                                      -0.5  1.0 0.0;
                                       0.0  0.0 3.0] )
            else

                error("Voigt_equivalent::EPD - implementar")

            end

        end

    end 

end
