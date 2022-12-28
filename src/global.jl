"""
Assembly the global stiffness matrix.

    Global_K(mesh::Mesh, x::Vector{Float64}, kparam::Function)

where kparam(x) can be, for example

function kparam(x,p=1.0)
    x^p
end

for a SIMP like material parametrization.
    
This function also considers entries :Stiffness in mesh.options
with [node dof value;]

"""
function Global_K(mesh::Mesh, x::Vector{T}, kparam::Function) where T

    # Alias
    ne = Get_ne(mesh)
    nn = Get_nn(mesh)

    # Basic check
    length(x)==ne || throw("Global_K::x must have the same dimensions as the number of elements")

    # Dimensão do problema
    dim = Get_dim(mesh)
    
    # Tipo de elemento
    etype = Get_etype(mesh)

    # Options
    options = mesh.options

    # Primeira coisa é alocar a matriz K 
    ng = dim*nn

    # Vamos usar um sizehint de 1% de esparsividade
    # hint = round(Int64,0.01*(ng^2))

    # Aloca arrays para usar o sparse
    VI = Int64[]; #sizehint!(I,hint)
    VJ = Int64[]; #sizehint!(J,hint)
    VV = Float64[]; #sizehint!(V,hint)

    # Chama dofs uma vez para depois reaproveitar 
    # o acesso de memória
    gls = DOFs(mesh,1) 

    # Comprimento do vetor gls
    s_gls = length(gls)

    # Chama o primeiro elemento aqui para podermos
    # reaproveitar nos loops
    Ke = Local_K(mesh,1) 

    # Mesma coisa
    Keg = similar(Ke)

    # Experimental!!
    # If this key is set, than 
    # use Keg from element 1 ONLY
    if haskey(options,:IS_TOPO) && Get_eclass(mesh)==:solid

       @inbounds for ele in mesh
            # Determina quais são os gls GLOBAIS que são "acessados"
            # por esse elemento
            gls .= DOFs(mesh,ele) 

            # Parametrização para esse elemento
            kx = kparam(x[ele])

            # Adiciona a matriz do elemento à matriz Global
            @inbounds for i=1:s_gls
                gi = gls[i]
                @inbounds for j=1:s_gls
                    gj = gls[j]
                    push!(VI,gi)
                    push!(VJ,gj)
                    push!(VV, Ke[i,j]*kx)
                end #j
            end #i       
        end #ele

    else

        # Loop pelos elementos, calculando a matriz local Ke de cada um
        # e posicionando na K
        @inbounds for ele in mesh
        
            # Local stiffness matrix
            Ke .= Local_K(mesh,ele) 

            # Determina quais são os gls GLOBAIS que são "acessados"
            # por esse elemento
            gls .= DOFs(mesh,ele) 

            # If needed, convert to global reference
            Keg .= To_global(Ke,mesh,ele)
            
            # Adiciona a matriz do elemento (rotacionada) à matriz Global
            # Parametrização para esse elemento
            kx = kparam(x[ele])

            # Adiciona a matriz do elemento à matriz Global
            @inbounds for i=1:s_gls
                gi = gls[i]
                @inbounds for j=1:s_gls
                    gj = gls[j]
                    push!(VI,gi)
                    push!(VJ,gj)
                    push!(VV, Keg[i,j]*kx)
                end #j
            end #i

        end #ele

    end # :IS_TOPO

    # Add options:: :Stiffness
    # If there are lumped stiffness, we add here
   
    # If :Stiffness is defined
    if haskey(options,:Stiffness)

       # Alias 
       stiffness = options[:Stiffness]

       # Dimensions
       nstiff,ncol = size(stiffness)

       # Check if ncols==3
       # node gl value
       ncol==3 || throw("Global_K:: :Stiffness must have 3 columns")

       # Add stiffness
       @inbounds for s=1:nstiff
 
           # Recover data for this stiffness
           node   = Int(stiffness[s,1])
           dof    = Int(stiffness[s,2])
           gl     = dim*(node-1)+dof

           value  = stiffness[s,3]
           
           # Assert if valid gl
           0<node<=nn   || throw("Global_K:: :Stiffness :: invalid node")
           0<dof<=dim   || throw("Global_K:: :Stiffness :: invalid dof")
           value >= 0.0 || throw("Global_K:: :Stiffness :: invalid value")

           push!(VI,gl)
           push!(VJ,gl)
           push!(VV,value)

       end #Stiff

    end # if :Stiffness

    # Generate the sparse matrix
    K = sparse(VI,VJ,VV)
    dropzeros!(K)

    # Return the global matrix
    return Symmetric(K)

end

#
# Version without parametrization
#
"""
Assembly the global stiffness matrix.

    Global_K(mesh::Mesh)

        
This function also considers entries :Stiffness in mesh.options
with [node dof value;]
    
"""
function Global_K(mesh::Mesh)

    # Dummy function 
    dummy_f(x)=1.0

    # x is ones
    ne = Get_ne(mesh)
    x = ones(ne)

    # Call Global_K
    Global_K(mesh, x, dummy_f)

end



"""
Assembly the global mass matrix.

     Global_M(mesh::Mesh, x=Float64[], mparam::Function)

where mparam(x) can be, for example

function param(x,p=2.0,cut=0.1)
    s = x
    if x<cut
       s = x^p
    end
    return s
end

for a SIMP like material parametrization.
        
This function also considers entries :Mass in mesh.options
with [node dof value;]
    
"""
 function Global_M(mesh::Mesh, x::Vector{T}, mparam::Function) where T

    # Alias
    ne = Get_ne(mesh)
    nn = Get_nn(mesh)
    
    # Options
    options = mesh.options

    # Basic check
    length(x)==ne || throw("Stresses::x must have the same dimensions as the number of elements")

    # Dimensão do problema
    dim = Get_dim(mesh)

    # Tipo de elemento
    etype = Get_etype(mesh)

    # Primeira coisa é alocar a matriz M 
    ng = dim*nn
    
    # Vamos usar um sizehint de 1% de esparsividade
    # hint = round(Int64,0.01*(ng^2))

    VI = Int64[]; #sizehint!(I,hint)
    VJ = Int64[]; #sizehint!(J,hint)
    VV = Float64[]; #sizehint!(V,hint)

    # Chama dofs uma vez para depois reaproveitar 
    # o acesso de memória
    gls = DOFs(mesh,1) 
    s_gls = length(gls)
   
    # Aloca M aqui para aproveitar depois
    Me = Local_M(mesh,1) 
    Meg = similar(Me)

    # Experimental!!
    # If this key is set, than 
    # use Meg from element 1 ONLY
    if haskey(options,:IS_TOPO) && Get_eclass(mesh)==:solid

        for ele in mesh
             # Determina quais são os gls GLOBAIS que são "acessados"
             # por esse elemento
             gls .= DOFs(mesh,ele) 
 
             # Parametrização do material
             mx = mparam(x[ele])

             # Adiciona a matriz do elemento (rotacionada) à matriz Global
             @inbounds for i=1:s_gls
                gi = gls[i]
                @inbounds for j=1:s_gls
                    gj = gls[j]
                    push!(VI,gi)
                    push!(VJ,gj)
                    push!(VV, Me[i,j]*mx)
                end #j
            end #i

        end
 
    else
 
        # Loop pelos elementos, calculando a matriz local Me de cada um
        # e posicionando na M
        for ele in mesh

            # Local mass matrix
            Me .= Local_M(mesh,ele)
            
            # Determina quais são os gls GLOBAIS que são "acessados"
            # por esse elemento
            gls .= DOFs(mesh,ele) 

            # If needed, convert to global reference
            Meg .= To_global(Me,mesh,ele)
            
            # Adiciona a matriz do elemento (rotacionada) a matriz Global
            # Parametrização do material
            mx = mparam(x[ele])

            # Adiciona a matriz do elemento (rotacionada) à matriz Global
            @inbounds for i=1:s_gls
               gi = gls[i]
               @inbounds for j=1:s_gls
                   gj = gls[j]
                   push!(VI,gi)
                   push!(VJ,gj)
                   push!(VV, Meg[i,j]*mx)
               end #j
           end #i


        end #ele

    end #IS_TOPO
        
    # Add options:: :Mass
    # If there are lumped mass, we add here
    options = mesh.options

    # If :Mass is defined
    if haskey(options,:Mass)

       # Alias 
       mass = options[:Mass]

       # Dimensions
       nmass,ncol = size(mass)

       # Check if ncols==3
       # node gl value
       ncol==3 || throw("Global_M:: :Mass must have 3 columns")

       # Add mass
       @inbounds for m=1:nmass
 
           # Recover data for this stiffness
           node   = Int(mass[m,1])
           dof    = Int(mass[m,2])
           gl     = dim*(node-1)+dof
           value  = mass[m,3]
           
           # Assert if valid gl
           0<node<=nn   || throw("Global_M:: :Mass :: invalid node")
           0<dof<=dim   || throw("Global_M:: :Mass :: invalid dof")
           value >= 0.0 || throw("Global_M:: :Mass :: invalid value")

           push!(VI,gl)
           push!(VJ,gl)
           push!(VV,value)

       end # Mass

    end # if :Mass

    # Generate the sparse matrix
    M = sparse(VI,VJ,VV)
    dropzeros!(M)

    # Return the global mass matrix
    return Symmetric(M)

end

"""
Assembly the global mass matrix.

     Global_M(mesh::Mesh)
         
This function also considers entries :Stiffness in mesh.options
with [node dof value;]
    
"""
function Global_M(mesh::Mesh)

    # Dummy function
    dummy_f(x)=1.0

    # x is ones
    ne = Get_ne(mesh)
    x = ones(ne)

    # Call function
    Global_M(mesh, x, dummy_f)
end


"""
Assembly the global damping matrix  C = α_c M + β_c K

     Global_C(M,K,mesh::Mesh,α_c=0.0,β_c=0.0)

This function also considers entries :Damper in mesh.options
with [node dof value;]
    
"""
function Global_C(M::AbstractMatrix{T},K::AbstractMatrix{T},mesh::Mesh,α_c=0.0,β_c=1E-6) where T

    # Damping matrix
    C = α_c*M .+ β_c*K

    # Add lumped dampers
    # Add options:: :Damper
    # If there are lumped mass, we add here
    options = mesh.options

    # If :Damper is defined
    if haskey(options,:Damper)

       # Alias
       dim = Get_dim(mesh) 
       nn = Get_nn(mesh)

       # Alias 
       damper = options[:Damper]

       # Dimensions
       ndamper,ncol = size(damper)

       # Check if ncols==3
       # node gl value
       ncol==3 || throw("Global_C:: :Damper must have 3 columns")

       # Add damper
       @inbounds for m=1:ndamper
 
           # Recover data for this damper
           node   = Int(damper[m,1])
           dof    = Int(damper[m,2])
           gl     = dim*(node-1)+dof
           value  = damper[m,3]
           
           # Assert if valid gl
           0<node<=nn   || throw("Global_C:: :Damper :: invalid node")
           0<dof<=dim   || throw("Global_C:: :Damper :: invalid dof")
           value >= 0.0 || throw("Global_C:: :Damper :: invalid value")

           # Add damper
           C[gl,gl] += value

       end # damper

    end # if :Damper

    return Symmetric(C)

end


"""
Assembly the global geometric stiffness matrix.

    Global_Ks(mesh::Mesh, stress::Array{Float64})

"""
function Global_Ks(mesh::Mesh, stress::Array{T}) where T

    # Alias
    ne = Get_ne(mesh)
    nn = Get_nn(mesh)

    # Dimensão do problema
    dim = Get_dim(mesh)
    
    # Tipo de elemento
    etype = Get_etype(mesh)

    # Options
    options = mesh.options

    # Primeira coisa é alocar a matriz Ks
    ng = dim*nn

    # Vamos usar um sizehint de 1% de esparsividade
    # hint = round(Int64,0.01*(ng^2))

    VI = Int64[]; #sizehint!(I,hint)
    VJ = Int64[]; #sizehint!(J,hint)
    VV = Float64[]; #sizehint!(V,hint)

    # Chama dofs uma vez para depois reaproveitar 
    # o acesso de memória
    gls = DOFs(mesh,1) 
    s_gls = length(gls)

    # Pré-aloca
    sigma = stress[1,:]

    # Local geometric stiffness matrix
    Kse = Local_Ks(mesh,1,sigma) 
    Kseg = similar(Kse)

    # Loop pelos elementos, calculando a matriz local Ke de cada um
    # e posicionando na K
    @inbounds for ele in mesh
        
        # Stress for element ele (returns a vector)
        sigma .= stress[ele,:]

        # Local geometric stiffness matrix
        Kse .= Local_Ks(mesh,ele,sigma) 

        # Determina quais são os gls GLOBAIS que são "acessados"
        # por esse elemento
        gls .= DOFs(mesh,ele) 

        # If needed, convert to global reference
        Kseg .= To_global(Kse,mesh,ele)
            
        # Adiciona a matriz do elemento (rotacionada) à matriz Global
        @inbounds for i=1:s_gls
            gi = gls[i]
            @inbounds for j=1:s_gls
                gj = gls[j]
                push!(VI,gi)
                push!(VJ,gj)
                push!(VV, Kseg[i,j])
            end #j
        end #i

    end #ele

    # Generate the sparse matrix
    Ks = sparse(VI,VJ,VV)
    dropzeros!(Ks)

    # Return the global matrix 
    return Symmetric(Ks)

end




########################## Static Stresses ###################

#
# Evaluate stresses for the entire mesh (central points only)
#
"""
Return stresses for the entire mesh. This version evaluates only in the central point.

  Stresses(mesh::Mesh,U::Vector{T};x=Float64[],sparam::Function; center=true)

The output is a matrix ne x ncol, where ncol is 1 for :truss2D and 3D, 3 for :solid2D and 6 for :solid3D
if center = true and ncol = 3*4 for 2D and 6*8 for 3D if center = false since we return stresses for 
each superconvergent (Gauss) point in the incompatible elements

Function sparam(x) is used to parametrize stress. One possibility is 

function sparam(x,p=1.0,q=0.0) 
       x^(p=q)
end

"""
function Stresses(mesh::Mesh,U::Vector{T},x::Vector{T1},sparam::Function; center=true)  where {T,T1}

    # Number of elements
    ne = Get_ne(mesh)

    # Basic check
    length(x)==ne || throw("Stresses::x must have the same dimensions as the number of elements")

   # Alocate the output array. It depends on the number of stresses 
   # for each element type and the number of Gauss Points
   etype = Get_etype(mesh)

   ncol = 1
   if etype==:solid2D
      ncol= ifelse(center,3,3*4)
   elseif etype==:solid3D
      ncol = ifelse(center,6,6*8)
   end

   # Alocate
   stresses = Matrix{T}(undef,ne,ncol)

   # If center, solve and bail out 
   if center || ncol==1
        
        # Allocate before the loop 
        s = Stress(mesh,1,U)

        # Loop 
        @inbounds for ele in mesh
            s .= Stress(mesh,ele,U)
            stresses[ele,:].=s[:]*sparam(x[ele])
        end

    else

        # Solid 2D,four Gauss Points        
        if etype==:solid2D
       
            # Gauss Points
            G = Gauss_2D()

            # pré-aloca
            s = Stress(G[1,1],G[2,1],mesh,1,U)

            # For each element,also loop at the Gauss Points
            @inbounds for ele in mesh
                xp = sparam(x[ele])
                @inbounds for j=1:4
                    s .= Stress(G[1,j],G[2,j],mesh,ele,U)
                    p1 = 3*(j-1)+1
                    p2 = 3*(j-1)+3
                    stresses[ele,p1:p2].=s[:]*xp
                end
            end
        elseif etype==:solid3D

             # Gauss Points
             G = Gauss_3D()

             # Pre-aloca
             s = Stress(G[1,1],G[2,1],G[3,1],mesh,1,U)

             # For each element,also loop at the Gauss Points
             @inbounds for ele in mesh
                 xp = sparam(x[ele])
                 @inbounds for j=1:8
                     s .= Stress(G[1,j],G[2,j],G[3,j],mesh,ele,U)
                     p1 = 6*(j-1)+1
                     p2 = 6*(j-1)+6
                     stresses[ele,p1:p2].=s[:]*xp
                 end
             end
            
        else
            error("Should not happen")
        end
    end   

   # Return stresses
   return stresses

end

"""
Return stresses for the entire mesh. This version evaluates only in the central point.

  Stresses(mesh::Mesh,U::Vector{T}; center=true)

The output is a matrix ne x ncol, where ncol is 1 for :truss2D and 3D, 3 for :solid2D and 6 for :solid3D
if center = true and ncol = 3*4 for 2D and 6*8 for 3D if center = false since we return stresses for 
each superconvergent (Gauss) point in the incompatible elements
    

"""
function Stresses(mesh::Mesh,U::Vector{T};center=true)  where T

    # Dummy function
    dummy_f(x)=1.0

    # x is ones
    ne = Get_ne(mesh)
    x = ones(ne)

    # Call function
    Stresses(mesh,U,x,dummy_f,center=center)
end


########################## Harmonic Stresses ###################

#
# Evaluate stresses for the entire mesh (central points only)
#
"""
Return (harmonic) stresses for the entire mesh. This version evaluates only in the central point.

  Harmonic_stresses(mesh::Mesh,U::Vector{T}, w::Float64, β_c::Float64,
                    x=Float64[],sparam::Function; center=true)

where w is the angular frequency and β_c is the damping parameter.

The output is a matrix ne x ncol, where ncol is 1 for :truss2D and 3D, 3 for :solid2D and 6 for :solid3D
if center = true and ncol = 3*4 for 2D and 6*8 for 3D if center = false since we return stresses for 
each superconvergent (Gauss) point in the incompatible elements
    

Function sparam(x) is used to parametrize stress. One possibility is 

function sparam(x,p=1.0,q=0.0) 
       x^(p=q)
end

"""
function Harmonic_stresses(mesh::Mesh,U::Vector{T}, w::Float64, β_c::Float64,
                           x::Vector{T1},sparam::Function;center=true)  where {T,T1}


        Stresses(mesh,U,x,sparam,center=center).*(1+im*w*β_c)

end
"""
Return stresses for the entire mesh. This version evaluates only in the central point.

    Harmonic_stresses(mesh::Mesh,U::Vector{T}, w::Float64, β_c::Float64)


where w is the angular frequency and β_c is the damping parameter.

The output is a matrix ne x ncol, where ncol is 1 for :truss2D and 3D, 3 for :solid2D and 6 for :solid3D
if center = true and ncol = 3*4 for 2D and 6*8 for 3D if center = false since we return stresses for 
each superconvergent (Gauss) point in the incompatible elements
    

"""
function Harmonic_stresses(mesh::Mesh,U::Vector{T}, w::Float64, β_c::Float64; center=true)  where T

    # Dummy function
    dummy_f(x)=1.0

    # x is ones
    ne = Get_ne(mesh)
    x = ones(ne)

    # Call function
    Harmonic_stresses(mesh,U,w, β_c,x,dummy_f,center=center)
end
