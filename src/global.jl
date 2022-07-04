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
function Global_K(mesh::Mesh, xin::Vector{Float64}, kparam::Function)

    # Alias
    ne = Get_ne(mesh)
    nn = Get_nn(mesh)

     # Basic test for xe
    if isempty(xin)
        x = ones(ne)
    else
        x = copy(xin)
    end

    # Vamos reforçar a ideia de que x tem a dimensão ne
    @assert length(x)==ne "Global_K::x deve ter dimensão igual o número de elementos"

    # Dimensão do problema
    dim = Get_dim(mesh)
    
    # Tipo de elemento
    etype = Get_etype(mesh)

    # Options
    options = mesh.options

    # Primeira coisa é alocar a matriz K 
    ng = dim*nn
    K = spzeros(ng,ng)

    # Chama dofs uma vez para depois reaproveitar 
    # o acesso de memória
    gls = DOFs(mesh,1) 

    # Faz o mesmo para a matriz de rigidez
    Keg = Local_K(mesh,1) 
    Ke  = similar(Keg)
 
    # Experimental!!
    # If this key is set, than 
    # use Keg from element 1 ONLY
    if haskey(options,:IS_TOPO) && Get_eclass(mesh)==:solid
 
       for ele in mesh
            # Determina quais são os gls GLOBAIS que são "acessados"
            # por esse elemento
            gls .= DOFs(mesh,ele) 

            # Adiciona a matriz do elemento (rotacionada) à matriz Global
            @inbounds K[gls,gls] .= K[gls,gls] .+ Keg*kparam(x[ele])
       end

    else

        # Loop pelos elementos, calculando a matriz local Ke de cada um
        # e posicionando na K
        for ele in mesh
        
            # Local stiffness matrix
            Ke .= Local_K(mesh,ele) 

            # Determina quais são os gls GLOBAIS que são "acessados"
            # por esse elemento
            gls .= DOFs(mesh,ele) 

            # If needed, convert to global reference
            Keg .= To_global(Ke,mesh,ele)
            
            # Adiciona a matriz do elemento (rotacionada) à matriz Global
            @inbounds K[gls,gls] .= K[gls,gls] .+ Keg*kparam(x[ele])

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

           K[gl,gl] += value

       end #Stiff

    end # if :Stiffness

    # Retorna a matriz global
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

    # Call Global_K
    Global_K(mesh, Float64[], dummy_f)
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
 function Global_M(mesh::Mesh, xin::Vector{Float64}, mparam::Function)

    # Alias
    ne = Get_ne(mesh)
    nn = Get_nn(mesh)

    # Basic test for xe
    if isempty(xin)
        x = ones(ne)
    else
        x = copy(xin)
    end
    
    # Options
    options = mesh.options

    # Vamos reforçar a ideia de que x tem a dimensão ne
    @assert length(x)==ne "Global_M::x deve ter dimensão igual o número de elementos"

    # Dimensão do problema
    dim = Get_dim(mesh)

    # Tipo de elemento
    etype = Get_etype(mesh)

    # Primeira coisa é alocar a matriz M 
    ng = dim*nn
    M = spzeros(ng,ng)

    # Chama dofs uma vez para depois reaproveitar 
    # o acesso de memória
    gls = DOFs(mesh,1) 

    # Faz o mesmo para a matriz de massa
    Meg = Local_M(mesh,1) 
    Me  = similar(Meg)

    # Experimental!!
    # If this key is set, than 
    # use Meg from element 1 ONLY
    if haskey(options,:IS_TOPO) && Get_eclass(mesh)==:solid

        for ele in mesh
             # Determina quais são os gls GLOBAIS que são "acessados"
             # por esse elemento
             gls .= DOFs(mesh,ele) 
 
             # Adiciona a matriz do elemento (rotacionada) à matriz Global
             @inbounds M[gls,gls] .= M[gls,gls] .+ Meg*mparam(x[ele])
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
            @inbounds M[gls,gls] .= M[gls,gls] .+ Meg*mparam(x[ele])

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

           M[gl,gl] += value

       end # Mass

    end # if :Mass

    # Retorna a matriz global
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

    # Call function
    Global_M(mesh, Float64[], dummy_f)
end
#
# Evaluate stresses for the entire mesh (central points only)
#
"""
Return stresses for the entire mesh. This version evaluates only in the central point.

  Stresses(mesh::Mesh,U::Vector{T};x=Float64[],sparam::Function)

The output is a matrix ne x ncol, where ncol is 1 for :truss2D and 3D, 3 for :solid2D and 6 for :solid3D

Function sparam(x) is used to parametrize stress. One possibility is 

function sparam(x,p=1.0,q=0.0) 
       x^(p=q)
end

"""
function Stresses(mesh::Mesh,U::Vector{T},xin::Vector{Float64},sparam::Function)  where T

    # Number of elements
    ne = Get_ne(mesh)

    # Basic test for xe
    if isempty(xin)
        x = ones(ne)
    else
        x = copy(xin)
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
   stresses = Matrix{Float64}(undef,ne,ncol)

   # Loop 
   for ele in mesh
      s = Stress(mesh,ele,U)
      @inbounds stresses[ele,:].=s[:]*sparam(x[ele])
   end

   # Return stresses
   return stresses

end

"""
Return stresses for the entire mesh. This version evaluates only in the central point.

  Stresses(mesh::Mesh,U::Vector{T})

The output is a matrix ne x ncol, where ncol is 1 for :truss2D and 3D, 3 for :solid2D and 6 for :solid3D

"""
function Stresses(mesh::Mesh,U::Vector{T})  where T

    # Dummy function
    dummy_f(x)=1.0

    # Call function
    Stresses(mesh,U,Float64[],dummy_f)
end
