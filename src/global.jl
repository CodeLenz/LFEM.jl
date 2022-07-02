"""
Assembly the global stiffness matrix.

    Global_K(mesh::Mesh, x::Vector{Float64}, kparam::Function)

where kparam(x) can be, for example

function kparam(x,p=1.0)
    x^p
end

for a SIMP like material parametrization.
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

    # Primeira coisa é alocar a matriz K 
    ng = dim*nn
    K = spzeros(ng,ng)

    # Loop pelos elementos, calculando a matriz local Ke de cada um
    # e posicionando na K
    for ele in mesh
    
        # Local stiffness matrix
        Ke = Local_K(mesh,ele) 

        # Determina quais são os gls GLOBAIS que são "acessados"
        # por esse elemento
        gls = DOFs(mesh,ele) 

        # If needed, convert to global reference
        Keg = To_global(Ke,mesh,ele)
        
        # Adiciona a matriz do elemento (rotacionada) à matriz Global
        K[gls,gls] .= K[gls,gls] .+ Keg*kparam(x[ele])

    end #ele

    # Retorna a matriz global
    return Symmetric(K)

end

#
# Version without parametrization
#
"""
Assembly the global stiffness matrix.

    Global_K(mesh::Mesh)

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
    
    # Vamos reforçar a ideia de que x tem a dimensão ne
    @assert length(x)==ne "Global_M::x deve ter dimensão igual o número de elementos"

    # Dimensão do problema
    dim = Get_dim(mesh)

    # Tipo de elemento
    etype = Get_etype(mesh)

    # Primeira coisa é alocar a matriz M 
    ng = dim*nn
    M = spzeros(ng,ng)

    # Loop pelos elementos, calculando a matriz local Me de cada um
    # e posicionando na M
    for ele in mesh

        # Local mass matrix
        Me = Local_M(mesh,ele)
           
        # Determina quais são os gls GLOBAIS que são "acessados"
        # por esse elemento
        gls = DOFs(mesh,ele) 

        # If needed, convert to global reference
        Meg = To_global(Me,mesh,ele)
        
        # Adiciona a matriz do elemento (rotacionada) a matriz Global
        M[gls,gls] .= M[gls,gls] .+ Meg*mparam(x[ele])

    end #ele

    # Retorna a matriz global
    return Symmetric(M)

end

"""
Assembly the global mass matrix.

     Global_M(mesh::Mesh)

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
   stresses = zeros(ne,ncol)

   # Loop 
   for ele in mesh
      s = Stress(mesh,ele,U)
      stresses[ele,:].=s[:]*sparam(x[ele])
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
