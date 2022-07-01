function Global_K(mesh::Mesh; x=Float64[], p=1.0)

    # Alias
    bmesh = mesh.bmesh
    ne = bmesh.ne

    # Basic test for x
    if isempty(x)
        x = ones(ne)
    end
    
    # Vamos reforçar a ideia de que x tem a dimensão ne
    @assert length(x)==bmesh.ne "Global_K::x deve ter dimensão igual o número de elementos"

    p>=1.0 || throw("Global_K:::p must be larger or equal to 1.0")
    
    # Dimensão do problema
    dim = Get_dim(mesh)
    
    # Tipo de elemento
    etype = Get_etype(mesh)

    # Primeira coisa é alocar a matriz K 
    ng = dim*bmesh.nn
    K = spzeros(ng,ng)

    # Flag se for barras
    flag_truss = Get_eclass(mesh)==:truss

    # Loop pelos elementos, calculando a matriz local Ke de cada um
    # e posicionando na K
    for ele in mesh
    
        # Local stiffness matrix
        Ke = Local_K(mesh,ele) 

        # Determina quais são os gls GLOBAIS que são "acessados"
        # por esse elemento
        gls = DOFs(mesh,ele) 

        # If truss
        Keg = flag_truss ? To_global(Ke,mesh,ele) :  Ke
        
        # Adiciona a matriz do elemento (rotacionada) à matriz Global
        K[gls,gls] .= K[gls,gls] .+ Keg*(x[ele]^p)

    end #ele

    # Retorna a matriz global
    return Symmetric(K)

end



 function Global_M(mesh::Mesh; x=Float64[])

    # Alias
    bmesh = mesh.bmesh
    ne = bmesh.ne

    # Basic test for x
    if isempty(x)
        x = ones(ne)
    end
    
    # Vamos reforçar a ideia de que x tem a dimensão ne
    @assert length(x)==bmesh.ne "Global_M::x deve ter dimensão igual o número de elementos"

    # Dimensão do problema
    dim = Get_dim(mesh)

    # Tipo de elemento
    etype = Get_etype(mesh)

    # Primeira coisa é alocar a matriz M 
    ng = dim*bmesh.nn
    M = spzeros(ng,ng)

    # Flag se for barras
    flag_truss = Get_eclass(mesh)==:truss

    # Loop pelos elementos, calculando a matriz local Me de cada um
    # e posicionando na M
    for ele in mesh

        # Local mass matrix
        Me = Local_M(mesh,ele)
           
        # Determina quais são os gls GLOBAIS que são "acessados"
        # por esse elemento
        gls = DOFs(mesh,ele) 

        # If truss
        Meg = flag_truss ? To_global(Me,mesh,ele) :  Me
        
        # Adiciona a matriz do elemento (rotacionada) a matriz Global
        M[gls,gls] .= M[gls,gls] .+ Meg*(x[ele])

    end #ele

    # Retorna a matriz global
    return Symmetric(M)

end

