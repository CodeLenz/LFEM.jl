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
    dim = 2
    if isa(mesh,Mesh3D)
        dim=3
    end

    # Tipo de elemento
    etype = bmesh.etype

    # Primeira coisa é alocar a matriz K 
    ng = dim*bmesh.nn
    K = spzeros(ng,ng)

    # Flag se for barras
    flag_truss = contains(string(etype),"truss")

    # Loop pelos elementos, calculando a matriz local Ke de cada um
    # e posicionando na K
    for ele=1:ne

        if flag_truss
           Keg = Keg_truss(mesh,ele)
        else 
           Keg = Keg_solid(mesh,ele)
        end
           
        # Determina quais são os gls GLOBAIS que são "acessados"
        # por esse elemento
        gls = DOFs(bmesh,ele) 

        # Adiciona a matriz do elemento (rotacionada) a matriz Global
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
    dim = 2
    if isa(mesh,Mesh3D)
        dim=3
    end

    # Tipo de elemento
    etype = bmesh.etype

    # Primeira coisa é alocar a matriz M 
    ng = dim*bmesh.nn
    M = spzeros(ng,ng)

    # Flag se for barras
    flag_truss = contains(string(etype),"truss")

    # Loop pelos elementos, calculando a matriz local Me de cada um
    # e posicionando na M
    for ele=1:ne

        if flag_truss
           Meg = Meg_truss(mesh,ele)
        else 
           Meg = Meg_solid2D(mesh,ele)
        end
           
        # Determina quais são os gls GLOBAIS que são "acessados"
        # por esse elemento
        gls = DOFs(bmesh,ele) 

        # Adiciona a matriz do elemento (rotacionada) a matriz Global
        M[gls,gls] .= M[gls,gls] .+ Meg*(x[ele])

    end #ele

    # Retorna a matriz global
    return Symmetric(M)

end

