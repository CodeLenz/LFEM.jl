################################################################################
#####                           GMSH 2.1                                 ######
################################################################################
#
# Cria o cabecalho com informacoes da malha
# para posterior adicao de vistas com saidas
#
function Gmsh_init(nome_arquivo::String,mesh::Mesh)

    # Verifica se já existe o arquivo, se sim, remove
    if isfile(nome_arquivo); rm(nome_arquivo); end

    # Abre o arquivo para escrita
    saida = open(nome_arquivo,"a")

    dim = 2
    if isa(mesh,Mesh3D)
        dim = 3
    end

    # Cabecalho do gmsh
    println(saida,"\$MeshFormat")
    println(saida,"2.2 0 8")
    println(saida,"\$EndMeshFormat")

    # Nodes
    println(saida,"\$Nodes")
    println(saida,mesh.bmesh.nn)
    if dim==2
        for i=1:mesh.bmesh.nn
            println(saida,i," ",mesh.bmesh.coord[i,1]," ",mesh.bmesh.coord[i,2]," 0.0 ")
        end
    else
        for i=1:mesh.bmesh.nn
            println(saida,i," ",mesh.bmesh.coord[i,1]," ",mesh.bmesh.coord[i,2]," ",mesh.bmesh.coord[i,3])
        end
    end    
    println(saida,"\$EndNodes")

    # Conectividades
    if dim==2
       mesh.bmesh.etype==:truss2D || throw("Gmsh_init::invalid element type for Mesh2D")
    else
       mesh.bmesh.etype==:truss3D || throw("Gmsh_init::invalid element type for Mesh3D")
    end   
    tipo_elemento = 1

    println(saida,"\$Elements")
    println(saida,mesh.bmesh.ne)
    for i=1:mesh.bmesh.ne
        con = string(i)*" "*string(tipo_elemento)*" 0 "*string(mesh.bmesh.connect[i,1])
        for j=2:size(mesh.bmesh.connect,2)
            con = con * " " * string(mesh.bmesh.connect[i,j])
        end
        println(saida,con)
    end
    println(saida,"\$EndElements")

    # Fecha o arquivo ... por hora
    close(saida)


end # Gera_Malha_Gmsh


#
# Adiciona uma vista escalar nodal a um arquivo (que ja deve ter o cabecalho)
# O vetor que contem os valores nodais deve ter dimensao nnos
#
function Gmsh_nodal_scalar(mesh::Mesh,escalares::Vector,nome_arquivo::String,
                          nome_vista::String,tempo=0.0)


    # Tenta abrir o arquivo para append
    if isfile(nome_arquivo)==false
            error("ERROR::Adiciona_Vista_Nodal_Escalar_Gmsh:: Nao foi possivel acessar $nome_arquivo")
    end

    saida = open(nome_arquivo,"a")


    # Verifica se a dimensao esta correta
    if size(escalares,1)!=mesh.bmesh.nn
        error("ERROR::Adiciona_Vista_Nodal_Escalar_Gmsh:: vetor com escalares deve ter dimensao nnos")
    end

    #
    #
    println(saida,"\$NodeData")
    println(saida,"1")
    println(saida,"\" $nome_vista \"")
    println(saida,"1")
    println(saida,tempo)
    println(saida,"3")
    println(saida,"0")
    println(saida,"1")
    println(saida,mesh.bmesh.nn)
    for i=1:mesh.bmesh.nn
        println(saida,i," ",escalares[i])
    end
    println(saida,"\$EndNodeData")

    # Fecha por hora
    close(saida)

end

#
# Adiciona uma vista escalar a um arquivo (que ja deve ter o cabecalho)
# O vetor que contem os valores centroidais deve ter dimensao nelems
#
function Gmsh_element_scalar(mesh::Mesh,escalares::Vector,nome_arquivo::String,
                             nome_vista::String,tempo=0.0)


    # Alias
    nelems = mesh.bmesh.ne

    # Tenta abrir o arquivo para append
    saida = try
                open(nome_arquivo,"a")
    catch
        error("ERROR::Adiciona_Vista_Nodal_Escalar_Gmsh:: Nao foi possivel acessar $nome_arquivo")
    end


    # Verifica se a dimensao esta correta
    if size(escalares,1)!=nelems
        error("ERROR::Adiciona_Vista_Nodal_Escalar_Gmsh:: vetor com escalares deve ter dimensao nnos")
    end

    #
    #
    println(saida,"\$ElementData")
    println(saida,"1")
    println(saida,"\" $nome_vista \"")
    println(saida,"1")
    println(saida,tempo)
    println(saida,"3")
    println(saida,"0")
    println(saida,"1")
    println(saida,nelems)
    for i=1:nelems
        println(saida,i," ",escalares[i])#,digits=15))
    end
    println(saida,"\$EndElementData")

    # Fecha por hora
    close(saida)

end


#
# Adiciona uma vista vetorial nodal a um arquivo (que ja deve ter o cabecalho)
# O vetor que contem os valores nodais deve ser expandido, ou seja,
# deve ter sido criado por
#
function Gmsh_nodal_vector(mesh::Mesh,vetor::Vector,nome_arquivo::String,
                           nome_vista::String,tempo=0.0)


    # Tenta abrir o arquivo para append
    saida = try
                open(nome_arquivo,"a")
    catch
        error("ERROR::Adiciona_Vista_Nodal_Vetorial_Gmsh:: Nao foi possivel acessar $nome_arquivo")
    end

    dim = 2
    if isa(mesh,Mesh3D)
        dim=3
    end

    # Verifica se a dimensao esta correta
    dim_total = dim*mesh.bmesh.nn
    if size(vetor,1)!= dim_total
        error("ERROR::Adiciona_Vista_Nodal_Vetorial_Gmsh:: vetor com escalares deve ter dimensao $dim_total")
    end

    #
    #
    println(saida,"\$NodeData")
    println(saida,"1")
    println(saida,"\" $nome_vista\"")
    println(saida,"1")
    println(saida,tempo)
    println(saida,"3")
    println(saida,"0")
    println(saida,"3")
    println(saida,mesh.bmesh.nn)
    for no=1:mesh.bmesh.nn
        pos1 = dim*(no-1)+1; val1 = vetor[pos1]
        pos2 = dim*(no-1)+2; val2 = vetor[pos2]
        val3 = 0.0
        if dim==3
            pos3 = dim*(no-1)+3; val3 = vetor[pos3]
        end 
        println(saida,no," ",val1," ",val2," ",val3 )
    end
    println(saida,"\$EndNodeData")

    # Fecha por hora
    close(saida)

end
