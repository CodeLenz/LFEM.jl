################################################################################
#####                           GMSH 2.1                                 ######
################################################################################

"""
(Overloaded from BMesh)
Initialize a gmsh mesh for post processing.

    Gmsh_init(nome_arquivo::String,mesh::Mesh)

"""
import BMesh:Gmsh_init
function Gmsh_init(nome_arquivo::String,mesh::Mesh)
         BMesh.Gmsh_init(nome_arquivo,mesh.bmesh)
end



"""
Export a nodal scalar view to gmsh

    Gmsh_nodal_scalar(mesh::Mesh,escalares::Vector,nome_arquivo::String,
                      nome_vista::String,tempo=0.0)

"""
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

"""
Export an element (centroidal) scalar view to gmsh

    Gmsh_element_scalar(mesh::Mesh,escalares::Vector,nome_arquivo::String,
                        nome_vista::String,tempo=0.0)

"""
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
    for i in mesh
        println(saida,i," ",escalares[i])#,digits=15))
    end
    println(saida,"\$EndElementData")

    # Fecha por hora
    close(saida)

end


"""
Export an nodal vectorial view to gmsh

    Gmsh_nodal_vector(mesh::Mesh,vetor::Vector,nome_arquivo::String,
                      nome_vista::String,tempo=0.0)

"""
function Gmsh_nodal_vector(mesh::Mesh,vetor::Vector,nome_arquivo::String,
                           nome_vista::String,tempo=0.0)


    # Tenta abrir o arquivo para append
    saida = try
                open(nome_arquivo,"a")
    catch
        error("ERROR::Adiciona_Vista_Nodal_Vetorial_Gmsh:: Nao foi possivel acessar $nome_arquivo")
    end

    dim = Get_dim(mesh)

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


#
# Adiciona uma vista tensorial a um arquivo (que ja deve ter o cabecalho)
# O vetor que contem os valores centroidais deve ter dimensao nelems e o número de 
# colunas depende do tipo de elemento
#
"""
Export a tensorial element view to gmsh

    Gmsh_element_stress(mesh::Mesh,stress::Matrix,nome_arquivo::String,
                        nome_vista::String,tempo=0.0)

"""

function Gmsh_element_stress(mesh::Mesh,stress::Matrix,nome_arquivo::String,
                             nome_vista::String,tempo=0.0)


    # Alias
    nelems = Get_ne(mesh)

    # Tenta abrir o arquivo para append
    saida = try
    open(nome_arquivo,"a")
    catch
        error("ERROR::Gmsh_element_stress:: Nao foi possivel acessar $nome_arquivo")
    end


    # Verifica se a dimensao esta correta
    if size(stress,1)!=nelems
        error("ERROR::Gmesh_element_stress:: vetor com escalares deve ter dimensao nnos")
    end

    # Now we have to proccess the input. Each element type has different number of
    # stress components. The maximum number is 6
    # components are xx yy zz xz yz zy 
    outs = zeros(nelems,6)

    etype = Get_etype(mesh)
    if etype==:truss2D || etype==:truss3D

        outs[:,1].=stress[:]

    elseif etype==:solid2D    

        # Normal
        outs[:,1:2].=stress[:,1:2]

        # Shear xy
        outs[:,6].= stress[:,3]
    else
        outs = stress    
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
    println(saida,"9")
    println(saida,nelems)
    for i in mesh
        comp =outs[i,:]
        println(saida,i," ",comp[1], " ", comp[6], " ", comp[5], " ", comp[6], " ",comp[2]," ",comp[4], " ",comp[5], " ",comp[4]," ",comp[3])
    end
    println(saida,"\$EndElementData")

    # Fecha por hora
    close(saida)

end
