################################################################################
#####                           GMSH 2.1                                 ######
################################################################################

import BMesh:Lgmsh_export_init
"""
(Overloaded from BMesh)
Initialize a gmsh mesh for post processing.

    Lgmsh_export_init(nome_arquivo::String,mesh::Mesh)

"""
function Gmsh_init(nome_arquivo::String,mesh::Mesh)
         BMesh.Lgmsh_export_init(nome_arquivo,mesh.bmesh)
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
        error("ERROR::Adiciona_Vista_Nodal_Vetorial_Gmsh:: vetor deve ter dimensao $dim_total")
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
# O array que contem os valores nodais deve ter dimensao nnos e o número de 
# colunas depende do tipo de elemento
#
"""
Export a tensorial element view to gmsh -> to be used when center=true 

    Gmsh_nodal_stress(mesh::Mesh,stress::Matrix,nome_arquivo::String,
                      nome_vista::String,tempo=0.0)

"""
function Gmsh_nodal_stress(mesh::Mesh,stress::Matrix,nome_arquivo::String,
                           nome_vista::String,tempo=0.0)


    # Alias
    nn = Get_nn(mesh)

    # Tenta abrir o arquivo para append
    saida = try
    open(nome_arquivo,"a")
    catch
        error("ERROR::Gmsh_nodal_stress:: Cannot open $nome_arquivo")
    end


    # Verifica se a dimensao esta correta
    if size(stress,1)!=nn
        error("ERROR::Gmsh_nodal_stress:: number of rows must be nnos")
    end

    # Now we have to proccess the input. Each element type has different number of
    # stress components. The maximum number is 6
    # components are xx yy zz xz yz zy 
    outs = zeros(nn,6)

    etype = Get_etype(mesh)
    if etype===:truss2D || etype===:truss3D

        outs[:,1].=stress[:]

    elseif etype===:solid2D    

        # Normal
        outs[:,1:2].=stress[:,1:2]

        # Shear xy
        outs[:,6].= stress[:,3]
    else
        outs = stress    
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
    println(saida,"9")
    println(saida,nn)
    for i=1:nn
        comp =outs[i,:]
        println(saida,i," ",comp[1], " ", comp[6], " ", comp[5], " ", comp[6], " ",comp[2]," ",comp[4], " ",comp[5], " ",comp[4]," ",comp[3])
    end
    println(saida,"\$EndNodeData")

    # Fecha por hora
    close(saida)

end

#
# Adiciona uma vista tensorial a um arquivo (que ja deve ter o cabecalho)
# O array que contem os valores centroidais deve ter dimensao nelems e o número de 
# colunas depende do tipo de elemento
#
"""
Export a tensorial element view to gmsh -> to be used when center=true 

    Gmsh_element_stress(mesh::Mesh,stress::Matrix,nome_arquivo::String,
                        nome_vista::String,tempo=0.0)

"""

function Gmsh_element_stress(mesh::Mesh,stress::Matrix,nome_arquivo::String,
                             nome_vista::String,tempo=0.0)



    
    #
    # Hint for the user
    #
    if isa(mesh,Mesh2D) && size(stress,2)==4*3
        trhow("Gmesh_element_stress:: this subroutine considers that stresses were evaluated with center=true. Use Gmsh_element_stresses.")
    end

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
        error("ERROR::Gmesh_element_stress:: entrada deve ter dimensao nnos")
    end

    # Now we have to proccess the input. Each element type has different number of
    # stress components. The maximum number is 6
    # components are xx yy zz xz yz zy 
    outs = zeros(nelems,6)

    etype = Get_etype(mesh)
    if etype===:truss2D || etype===:truss3D

        outs[:,1].=stress[:]

    elseif etype===:solid2D    

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


########################################################################################################################################

#
# Tensoes é uma matrizinha com 4 linhas (pto de gauss) por 3 colunas (tensão)
#

function Map_stress2nodes_Quad(tensoes::Matrix,processado::Bool=false)

    # Esta matriz foi obtida no maxima e os cálculos estão
    # na pasta de documentação.

    c1 = 1.866025403784438
    c2 = -0.5
    c3 = 0.1339745962155612

    invA  = [   c1 c2 c3 c2 ;
                c2 c1 c2 c3 ;
                c3 c2 c1 c2 ;
                c2 c3 c2 c1 ]

    # e calcula os valores nodais - se processado for true, usa uma cópia direta das tensões
    S = copy(tensoes)
    if !processado
       S .= invA*tensoes
    end

    # Já que estamos nos ocupando da visualização,
    # vamos montar a string para o gmsh
    # Cada nó vai ter que conter uma sequência como a
    # comp[1], " ", comp[3], " ", 0.0, " ", comp[3], " ",comp[2]," ",0.0, " ",0.0, " ",0.0," ",0.0
    saida = " "
    for i=1:4
        saida = saida * " " * string(S[i,1]) * " " * string(S[i,3]) * " 0.0 " * string(S[i,3]) * " " * string(S[i,2]) *  " 0.0 " *  " 0.0 " * " 0.0 "  * " 0.0 "
    end

    return saida
end


# Se processado for true, então os valores já são os nodais. Do contrário,
# extrapolamos os valores para os nós intermente.
function Gmsh_element_stresses(mesh::Mesh,stresses::Matrix,nome_arquivo::String,
                               nome_vista::String,tempo=0.0,processado::Bool=true)


    isa(mesh,Mesh2D) || throw("Gmsh_element_stresses:: 3D not implemented yet..")                           

    # Alias
    ne = mesh.bmesh.ne

    # Tenta abrir o arquivo para append
    saida = try
                open(nome_arquivo,"a")
    catch
        error("ERROR::Gmsh_element_stresses:: $(nome_arquivo) not found")
    end


    # Verifica se a dimensao esta correta
    if size(stresses,1) != ne
        error("ERROR::Gmsh_element_stresses:: number of rows must be $ne")
    end

    # Verifica se a dimensao esta correta
    if size(stresses,2) != 4*3
        error("ERROR::Gmsh_element_stresses:: number of cols must be 12")
    end

    #
    #
    println(saida,"\$ElementNodeData")
    println(saida,"1")
    println(saida,"\" $nome_vista\"")
    println(saida,"1")
    println(saida,tempo)
    println(saida,"3")
    println(saida,"0")
    println(saida,"9")
    println(saida,ne)
    for i=1:ne

            
            tensoes_elemento = vcat(transpose(stresses[i,1:3]),transpose(stresses[i,4:6]),
                                  transpose(stresses[i,7:9]),transpose(stresses[i,10:12]))

            texto = Map_stress2nodes_Quad(tensoes_elemento,processado)
            println(saida,i," 4 ", texto)

    end
    println(saida,"\$EndElementNodeData")

    # Fecha por hora
    close(saida)

end