@testset "Options=:IS_TOPO" begin

    #                        Valid tests 2D 

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(10,10,:solid2D)
    
    # Set option
    mesh.options[:IS_TOPO] = zeros(1,1)

    # Solve using all elements
    U, _ = Solve_linear(mesh)

    # Set option
    mesh.options[:IS_TOPO] = ones(1,1)

    # Solve using only one element
    U2, _ = Solve_linear(mesh)

    @test all(isapprox(U2,U))
     

    #                        Valid tests 3D 

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported3D(4,4,4,:solid3D)
    
    # Unset option
    delete!(mesh.options,:IS_TOPO)

    # Solve using all elements
    U, _ = Solve_linear(mesh)

    # Set option
    mesh.options[:IS_TOPO] = ones(1,1)

    # Solve using only one element
    U2, _ = Solve_linear(mesh)

    @test all(isapprox(U2,U))
    

    #                        Valid tests 2D 

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(10,10,:solid2D)
    
    # Unset option
    delete!(mesh.options,:IS_TOPO)

    # Solve using all elements
    w, _ = Solve_modal(mesh)

    # Set option
    mesh.options[:IS_TOPO] = ones(1,1)

    # Solve using only one element
    w2, _ = Solve_modal(mesh)

    @test all(isapprox(w2,w))
     

    #                        Valid tests 3D 

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported3D(4,4,4,:solid3D)
    
    # Unset option
    delete!(mesh.options,:IS_TOPO)

    # Solve using all elements
    w, _ = Solve_modal(mesh,nev=1)

    # Set option
    mesh.options[:IS_TOPO] = ones(1,1)

    # Solve using only one element
    w2, _ = Solve_modal(mesh,nev=1)

    @test all(isapprox(w2,w))
    



end

