@testset "Options=:Mass" begin

    #                        Valid tests 2D 

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)
    
    # Mass
    mesh.options[:Mass] = [2 1 10.0;
                           3 2 10.0]

    # Solve the linear static equilibrium
    @test try  Solve_modal(mesh)
            true
    catch err  
            false
    end
    @isinferred Solve_modal(mesh)

    #                        Valid tests 3D 

    # Load Simply supported 3D from TMeshes
    mesh = Simply_supported3D(6,6,6)
    
    # Mass
    mesh.options[:Mass] = [2 1 10.0;
                           3 2 10.0;
                           4 3 10.0]

    # Solve the linear static equilibrium
    @test try  Solve_modal(mesh)
            true
    catch err  
            false
    end
    @isinferred Solve_modal(mesh)

    ############################ THROWS #####################

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)

    # Invalid entry
    mesh.options[:Mass] = [1 2 1E5 10.0 20.0]
    @test_throws String Solve_modal(mesh)

    # Invalid entry
    mesh.options[:Mass] = [1 2]
    @test_throws String Solve_modal(mesh)

    # Mass - invalid nodes
    mesh.options[:Mass] = [0 2 10.0]
    @test_throws String Solve_modal(mesh)

    # Mass - invalid nodes
    mesh.options[:Mass] = [-1 2 10.0]    
    @test_throws String Solve_modal(mesh)

    # Mass - invalid nodes
    mesh.options[:Mass] = [Get_nn(mesh)+1 2 10.0]
    @test_throws String Solve_modal(mesh)

    # Mass - invalid dof
    mesh.options[:Mass] = [3 0 10.0]
    @test_throws String Solve_modal(mesh)

    # Mass - invalid dof (its 2D)
    mesh.options[:Mass] = [3 3 10.0]
    @test_throws String Solve_modal(mesh)

    # Mass - invalid value
    mesh.options[:Mass] = [3 2 -10.0]
    @test_throws String Solve_modal(mesh)

    ################# 3D ################

    # Load Simply supported 3D from TMeshes
    mesh = Simply_supported3D(6,6,6)

    # Mass - invalid nodes
    mesh.options[:Mass] = [0 2 10.0]
    @test_throws String Solve_modal(mesh)
    
    # Mass - invalid nodes
    mesh.options[:Mass] = [-1 2 10.0]    
    @test_throws String Solve_modal(mesh)
    
    # Mass - invalid nodes
    mesh.options[:Mass] = [Get_nn(mesh)+1 2 10.0]
    @test_throws String Solve_modal(mesh)
    
    # Mass - invalid dof
    mesh.options[:Mass] = [3 0 10.0]
    @test_throws String Solve_modal(mesh)
    
    # Mass - invalid dof (its 3D)
    mesh.options[:Mass] = [3 4 10.0]
    @test_throws String Solve_modal(mesh)
    
    # Mass - invalid value
    mesh.options[:Mass] = [3 2 -1E3]
    @test_throws String Solve_modal(mesh)
    
end

