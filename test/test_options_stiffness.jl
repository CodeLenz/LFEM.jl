@testset "Options=:Stiffness" begin

    #                        Valid tests 2D 

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)
    
    # Stiffness
    mesh.options[:Stiffness] = [2 1 1E5;
                                3 2 1E5]

    # Solve the linear static equilibrium
    @test try  Solve_linear(mesh)
            true
    catch err  
            false
    end
    @isinferred Solve_linear(mesh)

    #                        Valid tests 3D 

    # Load Simply supported 3D from TMeshes
    mesh = Simply_supported3D(6,6,6)
    
    # Stiffness
    mesh.options[:Stiffness] = [2 1 1E5;
                                3 2 1E5;
                                4 3 1E5]

    # Solve the linear static equilibrium
    @test try  Solve_linear(mesh)
            true
    catch err  
            false
    end
    @isinferred Solve_linear(mesh)

    ############################ THROWS #####################

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)

    # Invalid entry
    mesh.options[:Stiffness] = [0 2 1E5 10.0 20.0]
    @test_throws String Solve_linear(mesh)

    # Stiffness - invalid nodes
    mesh.options[:Stiffness] = [0 2 1E5]
    @test_throws String Solve_linear(mesh)

    # Stiffness - invalid nodes
    mesh.options[:Stiffness] = [-1 2 1E5]    
    @test_throws String Solve_linear(mesh)

    # Stiffness - invalid nodes
    mesh.options[:Stiffness] = [Get_nn(mesh)+1 2 1E5]
    @test_throws String Solve_linear(mesh)

    # Stiffness - invalid dof
    mesh.options[:Stiffness] = [3 0 1E5]
    @test_throws String Solve_linear(mesh)

    # Stiffness - invalid dof (its 2D)
    mesh.options[:Stiffness] = [3 3 1E5]
    @test_throws String Solve_linear(mesh)

    # Stiffness - invalid value
    mesh.options[:Stiffness] = [3 2 -1E3]
    @test_throws String Solve_linear(mesh)

    ################# 3D ################

    # Load Simply supported 3D from TMeshes
    mesh = Simply_supported3D(6,6,6)

    # Stiffness - invalid nodes
    mesh.options[:Stiffness] = [0 2 1E5]
    @test_throws String Solve_linear(mesh)
    
    # Stiffness - invalid nodes
    mesh.options[:Stiffness] = [-1 2 1E5]    
    @test_throws String Solve_linear(mesh)
    
    # Stiffness - invalid nodes
    mesh.options[:Stiffness] = [Get_nn(mesh)+1 2 1E5]
    @test_throws String Solve_linear(mesh)
    
    # Stiffness - invalid dof
    mesh.options[:Stiffness] = [3 0 1E5]
    @test_throws String Solve_linear(mesh)
    
    # Stiffness - invalid dof (its 3D)
    mesh.options[:Stiffness] = [3 4 1E5]
    @test_throws String Solve_linear(mesh)
    
    # Stiffness - invalid value
    mesh.options[:Stiffness] = [3 2 -1E3]
    @test_throws String Solve_linear(mesh)
    
end