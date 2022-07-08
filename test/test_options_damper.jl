@testset "Options=:Damper" begin

    #                        Valid tests 2D 

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)

    # Angular frequency
    w  = 20.0

    # Structural damping
    α_c = 0.0
    β_c = 1E-6

    # Damper
    mesh.options[:Damper] = [2 1 10.0;
                             3 2 10.0]


    # Solve the harmonic problem
    @test try  Solve_harmonic(mesh,w,α_c,β_c)
            true
    catch err  
            false
    end
    @isinferred Solve_harmonic(mesh,w,α_c,β_c)

    #                        Valid tests 3D 

    # Load Simply supported 3D from TMeshes
    mesh = Simply_supported3D(6,6,6)
    
    # Damper
    mesh.options[:Damper] = [2 1 10.0;
                             3 2 10.0;
                             3 3 10.0]

    # Solve the harmonic problem
    @test try  Solve_harmonic(mesh,w,α_c,β_c)
        true
    catch err  
        false
    end
    @isinferred Solve_harmonic(mesh,w,α_c,β_c)

    ############################ THROWS #####################

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)

    # Damper # Invalid entries (more columns)
    mesh.options[:Damper] = [2 1 10.0 10 10 10]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)

    # Damper # Invalid entries (few columns)
    mesh.options[:Damper] = [2 1]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)

    # Damper - invalid nodes
    mesh.options[:Damper] = [0 2 10.0]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)

    # Damper - invalid nodes
    mesh.options[:Damper] = [-1 2 10.0]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)

    # Damper - invalid dof
    mesh.options[:Damper] = [1 0 10.0]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)
    
    # Damper - invalid dof
    mesh.options[:Damper] = [1 -1 10.0]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)
    
    # Damper - invalid value
    mesh.options[:Damper] = [1 2 -10.0]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)

    ################# 3D ################

    # Load Simply supported 3D from TMeshes
    mesh = Simply_supported3D(6,6,6)

    # Damper - invalid nodes
    mesh.options[:Damper] = [0 2 10.0]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)
    
    # Damper - invalid nodes
    mesh.options[:Damper] = [-1 2 10.0]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)

    # Damper - invalid dof
    mesh.options[:Damper] = [1 0 10.0]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)
    
    # Damper - invalid dof
    mesh.options[:Damper] = [0 -1 10.0]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)

    # Damper - invalid value
    mesh.options[:Damper] = [1 2 -10.0]
    @test_throws String Solve_harmonic(mesh,w,α_c,β_c)
    
end

