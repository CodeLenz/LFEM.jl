@testset "Solve_linear" begin
    
    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)

    # Solve the linear static equilibrium
    @test try  Solve_linear(mesh)
            true
    catch err  
            false
    end
    @isinferred LFEM.Solve_linear(mesh)


    # Load Simply supported 2D from TMeshes (solid)
    mesh = Simply_supported2D(6,6,:solid2D)

    # Solve the linear static equilibrium
    @test try  Solve_linear(mesh)
            true
    catch err  
            false
    end
    @isinferred Solve_linear(mesh)

    # Load Simply supported 3D from TMeshes (truss)
    mesh = Simply_supported3D(2,2,2)

    # Solve the linear static equilibrium
    @test try  Solve_linear(mesh)
        true
    catch err  
        false   
    end
    @isinferred Solve_linear(mesh)

    # Load Simply supported 3D from TMeshes (solid)
    mesh = Simply_supported3D(2,2,2,:solid3D)

    # Solve the linear static equilibrium
    @test try  Solve_linear(mesh)
        true
    catch err  
        false
    end
    @isinferred Solve_linear(mesh)

end

@testset "Solve_linear (parametrization)" begin


    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)

    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    kparam(xe::Float64,p=3.0)= xe^p
    
    # Solve the linear static equilibrium
    @test try  Solve_linear(mesh,x,kparam)
        true
    catch err  
        false
    end
    @isinferred Solve_linear(mesh,x,kparam)
    
    # Load Simply supported 2D from TMeshes (solid)
    mesh = Simply_supported2D(6,6,:solid2D)
    
    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    
    # Solve the linear static equilibrium
    @test try  Solve_linear(mesh,x,kparam)
        true
    catch err  
        false
    end
    @isinferred Solve_linear(mesh,x,kparam)
        
    # Load Simply supported 3D from TMeshes (truss)
    mesh = Simply_supported3D(2,2,2)
    
    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    
    # Solve the linear static equilibrium
    @test try  Solve_linear(mesh,x,kparam)
        true
    catch err  
        false
    end
    @isinferred Solve_linear(mesh,x,kparam)
    
    # Load Simply supported 3D from TMeshes (solid)
    mesh = Simply_supported3D(2,2,2,:solid3D)
    
    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    
    # Solve the linear static equilibrium
    @test try  Solve_linear(mesh,x,kparam)
        true
    catch err  
        false
    end
    @isinferred Solve_linear(mesh,x,kparam)

end