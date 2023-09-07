@testset "Solve_modal" begin
    
    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)

    # Solve the modal problem
    @test try  Solve_modal(mesh)
            true
    catch err  
            false
    end
    @isinferred Solve_modal(mesh)


    # Load Simply supported 2D from TMeshes (solid)
    mesh = Simply_supported2D(6,6,:solid2D)

    # Solve the modal problem
    @test try  Solve_modal(mesh,lumped=false)
            true
    catch err  
            false
    end
    @isinferred Solve_modal(mesh,lumped=false)

    # Load Simply supported 3D from TMeshes (truss)
    mesh = Simply_supported3D(2,2,2)

    # Solve the modal problem
    @test try  Solve_modal(mesh)
        true
    catch err  
        false   
    end
    @isinferred Solve_modal(mesh)

    # Load Simply supported 3D from TMeshes (solid)
    mesh = Simply_supported3D(2,2,2,:solid3D)

    # Solve the modal problem
    @test try  Solve_modal(mesh,lumped=false)
        true
    catch err  
        false
    end
    @isinferred Solve_modal(mesh,lumped=false)

end

@testset "Solve_modal (parametrization)" begin


    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)

    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    kparam(xe::Float64,p=3.0)= xe^p
    mparam(xe::Float64,p=2.0,cut=0.1)= ifelse(xe>=cut,xe,xe^p)
    
    # Solve the modal problem
    @test try  Solve_modal(mesh,x,kparam,mparam)
        true
    catch err  
        false
    end
    @isinferred Solve_modal(mesh,x,kparam,mparam)
    
    # Load Simply supported 2D from TMeshes (solid)
    mesh = Simply_supported2D(6,6,:solid2D)
    
    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    
    # Solve the modal problem
    @test try  Solve_modal(mesh,x,kparam,mparam,lumped=false)
        true
    catch err  
        false
    end
    @isinferred Solve_modal(mesh,x,kparam,mparam,lumped=false)
        
    # Load Simply supported 3D from TMeshes (truss)
    mesh = Simply_supported3D(2,2,2)
    
    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    
    # Solve the modal problem
    @test try  Solve_modal(mesh,x,kparam,mparam)
        true
    catch err  
        false
    end
    @isinferred Solve_modal(mesh,x,kparam,mparam)
    
    # Load Simply supported 3D from TMeshes (solid)
    mesh = Simply_supported3D(2,2,2,:solid3D)
    
    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    
    # Solve the modal problem
    @test try  Solve_modal(mesh,x,kparam,mparam,lumped=false)
        true
    catch err  
        false
    end
    @isinferred Solve_modal(mesh,x,kparam,mparam,lumped=false)

end