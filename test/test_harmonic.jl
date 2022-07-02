@testset "Solve_harmonic" begin

    # Angular frequency [rad/s]
    w = 20.0
    
    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)

    # Solve the harmonic problem
    @test try  Solve_harmonic(mesh,w)
            true
    catch err  
            false
    end
    @isinferred Solve_harmonic(mesh,w)


    # Load Simply supported 2D from TMeshes (solid)
    mesh = Simply_supported2D(6,6,:solid2D)

    # Solve the harmonic problem
    @test try  Solve_harmonic(mesh,w)
            true
    catch err  
            false
    end
    @isinferred Solve_harmonic(mesh,w)

    # Load Simply supported 3D from TMeshes (truss)
    mesh = Simply_supported3D(2,2,2)

    # Solve the harmonic problem
    @test try  Solve_harmonic(mesh,w)
        true
    catch err  
        false   
    end
    @isinferred Solve_harmonic(mesh,w)

    # Load Simply supported 3D from TMeshes (solid)
    mesh = Simply_supported3D(2,2,2,:solid3D)

    # Solve the harmonic problem
    @test try  Solve_harmonic(mesh,w)
        true
    catch err  
        false
    end
    @isinferred Solve_harmonic(mesh,w)

end

@testset "Solve_harmonic (parametrization)" begin

    # Angular frequency [rad/s]
    w = 20.0

    # Load Simply supported 2D from TMeshes
    mesh = Simply_supported2D(6,6)

    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    kparam(xe::Float64,p=3.0)= xe^p
    mparam(xe::Float64,p=2.0,cut=0.1)= ifelse(xe>=cut,xe,xe^p)
    
    # Solve the harmonic problem
    @test try  Solve_harmonic(mesh,w,x,kparam,mparam)
        true
    catch err  
        false
    end
    @isinferred Solve_harmonic(mesh,w,x,kparam,mparam)
    
    # Load Simply supported 2D from TMeshes (solid)
    mesh = Simply_supported2D(6,6,:solid2D)
    
    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    
    # Solve the harmonic problem
    @test try  Solve_harmonic(mesh,w,x,kparam,mparam)
        true
    catch err  
        false
    end
    @isinferred Solve_harmonic(mesh,w,x,kparam,mparam)
        
    # Load Simply supported 3D from TMeshes (truss)
    mesh = Simply_supported3D(2,2,2)
    
    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    
    # Solve the harmonic problem
    @test try  Solve_harmonic(mesh,w,x,kparam,mparam)
        true
    catch err  
        false
    end
    @isinferred Solve_harmonic(mesh,w,x,kparam,mparam)
    
    # Load Simply supported 3D from TMeshes (solid)
    mesh = Simply_supported3D(2,2,2,:solid3D)
    
    # Define x, kparam and mparam
    x = ones(Get_ne(mesh))
    
    # Solve the harmonic problem
    @test try  Solve_harmonic(mesh,w,x,kparam,mparam)
        true
    catch err  
        false
    end
    @isinferred Solve_harmonic(mesh,w,x,kparam,mparam)

end