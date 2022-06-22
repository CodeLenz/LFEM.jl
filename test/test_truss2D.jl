@testset "Truss2D" begin

    ####################  2D BMesh ######################
    # Single bar in tension (L=0.1m F=100.0 A=1.0 E=100.0)
    #
    #    Valid inputs (no error) 
    #
    nn = 2
    ne = 1
    coord = [0.0 0.0 ;
             0.1 0.0 ]
    connect = [1 2 ]
    Lx = 1.0
    Ly = 1.0
    nx = 1
    ny = 1
    etype = :truss2D

    # truss2D
    b2 = Bmesh2D(etype,nn,ne,coord,connect,Lx,Ly,nx,ny)

    # Essential boundary conditions
    ebc = [1 1 0.0 ; 
           1 2 0.0 ; 
           2 2 0.0 ]

    # Natural boundary conditions
    nbc = [2 1 100.0]

    # Material and geometry
    materials = [Material(Ex=100.0)]
    geometries = [Geometry(A=1.0)]
           
    # Mesh
    m2 = Mesh2D(b2,materials,geometries,ebc,nbc)


    # FIRST TEST 
    # Local stiffness matrix
    K = K_truss2D(100.0,1.0,0.1)

    # Reference
    cte = 100.0*1.0/0.1
    refer = cte*[1.0 0.0 -1.0 0.0 ;
                 0.0 0.0  0.0 0.0 ;
                -1.0 0.0  1.0 0.0 ;
                 0.0 0.0  0.0 0.0]

    @test all(K.==refer)

    # SECOND TEST 
    # B matrix
    B = B_truss2D(0.1)

    # Rerefence
    refer =[-1/0.1 0.0 1/0.1 0.0]

    @test all(B.==refer)

    # THIRD TEST
    # stress
    U = [0.0;0.0;0.1;0.0]

    s = Stress_truss2D(m2,1,U)[1]

    @test s==100.0

   # Test the driver 
   s = Stress(m2,1,U) 
    
   @test s==[100.0] 

end
