@testset "Truss3D" begin

    ####################   3D BMesh ######################
    #
    # Single bar in tension (L=0.1m F=100.0 A=1.0 E=100.0)
    #
    #    Valid inputs (no error) 
    #
    nn = 2
    ne = 1
    coord = [0.0 0.0 0.0;
             0.1 0.0 0.0 ]
    connect = [1 2 ]
    Lx = 1.0
    Ly = 1.0
    Lz = 1.0
    nx = 1
    ny = 1
    nz = 1
    etype = :truss3D

    # truss3D
    b3 = Bmesh3D(etype,nn,ne,coord,connect,Lx,Ly,Lz,nx,ny,nz)

    # Essential boundary conditions
    nebc = 5
    ebc = [1 1 0.0 ; 
           1 2 0.0 ;  
           1 3 0.0 ;
           2 2 0.0 ; 
           2 3 0.0 ]

    # Natural boundary conditions
    nbc = [2 1 100.0]

    # Material and geometry
    materials = [Material(Ex=100.0,density=1.0)]
    geometries = [Geometry(A=1.0)]
           
    # Mesh
    m3 = Mesh3D(b3,materials,geometries,ebc,nbc)

    # FIRST TEST 
    # Local stiffness matrix
    K = Local_K(m3,1)

    # Reference
    cte = 100.0*1.0/0.1
    refer = cte*[1.0 0.0 0.0 -1.0 0.0 0.0;
                 0.0 0.0 0.0  0.0 0.0 0.0;
                 0.0 0.0 0.0  0.0 0.0 0.0
                -1.0 0.0 0.0  1.0 0.0 0.0;
                 0.0 0.0 0.0  0.0 0.0 0.0;
                 0.0 0.0 0.0  0.0 0.0 0.0]

    @test all(K.==refer)

    # SECOND TEST 
    # B matrix
    B = B_truss3D(m3,1)

    # Rerefence
    refer =[-1/0.1 0.0 0.0 1/0.1 0.0 0.0]

    @test all(B.==refer)

    # THIRD TEST
    # stress
    U = [0.0;0.0;0.0;0.1;0.0;0.0]

    s = Stress_truss3D(m3,1,U)[1]

    @test s==100.0
    
    # Test the driver 
    s = Stress(m3,1,U) 
    @test all(s.==[100.0])

    
    # Test local mass matrix
    M = Local_M(m3,1)

    # Reference
    cte = (1.0*1.0*0.1)/2
    refer = cte*[1.0 0.0 0.0  0.0 0.0 0.0;
                 0.0 1.0 0.0  0.0 0.0 0.0;
                 0.0 0.0 1.0  0.0 0.0 0.0;
                 0.0 0.0 0.0  1.0 0.0 0.0;
                 0.0 0.0 0.0  0.0 1.0 0.0;
                 0.0 0.0 0.0  0.0 0.0 1.0]

    @test all(M.==refer)


end
