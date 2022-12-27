@testset "Point load (3D)" begin

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
    nhebc = 5
    hebc = [1 1 ; 
           1 2  ;  
           1 3  ;
           2 2  ; 
           2 3  ]

    # Natural boundary conditions
    nbc = [2 1 100.0]

    # Material and geometry
    materials = [Material(Ex=100.0)]
    geometries = [Geometry(A=1.0)]
           
    # Mesh
    m3 = Mesh3D(b3,materials,geometries,hebc,nbc)

                 
    # Point load
    F = Point_load(m3)

    # Reference
    refer = [0.0;0.0;0.0;100.0;0.0;0.0]
    
    @test all(F.==refer)


end
