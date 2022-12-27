@testset "Solid2D" begin

    #    Valid inputs (no error) 
    #
    #  Bar in tension (EPT)
    #
    # >------------ -> F/2
    #  |         2| 
    # X------------ -> F/2
    #       5     o
    #            ---
    nn = 4
    ne = 1
    coord = [0.0 0.0 ;
             5.0 0.0 ;
             5.0 2.0 ;
             0.0 2.0]
    connect = [1 2 3 4]
    Lx = 1.0
    Ly = 1.0
    nx = 1
    ny = 1
    etype = :solid2D

    # Bmesh
    b2 = Bmesh2D(etype,nn,ne,coord,connect,Lx,Ly,nx,ny)

    # Essential boundary conditions
    hebc = [1 1  ; 
           1 2  ; 
           2 2  ;
           4 1 ]

    # Natural boundary conditions
    nbc = [2 1 50.0 ; 
           3 1 50.0 ]

    # Material and geometry
    materials = [Material(Ex=100.0,νxy=0.3,density=10.0)]
    geometries = [Geometry(thickness=0.1)]
           
    # Mesh
    m2 = Mesh2D(b2,materials,geometries,hebc,nbc)

    # Total mass is density*volume
    mass = 10.0*5*2*0.1

    # Mass matrix
    M = Local_M(m2,1)
    
    # sum should be equal 2*mass
    @test isapprox(sum(M), 2*mass)

    # Displacement
    U,_ = Solve_linear(m2)

    # reference is u = FL/EA, where A = height*thickness
    refu = 100*5/(100.0 * 2 * 0.1)
    @test isapprox(U[3],refu)
    @test isapprox(U[5],refu)
   
    # vertical displacement is proportional to Poisson's ratio
    # exx = F/EA
    exx = 100/(100.0*2*0.1)
    eyy = -0.3*exx
    uy = eyy*2.0
    @test isapprox(U[6],uy)
    @test isapprox(U[8],uy)

end
