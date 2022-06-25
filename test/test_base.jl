@testset "Base" begin
  
  # Test Expand_vector!
  
  U = zeros(10)
  u = ones(5)
  pos = [1;3;4;8;10]
  Expand_vector!(U,u,pos)
  @test all(U.==[1.0;0.0;1.0;1.0;0.0;0.0;0.0;1.0;0.0;1.0])
  
  # Should throw - U2 is smaller than u2
  U2 = zeros(2)
  u2 = zeros(5)
  @test_throws String Expand_vector!(U2,u2,pos)
  
  # Should throw - length of u and pos2 are different
  pos2 = [1;2]
  @test_throws String Expand_vector!(U,u,pos2)
  
  
  # Test Expand_vector
  
  u = ones(5)
  pos = [1;3;4;8;10]
  @test all(Expand_vector(u,10,pos).==[1.0;0.0;1.0;1.0;0.0;0.0;0.0;1.0;0.0;1.0])
  
  # Should throw - length of u and pos2 are different
  pos2 = [1;2]
  @test_throws String Expand_vector(u,10,pos2)
  
  #
  # Nodal coordinates
  #
  # Nodal_coordinates(m::Mesh2D,ele)

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
    materials = [Material(Ex=100.0,density=1.0)]
    geometries = [Geometry(A=1.0)]
           
    # Mesh
    m2 = Mesh2D(b2,materials,geometries,ebc,nbc)

    # Nodal coordinates of element 1 
    x,y=Nodal_coordinates(m2,1)
    xref = [0.0;0.1]
    yref = [0.0;0.0]
    @assert all(x.==xref)
    @assert all(y.==yref)


    ####################  2D BMesh ######################
    # Square
    #
    #    Valid inputs (no error) 
    #
    nn = 4
    ne = 1
    coord = [0.0 0.0 ;
             1.0 0.0 ;
             1.1 1.1 ;
             0.0 1.0]
    connect = [1 2 3 4]
    Lx = 1.0
    Ly = 1.0
    nx = 1
    ny = 1
    etype = :solid2D

    # truss2D
    b2 = Bmesh2D(etype,nn,ne,coord,connect,Lx,Ly,nx,ny)

    # Essential boundary conditions
    ebc = [1 1 0.0 ; 
           1 2 0.0 ; 
           2 1 0.0 ]

    # Natural boundary conditions
    nbc = [3 1 100.0]

    # Material and geometry
    materials = [Material(Ex=100.0,density=1.0)]
    geometries = [Geometry(thickness=1.0)]
           
    # Mesh
    m2 = Mesh2D(b2,materials,geometries,ebc,nbc)

    # Nodal coordinates of element 1 
    x,y=Nodal_coordinates(m2,1)
    xref = [0.0;1.0;1.1;0.0]
    yref = [0.0;0.0;1.1;1.0]
    @assert all(x.==xref)
    @assert all(y.==yref)

  
end
