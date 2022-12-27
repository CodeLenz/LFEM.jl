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
  
  
  # Change of reference

  ####################  2D BMesh ######################
  #    Valid inputs (no error) 
  #
  #
  #   Element 1 -> horizontal
  #   Element 2 -> vertical
  #
    nn = 3
    ne = 2
    coord = [0.0 0.0 ;
             1.0 0.0 ;
             0.0 1.0]
    connect = [1 2 ;
               1 3 ]
    Lx = 1.0
    Ly = 1.0
    nx = 1
    ny = 1
    etype = :truss2D

    # truss2D
    b2 = Bmesh2D(etype,nn,ne,coord,connect,Lx,Ly,nx,ny)

    # Essential boundary conditions
    hebc = [1 1  ; 
           1 2  ; 
           2 2  ;
           3 1 ]

    # Natural boundary conditions
    nbc = [2 1 100.0;
           3 2 100.0]

    # Material and geometry
    materials = [Material(Ex=100.0)]
    geometries = [Geometry(A=1.0)]
           
    # Mesh
    m2 = Mesh2D(b2,materials,geometries,hebc,nbc)

    # Local Matrices
    K1 = Local_K(m2,1)
    K2 = Local_K(m2,2)

    # They should be equal
    @test all(isapprox(K1,K2))

    # Change K1 to global (must be the same)
    K1g = To_global(K1,m2,1)
    @test all(isapprox(K1,K1g)) 

    # Change K2 to global
    K2g = To_global(K2,m2,2)
    
    cte = 100.0*1.0/1.0
    refer = cte*[0.0  0.0  0.0  0.0 ;
                 0.0  1.0  0.0 -1.0 ;
                 0.0  0.0  0.0  0.0 ;
                 0.0 -1.0  0.0  1.0]

    @test  all(isapprox(K2g,refer))             

    # Bring K2g back to local
    K2l = To_local(K2g,m2,2)
    @test all(isapprox(K2l,K2)) 



end
