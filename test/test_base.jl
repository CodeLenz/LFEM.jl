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
  
  
