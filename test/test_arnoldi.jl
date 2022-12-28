@testset "Arnoldi" begin

   # Define matrices A and B
   A = [20.0 0.0 ; 0.0 40.0]
   B = [30.0 0.0 ; 0.0 40.0]

   # Solve 
   av, AV = Solve_Eigen_(A,B,2)

   # Test eigenvalues
   @test isapprox(av[1],2/3)
   @test isapprox(av[2],1.0)

   # Test eigenvectors
   @test isapprox(AV[:,1],[1.0 ; 0.0])
   @test isapprox(AV[:,2],[0.0 ; -1.0])
   

    # Randon, larger problem
    A = rand(100,100); A = A .+ A' #+ 100*I(100)
    B = rand(100,100); B = B .+ B' #+ 200*I(100)
    av, AV = Solve_Eigen_(A,B,4)

    @test isapprox(norm(A*AV - B*AV*Diagonal(av)),0.0)


end