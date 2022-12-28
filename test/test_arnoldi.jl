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
   

    # Larger problem
    A = diagm(collect(10:10:100)) .+ rand(10,10)
    B = diagm(collect(1.0:-0.1:0.1)) 

    av, AV = Solve_Eigen_(A,B,6,false)

    @test isapprox(norm(A*AV .- B*AV*Diagonal(av)),0.0,atol=1E-6)

    # Test for ortogonality
    #orto = Array(qr(AV).Q)
    #for i=1:20
    #    for j=1:20
    #        if i!=j
    #           @test isapprox(dot(orto[:,i],orto[:,j]),0.0,atol=1E-6)
    #        end
    #    end
    #end


end