@testset "Arnoldi" begin

   # Define matrices A and B
   A = [20.0 0.0 ; 0.0 40.0]
   B = [30.0 0.0 ; 0.0 40.0]

   # Solve with our method
   flag, av, AV = Solve_Eigen_(A,B,2)

   # Test flag
   @test flag==1

   # Test eigenvalues
   @test isapprox(av[1],2/3)
   @test isapprox(av[2],1.0)

   # Test eigenvectors
   @test isapprox(AV[:,1],[1.0 ; 0.0])
   @test isapprox(AV[:,2],[0.0 ; 1.0])
   

    # Larger problem 
    A = diagm(collect(10:10:100)) .+ rand(10,10)
    B = diagm(collect(1.0:-0.1:0.1)) 

    # Solve with our method
    nev = 6
    flag, av, AV = Solve_Eigen_(A,B,nev,positive=false)

    # test flag
    @test flag==1

    @test isapprox(norm(A*AV .- B*AV*Diagonal(av)),0.0,atol=1E-6)

    # Solve with eigen
    aref, vref = eigen(A,B)

    # Test for eigenvalues
    for i=1:nev
        @test isapprox(aref[i],av[i])
    end


    # Larger problem - just positive eigenvalues
    A = diagm(collect(10:10:100)) .+ rand(10,10)
    B = diagm(collect(1.0:-0.1:0.1)) 

    # Solve with our method
    nev = 2
    flag, av, AV = Solve_Eigen_(A,B,nev,positive=true)

    # test flag
    @test flag==1

    @test isapprox(norm(A*AV .- B*AV*Diagonal(av)),0.0,atol=1E-6)

    # Solve with eigen
    aref, vref = eigen(A,B)

    # Test for eigenvalues
    for i=1:nev
        @test isapprox(aref[i],av[i])
    end


    # interesting test case - Arnoldi does not solve it 
    # since eigenvalues are too close.
    A = 2.0.*I(10)
    B = 5.0.*I(10)
    nev = 4
    flag, av, AV = Solve_Eigen_(A,B,nev,ortho=false)

    @test flag == -1
    @test av == zeros(1)
    @test AV == zeros(1,1)

    # Thus, we revert to 
    af, Vf = Failed_Arnoldi(A, B, nev, positive=true)

    @test isapprox(af, 0.4*ones(length(af)))
    #for i=1:nev
    #   rf = zeros(10)
    #   rf[i] = 1.0
    #   @test isapprox(Vf[:,i],rf)
    #end

end