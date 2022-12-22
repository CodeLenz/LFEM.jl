
#
# Auxiliary functions defined in 
# https://julialinearalgebra.github.io/ArnoldiMethod.jl/stable/usage/02_spectral_transformations.html
#
struct ShiftAndInvert{TA,TB,TT}
    A_lu::TA
    B::TB
    temp::TT
end

function (M::ShiftAndInvert)(y,x)
    mul!(M.temp, M.B, x)
    ldiv!(y, M.A_lu, M.temp)
end

function construct_linear_map(A,B)
    a = ShiftAndInvert(lu(A),B,Vector{eltype(A)}(undef, size(A,1)))
    LinearMap{eltype(A)}(a, size(A,1), ismutating=true)
end

#
# Solve the generalized eigenvalue problem
#
"""
Evaluates the smaller nev eigenvalues for the generalized problem (A - Bλ)X = 0

    Solve_Eigen_(A::AbstractMatrix, B::AbstractMatrix, nev=4, positive=true, order=true, tol=1E-5, restarts=500)

A, B -> Abstract matrices
nev  -> Number of eigenvalues (2*nev are evaluated to account for possible negative values)   
positive -> only positive values are returned
order    -> eigenvalues are sorted in ascending order
tol and restarts are defined in ArnoldiMethod

Returns

eigenvalues (real vector)
eigenvectors (real matrix)

"""
function Solve_Eigen_(A::AbstractMatrix, B::AbstractMatrix, nev=4, positive=true, order=true, tol=1E-5, restarts=500)

    # Target the largest eigenvalues of the inverted problem
    decomp, history  = partialschur(construct_linear_map(A, B), nev=2*nev, tol=tol, restarts=restarts, which=LM())
    λs_inv, X = partialeigen(decomp)

    # Eigenvalues have to be inverted to find the smallest eigenvalues of the non-inverted problem.
    λs = real.(1 ./ λs_inv)

    # There is no guarantee that eigenvalues are positive and are in crescent order. Thus, we 
    # can enforce this situation
    if positive
        
        # Extract just the positive values
        λp = λs[λs.>0.0] 

        # Verify if we have at least nev positive eigenvalues
        length(λp)>=nev || throw("Solve_Eigen_:: there are not enough positive eigenvalues $(length(λp)) of $nev")
    else
        λp = λs
    end

    if order
       p = sortperm(λp) 
    else
       p = collect(1:length(λp))
    end

    # Return both the eigenvalues and the eigenvectors
    return λp[p[1:nev]], real.(X[:,p[1:nev]])
end
