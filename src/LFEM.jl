module LFEM

using LinearAlgebra, SparseArrays, BMesh, LMesh

include("load.jl")
include("truss2D.jl")
include("truss3D.jl")
include("global.jl")

export Point_load
export K_truss2D
export K_truss3D
export Global_K, Solve_KU

end # module
