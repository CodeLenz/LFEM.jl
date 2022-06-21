module LFEM

using LinearAlgebra, SparseArrays, BMesh, LMesh, Plots

include("load.jl")
include("truss2D.jl")
include("truss3D.jl")
include("global.jl")
include("drivers.jl")
include("show.jl")

export Point_load
export K_truss2D, B_truss2D, Stress_truss2D
export K_truss3D, B_truss3D, Stress_truss3D
export Global_K, Solve_KU
export Stress
export plot

end # module
