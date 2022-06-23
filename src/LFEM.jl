module LFEM

using LinearAlgebra, SparseArrays, BMesh, LMesh, Plots

include("load.jl")
include("elements/truss2D.jl")
include("elements/truss3D.jl")
include("global.jl")
include("drivers.jl")
include("analysis/linear.jl")
include("show.jl")
include("gmsh.jl")

export Point_load
export K_truss2D, B_truss2D, Stress_truss2D
export K_truss3D, B_truss3D, Stress_truss3D
export Global_K, Global_M
export Solve_KU
export Stress, Keg_truss, Meg_truss
export plot
export Gmsh_init, Gmsh_nodal_scalar, Gmsh_element_scalar, Gmsh_nodal_vector

end # module
