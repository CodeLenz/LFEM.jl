using Test
using BMesh, LMesh, TMeshes
using LFEM
using LinearAlgebra

# base
include("test_base.jl")

# truss2D
include("test_truss2D.jl")

# truss3D
include("test_truss3D.jl")

# load.jl
include("test_point_load2D.jl")

# truss3D
include("test_point_load3D.jl")
