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

# solid2D
include("test_solid2D.jl")

# solid3D
include("test_solid3D.jl")

# load.jl
include("test_point_load2D.jl")

# truss3D
include("test_point_load3D.jl")

# Solve_linear
include("test_linear.jl")

# Solve_Eigen_
include("test_arnoldi.jl")

# Solve_modal
include("test_modal.jl")

# Solve_harmonic
include("test_harmonic.jl")

# Solve_newmark
include("test_newmark.jl")

# Options :Stiffness
include("test_options_stiffness.jl")

# Options :Mass
include("test_options_mass.jl")

# Options :IS_TOPO
include("test_options_istopo.jl")

# Options :Damper
include("test_options_damper.jl")
