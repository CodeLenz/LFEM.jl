module LFEM

using LinearAlgebra, SparseArrays
using StaticArrays
using Plots, Arpack
using BMesh, LMesh


include("load.jl")
include("base.jl")
include("material.jl")
include("elements/truss2D.jl")
include("elements/truss3D.jl")
include("elements/solid2D.jl")
include("elements/solid3D.jl")
include("global.jl")
include("drivers.jl")
include("analysis/linear.jl")
include("analysis/modal.jl")
include("analysis/harmonic.jl")
include("analysis/newmark.jl")
include("show.jl")
include("gmsh.jl")


export Point_load
export Expand_vector, Expand_vector!, To_global, To_local
export Constitutive
export K_truss2D, M_truss2D, B_truss2D, Stress_truss2D
export K_truss3D, M_truss3D, B_truss3D, Stress_truss3D
export K_solid2D, M_solid2D, B_solid2D, Stress_solid2D, Volume_solid2D
export K_solid3D, M_solid3D, B_solid3D, Stress_solid3D, Volume_solid3D
export Global_K, Global_M
export Solve_linear, Solve_modal, Solve_harmonic, Solve_newmark
export Stress, Local_K, Local_M
export Stresses, Volume_truss, Volume_element
export plot
export Gmsh_init, Gmsh_nodal_scalar, Gmsh_element_scalar, Gmsh_nodal_vector
export Gmsh_element_stress

# Precompile some important functions
precompile(dN_solid2D,(Float64,Float64))
precompile(Jacobian_solid2D,(Vector{Float64},Vector{Float64},Matrix{Float64}))
precompile(B_solid2D,(Float64,Float64,Vector{Float64},Vector{Float64}))
precompile(N_solid2D,(Float64,Float64))

precompile(dN_solid3D,(Float64,Float64,Float64))
precompile(Jacobian_solid3D,(Vector{Float64},Vector{Float64},Vector{Float64},Matrix{Float64}))
precompile(B_solid3D,(Float64,Float64,Float64,Vector{Float64},Vector{Float64},Vector{Float64}))
precompile(N_solid3D,(Float64,Float64,Float64))

precompile(Expand_vector!,(Vector{Float64},Vector{Float64},Vector{Int64}))
precompile(Expand_vector,(Vector{Float64},Int64,Vector{Int64}))

precompile(To_local,(SMatrix,Mesh2D,Int64))
precompile(To_local,(SMatrix,Mesh3D,Int64))
precompile(To_local,(MMatrix,Mesh2D,Int64))
precompile(To_local,(MMatrix,Mesh3D,Int64))

precompile(To_global,(SMatrix,Mesh2D,Int64))
precompile(To_global,(SMatrix,Mesh3D,Int64))
precompile(To_global,(MMatrix,Mesh2D,Int64))
precompile(To_global,(MMatrix,Mesh3D,Int64))

precompile(To_local,(SVector,Mesh2D,Int64))
precompile(To_local,(SVector,Mesh3D,Int64))
precompile(To_local,(MVector,Mesh2D,Int64))
precompile(To_local,(MVector,Mesh3D,Int64))

precompile(To_global,(SVector,Mesh2D,Int64))
precompile(To_global,(SVector,Mesh3D,Int64))
precompile(To_global,(MVector,Mesh2D,Int64))
precompile(To_global,(MVector,Mesh3D,Int64))

precompile(Volume_truss,(Mesh2D,Int64))
precompile(Volume_truss,(Mesh3D,Int64))


precompile(Volume_element,(Mesh2D,Int64))
precompile(Volume_element,(Mesh3D,Int64))

end # module
