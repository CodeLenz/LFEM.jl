module LFEM

       using LinearAlgebra, SparseArrays
       using ArnoldiMethod, LinearMaps
       using LinearSolve
       using StaticArrays
       using BMesh, LMesh

       # From version 0.4 of Arnoldy methods, one has to explicitly 
       # import the symbols
       using ArnoldiMethod: LM

       # If possible, set optimization to 3
       if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
              @eval Base.Experimental.@optlevel 3
       end

       include("load.jl")
       include("base.jl")
       include("material.jl")
       include("solve_arnoldi.jl")
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
       include("gmsh.jl")


       export Solve_Eigen_, Failed_Arnoldi, Residue_Eigenpair
       export Point_load
       export Expand_vector, Expand_vector!, To_global, To_local
       export Constitutive, Equivalent_stress
       export Gauss_2D, Gauss_3D
       export K_truss2D, Ks_truss2D, M_truss2D, B_truss2D, Stress_truss2D
       export K_truss3D, Ks_truss3D, M_truss3D, B_truss3D, Stress_truss3D
       export K_solid2D, M_solid2D, B_solid2D, Stress_solid2D, Volume_solid2D
       export K_solid3D, M_solid3D, B_solid3D, Stress_solid3D, Volume_solid3D
       export Global_K, Global_Ks, Global_M, Global_C
       export Solve_linear, Solve_modal, Solve_harmonic, Solve_newmark
       export Organize_Eigen
       export Stress, Local_K, Local_M, Local_Ks
       export Stresses, Volume_truss, Volume_element
       export Gmsh_init, Gmsh_nodal_scalar, Gmsh_element_scalar, Gmsh_nodal_vector
       export Gmsh_element_stress, Map_stress2nodes_Quad, Gmsh_element_stresses
       export Harmonic_stresses, Voigt_equivalent
       export B_element

end # module
