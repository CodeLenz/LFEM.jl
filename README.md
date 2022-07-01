# LFEM
Basic routines for FEM

Linear elastic analysis in 2D - truss
```julia

using BMesh, LMesh, TMeshes, Plots
using LFEM

# Load Simply supported 2D from TMeshes
m = Simply_supported2D(6,6)

# Solve the linear static equilibrium
U, F, Chol = Solve_linear(m)

# Show displacements
plot(m;U=U)

```

Linear elastic analysis in 2D - Plane Stress
```julia

using BMesh, LMesh, TMeshes, Plots
using LFEM

# Load Simply supported 2D from TMeshes
m = Simply_supported2D(6,6,:solid2D)

# Solve the linear static equilibrium
U, F, Chol = Solve_linear(m)

# It is also possible to export the results to
# gmsh (https://gmsh.info/)

# Initilize an output file
name = "output.pos"
Gmsh_init(name,m)

# Export displacements
Gmsh_nodal_vector(m,U,name,"Displacement [m]")

```

Linear elastic analysis in 3D (truss)
```julia

using BMesh, LMesh, TMeshes, Plots
using LFEM

# Load Simply supported 3D from TMeshes
m = Simply_supported3D(6,6,6)

# Solve the linear static equilibrium
U, F, Chol = Solve_linear(m)

# Show displacements
plot(m;U=U)

```

Linear elastic analysis in 3D
```julia

using BMesh, LMesh, TMeshes, Plots
using LFEM

# Load Simply supported 2D from TMeshes
m = Simply_supported3D(6,6,6,:solid3D)

# Solve the linear static equilibrium
U, F, Chol = Solve_linear(m)

# It is also possible to export the results to
# gmsh (https://gmsh.info/)

# Initilize an output file
name = "output.pos"
Gmsh_init(name,m)

# Export displacements
Gmsh_nodal_vector(m,U,name,"Displacement [m]")

```



Stresses
```julia

using BMesh, LMesh, TMeshes, Plots
using LFEM

# Load Simply supported 2D from TMeshes
m = Simply_supported2D(6,6)

# Solve the linear static equilibrium
U, F, Chol = Solve_linear(m)

# Array with stresses
sigma = Stresses(m,U)

# Show displacements and stresses
plot(m;U=U,N=sigma)

# It is also possible to export the results to
# gmsh (https://gmsh.info/)

# Initilize an output file
name = "output.pos"
Gmsh_init(name,m)

# Export stresses
Gmsh_element_scalar(m,sigma,name,"Stress [Pa]")

# Export displacements
Gmsh_nodal_vector(m,U,name,"Displacement [m]")

```

Modal analysis
```julia

using BMesh, LMesh, TMeshes, Plots
using LFEM

# Load Simply supported 3D from TMeshes
m = Simply_supported3D(6,6,6)

# Solve the modal problem with default parameters
λ, ϕ = Solve_modal(m)

# It is also possible to export the results to
# gmsh (https://gmsh.info/)

# Initilize an output file
name = "output.pos"
Gmsh_init(name,m)

# Export first mode
w1 = sqrt(real(λ[1]))
Gmsh_nodal_vector(m,vec(ϕ[:,1]),name,"First mode, w=$w1 [rad/s]")

```

Harmonic analysis
```julia

using BMesh, LMesh, TMeshes, Plots
using LFEM

# Load Simply supported 2D from TMeshes
m = Simply_supported2D(6,6)

# Angular frequency [rad/s]
w = 20.0

# Solve the harmonic problem
Ud,LU = Solve_harmonic(m,w)

# Initilize an output file
name = "output.pos"
Gmsh_init(name,m)

# Export to gmsh
Gmsh_nodal_vector(m,real.(Ud),name,"Harmonic displacement - real part")
Gmsh_nodal_vector(m,imag.(Ud),name,"Harmonic displacement - imag part")
Gmsh_nodal_vector(m,abs.(Ud),name,"Harmonic displacement - abs")

```

Transient  analysis - Newmark method
```julia

using BMesh, LMesh, TMeshes, Plots
using LFEM

# Load Simply supported 2D from TMeshes
m = Simply_supported2D(6,6)

# We need to define a function that modifies a force
# vector according to time t. Lets use the same point
# load as the static example, but with a cos(2*t) 
function f!(t,F,m)
         P  = Point_load(m)
         F .= cos(2*t)*P
end
 
# And a list of nodes/dofs to monitor. Lets monitor 
# the same DOF as the load
node = Int(m.nbc[1,1])
dof  = Int(m.nbc[1,2])
monitor = [node dof]

# timespan [s]
tspan = (0.0,5.0)

# Interval
dt = 1E-2

# Solve the transisent problem
U,V,A,T = Solve_newmark(m,f!,monitor,tspan,dt)

# Plot displacement 
plot(T,U,xlabel="time [s]",ylabel="Displacement [m]", label="")

```

