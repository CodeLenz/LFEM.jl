# LFEM
Basic routines for FEM

Docs: https://codelenz.github.io/LFEM.jl/

This package depends on

https://github.com/CodeLenz/BMesh.jl  
https://github.com/CodeLenz/LMesh.jl  
https://github.com/CodeLenz/TMeshes.jl  

It is recomended to install using the following sequence
```julia
using Pkg
Pkg.add(url="https://github.com/CodeLenz/BMesh.jl.git#main#main")
Pkg.add(url="https://github.com/CodeLenz/LMesh.jl.git#main#main")
Pkg.add(url="https://github.com/CodeLenz/TMeshes.git#main#main")
Pkg.add(url="https://github.com/CodeLenz/LFEM.git#main#main")
```

There are currently four types of elements: 2D and 3D bar (truss) elements and 
2D (four node bilinear isoparametric elements) and 3D (eigth node trilinear isoparametric elements).
The elasticity elements are compatible, without additional incompatible (or bubble) modes. To turn on
the elasticity elements into their incompatible couterparts, use
```julia
mesh.options[:INCOMPATIBLE]=[1.0 1.0]
```


# Examples 

## Linear elastic analysis in 2D - truss 
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported2D(6,6);

# Solve the linear static equilibrium
# using the traditional 4 node bilinear isop elements
U, F, linsolve = Solve_linear(mesh);

# Turn the incompatible mode on
mesh.options[:INCOMPATIBLE]=[1.0 1.0]

# Solve the linear static equilibrium
# using incompatible elements
U, F, linsolve = Solve_linear(mesh);

```

## Linear elastic analysis in 2D - Plane Stress 
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported2D(6,6,:solid2D);

# Turn the incompatible mode on
mesh.options[:INCOMPATIBLE]=[1.0 1.0]

# Solve the linear static equilibrium
U, F, linsolve = Solve_linear(mesh);

# It is also possible to export the results to
# gmsh (https://gmsh.info/)

# Initilize an output file
name = "output.pos";
Gmsh_init(name,mesh);

# Export displacements
Gmsh_nodal_vector(mesh,U,name,"Displacement [m]");

```

## Linear elastic analysis in 3D (truss) 
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 3D from TMeshes
mesh = Simply_supported3D(6,6,6);

# Solve the linear static equilibrium
U, F, linsolve = Solve_linear(mesh);

# Show displacements
plot(mesh;U=U);

```

## Linear elastic analysis in 3D 
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported3D(6,6,6,:solid3D);

# Turn the incompatible mode on
mesh.options[:INCOMPATIBLE]=[1.0 1.0]

# Solve the linear static equilibrium
U, F, linsolve = Solve_linear(mesh);

# Initilize an output file
name = "output.pos";
Gmsh_init(name,mesh);

# Export displacements
Gmsh_nodal_vector(mesh,U,name,"Displacement [m]");

```

## Stresses 

Evaluation and gmsh export are automatic for each element type
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported2D(6,6);

# Turn the incompatible mode on
mesh.options[:INCOMPATIBLE]=[1.0 1.0]

# Solve the linear static equilibrium
U, F, linsolve = Solve_linear(mesh);

# Array with stresses
sigma = Stresses(mesh,U);

# Initilize an output file
name = "output.pos";
Gmsh_init(name,mesh);

# Export stresses
Gmsh_element_stress(mesh,sigma,name,"Stress [Pa]");

# Export displacements
Gmsh_nodal_vector(mesh,U,name,"Displacement [m]");

```

## Stresses (solid) 

Evaluation and gmsh export are automatic for each element type
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported3D(6,6,6,:solid3D);

# Turn the incompatible mode on
mesh.options[:INCOMPATIBLE]=[1.0 1.0]

# Solve the linear static equilibrium
U, F, linsolve = Solve_linear(mesh);

# Array with stresses - Default is at the center of the element
sigma = Stresses(mesh,U);

# Initilize an output file
name = "output.pos";
Gmsh_init(name,mesh);

# Export stresses
Gmsh_element_stress(mesh,sigma,name,"Stress [Pa]");

# Export displacements
Gmsh_nodal_vector(mesh,U,name,"Displacement [m]");

```


## Stresses (solid) 

Evaluation and gmsh export are automatic for each element type
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported2D(6,6,:solid2D);

# Turn the incompatible mode on
mesh.options[:INCOMPATIBLE]=[1.0 1.0]

# Solve the linear static equilibrium
U, F, linsolve = Solve_linear(mesh);

# Array with stresses - Default is at the center of the element
sigma = Stresses(mesh,U);

# Initilize an output file
name = "output.pos";
Gmsh_init(name,mesh);

# Export stresses
Gmsh_element_stress(mesh,sigma,name,"Centroidal: Stress [Pa]");

# Export displacements
Gmsh_nodal_vector(mesh,U,name,"Displacement [m]");

# It is also possible to evaluate stresses at the superconvergent
# points (since 2D and 3D elasticity elements are nonconforming)
sigma_g = Stresses(mesh,U,center=false);

# Export stresses. In this case, stresses are interpolated to the 
# nodes to be complatible with gmsh.
# Note that the function has a different name!
Gmsh_element_stresses(mesh,sigma_g,name,"Gauss Points: Stress [Pa]");

# One can also apply nodal smoothing to stresses. There are some
# approaches avaliable. For example, the simple nodal average (scaled
# by the element volume)
nodal_smooth_sigma = Nodal_stress_smooth(mesh,sigma)

# Export the nodal stresses to gmsh
Gmsh_nodal_stress(mesh,nodal_smooth_sigma,name,"Nodal smooth: Stress [Pa]");

# Global smoothing
global_smooth_sigma = Global_stress_smooth(mesh,sigma)

# Export the nodal stresses to gmsh
Gmsh_nodal_stress(mesh,global_smooth_sigma,name,"Global smooth: Stress [Pa]");

# Patch (element centered) smoothing
patch_smooth_sigma = Patch_stress_smooth(mesh, sigma)

# Export the nodal stresses to gmsh
Gmsh_nodal_stress(mesh,patch_smooth_sigma,name,"Patch smooth: Stress [Pa]");


```

## Modal analysis 
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 3D from TMeshes
mesh = Simply_supported3D(6,6,6);

# Solve the modal problem with default parameters
λ, ϕ = Solve_modal(mesh);

# Initilize an output file
name = "output.pos";
Gmsh_init(name,mesh);

# Export first mode
w1 = sqrt(real(λ[1]));
Gmsh_nodal_vector(mesh,vec(ϕ[:,1]),name,"First mode, w=$w1 [rad/s]");

```

## Harmonic analysis 
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported2D(6,6);

# Angular frequency [rad/s]
w = 20.0

# Structural damping
α_c = 0.0
β_c = 1E-6

# Solve the harmonic problem
Ud, linsolve = Solve_harmonic(mesh,w,α_c,β_c);

# Initilize an output file
name = "output.pos";
Gmsh_init(name,mesh);

# Export to gmsh
Gmsh_nodal_vector(mesh,real.(Ud),name,"Harmonic displacement - real part");
Gmsh_nodal_vector(mesh,imag.(Ud),name,"Harmonic displacement - imag part");
Gmsh_nodal_vector(mesh,abs.(Ud),name,"Harmonic displacement - abs");

```

## Harmonic analysis with stress
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported2D(6,6);

# Angular frequency [rad/s]
w = 20.0

# Structural damping
α_c = 0.0
β_c = 1E-6

# Solve the harmonic problem
Ud, linsolve = Solve_harmonic(mesh,w,α_c,β_c);

# Harmonic stresses
sigma_h = Harmonic_stresses(mesh,Ud,w,β_c);

# Initilize an output file
name = "output.pos";
Gmsh_init(name,mesh);

# Export to gmsh
Gmsh_nodal_vector(mesh,real.(Ud),name,"Harmonic displacement [m] - real part");
Gmsh_nodal_vector(mesh,imag.(Ud),name,"Harmonic displacement [m] - imag part");
Gmsh_nodal_vector(mesh,abs.(Ud),name,"Harmonic displacement [m] - abs");

Gmsh_element_stress(mesh,real.(sigma_h),name,"Harmonic stress [Pa] - real part");
Gmsh_element_stress(mesh,imag.(sigma_h),name,"Harmonic stress [Pa] - imag part");
Gmsh_element_stress(mesh,abs.(sigma_h),name,"Harmonic stress [Pa] - abs");

```


## Transient  analysis - Newmark method 
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported2D(6,6);

# We need to define a function that modifies a force
# vector according to time t. Lets use the same point
# load as the static example, but with a cos(2*t) 
function f!(t,F,mesh,load=1)
         P  = Point_load(mesh)
         F .= cos(2*t)*P
end
 
# And a list of nodes/dofs to monitor. Lets monitor 
# the same DOF as the load
node = Int(mesh.nbc[1,1])
dof  = Int(mesh.nbc[1,2])
monitor = [node dof]

# timespan [s]
tspan = (0.0,5.0)

# Interval
dt = 1E-2

# Solve the transient problem
U,V,A,T,dofs = Solve_newmark(mesh,f!,monitor,tspan,dt);

# Plot displacement 
plot(T,U,xlabel="time [s]",ylabel="Displacement [m]", label="");

```

# Material Parametrizations 

It is possible to pass additional arguments to each one of the previous 
examples.

A vector of design variables
```julia
x = ones(Get_ne(mesh))
```
and three material parametrizations (Like SIMP). The first is for K(x)
```julia
kparam(xe::Float64,p=3.0)= xe^p
```
the second is for M(x)
```julia
mparam(xe::Float64,p=2.0,cut=0.1)= ifelse(xe>=cut,xe,xe^p)
```
and the third for stresses
```julia
sparam(xe::Float64,p=3.0,q=2.5)= xe^(p-q)
```
It is important that each one of those parametrizations should be 
function of xe (scalar) only.

## Linear Elastic Analysis with material parametrization
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported2D(6,6)

# Define x and kparam
x = ones(Get_ne(mesh))
kparam(xe::Float64,p=3.0)= xe^p

# Solve the linear static equilibrium passing x and kparam
U, F, linsolve = Solve_linear(mesh,x,kparam);

# Initilize an output file
name = "output.pos"
Gmsh_init(name,mesh)

# Export displacements
Gmsh_nodal_vector(mesh,U,name,"Displacement [m]")

# Export x (just for example)
Gmsh_element_scalar(mesh,x,name,"Design variables")

```

## Stresses with material parametrization
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported3D(6,6,6,:solid3D)

# Define x, kparam and sparam
x = ones(Get_ne(mesh))
kparam(xe::Float64,p=3.0)= xe^p
sparam(xe::Float64,p=3.0,q=2.5)= xe^(p-q)


# Solve the linear static equilibrium
U, F, linsolve = Solve_linear(mesh,x,kparam);

# Array with stresses
sigma = Stresses(mesh,U,x,sparam);

# Initilize an output file
name = "output.pos"
Gmsh_init(name,mesh)

# Export stresses
Gmsh_element_stress(mesh,sigma,name,"Stress [Pa]")

# Export displacements
Gmsh_nodal_vector(mesh,U,name,"Displacement [m]")

```

## Modal analysis with material parametrization
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 3D from TMeshes
mesh = Simply_supported3D(6,6,6);

# Define x, kparam and mparam
x = ones(Get_ne(mesh))
kparam(xe::Float64,p=3.0)= xe^p
mparam(xe::Float64,p=2.0,cut=0.1)= ifelse(xe>=cut,xe,xe^p)

# Solve the modal problem with default parameters
λ, ϕ = Solve_modal(mesh,x,kparam,mparam);

# Initilize an output file
name = "output.pos"
Gmsh_init(name,mesh)

# Export first mode
w1 = sqrt(real(λ[1]))
Gmsh_nodal_vector(mesh,vec(ϕ[:,1]),name,"First mode, w=$w1 [rad/s]")

```

## Harmonic analysis with material parametrization
```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
mesh = Simply_supported2D(6,6)

# Angular frequency [rad/s]
w = 20.0

# Structural damping
α_c = 0.0
β_c = 1E-6

# Define x, kparam and mparam
x = ones(Get_ne(mesh))
kparam(xe::Float64,p=3.0)= xe^p
mparam(xe::Float64,p=2.0,cut=0.1)= ifelse(xe>=cut,xe,xe^p)

# Solve the harmonic problem
Ud, linsolve = Solve_harmonic(mesh,w,α_c,β_c,x,kparam,mparam);

# Array with harmonic stresses
sigma_h = Harmonic_stresses(mesh,U,w,β_c,x,sparam);

# Initilize an output file
name = "output.pos"
Gmsh_init(name,mesh)

# Export to gmsh
Gmsh_nodal_vector(mesh,real.(Ud),name,"Harmonic displacement - real part")
Gmsh_nodal_vector(mesh,imag.(Ud),name,"Harmonic displacement - imag part")
Gmsh_nodal_vector(mesh,abs.(Ud),name,"Harmonic displacement - abs")

Gmsh_element_stress(mesh,real.(sigma_h),name,"Harmonic stress [Pa] - real part")
Gmsh_element_stress(mesh,imag.(sigma_h),name,"Harmonic stress [Pa] - imag part")
Gmsh_element_stress(mesh,abs.(sigma_h),name,"Harmonic stress [Pa] - abs")


```
