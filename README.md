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
sigma_h = Harmonic_stresses(mesh,Ud,w,β_c,x,sparam);

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


## Stress smoothing and error estimates (with refinement indicator)
```julia
#
# Cantilever beam with different stress computation strategies and
# stress smoothing, error estimates and refinement estimate
#
#
using BMesh, LMesh, TMeshes
using LFEM


# Data

# Length and heigt
Lx = 1.0
Ly = 0.1

# Thickness
esp = 0.1

# Number of elements in each direction
nx = 20*4
ny = 2*4

# Maximimum admissible error (for stress)
# in %
admissible = 5/100

# For some problems, stress can tend to infinity as the mesh is 
# refined (point loads and right angle corners, for examploe).
# Thus, we can use a skiplist to avoid compute errors
# in these elements.
#
# Avoid the elements in the corner of the clamped left side
# and the element with point load
# 
skiplist = [1;2;nx;
            nx+1;nx+2;
            nx*(ny-2)+1;nx*(ny-2)+2; nx*(ny-1)+1;nx*(ny-1)+2]

# Alternativelly, on can consider all elements
#skiplist = Int64[]

# Elements to be considered
elist = setdiff(collect(1:(nx*ny)),skiplist)

#
#
# BEAM 1 x 0.1 x 0.1, E = 1GPa e F=1kN
#
# Maximum vertical displacement is : FL³ / (3EI) = 0.04 m 
# 
# Maximum normal stress is : 6(FL)/(bh^2) = 6 Mpa
#
# Maximum shear stress is : (4/3)(F/bh) = 133.3 kPa
#
#

# Load Simply supported 2D from TMeshes
# We are using Poissons ratio = 0 to compare with the beam theory.
# 
mesh = Cantilever_beam_bottom2D(nx,ny,:solid2D,Lx=Lx, Ly=Ly, force=1000.0,
                                Ex=1E9,νxy=0.0,thickness=esp)

# Turn the incompatible mode on
mesh.options[:INCOMPATIBLE]=[1.0 1.0]

# Turn the incomplatible mode of
#delete!(mesh.options,:INCOMPATIBLE)

# Solve the linear static equilibrium
U, F, linsolve = Solve_linear(mesh);

# Reference displacement (beam theory)
println("Displacement is  ", U[end], " analytical solution is 0.04 m ")

# Array with stresses - Default is at the center of the element
# xx yy xy para cada elemento (linha)
sigma = Stresses(mesh,U);

# Initilize an output file
name = "output.pos";
Gmsh_init(name,mesh);

# Export displacements
Gmsh_nodal_vector(mesh,U,name,"Displacement [m]");

# Export stresses
Gmsh_element_stress(mesh,sigma,name,"Centroidal: Stress [Pa]");

# It is also possible to evaluate stresses at the superconvergent
# points (since 2D and 3D elasticity elements are nonconforming)
# xx yy xy | xx yy xy | ....
sigma_g = Stresses(mesh,U,center=false);

# Export stresses. In this case, stresses are interpolated to the 
# nodes to be compatible with gmsh. Thus, there is already some
# smoothing due to interpolation.
#
# Note that the function has a different name than before !
Gmsh_element_stresses(mesh,sigma_g,name,"Gauss Points: Stress [Pa]");

#
# Smooth, error and adaptive refinement
#

# One can also apply nodal smoothing to stresses. There are some
# approaches avaliable. 

# For example, the simple nodal average (scaled by the element volume)
nodal_smooth_sigma = Nodal_stress_smooth(mesh,sigma)

# Export the nodal stresses to gmsh
Gmsh_nodal_stress(mesh,nodal_smooth_sigma,name,"Nodal smooth: Stress [Pa]");

# Global smoothing
global_smooth_sigma = Global_stress_smooth(mesh,sigma)

# Export the nodal stresses to gmsh
Gmsh_nodal_stress(mesh,global_smooth_sigma,name,"Global smooth: Stress [Pa]");

# Patch (element centered) smoothing. Patches are bilinear by now
patch_smooth_sigma = Patch_stress_smooth(mesh, sigma)

# Export the nodal stresses to gmsh
Gmsh_nodal_stress(mesh,patch_smooth_sigma,name,"Patch smooth: Stress [Pa]");

# Evaluate the element "error" L2 for each element
errors_nodal = [Element_error_stress(mesh,ele,sigma,nodal_smooth_sigma) for ele=1:mesh.bmesh.ne]
errors_global = [Element_error_stress(mesh,ele,sigma,global_smooth_sigma) for ele=1:mesh.bmesh.ne]
errors_patch = [Element_error_stress(mesh,ele,sigma,patch_smooth_sigma) for ele=1:mesh.bmesh.ne]

# Expor to gmsh
Gmsh_element_scalar(mesh,errors_nodal,name,"e nodal")
Gmsh_element_scalar(mesh,errors_global,name,"e global")
Gmsh_element_scalar(mesh,errors_patch,name,"e patch")

# Zero out the elements in the skiplist
errors_nodal[skiplist] .= 0 
errors_global[skiplist] .= 0 
errors_patch[skiplist] .= 0 

# Global errors (sum)
error_nodal = sum(errors_nodal)
error_global = sum(errors_global)
error_patch = sum(errors_patch)

# "Quality" of each element (% of error)
quality_nodal = 100*(errors_nodal./(error_nodal))
quality_global = 100*(errors_global./(error_global))
quality_patch = 100*(errors_patch./(error_patch))

# Export to gmsh
Gmsh_element_scalar(mesh,quality_nodal,name,"Quality nodal %")
Gmsh_element_scalar(mesh,quality_global,name,"Quality global %")
Gmsh_element_scalar(mesh,quality_patch,name,"Quality patch %")

# Effective number of elements
nee = length(elist)

# Evaluate the "indicator". If the element error is larger than
# admissible*total error, the indicator is larger than 1. Thus
# this element must be refined
indicator_nodal = sqrt(nee)*errors_nodal/(admissible*error_nodal)
indicator_global = sqrt(nee)*errors_global/(admissible*error_global)
indicator_patch = sqrt(nee)*errors_patch/(admissible*error_patch)

# Export to gmsh
Gmsh_element_scalar(mesh,indicator_nodal,name,"Indicator: nodal smooth")
Gmsh_element_scalar(mesh,indicator_global,name,"Indicator: global smooth")
Gmsh_element_scalar(mesh,indicator_patch,name,"Indicator: patch smooth")

#
# h "refinement" (just the estimates, not the refinement)
#

# Element size (the mesh is regular and all elements are assumd to be 
# the same size)
h = min(Lx/nx,Ly/ny)
println("Element size " , h)

# Error estimate - nodal smooth
estimate_h_nodal = h ./ sqrt.(indicator_nodal.+1E-12)

# Error estimate - global smooth
estimate_h_global = h ./ sqrt.(indicator_global.+1E-12)

# Error estimate - element centered superconvergent patch recovery
estimate_h_patch = h ./ sqrt.(indicator_patch.+1E-12)

# Zero ou elements in the skiplist
estimate_h_nodal[skiplist] .= 0
estimate_h_global[skiplist] .= 0
estimate_h_patch[skiplist] .= 0

# Export new element sizes (estimates)
Gmsh_element_scalar(mesh,estimate_h_nodal,name,"h nodal (original is $(h))")
Gmsh_element_scalar(mesh,estimate_h_global,name,"h global (original is $(h))")
Gmsh_element_scalar(mesh,estimate_h_patch,name,"h patch (original is $(h))")

```
