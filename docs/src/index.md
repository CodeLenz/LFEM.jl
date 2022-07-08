# LFEM.jl

* Basic routines for FEM *

## Solvers

```@docs
Solve_linear
```

```@docs
Solve_modal
```

```@docs
Solve_harmonic
```

```@docs
Solve_newmark
```
## Global

```@docs
Global_K
```

```@docs
Global_M
```

```@docs
Global_C
```

```@docs
Stresses
```

```@docs
Harmonic_stresses
```

## Truss2D

```@docs
K_truss2D
```

```@docs
M_truss2D
```

```@docs
B_truss2D
```

```@docs
B_truss2D
```
```@docs
Stress_truss2D
```

## Truss3D

```@docs
K_truss3D
```

```@docs
M_truss3D
```

```@docs
B_truss3D
```

```@docs
B_truss3D
```

```@docs
Stress_truss3D
```

## Solid2D

```@docs
dN_solid2D
```

```@docs
Jacobian_solid2D
```

```@docs
B_solid2D
```

```@docs
K_solid2D
```

```@docs
N_solid2D
```

```@docs
M_solid2D
```

```@docs
Stress_solid2D
```

```@docs
Volume_solid2D
```

## Solid3D

```@docs
dN_solid3D
```

```@docs
Jacobian_solid3D
```

```@docs
B_solid3D
```

```@docs
K_solid3D
```

```@docs
N_solid3D
```

```@docs
M_solid3D
```

```@docs
Stress_solid3D
```

```@docs
Volume_solid3D
```

## Gmsh

```@docs
Gmsh_init
```

```@docs
Gmsh_nodal_scalar
```

```@docs
Gmsh_element_scalar
```

```@docs
Gmsh_nodal_vector
```

```@docs
Gmsh_element_stress
```