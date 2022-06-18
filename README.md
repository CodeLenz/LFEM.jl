# LFEM
Basic routines for FEM

```julia

using BMesh, LMesh, TMeshes
using LFEM

# Load Simply supported 2D from TMeshes
m = Simply_supported2D(6,6)

# Solve the linear system 
U, F, Chol = Solve_KU(m)

# Add a vector of design variables and use a SIMP
# exponent p=3.0
x = ones(length(m.bmesh.ne))
U, F, Chol = Solve_KU(m; x=x, p=3.0)
