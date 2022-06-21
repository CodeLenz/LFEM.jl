# LFEM
Basic routines for FEM

```julia

using BMesh, LMesh, TMeshes, Plots
using LFEM

# Load Simply supported 2D from TMeshes
m = Simply_supported2D(6,6)

# Solve the linear system 
U, F, Chol = Solve_KU(m)

# Stress for element 6
Stress(m,6,U)

# Vector with stresses
ne = m.bmesh.ne
sigma = zeros(ne)
for i=1:ne
   sigma[i] = Stress(m,i,U)[1] 
end

# Show displacements and stresses
plot(m;U=U,N=sigma)

# It is also possible to export the results to
# gmsh

# Initilize an output file
name = "output.pos"
Gmsh_init(name,m)

# Export stresses
Gmsh_element_scalar(m,sigma,name,"Stress")

# Export displacements
Gmsh_nodal_vector(m,U,name,"Displacement")

# Add a vector of design variables and use a SIMP
# exponent p=3.0
x = ones(m.bmesh.ne)
U, F, Chol = Solve_KU(m; x=x, p=3.0)

# qp parametrization
Stress(m,6,U;xe=0.4,p=3.0,q=2.5)
