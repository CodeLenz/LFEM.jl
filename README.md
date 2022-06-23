# LFEM
Basic routines for FEM

```julia

using BMesh, LMesh, TMeshes, Plots
using LFEM

# Load Simply supported 2D from TMeshes
m = Simply_supported2D(6,6)

# Solve the linear static equilibrium
U, F, Chol = Solve_linear(m)

# Solve the modal problem
λ, ϕ = Solve_modal(m)

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

# Export modes
Gmsh_nodal_vector(m,vec(ϕ[:,1]),name,"First mode, λ=$(λ[1])")

# Add a vector of design variables and use a SIMP
# exponent p=3.0
x = ones(m.bmesh.ne)
U, F, Chol = Solve_linear(m; x=x, p=3.0)

# qp parametrization
Stress(m,6,U;xe=0.4,p=3.0,q=2.5)
