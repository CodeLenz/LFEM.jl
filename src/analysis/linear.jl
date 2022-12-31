"""
Solve the linear system K(x)U(x)=F

    Solve_linear(mesh::Mesh, x::Vector{Float64}, kparam::Function; loadcase=1)

where 

    x is a ne x 1 vector of design varibles 
    kparam(xe): R->R is the material parametrization (SIMP like)
    loadcase is the loadcase

Returns:

    U = displacement vector (dim*nn x 1)
    F = Force vector (dim*nn x 1)
    linsolve = LinearSolve object with factored linear problem
"""
function Solve_linear(mesh::Mesh, x::Vector{Float64}, kparam::Function; loadcase::Int64=1)
  
    # Basic assertions
    length(x)==Get_ne(mesh) || throw("Solve_linear:: length of x must be ne")

    0<=loadcase<=mesh.nload || throw("Solve_linear:: invalid loadcase")

    # Assembly
    K = Global_K(mesh,x,kparam)
    F = Point_load(mesh,loadcase)

    # Free dofs
    free_dofs = mesh.free_dofs[loadcase]
    
    # View
    K =  K[free_dofs,free_dofs]

    # Create LinearSolve problem
    prob = LinearProblem(Symmetric(K),F[free_dofs])
    linsolve = init(prob)

    # Solve
    Ul = solve(linsolve)

    # Expand homogeneous ebc
    Us  = zeros(length(F))
    Expand_vector!(Us,Ul.u,free_dofs)
    
    return Us, F, linsolve
    
 end

 """
 Solve the linear system KU=F
 
     Solve_linear(mesh::Mesh; loadcase=1)
 
 Returns:
 
     U = displacement vector (dim*nn x 1)
     F = Force vector (dim*nn x 1)
     linsolve = LinearSolve object with factored linear problem
 """ 
function Solve_linear(mesh::Mesh;loadcase::Int64=1)

    # x->1.0 mapping
    dummy_f(x)=1.0

    # x is not used
    x = Vector{Float64}(undef,Get_ne(mesh))

    # Call Solve_linear
    Solve_linear(mesh,x,dummy_f,loadcase=loadcase)
    
end
  