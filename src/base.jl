#
# Copy v to V[free_dofs] inplace
#
function Expand_vector!(V::Vector{T},v::Vector{T},free_dofs::Vector{Int64}) where {T}

    # Basic assertions
    length(v)==length(free_dofs)||throw("Expand_vector!:: v and free_dofs must have the same size")
    length(V)>=length(free_dofs)||throw("Expand_vector!:: length of V must be larger than the length of v")

    # Copy
    V[free_dofs].=v
end


function Expand_vector(v::Vector{T},dim::Int64,free_dofs::Vector{Int64}) where {T}

    V = zeros(T,dim)
    Expand_vector!(V,v,free_dofs)
    return V

end

#
# Find the coordinates x and y of a given element
#
function Nodal_coordinates(m::Mesh2D,ele)

    # Alias
    bm = m.bmesh

    # Nodes
    nodes = Conect(bm,ele)

    # Number of nodes
    nn = length(nodes)

    # Coordinates 
    x = Vector{Float64}(undef,nn)
    y = Vector{Float64}(undef,nn)
    for i=1:nn

        # Local coordinates
        x[i],y[i] = Coord(bm,nodes[i])

    end

   return x, y
end


#
# Find the coordinates x and y of a given element
#
function Nodal_coordinates(m::Mesh3D,ele)

    # Alias
    bm = m.bmesh

    # Nodes
    nodes = Conect(bm,ele)

    # Number of nodes
    nn = length(nodes)

    # Coordinates 
    x = Vector{Float64}(undef,nn)
    y = Vector{Float64}(undef,nn)
    z = Vector{Float64}(undef,nn)
    
    for i=1:nn

        # Local coordinates
        x[i],y[i],z[i] = Coord(bm,nodes[i])

    end

   return x, y, z
end