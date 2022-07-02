push!(LOAD_PATH,"../src/")
using BMesh
using LMesh
using LFEM
using Documenter

makedocs(
         sitename = "LFEM",
         modules  = [LFEM],
         pages=[
                "Home" => "index.md"
               ])
               
deploydocs(;
    versions = nothing,
    repo="github.com/CodeLenz/LFEM.jl",
)
