module GridapGmshTests

using GridapGmsh
using Test

@testset "gmsh" begin
  include("gmshTests.jl")
end

@testset "GmshDiscreteModel" begin
 include("GmshDiscreteModelsTests.jl")
end

@testset "LagrangianFESpaces" begin
 include("LagrangianFESpacesTests.jl")
end

end # module
