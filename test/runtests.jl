module GridapGmshTests

using GridapGmsh
using Test

if GridapGmsh.GMSH_FOUND
  @time @testset "gmsh" begin
    include("gmshTests.jl")
  end

  @time @testset "ElementTypes" begin
    include("ElementTypesTests.jl")
  end

  @time @testset "GmshDiscreteModel" begin
   include("GmshDiscreteModelsTests.jl")
  end

  @time @testset "Distributed" begin
   include("DistributedTests.jl")
  end
else
  @warn "GridapGmsh is not loaded or installed properly."
end

end # module
