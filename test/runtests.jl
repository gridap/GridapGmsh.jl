module GridapGmshTests

using GridapGmsh
using Test

if GridapGmsh.GMSH_FOUND
  @testset "gmsh" begin
    include("gmshTests.jl")
  end

  @testset "GmshDiscreteModel" begin
   include("GmshDiscreteModelsTests.jl")
  end

  @testset "Distributed" begin
   include("DistributedTests.jl")
  end
else
  @warn "GridapGmsh is not loaded or installed properly."
end

end # module
