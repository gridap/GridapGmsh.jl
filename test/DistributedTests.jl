module DistributedTests

using Gridap
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using Test

function main(parts)
  mshfile = joinpath(@__DIR__,"..","demo","demo.msh")
  model = GmshDiscreteModel(parts,mshfile)
  k = 2
  Ω = Interior(model)
  reffe = ReferenceFE(lagrangian,Float64,k)
  V = FESpace(model,reffe)
  u(x) = sum(x)
  uh = interpolate(u,V)
  eh = u - uh
  dΩ = Measure(Ω,2*k)
  @test sum(∫( eh*eh )dΩ) < 1.0e-9
  writevtk(Ω,"Ω",cellfields=["uh"=>uh,"eh"=>eh])
end

prun(main,sequential,6)
prun(main,sequential,1)
prun(main,mpi,1)

end # module
