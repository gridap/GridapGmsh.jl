module GmshDiscreteModelsTests

using Gridap
using GridapGmsh
using Gridap.Geometry
using Gridap.Helpers
using Test

function check_interpolation(model)
  order = 2
  reffe = ReferenceFE(:Lagrangian,Float64,order)
  V = FESpace(model,reffe)
  if num_dims(model) == 2
    u = x -> x[1]^2 + 2*x[2]^2
  elseif num_dims(model) == 3
    u = x -> x[1]^2 + 2*x[2]^2 + 3*x[3]^2
  else
    @unreachable
  end
  uh = interpolate(u,V)
  e = u - uh
  Ω = Triangulation(model)
  dΩ = Measure(Ω,2*order)
  @test sum(∫( e*e )dΩ) < 1.0e-9
end

mshfile = joinpath(@__DIR__,"t1.msh")
model = GmshDiscreteModel(mshfile)
test_discrete_model(model)
writevtk(model,"model")
check_interpolation(model)

mshfile = joinpath(@__DIR__,"twoTetraeder.msh")
model = GmshDiscreteModel(mshfile; renumber=true)
test_discrete_model(model)
check_interpolation(model)

mshfile = joinpath(@__DIR__,"..","demo","demo.msh")
model = GmshDiscreteModel(mshfile)
test_discrete_model(model)
check_interpolation(model)

mshfile = joinpath(@__DIR__,"square.msh")
model = GmshDiscreteModel(mshfile)
test_discrete_model(model)
check_interpolation(model)

mshfile = joinpath(@__DIR__,"cube.msh")
model = GmshDiscreteModel(mshfile)
test_discrete_model(model)
check_interpolation(model)

end # module
