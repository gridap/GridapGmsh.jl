module GmshDiscreteModelsTests

using GridapGmsh
using Gridap.Geometry

mshfile = joinpath(@__DIR__,"t1.msh")
model = GmshDiscreteModel(mshfile)
test_discrete_model(model)

mshfile = joinpath(@__DIR__,"twoTetraeder.msh")
model = GmshDiscreteModel(mshfile; renumber=true)
test_discrete_model(model)

mshfile = joinpath(@__DIR__,"..","demo","demo.msh")
model = GmshDiscreteModel(mshfile)
test_discrete_model(model)

mshfile = joinpath(@__DIR__,"square.msh")
model = GmshDiscreteModel(mshfile)
test_discrete_model(model)

mshfile = joinpath(@__DIR__,"cube.msh")
model = GmshDiscreteModel(mshfile)
test_discrete_model(model)

end # module
