module GmshDiscreteModelsTests

using Test
using GridapGmsh
using Gridap

mshfile = joinpath(@__DIR__,"t1.msh")

model = GmshDiscreteModel(mshfile)

d = mktempdir()
f = joinpath(d,"model")
#f = "model"

writevtk(model,f)

@test filesize("$(f)_0.vtu") > 0
@test filesize("$(f)_1.vtu") > 0
@test filesize("$(f)_2.vtu") > 0

rm(d,recursive=true)

end # module
