module GmshDiscreteModelsTests

using Test
using GridapGmsh
using Gridap

d = mktempdir()
f = joinpath(d,"model")
f = "model"

mshfile = joinpath(@__DIR__,"t1.msh")
model = GmshDiscreteModel(mshfile)
writevtk(model,f)

@test filesize("$(f)_0.vtu") > 0
@test filesize("$(f)_1.vtu") > 0
@test filesize("$(f)_2.vtu") > 0

mshfile = joinpath(@__DIR__,"../demo/","demo.msh")
model = GmshDiscreteModel(mshfile)
writevtk(model,f)

@test filesize("$(f)_0.vtu") > 0
@test filesize("$(f)_1.vtu") > 0
@test filesize("$(f)_2.vtu") > 0
@test filesize("$(f)_3.vtu") > 0

rm(d,recursive=true)

end # module
