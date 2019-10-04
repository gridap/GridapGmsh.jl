module GmshDiscreteModelsTests

using Test
using GridapGmsh
using Gridap

d = mktempdir()
f = joinpath(d,"model")
#f = "model"

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

mshfile = joinpath(@__DIR__,"../demo/","demo.msh")
model = GmshDiscreteModel(mshfile)

order = 3

fespace = FESpace(
  reffe=:Lagrangian, order=order, valuetype=Float64,
  conformity=:H1, model=model)

U = TrialFESpace(fespace)

u(x) = x[1]*x[2]*x[3]

uh = interpolate(U,u)

trian = Triangulation(model)
quad = CellQuadrature(trian,degree=3*order)

e = u - uh

ce = integrate( inner(e,e), trian, quad)

#writevtk(trian,"trian",celldata=["ce"=>ce],cellfields=["e"=>e,"uh"=>uh,"u"=>CellField(trian,u)])

tol = 1.e-8
r = sum( ce  )
@test r < tol

end # module
