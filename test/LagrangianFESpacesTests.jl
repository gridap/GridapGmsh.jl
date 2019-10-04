module LagrangianFESpacesTests

using Test
using Gridap
using GridapGmsh

mshfile = joinpath(@__DIR__,"../demo/","demo.msh")

T = VectorValue{3,Float64}

dirinames = ["boundary1", "boundary2"]
dirimasks = [(true,false,true), (true,true,false)]

ufun(x) = x

a(v,u) = inner(e,e)

# CLagrangianFESpace

fespace = GmshCLagrangianFESpace(T, mshfile, dirinames, dirimasks)

grid = fespace.grid

trian = Triangulation(grid)
quad = CellQuadrature(trian,degree=2)

uh = interpolate(fespace,ufun)

u = CellField(trian,ufun)

e = u - uh

tol = 1.e-8
@test sum( integrate( a(e,e), trian, quad) ) < tol

# DLagrangianFESpace

fespace = GmshDLagrangianFESpace(T, mshfile, dirinames, dirimasks)

grid = fespace.grid

trian = Triangulation(grid)
quad = CellQuadrature(trian,degree=2)

uh = interpolate(fespace,ufun)

u = CellField(trian,ufun)

e = u - uh

tol = 1.e-8
@test sum( integrate( a(e,e), trian, quad) ) < tol

end # module
