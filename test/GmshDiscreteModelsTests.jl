module GmshDiscreteModelsTests

using Gridap
using GridapGmsh
using Gridap.Geometry
using Gridap.Helpers
using Test

function check_interpolation(model)
  order = 2
  reffe = ReferenceFE(lagrangian,Float64,order)
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

mshfile = joinpath(@__DIR__,"plane.msh")
model = GmshDiscreteModel(mshfile)
grid = get_grid(model)
@test num_cell_dims(grid) == 2
@test num_point_dims(grid) == 3
Ω = Triangulation(model)
n_Ω = get_normal_vector(Ω)
x_Ω = get_cell_points(Ω)
collect(n_Ω(x_Ω))
test_discrete_model(model)
check_interpolation(model)

mshfile = joinpath(@__DIR__,"periodic.msh")
model = GmshDiscreteModel(mshfile)
test_discrete_model(model)
u(x) = 2*x[1]
reffe = ReferenceFE(lagrangian,Float64,1)
V = TestFESpace(model,reffe,dirichlet_tags=["L","PL"])
U = TrialFESpace(V,u)
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model,tags="R")
n_Γ = get_normal_vector(Γ)
dΩ = Measure(Ω,2)
dΓ = Measure(Γ,2)
a(u,v) = ∫( ∇(u)⋅∇(v) )dΩ
l(v) = ∫( n_Γ⋅(v*∇(u)) )dΓ
op = AffineFEOperator(a,l,U,V)
uh = solve(op)
#writevtk(Ω,"sol",cellfields=["uh"=>uh,"u"=>u,"e"=>u-uh,"zh"=>zero(U)])
@test sum( ∫(abs2(u - uh))dΩ ) < 1.0e-9

mshfile = joinpath(@__DIR__,"full_periodic.msh")
model = GmshDiscreteModel(mshfile)
test_discrete_model(model)
model2 = Gridap.DiscreteModelFromFile(mshfile)
test_discrete_model(model2)

end # module
