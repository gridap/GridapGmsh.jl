using Gridap
using GridapGmsh

model = GmshDiscreteModel("demo/demo.msh")

order = 1
dirichlet_tags = ["boundary1","boundary2"]
u_boundary1(x) = 0.0
u_boundary2(x) = 1.0

reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,dirichlet_tags=dirichlet_tags)
U = TrialFESpace(V,[u_boundary1,u_boundary2])

Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)

a(u,v) = ∫( ∇(u)⋅∇(v) )dΩ
l(v) = 0

op = AffineFEOperator(a,l,U,V)
uh = solve(op)

writevtk(Ω,"demo",cellfields=["uh"=>uh])
