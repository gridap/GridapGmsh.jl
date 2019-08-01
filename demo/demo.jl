
using Gridap
using GridapGmsh

model = GmshDiscreteModel("demo/demo.msh")

order = 1
diritags = ["boundary1", "boundary2"]
fespace = ConformingFESpace(Float64,model,order,diritags)

ufun1(x) = 0.0
ufun2(x) = 1.0
V = TestFESpace(fespace)
U = TrialFESpace(fespace,[ufun1,ufun2])

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

a(v,u) = inner(∇(v), ∇(u))
t_Ω = LinearFETerm(a,trian,quad)

assem = SparseMatrixAssembler(V,U)

op = LinearFEOperator(V,U,assem,t_Ω)

ls = LUSolver()
solver = LinearFESolver(ls)

uh = solve(solver,op)
writevtk(trian,"demo",cellfields=["uh"=>uh])
