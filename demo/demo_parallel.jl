using Gridap
using GridapGmsh
using GridapPETSc
using GridapDistributed
using PartitionedArrays
n = 6
with_mpi() do distribute 
  ranks=distribute(LinearIndices((n,)))
  options = "-ksp_type cg -pc_type gamg -ksp_monitor"
  GridapPETSc.with(args=split(options)) do
    model = GmshDiscreteModel(ranks,"demo/demo.msh")
    order = 1
    dirichlet_tags = ["boundary1","boundary2"]
    u_boundary1(x) = 0.0
    u_boundary2(x) = 1.0
    reffe = ReferenceFE(lagrangian,Float64,order)
    V = TestFESpace(model,reffe,dirichlet_tags=dirichlet_tags)
    U = TrialFESpace(V,[u_boundary1,u_boundary2])
    Ω = Interior(model)
    dΩ = Measure(Ω,2*order)
    a(u,v) = ∫( ∇(u)⋅∇(v) )dΩ
    l(v) = 0
    op = AffineFEOperator(a,l,U,V)
    solver = PETScLinearSolver()
    uh = solve(solver,op)
    writevtk(Ω,"demo",cellfields=["uh"=>uh])
  end
end
