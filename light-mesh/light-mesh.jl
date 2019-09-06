
using Gridap
using GridapGmsh
using GridapGmsh: gmsh

model = GmshDiscreteModel("light-mesh/light-mesh.msh")

mshfile = "light-mesh/light-mesh.msh"
if !isfile(mshfile)
  error("Msh file not found: $mshfile")
end

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.open(mshfile)

D = _setup_dim(gmsh)

node_to_coords = _setup_node_coords(gmsh,D)

grid, cell_to_entity = _setup_grid(gmsh,D,node_to_coords)

dim_to_gface_to_nodes, dim_gface_to_entity = _setup_faces(gmsh,D)

push!(dim_to_gface_to_nodes,cells(grid))
push!(dim_gface_to_entity,cell_to_entity)

dim_to_group_to_entities = _setup_dim_to_group_to_entities(gmsh)

dim_to_group_to_name = _setup_dim_to_group_to_name(gmsh)

dim_to_offset = _setup_dim_to_offset(gmsh)

nnodes = length(node_to_coords)
node_to_label = _setup_node_to_label(gmsh,dim_to_offset,nnodes)

gmsh.finalize()

graph = FullGridGraph(grid)

facelabels = _setup_face_labels(
graph,
dim_to_gface_to_nodes,
dim_gface_to_entity,
dim_to_offset,
dim_to_group_to_entities,
dim_to_group_to_name,
node_to_label)

grids = _setup_model_grids(graph,grid,D)

DiscreteModelFromData( grids, graph, facelabels)

model.grids[2]

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



function _orient_simplex_connectivities(cell_to_nodes_data, cell_to_nodes_ptrs)
  cell_nodes = Gridap.CellValuesGallery.CellVectorFromDataAndPtrs(cell_to_nodes_data,cell_to_nodes_ptrs)
  oriented = Int64[]
  for i in cell_nodes
    push!(oriented,sort(i)...)
  end
  return oriented
end
