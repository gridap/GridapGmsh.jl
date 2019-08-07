
function _read_for_lagrangian_fe_space(mshfile)

  if !isfile(mshfile)
    error("Msh file not found: $mshfile")
  end

  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 1)
  gmsh.open(mshfile)

  D = _setup_dim(gmsh)

  node_to_coords = _setup_node_coords(gmsh,D)

  grid, cell_to_entity = _setup_grid(gmsh,D,node_to_coords)

  dim_to_offset = _setup_dim_to_offset(gmsh)

  nnodes = length(node_to_coords)
  node_to_label = _setup_node_to_label(gmsh,dim_to_offset,nnodes)

  dim_to_group_to_entities = _setup_dim_to_group_to_entities(gmsh)

  dim_to_group_to_name = _setup_dim_to_group_to_name(gmsh)

  tag_to_name, tag_to_groups, name_to_tag =
  
  _setup_tag_to_name(dim_to_group_to_name)

  tag_to_labels = _setup_tag_to_labels(
    tag_to_groups,dim_to_group_to_entities,dim_to_offset)

  gmsh.finalize()

  (grid, node_to_label, tag_to_labels, name_to_tag)

end

function GmshCLagrangianFESpace(
  ::Type{T},
  mshfile::String,
  dirinames::Vector{String},
  dirimasks::Vector{B}) where {D,Z,T,B}

  grid, node_to_label, tag_to_labels, name_to_tag =
    _read_for_lagrangian_fe_space(mshfile)

  diritags = [name_to_tag[name] for name in dirinames]

  CLagrangianFESpace(
    T, grid, node_to_label, tag_to_labels, diritags, dirimasks)

end

function GmshDLagrangianFESpace(
  ::Type{T},
  mshfile::String,
  dirinames::Vector{String},
  dirimasks::Vector{B}) where {D,Z,T,B}

  grid, node_to_label, tag_to_labels, name_to_tag =
    _read_for_lagrangian_fe_space(mshfile)

  diritags = [name_to_tag[name] for name in dirinames]

  DLagrangianFESpace(
    T, grid, node_to_label, tag_to_labels, diritags, dirimasks)

end

