
const D3=3
const POINT=15
const UNSET = 0

function GmshDiscreteModel(mshfile; renumber=true)
  @check_if_loaded
  if !isfile(mshfile)
    error("Msh file not found: $mshfile")
  end

  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 1)
  gmsh.option.setNumber("Mesh.SaveAll", 1)
  gmsh.option.setNumber("Mesh.MedImportGroupsOfNodes", 1)
  gmsh.open(mshfile)

  renumber && gmsh.model.mesh.renumberNodes()
  renumber && gmsh.model.mesh.renumberElements()

  Dc = _setup_cell_dim(gmsh)
  Dp = _setup_point_dim(gmsh,Dc)
  node_to_coords = _setup_node_coords(gmsh,Dp)
  nnodes = length(node_to_coords)
  vertex_to_node, node_to_vertex = _setup_nodes_and_vertices(gmsh,node_to_coords)
  grid, cell_to_entity = _setup_grid(gmsh,Dc,Dp,node_to_coords,node_to_vertex)
  cell_to_vertices = _setup_cell_to_vertices(get_cell_node_ids(grid),node_to_vertex,nnodes)
  grid_topology = UnstructuredGridTopology(grid,cell_to_vertices,vertex_to_node)
  labeling = _setup_labeling(gmsh,grid,grid_topology,cell_to_entity,vertex_to_node,node_to_vertex)
  gmsh.finalize()

  UnstructuredDiscreteModel(grid,grid_topology,labeling)
end

function  _setup_nodes_and_vertices(gmsh,node_to_coords)
  nnodes = length(node_to_coords)
  dimTags = gmsh.model.getEntities()
  if _has_periodic_bcs(gmsh,dimTags)
    dimTags = gmsh.model.getEntities()
    vertex_to_node, node_to_vertex = _setup_nodes_and_vertices_periodic(gmsh,dimTags,nnodes)
  else
    vertex_to_node = 1:nnodes
    node_to_vertex = vertex_to_node
  end
  vertex_to_node, node_to_vertex
end

function _setup_nodes_and_vertices_periodic(gmsh,dimTags,nnodes)
  # Assumes linear grid
  node_to_node_master = fill(UNSET,nnodes)
  _node_to_node_master!(node_to_node_master,gmsh,dimTags)
  slave_to_node_slave = findall(node_to_node_master .!= UNSET)
  slave_to_node_master = node_to_node_master[slave_to_node_slave]
  node_to_vertex = fill(UNSET,nnodes)
  vertex_to_node = findall(node_to_node_master .== UNSET)
  node_to_vertex[vertex_to_node] = 1:length(vertex_to_node)
  nmax = 20
  for i in 1:nmax
    node_to_vertex[slave_to_node_slave] = node_to_vertex[slave_to_node_master]
    if all(j->j!=0,node_to_vertex)
      break
    end
    if i == nmax
      @unreachable
    end
  end
  vertex_to_node, node_to_vertex
end

function _setup_grid(gmsh,Dc,Dp,node_to_coords,node_to_vertex)

  if Dp == 3 && Dc == 2
    orient_if_simplex = false
  else
    orient_if_simplex = true
  end

  cell_to_nodes, nminD = _setup_connectivity(gmsh,Dc,node_to_vertex,orient_if_simplex)
  cell_to_type, reffes, orientation = _setup_reffes(gmsh,Dc,orient_if_simplex)
  cell_to_entity = _setup_cell_to_entity(
    gmsh,Dc,length(cell_to_nodes),nminD)

  if Dp == 3 && Dc == 2
    cell_coords = lazy_map(Broadcasting(Reindex(node_to_coords)),cell_to_nodes)
    ctype_shapefuns = map(get_shapefuns,reffes)
    cell_shapefuns = expand_cell_data(ctype_shapefuns,cell_to_type)
    cell_map = lazy_map(linear_combination,cell_coords,cell_shapefuns)
    ctype_x = fill(zero(VectorValue{Dc,Float64}),length(ctype_shapefuns))
    cell_x = expand_cell_data(ctype_x,cell_to_type)
    cell_Jt = lazy_map(∇,cell_map)
    cell_n = lazy_map(Operation(_unit_outward_normal),cell_Jt)
    cell_nx = lazy_map(evaluate,cell_n,cell_x) |> collect
    facet_normal = lazy_map(constant_field,cell_nx)
  else
    facet_normal = nothing
  end

  grid = UnstructuredGrid(
    node_to_coords,
    cell_to_nodes,
    reffes,
    cell_to_type,
    orientation,
    facet_normal)

  (grid, cell_to_entity)

end

function _unit_outward_normal(v::MultiValue{Tuple{2,3}})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1,n2,n3)
  n/norm(n)
end

function _setup_cell_to_vertices(cell_to_nodes,node_to_vertex,nnodes)
  if isa(node_to_vertex,AbstractVector)
    cell_to_vertices = Table(lazy_map(Broadcasting(Reindex(node_to_vertex)),cell_to_nodes))
  else
    @assert node_to_vertex == 1:nnodes
    cell_to_vertices = cell_to_nodes
  end
  cell_to_vertices
end

function _has_periodic_bcs(gmsh,dimTags)
  for (dim,tag) in dimTags
    tagMaster, nodeTags, = gmsh.model.mesh.getPeriodicNodes(dim,tag)
    if length(nodeTags) > 0
      return true
    end
  end
  return false
end

function _node_to_node_master!(node_to_node_master,gmsh,dimTags)
  for (dim,tag) in dimTags
    tagMaster, nodeTags, nodeTagsMaster, = gmsh.model.mesh.getPeriodicNodes(dim,tag)
    if length(nodeTags) > 0
      node_to_node_master[nodeTags] .= nodeTagsMaster
    end
  end
end

function _setup_cell_dim(gmsh)
  entities = gmsh.model.getEntities()
  D = -1
  for e in entities
    D = max(D,e[1])
  end
  if D == -1
    gmsh.finalize()
    error("No entities in the msh file.")
  end
  D
end

function _setup_point_dim(gmsh,Dc)
  if Dc == D3
    return Dc
  end
  nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()
  for node in nodeTags
    j = D3
    k = (node-1)*D3 + j
    xj = coord[k]
    if !(xj + 1 ≈ 1)
      return D3
    end
  end
  Dc
end

function _setup_node_coords(gmsh,D)
  nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()
  nmin = minimum(nodeTags)
  nmax = maximum(nodeTags)
  nnodes = length(nodeTags)
  if !(nmax == nnodes && nmin == 1)
    gmsh.finalize()
    error("Only consecutive node tags allowed.")
  end
  node_to_coords = zeros(Point{D,Float64},nnodes)
  _fill_node_coords!(node_to_coords,nodeTags,coord,D)
  node_to_coords
end

function _fill_node_coords!(node_to_coords,nodeTags,coord,D)
  m = zero(Mutable(Point{D,Float64}))
  for node in nodeTags
    for j in 1:D
      k = (node-1)*D3 + j
      xj = coord[k]
      m[j] = xj
    end
    node_to_coords[node] = m
  end
end

function _setup_connectivity(gmsh,d,node_to_vertex,orient_if_simplex)

  elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)

  if length(elemTypes) == 0
    ncells = 0
    ndata = 0
    nmin = 1
    cell_to_nodes_prts = zeros(Int,ncells+1)
    cell_to_nodes_data = zeros(Int32,ndata)
    cell_to_nodes = Table(cell_to_nodes_data,cell_to_nodes_prts)
    return (cell_to_nodes, nmin)
  end

  ncells, nmin, nmax = _check_cell_tags(elemTags)

  etype_to_nlnodes = _setup_etype_to_nlnodes(elemTypes,gmsh)

  ndata = sum([length(t) for t in nodeTags])

  cell_to_nodes_data = zeros(Int,ndata)
  cell_to_nodes_prts = zeros(Int32,ncells+1)

  _fill_connectivity!(
    cell_to_nodes_data,
    cell_to_nodes_prts,
    nmin-1,
    etype_to_nlnodes,
    elemTypes,
    elemTags,
    nodeTags,
    d,
    node_to_vertex,
    orient_if_simplex)

  cell_to_nodes = Table(cell_to_nodes_data,cell_to_nodes_prts)

  (cell_to_nodes, nmin)

end

function _check_cell_tags(elemTags)
  nmin::Int = minimum([minimum(t) for t in elemTags])
  nmax::Int = maximum([maximum(t) for t in elemTags])
  ncells = sum([length(t) for t in elemTags])
  if !( (nmax-nmin+1) == ncells)
    gmsh.finalize()
    error("Only consecutive elem tags allowed.")
  end
  (ncells,nmin,nmax)
end

function _setup_etype_to_nlnodes(elemTypes,gmsh)
  netypes = maximum(elemTypes)
  etype_to_nlnodes = zeros(Int,netypes)
  for etype in elemTypes
    name, dim, order, numv, parv =
      gmsh.model.mesh.getElementProperties(etype)
    etype_to_nlnodes[etype] = numv
  end
  etype_to_nlnodes
end

function  _fill_connectivity!(
    cell_to_nodes_data,
    cell_to_nodes_prts,
    o,
    etype_to_nlnodes,
    elemTypes,
    elemTags,
    nodeTags,
    d,
    node_to_vertex,
    orient_if_simplex)

  for (j,etype) in enumerate(elemTypes)
    nlnodes = etype_to_nlnodes[etype]
    i_to_cell = elemTags[j]
    for cell in i_to_cell
      cell_to_nodes_prts[cell+1-o] = nlnodes
    end
  end

  length_to_ptrs!(cell_to_nodes_prts)

  for (j,etype) in enumerate(elemTypes)
    nlnodes = etype_to_nlnodes[etype]
    i_to_cell = elemTags[j]
    i_lnode_to_node = nodeTags[j]
    if (nlnodes == d+1) && orient_if_simplex
      # what we do here has to match with the OrientationStyle we
      # use when building the UnstructuredGrid
      _orient_simplex_connectivities!(nlnodes,i_lnode_to_node,node_to_vertex)
    elseif (nlnodes == 4)
      _sort_quad_connectivites!(nlnodes,i_lnode_to_node)
    elseif (nlnodes == 8)
      _sort_hex_connectivites!(nlnodes,i_lnode_to_node)
    end
    for (i,cell) in enumerate(i_to_cell)
      a = cell_to_nodes_prts[cell-o]-1
      for lnode in 1:nlnodes
        node = i_lnode_to_node[(i-1)*nlnodes+lnode]
        cell_to_nodes_data[a+lnode] = node
      end
    end
  end

end

function _orient_simplex_connectivities!(nlnodes,i_lnode_to_node,node_to_vertex)
  aux = zeros(eltype(i_lnode_to_node),nlnodes)
  offset = nlnodes-1
  for i in 1:nlnodes:length(i_lnode_to_node)
    nodes = i_lnode_to_node[i:i+offset]
    vertices = view(node_to_vertex,nodes)
    perm = sortperm(vertices)
    i_lnode_to_node[i:i+offset] = nodes[perm]
  end
end

function _sort_quad_connectivites!(nlnodes,i_lnode_to_node)
  aux = zero(eltype(i_lnode_to_node))
  offset = nlnodes-1
  for i in 1:nlnodes:length(i_lnode_to_node)
    aux = i_lnode_to_node[i+offset-1]
    i_lnode_to_node[i+offset-1] = i_lnode_to_node[i+offset]
    i_lnode_to_node[i+offset] = aux
  end
end

function _sort_hex_connectivites!(nlnodes,i_lnode_to_node)
  perm = [1, 2, 4, 3, 5, 6, 8, 7]
  aux = zeros(eltype(i_lnode_to_node),nlnodes)
  offset = nlnodes-1
  for i in 1:nlnodes:length(i_lnode_to_node)
    aux = i_lnode_to_node[i:i+offset]
    i_lnode_to_node[i:i+offset] = aux[perm]
  end
end

function _setup_reffes(gmsh,d,orient_if_simplex)

  elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)

  if !(length(elemTypes)==1)
    gmsh.finalize()
    s = """
    Only one element type per dimension allowed for the moment.
    Dimension $d has $(length(elemTypes)) different element types
    """
    error(s)
  end

  ncells, nmin, nmax = _check_cell_tags(elemTags)
  cell_to_type = fill(Int8(1),ncells)

  etype::Int = first(elemTypes)
  name, dim, order::Int, numv, parv = gmsh.model.mesh.getElementProperties(etype)

  if order == 0 && etype == POINT
    order = 1
  end

  if order != 1
    gmsh.finalize()
    error("For the moment only for first-order elements")
  end

  reffe = _reffe_from_etype(etype)
  reffes = [reffe,]

  boo = is_simplex(get_polytope(reffe)) && orient_if_simplex

  orientation = boo ? Oriented() : NonOriented()

  (cell_to_type, reffes, orientation)
end

function _reffe_from_etype(eltype)
  if eltype == 1
    SEG2
  elseif eltype == 2
    TRI3
  elseif eltype == 3
    QUAD4
  elseif eltype == 4
    TET4
  elseif eltype == 5
    HEX8
  elseif eltype == POINT
    VERTEX1
  else
    gmsh.finalize()
    error("Unsupported element type. elemType: $etype")
  end
end

function _setup_cell_to_entity(gmsh,d,ncells,nmin)
  cell_to_entity = fill(UNSET,ncells)
  entities = gmsh.model.getEntities(d)
  for e in entities
    _, elemTags, _ = gmsh.model.mesh.getElements(e[1], e[2])
    _fill_cell_to_entity!(cell_to_entity,elemTags,e[2],nmin-1)
  end
  cell_to_entity
end

function _fill_cell_to_entity!(cell_to_entity,elemTags,entity,o)
  for i_to_cell in elemTags
    for cell in i_to_cell
      cell_to_entity[cell-o] = entity
    end
  end
end

function _setup_labeling(gmsh,grid,grid_topology,cell_to_entity,vertex_to_node,node_to_vertex)

  D = num_cell_dims(grid)
  dim_to_gface_to_nodes, dim_gface_to_entity = _setup_faces(gmsh,D,node_to_vertex)
  push!(dim_to_gface_to_nodes,get_cell_node_ids(grid))
  push!(dim_gface_to_entity,cell_to_entity)

  dim_to_group_to_entities = _setup_dim_to_group_to_entities(gmsh)
  dim_to_group_to_name = _setup_dim_to_group_to_name(gmsh)
  dim_to_offset = _setup_dim_to_offset(gmsh)

  nnodes = num_nodes(grid)
  node_to_label = _setup_node_to_label(gmsh,dim_to_offset,nnodes)

  labeling = _setup_face_labels(
    vertex_to_node,
    node_to_vertex,
    grid_topology,
    dim_to_gface_to_nodes,
    dim_gface_to_entity,
    dim_to_offset,
    dim_to_group_to_entities,
    dim_to_group_to_name,
    node_to_label)

  labeling

end

function _setup_faces(gmsh,D,node_to_vertex)
  dim_to_gface_to_nodes = []
  dim_gface_to_entity = []
  for d in 0:(D-1)
    orient_if_simplex = true
    face_to_nodes, nmin  = _setup_connectivity(gmsh,d,node_to_vertex,orient_if_simplex)
    face_to_entity = _setup_cell_to_entity(gmsh,d,length(face_to_nodes),nmin)
    push!(dim_to_gface_to_nodes,face_to_nodes)
    push!(dim_gface_to_entity,face_to_entity)
  end
  (dim_to_gface_to_nodes, dim_gface_to_entity)
end

function _setup_dim_to_group_to_entities(gmsh)
  dim_to_group_to_entities = Vector{Vector{Int}}[]
  for d in 0:D3
    pgs = gmsh.model.getPhysicalGroups(d)
    n = length(pgs)
    group_to_entities = Vector{Int}[]
    for pg in pgs
      es = gmsh.model.getEntitiesForPhysicalGroup(pg...)
      entities = Int[]
      for e in es
        push!(entities,e)
      end
      push!(group_to_entities,entities)
    end
    push!(dim_to_group_to_entities,group_to_entities)
  end
  dim_to_group_to_entities
end

function _setup_dim_to_group_to_name(gmsh)
  u = 1
  dim_to_group_to_name = Vector{String}[]
  for d in 0:D3
    pgs = gmsh.model.getPhysicalGroups(d)
    n = length(pgs)
    names = String[]
    for pg in pgs
      _name = gmsh.model.getPhysicalName(pg...)
      if _name == ""
        name = "untitled$u"
        u += 1
      else
        name = _name
      end
      push!(names,name)
    end
    push!(dim_to_group_to_name,names)
  end
  dim_to_group_to_name
end

function _setup_dim_to_offset(gmsh)
  entities = gmsh.model.getEntities()
  dim_to_nentities = zeros(Int,D3+1)
  for e in entities
    d = e[1]
    id = e[2]
    _id = dim_to_nentities[d+1]
    dim_to_nentities[d+1] = max(id,_id)
  end
  dim_to_offset = zeros(Int,D3+1)
  for d=1:D3
    dim_to_offset[d+1] += dim_to_nentities[d-1+1] + dim_to_offset[d-1+1]
  end
  dim_to_offset
end

function _setup_node_to_label(gmsh,dim_to_offset,nnodes)
  node_to_label = zeros(Int,nnodes)
  entities = gmsh.model.getEntities()
  for (dim, entity) in entities
    nodes, _, _ = gmsh.model.mesh.getNodes(dim, entity)
    node_to_label[nodes] .= (entity + dim_to_offset[dim+1])
  end
  node_to_label
end

function _setup_face_labels(
  vertex_to_node,
  node_to_vertex,
  grid_topology,
  dim_to_gface_to_nodes,
  dim_gface_to_entity,
  dim_to_offset,
  dim_to_group_to_entities,
  dim_to_group_to_name,
  node_to_label)

  vertex_to_label = node_to_label[vertex_to_node]

  D = num_cell_dims(grid_topology)

  dim_to_face_to_label = [vertex_to_label,]

  for d in 1:D
    z = fill(UNSET,num_faces(grid_topology,d))
    push!(dim_to_face_to_label,z)
  end

  _fill_dim_to_face_to_label!(
    node_to_vertex,
    dim_to_face_to_label,
    grid_topology,
    dim_to_gface_to_nodes,
    dim_gface_to_entity,
    dim_to_offset)

  tag_to_name, tag_to_groups, _ = _setup_tag_to_name(dim_to_group_to_name)

  tag_to_labels = _setup_tag_to_labels(
    tag_to_groups,dim_to_group_to_entities,dim_to_offset)

  FaceLabeling(dim_to_face_to_label,tag_to_labels,tag_to_name)

end

function _fill_dim_to_face_to_label!(
  node_to_vertex,
  dim_to_face_to_label,
  grid_topology,
  dim_to_gface_to_nodes,
  dim_gface_to_entity,
  dim_to_offset)

  D = length(dim_to_face_to_label)-1
  d = D
  cell_to_entity =  dim_gface_to_entity[d+1]
  cell_to_label = dim_to_face_to_label[d+1]
  offset = dim_to_offset[d+1]
  _apply_offset_for_cells!(cell_to_label,cell_to_entity,offset)

  for d in 0:(D-1)

    gface_to_nodes = dim_to_gface_to_nodes[d+1]
    face_to_nodes = get_faces(grid_topology,d,0)
    node_to_faces = get_faces(grid_topology,0,d)
    gface_to_face = _setup_gface_to_face(
      node_to_vertex,
      face_to_nodes,
      node_to_faces,
      gface_to_nodes)
    face_to_label = dim_to_face_to_label[d+1]
    gface_to_entity = dim_gface_to_entity[d+1]
    offset = dim_to_offset[d+1]
    _apply_offset_for_faces!(face_to_label,gface_to_entity,gface_to_face,offset)

  end

  for d = 0:(D-1)
    for j in (d+1):D
      dface_to_jfaces = get_faces(grid_topology,d,j)
      dface_to_label = dim_to_face_to_label[d+1]
      jface_to_label = dim_to_face_to_label[j+1]
      _fix_dface_to_label!(dface_to_label,jface_to_label,dface_to_jfaces)
    end
  end

end

function _apply_offset_for_cells!(cell_to_label,cell_to_entity,offset)
  ncells = length(cell_to_label)
  for cell in 1:ncells
    entity = cell_to_entity[cell]
    cell_to_label[cell] = entity+offset
  end
end

function _apply_offset_for_faces!(face_to_label,gface_to_entity,gface_to_face,offset)
  ngfaces = length(gface_to_face)
  for gface in 1:ngfaces
    entity = gface_to_entity[gface]
    face = gface_to_face[gface]
    face_to_label[face] = entity+offset
  end
end

function _setup_gface_to_face(
  node_to_vertex,
  face_to_nodes,
  node_to_faces,
  gface_to_nodes)

  gface_to_face = _find_gface_to_face(
    node_to_vertex,
    face_to_nodes.data,
    face_to_nodes.ptrs,
    node_to_faces.data,
    node_to_faces.ptrs,
    gface_to_nodes.data,
    gface_to_nodes.ptrs)

  gface_to_face

end

function _find_gface_to_face(
  node_to_vertex,
  face_to_nodes_data,
  face_to_nodes_ptrs,
  node_to_faces_data::AbstractVector{T},
  node_to_faces_ptrs,
  gface_to_nodes_data,
  gface_to_nodes_ptrs) where T

  ngfaces = length(gface_to_nodes_ptrs) - 1
  gface_to_face = fill(T(UNSET),ngfaces)
  n = max_cells_arround_vertex(node_to_faces_ptrs)
  faces_around = fill(UNSET,n)
  faces_around_scratch = fill(UNSET,n)

  _fill_gface_to_face!(
    gface_to_face,
    node_to_vertex,
    face_to_nodes_data,
    face_to_nodes_ptrs,
    node_to_faces_data,
    node_to_faces_ptrs,
    gface_to_nodes_data,
    gface_to_nodes_ptrs,
    faces_around,
    faces_around_scratch)

  gface_to_face

end

function  _fill_gface_to_face!(
  gface_to_face,
  node_to_vertex,
  face_to_nodes_data,
  face_to_nodes_ptrs,
  node_to_faces_data,
  node_to_faces_ptrs,
  gface_to_nodes_data,
  gface_to_nodes_ptrs,
  faces_around,
  faces_around_scratch)

  ngfaces = length(gface_to_nodes_ptrs) - 1

  nfaces_around = UNSET
  nfaces_around_scratch = UNSET

  for gface in 1:ngfaces

    a = gface_to_nodes_ptrs[gface]-1
    b = gface_to_nodes_ptrs[gface+1]
    nlnodes = b-(a+1)

    for lnode in 1:nlnodes
      node = gface_to_nodes_data[lnode+a]
      vertex = node_to_vertex[node]
      if lnode == 1
        nfaces_around = _fill_cells_around_scratch!(
          faces_around,
          vertex,
          node_to_faces_data,
          node_to_faces_ptrs)
      else
        nfaces_around_scratch = _fill_cells_around_scratch!(
          faces_around_scratch,
          vertex,
          node_to_faces_data,
          node_to_faces_ptrs)
        _set_intersection!(
          faces_around,faces_around_scratch,
          nfaces_around,nfaces_around_scratch)
      end
    end

    for face in faces_around
      if face != UNSET
        gface_to_face[gface] = face
        break
      end
    end

  end

end

function _fix_dface_to_label!(dface_to_label,jface_to_label,dface_to_jfaces)

  ndfaces = length(dface_to_label)
  @assert ndfaces == length(dface_to_jfaces)

  for dface in 1:ndfaces

    dlabel = dface_to_label[dface]
    if dlabel != UNSET
      continue
    end

    jfaces = dface_to_jfaces[dface]
    for jface in jfaces
      jlabel = jface_to_label[jface]
      if jlabel != UNSET
        dface_to_label[dface] = jlabel
        break
      end
    end

  end

end

function _setup_tag_to_name(dim_to_group_to_name)
  tag_to_name = String[]
  tag_to_groups = Vector{Tuple{Int,Int}}[]
  name_to_tag = Dict{String,Int}()
  tag = 1
  for (dim, group_to_name) in enumerate(dim_to_group_to_name)
    for (id,name) in enumerate(group_to_name)
      group = (dim,id)
      if !haskey(name_to_tag,name)
        push!(tag_to_name,name)
        groups = [group,]
        push!(tag_to_groups,groups)
        name_to_tag[name] = tag
        tag += 1
      else
        _tag = name_to_tag[name]
        groups = tag_to_groups[_tag]
        push!(groups,group)
      end
    end
  end
  (tag_to_name, tag_to_groups, name_to_tag)
end

function _setup_tag_to_labels(
  tag_to_groups,dim_to_group_to_entities,dim_to_offset)

  tag_to_labels = Vector{Int}[]
  for groups in tag_to_groups

    labels = Int[]
    for group in groups
      dim, id = group
      offset = dim_to_offset[dim]
      entities = dim_to_group_to_entities[dim][id]
      for entity in entities
        label = entity + offset
        push!(labels,label)
      end
    end
    push!(tag_to_labels,labels)

  end

  tag_to_labels

end

## Parallel related

function GmshDiscreteModel(parts::PArrays.AbstractPData, args...;kwargs...)
  GmshDiscreteModel(parts,args...;kwargs...) do g,np
    if np == 1
      fill(Int32(1),size(g,1))
    else
      Metis.partition(g,np)
    end
  end
end

function GmshDiscreteModel(
  do_partition,
  parts::PArrays.AbstractPData,
  args...;kwargs...)

  smodel = GmshDiscreteModel(args...;kwargs...)
  g = GridapDistributed.compute_cell_graph(smodel)
  np = length(parts)
  cell_to_part = do_partition(g,np)
  DiscreteModel(parts,smodel,cell_to_part,g)
end

# Native serial Gridap if not parts given
function GmshDiscreteModel(parts::Nothing, args...;kwargs...)
  GmshDiscreteModel(args...;kwargs...)
end
