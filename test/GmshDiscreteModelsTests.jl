module GmshDiscreteModelsTests

using GridapGmsh
using GridapGmsh: gmsh
using Gridap
using Gridap.CellValuesGallery
using StaticArrays

const D3=3

function explore(mshfile)

  if !isfile(mshfile)
    error("Msh file not found: $mshfile")
  end

  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 1)
  gmsh.open(mshfile)

  D = 2
  node_to_coords = _setup_node_coords(gmsh,D)
  cell_to_nodes, cell_to_extrusion, cell_to_order = _setup_connectivity(gmsh,D)
  facet_to_nodes, facet_to_extrusion, facet_to_order = _setup_connectivity(gmsh,1)

  grid2 = UnstructuredGrid(
    node_to_coords,
    cell_to_nodes.data,
    cell_to_nodes.ptrs,
    cell_to_extrusion,
    cell_to_order)

  grid1 = UnstructuredGrid(
    node_to_coords,
    facet_to_nodes.data,
    facet_to_nodes.ptrs,
    facet_to_extrusion,
    facet_to_order)

  writevtk(grid2,"grid2")
  writevtk(grid1,"grid1")
  
  ## get all elementary entities in the model
  #entities = gmsh.model.getEntities()
  #
  #for e in entities
  #    # get the mesh nodes for each elementary entity
  #    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(e[1], e[2])
  #    # get the mesh elements for each elementary entity
  #    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
  #    # report some statistics
  #    numElem = sum([length(t) for t in elemTags])
  #    println(length(nodeTags), " mesh nodes and ", numElem, " mesh elements ",
  #            "on entity ", e, " of type ", gmsh.model.getType(e[1], e[2]))
  #    partitions = gmsh.model.getPartitions(e[1], e[2])
  #    if length(partitions) > 0
  #        println(" - Partition tag(s): ", partitions, " - parent entity",
  #                gmsh.model.getParent(e[1], e[2]))
  #    end
  #    for t in elemTypes
  #        name, dim, order, numv, parv = gmsh.model.mesh.getElementProperties(t)
  #        println(" - Element type: ", name, ", order ", order)
  #        println("   with ", numv, " nodes in param coord: ", parv)
  #    end
  #end
  
  gmsh.finalize()

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

  m = zero(MVector{D,Float64})

  for node in nodeTags
    for j in 1:D
      k = (node-1)*D3 + j
      xj = coord[k]
      m[j] = xj
    end
    node_to_coords[node] = m
  end

end

function _setup_connectivity(gmsh,d)

  elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)

  nmin::Int = minimum([minimum(t) for t in elemTags])
  nmax::Int = maximum([maximum(t) for t in elemTags])
  ncells = sum([length(t) for t in elemTags])
  if !( (nmax-nmin+1) == ncells)
    gmsh.finalize()
    error("Only consecutive elem tags allowed.")
  end
  ndata = sum([length(t) for t in nodeTags])
  etype_to_nlnodes = _setup_etype_to_nlnodes(elemTypes,gmsh)

  cell_to_nodes_prts = zeros(Int,ncells+1)
  cell_to_nodes_data = zeros(Int,ndata)

  _fill_connectivity!(
    cell_to_nodes_data,
    cell_to_nodes_prts,
    nmin-1,
    etype_to_nlnodes,
    elemTypes,
    elemTags,
    nodeTags)

  cell_to_nodes = CellVectorFromDataAndPtrs(cell_to_nodes_data,cell_to_nodes_prts)

  if !(length(elemTypes)==1)
    gmsh.finalize()
    s = """
    Only one element type per dimension allowed for the moment.
    Dimension $d has $(length(elemTypes)) different element types
    """
    error(s)
  end

  etype::Int = elemTypes[1]
  name, dim, order::Int, numv, parv = gmsh.model.mesh.getElementProperties(etype)
  extrusion = _extrussion_from_etype(Val(etype))
  if extrusion == nothing
    gmsh.finalize()
    error("Unsupported element type. elemType: $etype")
  end

  cell_to_extrusion = ConstantCellValue(extrusion,ncells)
  cell_to_order = ConstantCellValue(order,ncells)

  (cell_to_nodes, cell_to_extrusion, cell_to_order)

end

const t = TET_AXIS
const h = HEX_AXIS
_extrussion_from_etype(::Val{1}) = (h,)
_extrussion_from_etype(::Val{2}) = (t,t)
_extrussion_from_etype(::Val{3}) = (h,h)
_extrussion_from_etype(::Val{4}) = (t,t,t)
_extrussion_from_etype(::Val{5}) = (h,h,h)
_extrussion_from_etype(v) = nothing

function  _fill_connectivity!(
    cell_to_nodes_data,
    cell_to_nodes_prts,
    o,
    etype_to_nlnodes,
    elemTypes,
    elemTags,
    nodeTags)

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
    for (i,cell) in enumerate(i_to_cell)
      a = cell_to_nodes_prts[cell-o]-1
      for lnode in 1:nlnodes
        node = i_lnode_to_node[(i-1)*nlnodes+lnode]
        cell_to_nodes_data[a+lnode] = node
      end
    end
  end

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

mshfile = "test/t1.msh"

explore(mshfile)

end # module
