module ElementTypesTest

# gmsh element types from:
# https://gitlab.onelab.info/gmsh/gmsh/-/blob/165ded1ea7075b56480411d71ea7fe0f2c757504/src/common/GmshDefines.h

using Test
using Gridap.ReferenceFEs
using GridapGmsh: gmsh
using GridapGmsh: _reffe_from_etype
using GridapGmsh: _get_lnode_to_glnode

max_order = 10
orders = 1:max_order
lin_eltypes = [ 1, 8, 26, 27, 28, 62, 63, 64, 65, 66 ]
tri_eltypes = [ 2, 9, 21, 23, 25, 42, 43, 44, 45, 46 ]
tet_eltypes = [ 4, 11, 29, 30, 31, 71, 72, 73, 74, 75 ]
qua_eltypes = [ 3, 10, 36, 37, 38, 47, 48, 49, 50, 51 ]
hex_eltypes = [ 5, 12, 92, 93, 94, 95, 96, 97, 98 ]

gmsh.initialize()

for (i,et) in enumerate(lin_eltypes)
  reffe = _reffe_from_etype(gmsh,et)
  _get_lnode_to_glnode(gmsh,et)
  order = get_order(reffe)
  p = get_polytope(reffe)
  @test p == SEGMENT
  @test order == orders[i]
end

for (i,et) in enumerate(tri_eltypes)
  reffe = _reffe_from_etype(gmsh,et)
  _get_lnode_to_glnode(gmsh,et)
  order = get_order(reffe)
  p = get_polytope(reffe)
  @test p == TRI
  @test order == orders[i]
end

for (i,et) in enumerate(tet_eltypes)
  reffe = _reffe_from_etype(gmsh,et)
  _get_lnode_to_glnode(gmsh,et)
  order = get_order(reffe)
  p = get_polytope(reffe)
  @test p == TET
  @test order == orders[i]
end

for (i,et) in enumerate(qua_eltypes)
  reffe = _reffe_from_etype(gmsh,et)
  _get_lnode_to_glnode(gmsh,et)
  order = get_order(reffe)
  p = get_polytope(reffe)
  @test p == QUAD
  @test order == orders[i]
end

for (i,et) in enumerate(hex_eltypes)
  reffe = _reffe_from_etype(gmsh,et)
  _get_lnode_to_glnode(gmsh,et)
  order = get_order(reffe)
  p = get_polytope(reffe)
  @test p == HEX
  @test order == orders[i]
end

gmsh.finalize()

end # module
