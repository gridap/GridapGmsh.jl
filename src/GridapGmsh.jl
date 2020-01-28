module GridapGmsh

using Libdl

using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: max_cells_arround_vertex
using Gridap.Geometry: _fill_cells_around_scratch!
using Gridap.Geometry: _set_intersection!

export GmshDiscreteModel


deps_jl = joinpath(@__DIR__, "..", "deps", "deps.jl")
if !isfile(deps_jl)
  s = """
  Package GridapGmsh not installed properly.
  Run Pkg.build(\"GridapGmsh\"), restart Julia and try again
  """
  error(s)
end

include(deps_jl)

include(gmsh_jl)

# Hack taken from MPI.jl
function __init__()
    @static if Sys.isunix()
        Libdl.dlopen(gmsh.lib, Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
    end
end

include("GmshDiscreteModels.jl")

end # module
