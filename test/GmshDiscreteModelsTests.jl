module GmshDiscreteModelsTests

using GridapGmsh
using GridapGmsh: gmsh
using Gridap
using Gridap.CellValuesGallery
using Gridap.DiscreteModels: DiscreteModelFromData
using StaticArrays
using UnstructuredGrids.Kernels

include("../src/GmshDiscreteModels.jl")

mshfile = "test/t1.msh"

model = GmshDiscreteModel(mshfile)

writevtk(model,"model")

end # module
