module GmshDiscreteModelsTests

using GridapGmsh
using Gridap

mshfile = "test/t1.msh"

model = GmshDiscreteModel(mshfile)

writevtk(model,"model")

end # module
