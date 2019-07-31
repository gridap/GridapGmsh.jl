module GridapGmsh

using Libdl

include(joinpath(@__DIR__, "..", "deps", "deps.jl"))

include(gmsh_jl)

# Hack taken from https://github.com/JuliaParallel/MPI.jl/blob/master/src/MPI.jl
function __init__()
    @static if Sys.isunix()
        Libdl.dlopen(gmsh.lib, Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
    end
end

end # module
