# removeexpressions.jl
# This script deletes the "visit_expressions" dataset of a pixie3D HDF5 file so that it can be opened by Paraview cleanly
# Run this script from a shell via: julia removeexpressions.jl Relative/Path/To/MyPixie3dFile.h5
using HDF5
filename = pwd()*string(/)*ARGS[1]
h5open(filename,"cw") do f
	if haskey(f,"visit_expressions")
		delete_object(f,"visit_expressions")
	end
end
