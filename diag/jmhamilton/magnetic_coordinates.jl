# magnetic_coordinates.jl
using Interpolations
using Dierckx
using MPI
using NPZ
using Statistics
using HDF5
using PyPlot
using Plots
using CairoMakie
using FFTW
#=Parallel code that calculates the magnetic perturbation projected to straight field line coordinates. It produces a series 
of arrays (b_hat_rho is the most important) with quantities defined in a psi,theta_magnetic,phi coordinate system.

RUN INSTRUCTIONS: The program is called inside an allocation as "srun -n 2 julia magnetic_coordinates.jl [user input options]"
Use mpiexec instead of srun if that slurm command is not available in your environment. Julia and MPI modules must be loaded.

Change the number of processes from 2 to some other value as long as there is at least 1 timestep per process.
Parallelization is done by assigning different timesteps to different processes. Little to no communication is required until the end when output data is written.
See below for which "user input options" can be passed at execution (order doesn't matter, format does: "option=value").

Default user options:
filepath=$(pwd)
filename=/pixie3d.h5
tstart=0
tend=nothing
q_only=false
plot_surfaces=false
plot_q=true
plot_grid=false
plot_pert=false
closed_boundary=false

Authored by Jason Hamilton
Contact: jmhamilton@lanl.gov

Latest update: 6th August 2023

=#

# Set up of parallel communicator
MPI.Init()
comm = MPI.COMM_WORLD
size_comm = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)
global error_code = 0
if rank == 0 ; start = time() ; end

# User Inputs (defaults) -------------------------------------------------------------------------------

filepath = pwd(); # Directory of .h5 file
filename = "/pixie3d.h5"; # This is the default filename from pixplot. Note the leading slash.
tstart = 0; # tstart and tend specify the range of timestep indices that will be processed by this script: [tstart+1:tend]
tend = nothing; # Will self-limit itself to the total number of time steps in the HDF5 file!
q_only = false; # Set to true to only calculate q(psi,t) and then exit the script
plot_surfaces = false; # When true, the flux surfaces will be plotted and saved, and then the script will exit
plot_q = true;  # When false, plots/animations of the q profile won't be made
plot_grid = false; # When true, the script will make a high-res mp4 of the closed surfaces and the constructed magnetic grid on top of it
plot_pert = false; # When true, 2D plots of the magnetic perturbations will be made
closed_boundary = false; # Tells the normalization routine to use the boundary value of psi instead of locating an x-point

# ------------------------------------------------------------------------------------------------------

# Read input options to override defaults
for option in ARGS
    if startswith(option,"filepath=")
        global filepath = option[length("filepath=")+1:end]
    elseif startswith(option,"filename=")
        global filename = option[length("filename=")+1:end]
    elseif startswith(option,"closed_boundary=")
        global closed_boundary = ( option[length("closed_boundary=")+1:end] == "true" )
    elseif startswith(option,"q_only=")
        global q_only = ( option[length("q_only=")+1:end] == "true" )
    elseif startswith(option,"plot_surfaces=")
        global plot_surfaces = ( option[length("plot_surfaces=")+1:end] == "true" )
    elseif startswith(option,"plot_q=")
        global plot_q = ( option[length("plot_q=")+1:end] == "true" )
    elseif startswith(option,"plot_grid=")
        global plot_grid = ( option[length("plot_grid=")+1:end] == "true" )
    elseif startswith(option,"plot_pert=")
	global plot_pert = ( option[length("plot_pert=")+1:end] == "true" )
    elseif startswith(option,"tstart=")
        global tstart = parse(Int,option[length("tstart=")+1:end])
    elseif startswith(option,"tend=")
        global tend = parse(Int,option[length("tend=")+1:end])
    else
        if rank == 0 ; println("Warning! Your option: "*option*" , was not a valid option.") end
    end
end

##################################################################   FUNCTIONS   #####################################################################

function project(arr_eint::Interpolations.Extrapolation,N::Int)
    """
    Returns the provided N-dimensional interpolated array, projected onto the new magnetic grid (psi,theta_m,phi,t).
    """
    if N == 4
    ArrP = Array{Float64,4}(undef,psi_range+1,ufdim,fidim,tdim)
    psilist = pushfirst!(collect(pn[psi_min_ind:psi_max_ind]),0.)
    for t in 1:tdim
	if Xpnts[t] != [0,0]
        for (i,p) in enumerate(psilist)
            for (j,uf) in enumerate(ufn)
                for (k,phi) in enumerate(phin)
                    if i == 1 # Adding magnetic axis point
			r = raxisarray[t]
			u = uaxisarray[t]
                        ArrP[1,j,k,t] = arr_eint(r,u,phi,t)
                    else
                        r = GD[t,round(p,digits=2),Float64(uf)][1]
                        u = GD[t,round(p,digits=2),Float64(uf)][2]
                        ArrP[i,j,k,t] = arr_eint(r,u,phi,t)
                    end
                end
            end
        end
	else
	    ArrP[:,:,:,t] .= 0. # In the case of no closed surfaces, use placeholder values
	end
    end

    elseif N == 3 # Toroidally symmetric arrays
    ArrP = Array{Float64,3}(undef,psi_range+1,ufdim,tdim)
    psilist = pushfirst!(collect(pn[psi_min_ind:psi_max_ind]),0.)
    for t in 1:tdim
        if Xpnts[t] != [0,0]
        for (i,p) in enumerate(psilist)
            for (j,uf) in enumerate(ufn)
                    if i == 1 # Adding magnetic axis point
                        r = raxisarray[t]
                        u = uaxisarray[t]
                        ArrP[1,j,t] = arr_eint(r,u,t)
                    else
                        r = GD[t,round(p,digits=2),Float64(uf)][1]
                        u = GD[t,round(p,digits=2),Float64(uf)][2]
                        ArrP[i,j,t] = arr_eint(r,u,t)
                    end
            end
        end 
        else 
            ArrP[:,:,t] .= 0. # In the case of no closed surfaces, use placeholder values
        end
    end

    else # Otherwise, assume the provided interpolant is for a Cartesian map i.e. X,Y,Z
    ArrP = Array{Float64,4}(undef,psi_range+1,ufdim,fidim,tdim)
    psilist = pushfirst!(collect(pn[psi_min_ind:psi_max_ind]),0.)
    for t in 1:tdim
        if Xpnts[t] != [0,0]
        for (i,p) in enumerate(psilist)
            for (j,uf) in enumerate(ufn)
                for (k,phi) in enumerate(phin)
                    if i == 1 # Adding magnetic axis point
                        r = raxisarray[t]
                        u = uaxisarray[t]
                        ArrP[1,j,k,t] = arr_eint(r,u,phi)
                    else
                        r = GD[t,round(p,digits=2),Float64(uf)][1]
                        u = GD[t,round(p,digits=2),Float64(uf)][2]
                        ArrP[i,j,k,t] = arr_eint(r,u,phi)
                    end
                end
            end
        end
        else
            ArrP[:,:,:,t] .= 0. # In the case of no closed surfaces, use placeholder values
        end
    end

    end

    return ArrP
end

function getTimesteps(fileid::HDF5.File)
    """
    For a specified HDF5 'fileid', returns a list of all Timestep groups and their timestep values (in Alfven times)
    """
    timesteps = filter(contains(r"Timestep"),keys(fileid)) # String array of all timestep keys (not visit_expressions) in h5 file
    time_stamps = Array{Int,1}(undef, length(timesteps))
    timesteps2 = Array{Float64,1}(undef, length(timesteps))
    for i in 1:length(timesteps)
        time_stamps[i] = parse(Int64,split(timesteps[i],"_")[2]) # Integer array of timestep values (so that they can be sorted)
        attrlist = HDF5.attributes(fileid[timesteps[i]])
        timesteps2[i] = round(read(attrlist["Time"])[1], digits = 3)
    end
    return timesteps, sort(timesteps2), sort(time_stamps)   
end

function loadh5(fileid::HDF5.File,vartype::String,myvariable::String,tfirst::Integer,tlast::Integer)
    """
    For a specified HDF5 'fileid', returns the dataset specified by strings 'vartype', 'myvariable' as an array.
    'tfirst' and 'tlast' specify the number of timesteps to be included.
    """
    mytimesteps = ["Timestep_" * string(x) for x in time_stamps[tfirst+1:tlast]] # String array of only requested timestep keys
    testarray = read(fileid["Timestep_0"][vartype][myvariable]) # Used to probe dataset dimensions
    myarray = Array{Float64,ndims(testarray)+1}(undef, (size(testarray)...,length(mytimesteps))) # An example array structure is: V(r,theta,phi,t)
    for x in mytimesteps # Return the specified variables at specified timesteps as an array
        myarray[:,:,:,findfirst(i -> i == x, mytimesteps)] = read(fileid[x][vartype][myvariable])
    end
    return myarray
end

function NormalizePsi(psitarray::Array{Float64,3},B_pol::Interpolations.Extrapolation)
    """
    Calculates the minimum value and the X-point value of the poloidal flux to normalize all toroidally-averaged psi values at a given time t.
    """
    Xpnts = Array{Any}[] # List of grid coordinates for all chosen X-points
    raxisarray = Array{Float64}(undef, size(psitarray)[3]) # Logical coordinates for the grid's minimum psi for each time step
    uaxisarray = Array{Float64}(undef, size(psitarray)[3])
    rxpntarray = Array{Float64}(undef, size(psitarray)[3]) # Logical coordinates for the chosen X-points for each time step
    uxpntarray = Array{Float64}(undef, size(psitarray)[3])
    raxisindex, uaxisindex = LocateGridAxis(psitarray) # Obtain closest grid location to the magnetic axis at all time steps
    psitnorm = Array{Float64,3}(undef,size(psitarray)...)

    Xpnt_found_t1 = true
    if rank == 0 && tstart == 0 ; Xpnt_found_t1, Xpnt_t1 = LocateGridXpoint(psitarray[:,:,1],1) end
    MPI.Barrier(comm)
    Xpnt_found_t1 = MPI.bcast(Xpnt_found_t1,0,comm) # Tell all processes if an X-point exists at the first timestep
    if Xpnt_found_t1 == false
    	if rank == 0 ; println("Using boundary values to normalize psi.") end
	global closed_boundary = true # The boundary condition is a closed surface
    end

    for t in range(1,stop=size(psitarray)[3])
	
	XPnt_found, Xpnt = LocateGridXpoint(psitarray[:,:,t],t) # Obtain closest grid location to the x-point at a given timestep       

	if XPnt_found == false && closed_boundary == false
	    println("No X-point found at t = $(rank_tstart+t). Assuming there are no closed surfaces; skipping this timestep.")
	end
	push!(Xpnts, Xpnt) # Collect the grid saddle point location for each time step
        if Xpnt != [0,0] && closed_boundary == false
	    rxpntarray[t], uxpntarray[t], psit_xpnt = LocateTrueXpoint(psitarray[:,:,t],B_pol,Xpnt,t) # Use the interpolated B_pol to get an accurate location
            raxisarray[t], uaxisarray[t], psit_min = LocateTrueAxis(psitarray[:,:,t],raxisindex[t],uaxisindex[t])
	elseif closed_boundary == true
	    rxpntarray[t], uxpntarray[t], psit_xpnt = 1., argmin(abs.(psitarray[end,:,t]))*2pi/udim, minimum(abs.(psitarray[end,:,t])) # Use boundary value
	    raxisarray[t], uaxisarray[t], psit_min = LocateTrueAxis(psitarray[:,:,t],raxisindex[t],uaxisindex[t])
	else	    
	    rxpntarray[t], uxpntarray[t], psit_xpnt = 0., 0., maximum(psitarray[:,:,t]) # These timesteps will be skipped
	    raxisarray[t], uaxisarray[t], psit_min = 0., 0., minimum(psitarray[:,:,t])
	end
        psitnorm[:,:,t] = (psitarray[:,:,t] .- psit_min)./(psit_xpnt - psit_min) # Normalize between values at the axis and values at the x-point(or boundary)
    end
    MPI.Barrier(comm)
    if rank == 0 ; println("Normalization of psi is complete.") end

    return psitnorm, Xpnts, raxisindex, uaxisindex, raxisarray, uaxisarray, rxpntarray, uxpntarray
end

function LocateGridXpoint(psitarray::Array{Float64,2},t::Int)
    """
    Returns the logical coordinate (r,theta) INDICES of the location of the separatrix's x-point ON THE GRID.
    The fastest way to find a saddle point is to loop through the grid and check for a max/min in r/theta respectively.
    """
    x0 = Array{Any}[] # Declare empty list of x-point locations
    psi0 = Vector{Float64}() # Declare empty list of psi at x-point locations
    XPnt_found = false
    for i in range(3,stop=size(psitarray)[1]-3) # Ignore locations close to the origin
        for j in range(1,stop=size(psitarray)[2]) # Find saddle point (maximum in 'r', minimum in 'theta'), ignore r=0 and r=1
            if (psitarray[i,j] >= psitarray[i+1,j] && psitarray[i,j] >= psitarray[i-1,j]) || (psitarray[i,j] >= psitarray[i+2,j] && psitarray[i,j] >= psitarray[i-2,j])
                if j == 1
                    if (psitarray[i,1] < psitarray[i,2] && psitarray[i,1] < psitarray[i,end-1]) || (psitarray[i,1] < psitarray[i,3] && psitarray[i,1] < psitarray[i,end-2])
                        push!(x0, [i,j]) ; push!(psi0, psitarray[i,j]) ; XPnt_found = true
                    end
                elseif j == size(psitarray)[2]
                    if (psitarray[i,end] < psitarray[i,2] && psitarray[i,end] < psitarray[i,end-1]) || (psitarray[i,end] < psitarray[i,3] && psitarray[i,end] < psitarray[i,end-2])
                        push!(x0, [i,j]) ; push!(psi0, psitarray[i,j]) ; XPnt_found = true
                    end
                elseif (psitarray[i,j] < psitarray[i,j+1] && psitarray[i,j] < psitarray[i,j-1]) || (psitarray[i,j] < psitarray[i,j+2] && psitarray[i,j] < psitarray[i,j-2])
                        push!(x0, [i,j]) ; push!(psi0, psitarray[i,j]) ; XPnt_found = true
                end
            end
        end
    end

    if XPnt_found == true && closed_boundary == false
	if length(x0) < 2 ; println("Missing XPnt at t = $(rank_tstart+t) . # of Xpnts found = $(length(x0))") end
	Xpnt = x0[argmin(psi0)] # Use the innermost saddle point as the normalization
    elseif closed_boundary == true
	Xpnt = [rdim,argmin(psitarray[rdim,:])] # If the boundary condition is a closed surface there is no x-point, so use the boundary value of psi
    else
	Xpnt = [0,0] # Placeholder value
    end

    return XPnt_found, Xpnt
end

function LocateGridAxis(psitarray::Array{Float64,3})
    """
    Returns the logical coordinate (r,theta) INDICES of the location of the magnetic axis ON THE GRID.
    """
    raxis = Array{Int}(undef, size(psitarray)[3])
    uaxis = Array{Int}(undef, size(psitarray)[3])
    lastindex = convert(Int,round(0.7*size(psitarray)[1])) # Highest reasonable radial index for axis (avoid points near the boundary)    
    for t in range(1,stop=size(psitarray)[3])
        axisloc = argmin(psitarray[1:lastindex,:,t])
        raxis[t] = axisloc[1]
        uaxis[t] = axisloc[2]
	if raxis[t] == 1 ; uaxis[t] = 1 end # theta is undefined at r = 0 so use theta = 0
	if uaxis[t] == size(psitarray)[2] ; uaxis[t] = 1 end # prefer theta = 0 over theta = 2pi
    end
    return raxis, uaxis
end

function LocateTrueAxis(psitarray::Array{Float64,2},raxisindex::Integer,uaxisindex::Integer)
    """
    Interpolates over psi_t and finds the minimum at a specific time step. More accurate than finding the minimum on the grid.
    Returns the logical coordinates (r,theta) of the axis as well as the poloidal flux at the axis.
    """
    r_range = LinRange(rn[convert(Int,max(raxisindex-2,1))], rn[convert(Int,min(raxisindex+2,rdim))], 100) # Construct a sample of points located around the grid minimum
    u_range = LinRange(un[convert(Int,max(uaxisindex-2,1))], un[convert(Int,min(uaxisindex+2,udim))], 100)
    psi_eint = extrapolate(scale(Interpolations.interpolate(psitarray,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))))),rn,un), (Line(),Periodic()))
    psi_neighbourhood = psi_eint(r_range,u_range) # Values of psi in the neighbourhood of the grid minimum
    trueaxis_r = r_range[argmin(psi_neighbourhood)[1]] # Location of the axis
    trueaxis_u = u_range[argmin(psi_neighbourhood)[2]]

    return trueaxis_r, trueaxis_u, psi_eint(trueaxis_r,trueaxis_u)
end

function LocateTrueXpoint(psitarray::Array{Float64,2},B_pol::Interpolations.Extrapolation,Xpnt_ind::Vector{Any},t::Integer)
    """
    Interpolates over the poloidal magnetic field and finds the minimum at a specific time step. More accurate than finding the minimum on the grid.
    Returns the logical coordinates (r,theta) of the X-point as well as the poloidal flux at this point.
    """
    r_range = LinRange(rn[convert(Int,max(Xpnt_ind[1]-2,1))], rn[convert(Int,min(Xpnt_ind[1]+2,rdim))], 100) # Construct a sample of points located around the grid minimum
    u_range = LinRange(un[convert(Int,max(Xpnt_ind[2]-2,1))], un[convert(Int,min(Xpnt_ind[2]+2,udim))], 100)
    psi_eint = extrapolate(scale(Interpolations.interpolate(psitarray,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))))),rn,un), (Line(),Periodic()))
    B_pol_neighbourhood = B_pol(r_range,u_range,t) # Values of B_pol in the neighbourhood of the grid minimum
    truexpnt_r = r_range[argmin(B_pol_neighbourhood)[1]] # Location of the X-point
    truexpnt_u = u_range[argmin(B_pol_neighbourhood)[2]]

    return truexpnt_r, truexpnt_u, psi_eint(truexpnt_r,truexpnt_u)
end

function ConvertCell2Node(cellarray::Array{Float64,4})
    """
    Converts cell-centered quantities to node-based.
    Psi needs to be node-based, but is stored as cell-centered if car_diag_plots = f is set in the pixie3d input deck.
    """
    var_eint_cell = extrapolate(scale(Interpolations.interpolate(cellarray,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),NoInterp())), rc,uc,phic,tn), (Line(),Periodic(),Periodic(),Interpolations.Throw()))
    nodearray = var_eint_cell(rn,un,phin,tn)

    return nodearray
end

function NodeInterpolation(nodearray::Array{Float64})
    """
    Returns an interpolation object of the provided array over a node-based (r,theta,phi) or (r,theta) grid.
    """
    if ndims(nodearray) == 3
	node_eint = extrapolate(scale(Interpolations.interpolate(nodearray,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),NoInterp())),rn,un,tn), (Line(),Periodic(),Interpolations.Throw()));
    elseif ndims(nodearray) == 4
	node_eint = extrapolate(scale(Interpolations.interpolate(nodearray,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),NoInterp())),rn,un,phin,tn), (Line(),Periodic(),Periodic(),Interpolations.Throw()));
    else ; node_eint = nothing ; global error_code = 3
    end
    ErrorCheck()

    return node_eint
end

function PsiContour(psitarray::Array{Float64,3})
    """
    Computes contours of psi at specified isosurfaces (1%,2%,etc) and then returns the coordinates of the contour segments, for each time step.
    Example: [r,theta] of kth vertex on pth isosurface at time t is given by surfaces[p,t][:,k]
    Example: list of all radius values on pth isosurface at time t is given by surfaces[p,t][1,:] (theta values are at index 2)
    """
    surfaces = Array{Any,2}(undef,psi_range,tdim)
    psilevels = 0.01 .* [psi_min_ind-1:psi_max_ind-1;] # Create a vector of the desired psi isosurfaces
    rarray = LinRange(0.,1.,512) # Use higher resolution logical coordinate arrays to get smoother contours
    uarray = LinRange(0.,2pi,512)
    for t in range(1, stop=tdim)
	if Xpnts[t] != [0,0] # Don't attempt to contour isosurfaces at timesteps with no separatrix
	    psit_inside = psitarray[:,:,t] # To avoid contouring psi outside the separatrix, limit the array to a radius inside psi = 1
            for i in range(1, stop=udim)
                seploc = findfirst(x -> x >= 0.99, psit_inside[:,i]) # Determine the location of the separatrix
	        if seploc == 1 # If the separatrix does not enclose the origin, locate the range of r's that have psi < 1
	   	    rmin = findfirst(x -> x < 0.99, psit_inside[:,i])
		    if rmin != nothing
		        seploc = findfirst(x -> x >= 0.99, psit_inside[rmin:end,i])
		        if seploc != nothing ; seploc += rmin - 1 end # Use the r's inside the region of closed flux surfaces
		        if seploc == nothing ; seploc = 1 end # In this case, the regions with psi < 1 are open surfaces, so don't use them
		        psit_inside[1:rmin-1,i] .= 1.0
		    end
	        end
                if seploc == nothing ; seploc = rdim + 1 end # If there is no separatrix, use the full radius range
                psit_inside[seploc:rdim,i] .= 1.0 # Set these unwanted locations to some arbitrary value outside the range (0,1)
            end

   	    # psi_interp should be transposed before contours are taken for correct r,theta indexing
	    psi_interp = extrapolate(scale(Interpolations.interpolate(transpose(psit_inside[:,:]),(BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid()))))),un,rn), (Periodic(),Line()));

	    plt.clf()
	    CS = plt.contour(rarray,uarray,psi_interp(uarray,rarray),levels=psilevels)
            # Grab the contour segments at each isosurface, CS.allsegs is an array that contains the coordinates of all vertices of each segment
	    #if plot_surfaces ; plt.savefig("logicalsurface$(lpad(rank_tstart+t,4,"0")).png") end

            for p in range(1, stop=psi_range)
                tmp = []
	        for i in range(1, length(CS.allsegs[p,:]))
		    tmp = cat(tmp, CS.allsegs[p,i], dims=1) # Concatenate all segments into a single continuous array
	        end
                surfaces[p,t] = transpose(tmp) # List of all logical (r,theta) values for all segments
                if size(surfaces[p,t])[1] == 1 && length(surfaces[p,t]) != 0 # In the case of multiple segments for a given contour, CS.allsegs requires different indexing
		    tmp = surfaces[p,t][1]
		    for j in range(2, length(surfaces[p,t]))
		        tmp = cat(tmp,surfaces[p,t][j],dims=2)
		    end
		    surfaces[p,t] = tmp
	        end
		if length(surfaces[p,t]) == 0 ; surfaces[p,t] = [0.,0.] end # If the closed flux surface is very small, it might not be resolved
            end
	else
	    for p in range(1, stop=psi_range)
		surfaces[p,t] = [0.,0.] # Timesteps that do not have a separatrix will have no closed flux surfaces
	    end
	end
    end

    return surfaces
end

function PlotSurfaces(surfaces::Array{Any,2})
    """
    Plots the flux surfaces and saves them to the current directory.
    """
    psi_plotted = [5,20,40,60,80,95,psi_max_ind-1] .- (psi_min_ind - 2) # The list of psi isosurfaces to be plotted (in addition to psi_min_ind - 1)
    minY = minimum(H)
    maxY = maximum(H)
    avgX = 0.5*(maximum(R) + minimum(R))
    minX = avgX - 0.5*(maxY-minY)
    maxX = avgX + 0.5*(maxY-minY)
    ENV["GKSwstype"]="nul"; # This tells the plotting backend to not display plot objects
    for t in range(1, stop=tdim)
        fsplot = nothing # Clear plot
        xlist = R_eint.(surfaces[1,t][1,:],surfaces[1,t][2,:])
        ylist = H_eint.(surfaces[1,t][1,:],surfaces[1,t][2,:])
        fsplot = Plots.plot(xlist,ylist, 
#		annotate = (xlist[rand(1:length(xlist))],ylist[rand(1:length(ylist))],Plots.text("ψ = $(round(0.01*(psi_min_ind-1),digits=2))",pointsize=8)), 
		title = " t = $(rank_tstart+t) ", 
		xlims=(minX,maxX), ylims=(minY,maxY),
		dpi=600, aspect_ratio = :equal,
		xlabel="Radius (m)" , ylabel="Height (m)", linewidth=2, xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22, legend=false)
        for p in psi_plotted # Add all flux surfaces at time t to the same plot
            xlist = R_eint.(surfaces[p,t][1,:],surfaces[p,t][2,:])
            ylist = H_eint.(surfaces[p,t][1,:],surfaces[p,t][2,:])
            Plots.plot!(fsplot, xlist, ylist, 
#		annotate = (xlist[rand(1:length(xlist))],ylist[rand(1:length(ylist))],Plots.text("ψ = $(round(0.01*(p+psi_min_ind-2),digits=2))",pointsize=8)),
		linewidth=2)
        end
        Plots.savefig(fsplot,"surface$(lpad(rank_tstart+t-tstart,4,"0")).png")
    end

    for t in range(1, stop=tdim) # In addition to the Cartesian plot, also plot the surfaces in logical coordinates
	fsplot = nothing # Clear plot
	fsplot = Plots.plot(surfaces[1,t][1,:],surfaces[1,t][2,:],
		title = " t = $(rank_tstart+t) ",
		xlims = (0,1), ylims = (0,2pi),
		dpi=600,
		xlabel="r", ylabel="theta", linewidth=2, xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22, legend=false)
	for p in 2:psi_range
		Plots.plot!(fsplot,surfaces[p,t][1,:],surfaces[p,t][2,:], linewidth=2)
	end
	Plots.savefig(fsplot,"logicalsurface$(lpad(rank_tstart+t-tstart,4,"0")).png")
    end
        
    return nothing
end

function SurfaceAnimation()
    """
    Turns all surface plots in the current directory into an avi animation. Deletes the individual plots.
    """
    run(`/usr/bin/bash -c "ffmpeg -hide_banner -loglevel error -y -r 20 -f image2 -i surface%04d.png -vcodec mjpeg -q:v 2 -pix_fmt yuvj420p surfaces.mp4"`)
    run(`/usr/bin/bash -c "rm surface*.png"`)
    run(`/usr/bin/bash -c "ffmpeg -hide_banner -loglevel error -y -r 20 -f image2 -i logicalsurface%04d.png -vcodec mjpeg -q:v 2 -pix_fmt yuvj420p logicalsurfaces.mp4"`)
    run(`/usr/bin/bash -c "rm logicalsurface*.png"`)

    return nothing
end

function GridAnimation()
    """
    Similar to the above commands, but for the case with the magnetic grid points plotted as well.
    """
    run(`/usr/bin/bash -c "ffmpeg -hide_banner -loglevel error -y -r 20 -f image2 -i grid%04d.png -vcodec mjpeg -q:v 2 -pix_fmt yuvj420p grid.mp4"`)
    run(`/usr/bin/bash -c "rm grid*.png"`)

    println("Grid animation is completed.")

    return nothing
end

function q_profile(surfaces::Array{Any,2})
    """ 
    Calculates the q profile and magnetic angles by integrating around the pre-determined flux surfaces.
    """ 
    q = Array{Float64,2}(undef,psi_range,tdim)
    theta_m = Array{Any,2}(undef,psi_range,tdim)
    pathlength = Array{Any,2}(undef,psi_range,tdim)
    for t in 1:tdim
	if Xpnts[t] != [0,0]
	    for p in 1:psi_range
	        path = convert(Array{Float64,2},surfaces[p,t][:,:])
	        q[p,t], theta_m[p,t], pathlength[p,t] = flux_integrator(path,t)
	    end
        else
	    q[:,t] .= 0. # The equation for q used here is ill-posed if there are no closed flux surfaces, q=0 is a placeholder value
	    theta_m[:,t] .= 0.
	    pathlength[:,t] .= 0.
	end
    end
    MPI.Barrier(comm)
    if rank == 0 ; println("The q profiles are finished!") end

    return q, theta_m, pathlength
end

function flux_integrator(path::Array{Float64,2},t::Integer)
    """
    Uses a simple midpoint rule to numerically evaluate the integral of B_tor / B_pol around a path defined by an isosurface of psi.
    """

    integrand = B3t_eint.(path[1,:],path[2,:],t)./JBpN_eint.(path[1,:],path[2,:],t) # Ratio of B_tor / B_pol gives q after a full circuit around the path
    bad_indices = Vector{Int}()
    for i in 1:length(path[2,:])
        if Jac_eint(path[1,i],path[2,i],t) > 1e-4
            integrand[i] = integrand[i]/sqrt(Jac_eint(path[1,i],path[2,i],t)) # Geometric factor (avoid dividing by zero)
	else
            push!(bad_indices,i) # Mark the locations where the path crosses r=0, where the Jacobian is close to zero
        end
    end

    for i in bad_indices # For these locations, average neighbouring values on the path instead of dividing by a small Jacobian factor (avoid a large error)
        integrand[i] = 0.5*(circshift(integrand,1)[i] + circshift(integrand,-1)[i])
    end

    # Midpoint rule - Cartesian
    rlist = R_eint.(path[1,:],path[2,:]) # List of all cylindrical R values along the path
    hlist = H_eint.(path[1,:],path[2,:]) # List of all Z/height values along the path
    dl = sqrt.((rlist .- circshift(rlist,-1)).^2 .+ (hlist .- circshift(hlist,-1)).^2) # Length of line element connecting each vertex along the path
    dl_mp = 0.5*(dl + circshift(dl,1))
    q_m = sum(dl_mp.*integrand)/(2pi)

    # Evaluate the magnetic angle theta_m at each point along the path, along with the total path length up to each point
    theta_m = Array{Float64}(undef, length(path[2,:]))
    pathlength = Array{Float64}(undef, length(path[2,:]))
    for i in 1:length(path[2,:])
	theta_m[i] = sum(dl_mp[1:i].*integrand[1:i])/q_m
	if abs(theta_m[i]) < 1e-4 ; theta_m[i] = 0. + i*1e-7 end # Small negative angles are possible at r=0. This ensures positive monotonicity
	if abs(theta_m[i] - 2pi) < 1e-4 ; theta_m[i] = 2pi - 1e-4 + i*1e-7 end # This again ensures monotonicity at theta_m = 2pi
	pathlength[i] = sum(dl_mp[1:i])
    end

    return q_m, theta_m, pathlength
end

function stretch_q(q_array::Array{Float64,2})
    """
    Linearly extrapolates q to the magnetic axis and the separatrix.
    """ 
    psin_list = LinRange(0.0,1.0,psidim)
    small_grid = LinRange(psin_list[psi_min_ind],psin_list[psi_max_ind],psi_range)
    q_str = Array{Float64,2}(undef,psidim,tdim)
    for t in 1:size(q_array)[2] 
        q = Float64.(q_array[:,t])
        q_eint = Dierckx.Spline1D(small_grid,q; k=1, bc="extrapolate") # The polynomial order of this extrapolation impacts q(0), linear works best
        q_str[:,t] = q_eint(psin_list)
    end                     
                        
    return q_str        
end 

function Plotq(q_array::Array{Float64,2})
    """
    Plot the q(psi) profile over all time steps. Label the integer values of q.
    Also plot q(0) over time.
    """
    ENV["GKSwstype"]="nul"; # This tells the plotting backend to not display plot objects
    pn = LinRange(0.,1.,psidim)

    anim_q = @animate for t in 1:size(q_array)[2]
	Plots.plot(pn,q_array[:,t], title ="q(ψ) profile", linewidth=3, xlims = (0,1), xlabel="Poloidal Flux ψ/ψsep", ylabel="q", xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22, legend=false);
	Plots.annotate!(0.25,maximum(q_array)*0.9,"t = $(timesteps2[tstart+t]) τA ") # Label the timestep value        

	q_psi_eint = extrapolate(scale(Interpolations.interpolate(q_array[:,t],(BSpline(Cubic(Line(OnGrid()))))),pn), Line()) # Interpolate q(psi)
	max_q_integer = min(convert(Int,floor(q_array[end,t])),10) # Plot intersections up to q=10
	min_q_integer = min(convert(Int,ceil(q_array[begin,t])),11)
        if min_q_integer != 0
	    for i in min_q_integer:max_q_integer # These are the values of q that will be labelled with vertical lines
                p = findfirst(x -> x >= i,q_array[:,t]) # Intercept (q(p0) = i) is in between p-1 and p
		qi = q_array[p,t]
		pu = pn[p] # Upper limit
		pl = pn[p-1] # Lower limit
		p0 = pu
		count = 0
		while abs(qi - i) > 1e-3 && count < 100 # This is the error tolerance. Limit loop to 100 iterations.
			p0 = 0.5*(pu + pl) # Next guess
			qi = q_psi_eint(p0)
			if qi > i
				pu = p0
			else
				pl = p0
			end
			count += 1
		end
        	Plots.vline!([p0]) # Add a vertical line where q(p0) = i
	    end
	end
    end

    Plots.mp4(anim_q,"q_profile.mp4",fps=20,loop=1) # Set the framerate equal to the setting used in diagnosticplots.jl , to be consistent with those movies

    AxisPlot = Plots.plot(timesteps2[tstart+1:tend],q_array[1,:], title="q(ψ=0,t) time series", linewidth=3, xlabel="Alfven times", ylabel="q", xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22, formatter = :plain);

    Plots.savefig(AxisPlot,"q_axis.png")

    return nothing
end

function new_grid(surfaces::Array{Any,2},theta_m::Array{Any,2},pathlength::Array{Any,2})
    """
    Creates a dictionary that translates the new grid (t,psi,theta_mag) -> (r,theta) for projection of arrays.
    """
    local l_of_uf_rsc
    Grid_Dict = Dict{Tuple{Int,Float64,Float64},Tuple{Float64,Float64}}() # Returns (r,theta) for each magnetic grid point
    for t in 1:tdim
        if Xpnts[t] != [0,0]
            for p in 1:psi_range
	    if length(theta_m[p,t]) > 1 # Otherwise, the flux surface is singular/underresolved, and should be skipped
                path = surfaces[p,t][:,:] # Set of all 'path points' used to contour the flux surface
                uf_rsc = theta_m[p,t][:] # Values of theta_m evaluated at every path point
                l_path = pathlength[p,t][:] # Interpolate pathlength over theta_m 
                l_path_tmp = copy(l_path)
                try
		    l_of_uf_rsc = Spline1D(pushfirst!(uf_rsc[1:end],0.),pushfirst!(l_path_tmp[1:end],0.),k=1,bc="extrapolate")
	        catch e
		    println("Non-monotonic theta_m values at t=$(rank_tstart+t), psi=$(p+psi_min_ind-1). Likely r=0 issues. Exiting.")
		    global error_code = 6 ; @goto errorCleanUp
		end
		l_grid = l_of_uf_rsc(ufn) # Evaluate pathlength(theta_m) at the desired 'grid points' of theta_m between 0,2pi
                for j in 2:ufdim-1
                    k = findfirst(x -> x >= l_grid[j], l_path) # Find nearest path points where each desired grid point occurs
		    if k == nothing ; println("t=$(rank_tstart+t), p=$(p), j=$(j)"); break; end
		    r_kp = 0.5*(path[1,k] + circshift(path[1,:],-1)[k]) # Coordinates of path point ahead of grid point
                    th_kp = AverageTheta(path[2,k],circshift(path[2,:],-1)[k])
                    r_km = 0.5*(circshift(path[1,:],1)[k] + path[1,k]) # Coordinates of path point behind grid point
                    th_km = AverageTheta(path[2,k],circshift(path[2,:],1)[k])
                    dl_k = (l_grid[j] - circshift(l_path,1)[k]) / (l_path[k] - circshift(l_path,1)[k])  # Path length difference between grid point and path point
                    r_grid = r_km + dl_k*(r_kp - r_km) # A weighted average of the nearby path points gives the location of the grid point
                    th_grid = th_km + dl_k*(th_kp - th_km)
                    # Equivalance between each magnetic grid point and the logical grid coordinates
                    Grid_Dict[(Int64(round(t)),round(pn[p+psi_min_ind-1],digits=2),ufn[j])] = (r_grid,th_grid)
                end
                Grid_Dict[(Int64(round(t)),round(pn[p+psi_min_ind-1],digits=2),ufn[1])] = (path[1,1],path[2,1]) # First point is always theta_m = 0.
                Grid_Dict[(Int64(round(t)),round(pn[p+psi_min_ind-1],digits=2),ufn[end])] = (path[1,end],path[2,end]) # Last point is always theta_m = 2π
	    else
                for j in 1:ufdim
		    Grid_Dict[(Int64(round(t)),round(pn[p+psi_min_ind-1],digits=2),ufn[j])] = (0.,0.) # Placeholder values for singular/underresolved surfaces
                end
            end
    	    end
        else
            for p in 1:psi_range
                for j in 1:ufdim
                    Grid_Dict[(Int64(round(t)),round(pn[p+psi_min_ind-1],digits=2),ufn[j])] = (0.,0.) # Placeholder values for non-processed timesteps
                end
            end
        end
    end
    @label errorCleanUp # This label simply allows the script to cleanly leave the nested loop above without disrupting the function

    return Grid_Dict
end

function AverageTheta(theta1::Real,theta2::Real)
    """
    Takes two angles (radians) and finds their average. Avoids issues of 2π discontinuities.
    """
    if abs(theta1 - theta2) >= pi
        theta1 = mod(theta1,2pi)
        theta2 = mod(theta2,2pi)
        thetamin = min(theta1,theta2)
        thetamax = max(theta1,theta2)
        if thetamax - thetamin >= pi ; thetamax -= 2pi end
        Avg = 0.5*(thetamin + thetamax)

    else ; Avg = 0.5*(theta1 + theta2) end

    return Avg
end

function monotonic(arr)
    """
    Checks if an array is monotonically increasing. Used for de-bugging.
    """
    check = true
    for (i,f) in enumerate(arr[2:end])
        if f <= arr[i] ; check = false end
    end
    return check
end

function PlotGrid(GD::Dict,surfaces::Array{Any,2})
    """
    Plots the magnetic grid points along with the flux surfaces.
    """
    ENV["GKSwstype"]="nul"; # This tells the plotting backend to not display plot objects
    for t in 1:tdim
        gplot = nothing # Clear plot
        gplot = Plots.plot(surfaces[1,t][1,:],surfaces[1,t][2,:],
                title = " t = $(rank_tstart+t) ",
                xlims = (0,1), ylims = (0,2pi),
                dpi=600,
                xlabel="r", ylabel="theta", linewidth=1, xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22, legend=false)
	
	for p in 1:psi_range
	    rpoints = Array{Any}(undef,ufdim)
	    thpoints = Array{Any}(undef,ufdim)
	    for j in 1:ufdim
		rpoints[j] = GD[t,round(pn[p+psi_min_ind-1],digits=2),ufn[j]][1] 
		thpoints[j] = GD[t,round(pn[p+psi_min_ind-1],digits=2),ufn[j]][2]
	    end
	    if p > 1 ; Plots.plot!(gplot,surfaces[p,t][1,:],surfaces[p,t][2,:], linewidth=1) end
	    Plots.scatter!(gplot,rpoints,thpoints, mc=:black, ms=3, shape=:x)
	end
	Plots.savefig(gplot,"grid$(lpad(rank_tstart+t-tstart,4,"0")).png")
    end

    return nothing
end

function PlotdB(dBarray::Array{Float64,4},component::Int,toroidalAngleIndex::Int)
    """
    Uses CairoMakie to create 2D heat plots of the provided magnetic perturbation array.
    """
    ENV["GKSwstype"]="nul"; # This tells the plotting backend to not display plot objects

    if component == 1 ; mytitle = "Radial bψ perturbation" ; myfile = "psipert.mp4"
    elseif component == 3 ; mytitle = "Toroidal bϕ perturbation" ; myfile = "phipert.mp4"
    else ; mytitle = "Unknown perturbation" ; myfile = "pert.mp4" ; end

    tA = 1.8444743E-04 # Alfven time in ms
    tA = tA*sqrt(2) # For Deuterium
    tA = tA*sqrt(300.) # For the 300x density simulation
    timesteps_ms = round.(timesteps2.*tA, digits=2)

    k = min(toroidalAngleIndex,size(dBarray)[3]) # This determine at what toroidal angle to take the 2D slice
    k = max(1,k)

    # We first generate the first timestep alone
    plot2d, myaxis, mydata = CairoMakie.surface(Xproj[:,:,k,1], Zproj[:,:,k,1], dBarray[:,:,k,1],
				figure=(; resolution=(900,1200), fontsize=30),
				axis=(type=Axis3, aspect= :data, azimuth = 1.5pi, elevation = 0.5pi,
					zlabelvisible=false, zticksvisible=false, zticklabelsvisible=false,
					title=mytitle, titlesize=30,
					ylabel="Z", xlabel="R", ylabelrotation=0.5*pi,
					xlabelvisible=false, ylabelvisible=false,
					xticklabelsvisible=false, yticklabelsvisible=false,
					xticksvisible=false, yticksvisible=false),
				interpolation=true,
				colormap= :gist_rainbow,
				colorrange=(0.,maximum(dBarray)),
				shading=false)
    mytext = CairoMakie.text!(myaxis,0.6,0.9,0.; text="t = " * string(timesteps2[tstart+1]) * " τA = "*string(timesteps_ms[tstart+1])*" ms",
			      fontsize = 30, space = :relative)
    mytext.fontsize = 30

    CairoMakie.Colorbar(plot2d[1,2], mydata, height=Relative(1.), width = 40, ticklabelsize=24) # Add the colour bar to the right
 
    # Now we start to iterate over all timesteps of all data sent to the root process from all other processes
    myiterator = range(1,time_interval,step=1)

    CairoMakie.record(plot2d, myfile, myiterator; framerate = 20) do t
	mydata.x = Xproj[:,:,k,t]
	mydata.y = Zproj[:,:,k,t]
        mydata.z = dBarray[:,:,k,t]
	mytext.text = "t = " * string(timesteps2[tstart+t]) * " τA = "*string(timesteps_ms[tstart+t])*" ms"
    end

    return nothing
end

function ElapsedTime(start::Float64)
    """
    Calculates the elapsed runtime of the script in the appropriate units.
    """
    elapsed = time() - start
    minutes = elapsed/60.
    hours = minutes/60.
    days = hours/24.

    if days >= 1. ; units = "days" ; runtime = round(days,digits=2)
    elseif hours >= 1. ; units = "hours" ; runtime = round(hours,digits=2)
    elseif minutes >= 1. ; units = "minutes" ; runtime = round(minutes,digits=2)
    else ; units = "seconds" ; runtime = round(elapsed,digits=2)
    end

    return runtime, units
end

function ErrorCheck()
    """
    Performs error handling and cleanly exits all processes.
    """
    for i in 0:size_comm-1 # Search all processes for potential non-zero error code
	share_error = MPI.bcast(error_code,i,comm)
	if share_error != 0 ; global error_code = share_error end # Send non-zero error_code to all processes
    end
    MPI.Barrier(comm)
    if error_code != 0
        if rank == 0 && error_code == 1 ; println("\nUser Error! You've requested more processes than time steps! Exiting.")
        elseif rank == 0 && error_code == 2 ; println("\nError! Number of selected saddle points doesn't match the number of timesteps! Exiting.")
        elseif rank == 0 && error_code == 3 ; println("\nError! Arrays must either be 2D toroidally averaged or 3D. Bad data. Exiting.")
	elseif rank == 0 && error_code == 4 ; println("\nError! Magnetic axis is located on the boundary/past the separatrix! Bad data! Exiting.")
	elseif rank == 0 && error_code == 5 ; println("\nError! Dimensions of q do not match ($(psidim), $(tend - tstart)). Exiting.")
	elseif rank == 0 && error_code == 6 ; println("\nError! Failed to construct magnetic grid.")
	elseif rank == 0 ; println("\nUnknown error code. Exiting.") end
        MPI.Finalize()
        exit(error_code)
    else
        return nothing
    end
end

##################################################################   EXECUTION   #####################################################################

# Load pixie3d.h5
pixiefile = h5open(filepath * filename,"r"); # Load the .h5 file in Julia (file structure is retained)
timesteps, timesteps2, time_stamps = getTimesteps(pixiefile); # Load the timestep values

# Error handling
if tend != nothing
    if tend > length(timesteps)
        if rank == 0 ; println("\nWarning! Your specified 'tend' exceeds the final timestep in the HDF5 file! Reducing 'tend' to this final timestep.") end
        tend = length(timesteps);
    end
else ; tend = length(timesteps); end

if (tstart+1) > tend
    if rank == 0 ; println("\nWarning! Your specified 'tstart' is equal to or exceeds 'tend'! Reducing 'tstart' to 'tend - 1' .") end
    tstart = tend - 1;
end

tstart = max(0,tstart);
tend = max(tstart+1,tend);
if size_comm > (tend-tstart) ; global error_code = 1; end
ErrorCheck()

# Print set-up variables
if rank == 0
    println("Runnning magnetic_coordinates.jl with "*string(size_comm)*" parallel process(es).\n")
    println("Processing the HDF5 file located at: "*filepath*filename) 
    println("Timesteps numbered between "*string(tstart+1)*" and "*string(tend)*" will be processed.")
    if plot_surfaces == true && plot_q == false && plot_grid == false ; println("Flux surfaces will be plotted.\n") end
    if plot_surfaces == true && plot_q == true && plot_grid == false ; println("Flux surfaces and the q profile will be plotted.\n") end
    if plot_surfaces == true && plot_q == false && plot_grid == true ; println("Flux surfaces and the magnetic grid will be plotted.\n") end
    if plot_surfaces == true && plot_q == true && plot_grid == true ; println("Flux surfaces, q profile, and the magnetic grid will be plotted.\n") end
    if plot_surfaces == false && plot_q == true && plot_grid == false ; println("The q profile will be plotted.\n") end
    if plot_surfaces == false && plot_q == true && plot_grid == true ; println("The q profile and the magnetic grid will be plotted.\n") end
    if plot_surfaces == false && plot_q == false && plot_grid == true ; println("The magnetic grid will be plotted.\n") end
    if plot_surfaces == false && q_only == false && plot_grid == false ; println("The magnetic perturbations will be computed.\n") end
end

# Construct parallel partitions
time_interval = max(1,tend-tstart); 
chunk = max(1.,floor(time_interval/size_comm)); # Loading chunk: each process loads different times
rank_tstart = convert(Int,tstart+rank*chunk); # Starting and ending points of each process
rank_tend = convert(Int,tstart+(rank+1)*chunk);
if rank == (size_comm-1) ; rank_tend = tend; end
final_chunk_size = convert(Int,max(1.,time_interval-(size_comm-1)*chunk)); # The last partition may have a different number of timesteps

# Load variable data
psi = loadh5(pixiefile,"Diagnostics","Poloidal flux",rank_tstart,rank_tend); # Poloidal flux
B1 = loadh5(pixiefile,"Cnv_variables","B^1",rank_tstart,rank_tend); # Contravariant components
B2 = loadh5(pixiefile,"Cnv_variables","B^2",rank_tstart,rank_tend);
B3 = loadh5(pixiefile,"Cnv_variables","B^3",rank_tstart,rank_tend);
B_1 = loadh5(pixiefile,"Cov_variables","B_1",rank_tstart,rank_tend); # Covariant components
B_2 = loadh5(pixiefile,"Cov_variables","B_2",rank_tstart,rank_tend);
B_3 = loadh5(pixiefile,"Cov_variables","B_3",rank_tstart,rank_tend);
Bx = loadh5(pixiefile,"Car_variables","Bx",rank_tstart,rank_tend); # Cartesian vector components
By = loadh5(pixiefile,"Car_variables","By",rank_tstart,rank_tend);
Bz = loadh5(pixiefile,"Car_variables","Bz",rank_tstart,rank_tend);

# Load coordinate data (use any timestep to do this)
X = read(pixiefile[timesteps[1]]["nodes"]["X"]);
Y = read(pixiefile[timesteps[1]]["nodes"]["Y"]);
Z = read(pixiefile[timesteps[1]]["nodes"]["Z"]);

# Close the .h5 file
close(pixiefile)

# Ensure psi is monotonically increasing with radius
psi_sign = 1.
if psi[begin,1,1,1] > psi[end,1,1,1]
    psi .= -psi;
    psi_sign = -1.; # This keeps track of if this reversal of psi took place
    if rank == 0 ; println("\nSign of psi was flipped to ensure it is monotonically increasing.\n") end
end

# Special handling for the 2D case only: extrude to a second toroidal point to enable proper extrapolation boundary conditions
if size(B3)[3] == 1
    if size(psi)[3] == 2 ; psi = cat(dims=3,psi,psi)[:,:,1:3,:]; end; # If node-based, extrude psi to a 3rd grid point at phi=pi
    if size(psi)[3] == 1 ; psi = cat(dims=3,psi,psi); end; # If cell-based, extrude psi like the other variables
    B1 = cat(dims=3,B1,B1);
    B2 = cat(dims=3,B2,B2);
    B3 = cat(dims=3,B3,B3);
    B_1 = cat(dims=3,B_1,B_1);
    B_2 = cat(dims=3,B_2,B_2);
    B_3 = cat(dims=3,B_3,B_3);
    Bx = cat(dims=3,Bx,Bx)[:,:,1:3,:];
    By = cat(dims=3,By,By)[:,:,1:3,:];
    Bz = cat(dims=3,Bz,Bz)[:,:,1:3,:];
    X = cat(dims=3,X,X)[:,:,1:3];
    Y = cat(dims=3,Y,Y)[:,:,1:3];
    Z = cat(dims=3,Z,Z)[:,:,1:3];
end

# Definitions of cell grid
num_r_cells = size(B3)[1];
num_u_cells = size(B3)[2];
num_phi_cells = size(B3)[3];
dn_r = (1.0/num_r_cells); # size of each cell
dn_u = ((2.0*pi)/num_u_cells);
dn_p = ((2.0*pi)/num_phi_cells);

# Cell-based grid (cell-center locations)
rc = LinRange(0.0+(dn_r/2.0),1.0-(dn_r/2.0),num_r_cells);
uc = LinRange(0.0+(dn_u/2.0),2.0*pi-(dn_u/2.0),num_u_cells);
phic = LinRange(0.0+(dn_p/2.0),2.0*pi-(dn_p/2.0),num_phi_cells);

# Node-based grid
rn = LinRange(0.0,1.0,(num_r_cells+1));
un = LinRange(0.0,2.0*pi,(num_u_cells+1));
phin = LinRange(0.0,2.0*pi,(num_phi_cells+1));
tn = range(1,stop=size(B3)[4]);

if size(psi)[1] == num_r_cells # Check if psi is cell-centered, and convert to node-based
    if rank == 0; println("Psi is cell-centered... converting to node-based.") end
    psi = ConvertCell2Node(psi);
end

# Node-based grid dimensions
rdim = size(psi)[1];
udim = size(psi)[2];
fidim = size(psi)[3];
tdim = size(psi)[4];

# Convert magnetic field from Cell to Grid
B1 = ConvertCell2Node(B1);
B2 = ConvertCell2Node(B2);
B3 = ConvertCell2Node(B3);
B_1 = ConvertCell2Node(B_1);
B_2 = ConvertCell2Node(B_2);
B_3 = ConvertCell2Node(B_3);

# Ensure r=0 regularity by setting this first node point as the average over all theta values
for t in range(1,stop=tdim)
    for k in range(1,stop=fidim)
        B1[1,:,k,t] .= mean(B1[1,:,k,t]);
        B2[1,:,k,t] .= mean(B2[1,:,k,t]);
        B3[1,:,k,t] .= mean(B3[1,:,k,t]);
        B_1[1,:,k,t] .= mean(B_1[1,:,k,t]);
        B_2[1,:,k,t] .= mean(B_2[1,:,k,t]);
	B_3[1,:,k,t] .= mean(B_3[1,:,k,t]);
	Bx[1,:,k,t] .= mean(Bx[1,:,k,t]);
	By[1,:,k,t] .= mean(By[1,:,k,t]);
	Bz[1,:,k,t] .= mean(Bz[1,:,k,t]);
	psi[1,:,k,t] .= mean(psi[1,:,k,t]);
    end
end
for k in range(1,stop=fidim)
    X[1,:,k] .= mean(X[1,:,k]);
    Y[1,:,k] .= mean(Y[1,:,k]);
    Z[1,:,k] .= mean(Z[1,:,k]);
end

# Calculate the determinant of the Jacobian, avoid dividing by zero
Jac = @. (B1*B_1 + B2*B_2 + B3*B_3)/(Bx^2 + By^2 + Bz^2 + 1e-12)

# Here do the averaged fields corresponding to the node-projected arrays - otherwise, the removal of the n=0 component from the perturbation is not exact
B1t = dropdims(mean(B1,dims=3),dims=3);
B2t = dropdims(mean(B2,dims=3),dims=3);
B3t = dropdims(mean(B3,dims=3),dims=3);
B_1t = dropdims(mean(B_1,dims=3),dims=3);
B_2t = dropdims(mean(B_2,dims=3),dims=3);
psit = dropdims(mean(psi,dims=3),dims=3);
Jact = dropdims(mean(Jac,dims=3),dims=3);

# Calculate the node-based poloidal field (with jacobian factor): |J|B_p^2
JBpsqN = (B1t.*B_1t) .+ (B2t.*B_2t); # Sometimes the r=0 value may be slightly < 0 due to small errors when extrapolating from the first cell
JBpN = sqrt.(abs.(JBpsqN));

# Construct cylindrical coordinate arrays (H is the cylindrical Z toroidally averaged)
R = dropdims(mean(sqrt.(X.*X .+ Y.*Y),dims=3),dims=3);
H = dropdims(mean(Z,dims=3),dims=3);

# Interpolate field quantities on node-based grid
B1_eint = NodeInterpolation(B1);
B2_eint = NodeInterpolation(B2);
B3_eint = NodeInterpolation(B3);
B_1_eint = NodeInterpolation(B_1);
B_2_eint = NodeInterpolation(B_2);
B1t_eint = NodeInterpolation(B1t);
B2t_eint = NodeInterpolation(B2t);
B3t_eint = NodeInterpolation(B3t);
B_1t_eint = NodeInterpolation(B_1t);
B_2t_eint = NodeInterpolation(B_2t);
JBpsqN_eint = NodeInterpolation(JBpsqN);
JBpN_eint = NodeInterpolation(JBpN);
Jac_eint = NodeInterpolation(Jact);

# Interpolate Cartesian + Cylindrical grids
X_eint = extrapolate(scale(Interpolations.interpolate(X,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))))),rn,un,phin), (Line(),Periodic(),Periodic()));
Z_eint = extrapolate(scale(Interpolations.interpolate(Z,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))))),rn,un,phin), (Line(),Periodic(),Periodic()));
R_eint = extrapolate(scale(Interpolations.interpolate(R,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))))),rn,un), (Line(),Periodic()));
H_eint = extrapolate(scale(Interpolations.interpolate(H,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))))),rn,un), (Line(),Periodic()));

# Normalize psi_t by finding the magnetic axis and an X-point at each time step
psitnorm, Xpnts, raxisind, uaxisind, raxisarray, uaxisarray, rxpntarray, uxpntarray = NormalizePsi(psit,JBpsqN_eint);
if length(Xpnts) != tdim; global error_code = 2; end
ErrorCheck()

XpntsCar = Tuple{Float64, Float64}[]; # Cartesian grid locations of all selected x-points
AxisCar = Tuple{Float64, Float64}[]; # Cartesian grid coordinates of the magnetic axis

psit_eint = NodeInterpolation(psitnorm);

for (ip, p) in enumerate(Xpnts) # Confirm that the normalization was successful
    if abs(psit_eint(raxisarray[ip],uaxisarray[ip],ip)) > 1e-4 && Xpnts[ip] != [0,0]
        println("\nWarning! Poor Normalization of psi at t = "*string(ip)*" with psi_axis = $(psit_eint(raxisarray[ip],uaxisarray[ip],ip))")
    end
    if abs(psit_eint(rxpntarray[ip],uxpntarray[ip],ip) - 1.) > 1e-4 && Xpnts[ip] != [0,0]
        println("\nWarning! Poor Normalization of psi at t = "*string(ip)*" with psi_xpnt = $(psit_eint(rxpntarray[ip],uxpntarray[ip],ip))")
    end
    push!(XpntsCar, (R_eint(rxpntarray[ip],uxpntarray[ip]), H_eint(rxpntarray[ip],uxpntarray[ip])));
    push!(AxisCar, (R_eint(raxisarray[ip],uaxisarray[ip]), H_eint(raxisarray[ip],uaxisarray[ip])));
end

# Set up magnetic coordinate radial indices
psi_max_ind = 98; # Avoid separatrix by 3 cells
psi_min_ind = 3; # Avoid magnetic axis by 1 cell
psi_range = psi_max_ind-psi_min_ind+1; # number of psi grid points used for analysis

# Psi-magnetic angle grid dimensions
psidim = 101;
ufdim = udim; # ufdim same as udim to avoid aliasing

# Straight field line coordinate grid
pn = LinRange(0.0,1.0,psidim); # Grid for the radial psi coordinate with points at psi = 1%, 2%, etc.
ufn = LinRange(0.0,2.0*pi,ufdim);

# Obtain logical coordinates of all flux surfaces to construct line elements for the q integral
MPI.Barrier(comm)
if rank == 0 ; println("Starting to contour psi isosurfaces.") end
surfaces = PsiContour(psitnorm);

if plot_surfaces # Dump plots of the flux surfaces from the integrator
    MPI.Barrier(comm)
    if rank == 0 ; println("Starting to plot flux surfaces.") end
    PlotSurfaces(surfaces);
    MPI.Barrier(comm)
    if rank == 0 ; SurfaceAnimation() end
    if rank == 0 && q_only == false && plot_grid == false && plot_q == false; runtime, units = ElapsedTime(start) ; println("All flux surfaces have been plotted! Runtime = $(runtime) "*units*". Exiting") end
    if q_only == false && plot_grid == false && plot_q == false
	MPI.Finalize()
	exit(0)
    end
end
MPI.Barrier(comm)

# Calculate q-profile and the magnetic angles by integrating around the flux surfaces
q, theta_m, pathlength = q_profile(surfaces);
q_sign = 1.
if q[1,1] < 0 ; q .= -q ; q_sign = -1. end # Flip q to positive if needed
q_str = stretch_q(q);

# Diagnostic for theta_m
for t in 1:tdim
    for p in 1:psi_range
	if abs(theta_m[p,t][end] - 2pi) > 1e-3 && Xpnts[t] != [0,0]
	    println("\nError! theta_m != 2pi after a full circuit at t = $(rank_tstart+t), psi = $(0.01*(p+psi_min_ind-1)) : theta_m = $(theta_m[p,t][end])")
	end
	if maximum(pathlength[p,t]) != pathlength[p,t][end] || minimum(pathlength[p,t]) != pathlength[p,t][1]
	    println("\nError! pathlength is not monotonic at t = $(rank_tstart+t), psi = $(0.01*(p+psi_min_ind-1))")
	end
    end
end

MPI.Barrier(comm)
if rank != 0
    MPI.Send(q_str,0,rank,comm)
elseif rank == 0 # Write aggregate data from all processes to a single output NumPy file
    q_dim1 = size(q_str)[1];
    q_dim2 = size(q_str)[2];
    q_composite = Array{Float64,2}(undef,q_dim1,0);
    q_composite = cat(dims=2,q_composite,q_str);
    for r in range(1,stop=size_comm-2)
        q_recv = Array{Float64,2}(undef,q_dim1,q_dim2);
        MPI.Recv!(q_recv,r,r,comm)
        global q_composite = cat(dims=2,q_composite,q_recv);
    end
    if size_comm > 1
        q_recv_fin = Array{Float64,2}(undef,q_dim1,final_chunk_size);
        MPI.Recv!(q_recv_fin,size_comm-1,size_comm-1,comm)
        q_composite = cat(dims=2,q_composite,q_recv_fin);
    end
    q_file = "/q_"*string(tstart)*"_"*string(tend)*".npy";
    npzwrite(filepath * q_file ,q_composite);
end
MPI.Barrier(comm)

if rank == 0 && size(q_composite) != (101,tend-tstart) ; global error_code = 5 end
ErrorCheck()
if rank == 0 && plot_q == true ; println("Now plotting q profiles.") ; Plotq(q_composite) end

if q_only && plot_grid == false # Exit here if only q(psi,t) is desired
    if rank == 0 ; runtime, units = ElapsedTime(start) ; println("All done! Runtime = $(runtime) "*units*". Exiting.") end
    MPI.Finalize()
    exit(0)
end

# Create the new magnetic field-aligned grid (psi,theta_m,phi)
GD = new_grid(surfaces,theta_m,pathlength);
ErrorCheck()
MPI.Barrier(comm)
if rank == 0 ; println("Magnetic grid has been created.") end

if plot_grid
    PlotGrid(GD,surfaces)
    MPI.Barrier(comm)
    if rank == 0 ; GridAnimation() end 
    if rank == 0 ; runtime, units = ElapsedTime(start) ; println("All done! Runtime = $(runtime) "*units*". Exiting.") end
    MPI.Finalize()
    exit(0)
end

# Calculate perturbations in the magnetic field components
na = [CartesianIndex()]; # this provides padding for a symmetric dimension 
dB1 = B1.-B1t[:,:,na,:];
dB2 = B2.-B2t[:,:,na,:];    
dB3 = B3.-B3t[:,:,na,:];

# Normalized radial perturbation in the logicial coordinates (r,theta,phi)
b_hat_rho = Jac.*(dB1.*B2t[:,:,na,:] .- dB2.*B1t[:,:,na,:])./B3t[:,:,na,:];

# Normalized toroidal perturbation in the logicial coordinates (r,theta,phi)
db3 = dB3./B3t[:,:,na,:];

# Interpolate b_hat_rho avoiding problematic r=0 point
b_hat_rho_eint = extrapolate(scale(Interpolations.interpolate(b_hat_rho[2:end,:,:,:],(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),NoInterp())),rn[2:end],un,phin,tn), (Line(),Periodic(),Periodic(),Interpolations.Throw()));

# Interpolate dB3 
db3_eint = extrapolate(scale(Interpolations.interpolate(db3[2:end,:,:,:],(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),NoInterp())),rn[2:end],un,phin,tn), (Line(),Periodic(),Periodic(),Interpolations.Throw()));

# Project in new coordinates 
b_hat_rho_P = project(b_hat_rho_eint,4);
db3_P = project(db3_eint,4);
B3t_P = project(B3t_eint,3);

# Project Cartesian Maps
X_P = project(X_eint,0);
Z_P = project(Z_eint,0);

if rank == 0 ; println("Finished projection\n") end

# Interpolate b_hat_rho_P, avoid the axis so that the grid is evenly spaced
b_hat_rho_P_eint = extrapolate(scale(Interpolations.interpolate(b_hat_rho_P[2:end,:,:,:],(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),NoInterp())),pn[psi_min_ind:psi_max_ind],ufn,phin,tn), (Line(),Periodic(),Periodic(),Interpolations.Throw()));

# Interpolate dB3_P 
db3_P_eint = extrapolate(scale(Interpolations.interpolate(db3_P[2:end,:,:,:],(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),NoInterp())),pn[psi_min_ind:psi_max_ind],ufn,phin,tn), (Line(),Periodic(),Periodic(),Interpolations.Throw()));

# Interpolate B3t_P
B3t_P_eint = extrapolate(scale(Interpolations.interpolate(B3t_P[2:end,:,:],(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),NoInterp())),pn[psi_min_ind:psi_max_ind],ufn,tn), (Line(),Periodic(),Interpolations.Throw()));

# Interpolate Cartesian Maps
X_P_eint = extrapolate(scale(Interpolations.interpolate(X_P[2:end,:,:,:], (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),NoInterp())),pn[psi_min_ind:psi_max_ind],ufn,phin,tn), (Line(),Periodic(),Periodic(),Interpolations.Throw()));
Z_P_eint = extrapolate(scale(Interpolations.interpolate(Z_P[2:end,:,:,:], (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),NoInterp())),pn[psi_min_ind:psi_max_ind],ufn,phin,tn), (Line(),Periodic(),Periodic(),Interpolations.Throw()));

# Extrapolate in whole psi range and multiply with q(psi,uf,phi,t)
db3_P_ext = db3_P_eint(pn,ufn,phin,tn);
b_hat_rho_P_ext = b_hat_rho_P_eint(pn,ufn,phin,tn);
B3t_P_ext = B3t_P_eint(pn,ufn,tn);
X_P_ext = X_P_eint(pn,ufn,phin,tn);
Z_P_ext = Z_P_eint(pn,ufn,phin,tn);

b_hat_rho_ncq = q_str[:,na,na,:].*q_sign.*psi_sign.*b_hat_rho_P_ext; # Watch out for this sign, depends on original psi monotonicity.
db3_hat_P_ext = q_str[:,na,na,:].*q_sign.*db3_P_ext;

# Output files
X_file = "/X.npy";
Z_file = "/Z.npy";
db3_file = "/db3_"*string(tstart)*"_"*string(tend)*".npy";
B3t_file = "/B3t_"*string(tstart)*"_"*string(tend)*".npy";
b_hat_rho_file = "/b_hat_rho_"*string(tstart)*"_"*string(tend)*".npy";

# Send the array pieces to root and save
MPI.Barrier(comm)
if rank != 0 # Send to root from all other processes with unique message
    MPI.Send(b_hat_rho_ncq,0,rank+1001,comm)
    MPI.Send(X_P_ext,0,rank+1002,comm)
    MPI.Send(Z_P_ext,0,rank+1003,comm)
    MPI.Send(db3_hat_P_ext,0,rank+1004,comm)
    MPI.Send(B3t_P_ext,0,rank+1005,comm)
elseif rank == 0 # Receive to root process
    # dimension initialization
    dim1 = size(b_hat_rho_ncq)[1]
    dim2 = size(b_hat_rho_ncq)[2]
    dim3 = size(b_hat_rho_ncq)[3]
    dim4 = size(b_hat_rho_ncq)[4]    
    
    # Array initialization
    bhat_composite = Array{Float64,4}(undef,dim1,dim2,dim3,0)
    X_composite = Array{Float64,4}(undef,dim1,dim2,dim3,0)
    Z_composite = Array{Float64,4}(undef,dim1,dim2,dim3,0)
    db3_composite = Array{Float64,4}(undef,dim1,dim2,dim3,0)
    B3t_composite = Array{Float64,3}(undef,dim1,dim2,0)
    
    # Concatenate with element in the root process along the time dimension
    bhat_composite = cat(dims=4,bhat_composite,b_hat_rho_ncq)
    X_composite = cat(dims=4,X_composite,X_P_ext)
    Z_composite = cat(dims=4,Z_composite,Z_P_ext)
    db3_composite = cat(dims=4,db3_composite,db3_hat_P_ext)
    B3t_composite = cat(dims=3,B3t_composite,B3t_P_ext)
    
    # Receive from all other processes
    for r in range(1,stop=size_comm-2)
        # Initialize receiver buffer
        bhat_recv = Array{Float64,4}(undef,dim1,dim2,dim3,dim4)
        X_recv = Array{Float64,4}(undef,dim1,dim2,dim3,dim4)
        Z_recv = Array{Float64,4}(undef,dim1,dim2,dim3,dim4)
        db3_recv = Array{Float64,4}(undef,dim1,dim2,dim3,dim4)
        B3t_recv = Array{Float64,3}(undef,dim1,dim2,dim4)
        # Receive
        MPI.Recv!(bhat_recv,r,r+1001,comm)
        MPI.Recv!(X_recv,r,r+1002,comm)
        MPI.Recv!(Z_recv,r,r+1003,comm)
        MPI.Recv!(db3_recv,r,r+1004,comm)
        MPI.Recv!(B3t_recv,r,r+1005,comm)
        # Concatenate received array to global one
        global bhat_composite = cat(dims=4,bhat_composite,bhat_recv)
        global X_composite = cat(dims=4,X_composite,X_recv)
        global Z_composite = cat(dims=4,Z_composite,Z_recv)
        global db3_composite = cat(dims=4,db3_composite,db3_recv)
        global B3t_composite = cat(dims=3,B3t_composite,B3t_recv)
    end
    if size_comm > 1
        # Receive from final process which has variable number of elements
        # Initialize the buffer
        bhat_recv_fin = Array{Float64,4}(undef,dim1,dim2,dim3,final_chunk_size)
        X_recv_fin = Array{Float64,4}(undef,dim1,dim2,dim3,final_chunk_size)
        Z_recv_fin = Array{Float64,4}(undef,dim1,dim2,dim3,final_chunk_size)
        db3_recv_fin = Array{Float64,4}(undef,dim1,dim2,dim3,final_chunk_size)
        B3t_recv_fin = Array{Float64,3}(undef,dim1,dim2,final_chunk_size)
        # Receive
        MPI.Recv!(bhat_recv_fin,size_comm-1,size_comm-1+1001,comm)
        MPI.Recv!(X_recv_fin,size_comm-1,size_comm-1+1002,comm)
        MPI.Recv!(Z_recv_fin,size_comm-1,size_comm-1+1003,comm)
        MPI.Recv!(db3_recv_fin,size_comm-1,size_comm-1+1004,comm)
        MPI.Recv!(B3t_recv_fin,size_comm-1,size_comm-1+1005,comm)
        # Attach final chunk to the global array
        bhat_composite = cat(dims=4,bhat_composite,bhat_recv_fin)
        X_composite = cat(dims=4,X_composite,X_recv_fin)
        Z_composite = cat(dims=4,Z_composite,Z_recv_fin)
        db3_composite = cat(dims=4,db3_composite,db3_recv_fin)
        B3t_composite = cat(dims=3,B3t_composite,B3t_recv_fin)
    end
    bhatrho = Float64.(bhat_composite);
    Xproj = Float64.(X_composite);
    Zproj = Float64.(Z_composite);
    db3proj = Float64.(db3_composite);
    B3tproj = Float64.(B3t_composite);
    
    npzwrite(filepath * X_file, Xproj)
    npzwrite(filepath * Z_file, Zproj)
    npzwrite(filepath * db3_file, db3proj)
    npzwrite(filepath * B3t_file, B3tproj)
    npzwrite(filepath * b_hat_rho_file, bhatrho)
end    

if rank == 0 && plot_pert ; println("A 2D animation of the magnetic perturbations will now be made.\n") end
if rank == 0 && plot_pert ; PlotdB(bhatrho,1,1) ; PlotdB(db3proj,3,1) ; end

if rank == 0 ; runtime, units = ElapsedTime(start) ; println("All done! Runtime = $(runtime) "*units*". Exiting.") end
MPI.Finalize() # End of script.
exit(0)

result = Array{Float64}(undef,time_interval)
for t in 1:time_interval
	result[t] = imag(rfft(bhatrho[6,:,:,t])[2,2])
end
