# magnetic_coordinates.jl
using Interpolations
using Dierckx
using DifferentialEquations
using MPI
using NPZ
using Statistics
using FileIO
using HDF5
using EllipsisNotation
using PyPlot
using Plots
            
#=
This module contains the functions needed to do the integration of theta_M around a flux surface.
This method proved to be not robust enough, as the system of ODE's is stiff near the axis.
The contour method used in magnetic_coordinates.jl is both faster and more reliable.
This module is kept for posterity.

Authored by Jason Hamilton
Contact: jmhamilton@lanl.gov
        
Latest update: 26th July 2023

=#  

# Set up of parallel communicator
MPI.Init()
comm = MPI.COMM_WORLD
size_comm = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)
error_code = 0
if rank == 0 ; start = time() ; end

# User Inputs (defaults) -------------------------------------------------------------------------------

filepath = pwd(); # Directory of .h5 file
filename = "/pixie3d.h5"; # This is the default filename from pixplot. Note the leading slash.
q_only = false; # Set to true to only calculate q(psi,t) and then exit the script
plot_surfaces = false; # When true, the flux surfaces will be plotted and saved, and then the script will exit
tstart = 0; # tstart and tend specify the range of timestep indices that will be processed by this script: [tstart+1:tend]
tend = nothing; # Will self-limit itself to the total number of time steps in the HDF5 file!

# ------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------

# Read input options to override defaults
for option in ARGS
    if startswith(option,"filepath=")
        global filepath = option[length("filepath=")+1:end]
    elseif startswith(option,"filename=")
        global filename = option[length("filename=")+1:end]
    elseif startswith(option,"isShaped=")
        global isShaped = ( option[length("isShaped=")+1:end] == "true" )
    elseif startswith(option,"q_only=")
        global q_only = ( option[length("q_only=")+1:end] == "true" )
    elseif startswith(option,"plot_surfaces=")
        global plot_surfaces = ( option[length("plot_surfaces=")+1:end] == "true" )
    elseif startswith(option,"tstart=")
        global tstart = parse(Int,option[length("tstart=")+1:end])
    elseif startswith(option,"tend=")
        global tend = parse(Int,option[length("tend=")+1:end])
    elseif startswith(option,"SplineSensitivity=")
        global SplineSensitivity = parse(Float64,option[length("SplineSensitivity=")+1:end])
    else
        if rank == 0 ; println("Warning! Your option: "*option*" , was not a valid option.") end
    end
end

##################################################################   FUNCTIONS   #####################################################################

function Flux_surface!(du,u,p,t)
    """
    Defines the flux surface system of equations du/dt = f(u,t,p) .
    Here, p[1] = time , u = [r_surface,u_surface,theta_mag*q], t = line element of flux surface.
    """
    du[1] = B1t_eint(u[1],u[2],p[1])/JBpN_eint(u[1],u[2],p[1])
    du[2] = B2t_eint(u[1],u[2],p[1])/JBpN_eint(u[1],u[2],p[1])
    du[3] = B3t_eint(u[1],u[2],p[1])/JBpN_eint(u[1],u[2],p[1])
end

function condition_pi(u,t,integrator)
    """
    Callback: checks when the integrator has traversed to the opposite side of the origin.
    Preserve the sign of the argument so that the integrator can detect when the condition crosses zero.
    """
    mod(mod(u[2],2pi) - integrator.p[3] - pi, 2pi)*sign(mod(u[2],2pi) - integrator.p[3] - pi)
end

function condition_zero(u,t,integrator)
    """
    Callback: checks when the integrator has crossed the starting angle.
    Preserve the sign of the argument so that the integrator can detect when the condition crosses zero.
    """
    mod(mod(u[2],2pi) - integrator.p[3], 2pi)*sign(mod(u[2],2pi) - integrator.p[3])
end

function condition_stop(u,t,integrator)
    """
    Callback: checks if the above conditions have been met twice and the integrator has returned to the starting position (within the resolution of the grid).
    """
    integrator.p[4] >= 2 && abs(u[1] - integrator.p[2]) < 0.5*rn[2] && (abs(mod(u[2],2pi) - integrator.p[3]) < 0.5*un[2] || abs(abs(mod(u[2],2pi) - integrator.p[3]) - 2pi) < 0.5*un[2])
end

function condition_r_reg(u,t,integrator)
    """
    Callback: checks if 'r' gets out of bounds i.e. r < 0 .
    """
    u[1] < 1e-4
end

function condition_theta_reg1(u,t,integrator)
    """
    Callback: checks if 'theta' gets out of bounds i.e. theta < 0 .
    """
    u[2] < 0
end

function condition_theta_reg2(u,t,integrator)
    """
    Callback: checks if 'theta' gets out of bounds i.e. theta > 2pi .
    """
    u[2] > 2pi
end

function affect_pi!(integrator)
    """
    Counts the number of times the integrator has met the above conditions.
    This takes into account small flux surfaces that do not encompass the origin (shifted axis).
    """
    integrator.p[4] = integrator.p[4]+1
    integrator.p[5] = integrator.p[5]+1
end

function affect_zero!(integrator)
    """
    Counts the number of times the integrator has met the above conditions.
    This takes into account small flux surfaces that do not encompass the origin (shifted axis).
    """
    integrator.p[4] = integrator.p[4]+1
end

function affect_stop!(integrator)
    """
    Tells the integrator to stop.
    """
    terminate!(integrator)
end

function affect_r_reg!(integrator)
    """
    Flips negative 'r' to positive 'r'.
    """
    integrator.u[1] = max(integrator.u[1],1e-4)
end

function affect_theta_reg!(integrator)
    """
    Keeps 'theta' in the range (0,2pi).
    """
    integrator.u[2] = mod(integrator.u[2],2pi)
end

function fs_integration(rs::Float64,us::Float64,time::Int,index::Int)
    """
    Integrates around the flux surface to solve for q. Requires the DifferentialEquations package.
    The starting angle 'us' will correspond to theta_mag = 0.
    """
    tolerance = 1.e-12 # Set the error tolerance for the integrator (smaller is better)
    u0 = [rs,us,0.0] # Starting position for integration, phi = 0 is chosen arbitrarily
    p = [time,rs,us,0,0] # Parameter set: [time step value, starting radius, starting angle, number of times the integrator has crossed theta = 0,pi]
    tspan = (0.0,50.0) # Here, 'dt' is the integrator's line element, and tspan is its domain (not to be confused with the time step)
    cb1 = DiscreteCallback(condition_r_reg,affect_r_reg!)
    cb2 = DiscreteCallback(condition_theta_reg1,affect_theta_reg!)
    cb3 = DiscreteCallback(condition_theta_reg2,affect_theta_reg!)
    cb4 = ContinuousCallback(condition_pi,affect_pi!,interp_points=10) # If the integrator completes a half-circuit, count the event
    cb5 = ContinuousCallback(condition_zero,affect_zero!,interp_points=10) # If the integrator crosses the starting angle, count the event
    cb6 = DiscreteCallback(condition_stop,affect_stop!) # If the integrator completes a full 2pi circuit, stop!
    prob = ODEProblem(Flux_surface!,u0,tspan,p) # Define the system of equations to be solved
    cbs = CallbackSet(cb1,cb2,cb3,cb4,cb5,cb6) # Define what conditions end or change the integration
    integrator = init(prob,Tsit5(),callback=cbs,reltol=tolerance,abstol=tolerance) # This holds the integrator as an object that can be manipulated during the solve
    sol = solve(prob,Tsit5(),callback=cbs,reltol=tolerance,abstol=tolerance) # Perform the solve
    if integrator.p[4] < 2 # Debugging when the solver can't complete a full surface
        println("Integrator failed at: t = $(rank_tstart+time), r0 = $(rs), r = $(last(sol[1,length(sol.t)])), theta= $(last(sol[2,length(sol.t)])), psi = $(psi_min_ind+index-2), rank = $(rank), half-circuits = $(integrator.p[4]), theta = pi crossings = $(integrator.p[5]), Bpol = $(JBpN_eint(last(sol[1,length(sol.t)]),last(sol[2,length(sol.t)]),time))")
    end
    if integrator.p[5] > 0 && psi_min_ind+index-2 < psit_eint(0.,0.,time)*100.
        println("Integrator incorrectly went around the origin at: t = $(rank_tstart+time), r0 = $(rs), r = $(last(sol[1,length(sol.t)])), theta= $(last(sol[2,length(sol.t)])), psi = $(psi_min_ind+index-2), rank = $(rank), half-circuits = $(integrator.p[4]), theta = pi crossings = $(integrator.p[5]), Bpol = $(JBpN_eint(last(sol[1,length(sol.t)]),last(sol[2,length(sol.t)]),time))")
    end
    return sol
end

function q_profile2(psit_eint)
    """
    Calculates the q profile by integrating around the flux surfaces intersecting (r(psi),theta(axis)) at the specified psi values.
    The 3rd component of the solution corresponds to theta_mag * q, and the last element of the integration is when theta_mag = 2pi.
    """
    if rank == 0 ; println("Starting to calculate the q profile.") end
    q = Array{Float64,2}(undef,psi_range,tdim)

    for time in range(1, stop=tdim)
        for (i,r) in enumerate(r_of_psi_array[psi_min_ind:psi_max_ind,time])
            sol = fs_integration(r,mod(un[uaxisind[time]],2pi),time,i) # Performs the integration starting from the point (r,theta_axis), returns an array
            q[i,time] = last(sol[3,length(sol.t)])/(2.0*pi) # After a complete circuit around the flux surface, theta_mag = 2pi and this solution corresponds to 2pi*q
            if abs(psit_eint(last(sol[1,length(sol.t)]),last(sol[2,length(sol.t)]),time)*100. - psi_min_ind - i + 2) > 1e-1 # Check that the integrator stayed on the psi contour
                println("The contour is inaccurate at: t = $(rank_tstart+time), r0 = $(r), r = $(last(sol[1,length(sol.t)])), theta = $(mod(last(sol[2,length(sol.t)]),2pi)) , psi = $(psi_min_ind+i-2), psit_eint = $(psit_eint(last(sol[1,length(sol.t)]),last(sol[2,length(sol.t)]),time)*100.), rank = $(rank)")
            end
            if abs(last(sol[1,length(sol.t)]) - r) > 0.5*rn[2] || (abs(mod(last(sol[2,length(sol.t)]),2pi) - mod(un[uaxisind[time]],2pi)) > 0.5*un[2] && abs(abs(mod(last(sol[2,length(sol.t)]),2pi) - mod(un[uaxisind[time]],2pi)) - 2pi) > 0.5*un[2]) # Check that the integrator stopped at the same place it began
                println("The integrator did not return to the starting position at: t = $(rank_tstart+time), r0 = $(r), psi = $(psi_min_ind+i-2), r = $(last(sol[1,length(sol.t)])), theta = $(mod(last(sol[2,length(sol.t)]),2pi)), rank = $(rank)")
            end
        end
        println("The q profile is finished at t = $(rank_tstart+time)!")
    end

    if rank == 0 ; println("The q profiles are finished!") end

    return q
end

# Creates dictionary of grid in the form {t,psi,theta_mag} -> (r,theta)
function new_grid2()
    Grid_Dict = Dict{Tuple{Int64,Float64,Float64},Tuple{Float64,Float64}}()
    for t in range(0, stop=tdim-1)
        for (i,r) in enumerate(r_of_psi_array[psi_min_ind:psi_max_ind,t+1])
            sol = fs_integration(r,0.0,t,i)
            uf_sol = sol[3,:] # Solution to integral unnormalized at all points 't' along the path spanning 't'
            q_redef = last(sol[3,length(sol.t)])/(2.0*pi)
            uf_rsc = uf_sol/q_redef # Normalized solution
            l = sol.t # Parameter denoting position along flux surface/integration path
            l_of_uf_rsc = Spline1D(sort!(uf_rsc[1:end-1]),l[1:end-1],bc="extrapolate",s=1e-3) # Obtain l(theta_m)
            ls = l_of_uf_rsc(ufn) # evaluate l(theta_m) at the desired resolution of theta_m between 0,2pi
            polar_coords = sol(ls,idxs=1:2) #obtain the r,theta coordinates of the corresponding positions along the path
            psi_val = pn[i+psi_min_ind-1]
            for j in 1:udim # Equivalance between each magnetic grid point and the logical grid coordinates
                Grid_Dict[(Int64(round(t)),round(psi_val,digits=2),Float64(ufn[j]))] = (polar_coords[j][1],polar_coords[j][2])
            end
        end
    end
    return Grid_Dict
end

# All project functions work with the new grid calculated in the above function
function project(arr_eint)
    ArrP = []
    for t in tn
        for p in pn[1:psi_max_ind]
            for uf in ufn
                for phi in phin
                    if p == 0.0 # Adding magnetic axis point
                        r = rmaxis[Int64(round(t+1))]
                        u = umaxis[Int64(round(t+1))]
                        if tdim > 1
                            ArrVal = arr_eint(r,u,phi,t)
                        else
                            ArrVal = arr_eint(r,u,phi)
                        end
                        append!(ArrP,ArrVal)
                    else
                        r = GD[Int64(round(t)),round(p,digits=2),Float64(uf)][1]
                        u = GD[Int64(round(t)),round(p,digits=2),Float64(uf)][2]
                        if tdim > 1
                            ArrVal = arr_eint(r,u,phi,t)
                        else
                            ArrVal = arr_eint(r,u,phi)
                        end
                        append!(ArrP,ArrVal)
                    end
                end
            end
        end
    end
    Arr = permutedims(reshape(ArrP,size(phin)[1],size(ufn)[1],psi_range+1,tdim),(3,2,1,4))
    return Arr
end

function project_toroidally_symmetric(arr_eint)
    ArrP = []
    for t in tn
        for p in pn[1:psi_max_ind]
            for uf in ufn
                if p == 0.0 # Adding magnetic axis point
                    r = rmaxis[Int64(round(t+1))]
                    u = umaxis[Int64(round(t+1))]
                    if tdim > 1
                        ArrVal = arr_eint(r,u,t)
                    else
                        ArrVal = arr_eint(r,u)
                    end
                    append!(ArrP,ArrVal)
                else
                    r = GD[Int64(round(t)),round(p,digits=2),Float64(uf)][1]
                    u = GD[Int64(round(t)),round(p,digits=2),Float64(uf)][2]
                    if tdim > 1
                        ArrVal = arr_eint(r,u,t)
                    else
                        ArrVal = arr_eint(r,u)
                    end
                    append!(ArrP,ArrVal)
                end
            end
        end
    end
    Arr = permutedims(reshape(ArrP,size(ufn)[1],psi_range+1,tdim),(2,1,3))
    return Arr
end

function project_Cartesian_Map(arr_eint)
    ArrP = []
    for t in tn
        for p in pn[1:psi_max_ind]
            for uf in ufn
                for phi in phin
                    if p == 0.0 # Adding magnetic axis point
                        r = rmaxis[Int64(round(t+1))]
                        u = umaxis[Int64(round(t+1))]
                        ArrVal = arr_eint(r,u,phi)
                        append!(ArrP,ArrVal)
                    else
                        r = GD[Int64(round(t)),round(p,digits=2),Float64(uf)][1]
                        u = GD[Int64(round(t)),round(p,digits=2),Float64(uf)][2]
                        ArrVal = arr_eint(r,u,phi)
                        append!(ArrP,ArrVal)
                    end
                end
            end
        end
    end
    Arr = permutedims(reshape(ArrP,size(phin)[1],size(ufn)[1],psi_range+1,tdim),(3,2,1,4))
    return Arr
end


# Interpolate across psi to obtain r(psi)
r_of_psi_array = Array{Float64, 2}(undef, 101, tdim)

for t in range(1,stop=tdim)
    psi_sample = Array{Float64}(undef, rdim + 1 - raxisind[t]) # Initialize arrays with number of node points between the axis and wall along theta_axis
    r_sample = Array{Float64}(undef, rdim + 1 - raxisind[t])
    psi_sample[1] = 0. # Poloidal flux at the magnetic axis is zero by definition
    r_sample[1] = raxisarray[t] # The true location of the axis based on the interpolated psitnorm function
    psi_sample[2:end] = psitnorm[raxisind[t] + 1 : end, uaxisind[t], t] # Values of psi_t along chord from outside the axis to the boundary
    r_sample[2:end] = rn[raxisind[t] + 1 : end] # Values of logical coordinate 'r' along the same chord at the same points
    lastindex = findfirst(x -> x > 1.0, psi_sample) # Find index of when psi exceeds 1.0, don't interpolate past this index

    # Interpolate over r(psi) - must use Dierckx package to get a cubic interpolation over the irregular psi grid
    r_psi_spline = Dierckx.Spline1D(psi_sample[begin:lastindex],r_sample[begin:lastindex],bc="extrapolate")
    r_of_psi_array[:,t] = r_psi_spline(pn) 
end

r_of_psi_array[r_of_psi_array .< 0.] .= 0. # Again, extrapolation to r=0 may cause small negative numbers


