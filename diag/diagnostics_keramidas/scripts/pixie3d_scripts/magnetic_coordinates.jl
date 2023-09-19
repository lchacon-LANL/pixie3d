using Interpolations
using Dierckx
using DifferentialEquations
using PyCall
using MPI
using NPZ
using Statistics
using QuadGK
using FileIO
using JLD2
MPI.Init()

#=Parallel code that calculates the magnetic perturbation projected to straight field line coordinates. It produces a series 
of arrays (b_hat_rho is the most important) with quantities defined in a psi,theta_magnetic,phi coordinate system.

For this programm to run correctly when submitted one needs to uncomment the line "matplotlib.use('Agg')" in the "pixie_read_st" module. It requires an initial guess for the X-point that is fed in the python module.  

RUN INSTRUCTIONS: The program is called inside an allocation as "mpiexec -np 4 julia magnetic_coordinates.jl"

Inside a single node, the program runs into segmentation faults if more than 200 time steps (here we mean time steps of the 
pixie3d.h5 file so 200 time steps of it might correspond to 2000 simulation time steps, if results are dumped every 10 time steps) are processed. Therefore, it is useful to specify the duration in the names of the output files.
=#

# Set up of parallel communicator
comm = MPI.COMM_WORLD
size_comm = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)

# Inputs

#filepath = "/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/sawtooth2.scratch/"
#filepath = "/net/scratch3/giannis_kx/pixie3d/tests/sawtooth/sawtooth.scratch/"
#filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_tear/dt_sh_m3_n2.scratch/"
#filepath = "/net/scratch4/giannis_kx/pixie3d/iter/int_kink/11/11_new_visc.scratch/"
#filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_new/dt2.scratch./"
#filepath = "/net/scratch3/chacon/pixie3d/EFIT/ITER/ITER3-chipar/3d/pixie3d-iter-SN-fr_11-refined.scratch/"
#filepath = "/net/scratch4/giannis_kx/pixie3d/iter/int_kink/11/11_visc_old_nodiff/11_visc_old_nodiff.scratch/"
#filepath = "/net/scratch3/chacon/pixie3d/EFIT/ITER/ITER3-chipar/3d/pixie3d-iter-SN-fr_11-n0=300x.scratch/"
#filepath = "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/300x/t10000/pixie3d-n0=300x.scratch/"
#filepath = "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/db_tear/dt-chi-par/"
filepath = "/lustre/scratch5/.mdt1/jmhamilton/2022/pixie3d-ppcf17/pixie3d-ppcf17.scratch/magneticcoord/"

#metricpath = "/net/scratch3/giannis_kx/FTLE/11/shaped_metric_coeff6.npz"

tstart = 0 #800;
tend = 200 #980;
psi_max_ind = 97; # Avoid separatrix
psi_min_ind = 2; # Avoid magnetic axis
psi_range = psi_max_ind-psi_min_ind+1;
x0 = (0.9,0.44); # Initial guess for the X-point

#Output files
q_file = "q_mpi_0_200.npy";
X_file = "X.npy"
Z_file = "Z.npy"
db3_file = "db3_mpi_0_200.npy"
B3t_file = "B3t_mpi_0_200.npy"
b_hat_rho_file = "b_hat_rho_mpi_0_200.npy"

function dAdr(arr::Array{Float64,4},theta::Int,phi::Int,t::Int)
    dr = 1.0 / size(arr,1)
    p = arr[:,theta,phi,t]
    theta_mid = convert(Int,ceil(udim/2))
    if theta<theta_mid
        mid = theta + convert(Int,(size(arr,2)-1)/2)
    else
        mid = theta - convert(Int,(size(arr,2)-1)/2)
    end
    p_pi = arr[:,mid,phi,t]
    L = length(p)
    dpsidr = zeros(L)
    for i in range(1,stop=L)
        if i==1
            dpsidr[i] = (-p[3]+8*p[2]-8*p_pi[2]+p_pi[3])/(12.0*dr)
        elseif i==2
            dpsidr[i] = (-p[4]+8*p[3]-8*p[1]+p_pi[2])/(12.0*dr)
        elseif i==L-1
            dpsidr[i] = (p[L]-p[L-2])/(2.0*dr)
        elseif i==L-2 
            dpsidr[i] = (-p[L]+8*p[L-1]-8*p[L-3]+p[L-4])/(12.0*dr)
        elseif i==L
            dpsidr[i] = (p[L]-p[L-1])/dr
        else
            dpsidr[i] = (-p[i+2]+8*p[i+1]-8*p[i-1]+p[i-2])/(12.0*dr)
        end
    end
    return dpsidr
end

function dAdt(arr::Array{Float64,4},r::Int,phi::Int,t::Int)
    dtheta = 2.0*pi / size(arr,2)
    p = arr[r,:,phi,t]
    L = length(p)
    dpsidtheta = zeros(L)
    for i in range(1,stop=L)
        if i==1
            dpsidtheta[i] = (-p[3]+8*p[2]-8*p[L-1]+p[L-2])/(12.0*dtheta)
        elseif i==2
            dpsidtheta[i] = (-p[4]+8*p[3]-8*p[1]+p[L-1])/(12.0*dtheta)
        elseif i==L-1
            dpsidtheta[L-1] = (-p[2]+8*p[1]-8*p[L-2]+p[L-3])/(12.0*dtheta)
        elseif i==L-2
            dpsidtheta[L-2] = (-p[1]+8*p[L-1]-8*p[L-3]+p[L-4])/(12.0*dtheta)
        elseif i==L
            dpsidtheta[L] = dpsidtheta[1]
        else
            dpsidtheta[i] = (-p[i+2]+8*p[i+1]-8*p[i-1]+p[i-2])/(12.0*dtheta)
        end
    end
    return dpsidtheta
end 
    
function dpsi_dr(arr::Array{Float64,4})
    r_deriv = Float64[] 
    for theta in range(1,stop=size(arr,2))
        for phi in range(1,stop=size(arr,3))
            for t in range(1,stop=size(arr,4))
                r_deriv = append!(r_deriv,dAdr(arr,theta,phi,t))
            end
        end
    end
     # In Julia the convention in reshape is the opposite of Python
    return permutedims(reshape(r_deriv,(size(arr,1),size(arr,4),size(arr,3),size(arr,2))),(1,4,3,2))
end

function dpsi_dtheta(arr::Array{Float64,4})
    t_deriv = Float64[] 
    for r in range(1,stop=size(arr,1))
        for phi in range(1,stop=size(arr,3))
            for t in range(1,stop=size(arr,4))
                t_deriv = append!(t_deriv,dAdt(arr,r,phi,t))
            end
        end
    end
    return permutedims(reshape(t_deriv,(size(arr,2),size(arr,4),size(arr,3),size(arr,1))),(4,1,3,2))
end

function dAdr(arr::Array{Float64,3},theta::Int,t::Int)
    dr = 1.0 / size(arr,1)
    p = arr[:,theta,t]
    theta_mid = convert(Int,ceil(udim/2))
    if theta<theta_mid
        mid = theta + convert(Int,(size(arr,2)-1)/2)
    else
        mid = theta - convert(Int,(size(arr,2)-1)/2)
    end
    p_pi = arr[:,mid,t]
    L = length(p)
    dpsidr = zeros(L)
    for i in range(1,stop=L)
        if i==1
            dpsidr[i] = (-p[3]+8*p[2]-8*p_pi[2]+p_pi[3])/(12.0*dr)
        elseif i==2
            dpsidr[i] = (-p[4]+8*p[3]-8*p[1]+p_pi[2])/(12.0*dr)
        elseif i==L-1
            dpsidr[i] = (p[L]-p[L-2])/(2.0*dr)
        elseif i==L-2 
            dpsidr[i] = (-p[L]+8*p[L-1]-8*p[L-3]+p[L-4])/(12.0*dr)
        elseif i==L
            dpsidr[i] = (p[L]-p[L-1])/dr
        else
            dpsidr[i] = (-p[i+2]+8*p[i+1]-8*p[i-1]+p[i-2])/(12.0*dr)
        end
    end
    return dpsidr
end

function dAdt(arr::Array{Float64,3},r::Int,t::Int)
    dtheta = 2.0*pi / size(arr,2)
    p = arr[r,:,t]
    L = length(p)
    dpsidtheta = zeros(L)
    for i in range(1,stop=L)
        if i==1
            dpsidtheta[i] = (-p[3]+8*p[2]-8*p[L-1]+p[L-2])/(12.0*dtheta)
        elseif i==2
            dpsidtheta[i] = (-p[4]+8*p[3]-8*p[1]+p[L-1])/(12.0*dtheta)
        elseif i==L-1
            dpsidtheta[L-1] = (-p[2]+8*p[1]-8*p[L-2]+p[L-3])/(12.0*dtheta)
        elseif i==L-2
            dpsidtheta[L-2] = (-p[1]+8*p[L-1]-8*p[L-3]+p[L-4])/(12.0*dtheta)
        elseif i==L
            dpsidtheta[L] = dpsidtheta[1]
        else
            dpsidtheta[i] = (-p[i+2]+8*p[i+1]-8*p[i-1]+p[i-2])/(12.0*dtheta)
        end
    end
    return dpsidtheta
end 
    
function dpsi_dr(arr::Array{Float64,3})
    r_deriv = Float64[] 
    for theta in range(1,stop=size(arr,2))
        for t in range(1,stop=size(arr,3))
            r_deriv = append!(r_deriv,dAdr(arr,theta,t))
        end
    end
     # In Julia the convention in reshape is the opposite of Python
    return permutedims(reshape(r_deriv,(size(arr,1),size(arr,3),size(arr,2))),(1,3,2))
end

function dpsi_dtheta(arr::Array{Float64,3})
    t_deriv = Float64[] 
    for r in range(1,stop=size(arr,1))
        for t in range(1,stop=size(arr,3))
            t_deriv = append!(t_deriv,dAdt(arr,r,t))
        end
    end
    return permutedims(reshape(t_deriv,(size(arr,2),size(arr,3),size(arr,1))),(3,1,2))
end

function quad_r_int(B_eint,xmin,xmax,yo,to)
    res, err = quadgk(x -> B_eint(x,yo,to),xmin,xmax)
    return res
end

function quad_r_int(B_eint,xmin,xmax,yo,zo,to)
    res, err = quadgk(x -> B_eint(x,yo,zo,to),xmin,xmax)
    return res
end

function quad_u_int(B_eint,ymin,ymax,xo,zo,to)
    res, err = quadgk(y -> B_eint(xo,y,zo,to),ymin,ymax)
    return res
end

function quad_u_int(B_eint,ymin,ymax,xo,to)
    res, err = quadgk(y -> B_eint(xo,y,to),ymin,ymax)
    return res
end

function Au(B3_eint,ro,uo,rMA,uMA,to)
    A = quad_r_int(B3_eint,rMA,ro,uo,to)
    return A
end

function Au(B3_eint,ro,uo,rMA,uMA,zo,to)
    A = quad_r_int(B3_eint,rMA,ro,uo,zo,to)
    return A
end

function Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA,to)
    A = quad_u_int(B1_eint,uMA,uo,ro,to) - quad_r_int(B2_eint,rMA,ro,uMA,to)
    return A
end

function Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA,zo,to)
    A = quad_u_int(B1_eint,uMA,uo,ro,zo,to) - quad_r_int(B2_eint,rMA,ro,uMA,zo,to)
    return A
end

function grid_Au(B3_eint,Ndims::Int)
    if Ndims == 3
        Au_arr = []
        for ro in rn
            for uo in un
                for fo in phin
                    for to in tn
                        rMA = rmaxis[convert(Int,to)+1]
                        uMA = umaxis[convert(Int,to)+1]
                        append!(Au_arr,Au(B3_eint,ro,uo,rMA,uMA,fo,to))
                    end
                end
            end
        end
        Au_arr = permutedims(reshape(Au_arr, size(tn)[1], size(phin)[1], size(un)[1],size(rn)[1]),(4,3,2,1))
    elseif Ndims == 2
        Au_arr = []
        for ro in rn
            for uo in un
                for to in tn
                    rMA = rmaxis[convert(Int,to)+1]
                    uMA = umaxis[convert(Int,to)+1]
                    append!(Au_arr,Au(B3_eint,ro,uo,rMA,uMA,to))
                end
            end
        end
        Au_arr = permutedims(reshape(Au_arr, size(tn)[1], size(un)[1],size(rn)[1]),(3,2,1))
    end
    return Au_arr
end

function grid_Aphi(B1_eint,B2_eint,Ndims::Int)
    if Ndims == 3
        Aphi_arr = []
        for ro in rn
            for uo in un
                for fo in phin
                    for to in tn
                        rMA = rmaxis[convert(Int,to)+1]
                        uMA = umaxis[convert(Int,to)+1]
                        append!(Aphi_arr,Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA,fo,to))
                    end
                end
            end
        end
        Aphi_arr = permutedims(reshape(Aphi_arr, size(tn)[1],size(phin)[1],size(un)[1],size(rn)[1]),(4,3,2,1))
    elseif Ndims == 2
        Aphi_arr = []
        for ro in rn
            for uo in un
                for to in tn
                    rMA = rmaxis[convert(Int,to)+1]
                    uMA = umaxis[convert(Int,to)+1]
                    append!(Aphi_arr,Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA,to))
                end
            end
        end
        Aphi_arr = permutedims(reshape(Aphi_arr, size(tn)[1],size(un)[1],size(rn)[1]),(3,2,1))
    end
    return Aphi_arr
end

function b1_an(GAp_eint, Ndims::Int)
    if Ndims == 2
        b1_an = []
        for r in rn
            for u in un
                for t in tn
                    append!(b1_an, Interpolations.gradient(GAp_eint,r,u,t)[2])
                end
            end
        end
        b1_an = permutedims(reshape(b1_an,size(tn)[1],size(un)[1],size(rn)[1]),(3,2,1))
    elseif Ndims == 3
        b1_an = []
        for r in rn
            for u in un
                for f in phin
                    for t in tn
                        append!(b1_an, Interpolations.gradient(GAp_eint,r,u,f,t)[2])
                    end
                end
            end
        end
        b1_an = permutedims(reshape(b1_an,size(tn)[1],size(phin)[1],size(un)[1],size(rn)[1]),(4,3,2,1))
    end
    return b1_an
end
    
function b2_an(GAp_eint,Ndims::Int)
    if Ndims == 2
        b2_an = []
        for r in rn
            for u in un
                for t in tn
                    append!(b2_an, -Interpolations.gradient(GAp_eint,r,u,t)[1])
                end
            end
        end
        b2_an = permutedims(reshape(b2_an,size(tn)[1],size(un)[1],size(rn)[1]),(3,2,1))
    elseif Ndims == 3 
        b2_an = []
        for r in rn
            for u in un
                for f in phin
                    for t in tn
                        append!(b2_an, -Interpolations.gradient(GAp_eint,r,u,f,t)[1])
                    end
                end
            end
        end
        b2_an = permutedims(reshape(b2_an,size(tn)[1],size(phin)[1],size(un)[1],size(rn)[1]),(4,3,2,1))
    end
    return b2_an
end 

function b3_an(GAu_eint,Ndims::Int)
    if Ndims == 2
        b3_an = []
        for r in rn
            for u in un
                for t in tn
                    append!(b3_an, Interpolations.gradient(GAu_eint,r,u,t)[1] -Interpolations.gradient(GAu_eint,r,u,t)[2])
                end
            end
        end
        b3_an = permutedims(reshape(b3_an,size(tn)[1],size(un)[1],size(rn)[1]),(3,2,1))
    elseif Ndims == 3
        b3_an = []
        for r in rn
            for u in un
                for f in phin
                    for t in tn
                        append!(b3_an, Interpolations.gradient(GAu_eint,r,u,f,t)[1] -Interpolations.gradient(GAu_eint,r,u,f,t)[2])
                    end
                end
            end
        end
        b3_an = permutedims(reshape(b3_an,size(tn)[1],size(phin)[1],size(un)[1],size(rn)[1]),(4,3,2,1))
    end
    return b3_an
end

# Definition of flux surface system of equations
function Flux_surface!(du,u,p,t)
    du[1] = B1t_eint(u[1],u[2],p[1])/JBpN_eint(u[1],u[2],p[1]) # use of n=o components          
    du[2] = B2t_eint(u[1],u[2],p[1])/JBpN_eint(u[1],u[2],p[1])
    du[3] = B3t_eint(u[1],u[2],p[1])/JBpN_eint(u[1],u[2],p[1])
    end

# callback for crossing 2pi and finish integration
function condition(u,t,integrator)
    abs(u[2]) - 2*pi 
end

# callback for crossing zero for flux surfaces close to shifted magnetic axis
function condition_zero(u,t,integrator)
    u[2] - 0.0
end

# check if zero has been crossed twice
function condition_crossing(u,t,integrator)
    integrator.p[2] == 2.0
end 

function affect_cross!(integrator)
    integrator.p[2] = integrator.p[2]+1
end


function affect!(integrator)
    terminate!(integrator)
end

cb1 = ContinuousCallback(condition,affect!,rootfind = true)
cb2 = ContinuousCallback(condition_zero,affect_cross!,rootfind = true)
cb3 = DiscreteCallback(condition_crossing,affect!)

# Integrator
function fs_integration(rs::Float64,us::Float64,time::Int)
    u0 = [rs,us,0.0]
    p = [time,0]
    tspan = (0.0,180.0)
    prob = ODEProblem(Flux_surface!,u0,tspan,p)
    cbs = CallbackSet(cb1,cb2,cb3)
    integrator = init(prob,Vern9(),callback=cbs,reltol=1.e-10,abstol=1.e-10)
    sol = solve(prob,Vern9(),callback=cbs,reltol=1.e-10,abstol=1.e-10)
    return sol
end

function q_profile()
    q = []
    for t in range(0, stop=tdim-1)
        for r in r_of_psi_array[t+1,psi_min_ind:psi_max_ind]
            sol = fs_integration(r,0.0,t)
            append!(q, last(sol[3,length(sol.t)])/(2.0*pi))
        end
    end
    println("q finished")
    q_fun = reshape(q,psi_range,tdim)
    return q_fun
end

function resize_q(arr::Array{Float64,2})
    big_arr = []
    for psi_ind in 1:psidim
        for uf_ind in 1:ufdim
            for t_ind in 1:tdim
                for fi_ind in 1:fidim
                    append!(big_arr, arr[psi_ind,t_ind])
                end
            end
        end
    end
    
    aug_arr = permutedims(reshape(big_arr,fidim,tdim,ufdim,psidim),(4,3,1,2))
    return aug_arr
end

function stretch_q(q_array::Array{Float64,2})
    psin_list = LinRange(0.0,1.0,psidim)
    small_grid = LinRange(psin_list[psi_min_ind],psin_list[psi_max_ind],psi_range)
    q_str = []
    for t in 1:size(q_array)[2]
        q = Float64.(q_array[:,t])
        q_int = Interpolations.interpolate(q,(BSpline(Cubic(Line(OnGrid())))))
        q_sint = scale(q_int,small_grid)
        q_eint = extrapolate(q_sint,Line())
        for i in psin_list
            append!(q_str,q_eint(i))
        end
    end
    q_str = reshape(q_str,(length(psin_list),size(q_array,2)))    
    return q_str
end

# Creates dictionary of grid in the form {t,psi,theta_mag} -> (r,theta)
function new_grid()
    Grid_Dict = Dict{Tuple{Int64,Float64,Float64},Tuple{Float64,Float64}}()
    for t in range(0, stop=tdim-1)
        for (i,r) in enumerate(r_of_psi_array[t+1,psi_min_ind:psi_max_ind])
            sol = fs_integration(r,0.0,t)
            uf_sol = sol[3,:]
            q_redef = last(sol[3,length(sol.t)])/(2.0*pi)
            uf_rsc = uf_sol/q_redef
            l = sol.t
            l_of_uf_rsc = Spline1D(sort!(uf_rsc[1:end-1]),l[1:end-1],bc="extrapolate",s=1e-4)
            ls = l_of_uf_rsc(ufn)
            polar_coords = sol(ls,idxs=1:2)
            psi_val = pn[i+psi_min_ind-1]
            for j in 1:udim
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
                        ArrVal = arr_eint(r,u,phi,t)
                        append!(ArrP,ArrVal)
                    else
                        r = GD[Int64(round(t)),round(p,digits=2),Float64(uf)][1]
                        u = GD[Int64(round(t)),round(p,digits=2),Float64(uf)][2]
                        ArrVal = arr_eint(r,u,phi,t)
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
                    ArrVal = arr_eint(r,u,t)
                    append!(ArrP,ArrVal)
                else
                    r = GD[Int64(round(t)),round(p,digits=2),Float64(uf)][1]
                    u = GD[Int64(round(t)),round(p,digits=2),Float64(uf)][2]
                    ArrVal = arr_eint(r,u,t)
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


########################################################################################################################################################################################################################################################################################################################################## EXECUTION ############################################################################
######################################################################################################################################################################################################################################################################################
na = CartesianIndex();

# Simulation loading
time_interval = tend-tstart; 
chunk = floor(time_interval/size_comm); # loading chunk
rank_tstart = convert(Int,tstart+rank*chunk); # Starting and ending points of each proccess
rank_tend = convert(Int,tstart+(rank+1)*chunk);
final_chunk_size = convert(Int,tend-(tstart+(size_comm-1)*chunk));


pxr = pyimport("pixie_read_st")
pxr.pixieload(filepath * "pixie3d.h5")

# Loading of data
if rank == size_comm-1
    psi = pxr.load_array(3,4,rank_tstart,tend);
    B1 = pxr.load_array(1,0,rank_tstart,tend); # Contravariant components
    B2 = pxr.load_array(1,1,rank_tstart,tend);
    B3 = pxr.load_array(1,2,rank_tstart,tend);
    B_1 = pxr.load_array(2,0,rank_tstart,tend); # Covariant components
    B_2 = pxr.load_array(2,1,rank_tstart,tend);
else
    psi = pxr.load_array(3,4,rank_tstart,rank_tend);
    B1 = pxr.load_array(1,0,rank_tstart,rank_tend);
    B2 = pxr.load_array(1,1,rank_tstart,rank_tend);
    B3 = pxr.load_array(1,2,rank_tstart,rank_tend);
    B_1 = pxr.load_array(2,0,rank_tstart,rank_tend);
    B_2 = pxr.load_array(2,1,rank_tstart,rank_tend);
end

#psi=-psi # Make psi be monotonically increasing. Inspect if needed. Fix last sign of Q accordingly. Needed for sawtooth.

psit = dropdims(mean(psi,dims=3),dims=3);
B1tc = dropdims(mean(B1,dims=3),dims=3);
B2tc = dropdims(mean(B2,dims=3),dims=3);
B_1tc = dropdims(mean(B_1,dims=3),dims=3);
B_2tc = dropdims(mean(B_2,dims=3),dims=3);

# Calculate JBp^2 term with cell components 
JBpsqc = (B1tc.*B_1tc) .+ (B2tc.*B_2tc);
JBpc = sqrt.((B1tc.*B_1tc) .+ (B2tc.*B_2tc));


# Make B-arrays into node-based quantities
# definitions of cell grid
num_r_cells = size(B3)[1];
num_u_cells = size(B3)[2];
num_phi_cells = size(B3)[3];
dn_r = (1.0/num_r_cells);
dn_u = ((2.0*pi)/num_u_cells);

# Cell-based grid
rc = LinRange(0.0+(dn_r/2.0),1.0-(dn_r/2.0),num_r_cells);
uc = LinRange(0.0+(dn_u/2.0),2.0*pi-(dn_u/2.0),num_u_cells);
phic = LinRange(0.0+(dn_u/2.0),2.0*pi-(dn_u/2.0),num_phi_cells);
tn = LinRange(0, size(B3)[4]-1,size(B3)[4]);

# Node-based grid
rn = LinRange(0.0,1.0,(num_r_cells+1));
un = LinRange(0.0,2.0*pi,(num_u_cells+1));
phin = LinRange(0.0,2.0*pi,(num_phi_cells+1));

# Node-based grid dimensions
rdim = size(psi)[1];
udim = size(psi)[2];
fidim = size(psi)[3];
tdim = size(psi)[4];

# Psi-magnetic angle grid dimensions
psidim = 101;
ufdim = udim; # ufdim same as udim to avoid aliasing

# Straight field line coordinate grid
pn = LinRange(0.0,1.0,psidim);
ufn = LinRange(0.0,2.0*pi,ufdim);

# metric tensor arrays in logical coordinates
#grr_l = grr_eint(rn,un);
#gtt_l = gtt_eint(rn,un);

# Interpolate on cell-based grid
B1_int_cell = Interpolations.interpolate(B1,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Line(OnGrid())))));
B2_int_cell = Interpolations.interpolate(B2,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Line(OnGrid())))));
B3_int_cell = Interpolations.interpolate(B3,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Line(OnGrid())))));
B_1_int_cell = Interpolations.interpolate(B_1,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Line(OnGrid())))));
B_2_int_cell = Interpolations.interpolate(B_2,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Line(OnGrid())))));
JBpsq_int_cell = Interpolations.interpolate(JBpsqc,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Line(OnGrid())))));
JBpc_int_cell = Interpolations.interpolate(JBpc,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Line(OnGrid())))));



B1_sint_cell = scale(B1_int_cell,rc,uc,phic,tn);
B2_sint_cell = scale(B2_int_cell,rc,uc,phic,tn);
B3_sint_cell = scale(B3_int_cell,rc,uc,phic,tn);
B_1_sint_cell = scale(B_1_int_cell,rc,uc,phic,tn);
B_2_sint_cell = scale(B_2_int_cell,rc,uc,phic,tn);
JBpsq_sint_cell = scale(JBpsq_int_cell,rc,uc,tn);
JBpc_sint_cell = scale(JBpc_int_cell,rc,uc,tn);

B1_eint_cell = extrapolate(B1_sint_cell, (Line(),Periodic(),Periodic(),Line()));
B2_eint_cell = extrapolate(B2_sint_cell, (Line(),Periodic(),Periodic(),Line()));
B3_eint_cell = extrapolate(B3_sint_cell, (Line(),Periodic(),Periodic(),Line()));
B_1_eint_cell = extrapolate(B_1_sint_cell, (Line(),Periodic(),Periodic(),Line()));
B_2_eint_cell = extrapolate(B_2_sint_cell, (Line(),Periodic(),Periodic(),Line()));
JBpsq_eint_cell = extrapolate(JBpsq_sint_cell, (Line(),Periodic(),Line()));
JBpc_eint_cell = extrapolate(JBpc_sint_cell, (Line(),Periodic(),Line()));

# Evaluate B on node grid
B1 = B1_eint_cell(rn,un,phin,tn);
B2 = B2_eint_cell(rn,un,phin,tn);
B3 = B3_eint_cell(rn,un,phin,tn);
B_1 = B_1_eint_cell(rn,un,phin,tn);
B_2 = B_2_eint_cell(rn,un,phin,tn);
JBpsqN = JBpsq_eint_cell(rn,un,tn);
JBpN = JBpc_eint_cell(rn,un,tn);

# Here do the averaged fields corresponding to the node-projected ones-otherwise, the removal of the n=0 component from the perturbation is not exact
B1t = dropdims(mean(B1,dims=3),dims=3);
B2t = dropdims(mean(B2,dims=3),dims=3);
B3t = dropdims(mean(B3,dims=3),dims=3);
B_1t = dropdims(mean(B_1,dims=3),dims=3);
B_2t = dropdims(mean(B_2,dims=3),dims=3);


MPI.Barrier(comm)

# Interpolating on grid nodes
#psi_int = Interpolations.interpolate(psi,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B3_int = Interpolations.interpolate(B3,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B1_int = Interpolations.interpolate(B1,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B2_int = Interpolations.interpolate(B2,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B3t_int = Interpolations.interpolate(B3t,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B1t_int = Interpolations.interpolate(B1t,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B2t_int = Interpolations.interpolate(B2t,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B_1_int = Interpolations.interpolate(B_1,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B_2_int = Interpolations.interpolate(B_2,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B_1t_int = Interpolations.interpolate(B_1t,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B_2t_int = Interpolations.interpolate(B_2t,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
JBpsqN_int = Interpolations.interpolate(JBpsqN,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
JBpN_int = Interpolations.interpolate(JBpN,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));

# Splining Cartesian Maps
X_int = Interpolations.interpolate(pxr.X,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
Z_int = Interpolations.interpolate(pxr.Z,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));


# Rescaling interpolation intervals
#psi_sint = scale(psi_int, rn,un,phin,tn);
B3_sint = scale(B3_int, rn,un,phin,tn);
B1_sint = scale(B1_int, rn,un,phin,tn);
B2_sint = scale(B2_int, rn,un,phin,tn);
B3t_sint = scale(B3t_int, rn,un,tn);
B1t_sint = scale(B1t_int, rn,un,tn);
B2t_sint = scale(B2t_int, rn,un,tn);
B_1_sint = scale(B_1_int, rn,un,phin,tn);
B_2_sint = scale(B_2_int, rn,un,phin,tn);
B_1t_sint = scale(B_1t_int, rn,un,tn);
B_2t_sint = scale(B_2t_int, rn,un,tn);
JBpsqN_sint = scale(JBpsqN_int, rn,un,tn);
JBpN_sint = scale(JBpN_int, rn,un,tn);


X_sint = scale(X_int, rn,un,phin);
Z_sint = scale(Z_int, rn,un,phin);


# Extrapolations of rescaled functions
#psi_eint = extrapolate(psi_sint, (Line(),Periodic(),Periodic(),Line()));
B3_eint = extrapolate(B3_sint, (Line(),Periodic(),Periodic(),Line()));
B1_eint = extrapolate(B1_sint, (Line(),Periodic(),Periodic(),Line()));
B2_eint = extrapolate(B2_sint, (Line(),Periodic(),Periodic(),Line()));
B3t_eint = extrapolate(B3t_sint, (Line(),Periodic(),Line()));
B1t_eint = extrapolate(B1t_sint, (Line(),Periodic(),Line()));
B2t_eint = extrapolate(B2t_sint, (Line(),Periodic(),Line()));
B_1_eint = extrapolate(B_1_sint, (Line(),Periodic(),Periodic(),Line()));
B_2_eint = extrapolate(B_2_sint, (Line(),Periodic(),Periodic(),Line()));
B_1t_eint = extrapolate(B_1t_sint, (Line(),Periodic(),Line()));
B_2t_eint = extrapolate(B_2t_sint, (Line(),Periodic(),Line()));
JBpsqN_eint = extrapolate(JBpsqN_sint, (Line(),Periodic(),Line()));
JBpN_eint = extrapolate(JBpN_sint, (Line(),Periodic(),Line()));

X_eint = extrapolate(X_sint, (Line(),Periodic(),Periodic()));
Z_eint = extrapolate(Z_sint, (Line(),Periodic(),Periodic()));


# Preparations of the Python module

pxr.Calculation_of_Units_and_Sizes()
pxr.Axes_of_Interpolation(B3)
psi_min,norm = pxr.Normalization_numbers(psit,JBpsqN,x0);
pythonresult = pxr.create_r_psi_list(psit,JBpsqN,x0);
r_of_psi_array = pythonresult[1]; # Pick python outputs
rmaxis = pythonresult[2];
umaxis = pythonresult[3];

q_tmp = q_profile()
q = Float64.(q_tmp) # 2D q array
q_str_tmp = stretch_q(q)
q_str = Float64.(q_str_tmp)
q_rsz = resize_q(q_str)
Q = Float64.(q_rsz) # 4D q array

MPI.Barrier(comm)
if rank != 0
    MPI.Send(q_str,0,rank,comm)
elseif rank == 0
    q_dim1 = size(q_str)[1]
    q_dim2 = size(q_str)[2]
    q_composite = Array{Float64,2}(undef,q_dim1,0)
    q_composite = cat(dims=2,q_composite,q_str)
    for r in range(1,stop=size_comm-2)
        q_recv = Array{Float64,2}(undef,q_dim1,q_dim2)
        MPI.Recv!(q_recv,r,r,comm)
        global q_composite = cat(dims=2,q_composite,q_recv)
    end
    q_recv_fin = Array{Float64,2}(undef,q_dim1,final_chunk_size)
    MPI.Recv!(q_recv_fin,size_comm-1,size_comm-1,comm)
    q_composite = cat(dims=2,q_composite,q_recv_fin)
    
    npzwrite(filepath * q_file ,q_composite)
end
MPI.Barrier(comm)
   
# Grid Dictionary creation
GD = new_grid()

na = [CartesianIndex()];

# dB calculation
dB1 = B1.-B1t[:,:,na,:];
dB2 = B2.-B2t[:,:,na,:];    
dB3 = B3.-B3t[:,:,na,:];

# b-hat-rho in r,theta,phi
brho = dB1.*B2t[:,:,na,:] .- dB2.*B1t[:,:,na,:];
b_hat_rho = brho./B3t[:,:,na,:];
db3 = dB3./B3t[:,:,na,:]; # scaling db3 - first part

# Interpolate b_hat_rho avoiding problematic r=0 point
b_hat_rho_int = Interpolations.interpolate(b_hat_rho[2:end,:,:,:],(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
b_hat_rho_sint = scale(b_hat_rho_int, rn[2:end],un,phin,tn);
b_hat_rho_eint = extrapolate(b_hat_rho_sint, (Line(),Periodic(),Periodic(),Line()));

# Interpolate dB3 
db3_int = Interpolations.interpolate(db3[2:end,:,:,:],(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
db3_sint = scale(db3_int, rn[2:end],un,phin,tn);
db3_eint = extrapolate(db3_sint, (Line(),Periodic(),Periodic(),Line()));

# Project in new coordinates 
b_hat_rho_P = project(b_hat_rho_eint);
db3_P = project(db3_eint);
println("Finished projection")

# Project Cartesian Maps
X_P = project_Cartesian_Map(X_eint);
Z_P = project_Cartesian_Map(Z_eint);
B3t_P = project_toroidally_symmetric(B3t_eint);

# Interpolate b_hat_rho_P 
b_hat_rho_P_int = Interpolations.interpolate(Float64.(b_hat_rho_P),(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
b_hat_rho_P_sint = scale(b_hat_rho_P_int, pn[1:psi_max_ind],ufn,phin,tn);
b_hat_rho_P_eint = extrapolate(b_hat_rho_P_sint, (Line(),Periodic(),Periodic(),Line()));

# Interpolate dB3_P 
db3_P_int = Interpolations.interpolate(Float64.(db3_P),(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
db3_P_sint = scale(db3_P_int, pn[1:psi_max_ind],ufn,phin,tn);
db3_P_eint = extrapolate(db3_P_sint, (Line(),Periodic(),Periodic(),Line()));


# Extrapolate in whole psi range and multiply with q(psi,uf,phi,t)
db3_P_ext = db3_P_eint(pn,ufn,phin,tn);
b_hat_rho_P_ext = b_hat_rho_P_eint(pn,ufn,phin,tn);
b_hat_rho_ncq = Q.*b_hat_rho_P_ext; # Watch out for this sign, depends on original psi monotonicity.
db3_hat_P_ext = Q.*db3_P_ext; # scaling db3 - second part.

# Interpolate Cartesian Maps
X_P_int = Interpolations.interpolate(Float64.(X_P), (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
Z_P_int = Interpolations.interpolate(Float64.(Z_P), (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));

X_P_sint = scale(X_P_int, pn[1:psi_max_ind],ufn,phin,tn);
Z_P_sint = scale(Z_P_int, pn[1:psi_max_ind],ufn,phin,tn);

X_P_eint = extrapolate(X_P_sint, (Line(),Periodic(),Periodic(),Line()));
Z_P_eint = extrapolate(Z_P_sint, (Line(),Periodic(),Periodic(),Line()));

# Extrapolate in whole psi range
X_P_ext = X_P_eint(pn,ufn,phin,tn);
Z_P_ext = Z_P_eint(pn,ufn,phin,tn);

# Interpolate B3t
B3t_P_int = Interpolations.interpolate(Float64.(B3t_P),(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
B3t_P_sint = scale(B3t_P_int, pn[1:psi_max_ind],ufn,tn);
B3t_P_eint = extrapolate(B3t_P_sint, (Line(),Periodic(),Line()));
B3t_P_ext = B3t_P_eint(pn,ufn,tn);

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
    
    bhatrho = Float64.(bhat_composite)
    Xproj = Float64.(X_composite)
    Zproj = Float64.(Z_composite)
    db3proj = Float64.(db3_composite)
    B3tproj = Float64.(B3t_composite)
    
    npzwrite(filepath * X_file, Xproj)
    npzwrite(filepath * Z_file, Zproj)
    npzwrite(filepath * db3_file, db3proj)
    npzwrite(filepath * B3t_file, B3tproj)
    npzwrite(filepath * b_hat_rho_file, bhatrho)
end    
    


###################################################################################################################################################################################################################################################################################################################################################################################################

#function Flux_surface!(du,u,p,t)
#    alpha = 2.24
    #du[1] = -dpsidthetat_eint(u[1],u[2],p[1])/(sqrt((dpsidthetat_eint(u[1],u[2],p[1]))^2+((u[1])^2)*(dpsidrt_eint(u[1],u[2],p[1]))^2))
    #du[2] = dpsidrt_eint(u[1],u[2],p[1])/(sqrt((dpsidthetat_eint(u[1],u[2],p[1]))^2+((u[1])^2)*(dpsidrt_eint(u[1],u[2],p[1]))^2))
    #du[3] = -B3t_eint(u[1],u[2],p[1])/(sqrt((dpsidthetat_eint(u[1],u[2],p[1]))^2+(u[1]^2)*(dpsidrt_eint(u[1],u[2],p[1]))^2))
#    du[1] = -B1t_eint(u[1],u[2],p[1])/(sqrt(abs(grr_eint(u[1]*alpha,u[2])*((B1t_eint(u[1],u[2],p[1]))^2) + 2*grt_eint(u[1]*alpha,u[2])*B1t_eint(u[1],u[2],p[1])*B2t_eint(u[1],u[2],p[1]) + gtt_eint(u[1]*alpha,u[2])*((B2t_eint(u[1],u[2],p[1]))^2))))           
#    du[2] = -B2t_eint(u[1],u[2],p[1])/(sqrt(abs(grr_eint(u[1]*alpha,u[2])*((B1t_eint(u[1],u[2],p[1]))^2) + 2*grt_eint(u[1]*alpha,u[2])*B1t_eint(u[1],u[2],p[1])*B2t_eint(u[1],u[2],p[1]) + gtt_eint(u[1]*alpha,u[2])*((B2t_eint(u[1],u[2],p[1]))^2))))
#    du[3] = B3t_eint(u[1],u[2],p[1])/(sqrt(abs(grr_eint(u[1]*alpha,u[2])*((B1t_eint(u[1],u[2],p[1]))^2) + 2*grt_eint(u[1]*alpha,u[2])*B1t_eint(u[1],u[2],p[1])*B2t_eint(u[1],u[2],p[1]) + gtt_eint(u[1]*alpha,u[2])*((B2t_eint(u[1],u[2],p[1]))^2))))
#    end

# read metric
#metric = npzread(metricpath);

# load covariant components
#grr = metric["grr_do"];
#gtt = metric["gtt_do"];
#grt = metric["grt_do"];
#Jac = metric["Jacobian"];
#radial_grid = metric["radial_grid"];
#pol_grid = metric["pol_grid"];

# radial and poloidal grid for metric interpolation
#rm = LinRange(radial_grid[1],radial_grid[end],size(radial_grid)[1]);
#um = LinRange(pol_grid[1],pol_grid[end],size(pol_grid)[1]);

# metric interpolation
#grr_int = Interpolations.interpolate(grr, (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
#gtt_int = Interpolations.interpolate(gtt, (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
#grt_int = Interpolations.interpolate(grt, (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
#Jac_int = Interpolations.interpolate(Jac, (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));

#grr_sint = scale(grr_int,rm,um);
#gtt_sint = scale(gtt_int,rm,um);
#grt_sint = scale(grt_int,rm,um);
#Jac_sint = scale(Jac_int,rm,um);

#grr_eint = extrapolate(grr_sint,(Line(),Line()));
#gtt_eint = extrapolate(gtt_sint,(Line(),Line()));
#grt_eint = extrapolate(grt_sint,(Line(),Line()));
#Jac_eint = extrapolate(Jac_sint,(Line(),Line()));



# Divergence Cleanup
# Vector potential 
#GAu = grid_Au(B3_eint,3);
#GAp = grid_Aphi(B1_eint,B2_eint,3);

#GAu = Float64.(GAu);
#GAp = Float64.(GAp);

# Spline the vector potential
#GAu_int = Interpolations.interpolate(GAu,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
#GAp_int = Interpolations.interpolate(GAp,(BSpline(Quadratic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));

#GAu_sint = scale(GAu_int,rn,un,phin,tn);
#GAp_sint = scale(GAp_int,rn,un,phin,tn);

#GAu_eint = extrapolate(GAu_sint, (Line(),Periodic(),Periodic(),Line()));
#GAp_eint = extrapolate(GAp_sint, (Line(),Periodic(),Periodic(),Line()));

# Create divergence free magnetic field
#B1df = b1_an(GAp_eint,3);
#B2df = b2_an(GAp_eint,3);
#B3df = b3_an(GAu_eint,3);

#B1df = Float64.(B1df);
#B2df = Float64.(B2df);
#B3df = Float64.(B3df);

# JBp with divergence free field
#JBp = sqrt.((B1df.*grr_l[:,:,na,na].*B1df) .+ (B2df.*gtt_l[:,:,na,na].*B2df));

# Toroidally average
#B1dft = dropdims(mean(B1df,dims=3),dims=3);
#B2dft = dropdims(mean(B2df,dims=3),dims=3);
#B3dft = dropdims(mean(B3df,dims=3),dims=3);
#JBpt = dropdims(mean(JBp,dims=3),dims=3);

# Spline divergence free magnetic field
#B1df_int = Interpolations.interpolate(B1df,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
#B2df_int = Interpolations.interpolate(B2df,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
#B3df_int = Interpolations.interpolate(B3df,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
#B1dft_int = Interpolations.interpolate(B1dft,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
#B2dft_int = Interpolations.interpolate(B2dft,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
#B3dft_int = Interpolations.interpolate(B3dft,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
#JBpt_int = Interpolations.interpolate(JBpt,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));


# Scale interpolations
#B1df_sint = scale(B1df_int,rn,un,phin,tn);
#B2df_sint = scale(B2df_int,rn,un,phin,tn);
#B3df_sint = scale(B3df_int,rn,un,phin,tn);
#B1dft_sint = scale(B1dft_int,rn,un,tn);
#B2dft_sint = scale(B2dft_int,rn,un,tn);
#B3dft_sint = scale(B3dft_int,rn,un,tn);
#JBpt_sint = scale(JBpt_int,rn,un,tn);

# Extrapolate interpolations
#B1df_eint = extrapolate(B1df_sint, (Line(),Periodic(),Periodic(),Line()));
#B2df_eint = extrapolate(B2df_sint, (Line(),Periodic(),Periodic(),Line()));
#B3df_eint = extrapolate(B3df_sint, (Line(),Periodic(),Periodic(),Line()));
#B1dft_eint = extrapolate(B1dft_sint, (Line(),Periodic(),Line()));
#B2dft_eint = extrapolate(B2dft_sint, (Line(),Periodic(),Line()));
#B3dft_eint = extrapolate(B3dft_sint, (Line(),Periodic(),Line()));
#JBpt_eint = extrapolate(JBpt_sint, (Line(),Periodic(),Line()));

