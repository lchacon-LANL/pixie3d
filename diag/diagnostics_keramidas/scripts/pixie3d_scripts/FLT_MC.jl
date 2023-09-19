using Interpolations
using DifferentialEquations
using PyCall
using NPZ
using Statistics
using QuadGK
using Plots
using LaTeXStrings
using Printf
using FFTW

#filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_tear/dt_sh_m3_n2.scratch/"
#filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_new/dt2.scratch./"
filepath = "/net/scratch4/giannis_kx/pixie3d/iter/int_kink/11/11_visc_old_nodiff/11_visc_old_nodiff.scratch/"

timestep = 80;
timestep = 21;
Npsi = 100; # Number of initial conditions
Nuf = 3;
Lmax = 20000;
savestring = @sprintf("PP_MC_11_visc_old_nodiff@%s.png",(timestep))

# Loading fields
bpsi = npzread(filepath *  "b_hat_rho_mpi.npy")
bpsi = 2*pi*bpsi # Rescale the perturbation by 2pi
q = npzread(filepath * "q_mpi.npy")
Nm = size(bpsi)[2]; # number of m-modes
Nn = size(bpsi)[3]; # number of n-modes

function Field_line!(du,u,p,t)
    du[1] = bpsi_eint(u[1],u[2],u[3],p[1])
    du[2] = 1 + buf_eint(u[1],u[2],u[3])
    du[3] = q_eint(u[1],p[1])
end

function condition_cross(u,t,integrator)
    u[3] - integrator.p[2]*2.0*pi
end
function count_cross!(integrator)
    integrator.p[2] = integrator.p[2]+1
end

cb1 = ContinuousCallback(condition_cross,count_cross!,rootfind = true,save_positions=(false,true))

function field_line_integration(psis::Float64,ufs::Float64,time::Int)
    u0 = [psis,ufs,0.0]
    p = [time,1]
    tspan = (0.0,Lmax)
    prob = ODEProblem(Field_line!,u0,tspan,p)
    sol = solve(prob,BS5(),callback = cb1, reltol=1.e-10,abstol=1.e-10,save_everystep=false,save_start=false, save_end=false)
    return sol
end

function Poincare_Map(Npsi::Int,Nuf::Int,timestep::Int)
    icpsi = LinRange(0.01,0.9,Npsi)
    icuf = LinRange(0.0,2*pi,Nuf)
    ic_tot = Npsi*Nuf

    psi = [[] for i in 1:ic_tot]
    theta_f = [[] for i in 1:ic_tot]
    phi = [[] for i in 1:ic_tot]
    for i in 1:Npsi
        for j in 1:Nuf
            sol = field_line_integration(icpsi[i],icuf[j],timestep)
            psi[i*j] = sol[1,:]
            theta_f[i*j] = sol[2,:]
            phi[i*j] = sol[3,:]
        end
    end
    return psi, theta_f
end

function n_mode_num(n)
    if n == 1
        mode = 0
    end
    if n <= (Nn-1)/2 && n > 1
        mode =  n-1
    end
    if n > (Nn-1)/2
        mode = n-(Nn+1)
    end
    return mode
end

function m_mode_num(m)
    if m == 1
        mode = 0
    end
    if m <= (Nm-1)/2 && m > 1
        mode = -(m-1) # m's have flipped frequencies
    end
    if m > (Nm-1)/2
        mode = Nm-m+1
    end
    return mode
end



############################################################################################################################################################################################################################################################################################################# EXECUTION ########################################################################
################################################################################################################################################################################################################################################################

# Import the eqdsk info
eqdsk = pyimport("poinc") 
eqdsk.eqdsk_info()


# Node-based grid
psin = LinRange(0.0,1.0,size(bpsi)[1]);
ufn = LinRange(0.0,2.0*pi,size(bpsi)[2]);
phin = LinRange(0.0,2.0*pi,size(bpsi)[3]);
tn = LinRange(0, size(bpsi)[4]-1,size(bpsi)[4]);

# Spline the field
bpsi_int = Interpolations.interpolate(bpsi,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));
q_int = Interpolations.interpolate(q,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));

bpsi_sint = scale(bpsi_int,psin,ufn,phin,tn);
q_sint = scale(q_int,psin,tn);

bpsi_eint = extrapolate(bpsi_sint, (Line(),Periodic(),Periodic(),Line()));
q_eint = extrapolate(q_sint, (Line(),Line()));

# Calculate perturbation along magnetic angle (uf) with Fourier components
# Fourier transform to get bpsi_mn
bpsi_t = bpsi[:,:,:,timestep];
bpsi_mn = fft(bpsi_t,[2,3]);

# Create lists to scale the interpolations for bpsi_mn
psin_list = LinRange(0,1.0,size(bpsi_mn)[1]);
m_index_list = LinRange(1,size(bpsi_mn)[2],size(bpsi_mn)[2]);
n_index_list = LinRange(1,size(bpsi_mn)[3],size(bpsi_mn)[3]);

# Spline the Fourier components of bpsi 
bpsi_mn_int = Interpolations.interpolate(bpsi_mn,(BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Periodic(OnGrid()))), BSpline(Cubic(Periodic(OnGrid())))));
bpsi_mn_sint = scale(bpsi_mn_int,psin_list,m_index_list,n_index_list);
bpsi_mn_eint = extrapolate(bpsi_mn_sint,(Line(),Periodic(),Periodic()));

# Create the perturbation along uf based on buf_mn = -i/m d(bpsi_mn)/dpsi
buf_mn = []
for p in psin_list
    for m in m_index_list
        for n in n_index_list
            if m !=1 # avoid zero mode
                mode = m_mode_num(m)
                b = (-1im/mode)*Interpolations.gradient(bpsi_mn_eint,p,m,n)[1]
                append!(buf_mn,b)
            else
                append!(buf_mn,0) # zero-out the zero mode
            end
        end
    end
end
buf_mn = permutedims(reshape(buf_mn,length(n_index_list),length(m_index_list),length(psin_list)),(3,2,1));

buf_mn = Complex{Float64}.(buf_mn);

# Inverse Fourier transform to get perturbation along uf
buf = ifft(buf_mn,[2,3]);
buf = real.(buf);

# Spline the calculated field so it can be used by integrator
buf_int = Interpolations.interpolate(buf,(BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Periodic(OnGrid()))), BSpline(Cubic(Periodic(OnGrid())))));
buf_sint = scale(buf_int,psin,ufn,phin);
buf_eint = extrapolate(buf_sint, (Line(),Periodic(),Periodic()));

# Create Poincare maps
P,U = Poincare_Map(Npsi,Nuf,timestep)

# Plotting
gr(size=(1000,1000));
Plots.plot();
for i in 1:(Npsi*Nuf)
    Plots.scatter!(P[i,:],[mod2pi(x) for x in U[i][:]], markerstrokewidth=0,markersize=1.0,legend=false,xlims=(0.0,1.0),ylims=(0.0,2*pi));
end
xlabel!(L"\Psi_N")
ylabel!(L"\theta_f")
poinc_plot = current();

savefig(poinc_plot,savestring);




