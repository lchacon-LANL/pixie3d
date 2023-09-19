using Interpolations
using DifferentialEquations
using PyCall
using NPZ
using Statistics
using QuadGK
using Plots
using FFTW
using Printf
using LaTeXStrings

#filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_tear/dt_sh_m3_n2.scratch/"
filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_new/dt.scratch/"

# Loading fields
bpert = npzread(filepath *  "b_hat_rho_mpi_nm.npy")
q = npzread(filepath * "q_mpi_nm.npy")

# Initialization
#NmodesInPol = 10; # number of modes to discard in poloidal direction
#NmodesInTor = 5; # number of modes to discard in toroidal direction
NmodesToThrow = 1000; # number of modes to discard
timestep = 110;
Npsi = 120; # Number of initial conditions
Nuf = 3;
Lmax = 20000;

function Field_line!(du,u,p,t)
    du[1] = bpsi_eint(u[1],u[2],u[3])
    du[2] = 1
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
    icpsi = LinRange(0.01,0.95,Npsi)
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


############################################################################################################################################################################################################################################################################################################# EXECUTION ########################################################################
################################################################################################################################################################################################################################################################

# Import the eqdsk info
eqdsk = pyimport("poinc") 
eqdsk.eqdsk_info()


# Node-based grid
psin = LinRange(0.0,1.0,size(bpert)[1]);
ufn = LinRange(0.0,2.0*pi,size(bpert)[2]);
phin = LinRange(0.0,2.0*pi,size(bpert)[3]);
tn = LinRange(0, size(bpert)[4]-1,size(bpert)[4]);

# Truncation of original field using Fourier representation
# Passing in Fourier space
bp = bpert[:,:,:,timestep];
bhat = fft(bp,[2,3]);

Nm = size(bhat)[2]; # Poloidal mode range
Nn = size(bhat)[3]; # Toroidal mode range

# Storing of coefficients by maximum value along psi dimension
coeffs = []
for m in 1:Nm
    for n in 1:Nn
        max = maximum(abs.(bhat[:,m,n]))
        append!(coeffs,max)
    end
end

# Ordering of coefficients by size
scoeffs = sortperm(coeffs,rev=true);

# Extracting m and n mode numbers corresponding to indices of the coeeficient vector
# Very important not to mix up the order of the previous nested for-loop!
nvec = []
mvec = []
for i in 1:length(scoeffs)
    m = div(scoeffs[i],Nn)
    if m == Nm
        append!(mvec,m)
    else
        append!(mvec,m+1)
    end
    n = mod(scoeffs[i],Nn)
    if n == 0
        append!(nvec,Nn)
    else
        append!(nvec,n)
    end
end

# Zeroing of the smallest modes
for i in 0:(NmodesToThrow-1)
    bhat[:,mvec[end-i],nvec[end-i]] .= 0;
end

##################################################################################################################
# Vector of averaged poloidal modes
#mvec = [];
#for m in 1:size(bhat)[2]
#    m_mode = mean(bhat[:,m,:])
#    append!(mvec,abs(imag(m_mode)))
#end

# Vector of averaged toroidal modes
#nvec = [];
#for n in 1:size(bhat)[3]
#    n_mode = mean(bhat[:,:,n])
#    append!(nvec,abs(imag(n_mode)))
#end

# Discarding poloidal modes
#smvec = sortperm(mvec,rev=true);
#keep_m_ind = smvec[1:end-NmodesInPol];
#disc_m_ind = smvec[end-NmodesInPol+1:end];

#bhat[:,disc_m_ind,:] .= 0;

# Discarding toroidal modes
#snvec = sortperm(nvec,rev=true);
#keep_n_ind = snvec[1:end-NmodesInTor];
#disc_n_ind = snvec[end-NmodesInTor+1:end];

#bhat[:,:,disc_n_ind] .= 0;

####################################################################################################################
# Shift the frequencies so you can throw away the edges
#bhat_shift = fftshift(bhat,[2,3]);

# Truncate the fft by throwing out high frequencies
#Pol_cutoff = convert(Int64,ceil(NmodesInPol/2));
#Tor_cutoff = convert(Int64,ceil(NmodesInTor/2));
#bhat_trunc = bhat_shift[:,Pol_cutoff:end-Pol_cutoff+1,Tor_cutoff:end-Tor_cutoff+1];

# Shift the frequencies back so inverse works properly
#bhat_trunc_ishift = ifftshift(bhat_trunc,[2,3]);

# Invert truncated fft
#bpsi = ifft(bhat_trunc_ishift,[2,3]);

# Reduced dimensionality grids
#ufn = LinRange(0.0,2.0*pi,size(bpsi)[2]);
#phin = LinRange(0.0,2.0*pi,size(bpsi)[3]);
#######################################################################################################################

# Invert truncated fft
bpsi = ifft(bhat,[2,3]);
bpsi = real.(bpsi);

# Spline the field
bpsi_int = Interpolations.interpolate(bpsi,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
q_int = Interpolations.interpolate(q,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));

bpsi_sint = scale(bpsi_int,psin,ufn,phin);
q_sint = scale(q_int,psin,tn);

bpsi_eint = extrapolate(bpsi_sint, (Line(),Periodic(),Periodic()));
q_eint = extrapolate(q_sint, (Line(),Line()));

P,U = Poincare_Map(Npsi,Nuf,timestep)

gr(size=(1000,1000));
Plots.plot();
for i in 1:(Npsi*Nuf)
    Plots.scatter!(P[i,:],[mod2pi(x) for x in U[i][:]], markerstrokewidth=0,markersize=1.0,legend=false,xlims=(0.0,1.0),ylims=(0.0,2*pi));
end
xlabel!(L"\Psi_N")
ylabel!(L"\theta_f")
poinc_plot = current();

savestring = @sprintf("PP_MC_red_50pc@%s.png",(timestep))
savefig(poinc_plot,savestring);





