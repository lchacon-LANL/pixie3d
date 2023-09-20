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
include("fftUtils.jl")
using .fftUtils

# Paths
filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_new/dt2.scratch./"
#filepath = "/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/sawtooth2.scratch/";

# Inputs
timestep = 143;
Npsi = 120; # Number of initial conditions
Nuf = 3;
Lmax = 10000;
savestring = @sprintf("PP_MC2CAR_db_new@%s.png",(timestep)) # saved plot name

# Loading fields
bpsi_load = npzread(filepath *  "b_hat_rho_mpi.npy");
q = npzread(filepath * "q_mpi.npy");
db3_load = npzread(filepath * "db3.npy");
Nm = size(bpsi_load)[2]; # number of m-modes
Nn = size(bpsi_load)[3]; # number of n-modes

# Loading Cartesian maps
X_load = npzread(filepath * "X.npy");
Z_load = npzread(filepath * "Z.npy");


# Integrator functions for the problem
function Field_line!(du,u,p,t)
    du[1] = bpsi_eint(u[1],u[2],u[3])
    du[2] = 1 + buf_eint(u[1],u[2],u[3])
    du[3] = q_eint(u[1],p[1]) + db3_eint(u[1],u[2],u[3])
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
    #icpsi = sqrt.(LinRange(0.01,0.95,Npsi))
    #icpsi = p_points
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
    xx = [[] for i in 1:ic_tot]
    zz = [[] for i in 1:ic_tot]
    for i in 1:ic_tot
        xcoord = [X_eint(psi[i][j],mod2pi(theta_f[i][j]),0) for j in 1:length(psi[i])]
        zcoord = [Z_eint(psi[i][j],mod2pi(theta_f[i][j]),0) for j in 1:length(psi[i])]
        xx[i] = xcoord
        zz[i] = zcoord
    end
    return xx,zz
end

function div_free_buf(psin_list,m_index_list,n_index_list,bpsi_mn_eint)
    """Create the perturbation along uf based on buf_mn = -i/m d(bpsi_mn)/dpsi"""
    buf_mn = []
    Nm = length(m_index_list)
    for p in psin_list
        for m in m_index_list
            for n in n_index_list
                if m !=1 # avoid zero mode
                    mode = m_mode_num(m,Nm)
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
    return buf
end
###########################################################################################################################################################################################EXECUTION##################################################################################################################################################################################################

# Import the eqdsk info
eqdsk = pyimport("poinc") 
eqdsk.eqdsk_info()

# Read p-points
#p_points = npzread("/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_new/dt.scratch/p_points.npy")

# Select timestep
bpsi = 2*pi*bpsi_load[:,:,:,timestep]; # Rescale the perturbation by 2pi
X = X_load[:,:,:,timestep];
Z = Z_load[:,:,:,timestep];
db3 = db3_load[:,:,:,timestep];


# Node-based grid
psin = LinRange(0.0,1.0,size(bpsi)[1]);
ufn = LinRange(0.0,2.0*pi,size(bpsi)[2]);
phin = LinRange(0.0,2.0*pi,size(bpsi)[3]);
tn = LinRange(0, size(bpsi_load)[4]-1,size(bpsi_load)[4]);

# Spline the field
bpsi_int = Interpolations.interpolate(bpsi,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
db3_int = Interpolations.interpolate(db3,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
q_int = Interpolations.interpolate(q,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));

bpsi_sint = scale(bpsi_int,psin,ufn,phin);
db3_sint = scale(db3_int,psin,ufn,phin);
q_sint = scale(q_int,psin,tn);

bpsi_eint = extrapolate(bpsi_sint, (Line(),Periodic(),Periodic()));
db3_eint = extrapolate(db3_sint, (Line(),Periodic(),Periodic()));
q_eint = extrapolate(q_sint, (Line(),Line()));

# Spline the maps
X_int = Interpolations.interpolate(X, (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
Z_int = Interpolations.interpolate(Z, (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));

X_sint = scale(X_int, psin,ufn,phin);
Z_sint = scale(Z_int, psin,ufn,phin);

X_eint = extrapolate(X_sint, (Line(),Periodic(),Periodic()));
Z_eint = extrapolate(Z_sint, (Line(),Periodic(),Periodic()));

# Calculate perturbation along magnetic angle (uf) with Fourier components
# Fourier transform to get bpsi_mn
bpsi_mn = fft(bpsi,[2,3]);

# Create lists to scale the interpolations for bpsi_mn
psin_list = LinRange(0,1.0,size(bpsi_mn)[1]);
m_index_list = LinRange(1,size(bpsi_mn)[2],size(bpsi_mn)[2]);
n_index_list = LinRange(1,size(bpsi_mn)[3],size(bpsi_mn)[3]);

# Spline the Fourier components of bpsi 
bpsi_mn_int = Interpolations.interpolate(bpsi_mn,(BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Periodic(OnGrid()))), BSpline(Cubic(Periodic(OnGrid())))));
bpsi_mn_sint = scale(bpsi_mn_int, psin_list,m_index_list,n_index_list);
bpsi_mn_eint = extrapolate(bpsi_mn_sint, (Line(),Periodic(),Periodic()));

# Calculate divergence free buf
buf = div_free_buf(psin_list,m_index_list,n_index_list,bpsi_mn_eint);

# Spline the calculated field so it can be used by integrator
buf_int = Interpolations.interpolate(buf,(BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Periodic(OnGrid()))), BSpline(Cubic(Periodic(OnGrid())))));
buf_sint = scale(buf_int,psin,ufn,phin);
buf_eint = extrapolate(buf_sint, (Line(),Periodic(),Periodic()));

# Create Poincare maps
xx,zz = Poincare_Map(Npsi,Nuf,timestep)

# Plot sizes
x_r = maximum(eqdsk.m.DS.rlim) + 0.1;
x_l = minimum(eqdsk.m.DS.rlim) - 0.1;
z_u = maximum(eqdsk.m.DS.zlim) + 0.1;
z_d = minimum(eqdsk.m.DS.zlim) - 0.1;

# Plotting
gr(size=(2000,2000));
Plots.plot();
Plots.plot!(eqdsk.m.DS.rlim,eqdsk.m.DS.zlim,linecolor=:black);
Plots.plot!(eqdsk.m.DS.rbbbs,eqdsk.m.DS.zbbbs,linecolor = :red);
for i in 1:(Npsi*Nuf)
    Plots.scatter!(xx[i]*eqdsk.e.a,zz[i]*eqdsk.e.a,markerstrokewidth=0,markersize=1.0,legend=false,aspect_ratio=:equal,xlims=(x_l,x_r),ylims=(z_d,z_u), dpi=500);
    #Plots.scatter!(xx[i],zz[i],markerstrokewidth=0,markersize=1.0,legend=false,aspect_ratio=:equal,xlims=(1.0,2.25),ylims=(-1.0,1.0), dpi=500);
end
xlabel!("R (m)")
ylabel!("Z (m)")

poinc_plot = current();

savefig(poinc_plot,savestring);



