using Interpolations
using Dierckx
using DifferentialEquations
using PyCall
using NPZ
using Statistics
using QuadGK
using Plots
using Printf
using LaTeXStrings

filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_new/dt2.scratch./pixie3d.h5"
#filepath = "/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/sawtooth2.scratch/pixie3d.h5"

timestep = 143;
Nr = 100; # Number of initial conditions
Np = 3;
Lmax = 15000;
savestring = @sprintf("PP_db_new@%s.png",(timestep))

function quad_r_int(B_eint,xmin,xmax,yo)
    res, err = quadgk(x -> B_eint(x,yo),xmin,xmax)
    return res
end

function quad_r_int(B_eint,xmin,xmax,yo,zo)
    res, err = quadgk(x -> B_eint(x,yo,zo),xmin,xmax)
    return res
end

function quad_u_int(B_eint,ymin,ymax,xo,zo)
    res, err = quadgk(y -> B_eint(xo,y,zo),ymin,ymax)
    return res
end

function quad_u_int(B_eint,ymin,ymax,xo)
    res, err = quadgk(y -> B_eint(xo,y),ymin,ymax)
    return res
end

function Au(B3_eint,ro,uo,rMA,uMA)
    A = quad_r_int(B3_eint,rMA,ro,uo)
    return A
end

function Au(B3_eint,ro,uo,rMA,uMA,zo)
    A = quad_r_int(B3_eint,rMA,ro,uo,zo)
    return A
end

function Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA)
    A = quad_u_int(B1_eint,uMA,uo,ro) - quad_r_int(B2_eint,rMA,ro,uMA)
    return A
end

function Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA,zo)
    A = quad_u_int(B1_eint,uMA,uo,ro,zo) - quad_r_int(B2_eint,rMA,ro,uMA,zo)
    return A
end

function grid_Au(B3_eint,Ndims::Int)
    if Ndims == 3
        Au_arr = []
        for ro in rn
            for uo in un
                for fo in phin
                    rMA = rmaxis[1]
                    uMA = umaxis[1]
                    append!(Au_arr,Au(B3_eint,ro,uo,rMA,uMA,fo))
                end
            end
        end
        Au_arr = permutedims(reshape(Au_arr, size(phin)[1], size(un)[1],size(rn)[1]),(3,2,1))
    elseif Ndims == 2
        Au_arr = []
        for ro in rn
            for uo in un
                rMA = rmaxis[1]
                uMA = umaxis[1]
                append!(Au_arr,Au(B3_eint,ro,uo,rMA,uMA))
            end
        end
        Au_arr = permutedims(reshape(Au_arr, size(un)[1],size(rn)[1]),(2,1))
    end
    return Au_arr
end

function grid_Aphi(B1_eint,B2_eint,Ndims::Int)
    if Ndims == 3
        Aphi_arr = []
        for ro in rn
            for uo in un
                for fo in phin
                    rMA = rmaxis[1]
                    uMA = umaxis[1]
                    append!(Aphi_arr,Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA,fo))
                end
            end
        end
        Aphi_arr = permutedims(reshape(Aphi_arr, size(phin)[1],size(un)[1],size(rn)[1]),(3,2,1))
    elseif Ndims == 2
        Aphi_arr = []
        for ro in rn
            for uo in un
                rMA = rmaxis[1]
                uMA = umaxis[1]
                append!(Aphi_arr,Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA))
            end
        end
        Aphi_arr = permutedims(reshape(Aphi_arr, size(un)[1],size(rn)[1]),(2,1))
    end
    return Aphi_arr
end

function b1_an(GAp_eint, Ndims::Int)
    if Ndims == 2
        b1_an = []
        for r in rn
            for u in un
                append!(b1_an, Interpolations.gradient(GAp_eint,r,u)[2])
            end
        end
        b1_an = permutedims(reshape(b1_an, size(un)[1],size(rn)[1]),(2,1))
    elseif Ndims == 3
        b1_an = []
        for r in rn
            for u in un
                for f in phin
                    append!(b1_an, Interpolations.gradient(GAp_eint,r,u,f)[2])
                end
            end
        end
        b1_an = permutedims(reshape(b1_an, size(phin)[1],size(un)[1],size(rn)[1]),(3,2,1))
    end
    return b1_an
end
    
function b2_an(GAp_eint,Ndims::Int)
    if Ndims == 2
        b2_an = []
        for r in rn
            for u in un
                append!(b2_an, -Interpolations.gradient(GAp_eint,r,u)[1])
            end
        end
        b2_an = permutedims(reshape(b2_an, size(un)[1],size(rn)[1]),(2,1))
    elseif Ndims == 3 
        b2_an = []
        for r in rn
            for u in un
                for f in phin
                    append!(b2_an, -Interpolations.gradient(GAp_eint,r,u,f)[1])
                end
            end
        end
        b2_an = permutedims(reshape(b2_an, size(phin)[1],size(un)[1],size(rn)[1]),(3,2,1))
    end
    return b2_an
end 

function b3_an(GAu_eint,Ndims::Int)
    if Ndims == 2
        b3_an = []
        for r in rn
            for u in un
                append!(b3_an, Interpolations.gradient(GAu_eint,r,u)[1])
            end
        end
        b3_an = permutedims(reshape(b3_an, size(un)[1],size(rn)[1]),(2,1))
    elseif Ndims == 3
        b3_an = []
        for r in rn
            for u in un
                for f in phin
                    append!(b3_an, Interpolations.gradient(GAu_eint,r,u,f)[1])
                end
            end
        end
        b3_an = permutedims(reshape(b3_an, size(phin)[1],size(un)[1],size(rn)[1]),(3,2,1))
    end
    return b3_an
end

pxr = pyimport("pixie_read_st")
pxr.pixieload(filepath)
eqdsk = pyimport("poinc") # Importing the eqdsk information
eqdsk.eqdsk_info()

# Load data
psi = pxr.load_array(3,4,timestep,timestep+1);
B1 = pxr.load_array(1,0,timestep,timestep+1); # Contravariant components
B2 = pxr.load_array(1,1,timestep,timestep+1);
B3 = pxr.load_array(1,2,timestep,timestep+1);

#psi = -psi; # inverting psi

# Drop singleton dimensions
psi = dropdims(psi,dims=4);
B1 = dropdims(B1,dims=4);
B2 = dropdims(B2,dims=4);
B3 = dropdims(B3,dims=4);

# Toroidally averaged fields
psit = dropdims(mean(psi,dims=3),dims=3);
B1t = dropdims(mean(B1,dims=3),dims=3);
B2t = dropdims(mean(B2,dims=3),dims=3);
B3t = dropdims(mean(B3,dims=3),dims=3);

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
#tn = LinRange(0, size(B3)[4]-1,size(B3)[4]);

# Node-based grid
rn = LinRange(0.0,1.0,(num_r_cells+1));
un = LinRange(0.0,2.0*pi,(num_u_cells+1));
phin = LinRange(0.0,2.0*pi,(num_phi_cells+1));

# Node-based grid dimensions
rdim = size(psi)[1];
udim = size(psi)[2];
fidim = size(psi)[3];
#tdim = size(psi)[4];

# Interpolate on cell-based grid
B1_int_cell = Interpolations.interpolate(B1,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell())))));
B2_int_cell = Interpolations.interpolate(B2,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell())))));
B3_int_cell = Interpolations.interpolate(B3,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell())))));
B1t_int_cell = Interpolations.interpolate(B1t,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell())))));
B2t_int_cell = Interpolations.interpolate(B2t,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell())))));
B3t_int_cell = Interpolations.interpolate(B3t,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell())))));

B1_sint_cell = scale(B1_int_cell,rc,uc,phic);
B2_sint_cell = scale(B2_int_cell,rc,uc,phic);
B3_sint_cell = scale(B3_int_cell,rc,uc,phic);
B1t_sint_cell = scale(B1t_int_cell,rc,uc);
B2t_sint_cell = scale(B2t_int_cell,rc,uc);
B3t_sint_cell = scale(B3t_int_cell,rc,uc);

B1_eint_cell = extrapolate(B1_sint_cell, (Line(),Periodic(),Periodic()));
B2_eint_cell = extrapolate(B2_sint_cell, (Line(),Periodic(),Periodic()));
B3_eint_cell = extrapolate(B3_sint_cell, (Line(),Periodic(),Periodic()));
B1t_eint_cell = extrapolate(B1t_sint_cell, (Line(),Periodic()));
B2t_eint_cell = extrapolate(B2t_sint_cell, (Line(),Periodic()));
B3t_eint_cell = extrapolate(B3t_sint_cell, (Line(),Periodic()));

# Evaluate B on node grid
B1 = B1_eint_cell(rn,un,phin);
B2 = B2_eint_cell(rn,un,phin);
B3 = B3_eint_cell(rn,un,phin);
B1t = B1t_eint_cell(rn,un);
B2t = B2t_eint_cell(rn,un);
B3t = B3t_eint_cell(rn,un);

B3_int = Interpolations.interpolate(B3,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
B1_int = Interpolations.interpolate(B1,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
B2_int = Interpolations.interpolate(B2,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
B3t_int = Interpolations.interpolate(B3t,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
B1t_int = Interpolations.interpolate(B1t,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
B2t_int = Interpolations.interpolate(B2t,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));

B3_sint = scale(B3_int, rn,un,phin);
B1_sint = scale(B1_int, rn,un,phin);
B2_sint = scale(B2_int, rn,un,phin);
B3t_sint = scale(B3t_int, rn,un);
B1t_sint = scale(B1t_int, rn,un);
B2t_sint = scale(B2t_int, rn,un);

B3_eint = extrapolate(B3_sint, (Line(),Periodic(),Periodic()));
B1_eint = extrapolate(B1_sint, (Line(),Periodic(),Periodic()));
B2_eint = extrapolate(B2_sint, (Line(),Periodic(),Periodic()));
B3t_eint = extrapolate(B3t_sint, (Line(),Periodic()));
B1t_eint = extrapolate(B1t_sint, (Line(),Periodic()));
B2t_eint = extrapolate(B2t_sint, (Line(),Periodic()));

# Cartesian maps
X_int = Interpolations.interpolate(pxr.X,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
Y_int = Interpolations.interpolate(pxr.Y,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
Z_int = Interpolations.interpolate(pxr.Z,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));

X_sint = scale(X_int,rn,un,phin);
Y_sint = scale(Y_int,rn,un,phin);
Z_sint = scale(Z_int,rn,un,phin);

X_eint = extrapolate(X_sint,(Line(),Periodic(),Periodic()));
Y_eint = extrapolate(Y_sint,(Line(),Periodic(),Periodic()));
Z_eint = extrapolate(Z_sint,(Line(),Periodic(),Periodic()));

# Preparations of the Python module
pxr.Axes_of_Interpolation(B3)
pxr.Grid_Interpolations(psit,B1,B2,B3)
pxr.Calculation_of_Units_and_Sizes()

psi_min,norm = pxr.Normalization_numbers(psit,B1t,B2t);
pythonresult = pxr.create_r_psi_list(psit,B1t,B2t)
r_of_psi_array = pythonresult[1]; # Pick python outputs
rmaxis = pythonresult[2]
umaxis = pythonresult[3]

# Divergence Cleanup
println("Entering divergence cleanup")
# Vector potential 
GAu = grid_Au(B3_eint,3);
GAp = grid_Aphi(B1_eint,B2_eint,3);

# Spline Vector Potential
GAu = Float64.(GAu);
GAp = Float64.(GAp);

GAu_int = Interpolations.interpolate(GAu,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));
GAp_int = Interpolations.interpolate(GAp,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Periodic(OnGrid())))));

GAu_sint = scale(GAu_int,rn,un,phin);
GAp_sint = scale(GAp_int,rn,un,phin);

GAu_eint = extrapolate(GAu_sint, (Line(),Periodic(),Periodic()));
GAp_eint = extrapolate(GAp_sint, (Line(),Periodic(),Periodic()));

function B1_df(r,u,phi)
    return Interpolations.gradient(GAp_eint,r,u,phi)[2]
end

function B2_df(r,u,phi)
    return -Interpolations.gradient(GAp_eint,r,u,phi)[1]
end

function B3_df(r,u,phi)
    return Interpolations.gradient(GAu_eint,r,u,phi)[1] 
end

function Field_line_df!(du,u,p,t)
    du[1] = B1_df(u[1],u[2],u[3])
    du[2] = B2_df(u[1],u[2],u[3])
    du[3] = B3_df(u[1],u[2],u[3])
end

function condition_cross(u,t,integrator)
    u[3] - integrator.p[1]*2.0*pi
end
function count_cross!(integrator)
    integrator.p[1] = integrator.p[1] + 1
end

cb1 = ContinuousCallback(condition_cross,count_cross!,rootfind = true,save_positions=(false,true))

function field_line_integration_df(rs::Float64,us::Float64)
    u0 = [rs,us,0.0]
    p = [1]
    tspan = (0.0,Lmax)
    prob = ODEProblem(Field_line_df!,u0,tspan,p)
    sol = solve(prob,BS5(),callback = cb1,reltol=1.e-10,abstol=1.e-12,save_everystep=false,save_start=false,save_end=false)
    return sol
end

function Poincare_Map(Nr::Int,Np::Int)
    icr = LinRange(0.0,0.8,Nr)
    icu = LinRange(0.0,2*pi,Np)
    ic_tot = Nr*Np

    r = [[] for i in 1:ic_tot]
    theta = [[] for i in 1:ic_tot]
    phi = [[] for i in 1:ic_tot]
    for i in 1:Nr
        for j in 1:Np
            sol = field_line_integration_df(icr[i],icu[j])
            r[i*j] = sol[1,:]
            theta[i*j] = sol[2,:]
            phi[i*j] = sol[3,:]
        end
    end
    xx = [[] for i in 1:ic_tot]
    zz = [[] for i in 1:ic_tot]
    for i in 1:ic_tot
        xcoord = [X_eint(r[i][j],theta[i][j],0) for j in 1:length(r[i])]
        zcoord = [Z_eint(r[i][j],theta[i][j],0) for j in 1:length(r[i])]
        xx[i] = xcoord
        zz[i] = zcoord
    end
    return xx,zz
end

xx,zz = Poincare_Map(Nr,Np);

# Creating r(psi) and R(r) interpolants for plotting with psi labels

# Lists to be interpolated:
#psi_list = LinRange(0,1,101);
#R_list = [pxr.X[i,1,1] for i in 1:size(pxr.X)[1]];
#r_list = [i/(size(pxr.X)[1]-1) for i in 0:size(pxr.X)[1]-1];

#p2r = Spline1D(psi_list,r_of_psi_array[1],bc="zero") # r(psi)
#r2p = Spline1D(r_of_psi_array[1],psi_list,bc="zero") # psi(r)
#r_points = LinRange(0.2,0.9,Nr);
#p_points = r2p(r_points)
#npzwrite("/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_new/dt.scratch/p_points.npy",p_points)
#r2R = Spline1D(r_list,R_list,bc="zero") # R(r)

# Preparing labels and ticks
#psi_labels=[0.0,0.025,0.1,0.3,0.6,0.8,1.0];
#x_tick_pos = r2R(p2r(psi_labels));
#x_tick_lab = [string(x) for x in psi_labels]

# Plot sizes
x_r = maximum(eqdsk.m.DS.rlim) + 0.1;
x_l = minimum(eqdsk.m.DS.rlim) - 0.1;
z_u = maximum(eqdsk.m.DS.zlim) + 0.1;
z_d = minimum(eqdsk.m.DS.zlim) - 0.1;

gr(size=(2000,2000));
Plots.plot();
Plots.plot!(eqdsk.m.DS.rlim,eqdsk.m.DS.zlim,linecolor=:black);
Plots.plot!(eqdsk.m.DS.rbbbs,eqdsk.m.DS.zbbbs,linecolor = :red);
for i in 1:(Nr*Np)
    Plots.scatter!(xx[i]*eqdsk.e.a,zz[i]*eqdsk.e.a,markerstrokewidth=0,markersize=1.0,legend=false,aspect_ratio=:equal,xlims=(x_l,x_r),ylims=(z_d,z_u),dpi=500);
    #Plots.scatter!(xx[i],zz[i],markerstrokewidth=0,markersize=1.0,legend=false,aspect_ratio=:equal,xlims=(1.0,2.25),ylims=(-1.0,1.0),dpi=500);
end
#xticks!(x_tick_pos*eqdsk.e.a,x_tick_lab)
#xlabel!(L"R (\Psi_N)")
xlabel!("R (m)")
ylabel!("Z (m)")


poinc_plot = current();

savefig(poinc_plot,savestring);







