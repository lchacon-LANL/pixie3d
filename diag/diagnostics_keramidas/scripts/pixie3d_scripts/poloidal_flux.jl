#poloidal_flux.jl

module poloidal_flux

"""Module that reads magnetic field arrays and calculates the poloidal magnetic flux, in the same coordinate system that the field components are given."""

using Interpolations
using Dierckx
using NPZ
using DifferentialEquations
using Plots
using PyPlot 
import PyPlot
using LaTeXStrings
using QuadGK
using Statistics

filepath = "/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/python_arrays/";

"""Functions that convert r and u indexes to (logical) r and u values."""
function CnvNumber2LogicalR(r_num)
    r_dim = size(X,1)
    r_log = r_num/r_dim
    return r_log
end

function CnvNumber2LogicalU(u_num)
    u_dim = size(X,2)
    u_log = (2*pi)*u_num/u_dim
    return u_log
end

"""Locates indices of magnetic axis."""
function magnetic_axis(t)
    Bp = sqrt.(B1t.^2 .+ B2t.^2)
    index = argmin(Bp[:,:,t])
    return index[1], index[2]
end

"""Function that does line integrals of field components to calculate poloidal flux at point r,u."""
function pol_flux(time,r_ind,u_ind)
    r_ma_index, u_ma_index = magnetic_axis(time) # find magnetic axis location
    
    r_ma = CnvNumber2LogicalR(r_ma_index) # Convert magnetic axis index in logical r,theta 
    u_ma = CnvNumber2LogicalU(u_ma_index)
    
    ro = CnvNumber2LogicalR(r_ind) # Convert integration indices to logical r,theta
    uo = CnvNumber2LogicalU(u_ind)
    
    # Integration path is (r_ma,theta_ma)-> (ro,theta_ma)-> (ro,thetao)
    Int1,_ = quadgk(r->B2_eint(r,u_ma,time),r_ma,ro)
    Int2,_ = quadgk(u->B1_eint(ro,u,time),u_ma,uo)
    
    pol_flux = Int1-Int2
    return pol_flux
end

function pol_flux2(time,r_ind,u_ind)
    r_ma_index, u_ma_index = magnetic_axis(time) # find magnetic axis location
    
    r_ma = CnvNumber2LogicalR(r_ma_index) # Convert magnetic axis index in logical r,theta 
    u_ma = CnvNumber2LogicalU(u_ma_index)
    
    ro = CnvNumber2LogicalR(r_ind) # Convert integration indices to logical r,theta
    uo = CnvNumber2LogicalU(u_ind)
    
    # Integration path is (r_ma,theta_ma)-> (ro,theta_ma)-> (ro,thetao)
    Int1,_ = quadgk(r->B2_eint(r,u_ma,time)/X_eint(r,u_ma,0),r_ma,ro)
    Int2,_ = quadgk(u->B1_eint(ro,u,time)/X_eint(ro,u,0),u_ma,uo)
    
    pol_flux = X_eint(ro,uo,0)*(Int1-Int2)
    return pol_flux
end

"""Calculates the psi array for a specific time."""
function psi(t)
    psi = Float64[]
    r_dim = size(X,1)
    u_dim = size(X,2)
    for r in 1:r_dim
        for u in 1:u_dim
            append!(psi,pol_flux2(t,r,u))
        end
    end
    return permutedims(reshape(psi,(u_dim,r_dim)),(2,1))
end

"""Calculates the psi array for all times."""
function psi_array()
    t_dim = size(B1,4)
    r_dim = size(B1,1)
    u_dim = size(B1,2)
    
    P = Array{Float64}(undef,(r_dim,u_dim,0))
    
    for t in 1:t_dim
        P_temp = psi(t)
        P = cat(dims=3,P,P_temp)
        if mod(t,5)==0
            println("t:",t)
        end
    end
    return P
end

################### RUNNING THE PROGRAM #####################
println("Starting...")
# Reading of arrays
X = npzread(filepath * "X.npy")
Z = npzread(filepath * "Z.npy")
B1 = npzread(filepath * "B1_first_line.npy");
B2 = npzread(filepath * "B2_first_line.npy");
B3 = npzread(filepath * "B3_first_line.npy");

println("Arrays read...")

# Toroidal average of fields
B1t = dropdims(mean(B1,dims=3),dims=3);
B2t = dropdims(mean(B2,dims=3),dims=3);
B3t = dropdims(mean(B3,dims=3),dims=3);

# Performing interpolations
X_int = interpolate(X,(BSpline(Quadratic(Periodic(OnGrid()))),BSpline(Quadratic(Periodic(OnGrid()))),BSpline(Quadratic(Periodic(OnGrid())))));
Z_int = interpolate(Z,(BSpline(Quadratic(Line(OnGrid()))),BSpline(Quadratic(Periodic(OnGrid()))),BSpline(Quadratic(Periodic(OnGrid())))));
B1_int = interpolate(B1t,(BSpline(Quadratic(Line(OnGrid()))),BSpline(Quadratic(Periodic(OnGrid()))),BSpline(Quadratic(Line(OnGrid())))));
B2_int = interpolate(B2t,(BSpline(Quadratic(Line(OnGrid()))),BSpline(Quadratic(Periodic(OnGrid()))),BSpline(Quadratic(Line(OnGrid())))));
B3_int = interpolate(B3t,(BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Periodic(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));

# Interpolation scale
r = range(0.0, length=size(B3)[1], stop=1.0);
u = range(0.0, length=size(B3)[2], stop=2.0*pi);
fi = range(0.0, length=size(B3)[3], stop=2.0*pi);
t = range(0, stop=size(B3)[4]-1);

# Rescaling interpolation intervals
X_sint = scale(X_int, r,u,fi);
Z_sint = scale(Z_int, r,u,fi);
B1_sint = scale(B1_int, r,u,t);
B2_sint = scale(B2_int, r,u,t);
B3_sint = scale(B3_int, r,u,t);

# Extrapolations on rescaled
X_eint = extrapolate(X_sint, Line());
Z_eint = extrapolate(Z_sint, Line());
B1_eint = extrapolate(B1_sint, Line());
B2_eint = extrapolate(B2_sint, Line());
B3_eint = extrapolate(B3_sint, Line());

println("Interpolations done...")
# Calculating the poloidal flux array
P = psi_array()

# Saving the result
npzwrite(filepath * "psi_first_line2.npy", P)
println("Done!")
end