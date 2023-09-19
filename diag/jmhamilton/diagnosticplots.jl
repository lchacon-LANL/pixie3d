# diagnosticplots.jl
using HDF5
using Plots
#using EllipsisNotation
using Statistics
using CairoMakie

filepath = pwd()*string(/); # Directory of .h5 file
filename = "pixie3d.h5"; # This is the default filename from pixplot

function getTimesteps(fileid)
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
    return sort(timesteps2), sort(time_stamps)
end 
    
function loadh5(fileid,vartype,myvariable,tfirst,tlast)
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

pixiefile = h5open(filepath * filename,"r"); # Load the .h5 file in Julia (file structure is retained)
timesteps, time_stamps = getTimesteps(pixiefile); # Load the timestep values
tA = 1.8444743E-04 # Alfven time in ms
tA = tA*sqrt(2) # For Deuterium
tA = tA*sqrt(300.) # For the 300x density simulation
timesteps_ms = round.(timesteps*tA, digits=2) # Convert to ms

P = loadh5(pixiefile,"Diagnostics","pi+pe",0,length(timesteps)); # Scalars
Te = loadh5(pixiefile,"Car_variables","Te",0,length(timesteps));
#Tecnv = loadh5(pixiefile,"Cnv_variables","JxTe",0,length(timesteps));
rho = loadh5(pixiefile,"Car_variables","rho",0,length(timesteps));

Bx = loadh5(pixiefile,"Car_variables","Bx",0,length(timesteps)); # Cartesian vector components
By = loadh5(pixiefile,"Car_variables","By",0,length(timesteps));
Bz = loadh5(pixiefile,"Car_variables","Bz",0,length(timesteps));
Jx = loadh5(pixiefile,"Car_variables","Jx",0,length(timesteps));
Jy = loadh5(pixiefile,"Car_variables","Jy",0,length(timesteps));
Jz = loadh5(pixiefile,"Car_variables","Jz",0,length(timesteps));

X = read(pixiefile["Timestep_" * string(time_stamps[1])]["nodes"]["X"]);
Y = read(pixiefile["Timestep_" * string(time_stamps[1])]["nodes"]["Y"]);
Z = read(pixiefile["Timestep_" * string(time_stamps[1])]["nodes"]["Z"]);
R = sqrt.(X.^2 + Y.^2); # Cylindrical radius
A = sqrt.((X[:,:,1] .- X[1,1,1]).^2 + (Z[:,:,1] .- Z[1,1,1]).^2); # Poloidal radius, assumes grid is axisymmetric
theta = mod.(atan.((Z[:,:,1] .- Z[1,1,1]),(R[:,:,1] .- R[1,1,1])),2pi); # Poloidal angle
#phi = mod.(atan.(-Y,X),2pi); # Toroidal angle
X = transpose(X[:,:,1]);
Z = transpose(Z[:,:,1]);

Xcell = read(pixiefile["Timestep_" * string(time_stamps[1])]["cells"]["X"]);
Ycell = read(pixiefile["Timestep_" * string(time_stamps[1])]["cells"]["Y"]);
Rcell = sqrt.(Xcell.^2 + Ycell.^2);
J_phi = loadh5(pixiefile,"Cov_variables","J_3",0,length(timesteps)); # Cell-based Jphi*R

close(pixiefile)

J_phi = J_phi./Rcell; # Remove Jacobian factor
if J_phi[1,1,1,1] < 0. ; J_phi .= - J_phi ; end

Br2 = Bx.^2 + By.^2;
Bmag2 = Bz.^2 + Br2;
Jpar = Jx.*Bx + Jy.*By + Jz.*Bz; # Node-based J*B
JparBi = Jpar./Bmag2;
#Jphi = Jy.*cos.(phi) + Jx.*sin.(phi); # Node-based

rn = LinRange(0.,1.,size(J_phi)[1]);
dr = rn[2]; 
Jphi_total = Array{Float64}(undef, size(J_phi)[4]);
Te_total = Array{Float64}(undef, size(Te)[4]);
for t in range(1,size(J_phi)[4]) # Cell-based iteration
    global Jtot = 0. ;
    for th in range(1,size(J_phi)[2]) # Sum up J_phi*r over all cells
	for r in range(1,size(J_phi)[1])
		dr1 = 0.25*(A[r+1,th] + A[r+1,th+1] - A[r,th] - A[r,th+1])^2; # Calculate area of cell (error ~ dx**2)
		dr2 = 0.25*(A[r+1,th+1] + A[r,th+1] - A[r+1,th] - A[r,th])^2;
		dth1 = 0.0625*(A[r+1,th] + A[r+1,th+1] + A[r,th] + A[r,th+1])^2*(mod(theta[r+1,th] + theta[r+1,th+1] - theta[r,th] - theta[r,th+1],2pi))^2;
		dth2 = 0.0625*(A[r+1,th] + A[r+1,th+1] + A[r,th] + A[r,th+1])^2*(mod(theta[r+1,th+1] + theta[r,th+1] - theta[r+1,th] - theta[r,th],2pi))^2;
		dAp = sqrt((dr1+dth1)*(dr+dth2)); # Area of poloidal cell (axisymmetric)
		global Jtot += J_phi[r,th,1,t]*dAp;
	end
    end
    Jphi_total[t] = Jtot;
end
for t in range(1,size(Te)[4]) # Node-based iteration
    global Ttot = Te[1,1,1,t];
    for th in range(2,size(Te)[2]) # Sum up Te at all other points
        for r in range(2,size(Te)[1]) # Do not count r=0 for all theta values (singular point), or double-count theta=0,2pi
                global Ttot += Te[r,th,1,t];
        end
    end
    Te_total[t] = Ttot;
end

Jphi_total = Jphi_total./Jphi_total[1]; # Normalize to initial value
Te_total = Te_total./Te_total[1];

ENV["GKSwstype"]="nul";

r = range(0,1,length=size(P)[1]);

P0 = P[1,1,1,1]
T0 = Te[1,1,1,1]

Ptor = dropdims(mean(P,dims=3),dims=3)./P0;
Tetor = dropdims(mean(Te,dims=3),dims=3)./T0;
#Tecnvtor = dropdims(mean(Tecnv,dims=3),dims=3)./T0;
rhotor = dropdims(mean(rho,dims=3),dims=3);
JpBtor = dropdims(mean(JparBi,dims=3),dims=3); 

iv = trunc(Int,size(Ptor)[2]*0.6875) # Index location of the line reaching the divertor from the axis
it = trunc(Int,size(Ptor)[3]*0.10) # Index location of the time 10% of the way through the simulation

plot1 = Plots.plot(r, [Ptor[:,1,1], Ptor[:,iv,it*5], Ptor[:,iv,it*9], Ptor[:,iv,end]], title ="Vertical cut n=0 Pressure profiles", labels=["t = 1" "t = $(it*5)" "t = $(it*9)" "t = $(it*10)"], linewidth=3, ylims=(0,maximum(Ptor)*1.05), xlabel="Radial Position at θ = 3π/2", ylabel="Pressure", xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22, legendfontsize=18)

Plots.savefig(plot1,"P_profile.png")

animP = @animate for i in 1:size(Ptor)[3]
    Plots.plot(r,Ptor[:,iv,i], title ="Vertical cut n=0 Pressure profile", label="P(r)", linewidth=3, ylims=(0,maximum(Ptor)*1.05), xlabel="Radial Position at θ = 3π/2", ylabel="Pressure", xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22, legendfontsize=18);
end

Plots.gif(animP,"P_profile.gif",fps=6)

animT = @animate for i in 1:size(Tetor)[3]
    Plots.plot(r,Tetor[:,iv,i], title ="Vertical cut n=0 Temperature profile", label="Te(r)", linewidth=3, ylims=(0,maximum(Tetor)*1.05), xlabel="Radial Position at θ = 3π/2", ylabel="Temperature", xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22, legendfontsize=18);
end

Plots.gif(animT,"T_profile.gif",fps=6)

animD = @animate for i in 1:size(rhotor)[3]
    Plots.plot(r,rhotor[:,iv,i], title ="Vertical cut n=0 Density profile", label="ρ(r)", linewidth=3, ylims=(0,maximum(rhotor)*1.05), xlabel="Radial Position at θ = 3π/2", ylabel="Density", xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22, legendfontsize=18);
end
   
Plots.gif(animD,"D_profile.gif",fps=6)

plotJ = Plots.scatter(timesteps_ms, Jphi_total, title = "Total Toroidal Current Density", linewidth=3, ylims=(0,maximum(Jphi_total)), xlabel="Time (ms)", ylabel="Jϕ" , xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22)

Plots.savefig(plotJ,"J(t).png")

plotE = Plots.scatter(timesteps_ms, Te_total, title = "Total Thermal Energy", linewidth=3, ylims=(0,maximum(Te_total)), xlabel="Time (ms)", ylabel="Thermal Energy" , xtickfontsize=18, ytickfontsize=18, guidefontsize=20, titlefontsize=22)

Plots.savefig(plotE,"E(t).png")

######################################
#anim2d = @animate for i in 1:size(Tetor)[3]
#    Plots.heatmap(range(0,2,size(Tecnvtor)[2]),range(0,1,size(Tecnvtor)[1]), Tecnvtor[:,:,i], c = :thermal, titlefontsize=22, title="Temperature", ylabel="r", xlabel="θ (π radians)", guidefontsize=22, tickfontsize=18, legendfontsize=16, clims=(0,maximum(Tecnvtor)*1.05))
#end

#Plots.gif(anim2d,"Temp_2D_cnv.gif",fps=6)

######################################

plot2d, myaxis, mydata = CairoMakie.surface(X[:,:,1], Z[:,:,1], transpose(Tetor[:,:,1]), figure=(; resolution=(900,1200), fontsize=30), axis=(type=Axis3, aspect= :data, azimuth = 1.5pi, elevation = 0.5pi, zlabelvisible=false, zticksvisible=false, zticklabelsvisible=false, title="Temperature", titlesize=30, ylabel="Z", xlabel="R", ylabelrotation=0.5*pi, xlabelvisible=false, ylabelvisible=false, xticklabelsvisible=false, yticklabelsvisible=false, xticksvisible=false, yticksvisible=false), interpolation=true, colormap= :gist_rainbow, colorrange=(0.,maximum(Tetor[:,:,:])), shading=false)

mytext = CairoMakie.text!(myaxis,0.6,0.9,0.; text="t = " * string(timesteps[1]) * " τA = "*string(timesteps_ms[1])*" ms", fontsize = 30, space = :relative)
mytext.fontsize = 30

CairoMakie.Colorbar(plot2d[1,2], mydata, height=Relative(1.), width = 40, ticklabelsize=24)

CairoMakie.save("Temp_2D_car.png",plot2d, px_per_unit = 1)

myiterator = range(1,size(Tetor)[3],step=1)

CairoMakie.record(plot2d,"Temp_2D_car.mp4", myiterator; framerate = 20) do i
	mydata.z = transpose(Tetor[:,:,i])
	mytext.text = "t = " * string(timesteps[i]) * " τA = "*string(timesteps_ms[i])*" ms"
end

plot2drho, myaxisrho, mydatarho = CairoMakie.surface(X[:,:,1], Z[:,:,1], transpose(rhotor[:,:,1]), figure=(; resolution=(900,1200), fontsize=30), axis=(type=Axis3, aspect= :data, azimuth = 1.5pi, elevation = 0.5pi, zlabelvisible=false, zticksvisible=false, zticklabelsvisible=false, title="Density", titlesize=30, ylabel="Z", xlabel="R", ylabelrotation=0.5*pi, xlabelvisible=false, ylabelvisible=false, xticklabelsvisible=false, yticklabelsvisible=false, xticksvisible=false, yticksvisible=false), interpolation=true, colormap= :dense, colorrange=(minimum(rhotor[:,:,:]),maximum(rhotor[:,:,:])), shading=false)

mytextrho = CairoMakie.text!(myaxisrho,0.6,0.9,0.; text="t = " * string(timesteps[1]) * " τA = "*string(timesteps_ms[1])*" ms", fontsize = 30, space = :relative)

CairoMakie.Colorbar(plot2drho[1,2], mydatarho, height=Relative(1.), width = 40, ticklabelsize=24)

CairoMakie.save("Den_2D_car.png",plot2drho, px_per_unit = 1)

myiterator = range(1,size(rhotor)[3],step=1)

CairoMakie.record(plot2drho,"Den_2D_car.mp4", myiterator; framerate = 20) do i
        mydatarho.z = transpose(rhotor[:,:,i])
	mytextrho.text = "t = " * string(timesteps[i]) * " τA = "*string(timesteps_ms[i])*" ms"
end

plot2dp, myaxisp, mydatap = CairoMakie.surface(X[:,:,1], Z[:,:,1], transpose(Ptor[:,:,1]), figure=(; resolution=(900,1200), fontsize=30), axis=(type=Axis3, aspect= :data, azimuth = 1.5pi, elevation = 0.5pi, zlabelvisible=false, zticksvisible=false, zticklabelsvisible=false, title="Pressure", titlesize=30, ylabel="Z", xlabel="R", ylabelrotation=0.5*pi, xlabelvisible=false, ylabelvisible=false, xticklabelsvisible=false, yticklabelsvisible=false, xticksvisible=false, yticksvisible=false), interpolation=true, colormap= :gist_rainbow, shading=false)

mytextp = CairoMakie.text!(myaxisp,0.6,0.9,0.; text="t = " * string(timesteps[1]) * " τA = "*string(timesteps_ms[1])*" ms", fontsize = 30, space = :relative)

CairoMakie.Colorbar(plot2dp[1,2], mydatap, height=Relative(1.), width = 40, ticklabelsize=24)

CairoMakie.save("P_2D_car.png",plot2dp, px_per_unit = 1)

myiterator = range(1,size(Ptor)[3],step=1)

CairoMakie.record(plot2dp,"P_2D_car.mp4", myiterator; framerate = 20) do i
        mydatap.z = transpose(Ptor[:,:,i])
	mytextp.text = "t = " * string(timesteps[i]) * " τA = "*string(timesteps_ms[i])*" ms"
end

plot2dJ, myaxisJ, mydataJ = CairoMakie.surface(X[:,:,1], Z[:,:,1], transpose(JpBtor[:,:,1]), figure=(; resolution=(900,1200), fontsize=30), axis=(type=Axis3, aspect= :data, azimuth = 1.5pi, elevation = 0.5pi, zlabelvisible=false, zticksvisible=false, zticklabelsvisible=false, title="J*B/|B|^2", titlesize=30, ylabel="Z", xlabel="R", ylabelrotation=0.5*pi, xlabelvisible=false, ylabelvisible=false, xticklabelsvisible=false, yticklabelsvisible=false, xticksvisible=false, yticksvisible=false), interpolation=true, colormap= :gist_rainbow, shading=false)

mytextJ = CairoMakie.text!(myaxisJ,0.6,0.9,0.; text="t = " * string(timesteps[1]) * " τA = "*string(timesteps_ms[1])*" ms", fontsize = 30, space = :relative)

CairoMakie.Colorbar(plot2dJ[1,2], mydataJ, height=Relative(1.), width = 40, ticklabelsize=24)

CairoMakie.save("J_2D_car.png",plot2dJ, px_per_unit = 1)

myiterator = range(1,size(JpBtor)[3],step=1)

CairoMakie.record(plot2dJ,"J_2D_car.mp4", myiterator; framerate = 20) do i
        mydataJ.z = transpose(JpBtor[:,:,i])
	mytextJ.text = "t = " * string(timesteps[i]) * " τA = "*string(timesteps_ms[i])*" ms"
end


