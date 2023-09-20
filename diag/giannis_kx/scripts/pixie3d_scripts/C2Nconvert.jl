using Interpolations
using PyCall
using MPI
using NPZ
using Statistics
MPI.Init()

# Set up of parallel communicator
comm = MPI.COMM_WORLD
size_comm = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)

# Inputs

filepath = "/net/scratch3/giannis_kx/pixie3d/tests/sawtooth/sawtooth.scratch/pixie3d.h5"
#filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_tear/dt_sh_m3_n2.scratch/pixie3d.h5"
#filepath = "/net/scratch4/giannis_kx/pixie3d/iter/int_kink/11/11_new_visc.scratch/pixie3d.h5"
#filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_new/dt.scratch/pixie3d.h5"
#metricpath = "/net/scratch3/giannis_kx/FTLE/11/shaped_metric_coeff6.npz"

tstart = 0;
tend = 1000;

# Simulation loading
time_interval = tend-tstart; 
chunk = floor(time_interval/size_comm); # loading chunk
rank_tstart = convert(Int,tstart+rank*chunk); # Starting and ending points of each proccess
rank_tend = convert(Int,tstart+(rank+1)*chunk);
final_chunk_size = convert(Int,tend-(tstart+(size_comm-1)*chunk));


pxr = pyimport("pixie_read_st")
pxr.pixieload(filepath)

# Loading of data
if rank == size_comm-1
    psi = pxr.load_array(3,4,rank_tstart,tend);
    B1 = pxr.load_array(1,0,rank_tstart,tend); # Contravariant components
    B2 = pxr.load_array(1,1,rank_tstart,tend);
    B3 = pxr.load_array(1,2,rank_tstart,tend);
else
    psi = pxr.load_array(3,4,rank_tstart,rank_tend);
    B1 = pxr.load_array(1,0,rank_tstart,rank_tend);
    B2 = pxr.load_array(1,1,rank_tstart,rank_tend);
    B3 = pxr.load_array(1,2,rank_tstart,rank_tend);
end


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


# Interpolate on cell-based grid
B1_int_cell = Interpolations.interpolate(B1,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Line(OnGrid())))));
B2_int_cell = Interpolations.interpolate(B2,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Line(OnGrid())))));
B3_int_cell = Interpolations.interpolate(B3,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Line(OnGrid())))));

B1_sint_cell = scale(B1_int_cell,rc,uc,phic,tn);
B2_sint_cell = scale(B2_int_cell,rc,uc,phic,tn);
B3_sint_cell = scale(B3_int_cell,rc,uc,phic,tn);

B1_eint_cell = extrapolate(B1_sint_cell, (Line(),Periodic(),Periodic(),Line()));
B2_eint_cell = extrapolate(B2_sint_cell, (Line(),Periodic(),Periodic(),Line()));
B3_eint_cell = extrapolate(B3_sint_cell, (Line(),Periodic(),Periodic(),Line()));

# Evaluate B on node grid
B1 = B1_eint_cell(rn,un,phin,tn);
B2 = B2_eint_cell(rn,un,phin,tn);
B3 = B3_eint_cell(rn,un,phin,tn);

MPI.Barrier(comm)
if rank != 0
    MPI.Send(B1,0,rank+1001,comm)
    MPI.Send(B2,0,rank+1002,comm)
    MPI.Send(B3,0,rank+1003,comm)
elseif rank == 0
    b_dim1 = size(B1)[1]
    b_dim2 = size(B1)[2]
    b_dim3 = size(B1)[3]
    b_dim4 = size(B1)[4]    
    B1_composite = Array{Float64,4}(undef,b_dim1,b_dim2,b_dim3,0)
    B1_composite = cat(dims=4,B1_composite,B1)
    B2_composite = Array{Float64,4}(undef,b_dim1,b_dim2,b_dim3,0)
    B2_composite = cat(dims=4,B2_composite,B2)
    B3_composite = Array{Float64,4}(undef,b_dim1,b_dim2,b_dim3,0)
    B3_composite = cat(dims=4,B3_composite,B3)
    for r in range(1,stop=size_comm-2)
        b1_recv = Array{Float64,4}(undef,b_dim1,b_dim2,b_dim3,b_dim4)
        b2_recv = Array{Float64,4}(undef,b_dim1,b_dim2,b_dim3,b_dim4)
        b3_recv = Array{Float64,4}(undef,b_dim1,b_dim2,b_dim3,b_dim4)
        MPI.Recv!(b1_recv,r,r+1001,comm)
        MPI.Recv!(b2_recv,r,r+1002,comm)
        MPI.Recv!(b3_recv,r,r+1003,comm)
        global B1_composite = cat(dims=4,B1_composite,b1_recv)
        global B2_composite = cat(dims=4,B2_composite,b2_recv)
        global B3_composite = cat(dims=4,B3_composite,b3_recv)
    end
    b1_recv_fin = Array{Float64,4}(undef,b_dim1,b_dim2,b_dim3,final_chunk_size)
    b2_recv_fin = Array{Float64,4}(undef,b_dim1,b_dim2,b_dim3,final_chunk_size)
    b3_recv_fin = Array{Float64,4}(undef,b_dim1,b_dim2,b_dim3,final_chunk_size)
    MPI.Recv!(b1_recv_fin,size_comm-1,size_comm-1+1001,comm)
    MPI.Recv!(b2_recv_fin,size_comm-1,size_comm-1+1002,comm)
    MPI.Recv!(b3_recv_fin,size_comm-1,size_comm-1+1003,comm)
    B1_composite = cat(dims=4,B1_composite,b1_recv_fin)
    B2_composite = cat(dims=4,B2_composite,b2_recv_fin)
    B3_composite = cat(dims=4,B3_composite,b3_recv_fin)
    
    B1N = Float64.(B1_composite)
    B2N = Float64.(B2_composite)
    B3N = Float64.(B3_composite)
    
    npzwrite("/net/scratch3/giannis_kx/pixie3d/tests/sawtooth/python_arrays/" * "B1.npy",Float64.(B1N))
    npzwrite("/net/scratch3/giannis_kx/pixie3d/tests/sawtooth/python_arrays/" * "B2.npy",Float64.(B2N))
    npzwrite("/net/scratch3/giannis_kx/pixie3d/tests/sawtooth/python_arrays/" * "B3.npy",Float64.(B3N))
end    
MPI.Barrier(comm)
