#Cell2Node.jl
module C2N

"""Functions that take magnetic field component arrays defined on cells and interpolates them on grid nodes"""

using Interpolations
using NPZ
using Plots
using PyPlot 
import PyPlot

function cell_node_constants(B::Array{Float64,4})
    # definitions of cell grid
    num_r_cells = size(B)[1];
    num_u_cells = size(B)[2];
    num_phi_cells = size(B)[3];
    dn_r = (1.0/num_r_cells);
    dn_u = ((2.0*pi)/num_u_cells);

    # Cell-based grid
    rc = LinRange(0.0+(dn_r/2.0),1.0-(dn_r/2.0),num_r_cells);
    uc = LinRange(0.0+(dn_u/2.0),2.0*pi-(dn_u/2.0),num_u_cells);
    phic = LinRange(0.0+(dn_u/2.0),2.0*pi-(dn_u/2.0),num_phi_cells);
    tn = LinRange(0, size(B)[4],size(B)[4]);
    
    # Node-based grid
    rn = LinRange(0.0,1.0,(num_r_cells+1));
    un = LinRange(0.0,2.0*pi,(num_u_cells+1));
    phin = LinRange(0.0,2.0*pi,(num_phi_cells+1));
    
    return rc,uc,phic,rn,un,phin,tn
end

function cell_node_constants(B::Array{Float64,3})
    # definitions of cell grid
    num_r_cells = size(B)[1];
    num_u_cells = size(B)[2];
    dn_r = (1.0/num_r_cells);
    dn_u = ((2.0*pi)/num_u_cells);

    # Cell-based grid
    rc = LinRange(0.0+(dn_r/2.0),1.0-(dn_r/2.0),num_r_cells);
    uc = LinRange(0.0+(dn_u/2.0),2.0*pi-(dn_u/2.0),num_u_cells);
    tn = LinRange(0, size(B)[3],size(B)[3]);
    
    # Node-based grid
    rn = LinRange(0.0,1.0,(num_r_cells+1));
    un = LinRange(0.0,2.0*pi,(num_u_cells+1));
    
    return rc,uc,rn,un,tn
end
    

# Interpolations with boundary conditions
function cell_interp_and_scale(B::Array{Float64,4},rc,uc,phic,tn)
    B_int = interpolate(B,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnGrid())))));
    B_sint = scale(B_int,rc,uc,phic,tn)
    return B_sint
end
        
function cell_interp_and_scale(B::Array{Float64,3},rc,uc,tn)
    B_int = interpolate(B,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnGrid())))));
    B_sint = scale(B_int,rc,uc,tn)
    return B_sint
end
        
function node_interp_and_scale(B::Array{Float64,4},rn,un,phin,tn)
    B_int = interpolate(B,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnGrid())))));
    B_sint = scale(B_int,rn,un,phin,tn)
    return B_sint
end
        
function node_interp_and_scale(B::Array{Float64,3},rn,un,tn)
    B_int = interpolate(B,(BSpline(Cubic(Line(OnCell()))),BSpline(Cubic(Periodic(OnCell()))),BSpline(Cubic(Periodic(OnGrid())))));
    B_sint = scale(B_int,rn,un,tn)
    return B_sint
end
        
function extrap4d(Bsint)
    return extrapolate(Bsint, (Line(),Periodic(),Periodic(),Line()))
end
        
function extrap3d(Bsint)
    return extrapolate(Bsint, (Line(),Periodic(),Line()))
end 
    
function cell2node(B::Array{Float64,4})
    rc,uc,phic,rn,un,phin,tn = cell_node_constants(B)
    Bc_sint = cell_interp_and_scale(B,rc,uc,phic,tn)
    Bc_eint = extrap4d(Bc_sint)
    Bn = Bc_eint(rn,un,phin,tn)
    Bn_sint = node_interp_and_scale(Bn,rn,un,phin,tn)
    Bn_eint = extrap4d(Bn_sint)
    return Bn, Bn_eint
end
    
function cell2node(B::Array{Float64,3})
    rc,uc,rn,un,tn = cell_node_constants(B)
    Bc_sint = cell_interp_and_scale(B,rc,uc,tn)
    Bc_eint = extrap3d(Bc_sint)
    Bn = Bc_eint(rn,un,tn)
    Bn_sint = node_interp_and_scale(Bn,rn,un,tn)
    Bn_eint = extrap3d(Bn_sint)
    return Bn, Bn_eint
end
    
    
end