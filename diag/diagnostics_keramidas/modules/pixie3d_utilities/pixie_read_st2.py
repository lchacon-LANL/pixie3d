import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata
from skimage import measure
import scipy.optimize
import scipy.integrate
from scipy.integrate import romb
from scipy.integrate import odeint
import warnings
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit
from matplotlib.tri import Triangulation, LinearTriInterpolator
import math
import pylab as py
import sys
import numpy.ma as ma
from decimal import Decimal
import os
from importlib import reload
import h5py
import pickle as pkl
import mpmath
from functools import wraps
from numba import jit, float64
#import solve_v30 as so # REMOVED BY JMHAMILTON
#import pyximport
#pyximport.install(reload_support=True)
#del sys.modules['solve']
#reload(so)
#import solve as so
mpmath.mp.dps=30 # floating point precision

"""Collection of functions for the processing of Pixie3d output files. These functions are mainly used by the 
   magnetic_coordinates.jl script by they are also useful by themselves. The main uses of this module is to:
   1) load variables in numpy arrays
   2) transform cell-centered data into node centered data
   3) locate position of X-point and use it to normalize poloidal flux
   4) provide a list of point on the theta = 0 line that lie between the magnetic axis and the separatrix and they are psi = 0.1 
   apart. Those points will be used by the magnetic_coordinates script to do flux surface integration for the projection to 
   magnetic coordinates.
   5) Contains a function that can initialize a blob of initial conditions so that it can be run with nemato.
   
   CAUTION: The functions below are sensitive to the following things:
   1) Positioning of different variable types inside the h5 tree. Functions affected: pixieload(), load_array() (e.g., the 
   covariant_variables dictionary should be the 3rd variable type)
   2) Positioning of different quantities within their variable type dictionary.
   3) Some quantities are defined on nodes and some are defined on cell centers. There exist functions that turn cell-centered
   quantities into node centered ones so that they can be plotted without a missing part at the edge.
   4) Certain functions require an initial guess because they contain non-linear solvers that need to be initialized. Those are:
   intersection_of_horizontal_line_with_flux_surface()-> requires guess for x-coordinate of q95 surface at theta=0.
   findXpoint(), magneticAxisAndXpoint(), Normalization_numbers() -> require initial guess for (r,theta) coordinate of X-point 
   with r in [0,1], theta in [0,2pi].
   pntCnvInGrid_shaped(),pntCnvInGrid_simple_toroidal() -> requires initial guess (r,theta) for conversion to logical grid.
   
   IMPORTANT: When one wants to run magnetic_coordinates.jl, one needs to uncomment the line matplotlib.use('Agg'). If one wants 
   to import the module and do some analysis then that line should be commented back. It changes the matplotlib backend required
   for plotting and disables it when the code runs in an allocation.
"""


def ignore_warnings(f):
    """Function wrapper for ignoring warnings."""
    @wraps(f)
    def inner(*args, **kwargs):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("ignore")
            response = f(*args, **kwargs)
        return response
    return inner

# Major radius
R0= 1.65

def timestamps(filepath):
    """Sorting the time steps by time value."""
    f = h5py.File(filepath,'r')
    timesteps = list(f.keys())
    
    tstamp_list = []
    for ts in timesteps[:-1]: # Sorting timestep list / avoiding Visit expressions
        tstamp = ts.split('_')[1]
        tstamp_list.append(int(tstamp))
    timesteps = list(['Timestep_'+str(x) for x in sorted(tstamp_list)])
    for t in timesteps:
        print(t+"\n")


def pixieload(filepath):
    """Function that reads in the pixie h5 output file and returns dictionary of terms.""" 
    global f, timesteps, variables, car_var, cnv_var, cov_var, diag_var, pert_var,cell_var, X, Y, Z, Xc, Yc, Zc
    
    f = h5py.File(filepath,'r')
    timesteps = list(f.keys())
    
    tstamp_list = []
    for ts in timesteps[:-1]: # Sorting timestep list / avoiding Visit expressions
        tstamp = ts.split('_')[1]
        tstamp_list.append(int(tstamp))
    timesteps = list(['Timestep_'+str(x) for x in sorted(tstamp_list)])

    variables = list(f[timesteps[0]].keys())
    car_var = list(f[timesteps[0]][variables[0]].keys())
    cnv_var = list(f[timesteps[0]][variables[1]].keys())
    cov_var = list(f[timesteps[0]][variables[2]].keys())
    diag_var = list(f[timesteps[0]][variables[3]].keys())
    pert_var = list(f[timesteps[0]][variables[4]].keys())
    cell_var = list(f[timesteps[0]][variables[5]].keys())
    node_var = list(f[timesteps[0]][variables[6]].keys())
    X = np.asarray(f[timesteps[0]][variables[6]][node_var[0]]) # Node coordinates
    Y = np.asarray(f[timesteps[0]][variables[6]][node_var[1]])
    Z = np.asarray(f[timesteps[0]][variables[6]][node_var[2]])
    X = np.swapaxes(X,0,2)
    Y = np.swapaxes(Y,0,2)
    Z = np.swapaxes(Z,0,2)
    Xc = np.asarray(f[timesteps[0]][variables[5]][cell_var[0]]) # Cell coordinates
    Yc = np.asarray(f[timesteps[0]][variables[5]][cell_var[1]])
    Zc = np.asarray(f[timesteps[0]][variables[5]][cell_var[2]])
    Xc = np.swapaxes(Xc,0,2)
    Yc = np.swapaxes(Yc,0,2)
    Zc = np.swapaxes(Zc,0,2)
    print("timesteps=",len(timesteps))
    print("Dictionary of terms:")
    print("Variables:", variables)
    print("Cartesian:", car_var)
    print("Contravariant:", cnv_var)
    print("Covariant:", cov_var)
    print("Diagnostic:", diag_var)
    print("Perturbations:", pert_var)
    print("Cell:",cell_var)
    print("Node:",node_var)

def load_array(var_type,var_num,t_start,t_stop):
    """"Loads individual arrays based in their position in the dictionary of terms of pixieload."""
    if t_start == None:
        t_start = 0
    if t_stop == None:
        t_stop = -1
    if var_type == 0:
        arr = np.asarray([f[ts][variables[0]][car_var[var_num]] for ts in timesteps[t_start:t_stop]])
    if var_type == 1:
        arr = np.asarray([f[ts][variables[1]][cnv_var[var_num]] for ts in timesteps[t_start:t_stop]])
    if var_type == 2:
        arr = np.asarray([f[ts][variables[2]][cov_var[var_num]] for ts in timesteps[t_start:t_stop]])
    if var_type == 3:
        arr = np.asarray([f[ts][variables[3]][diag_var[var_num]] for ts in timesteps[t_start:t_stop]])
    if var_type == 4:
        arr = np.asarray([f[ts][variables[4]][pert_var[var_num]] for ts in timesteps[t_start:t_stop]])
    if var_type == 6:
        arr = np.asarray(f[timesteps[0]][variables[6]][cell_var[var_num]])
    
    if var_type in range(0,5):
        arr = np.swapaxes(arr,0,3)
        arr = np.swapaxes(arr,1,2)
    else:
        arr = np.swapaxes(arr,0,2)
    return arr

def fft_in_phi(arr,t=None):
    dphi = 2*np.pi/arr.shape[2]
    N_p = arr.shape[2]
    dt_p = dphi
    T_p = dt_p*N_p
    df_p = 1/T_p
    dw_p = 2*np.pi/T_p
    freq_p = np.fft.fftfreq(N_p)

    arr_tilda = []
    if t == None:
        for r in range(arr.shape[0]):
            for theta in range(arr.shape[1]):
                for time in range(arr.shape[3]):
                    arr_tilda.append(np.fft.fft(arr[r,theta,:,time]))
        arr_tilda = np.reshape(np.asarray(arr_tilda),(arr.shape[0],arr.shape[1],arr.shape[3],arr.shape[2]))
        arr_tilda = np.swapaxes(arr_tilda,2,3)
    else:
        for r in range(arr.shape[0]):
            for theta in range(arr.shape[1]):
                arr_tilda.append(np.fft.fft(arr[r,theta,:,t]))
        arr_tilda = np.reshape(np.asarray(arr_tilda),(arr.shape[0],arr.shape[1],arr.shape[2]))
    return arr_tilda


def nm_array(b_hat_rho):
    b_hat_rho = (1/(2*np.pi))*b_hat_rho
    c = (1/(2*np.pi))**2 # Normalization factor
    bnt = [] # n-transform
    for r in range(b_hat_rho.shape[0]):
        for uf in range(b_hat_rho.shape[1]):
            for time in range(b_hat_rho.shape[3]):
                bnt.append(np.fft.fft(b_hat_rho[r,uf,:,time]))
    bnt = np.reshape(np.asarray(bnt),(b_hat_rho.shape[0],b_hat_rho.shape[1],b_hat_rho.shape[3],b_hat_rho.shape[2]))
    bnt = np.swapaxes(bnt,2,3)
    
    bnmt = [] # m-transform
    for r in range(b_hat_rho.shape[0]):
        for fn in range(b_hat_rho.shape[2]):
            for time in range(b_hat_rho.shape[3]):
                bnmt.append(np.fft.fft(bnt[r,:,fn,time]))
    bnmt = np.reshape(np.asarray(bnmt),(b_hat_rho.shape[0],b_hat_rho.shape[2],b_hat_rho.shape[3],b_hat_rho.shape[1]))
    bnmt = np.swapaxes(bnmt,2,3)
    bnmt = np.swapaxes(bnmt,1,2)
    bnmt = c*bnmt
    return bnmt

def intersection_of_horizontal_line_with_flux_surface(psit,psi_min,norm,X0,Z0,init_guess = 3.4):
    """Return the intersection point of a horizontal line from the magnetic axis with the q_95 surface."""
    #init_guess = 3.4 Double Tearing shaped simulation, 2.4 Sawtooth, 3.7 dt-nodiff-64nzd 
    cs1 = plt.contour(X[:,:,0],Z[:,:,0],(psit[:,:]-psi_min)/(norm-psi_min),levels=[0.95],colors="y")
    #for i in range(len(cs1.collections[0].get_paths())): # get the coordinates of q_95 surface
    #    p1 = cs1.collections[0].get_paths()[i]
    #    v1 = p1.vertices
    #    x_l = v1[:,0] 
    #    y_l = v1[:,1]
    v1 = cs1.allsegs[0][0]
    x_l = v1[:,0]
    y_l = v1[:,1]
    
    
    lbind = np.where(v1[:,0] < X0) # left branch indices
    lbt = int(lbind[0][0]) # top of left branch
    lbb = int(lbind[0][-1]) # bottom of left branch
    
     # Choose part of flux surface suitable for interpolation
    if Z0 > y_l[0]: # magnetic axis higher than 1st point of flux surface
        list_fs_x = list(x_l[0:lbt-2]) # pick upper branch. Leave two points to avoid endpoint problems
        list_fs_y = list(y_l[0:lbt-2])
    else:
        list_fs_x = list(x_l[lbb+2:-1]) # pick lower branch
        list_fs_y = list(y_l[lbb+2:-1])

    xmin = min(X[:,int(X.shape[1]/2),0])
    xmax = max(X[:,0,0])
    xhl_list = np.linspace(xmin,xmax,100) # list of points for horizontal line interpolation
    zhl_list = Z0*np.ones(100)

    # Interpolation of the two lines
    flux_surf_interp = scipy.interpolate.interp1d(list_fs_x,list_fs_y,fill_value="extrapolate")
    hline_interp = scipy.interpolate.interp1d(xhl_list,zhl_list,fill_value="extrapolate")

    # defines error that needs to be minimized to find crossing
    def difference(x):
        return np.abs(hline_interp(x)-flux_surf_interp(x))
    # find crossing point: NEEDS INITIAL GUESS
    Xx = scipy.optimize.fsolve(difference,x0=init_guess)
    
    return Xx

def r_psi_list(t,psit,psi_min,norm):
    """Finding the r(psi) relationship via interpolation."""
    psit_int = Grid_Interpolation_Single_Array(psit,tor=True) # Interpolator for psi_toroidal_average.
    if t_dim == 0: # works for single time step
        r_ma,u_ma = np.unravel_index(np.argmin(psit[:,:]),(psit.shape[0],psit.shape[1])) # find index of magnetic axis
        
        X0 = X[r_ma,u_ma,0] # x,z locations of magnetic /lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/300x/pixie3d-n0=300x.scratchaxis
        Z0 = Z[r_ma,u_ma,0]
        
        r_ma_log = CnvNumber2LogicalR(r_ma) # r,u logical coordinates of magnetic axis 
        u_ma_log = CnvNumber2LogicalU(u_ma)
        
        Xx = intersection_of_horizontal_line_with_flux_surface(psit,psi_min,norm,X0,Z0) # intersection of horizontal line with q95
        
        H_line = np.linspace(X0,Xx,100) # horizontal line that connects magnetic axis and q95 surface

        # convert X,Z0 pairs to logical coordinates in the shaped grid
        r_of_hline = [r_ma_log] # start with first position being the magnetic axis. Procedure leaves a systematic error of 2cm.
        u_of_hline = [u_ma_log]
        for x in H_line[1:]:
            r,u = pntCnvInGrid_shaped(x,Z0,r_of_hline[-1],u_of_hline[-1]) # last position is initial guess for new one
            #r,u = pntCnvInGrid_simple_toroidal(x,Z0,r_of_hline[-1],u_of_hline[-1])
            r_of_hline.append(r)
            u_of_hline.append(abs(u))

        # psin values of those points    
        psi_val_list = (psit_int((r_of_hline,u_of_hline))-psi_min)/(norm-psi_min)
        
       
        # R(psin) interpolation
        r_of_psin = interp1d(psi_val_list,r_of_hline, kind = "cubic", fill_value = "extrapolate")

        # Final r(psin) list
        psin_list = np.linspace(0.0,1.0,101)
        r_of_psin_list = r_of_psin(psin_list)
        
    else:# works for multiple time steps
        #print(t)
        r_ma,u_ma = np.unravel_index(np.argmin(psit[:,:,t]),(psit.shape[0],psit.shape[1])) # find index of magnetic axis
        
        X0 = X[r_ma,u_ma,0] # x,z locations of magnetic axis
        Z0 = Z[r_ma,u_ma,0]
        
        r_ma_log = CnvNumber2LogicalR(r_ma) # r,u logical coordinates of magnetic axis 
        u_ma_log = CnvNumber2LogicalU(u_ma)
        #print(r_ma_log)
        
        Xx = intersection_of_horizontal_line_with_flux_surface(psit[:,:,t],psi_min,norm,X0,Z0)
        
        H_line = np.linspace(X0,Xx,100) # horizontal line that connects magnetic axis and q95 surface
        #print(H_line)
        
        # convert X,Z0 pairs to logical coordinates in the shaped grid
        r_of_hline = [r_ma_log] # start with first position being the magnetic axis. Procedure leaves a systematic error of 2cm.
        u_of_hline = [u_ma_log]
        for x in H_line[1:]:
            r,u = pntCnvInGrid_shaped(x,Z0,r_of_hline[-1],u_of_hline[-1]) # last position is initial guess for new one
            #r,u = pntCnvInGrid_simple_toroidal(x,Z0,r_of_hline[-1],u_of_hline[-1])
            if r>0:
                r_of_hline.append(r)
                u_of_hline.append(abs(u) % (2*np.pi))
            else:
                pass

        # psin values of those points    
        psi_val_list = (psit_int((r_of_hline,u_of_hline,t))-psi_min)/(norm-psi_min)
        #psi_val_list = psit_int((r_of_hline,u_of_hline,t))
        #print("R:",r_of_hline)
        #print("U:",u_of_hline)
        #print("Psi:",psi_val_list)
        #plt.clf()
        #plt.figure()
        #plt.scatter(X_int((r_of_hline[5],0,0)),Z_int((r_of_hline[5],0,0)),"*")
        #plt.plot(psi_val_list, r_of_hline)
        #plt.show()
        #print(r_of_hline)
        # R(psin) interpolation
        r_of_psin = interp1d(psi_val_list,r_of_hline, kind = "cubic", fill_value = "extrapolate")
        
        # Final r(psin) list
        psin_list = np.linspace(0.0,1.0,101)
        r_of_psin_list = r_of_psin(psin_list)

    return r_of_psin_list, r_ma_log, u_ma_log

def create_r_psi_list(psit,Bpsq_tor,x0):
    """Functions that create the r(t,psi) arrays via the two above methods."""
    psi_min_list, norm_list = Normalization_numbers(psit,Bpsq_tor,x0)
    R_of_psi = []
    r_ma = []
    u_ma = []
    if t_dim == 0:
        Rs, rmaxis, umaxis = r_psi_list(0,psit,psi_min_list[0],norm_list[0])
        R_of_psi.append(Rs)
        r_ma.append(rmaxis)
        u_ma.append(umaxis)
    else:
        t_range = psit.shape[2]
        for t in range(t_range):
            Rs, rmaxis, umaxis = r_psi_list(t,psit,psi_min_list[t],norm_list[t])
            R_of_psi.append(Rs)
            r_ma.append(rmaxis)
            u_ma.append(umaxis)
        R_of_psi = np.reshape(np.asarray(R_of_psi),(t_range,101))
    return R_of_psi, r_ma, u_ma
            

def norm_Psi2Psi_Init(psin,psi_min,norm):
    """Function that takes a normalized psi value and returns the value of the original psi array."""
    return psin*(norm - psi_min) + psi_min

def CnvNumber2LogicalR(r_num):
    """Convert a point of the numerical grid to logical coordinates, i.e., r in [0,1]."""
    r_dim = float(X.shape[0]-1)
    r_log = np.asarray(r_num)/r_dim
    return r_log

def CnvNumber2LogicalU(u_num):
    """Convert a point of the numerical grid to logical coordinates, i.e., u in [0,2pi]."""
    u_dim = float(X.shape[1]-1)
    u_log = (2*np.pi)*(np.asarray(u_num)/u_dim)
    return u_log


def Axes_of_Interpolation(arr):    
    """Create axes of logical coordinates for intepolating arrays on grid nodes. Takes array for finding time dimension."""
    global r, theta, phi, t, t_dim
    
    r = np.linspace(0.0,1.0,num=X.shape[0]) 
    theta = np.linspace(0.0,2.0*np.pi,num=X.shape[1])
    phi = np.linspace(0.0,2.*np.pi,num=X.shape[2])
    
    if arr.ndim == 4:
        t_dim = arr.shape[3]
        t = np.linspace(0,arr.shape[3]-1,num=arr.shape[3])
    else:
        t_dim = 0
        t = 0

def Grid_for_Cell_Interpolation(arr):    
    """Create axes of logical coordinates for intepolating arrays on grid cells. Takes array defined on cells."""
    global rc, uc, phic 
    
    num_r_cells = arr.shape[0]
    num_u_cells = arr.shape[1]
    num_phi_cells = arr.shape[2]
    dn_r = (1.0/num_r_cells)
    dn_u = ((2.0*np.pi)/num_u_cells)
    
    # Cell-based grid
    rc = np.linspace(0.0+(dn_r/2.0),1.0-(dn_r/2.0),num_r_cells)
    uc = np.linspace(0.0+(dn_u/2.0),2.0*np.pi-(dn_u/2.0),num_u_cells)
    phic = np.linspace(0.0+(dn_u/2.0),2.0*np.pi-(dn_u/2.0),num_phi_cells)
    
def Coordinate_Maps_Interpolations():
    global X_int, Y_int, Z_int
    
    r = np.linspace(0.0,1.0,num=X.shape[0]) 
    theta = np.linspace(0.0,2.0*np.pi,num=X.shape[1])
    phi = np.linspace(0.0,2.*np.pi,num=X.shape[2])
    
    X_int = RegularGridInterpolator((r,theta,phi), X[:,:,:], method='linear', bounds_error=False, fill_value = 0)
    Y_int = RegularGridInterpolator((r,theta,phi), Y[:,:,:], method='linear', bounds_error=False, fill_value = 0)
    Z_int = RegularGridInterpolator((r,theta,phi), Z[:,:,:], method='linear', bounds_error=False, fill_value = 0)
    
def Grid_Interpolation_Single_Array(arr,tor=False):
    """Interpolate single array on logical coordinate grid nodes."""
    assert arr.shape[0] == len(r), f"r-node dimension is not the same as array's r-dimension."
    assert arr.shape[1] == len(theta), f"theta-node dimension is not the same as array's theta-dimension."
    if tor == True:
        try:
            if arr.ndim == 3:
                assert arr.shape[2] == len(t), f"time dimension of array is inconsistent with reference."
                arr_int = RegularGridInterpolator((r,theta,t), arr[:,:,:], method='linear', bounds_error=False, fill_value = 0)
            if arr.ndim == 2:
                arr_int = RegularGridInterpolator((r,theta), arr[:,:], method='linear', bounds_error=False, fill_value = 0)
        except:
            print("Array has prohibited number of dimensions")
    else:
        assert arr.shape[2] == len(phi), f"phi-node dimension is not the same as array's phi-dimension."
        try:
            if arr.ndim == 4:
                assert arr.shape[3] == len(t), f"time dimension of array is inconsistent with reference."
                arr_int = RegularGridInterpolator((r,theta,phi,t), arr[:,:,:,:], method='linear', bounds_error=False, fill_value = 0)
            if arr.ndim == 3:
                arr_int = RegularGridInterpolator((r,theta,phi), arr[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        except:
            print("Array has prohibited number of dimensions")
    
    return arr_int

def Grid_Cell_Interpolation_Single_Array(arr,tor=False):
    """Interpolate single array on logical coordinate grid cells."""
    assert arr.shape[0] == len(rc), f"r-cell dimension is not the same as array's r-dimension."
    assert arr.shape[1] == len(uc), f"theta-cell dimension is not the same as array's theta-dimension."
    if tor == True:
        try:
            if arr.ndim == 3:
                assert arr.shape[2] == len(t), f"time dimension of array is inconsistent with reference."
                arr_int = RegularGridInterpolator((rc,uc,t), arr[:,:,:], method='linear', bounds_error=False, fill_value = 0)
            if arr.ndim == 2:
                arr_int = RegularGridInterpolator((rc,uc), arr[:,:], method='linear', bounds_error=False, fill_value = 0)
        except:
            print("Array has prohibited number of dimensions")
    else:
        assert arr.shape[2] == len(phic), f"phi-cell dimension is not the same as array's phi-dimension."
        try:
            if arr.ndim == 4:
                assert arr.shape[3] == len(t), f"time dimension of array is inconsistent with reference."
                arr_int = RegularGridInterpolator((rc,uc,phic,t), arr[:,:,:,:], method='linear', bounds_error=False, fill_value = 0)
            if arr.ndim == 3:
                arr_int = RegularGridInterpolator((rc,uc,phic), arr[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        except:
            print("Array has prohibited number of dimensions")
    
    return arr_int               
        
def MeshGrids_Creation():    
    """Creating grids of logical coordinates for evaluation of interpolating arrays."""
    global RRRI, TTTI, PPPI,TIM, RI, TI, PI, RRI, TTI, TR, TU, TT
    
    # 2D-3D-4D grids
    RRRI,TTTI,PPPI,TIM = np.meshgrid(r,theta,phi,t,indexing='ij') # 4D
    RI,TI,PI = np.meshgrid(r,theta,phi,indexing='ij') # 3D
    TR,TU,TT = np.meshgrid(r,theta,t,indexing='ij')
    RRI,TTI = np.meshgrid(r,theta,indexing='ij') # 2D
    
def C2N_Evaluation(arr_int,arr,tor=False):
    """Converts cell-based arrays to node-based. If tor switch is turned on, we assume that in the array, the phi
    dimension is missing. Otherwise, we assume it's there."""
    if tor == False:
        if arr.ndim == 4:
            arrNode = arr_int((RRRI,TTTI,PPPI,TIM))
        if arr.ndim == 3:
            arrNode == arr_int((RI,TI,PI))
        if arr.ndim == 2:
            arrNode = arr_int((RRI,TTI))
    if tor == True:
        if arr.ndim == 3:
            arrNode = arr_int((TR,TU,TT))
        if arr.ndim == 2:
            arrNode = arr_int((RRI,TTI))
    return arrNode
        
def units():
    """Calculates units."""
    unitR = (-X[0,0,0]+X[-1,0,0])/(X.shape[0]) # Units of steps of minor radius
    unitU = 2*np.pi/(X.shape[1]-1) # In rad (first and last point is the same)
    
    return unitR, unitU

def pntCnvInGrid_simple_toroidal(x_co,z_co,initial_guess_r = 0.8,initial_guess_u = 4.2):
    """Convert an x,y point in the numerical grid. Works only for toroidal computational domain. 
    Default initial guess arguments work for the Bonfiglio sawtooth simulation X-point.
    """
    def shape_fun(x,*In):
        xx = In[0]
        zz = In[1]
        return [xx-1.65-x[0]*np.cos(x[1]), zz-x[0]*np.sin(x[1])]
    params = (x_co,z_co)
    initial_guess = [initial_guess_r,initial_guess_u]
    xcross = scipy.optimize.fsolve(shape_fun,x0=initial_guess,args=params)
    return xcross[0],xcross[1]

def pntCnvInGrid_shaped(R,Z,initial_guess_r=0.9,initial_guess_u=4.5):
    """Converts x,z points in the r,theta logical grid. Based on the plasma shaping function:
    R(r,u) = Ro + a*r*cos[u+sin^(-1)(delta*r^2*sinu)], Z(r,u)=Zo + a*k*r*sin[u+zeta*r^2*sin(2u)]. 
    Initial guess works for the X-point of the shaped double tearing simulations.
    """
    alpha = 2.24 # minor radius
    alpha_pixie = 2.18095#2.1895 # Pixie3d normalization parameter
    def shape_func(x,*In_par):
        R=In_par[0]
        Z=In_par[1]
        Ro = 6.219577546
        Zo = 0.5143555944
        k = 1.9
        delta = 0.6
        zeta = 0.06
        return [R-Ro-alpha*x[0]*np.cos(x[1]+np.arcsin(delta*(x[0]**2)*np.sin(x[1]))),Z-Zo-alpha*k*x[0]*np.sin(x[1]+zeta*(x[0]**2)*np.sin(2*x[1]))]
    initial_guess=[initial_guess_r,initial_guess_u] # visual inspection is needed to determind initial guess
    params = (R*alpha_pixie,Z*alpha_pixie) # shape function takes unnormalized values
    xcross = scipy.optimize.fsolve(shape_func,x0=initial_guess,args=params,xtol=1.0e-12)
    r_log = xcross[0] # convert r to [0,1] range
    u_log = xcross[1]
    return r_log,u_log
    

def Calculation_of_Units_and_Sizes():    
    """Units in R and theta."""
    global unitR, unitU, r_dim, u_dim, fi_dim, t_dim
    unitR, unitU = units()
    r_dim = X.shape[0]
    u_dim = X.shape[1]
    fi_dim = X.shape[2]
    print("Units and sizes calculated.")
        

def findXpoint(Bpsq_tor,t,x0):
    """Function that locates the X-point as the minimum of the Bp^2. Requires visual inspection for initial guess. 
    Initial guess is given as coordinates in logical space."""
    Bpsq_int = Grid_Interpolation_Single_Array(Bpsq_tor[:,:,t],tor=True) # Interpolator for each time step.
    Xpnt = scipy.optimize.fmin(Bpsq_int,x0)
    return Xpnt

def magneticAxisAndXpoint(psit,Bpsq_tor,t,x0):
    """Calculates values of flux at magnetic axis and X-point. Requires visual inspection for the X-point initial guess. 
    """
    psit_int = Grid_Interpolation_Single_Array(psit,tor=True) # Interpolator for psi_toroidal_average.
    if t_dim == 0:
        psi_min = np.amin(psit[:,:]) # find min of flux
        Xpnt = findXpoint(Bpsq_tor,0,x0) # location of X-point
        norm = psit_int((Xpnt[0],Xpnt[1])) # psi value at X-point (normalization constant)
    else:
        psi_min = np.amin(psit[:,:,t]) # find min of flux 
        Xpnt= findXpoint(Bpsq_tor,t,x0) # location of X-point
        norm = psit_int((Xpnt[0],Xpnt[1],t)) # psi value at X-point (normalization constant)
    return psi_min,norm # norm comes out as 1-dimensional array


def Normalization_numbers(psit,Bpsq_tor,x0):
    """Function that calculates Normalization numbers for each time step.
    dt-nodiff-64nzd: x0=(0.9,4.44)
    """
    psi_min_list = []
    norm_list = []
    if t_dim == 0:
        psi_min, norm = magneticAxisAndXpoint(psit,Bpsq_tor,0,x0)
        psi_min_list.append(psi_min)
        norm_list.append(norm)
    else:
        t_range = psit.shape[2]
        for t in range(t_range):
            psi_min, norm = magneticAxisAndXpoint(psit,Bpsq_tor,t,x0)
            psi_min_list.append(psi_min)
            norm_list.append(norm)
            #if t%10==0:
                #print("t=",t,"is done")
            
    return psi_min_list,norm_list

def x_array(rs,us):
    """Returns X-coordinates from logical ones."""
    X = []
    for r in rs:
        for u in us:
            X.append(X_int((r,u,0)))
    X = np.reshape(X,(rs.shape[0],us.shape[0]))
    return X

def z_array(rs,us):
    """Returns Z-coordinates from logical ones."""
    Z = []
    for r in rs:
        for u in us:
            Z.append(Z_int((r,u,0)))
    Z = np.reshape(Z,(rs.shape[0],us.shape[0]))
    return Z
    
def blob_coordinates(ro,uo,radius,N_radial,N_poloidal):
    """Function that creates logical coordinates of blob of initial conditions.
    ro,uo: center of blob coordinates in logical grid. Found by visual inspection based on where we want to
    place the blob.
    radius: radial extent of the blob in real space.
    N_radial: number of radial steps to go from center to edge of blob.
    N_poloidal: number of poloidal steps in each radial section.
    Return the radial and poloidal coordinate lists of each blob point."""
    dr = 1/X.shape[0] # dr-step that produces the dx
    dx,_ = units() # dx-step produced by dr
    steps = radius/dx # number of dx-steps to make the radius of the blob
    Ro = X_int((ro,uo,0)) # Blob's center in real space
    Zo = Z_int((ro,uo,0))
    r_sblob = [] # Initialize blob points in logical space
    u_sblob = []
    radius = steps*dr # radial extent of the blob in logical space
    rlens = np.arange(0,radius,radius/N_radial) # grid of blob points in logical space, centered at (0,0)
    ulens = np.arange(0,2*np.pi,2*np.pi/N_poloidal)
    for r in rlens: # lists of blob coordinates centered at (0,0)
        for u in ulens:
            r_sblob.append(r)
            u_sblob.append(u)
    ru_sblob = zip(r_sblob,u_sblob) # iterator (only evaluated once!)
    
    X_sblob = [] # Initialize blob points in real space
    Z_sblob = []
    for coor in ru_sblob: # Transport blob points to center them at (Ro,Zo)
        X_sblob.append(Ro + X_int((coor[0],coor[1],0)) - X[0,0,0])
        Z_sblob.append(Zo + Z_int((coor[0],coor[1],0)) - Z[0,0,0])
    XZ_sblob = zip(X_sblob,Z_sblob) # iterator
    
    r_transblob = [] # Initialize transported blob points in logical space
    u_transblob = []
    for coor in XZ_sblob: # transform real space points back to logical space
        r_tr, u_tr = pntCnvInGrid_shaped(coor[0],coor[1],ro,uo)
        r_transblob.append(r_tr)
        u_transblob.append(u_tr)
        
    return r_transblob, u_transblob    
################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################    
################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################        
################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################    
############################################################################################################################################################################################################################################################################################################################################################################################################################################   DEPRECATED FUNCTIONS  ###########################################################  ################################################################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################    
################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################    
################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################    

"""
# Calculating the psi for all times
A = A_phi()
A = -A

A_int = RegularGridInterpolator((r,theta,t), A[:,:,:], method='linear', bounds_error=False, fill_value = 0)

PN = norm_psi()

PN_int = RegularGridInterpolator((r,theta,t), PN[:,:,:], method='linear', bounds_error=False, fill_value = 0)
PN_grid = PN_int((TR,TU,TT))

# Derivatives and interpolation
dpdr_t, dpdt_t = derivative_calculation_t()
dpsidr_t_int = RegularGridInterpolator((r,theta,t), dpdr_t[:,:,:], method='linear', bounds_error=False, fill_value = 0)
dpsidtheta_t_int = RegularGridInterpolator((r,theta,t), dpdt_t[:,:,:], method='linear', bounds_error=False, fill_value = 0)
print("Psi derivatives over all times done...")
"""
"""
"Function that reconstructs later q-profiles."
def q_prof_t():
    
    q = []
    psi = 0.1
 
    for t in range(PN.shape[2]):
        r_fs, u_fs = fs_locator_t((psi),t)
        rm, u, uf, l = integrator_t(r_fs,u_fs,t)
        q_redef = uf[-1]/(2*np.pi)
        uf_rsc = uf/q_redef
        q.append(q_redef)
    
    return q

def q_at_zero():
    
    psi_list = np.linspace(0.001,1,100)
    q = []
    
    for psi in psi_list:
        r_fs, u_fs = fs_locator_t((psi),108)
        #print("R:", r_fs)
        rm, u, uf, l = integrator_t(r_fs,u_fs,108)
        q_redef = uf[-1]/(2*np.pi)
        uf_rsc = uf/q_redef
        q.append(q_redef)
    
    return q, psi_list

def q_cleaning():
    "Function that cleans noisy data from q-profile leaving only points on the sawtooth cycle."
    global q,t
    q_mag_axis = q[2,0,0,:]
    tolerance = 0.006
    # Extracting points on the sawtooth curve
    q_clean = []
    dim = len(q_mag_axis)
    for i in range(dim-1):
        qi = q_mag_axis[i]
        if i>0:
            qim = q_mag_axis[i-1]
        else:
            qim = qi
        if i<dim:
            qip = q_mag_axis[i+1]
        else:
            qip = qi
        if np.abs(qi-qim)<tolerance or np.abs(qi-qip)<tolerance:
            q_clean.append(qi)
    
    # Splining the new q
    reduced_dim = len(q_clean)
    time = np.linspace(0,len(t),reduced_dim)
    q_interp = interp1d(time,q_clean)
    q_new = q_interp(t)
    
    # Plotting the q(0)
    plt.plot(t,q_new)
    plt.xlabel(r'$\frac{t}{\tau_A}$')
    plt.ylabel(r'q(0)')
    plt.show()

def A_phi():
    "Function that returns the phi component of the vector potential for all times."    
    B1_tor = np.mean(B1,2)
    B2_tor = np.mean(B2,2)
    r_dim = B1.shape[0]
    u_dim = B1.shape[1]
    t_dim = B1.shape[3]
    A_u = []
  
    mRadius = [x*unitR for x in range(r_dim)] # Minor radius
    mROv2 = np.power(1./np.asarray(mRadius),0) # One over r^2
    
    MRadius = [] # Major radius
    for r in range(r_dim):
        for u in range(u_dim):
            MRadius.append(X[0,0,0] + (r*unitR)*np.cos(unitU*u))
    MRadius = np.reshape(np.asarray(MRadius),(r_dim,u_dim))
    MROv1 = 1./MRadius
    MROv1 = MROv1[:,:,np.newaxis]
    
    B1Int = np.multiply(MROv1,B1_tor)
    B2Int = np.multiply(MROv1,B2_tor)
    #B1Int = MROv1*B1_tor
    #B2Int = MROv1*B2_tor
    
    for t in range(t_dim):
        # Locate magnetic axis
        A = np.where(np.abs(B2_tor[0:10,0,t]) == np.amin(np.abs(B2_tor[0:10,0,t])))
        r_ma = A[0][0]
        for u in range(u_dim):
            for r in range(r_dim):
                if r == r_ma:
                #if r == 0:
                    A_u.append(0)
                else:
                    summa_u = np.sum(B1_tor[r,1:u-1,t])
                    #summa_r = np.sum(B2_tor[1:r-1,u,t])
                    summa_r = np.sum(B2_tor[1:r-1,0,t])
                    #summa_u = np.sum(MROv1[r,1:u-1]*B1_tor[r,1:u-1,t])
                    #summa_r = np.sum(MROv1[1:r-1,u]*B2_tor[1:r-1,u,t])
                    #integral = -(1./r_dim)*(((MROv1[1,u]*B2_tor[1,u,t])/2.) + summa_r + ((MROv1[r,u]*B2_tor[r,u,t])/2.)) + (2.*np.pi/u_dim)*((MROv1[r,0]*B1_tor[r,0,t]/2.) + summa_u + (MROv1[r,u]*B1_tor[r,u,t]/2.))
                    #A_u.append(-integral*2*np.pi*MRadius[r,u])
                    integral = -((X[-1,0,0]-X[0,0,0])/(r_dim))*(((B2_tor[1,0,t])/2.) + summa_r + ((B2_tor[r,0,t])/2.)) + (2.*np.pi/(u_dim))*((B1_tor[r,0,t]/2.) + summa_u + (B1_tor[r,u,t]/2.))
                    A_u.append(integral)
                    
    A_u = np.reshape(A_u,(t_dim,u_dim,r_dim))
    A_u = np.swapaxes(A_u,0,2)
    
    return A_u        

def Load_Presaved_Arrays():
    "Loads a pre-saved array. We treat this all time-dependent quantities."
    global B1,B_1,B2,B_2,B3,B_3,q
    B1 = np.load("/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/python_arrays/B1.npy")
    #B_1 = np.load("/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/python_arrays/B_1.npy")
    B2 = np.load("/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/python_arrays/B2.npy")
    #B_2 = np.load("/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/python_arrays/B_2.npy")
    B3 = np.load("/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/python_arrays/B3.npy")
    #B_3 = np.load("/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/python_arrays/B_3.npy")
    q = np.load("/net/scratch3/giannis_kx/pixie3d/tests/bonfiglio_div_tok/python_arrays/q.npy")

def Data_Load_And_Preparation():
    load_data()
    ax_rearrange()
    f2c_transformations()
    fix_periodicity()
    
def Grids_And_Interpolations():
    "Creates grids, interpolations and evaluations. To be executed in notebook."
    # Should be called only once. Can go wrong with the non-consistency of face centered and cell centered data. If it goes wrong
    # check which variable has the wrong dimension and fix it with the f2c function.
    Axes_of_Interpolation()
    Grid_Interpolations()
    Grids_Creation()
    Grid_Evaluations()

def fs_locator_t(psi,t):
    "Function that locates the flux surface coordinates for later times."    
    fs = np.asarray(measure.find_contours(PN_grid[:,:,t], psi))
    raxis_fs=[]
    uaxis_fs=[]
    if psi<0 or psi>1:
        raxis_fs = np.nan
        uaxis_fs = np.nan
    else:
        try:
            for i in range(fs.shape[1]):
                raxis_fs.append(fs[0][i][0])
                uaxis_fs.append(fs[0][i][1])
        except IndexError:
            try:
                for i in range(fs[0].shape[0]):
                    raxis_fs.append(fs[0][i][0])
                    uaxis_fs.append(fs[0][i][1])
            except IndexError:
                raxis_fs = np.nan
                uaxis_fs = np.nan
    warnings.resetwarnings()
    
    return raxis_fs,uaxis_fs

def dPsidr2prec(r,theta,phi):
    "Partial w.r.t r."    
    global Psi2_grid
    dr = 1./(2*Psi.shape[0])
    psi = Psi2_grid[:,theta,phi]
    
    psiplus = np.roll(psi,-1)
    psiminus = np.roll(psi,1)
    dpsidr = (psiplus-psiminus)/(2.*dr)
    dpsidr[0] = (psi[1]-psi[0])/dr
    dpsidr[-1] = (psi[-1]-psi[-2])/dr
    
    return dpsidr[r]


def dPsidr_t(r,theta,t):
    "Partial w.r.t r."    
    dr = 1./PN.shape[0]
    psi = PN[:,theta,t]
    
    psiplus = np.roll(psi,-1)
    psiminus = np.roll(psi,1)
    dpsidr = (psiplus-psiminus)/(2.*dr)
    dpsidr[0] = (psi[1]-psi[0])/dr
    dpsidr[-1] = (psi[-1]-psi[-2])/dr
    
    return dpsidr[r]

def dPsidtheta2prec(r,theta,phi):
    "Partial w.r.t theta (double precision)."    
    global Psi2_grid
    dtheta = (2.*np.pi)/(2*Psi.shape[1])
    psi = Psi2_grid[r,:,phi]
    
    psiplus = np.roll(psi,-1)
    psiminus = np.roll(psi,1)
    dpsidtheta = (psiplus-psiminus)/(2.*dtheta)
    dpsidtheta[0] = (psi[1]-psi[0])/dtheta
    dpsidtheta[-1] = (psi[-1]-psi[-2])/dtheta
    
    return dpsidtheta[theta]


def dPsidtheta_t(r,theta,t):
    "Partial w.r.t theta."    
    dtheta = (2.*np.pi)/PN.shape[1]
    psi = PN[r,:,t]
    
    psiplus = np.roll(psi,-1)
    psiminus = np.roll(psi,1)
    dpsidtheta = (psiplus-psiminus)/(2.*dtheta)
    dpsidtheta[0] = (psi[1]-psi[0])/dtheta
    dpsidtheta[-1] = (psi[-1]-psi[-2])/dtheta
    
    return dpsidtheta[theta]

def derivative_calculation_2_prec():
    "Function that calculates both derivatives (double precision)."    
    dpdr = []
    dpdt = []
    
    for r in range(2*Psi.shape[0]):
        for theta in range(2*Psi.shape[1]):
            for phi in range(2*Psi.shape[2]):
                dpdr.append(dPsidr2prec(r,theta,phi))
                dpdt.append(dPsidtheta2prec(r,theta,phi))
    
    return np.reshape(np.asarray(dpdr),(2*Psi_grid.shape[0],2*Psi_grid.shape[1],2*Psi_grid.shape[2])), np.reshape(np.asarray(dpdt),(2*Psi_grid.shape[0],2*Psi_grid.shape[1],2*Psi_grid.shape[2]))


def derivative_calculation_t():
    "Function that calculates both derivatives."    
    dpdr = []
    dpdt = []
    
    for r in range(PN.shape[0]):
        for theta in range(PN.shape[1]):
            for t in range(PN.shape[2]):
                dpdr.append(dPsidr_t(r,theta,t))
                dpdt.append(dPsidtheta_t(r,theta,t))
    
    return np.reshape(np.asarray(dpdr),(PN.shape[0],PN.shape[1],PN.shape[2])), np.reshape(np.asarray(dpdt),(PN.shape[0],PN.shape[1],PN.shape[2]))

def flux_surf_ang_2_prec(l,state_vec):
    "Differential equations for arc lenth of a flux surface and flux angle. Inputs are given in terms of logical coordinates. Using double precision derivatives. Even then, it is failing below psi=0.009."    
    r = state_vec[0]
    u = state_vec[1]
    uf = state_vec[2]
    phi = state_vec[3]
   
    drdl = -dpsidtheta_2_prec_int((r,u,phi))/(np.sqrt(np.power(dpsidtheta_2_prec_int((r,u,phi)),2) \
                                            + np.power(r,2)*np.power(dpsidr_2_prec_int((r,u,phi)),2)))
    dthetadl = dpsidr_2_prec_int((r,u,phi))/(np.sqrt(np.power(dpsidtheta_2_prec_int((r,u,phi)),2) \
                                            + np.power(r,2)*np.power(dpsidr_2_prec_int((r,u,phi)),2)))
    
    dufdl = B3_int((r,u,phi,0))/(np.sqrt(np.power(dpsidtheta_2_prec_int((r,u,phi)),2) \
                                            + np.power(r,2)*np.power(dpsidr_2_prec_int((r,u,phi)),2)))
    dphidl = 0
    
    return[drdl,dthetadl,dufdl,dphidl]



def flux_surf_ang_t(l,state_vec):
    "Differential equations for arc lenth of a flux surface and flux angle. Inputs are given in terms of logical coordinates."    
    r = state_vec[0]
    u = state_vec[1]
    uf = state_vec[2]
    time = state_vec[3]
    
    drdl = -dpsidtheta_t_int((r,u,time))/(np.sqrt(np.power(dpsidtheta_t_int((r,u,time)),2) \
                                            + np.power(r,2)*np.power(dpsidr_t_int((r,u,time)),2)))
    dthetadl = dpsidr_t_int((r,u,time))/(np.sqrt(np.power(dpsidtheta_t_int((r,u,time)),2) \
                                            + np.power(r,2)*np.power(dpsidr_t_int((r,u,time)),2)))
    
    dufdl = B3_tor_int((r,u,time))/(np.sqrt(np.power(dpsidtheta_t_int((r,u,time)),2) \
                                            + np.power(r,2)*np.power(dpsidr_t_int((r,u,time)),2)))
    dtimedl = 0
    
    return[drdl,dthetadl,dufdl, dtimedl]

def integrator_t(rs,us,time):
    "Integrating to find r(l), theta(l) around a flux surface. Returning r's, theta's, theta_f's and l's."    
    rs = np.asarray(rs)/float(PN.shape[0]-1)
    us = ((2*np.pi)/float(PN.shape[1]-1))*np.asarray(us)
    solver = scipy.integrate.ode(flux_surf_ang_t).set_integrator('dopri5')
    
    uf0 = 0.0
    dl = 0.01
            
    state_vec0,l0 = [rs[0],0.0,uf0,time],0.0
    solver.set_initial_value(state_vec0,l0)
    y,t = [], []
    warnings.filterwarnings("ignore")
    
    while solver.successful() and solver.y[1]<2*np.pi:
        solver.set_initial_value([solver.y[0],solver.y[1],solver.y[2],solver.y[3]],solver.t)
        solver.integrate(solver.t+dl)
        y.append(solver.y)
        t.append(solver.t)
    warnings.resetwarnings()
    y = np.array(y)
    t = np.array(t)
    
    return y[:,0],y[:,1],y[:,2],t

def norm_psi():
    "Function that normalizes the psi array for all times."
    B2_tor = np.mean(B2,2)
    PN = np.zeros((A.shape[0],A.shape[1],A.shape[2]), dtype = A.dtype)
    for t in range(A.shape[2]):
        psi_min = np.amin(A[:,:,t])
        # Finding X-point
        Xx, Xy = findXpoint(t)
        x_int, u_int = pntCnvInGrid(Xx,Xy)
        # Normalizing value from the X-point
        norm = A_int((x_int,u_int,t)) 
        
        pn = (A[:,:,t] - psi_min)/norm
        PN[:,:,t] = pn
    
    return PN   


def load_data():
    "Loading the data"
    global X, Y, Z, Psi
    X, Y, Z, Psi = pixieload(filepath)
    print("Data loaded...")
    
    
def ax_rearrange():
    "Before we start anything we are putting the arrays in the form r,theta,phi,t and we fix the periodicity."
    global X, Y, Z, B1, B2, B3, Bx, By, Bz, B_1, B_2, B_3, Psi, Pr, q
    Psi = np.swapaxes(Psi,0,2)
    X = np.swapaxes(X,0,2)
    Y = np.swapaxes(Y,0,2)
    Z = np.swapaxes(Z,0,2)
    B1 = np.swapaxes(B1,0,3)
    B2 = np.swapaxes(B2,0,3)
    B3 = np.swapaxes(B3,0,3)
    B1 = np.swapaxes(B1,1,2)
    B2 = np.swapaxes(B2,1,2)
    B3 = np.swapaxes(B3,1,2)
    #Pr = np.swapaxes(Pr,0,3)
    q = np.swapaxes(q,0,3)
    #Pr = np.swapaxes(Pr,1,2)
    q = np.swapaxes(q,1,2)
    #Bx = np.swapaxes(Bx,0,3)
    #By = np.swapaxes(By,0,3)
    #Bz = np.swapaxes(Bz,0,3)
    #Bx = np.swapaxes(Bx,1,2)
    #By = np.swapaxes(By,1,2)
    #Bz = np.swapaxes(Bz,1,2)
    #B_1 = np.swapaxes(B_1,0,3)
    #B_2 = np.swapaxes(B_2,0,3)
    #B_3 = np.swapaxes(B_3,0,3)
    #B_1 = np.swapaxes(B_1,1,2)
    #B_2 = np.swapaxes(B_2,1,2)
    #B_3 = np.swapaxes(B_3,1,2)
    print("Axes swapped...")

def f2c_transformations():   
    "Transforming face centered to cell centered data"
    global Bx, By, Bz, Psi, Pr, q
    #Bx = f2c(Bx)
    #By = f2c(By)
    #Bz = f2c(Bz)
    Psi = f2c(Psi)
    q = f2c(q)
    #Pr = f2c(Pr)
    print("f2c transormations done...")

def fix_periodicity():    
    "Fixing the periodicity of the arrays"
    global X, Y, Z, B1, B2, B3, Bx, By, Bz, B_1, B_2, B_3, Psi, Pr, q
    Psi = periodicity(Psi)
    X = periodicity(X)
    Y = periodicity(Y)
    Z = periodicity(Z)
    B1 = periodicity(B1)
    B2 = periodicity(B2)
    B3 = periodicity(B3)
    #B_1 = periodicity(B_1)
    #B_2 = periodicity(B_2)
    #B_3 = periodicity(B_3)
    #Pr = periodicity(Pr)
    q = periodicity(q)
    #Bx = periodicity(Bx)
    #By = periodicity(By)
    #Bz = periodicity(Bz)
    print("Periodicity of the arrays fixed...")

def Double_Precision_Evaluation_Grid():
    "Create grid of logical coordinates for evaluating psi array with double precision."
    global RR, THTH, PHPH, rr, thth, phph    
    rr = np.linspace(0.0,1.0,num=2*Psi.shape[0])
    thth = np.linspace(0.0,2.*np.pi,num=2*Psi.shape[1])
    phph = np.linspace(0.0,2.*np.pi,num=2*Psi.shape[2])
    RR,THTH,PHPH = np.meshgrid(rr,thth,phph,indexing='ij')
    
def Double_Precision_Grid_Evaluation():
    "Evaluating psi with double precision."
    global Psi2_grid
    Psi2_grid = Psi_int((RR,THTH,PHPH))

def f2c_angles(arr):
    "Function the reduces the array shape just in the angle coordinates."
    rp = np.roll(arr,1,axis=0)
    up = np.roll(arr,1,axis=1)
    pp = np.roll(arr,1,axis=2)
    
    c_arr = (rp+arr+up+arr+pp+arr)/6.
    
    if arr.ndim == 4:    
        c_arr = c_arr[:,:-1,:-1,:] 
    
    if arr.ndim ==3: 
        c_arr = c_arr[:,:-1,:-1]
    
    return c_arr

def Jacobian():
    global Jac
    R_maj = np.sqrt(np.power(X[:,:,0],2)+np.power(Y[:,:,0],2))
    r_min = np.sqrt(np.power((R_maj-R0),2) + np.power(Z[:,:,0],2))
    Jac = np.multiply(R_maj,r_min)

def dAdr(r,theta,A):
    "Partial of (2d) A w.r.t r."    
    dr = 1./(A.shape[0])
    A_r = A[:,theta]
    
    A_r_plus = np.roll(A_r,-1)
    A_r_minus = np.roll(A_r,1)
    dAdr = (A_r_plus-A_r_minus)/(2.*dr)
    dAdr[0] = (A_r[1]-A_r[0])/dr
    dAdr[-1] = (A_r[-1]-A_r[-2])/dr
    
    return dAdr[r]

def dAdtheta(r,theta,A):
    "Partial of (2d) A w.r.t theta."    
    dtheta = (2.*np.pi)/(A.shape[1])
    A_u = A[r,:]
    
    A_u_plus = np.roll(A_u,-1)
    A_u_minus = np.roll(A_u,1)
    dAdtheta = (A_u_plus-A_u_minus)/(2.*dtheta)
    dAdtheta[0] = (A_u[1]-A_u[0])/dtheta
    dAdtheta[-1] = (A_u[-1]-A_u[-2])/dtheta
    
    return dAdtheta[theta]

def derivative_calculation_2d(A):
    "Function that calculates both derivatives of 2d array, A."    
    dadr = []
    dadt = []
    
    for r in range(A.shape[0]):
        for theta in range(A.shape[1]):
            dadr.append(dAdr(r,theta,A))
            dadt.append(dAdtheta(r,theta,A))
    
    return np.reshape(np.asarray(dadr),(A.shape[0],A.shape[1])), np.reshape(np.asarray(dadt),(A.shape[0],A.shape[1]))

def derivatives_and_interpolations_2_prec():
    global dpdr2pr, dpdt2pr, dpsidr_2_prec_int, dpsidtheta_2_prec_int
    "Evaluating psi derivatives and their interpolation with double precision"
    dpdr2pr, dpdt2pr = derivative_calculation_2_prec()
    dpsidr_2_prec_int = RegularGridInterpolator((rr,thth,phph), dpdr2pr[:,:,:], method='linear', bounds_error=False, fill_value = 0)
    dpsidtheta_2_prec_int = RegularGridInterpolator((rr,thth,phph), dpdt2pr[:,:,:], method='linear', bounds_error=False, fill_value = 0)
    print("Double precision Psi derivatives done...")    
    
def Interpolating_Grid(arr):
    "Creates the interpolating grid for an array. Created for parallel program where each toroidal plane is given to a different
    processor and the arrays are in the form [r,theta,time]."
    r = np.linspace(0.0,1.0,num=arr.shape[0])
    theta = np.linspace(0.0,2.*np.pi,num=arr.shape[1])
    t = np.linspace(0,arr.shape[2],num=arr.shape[2])
    return r,theta,t
    
def Interpolate_Array(arr):
    "Interpolates an array on logical grid."
    r,theta,t = Interpolating_Grid(arr)
    arr_int = RegularGridInterpolator((r,theta,t), arr[:,:,:], method='linear', bounds_error=False, fill_value = 0)
    return arr_int  
    
def new_coords2(psi_dict, arr, time):
    psi_list = []
    uf_list = []
    arr_int_list = []
    arr_int = Interpolate_Array(arr)
    for key in list(psi_dict.keys())[:]:
        l_min = psi_dict[key[0],str(0)]['l_u'](0.001)
        l_max = psi_dict[key[0],str(0)]['l_u'](6.28)
        ls = np.linspace(l_min,l_max,30)
        for l in ls:
            r = psi_dict[key[0],str(0)]['r_l'](l)
            u = psi_dict[key[0],str(0)]['u_l'](l)
            uf = psi_dict[key[0],str(0)]['uf_l'](l)
            value = arr_int((r,u,time))
            psi_list.append((float(key[0])-psi_min)/(norm-psi_min))
            uf_list.append(uf)
            arr_int_list.append(value)
    
    triangObj = Triangulation(psi_list,uf_list)
    tcp = LinearTriInterpolator(triangObj, arr_int_list)
    
    return tcp

def pntCnvInGrid2(x,y):
    "Converts x,y point in the numerical grid. Works for any type of computational domain but sometimes it fails and then it
    returns an absurdly high value to alert the program not to use these numbers."
    tol = 0.0001
    XY1 = np.array([],dtype=np.int64)
    XY2 = np.array([],dtype=np.int64)
    try:
        while XY1.size==0 or XY2.size==0: # keep increasing tolerance until you find some points
            XY1 = np.where(np.abs(X[:,:,0]-x)<tol)
            XY2 = np.where(np.abs(Z[:,:,0]-y)<tol)
            tol = 10*tol
    except:
        while XY1[0].size==0 or XY2[0].size==0: # keep increasing tolerance until you find some points
            XY1 = np.where(np.abs(X[:,:,0]-x)<tol)
            XY2 = np.where(np.abs(Z[:,:,0]-y)<tol)
            tol = 10*tol
        
    r_grid_pnt = 1000 # Initialize to avoid not finding points
    u_grid_pnt = 1000
    
    for i in range(len(XY1[0])):
        for j in range(len(XY2[0])):
            if XY1[0][i]==XY2[0][j] and XY1[1][i]==XY2[1][j]:
                r_grid_pnt = XY1[0][i]
                u_grid_pnt = XY1[1][i]
    
    return r_grid_pnt,u_grid_pnt


def r_psi_list2(t,psi,psi_min,norm):
    "Finding the r(psi) relationship via contour identification."
    r = []
    psi_list = np.linspace(0.05,0.99,100)
    for psin in psi_list:
        psi_value = norm_Psi2Psi_Init(psin,psi_min,norm) 
        r_fs, _ = fs_locator(psi,(psi_value),t)
        try:
            r.append(CnvNumber2LogicalR(r_fs[0])) # convert r's to logical grid
        except:
            r.append(np.nan)
    return r

def create_r_psi_list2(psi,B1,B2):
    psi_min_list,norm_list = Normalization_numbers(psi,B1,B2)
    R_of_psi = []
    for t in range(t_dim):
        Rs = r_psi_list(t,psi,psi_min_list[t],norm_list[t])
        R_of_psi.append(Rs)
    R_of_psi = np.reshape(np.asarray(R_of_psi),(t_dim,100))
    return R_of_psi

def n2c(arr):
    "Function that transforms node-centered data to cell-centered data" 
    if arr.ndim == 4:
        cell_arr = np.zeros((arr.shape[0]-1,arr.shape[1]-1,arr.shape[2]-1,arr.shape[3]))
        for n in range(cell_arr.shape[0]):
            for m in range(cell_arr.shape[1]):
                for k in range(cell_arr.shape[2]):
                    for t in range(cell_arr.shape[3]):
                        cell_arr[n,m,k,t] = (arr[n,m,k,t] + arr[n+1,m,k,t] + arr[n,m+1,k,t] + arr[n+1,m+1,k,t] + arr[n,m,k+1,t] + arr[n+1,m,k+1,t] + arr[n,m+1,k+1,t] + arr[n+1,m+1,k+1,t])/8.0
    
    if arr.ndim == 3:
        cell_arr = np.zeros((arr.shape[0]-1,arr.shape[1]-1,arr.shape[2]-1))
        for n in range(cell_arr.shape[0]):
            for m in range(cell_arr.shape[1]):
                for k in range(cell_arr.shape[2]):
                    cell_arr[n,m,k] = (arr[n,m,k] + arr[n+1,m,k] + arr[n,m+1,k] + arr[n+1,m+1,k] + arr[n,m,k+1] + arr[n+1,m,k+1] + arr[n,m+1,k+1] + arr[n+1,m+1,k+1])/8.0
     
    return cell_arr                       
                                       
def periodicity(arr):
    "Function that takes an array and fixes the endpoints by enforcing periodicity in the theta and phi dimensions"    
    if arr.ndim == 3:
        #fixing periodicity in theta
        arr_up = np.zeros((arr.shape[0],arr.shape[1]+1,arr.shape[2]), dtype = arr.dtype)
        arr_up[0:arr.shape[0],0:arr.shape[1],0:arr.shape[2]] = arr
        arr_up[:,-1,:] = arr[:,0,:]
        #fixing periodicity in phi
        arr_new = np.zeros((arr.shape[0],arr.shape[1]+1,arr.shape[2]+1), dtype = arr.dtype)
        arr_new[0:arr_up.shape[0],0:arr_up.shape[1],0:arr_up.shape[2]] = arr_up
        arr_new[:,:,-1] = arr_up[:,:,-1]
    
    if arr.ndim == 4:
        #fixing periodicity in theta
        arr_up = np.zeros((arr.shape[0],arr.shape[1]+1,arr.shape[2],arr.shape[3]), dtype = arr.dtype)
        arr_up[0:arr.shape[0],0:arr.shape[1],0:arr.shape[2],0:arr.shape[3]] = arr
        arr_up[:,-1,:,:] = arr[:,0,:,:]
        #fixing periodicity in phi
        arr_new = np.zeros((arr.shape[0],arr.shape[1]+1,arr.shape[2]+1,arr.shape[3]), dtype = arr.dtype)
        arr_new[0:arr_up.shape[0],0:arr_up.shape[1],0:arr_up.shape[2],0:arr_up.shape[3]] = arr_up
        arr_new[:,:,-1,:] = arr_up[:,:,-1,:]
    
    return arr_new

@ignore_warnings
def fs_locator(psi,psi_value,t):
    "Function that returns r-theta pairs of a flux surface, given Psi and phi. Returns numerical grid locations which 
     are needed for the derivative functions."
     
    fs = np.asarray(measure.find_contours(psi[:,:,0,t], psi_value))
    raxis_fs=[]
    uaxis_fs=[]
    
    if psi_value<0 or psi_value>1:
        raxis_fs = np.nan
        uaxis_fs = np.nan
    else:
        num_of_curves = fs.shape[0]
        num_of_points_of_curves = []
        
        for i in range(num_of_curves):
            num_of_points_of_curves.append(fs[i].shape[0])
        
        ind_of_good_fs = num_of_points_of_curves.index(max(num_of_points_of_curves))
        
        try:
            for i in range(fs[ind_of_good_fs].shape[0]):
                raxis_fs.append(fs[ind_of_good_fs][i][0])
                uaxis_fs.append(fs[ind_of_good_fs][i][1])
        except IndexError:
            raxis_fs = np.nan
            uaxis_fs = np.nan
    
    return raxis_fs,uaxis_fs
    
def flux_surf_ang(l,state_vec):
    "Differential equations for arc lenth of a flux surface and flux angle. Inputs are given in terms of logical coordinates."    
    r = state_vec[0]
    u = state_vec[1]
    uf = state_vec[2]
    phi = state_vec[3]
   
    drdl = -dpsidtheta_int((r,u,phi))/(np.sqrt(np.power(dpsidtheta_int((r,u,phi)),2) \
                                            + np.power(r,2)*np.power(dpsidr_int((r,u,phi)),2)))
    dthetadl = dpsidr_int((r,u,phi))/(np.sqrt(np.power(dpsidtheta_int((r,u,phi)),2) \
                                            + np.power(r,2)*np.power(dpsidr_int((r,u,phi)),2)))
    
    dufdl = B3_int((r,u,phi,0))/(np.sqrt(np.power(dpsidtheta_int((r,u,phi)),2) \
                                            + np.power(r,2)*np.power(dpsidr_int((r,u,phi)),2)))
    dphidl = 0
    
    return[drdl,dthetadl,dufdl,dphidl]

def integrator(rs,us,phi):
    '''Integrating to find r(l), theta(l) around a flux surface. Returning r's, theta's, theta_f's and l's.'''    
    rs = CnvNumber2LogicalR(rs)
    solver = scipy.integrate.ode(flux_surf_ang).set_integrator('dopri5')
    
    uf0 = 0.0
    dl = 0.01
            
    state_vec0,l0 = [rs[0],0.0,uf0,phi],0.0
    solver.set_initial_value(state_vec0,l0)
    y,t = [], []
    warnings.filterwarnings("ignore")
    
    while solver.successful() and solver.y[1]<2*np.pi:
        solver.set_initial_value([solver.y[0],solver.y[1],solver.y[2],solver.y[3]],solver.t)
        solver.integrate(solver.t+dl)
        y.append(solver.y)
        t.append(solver.t)
    warnings.resetwarnings()
    y = np.array(y)
    t = np.array(t)
    
    return y[:,0],y[:,1],y[:,2],t

def dict_interp(dictionary, psi, phi, r, u, uf, l):
    "Function that takes lists of values of r, u, uf and finds interpolating functions. Then, it adds those in a psi, phi
    dictionary."    
    try:
        r_of_l = interp1d(l,r,kind='quadratic',fill_value='extrapolate')
    except ValueError:
        r_of_l = 0
    try:
        u_of_l = interp1d(l,u,kind='quadratic',fill_value='extrapolate')
    except ValueError:
        u_of_l = 0
    try:
        l_of_u = interp1d(u,l,kind='quadratic',fill_value='extrapolate')
    except ValueError:
        l_of_u = 0
    try:
        uf_of_l = interp1d(l,uf,kind='quadratic',fill_value='extrapolate')
    except ValueError:
        uf_of_l = 0
        
    dictionary[(str(psi), str(phi))] = {'r_l': r_of_l, 'u_l': u_of_l, 'l_u': l_of_u, 'uf_l': uf_of_l}
    
    return dictionary

def q_prof():
    "Function that reconstructs the initial q-profile. It also is a dictionary. Because the psi array is axisymmetric w.r.t the
    phi axis, it is only needed to create the dictionary for phi=0."    
    q = []
    psi_dict = {}
    psi_list = np.linspace(0.01,1.0,100)
    for psin in psi_list:
        psi = norm_Psi2Psi_Init(psin) 
        r_fs, u_fs = fs_locator((psi),0) 
        r, u, uf, l = so.integrator(np.asarray(r_fs),np.asarray(u_fs),0)
        q_redef = uf[-1]/(2*np.pi)
        uf_rsc = uf/q_redef
        # Dictionary creation
        psi_dict = dict_interp(psi_dict, psi, 0, r, u, uf_rsc, l)
        f = open('dict_bonfiglio.pkl','wb')
        pkl.dump(psi_dict,f)
        f.close()
        #if phi == 0:
        q.append(q_redef)
        print('psi = ', psin)
    return q, psi_list, psi_dict
def q_prof_no_dict():
    "Function that reconstructs the initial q-profile. It also creates a dictionary. Because the psi array is axisymmetric w.r.t
    the phi axis, it is only needed to create the dictionary for phi=0."    
    q = []
    psi_list = np.linspace(0.01,1.0,100)
    for psin in psi_list:
        psi = norm_Psi2Psi_Init(psin)
        r_fs, u_fs = fs_locator((psi),0) 
        r, u, uf, l = integrator(np.asarray(r_fs),np.asarray(u_fs),0)
        q_redef = uf[-1]/(2*np.pi)
        uf_rsc = uf/q_redef
        q.append(q_redef)
        print('psi = ', psin)
    return q, psi_list

def psi_proj(t,psi_min,norm):
    "Function that projects the psi array in the new grid."    
    psi_list = []
    uf_list = []
    psi_int_list = []
    for key in list(loaded.keys())[:]:
        l_min = loaded[key[0],str(0)]['l_u'](0.001)
        l_max = loaded[key[0],str(0)]['l_u'](6.28)
        ls = np.linspace(l_min,l_max,30)
        for l in ls:
            r = loaded[key[0],str(0)]['r_l'](l)
            u = loaded[key[0],str(0)]['u_l'](l)
            uf = loaded[key[0],str(0)]['uf_l'](l)
            p = (Psi_int((r,u,0.0,t))-psi_min[t])/(norm[t]-psi_min[t])
            psi_list.append((float(key[0])-psi_min[t])/(norm[t]-psi_min[t]))
            uf_list.append(uf)
            psi_int_list.append(p)
    triangObj = Triangulation(psi_list,uf_list)
    try:
        psi_int_list = np.asarray(psi_int_list)
    except:
        psi_int_list = psi_int_list[:,0] # Sometimes Triangulation gives an array. Uncomment
    tcp = LinearTriInterpolator(triangObj, psi_int_list)
    tcp_int = tcp(PS,TF)
    
    plt.contour(PS,TF,tcp_int,20)
    plt.colorbar()
    plt.xlabel(r'$\Psi$')
    plt.ylabel(r'$\theta_f$')
    plt.show()
    
def new_coords(psi_dict, arr_int, time):
    "Function that takes the dictionary that relates psi values and interpolating functions and returns an array on the new psi-
    theta_f mesh."    
    psi_list = []
    uf_list = []
    arr_int_list = []
    for key in list(psi_dict.keys())[:]:
        l_min = psi_dict[key[0],str(0)]['l_u'](0.001)
        l_max = psi_dict[key[0],str(0)]['l_u'](6.28)
        ls = np.linspace(l_min,l_max,30)
        for l in ls:
            r = psi_dict[key[0],str(0)]['r_l'](l)
            u = psi_dict[key[0],str(0)]['u_l'](l)
            uf = psi_dict[key[0],str(0)]['uf_l'](l)
            value = arr_int((r,u,0.0,time))
            psi_list.append((float(key[0])-psi_min)/(norm-psi_min))
            uf_list.append(uf)
            arr_int_list.append(value)
    
    triangObj = Triangulation(psi_list,uf_list)
    tcp = LinearTriInterpolator(triangObj, arr_int_list)
    
    return tcp

def flux_angle():
    "Returns an array where theta_f is evaluated at every r, theta point"
    r_list = []
    u_list = []
    uf_list = []
    for key in list(loaded.keys())[:]:
        l_min = loaded[key[0],str(0)]['l_u'](0.001)
        l_max = loaded[key[0],str(0)]['l_u'](6.28)
        ls = np.linspace(l_min,l_max,30)
        for l in ls:
            r = loaded[key[0],str(0)]['r_l'](l)
            u = loaded[key[0],str(0)]['u_l'](l)
            uf = loaded[key[0],str(0)]['uf_l'](l)
            r_list.append(r)
            u_list.append(u)
            uf_list.append(uf)
            
    triangObj = Triangulation(r_list,u_list)
    tcp = LinearTriInterpolator(triangObj, uf_list)
    fa = np.asarray(tcp(RRI,TTI))
    return fa

def dB():
    global dB1, dB2, dB3
    B10 = B1[:,:,:,0]
    dB1 = B1-B10[:,:,:,np.newaxis]
    B20 = B2[:,:,:,0]
    dB2 = B2-B20[:,:,:,np.newaxis]
    B30 = B3[:,:,:,0]
    dB3 = B3-B30[:,:,:,np.newaxis]
        
def Dictionary_Creation():
    "Creating the dictionary"
    global init_q_list, init_psi_list
    init_q_list, init_psi_list, psi_dict = q_prof()
    f = open('dict_11_visc.pkl','wb')
    pkl.dump(psi_dict,f)
    f.close()
    print("Dictionary has been created and saved.")

def Load_Dictionary():
    "Loading the dictionary"
    global loaded
    f = open("dict_bonfiglio.pkl", "rb")
    loaded = pkl.load(f)
    f.close()
    #loaded = pkl.load(open("psi_phi_dict.pkl", "rb"))
    #loaded = pkl.load(open("dict_11_visc.pkl", "rb"))
    print("Dictionary loaded.")
    
def test_psi_projection(t,psi_min,norm):
    "Test that should plot straight flux surfaces"
    global PS,TF
    # Creating the meshgrid
    #Ps = np.linspace(np.amin(P_N[:,:,0,0]),np.amax(P_N[:,:,0,0]),100)
    Ps = np.linspace(0,1,100)
    Tf = np.linspace(0.0,2.0*np.pi,100)
    (PS,TF) = np.meshgrid(Ps,Tf)
    # Projecting psi
    psi_proj(t,psi_min,norm)
    print("If field lines are straight, test passed.")
        
def test_locator(psi,psin,t,psi_min,norm):
    "Test of the flux surface locator"
    psi_min_t = psi_min[t]
    norm_t = norm[t]
    psi_actual = norm_Psi2Psi_Init(psin,psi_min_t,norm_t)
    rs, us = fs_locator(psi,psi_actual,t) # fs_locator works on original (flipped) Psi array
    r_log = CnvNumber2LogicalR(rs) # fs_locator gives values on numerical grid
    u_log = CnvNumber2LogicalU(us)
    #x_X,y_X = findXpoint(t)
    plt.plot(X_int((r_log,u_log,0)),Z_int((r_log,u_log,0)))
    #plt.plot(x_X,y_X,"x")
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.show()

def test_integrator(psin):
    psi_actual = norm_Psi2Psi_Init(psin)
    rs, us = fs_locator(psi_actual,0) # fs_locator works on original (flipped) Psi array
    r_log = CnvNumber2LogicalR(rs) # fs_locator gives values on numerical grid
    u_log = CnvNumber2LogicalU(us) # integrator takes values on logical grid
    r_Integr, u_Integr, uf, _ = so.integrator(np.asarray(rs),np.asarray(us),0)
    q_redef = uf[-1]/(2*np.pi)
    print("Q:",q_redef)
    plt.plot(X_int((r_log,u_log,0)),Z_int((r_log,u_log,0)),'x')
    plt.plot(X_int((r_Integr,u_Integr,0)),Z_int((r_Integr,u_Integr,0)),'*')
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.show()
    
    
def Grid_Interpolations(psit,B1,B2,B3,J3=False):
    "Find Interpolators on logical coordinate grid nodes."    
    global X_int, Z_int, B1tor_int, B2tor_int, B3_int, B1_int, B2_int, B1tor, B2tor, psit_int, J3_int
    
    B1tor = np.mean(B1,axis=2)
    B2tor = np.mean(B2,axis=2)
    
    if t_dim == 0:
        # Interpolators
        psit_int = RegularGridInterpolator((r,theta), psit[:,:], method='linear', bounds_error=False, fill_value = 0)
        X_int = RegularGridInterpolator((r,theta,phi), X[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        Z_int = RegularGridInterpolator((r,theta,phi), Z[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        B1_int = RegularGridInterpolator((r,theta,phi), B1[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        B2_int = RegularGridInterpolator((r,theta,phi), B2[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        B3_int = RegularGridInterpolator((r,theta,phi), B3[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        B1tor_int = RegularGridInterpolator((r,theta), B1tor[:,:], method='linear', bounds_error=False, fill_value = 0)
        B2tor_int = RegularGridInterpolator((r,theta), B2tor[:,:], method='linear', bounds_error=False, fill_value = 0)
        try:
            J3_int = RegularGridInterpolator((r,theta,phi), J3[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        except:
            J3_int=0
    else:
        psit_int = RegularGridInterpolator((r,theta,t), psit[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        X_int = RegularGridInterpolator((r,theta,phi), X[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        Z_int = RegularGridInterpolator((r,theta,phi), Z[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        B1_int = RegularGridInterpolator((r,theta,phi,t), B1[:,:,:,:], method='linear', bounds_error=False, fill_value = 0)
        B2_int = RegularGridInterpolator((r,theta,phi,t), B2[:,:,:,:], method='linear', bounds_error=False, fill_value = 0)
        B3_int = RegularGridInterpolator((r,theta,phi,t), B3[:,:,:,:], method='linear', bounds_error=False, fill_value = 0)
        B1tor_int = RegularGridInterpolator((r,theta,t), B1tor[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        B2tor_int = RegularGridInterpolator((r,theta,t), B2tor[:,:,:], method='linear', bounds_error=False, fill_value = 0)
        try:
            J3_int = RegularGridInterpolator((r,theta,phi,t), J3[:,:,:,:], method='linear', bounds_error=False, fill_value = 0)
        except:
            J3_int=0
    
    print("Logical Grid Interpolations done. All array (A) interpolators take A_int.")
    
def Grid_Cell_Interpolations(B1,B2,B3,psi,J3=False):
    "Find Interpolators on logical coordinate grid cells."    
    global B1c_int, B2c_int, B3c_int, J3c_int, psic_int
    
    if t_dim == 0:
        # Interpolators
        B1c_int = RegularGridInterpolator((rc,uc,phic), B1[:,:,:], method='linear', bounds_error=False, fill_value = None)
        B2c_int = RegularGridInterpolator((rc,uc,phic), B2[:,:,:], method='linear', bounds_error=False, fill_value = None)
        B3c_int = RegularGridInterpolator((rc,uc,phic), B3[:,:,:], method='linear', bounds_error=False, fill_value = None)
        psic_int = RegularGridInterpolator((rc,uc,phic), psi[:,:,:], method='linear', bounds_error=False, fill_value = None)
        try:
            J3c_int = RegularGridInterpolator((rc,uc,phic), J3[:,:,:], method='linear', bounds_error=False, fill_value = None)
        except:
            J3c_int = 0
        
    else:
        B1c_int = RegularGridInterpolator((rc,uc,phic,t), B1[:,:,:,:], method='linear', bounds_error=False, fill_value = None)
        B2c_int = RegularGridInterpolator((rc,uc,phic,t), B2[:,:,:,:], method='linear', bounds_error=False, fill_value = None)
        B3c_int = RegularGridInterpolator((rc,uc,phic,t), B3[:,:,:,:], method='linear', bounds_error=False, fill_value = None)
        psic_int = RegularGridInterpolator((rc,uc,phic,t), psi[:,:,:,:], method='linear', bounds_error=False, fill_value = None)
        try:
            J3c_int = RegularGridInterpolator((rc,uc,phic,t), J3[:,:,:,:], method='linear', bounds_error=False, fill_value = None)
        except:
            J3c_int = 0

def Grid_Cell_Interpolations_Covariant(B_1,B_2,B_3):
    "Find Interpolators on logical coordinate grid cells."    
    global B_1c_int, B_2c_int, B_3c_int
    
    if t_dim == 0:
        # Interpolators
        B_1c_int = RegularGridInterpolator((rc,uc,phic), B_1[:,:,:], method='linear', bounds_error=False, fill_value = None)
        B_2c_int = RegularGridInterpolator((rc,uc,phic), B_2[:,:,:], method='linear', bounds_error=False, fill_value = None)
        B_3c_int = RegularGridInterpolator((rc,uc,phic), B_3[:,:,:], method='linear', bounds_error=False, fill_value = None)
    else:
        B_1c_int = RegularGridInterpolator((rc,uc,phic,t), B_1[:,:,:,:], method='linear', bounds_error=False, fill_value = None)
        B_2c_int = RegularGridInterpolator((rc,uc,phic,t), B_2[:,:,:,:], method='linear', bounds_error=False, fill_value = None)
        B_3c_int = RegularGridInterpolator((rc,uc,phic,t), B_3[:,:,:,:], method='linear', bounds_error=False, fill_value = None)

def C2N_Evaluations():
    "Takes the cell interpolators and evaluates them in the node grid."    
    if t_dim == 0:
        B1 = B1c_int((RI,TI,PI))
        B2 = B2c_int((RI,TI,PI))
        B3 = B3c_int((RI,TI,PI))
        if J3c_int !=0:
            J3 = J3c_int((RI,TI,PI))
        else:
            J3=0
    else:
        B1 = B1c_int((RRRI,TTTI,PPPI,TIM))
        B2 = B2c_int((RRRI,TTTI,PPPI,TIM))
        B3 = B3c_int((RRRI,TTTI,PPPI,TIM))
        psi = psic_int((RRRI,TTTI,PPPI,TIM))
        if J3c_int !=0:
            J3 = J3c_int((RRRI,TTTI,PPPI,TIM))
        else:
            J3=0
    
    return B1, B2, B3, J3, psi

def C2N_Evaluations_Covariant():
    "Takes the cell interpolators and evaluates them in the node grid."    
    if t_dim == 0:
        B_1 = B_1c_int((RI,TI,PI))
        B_2 = B_2c_int((RI,TI,PI))
        B_3 = B_3c_int((RI,TI,PI))
    else:
        B_1 = B_1c_int((RRRI,TTTI,PPPI,TIM))
        B_2 = B_2c_int((RRRI,TTTI,PPPI,TIM))
        B_3 = B_3c_int((RRRI,TTTI,PPPI,TIM))
    
    return B_1, B_2, B_3
    
def findLevelsZero(B1t,B2t,t,lower_level=-0.2):
    "Function that returns the x,y points of the zero level of Br, Btheta fields. Their intersection is the X-point.
    Care should be taken to change the minimum level so that you avoid the points at the axis.
    "      
    if t_dim == 0:
        cs1 = plt.contour(X[:,:,0],Z[:,:,0],B1t[:,:],levels=[0])
        cs2 = plt.contour(X[:,:,0],Z[:,:,0],B2t[:,:],levels=[0])
    else:
        cs1 = plt.contour(X[:,:,0],Z[:,:,0],B1t[:,:,t],levels=[0])
        cs2 = plt.contour(X[:,:,0],Z[:,:,0],B2t[:,:,t],levels=[0])
    plt.close()
    
    for i in range(len(cs1.collections[0].get_paths())):
        p1 = cs1.collections[0].get_paths()[i]
        v1 = p1.vertices
        x1_int = v1[:,0]
        y1_int = v1[:,1]
        if min(y1_int) < lower_level: # Level to avoid points on axis and find lower X-point
            x1 = x1_int
            y1 = y1_int
        else:
            pass
        
    for i in range(len(cs2.collections[0].get_paths())):
        p2 = cs2.collections[0].get_paths()[i]
        v2 = p2.vertices
        x2_int = v2[:,0]
        y2_int = v2[:,1]
        if min(y2_int) < lower_level: # Level to avoid points on axis and find lower X-point
            x2 = x2_int
            y2 = y2_int
        else:
            pass
    plt.plot(x1,y1)
    plt.plot(x2,y2)
    plt.show()
    return x1, y1, x2, y2

def sort_line(x,y,lower_level):
    "Sometimes the zero level of B1 or B2 might be not be a function. By visual inspection, we can tweak the zero level curve
    to contain only points that make it a function.
    "
    new_x = []
    new_y = []
    for i in range(len(y)):
        if y[i] < lower_level:
            new_x.append(x[i])
            new_y.append(y[i])
        else:
            pass
    return new_x, new_y

def findXpoint(B1t,B2t,t):
    "Find the intersection of the zeros of Br, Btheta which is the X-point. Care should be taken to set the right x0 for the
    numerical method so it is not outside the interpolation range.
    "
    # Level to avoid points on axis and find lower X-point
    # -1.2 for double tearing, -0.2 for sawtooth
    lower_level = -1.0
    x1,y1,x2,y2 = findLevelsZero(B1t,B2t,t,lower_level)
    x1,y1 = sort_line(x1,y1,lower_level)
    x2,y2 = sort_line(x2,y2,lower_level)

    f1 = interp1d(x1, y1, kind = 'cubic')
    f2 = interp1d(x2, y2, kind = 'cubic')

    def difference(x):
        return np.abs(f1(x) - f2(x))
    
    # Sometimes visual inspection is needed to select sp
    #sp = findOverlap(x1,x2) # works well for sawtooth case
    sp = x2[math.ceil(len(x2)/2)] # works well for double tearing case
    
    xcross = None
    i = 0
    while xcross is None:
        try:
            xcross = scipy.optimize.fsolve(difference, x0=x1[i]) # Set the x0 appropriately to avoid being outside interpolation range
        except ValueError:
            i += 1
        
    #    plt.title("t=%s"%(t))
    #    plt.plot(x1,y1,label='B1 zero')
    #    plt.plot(x2,y2,label='B2 zero')
    #    plt.plot(xcross,f1(xcross),'kx',markersize=30)
    #    plt.legend()
    #    plt.show()

    return xcross, f1(xcross)

def findOverlap(x1,x2):
    "Find the overlap region between two ranges."
    x1L = x1[0]
    x1R = x1[-1]
    if x1L > x1R:
        x1R = x1[0]
        x1L = x1[-1]
    
    x2L = x2[0]
    x2R = x2[-1]
    if x2L > x2R:
        x2R = x2[0]
        x2L = x2[-1]
        
    if x1L < x2L:
        inner_left = x2L
    else:
        inner_left = x1L
        
    if x1R > x2R:
        inner_right = x2R
    else:
        inner_right = x1R
        
    starting_point = (inner_right-inner_left)/2. + inner_left
    
    return starting_point
    
    
def magneticAxisAndXpoint(psit,psit_int,B1t,B2t,t):
    "Calculates values of flux at magnetic axis and X-point. This needs at least the first norm value to be calculated
    correctly. If it fails on the first step, maybe visual inspection is needed to initialize it properly.
    "
    #print("t:",t)
    if t_dim == 0:
        psi_min = np.amin(psit[:,:]) # find min of flux
        x_X,y_X = findXpoint(B1t,B2t,0) # location of X-point
        #x_X_log_grid,y_X_log_grid = pntCnvInGrid_simple_toroidal(x_X[0],y_X[0]) # Convert X-point coords in logical grid
        x_X_log_grid,y_X_log_grid = pntCnvInGrid_shaped(x_X[0],y_X[0]) # Convert X-point coords in logical grid
        norm = psit_int((x_X_log_grid,y_X_log_grid)) # psi value at X-point (normalization constant)
    else:
        psi_min = np.amin(psit[:,:,t]) # find min of flux 
        #print("psi_min=",psi_min)
        x_X,y_X = findXpoint(B1t,B2t,t) # location of X-point
        #print("x_X,y_X=",x_X,y_X)
        x_X_log_grid,y_X_log_grid = pntCnvInGrid_shaped(x_X[0],y_X[0])
        #x_X_log_grid,y_X_log_grid = pntCnvInGrid_simple_toroidal(x_X[0],y_X[0]) # Convert X-point coords in logical grid
        #print("x_X_log_grid,y_X_log_grid=",x_X_log_grid,y_X_log_grid)
        norm = psit_int((x_X_log_grid,y_X_log_grid,t)) # psi value at X-point (normalization constant)
        #print("norm=",norm)
    return psi_min,norm # norm comes out as 1-dimensional array
    
    
"""
