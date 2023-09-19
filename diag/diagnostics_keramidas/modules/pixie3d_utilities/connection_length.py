import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import os
import h5py
import pixie_read_st as pxr

"""Contains functions meant to process nemato connection length files. 
   Main use is to read the nemato output and create an array that has the same dimensions as the initial condition grid.
   The output initially will not have the same dimensions because some initial conditions fail to finish. This is corrected by 
   the cnc_coord_tuples() and conn_array(). It also contains a function that bins connection lengths.
   
   CAUTION: Some functions require the filepath and file name convention (e.g., read_connection_lengths()."""

def store_coordinate_arrays(filepath):
    f = h5py.File(filepath+'pixie3d.h5','r')
    timesteps = list(f.keys())
    variables = list(f[timesteps[0]].keys())
    cell_var = list(f[timesteps[0]][variables[5]].keys())
    node_var = list(f[timesteps[0]][variables[6]].keys())
    diag_var = list(f[timesteps[0]][variables[3]].keys())
    X = np.asarray(f[timesteps[0]][variables[6]][node_var[0]])
    Y = np.asarray(f[timesteps[0]][variables[6]][node_var[1]])
    Z = np.asarray(f[timesteps[0]][variables[6]][node_var[2]])
    Psi = np.asarray(f[timesteps[0]][variables[3]][diag_var[4]])
    X = np.swapaxes(X,0,2)
    Z = np.swapaxes(Z,0,2)
    Psi = np.swapaxes(Psi,0,2)
    if not os.path.exists(filepath + 'cnc/'):
        os.makedirs(filepath + 'cnc/')
        np.save(filepath + 'cnc/'+'X.npy',X)
        np.save(filepath + 'cnc/'+'Z.npy',Z)
        np.save(filepath + 'cnc/'+'Psi.npy',Psi)
    
def read_coordinates(filepath):
    X = np.load(filepath + 'cnc/' + 'X.npy')
    Z = np.load(filepath + 'cnc/' + 'Z.npy')
    Psi = np.load(filepath + 'cnc/' + 'Psi.npy')
    return X, Z, Psi
    
def interpolate(X,Z,Psi):
    #Grid of interpolation
    r = np.linspace(0.0,1.0,num=Psi.shape[0])
    theta = np.linspace(0.0,2.*np.pi,num=Psi.shape[1])
    phi = np.linspace(0.0,2.*np.pi,num=Psi.shape[2])

    #Interpolation
    Psi_int = RegularGridInterpolator((r,theta,phi), Psi[:,:,:], method='linear', bounds_error=False, fill_value = 0)
    X_int = RegularGridInterpolator((r,theta,phi), X[:,:,:], method='linear', bounds_error=False, fill_value = 0)
    Z_int = RegularGridInterpolator((r,theta,phi), Z[:,:,:], method='linear', bounds_error=False, fill_value = 0)
    return X_int, Z_int, Psi_int

def read_connection_lengths(filepath,timestep,filenumber):
    xs = []
    zs = []
    conn = []
    
    for i in range(filenumber):
        file = open(filepath+"cnc-poinc_t=%s.txt_proc%s" %(str(timestep),str(i)),'r')
        
        for columns in (raw.strip().split() for raw in file):
            xs.append(float(columns[0]))
            zs.append(float(columns[1]))
            conn.append(float(columns[2]))
    return xs,zs,conn

def read_connection_lengths_logical(filepath,timestep,filenumber):
    """Reads cnc-files with dumped logical points"""
    rs = []
    us = []
    conn = []
    
    for i in range(filenumber):
        file = open(filepath+"cnc-poinc_t=%s.txt_proc%s" %(str(timestep),str(i)),'r')
        
        for columns in (raw.strip().split() for raw in file):
            rs.append(float(columns[0]))
            us.append(float(columns[1]))
            conn.append(float(columns[3]))
    return rs,us,conn

def cnc_coord_tuples(rs,us,conn):
    """Creates tuple of coordinates and coordinates+connection lenghts."""
    CNL = list(zip(rs,us,conn))
    Ctuples = [(item[0],item[1]) for item in CNL]
    
    return CNL,Ctuples

def conn_array(CNL, Ctuples, R_mg, U_mg, Nr, Nu):
    """Creates connection length array with the same dimensions as the initial grid"""
    conn_array = np.zeros((Nr,Nu))
    for i in range(Nr):
        for j in range(Nu):
            if (R_mg[i,j],U_mg[i,j]) in Ctuples:
                ind_of_tup = Ctuples.index((R_mg[i,j],U_mg[i,j]))
                conn_array[i,j] = CNL[ind_of_tup][2]
            else:
                conn_array[i,j] = np.nan
    return conn_array

def R_Bins2(rs, conn, conn_norm, conn_max, n_bins):
    """Splits the connection lengths in bins of minor radius and returns averages, mins and maxs."""
    R_bins = list(set(rs))
    R_bins.sort()
    chunk = int(math.floor(len(R_bins)/(n_bins+1)))
    
    bin_list = np.zeros(n_bins)
    
    num_list = np.zeros(n_bins)
    
    for i in range(len(rs)):
        if conn[i] == conn_max: # exclude initial conditions that do not terminate
            pass
        else:                            
            if rs[i]>R_bins[0] and rs[i]<R_bins[chunk]: # edge case, start of array 
                bin_list[0] += conn[i]
                num_list[0] += 1
            if rs[i]>R_bins[(n_bins-1)*chunk+1] and rs[i]<R_bins[-1]: # edge case, end of array
                bin_list[n_bins-1] += conn[i]
                num_list[n_bins-1] += 1
            for j in range(1,n_bins-1): # general case
                if rs[i]>R_bins[j*chunk +1] and rs[i]<R_bins[(j+1)*chunk]:
                    bin_list[j] += conn[i]
                    num_list[j] += 1
    
    num_list_wo_zero = [x if x>0 else x+1 for x in num_list] # correct for case where list is empty            
    
    avg_list_norm = [x/y for x,y in zip(bin_list,num_list_wo_zero)] # average
    avg_list = [x*(conn_norm/1000) for x in avg_list_norm] # average denormed in Kilometers
    
    r_list = [R_bins[j*chunk] for j in range(1,n_bins)]
    r_list = r_list + [R_bins[-1]] # edge case
    
    return avg_list, r_list


    
    
    
    