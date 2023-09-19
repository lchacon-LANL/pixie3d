import numpy as np
import scipy.interpolate as si
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
#from scipy.interpolate import RegularGridInterpolator # outdated module
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import math
import pylab as py
import master_read as m

def shaping_plot():
    Ro = 6.219577546
    Zo = 0.5143555944
    alpha = 2.24
    kappa = 1.9
    delta = 0.6
    zeta = 0.06

    theta = np.linspace(0,2*np.pi,100)
    Rs = Ro + alpha*np.cos(theta + np.arcsin(delta*np.sin(theta)))
    Zs = Zo + kappa*alpha*np.sin(theta + zeta*np.sin(2*theta))

    return Rs, Zs

def gridsNunits():    
    global psi_norm, psig, psig_gap, Rg, RR, ZZ, a, Rs, Zs, psi_grid_dim
    
    #R-Z grid
    RR, ZZ = np.meshgrid(m.rr,m.zz)
    Rs, Zs = shaping_plot()
    
    #psi grid
    psi_grid_dim = len(m.DS.qpsi)
    psig = np.linspace(m.DS.simag,m.DS.sibry,psi_grid_dim)
    psig_gap = psig[1]-psig[0]
    psi_norm = (psig-m.DS.simag)/(m.DS.sibry - m.DS.simag)

    #radial points on the flux grid (extends from magnetic axis to right boundary)
    #Normalization in Pixie3d units with a=(0.5)*(right_limiter-left_limiter)
    try:
        a = 0.5*(max(m.DS.rlim)-min(m.DS.rlim))
        Rg = np.linspace(m.DS.rmaxis/a, m.r_sep_right_intersection/a, psi_grid_dim)
        Rg_gap = Rg[1]-Rg[0]
    except ValueError:
        print("WARNING: No limiter. Units not normalized to Pixie units.")
        psi_grid_dim = len(m.DS.qpsi)
        Rg = np.linspace(m.DS.rmaxis, max(m.DS.rbbbs), psi_grid_dim)
        Rg_gap = Rg[1]-Rg[0]

def Profile_plots():
    '''Pressure Profile Plots'''

    plt.title(r'$\Psi$ contours')
    plt.contour(RR,ZZ,m.DS.psirz,80)
    plt.plot(m.DS.rmaxis,m.DS.zmaxis,'or')
    plt.plot(m.DS.rlim,m.DS.zlim,'k')
    plt.plot(m.DS.rbbbs,m.DS.zbbbs,'r')
    plt.plot(Rs,Zs,'b')
    plt.show()

    plt.title(r'Pressure vs $\Psi_N$')
    plt.plot(psi_norm,m.DS.pres)
    plt.show()
    
    plt.title(r'pprime vs $\Psi$')
    plt.plot(psig,m.DS.pprime)
    plt.show()

    plt.title(r'Ffprime vs $\Psi$')
    plt.plot(psig,m.DS.ffprim)
    plt.show()

    plt.title(r'pprime vs $\Psi_N$')
    plt.plot(psi_norm,m.DS.pprime)
    plt.show()

    plt.title(r'Ffprime vs $\Psi_N$')
    plt.plot(psi_norm,m.DS.ffprim)
    plt.show()


    dpress = (-np.roll(m.DS.pres,1)+np.roll(m.DS.pres,-1))/psig_gap
    plt.title(r'fin. diff. pprime vs $\Psi$')
    plt.plot(psig[1:-1], dpress[1:-1])
    plt.show()

    if m.DS.limitr > 0:
        plt.title(r'$\Psi$ along the midline')
        plt.plot(m.rr,m.DS.psirz[m.grid_mag_center_z_index,:])
        plt.plot((m.rr[m.r_sep_index_left_intersection],m.rr[m.r_sep_index_left_intersection]),(min(m.DS.psirz[m.grid_mag_center_z_index,:]),max(m.DS.psirz[m.grid_mag_center_z_index,:])),'r-')
        plt.plot((m.rr[m.r_sep_index_right_intersection],m.rr[m.r_sep_index_right_intersection]),(min(m.DS.psirz[m.grid_mag_center_z_index,:]),max(m.DS.psirz[m.grid_mag_center_z_index,:])),'r-')
        plt.plot((m.rr[m.r_lim_index_left_intersection],m.rr[m.r_lim_index_left_intersection]),(min(m.DS.psirz[m.grid_mag_center_z_index,:]),max(m.DS.psirz[m.grid_mag_center_z_index,:])),'k-')
        plt.plot((m.rr[m.r_lim_index_right_intersection],m.rr[m.r_lim_index_right_intersection]),(min(m.DS.psirz[m.grid_mag_center_z_index,:]),max(m.DS.psirz[m.grid_mag_center_z_index,:])),'k-')
        plt.plot((m.rr[m.grid_mag_center_r_index],m.rr[m.grid_mag_center_r_index]),(min(m.DS.psirz[m.grid_mag_center_z_index,:]),max(m.DS.psirz[m.grid_mag_center_z_index,:])),'k-')
        plt.show()
    else:
        pass

    pres_int = interp1d(psig,m.DS.pres, kind='quadratic',bounds_error = False, fill_value=np.nan)

def int_arr(interpolant,psi_norm):
    arr = np.zeros((m.DS.psirz.shape[0],m.DS.psirz.shape[1]))
    for i in range(0,m.DS.psirz.shape[0]):
        for j in range(0,m.DS.psirz.shape[1]):
            arr[i,j]=float(interpolant(m.DS.psirz[i,j]))
    return arr

def pressure_contours():
    pres_arr = int_arr(pres_int, psi_norm)
    
    plt.title('Pressure Contours')
    plt.contour(RR,ZZ,pres_arr,80)
    plt.plot(m.DS.rmaxis,m.DS.zmaxis,'or')
    plt.plot(m.DS.rlim,m.DS.zlim,'k')
    plt.plot(m.DS.rbbbs,m.DS.zbbbs,'r')
    plt.show()

    if m.DS.limitr > 0:
        
        l_index = min(m.r_lim_index_left_intersection , m.r_lim_index_right_intersection)
        r_index = max(m.r_lim_index_left_intersection , m.r_lim_index_right_intersection)

        plt.title('Pressure profile in device on zmaxis line')
        plt.plot(m.rr[l_index:r_index], pres_arr[m.grid_mag_center_z_index,l_index:r_index,])
#plt.plot((m.rr[m.r_sep_index_left_intersection],m.rr[m.r_sep_index_left_intersection(np.nanmin(pres_arr[m.grid_mag_center_z_index,l_index:r_index]),np.nanmax(pres_arr[m.grid_mag_center_z_index,l_index:r_index])),'r-')
        plt.plot((m.rr[m.r_sep_index_right_intersection],m.rr[m.r_sep_index_right_intersection]),(np.nanmin(pres_arr[m.grid_mag_center_z_index,l_index:r_index]),np.nanmax(pres_arr[m.grid_mag_center_z_index,l_index:r_index])),'r-')
        plt.show()

    else:
        pass

    
def interp2d_pairs(*args,**kwargs):
    """ Same interface as interp2d but the returned interpolant will evaluate its inputs as pairs of values.
    """
    # Internal function, that evaluates pairs of values, output has the same shape as input
    def interpolant(x,y,f):
        x,y = np.asarray(x), np.asarray(y)
        return (si.dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], x.ravel(), y.ravel())[0]).reshape(x.shape)
    # Wrapping the scipy interp2 function to call out interpolant instead
    return lambda x,y: interpolant(x,y,si.interp2d(*args,**kwargs))

def Psi_on_sep():
    '''Scaled Contour Plots'''
    ys = np.linspace(0,1,len(m.DS.rbbbs))

    psi_int = interp2d_pairs(m.rr,m.zz,m.DS.psirz)
    psi_on_sep = psi_int(m.DS.rbbbs,m.DS.zbbbs)
    plt.title(r'$\Psi$ on separatrix')
    plt.plot(ys,psi_on_sep)
    plt.show()

    plt.title('R-Z coordinates of separatrix')
    plt.plot(ys, m.DS.rbbbs,'r')
    plt.plot(ys, m.DS.zbbbs,'b')
    plt.show()

def intersections():
    global idx32, idx43, idx2, idx54, idx65, idx76, idx87, idx4, idx3, idx83, idx52,\
           idx73, idx72, idx5, idx1, idx74, idx85, idx53
    
    if m.DS.qpsi[-1]<=0:
        qpsi = [-x for x in m.DS.qpsi]
    else:
        qpsi = m.DS.qpsi     
 
    idx1 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(1./1.)))).flatten()
    idx32 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(3./2.)))).flatten()
    idx43 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(4./3.)))).flatten()
    idx2 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(2./1.)))).flatten()
    idx54 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(5./4.)))).flatten()
    idx65 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(6./5.)))).flatten()
    idx76 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(7./6.)))).flatten()
    idx87 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(8./7.)))).flatten()
    idx4 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(4./1.)))).flatten()
    idx3 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(3./1.)))).flatten()
    idx83 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(8./3.)))).flatten()
    idx52 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(5./2.)))).flatten()
    idx73 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(7./3.)))).flatten()
    idx72 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(7./2.)))).flatten()
    idx5 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(5./1.)))).flatten()
    idx74 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(7./4.)))).flatten()
    idx85 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(8./5.)))).flatten()
    idx53 = np.argwhere(np.diff(np.sign(np.asarray(qpsi)-(5./3.)))).flatten()
    
    print('1/1:',Rg[idx1])
    #print('3/2:',Rg[idx32])
    print('4/3:',Rg[idx43])
    print('2/1:',Rg[idx2])
    #print('5/4:',Rg[idx54])
    #print('6/5:',Rg[idx65])
    #print('7/6:',Rg[idx76])
    #print('8/7:',Rg[idx87])
    #print('4/1:',Rg[idx4])
    #print('3/1:',Rg[idx3])
    #print('8/3:',Rg[idx83])
    #print('5/2:',Rg[idx52])
    #print('7/3:',Rg[idx73])
    #print('7/2:',Rg[idx72])
    #print('5/1:',Rg[idx5])

"""Transformations R(psi) and psi(R)."""    
def R2psi(x):
    psi = m.DS.psirz[m.grid_mag_center_z_index,m.grid_mag_center_r_index:m.r_sep_index_right_intersection]
    psin = (psi-m.DS.simag)/(m.DS.sibry - m.DS.simag)
    R = m.rr[m.grid_mag_center_r_index:m.r_sep_index_right_intersection]
    psi_r = interp1d(R, psin, bounds_error=False, fill_value='extrapolate')
    return psi_r(x)

def psin2R(x):
    psi = m.DS.psirz[m.grid_mag_center_z_index,m.grid_mag_center_r_index:m.r_sep_index_right_intersection]
    psin = (psi-m.DS.simag)/(m.DS.psirz[m.grid_mag_center_z_index,m.r_sep_index_right_intersection]-m.DS.simag)
    R = m.rr[m.grid_mag_center_r_index:m.r_sep_index_right_intersection]
    r_psin = interp1d(psin, R, bounds_error=False, fill_value='extrapolate')
    return r_psin(x)


def psi2R(x):
    psi = m.DS.psirz[m.grid_mag_center_z_index,m.grid_mag_center_r_index:m.r_sep_index_right_intersection]
    psin = (psi-m.DS.simag)/(m.DS.sibry - m.DS.simag)
    R = m.rr[m.grid_mag_center_r_index:m.r_sep_index_right_intersection]
    r_psi = interp1d(psin, R, bounds_error=False, fill_value='extrapolate')
    return r_psi(x)


def compose2(f,g):
    """Returns composition of two functions"""
    return lambda x: f(g(x))

"""Functions that give transformation between x and R for figure coordinates."""
def x2R(x):
    xs = np.linspace(0,1,psi_grid_dim)
    R_x = interp1d(xs, Rg, bounds_error=False, fill_value='extrapolate')
    return R_x(x)

def R2x(x):
    xs = np.linspace(0,1,psi_grid_dim)
    X_R = interp1d(Rg, xs, bounds_error=False, fill_value='extrapolate')
    return X_R(x)


def q_plots():
    '''Q Check plots'''
    plt.title(r'q vs. $\Psi_N$')
    plt.plot(psi_norm,m.DS.qpsi)
    plt.axhline(y=1., color='r', linestyle='--')
    plt.axhline(y=3./2., color='y', linestyle='-.')
    plt.show()
    

    fig, ax = plt.subplots()
    if m.DS.qpsi[-1]>=0:
        ax.plot(psi_norm,m.DS.qpsi)
    else:
        ax.plot(psi_norm,[-x for x in m.DS.qpsi])
    ax.set_ylabel(r'q')
    ax.set_xlabel(r'$\Psi_N$')
    ax.set_xlim(0,1)
    secax = ax.twiny()
    newlabel = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    newlabel = [psi2R(x) for x in newlabel]
    newlabel = [round(float(x),2) for x in newlabel]
    newpos = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    #newpos = [R2psi(p) for p in newlabel]
    ax.set_xticks(newpos)
    secax.set_xticks(newpos) 
    secax.set_xticklabels(newlabel)
    secax.set_xlabel(r'R')
    plt.show()
    
        

    plt.title(r'$I_p$ vs. $\Psi_N$')
    plt.plot(psi_norm, -2*np.pi*np.array(m.DS.fpol))
    plt.show()
    

def norm_shape():
    plt.plot(Rs/a,Zs/a)
    plt.show()
