import matplotlib.pyplot as plt
import pylab as py
import numpy as np
from imp import reload
import pixie_read_st as pxr
from matplotlib import animation, rc
import matplotlib.animation as animation
import types
from IPython.display import HTML
import master_read as m
import efit_plots as e
import os
import sys

plt.style.use('ggplot')
plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'heavy'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12
plt.rcParams['text.usetex']=True

#eqdsk_file = "/net/scratch3/chacon/pixie3d/EFIT/ITER/ITER3-chipar/3d/SN_fr_11_sh.geqdsk"
#filepath = "/net/scratch4/giannis_kx/pixie3d/iter/int_kink/11/11_visc_old_nodiff/11_visc_old_nodiff.scratch/"
eqdsk_file = "/users/giannis_kx/eqdsks/eqdsk_9MA_SS.Gpolevoa"
filepath = "/net/scratch3/giannis_kx/pixie3d/iter/iter3d/db_new/dt2.scratch./"
t_end = 204

def eqdsk_info():
    """Extracts parameters from the eqdsk file."""
    global a
    sys.stdout = open(os.devnull, 'w')
    m.read_geqdsk(eqdsk_file)
    m.struct_hor_ax_det()
    e.gridsNunits()
    e.intersections()
    sys.stdout = sys.__stdout__

# Bug fix for animating contour plots    
def setvisible(self,vis):
    for c in self.collections: c.set_visible(vis)
def setanimated(self,ani):
    for c in self.collections: c.set_animated(ani)
        
            
def main():
    eqdsk_info()
    
    # Extract the separatrix
    r_sep = [x for x in m.DS.rbbbs]
    z_sep = [x for x in m.DS.zbbbs]
    
    # Data loading
    pxr.pixieload(filepath + "pixie3d.h5")
    Te = pxr.load_array(0,6,0,t_end)
    psi_pol = pxr.load_array(3,4,0,1)
    
    # Toroidal averaging
    psi_pol_n0 = np.mean(psi_pol,axis=2)
    Te_n0 = np.mean(Te,axis=2)
    
    # Extracting magnetic axis and peak temperature
    MA = np.unravel_index(np.argmin(psi_pol_n0[:,:,0]),(psi_pol.shape[0],psi_pol.shape[1]))
    Tpeak = Te_n0[MA[0],MA[1],0]
    Tnorm = 1/Tpeak
    
    # Animation
    fig = plt.figure(figsize=(10,10))
    plt.plot(r_sep,z_sep,color='red',linewidth=3,linestyle='--')
    plt.axes().set_aspect("equal")
    plt.xlabel("R (m)")
    plt.ylabel("Z (m)")

    ims = []
    for i in range(86):
        im = plt.contourf(pxr.X[:,:,0]*e.a,pxr.Z[:,:,0]*e.a,np.log10(Te[:,:,0,i]*Tnorm),80,cmap="coolwarm")
        text = 't = '+str(i*10)+ r'$\;\tau_A$'
        an = plt.annotate(text, xy=(0.8, 0.94), xycoords='axes fraction',fontsize=14)

        #################################################################
        ## Bug fix for Quad Contour set not having attribute 'set_visible'
        im.set_visible = types.MethodType(setvisible,im)
        im.set_animated = types.MethodType(setanimated,im)
        im.axes = plt.gca()
        im.figure=fig
        ####################################################################

        ims.append(im.collections+[an])
    cbar = plt.colorbar()
    cbar.set_label(r"$\log_{10}\left(\frac{T_e}{T_o}\right)$",rotation=0,labelpad=42)

    ani = animation.ArtistAnimation(fig, ims, interval=200, blit=False,repeat_delay=100)
    
    return ani
        
if __name__ == "__main__":
    main()
