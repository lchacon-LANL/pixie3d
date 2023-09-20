import numpy as np
from scipy.io import FortranFile

filename = "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/11_visc_old_nodiff/ftle/t-416/ic_r400xu200_dr1e-8.bin"

Nr = 400 # r-points on main grid
Nu = 200 # u-points on main grid

rmin = 0.15 # don't put it to zero!
rmax = 0.8 # avoid exiting the confinement

dx = 1e-8 # auxilliary grid spacing

r_mg = np.linspace(0.15,0.82,Nr)
u_mg = np.linspace(0.001,2*np.pi,Nu)

icf = FortranFile(filename,"w")
for r in r_mg:
    for u in u_mg:
        icf.write_record(np.array([r, u, 0]))
        icf.write_record(np.array([r + dx, u, 0]))
        icf.write_record(np.array([r - dx, u, 0]))
        icf.write_record(np.array([r, u + dx, 0]))
        icf.write_record(np.array([r, u - dx, 0]))

icf.close()
