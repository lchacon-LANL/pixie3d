# README #

>0.0 Table of contents
>
>1.0 Introduction
>
>2.0 Requirements
>
> > 2.1 Compilers
>
> > 2.2 External libraries
>
> > 2.3 Graphics: Xdraw, Visit
>
>3.0 Installation instructions
>
>4.0 Compilation instructions
>
> > 4.1 Serial compilation
> >
> > 4.2 Parallel compilation
> >
> > 4.3 Special Makefile flags
>
> 5.0 Using PIXIE3D
>
> > 5.1 Input file
> >
> > 5.2 Serial execution
> >
> > 5.3 Parallel execution
> >
> > 5.4 Postprocessing
>
> 6.0 References
>
> Appendix 1: Default values for input parameters for "pixie3d.x"
>
> Appendix 2: Default values for input parameters for "pixplot.x" (graphics postprocessor)

----------------------------------------------------------------------

# Introduction

This readme is intended to describe the requirements to compile and
run the code PIXIE3D. PIXIE3D solves the extended magnetohydro-
dynamics (MHD) equations in 3D arbitrary geometries using fully
implicit Newton-Krylov algorithms. The code employs the open-source
PETSc library for parallelization, and the HDF5 library for data
management. It is primarily used to study plasma behavior in magnetic
confinement fusion devices. The spatial discretization of the MHD
equations in PIXIE3D is described in Ref. [1]. Details of the temporal
discretization can be found in Refs. [2-5].

# Requirements

We proceed to describe the minimum requirements to compile and run
PIXIE3D.

## Compilers

The serial version of PIXIE3D (pixie3d.x) only requires a Fortran
90/95 compiler. PIXIE3D has been compiled and tested with Absoft,
Lahey, and Intel f95 compilers. Others will be added in time.

The parallel version of PIXIE3D (pixie3d.petsc.x) requires an ANSI C
compiler in addition to the fortran compiler. In addition, certain
libraries need to be present. We proceed to discuss these next.

## External libraries

The serial version of PIXIE3D does not require any external
libraries. The parallel version currently uses PETSC 2.2.0 for
parallelization (http://www.anl??).

Advanced graphics postprocessing in 3D requires HDF5 libraries
(v. 1.8 works; http://www.hdfgroup.org/HDF5/). Such libraries allows
the user to generate an HDF5 format file, readable by the LLNL
graphics interpreter Visit (see next section).

All these libraries contain documentation as to how to compile them
and test them, and these are not reproduced here in the interest of
brevity.

## Graphics

PIXIE3D comes with its own postprocessor, called pixplot.x. This file
will take output from PIXIE3D and generate various files with time traces,
contour plots (viewable with Xdraw) and fully 3D plots (in HDF5 format,
which can be viewed with the free Visit software). We proceed to describe 
the visualization software below.

### Xdraw

Xdraw is a free, simple (yet powerful) 1D and 2D plotter (included in
this electronic distribution). The postprocessor "pixplot.x" will
always generate xdraw files for 2D contour plots, even in 3D runs. The
user can choose how to slice the solution among the main logical
coordinate axes (see Sec. 5.1 below).

Xdraw is invoked with "xdraw [keyword]", where "keyword" can be "car"
(for cartesian coordinate plots, "cnv" for contravariant component
plots, "cov" for covariant component plots, "diag" for diagnostic
plots, "pert" for perturbation plots, and "gamma" for time histories
of various diagnostics (described in Sec. 5.1 below). When invoked,
xdraw reads the file draw[keyword].in to read its configuration, and
plot accordingly. These draw*.in files are generated automatically by
pixplot.x.

### Visit

This is an open source, freely available visualization solfware
developed at LLNL. An interface has been developed for PIXIE3D to
render 2D and 3D plots generated by the latter. Version 1.5.2 or later
is required, and can be obtained from http://www.llnl.gov/VisIt/home.html

# Installation instructions

PIXIE3D is easy to install. Just execute:

> tar xzvf pixie3d.tar.gz

to generate the directory structure, and change directories into
"pixie3d".  From there, type:

> make setup

to setup all required compilation files.

# Compilation instructions

In order to compile PIXIE3D, the user will need to edit the Makefile
in the main distribution directory to setup some system-dependent
variables, such as the installation directory for PETSC (PETSC_DIR; if
a parallel compilation is required) and that for the HDF5 distribution
(HDF5_HOME; if an HDF5 file is to be generated after postprocessing).

PIXIE3D has been successfully compiled in various environments (Linux,
AIX) with several compilers (Absoft, Lahey, Intel, gfortran, PGI). Before
attempting compilation, one must manually setup the workstation configuration
in "common/make/make.mach.inc". There are various examples there one can
follow. A typical workstation configuration is as follows:

     ifeq ($(findstring sith,$(HOST)),sith)
       FC = pgf95

       BOPT = O
       PETSC_DIR ?= /sw/sith/petsc/3.0.0/centos5.5_pgi10.9_opt/petsc-3.0.0-p10
       PETSC_ARCH ?= sith-opt

       HDF5 = t
       HDF5_DIR = /sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.2
       HDF5_HOME = $(HDF5_DIR)
       HDF5_LIBS = $(HDF5_FLIB)
       HDF5_MOD  = $(HDF5_HOME)/lib -I$(HDF5_DIR)/include

       ADIOS = t
       ADIOS_HOME = $(ADIOS_DIR)
       ADIOS_LIBS = $(ADIOS_INC) $(ADIOS_FLIB)
       ADIOS_VER  = 1.3

       MACHINE = sith
     endif

Typical variables to be defined in the configuration file include:

   * MACHINE: name of machine
   * FC: fortran compiler
   * HDF5=t,f: whether to use HDF5.
   * BOPT: if defined, use PETSc.
   * ADIOS=t,f: whether to use ADIOS I/O library.
   * NETCDF=t,f: whether to use NETCDF
   * LIBS: LAPACK and BLAS libraries (in serial verion only).

## Serial compilation

For a serial compilation, type

> make distclean
> make

For testing serially, type:

> make tests-b

## Parallel compilation

Linking against the PETSC library requires one to set the variable
BOPT to either "g" (debugging) or "O" (optimized). For production
runs, "O" should normally be used, and one should type

> make BOPT=O distclean
> make BOPT=O all

to generate pixie3d.petsc.x and pixplot.petsc.x (the parallel
postprocessor).  For testing in parallel, type:

> make BOPT=O tests-b

## Options when invoking "make"

There are several specific options for the make command:

OPT: to set the optimization level. Can be set to "g" (debugging), "O"
(optimized), "p" (profiling), "s" (static; for Absoft only), and "v"
(verbose), and combinations (e.g. one can set OPT=Opv to combine
optimization with profiling and verbose options).

VECPOT: if set, the vector potential version of PIXIE3D will be
compiled, generating the executables "pixie3d_a.x" and "pixplot_a.x"

SAMR: if set, compiles AMR version of PIXIE3D using SAMRAI. Requires
SAMRAI and SAMRAI_UTILS. This is in an experimental stage.

# Using PIXIE3D

## Input file

Before execution, the user should set the input file, "pixie3d.in". 
Defaults are provided within PIXIE3D for all quantities,
which are explained in the readInput.F file. An excerpt of this file
is included in App. 1 (for pixie3d.x) and App. 2 (for pixplot.x).

## Serial execution

Serial execution requires the user simply to type:

> ./pixie3d.x (or ./pixie3d_a.x for the vector potential version).

After completion of the run, the user should invoke:

> ./pixplot.x (or ./pixplot_a.x)

for postprocessing. Plots can be viewed by invoking xdraw or visit (if
an HDF5 file has been created).

## Parallel execution

Parallel execution will be largely dependent on the MPI distribution
that has been used to compile PETSc, and on the availability of a
queue management system. With mpich-1.2.5.2, one would
type, for instance:

> mpiexec [mpi options] pixie3d.petsc.x [PIXIE3D options] [petsc options]

Several mpirun options are available according to the MPI
distribution. Similarly with PETSC. Specific PIXIE3D options include:

     OUTPUT CONTROL:                                     
        -ilevel <ilevel>: level of output info
        -test: whether we are in test mode (suppresses output)

     TIME STEPPING:                                      
        -nmax <nmax>: number of time steps     
        -tmax <tmax>: final time                         
                                                         
     PARALLEL DA:                                        
        -npx: processors in X direction                  
        -npy: processors in Y direction                  
        -npz: processors in Z direction                  

## Postprocessing

Postprocessing is done with pixplot.x (or any of the following
variants: pixplot_a.x for vector potential version; pixplot.petsc.x
for parallel version; and pixplot_a.petsc.x for parallel vector
potential version).

pixplot.x is configured via the namelist "graphdef", also located in
the file "pixie3d.in". The components of the namelist are described in
App. 2. 

The serial postprocessor can be used to postprocess serial or parallel
runs (provided the information fits in memory). However, the parallel
postprocessor cannot be used to postprocess serial runs (unless it
employs only one processor).

The parallel postprocessor does not require any runtime options. If
postprocessing is done in parallel, one needs to use the same number
of processors as were used with pixie3d.petsc.x. Furthermore, XDRAW
contour plots are not generated (although time histories still are),
and the HDF5 file pixie3d.h5 is generated instead. This file can be
interpreted using the freely available graphics software Visit.

# References

[1] L. Chac??n, "A non-staggered, conservative, , finite-volume scheme
for 3D implicit extended magnetohydrodynamics in curvilinear
geometries," Comput. Phys. Comm., 163 (3), 143-171 (2004).

[2] L. Chac??n, D. A. Knoll, J. M. Finn, "An implicit, nonlinear
reduced resistive MHD solver," J. Comput. Phys., 178, 15-36 (2002).

[3] L. Chac??n, D. A. Knoll, "A 2D high-beta Hall MHD implicit
nonlinear solver," J. Comput. Phys., 188 (2), 573-592 (2003).

[4] L. Chac??n, "A fully implicit 3D extended MHD algorithm," 32nd EPS
Conference on Plasma Physics, Tarragona, Spain, June 27-July 1, 2005.

[5] L. Chac??n, D. A. Knoll, "A fully implicit 3D extended MHD
algorithm," 33rd EPS Conference on Plasma Physics, Rome, Italy, June
19-23, 2006.

# Appendix 1: Default values for input parameters for "pixie3d.x"

      !General setup
      nxd      = 64            ! Mesh points in x-direction
      nyd      = 64            ! Mesh points in y-direction
      nzd      = 64            ! Mesh points in z-direction

      coords   = 'car'         ! Coordinate system (see grid_anal_mod.F)

      xmax     = 1d0           ! Length in x-direction
      ymax     = 1d0           ! Length in y-direction
      zmax     = 1d0           ! Length in z-direction

      xmin     = 0d0           ! Length in x-direction
      ymin     = 0d0           ! Length in y-direction
      zmin     = 0d0           ! Length in z-direction

      gparams  = 0d0           ! Array with additional grid parameters (grid-dependent)

      numerical_grid = .false. ! Whether grid metrics are calculated numerically (.true.)
                               !   or analytically.

      bcs      = (/ 'def','def','per','per','per','per' /) 
                               ! Defines topological boundary conditions:
                               !   + 'def' = default
                               !     (see applyBoundaryCondtions.F)
                               !   + 'per' = periodic
                               !   + 'spt' = singular point
                               !   + 'sym' = symmetry 
                               !     (homogeneous Neumann/Dirichlet
                               !      for even/odd quantities)

      !Time stepping
      dt       = 5.            ! Time step (if zero, dt is calculated in code)
      tmax     = -1d0          ! Target time, in Alfven times.
      numtime  = -1            ! Number of time steps
      ndstep   = 0             ! # time steps between plots (if zero,
                               !        ndstep is calculated in code)
      dstep    = 0.            ! Time interval between plots (if zero,
                               !        dstep is calculated in code)

      restart  = .false.       ! Restarting flag
      timecorr = .true.        ! Time adaptive algorithm (based on Newton convergence)
      postprocess = .false.    ! Whether to run solution postprocessing after time step

      eigen    = .false.       ! Whether to run PIXIE3D in eigensolver mode
      eig_dt   = 1d0           ! Time step for eigenvalue dynamical system

      cnfactor = 0.5           ! Crank-Nicolson factor (implicit if <= 0.5)
      sm_flag  = 0             ! Time smoothing flag:
                               !   0 -> Theta scheme
                               !   1 -> Rannacher time stepping
                               !   2 -> BDF2
      sm_pass  = 2             ! Number of initial smoother passes for Rannacher TS

      !NK parameters
      fpa      = .false.       ! Whether to use Carlson's Fixed-Point-Accelerator instead of JFNK
      tolgm    = 8d-1          ! Inexact Newton parameter
      rtol     = 1d-4          ! Newton relative convergence tolerance
      atol     = 0d0           ! Newton absolute convergence tolerance
      stol     = 0d0           ! Newton update convergence tolerance
      mf_eps   = 1d-6          ! Newtom matrix-free differencing parameter
      maxitnwt = 20            ! Maximum number of Newton its. (if zero, maxitnwt
                               !        is determined in code)
      maxksp   = 15            ! Maximum krylov subspace dimension
      maxitgm  = maxksp        ! Maximum GMRES iterations
      method   = 1             ! Inexact Newton method:
                               !   + 0: constant forcing
                               !   + other: adaptive (Eisenstat-Walker)

      global   = 0             ! Newton's globalization method:
                               !   + 0 --> no globalization
                               !   + 1 --> linesearch backtracking 
                               !   + 2 --> pseudo-transient

      damp     = 1d0           ! Damping parameter in Newton (see nk_mod.f doc))
      dt0      = 1d30          ! Initial pseudo-transient time step (" " ")
      iguess   = 0             ! Whether preconditioner is used to give initial
                               !   guess to GMRES (when =1) 
      !Physics parameters
      nu       = 1d-3          ! Reynolds number
      eta      = 1d-3          ! Lundquist number
      dd       = 1d-3          ! Particle diffusion
      chi      = 1d-3          ! Perpendicular thermal diffusivity
      chi_par  = 0d0           ! Parallel thermal diffusivity

      hrtmn    = 0d0           ! Hartmann number (1/sqrt(nu*eta))
      prndtl   = 0d0           ! Prandtl number (nu/eta)

      di       = 0d0           ! Ion inertia parameter
      de       = 0d0           ! Electron inertia parameter

      aa_eta   = 0d0           ! Coefficient for eta profile
      bb_eta   = 0d0           ! Coefficient for eta profile
      cc_eta   = 1d0           ! Coefficient for eta profile
      aa_nu    = 0d0           ! Coefficient for nu profile
      bb_nu    = 0d0           ! Coefficient for nu profile
      cc_nu    = 1d0           ! Coefficient for nu profile

      heta     = 0d0           ! Coefficient for hyper-resistivity

      gamma    = 5./3.         ! Polytropic constant of plasma

      adiabatic  = .false.     ! Whether we use adiabatic EoS or not

      temp_ratio = 1d0         ! Ion-electron temperature ratio, Ti/Te

      lagrangian = .false.     ! Whether we perform Lagrangian step for Te

      spitzer  = .false.       ! Whether we are using Spitzer resistivity

      E0       = 0d0           ! Boundary electric field (cov)
      B0       = 0d0           ! Boundary magnetic field (cnv)

      !Nonlinear function parameters
      nc_eom_jxb = .false.     ! Whether we use non-conservative form of jxB in EOM
      nc_eom_gp  = .false.     ! Whether we use non-conservative form of grad(p) in EOM
      nc_eom_f   = .false.     ! Implies both jxb and grad(p) in EOM
      nc_eom_v   = .false.     ! Whether we use non-conservative form of inertia in EOM
      solenoidal = .true.      ! Whether we use solenoidal discret. of Faraday's law
      solve_rho  = .true.      ! Whether we solve continuity equation or not
      solve_prs  = .false.     ! Whether we solve for electron pressure or temperature
      sym_st     = .false.     ! Whether we use the symmetric form of the viscous stress
                               !   tensor or not
      vol_wgt    = .true.      ! Whether residual is volume weighed or not

      advect     = 2           ! Type of advective scheme for scalars
                               ! Available options:
                               !    1 -> upwind,
                               !    2 -> ZIP (centered),
                               !    3 -> QUICK,
                               !    4 -> SMART,
                               !    5 -> smooth SMART,
                               !    6 -> centered 4th order

      v_advect   = 0           ! Type of advective scheme for vectors
                               !   (if zero, determine from advect below)

      ion_hall   = .false.     ! Whether to use new Hall formulation or not
      fake_ve    = .true.      ! Whether to use fake electron velocity in PIe in EOM

      slava      = .false.     ! Whether to use Slava Lukin's Hall MHD implementation
      noise      = .false.     ! Whether to add white noise to EOM
      noise_lev  = 0d0         ! Noise level

      !Preconditioner parameters
      pc_type    = 'id'          ! Type of preconditioner. Currently:
                               !   - 'id': identity (default)
                               !   - 's1': SI without flow
                               !   - 's2': SI with flow
      pc_iter    = 1           ! Number of SI iterations in preconditioner
      pc_tol     = 1d-3        ! PC solvers convergence tolerance

      pc_asm     = .false.     ! Whether we are doing additive Schwartz PC (in parallel)
      pc_debug   = .false.     ! PC debugging flag
      pc_divclean= .false.     ! Whether to perform divergence cleaning in PC or not

      mg_ratio = 2             ! MG coarsening ratio
      mg_vcyc  = 1             ! Maximum number of MG V-cycles
      mg_ores  = 0             ! Restriction order for MG
      mg_oprol = 2             ! Prolongation order for MG

      sm_iter  = 5             ! Number of MG smoothing passes
      sm_zebra_relax = .false. ! Whether to use ZEBRA relaxation

      !Initial condition
      equil    = ''            ! Type of equilibrium (see setEquilibrium.F)
      eq_params= 0d0           ! Equilibrium parameters (see      "       )
      dlambda  = .2            ! Characteristic equilibrium scale length
      rshear   = 1.            ! Ratio of magnetic to fluid sheet thicknesses
      vparflow = 0.            ! Maximum parallel fluid flow
      vperflow = 0.            ! Maximum perpendicular fluid flow
      source   = .true.        ! Impose source to provide equilibrium
      chk_src  = .false.       ! Whether to check source (e.g., to check equilibria)

      prho     = 0d0           ! Density perturbation
      pvx      = 0d0           ! Vx perturbation
      pvy      = 0d0           ! Vy perturbation
      pvz      = 0d0           ! Vz perturbation
      pbx      = 0d0           ! Bx perturbation
      pby      = 0d0           ! By perturbation
      pbz      = 0d0           ! Bz perturbation
      ptemp    = 0d0           ! Temperature perturbation
      nh1      = 0             ! Starting harmonic number for perturbation in axis 1
      nh2      = 0             ! Starting harmonic number for perturbation in axis 2
      nh3      = 0             ! Starting harmonic number for perturbation in axis 3
      npw1     = 1             ! Number of harmonics for perturbation in axis 1
      npw2     = 1             ! Number of harmonics for perturbation in axis 2
      npw3     = 1             ! Number of harmonics for perturbation in axis 3
      odd      = .false.       ! Symmetry of perturbation
      random   = .false.       ! Random initialization if true

      br_pert_bc = 0d0         ! Br perturbation magnitude at boundary (for pinch)
                               ! Perturbation modes defined in grid_params 1 and 2
                               !   for all geometries (hel, cyl, tor)

      !Logical grid configuration
      gp1%pack = .false.       ! Do not pack in X-direction
      gp1%xp   = 0d0
      gp1%dx0  = 0d0
      gp2%pack = .false.       ! Do not pack in Y-direction
      gp2%xp   = 0d0
      gp2%dx0  = 0d0
      gp3%pack = .false.       ! Do not pack in Z-direction
      gp3%xp   = 0d0
      gp3%dx0  = 0d0
                               ! To select packing, one needs to set the fields
                               ! of gp1, gp2, gp3 as follows:
                               !    gp* = pack,xp,dx0
     .                         ! where:
                               !   + pack (logical): whether to pack
                               !   + xp (real): where to pack
                               !   + dx0 (real): initial grid spacing (at xp)

      check_grid = .false.     ! Whether to dump grid info or not

      !I/O parameters
      ilevel     = 0           ! Level of output information:
                               !   -  0: time step      level
                               !   -  1: Newton solver  level (basic)
                               !   -  2: Newton solver  level (advanced)
                               !   -  3: Krylov solver  level
                               !   -  4: Preconditioner level (basic)
                               !   - >4: Preconditioner level (advanced)
                               ! Each level encompasses previous ones.

      test       = .false.     ! Whether we are performing tests (do not dump namelist)

      equ_file   ='pixie3d.equ'! Default equilibrium file (when needed)
      prt_file   ='pixie3d.eig'! Default perturbation file (when needed)
      map_file   ='pixie3d.map'! Default map          file (when needed)

      dcon       = .false.     ! If using VMEC input, whether to dump DCON output

# Appendix 2: Default values for input parameters for "pixplot.x" (graphics postprocessor):

      ndplot = 0          ! Postprocessing interval (# of time steps; integer)
      dplot  = 0d0        !       "          "      (time interval; real)
      hdf_plot =.false.   ! Whether an HDF5 file is to be created
      adios_plot =.false. ! Whether an ADIOS-BP file is to be created
      xdraw_plot = .true. ! Whether to dump XDRAW files

      !Local (point) diagnostic configuration (XDRAW)
      iplot = 1
      jplot = 1
      kplot = 1

      !Line profile configuration (XDRAW)
      prof_conf%line   = 1          ! Direction (1 -> x, 2 -> y, 3 -> z)
      prof_conf%label  ='x'         ! Label
      prof_conf%coords = (/1,1,1/)  ! Line coordinates 

      !Contour slice configuration (XDRAW)
      cont_conf%plane  = 3          ! Normal to cut plane (1 -> x, 2 -> y, 3 -> z)
      cont_conf%label  = (/'x','y'/)! Contour plot axes labels
      cont_conf%coords =(/1,1,1/)   ! Plane coordinates along normal

      sel_diag = 0   !Array of size "xdraw_cont_lim" (currently set to 16)
                     !  indicating time histories to be dumped (see
                     !  drawgamma.in --generated after postprocessing--
                     !  or diagnostics.F for information on available variables)

      sel_graph = (/ (i,i=1,xdraw_cont_lim) /)  !(obsolete)
                     !Selects graphics to be shown with XDRAW contour plotter.