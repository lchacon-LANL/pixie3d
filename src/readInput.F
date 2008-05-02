c readInput
c######################################################################
      subroutine readInput

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use equilibrium

      use grid

      use timeStepping

      use newtongm

      use nlfunction_setup

      use precond_setup

      use iosetup

      use constants

      use transport_params

      use icond

      implicit none

c Call variables

c Local variables

      integer       :: dim(1),loc(1),ioerr
      real(8)       :: prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
      real(8)       :: prndtl,hrtmn,temp_ratio
      character*(3) :: bcs(6)
      type(dim_pack):: gp1,gp2,gp3

      logical       :: nc_eom_f

c Namelists

      !Old namelist, for backward compatibility
      namelist /datin/ neqd,nxd,nyd,nzd,coords,bcs,xmax,ymax,zmax
     .                   ,xmin,ymin,zmin,gparams,mg_ratio,numerical_grid
     .                   ,npx,npy,npz
     .                ,ilevel,debug,debug_it
     .                ,nu,eta,dd,chi,gamma,prndtl,hrtmn,di,temp_ratio
     .                   ,nn_eta,aa_eta,c_hyper
     .                ,rtol,atol,stol,mf_eps,maxitnwt,tolgm,maxksp
     .                   ,iguess,maxitgm,global,method,damp,dt0,vol_wgt
     .                ,equil,dlambda,rshear,vparflow,vperflow,source
     .                   ,chk_src
     .                ,nh1,nh2,nh3,prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
     .                   ,odd,random
     .                ,precon,maxvcyc,nsweep,mgtol,iguess,precpass
     .                   ,asm_PC
     .                ,dt,cnfactor,tmax,dstep,timecorr,numtime,restart
     .                   ,ndstep,sm_flag,sm_pass,predictor,postprocess
     .                ,gp1,gp2,gp3,check_grid
     .                ,nc_eom_f,nc_eom_jxb,nc_eom_gp,nc_eom_v,solenoidal
     .                   ,k_si,solve_rho,advect,sym_st,adiabatic
     .                ,inputfile,recordfile,equ_file,dcon

      !Problem setup
      namelist /setup/ neqd,nxd,nyd,nzd,coords,bcs,xmax,ymax,zmax
     .                   ,xmin,ymin,zmin,gparams,mg_ratio,numerical_grid
     .                   ,npx,npy,npz,gp1,gp2,gp3,equ_file

      !Physics setup
      namelist /phys/  nu,eta,dd,chi,gamma,prndtl,hrtmn,di,temp_ratio
     .                   ,nn_eta,aa_eta,c_hyper
     .                   ,equil,dlambda,rshear,vparflow,vperflow

      !Spatial discretization
      namelist /fv/    nc_eom_f,nc_eom_jxb,nc_eom_gp,nc_eom_v,solenoidal
     .                   ,solve_rho,advect,sym_st,adiabatic

      !Temporal discretization
      namelist /ts/    dt,cnfactor,tmax,dstep,timecorr,numtime,restart
     .                   ,ndstep,sm_flag,sm_pass,predictor,postprocess
     .                   ,nh1,nh2,nh3,prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
     .                   ,odd,random,source,k_si

      !Newton-Krylov
      namelist /nks/   rtol,atol,stol,mf_eps,maxitnwt,tolgm,maxksp
     .                   ,maxitgm,global,method,damp,dt0,vol_wgt

      !Preconditioner
      namelist /pc/    precon,maxvcyc,nsweep,mgtol,iguess,precpass
     .                   ,asm_PC,iguess

      !I/O
      namelist /io/    ilevel,debug,debug_it,check_grid
     .                ,recordfile,dcon,chk_src

c ******************************************************************

c Begin program

c Set defaults

      !General setup
      neqd     = 8             ! Number of degrees of freedom (equations)

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

      gparams  = 0d0           ! Array with additional grid parameters
                               ! (grid-dependent; see setEquilibrium.F)

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
      tmax     = 0d0           ! Target time, in Alfven times.
      numtime  = -1            ! Number of time steps
      ndstep   = 0             ! # time steps between plots (if zero,
                               !        ndstep is calculated in code)
      dstep    = 0.            ! Time interval between plots (if zero,
                               !        dstep is calculated in code)

      restart  = .false.       ! Restarting flag
      timecorr = .true.        ! Time adaptive algorithm (based on Newton convergence)

      cnfactor = 0.5           ! Crank-Nicolson factor (implicit if <= 0.5)
      sm_flag  = 0             ! Time smoothing flag:
                               !   0 -> Theta scheme
                               !   1 -> Rannacher time stepping
                               !   2 -> BDF2
      sm_pass  = 2             ! Number of initial smoother passes for Rannacher TS

      !NK parameters
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

      !Preconditioner parameters
      mg_ratio = 2             ! MG coarsening ratio
      nsweep   = 5             ! Number of MG smoothing passes
      maxvcyc  = 1             ! Maximum number of MG V-cycles
      mgtol    = 1d-3          ! MG convergence tolerance

      precon   = 'id'          ! Type of preconditioner. Currently:
                               !   - 'id': identity (default)
                               !   - 's1': SI without flow
                               !   - 's2': SI with flow
      precpass = 1             ! Number of SI iterations in preconditioner
      asm_PC   = .false.       ! Whether we are doing additive Schwartz PC (in parallel)

      !Physics parameters
      nu       = 1d-3          ! Reynolds number
      eta      = 1d-3          ! Lundquist number
      dd       = 1d-3          ! Particle diffusion
      chi      = 1d-3          ! Thermal diffusivity

      hrtmn    = 0d0           ! Hartmann number (1/sqrt(nu*eta))
      prndtl   = 0d0           ! Prandtl number (nu/eta)

      nn_eta   = 0             ! Exponent for eta profile
      aa_eta   = 0d0           ! Coefficient for eta profile

      c_hyper  = 0d0           ! Coefficient for hyper-resistivity

      gamma    = 5./3.         ! Polytropic constant of plasma

      di       = 0d0           ! Hall parameter

      temp_ratio = 1d0         ! Ion-electron temperature ratio

      !Nonlinear function parameters
      k_si       = 0d0         ! SI constant

      nc_eom_jxb = .false.     ! Whether we use non-conservative form of jxB in EOM
      nc_eom_gp  = .false.     ! Whether we use non-conservative form of grad(p) in EOM
      nc_eom_f   = .false.     ! Implies both jxb and grad(p) in EOM
      nc_eom_v   = .false.     ! Whether we use non-conservative form of inertia in EOM
      solenoidal = .true.      ! Whether we use solenoidal discret. of Faraday's law
      solve_rho  = .true.      ! Whether we solver continuity equation or not
      sym_st     = .false.     ! Whether we use the symmetric form of the viscous stress
                               !   tensor or not
      adiabatic  = .false.     ! Whether we use adiabatic EoS or not
      vol_wgt    = .true.      ! Whether residual is volume weighed or not

      advect     = 2           ! Type of advective scheme (currently available only
                               ! for scalars):
                               !    1 -> upwind,
                               !    2 -> centered (ZIP),
                               !    3 -> QUICK,
                               !    4 -> SMART,
                               !    5 -> smooth SMART,
                               !    6 -> centered 4th order

      !Initial condition
      equil    = ''            ! Type of equilibrium (see setEquilibrium.F)
cc      eq_params= 0d0           ! Equilibrium parameters (see      "       )
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
      nh1      = 1             ! Harmonic number for perturbation in axis 1
      nh2      = 1             ! Harmonic number for perturbation in axis 2
      nh3      = 1             ! Harmonic number for perturbation in axis 3
      odd      = .false.       ! Symmetry of perturbation
      random   = .false.       ! Random initialization if true

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

      debug      = .false.     ! PC debugging flag
      debug_it   = 1           ! Newton iteration level for debugging

      recordfile ='record.bin' ! Default output file
      equ_file   ='pixie3d.equ'! Default equilibrium file (when needed)
      prt_file   ='pixie3d.eig'! Default perturbation file (when needed)

      dcon       = .false.     ! If using VMEC input, whether to dump DCON output

c Read initialization parameters

      inputfile  ='pixie3d.in' ! Default input file

      open(unit=uinput,file=inputfile,status='old')

      read(uinput,datin,iostat=ioerr)    !Compatibility mode

      if (ioerr /= 0) then
        !Reopen input file
        close(unit=uinput)
        open(unit=uinput,file=inputfile,status='old')

        ioerr = 0
        read(uinput,setup,iostat=ioerr)
        if (ioerr/=0) call pstop('readInput','Problem reading "setup"')

        read(uinput,phys,iostat=ioerr)
        if (ioerr/=0) call pstop('readInput','Problem reading "phys"')

        read(uinput,fv,iostat=ioerr)
        if (ioerr/=0) call pstop('readInput','Problem reading "fv"')

        read(uinput,ts,iostat=ioerr)
        if (ioerr/=0) call pstop('readInput','Problem reading "ts"')

        read(uinput,nks,iostat=ioerr)
        if (ioerr/=0) call pstop('readInput','Problem reading "nks"')

        read(uinput,pc,iostat=ioerr)
        if (ioerr/=0) call pstop('readInput','Problem reading "pc"')

        read(uinput,io,iostat=ioerr)
        if (ioerr/=0) call pstop('readInput','Problem reading "io"')
      endif

      close(unit=uinput)

c Obtain eta, nu, dd from Prandtl, Hartmann

      if (hrtmn.gt.0d0.and.prndtl.gt.0d0) then
        nu  = sqrt(prndtl)/hrtmn
        eta = 1d0/hrtmn/sqrt(prndtl)
      endif

c Find hyperresistivity

      if (di > 0d0) heta = c_hyper

c Find alpha_p

      a_p = 1d0+temp_ratio

c Initialize grid packing structure

      g_pack%dim(1) = gp1
      g_pack%dim(2) = gp2
      g_pack%dim(3) = gp3

c Consistency check

cc      if (precon == 'id') iguess = 0

      !Adequate BCs for collapsed dimension
      if (nxd == 1) bcs(1:2) = 'per'
      if (nyd == 1) bcs(3:4) = 'per'
      if (nzd == 1) bcs(5:6) = 'per'

      !Non-conservative EOM
      if (nc_eom_f) then
        nc_eom_jxb = .true.
        nc_eom_gp  = .true.
      endif

c Map perturbations

      pert(IRHO)= prho
      pert(IVX) = pvx
      pert(IVY) = pvy
      pert(IVZ) = pvz
#if !defined(vec_pot)
      pert(IBX) = pbx
      pert(IBY) = pby
      pert(IBZ) = pbz
#else
      pert(IAX) = pbx
      pert(IAY) = pby
      pert(IAZ) = pbz
#endif
      pert(ITMP)= ptemp

c Translate boundary conditions

      bcond = -1

      where (bcs == 'def')
        bcond = DEF
      elsewhere (bcs == 'per')
        bcond = PER
      elsewhere (bcs == 'spt')
        bcond = SP
      elsewhere (bcs == 'sym')
        bcond = SYM
      elsewhere (bcs == 'equ')
        bcond = EQU
      end where

      if (minval(bcond) < 0) then
        loc = 1 - mod(minloc(bcond),2)
        dim = 1+(minloc(bcond) - (1+loc))/2
        write (*,*) 'Error in defining boundary conditions'
        write (*,*) 'Undefined boundary condition in axis',dim,
     .              ', location',loc
        write (*,*) 'Aborting'
        write (*,*) bcond
        stop
      endif

c End program

      end subroutine readInput

c readGraphicsInput
c######################################################################
      subroutine readGraphicsInput

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use iosetup

      use graphics_variables

      use transport_params

      implicit none

c Call variables

c Local variables

      integer    :: i,ioerr

      !For compatibility
      namelist /graphdef/ sel_diag,sel_graph,ndplot,dplot,hdf_plot
     .                   ,prof_conf,cont_conf,clean,E0,B0
     .                   ,iplot,jplot,kplot


      namelist /plots/    sel_diag,sel_graph,ndplot,dplot,hdf_plot
     .                   ,prof_conf,cont_conf,clean,E0,B0
     .                   ,iplot,jplot,kplot

c Begin program

c Read computation initializations (external)

      call readInput

c Graphics defaults

      ndplot = 0          ! Postprocessing interval (# of time steps; integer)
      dplot  = 0d0        !       "          "      (time interval; real)
      hdf_plot =.false.   ! Whether an HDF5 file is to be created

      !Line profile configuration (XDRAW)
      prof_conf%line   = 1          ! Direction (1 -> x, 2 -> y, 3 -> z)
      prof_conf%label  ='x'         ! Label
      prof_conf%coords = (/1,1,1/)  ! Line coordinates 

      !Contour slice configuration (XDRAW)
      cont_conf%plane  = 3          ! Normal to cut plane (1 -> x, 2 -> y, 3 -> z)
      cont_conf%label  = (/'x','y'/)! Contour plot axes labels
      cont_conf%coords =(/1,1,1/)   ! Plane coordinates along normal

      sel_diag = 0   !Array of size "xdraw_cont_lim (currently set to 16)
                     !  indicating time histories to be dumped (see
                     !  drawgamma.in --generated after postprocessing--
                     !  or diagnostics.F for information on available variables)

      sel_graph = (/ (i,i=1,xdraw_cont_lim) /)  !(obsolete)
                     !Selects graphics to be shown with XDRAW contour plotter.


c Read graphics initialization parameters

      open(unit=uinput,file=inputfile,status='old')

      read(uinput,graphdef,iostat=ioerr)
      if (ioerr /= 0) then
        close(uinput)
        open(unit=uinput,file=inputfile,status='old')
        read(uinput,plots,iostat=ioerr)
        if (ioerr /= 0) then
          write (*,*) 'Problem reading "plots"' ; stop
        endif
      endif

      close(unit=uinput)

c End program

      end subroutine readGraphicsInput
