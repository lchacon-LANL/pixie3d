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

      use nk_setup

      use nlfunction_setup

      use precond_setup

      use app_iosetup

      use equilibrium

      implicit none

c Call variables

c Local variables

      integer       :: dim(1),loc(1),ioerr
      real(8)       :: prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
      real(8)       :: prndtl,hrtmn
      character*(3) :: bcs(6)
      type(grid_pack):: gp1,gp2,gp3

      logical       :: nc_eom_f

      !Legacy parameters
      character(2)  :: precon
      integer       :: maxvcyc,precpass,nsweep,ores,oprol,debug_it
     .                ,pc_debug_it,pc_max_iter
      logical       :: zebra_relax,hall_new,lag_pinch_bc,debug
     .                ,mg_zebra_relax
      real(8)       :: mgtol

c Namelists

      !Old namelist, for backward compatibility
      namelist /datin/ nxd,nyd,nzd,coords,bcs,xmax,ymax,zmax
     .                   ,xmin,ymin,zmin,gparams,mg_ratio,numerical_grid
     .                   ,npx,npy,npz,gp1,gp2,gp3,check_grid
     .                ,ilevel,debug,test
     .                ,nu,eta,dd,chi,chi_par,gamma,prndtl,hrtmn,di,de
     .                   ,temp_ratio,aa_eta,bb_eta,cc_eta,E0,B0
     .                   ,aa_nu,bb_nu,cc_nu,heta,spitzer,br_pert_bc
     .                ,rtol,atol,stol,mf_eps,maxitnwt,tolgm,maxksp,fpa
     .                   ,iguess,maxitgm,global,method,damp,dt0,vol_wgt
     .                ,equil,dlambda,rshear,vparflow,vperflow,source
     .                   ,chk_src,eq_params
     .                ,nh1,nh2,nh3,npw1,npw2,npw3
     .                   ,odd,random,prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
     .                ,pc_type,pc_tol,pc_iter,pc_asm,pc_debug
     .                   ,pc_divclean,mg_vcyc,sm_iter
     .                   ,sm_zebra_relax,mg_ores,mg_oprol
     .                ,dt,cnfactor,tmax,dstep,timecorr,numtime,restart
     .                   ,ndstep,sm_flag,sm_pass,predictor,postprocess
     .                   ,eigen,eig_dt
     .                ,nc_eom_f,nc_eom_jxb,nc_eom_gp,nc_eom_v,solenoidal
     .                   ,solve_rho,advect,v_advect,sym_st,adiabatic
     .                   ,solve_prs,ion_hall,fake_ve,slava,noise
     .                   ,noise_lev,lagrangian
     .                ,equ_file,map_file,dcon
     .           !Legacy parameters
     .                ,precon,mgtol,maxvcyc,precpass,nsweep,hall_new
     .                   ,zebra_relax,ores,oprol,debug,debug_it
     .                   ,lag_pinch_bc,recordfile,pc_debug_it
     .                   ,mg_zebra_relax,pc_max_iter

      !Problem setup
      namelist /setup/ nxd,nyd,nzd,coords,bcs,xmax,ymax,zmax
     .                   ,xmin,ymin,zmin,gparams,mg_ratio,numerical_grid
     .                   ,npx,npy,npz,gp1,gp2,gp3

      !Physics setup
      namelist /phys/  nu,eta,dd,chi,chi_par,gamma,prndtl,hrtmn,di,de
     .                   ,temp_ratio,aa_eta,bb_eta,cc_eta,E0,B0
     .                   ,aa_nu,bb_nu,cc_nu,heta,adiabatic
     .                   ,equil,dlambda,rshear,vparflow,vperflow
     .                   ,eq_params,lagrangian,spitzer,br_pert_bc

      !Spatial discretization
      namelist /fv/    nc_eom_f,nc_eom_jxb,nc_eom_gp,nc_eom_v,solenoidal
     .                   ,solve_rho,solve_prs,advect,sym_st,ion_hall
     .                   ,fake_ve,slava,noise,noise_lev,v_advect

      !Temporal discretization
      namelist /ts/    dt,cnfactor,tmax,dstep,timecorr,numtime,restart
     .                   ,ndstep,sm_flag,sm_pass,predictor,postprocess
     .                   ,nh1,nh2,nh3,npw1,npw2,npw3
     .                   ,prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
     .                   ,odd,random,source,eigen,eig_dt,iguess

      !Newton-Krylov
      namelist /nks/   rtol,atol,stol,mf_eps,maxitnwt,tolgm,maxksp
     .                   ,maxitgm,global,method,damp,dt0,vol_wgt,fpa

      !Preconditioner
      namelist /pc/    pc_type,mg_vcyc,sm_iter,pc_tol,pc_iter
     .                   ,pc_asm,sm_zebra_relax,pc_debug
     .                   ,pc_divclean

      !I/O
      namelist /i_o/   ilevel,debug,check_grid,equ_file,map_file
     .                   ,dcon,chk_src,test

c ******************************************************************

c Begin program

#if defined(petsc)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
#endif

c Set fundamental problem parameters (number of dep. variables, aux. variables)

      call setDepVarDims(2,2)   !2 scalars, 2 vectors (vectors not operational yet)

#if !defined(vec_pot)
      call setAuxVarDims(8,22)  !8 scalars, 22 vectors
#else
      call setAuxVarDims(8,25)  !8 scalars, 25 vectors
#endif

c Set defaults

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
cc      lag_pinch_bc = .false.   ! Whether to lag the paramagnetic pinch BC

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

      !Legacy parameters
      precon   = ''
      mgtol    = 0d0
      maxvcyc  = 0
      precpass = 0
      nsweep   = 0
      zebra_relax = .false.
      mg_zebra_relax = .false.
      ores     = 0
      oprol    = 0
      debug    = .false.
      debug_it = -1
      hall_new = .false.
      pc_max_iter = 0

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

        read(uinput,i_o,iostat=ioerr)
        if (ioerr/=0) call pstop('readInput','Problem reading "i_o"')
      endif

      close(unit=uinput)

c Check legacy parameters

      if (hall_new     )  ion_hall = .true.

      if (precon /= '' )  pc_type = precon
      if (mgtol  /= 0d0)  pc_tol  = mgtol
      if (precpass/=0  )  pc_iter = precpass
      if (debug        )  pc_debug = debug

      if (maxvcyc/= 0  )  mg_vcyc = maxvcyc
      if (zebra_relax.or.mg_zebra_relax  )  sm_zebra_relax = .true.
      if (ores   /= 0  )  mg_ores = ores
      if (oprol  /= 0  )  mg_oprol= oprol

      if (nsweep /= 0  )  sm_iter = nsweep

      if (pc_max_iter>0)  pc_iter = pc_max_iter

c Obtain eta, nu, dd from Prandtl, Hartmann

      if (hrtmn.gt.0d0.and.prndtl.gt.0d0) then
        nu  = sqrt(prndtl)/hrtmn
        eta = 1d0/hrtmn/sqrt(prndtl)
      endif

c Find alpha_p=p_t/p_e; p_t=p_i+p_e

      a_p = 1d0+temp_ratio

c Initialize grid packing structure

      g_pack(1) = gp1
      g_pack(2) = gp2
      g_pack(3) = gp3

c Consistency check

cc      if (di == 0d0) heta = 0d0

      !Adequate BCs for collapsed dimension
      if (nxd == 1) bcs(1:2) = 'per'
      if (nyd == 1) bcs(3:4) = 'per'
      if (nzd == 1) bcs(5:6) = 'per'

      !NC EOM for non-cartesian geometry
      if (coords /= 'car') nc_eom_f = .true.

      !Non-conservative EOM
      if (nc_eom_f) then
        nc_eom_jxb = .true.
        nc_eom_gp  = .true.
      endif

      if (cnfactor == 1d0) ion_hall = .false.

      if (.not.ion_hall) fake_ve = .false.

      if (ion_hall) slava = .false.

      if (.not.solenoidal) nc_eom_jxb = .true. !To avoid parallel forces

      if (eigen) then
        vol_wgt = .false.
        source  = .false.
        cnfactor = 0d0
        sm_flag  = 0
      endif

      if (tst_flg) test = .true.

      if (restart) dcon = .false.

cc      if (adiabatic.and.gamma /= 1d0.and.(.not.solve_rho)) then
cccc        solve_rho = .false.
cc        postprocess = .true.
cc      endif

      if (coords == 'hel') then   !Special case for 2D helical equilibria
         nh2 = gparams(1)
         nh1 = nh2
      endif

      if ((nh1 == 0 .and. nh2 == 0 .and. nh3 == 0)
     .   .and.(    equil /= 'ppn3d'
     .        .and.equil /= 'p3nsl'
     .        .and.equil /= 'adv1'
     .        .and.equil /= 'adv2')
     .        .and.equil /= 'sod' ) then
        call pstop('readInput','No sinusoidal perturbation specified')
      endif

      if (chi_par == 0d0) chi_par = chi

      if (di == 0d0) de = 0d0

      if (de > di) then
        call pstop('readInput','Unphysical value for di and/or de')
      endif

      !Performs divergence cleaning at postprocessing stage
      if (solenoidal.and.(.not.pc_divclean)) postprocess = .true.

      !Vector advection
      if (v_advect == 0) then
        if (nc_eom_v) then
          v_advect = min(advect,2) !Non-conservative EOM
        else
          v_advect = advect
        endif
      endif

      !Perform Lagrangian step at postprocessing time
      if (lagrangian) postprocess = .true.

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
      pert(ITMP)= a_p*ptemp  !Transform to total pressure

c Translate boundary conditions

      bcond = -1

      where (bcs == 'def')
        bcond = DEF
      elsewhere (bcs == 'per')
        bcond = PER
      elsewhere (bcs == 'spt')
        bcond = SP
      elsewhere (bcs == 'sym')
        bcond = FSYM  !Symmetry at faces (for conservation)
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

c Write initialization parameters to standard output

      if (my_rank == 0 .and. (.not. test)) then
        write (*,*) '-------------------'
        write (*,*) 'SETUP configuration'
        write (*,*) '-------------------'
        write (*,setup)

        write (*,*) '------------------'
        write (*,*) 'PHYS configuration'
        write (*,*) '------------------'
        write (*,phys)

        write (*,*) '----------------- '
        write (*,*) 'FV  configuration'
        write (*,*) '----------------- '
        write (*,fv)

        write (*,*) '----------------- '
        write (*,*) 'TS  configuration'
        write (*,*) '----------------- '
        write (*,ts)

        write (*,*) '----------------- '
        write (*,*) 'NKS configuration'
        write (*,*) '----------------- '
        write (*,nks)

        write (*,*) '----------------- '
        write (*,*) 'PC  configuration'
        write (*,*) '----------------- '
        write (*,pc)

        write (*,*) '----------------- '
        write (*,*) 'IO  configuration'
        write (*,*) '----------------- '
        write (*,i_o)
      endif

c End program

      end subroutine readInput

c readGraphicsInput
c######################################################################
      subroutine readGraphicsInput

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use app_iosetup

      use graphics_io

      use equilibrium

      implicit none

c Call variables

c Local variables

      integer :: i,ioerr

      !For compatibility
      namelist /graphdef/ sel_diag,sel_graph,ndplot,dplot,hdf_plot
     .                   ,prof_conf,cont_conf,clean
     .                   ,iplot,jplot,kplot,poincare,xdraw_plot
     .                   ,adios_plot


      namelist /plots/    sel_diag,sel_graph,ndplot,dplot,hdf_plot
     .                   ,prof_conf,cont_conf,clean
     .                   ,iplot,jplot,kplot,poincare,xdraw_plot
     .                   ,adios_plot

c Begin program

c Read computation initializations (external)

      call readInput

c Graphics defaults

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

c Consistency check

      if (iplot > nxd .or. jplot > nyd .or. kplot > nzd) then
        call pstop('readGraphicsInput'
     .            ,'Plot indices out of range')
      endif

c End program

      end subroutine readGraphicsInput
