c readInput
c######################################################################
      subroutine readInput

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

#if defined(plotting) || defined(_equ)
      use ts_setup
#else
      use timeStepping
#endif

      use nlfunction_setup

      use PC_setup

      use local_BCS_variables, ONLY: normal_flow_bc,sheath_flow_bc
     $                              ,v00_ExB_bc,res_wall_bc

      use rw_bc, ONLY: rw_fft_tupdate,rw_implicit_bc,rw_exp_ssp_bc
     $                ,rw_exp_pc_bc,rw_exp_fix,tau_w,r_v,rw_bcs

      use app_iosetup

#if !defined(_equ)
      use re, ONLY: re_diag
#endif

      implicit none

c Call variables

c Local variables

#if defined(_equ)
      logical :: re_diag
#endif

      integer       :: dim(1),loc(1),ioerr,mgcf
      real(8)       :: prho,pvx,pvy,pvz,pbx,pby,pbz,pax,pay,paz,ptmp
      real(8)       :: prndtl,hrtmn
      character*(3) :: bcs(6)
      type(grid_pack):: gp1,gp2,gp3

      logical       :: nc_eom_f,mg_coarse_solve

      integer       :: mg_coarse_size

      !Legacy parameters
      character(2)  :: precon
      integer       :: maxvcyc,precpass,nsweep,ores,oprol,debug_it
     .                ,pc_debug_it,pc_max_iter
      logical       :: zebra_relax,hall_new,lag_pinch_bc,debug
     .                ,mg_zebra_relax,restart,timecorr,fake_ve
     $                ,zero_ExB_tang,regularize_A,zero_vpar_bc
     $                ,sheath_bc
      real(8)       :: mgtol

c Namelists

      !Old namelist, for backward compatibility
      namelist /datin/ nxd,nyd,nzd,coords,bcs,xmax,ymax,zmax
     .                   ,xmin,ymin,zmin,gparams,mg_ratio,numerical_grid
     .                   ,npx,npy,npz,gp1,gp2,gp3,check_grid,equil
     .                   ,efit_coords,gstor_I_drive
     .                ,ilevel,debug,test
     .                ,nu,eta,dd,chi,chi_par,gamma,prndtl,hrtmn,di,de
     .                   ,temp_ratio,aa_eta,bb_eta,cc_eta,E0,B0,M0
     .                   ,aa_nu,bb_nu,cc_nu,heta,irrot_br_pert,eig_is_T
     .                   ,spitzer,chi_bgnskii,n_0,mi_mp,log_lamb,min_te
     .                   ,br_pert_bc,br_pert_phase,br_pert_freq
     .                   ,br_pert_bc_old,br_pert_phase_old,xport_order
     .                   ,r_v,tau_w,normal_flow_bc,d_rw,eta_rw
     .                   ,re_pkg,c_va,re_diag
     .                ,rtol,atol,stol,mf_eps,maxitnwt,tolgm,maxksp
     .                   ,iguess,maxitgm,global,method,damp,dt0,vol_wgt
     .                   ,pc_as_solver,fpa,fpa_it
     .                ,dlambda,rshear,vparflow,vperflow,source
     .                   ,chk_src,eq_params
     .                ,nh1,nh2,nh3,npw1,npw2,npw3
     .                   ,odd,random,prho,pvx,pvy,pvz,pbx,pby,pbz,ptmp
     .                   ,pax,pay,paz
     .                ,pc_type,pc_tol,pc_iter,pc_asm,pc_debug,pc_B_solve
     .                   ,pc_divclean,mg_vcyc,sm_iter,sm_type,sm_omega
     .                   ,sm_zebra_relax,sm_zebra_it,sm_line_relax
     .                   ,mg_ores,mg_oprol,mg_coarse_solve,mgcf
     .                   ,pc_gm_driver,mg_coarse_size
     .                ,dt,cnfactor,tmax,dstep,numtime,ts_restart
     .                   ,trestart
     .                   ,ndstep,sm_flag,sm_pass,predictor,postprocess
     .                   ,eigen,eig_dt,perturb_rstrt,ts_adapt_dt
     .                   ,limit_dt_flow_cfl,ts_ramp_up_dt
     .                   ,ts_ramp_up_dt_fraction
     .                   ,ts_ramp_up_dt_recovery,refine_eq
     .                ,nc_eom_f,nc_eom_jxb,nc_eom_gp,nc_eom_v,solenoidal
     .                   ,solve_rho,advect,v_advect,sym_st,adiabatic
     .                   ,adiabatic_rho,solve_prs,no_visc_heating
     .                   ,use_p_eom,use_p_diff,ion_hall,fake_ve,slava
     .                   ,noise,noise_lev,lagrangian,no_eom_divpe
     .                   ,vlap_etaj,nu_is_rhonu,heta_is_rhoheta
     .                   ,subtract_E0,e_cons_vstar,post_smooth_B
     .                   ,use_p_classic_eq,post_filter_hel
     .                   ,binom_alph,v00_ExB_bc,sheath_flow_bc
     .                   ,B_angle_threshold,brem_rad_sink
     .                   ,rw_fft_tupdate,rw_implicit_bc,rw_exp_ssp_bc
     .                   ,rw_exp_pc_bc,rw_exp_fix
     .                ,prt_file,equ_file,map_file,dcon,min_sheath_flow
     .                !Legacy parameters
     .                ,precon,mgtol,maxvcyc,precpass,nsweep,hall_new
     .                   ,zebra_relax,ores,oprol,debug,debug_it
     .                   ,lag_pinch_bc,recordfile,pc_debug_it
     .                   ,mg_zebra_relax,pc_max_iter,timecorr,restart
     .                   ,aa_chi,bb_chi ! Daniele, 02/07/2014
     .                   ,cc_chi,dd_chi ! Daniele, 02/17/2014
     .                   ,zero_ExB_tang,regularize_A,zero_vpar_bc
     .                   ,sheath_bc
      
      !Problem setup
      namelist /setup/ nxd,nyd,nzd,coords,bcs,xmax,ymax,zmax
     .                   ,xmin,ymin,zmin,gparams,mg_ratio,numerical_grid
     .                   ,npx,npy,npz,gp1,gp2,gp3,efit_coords

      !Physics setup
      namelist /phys/  nu,eta,dd,chi,chi_par,gamma,prndtl,hrtmn,di,de
     .                   ,temp_ratio,aa_eta,bb_eta,cc_eta,E0,B0,M0
     .                   ,aa_nu,bb_nu,cc_nu,heta,adiabatic,adiabatic_rho
     .                   ,dlambda,rshear,vparflow,vperflow,xport_order
     .                   ,eq_params,lagrangian,irrot_br_pert,eig_is_T
     .                   ,spitzer,chi_bgnskii,nu_is_rhonu,min_te
     .                   ,br_pert_bc,br_pert_phase,br_pert_freq
     .                   ,br_pert_bc_old,br_pert_phase_old
     .                   ,nh1,nh2,nh3,npw1,npw2,npw3,random,odd
     .                   ,prho,pvx,pvy,pvz,pbx,pby,pbz,ptmp,pax,pay,paz
     .                   ,equil,gstor_I_drive,refine_eq
     .                   ,aa_chi,bb_chi ! Daniele, 02/07/2014
     .                   ,cc_chi,dd_chi ! Daniele, 02/17/2014
     .                   ,aa_dff,bb_dff,cc_dff
     .                   ,r_v,tau_w,zero_vpar_bc,d_rw,eta_rw
     .                   ,sheath_flow_bc,B_angle_threshold
     .                   ,min_sheath_flow
     .                   ,brem_rad_sink,brem_rad
     .                   ,re_pkg,c_va,re_diag
     .                   ,n_0,mi_mp,log_lamb
     $                   ,no_visc_heating

      !Spatial discretization
      namelist /fv/    nc_eom_f,nc_eom_jxb,nc_eom_gp,nc_eom_v,solenoidal
     .                   ,solve_rho,use_p_eom,use_p_diff,advect,sym_st
     .                   ,ion_hall,slava,noise,noise_lev,post_smooth_B
     .                   ,v_advect,no_eom_divpe,vlap_etaj,nu_is_rhonu
     .                   ,heta_is_rhoheta,subtract_E0,e_cons_vstar
     .                   ,use_p_classic_eq,post_filter_hel
     .                   ,solve_prs,binom_alph,xport_order
     .                   ,rw_implicit_bc,rw_exp_ssp_bc,rw_exp_pc_bc
     .                   ,rw_exp_fix,rw_fft_tupdate,v00_ExB_bc
      
      !Temporal discretization
      namelist /ts/    dt,cnfactor,tmax,dstep,ts_adapt_dt,numtime
     .                   ,ndstep,sm_flag,sm_pass,predictor,postprocess
     .                   ,source,eigen,eig_dt,perturb_rstrt,ts_restart
     .                   ,restart,trestart,limit_dt_flow_cfl
     .                   ,ts_ramp_up_dt
     .                   ,ts_ramp_up_dt_fraction,ts_ramp_up_dt_recovery

      !Newton-Krylov
      namelist /nks/   rtol,atol,stol,mf_eps,maxitnwt,tolgm,maxksp
     .                   ,maxitgm,global,method,damp,dt0,vol_wgt
     .                   ,iguess,pc_as_solver,fpa,fpa_it

      !Preconditioner
      namelist /pc/    pc_type,pc_tol,pc_iter,pc_debug,pc_asm,mg_vcyc
     .                   ,mg_ores,mg_oprol,sm_iter,sm_zebra_relax
     .                   ,sm_zebra_it,sm_line_relax,sm_type,sm_omega
     .                   ,mg_coarse_solve,mgcf,pc_gm_driver
     .                   ,mg_coarse_size,pc_divclean,pc_B_solve

      !I/O
      namelist /i_o/   ilevel,check_grid,equ_file,map_file
     .                   ,dcon,chk_src,test,prt_file

c ******************************************************************

c Begin program

cc      call set_bc_debug(.true.)

#if defined(petsc)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
#endif

c Set fundamental problem parameters (number of dep. variables, aux. variables)

      call setDepVarDims(2,2)   !2 scalars, 2 vectors (vectors not operational yet)

#if !defined(vec_pot)
      call setAuxVarDims(13,22)  !13 scalars, 22 vectors
#else
      call setAuxVarDims(13,26)  !13 scalars, 26 vectors
#endif

c Set defaults

      call set_defaults

c Read initialization parameters

      call read_namelists

c Dump initial IO headers
      
      if (.not.test) call output_header

c Obtain eta, nu from Prandtl, Hartmann

      if (hrtmn.gt.0d0.and.prndtl.gt.0d0) then
        nu  = sqrt(prndtl)/hrtmn
        eta = 1d0/hrtmn/sqrt(prndtl)
      endif

c Find alpha_p=T_t/T_e; T_t=T_i+T_e

      a_p = 1d0+temp_ratio

c Initialize MG parameters

#if defined(petsc) && defined(coarse_MG)
      if (mgcf > 0) call MG_set_proc_coarsen_factor(mgcf)
      call MG_set_coarse_solve(mg_coarse_solve)
#endif

      call MG_set_coarse_size(mg_coarse_size)

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
      pert(IAX) = pax
      pert(IAY) = pay
      pert(IAZ) = paz
#endif
      pert(IPRS)= a_p*ptmp  !Transform to total pressure

c Consistency check

      call consistency_chk

c Init boundary conditions

      ioerr = initBCs(bcs)

      if (ioerr /= 0)
     .     call critical('readInput','Error in initialization of BCs')

c Write initialization parameters to standard output

cc      if (.not.test) call write_namelists

c Initialize graphics

      ioerr=readGraphicsConfig(in_file)

      if (ioerr /= 0.and.my_rank == 0.and.(.not.test)) then
        call warning('  Problem reading "graphics" namelist')
      endif

c Debug profile stage

#if defined(petsc)
cc      iDBG=add_profile_stage("Debug")
#endif

c End program

      contains

c     read_namelists
c     #####################################################################
      subroutine read_namelists

      implicit none

      open(unit=uinput,file=in_file,status='old')
      read(uinput,datin,iostat=ioerr)    !Compatibility mode
      close(unit=uinput)

      if (ioerr /= 0) then
        !Reopen input file

        ioerr = 0
        open(unit=uinput,file=in_file,status='old')
        read(uinput,setup,iostat=ioerr)
        close(unit=uinput)
        if (ioerr/=0)
     $       call critical('readInput','Problem reading "setup"')

        open(unit=uinput,file=in_file,status='old')
        read(uinput,phys,iostat=ioerr)
        close(unit=uinput)
        if (ioerr/=0)
     $       call critical('readInput','Problem reading "phys"')

        open(unit=uinput,file=in_file,status='old')
        read(uinput,fv,iostat=ioerr)
        close(unit=uinput)
        if (ioerr/=0) call critical('readInput','Problem reading "fv"')

        ioerr = readTSInput(in_file,.false.)
        if (ioerr/=0) then
          open(unit=uinput,file=in_file,status='old')
          read(uinput,ts,iostat=ioerr)
          close(unit=uinput)
        endif
        if (ioerr/=0) call critical('readInput','Problem reading "ts"')

        open(unit=uinput,file=in_file,status='old')
        read(uinput,nks,iostat=ioerr)
        close(unit=uinput)
        if (ioerr/=0) call critical('readInput','Problem reading "nks"')

        open(unit=uinput,file=in_file,status='old')
        read(uinput,pc,iostat=ioerr)
        close(unit=uinput)
        if (ioerr/=0) call critical('readInput','Problem reading "pc"')

        open(unit=uinput,file=in_file,status='old')
        read(uinput,i_o,iostat=ioerr)
        close(unit=uinput)
        if (ioerr/=0) call critical('readInput','Problem reading "i_o"')
      endif

#if !defined(plotting) && !defined(samrai) && !defined(_equ)
      !Read projective integration configuration
      call read_PrjInt_Input(in_file)
#endif

      call set_grid_pack_config(gp1,gp2,gp3)
      
      end subroutine read_namelists

c     write_namelists
c     #####################################################################
      subroutine write_namelists

      implicit none

      if (my_rank == 0) then
        write (*,*) '------------------'
        write (*,*) 'GRID configuration'
        write (*,*) '------------------'
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

      end subroutine write_namelists

c     set_defaults
c     #####################################################################
      subroutine set_defaults

      implicit none

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
                               !   + 'fsm' = symmetry (at faces)
                               !     (homogeneous Neumann/Dirichlet
                               !      for even/odd quantities)

      mg_ratio = 2             ! MG coarsening ratio

      !Time stepping
      dt       = 5.            ! Time step (if zero, dt is calculated in code)
      tmax     = -1d0          ! Target time, in Alfven times.
      numtime  = -1            ! Number of time steps
      ndstep   = 0             ! # time steps between plots (if zero,
                               !        ndstep is calculated in code)
      dstep    = 0.            ! Time interval between plots (if zero,
                               !        dstep is calculated in code)

      ts_restart     = .false. ! Restarting flag
      ts_adapt_dt    = .true.  ! Time adaptive algorithm (based on Newton convergence)
      perturb_rstrt  = .false. ! Whether to apply perturbation after restart
      postprocess    = .false. ! Whether to run solution postprocessing after time step
      source         = .true.  ! Impose source to provide equilibrium

      eigen    = .false.       ! Whether to run PIXIE3D in eigensolver mode
      eig_dt   = 1d0           ! Time step for eigenvalue dynamical system

      cnfactor = 0.5           ! Crank-Nicolson factor (implicit if <= 0.5)
      sm_flag  = 0             ! Time smoothing flag:
                               !   0 -> Theta scheme
                               !   1 -> Rannacher time stepping
                               !   2 -> BDF2
      sm_pass  = 2             ! Number of initial smoother passes for Rannacher TS

      limit_dt_flow_cfl=.true. ! Whether to limit Dt according to flow CFL

      refine_eq = .false.      ! Whether we run PIXIE3D in equilibrium refinement mode

      !NK parameters
      fpa      = .false.       ! Whether to use Carlson's Fixed-Point-Accelerator instead of JFNK
      tolgm    = 8d-1          ! Inexact Newton parameter
      rtol     = 1d-4          ! Newton relative convergence tolerance
      atol     = 0d0           ! Newton absolute convergence tolerance
      stol     = 0d0           ! Newton update convergence tolerance
      mf_eps   = 1d-6          ! Newtom matrix-free differencing parameter
      maxitnwt = 10            ! Maximum number of Newton its. (if zero, maxitnwt
                               !        is determined in code)
      maxksp   = 15            ! Maximum krylov subspace dimension
      maxitgm  = maxksp        ! Maximum GMRES iterations
      method   = 1             ! Inexact Newton method:
                               !   + 0: constant forcing
                               !   + other: adaptive (Eisenstat-Walker)

      global   = 1             ! Newton's globalization method:
                               !   + 0 --> no globalization
                               !   + 1 --> linesearch backtracking 
                               !   + 2 --> pseudo-transient

      damp     = 1d0           ! Damping parameter in Newton (see nk_mod.f doc))
      dt0      = 1d30          ! Initial pseudo-transient time step (" " ")
      iguess   = 0             ! Whether preconditioner is used to give initial
                               !   guess to GMRES (when =1) 

      pc_as_solver= .false.    ! Whether we use PC as solver, instead of preconditioned GMRES

      vol_wgt  = .true.        ! Whether residual is volume weighed or not

      !Physics parameters
      n_0      = 1d20          ! Number density, SI units (m^(-3))
      mi_mp    = 1d0           ! Ion mass in proton units
      log_lamb = 1d1           ! Coulomb log

      nu       = 1d-3           ! Reynolds number
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
! Daniele, 02/07/2014
      aa_chi   = 1d0           ! Coefficient for chi profile
      bb_chi   = 1d0           ! Coefficient for chi profile
! Daniele, 02/17/2014
      cc_chi   = 1d0           ! Coefficient for chi profile
      dd_chi   = 0d0           ! Coefficient for chi profile
      aa_dff   = 0d0           ! Coefficient for diffusivity profile
      bb_dff   = 0d0           ! Coefficient for diffusivity profile
      cc_dff   = 1d0           ! Coefficient for diffusivity profile

      heta     = 0d0           ! Coefficient for hyper-resistivity (if negative, it is determined by mesh)

      spitzer  = .false.       ! Whether we are using Spitzer resistivity
      chi_bgnskii = .false.    ! Whether to use Braginskii conductivities

      no_visc_heating = .false.! Whether we allow for viscous heating or not
      
      gamma    = 5./3.         ! Polytropic constant of plasma

      temp_ratio = 1d0         ! Ion-electron temperature ratio, Ti/Te

      lagrangian = .false.     ! Whether we perform Lagrangian step for Te

      E0       = 0d0           ! Boundary electric field (cov)
      B0       = 0d0           ! Boundary magnetic field (cnv)
      M0       = 0d0           ! Momentum source

      tau_w    = 1d2           ! Resistive wall thin-shell time
      r_v      = 1.5           ! Radius of outer conducting shell
      normal_flow_bc = .true.  ! Whether we zero out parallel flow at bdry
      sheath_flow_bc = .false. ! Whether we apply sheath BC in v_parallel
      B_angle_threshold = 5d-2 ! Ratio of |B^1|/|B| at the wall that allows
                               ! the full Bohm flow to escape (sheath_flow_bc)
      brem_rad_sink = .false.  ! Whether Bremsstrahlung radiation sink is turned on
      brem_rad = 1.6d-7        ! Typical value of Bremsstrahlung coefficient
      min_sheath_flow = 0d0    ! Temporary coefficient of minimum sheath outflow (sheath_flow_bc = true)

      rw_fft_tupdate = .false. ! Whether we integrate RW eqs in Fourier space
      rw_implicit_bc = .false. ! Whether we use implicit RW BCs
      rw_exp_ssp_bc  = .true.  ! Whether we use explicit SSP RW BCs
      rw_exp_pc_bc   = .false. ! Whether we use explicit predictor-corrector BCs
      rw_exp_fix     = .true.  ! Whether we use instability fix for small tau_w

      re_pkg         = .false. ! Whether we activate the RE package
      c_va           = 1d2     ! Ratio speed of light / Alfven speed
      min_te         = 0d0     ! Reference min temperature (eV)
      
      !Nonlinear function parameters
      adiabatic  = .false.     ! Whether we use adiabatic EoS or not
      adiabatic_rho = .false.  ! Whether we use adiabatic rho formula
      nc_eom_jxb = .false.     ! Whether we use non-conservative form of jxB in EOM
      nc_eom_gp  = .false.     ! Whether we use non-conservative form of grad(p) in EOM
      nc_eom_f   = .false.     ! Implies both jxb and grad(p) in EOM
      nc_eom_v   = .false.     ! Whether we use non-conservative form of inertia in EOM
      no_eom_divpe = .true.    ! Whether we include div(Pe) in EOM
      solenoidal = .true.      ! Whether we use solenoidal discret. of Faraday's law
      solve_rho  = .true.      ! Whether we solve continuity equation or not
      solve_prs  = .true.      ! Whether we solve for pressure or temperature
      use_p_eom  = .false.     ! Whether we drive EOM with pressure or n*T
      use_p_diff = .false.     ! Whether we use pressure in diffusion term in pressure eq.
      use_p_classic_eq=.false. ! Whether we use classic form of prs/tmp eq (p.div(v), T.div(v))
      sym_st     = .false.     ! Whether we use the symmetric form of the viscous stress
                               !   tensor or not
      vlap_etaj  = .false.     ! Whether to use eta*lap(A) instead of eta*j in resistive term

      binom_alph = 0d0         ! Binomial smoothing constant of v in Ohm's law (0.5 to kill Nyquist mode)
      
      nu_is_rhonu= .false.     ! Whether we use rho*nu or just nu in ion viscosity term
      heta_is_rhoheta= .true.  ! Whether we use rho*heta or just heta in e viscosity term
      subtract_E0= .false.     ! Whether we subtract initial E in Ohm's law

      advect     = 2           ! Type of advective scheme for scalars
                               ! Available options:
                               !    1 -> upwind,
                               !    2 -> ZIP (centered),
                               !    3 -> QUICK,
                               !    4 -> SMART,
                               !    5 -> smooth SMART,
                               !    6 -> centered 4th order

      v_advect   = 0           ! Type of advective scheme for vectors (if zero, set to advect)

      ion_hall   = .false.     ! Whether to use new Hall formulation or not

      slava      = .false.     ! Whether to use Slava Lukin's Hall MHD implementation
      noise      = .false.     ! Whether to add white noise to EOM
      noise_lev  = 0d0         ! Noise level

      post_smooth_B = .false.  ! Whether to smooth B at postprocessing stage

      post_filter_hel = .false.! Whether to Fourier filter helical geometry sim

      xport_order   = 4        ! Order of anisotropic transport operator
                               ! Available options:
                               !    2 -> second order,
                               !    4 -> 4th order, 7-pt stencil,
                               !   -4 -> 4th order, 5-pt stencil

      !Preconditioner parameters
      pc_type    = 'id'        ! Type of preconditioner. Currently:
                               !   - 'id': identity (default)
                               !   - 's1': SI without flow
                               !   - 's2': SI with flow (iterative)
                               !   - 's3': SI with flow (commuted)
                               !   - 's4': EMHD split PC
      pc_iter    = 1           ! Number of SI iterations in preconditioner
      pc_tol     = 1d-3        ! PC solvers convergence tolerance

      pc_asm     = .false.     ! Whether we are doing additive Schwartz PC (in parallel)
      pc_debug   = .false.     ! PC debugging flag
      pc_divclean= .false.     ! Whether to perform divergence cleaning in PC or not
      pc_B_solve = .true.      ! Whether to perform B-solve in PC

      pc_gm_driver = .false.   ! Whether GMRES drives PC solve
      mg_vcyc  = 1             ! Maximum number of MG V-cycles
      mg_ores  = 0             ! Restriction order for MG
      mg_oprol = 1             ! Prolongation order for MG

      mg_coarse_solve = .true. ! Whether to perform MG solve beyond proc mesh
      mgcf = 0                 ! Default processor coarsening ratio per dimension in MG 
                               !   (if zero, leave default framework value unchanged)

      mg_coarse_size = 2       ! MG minimum mesh coarse size (integer exponent of coarsening ratio)

      sm_type  = 'jb'          ! Type of smoother ('jg','gs','gm')
      sm_omega = 7d-1          ! Damping for MG smoother (only when 'jb' is chosen)
      sm_iter  = 5             ! Number of MG smoothing passes
      sm_zebra_relax = .false. ! Whether to use ZEBRA relaxation
      sm_zebra_it    = 1       ! Number of ZEBRA relaxation passes
      sm_line_relax  = .false. ! Whether to perform line relaxation in MG

      !Initial condition
      equil    = ''            ! Type of equilibrium (see setEquilibrium.F)
      eq_params= 0d0           ! Equilibrium parameters (see      "       )
      dlambda  = .2            ! Characteristic equilibrium scale length
      rshear   = 1.            ! Ratio of magnetic to fluid sheet thicknesses
      vparflow = 0.            ! Maximum parallel fluid flow
      vperflow = 0.            ! Maximum perpendicular fluid flow
      chk_src  = .false.       ! Whether to check source (e.g., to check equilibria)

      gstor_I_drive = .true.   ! Whether we drive GS toroidal solver with I=Bphi*R or q-profile
      
      prho     = 0d0           ! Density perturbation
      pvx      = 0d0           ! Vx perturbation
      pvy      = 0d0           ! Vy perturbation
      pvz      = 0d0           ! Vz perturbation
      pbx      = 0d0           ! Bx perturbation
      pby      = 0d0           ! By perturbation
      pbz      = 0d0           ! Bz perturbation
      pax      = 0d0           ! Ax perturbation
      pay      = 0d0           ! Ay perturbation
      paz      = 0d0           ! Az perturbation
      ptmp     = 0d0           ! Temperature perturbation
      nh1      = 0             ! Starting harmonic number for perturbation in axis 1
      nh2      = 0             ! Starting harmonic number for perturbation in axis 2
      nh3      = 0             ! Starting harmonic number for perturbation in axis 3
      npw1     = 1             ! Number of harmonics for perturbation in axis 1
      npw2     = 1             ! Number of harmonics for perturbation in axis 2
      npw3     = 1             ! Number of harmonics for perturbation in axis 3
      odd      = .false.       ! Symmetry of perturbation
      random   = .false.       ! Random initialization if true

      br_pert_bc = 0d0         ! Br perturbation magnitude at boundary (for pinch)
                               ! Perturbation modes defined in nh2, nh3
                               !   for all geometries (hel, cyl, tor)

      br_pert_phase = 0d0      ! Phase of Br perturbation

      br_pert_freq  = 0d0      ! Frequency of Br perturbation

      irrot_br_pert = .false.  ! Whether to impose irrotational magnetic perturbation

      eig_is_T = .true.        ! Whether eigenvalue pert gives T (instead of P)
      
      !Logical grid configuration
      call init_grid_pack_config(gp1) ! To select packing, one needs to set the fields
      call init_grid_pack_config(gp2) ! of gp1, gp2, gp3 as follows:
      call init_grid_pack_config(gp3) !    gp* = pack,xp,dx0                          
                                      ! where:
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
      timecorr = .true.
      restart  = .false.
      zero_ExB_tang = .true.
      zero_vpar_bc = .false.
      sheath_bc = .false.

      end subroutine set_defaults

c     consistency_chk
c     #####################################################################
      subroutine consistency_chk

      implicit none

      integer :: i
      real(8) :: kk,rr

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !Legacy parameters checks!
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      
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

      if (.not.timecorr) ts_adapt_dt = timecorr

      if (restart) ts_restart = restart

      if (.not.zero_ExB_tang) normal_flow_bc = .false.
      if (zero_vpar_bc) normal_flow_bc = .false.

      if (sheath_bc) sheath_flow_bc = .true.
      
      !!!!!!!!!!!!!!!!!!!!!!
      !Time stepping checks!
      !!!!!!!!!!!!!!!!!!!!!!

      if (sm_flag == 2) cnfactor = 0d0  !Force BDF2

      if (tst_flg) test = .true.

      if (refine_eq) then
        ts_restart = .false.
        aa_nu = 0d0 ; bb_nu = 0d0 ; cc_nu = 0d0 !Flat nu profile
        dd  = 0d0         !To avoid dens decay
        chi = 0d0         !To avoid temp decay due to perpendicular xport
        nu  = 1d0         !To damp flow
        eta = 0d0         !To avoid resistive diffusion
        source = .false.  !To let equilibrium float
        pert = 0d0        !Zero out perturbations
        chk_src  = .true.
        sheath_flow_bc = .false.!Do not impose sheath flow BCs
      endif

      if (ts_restart) dcon = .false.

      !!!!!!!!!!!!!!!!!
      !BC setup checks!
      !!!!!!!!!!!!!!!!!

      !Face symmetric if "sym"
      where (bcs == 'sym') bcs = 'fsm'

      !Adequate BCs for collapsed dimension
      if (nxd == 1) bcs(1:2) = 'per'
      if (nyd == 1) bcs(3:4) = 'per'
      if (nzd == 1) bcs(5:6) = 'per'

      !Setup resistive wall logic
      where (bcs == 'rsw')
        rw_bcs = .true.
        bcs = 'def'
      end where

      res_wall_bc = .false.
      do i=1,size(rw_bcs)
        res_wall_bc = res_wall_bc.or.rw_bcs(i)
      enddo

      if (res_wall_bc.and.(.not.rw_implicit_bc)) then
        if (dt > 0.5*tau_w/nxd) then
          dt = 0.5*tau_w/nxd
          write (messg,*) " Adjusted dt for RW BC stability"
          call warning(messg)
        endif
      endif

      if (sheath_flow_bc) normal_flow_bc = .false. !Consider full ExB flow

      !!!!!!!!!!!!!!!!!!!!!!
      !Physics setup checks!
      !!!!!!!!!!!!!!!!!!!!!!

      !Non-conservative EOM
      if (nc_eom_f) then
        nc_eom_jxb = .true.
        nc_eom_gp  = .true.
      endif

      !Parallel transport checks
      if (chi_par == 0d0) chi_par = chi

      if (chi_par /= chi) solve_prs = .false.  !Solve for Te

      !Hall MHD checks
      if (di == 0d0) then
        de = 0d0
        ion_hall = .false.
      endif

cc      if (di == 0d0) heta = 0d0

      if (de > di) then
        call critical('readInput','Unphysical value for di and/or de')
      endif

      if (cnfactor==1d0.and.ion_hall) then
        call critical('readInput','Cannot do ion Hall explicitly')
      endif

      if (de>0d0.and.ion_hall) then
        call critical('readInput','Cannot do ion Hall with e inertia')
      endif

      if ((de>0d0).and.(cnfactor==0d0.and.sm_flag==2)) then
        call critical('readInput','Cannot do BDF2 with e inertia')
      endif

c$$$      if (ion_hall.and.(cnfactor==0d0.and.sm_flag==2)) then
c$$$        if (my_rank == 0)
c$$$     .       write (*,*) 'Doing NC velocity update for ion Hall'
c$$$        nc_eom_v = .true.
c$$$      endif

      electron_hall = (di > 0d0).and.(.not.ion_hall)

      !Setup for comparison with HiFi (S. Lukin)
      if (slava) then
        nu_is_rhonu = .true.
        heta_is_rhoheta = .true.
        e_cons_vstar = .false.
      else
#if defined(RFX)
        nu_is_rhonu = .false.
#endif
      endif

      if (ion_hall) no_eom_divpe = .false.       !Include divpe in ion EOM

cc      if (.not.solenoidal) nc_eom_jxb = .true. !To avoid parallel forces

      if (.not.adiabatic) adiabatic_rho = .false.

      if (adiabatic_rho) solve_rho = .false.

      if (.not.solve_prs) use_p_diff = .false.

      if (v_advect == 0) v_advect = advect

      if (br_pert_freq /= 0d0) irrot_br_pert = .true.
      if (br_pert_bc   /= 0d0) irrot_br_pert = .true.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Coordinate representation checks!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
cc      if (bcs(1) == 'spt') vlap_etaj = .false.

      !Adequate size for collapsed dimension
      if (coords /= 'car') then
        if (nxd == 1) then
          xmax = LMIN_NEQ_ONE ; xmin = 0d0
        endif
        if (nyd == 1) then
          ymax = LMIN_NEQ_ONE ; ymin = 0d0
        endif
        if (nzd == 1) then
          zmax = LMIN_NEQ_ONE ; zmin = 0d0
        endif
      endif

      !Helical geometry checks
      if (coords == 'hel') then  !Special case for 2D helical geometry
        if (nh2 == 0) nh2 = nint(gparams(1))
        nh1 = nh2
        
        kk = gparams(2)
        rr = gparams(3)
        if (abs(kk*rr-1d0*nint(kk*rr)) > 1d-13) then
          call critical("readInput"
     $                 ," Inconsistent helical input: k*R not integer")
        endif

        if (.not.test)
     $     call warning("  Hardwired nh1 for this coordinate system")
      endif

      !Cyl-tor geometries checks
      if (   coords == 'cyl'
     $   .or.coords == 'tor'
     $   .or.trim(equil) == '3dkai'
     $   .or.trim(equil) == 'ppn3d'
     $   .or.trim(equil) == 'p3nsl'
     $   .or.trim(equil) == '3drfp') then 
        if (nh2 == 0) nh2 = nint(gparams(1))
        if (nh3 == 0) nh3 = nint(gparams(2)*gparams(3))
        nh1 = nh2                
        if (.not.test)
     $     call warning("  Hardwired nh1 for this equilibrium")
      endif

      !!!!!!!!!!!!!!!!!!!!
      !Equilibrium checks!
      !!!!!!!!!!!!!!!!!!!!
      
      if (   trim(equil) == 'gstor'
     $   .or.trim(equil) == 'ohtor'
     $   .or.trim(equil) == 'stell') then
        if (nh2 == 0) nh2 =  nint(gparams(4))
        if (nh3 == 0) nh3 = -nint(gparams(5))
        nh1 = nh2
        if (.not.test)
     $     call warning("  Hardwired nh1 for this equilibrium")
      endif
      
      if (   trim(equil) == 'vmec'
     $   .or.trim(equil) == 'efit'
     $   .or.trim(equil) == '3drfp') then
        nh1 = nh2
        if (.not.test)
     $     call warning("  Hardwired nh1 for this equilibrium")
      endif

      if ((nh1 == 0 .and. nh2 == 0 .and. nh3 == 0)
     .   .and.(    equil /= 'ppn3d'
     .        .and.equil /= 'p3nsl'
     .        .and.equil(1:4) /= 'vela'
     .        .and.equil(1:4) /= 'tbbr'
     .        .and.equil /= 'rfp83'
     .        .and.equil /= 'bpnch'
     .        .and.equil /= 'bpchT'
     .        .and.equil /= 'adv1'
     .        .and.equil /= 'adv2')
     .        .and.equil /= 'sod' ) then
        call critical('readInput'
     $               ,'No sinusoidal perturbation specified')
      endif

      !!!!!!!!!!!!!!!!!!!!!!!
      !Postprocessing checks!
      !!!!!!!!!!!!!!!!!!!!!!!
      
#if defined(vec_pot)
      post_divclean = .false.
#else
      post_divclean = .not.((.not.pc_B_solve)
     .                .or.   pc_divclean
     .                .or.  (pc_type == 'id'.and.solenoidal)
     .                .or.  (nxd==1.and.nyd==1)
     .                .or.  (nxd==1.and.nzd==1)
     .                .or.  (nyd==1.and.nzd==1)
     .                .or.  (cnfactor == 1d0))
      post_divclean = post_divclean.or.post_smooth_B
#endif

      if (.not.postprocess)
     .     postprocess =    noise
     .                  .or.lagrangian
     .                  .or.post_divclean
     .                  .or.post_smooth_B
     .                  .or.(br_pert_freq /= 0d0)
     .                  .or.refine_eq
     .                  .or.post_filter_hel
     .                  .or.res_wall_bc

#if defined(pit)
      postprocess = .true.
#endif

      !!!!!!!!!!!!!!!
      !Solver checks!
      !!!!!!!!!!!!!!!

      !Vector potential exception
#if defined(vec_pot)
      if (.not.test) rtol = rtol*5d0
#endif

      !Use PC as linear solver
      if (pc_as_solver) then
         iguess = 1
         maxitgm = 0
cc         maxitnwt = 10
         rtol = 1d-1
      endif

      !PC
      if (pc_type == 's3') pc_B_solve = .true.   !Invert B for large-flow PC

      if (bcs(1) == 'spt') mg_oprol = 1

      if (test.and.
     .    electron_hall.and.pc_type /='id'
     .                 .and.pc_type /='s1'
     .                 .and.pc_type /='s4') pc_type = 's4'!Exception for tests

      if (test.and.
     .    ion_hall.and.pc_type /='id'
     .            .and.pc_type /='s1'
     .            .and.pc_type /='s2') pc_type = 's1'     !Exception for tests

cc      if (pc_type == 's1' .and. (electron_hall) .and. (di > 0d0))
cc     .                     pc_B_solve = .false.   !Do not invert B (or A) when doing Hall

      !Eigenvalue computation
      if (eigen) then
        vol_wgt = .false.
        source  = .false.
        cnfactor = 0d0
        sm_flag  = 0
      endif

      !Parallelization
#if defined(petsc)
      if (npx*npy*npz > np) then
         call critical("readInput","Too many procs specified")
      endif
#else
      npx = 1 ; npy = 1 ; npz = 1
#endif

      if (my_rank == 0) write (*,*)
      
      end subroutine consistency_chk

      end subroutine readInput
