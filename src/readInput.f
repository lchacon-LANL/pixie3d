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

      integer(4)    :: dim(1),loc(1)
      real(8)       :: prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
      real(8)       :: prndtl,hrtmn
      character*(3) :: bcs(6)

c Namelist

      namelist /datin/ neqd,nxd,nyd,nzd,coords,bcs,xmax,ymax,zmax
     .                   ,xmin,ymin,zmin,gparams,mg_ratio,numerical_grid
     .                ,ilevel
     .                ,nu,eta,dd,chi,gamma,prndtl,hrtmn
     .                ,rtol,atol,maxitnwt,tolgm,maxksp,iguess,maxitgm
     .                   ,global,method,damp,dt0
     .                ,equil,dlambda,rshear,vparflow,vperflow,source
     .                ,nh1,nh2,nh3,prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
     .                   ,odd,random
     .                ,precon,maxvcyc,nsweep,mgtol,iguess,precpass
     .                ,dt,cnfactor,tmax,dstep,timecorr,numtime,restart
     .                   ,ndstep,sm_pass

c ******************************************************************

c Begin program

c Set defaults

      !General setup
      neqd     = 8             ! Number of degrees of freedom (equations)
      nxd      = 64            ! Mesh points in x-direction
      nyd      = 64            ! Mesh points in y-direction
      nzd      = 64            ! Mesh points in z-direction

      coords   = 'car'         ! Coordinate system (see grid_mod.f)

      xmax     = 1d0           ! Length in x-direction
      ymax     = 1d0           ! Length in y-direction
      zmax     = 1d0           ! Length in z-direction

      xmin     = 0d0           ! Length in x-direction
      ymin     = 0d0           ! Length in y-direction
      zmin     = 0d0           ! Length in z-direction

      gparams  = 0d0           ! Array with additional grid parameters (grid-dependent)

      numerical_grid = .false. ! Whether grid metrics are calculated numerically (.true.)
                               !    or analytically.

      mg_ratio = 2             ! MG coarsening ratio

      bcs      = (/ 'def','def','per','per','per','per' /) 
                               ! Defines topological constraints. Convention:
                               !   + 'def' = default (set in code)
                               !   + 'per' = periodic
                               !   + 'spt' = singular point
                               !   + 'sym' = symmetry (homogeneous Neumann)

      !Time stepping
      dt       = 5.            ! Time step (if zero, dt is calculated in code)
      tmax     = 0d0           ! Target time, in Alfven times.
      numtime  = -1            ! Number of time steps
      ndstep   = 0             ! # time steps between plots (if zero,
                               !        ndstep is calculated in code)
      dstep    = 0.            ! Time interval between plots (if zero,
                               !        dstep is calculated in code)
      timecorr = .true.        ! Time adaptiveness based on Newton convergence
      cnfactor = -.48          ! Crank-Nicolson factor
      sm_pass  = 2             ! Number of initial smoother passes for CN

      !NK parameters
      tolgm    = 5.0d-2        ! Inexact Newton parameter (GMRES conv. tolerance)
      rtol     = 1.0d-4        ! Newton relative convergence tolerance
      atol     = 0d0           ! Newton absolute convergence tolerance
      maxitnwt = 0             ! Maximum number of Newton its. (if zero, maxitnwt
                               !        is determined in code)
      maxksp   = 15            ! Maximum krylov subspace dimension
      maxitgm  = maxksp        ! Maximum GMRES iterations
      method   = 0             ! Constant forcing parameter for inexact Newton
                               !        (see etak_meth in NewtonGmres doc)
      global   = 0             ! Do not use globalization (see global in NewtonGmres doc)
      damp     = 1d0           ! Damping parameter in Newton (see NewtonGmres doc))
      dt0      = 1d30          ! Initial pseudo-transient time step (" " ")
      iguess   = 1             ! Whether preconditioner is used to give initial
                               !        guess to GMRES (1)

      !Preconditioner parameters
      nsweep   = 5             ! Number of SGS sweeps
      maxvcyc  = 1             ! Maximum number of V-cycles
      mgtol    = 1d-3          ! MG convergence tolerance
      precon   = 'si'          ! Type of preconditioner
      precpass = 1             ! Number of SI iterations in preconditioner

      !Physics parameters
      nu       = 1d-3          ! Reynolds number
      eta      = 1d-3          ! Lundquist number
      dd       = 1d-3          ! Particle diffusion
      chi      = 1d-3          ! Thermal diffusivity

      hrtmn    = 0d0           ! Hartmann number (1/sqrt(nu*eta))
      prndtl   = 0d0           ! Prandtl number (nu/eta)

      gamma    = 5./3.         ! Polytropic constant of plasma

      !Initial condition
      equil    = ''            ! Type of equilibrium
      dlambda  = .2            ! Characteristic width of Harris sheet
      rshear   = 1.            ! Ratio of magnetic to fluid sheet thicknesses
      vparflow = 0.            ! Maximum parallel fluid flow
      vperflow = 0.            ! Maximum perpendicular fluid flow
      source   = .true.        ! Impose source to provide equilibrium

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

      !I/O parameters
      restart  = .false.       ! Restarting flag
      ilevel   = 0             ! Level of solver output information

c Read initialization parameters

      inputfile='pixie3d.in'

      open(unit=25,file=inputfile,status='old')
      read(25,datin)
      close(unit=25)

c Obtain eta, nu, dd from Prandtl, Hartmann

      if (hrtmn.gt.0d0.and.prndtl.gt.0d0) then
        nu  = sqrt(prndtl)/hrtmn
        eta = 1d0/hrtmn/sqrt(prndtl)
      endif

c Map perturbations

      pert(IRHO)= prho
      pert(IVX) = pvx
      pert(IVY) = pvy
      pert(IVZ) = pvz
      pert(IBX) = pbx
      pert(IBY) = pby
      pert(IBZ) = pbz
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
        bcond = NEU
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
        stop
      endif

c End program

      end subroutine readInput

c readGraphicsInput
c######################################################################
      subroutine readGraphicsInput(sel_diag,sel_graph,ndplot,dplot)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use iosetup

      implicit none

c Call variables

      integer(4)    :: sel_diag(9),sel_graph(9),ndplot
      real(8)       :: dplot

c Local variables

      namelist /graphdef/ sel_diag,sel_graph,ndplot,dplot

c Begin program

c Read computation initializations (external)

      call readInput

c Graphics defaults

      ndplot = 0
      dplot  = 0d0

c Read graphics initialization parameters

      open(unit=25,file=inputfile,status='old')
      read(25,graphdef)
      close(unit=25)

c End program

      end subroutine readGraphicsInput
