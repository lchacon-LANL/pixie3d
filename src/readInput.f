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

      use graphics

      use icond

      implicit none

c Call variables

c Local variables

      integer(4) :: dim(1),loc(1)
      real(8)    :: prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
      real(8)    :: prndtl,hrtmn
      character*(3) :: bcs(6)

c Namelist

      namelist /datin/ neqd,nxd,nyd,nzd,coords,bcs,xmax,ymax,zmax
     .                   ,xmin,ymin,zmin,gparams,mg_ratio
     .                ,plot,ilevel,debug
     .                ,nu,eta,dd,chi,gamma,prndtl,hrtmn
     .                ,rtol,atol,maxitnwt,tolgm,maxksp,iguess,maxitgm
     .                   ,global,method,damp,dt0
     .                ,equil,dlambda,rshear,vparflow,vperflow,source
     .                ,nh1,nh2,nh3,prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
     .                   ,odd,random
     .                ,precon,maxvcyc,nsweep,precpass,iguess
     .                ,dt,cnfactor,tmax,dstep,timecorr,numtime,restart
     .                   ,ndstep,sm_pass
     .                ,sel_diag,sel_graph

c ******************************************************************

c Begin program

c Set defaults

      !General setup
      neqd     = 8             ! Number of degrees of freedom (equations)
      nxd      = 64            ! Mesh points in x-direction
      nyd      = 64            ! Mesh points in y-direction
      nzd      = 64            ! Mesh points in z-direction

      coords   = 'car'         ! Coordinate system (car,cyl,tor)

      xmax     = 1d0           ! Length in x-direction
      ymax     = 1d0           ! Length in y-direction
      zmax     = 1d0           ! Length in z-direction

      xmin     = 0d0           ! Length in x-direction
      ymin     = 0d0           ! Length in y-direction
      zmin     = 0d0           ! Length in z-direction

      gparams  = 0d0           ! Array with additional grid parameters (grid-dependent)

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
      maxitnwt  = 0            ! Maximum number of Newton its. (if zero, maxitnwt
                               !        is determined in code)
      maxksp   = 15            ! Maximum krylov subspace dimension
      maxitgm  = maxksp        ! Maximum GMRES iterations
      method   = 0             ! Constant forcing parameter for inexact Newton
                               !        (see NewtonGmres doc)
      global   = 0             ! Do not use globalization (see NewtonGmres doc)
      damp     = 1d0           ! Damping parameter in Newton (see NewtonGmres doc))
      dt0      = 1d30          ! Initial pseudo-transient time step (" " ")
      iguess   = 1             ! Whether preconditioner is used to give initial
                               !        guess to GMRES (1)

      !Preconditioner parameters
      nsweep   = 4             ! Number of SGS sweeps
      maxvcyc  = 1             ! Maximum number of V-cycles
      precpass = 4             ! Number of jacobi passes in MHD preconditioner
      precon   = '4o'          ! Type of preconditioner

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
      debug    = .false.       ! Debugging flag
      plot     = .true.        ! Plots ouput
      ilevel   = 0             ! Level of solver output information

      sel_diag = (/ 2,5,15,16,17,18,19,12,20/) 
                               ! Selects diagnostics for xdraw output
                               ! Currently:
                               !    1 -> 'ln(drho)'       
                               !    2 -> 'ln(dvx)'        
                               !    3 -> 'ln(dvy)'        
                               !    4 -> 'ln(dvz)'        
                               !    5 -> 'ln(dbx)'        
                               !    6 -> 'ln(dby)'        
                               !    7 -> 'ln(dbz)'        
                               !    8 -> 'ln(dtmp)'       
                               !    9 -> 'Magnetic energy'
                               !    10-> 'Kinetic energy' 
                               !    11-> 'Thermal energy' 
                               !    12-> 'Total energy'   
                               !    13-> 'Time step'      
                               !    14-> 'Growth rate'    
                               !    15-> 'div(B)'
                               !    16-> 'Conservation of flux'
                               !    17-> 'Total particles'
                               !    18-> 'Total X momentum'
                               !    19-> 'Total Y momentum'
                               !    20-> 'Total Z momentum'
                               !    21-> 'Flow flux at boundaries'
                               !    22-> 'Toroidal current Iz'
                               !    23-> 'Toroidal flux'

      sel_graph = (/ 1,-15,-18,-9,11,14,0,0,0 /) 
                               ! Selects diagnostics for xdraw output
                               ! A negative value -i indicates to construct
                               !    a vector plot with the i and i+1 arrays.
                               ! Currently:
                               !    1 -> rho
                               !    2 -> Px
                               !    3 -> Py
                               !    4 -> Pz
                               !    5 -> Bx
                               !    6 -> By
                               !    7 -> Bz
                               !    8 -> Temp
                               !    9 -> jx
                               !    10-> jy
                               !    11-> jz

c Open input file

      open(unit=25,file='3dmhd.in',status='old')
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
