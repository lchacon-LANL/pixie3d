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

      real(8)    :: prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
      real(8)    :: prndtl,hrtmn

c Namelist

      namelist /datin/ neqd,nxd,nyd,nzd,coords,xlength,ylength,zlength
     .                   ,major_r,mg_ratio
     .                ,plot,ilevel,debug
     .                ,nu,eta,dd,chi,gamma,prndtl,hrtmn
     .                ,tolnewt,maxitnwt,tolgm,maxksp,iguess,maxitgm
     .                   ,global,method
     .                ,equil,dlambda,rshear,vparflow,vperflow
     .                   ,source,nh,prho,pvx,pvy,pvz,pbx,pby,pbz,ptemp
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

      xlength  = 1d0           ! Length in x-direction
      ylength  = 1d0           ! Length in y-direction
      zlength  = 1d0           ! Length in z-direction

      major_r  = 0d0           ! Major radius in toroidal coords.

      mg_ratio = 2             ! MG coarsening ratio

      !Time stepping
      dt       = 5.            ! Time step (if zero, dt is calculated in code)
      tmax     = 0d0           ! Target time, in Alfven times.
      numtime  = ntimemax      ! Number of time steps
      ndstep   = 0             ! # time steps between plots (if zero,
                               !        ndstep is calculated in code)
      dstep    = 0.            ! Time interval between plots (if zero,
                               !        dstep is calculated in code)
      timecorr = .true.        ! Time adaptiveness based on Newton convergence
      cnfactor = -.48          ! Crank-Nicolson factor
      sm_pass  = 2             ! Number of initial smoother passes for CN

      !NK parameters
      tolgm    = 5.0d-2        ! Inexact Newton parameter (GMRES conv. tolerance)
      tolnewt  = 1.0d-4        ! Newton convergence tolerance
      maxitnwt  = 0            ! Maximum number of Newton its. (if zero, maxitnwt
                               !        is determined in code)
      maxksp   = 15            ! Maximum krylov subspace dimension
      maxitgm  = maxksp        ! Maximum GMRES iterations
      method   = 0             ! Constant forcing parameter for inexact Newton
                               !        (see NewtonGmres doc)
      global   = 0             ! Do not use globalization (see NewtonGmres doc)
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
      nh       = 1             ! Harmonic number for perturbation

      !I/O parameters
      restart  = .false.       ! Restarting flag
      debug    = .false.       ! Debugging flag
      plot     = .true.        ! Plots ouput
      ilevel   = 0             ! Level of solver output information

      sel_diag = (/ 2,5,15,16,17,18,19,12,0/) 
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
                               !    16-> 'Total particles'
                               !    17-> 'Total X momentum'
                               !    18-> 'Total Y momentum'
                               !    19-> 'Total Z momentum'


      sel_graph = (/ 1,-15,-18,-9,11,14,0,0,0 /) 
                               ! Selects diagnostics for xdraw output
                               ! A negative value -i indicates to construct
                               !    a vector plot with the i and i+1 arrays.
                               ! Currently:
                               !    1 -> rho
                               !    2 -> rho*Vx
                               !    3 -> rho*Vy
                               !    4 -> rho*Vz
                               !    5 -> Bx
                               !    6 -> By
                               !    7 -> Bz
                               !    8 -> Temp
                               !    9 -> jx
                               !    10-> jy
                               !    11-> jz
                               !    12-> Nu
                               !    13-> Eta
                               !    14-> Divergence B
                               !    15-> Vx cartesian
                               !    16-> Vy cartesian
                               !    17-> Vz cartesian
                               !    18-> Bx cartesian
                               !    19-> By cartesian
                               !    20-> Bz cartesian

c Open input file

      open(unit=25,file='3dmhd.in',status='old')
      read(25,datin)
      close(unit=25)

c Obtain eta, nu, dd from Prandtl, Hartmann

      if (hrtmn.gt.0d0.and.prndtl.gt.0d0) then
        nu  = sqrt(prndtl)/hrtmn
        eta = 1d0/hrtmn/sqrt(prndtl)
      endif

c Check for autoinitializing parameters

      if (maxitnwt.eq.0) 
     .      maxitnwt = max(floor(1.5*log(tolnewt)/log(tolgm)),10)

      alpha = 1. - cnfactor

      dtbase = dt   

      sm_flag= 0
      if (cnfactor.lt.0d0) sm_flag= 1

      cnf_d = 1d0

c Map perturbations

      pert(IRHO)= prho
      pert(IVX) = pvx
      pert(IVY) = pvy
      pert(IVZ) = pvz
      pert(IBX) = pbx
      pert(IBY) = pby
      pert(IBZ) = pbz
      pert(ITMP)= ptemp

c Define topological constraints (i.e., singular points, periodic)

      select case (coords)
      case ('car','scl')
        bcond(1) = PER          !x1
        bcond(2) = PER          !x2
        bcond(3) = OTH          !y1
        bcond(4) = OTH          !y2
        bcond(5) = PER          !z1
        bcond(6) = PER          !z2
      case ('cyl')
        bcond(1) = SP           !x1
        bcond(2) = OTH          !x2
        bcond(3) = PER          !y1
        bcond(4) = PER          !y2
        bcond(5) = PER          !z1
        bcond(6) = PER          !z2
      case ('tor')
        bcond(1) = SP           !x1
        bcond(2) = OTH          !x2
        bcond(3) = PER          !y1
        bcond(4) = PER          !y2
        bcond(5) = PER          !z1
        bcond(6) = PER          !z2
      case default
        write (*,*) 'Could not define topological constraints for grid'
        write (*,*) 'Aborting...'
        stop
      end select

c End program

      end subroutine readInput
