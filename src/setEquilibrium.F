c setEquilibrium
c####################################################################
      subroutine setEquilibrium(iigx,iigy,iigz,var,bcs,label)

c--------------------------------------------------------------------
c     Set initial conditions for initial value calculation, define
c     boundary conditions, and label physical quantities.
c
c     On call:
c       * var (output): array with initial conditions for all variables
c       * bcs (output): array with boundary conditions "   "     "
c       * label (")   : array with labels              "   "     "
c
c     Boundary conditions specification in bcs is defined in routine
c     applyBoundaryConditions.f. 
c
c     Equilibrium quantities are:
c       * ieq=IRHO: density
c       * ieq=IVX:IVZ: momentum if conservative, velocity otherwise
c       * ieq=IBX:IBZ: magnetic field
c       * ieq=IAX:IAZ: vector potential
c       * ieq=ITMP: electron temperature (adiabatic=.false.) or
c                   electron pressure (adiabatic=.true.). Total
c                   pressure is found as p=n*(Ti+Te)=a_p*n*Te (or
c                   p=a_p*Te if adiabatic).
c--------------------------------------------------------------------

      use parameters

      use equilibrium

      use grid

      use timeStepping

      use constants

      use icond

      use nlfunction_setup

      use error

      use oned_int

      use local_BCS_variables

      use app_iosetup

      implicit none

c Call variables

      integer       :: iigx,iigy,iigz

      real(8)       :: var(ilom:ihip,jlom:jhip,klom:khip,neqd)

      character(20) :: label(neqd)

      integer       :: bcs(6,neqd)

c Local variables

      integer       :: i,j,k,ig,jg,kg,ieq,nmax,iglobal,jglobal
     .                ,ierr,sh,it,icomp,bcnd(6,3),brank
     .                ,igx,igy,igz,nx,ny,nz

      real(8)       :: ldaflow,x1,y1,z1,cx,cy,cz

      logical       :: covariant,to_cartsn,cartsn,a_bc

      real(8)       :: r,jac1,bnorm,aaa,bbb,ccc,qq,qqp,q0,RR,ff,nn
     .                ,aspect_ratio,mm,kk,Iz,bz0,btheta,vr,rrho,r1,bb
     .                ,r0(1),aa,dummy1,dummy2,dummy3,eta_in,err,bz_avg,dr
     .                ,bzh,bth,rh,e0_over_eta,theta_ppnch,p0,eps
     .                ,phi,rho_norm,beta_tor,Rtor,sigma,gsub(3,3)

      real(8)       :: rint(0:nxdp),dh,offset(1),bth0(0:nxdp)
     .                ,bzz0(0:nxdp),vr0(0:nxdp),rho0(0:nxdp)
     .                ,eth0(0:nxdp),ez0(0:nxdp),jth0(0:nxdp)
     .                ,jz0 (0:nxdp),bth0_old(0:nxdp),p00(0:nxdp)
     .                ,bzz0_old(0:nxdp),dummy(0:nxdp),bmag2(0:nxdp)
     .                ,q00(0:nxdp),dum (ilom:ihip,2)
     .                ,ldum(ilom:ihip,2)
     .                ,r2  (ilom:ihip,jlom:jhip)
     .                ,br02d (0:nxdp,0:nydp),bth02d(0:nxdp,0:nydp)
     .                ,bzz02d(0:nxdp,0:nydp)
     .                ,vth02d(0:nxdp,0:nydp),vr02d(0:nxdp,0:nydp)

      real(8)       :: a (ilom:ihip,jlom:jhip,klom:khip,3)
     .                ,b1(ilom:ihip,jlom:jhip)
     .                ,b2(ilom:ihip,jlom:jhip)
     .                ,b3(ilom:ihip,jlom:jhip)
     .                ,v1(ilom:ihip,jlom:jhip)
     .                ,v2(ilom:ihip,jlom:jhip)
     .                ,pp(ilom:ihip,jlom:jhip)

      character(200):: command

#if defined(vec_pot)
      real(8)       :: b(0:nxdp,0:nydp,0:nzdp,3)
#endif

c Functions

      real(8) :: qprof,pprof
      external   qprof,pprof
      
cc      !q-profile
cc      qq (r1) = 0.6125*(1 - 1.8748*r1**2 + 0.8323*r1**4)
cc      qqp(r1) = 0.6125*( -2*1.8748*r1   +4*0.8323*r1**3)
cccc      qq(r1)  = 0.3*(1 - 1.8748*r1**2 + 0.8323*r1**4)
cccc      qqp(r1) = 0.3*( -2*1.8748*r1   +4*0.8323*r1**3)
cc
cc      ff(r1) = r1**2 + qq(r1)**2

c Begin program

      ldaflow = dlambda/rshear

      var = 0d0

      E0 = 0d0 ; B0 = 0d0

#if !defined(vec_pot) 
      a_bc = .false.
#else
      a_bc = .true.
#endif

c Initialize required local grid information

      igx = iigx
      igy = iigy
      igz = iigz

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Label variables
     
      label(IRHO) = 'Rho'
#if !defined(vec_pot)
      label(IBX)  = 'B^1'
      label(IBY)  = 'B^2'
      label(IBZ)  = 'B^3'
#else
      label(IAX)  = 'A_1'
      label(IAY)  = 'A_2'
      label(IAZ)  = 'A_3'
#endif
      label(IVX)  = 'P^1'
      label(IVY)  = 'P^2'
      label(IVZ)  = 'P^3'
      if (adiabatic) then
        label(ITMP) = 'P_e'        
      else
        label(ITMP) = 'T_e'
      endif

c Define boundary conditions

      call defineBoundaryConditions(neqd,bcs)

c Set initial guess

      select case (trim(equil))

c     --------------------------------------------
c     ---------------  WAVE TESTS ----------------
c     --------------------------------------------

      case ('msw') !Magnetosonic wave

        gamma = 1d0

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

#if !defined(vec_pot)
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = 1d0
            enddo
          enddo
        enddo
#else
        B0(3) = 1d0

        var(:,:,:,IAX:IAZ) = 0d0
#endif

      case ('3wv') !3-wave MHD test

cc        gamma = 5./3.
cc        nh2   = nh1                     !This is required by setup(???)

        phi   = grid_params%params(1)   !Angle of B with x-axis, in degrees
        beta  = grid_params%params(2)   !plasma beta

        phi   = phi/180.*pi             !Phi in radians

#if !defined(PER_BC_SYNC)
        if (di > 0d0) then
          call pstop('setEquilibrium','Compile with PER_BC_SYNC=t')
        endif
#endif

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 0.5*beta/a_p

        B0(1) = cos(phi)
        B0(2) = sin(phi)
        B0(3) = 0d0

#if !defined(vec_pot)
        var(:,:,:,IBX)  = B0(1)
        var(:,:,:,IBY)  = B0(2)
        var(:,:,:,IBZ)  = B0(3)
#else
        var(:,:,:,IAX:IAZ) = 0d0
#endif

      case ('mswsn') !Magnetosonic wave in sinusoidal geometry

        gamma = 1d0

c     Check coordinates

        if (coords /= 'sin') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        var(:,:,:,IVX:IVZ) = 0d0
        var(:,:,:,IRHO) = 1d0
        var(:,:,:,ITMP) = 1d0

#if !defined(vec_pot)
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

              jac1 = gmetric%grid(igx)%jac(i,j,k)

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = jac1*1d0
            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')
cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc
cc              jac1 = gmetric%grid(igx)%jac(i,j,k)  !Not defined in the global domain
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = 0d0
cc              b(i,j,k,3) = jac1*1d0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('whslr') !Whistler wave

        gamma = 1d0

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

#if !defined(vec_pot)
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

              var(i,j,k,IBX)  = 1d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = 0d0
            enddo
          enddo
        enddo
#else
        B0(1) = 1d0

        var(:,:,:,IAX:IAZ) = 0d0
#endif

      case ('kaw') !KAW wave

        gamma = 1d0

        bb = 1d4

        di = (xmax-xmin)/2/pi/nh1*sqrt(bb/(bb-2))  !di k_par = 1

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 0.5*bb

#if !defined(vec_pot)
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

              var(i,j,k,IBX)  = 1d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = 0d0
            enddo
          enddo
        enddo
#else
        B0(1) = 1d0

        var(:,:,:,IAX:IAZ) = 0d0

#endif

c     --------------------------------------------
c     ------------  ADVECTION TESTS --------------
c     --------------------------------------------

      case ('adv1') !Advection test

        gamma = 2d0

        eps   = grid_params%params(1)
        sigma = grid_params%params(2)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Opposite Gaussian perturbations in rho, T, uniform V, no B

        var(:,:,:,IVX)  = 1d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              var(i,j,k,IRHO) = 1d0 + eps*exp(-0.5*x1**2/sigma**2)
              var(i,j,k,ITMP) = 1d0 - eps*exp(-0.5*x1**2/sigma**2)

            enddo
          enddo
        enddo

#if !defined(vec_pot)
        var(:,:,:,IBX:IBZ) = 0d0
#else
        var(:,:,:,IAX:IAZ) = 0d0
#endif

c     --------------------------------------------------
c     ---------------  KELVIN-HELMHOLTZ ----------------
c     --------------------------------------------------

      case ('khcar')

c     Define vector potential (in curvilinear coordinates) for initialization

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Kelvin-Helmholtz with constant magnetic field in cartesian coordinates

        var(:,:,:,IRHO) = 1d0
        var(:,:,:,ITMP) = 1d0

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)
              var(i,j,k,IVY)=vperflow
     .                      *tanh((x1-0.5d0*(xmax+xmin))/dlambda)
            enddo
          enddo
        enddo

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVZ)  = 0d0

#if !defined(vec_pot)
        var(:,:,:,IBX)  = 0d0
        var(:,:,:,IBY)  = 0d0
        var(:,:,:,IBZ)  = 1d0
#else
        call fillVectorPotential(a,igx,igy,igz,a_bc)

        var(:,:,:,IAX:IAZ) = a
#endif

      case ('3dkh')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a,igx,igy,igz,.false.)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (adiabatic.or.(a_p /= 2d0).or. (.not.solve_rho)) then
          messg = 'Wrong input for this equilibrium. '
     .          //'Check adiab, Ti/Te, solve_rho'
          call pstop('setEquilibrium',messg)
        endif

c     Tearing mode in cartesian coordinates

        var(:,:,:,ITMP) = 0.25

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

#if !defined(vec_pot)
             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),1)
             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),2)
             var(i,j,k,IBZ)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),3)

             var(i,j,k,IRHO)= (2d0 - var(i,j,k,IBZ)**2)
#else
             bz0 = curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                  ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),3)

             var(i,j,k,IRHO) = 2 - bz0**2
#endif

             var(i,j,k,IVX) = vperflow*tanh((y1-.5*(ymax+ymin))/ldaflow)
     .                                *sin (2*pi*z1/(zmax-zmin))
             var(i,j,k,IVY) = 0d0
             var(i,j,k,IVZ) = 0d0

            enddo
          enddo
        enddo

#if defined(vec_pot)

        call fillVectorPotential(a,igx,igy,igz,.true.)

        var(:,:,:,IAX:IAZ) = a

#endif

c     ----------------------------------------------
c     ---------------  TEARING MODE ----------------
c     ----------------------------------------------

      case ('tmcar')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a,igx,igy,igz,a_bc)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Tearing mode in cartesian coordinates

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

        bz0 = 1d0

#if !defined(vec_pot)
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),1)
             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),2)
             var(i,j,k,IBZ)=sqrt(bz0**2 - var(i,j,k,IBY)**2)
            enddo
          enddo
        enddo
#else
        var(:,:,:,IAX:IAZ) = a

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              b(i,j,k,1)  = 0d0                  
cc              b(i,j,k,2)  = tanh((x1-0.5d0*(xmax+xmin))/dlambda)
cc              b(i,j,k,3)  = sqrt(bz0**2 - b(i,j,k,2)**2)
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('tmsin')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a,igx,igy,igz,a_bc)

c     Check coordinates

        if (coords /= 'sin') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Tearing mode in sinusoidal coordinates

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

        bz0 = 1d0

#if !defined(vec_pot)
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

              !X-Y equilibrium
             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),1)
             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),2)
             var(i,j,k,IBZ)=0d0

             gsub = gmetric%grid(igx)%gsub(i,j,k,:,:)
             jac1 = gmetric%grid(igx)%jac (i,j,k)

             bnorm= vectorNorm(i,j,k,igx,igy,igz
     .                         ,var(i,j,k,IBX)
     .                         ,var(i,j,k,IBY)
     .                         ,var(i,j,k,IBZ),.false.)

             ccc = jac1*(bz0**2 - bnorm)
             bbb = gsub(3,2)*var(i,j,k,IBY) + gsub(3,1)*var(i,j,k,IBX)
             aaa = gsub(3,3)

             var(i,j,k,IBZ)  = (-bbb+sqrt(bbb**2+4*aaa*ccc))/2./aaa

             if (bbb**2+4*aaa*ccc < 0d0) then
               write (*,*) var(i,j,k,IBZ),bnorm,aaa,bbb**2+4*aaa*ccc
             endif
            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,x1,y1,z1)
cc
cc              cx = 0d0
cc              cy = tanh((x1-0.5d0*(xmax-xmin))/dlambda)
cc              cz = sqrt(bz0**2 - cy**2)  
cc
cc              call transformVectorToCurvilinear(i,j,k,igx,igy,igz
cc     .                                       ,cx,cy,cz,.false.
cc     .                                       ,b(i,j,k,1)
cc     .                                       ,b(i,j,k,2)
cc     .                                       ,b(i,j,k,3))          
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('gem')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a,igx,igy,igz,a_bc)

c     Check coordinates

        if (coords /= 'car' .and. coords /= 'scl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (adiabatic.or.(a_p /= 6d0).or. (.not.solve_rho)) then
          messg = 'Wrong input for this equilibrium. '
     .          //'Check adiab, Ti/Te, solve_rho'
          call pstop('setEquilibrium',messg)
        endif

c     GEM challenge

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,ITMP) = 1d0/12. !Required for pressure balance

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

             var(i,j,k,IRHO)=1d0/cosh((x1-0.5d0*(xmax+xmin))/dlambda)**2
     .                      +0.2
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),1)
             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),2)
            enddo
          enddo
        enddo
#else
        var(:,:,:,IAX:IAZ) = a
#endif

      case ('kai1')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a,igx,igy,igz,a_bc)

c     Check coordinates

        if (coords /= 'car' .and. coords /= 'scl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (adiabatic.or. (.not.solve_rho)) then
          messg = 'Wrong input for this equilibrium. '
     .          //'Check adiab, Ti/Te, solve_rho'
          call pstop('setEquilibrium',messg)
        endif

c     Kai TM test

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,ITMP) = 1d0/a_p !Required for pressure balance

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

             var(i,j,k,IRHO)=0.5/cosh((x1-0.5d0*(xmax+xmin))/dlambda)**2
     .                      +1d0
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),1)
             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),2)
            enddo
          enddo
        enddo
#else
        var(:,:,:,IAX:IAZ) = a
#endif

c     -------------------------------------------------
c     ---------------  ISLAND COALESCENCE -------------
c     -------------------------------------------------

      case ('ic')

        p0     = grid_params%params(1)
        eps    = grid_params%params(2)
        dlambda = (ymax-ymin)/4./pi

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a,igx,igy,igz,a_bc)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Tearing mode in cartesian coordinates

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        bz0 = 1d0

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

             var(i,j,k,ITMP)=(p0 + 0.5*(1-eps**2)
     $                       /(    cosh(x1/dlambda)
     $                        +eps*cos (y1/dlambda))**2)
     $                       /a_p/var(i,j,k,IRHO)
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),1)
             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,a(:,:,:,1),a(:,:,:,2),a(:,:,:,3),2)
             var(i,j,k,IBZ)=0d0
            enddo
          enddo
        enddo
#else
        var(:,:,:,IAX:IAZ) = a
#endif

c     ------------------------------------------------
c     ---------------  3D CYLINDRICAL ----------------
c     ------------------------------------------------

      case ('3dkai')

cc        mm = grid_params%params(1)
cc        kk = grid_params%params(2)
        RR = grid_params%params(3)

c     Check coordinates

        if (coords /= 'cyl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (abs(RR-zmax/2./pi) > 1d-5) then
          messg = 'Major radius and zmax do not agree'
          call pstop('setEquilibrium',messg)
        endif

        if ((a_p /= 2d0).or. (.not.solve_rho)) then
          messg = 'Wrong input for this equilibrium. '
     .          //'Check adiab, Ti/Te, solve_rho'
          call pstop('setEquilibrium',messg)
        endif

c     Tokamak equilibria (analytical)

        Iz = 2.22
        bb = 0.5*Iz

cc        nh2 = mm       !To set the right perturbation wavelength
cc        nh3 = kk*RR    !To set the right perturbation wavelength

        nh1 = nh2      !For perturbations to be of the right order
                        !(nh2,nh3 set by input file)

        if (nh2 == 0 .or. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*x1/(1+(x1/dlambda)**2)
              bz0    = RR

#if !defined(vec_pot)
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta
              var(i,j,k,IBZ)  = x1*bz0
#endif
              var(i,j,k,ITMP) = 1d-5
     .                     + 1/16.*(Iz*dlambda/(1+(x1/dlambda)**2))**2

            enddo
          enddo
        enddo

#if defined(vec_pot)
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              btheta = bb*x1/(1+(x1/dlambda)**2)
cc              bz0    = RR
cc
cc              b(i,j,k,1)  = 0d0                  
cc              b(i,j,k,2)  = btheta
cc              b(i,j,k,3)  = x1*bz0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('3dpar')  !EQ for parallel transport test

cc        mm = grid_params%params(1)
cc        kk = grid_params%params(2)
        RR = grid_params%params(3)
        q0 = grid_params%params(4)

c     Check coordinates

        if (coords /= 'cyl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (abs(RR-zmax/2./pi) > 1d-5) then
          messg = 'Major radius and zmax do not agree'
          call pstop('setEquilibrium',messg)
        endif

c     Simple screw pinch equilibria

cc        nh2 = mm       !To set the right perturbation wavelength
cc        nh3 = kk*RR    !To set the right perturbation wavelength

        nh1 = nh2      !For perturbations to be of the right order
                       !(nh2,nh3 set by input file)

        if (nh2 == 0 .or. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

        !Build equilibrium

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d-5

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

cc              btheta = x1  !Green's function orbit test
cc              bz0 = 1d0
cc
cc              btheta = 2*x1  !CosT orbit test
cc              bz0 = 1d0

              bz0 = 1d0
              bb  = q0*bz0*dlambda/RR
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)

              var(i,j,k,IBX) = 0d0
              var(i,j,k,IBY) = btheta
              var(i,j,k,IBZ) = x1*bz0

            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              if (coords == 'hl2') x1 = sign(sqrt(2*abs(x1)),x1)
cc
cc              bb     = (dlambda**2+aa)
cc     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              if (coords == 'hel') then
cc                b(i,j,k,1)  = 0d0
cc                b(i,j,k,2)  = btheta + kk*x1/mm*bz0
cc                b(i,j,k,3)  = x1*bz0
cc              else
cc                b(i,j,k,1)  = 0d0
cc                b(i,j,k,2)  = btheta/x1 + kk/mm*bz0
cc                b(i,j,k,3)  = bz0
cc              endif
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('3drfp')

cc        mm = grid_params%params(1)
cc        kk = grid_params%params(2)
        RR = grid_params%params(3)
        aa = grid_params%params(4)

c     Check coordinates

        if (coords /= 'cyl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (abs(RR-zmax/2./pi) > 1d-5) then
          messg = 'Major radius and zmax do not agree'
          call pstop('setEquilibrium',messg)
        endif

c     RFP equilibria (analytical)

cc        nh2 = mm       !To set the right perturbation wavelength
cc        nh3 = kk*RR    !To set the right perturbation wavelength

        nh1 = nh2      !For radial perturbations to be of the right order
                       !(nh2,nh3 set by input file)

        if (nh2 == 0 .or. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

        !Build equilibrium

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d-5

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta
              var(i,j,k,IBZ)  = x1*bz0
            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              bb     = (dlambda**2+aa)
cc     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = btheta
cc              b(i,j,k,3) = x1*bz0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

c     ----------------------------------------------------
c     ---------------  HELICAL EQUILIBRIA ----------------
c     ----------------------------------------------------

      case ('spnch')  !SCREW PINCH EQ.

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     RFP equilibria (analytical)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        aa = grid_params%params(4)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        !Build equilibrium

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d-5

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              if (coords == 'hel') then
                var(i,j,k,IBX)  = 0d0
                var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
                var(i,j,k,IBZ)  = x1*bz0
              endif
            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              if (coords == 'hl2') x1 = sign(sqrt(2*abs(x1)),x1)
cc
cc              bb     = (dlambda**2+aa)
cc     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              if (coords == 'hel') then
cc                b(i,j,k,1)  = 0d0
cc                b(i,j,k,2)  = btheta + kk*x1/mm*bz0
cc                b(i,j,k,3)  = x1*bz0
cc              else
cc                b(i,j,k,1)  = 0d0
cc                b(i,j,k,2)  = btheta/x1 + kk/mm*bz0
cc                b(i,j,k,3)  = bz0
cc              endif
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('2dtok')

c     Cylindrical Tokamak analytical equilibrium (Furth et al., 1973)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        q0 = grid_params%params(4)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        bb = dlambda/RR/q0

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d-3

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              !X-Y equilibrium
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
              var(i,j,k,IBZ)  = x1*bz0
            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = btheta + kk*x1/mm*bz0
cc              b(i,j,k,3) = x1*bz0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('2dkai')

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if ((a_p /= 2d0).or. (.not.solve_rho)) then
          messg = 'Wrong input for this equilibrium. '
     .          //'Check adiab, Ti/Te, solve_rho'
          call pstop('setEquilibrium',messg)
        endif

c     Tokamak equilibria (analytical)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)

        Iz           = 2.22
        bb           = 0.5*Iz

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*x1/(1+(x1/dlambda)**2)
              bz0    = RR

#if !defined(vec_pot)
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
              var(i,j,k,IBZ)  = x1*bz0
#endif
              var(i,j,k,ITMP) = 1d0
     .                    + 1./16.*(Iz*dlambda/(1+(x1/dlambda)**2))**2
            enddo
          enddo
        enddo

#if defined(vec_pot)
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              btheta = bb*x1/(1+(x1/dlambda)**2)
cc              bz0    = RR
cc
cc              b(i,j,k,1)  = 0d0                  
cc              b(i,j,k,2)  = btheta + kk*x1/mm*bz0
cc              b(i,j,k,3)  = x1*bz0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('rfp83')

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     RFP EQUILIBRIA (Caramana et al, PoP, 1983)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        !Find pressure and q-profiles
        bz0      = 1d0
        beta_tor = 0d0

        do i=0,nxdp
          x1  = grid_params%xg(i)
          p00(i) = 2d-5
          q00(i) = RR*qprof(x1)
        enddo

        !Initialize toroidal field
        bzz0 = bz0  

        !Solve 1D Grad-Shafranov equation
        call GS_equil_1d(q00,p00,bth0,bzz0)

        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = bth0(iglobal) + kk*x1/mm*bzz0(iglobal)
              var(i,j,k,IBZ)  = x1*bzz0(iglobal)

              var(i,j,k,ITMP) = p00(iglobal)/a_p/var(i,j,k,IRHO)

            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = bth0(i) + kk*x1/mm*bzz0(i)
cc              b(i,j,k,3) = x1*bzz0(i)
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('ntm2d')

        bootstrap= .true.

        mm       = grid_params%params(1) !Poloidal mode number
        kk       = grid_params%params(2) !Helical paramter (-nn/RR, nn-> toroidal mode number)
        RR       = grid_params%params(3) !Major radius (determines aspect ratio)
        bz0      = grid_params%params(4) !Toroidal field at magnetic axis
        beta_tor = grid_params%params(5) !Beta toroidal

        nh2 = mm  !To set the right perturbation wavelength
        nh1 = nh2 !For perturbations to be of the right order

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Equilibria (numerical)

        !Find pressure and q-profiles
        do i=0,nxdp
          x1  = grid_params%xg(i)
          p00(i) = 0.5*beta_tor*bz0**2*pprof(x1,0)/pprof(0d0,0)
          q00(i) = RR*qprof(x1)
        enddo

        !Initialize toroidal field
        bzz0 = bz0  

        !Solve 1D Grad-Shafranov equation
        call GS_equil_1d(q00,p00,bth0,bzz0)

        !Build equilibrium

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1

#if !defined(vec_pot)
              var(i,j,k,IBX) = 0d0
              var(i,j,k,IBY) = bth0(iglobal) + kk*x1/nh2*bzz0(iglobal)
              var(i,j,k,IBZ) = x1*bzz0(iglobal)
#endif

              var(i,j,k,ITMP) = 0.5*p00(iglobal)/var(i,j,k,IRHO)

            enddo
          enddo
        enddo

#if defined(vec_pot)
cc        call pstop('setEquilibrium','Equilibrium not working')

        do k = 0,nzdp
          do j = 0,nydp
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              b(i,j,k,1) = 0d0
              b(i,j,k,2) = bth0(i) + kk*x1/nh2*bzz0(i)
              b(i,j,k,3) = x1*bzz0(i)
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
#endif

cc      case ('3dtok')
cc
ccc     Cylindrical Tokamak analytical equilibrium (Furth et al., 1973)
cc
cc        RR = grid_params%params(1)
cccc        mm = grid_params%params(4)
cccc        nn = grid_params%params(5)
cc        q0 = grid_params%params(6)
cc
cc        kk = 0d0   !To avoid affecting q-profile
cc
cc        bb = dlambda/RR/q0
cc
ccc     Check coordinates
cc
cc        if (coords /= 'tor') then
cc          messg = 'Need "tor" coordinates for equilibrium '//equil
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
ccc     Failsafes
cc
cccc        nh2 =  mm       !To set the right perturbation wavelength
cccc        nh3 = -nn       !To set the right perturbation wavelength
cc
cc        nh1 = nh2       !For perturbations to be of the right order
cc                        !(nh2, nh3 set in input file)
cc
cc        if (nh2 == 0 .or. nh3 == 0) then
cc          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
cc     .          //' equilibrium not set'
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
ccc     Build equilibrium
cc
cc        var(:,:,:,IVX)  = 0d0
cc        var(:,:,:,IVY)  = 0d0
cc        var(:,:,:,IVZ)  = 0d0
cc
cc        var(:,:,:,IRHO) = 1d0
cc
cc        var(:,:,:,ITMP) = 1d-3
cc
cc#if !defined(vec_pot)
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc
cc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)
cc
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,x1,y1,z1)
cc
cc              !X-Y equilibrium
cc              var(i,j,k,IBX)  = 0d0
cc              var(i,j,k,IBY)  = btheta
cc              var(i,j,k,IBZ)  = gmetric%grid(igx)%gsup(i,j,k,3,3)*bz0
cc            enddo
cc          enddo
cc        enddo
cc#else
cc        call pstop('setEquilibrium','Equilibrium not working')
cc
cccc        do k = klom,khip
cccc          do j = jlom,jhip
cccc            do i = ilom,ihip
cccc
cccc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cccc     .                                      ,x1,y1,z1)
cccc
cccc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cccc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cccc
cccc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cccc     .                                    ,x1,y1,z1)
cccc
cccc              !X-Y equilibrium
cccc              b(i,j,k,1) = 0d0
cccc              b(i,j,k,2) = btheta
cccc              b(i,j,k,3) = gmetric%grid(igx)%gsup(i,j,k,3,3)*bz0
cccc            enddo
cccc          enddo
cccc        enddo
cccc
cccc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cccc     .                          ,b
cccc     .                          ,var(:,:,:,IAX:IAZ))
cc#endif

#if defined(vmec)

c     --------------------------------------------------
c     ---------------  VMEC EQUILIBRIUM ----------------
c     --------------------------------------------------

      case ('vmec')    

        nh1 = nh2       !For perturbations to be of the right order in radius
                        !(nh2,nh3 set by input file)

        if (nh2 == 0 .or. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

c     Check coordinates

        if (coords /= 'vmc') then
          messg = 'Need "vmc" coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Read VMEC equilibrium

        !Build equilibrium

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

#if !defined(vec_pot)
        call vmec_equ(igx,nx,ny,nz,var(:,:,:,IBX),var(:,:,:,IBY)
     .                            ,var(:,:,:,IBZ),var(:,:,:,ITMP))

c     Clean up flux surface quantities

        k = klo   !Fix k slice

        !Find R^2
        do i = ilom,ihip
          do j = jlom,jhip
            call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)
            r2(i,j) = x1**2+y1**2
          enddo
        enddo

        !Average in j as a function of i (for a single k slice per processor)
        dum = 0
        do j = jlo,jhi
          dum(:,2) = dum(:,2) + var(:,j,k,IBZ)*r2(:,j)
     .                         /gmetric%grid(igx)%jac(:,j,k)
          dum(:,1) = dum(:,1) + var(:,j,k,IBY)
        enddo

#if defined(petsc)
        ldum = dum
        call MPI_Allreduce(ldum,dum,2*(nx+2),MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_POL,mpierr)
#endif
        dum = dum/nyd           !Average in j

        !Fill flux (j,k) surfaces
        do k = klom,khip
           do j = jlom,jhip
              var(:,j,k,IBZ) = dum(:,2)*gmetric%grid(igx)%jac(:,j,k)
     $                                 /r2(:,j)
              var(:,j,k,IBY) = dum(:,1)
           enddo
        enddo

c$$$        !Transform to covariant
c$$$cc        call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
c$$$cc     .                    ,var(:,:,:,IBX),var(:,:,:,IBY),var(:,:,:,IBZ)
c$$$cc     .                    ,'cnv','cov')
c$$$
c$$$        !Perform flux average of B_3, B^2 (both flux functions in Tokamaks)
c$$$        !(instead of B_3, we can also do I(phi)=R^2 B^3/J)
c$$$        do i = ilom,ihip
c$$$          dum = 0
c$$$          do j = jlo,jhi
c$$$
c$$$cc            dum = dum + var(i,j,:,IBZ)  !For covariant B_3
c$$$
c$$$            call getCartesianCoordinates(i,j,1,igx,igy,igz,ig,jg,kg
c$$$     .                                    ,x1,y1,z1)
c$$$            dum(:,2) = dum(:,2) + var(i,j,:,IBZ)*(x1**2+y1**2)
c$$$     .                           /gmetric%grid(igx)%jac(i,j,:)
c$$$            dum(:,1) = dum(:,1) + var(i,j,:,IBY)
c$$$          enddo
c$$$
c$$$#if defined(petsc)
c$$$          ldum = dum
c$$$          call MPI_Allreduce(ldum,dum,2*(nz+2),MPI_DOUBLE_PRECISION
c$$$     .                      ,MPI_SUM,MPI_COMM_POL,mpierr)
c$$$#endif
c$$$
c$$$          dum = dum/nyd         !Average
c$$$
c$$$          do j = jlo,jhi
c$$$            call getCartesianCoordinates(i,j,1,igx,igy,igz,ig,jg,kg
c$$$     .                                    ,x1,y1,z1)
c$$$cc            var(i,j,:,IBZ) = dum  !For covariant B_3
c$$$            var(i,j,:,IBZ) = dum(:,2)*gmetric%grid(igx)%jac(i,j,:)
c$$$     .                              /(x1**2+y1**2)
c$$$            var(i,j,:,IBY) = dum(:,1)
c$$$          enddo
c$$$        enddo
c$$$
c$$$cc        call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
c$$$cc     .                    ,var(:,:,:,IBX),var(:,:,:,IBY),var(:,:,:,IBZ)
c$$$cc     .                    ,'cov','cnv')

#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        call vmec_equ(igx,nx,ny,nz,b(:,:,:,1),b(:,:,:,2)
cc     .                            ,b(:,:,:,3),var(:,:,:,ITMP))
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

        var(:,:,:,IRHO) = abs(var(:,:,:,ITMP))**(1d0/gamma)   !Find density from pressure

cc        if (solve_rho) then
          rho_norm = var(1,1,1,IRHO)
#if defined(petsc)
          !Broadcast value of rho at SP
          if (      grid_params%ilog == 1
     .        .and. grid_params%jlog == 1
     .        .and. grid_params%klog == 1) brank = my_rank
          call MPI_Bcast(rho_norm,1,MPI_DOUBLE_PRECISION,brank
     .                  ,MPI_COMM_WORLD,mpierr)
#endif
          var(:,:,:,IRHO) = var(:,:,:,IRHO)/rho_norm !Normalize density to 1 at r=0
cc        else
cc          var(:,:,:,IRHO) = 1d0
cc        endif

        if (adiabatic) then
          var(:,:,:,ITMP) = var(:,:,:,ITMP)/a_p !Transform to electron pressure
        else
          var(:,:,:,ITMP) = var(:,:,:,ITMP)/var(:,:,:,IRHO)/a_p !Transform to temperature
        endif
#endif

c     --------------------------------------------------
c     ---------------  OHMIC EQUILIBRIA ----------------
c     --------------------------------------------------

      case ('ppnch','ppnsl','ppnst')

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        theta_ppnch = grid_params%params(4)
        e0_over_eta = grid_params%params(5)
        beta_tor    = grid_params%params(6)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Equilibria (numerical)

        if (e0_over_eta == 0d0) then     !External file equilibrium

          if (my_rank == 0)
     .      write (*,*) 'Reading numerical equilibrium ',equil,'...'

          !Read file
          open (unit=1000,file=trim(equ_file),status='old')

          read(1000,'(a)')
          do i=0,nxdp
            read(1000,'(9f12.8)') dummy1,bth0(i),bzz0(i),jth0(i),jz0(i)
     .                  ,vr0(i),ez0(i),eth0(i),rho0(i)
cc            write (*,*)
cc            write (*,*) dummy1,bth0(i),bzz0(i),jth0(i),jz0(i)
cc     .                  ,vr0(i),ez0(i),eth0(i),rho0(i)
          enddo

          read(1000,'(a)')
          read(1000,'(a)')
          read(1000,'(f7.1)') eta_in

          close(1000)

          if (my_rank ==0) write (*,*) 'Done!'

          !Find electric field at face r=1d0
          E0(2) = 0.5*(eth0(nxd)+eth0(nxdp))
          E0(3) = 0.5*(ez0 (nxd)+ez0 (nxdp))
          E0 = E0*eta/eta_in

        else   !Calculate equilibrium

          if (my_rank == 0)
     .       write (*,*) 'Calculating equilibrium ',equil,'...'

          bz0 = 1d0

          !Find pressure and q-profiles
          do i=0,nxdp
            x1  = grid_params%xg(i)
            p00(i) = 0.5*beta_tor*bz0**2*pprof(x1,0)
          enddo

	  !Initialize toroidal field
          bzz0 = bz0

          call Ohm_equil_1d(e0_over_eta,theta_ppnch
     .                     ,p00,bth0,bzz0,vr0)

	  !Find electric field
          E0(1:2) = 0d0
          E0(3) = e0_over_eta*eta
          eta_in = eta

          !Adjust radial flow to actual resistivity
          vr0 = eta*vr0

	  !Find density
          rho0 = 1d0

        endif

#if defined(RFX)
        if (my_rank == 0) then
          write (*,*) 'r,bth0,bzz0,q0,vr0'
          do i=0,nxdp
            x1 = grid_params%xg(i) !Global coordinate
            write(*,10) x1,bth0(i),bzz0(i)
     .             ,(x1*bzz0(i))/(RR*bth0(i)),vr0(i)
          enddo
        endif
#else
cc        write (*,*) 'r,bth0,bzz0,vr0'
cc        do i=0,nxdp
cc          x1 = grid_params%xg(i) !Global coordinate
cc          write(*,*) x1,bth0(i),bzz0(i),vr0(i)
cc        enddo
cc        stop
#endif
        !Convert current to covariant at r=1d0
        E0(3) = -kk/mm*E0(2) + E0(3)

        !Build equilibrium
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1

              !X-Y equilibrium
#if !defined(vec_pot)
              var(i,j,k,IBX) = 0d0
              var(i,j,k,IBY) = bth0(iglobal) + kk*x1/mm*bzz0(iglobal)
              var(i,j,k,IBZ) = x1*bzz0(iglobal)
#endif
              var(i,j,k,IVX) = x1*vr0(iglobal)*eta/eta_in
              var(i,j,k,IVY) = 0d0
              var(i,j,k,IVZ) = 0d0

              if (solve_rho) then
                var(i,j,k,IRHO) = rho0(iglobal)
              else
                var(i,j,k,IRHO) = 1d0
              endif

              var(i,j,k,ITMP) = p00(iglobal)/a_p/var(i,j,k,IRHO)

            enddo
          enddo
        enddo

#if defined(vec_pot)
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = bth0(i) + kk*x1/mm*bzz0(i)
cc              b(i,j,k,3) = x1*bzz0(i)
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

        !Shift vr for EQU BC
        if (ihip + grid_params%ilo(igx)-1 == nxdp) then
          var(ihi+1,:,:,IVX) = 0.25*(var(ihi  ,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi,:,:)
     .                              +var(ihi+1,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi+1,:,:))
     .                         *(gmetric%grid(igx)%jac(ihi  ,:,:)
     .                          +gmetric%grid(igx)%jac(ihi+1,:,:))
        endif

      case ('2dhel')

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        theta_ppnch = grid_params%params(4)
        e0_over_eta = grid_params%params(5)
        bz0         = grid_params%params(6)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Equilibria (numerical)

        if (my_rank == 0)
     .     write (*,*) 'Calculating equilibrium ',equil,'...'

	!Initialize toroidal field
        bzz02d = 1d0

        call Ohm_equil_2d(-kk,nh2,e0_over_eta,theta_ppnch,bth02d,bzz02d)

	!Find electric field
        E0(1:2) = 0d0
        E0(3) = e0_over_eta*eta
        eta_in = eta

        !Find radial flow
        vr02d = -(E0(3)*bth02d-E0(2)*bzz02d)/(bzz02d**2+bth02d**2)

	!Find density
        rho0 = 1d0

        !Convert current to covariant at r=1d0
        E0(3) = -kk/mm*E0(2) + E0(3)

        !Build equilibrium
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1
              jglobal = j + grid_params%jlo(igy) - 1

              !X-Y equilibrium
#if !defined(vec_pot)
              var(i,j,k,IBX) = 0d0
              var(i,j,k,IBY) = bth02d(iglobal,jglobal)
     .                       + kk*x1/mm*bzz02d(iglobal,jglobal)
              var(i,j,k,IBZ) = x1*bzz02d(iglobal,jglobal)
#endif
              var(i,j,k,IVX) = x1*vr02d(iglobal,jglobal)*eta/eta_in
              var(i,j,k,IVY) = 0d0
              var(i,j,k,IVZ) = 0d0

              if (solve_rho) then
                var(i,j,k,IRHO) = rho0(iglobal)
              else
                var(i,j,k,IRHO) = 1d0
              endif

              var(i,j,k,ITMP) = 1d-5

            enddo
          enddo
        enddo

#if defined(vec_pot)
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = bth0(i) + kk*x1/mm*bzz0(i)
cc              b(i,j,k,3) = x1*bzz0(i)
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

        !Shift vr for EQU BC
        if (ihip + grid_params%ilo(igx)-1 == nxdp) then
          var(ihi+1,:,:,IVX) = 0.25*(var(ihi  ,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi,:,:)
     .                              +var(ihi+1,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi+1,:,:))
     .                         *(gmetric%grid(igx)%jac(ihi  ,:,:)
     .                          +gmetric%grid(igx)%jac(ihi+1,:,:))
        endif

      case ('ppn3d','p3nsl')

cc        mm = grid_params%params(1)
cc        kk = grid_params%params(2)
        RR = grid_params%params(3)
        theta_ppnch = grid_params%params(4)
        e0_over_eta = grid_params%params(5)
        bz0         = grid_params%params(6)
        beta_tor    = grid_params%params(7)

cc        nh2 = mm       !To set the right perturbation wavelength
cc        nh3 = kk*RR    !To set the right perturbation wavelength

        nh1 = nh2      !For perturbations to be of the right order
                       !(nh2,nh3 set by input file)

        if (nh2 == 0 .or. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

c     Check coordinates

        if (coords /= 'cyl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (abs(RR-zmax/2./pi) > 1d-5) then
          messg = 'Major radius and zmax do not agree'
          call pstop('setEquilibrium',messg)
        endif

c     Equilibria (numerical)

        if (e0_over_eta == 0d0) then     !External file equilibrium

          if (my_rank == 0)
     .      write (*,*) 'Reading numerical equilibrium ',equil,'...'

          !Read file
          open (unit=1000,file=trim(equ_file),status='old')

          read(1000,'(a)')
          do i=0,nxdp
            read(1000,'(9f12.8)') dummy1,bth0(i),bzz0(i),jth0(i),jz0(i)
     .                  ,vr0(i),ez0(i),eth0(i),rho0(i)
          enddo

          read(1000,'(a)')
          read(1000,'(a)')
          read(1000,'(f7.1)') eta_in

          close(1000)

          if (my_rank ==0) write (*,*) 'Done!'

          !Find electric field at face r=1d0
          E0(2) = 0.5*(eth0(nxd)+eth0(nxdp))
          E0(3) = 0.5*(ez0 (nxd)+ez0 (nxdp))
          E0 = E0*eta/eta_in

        else   !Calculate equilibrium

          if (my_rank == 0)
     .       write (*,*) 'Calculating equilibrium ',equil,'...'

          bz0 = 1d0

          !Find pressure and q-profiles
          do i=0,nxdp
            x1  = grid_params%xg(i)
            p00(i) = 0.5*beta_tor*bz0**2*pprof(x1,0)
          enddo

	  !Initialize toroidal field
          bzz0 = bz0

          call Ohm_equil_1d(e0_over_eta,theta_ppnch
     .                     ,p00,bth0,bzz0,vr0)

          !Find electric field
          E0(1:2) = 0d0
          E0(3) = e0_over_eta*eta
          eta_in = eta

          !Adjust radial flow to actual resistivity
          vr0 = eta*vr0

          !Find density
          rho0 = 1d0

        endif

#if defined(RFX)
        if (my_rank == 0) then
          write (*,*) 'r,bth0,bzz0,q0,vr0'
          do i=0,nxdp
            x1 = grid_params%xg(i) !Global coordinate
            write(*,10) x1,bth0(i),bzz0(i)
     .             ,(x1*bzz0(i))/(RR*bth0(i)),vr0(i)
          enddo
        endif
#else
cc        write (*,*) 'r,bth0,bzz0,jth0,jz0,vr0'
cc        do i=0,nxdp
cc          call getCurvilinearCoordinates(i,1,1,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,y1,z1)
cc          write(*,*) x1,bth0(i),bzz0(i),vr0(i)
cc        enddo
cc        stop
#endif

        !Build equilibrium
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1

#if !defined(vec_pot)
              var(i,j,k,IBX) = 0d0
              var(i,j,k,IBY) = bth0(iglobal)
              var(i,j,k,IBZ) = x1*bzz0(iglobal)
#endif
              var(i,j,k,IVX) = x1*vr0(iglobal)*eta/eta_in
              var(i,j,k,IVY) = 0d0
              var(i,j,k,IVZ) = 0d0

              if (solve_rho) then
                var(i,j,k,IRHO) = rho0(iglobal)
              else
                var(i,j,k,IRHO) = 1d0
              endif

              var(i,j,k,ITMP) = 1d-5

            enddo
          enddo
        enddo

#if defined(vec_pot)
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = bth0(i)
cc              b(i,j,k,3) = x1*bzz0(i)
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

        !Shift vr for EQU BC
        if (ihip + grid_params%ilo(igx)-1 == nxdp) then
          var(ihi+1,:,:,IVX) = 0.25*(var(ihi  ,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi,:,:)
     .                              +var(ihi+1,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi+1,:,:))
     .                         *(gmetric%grid(igx)%jac(ihi  ,:,:)
     .                          +gmetric%grid(igx)%jac(ihi+1,:,:))
        endif

      case ('ohtor')

        !Parameters 2 and 3 define elliptical shape; used in grid_anal_map
        RR = grid_params%params(1)              !Major radius
cc        mm = grid_params%params(4)              !m mode
cc        nn = grid_params%params(5)              !n mode
        e0_over_eta = grid_params%params(6)*RR  !Covariant component is constant, needs RR factor
        beta_tor    = grid_params%params(7)     !Toroidal beta

        kk = 0d0   !To avoid affecting q-profile

cc        nh2 =  mm       !To set the right perturbation wavelength
cc        nh3 = -nn       !To set the right perturbation wavelength

        nh1 = nh2       !For perturbations to be of the right order
                        !(nh2, nh3 set in input file)

        if (nh2 == 0 .or. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

c     Check coordinates

        if (coords /= 'tor') then
          messg = 'Need "tor" coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Ohmic GS equilibrium

        !Build equilibrium
        call Ohm_GS_equil_2d(RR,e0_over_eta,beta_tor,nx,ny,nz,b1,b2,b3
     .                      ,v1,v2,pp)

	!Find electric field
        E0(1:2) = 0d0
        E0(3)   = e0_over_eta*eta
        eta_in  = eta

        !Adjust pinch flow to current resistivity
        v1 = v1*eta
        v2 = v2*eta

        !Fill arrays
#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              var(i,j,k,IRHO) = 1d0

              !R-Z equilibrium
              var(i,j,k,IBX)  = b1(i,j)
              var(i,j,k,IBY)  = b2(i,j)
              var(i,j,k,IBZ)  = b3(i,j)

              var(i,j,k,IVX)  = v1(i,j)
              var(i,j,k,IVY)  = v2(i,j)
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,ITMP) = pp(i,j)/a_p/var(i,j,k,IRHO)
            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc
cc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)
cc
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,x1,y1,z1)
cc
cc              !X-Y equilibrium
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = btheta
cc              b(i,j,k,3) = gmetric%grid(igx)%gsup(i,j,k,3,3)*bz0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

        !Shift vr for EQU BC
        if (ihip + grid_params%ilo(igx)-1 == nxdp) then
          var(ihi+1,:,:,IVX) = 0.25*(var(ihi  ,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi,:,:)
     .                              +var(ihi+1,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi+1,:,:))
     .                         *(gmetric%grid(igx)%jac(ihi  ,:,:)
     .                          +gmetric%grid(igx)%jac(ihi+1,:,:))
          var(ihi+1,:,:,IVY) = 0.5*(var(ihi  ,:,:,IVY)
     .                             +var(ihi+1,:,:,IVY))
        endif

      case default

        messg = 'Equilibrium '//trim(equil)//' undefined'
        call pstop('setEquilibrium',messg)

      end select

c Find auxiliary equilibrium quantities

cc#if defined(vec_pot)
cc      call findEquilJ(igx,igy,igz,nx,ny,nz,var(:,:,:,IAX:IAZ))
cc#else
cc      call findEquilJ(igx,igy,igz,nx,ny,nz,var(:,:,:,IBX:IBZ))
cc#endif

c Find momentum components

      if (.not.nc_eom_v) then
        var(:,:,:,IVX) = var(:,:,:,IRHO)*var(:,:,:,IVX)
        var(:,:,:,IVY) = var(:,:,:,IRHO)*var(:,:,:,IVY)
        var(:,:,:,IVZ) = var(:,:,:,IRHO)*var(:,:,:,IVZ)
      endif

c Transfer electric field info to postprocessor via input file

      if (sum(E0) /= 0d0 .and. my_rank == 0) then

        command = 'grep E0 ' // trim(inputfile) // ' > /dev/null'
        ierr = system(trim(command))

        write (command,*) '   E0 =',real(E0(1),4)
     .                         ,',',real(E0(2),4)
     .                         ,',',real(E0(3),4)


        if (ierr /= 0) then  !Replace E0
          command=
     .       'sed -i --follow-symlinks -e"/graphdef/{" -e' // achar(39)
#if defined(pgf90) || defined(xlf)
     .       // 'a\\' // achar(39) //' -e"'
#else
     .       // 'a\' // achar(39) //' -e"'
#endif
     .       // trim(command)   // '" -e"}" '
     .       // trim(inputfile)
        else                 !Add E0
          command='sed -i --follow-symlinks -e"/E0/{" -e' // achar(39)
#if defined(pgf90) || defined(xlf)
     .       // 'c\\' // achar(39) //' -e"'
#else
     .       // 'c\' // achar(39) //' -e"'
#endif
     .       // trim(command)   // '" -e"}" '
     .       // trim(inputfile)
        endif

        ierr = system(trim(command))
      endif

c Transfer magnetic field info to postprocessor via input file

      if (sum(B0) /= 0d0 .and. my_rank == 0) then

        command = 'grep B0 ' // trim(inputfile) // ' > /dev/null'
        ierr = system(trim(command))

        write (command,*) '   B0 =',real(B0(1),4)
     .                         ,',',real(B0(2),4)
     .                         ,',',real(B0(3),4)


        if (ierr /= 0) then  !Replace B0
          command=
     .       'sed -i --follow-symlinks -e"/graphdef/{" -e' // achar(39)
#if defined(pgf90) || defined(xlf)
     .       // 'a\\' // achar(39) //' -e"'
#else
     .       // 'a\' // achar(39) //' -e"'
#endif
     .       // trim(command)   // '" -e"}" '
     .       // trim(inputfile)
        else                 !Add B0
          command='sed -i --follow-symlinks -e"/B0/{" -e' // achar(39)
#if defined(pgf90) || defined(xlf)
     .       // 'c\\' // achar(39) //' -e"'
#else
     .       // 'c\' // achar(39) //' -e"'
#endif
     .       // trim(command)   // '" -e"}" '
     .       // trim(inputfile)
        endif

        ierr = system(trim(command))
      endif

c End program

10    format (4f7.3,e12.3)

      end subroutine setEquilibrium

c     GS_equil_1d
c     ###############################################################
      subroutine GS_equil_1d(Rxq0,p0,bth0,bzz0)

c     ---------------------------------------------------------------
c     Solves 1D Grad-Shafranov equation to find cylindrical components
c     of magnetic field B_theta, B_z from specified q-profile and
c     pressure profile. On input:
c       * RR: major radius (aspect ratio -> 1/RR)
c       * Rxq0,p0: normalized q-profile (RR x q) and pressure profile
c       * bth0,bzz0: output magnetic field components (bzz0 contains
c                    initial guess, and bzz0(0) contains value at axis)
c     ---------------------------------------------------------------

        use parameters

        use equilibrium

        use grid

        use app_iosetup

        implicit none

c     Call variables

        real(8)    :: bth0(0:nxdp),bzz0(0:nxdp)
     .               ,p0(0:nxdp),Rxq0(0:nxdp)

c     Local variables

        integer    :: i,j,k,ig,jg,kg,it
        real(8)    :: bz_avg,x1,bzh,bth,gs_err,rh,dh

        real(8)    :: bz0,rint(0:nxdp),dummy(0:nxdp)
     .               ,p0h(0:nxdp),q0h(0:nxdp)
     .               ,xh(0:nxd),dr(0:nxd)

c     Begin program

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '************************************'
          write (*,*) '1D Grad-Shafranov equilibrium solver'
          write (*,*) '************************************'
          write (*,*)
        endif

c       Initialize iteration

        bz0 = bzz0(0)

        !Find half-mesh positions
        do i=0,nxd
          xh(i) = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))
          dr(i) =      grid_params%xg(i+1)-grid_params%xg(i)
        enddo

        !Average pressure and q-profiles to half mesh
        do i=0,nxd
          p0h(i) = 0.5*(p0(i)  +p0  (i+1)) !Average pressure  to half mesh
          q0h(i) = 0.5*(Rxq0(i)+Rxq0(i+1)) !Average q-profile to half mesh
        enddo

c       Iteration

        do it=1,100

cc          !Perform integral of B_z/qprof on half-mesh
cc          dummy(0) = 0d0
cc          bzz0 (0) = bz0
cc          do i=1,nxd
cc            rh  = grid_params%xg(i)
cc
cc            bzh = 0.5*(bzz0(i)+bzz0(i-1))
cc
cc            dummy(i) = dummy(i-1) + 2*(bzh/Rxq0(i))**2*rh*dr(i)
cc          enddo

          !Perform integral of B_th^2/r on half-mesh
          dummy(0) = 0d0
          do i=1,nxd
            rh  = grid_params%xg(i)
            bth = 0.5*(bth0(i)+bth0(i-1))
            dummy(i) = dummy(i-1) + 2*bth**2/rh*dr(i)
          enddo

          !Safeguard for solvability
          where (dummy>bz0**2+2*(p0h(0)-p0h))
     .           dummy=bz0**2+2*(p0h(0)-p0h)

          !Calculate Bz, Btheta on half-mesh
          do i=0,nxd
            x1  = xh(i)
            bth0(i) = x1*bzz0(i)/q0h(i)
            bzz0(i) = sign(1d0,q0h(i))
     .               *sqrt((bz0**2+2*(p0h(0)-p0h(i))-dummy(i))
     .                    /((x1/q0h(i))**2+1))
          enddo

          !Check convergence (force balance)
          gs_err = 0d0
          do i=1,nxd
            rh  = grid_params%xg(i)

            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))
            gs_err = gs_err
     .                +(bzh*(bzz0(i)-bzz0(i-1))/dr(i)
     .                     +(p0h (i)-p0h (i-1))/dr(i)
cc     $                 +bth/rh*(bth0(i)*i-bth0(i-1)*(i-1)))**2
     $                 +0.5*(bth0(i)**2-bth0(i-1)**2)/dr(i)
     .                 +bth**2/rh)**2
          enddo
          gs_err = sqrt(gs_err/nxd)

          if (ilevel > 0 .and. my_rank == 0) then
            write (*,*) '1D Grad-Shafranov equilibrium iter =',it
     .                 ,' ; Force balance rror=',gs_err
          endif

          dh = xmax/nxd
          if (gs_err < 0.01*dh**2) exit

        enddo

        if (ilevel == 0 .and. my_rank == 0) then
           write (*,'(a,i3,a,1pe10.2)')
     .          ' 1D equilibrium converged in '
     .          ,it,' iterations with force-balance error =',gs_err
        endif

cc        !Check force balance
cc        do i=1,nxd
cc          rh  = grid_params%xg(i)
cc
cc          bzh = 0.5*(bzz0(i)+bzz0(i-1))
cc          bth = 0.5*(bth0(i)+bth0(i-1))
cc          dummy(i) =  bzh*(bzz0(i)-bzz0(i-1))/dr(i)
cc     .                   +(p0h (i)-p0h (i-1))/dr(i)
cc     $              +bth/rh*(bth0(i)*i-bth0(i-1)*(i-1))
cc        enddo
cc
cc        if (my_rank == 0) then
cc          write(*,*)
cc          write(*,'(a,1pe14.7)') ' Force-balance error ='
cc     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
cccc          write (*,*) dummy(1:nxd)
cc        endif

        !Calculate equilibrium properties
        bz_avg = 0d0
        do i=1,nxd
          x1  = xh(i)
          bz_avg = bz_avg + 2*bzz0(i)*x1*dr(i)
        enddo

        if (my_rank == 0) then
          write (*,*) 
          write (*,999) ' 1D Grad-Shafranov equilibrium properties:'
          write (*,999) '   Theta  =',bth0(nxd)/bz_avg
          write (*,999) '   Bz(r=0)=',bzz0(0)
          write (*,999) '   I0     =',2*pi*xmax*bth0(nxd)
          write (*,999) '   beta_t =',2*p0h(0)/(bzz0(0)**2)
          write (*,999) '   beta_p =',2*p0h(0)/(bth0(nxd)**2)
 999      format (a,f10.5)
        endif

        !Save some equilibrium properties (to calculate bootstrap current)
        eq_params(1) = 1d0/grid_params%params(3)       !Aspect ratio
        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
        eq_params(3) = q0h(nxd)/grid_params%params(3)  !q(r=a)

        !Transfer to integer mesh
        call IntDriver1d(nxdp,xh(0:nxd),bth0(0:nxd)
     .                  ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp)
     .                  ,2,0)
        bth0 = dummy

        call IntDriver1d(nxdp,xh(0:nxd),bzz0(0:nxd)
     .                  ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp)
     .                  ,2,0)
        bzz0 = dummy

        !Extrapolate linearly in radius @ r=0
        bth0(0) = -bth0(1)
        bzz0(0) = bzz0(1)

        !Dump equilibrium for eigen code

        if (my_rank == 0 .and. equil == 'ntm2d') then
          open(unit=123,file='ntm2d.txt',status='unknown')
          write (123,*) nxdp+1
          write (123,*) grid_params%xg(0:nxdp)
          write (123,*) bth0
          write (123,*) bzz0
          write (123,*) 0.5*p0
          write (123,*) grid_params%params(3),eq_params(2),eq_params(3)
          close(123)
cc          stop
        endif

      end subroutine GS_equil_1d

c     Ohm_equil_1d
c     ###############################################################
      subroutine Ohm_equil_1d(Eoeta,theta,p0,bth0,bzz0,vr0)

c     ---------------------------------------------------------------
c     Solves 1D Ohmic eq. equations to find cylindrical components
c     of magnetic field B_theta, B_z and radial velocity vr from 
c     imposed electric field, theta parameter, and pressure profile.
c     On input:
c       * Eoeta: electric field over eta
c       * theta: RFP theta parameter
c       * b_tor: toroidal beta
c       * p0: input presure profile
c       * bth0,bzz0: output magnetic field components (bzz0 contains
c                    initial guess, and bzz0(0) contains value at axis)
c       * vr0: output radial pinch flow
c     ---------------------------------------------------------------

        use parameters

        use equilibrium

        use transport_params

        use grid

        use app_iosetup

        implicit none

c     Call variables

        real(8)    :: bth0(0:nxdp),bzz0(0:nxdp),vr0(0:nxdp)
     .               ,p0(0:nxdp),Eoeta,theta,b_tor

c     Local variables

        integer    :: i,j,k,ig,jg,kg,it
        real(8)    :: bz_avg,x1,bzh,bth,err,rh,dh,etah,p00

        real(8)    :: bz0,rint(0:nxdp),dummy(0:nxdp),bth0_old(0:nxdp)
     .               ,bzz0_old(0:nxdp),xh(0:nxd),dr(0:nxd)

        real(8)    :: pprof
        external      pprof

c     Begin program

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '************************************'
          write (*,*) '1D Ohmic equilibrium solver'
          write (*,*) '************************************'
          write (*,*)
        endif

c       Initialize iteration

        bz0 = bzz0(0)

cc        p00 = 0.5*b_tor*bz0**2

        !Find half-mesh positions
        do i=0,nxd
          xh(i) = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))
          dr(i) =      grid_params%xg(i+1)-grid_params%xg(i)
        enddo

        !Initial guess for B_theta
        do i=0,nxd
          bth0(i) = 0.5*xh(i)*Eoeta
        enddo

c       Iteration

        do it=1,100

          bth0_old = bth0
          bzz0_old = bzz0

          !Perform integral of B_theta
          dummy(0) = 0d0
          bth0 (0) = 0d0
          do i=1,nxd
            x1 = xh(i)
            rh = grid_params%xg(i)

            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))
cc            etah = (1.+aa_eta*rh**bb_eta)**cc_eta
            etah =res(i,1,1,nxd,nyd,nzd,1,1,1,rrr=rh)/eta

            dummy(i)=dummy(i-1)
cc     .               + (Eoeta/etah*bzh**2-bth*p00*pprof(rh,1))
     .              +(Eoeta/etah*bzh**2-bth*(p0(i+1)-p0(i-1))/2./dr(i))
     .                /(bzh**2+bth**2)*rh*dr(i)

            bth0(i) = dummy(i)/x1
          enddo

          !Perform integral of Bz
          dummy(0) = 0d0
          do i=1,nxd
            x1 = xh(i)
            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))
            rh  = grid_params%xg(i)
cc            etah= (1.+aa_eta*rh**bb_eta)**cc_eta
            etah =res(i,1,1,nxd,nyd,nzd,1,1,1,rrr=rh)/eta

            dummy(i) = dummy(i-1)
cc     .                +(Eoeta/etah*bth+p00*pprof(rh,1))
     .                +(Eoeta/etah*bth+(p0(i+1)-p0(i-1))/2./dr(i))
     .                /(bzh**2+bth**2)*dr(i)
          enddo

          bzz0 = exp(-dummy)

          if (theta /= 0d0) then
            bz_avg = 0d0
            do i=1,nxd
               bz_avg = bz_avg + (bzz0(i)*i + bzz0(i-1)*(i-1))*dr(i)**2
            enddo
            bz_avg = bz_avg/xmax**2
            bzz0 = bth0(nxd)/theta/bz_avg*bzz0
          elseif (bz0 /= 0d0) then
            bzz0 = bz0*bzz0
          else
            call pstop('setEquilibrium'
     $                ,'Cannot find PPNCH equilibrium:'
     $                 //' Both Bz0 and theta are zero')
          endif

          !Check convergence
          err = sqrt(0.5*(sum((bth0-bth0_old)**2)
     .                   +sum((bzz0-bzz0_old)**2))/nxdp)

          if (ilevel > 0 .and. my_rank == 0) then
            write (*,*) '1D Ohmic equilibrium iter =',it
     .                 ,' ; Error=',err
          endif

          dh = xmax/nxd
          if (err < 0.1*dh**2) exit

        enddo

        if (ilevel == 0 .and. my_rank == 0) then
           write (*,'(a,i3,a,1pe10.2)')
     .          ' 1D equilibrium converged in '
     .          ,it,' iterations with error =',err
        endif

        !Check force balance
        do i=1,nxd
          rh  = grid_params%xg(i)

          bzh = 0.5*(bzz0(i)+bzz0(i-1))
          bth = 0.5*(bth0(i)+bth0(i-1))
          dummy(i) = bzh*(bzz0(i)-bzz0(i-1))/dr(i)
     $              +bth/rh*(bth0(i)*i-bth0(i-1)*(i-1))
cc     .              +p00*pprof(rh,1)
     .              +(p0(i+1)-p0(i-1))/2./dr(i)
        enddo

        if (my_rank == 0) then
          write(*,*)
          write(*,'(a,1pe14.7)') ' Force-balance error ='
     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
        endif

        !Check Ohm's law
        do i=1,nxd
          bzh = 0.5*(bzz0(i)+bzz0(i-1))
          bth = 0.5*(bth0(i)+bth0(i-1))
          rh  = grid_params%xg(i)
cc          etah= (1.+aa_eta*rh**bb_eta)**cc_eta
          etah =res(i,1,1,nxd,nyd,nzd,1,1,1,rrr=rh)/eta

          dummy(i) = Eoeta*bzh/etah
     $              -bzh/rh*(bth0(i)*i-bth0(i-1)*(i-1))
     $              +bth*(bzz0(i)-bzz0(i-1))/dr(i)
        enddo

        if (my_rank == 0) then
          write(*,'(a,1pe14.7)') ' Ohms law      error ='
     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
        endif

        !Calculate pinch equilibrium features (theta parameter, Bz(r=0), I0)
        bz_avg = 0d0
        do i=1,nxd
          bz_avg = bz_avg + (bzz0(i)*i + bzz0(i-1)*(i-1))*dr(i)**2
        enddo

        if (my_rank == 0) then
          write (*,*) 
          write (*,999) ' Pinch equilibrium features:'
          write (*,999) '   Theta  =',bth0(nxd)/bz_avg
          write (*,999) '   Bz(r=0)=',bzz0(0)
          write (*,999) '   I0     =',2*pi*xmax*bth0(nxd)
          write (*,999) '   E0/eta =',Eoeta
          write (*,999) '   T_flux =',pi*xmax**2*bz_avg
 999      format (a,f10.7)
        endif

        !Save some equilibrium properties (to calculate bootstrap current)
cc        eq_params(1) = 1d0/grid_params%params(3)       !Aspect ratio
cc        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
cc        eq_params(3) = q0h(nxd)/grid_params%params(3)  !q(r=a)

        !Transfer to integer mesh
        call IntDriver1d(nxdp,xh(0:nxd),bth0(0:nxd)
     .                  ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp)
     .                  ,2,0)
        bth0 = dummy

        call IntDriver1d(nxdp,xh(0:nxd),bzz0(0:nxd)
     .                  ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp)
     .                  ,2,0)
        bzz0 = dummy

        !Extrapolate linearly in radius @ r=0
        bth0(0) = -bth0(1)
        bzz0(0) = bzz0(1)

        !Find pinch flow (normalized to nominal resistivity)
        do i=1,nxd
          rh  = grid_params%xg(i)
cc          etah= (1.+aa_eta*rh**bb_eta)**cc_eta
          etah =res(i,1,1,nxd,nyd,nzd,1,1,1,rrr=rh)/eta

cc          vr0(i) = -(Eoeta*bth0(i)+etah*p00*pprof(rh,1))
          vr0(i) = -(Eoeta*bth0(i)+etah*(p0(i+1)-p0(i-1))/2./dr(i))
     .             /(bzz0(i)**2+bth0(i)**2)
        enddo

        !Find pressure
cc        do i=1,nxd
cc          rh  = grid_params%xg(i)
cc          p0(i) = p00*pprof(rh,0)
cc        enddo

        !Dump equilibrium for eigen code
cc        if (my_rank == 0) then
cc          open(unit=123,file='ntm2d.txt',status='unknown')
cc          write (123,*) nxdp+1
cc          write (123,*) grid_params%xg(0:nxdp)
cc          write (123,*) bth0
cc          write (123,*) bzz0
cc          write (123,*) 0.5*p0
cc          write (123,*) grid_params%params(3),eq_params(2),eq_params(3)
cc          close(123)
cccc          stop
cc        endif

      end subroutine Ohm_equil_1d

c     Ohm_equil_2d
c     ###############################################################
      subroutine Ohm_equil_2d(kk,mm,Eoeta,theta,bth0,bzz0)

c     ---------------------------------------------------------------
c     Solves 2D helical Ohmic equil. equations to find cylindrical
c     components of magnetic field B_theta, B_z. On input:
c       * kk: helical parameter (n/m/R) (defined as theta-kk*z)
c       * mm: poloidal perturbation number
c       * Eoeta: electric field over eta
c       * theta: RFP theta parameter
c       * bth0,bzz0: output magnetic field components (on input, bzz0
c                    contains initial guess, and bzz0(0) contains
c                    value at axis)
c
c     In the calculation below, Psi is at cell centers, and magnetic
c     fields for iteration are at faces. They are transferred to
c     cell centers at end of computation.
c     ---------------------------------------------------------------

        use parameters

        use equilibrium

        use grid

        use app_iosetup

        use mg_solver

        use mg_internal, ONLY: mapMGVectorToArray

        use imposeBCinterface

        use ts_setup, ONLY: vol_wgt

        use orbit

        implicit none

c     Call variables

        integer    :: mm
        real(8)    :: bth0(0:nxdp,0:nydp)
     .               ,bzz0(0:nxdp,0:nydp)
     .               ,Eoeta,theta,kk

c     Local variables

        integer    :: i,j,k,ig,jg,kg,ii,jj,it,bcnd(6),guess
     .               ,ip,im,jp,jm,kp,km
        real(8)    :: bz_avg,bth_avg,x1,err,rh,dh,I0,th
     .               ,bth(0:nydp),bzh(0:nydp),dbz(0:nydp),xih(0:nydp)

        real(8)    :: bz0,rint(0:nxdp,0:nydp),dummy(0:nxdp,0:nydp)
     .               ,bth0_old(0:nxdp,0:nydp),bzz0_old(0:nxdp,0:nydp)
     .               ,xh(0:nxd),dr(0:nxd),psi(0:nxdp,0:nydp,0:nzdp)
     .               ,psi_rhs(nxd*nyd),psiv(nxd*nyd)
     .               ,ffpsi(0:nxdp,0:nydp)

        external   cyl_lap

        !Interfaces
        INTERFACE
          function hel_flux_avg(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: hel_flux_avg,fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function hel_flux_avg
        END INTERFACE

c     Begin program

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '************************************'
          write (*,*) '2D Ohmic equilibrium solver'
          write (*,*) '************************************'
          write (*,*)
        endif

        if (Eoeta == 0d0) then
          messg = 'E0/eta = 0!!'
          call pstop('Ohm_equil_2d',messg)
        endif

c     Initialize iteration (Psi -> integer mesh; bth, bz -> half mesh)

        bz0 = bzz0(0,1)

        !Find half-mesh positions
        do i=0,nxd
          xh(i) = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))
          dr(i) =      grid_params%xg(i+1)-grid_params%xg(i)
        enddo

        !Initial guess for B_theta and psi=integral_r(Bth-k*r*Bz)
        do k=0,nzd+1
          do j=0,nyd+1
            psi (:,j,k) = 0.25  *grid_params%xg(0:nxdp)**2*Eoeta
     .                  - 0.5*kk*grid_params%xg(0:nxdp)**2*bz0
cc            psi (:,j,k) = psi(:,j,k)
cc     .         +0.05*(sin(pi*grid_params%xg(0:nxdp)/xmax)**(mm+1))
cccc     .         +0.05*(0.1+sin(pi*grid_params%xg(0:nxdp)/xmax)**(mm+1))
cc     .                 *sin(grid_params%yg(j)*mm)
cccc            psi (:,j,k) = sin(pi*grid_params%xg(0:nxdp)/xmax)**(mm+1)
cccc     .                   *sin(grid_params%yg(j)*mm)
            bth0(1:nxd,j) = (psi(2:nxdp,j,k)-psi(0:nxd-1,j,k))
     .                    /(dr(1:nxd)+dr(0:nxd-1))
     .                    + kk*grid_params%xg(1:nxd)*bzz0(1:nxd,j)
          enddo
        enddo

        call mapArrayToMGVector(1,nxd,nyd,nzd,psi,psiv,1)

c     Iteration

        do it=1,100

          bth0_old = bth0
          bzz0_old = bzz0

#if 1
          call flux_average(psi,hel_flux_avg)
#endif

          !Find GS equation RHS (at cell centers)
          do j=1,nyd
            do i=1,nxd
              ii  = i + nxd*(j-1)

              x1 = xh(i)
              rh = grid_params%xg(i)
              dh = 0.5*(dr(i)+dr(i-1))

              bzh(j) = 0.5*(bzz0(i,j)+bzz0(i-1,j))

              dbz(j) = kk*(xh(i  )**2*bzz0(i  ,j)
     .                    -xh(i-1)**2*bzz0(i-1,j))/dh/rh
#if 0
              xih(j) = (psi(i+1,j,1)-psi(i-1,j,1))/(dr(i)+dr(i-1))
              bth(j) = xih(j) + kk*rh*bzh(j)
              
              psi_rhs(ii)= Eoeta*bzh(j)**2/(bzh(j)**2+bth(j)**2)-dbz(j)
#else
              psi_rhs(ii)= ffpsi(i,j)*bzh(j) - dbz(j)
#endif
              if (vol_wgt) psi_rhs(ii)=psi_rhs(ii)*rh*dh*2*pi/nyd
            enddo
          enddo

          !Solve GS Poisson eq (find psi at cell centers)
          !Define BCs
          bcnd = bcond
          where (bcnd == DEF) bcnd = DIR

          guess = 1
          call cSolver(1,nxd*nyd,psi_rhs,psiv,bcnd
     .                ,1,ilevel,guess,cyl_lap,vol_wgt
     .                ,tol=1d-2
     .                ,gm_driver=.true.
     .                ,gm_it=10
     .                ,mg_order_res  = 0
     .                ,mg_order_prol = 0
     .                ,mg_smooth     = 'jb'
     .                ,mg_gm_coarse_solve  =.true.
     .                ,mg_vcyc       = 1
     .                ,sm_ncolors    = 2  !Star stencil
     .                ,sm_zebra_relax=.true.
     .                )

          !Map to array (w BCs)
          call mapMGVectorToArray(0,1,psiv,nxd,nyd,nzd,psi,1,.false.)

          call setBC(IRHO,nxd,nyd,nzd,psi,zeros,bcnd,1,1,1,iorder=1)

          !Perform integral of Bz (at faces)
          dummy(0,:) = 0d0
          do i=1,nxd
            xih = (psi (i+1,:,1)-psi (i-1,:,1))/(dr(i)+dr(i-1))
#if 0
            rh  = grid_params%xg(i)
            bzh = 0.5*(bzz0(i,:)+bzz0(i-1,:))
            bth = xih + kk*rh*bzh

            dummy(i,:)=dummy(i-1,:)+Eoeta*bzh*xih/(bzh**2+bth**2)*dr(i)
#else
            dummy(i,:)=dummy(i-1,:) + ffpsi(i,:)*xih*dr(i)
#endif
          enddo

          do i=0,nxd
            bzz0(i,:)=( bz0
     .                 -kk*xh(i)*(psi(i+1,:,1)-psi(i,:,1))/dr(i)
     .                 -dummy(i,:) )/(1+(kk*xh(i))**2)
          enddo

          !Find B_theta (at faces)
          do i=0,nxd
            bth0(i,:) = (psi(i+1,:,1)-psi(i,:,1))/dr(i)
     .                 + kk*xh(i)*bzz0(i,:)
          enddo

          !Check convergence
          err = sqrt(0.5*(sum((bth0-bth0_old)**2)
     .                   +sum((bzz0-bzz0_old)**2))/nxdp/nydp)

          if (ilevel > 0 .and. my_rank == 0) then
            write (*,*) '2D Ohmic equilibrium iter =',it
     .                 ,' ; Error=',err
          endif

          dh = (xmax-xmin)/nxd
          if (err < 0.1*dh**2) exit

        enddo

        if (ilevel == 0 .and. my_rank == 0) then
           write (*,'(a,i3,a,1pe10.2)')
     .          ' 2D equilibrium converged in '
     .          ,it,' iterations with error =',err
        endif

c     Check force balance

        do i=1,nxd
          rh  = grid_params%xg(i)

          bzh = 0.5*(bzz0(i,:)+bzz0(i-1,:))
          bth = 0.5*(bth0(i,:)+bth0(i-1,:))
          dummy(i,:) =  bzh*(bzz0(i,:)-bzz0(i-1,:))/dr(i)
     $                 +bth/rh*(bth0(i,:)*i-bth0(i-1,:)*(i-1))
        enddo

        if (my_rank == 0) then
          write(*,*)
          write(*,'(a,1pe14.7)') ' Force-balance error ='
     .                        ,sqrt(sum(dummy(1:nxd,1:nyd)**2)/nxd/nyd)
        endif

c     Check Ohm's law

        do i=1,nxd
          bzh = 0.5*(bzz0(i,:)+bzz0(i-1,:))
          bth = 0.5*(bth0(i,:)+bth0(i-1,:))
          rh  = grid_params%xg(i)
          dummy(i,:) = Eoeta*bzh
     $              -bzh/rh*(bth0(i,:)*i-bth0(i-1,:)*(i-1))
     $              +bth*(bzz0(i,:)-bzz0(i-1,:))/dr(i)
        enddo

        if (my_rank == 0) then
          write(*,'(a,1pe14.7)') ' Ohms law      error ='
     .                        ,sqrt(sum(dummy(1:nxd,1:nyd)**2)/nxd/nyd)
        endif

c     Calculate pinch equilibrium features (theta parameter, Bz(r=0), I0)

        bz_avg = 0d0
        do j=1,nyd
          do i=1,nxd
            bz_avg = bz_avg
     .             + (bzz0(i,j)*i + bzz0(i-1,j)*(i-1))*dr(i)**2*2*pi/nyd
          enddo
        enddo
        bz_avg = 0.5*bz_avg/pi/xmax**2

        I0 = 0d0
        do j=1,nyd
          I0 = I0 + bth0(nxd,j)*2*pi*xmax/nyd
        enddo
        bth_avg = I0/2/pi/xmax

        if (my_rank == 0) then
          write (*,*) 
          write (*,999) ' Pinch equilibrium features:'
          write (*,999) '   Theta  =',bth_avg/bz_avg
          write (*,999) '   Bz(r=0)=',bzz0(0,1)
          write (*,999) '   I0     =',I0
          write (*,999) '   E0/eta =',Eoeta
          write (*,999) '   T_flux =',pi*xmax**2*bz_avg
 999      format (a,f10.7)
        endif

        !Save some equilibrium properties (to calculate bootstrap current)
cc        eq_params(1) = 1d0/grid_params%params(3)       !Aspect ratio
cc        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
cc        eq_params(3) = q0h(nxd)/grid_params%params(3)  !q(r=a)

c       Transfer magnetic field components to integer mesh

        do j=0,nydp
          call IntDriver1d(nxdp,xh(0:nxd),bth0(0:nxd,j)
     .                  ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp,j)
     .                  ,2,0)
        enddo
        bth0 = dummy

        do j=0,nydp
          call IntDriver1d(nxdp,xh(0:nxd),bzz0(0:nxd,j)
     .                  ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp,j)
     .                  ,2,0)
        enddo
        bzz0 = dummy

        !Extrapolate linearly in radius @ r=0
        bth0(0,:) = -bth0(1,:)
        bzz0(0,:) =  bzz0(1,:)

        !Dump equilibrium for eigen code
cc        if (my_rank == 0) then
cc          open(unit=123,file='ntm2d.txt',status='unknown')
cc          write (123,*) nxdp+1
cc          write (123,*) grid_params%xg(0:nxdp)
cc          write (123,*) bth0
cc          write (123,*) bzz0
cc          write (123,*) 0.5*p0
cc          write (123,*) grid_params%params(3),eq_params(2),eq_params(3)
cc          close(123)
cccc          stop
cc        endif

      contains

c     flux_average
c     ##############################################################
      subroutine flux_average(psi,flux_avg)

c     --------------------------------------------------------------
c     Performs flux averages according to psi and flux_avg routine.
c     --------------------------------------------------------------

c     Call variables

        real(8) :: psi(0:nxdp,0:nydp,0:nzdp)

        INTERFACE
          function flux_avg(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: flux_avg,fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function flux_avg
        END INTERFACE

c     Local variables

        real(8),allocatable,dimension(:,:,:) :: bx,by,bz
     .                                         ,bxcar,bycar,bzcar,jac
        real(8),allocatable,dimension(:,:,:,:) :: xmap

c     Begin program

c     Initialize orbit integrator work space

        allocate(bx   (0:nxd+1,0:nyd+1,0:nzd+1)
     .          ,by   (0:nxd+1,0:nyd+1,0:nzd+1)
     .          ,bz   (0:nxd+1,0:nyd+1,0:nzd+1)
     .          ,bxcar(0:nxd+1,0:nyd+1,0:nzd+1)
     .          ,bycar(0:nxd+1,0:nyd+1,0:nzd+1)
     .          ,bzcar(0:nxd+1,0:nyd+1,0:nzd+1)
     .          ,jac  (0:nxd+1,0:nyd+1,0:nzd+1)
     .          ,xmap (0:nxd+1,0:nyd+1,0:nzd+1,3))

c diag ****
cc        open(unit=110,file='debug.bin',form='unformatted'
cc     .      ,status='replace')
cc        call contour(psi(0:nxd+1,0:nyd,1),nxd+2,nyd+1
cc     .              ,0d0,xmax,0d0,ymax,0,110)
cc        close(110)
cccc        stop
c diag ****

c     Orbit integrator setup: cnv components (with BCs)

        do k=0,nzd+1
          do j=0,nyd+1
            do i=0,nxd+1
              rh = grid_params%xg(i)

              ip = min(i+1,nxd+1)
              im = max(i-1,0)
              jp = min(j+1,nyd+1)
              jm = max(j-1,0)
              kp = min(k+1,nzd+1)
              km = max(k-1,0)

              !Cell center cnv components
              bx(i,j,k) =-(psi(i,jp,k)-psi(i,jm,k))
     .                   /(grid_params%yg(jp)-grid_params%yg(jm))
              by(i,j,k) = (psi(ip,j,k)-psi(im,j,k))
     .                   /(grid_params%xg(ip)-grid_params%xg(im))
              bz(i,j,k) = rh*0.5*(bzz0(i,j)+bzz0(im,j))

            enddo
          enddo
        enddo

        !Periodic BCs
        bx(:,0    ,:)=bx(:,nyd,:)
        bx(:,nyd+1,:)=bx(:,1  ,:)
        by(:,0    ,:)=by(:,nyd,:)
        by(:,nyd+1,:)=by(:,1  ,:)
        bz(:,0    ,:)=bz(:,nyd,:)
        bz(:,nyd+1,:)=bz(:,1  ,:)

        bx(:,:,0    )=bx(:,:,nzd)
        bx(:,:,nzd+1)=bx(:,:,1  )
        by(:,:,0    )=by(:,:,nzd)
        by(:,:,nzd+1)=by(:,:,1  )
        bz(:,:,0    )=bz(:,:,nzd)
        bz(:,:,nzd+1)=bz(:,:,1  )

        !SP BCs
        do k=0,nzd+1
          do j=0,nyd+1
            jj = mod(j+nyd/2,nyd)
            if (jj == 0) jj = nyd

            bx(0,j,k) = bx(1,jj,k)
            by(0,j,k) =-by(1,jj,k)
            bz(0,j,k) =-bz(1,jj,k)
          enddo
        enddo

        !Cartesian components, maps
        do k=0,nzd+1
          do j=0,nyd+1
            do i=0,nxd+1
              rh = grid_params%xg(i)

              !Cartesian components
              th = grid_params%yg(j)

              bxcar(i,j,k) = bx(i,j,k)/rh*cos(th)
     .                     -(by(i,j,k)+kk*bz(i,j,k))*sin(th)
              bycar(i,j,k) = bx(i,j,k)/rh*sin(th)
     .                     +(by(i,j,k)+kk*bz(i,j,k))*cos(th)
              bzcar(i,j,k) = bz(i,j,k)/rh

              !Map (cylindrical)
              jac(i,j,k)   = rh
              xmap(i,j,k,1)= rh*cos(th)
              xmap(i,j,k,2)= rh*sin(th)
              xmap(i,j,k,3)= grid_params%zg(k)
            enddo
          enddo
        enddo

c     Orbit setup (input Az directly)

        bx = 0d0
        by = 0d0
        bz = -psi
        call orbit_setup(nxd+2,nyd+2,nzd+2
     .                  ,grid_params%xg
     .                  ,grid_params%yg
     .                  ,grid_params%zg
     .                  ,bx,by,bz,bxcar,bycar,bzcar,jac,xmap
     .                  ,dtime=1d-1,solen=.false.,B_input_is_A=.true.
     .                  ,r_min=1d-3)

        !Orbit average (at cell centers)
        if (my_rank == 0) then
          if (ilevel == 1) then
            write (*,FMT="(a)",ADVANCE="NO") 'Orbit averaging: '
          elseif (ilevel == 0) then
            write (*,*) 'Orbit averaging...'
          endif
        endif

c     Orbit integration

        k = 1
        do j=1,nyd
          do i=1,nxd

            x0 = grid_params%xg(i)
            y0 = grid_params%yg(j)
            z0 = grid_params%zg(k)

            bzz_avg = 0d0
            b2_avg  = 0d0
            ssmax    = 0d0

            call orbit_find(ilevel-2,x0,y0,z0,.true.
     .                     ,line_int=flux_avg)

            if ((my_rank == 0).and.(ilevel == 1))
     .                  write (*,FMT="(a)",ADVANCE="NO") '+'

            ffpsi(i,j) = Eoeta*bzz_avg/b2_avg
cc            ffpsi(i,j) = ssmax  !diag

          enddo
        enddo

        ffpsi(:,0    )=ffpsi(:,nyd)
        ffpsi(:,nyd+1)=ffpsi(:,1  )

c diag ****
cc          open(unit=110,file='debug.bin',form='unformatted'
cc     .        ,status='replace')
cc          call contour(psi(0:nxd+1,0:nyd,1),nxd+2,nyd+1
cc     .                ,0d0,xmax,0d0,ymax,0,110)
cc          close(110)
cccc          stop
cc
cccc          write (*,*) 'here'
cc          open(unit=110,file='debug2.bin',form='unformatted'
cc     .        ,status='replace')
cc          call contour(ffpsi(1:nxd,1:nyd),nxd,nyd
cc     .                ,0d0,xmax,0d0,ymax,0,110)
cc          close(110)
cc          stop
c diag ****

c     Free up orbit integral work space

        if ((my_rank == 0).and.(ilevel == 1)) write (*,*)

        call orbit_destroy

        deallocate(bx,by,bz,bxcar,bycar,bzcar,jac,xmap)

      end subroutine flux_average

      end subroutine Ohm_equil_2d

c     hel_flux_avg
c     ###################################################################
      function hel_flux_avg(fwd,xold,x,yold,y,zold,z,sold,ss)
     .         result(exloop)

c     -------------------------------------------------------------------
c     Performs flux average of various quantities for 2D helical
c     equilibrium
c     -------------------------------------------------------------------

        use orbit, xo => xold, yo => yold, zo => zold

        use equilibrium

        use constants

        implicit none

c     Call variables

        logical :: exloop,fwd
        real(8) :: xold,x,yold,y,zold,z,sold,ss

c     Local variables

        real(8) :: xh,yh,zh,sh,bb,bx,by,bz
     .            ,xs,ys,zs,ts,y00,dy0
        save dy0

c     Begin program

c     Check whether orbit should terminate

        if (yold == y0) dy0 = y-y0

        call ss_init(1,.false.,0d0,.true.,y0,.false.,0d0)

        ts = 0d0
        call sschk(1,0d0,sold,ss,0d0,xold,x,0d0,yold,y,0d0,zold,z
     .            ,xs,ys,zs,ts)

c     Update current position

        if ( ts /= 0d0 .and. (y-yold)*dy0 > 0d0 ) then
          x = xs ; y = ys ; z = zs
          exloop = .true.
        else
          ts = ss
          exloop = .false.
        endif

c     Find middle points

        xh = 0.5*(x +xold)
        yh = 0.5*(y +yold)
        zh = 0.5*(z +zold)
        sh = 0.5*(ts+sold)

c     Find magnetic field Cartesian components

        call getB(xh,yh,zh,bx,by,bz,.false.,.true.)

        bb = sqrt(bx**2+by**2+bz**2)

c     Perform flux averages

        bzz_avg = bzz_avg + bz*(ts-sold)/bb
        b2_avg  = b2_avg  + bb*(ts-sold)
        ssmax   = ssmax   + (ts-sold)/bb

      end function hel_flux_avg

c     Ohm_GS_equil_2d
c     ###############################################################
      subroutine Ohm_GS_equil_2d(R0,Eoeta,b_tor,nx,ny,nz,b1,b2,b3
     .                          ,v1,v2,pp)

c     ---------------------------------------------------------------
c     Solves 2D Ohmic Grad-Shafranov equilibrium equations to find 
c     contravariant components of magnetic field. On input:
c       * R0: Major radius
c       * Eoeta: electric field over eta
c       * b_tor: toroidal beta
c       * nx,ny,nz: local (processor) dimensions
c       * b1,b2,b3: output magnetic field contravariant components 
c       * v1,v2: output velocity contravariant components
c
c     In the calculation below, psi is at cell centers, and B_3=I(psi)
c     is defined at faces.
c     ---------------------------------------------------------------

        use parameters

        use equilibrium

        use grid

        use grid_debug

        use app_iosetup

        use mg_solver

        use mg_internal, ONLY: mapMGVectorToArray

        use imposeBCinterface

        use ts_setup, ONLY: vol_wgt

        use operators, ONLY: laplacian

        use transport_params

        implicit none

c     Call variables

        integer    :: nx,ny,nz
        real(8)    :: b1(0:nx+1,0:ny+1)
     .               ,b2(0:nx+1,0:ny+1)
     .               ,b3(0:nx+1,0:ny+1)
     .               ,v1(0:nx+1,0:ny+1)
     .               ,v2(0:nx+1,0:ny+1)
     .               ,pp(0:nx+1,0:ny+1),Eoeta,R0,b_tor

c     Local variables

        integer    :: i,j,k,ig,jg,kg,ii,it,bcnd(6),igrid,guess
     .               ,iglobal,jglobal,pos(2)

        real(8)    :: bz_avg,bth_avg,gserr,oherr,rh,dh,dt,I0,bnorm2
     .               ,psir,psith,x1,y1,z1,IIh,spsi,psi0,psia,gpsi2
     .               ,cc,p00,etah,lpsi0,dsdpsi,dpdpsi,spsip,spsim
     .               ,dpdr,dpdth,f_1,f_2,f_3,j1,j2,j3,jac,e1,e2,e3
     .               ,psip,psim,q0,qa,tor_flx,bz2_avg,p_avg,b_pol
     .               ,cov(3),e4,e5,e6

        real(8)    :: xh(0:nxd),dr(0:nxd),dth(0:nyd)

        real(8)    :: psi    (0:nx+1,0:ny+1,0:nz+1)
     .               ,iR2    (0:nx+1,0:ny+1,0:nz+1)
     .               ,qprof  (0:nx+1,0:ny+1)
     .               ,dummy  (0:nxdp,0:ny+1)
     .               ,psig   (0:nxdp,0:ny+1)
     .               ,III    (0:nxdp,0:ny+1)
     .               ,iR2psi (0:nxdp,0:nydp)
     .               ,B2psi  (0:nxdp,0:nydp)
     .               ,RR     (0:nxdp,0:nydp)

        real(8),allocatable,dimension(:) :: psi_rhs,psiv,mpivec,errv

        real(8)    :: pprof
        external   del_star,pprof

        !Interfaces
        INTERFACE
          function tor_flux_avg(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: tor_flux_avg,fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function tor_flux_avg
        END INTERFACE

c     Begin program

        igrid = 1

        allocate(psi_rhs(nx*ny*nz),psiv(nx*ny*nz),errv(nx*ny*nz))

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '************************************'
          write (*,*) '2D Ohmic GS equilibrium solver'
          write (*,*) '************************************'
          write (*,*)
        endif

        if (Eoeta == 0d0) then
          messg = 'E0/eta = 0!!'
          call pstop('Ohm_GS_equil_2d',messg)
        endif

c       Initialize iteration

        !Find half-mesh positions
        do i=0,nxd
          xh(i) = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))
          dr(i) =      grid_params%xg(i+1)-grid_params%xg(i)
        enddo

        do j=0,nyd
          dth(j)=      grid_params%yg(j+1)-grid_params%yg(j)
        enddo

        !Initial guess for psi
        psi = 0d0

        !Initalize toroidal radius (R=R0+r*cos(theta)) and 1/R^2
        do j=0,nydp
          do i=0,nxdp
cc            call getCartesianCoordinates(i,j,1,igrid,igrid,igrid
cc     .                                  ,ig,jg,kg,x1,y1,z1)
cc            RR(i,j) = sqrt(x1**2+y1**2)

            RR(i,j) = R0 +  grid_params%xg(i)*cos(grid_params%yg(j))
          enddo
        enddo

        if (minval(RR) < 0d0) then
          call pstop('Ohm_GS_equil_2d','r > R in torus')
        endif

        do j=0,ny+1
          do i=0,nx+1
            iglobal = i + grid_params%ilo(igrid) - 1
            jglobal = j + grid_params%jlo(igrid) - 1

            iR2(i,j,:) = 1d0/RR(iglobal,jglobal)**2
          enddo
        enddo

        !Initial guess for III=B_phi*R
        do j=0,ny+1
          do iglobal=0,nxdp
            jglobal = j + grid_params%jlo(igrid) - 1

            III(iglobal,j) = RR(iglobal,jglobal)   !B_phi=1 uniform
          enddo
        enddo

        !Maximum pressure
        p00 = 0.5*b_tor

c       Iteration

        psig = 0d0

        psi0 =-1d0
        psia = 0d0
        oherr = 0d0

        do it=1,100

cc          if (it > 1) call flux_average(psig,tor_flux_avg)

          !Find GS equation RHS (at cell centers)
          do k=1,nz
            do j=1,ny
              do i=1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              iglobal = i + grid_params%ilo(igrid) - 1
              jglobal = j + grid_params%jlo(igrid) - 1

              x1 = xh(iglobal)
              rh = grid_params%xg(iglobal)
              dh = 0.5*(dr (iglobal)+dr (iglobal-1))
              dt = 0.5*(dth(jglobal)+dth(jglobal-1))

              IIh = 0.5*(III(iglobal,j)+III(iglobal-1,j))

cc              psir   = (psi(i+1,j,k)-psi(i-1,j,k))
cc     .                /(dr (iglobal)+dr (iglobal-1))
cc              psith  = (psi(i,j+1,k)-psi(i,j-1,k))
cc     .                /(dth(jglobal)+dth(jglobal-1))/rh
cc              bnorm2 = IIh**2 + psir**2 + psith**2
              psir   = (psi(i+1,j,1)-psi(i-1,j,1))/2./dh
              psith  = (psi(i,j+1,1)-psi(i,j-1,1))/2./dt/rh
              gpsi2  = psir**2 + psith**2
cc              bnorm2 = IIh**2 + gpsi2

              etah =res(iglobal,1,1,nxd,nyd,nzd,1,1,1,rrr=rh)/eta

              spsi = sqrt(1d0-psi(i,j,1)/psi0)+1d-8
              dsdpsi =-0.5/spsi/psi0
              dpdpsi = dsdpsi*p00*pprof(spsi,1)

              spsip = sqrt(1d0-psi(i+1,j,1)/psi0)
              spsim = sqrt(1d0-psi(i-1,j,1)/psi0)
              dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dh
cc              psip  = 0.5*(psi(i+1,j,1)+psi(i,j,1))
cc              psim  = 0.5*(psi(i-1,j,1)+psi(i,j,1))
cc              spsip = sqrt(1d0-psip/psi0)
cc              spsim = sqrt(1d0-psim/psi0)
cc              dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/dh

              spsip = sqrt(1d0-psi(i,j+1,1)/psi0)
              spsim = sqrt(1d0-psi(i,j-1,1)/psi0)
              dpdth = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dt/rh
cc              psip  = 0.5*(psi(i,j+1,1)+psi(i,j,1))
cc              psim  = 0.5*(psi(i,j-1,1)+psi(i,j,1))
cc              spsip = sqrt(1d0-psip/psi0)
cc              spsim = sqrt(1d0-psim/psi0)
cc              dpdth = p00*(pprof(spsip,0)-pprof(spsim,0))/dt/rh

              psi_rhs(ii)= (Eoeta/etah*(IIh/RR(iglobal,jglobal))**2
     .                     *iR2psi(iglobal,jglobal)
     .                      /B2psi(iglobal,jglobal))
     .                    +dpdpsi*((IIh/RR(iglobal,jglobal))**2
     .                             /B2psi(iglobal,jglobal) - 1)

              if (vol_wgt) psi_rhs(ii)=psi_rhs(ii)
     .                                *gmetric%grid(igrid)%dvol(i,j,k)

              errv(ii) = (laplacian(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                             ,psi,dff=iR2,vol=vol_wgt)
     .                   -psi_rhs(ii))*RR(iglobal,jglobal)  !Multiply by R to ensure proper limit R->inf

              oherr = oherr
     .               +(Eoeta/etah*IIh/RR(iglobal,jglobal)**2*psir
     .                +IIh*dpdr
     .                +bnorm2/RR(iglobal,jglobal)**2
     .                *(III(iglobal,j)-III(iglobal-1,j))/dh)**2
              enddo
            enddo
          enddo

          !Calculate L2 norm of error
          gserr = sqrt(dot(nx*ny*nz,errv,errv))
          oherr = sqrt(oherr)

          !Solve GS Poisson eq (find psi at cell centers)
          !Define BCs
          bcnd = bcond
          where (bcnd == DEF) bcnd = DIR

          guess = 1

          call cSolver(1,nx*ny*nz,psi_rhs,psiv,bcnd
     .              ,1,ilevel,guess,del_star,vol_wgt
     .              ,tol=5d-1
     .              ,gm_driver=.true.
     .              ,mg_order_res  = 0
     .              ,mg_order_prol = 1
     .              ,mg_gm_coarse_solve=.true.
     .              ,mg_smooth     = 'jb'
     .              ,mg_vcyc       = 10
     .              ,sm_ncolors    = 2  !Star stencil
     .              )
cc     .               ,ncolors=2,order_res=0,order_prol=1
cccc     .               ,smooth='gm',gm_coarse_solve=.true.
cc     .               ,smooth='jb',gm_coarse_solve=.true.
cc     .               ,mvcyc=10,tol=5d-1,sm_it=5)
          !Map to array (w BCs)
          call mapMGVectorToArray(0,1,psiv,nx,ny,nz,psi,igrid
     .                           ,.false.)

          call setBC(IRHO,nx,ny,nz,psi,zeros,bcnd,igrid,igrid,igrid
     .              ,iorder=2)

          !Gather global psi in radius
#if defined(petsc)
          if (npx > 1) then
            allocate(mpivec(npx*(nx+2)*(ny+2)))
            call MPI_Allgather(psi(0:nx+1,0:ny+1,1)
     .                             ,(nx+2)*(ny+2),MPI_DOUBLE_PRECISION
     .                      ,mpivec,(nx+2)*(ny+2),MPI_DOUBLE_PRECISION
     .                      ,MPI_COMM_RAD,mpierr)

            do j=1,npx
              psig(nx*(j-1):nx*j+1,:)
     .         = reshape(mpivec(1+(nx+2)*(ny+2)*(j-1):(nx+2)*(ny+2)*j)
     .                   ,(/nx+2,ny+2/))
            enddo
            deallocate(mpivec)
cc            allocate(mpivec(npx*(nx+2)))
cc            do j=0,ny+1
cc              call MPI_Allgather(psi(0:nx+1,j,1)
cc     .                                 ,(nx+2),MPI_DOUBLE_PRECISION
cc     .                          ,mpivec
cc     .                                 ,(nx+2),MPI_DOUBLE_PRECISION
cc     .                          ,MPI_COMM_RAD,mpierr)
cc
cc              do i=1,npx
cc                psig(nx*(i-1):nx*i+1,j)= mpivec(1+(nx+2)*(i-1):(nx+2)*i)
cc              enddo
cc            enddo
cc            deallocate(mpivec)
          else
            psig = psi(:,:,1)
          endif
#else
          psig = psi(:,:,1)
#endif
          !Perform integral of III (at half-mesh on global grid in radius)

          psi0 = minval(psi)
#if defined (petsc)
          lpsi0 = psi0
          call MPI_Allreduce(lpsi0,psi0,1,MPI_DOUBLE_PRECISION
     .                      ,MPI_MIN,MPI_COMM_WORLD,mpierr)
#endif

          do j=1,ny
            jglobal = j + grid_params%jlo(igrid) - 1

            dummy(0,j) = 0d0
            do iglobal=1,nxd
              IIh = 0.5*(III(iglobal,j)+III(iglobal-1,j))

              rh = grid_params%xg(iglobal)
              dh = 0.5*(dr (iglobal)+dr (iglobal-1))
              dt = 0.5*(dth(jglobal)+dth(jglobal-1))

              psir   = (psig(iglobal+1,j)-psig(iglobal-1,j))/2./dh
              psith  = (psig(iglobal,j+1)-psig(iglobal,j-1))/2./dt/rh
              bnorm2 = IIh**2 + psir**2 + psith**2

cc              spsi = sqrt(1d0-psig(iglobal,j)/psi0)+1d-8
cc              dsdpsi =-0.5/spsi/psi0
cc              dpdpsi = dsdpsi*p00*pprof(spsi,1)

              psip  = 0.5*(psig(iglobal+1,j)+psig(iglobal,j))
              psim  = 0.5*(psig(iglobal-1,j)+psig(iglobal,j))
              spsip = sqrt(1d0-psip/psi0)
              spsim = sqrt(1d0-psim/psi0)
              dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/dh

              etah  = res(iglobal,1,1,nxd,nyd,nzd,1,1,1,rrr=rh)/eta

              dummy(iglobal,j) = dummy(iglobal-1,j)
     .                         - RR(iglobal,jglobal)**2/bnorm2
     .                          *(Eoeta*psir/etah/RR(iglobal,jglobal)**2
     .                          + dpdr)*dr(iglobal)
cc     .                         - RR(iglobal,jglobal)**2*psir/bnorm2
cc     .                          *(Eoeta/etah/RR(iglobal,jglobal)**2
cc     .                          + dpdpsi)*dr(iglobal)
            enddo
          enddo

          III = R0*exp(dummy)

          !Check convergence
          if (ilevel > 0 .and. my_rank == 0) then
            write (*,*) '2D Ohmic equilibrium iter =',it
     .                 ,' ; GS error=',gserr,' Ohm error=',oherr
            if (ilevel > 1) write (*,*)
          endif

          dh = xmax/nxd
          if (gserr < 0.001*dh**2) exit

        enddo

        if (ilevel == 0 .and. my_rank == 0) then
           write (*,'(a,i3,a,1pe10.2)')
     .          ' 2D equilibrium converged in '
     .          ,it,' iterations with error =',gserr
        endif

c       Find vectors and scalars

        do j=1,ny
          do i=1,nx
            iglobal = i + grid_params%ilo(igrid) - 1
            jglobal = j + grid_params%jlo(igrid) - 1

            rh = grid_params%xg(iglobal)
            dh = 0.5*(dr (iglobal)+dr (iglobal-1))
            dt = 0.5*(dth(jglobal)+dth(jglobal-1))

            IIh    = 0.5*(III(iglobal,j)+III(iglobal-1,j))
            psir   = (psi(i+1,j,1)-psi(i-1,j,1))/2./dh
            psith  = (psi(i,j+1,1)-psi(i,j-1,1))/2./dt/rh
            bnorm2 = (IIh**2+psir**2+psith**2)/RR(iglobal,jglobal)**2

            !Find contravariant magnetic field components in integer mesh
            b1(i,j) =-psith*rh
            b2(i,j) = psir
            b3(i,j) = IIh*rh/RR(iglobal,jglobal)

            !Find contravariant velocity components
            etah =res(iglobal,jglobal,1,nxd,nyd,nzd,1,1,1,rrr=rh)/eta

cc            spsi = sqrt(1d0-psi(i,j,1)/psi0)+1d-8
cc            dsdpsi =-0.5/spsi/psi0
cc            dpdpsi = dsdpsi*p00*pprof(spsi,1)
cc
cc            cc = (Eoeta/RR(iglobal,jglobal)**2+etah*dpdpsi)
cc     .           *RR(iglobal,jglobal)/bnorm2
cc
cc            v1(i,j) =-cc*rh*psir
cc            v2(i,j) =-cc   *psith

cc            psip  = 0.5*(psi(i+1,j,1)+psi(i,j,1))
cc            psim  = 0.5*(psi(i-1,j,1)+psi(i,j,1))
cc            spsip = sqrt(1d0-psip/psi0)
cc            spsim = sqrt(1d0-psim/psi0)
cc            dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/dh
            spsip = sqrt(1d0-psi(i+1,j,1)/psi0)
            spsim = sqrt(1d0-psi(i-1,j,1)/psi0)
            dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dh

cc            psip  = 0.5*(psi(i,j+1,1)+psi(i,j,1))
cc            psim  = 0.5*(psi(i,j-1,1)+psi(i,j,1))
cc            spsip = sqrt(1d0-psip/psi0)
cc            spsim = sqrt(1d0-psim/psi0)
cc            dpdth = p00*(pprof(spsip,0)-pprof(spsim,0))/dt/rh
            spsip = sqrt(1d0-psi(i,j+1,1)/psi0)
            spsim = sqrt(1d0-psi(i,j-1,1)/psi0)
            dpdth = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dt/rh

            v1(i,j) =-(Eoeta*psir /RR(iglobal,jglobal)**2+etah*dpdr )
     .                *RR(iglobal,jglobal)*rh/bnorm2
            v2(i,j) =-(Eoeta*psith/RR(iglobal,jglobal)**2+etah*dpdth)
     .                *RR(iglobal,jglobal)/bnorm2

            !Find pressure
            spsi = sqrt(1d0-psi(i,j,1)/psi0)
            pp(i,j) = p00*pprof(spsi,0)
          enddo
        enddo

c       Calculate GS equilibrium properties

        !Convergence properties
        e1 = 0d0
        e2 = 0d0
        e3 = 0d0
        e4 = 0d0
        e5 = 0d0
        e6 = 0d0

        tor_flx = 0d0
        p_avg   = 0d0
        bz2_avg = 0d0

        do j=2,ny-1
          do i=2,nx-1
            iglobal = i + grid_params%ilo(igrid) - 1
            jglobal = j + grid_params%jlo(igrid) - 1

            rh = grid_params%xg(iglobal)
            dh = 0.5*(dr (iglobal)+dr (iglobal-1))
            dt = 0.5*(dth(jglobal)+dth(jglobal-1))

            IIh    = 0.5*(III(iglobal,j)+III(iglobal-1,j))
            psir   = (psi(i+1,j,1)-psi(i-1,j,1))/2./dh
            psith  = (psi(i,j+1,1)-psi(i,j-1,1))/2./dt/rh
            bnorm2 = (IIh**2+psir**2+psith**2)/RR(iglobal,jglobal)**2

            jac = rh*RR(iglobal,jglobal)

            !Find pressure derivatives
            spsip = sqrt(1d0-psi(i+1,j,1)/psi0)
            spsim = sqrt(1d0-psi(i-1,j,1)/psi0)
            dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dh
cc            psip  = 0.5*(psi(i+1,j,1)+psi(i,j,1))
cc            psim  = 0.5*(psi(i-1,j,1)+psi(i,j,1))
cc            spsip = sqrt(1d0-psip/psi0)
cc            spsim = sqrt(1d0-psim/psi0)
cc            dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/dh

            spsip = sqrt(1d0-psi(i,j+1,1)/psi0)
            spsim = sqrt(1d0-psi(i,j-1,1)/psi0)
            dpdth = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dt
cc            psip  = 0.5*(psi(i,j+1,1)+psi(i,j,1))
cc            psim  = 0.5*(psi(i,j-1,1)+psi(i,j,1))
cc            spsip = sqrt(1d0-psip/psi0)
cc            spsim = sqrt(1d0-psim/psi0)
cc            dpdth = p00*(pprof(spsip,0)-pprof(spsim,0))/dt

            !Find contravariant current components
            j1 = (III(iglobal,j+1)+III(iglobal-1,j+1)
     .           -III(iglobal,j-1)-III(iglobal-1,j-1))/4./dt
            j2 =-(III(iglobal,j)  -III(iglobal-1,j))/dh
            j3 = jac*laplacian(i,j,1,nx,ny,nz,igrid,igrid,igrid
     .                        ,psi,dff=iR2,vol=.false.)

            !Find covariant EOM force components
            f_1 = (j2*b3(i,j)-j3*b2(i,j))/jac - dpdr
            f_2 = (j3*b1(i,j)-j1*b3(i,j))/jac - dpdth
            f_3 = (j1*b2(i,j)-j2*b1(i,j))/jac

cc            write (*,*) 'jxb equil',i,j,f_1+dpdr,f_2+dpdth,f_3

            !Accumulate errors
            e1 = e1 + f_1**2
            e2 = e2 + f_2**2
            e3 = e3 + f_3**2

            !Find covariant Ohm's law components
            etah =res(iglobal,1,1,nxd,nyd,nzd,1,1,1,rrr=rh)/eta

            cov = matmul(gmetric%grid(1)%gsub(i,j,1,:,:),(/ j1,j2,j3/))

            f_1 =      (v2(i,j)*b3(i,j)                )/jac-etah*cov(1)
            f_2 =      (               -v1(i,j)*b3(i,j))/jac-etah*cov(2)
            f_3 =Eoeta+(v1(i,j)*b2(i,j)-v2(i,j)*b1(i,j))/jac-etah*cov(3)

            !Accumulate errors
            e4 = e4 + f_1**2
            e5 = e5 + f_2**2
            e6 = e6 + f_3**2

            !Q-profile
            qprof(i,j) = b3(i,j)/b2(i,j)

            !Toroidal flux
            tor_flx = tor_flx + b3(i,j)*dr(iglobal)*dth(jglobal)

            !Bz^2 average
            bz2_avg = bz2_avg + (b3(i,j))**2/rh*dr(iglobal)*dth(jglobal)

            !Pressure average
            p_avg = p_avg + pp(i,j)*rh*dr(iglobal)*dth(jglobal)

          enddo
        enddo

#if defined (petsc)
        lpsi0 = e1
        call MPI_Allreduce(lpsi0,e1,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)

        lpsi0 = e2
        call MPI_Allreduce(lpsi0,e2,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)

        lpsi0 = e3
        call MPI_Allreduce(lpsi0,e3,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)

        lpsi0 = e4
        call MPI_Allreduce(lpsi0,e4,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)

        lpsi0 = e5
        call MPI_Allreduce(lpsi0,e5,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)

        lpsi0 = e6
        call MPI_Allreduce(lpsi0,e6,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
        gserr = e1 + e2 + e3

        e1 = sqrt(e1)/(nxd-2)/(nyd-2)
        e2 = sqrt(e2)/(nxd-2)/(nyd-2)
        e3 = sqrt(e3)/(nxd-2)/(nyd-2)
        gserr = sqrt(gserr)/(nxd-2)/(nyd-2)

        oherr = e4 + e5 + e6

        e4 = sqrt(e4)/(nxd-2)/(nyd-2)
        e5 = sqrt(e5)/(nxd-2)/(nyd-2)
        e6 = sqrt(e6)/(nxd-2)/(nyd-2)
        oherr = sqrt(oherr)/(nxd-2)/(nyd-2)

        if (my_rank == 0) then
          write (*,*) 'EOM force residual comp 1=',e1
          write (*,*) 'EOM force residual comp 2=',e2
          write (*,*) 'EOM force residual comp 3=',e3
          write (*,*) 'EOM force residual total =',gserr
          write (*,*)
          write (*,*) 'Ohms law residual comp 1=',e4
          write (*,*) 'Ohms law residual comp 2=',e5
          write (*,*) 'Ohms law residual comp 3=',e6
          write (*,*) 'Ohms law residual total =',oherr
        endif

        !Q-profile diagnostics (not working)
cc        q0 = mval(qprof)
cc#if defined (petsc)
cc        lpsi0 = q0
cc        call MPI_Allreduce(lpsi0,q0,1,MPI_DOUBLE_PRECISION
cc     .                    ,MPI_MIN,MPI_COMM_WORLD,mpierr)
cc#endif
cc        qa = maxval(qprof)
cc#if defined (petsc)
cc        lpsi0 = qa
cc        call MPI_Allreduce(lpsi0,qa,1,MPI_DOUBLE_PRECISION
cc     .                    ,MPI_MAX,MPI_COMM_WORLD,mpierr)
cc#endif

        !Poloidal beta (2*p(0)/Bth(r=a)^2)
        I0 = 0d0
        do j=1,ny
          jglobal = j + grid_params%jlo(igrid) - 1
          bth_avg = (psig(nxd+1,j)-psig(nxd,j))/dr(iglobal)
     .              /RR(iglobal,jglobal)
          I0 = I0 + bth_avg*xmax*dth(jglobal)
        enddo

#if defined (petsc)
        lpsi0 = I0
        call MPI_Allreduce(lpsi0,I0,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_POL,mpierr)
#endif

        bth_avg = I0/2/pi/xmax

        b_pol = 2*p00*pprof(0,0)/bth_avg**2

        !Toroidal beta
#if defined (petsc)
        lpsi0 = bz2_avg
        call MPI_Allreduce(lpsi0,bz2_avg,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
        lpsi0 = p_avg
        call MPI_Allreduce(lpsi0,p_avg,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

        b_tor = 2*p_avg/bz2_avg

        !Toroidal flux
#if defined (petsc)
        lpsi0 = tor_flx
        call MPI_Allreduce(lpsi0,tor_flx,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

        bz_avg = 0.5*tor_flx/pi/xmax**2

        if (my_rank == 0) then
          write (*,*) 
          write (*,999) ' Toroidal Ohmic equilibrium features:'
          write (*,999) '   Theta   =',bth_avg/bz_avg
          write (*,999) '   Bz(r=0) =',1d0
          write (*,999) '   I0      =',I0
          write (*,999) '   E0/eta  =',Eoeta/R0
          write (*,999) '   Tor_flx =',tor_flx
cc          write (*,999) '   q(psi=0)=',q0
cc          write (*,999) '   q(psi=1)=',qa
          write (*,999) '   beta_pol=',b_pol
          write (*,999) '   beta_tor=',b_tor
 999      format (a,f10.7)
        endif

        !Save some equilibrium properties (to calculate bootstrap current)
cc        eq_params(1) = 1d0/grid_params%params(3)       !Aspect ratio
cc        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
cc        eq_params(3) = q0h(nxd)/grid_params%params(3)  !q(r=a)

        !Dump equilibrium for eigen code
cc        if (my_rank == 0) then
cc          open(unit=123,file='ntm2d.txt',status='unknown')
cc          write (123,*) nxdp+1
cc          write (123,*) grid_params%xg(0:nxdp)
cc          write (123,*) bth0
cc          write (123,*) bzz0
cc          write (123,*) 0.5*p0
cc          write (123,*) grid_params%params(3),eq_params(2),eq_params(3)
cc          close(123)
cccc          stop
cc        endif

        deallocate(psi_rhs,psiv,errv)

cc      contains
cc
ccc     flux_average
ccc     ##############################################################
cc      subroutine flux_average(psi,flux_avg)
cc
ccc     --------------------------------------------------------------
ccc     Performs flux averages according to psi and flux_avg routine.
ccc     --------------------------------------------------------------
cc
ccc     Call variables
cc
cc        real(8) :: psi(0:nxdp,0:nydp,0:nzdp)
cc
cc        INTERFACE
cc          function flux_avg(fwd,xo,x,yo,y,zo,z,so,ss)
cc          logical :: flux_avg,fwd
cc          real(8) :: xo,x,yo,y,zo,z,so,ss
cc          end function flux_avg
cc        END INTERFACE
cc
ccc     Local variables
cc
cc        integer :: ip,im,jp,jm,kp,km
cc        real(8) :: th
cc
cc        real(8),allocatable,dimension(:,:,:) :: bx,by,bz
cc     .                                         ,bxcar,bycar,bzcar,jac
cc        real(8),allocatable,dimension(:,:,:,:) :: xmap
cc
ccc     Begin program
cc
ccc     Initialize orbit integrator work space
cc
cc        allocate(bx   (0:nxd+1,0:nyd+1,0:nzd+1)
cc     .          ,by   (0:nxd+1,0:nyd+1,0:nzd+1)
cc     .          ,bz   (0:nxd+1,0:nyd+1,0:nzd+1)
cc     .          ,bxcar(0:nxd+1,0:nyd+1,0:nzd+1)
cc     .          ,bycar(0:nxd+1,0:nyd+1,0:nzd+1)
cc     .          ,bzcar(0:nxd+1,0:nyd+1,0:nzd+1)
cc     .          ,jac  (0:nxd+1,0:nyd+1,0:nzd+1)
cc     .          ,xmap (0:nxd+1,0:nyd+1,0:nzd+1,3))
cc
ccc diag ****
cccc        open(unit=110,file='debug.bin',form='unformatted'
cccc     .      ,status='replace')
cccc        call contour(psi(0:nxd+1,0:nyd,1),nxd+2,nyd+1
cccc     .              ,0d0,xmax,0d0,ymax,0,110)
cccc        close(110)
cccccc        stop
ccc diag ****
cc
ccc     Orbit integrator setup: cnv components (with BCs)
cc
cc        do k=0,nzd+1
cc          do j=0,nyd+1
cc            do i=0,nxd+1
cc              rh = grid_params%xg(i)
cc
cc              ip = min(i+1,nxd+1)
cc              im = max(i-1,0)
cc              jp = min(j+1,nyd+1)
cc              jm = max(j-1,0)
cc              kp = min(k+1,nzd+1)
cc              km = max(k-1,0)
cc
cc              !Cell center cnv components
cc              bx(i,j,k) =-(psi(i,jp,k)-psi(i,jm,k))
cc     .                   /(grid_params%yg(jp)-grid_params%yg(jm))
cc              by(i,j,k) = (psi(ip,j,k)-psi(im,j,k))
cc     .                   /(grid_params%xg(ip)-grid_params%xg(im))
cc              bz(i,j,k) = rh*0.5*(bzz0(i,j)+bzz0(im,j))
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc        !Periodic BCs
cc        bx(:,0    ,:)=bx(:,nyd,:)
cc        bx(:,nyd+1,:)=bx(:,1  ,:)
cc        by(:,0    ,:)=by(:,nyd,:)
cc        by(:,nyd+1,:)=by(:,1  ,:)
cc        bz(:,0    ,:)=bz(:,nyd,:)
cc        bz(:,nyd+1,:)=bz(:,1  ,:)
cc
cc        bx(:,:,0    )=bx(:,:,nzd)
cc        bx(:,:,nzd+1)=bx(:,:,1  )
cc        by(:,:,0    )=by(:,:,nzd)
cc        by(:,:,nzd+1)=by(:,:,1  )
cc        bz(:,:,0    )=bz(:,:,nzd)
cc        bz(:,:,nzd+1)=bz(:,:,1  )
cc
cc        !SP BCs
cc        do k=0,nzd+1
cc          do j=0,nyd+1
cc            jj = mod(j+nyd/2,nyd)
cc            if (jj == 0) jj = nyd
cc
cc            bx(0,j,k) = bx(1,jj,k)
cc            by(0,j,k) =-by(1,jj,k)
cc            bz(0,j,k) =-bz(1,jj,k)
cc          enddo
cc        enddo
cc
cc        !Cartesian components, maps
cc        do k=0,nzd+1
cc          do j=0,nyd+1
cc            do i=0,nxd+1
cc              rh = grid_params%xg(i)
cc
cc              !Cartesian components
cc              th = grid_params%yg(j)
cc
cc              bxcar(i,j,k) = bx(i,j,k)/rh*cos(th)
cc     .                     -(by(i,j,k)+kk*bz(i,j,k))*sin(th)
cc              bycar(i,j,k) = bx(i,j,k)/rh*sin(th)
cc     .                     +(by(i,j,k)+kk*bz(i,j,k))*cos(th)
cc              bzcar(i,j,k) = bz(i,j,k)/rh
cc
cc              !Map (cylindrical)
cc              jac(i,j,k)   = rh
cc              xmap(i,j,k,1)= rh*cos(th)
cc              xmap(i,j,k,2)= rh*sin(th)
cc              xmap(i,j,k,3)= grid_params%zg(k)
cc            enddo
cc          enddo
cc        enddo
cc
ccc     Orbit setup (input Az directly)
cc
cc        bx = 0d0
cc        by = 0d0
cc        bz = -psi
cc        call orbit_setup(nxd+2,nyd+2,nzd+2
cc     .                  ,grid_params%xg
cc     .                  ,grid_params%yg
cc     .                  ,grid_params%zg
cc     .                  ,bx,by,bz,bxcar,bycar,bzcar,jac,xmap
cc     .                  ,dtime=1d-1,solen=.false.,B_input_is_A=.true.)
cc
cc        !Orbit average (at cell centers)
cc        if (my_rank == 0) then
cc          if (ilevel == 1) then
cc            write (*,FMT="(a)",ADVANCE="NO") 'Orbit averaging: '
cc          elseif (ilevel == 0) then
cc            write (*,*) 'Orbit averaging...'
cc          endif
cc        endif
cc
ccc     Orbit integration
cc
cc        k = 1
cc        do j=1,nyd
cc          do i=1,nxd
cc
cc            x0 = grid_params%xg(i)
cc            y0 = grid_params%yg(j)
cc            z0 = grid_params%zg(k)
cc
cc            bzz_avg = 0d0
cc            b2_avg  = 0d0
cc            ssmax   = 0d0
cc
cc            call orbit_find(ilevel-2,x0,y0,z0,.true.
cc     .                     ,line_int=flux_avg)
cc
cc            if ((my_rank == 0).and.(ilevel == 1))
cc     .                  write (*,FMT="(a)",ADVANCE="NO") '+'
cc
cc            iR2psi(i,j) = iR2
cc            B2psi (i,j) = b2_avg
cc
cc          enddo
cc        enddo
cc
cc        B2psi(:,0    )=B2psi(:,nyd)
cc        B2psi(:,nyd+1)=B2psi(:,1  )
cc
cc        iR2psi(:,0    )=iR2psi(:,nyd)
cc        iR2psi(:,nyd+1)=iR2psi(:,1  )
cc
ccc diag ****
cccc          open(unit=110,file='debug.bin',form='unformatted'
cccc     .        ,status='replace')
cccc          call contour(psi(0:nxd+1,0:nyd,1),nxd+2,nyd+1
cccc     .                ,0d0,xmax,0d0,ymax,0,110)
cccc          close(110)
cccccc          stop
cccc
cccccc          write (*,*) 'here'
cccc          open(unit=110,file='debug2.bin',form='unformatted'
cccc     .        ,status='replace')
cccc          call contour(ffpsi(1:nxd,1:nyd),nxd,nyd
cccc     .                ,0d0,xmax,0d0,ymax,0,110)
cccc          close(110)
cccc          stop
ccc diag ****
cc
ccc     Free up orbit integral work space
cc
cc        if ((my_rank == 0).and.(ilevel == 1)) write (*,*)
cc
cc        call orbit_destroy
cc
cc        deallocate(bx,by,bz,bxcar,bycar,bzcar,jac,xmap)
cc
cc      end subroutine flux_average

      end subroutine Ohm_GS_equil_2d

c     qprof
c     ###############################################################
      function qprof(r)

c     ---------------------------------------------------------------
c     Finds cylindrical components of magnetic field B_theta, B_z
c     from q profile.
c     ---------------------------------------------------------------

        use equilibrium

        use app_iosetup

        implicit none

c     Call variables

        real(8)    :: qprof,r

c     Local variables

        integer :: p
        real(8) :: q0,qa

c     Begin program

        select case(equil)
        case('ntm2d')

          q0 = 1.02
          qa = 4.
          p  = 1

          qprof = q0*(1+r**(2*p)*( (qa/q0)**p - 1 ))**(1./p)

        case('rfp83')

          qprof = 0.6125*(1 - 1.8748*r**2 + 0.8323*r**4)

        case default

          call pstop('qprof','Undefined q profile')

        end select

c     End program

      end function qprof

c     pprof
c     ###############################################################
      function pprof(s,der)

c     ---------------------------------------------------------------
c     Specifies pressure profile (der = 0) and its derivative (der=1)
c     as a function of poloidal flux.
c     ---------------------------------------------------------------

        use equilibrium

        implicit none

c     Call variables

        integer    :: der
        real(8)    :: pprof,s

c     Local variables

c     Begin program

        select case(equil)
        case('ntm2d','ppnsl')

          if (der == 0) then
            pprof = (1-s**2)**4
          else
            pprof =-8*s*(1-s**2)**3
          endif

        case('ohtor')

          if (der == 0) then
            pprof = (1-s**2)**2
          else
            pprof =-4*s*(1-s**2)
          endif

        case default

          pprof = 0d0

        end select

c     End program

      end function pprof

c fillVectorPotential
c####################################################################
      subroutine fillVectorPotential(a,iigx,iigy,iigz,bc)

c--------------------------------------------------------------------
c     Defines COVARIANT vector potential for initialization of 
c     equilibrium.
c--------------------------------------------------------------------

      use parameters

      use grid

      use equilibrium

      use imposeBCinterface

      implicit none

c Call variables

      integer    :: iigx,iigy,iigz
      real(8)    :: a(ilom:ihip,jlom:jhip,klom:khip,3)
      logical    :: bc

c Local variables

      integer    :: ig,jg,kg,i,j,k,bcnd(6,3),nx,ny,nz
      real(8)    :: x1,y1,z1,eps,car(3)
      real(8)    :: a_cnv(ilom:ihip,jlom:jhip,klom:khip,3)
      logical    :: cartsn

c Begin program

      nx = grid_params%nxv(iigx)
      ny = grid_params%nyv(iigy)
      nz = grid_params%nzv(iigz)

      select case (trim(equil))

      case ('khcar')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,iigx,iigy,iigz,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = x1
              a(i,j,k,3) = 0d0
            enddo
          enddo
        enddo

      case ('gem','kai1')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,iigx,iigy,iigz,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = -dlambda
     .             *dlog(cosh((x1-0.5d0*(xmax+xmin))/dlambda)) 

            enddo
          enddo
        enddo

      case ('tmcar','tmsin')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,iigx,iigy,iigz,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = 2*dlambda
     .               *atan(exp((x1-0.5d0*(xmax+xmin))/dlambda))
              a(i,j,k,3) = -dlambda
     .               *log(cosh((x1-0.5d0*(xmax+xmin))/dlambda)) 
            enddo
          enddo
        enddo

      case ('ic')

        eps = grid_params%params(2)
        dlambda = (ymax-ymin)/2./pi/2.

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,iigx,iigy,iigz,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = -dlambda*dlog(    cosh(x1/dlambda)
     $                                  +eps*cos (y1/dlambda)) 
            enddo
          enddo
        enddo

      case ('tmcyz','3dkh')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,iigx,iigy,iigz,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = dlambda
     .             *dlog(cosh((y1-0.5d0*(ymax+ymin))/dlambda)) 
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = 0d0
            enddo
          enddo
        enddo

      end select

c Set boundary values to face values for EQU BC

      if (bc) then
        bcnd(:,1) = bcond
        bcnd(:,2) = bcond
        bcnd(:,3) = bcond
        where (bcnd == DEF) bcnd = -IFC

        call setBC(IAX,3,nx,ny,nz,a_cnv,a,a,bcnd,iigx,iigy,iigz
     .          ,is_cnv=.false.)
      endif

c End program

      end subroutine fillVectorPotential

c findVectorPotential
c####################################################################
      subroutine findVectorPotential(v,b,a_cov)

c--------------------------------------------------------------------
c     Finds COVARIANT vector potential from CONTRAVARIANT magnetic
c     field components with gauge:
c                grad(phi + eta div(a)) = 0
c
c     On input:
c       * v -> cnv velocity
c       * b -> cnv magnetic field
c     On output:
c       * a -> cov vector potential
c--------------------------------------------------------------------

      use parameters

      use grid_aliases

      use equilibrium

cc      use precond_variables, ONLY: poisson
cc
cc      use  mg_internal, ONLY: mapMGVectorToArray

      use operators

      use imposeBCinterface

      use timeStepping

      implicit none

c Call variables

      real(8)    :: a_cov(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,v(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,b(0:nxdp,0:nydp,0:nzdp,3)

c Local variables

      integer    :: i,j,k,igrid,ii,bcnd(6,3),nx,ny,nz
      real(8)    :: cm,c0,cp,dxx,dyy,dvol,cov(3)
     .             ,cai_rhs(ntotdp),caiv(ntotdp)
     .             ,phi(ilom:ihip,jlom:jhip,klom:khip)
     .             ,cai(ilom:ihip,jlom:jhip,klom:khip)
     .             ,a_cnv(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,a(0:nxdp,0:nydp,0:nzdp,3)

      logical    :: spoint

c Begin program

      igrid = 1

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

c Find equilibrium (d_t A = 0) scalar potential (vxb electric field)

cc      call findScalarPotential(v
cc     .                        ,b(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)
cc     .                        ,phi)

c Find vector potential with gauge Ax = 0 in GLOBAL domain

      call curl_inv(nxd,nyd,nzd,igrid,igrid,igrid,b,a,global=.true.)

c Restrict solution to current domain

      a_cov(ilom:ihip,jlom:jhip,klom:khip,:)
     .     = a(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)

ccc Impose **topological constraints** and find a_cnv
cc
cc      bcnd(:,1) = bcond
cc      bcnd(:,2) = bcond
cc      bcnd(:,3) = bcond
cccc      where (bcnd == DEF) bcnd = EXT
cc
cc      call setBC(IAX,3,nx,ny,nz,a_cnv,a_cov,vzeros,bcnd
cc     .          ,igrid,igrid,igrid,is_cnv=.false.,iorder=2)
cc
ccc Transform A to gauge div(A)=0
cc
cc      !Find gauge transformation potential cai
cc      do k=klo,khi
cc        do j=jlo,jhi
cc          do i=ilo,ihi
cc            ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc            if (vol_wgt) then
cc              dvol = gmetric%grid(igrid)%dvol(i,j,k)
cc            else
cc              dvol = 1d0
cc            endif
cc
cc            cai_rhs(ii)=-dvol*(div(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                            ,a_cnv(:,:,:,1)
cc     .                            ,a_cnv(:,:,:,2)
cc     .                            ,a_cnv(:,:,:,3)))
cc          enddo
cc        enddo
cc      enddo
cc
cc      !Transform A (if rhs /= 0d0)
cc      if (sqrt(sum(cai_rhs*cai_rhs)) > 0d0) then
cc
cc        !Solve Poisson problem to find Gauge potential cai
cc#if defined(petsc)
cc        call pstop('findVectorPotential'
cc     .          ,'Poisson solver not implemented in parallel')
cc#else
cc        write (*,*) 'Solving for Gauge potential'
cc
cc        !Define BCs
cc        bcnd(:,1) = bcond
cc        where (bcnd == DEF) bcnd = DIR
cc
cc        !Solve Poisson eq
cc        call poisson(ntotdp,cai_rhs,caiv,30,1d-6,bcnd(:,1),ilevel,0)
cc
cc        !Map to array (w BCs)
cc        call mapMGVectorToArray(0,1,caiv,nx,ny,nz,cai,igrid,.false.)
cc
cc        call setBC(IRHO,nx,ny,nz,cai,zeros,bcnd(:,1),igrid,igrid,igrid
cc     .            ,iorder=2)
cc#endif
cc
cc        !Gauge transformation A' = A + grad(cai)
cc        do k=klom,khip
cc          do j=jlom,jhip
cc            do i=ilom,ihip
cc              call grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,cai
cc     .                 ,cov(1),cov(2),cov(3))
cc              a_cov(i,j,k,:) = a_cov(i,j,k,:) + cov
cc            enddo
cc          enddo
cc        enddo
cc
cc      endif

c Set boundary values to face values for EQU BC

      bcnd(:,1) = bcond
      bcnd(:,2) = bcond
      bcnd(:,3) = bcond
      where (bcnd == DEF) bcnd = -IFC
cc      where (bcnd == DEF) bcnd = -EXT

      call setBC(IAX,3,nx,ny,nz,a_cnv,a_cov,vzeros,bcnd
     .          ,igrid,igrid,igrid,is_cnv=.false.)

c End program

      end subroutine findVectorPotential

ccc findScalarPotential
ccc####################################################################
cc      subroutine findScalarPotential(v,b,phi)
cc
ccc--------------------------------------------------------------------
ccc     Finds scalar equilibrium potential from equilibrium velocity and 
ccc     magnetic fields, obeying the relation:
ccc     
ccc         lap(phi_0) = div(v0 x B0)
ccc
ccc     This automatically takes into account ExB flows at the boundary.
ccc--------------------------------------------------------------------
cc
cc      use parameters
cc
cc      use grid
cc
cc      use equilibrium
cc
cc      use precond_variables
cc
cc      implicit none
cc
ccc Call variables
cc
cc      real(8)    :: v(ilom:ihip,jlom:jhip,klom:khip,3)
cc     .             ,b(ilom:ihip,jlom:jhip,klom:khip,3)
cc     .             ,phi(ilom:ihip,jlom:jhip,klom:khip)
cc
ccc Local variables
cc
cc      integer    :: igrid,ii,bcnd(6,1)
cc      real(8)    :: phi_rhs(ntotdp),phiv(ntotdp),cov(3)
cc     .             ,dum(ilom:ihip,jlom:jhip,klom:khip,3),dvol
cc
ccc Begin program
cc
cc      igrid = 1
cc
cc      nx = grid_params%nxv(igrid)
cc      ny = grid_params%nyv(igrid)
cc      nz = grid_params%nzv(igrid)
cc
ccc Find rhs
cc
cc      do k=klom,khip
cc        do j=jlom,jhip
cc          do i=ilom,ihip
cc            jac = gmetric%grid(igrid)%jac(i,j,k)
cc
cc            cov(1) = v(i,j,k,2)*b(i,j,k,3)
cc     .             - v(i,j,k,3)*b(i,j,k,2)
cc
cc            cov(2) = v(i,j,k,3)*b(i,j,k,1)
cc     .             - v(i,j,k,1)*b(i,j,k,3)
cc
cc            cov(3) = v(i,j,k,1)*b(i,j,k,2)
cc     .             - v(i,j,k,2)*b(i,j,k,1)
cc
cc            cov = cov/jac
cc
cc            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
cc     .                                  ,cov(1),cov(2),cov(3)
cc     .                                  ,dum(i,j,k,1)
cc     .                                  ,dum(i,j,k,2)
cc     .                                  ,dum(i,j,k,3),.true.)
cc
cc          enddo
cc        enddo
cc      enddo
cc
cc      do k=klo,khi
cc        do j=jlo,jhi
cc          do i=ilo,ihi
cc            ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc            if (vol_wgt) then
cc              dvol = volume(i,j,k,igrid,igrid,igrid)
cc            else
cc              dvol = 1d0
cc            endif
cc
cc            phi_rhs(ii) = dvol*div(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                        ,dum(:,:,:,1),dum(:,:,:,2),dum(:,:,:,3))
cc          enddo
cc        enddo
cc      enddo
cc
ccc Solve for potential
cc
cc      if (sqrt(sum(phi_rhs*phi_rhs))>0d0) then
cc
cc#if defined(petsc)
cc
cc        call pstop('findScalarPotential'
cc     .          ,'Poisson solver not implemented in parallel')
cc
cc#else
cc
cc        !Define BCs
cc        bcnd(:,1) = bcond
cc        where (bcnd == DEF) bcnd = DIR
cc
cc        !Solve Poisson eq
cc        call poisson(ntotdp,phi_rhs,phiv,30,1d-6,bcnd,ilevel,0)
cc
cc        !Map to array (w/o BCs)
cc        call mapMGVectorToArray(0,1,phiv,nx,ny,nz,phi,igrid,.false.)
cc#endif
cc
cc      else
cc
cc        phi = 0d0
cc
cc      endif
cc
ccc End program
cc
cc      end subroutine findScalarPotential
