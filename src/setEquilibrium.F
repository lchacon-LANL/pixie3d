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

      use imposeBCinterface

      use app_iosetup

      implicit none

c Call variables

      integer       :: iigx,iigy,iigz

      real(8)       :: var(ilom:ihip,jlom:jhip,klom:khip,neqd)

      character(20) :: label(neqd)

      integer       :: bcs(6,neqd)

c Local variables

      integer       :: i,j,k,ig,jg,kg,ieq,nmax,iglobal,ierr,sh,it
     .                ,icomp,bcnd(6,3),brank

      real(8)       :: ldaflow,x1,y1,z1,cx,cy,cz

      logical       :: covariant,to_cartsn,cartsn,a_bc

      real(8)       :: r,jac1,bnorm,aaa,bbb,ccc,qq,qqp,q0,RR,ff,nn
     .                ,aspect_ratio,mm,kk,Iz,bz0,btheta,vr,rrho,r1,bb
     .                ,r0(1),aa,dummy1,dummy2,dummy3,ieta,err,bz_avg,dr
     .                ,bzh,bth,rh,e0_over_eta,theta_ppnch,p0,eps
     .                ,phi,rho_norm,beta_tor

      real(8)       :: rint(0:nxdp),dh,offset(1),bth0(0:nxdp)
     .                ,bzz0(0:nxdp),vr0(0:nxdp),rho0(0:nxdp)
     .                ,eth0(0:nxdp),ez0(0:nxdp),jth0(0:nxdp)
     .                ,jz0 (0:nxdp),bth0_old(0:nxdp),p00(0:nxdp)
     .                ,bzz0_old(0:nxdp),dummy(0:nxdp),bmag2(0:nxdp)
     .                ,q00(0:nxdp),dum (ilom:ihip,2)
     .                ,ldum(ilom:ihip,2)
     .                ,r2  (ilom:ihip,jlom:jhip)

      real(8)       :: a(ilom:ihip,jlom:jhip,klom:khip,3)

      character(200):: command

#if defined(vec_pot)
      real(8)       :: b(0:nxdp,0:nydp,0:nzdp,3)
#endif

c Functions

      real(8) :: qprof,pprof
      external   qprof,pprof

      !q-profile
      qq (r1) = 0.6125*(1 - 1.8748*r1**2 + 0.8323*r1**4)
      qqp(r1) = 0.6125*( -2*1.8748*r1   +4*0.8323*r1**3)
cc      qq(r1)  = 0.3*(1 - 1.8748*r1**2 + 0.8323*r1**4)
cc      qqp(r1) = 0.3*( -2*1.8748*r1   +4*0.8323*r1**3)

      ff(r1) = r1**2 + qq(r1)**2

c Begin program

      ldaflow = dlambda/rshear

      var = 0d0

      jcov_0 = 0d0

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
      label(ITMP) = 'Temp'

c Define boundary conditions

      call defineBoundaryConditions(neqd,bcs)

c Set initial guess

      select case (trim(equil))

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

cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,x1,y1,z1)
cc
cc              var(i,j,k,IAX)  = 0d0
cc              var(i,j,k,IAY)  = x1
cc              var(i,j,k,IAZ)  = 0d0
cc            enddo
cc          enddo
cc        enddo
#endif

      case ('3wv') !3-wave MHD test

        gamma = 5./3.

        phi   = grid_params%params(1)   !Angle of B with x-axis, in degrees
        beta  = grid_params%params(2)   !plasma beta

        phi   = phi/180.*pi             !Phi in radians

        nh2   = nh1                     !This is required by setup

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

        var(:,:,:,ITMP) = bb/2.

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

        var(:,:,:,ITMP) = 0.25

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

cc      case ('tmcxz')
cc
ccc     Define vector potential (in curvilinear coordinates) for initialization
cc
cc        call fillVectorPotential(a,igx,igy,igz,a_bc)
cc
ccc     Check coordinates
cc
cc        if (coords /= 'car') then
cc          messg = 'Wrong coordinates for equilibrium '//equil
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
ccc     Tearing mode in cartesian coordinates
cc
cc        bz0 = 1d0
cc
cc        do k = klo,khi
cc          do j = jlo,jhi
cc            do i = ilo,ihi
cc             !X-Z equilibrium
cc             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,1)
cc             var(i,j,k,IBZ)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,2)
cc             var(i,j,k,IBY)=sqrt(bz0**2 - var(i,j,k,IBZ)**2)
cc            enddo
cc          enddo
cc        enddo
cc
cc        var(:,:,:,IVX)  = 0d0
cc        var(:,:,:,IVY)  = 0d0
cc        var(:,:,:,IVZ)  = 0d0
cc
cc        var(:,:,:,IRHO) = 1d0
cc
cc        var(:,:,:,ITMP) = 1d0
cc
cc      case ('tmcyz')
cc
ccc     Define vector potential (in curvilinear coordinates) for initialization
cc
cc        call fillVectorPotential(a,igx,igy,igz,a_bc)
cc
ccc     Check coordinates
cc
cc        if (coords /= 'car') then
cc          messg = 'Wrong coordinates for equilibrium '//equil
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
ccc     Tearing mode in cartesian coordinates
cc
cc        bz0 = 1d0
cc
cc        do k = klo,khi
cc          do j = jlo,jhi
cc            do i = ilo,ihi
cc             !X-Z equilibrium
cc             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,2)
cc             var(i,j,k,IBZ)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,3)
cc             var(i,j,k,IBX)=sqrt(bz0**2 - var(i,j,k,IBZ)**2)
cc            enddo
cc          enddo
cc        enddo
cc
cc        var(:,:,:,IVX)  = 0d0
cc        var(:,:,:,IVY)  = 0d0
cc        var(:,:,:,IVZ)  = 0d0
cc
cc        var(:,:,:,IRHO) = 1d0
cc
cc        var(:,:,:,ITMP) = 1d0

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

      case ('kaitm')

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

      case ('kai3d')

        mm = grid_params%params(1)
        kk = grid_params%params(2)
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

        nh2 = mm       !To set the right perturbation wavelength
        nh3 = kk*RR    !To set the right perturbation wavelength

        nh1 = nh2      !For perturbations to be of the right order

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

      case ('ppnch','ppnsl','ppnst')

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

cc        command = '[ -f '//trim(equ_file)//' ]'

cc        ierr = system(trim(command))

cc        if (ierr == 0) then     !External file equilibrium
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
          read(1000,'(f7.1)') ieta

          close(1000)

          if (my_rank ==0) write (*,*) 'Done!'

          !Find electric field at face r=1d0
          E0(2) = 0.5*(eth0(nxd)+eth0(nxdp))
          E0(3) = 0.5*(ez0 (nxd)+ez0 (nxdp))
          E0 = E0*eta/ieta

        else   !Calculate equilibrium

          if (my_rank == 0)
     .       write (*,*) 'Calculating equilibrium ',equil,'...'

          do i=0,nxdp
            x1 = i*xmax/nxd
            bth0(i) = 0.5*x1*e0_over_eta
          enddo

	  !Initialize iteration
          bzz0 = 1d0

	  !Iteration
          dr = xmax/nxd

          do it=1,100

            bth0_old = bth0
            bzz0_old = bzz0

	    !Perform integral of B_theta
            dummy(0) = 0d0
            bth0 (0) = 0d0
            do i=1,nxdp
              x1 = i*dr
              bzh = 0.5*(bzz0(i)+bzz0(i-1))
              bth = 0.5*(bth0(i)+bth0(i-1))
              rh  = 0.5*(2*i-1)*dr

              dummy(i) = dummy(i-1) + bzh**2/(bzh**2+bth**2)*rh*dr

              bth0(i) = e0_over_eta*dummy(i)/x1
            enddo

            !Perform integral of Bz
            dummy(0) = 0d0
            do i=1,nxdp
              x1 = i*xmax/nxd
              bzh = 0.5*(bzz0(i)+bzz0(i-1))
              bth = 0.5*(bth0(i)+bth0(i-1))
              rh  = 0.5*(2*i-1)*dr
              dummy(i) = dummy(i-1) + bth/(bzh**2+bth**2)*dr
            enddo

            bzz0 = exp(-e0_over_eta*dummy)

            if (theta_ppnch /= 0d0) then
              bz_avg = 0d0
              do i=1,nxd
                 bz_avg = bz_avg + (bzz0(i)*i + bzz0(i-1)*(i-1))*dr**2
              enddo
              bz_avg = bz_avg/xmax**2
              bzz0 = bth0(nxd)/theta_ppnch/bz_avg*bzz0
            elseif (bz0 /= 0d0) then
              bzz0 = bz0*bzz0
            else
              call pstop('setEquilibrium'
     $                  ,'Can not find PPNCH equilibrium:'
     $                   //'Both Bz0 and theta are zero')
            endif

	    !Check convergence
            err = sqrt(0.5*(sum((bth0-bth0_old)**2)
     .                     +sum((bzz0-bzz0_old)**2))/nxdp)

            if (ilevel > 0 .and. my_rank == 0) then
              write (*,*) 'Equilibrium iter =',it,' ; Error=',err
            endif

            if (err < 0.1*dr**2) exit

          enddo

          if (ilevel == 0 .and. my_rank == 0) then
             write (*,*)
             write (*,'(a,i3,a,1pe10.2)') ' Equilibrium converged in '
     .            ,it,' iterations with error =',err
          endif

          !Check force free diagnostic
          do i=1,nxd
            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))
            rh  = 0.5*(2*i-1)*dr
            dummy(i) =  bzh*(bzz0(i)-bzz0(i-1))/dr
     $                 +bth/rh*(bth0(i)*i-bth0(i-1)*(i-1))
          enddo

          if (my_rank == 0) then
            write(*,*)
            write(*,'(a,1pe14.7)') ' Force-free error ='
     .                          ,sqrt(sum(dummy(1:nxd)**2)/nxd)
          endif

          !Check Ohm's law
          do i=1,nxd
            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))
            rh  = 0.5*(2*i-1)*dr
            dummy(i) = e0_over_eta*bzh
     $                -bzh/rh*(bth0(i)*i-bth0(i-1)*(i-1))
     $                +bth*(bzz0(i)-bzz0(i-1))/dr
          enddo

          if (my_rank == 0) then
            write(*,'(a,1pe14.7)') ' Ohms law   error ='
     .                          ,sqrt(sum(dummy(1:nxd)**2)/nxd)
          endif

          !Calculate pinch equilibrium features (theta parameter, Bz(r=0), I0)
          bz_avg = 0d0
          do i=1,nxd
            x1 = i*xmax/nxd
            bz_avg = bz_avg + (bzz0(i)*i + bzz0(i-1)*(i-1))*dr**2
          enddo

          if (my_rank == 0) then
            write (*,*) 
            write (*,999) ' Pinch equilibrium features:'
            write (*,999) '   Theta  =',bth0(nxd)/bz_avg
            write (*,999) '   Bz(r=0)=',bzz0(0)
            write (*,999) '   I0     =',2*pi*xmax*bth0(nxd)
            write (*,999) '   E0/eta =',e0_over_eta
            write (*,999) '   T_flux =',pi*xmax**2*bz_avg
 999        format (a,f10.7)
          endif

	  !Transfer to integer grid
          bth0(1:nxdp) = 0.5*(bth0(0:nxd)+bth0(1:nxdp))
          bzz0(1:nxdp) = 0.5*(bzz0(0:nxd)+bzz0(1:nxdp))

	  !Find electric field
          E0(1:2) = 0d0
          E0(3) = e0_over_eta*eta
          ieta = eta

          !Find radial flow
          vr0 = -(E0(3)*bth0-E0(2)*bzz0)/(bzz0**2+bth0**2)

	  !Find density
          rho0 = 1d0

        endif

c$$$        write (*,*) 'r,bth0,bzz0,jth0,jz0,vr0'
c$$$        do i=0,nxdp
c$$$          call getCurvilinearCoordinates(i,1,1,igx,igy,igz,ig,jg,kg
c$$$     .                                  ,x1,y1,z1)
c$$$          write(*,*) x1,bth0(i),bzz0(i),jth0(i),jz0(i),vr0(i)
c$$$        enddo

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
              var(i,j,k,IVX) = x1*vr0(iglobal)*eta/ieta
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
cc          var(ihi+1,:,:,IVX) = 0.5*(var(ihi,:,:,IVX)+var(ihi+1,:,:,IVX))
          var(ihi+1,:,:,IVX) = 0.25*(var(ihi  ,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi,:,:)
     .                              +var(ihi+1,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi+1,:,:))
     .                         *(gmetric%grid(igx)%jac(ihi  ,:,:)
     .                          +gmetric%grid(igx)%jac(ihi+1,:,:))
        endif

      case ('ppn3d')

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)

        nh2 = mm       !To set the right perturbation wavelength
        nh3 = kk*RR    !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

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

        write (*,*) 'Reading numerical equilibrium ',equil,'...'

        !Read file
        open (unit=1000,file=trim(equ_file),status='old')

        read(1000,'(a)')
        do i=0,nxdp
          read(1000,'(9f12.8)') dummy1,bth0(i),bzz0(i),jth0(i),jz0(i)
     .                ,vr0(i),ez0(i),eth0(i),rho0(i)
        enddo

        read(1000)
        read(1000,'(a)')
        read(1000,*) ieta

        close(1000)

        write (*,*) 'Done!'

        !Find electric field at face r=1d0
        E0(2) = 0.5*(eth0(nxd)+eth0(nxdp))
        E0(3) = 0.5*(ez0 (nxd)+ez0 (nxdp))
        E0 = E0*eta/ieta

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
cc          var(ihi+1,:,:,IVX) = 0.5*(var(ihi,:,:,IVX)+var(ihi+1,:,:,IVX))
          var(ihi+1,:,:,IVX) = 0.25*(var(ihi  ,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi,:,:)
     .                              +var(ihi+1,:,:,IVX)
     .                              /gmetric%grid(igx)%jac(ihi+1,:,:))
     .                         *(gmetric%grid(igx)%jac(ihi  ,:,:)
     .                          +gmetric%grid(igx)%jac(ihi+1,:,:))
        endif

      case ('rfp1')      !RFP EQUILIBRIA (Caramana et al, PoP, 1983)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     RFP equilibria (Caramana et al, PoP, 1983)

        !Integral in r
        rint(nxdp) = 0d0  !First guess to rint(r=1d0) = 0d0
        do i=nxd,1,-1
          x1 = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))  !Position in global grid
          dh = grid_params%xg(i+1)-grid_params%xg(i)
          rint(i) = rint(i+1) + dh*qq(x1)*qqp(x1)/ff(x1)
        enddo
        rint(0) = rint(1)  !Neumann BC for integral at r=0d0

        !Correction to ensure rint(r=1d0) = 0d0
        r0 = 1d0
        call IntDriver1d(nxd+2,grid_params%xg(0:nxdp)
     .                  ,rint,1,r0,offset,2)

        rint = rint - offset(1)

        !Determine Iz so that Bz(r=0)=1
        r0 = 0d0
        call IntDriver1d(nxd+2,grid_params%xg(0:nxdp)
     .                  ,rint,1,r0,offset,2)

cc        btheta = 1./2./pi*ff(1d0)/ff(0d0)*exp(-rint(0))
        btheta = 1./2./pi*ff(1d0)/ff(0d0)*exp(-offset(1))
        Iz = 1./btheta/qq(0d0)

cc        call B_from_q(bth0,bzz0)

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

              iglobal = i + grid_params%ilo(igx) - 1

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = bth0(iglobal) + kk*x1/mm*bzz0(iglobal)
              var(i,j,k,IBZ)  = x1*bzz0(iglobal)
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

      case ('rfp2','spnch')  !SCREW PINCH EQ.

c     Check coordinates

        if (coords /= 'hel' .and. coords /= 'hl2') then
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

              if (coords == 'hl2') x1 = sign(sqrt(2*abs(x1)),x1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              if (coords == 'hel') then
                var(i,j,k,IBX)  = 0d0
                var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
                var(i,j,k,IBZ)  = x1*bz0
              else
                var(i,j,k,IBX)  = 0d0
                var(i,j,k,IBY)  = btheta/x1 + kk/mm*bz0
                var(i,j,k,IBZ)  = bz0
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

      case ('rfp3d')

        mm = grid_params%params(1)
        kk = grid_params%params(2)
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

        nh2 = mm       !To set the right perturbation wavelength
        nh3 = kk*RR    !To set the right perturbation wavelength

        nh1 = nh2      !For radial perturbations to be of the right order

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

      case ('tok')

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

#if defined(vmec)

      case ('vmec')    

        nh1 = nh2       !For perturbations to be of the right order in radius
                        !(nh2,nh3 set by input file)

        if (nh2 == 0 .or. nh3 == 0) then
          messg='Perturbation harmonics for "vmec" equilibrium not set'
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

        var(:,:,:,IRHO) = var(:,:,:,ITMP)**(1d0/gamma)   !Find density from pressure

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

      case ('tor')

        RR = grid_params%params(1)
        mm = grid_params%params(4)
        nn = grid_params%params(5)
        q0 = grid_params%params(6)

        kk = 0d0   !To avoid affecting q-profile

        bb = dlambda/RR/q0

c     Check coordinates

        if (coords /= 'tor') then
          messg = 'Need "tor" coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Read VMEC equilibrium

        nh2 =  mm       !To set the right perturbation wavelength
        nh3 = -nn       !To set the right perturbation wavelength

        nh1 = nh2       !For perturbations to be of the right order

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

              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              !X-Y equilibrium
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta
              var(i,j,k,IBZ)  = gmetric%grid(igx)%gsup(i,j,k,3,3)*bz0
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
          p00(i) = 0.5*beta_tor*bz0**2*pprof(x1)/pprof(0d0)
          q00(i) = RR*qprof(x1)
        enddo

        !Initialize toroidal field
        bzz0 = bz0  

        !Solve 1D Grad-Shafranov equation
        call GS_equil_1d(q00,p00,bth0,bzz0)

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
              var(i,j,k,IBY) = bth0(iglobal) + kk*x1/nh2*bzz0(iglobal)
              var(i,j,k,IBZ) = x1*bzz0(iglobal)
#endif
              var(i,j,k,IVX) = 0d0
              var(i,j,k,IVY) = 0d0
              var(i,j,k,IVZ) = 0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,ITMP) = 0.5*p00(iglobal)/var(i,j,k,IRHO)

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

cc        !Shift vr for EQU BC
cc        if (ihip + grid_params%ilo(igx)-1 == nxdp) then
cccc          var(ihi+1,:,:,IVX) = 0.5*(var(ihi,:,:,IVX)+var(ihi+1,:,:,IVX))
cc          var(ihi+1,:,:,IVX) = 0.25*(var(ihi  ,:,:,IVX)
cc     .                              /gmetric%grid(igx)%jac(ihi,:,:)
cc     .                              +var(ihi+1,:,:,IVX)
cc     .                              /gmetric%grid(igx)%jac(ihi+1,:,:))
cc     .                         *(gmetric%grid(igx)%jac(ihi  ,:,:)
cc     .                          +gmetric%grid(igx)%jac(ihi+1,:,:))
cc        endif

      case default

        messg = 'Equilibrium '//trim(equil)//' undefined'
        call pstop('setEquilibrium',messg)

      end select

c Find momentum components

      if (.not.nc_eom_v) then
        var(:,:,:,IVX) = var(:,:,:,IRHO)*var(:,:,:,IVX)
        var(:,:,:,IVY) = var(:,:,:,IRHO)*var(:,:,:,IVY)
        var(:,:,:,IVZ) = var(:,:,:,IRHO)*var(:,:,:,IVZ)
      endif

c Find equilibrium current (-> jcov_0)

#if defined(vec_pot)
      acov = var(:,:,:,IAX:IAZ)
#else
      bcnv = var(:,:,:,IBX:IBZ)
#endif

      call findEquilJ

c Transfer electric field info to postprocessor via input file

      if (sum(E0) /= 0d0 .and. my_rank == 0) then

        command = 'grep E0 ' // trim(inputfile) // ' > /dev/null'
        ierr = system(trim(command))

        write (command,*) '   E0 =',real(E0(1),4)
     .                         ,',',real(E0(2),4)
     .                         ,',',real(E0(3),4)


        if (ierr /= 0) then  !Replace E0
          command='sed -i -e"/graphdef/{" -e' // achar(39)
#if defined(pgf90) || defined(xlf)
     .       // 'a\\' // achar(39) //' -e"'
#else
     .       // 'a\' // achar(39) //' -e"'
#endif
     .       // trim(command)   // '" -e"}" '
     .       // trim(inputfile)
        else                 !Add E0
          command='sed -i -e"/E0/{" -e' // achar(39)
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
          command='sed -i -e"/graphdef/{" -e' // achar(39)
#if defined(pgf90) || defined(xlf)
     .       // 'a\\' // achar(39) //' -e"'
#else
     .       // 'a\' // achar(39) //' -e"'
#endif
     .       // trim(command)   // '" -e"}" '
     .       // trim(inputfile)
        else                 !Add B0
          command='sed -i -e"/B0/{" -e' // achar(39)
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

      contains

c     B_from_q
c     ###############################################################
      subroutine B_from_q(bth0,bz0)

c     ---------------------------------------------------------------
c     Finds cylindrical components of magnetic field B_theta, B_z
c     from q profile.
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        real(8)    :: bth0(0:nxdp),bz0(0:nxdp)

c     Local variables

        integer    :: i,j,k,ig,jg,kg
        real(8)    :: qq,qqp,ff,x1,dh,offset(1),r1,btheta,Iz,r0(1)

        real(8)    :: rint(0:nxdp)

c     Functions

        !q-profile
        qq (r1) = 0.6125*(1 - 1.8748*r1**2 + 0.8323*r1**4)
        qqp(r1) = 0.6125*( -2*1.8748*r1   +4*0.8323*r1**3)
cc        qq(r1)  = 0.3*(1 - 1.8748*r1**2 + 0.8323*r1**4)
cc        qqp(r1) = 0.3*( -2*1.8748*r1   +4*0.8323*r1**3)

        ff(r1) = r1**2 + qq(r1)**2

c     Begin program

c     Find integrals

        !Integral in r
        rint(nxdp) = 0d0  !First guess to rint(r=1d0) = 0d0
        do i=nxd,1,-1
          x1 = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))  !Position in global grid
          dh = grid_params%xg(i+1)-grid_params%xg(i)
          rint(i) = rint(i+1) + dh*qq(x1)*qqp(x1)/ff(x1)
        enddo
        rint(0) = rint(1)  !Neumann BC for integral at r=0d0

        !Correction to ensure rint(r=1d0) = 0d0
        r0 = 1d0
        call IntDriver1d(nxd+2,grid_params%xg(0:nxdp)
     .                  ,rint,1,r0,offset,2)

        rint = rint - offset(1)

        !Determine Iz so that Bz(r=0)=1
        r0 = 0d0
        call IntDriver1d(nxd+2,grid_params%xg(0:nxdp)
     .                  ,rint,1,r0,offset,2)

cc        btheta = 1./2./pi*ff(1d0)/ff(0d0)*exp(-rint(0))
        btheta = 1./2./pi*ff(1d0)/ff(0d0)*exp(-offset(1))
        Iz = 1./btheta/qq(0d0)

c     Evaluate magnetic field

        do i = 0,nxdp

          x1 = grid_params%xg(i)  !Position in global grid

          bth0(i) = Iz/2./pi*x1*ff(1d0)/ff(x1)*exp(-rint(i))
          bz0 (i) = qq(x1)*bth0(i)/x1

        enddo

      end subroutine B_from_q

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
        real(8)    :: dr,bz_avg,x1,bzh,bth,err,rh

        real(8)    :: bz0,rint(0:nxdp),dummy(0:nxdp),bth0_old(0:nxdp)
     .               ,bzz0_old(0:nxdp),p0h(0:nxdp),q0h(0:nxdp)

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

        !Average pressure and q-profiles to half mesh
        do i=0,nxd
          p0h(i) = 0.5*(p0(i)+p0(i+1)) !Average pressure to half mesh
          q0h(i) = 0.5*(Rxq0(i)+Rxq0(i+1)) !Average pressure to half mesh
        enddo

c       Iteration

        do it=1,100

          bth0_old = bth0
          bzz0_old = bzz0

          !Perform integral of B_z/qprof on half-mesh
          dummy(0) = 0d0
          bzz0 (0) = bz0
          do i=1,nxd
cc            x1 = i*dr
cc            rh  = 0.5*(2*i-1)*dr

            x1  = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))
            dr  =      grid_params%xg(i+1)-grid_params%xg(i)
            rh  = grid_params%xg(i)

            bzh = 0.5*(bzz0(i)+bzz0(i-1))

            dummy(i) = dummy(i-1) + 2*(bzh/Rxq0(i))**2*rh*dr
          enddo

          !Calculate Bz, Btheta on half-mesh
          do i=0,nxd
cc            x1 = i*dr
            x1  = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))

cc            p0  (i) = p00*pprof(x1)
            bzz0(i) = sqrt((bz0**2+2*(p0h(0)-p0h(i))
     .                 -dummy(i))/((x1/q0h(i))**2+1))
            bth0(i) = x1*bzz0(i)/q0h(i)
          enddo

          !Check convergence
          err = sqrt(0.5*(sum((bth0-bth0_old)**2)
     .                   +sum((bzz0-bzz0_old)**2))/nxdp)

          if (ilevel > 0 .and. my_rank == 0) then
            write (*,*) '1D Grad-Shafranov equilibrium iter =',it
     .                 ,' ; Error=',err
          endif

          dr = xmax/nxd
          if (err < 0.01*dr**2) exit

        enddo

        if (ilevel == 0 .and. my_rank == 0) then
           write (*,'(a,i3,a,1pe10.2)')
     .          ' 1D equilibrium converged in '
     .          ,it,' iterations with error =',err
        endif

        !Check force balance
        do i=1,nxd
cc          rh  = 0.5*(2*i-1)*dr
          dr  = grid_params%xg(i+1)-grid_params%xg(i)
          rh  = grid_params%xg(i)

          bzh = 0.5*(bzz0(i)+bzz0(i-1))
          bth = 0.5*(bth0(i)+bth0(i-1))
          dummy(i) =  bzh*(bzz0(i)-bzz0(i-1))/dr
     .                   +(p0h (i)-p0h (i-1))/dr
     $              +bth/rh*(bth0(i)*i-bth0(i-1)*(i-1))
        enddo

        if (my_rank == 0) then
          write(*,*)
          write(*,'(a,1pe14.7)') ' Force-balance error ='
     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
        endif

        !Calculate GS equilibrium properties
        bz_avg = 0d0
        do i=1,nxd
          x1  = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))
          dr  =      grid_params%xg(i+1)-grid_params%xg(i)
          bz_avg = bz_avg + 2*bzz0(i)*x1*dr
        enddo

        if (my_rank == 0) then
          write (*,*) 
          write (*,999) ' 1D GS equilibrium properties:'
          write (*,999) '   Theta  =',bth0(nxd)/bz_avg
          write (*,999) '   Bz(r=0)=',bzz0(0)
          write (*,999) '   I0     =',2*pi*xmax*bth0(nxd)
          write (*,999) '   beta_t =',2*p0h(0)/(bzz0(0)**2)
          write (*,999) '   beta_p =',2*p0h(0)/(bth0(nxd)**2)
 999      format (a,f10.7)
        endif

        !Save some equilibrium properties (to calculate bootstrap current)
        eq_params(1) = 1d0/grid_params%params(3)       !Aspect ratio
        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
        eq_params(3) = q0h(nxd)/grid_params%params(3)  !q(r=a)

        !Transfer to integer mesh
        bth0(1:nxd) = 0.5*(bth0(0:nxd-1)+bth0(1:nxd))
        bzz0(1:nxd) = 0.5*(bzz0(0:nxd-1)+bzz0(1:nxd))

      end subroutine GS_equil_1d

c     qprof
c     ###############################################################
      function qprof(r)

c     ---------------------------------------------------------------
c     Finds cylindrical components of magnetic field B_theta, B_z
c     from q profile.
c     ---------------------------------------------------------------

        use equilibrium

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

        end select

c     End program

      end function qprof

c     pprof
c     ###############################################################
      function pprof(r)

c     ---------------------------------------------------------------
c     Finds cylindrical components of magnetic field B_theta, B_z
c     from q profile.
c     ---------------------------------------------------------------

        use equilibrium

        implicit none

c     Call variables

        real(8)    :: pprof,r

c     Local variables

c     Begin program

        select case(equil)
        case('ntm2d')

          pprof = (1-r**2)**4

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

      case ('gem')

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

ccc findVectorPotential
ccc####################################################################
cc      subroutine findVectorPotential(v,b,a_cov)
cc
ccc--------------------------------------------------------------------
ccc     Finds COVARIANT vector potential from CONTRAVARIANT magnetic
ccc     field components with gauge:
ccc                grad(phi + eta div(a)) = 0
ccc
ccc     On input:
ccc       * v -> cnv velocity
ccc       * b -> cnv magnetic field
ccc     On output:
ccc       * a -> cov vector potential
ccc--------------------------------------------------------------------
cc
cc      use parameters
cc
cc      use grid_aliases
cc
cc      use equilibrium
cc
cc      use precond_variables, ONLY: poisson
cc
cc      use  mg_internal, ONLY: mapMGVectorToArray
cc
cc      use operators
cc
cc      use imposeBCinterface
cc
cc      use timeStepping
cc
cc      implicit none
cc
ccc Call variables
cc
cc      real(8)    :: a_cov(ilom:ihip,jlom:jhip,klom:khip,3)
cc     .             ,v(ilom:ihip,jlom:jhip,klom:khip,3)
cc     .             ,b(0:nxdp,0:nydp,0:nzdp,3)
cc
ccc Local variables
cc
cc      integer    :: i,j,k,igrid,ii,bcnd(6,3)
cc      real(8)    :: cm,c0,cp,dxx,dyy,dvol,cov(3)
cc     .             ,cai_rhs(ntotdp),caiv(ntotdp)
cc     .             ,phi(ilom:ihip,jlom:jhip,klom:khip)
cc     .             ,cai(ilom:ihip,jlom:jhip,klom:khip)
cc     .             ,a_cnv(ilom:ihip,jlom:jhip,klom:khip,3)
cc     .             ,a(0:nxdp,0:nydp,0:nzdp,3)
cc
cc      logical    :: spoint
cc
ccc Begin program
cc
cc      igrid = 1
cc
cc      nx = grid_params%nxv(igrid)
cc      ny = grid_params%nyv(igrid)
cc      nz = grid_params%nzv(igrid)
cc
ccc Find equilibrium (d_t A = 0) scalar potential (vxb electric field)
cc
cccc      call findScalarPotential(v
cccc     .                        ,b(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)
cccc     .                        ,phi)
cc
ccc Find vector potential with gauge Ax = 0 in GLOBAL domain
cc
cc      call curl_inv(nxd,nyd,nzd,igrid,igrid,igrid,b,a,global=.true.)
cc
ccc Restrict solution to current domain
cc
cc      a_cov(ilom:ihip,jlom:jhip,klom:khip,:)
cc     .     = a(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)
cc
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
cc
ccc Set boundary values to face values for EQU BC
cc
cc      bcnd(:,1) = bcond
cc      bcnd(:,2) = bcond
cc      bcnd(:,3) = bcond
cc      where (bcnd == DEF) bcnd = -IFC
cc
cc      call setBC(IAX,3,nx,ny,nz,a_cnv,a_cov,a,bcnd,igrid,igrid,igrid
cc     .          ,is_cnv=.false.)
cc
ccc End program
cc
cc      end subroutine findVectorPotential

c findEquilJ
c####################################################################
      subroutine findEquilJ

c--------------------------------------------------------------------
c     Finds equilibrium current from vector vv (either acov, or bcnv,
c     depending on whether we are using B or A).
c--------------------------------------------------------------------

      use auxiliaryVariables

      use operators

      implicit none

c Call variables

c Local variables

      integer :: i,j,k,icomp,bcs(6,neqd),bcnd(6,3)

c Begin program

cc      jcov_0 = 0d0
cc      return

      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Check for imposed electric field

      if (sum(E0) /= 0d0) then
        jcov_0(:,:,:,1) = E0(1)/eta
        jcov_0(:,:,:,2) = E0(2)/eta
        jcov_0(:,:,:,3) = E0(3)/eta
        return
      endif

c Define extrapolation BCs

      call defineBoundaryConditions(neqd,bcs)

cc      bcnd(:,1) = bcond
cc      bcnd(:,2) = bcond
cc      bcnd(:,3) = bcond
cc      where (bcnd == DEF) bcnd = -EXT !Extrapolate covariant components

c Find magnetic field (A version)

#if defined(vec_pot)

      bcnd = bcs(:,IAX:IAZ)
      where (bcnd == -EQU)
        bcnd = -EXT             !Extrapolate covariant components
      end where

      call setBC(IAX,3,nx,ny,nz,acnv,acov,acov,bcnd,igx,igy,igz
     .          ,is_cnv=.false.,iorder=2)

      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            do icomp=1,3
              bcnv(i,j,k,icomp)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                                ,acov(:,:,:,1)
     .                                ,acov(:,:,:,2)
     .                                ,acov(:,:,:,3),icomp)
            enddo
          enddo
        enddo
      enddo
#else

      bcnd = bcs(:,IBX:IBZ)
      where (bcnd == -NEU)
        bcnd = -EXT             !Extrapolate covariant components
      end where

#endif

c     Extrapolate covariant magnetic field to ghost cells

      call setBC(IBX,3,nx,ny,nz,bcnv,bcov,vzeros,bcnd,igx,igy,igz
#if defined(vec_pot)
     .          ,is_cnv=.true.,iorder=1)
#else
     .          ,is_cnv=.true.,iorder=3)
#endif
cc      do k = 0,nz+1
cc        do j = 0,ny+1
cc          do i = 0,nx+1
cc            bcov(i,j,k,:)=matmul(gmetric%grid(igx)%gsub(i,j,k,:,:)
cc     .                          ,bcnv(i,j,k,:))
cc          enddo
cc        enddo
cc      enddo

      !Find equilibrium current (with BCs)
      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            do icomp=1,3
              jcnv(i,j,k,icomp)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                                ,bcov(:,:,:,1)
     .                                ,bcov(:,:,:,2)
     .                                ,bcov(:,:,:,3),icomp)
            enddo

cc            jcov_0(i,j,k,:) =
cc     .             matmul(gmetric%grid(igx)%gsub(i,j,k,:,:)
cc     .                   ,jcnv(i,j,k,:))
          enddo
        enddo
      enddo

      call setBC(IJX,3,nx,ny,nz,jcnv,jcov_0,vzeros,bcnd,igx,igy,igz
#if defined(vec_pot)
     .          ,is_cnv=.true.,iorder=1)
#else
     .          ,is_cnv=.true.,iorder=2)
#endif
cc      write (*,*) jcov_0(:,10,1,1)
cc      write (*,*)
cc      write (*,*) jcov_0(:,10,1,2)
cc      write (*,*)
cc      write (*,*) jcov_0(:,10,1,3)

c End program

      end subroutine findEquilJ

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


