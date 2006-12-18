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

      implicit none

c Call variables

      integer(4)    :: iigx,iigy,iigz

      real(8)       :: var(ilom:ihip,jlom:jhip,klom:khip,neqd)

      character(20) :: label(neqd)

      integer(4)    :: bcs(6,neqd)

c Local variables

      integer(4)    :: i,j,k,ig,jg,kg,ieq,nmax,iglobal,ierr,system,sh

      real(8)       :: ldaflow,x1,y1,z1,cx,cy,cz

      logical       :: covariant,to_cartsn,cartsn,a_bc

      real(8)       :: r,jac1,bnorm,aaa,bbb,ccc,qq,qqp,q0,RR,ff,nn
     .                ,aspect_ratio,mm,kk,Iz,bz0,btheta,vr,rrho,r1,bb
     .                ,r0(1),aa,dummy,dummy1,dummy2,dummy3,ieta

      real(8)       :: rint(0:nxdp),dh,offset(1),bth0(0:nxdp)
     .                ,bzz0(0:nxdp),vr0(0:nxdp),rho0(0:nxdp)
     .                ,eth0(0:nxdp),ez0(0:nxdp),jth0(0:nxdp)
     .                ,jz0 (0:nxdp)

      real(8)       :: a(ilom:ihip,jlom:jhip,klom:khip,3)

      character(200):: command

#if defined(vec_pot)
      real(8)       :: b(0:nxdp,0:nydp,0:nzdp,3)
#endif

c Functions

      !q-profile
      qq (r1) = 0.6125*(1 - 1.8748*r1**2 + 0.8323*r1**4)
      qqp(r1) = 0.6125*( -2*1.8748*r1   +4*0.8323*r1**3)
cc      qq(r1)  = 0.3*(1 - 1.8748*r1**2 + 0.8323*r1**4)
cc      qqp(r1) = 0.3*( -2*1.8748*r1   +4*0.8323*r1**3)

      ff(r1) = r1**2 + qq(r1)**2

c Begin program

      ldaflow = dlambda/rshear

      var = 0d0

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
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              var(i,j,k,IAX)  = 0d0
              var(i,j,k,IAY)  = x1
              var(i,j,k,IAZ)  = 0d0
            enddo
          enddo
        enddo
#endif

      case ('mswsn') !Magnetosonic wave

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
              var(i,j,k,IVX:IVZ)=vperflow
     .                          *tanh((x1-0.5d0*(xmax-xmin))/dlambda)
            enddo
          enddo
        enddo

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

             var(i,j,k,IVX) = vperflow*tanh((y1-.5*(ymax-ymin))/ldaflow)
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
cc        var(:,:,:,ITMP) = 0d0

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
cc        var(:,:,:,IAX:IAZ) = a
        do k = 0,nzdp
          do j = 0,nydp
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              b(i,j,k,1)  = 0d0                  
              b(i,j,k,2)  = tanh((x1-0.5d0*(xmax-xmin))/dlambda)
              b(i,j,k,3)  = sqrt(bz0**2 - b(i,j,k,2)**2)
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
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

             var(i,j,k,IRHO)=1d0/cosh((x1-0.5d0*(xmax-xmin))/dlambda)**2
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

c     Tokamak equilibria (analytical)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)

        Iz           = 2.22
        bb           = 0.5*Iz

        nh2 = mm !To set the right perturbation wavelength

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
     .                        + 1/16./var(i,j,k,IRHO)
     .                         *(Iz*dlambda/(1+(x1/dlambda)**2))**2
cc              var(i,j,k,ITMP) = 1d-0
            enddo
          enddo
        enddo

#if defined(vec_pot)
        do k = 0,nzdp
          do j = 0,nydp
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              btheta = bb*x1/(1+(x1/dlambda)**2)
              bz0    = RR

              b(i,j,k,1)  = 0d0                  
              b(i,j,k,2)  = btheta + kk*x1/mm*bz0
              b(i,j,k,3)  = x1*bz0
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
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

c     Tokamak equilibria (analytical)

        Iz = 2.22
        bb = 0.5*Iz

        nh2 = mm !To set the right perturbation wavelength
        nh3 = kk*RR    !To set the right perturbation wavelength

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
     .                        + 1/16./var(i,j,k,IRHO)
     .                         *(Iz*dlambda/(1+(x1/dlambda)**2))**2

            enddo
          enddo
        enddo

#if defined(vec_pot)
        do k = 0,nzdp
          do j = 0,nydp
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              btheta = bb*x1/(1+(x1/dlambda)**2)
              bz0    = RR

              b(i,j,k,1)  = 0d0                  
              b(i,j,k,2)  = btheta
              b(i,j,k,3)  = x1*bz0
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('ppnch')

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)

        nh2 = mm !To set the right perturbation wavelength

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
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
cc          write (*,*)
cc          write (*,*) dummy1,bth0(i),bzz0(i),dummy2,dummy3,vr0(i)
cc     .                ,dummy4,dummy5,rho0(i)
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
        do k = 0,nzdp
          do j = 0,nydp
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              b(i,j,k,1) = 0d0
              b(i,j,k,2) = bth0(i) + kk*x1/mm*bzz0(i)
              b(i,j,k,3) = x1*bzz0(i)
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
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
        do k = 0,nzdp
          do j = 0,nydp
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              b(i,j,k,1) = 0d0
              b(i,j,k,2) = bth0(i)
              b(i,j,k,3) = x1*bzz0(i)
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
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

        nh2 = mm !To set the right perturbation wavelength

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
        do k = 0,nzdp
          do j = 0,nydp
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              b(i,j,k,1) = 0d0
              b(i,j,k,2) = bth0(i) + kk*x1/mm*bzz0(i)
              b(i,j,k,3) = x1*bzz0(i)
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
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

        nh2 = mm !To set the right perturbation wavelength

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
        do k = 0,nzdp
          do j = 0,nydp
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              if (coords == 'hl2') x1 = sign(sqrt(2*abs(x1)),x1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              if (coords == 'hel') then
                b(i,j,k,1)  = 0d0
                b(i,j,k,2)  = btheta + kk*x1/mm*bz0
                b(i,j,k,3)  = x1*bz0
              else
                b(i,j,k,1)  = 0d0
                b(i,j,k,2)  = btheta/x1 + kk/mm*bz0
                b(i,j,k,3)  = bz0
              endif
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
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
        do k = 0,nzdp
          do j = 0,nydp
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              b(i,j,k,1) = 0d0
              b(i,j,k,2) = btheta
              b(i,j,k,3) = x1*bz0
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('tok')

c     Cylindrical Tokamak analytical equilibrium (Furth et al., 1973)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        q0 = grid_params%params(4)

        nh2 = mm !To set the right perturbation wavelength

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
        do k = 0,nzdp
          do j = 0,nydp
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              b(i,j,k,1) = 0d0
              b(i,j,k,2) = btheta + kk*x1/mm*bz0
              b(i,j,k,3) = x1*bz0
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
#endif

#if defined(vmec)

      case ('vmec')

        mm = grid_params%params(1)
        nn = grid_params%params(2)

        RR = 1d0   !To avoid affecting q-profile

        kk = 0d0   !To avoid affecting q-profile

c     Check coordinates

        if (coords /= 'vmc') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Read VMEC equilibrium

        nh2 =  mm       !To set the right perturbation wavelength
        nh3 = -nn       !To set the right perturbation wavelength

        !Build equilibrium

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

#if !defined(vec_pot)
        call vmec_equ(igx,nx,ny,nz,var(:,:,:,IBX),var(:,:,:,IBY)
     .                            ,var(:,:,:,IBZ),var(:,:,:,ITMP))
#else
        call vmec_equ(igx,nx,ny,nz,b(:,:,:,1),b(:,:,:,2)
     .                            ,b(:,:,:,3),var(:,:,:,ITMP))

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))
#endif

#endif
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

c Transfer electric field info to postprocessor via input file

      if (sum(E0) /= 0d0 .and. my_rank == 0) then

        command = 'grep E0 ' // trim(inputfile) // ' > /dev/null'
        ierr = system(trim(command))

        write (command,*) '   E0 =',real(E0(1),4)
     .                         ,',',real(E0(2),4)
     .                         ,',',real(E0(3),4)


        if (ierr /= 0) then  !Replace E0
          command='sed -e"/graphdef/{" -e' // achar(39)
     .       // 'a\' // achar(39) //' -e"'
     .       // trim(command)   // '" -e"}" '
     .       // trim(inputfile) // ' > temp'
        else                 !Add E0
          command='sed -e"/E0/{" -e' // achar(39)
     .       // 'c\' // achar(39) //' -e"'
     .       // trim(command)   // '" -e"}" '
     .       // trim(inputfile) // ' > temp'
        endif

        ierr = system(trim(command))

        command = 'mv temp ' // trim(inputfile)
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

        integer(4) :: i,j,k,ig,jg,kg
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

      integer(4) :: iigx,iigy,iigz
      real(8)    :: a(ilom:ihip,jlom:jhip,klom:khip,3)
      logical    :: bc

c Local variables

      integer(4) :: ig,jg,kg,i,j,k,bcnd(6,3),nx,ny,nz
      real(8)    :: x1,y1,z1
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
     .                                      ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = -dlambda
     .             *dlog(dcosh((x1-0.5d0*(xmax-xmin))/dlambda)) 
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
     .               *atan(exp((x1-0.5d0*(xmax-xmin))/dlambda))
              a(i,j,k,3) = -dlambda
     .              *dlog(dcosh((x1-0.5d0*(xmax-xmin))/dlambda)) 
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
     .             *dlog(dcosh((y1-0.5d0*(ymax-ymin))/dlambda)) 
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

      use grid

      use equilibrium

      use precond_variables

      use matvec

      use operators

      use imposeBCinterface

      implicit none

c Call variables

      real(8)    :: a_cov(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,v(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,b(0:nxdp,0:nydp,0:nzdp,3)

c Local variables

      integer(4) :: igrid,ii,bcnd(6,3)
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

cc      a(:,:,:,1)=0d0  !Gauge
cc
cc      a(0,0,:,3)=0d0
cc
cc      i = 0
cc
cccc      spoint = isSP(1,1,1,igx,igy,igz)  !Check whether current domain contains SP
cc      spoint = .false.
cc
cc      do j=1,nyd+1
cc        dyy = grid_params%yg(j)-grid_params%yg(j-1)
cc
cc        if (spoint) then  !r=0 face
cc          a(i,j,:,3)=a(i,j-1,:,3) + dyy*2.5d-1*(b(i  ,j-1,:,1)
cc     .                                         +b(i  ,j  ,:,1)
cc     .                                         +b(i+1,j-1,:,1)
cc     .                                         +b(i+1,j  ,:,1))
cc        else
cc          a(i,j,:,3)=a(i,j-1,:,3) + dyy*  5d-1*(b(i,j-1,:,1)
cc     .                                         +b(i,j  ,:,1))
cc        endif
cc      enddo
cc
cc      a(0,:,:,2)=0d0
cc
cc      do i=1,nxd+1
cc        if (spoint) then
cc
cc          dxx = 0.5*(grid_params%xg(i+1)-grid_params%xg(i-1))
cc
cc          a(i,:,:,3) = a(i-1,:,:,3) - dxx*b(i,:,:,2)
cc
cc          if (isSP(i,1,1,igx,igy,igz)) then
cc            a(i,:,:,2) = a(i-1,:,:,2) + 5d-1*dxx*b(i,:,:,3)  !Factor of 1/2 due to geometry
cc          else
cc            a(i,:,:,2) = a(i-1,:,:,2) +      dxx*b(i,:,:,3)
cc          endif
cc
cc        else
cc
cc          dxx = grid_params%xg(i)-grid_params%xg(i-1)
cc
cc          a(i,:,:,3) = a(i-1,:,:,3) - dxx*5d-1*(b(i-1,:,:,2)
cc     .                                         +b(i  ,:,:,2))
cc          a(i,:,:,2) = a(i-1,:,:,2) + dxx*5d-1*(b(i-1,:,:,3)
cc     .                                         +b(i  ,:,:,3))
cc
cc        endif
cc      enddo
cc
cc      !Average from radial faces to nodes in SP coordinate systems
cc      if (spoint) a(1:nxd+1,:,:,2:3) = 0.5*(a(1:nxd+1,:,:,2:3)
cc     .                                     +a(0:nxd  ,:,:,2:3))

      call curl_inv(nxd,nyd,nzd,igrid,igrid,igrid,b,a,global=.true.)

c Restrict solution to current domain

      a_cov(ilom:ihip,jlom:jhip,klom:khip,:)
     .     = a(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)

c Impose topological constraints and find a_cnv

      bcnd(:,1) = bcond
      bcnd(:,2) = bcond
      bcnd(:,3) = bcond
cc      where (bcnd == DEF) bcnd = EXT

      call setBC(IAX,3,nx,ny,nz,a_cnv,a_cov,vzeros,bcnd
     .          ,igrid,igrid,igrid,is_cnv=.false.,iorder=2)

c Transform A to gauge div(A)=0

      !Find gauge transformation potential cai
      do k=klo,khi
        do j=jlo,jhi
          do i=ilo,ihi
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            if (vol_wgt) then
              dvol = volume(i,j,k,igrid,igrid,igrid)
            else
              dvol = 1d0
            endif

            cai_rhs(ii)=-dvol*(div(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                            ,a_cnv(:,:,:,1)
     .                            ,a_cnv(:,:,:,2)
     .                            ,a_cnv(:,:,:,3)))
          enddo
        enddo
      enddo

      !Transform A (if rhs /= 0d0)
      if (sqrt(sum(cai_rhs*cai_rhs)) > 0d0) then

        !Solve Poisson problem to find Gauge potential cai
#if defined(petsc)
        call pstop('findVectorPotential'
     .          ,'Poisson solver not implemented in parallel')
#else
        write (*,*) 'Solving for Gauge potential'

        !Define BCs
        bcnd(:,1) = bcond
        where (bcnd == DEF) bcnd = DIR

        !Solve Poisson eq
        call poisson(ntotdp,cai_rhs,caiv,30,1d-6,bcnd(:,1),ilevel,0)

        !Map to array (w BCs)
        call mapMGVectorToArray(0,1,caiv,nx,ny,nz,cai,igrid,.false.)

        call setBC(IRHO,nx,ny,nz,cai,zeros,bcnd(:,1),igrid,igrid,igrid
     .            ,iorder=2)
#endif

        !Gauge transformation A' = A + grad(cai)
        do k=klom,khip
          do j=jlom,jhip
            do i=ilom,ihip
              call grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,cai
     .                 ,cov(1),cov(2),cov(3))
              a_cov(i,j,k,:) = a_cov(i,j,k,:) + cov
            enddo
          enddo
        enddo

      endif

c Set boundary values to face values for EQU BC

      bcnd(:,1) = bcond
      bcnd(:,2) = bcond
      bcnd(:,3) = bcond
      where (bcnd == DEF) bcnd = -IFC

      call setBC(IAX,3,nx,ny,nz,a_cnv,a_cov,a,bcnd,igrid,igrid,igrid
     .          ,is_cnv=.false.)

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
cc      integer(4) :: igrid,ii,bcnd(6,1)
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


