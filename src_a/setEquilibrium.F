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

      character*(20):: label(neqd)

      integer(4)    :: bcs(6,neqd)

c Local variables

      integer(4)    :: i,j,k,ig,jg,kg,ieq,nmax,iglobal
      real(8)       :: ldaflow,x1,y1,z1

      logical       :: covariant,to_cartsn,cartsn

      real(8)       :: r,jac1,bnorm,aaa,bbb,ccc,qq,qqp,q0,RR,ff,nn
     .                ,aspect_ratio,mm,kk,Iz,bz0,btheta,vr,rrho,r1,bb
     .                ,r0(1),aa,dummy,dummy1,dummy2,dummy3,dummy4,dummy5
     .                ,cx,cy,cz

      real(8)       :: rint(0:nxd+1),dh,offset(1),bth0(0:nxd+1)
     .                ,bzz0(0:nxd+1),vr0(0:nxd+1),rho0(0:nxd+1)

      real(8)       :: b(ilom:ihip,jlom:jhip,klom:khip,3)

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

c Initialize required local grid information

      igx = iigx
      igy = iigy
      igz = iigz

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Label variables
     
      label(IRHO) = 'Rho'
      label(IAX)  = 'A_1'
      label(IAY)  = 'A_2'
      label(IAZ)  = 'A_3'
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

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

              var(i,j,k,IRHO) = 1d0

cc              var(i,j,k,IBX)  = 0d0
cc              var(i,j,k,IBY)  = 0d0
cc              var(i,j,k,IBZ)  = 1d0
              var(i,j,k,IAX)  = -y1
              var(i,j,k,IAY)  = 0d0
              var(i,j,k,IAZ)  = 0d0

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,ITMP) = 1d0
            enddo
          enddo
        enddo

      case ('mswsn') !Magnetosonic wave

        gamma = 1d0

c     Check coordinates

        if (coords /= 'sin') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,ITMP) = 1d0
            enddo
          enddo
        enddo

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

              jac1 = gmetric%grid(igx)%jac(i,j,k)

              b(i,j,k,1) = 0d0
              b(i,j,k,2) = 0d0
              b(i,j,k,3) = jac1*1d0
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))

      case ('khcar')

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Kelvin-Helmholtz with constant magnetic field in cartesian coordinates

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

             var(i,j,k,IVX)=vperflow*tanh((y1-.5*(ymax-ymin))/ldaflow)
             var(i,j,k,IVY)=0d0
             var(i,j,k,IVZ)=0d0
            enddo
          enddo
        enddo

        var(:,:,:,IRHO) = 1d0
        var(:,:,:,IBX)  = 0d0
        var(:,:,:,IBY)  = 0d0
        var(:,:,:,IBZ)  = 1d0
        var(:,:,:,ITMP) = 1d0

      case ('3dkh')

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

             !X-Y equilibrium
             var(i,j,k,IAX) = dlambda
     .             *dlog(dcosh((y1-0.5d0*(ymax-ymin))/dlambda))
             var(i,j,k,IAY) = 0d0
             var(i,j,k,IAZ) = 0d0

             var(i,j,k,IRHO)= (2d0 - var(i,j,k,IBZ)**2)

             var(i,j,k,IVX) = vperflow*tanh((y1-.5*(ymax-ymin))/ldaflow)
     .                                *sin (2*pi*z1/(zmax-zmin))
             var(i,j,k,IVY) = 0d0
             var(i,j,k,IVZ) = 0d0

            enddo
          enddo
        enddo

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

             var(i,j,k,IRHO)= (2d0
     .                       - curl2(i,j,k,nxd,nyd,nzd,igx,igy,igz
     .                              ,var(:,:,:,IAX)
     .                              ,var(:,:,:,IAY)
     .                              ,var(:,:,:,IAZ),3)**2)

            enddo
          enddo
        enddo

        var(:,:,:,ITMP) = 0.25

      case ('tmcar')

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

        var(:,:,:,ITMP) = 1d-0

        bz0 = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              b(i,j,k,1) = 0d0
              b(i,j,k,2) = tanh((x1-0.5d0*(xmax-xmin))/dlambda)
              b(i,j,k,3) = sqrt(bz0**2 - b(i,j,k,2)**2)            
            enddo
          enddo
        enddo

c diag *****
cc        b(:,:,:,1) = 0d0
cc        b(:,:,:,2) = 1d0
cc        b(:,:,:,3) = 0d0
c diag *****

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))

cc      case ('tmcxz')
cc
ccc     Define vector potential (in curvilinear coordinates) for initialization
cc
cc        call fillVectorPotential(a1,a2,a3,igx,igy,igz)
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
cc        call fillVectorPotential(a1,a2,a3,igx,igy,igz)
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

c     Check coordinates

        if (coords /= 'sin') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Tearing mode in sinusoidal coordinates

cc        do k = klo,khi
cc          do j = jlo,jhi
cc            do i = ilo,ihi
cc
cc              !X-Y equilibrium
cc             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,1)
cc             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,2)
cc             var(i,j,k,IBZ)=0d0
cc
cc             gsub = gmetric%grid(igx)%gsub(i,j,k,:,:)
cc             jac1 = gmetric%grid(igx)%jac (i,j,k)
cc
cc             bnorm= vectorNorm(i,j,k,igx,igy,igz
cc     .                         ,var(i,j,k,IBX)
cc     .                         ,var(i,j,k,IBY)
cc     .                         ,var(i,j,k,IBZ),.false.)
cc
cc             ccc = jac1*(bz0**2 - bnorm)
cc             bbb = gsub(3,2)*var(i,j,k,IBY) + gsub(3,1)*var(i,j,k,IBX)
cc             aaa = gsub(3,3)
cc
cc             var(i,j,k,IBZ)  = (-bbb+sqrt(bbb**2+4*aaa*ccc))/2./aaa
cc
ccc diag ****
cccc             var(i,j,k,IBX)  = 0d0
cccc             var(i,j,k,IBY)  = 0d0
cccc             var(i,j,k,IBZ)  = 0d0
ccc diag ****
cc
cc
cccc             if (bbb**2+4*aaa*ccc < 0d0) then
cccc               write (*,*) var(i,j,k,IBZ),bnorm,aaa,bbb**2+4*aaa*ccc
cccc             endif
cc            enddo
cc          enddo
cc        enddo

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

        bz0 = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              cx = 0d0
              cy = tanh((x1-0.5d0*(xmax-xmin))/dlambda)
              cz = sqrt(bz0**2 - cy**2)  

              call transformVectorToCurvilinear(i,j,k,igx,igy,igz
     .                                       ,cx,cy,cz,.false.
     .                                       ,b(i,j,k,1)
     .                                       ,b(i,j,k,2)
     .                                       ,b(i,j,k,3))          
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))

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

cc        aspect_ratio = 0.3
cc        mm           = 1d0
cc        RR           = 1d0/aspect_ratio
cc        kk           = -mm/RR

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

              b(i,j,k,1)  = 0d0                  
              b(i,j,k,2)  = btheta + kk*x1/mm*bz0
              b(i,j,k,3)  = x1*bz0               
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

cc              btheta = bb*x1/(1+(x1/dlambda)**2)
cc              bz0    = RR
cc
cc              var(i,j,k,IBX)  = 0d0
cc              var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
cc              var(i,j,k,IBZ)  = x1*bz0               

              var(i,j,k,ITMP) = 1d-5
     .                        + 1/16./var(i,j,k,IRHO)
     .                         *(Iz*dlambda/(1+(x1/dlambda)**2))**2

            enddo
          enddo
        enddo

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

              b(i,j,k,1)  = 0d0                  
              b(i,j,k,2)  = btheta
              b(i,j,k,3)  = x1*bz0               
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

cc              btheta = bb*x1/(1+(x1/dlambda)**2)
cc              bz0    = RR
cc
cc              var(i,j,k,IBX)  = 0d0
cc              var(i,j,k,IBY)  = btheta
cc              var(i,j,k,IBZ)  = x1*bz0

              var(i,j,k,ITMP) = 1d-5
     .                        + 1/16./var(i,j,k,IRHO)
     .                         *(Iz*dlambda/(1+(x1/dlambda)**2))**2

            enddo
          enddo
        enddo

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
        open (unit=1000,file='par_pinch_eq.txt',status='old')

        read(1000,'(a)')
        do i=0,nxd+1
          read(1000,'(9f12.8)') dummy1,bth0(i),bzz0(i),dummy2,dummy3
     .                ,vr0(i),dummy4,dummy5,rho0(i)
cc          write (*,*)
cc          write (*,*) dummy1,bth0(i),bzz0(i),dummy2,dummy3,vr0(i)
cc     .                ,dummy4,dummy5,rho0(i)
        enddo

        close(1000)

        write (*,*) 'Done!'

        !Build equilibrium
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1
              vr      = vr0 (iglobal)
              rrho    = rho0(iglobal)

              if (nc_eom_v) then
                var(i,j,k,IVX) = x1*vr
              else
                var(i,j,k,IVX) = x1*vr*rrho
              endif
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = rrho

              var(i,j,k,ITMP) = 1d-5

            enddo
          enddo
        enddo

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1
              btheta  = bth0(iglobal)
              bz0     = bzz0(iglobal)

              !X-Y equilibrium
              b(i,j,k,1) = 0d0
              b(i,j,k,2) = btheta + kk*x1/mm*bz0
              b(i,j,k,3) = x1*bz0
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))

      case ('rfp1')

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
        rint(nxd+1) = 0d0  !First guess to rint(r=1d0) = 0d0
        do i=nxd,1,-1
          x1 = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))  !Position in global grid
          dh = grid_params%xg(i+1)-grid_params%xg(i)
          rint(i) = rint(i+1) + dh*qq(x1)*qqp(x1)/ff(x1)
        enddo
        rint(0) = rint(1)  !Neumann BC for integral at r=0d0

        !Correction to ensure rint(r=1d0) = 0d0
        r0 = 1d0
        call IntDriver1d(nxd+2,grid_params%xg(0:nxd+1)
     .                  ,rint,1,r0,offset,2)

        rint = rint - offset(1)

        !Determine Iz so that Bz(r=0)=1
        r0 = 0d0
        call IntDriver1d(nxd+2,grid_params%xg(0:nxd+1)
     .                  ,rint,1,r0,offset,2)

cc        btheta = 1./2./pi*ff(1d0)/ff(0d0)*exp(-rint(0))
        btheta = 1./2./pi*ff(1d0)/ff(0d0)*exp(-offset(1))
        Iz = 1./btheta/qq(0d0)

        !Build equilibrium

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d-5

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1
              btheta  = Iz/2./pi*x1*ff(1d0)/ff(x1)*exp(-rint(iglobal))
              bz0     = qq(x1)*btheta/x1

              !X-Y equilibrium
              b(i,j,k,1) = 0d0
              b(i,j,k,2) = btheta + kk*x1/mm*bz0
              b(i,j,k,3) = x1*bz0
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))

      case ('rfp2')

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

        nh2 = mm !To set the right perturbation wavelength

        !Build equilibrium

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d-5

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              b(i,j,k,1)  = 0d0
              b(i,j,k,2)  = btheta + kk*x1/mm*bz0
              b(i,j,k,3)  = x1*bz0
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))

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

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

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

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              !X-Y equilibrium
              b(i,j,k,1) = 0d0
              b(i,j,k,2) = btheta + kk*x1/mm*bz0
              b(i,j,k,3) = x1*bz0
            enddo
          enddo
        enddo

        call findVectorPotential(var(:,:,:,IVX:IVZ)
     .                          ,b
     .                          ,var(:,:,:,IAX:IAZ))

      case ('rhots')

c     Cylindrical Tokamak analytical equilibrium (Furth et al., 1973)

        if (coords /= 'cyl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        !Build equilibrium
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              !X-Y equilibrium
              var(i,j,k,IAX)  = 0d0
              var(i,j,k,IAY)  = 0d0
              var(i,j,k,IAZ)  = 0d0

              var(i,j,k,IVX)  = x1*cos(y1)
              var(i,j,k,IVY)  = -sin(y1)
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0
     .                      +1d-1*x1*exp(-(x1-0.3)**2/dlambda**2)

              var(i,j,k,ITMP) = 0d0

            enddo
          enddo
        enddo

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

c End program

      end subroutine setEquilibrium

ccc fillVectorPotential
ccc####################################################################
cc      subroutine fillVectorPotential(a1,a2,a3,igx,igy,igz)
cc
ccc--------------------------------------------------------------------
ccc     Defines COVARIANT vector potential for initialization of 
ccc     equilibrium.
ccc--------------------------------------------------------------------
cc
cc      use parameters
cc
cc      use grid
cc
cc      use equilibrium
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer(4) :: igx,igy,igz
cc      real(8)    :: a1(ilom:ihip,jlom:jhip,klom:khip)
cc     .             ,a2(ilom:ihip,jlom:jhip,klom:khip)
cc     .             ,a3(ilom:ihip,jlom:jhip,klom:khip)
cc
ccc Local variables
cc
cc      integer(4) :: ig,jg,kg,i,j,k
cc      real(8)    :: xx,yy,zz
cc      logical    :: cartsn
cc
ccc Begin program
cc
cc      select case (trim(equil))
cc
cc      case ('khcar','tmcar','tmsin','tmcxz')
cc
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,xx,yy,zz)
cc              a1(i,j,k) = 0d0
cc              a2(i,j,k) = 0d0
cc              a3(i,j,k) = dlambda
cc     .             *dlog(dcosh((xx-0.5d0*(xmax-xmin))/dlambda)) 
cc            enddo
cc          enddo
cc        enddo
cc
cc      case ('tmcyz','3dkh')
cc
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,xx,yy,zz)
cc              a1(i,j,k) = dlambda
cc     .             *dlog(dcosh((yy-0.5d0*(ymax-ymin))/dlambda)) 
cc              a2(i,j,k) = 0d0
cc              a3(i,j,k) = 0d0
cc            enddo
cc          enddo
cc        enddo
cc
cc      end select
cc
ccc End program
cc
cc      end subroutine fillVectorPotential

c findVectorPotential
c####################################################################
      subroutine findVectorPotential(v,b,a)

c--------------------------------------------------------------------
c     Finds COVARIANT vector potential from CONTRAVARIANT magnetic
c     field components for initialization of equilibrium.
c     
c     Integrals done by Numerical Recipes p. 128, 4.1.12 -- actually
c     one bit better -- extrapolates quadratically to f(-1) and uses
c     the internal formula, rather than using midpoint for the first
c     point. Thus even the first point has third order accuracy.
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

      real(8)    :: a(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,v(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,b(ilom:ihip,jlom:jhip,klom:khip,3)

c Local variables

      integer(4) :: igrid,ii,bcnd(6,3)
      real(8)    :: cm,c0,cp,dxx,dyy,dvol,cov(3)
     .             ,cai_rhs(ntotdp),caiv(ntotdp)
     .             ,phi(ilom:ihip,jlom:jhip,klom:khip)
     .             ,cai(ilom:ihip,jlom:jhip,klom:khip)
     .             ,a_cnv(ilom:ihip,jlom:jhip,klom:khip,3)

c Begin program

      igrid = 1

      nx = nxd
      ny = nyd
      nz = nzd

c Find scalar potential

      call findScalarPotential(v,b,phi)

c Find vector potential with gauge ax = 0

cc      cm=-0.0833333333333
cc      c0= 0.6666666666667
cc      cp= 0.4166666666667
cc
cc      do k=1,nz
cc        a(1,1,k,3)=0d0
cc      enddo
cc
cc      do k=1,nz
cc        do j=2,ny
cc          call getCartesianCoordinates(i,j,k,igrid,igrid,igrid,ig,jg,kg
cc     .                                ,x1,y1,z1)
cc
cc          dyy = grid_params%dyh(jg)
cc
cc          if(j.eq.2) then
cc            a(1,j,k,3)=a(1,j-1,k,3)+
cc     $         dyy*(cp*b(1,j-1,k,1)+c0*b(1,j,k,1)+cm*b(1,j+1,k,1))
cc          else
cc            a(1,j,k,3)=a(1,j-1,k,3)+dyy*(cm*b(1,j-2,k,1)+
cc     $                  c0*b(1,j-1,k,1)+cp*b(1,j,k,1))
cc          endif
cc        enddo
cc      enddo
cc
cc      do k=1,nz
cc        do j=1,ny
cc          a(1,j,k,2)=0d0
cc        enddo
cc      enddo
cc
cc      do k=1,nz
cc        do j=1,ny
cc          do i=2,nx
cc
cc            call getCartesianCoordinates(i,j,k,igrid,igrid,igrid
cc     .                                  ,ig,jg,kg,x1,y1,z1)
cc
cc            dxx = grid_params%dxh(ig)
cc
cc            if(i.eq.2) then
cc              a(i,j,k,3)=a(i-1,j,k,3)-
cc     .           dxx*(cp*b(i-1,j,k,2)+c0*b(i,j,k,2)+cm*b(i+1,j,k,2))
cc              a(i,j,k,2)=a(i-1,j,k,2)+
cc     .           dxx*(cp*b(i-1,j,k,3)+c0*b(i,j,k,3)+cm*b(i+1,j,k,3))
cc            else
cc              a(i,j,k,3)=a(i-1,j,k,3)-dxx*(cm*b(i-2,j,k,2)+
cc     .                     c0*b(i-1,j,k,2)+cp*b(i,j,k,2))
cc              a(i,j,k,2)=a(i-1,j,k,2)+dxx*(cm*b(i-2,j,k,3)+
cc     .                     c0*b(i-1,j,k,3)+cp*b(i,j,k,3))
cc            endif
cc          enddo
cc        enddo
cc      enddo

c Find vector potential with gauge ax = 0
c (with SP, it integrates at faces in r, and then averages to nodes)

      a(0,0,:,3)=0d0

      i = 0
      do k=1,nz+1
        do j=1,ny+1
          call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

          dyy = grid_params%dy(jg-1)

cc          if(j.eq.2) then
cc            a(i,j,k,3)=a(i,j-1,k,3)
cc     $                +dyy*(cp*b(i,j-1,k,1)
cc     .                     +c0*b(i,j  ,k,1)
cc     .                     +cm*b(i,j+1,k,1))
cc          else
cc            a(i,j,k,3)=a(i,j-1,k,3)
cc     .                +dyy*(cm*b(i,j-2,k,1)
cc     $                     +c0*b(i,j-1,k,1)
cc     .                     +cp*b(i,j  ,k,1))
cc          endif

          if (bcond(1) == SP) then
            a(i,j,k,3)=a(i,j-1,k,3) + dyy*0.25*(b(i  ,j-1,k,1)
     .                                         +b(i  ,j  ,k,1)
     .                                         +b(i-1,j-1,k,1)
     .                                         +b(i-1,j  ,k,1))
          else
            a(i,j,k,3)=a(i,j-1,k,3) + dyy*0.5*(b(i,j-1,k,1)
     .                                        +b(i,j  ,k,1))
          endif

        enddo
      enddo

      a(0,:,:,2)=0d0

      do k=1,nz+1
        do j=1,ny+1
          do i=1,nx+1

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            dxx = grid_params%dx(ig-1)

cc            if(i.eq.2) then
cc              a(i,j,k,3)=a(i-1,j,k,3)
cc     .                  -dxx*(cp*b(i-1,j,k,2)
cc     .                       +c0*b(i  ,j,k,2)
cc     .                       +cm*b(i+1,j,k,2))
cc              a(i,j,k,2)=a(i-1,j,k,2)
cc     .                  +dxx*(cp*b(i-1,j,k,3)
cc     .                       +c0*b(i  ,j,k,3)
cc     .                       +cm*b(i+1,j,k,3))
cc            else
cc              a(i,j,k,3)=a(i-1,j,k,3)
cc     .                   -dxx*(cm*b(i-2,j,k,2)
cc     .                        +c0*b(i-1,j,k,2)
cc     .                        +cp*b(i  ,j,k,2))
cc              a(i,j,k,2)=a(i-1,j,k,2)
cc     .                   +dxx*(cm*b(i-2,j,k,3)
cc     .                        +c0*b(i-1,j,k,3)
cc     .                        +cp*b(i  ,j,k,3))
cc            endif

            if (bcond(1) == SP) then
              dxx = grid_params%dxh(ig-1)
              a(i,j,k,3)=a(i-1,j,k,3) -dxx*(b(i-1,j,k,2))

              if (isSP(i,j,k,igrid,igrid,igrid)) then
                a(i,j,k,2)=a(i-1,j,k,2) +dxx*0.5*(b(i-1,j,k,3))  !Factor of 1/2 due to geometry
              else
                a(i,j,k,2)=a(i-1,j,k,2) +dxx*(b(i-1,j,k,3))
              endif
            else
              dxx = grid_params%dx(ig-1)
              a(i,j,k,3)=a(i-1,j,k,3) -dxx*0.5*(b(i-1,j,k,2)
     .                                         +b(i  ,j,k,2))

              a(i,j,k,2)=a(i-1,j,k,2) +dxx*0.5*(b(i-1,j,k,3)
     .                                         +b(i  ,j,k,3))
            endif
          enddo
        enddo
      enddo

      !Average from radial faces to nodes in SP coordinate systems
      if (bcond(1) == SP) then
        do i=1,nx
          a(i,:,:,:) = 0.5*(a(i,:,:,:)+a(i+1,:,:,:))
        enddo
      endif

      !Impose boundary conditions (only topological constraints:PER,SP)
      bcnd(:,1) = bcond
      bcnd(:,2) = bcond
      bcnd(:,3) = bcond
cc      where (bcnd == DEF) bcnd = EXT

      call setBC(IAX,3,nx,ny,nz,a_cnv,a,vzeros,bcnd,igrid,igrid,igrid
     .          ,is_cnv=.false.,iorder=2)

c Find gauge transformation potential cai (to gauge phi+eta div(A)=0)

      do k=klo,khi
        do j=jlo,jhi
          do i=ilo,ihi
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            if (vol_wgt) then
              dvol = volume(i,j,k,igrid,igrid,igrid)
            else
              dvol = 1d0
            endif

            cai_rhs(ii)=-dvol*(phi(i,j,k)
     .                        +eta*div(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                                ,a_cnv(:,:,:,1)
     .                                ,a_cnv(:,:,:,2)
     .                                ,a_cnv(:,:,:,3)))
          enddo
        enddo
      enddo

      if (sqrt(sum(cai_rhs*cai_rhs)) > 0d0) then

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

        call setMGBC(0,1,nx,ny,nz,igrid,cai,bcnd,icomp=IRHO
     .            ,is_vec=.false.,iorder=2)
#endif

      else

        cai = 0d0

      endif

c Find gauge-transformed vector potential

      do k=klom,khip
        do j=jlom,jhip
          do i=ilom,ihip
            cov = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,cai)
            a(i,j,k,:) = a(i,j,k,:) + cov
          enddo
        enddo
      enddo

c Set boundary values to face values for EQU BC

      bcnd(:,1) = bcond
      bcnd(:,2) = bcond
      bcnd(:,3) = bcond
      where (bcnd == DEF) bcnd = -IFC

      call setBC(IAX,3,nx,ny,nz,a_cnv,a,a,bcnd,igrid,igrid,igrid
     .          ,is_cnv=.false.,iorder=2)

c End program

      end subroutine findVectorPotential

c findScalarPotential
c####################################################################
      subroutine findScalarPotential(v,b,phi)

c--------------------------------------------------------------------
c     Finds COVARIANT vector potential from CONTRAVARIANT magnetic
c     field components for initialization of equilibrium.
c--------------------------------------------------------------------

      use parameters

      use grid

      use equilibrium

      use precond_variables

      implicit none

c Call variables

      real(8)    :: v(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,b(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,phi(ilom:ihip,jlom:jhip,klom:khip)

c Local variables

      integer(4) :: igrid,ii,bcnd(6,1)
      real(8)    :: phi_rhs(ntotdp),phiv(ntotdp),cov(3)
     .             ,dum(ilom:ihip,jlom:jhip,klom:khip,3),dvol

c Begin program

      igrid = 1

      nx = nxd
      ny = nyd
      nz = nzd

c Find rhs

      do k=klom,khip
        do j=jlom,jhip
          do i=ilom,ihip
            jac = gmetric%grid(igrid)%jac(i,j,k)

            cov(1) = v(i,j,k,2)*b(i,j,k,3)
     .             - v(i,j,k,3)*b(i,j,k,2)

            cov(2) = v(i,j,k,3)*b(i,j,k,1)
     .             - v(i,j,k,1)*b(i,j,k,3)

            cov(3) = v(i,j,k,1)*b(i,j,k,2)
     .             - v(i,j,k,2)*b(i,j,k,1)

            cov = cov/jac

            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,dum(i,j,k,1)
     .                                  ,dum(i,j,k,2)
     .                                  ,dum(i,j,k,3),.true.)

          enddo
        enddo
      enddo

      do k=klo,khi
        do j=jlo,jhi
          do i=ilo,ihi
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            if (vol_wgt) then
              dvol = volume(i,j,k,igrid,igrid,igrid)
            else
              dvol = 1d0
            endif

            phi_rhs(ii) = dvol*div(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,dum(:,:,:,1),dum(:,:,:,2),dum(:,:,:,3))
          enddo
        enddo
      enddo

c Solve for potential

      if (sqrt(sum(phi_rhs*phi_rhs))>0d0) then

#if defined(petsc)

        call pstop('findScalarPotential'
     .          ,'Poisson solver not implemented in parallel')

#else

        !Define BCs
        bcnd(:,1) = bcond
        where (bcnd == DEF) bcnd = DIR

        !Solve Poisson eq
        call poisson(ntotdp,phi_rhs,phiv,30,1d-6,bcnd,ilevel,0)

        !Map to array (w/o BCs)
        call mapMGVectorToArray(0,1,phiv,nx,ny,nz,phi,igrid,.false.)
#endif

      else

        phi = 0d0

      endif

c End program

      end subroutine findScalarPotential

