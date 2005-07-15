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
      real(8)       :: ldaflow,x1,y1,z1

      logical       :: covariant,to_cartsn,cartsn

      real(8)       :: r,jac1,bnorm,aaa,bbb,ccc,qq,qqp,q0,RR,ff,nn
     .                ,aspect_ratio,mm,kk,Iz,bz0,btheta,vr,rrho,r1,bb
     .                ,r0(1),aa,dummy,dummy1,dummy2,dummy3,ieta

      real(8)       :: rint(0:nxd+1),dh,offset(1),bth0(0:nxd+1)
     .                ,bzz0(0:nxd+1),vr0(0:nxd+1),rho0(0:nxd+1)
     .                ,eth0(0:nxd+1),ez0(0:nxd+1)

      real(8)       :: a1(ilom:ihip,jlom:jhip,klom:khip)
     .                ,a2(ilom:ihip,jlom:jhip,klom:khip)
     .                ,a3(ilom:ihip,jlom:jhip,klom:khip)

      character(200) :: command

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
      label(IBX)  = 'B^1'
      label(IBY)  = 'B^2'
      label(IBZ)  = 'B^3'
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

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = 1d0

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

              jac1 = gmetric%grid(igx)%jac(i,j,k)

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = jac1*1d0

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,ITMP) = 1d0
            enddo
          enddo
        enddo

      case ('khcar')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a1,a2,a3,igx,igy,igz)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Kelvin-Helmholtz with constant magnetic field in cartesian coordinates

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
             var(i,j,k,IVX:IVZ)=vperflow
     .                      *curl(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3)
            enddo
          enddo
        enddo

        var(:,:,:,IRHO) = 1d0
        var(:,:,:,IBX)  = 0d0
        var(:,:,:,IBY)  = 0d0
        var(:,:,:,IBZ)  = 1d0
        var(:,:,:,ITMP) = 1d0

      case ('3dkh')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a1,a2,a3,igx,igy,igz)

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
             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,1)
             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,2)
             var(i,j,k,IBZ)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,3)

             var(i,j,k,IRHO)= (2d0 - var(i,j,k,IBZ)**2)

             var(i,j,k,IVX) = vperflow*tanh((y1-.5*(ymax-ymin))/ldaflow)
     .                                *sin (2*pi*z1/(zmax-zmin))
             var(i,j,k,IVY) = 0d0
             var(i,j,k,IVZ) = 0d0

            enddo
          enddo
        enddo

        var(:,:,:,ITMP) = 0.25

      case ('tmcar')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a1,a2,a3,igx,igy,igz)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Tearing mode in cartesian coordinates

        bz0 = 1d0

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
             !X-Y equilibrium
             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,1)
             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,2)
             var(i,j,k,IBZ)=sqrt(bz0**2 - var(i,j,k,IBY)**2)

cc             var(i,j,k,IBX)=0d0
cc             var(i,j,k,IBY)=1d0
cc             var(i,j,k,IBZ)=0d0
            enddo
          enddo
        enddo

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0
cc        var(:,:,:,ITMP) = 0d0

      case ('tmcxz')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a1,a2,a3,igx,igy,igz)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Tearing mode in cartesian coordinates

        bz0 = 1d0

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
             !X-Z equilibrium
             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,1)
             var(i,j,k,IBZ)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,2)
             var(i,j,k,IBY)=sqrt(bz0**2 - var(i,j,k,IBZ)**2)
            enddo
          enddo
        enddo

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

      case ('tmcyz')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a1,a2,a3,igx,igy,igz)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Tearing mode in cartesian coordinates

        bz0 = 1d0

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
             !X-Z equilibrium
             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,2)
             var(i,j,k,IBZ)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,3)
             var(i,j,k,IBX)=sqrt(bz0**2 - var(i,j,k,IBZ)**2)
            enddo
          enddo
        enddo

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

      case ('tmsin')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a1,a2,a3,igx,igy,igz)

c     Check coordinates

        if (coords /= 'sin') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Tearing mode in sinusoidal coordinates

        bz0 = 1d0

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

              !X-Y equilibrium
             var(i,j,k,IBX)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,1)
             var(i,j,k,IBY)=curl2(i,j,k,nx,ny,nz,igx,igy,igz,a1,a2,a3,2)
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

c diag ****
cc             var(i,j,k,IBX)  = 0d0
cc             var(i,j,k,IBY)  = 0d0
cc             var(i,j,k,IBZ)  = 0d0
c diag ****

             var(i,j,k,IVX)  = 0d0
             var(i,j,k,IVY)  = 0d0
             var(i,j,k,IVZ)  = 0d0

             var(i,j,k,IRHO) = 1d0

             var(i,j,k,ITMP) = 1d0

             if (bbb**2+4*aaa*ccc < 0d0) then
               write (*,*) var(i,j,k,IBZ),bnorm,aaa,bbb**2+4*aaa*ccc
             endif
            enddo
          enddo
        enddo

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
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*x1/(1+(x1/dlambda)**2)
              bz0    = RR

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
              var(i,j,k,IBZ)  = x1*bz0
cc              var(i,j,k,IBX)  = 0d0
cc              var(i,j,k,IBY)  = 0d0
cc              var(i,j,k,IBZ)  = 0d0

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,ITMP) = 1d0
     .                        + 1/16./var(i,j,k,IRHO)
     .                         *(Iz*dlambda/(1+(x1/dlambda)**2))**2
cc              var(i,j,k,ITMP) = 1d-5
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
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*x1/(1+(x1/dlambda)**2)
              bz0    = RR

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta
              var(i,j,k,IBZ)  = x1*bz0

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

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
        open (unit=1000,file='pixie3d.equ',status='old')

        read(1000,'(a)')
        do i=0,nxd+1
          read(1000,'(9f12.8)') dummy1,bth0(i),bzz0(i),dummy2,dummy3
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
        J0(2) = 0.5*(eth0(nxd)+eth0(nxd+1))
        J0(3) = 0.5*(ez0 (nxd)+ez0 (nxd+1))

        !Convert electric field to covariant at r=1d0
        J0(3) = -kk/mm*J0(2) + J0(3)

        !Find current
        J0 = J0/ieta

        !Build equilibrium
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1

              !X-Y equilibrium
              var(i,j,k,IBX) = 0d0
              var(i,j,k,IBY) = bth0(iglobal) + kk*x1/mm*bzz0(iglobal)
              var(i,j,k,IBZ) = x1*bzz0(iglobal)

              if (nc_eom_v) then
                var(i,j,k,IVX) = x1*vr0(iglobal)
              else
                var(i,j,k,IVX) = x1*vr0(iglobal)*rho0(iglobal)
              endif
              var(i,j,k,IVX) = 0d0
              var(i,j,k,IVY) = 0d0
              var(i,j,k,IVZ) = 0d0

              var(i,j,k,IRHO) = rho0(iglobal)

              var(i,j,k,ITMP) = 1d-5

            enddo
          enddo
        enddo

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

cc        call B_from_q(bth0,bzz0)

        !Build equilibrium
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1

cc              btheta  = Iz/2./pi*x1*ff(1d0)/ff(x1)*exp(-rint(iglobal))
cc              bz0     = qq(x1)*btheta/x1
cc
cc              !X-Y equilibrium
cc              var(i,j,k,IBX)  = 0d0
cc              var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
cc              var(i,j,k,IBZ)  = x1*bz0

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = bth0(iglobal) + kk*x1/mm*bzz0(iglobal)
              var(i,j,k,IBZ)  = x1*bzz0(iglobal)

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,ITMP) = 1d-5

            enddo
          enddo
        enddo

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
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
              var(i,j,k,IBZ)  = x1*bz0

c diag ****
cc              var(i,j,k,IBX)  = 0d0
cc              var(i,j,k,IBY)  = 0d0
cc              var(i,j,k,IBZ)  = x1*0d0
c diag ****

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

c diag ****
cc              var(i,j,k,ITMP) = 1d0
cc              var(i,j,k,ITMP) = 0d0
c diag ****
              var(i,j,k,ITMP) = 1d-5

            enddo
          enddo
        enddo

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
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta
              var(i,j,k,IBZ)  = x1*bz0

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,ITMP) = 1d-5

            enddo
          enddo
        enddo

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
        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              !X-Y equilibrium
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
              var(i,j,k,IBZ)  = x1*bz0

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,ITMP) = 1d-3

            enddo
          enddo
        enddo

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
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = 0d0

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

c Transfer electric field info to postprocessor via input file

      if (sum(J0) /= 0d0) then

        command = 'grep J0 ' // trim(inputfile) // ' > /dev/null'
        ierr = system(trim(command))

        write (command,*) '   J0 =',real(J0(1),4)
     .                         ,',',real(J0(2),4)
     .                         ,',',real(J0(3),4)


        if (ierr /= 0) then  !Replace J0
          command='sed -e"/graphdef/{" -e' // achar(39)
     .       // 'a\' // achar(39) //' -e"'
     .       // trim(command)   // '" -e"}" '
     .       // trim(inputfile) // ' > temp'
        else                 !Add J0
          command='sed -e"/J0/{" -e' // achar(39)
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

        real(8)    :: bth0(0:nxd+1),bz0(0:nxd+1)

c     Local variables

        integer(4) :: i,j,k,ig,jg,kg
        real(8)    :: qq,qqp,ff,x1,dh,offset(1),r1,btheta,Iz,r0(1)

        real(8)    :: rint(0:nxd+1)

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

c     Evaluate magnetic field

        do i = 0,nxd+1

          x1 = grid_params%xg(i)  !Position in global grid

          bth0(i) = Iz/2./pi*x1*ff(1d0)/ff(x1)*exp(-rint(i))
          bz0 (i) = qq(x1)*bth0(i)/x1

        enddo

      end subroutine B_from_q

      end subroutine setEquilibrium

c fillVectorPotential
c####################################################################
      subroutine fillVectorPotential(a1,a2,a3,igx,igy,igz)

c--------------------------------------------------------------------
c     Defines COVARIANT vector potential for initialization of 
c     equilibrium.
c--------------------------------------------------------------------

      use parameters

      use grid

      use equilibrium

      implicit none

c Call variables

      integer(4) :: igx,igy,igz
      real(8)    :: a1(ilom:ihip,jlom:jhip,klom:khip)
     .             ,a2(ilom:ihip,jlom:jhip,klom:khip)
     .             ,a3(ilom:ihip,jlom:jhip,klom:khip)

c Local variables

      integer(4) :: ig,jg,kg,i,j,k
      real(8)    :: xx,yy,zz
      logical    :: cartsn

c Begin program

      select case (trim(equil))

      case ('khcar','tmcar','tmsin','tmcxz')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,xx,yy,zz)
              a1(i,j,k) = 0d0
              a2(i,j,k) = 0d0
              a3(i,j,k) = dlambda
     .             *dlog(dcosh((xx-0.5d0*(xmax-xmin))/dlambda)) 
            enddo
          enddo
        enddo

      case ('tmcyz','3dkh')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,xx,yy,zz)
              a1(i,j,k) = dlambda
     .             *dlog(dcosh((yy-0.5d0*(ymax-ymin))/dlambda)) 
              a2(i,j,k) = 0d0
              a3(i,j,k) = 0d0
            enddo
          enddo
        enddo

      end select

c End program

      end subroutine fillVectorPotential

