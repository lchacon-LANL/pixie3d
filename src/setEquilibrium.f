c setEquilibrium
c####################################################################
      subroutine setEquilibrium(var,bcs,label)

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

      implicit none

c Call variables

      real(8)       :: var(0:nxdp,0:nydp,0:nzdp,neqd)

      character*(20):: label(neqd)

      integer(4)    :: bcs(6,neqd)

c Local variables

      integer(4)    :: i,j,k,ig,jg,kg,ieq,nmax
      real(8)       :: ldaflow,x1,y1,z1

      logical       :: covariant,to_cartsn,cartsn

      real(8)       :: bz0,r,jac1,bnorm,aaa,bbb,ccc,qq,qqp,q0,rr,ff
     .                ,aspect_ratio,mm,kk,Iz,btheta,bzz,rint(0:nxd+1)
     .                ,r1,bb,aa,nn

      real(8)       :: a1(0:nxdp,0:nydp,0:nzdp)
     .                ,a2(0:nxdp,0:nydp,0:nzdp)
     .                ,a3(0:nxdp,0:nydp,0:nzdp)

c Functions

      !q-profile
      qq (rr) = 0.6125*(1 - 1.8748*rr**2 + 0.8323*rr**4)
      qqp(rr) = 0.6125*( -2*1.8748*rr   +4*0.8323*rr**3)
cc      qq(rr)  = 0.3*(1 - 1.8748*rr**2 + 0.8323*rr**4)
cc      qqp(rr) = 0.3*( -2*1.8748*rr   +4*0.8323*rr**3)

      ff(rr) = rr**2 + qq(rr)**2

c Begin program

      ldaflow = dlambda/rshear

      var = 0d0

c Initialize required grid information

      igx = 1
      igy = 1
      igz = 1

      nx = nxd
      ny = nyd
      nz = nzd

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

c Decide whether to use the alternative EOM formulation

      alt_eom = .false.
      if (bcond(1) == SP) alt_eom=.true.
      
c Set initial guess

      select case (trim(equil))

      case ('msw') !Magnetosonic wave

        gamma = 1d0

c     Check coordinates

        if (coords /= 'car') then
          write (*,*) 'Wrong coordinates for equilibrium ',equil
          write (*,*) 'Aborting...'
          stop
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd
cc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)

              var(i,j,k,IRHO) = 1d0
cc     .                         + 1d-3*cos(2*pi*(x1-xmin)/(xmax-xmin)
cc     .                                   +2*pi*(y1-ymin)/(ymax-ymin))

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = 1d0
cc     .                         + 1d-3*cos(2*pi*(x1-xmin)/(xmax-xmin)
cc     .                                   +2*pi*(y1-ymin)/(ymax-ymin))

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
          write (*,*) 'Wrong coordinates for equilibrium ',equil
          write (*,*) 'Aborting...'
          stop
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd
cc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)

              jac1 = gmetric%grid(igx)%jac(i,j,k)

              var(i,j,k,IRHO) = 1d0
cc     .                        + 1d-3*cos(2*pi*(x1-xmin)/(xmax-xmin)
cc     .                                  +2*pi*(y1-ymin)/(ymax-ymin))

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = jac1*1d0
cc     .                        + jac1*1d-3*cos(2*pi*(x1-xmin)/(xmax-xmin)
cc     .                                       +2*pi*(y1-ymin)/(ymax-ymin))

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
          write (*,*) 'Wrong coordinates for equilibrium ',equil
          write (*,*) 'Aborting...'
          stop
        endif

c     Kelvin-Helmholtz with constant magnetic field in cartesian coordinates

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd

              var(i,j,k,IVX)=vperflow*curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,1)
              var(i,j,k,IVY)=vperflow*curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,2)
              var(i,j,k,IVZ)=vperflow*curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,3)

            enddo
          enddo
        enddo

        var(:,:,:,IRHO) = 1d0
        var(:,:,:,IBX)  = 0d0
        var(:,:,:,IBY)  = 0d0
        var(:,:,:,IBZ)  = 1d0
        var(:,:,:,ITMP) = 1d0

      case ('tmcar')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a1,a2,a3,igx,igy,igz)

c     Check coordinates

        if (coords /= 'car') then
          write (*,*) 'Wrong coordinates for equilibrium ',equil
          write (*,*) 'Aborting...'
          stop
        endif

c     Tearing mode in cartesian coordinates

        bz0 = 1d0

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd

              !X-Y equilibrium
              var(i,j,k,IBX)  = curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,1)
              var(i,j,k,IBY)  = curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,2)
              var(i,j,k,IBZ)  = sqrt(bz0**2 - var(i,j,k,IBY)**2)

              !X-Z equilibrium
cc              var(i,j,k,IBX)  = curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,1)
cc              var(i,j,k,IBZ)  = curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,2)
cc              var(i,j,k,IBY)  = sqrt(bz0**2 - var(i,j,k,IBZ)**2)

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
          write (*,*) 'Wrong coordinates for equilibrium ',equil
          write (*,*) 'Aborting...'
          stop
        endif

c     Tearing mode in sinusoidal coordinates

        bz0 = 1d0

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd

              !X-Y equilibrium
              var(i,j,k,IBX)  = curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,1)
              var(i,j,k,IBY)  = curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,2)
              var(i,j,k,IBZ)  = 0d0

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

              !X-Z equilibrium
cc              var(i,j,k,IBX)  = curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,1)
cc              var(i,j,k,IBZ)  = curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,2)
cc              var(i,j,k,IBY)  = sqrt(bz0**2 - var(i,j,k,IBZ)**2)

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

      case ('dfcyl')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(a1,a2,a3,igx,igy,igz)

c     Check coordinates

        if (coords /= 'cyl') then
          write (*,*) 'Wrong coordinates for equilibrium ',equil
          write (*,*) 'Aborting...'
          stop
        endif

c     Dipolar flow in cylindrical coordinates

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
              r = grid_params%xx(ig)

              var(i,j,k,IRHO) = 1d0
cc              var(i,j,k,IRHO) = 1.+0.1*exp(-(r/0.25)**2)
cc              var(i,j,k,IRHO) = 1d0+0.5*sin(grid_params%xx(ig))

cc              var(i,j,k,IBX)=0d0
cc              var(i,j,k,IBY)=0d0
cc              var(i,j,k,IBZ)= r

              var(i,j,k,IBX)=curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,1)
              var(i,j,k,IBY)=curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,2)
              var(i,j,k,IBZ)=curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,3)

              var(i,j,k,IVX)=vperflow*curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,1)
              var(i,j,k,IVY)=vperflow*curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,2)
              var(i,j,k,IVZ)=vperflow*curl(i,j,k,nxd,nyd,nzd,a1,a2,a3,3)

              var(i,j,k,ITMP) = 1d0

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
          write (*,*) 'Wrong coordinates for equilibrium ',equil
          write (*,*) 'Aborting...'
          stop
        endif

c     RFP equilibria (Caramana et al, PoP, 1983)

        !Integral in r
        rint(nxd+1) = 0d0
        do i=nxd,1,-1
          call getCurvilinearCoordinates(i,1,1,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)
c1          rint(i) = rint(i+1) + dxh(ig)*x1/ff(x1)
          rint(i) = rint(i+1) + dxh(ig)*qq(x1)*qqp(x1)/ff(x1)
        enddo
        rint(0) = rint(1)

        !Determine Iz so that Bz(r=0)=1
c1        btheta = 1./2./pi*sqrt(ff(1d0)/ff(0d0))*exp(rint(0))
        btheta = 1./2./pi*ff(1d0)/ff(0d0)*exp(-rint(0))
        Iz = 1./btheta/qq(0d0)

        !Build equilibrium
        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd+1

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

c1              btheta = Iz/2./pi*x1*sqrt(ff(1d0)/ff(x1))*exp(rint(i))
              btheta = Iz/2./pi*x1*ff(1d0)/ff(x1)*exp(-rint(i))
              bzz    = qq(x1)*btheta/x1

              !X-Y equilibrium
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta + kk*x1/mm*bzz
              var(i,j,k,IBZ)  = x1*bzz

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,ITMP) = 1d-3

            enddo
          enddo
        enddo

      case ('rfp2')

c     Check coordinates

        if (coords /= 'hel') then
          write (*,*) 'Wrong coordinates for equilibrium ',equil
          write (*,*) 'Aborting...'
          stop
        endif

c     RFP equilibria (analytical)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        aa = grid_params%params(4)

        nh2 = mm !To set the right perturbation wavelength

        !Build equilibrium
        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd+1

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bzz    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta + kk*x1/mm*bzz
              var(i,j,k,IBZ)  = x1*bzz

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

cc              var(i,j,k,ITMP) = 1d0
              var(i,j,k,ITMP) = 1d-5

            enddo
          enddo
        enddo

      case ('rfp3d')

c     Check coordinates

        if (coords /= 'cyl') then
          write (*,*) 'Wrong coordinates for equilibrium ',equil
          write (*,*) 'Aborting...'
          stop
        endif

c     RFP equilibria (analytical)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        aa = grid_params%params(4)

        nh2 = mm       !To set the right perturbation wavelength
        nh3 = kk*RR    !To set the right perturbation wavelength

        !Build equilibrium
        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bzz    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta
              var(i,j,k,IBZ)  = x1*bzz

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
          write (*,*) 'Wrong coordinates for equilibrium ',equil
          write (*,*) 'Aborting...'
          stop
        endif

        !Build equilibrium
        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd+1

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bzz    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              !X-Y equilibrium
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta + kk*x1/mm*bzz
              var(i,j,k,IBZ)  = x1*bzz

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,ITMP) = 1d-3

            enddo
          enddo
        enddo

      case default

        write (*,*)
        write (*,*) 'Equilibrium ',trim(equil),' undefined'
        write (*,*) 'Aborting...'
        stop

      end select

c Transform vectors to curvilinear coordinates

cc      to_cartsn = .false.
cc      covariant    = .false.

      !Velocity

cc      call transformVector(igx,igy,igz,0,nxdp,0,nydp,0,nzdp
cc     .                    ,var(:,:,:,IVX)
cc     .                    ,var(:,:,:,IVY)
cc     .                    ,var(:,:,:,IVZ)
cc     .                    ,covariant,to_cartsn)

      !Magnetic field

cc      call transformVector(igx,igy,igz,0,nxdp,0,nydp,0,nzdp
cc     .                    ,var(:,:,:,IBX)
cc     .                    ,var(:,:,:,IBY)
cc     .                    ,var(:,:,:,IBZ)
cc     .                    ,covariant,to_cartsn)

c Find momentum components

      var(:,:,:,IVX) = var(:,:,:,IRHO)*var(:,:,:,IVX)
      var(:,:,:,IVY) = var(:,:,:,IRHO)*var(:,:,:,IVY)
      var(:,:,:,IVZ) = var(:,:,:,IRHO)*var(:,:,:,IVZ)

c End program

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
      real(8)    :: a1(0:nxdp,0:nydp,0:nzdp)
     .             ,a2(0:nxdp,0:nydp,0:nzdp)
     .             ,a3(0:nxdp,0:nydp,0:nzdp)

c Local variables

      integer(4) :: ig,jg,kg,i,j,k
      real(8)    :: xx,yy,zz
      logical    :: cartsn

c Begin program

      select case (trim(equil))

      case ('khcar','tmcar','tmsin')

        do k=0,nzdp
          do j=0,nydp
            do i=0,nxdp

              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,xx,yy,zz)
              a1(i,j,k) = 0d0
              a2(i,j,k) = 0d0
              a3(i,j,k) = dlambda*dlog(dcosh((xx-0.5d0)/dlambda)) 
            enddo
          enddo
        enddo

      case ('dfcyl')

        do k=0,nzdp
          do j=0,nydp
            do i=0,nxdp
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,xx,yy,zz)

              a1(i,j,k) = 0d0
              a2(i,j,k) = 0d0
cc              a3(i,j,k) = sin(yy)*(xx**2-1.)**3
              a3(i,j,k) = sin(yy)*(xx**2-1.)**3*xx**2
cc              a3(i,j,k) = -xx*sin(yy)
cc              a3(i,j,k) = (xx**2-1.)**3*xx**2
cc              a3(i,j,k) = xx*sin(yy)*(1.-1./xx**2)
cc              a3(i,j,k) = 0d0
            enddo
          enddo
        enddo

      end select

c End program

      end subroutine fillVectorPotential

