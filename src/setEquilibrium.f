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

      integer(4)    :: i,j,k,ig,jg,kg,ieq
      real(8)       :: ldaflow,x1,y1,z1

      logical       :: covariant,to_cartsn,cartsn

      real(8)       :: bz0,r,jac1,bnorm,aaa,bbb,ccc

      real(8)       :: a1(0:nxdp,0:nydp,0:nzdp)
     .                ,a2(0:nxdp,0:nydp,0:nzdp)
     .                ,a3(0:nxdp,0:nydp,0:nzdp)

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
      label(IBX)  = 'Bx (cnv)'
      label(IBY)  = 'By (cnv)'
      label(IBZ)  = 'Bz (cnv)'
      label(IVX)  = 'Px (cnv)'
      label(IVY)  = 'Py (cnv)'
      label(IVZ)  = 'Pz (cnv)'
      label(ITMP) = 'Temp'

c Define boundary conditions

      call defineBoundaryConditions(neqd,bcs)

c Define vector potential (in curvilinear coordinates) for initialization

      call fillVectorPotential(a1,a2,a3,igx,igy,igz)

c Set initial guess

      select case (trim(equil))

      case ('khcar')

        coords = 'car'

c     Kelvin-Helmholtz with constant magnetic field in cartesian coordinates

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd
              var(i,j,k,IRHO) = 1d0

              var(i,j,k,IBX)=0d0
              var(i,j,k,IBY)=0d0
              var(i,j,k,IBZ)=1d0

              var(i,j,k,IVX)=vperflow*curl(i,j,k,a1,a2,a3,1)
              var(i,j,k,IVY)=vperflow*curl(i,j,k,a1,a2,a3,2)
              var(i,j,k,IVZ)=vperflow*curl(i,j,k,a1,a2,a3,3)

              var(i,j,k,ITMP) = 1d0
            enddo
          enddo
        enddo

      case ('tmcar')

c     Tearing mode in cartesian coordinates

        coords = 'car'

        bz0 = 1d0

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd

              !X-Y equilibrium
              var(i,j,k,IBX)  = curl(i,j,k,a1,a2,a3,1)
              var(i,j,k,IBY)  = curl(i,j,k,a1,a2,a3,2)
              var(i,j,k,IBZ)  = sqrt(bz0**2 - var(i,j,k,IBY)**2)

              !X-Z equilibrium
cc              var(i,j,k,IBX)  = curl(i,j,k,a1,a2,a3,1)
cc              var(i,j,k,IBZ)  = curl(i,j,k,a1,a2,a3,2)
cc              var(i,j,k,IBY)  = sqrt(bz0**2 - var(i,j,k,IBZ)**2)

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,ITMP) = 1d0

            enddo
          enddo
        enddo

      case ('dfcyl')

c     Dipolar flow in cylindrical coordinates

        coords = 'cyl'

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

              var(i,j,k,IBX)=curl(i,j,k,a1,a2,a3,1)
              var(i,j,k,IBY)=curl(i,j,k,a1,a2,a3,2)
              var(i,j,k,IBZ)=curl(i,j,k,a1,a2,a3,3)

              var(i,j,k,IVX)=vperflow*curl(i,j,k,a1,a2,a3,1)
              var(i,j,k,IVY)=vperflow*curl(i,j,k,a1,a2,a3,2)
              var(i,j,k,IVZ)=vperflow*curl(i,j,k,a1,a2,a3,3)

              var(i,j,k,ITMP) = 1d0

            enddo
          enddo
        enddo

      case ('tmsin')

c     Tearing mode in sinusoidal coordinates

        coords = 'sin'

        bz0 = 1d0

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd

              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                           ,cartsn)

              !X-Y equilibrium
              var(i,j,k,IBX)  = curl(i,j,k,a1,a2,a3,1)
              var(i,j,k,IBY)  = curl(i,j,k,a1,a2,a3,2)
              var(i,j,k,IBZ)  = 0d0

              gsub = G_sub   (x1,y1,z1,cartsn)
              jac1 = jacobian(x1,y1,z1,cartsn)
              bnorm= vectorNorm(x1,y1,z1,var(i,j,k,IBX),var(i,j,k,IBY)
     .                         ,var(i,j,k,IBZ),.false.,cartsn)

              ccc = jac1*(bz0**2 - bnorm)
              bbb = gsub(3,2)*var(i,j,k,IBY) + gsub(3,1)*var(i,j,k,IBX)
              aaa = gsub(3,3)

              var(i,j,k,IBZ)  = (-bbb+sqrt(bbb**2+4*aaa*ccc))/2./aaa

              !X-Z equilibrium
cc              var(i,j,k,IBX)  = curl(i,j,k,a1,a2,a3,1)
cc              var(i,j,k,IBZ)  = curl(i,j,k,a1,a2,a3,2)
cc              var(i,j,k,IBY)  = sqrt(bz0**2 - var(i,j,k,IBZ)**2)

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,ITMP) = 1d0

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

              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,xx,yy,zz
     .                           ,cartsn)
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              xx = grid_params%xx(ig)
cc              yy = grid_params%yy(jg)
cc              zz = grid_params%zz(kg)

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
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              xx = grid_params%xx(ig)
              yy = grid_params%yy(jg)
              zz = grid_params%zz(kg)

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

      case default

        write (*,*)
        write (*,*) 'Equilibrium ',trim(equil),' undefined'
        write (*,*) 'Aborting...'
        stop

      end select

c End program

      end subroutine fillVectorPotential

c transformVector
c######################################################################
      subroutine transformVector(igx,igy,igz
     .                          ,imin,imax,jmin,jmax,kmin,kmax
     .                          ,arr1,arr2,arr3,covariant,to_cartsn)

c----------------------------------------------------------------------
c     Transforms vectors components in arrays arr1,arr2,arr3
c     from Cartesian to curvilinear (to_cartesian=.false.)
c     and viceversa, either with covariant (covariant=.true.) or 
c     contravariant curvilinear vectors.
c----------------------------------------------------------------------

      use grid

      implicit none

c     Input variables

        integer(4) :: imin,imax,jmin,jmax,kmin,kmax
        integer(4) :: igx,igy,igz
        logical    :: covariant,to_cartsn
        real(8)    :: arr1(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr2(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr3(imin:imax,jmin:jmax,kmin:kmax)

c     Local variables

        integer(4) :: i,j,k
        real(8)    :: vec(3)

c     Begin program

        if (to_cartsn) then

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call transformVectorToCartesian
     .               (i,j,k,igx,igy,igz
     .               ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .               ,covariant
     .               ,vec(1),vec(2),vec(3))

                arr1(i,j,k) = vec(1)
                arr2(i,j,k) = vec(2)
                arr3(i,j,k) = vec(3)
                
              enddo
            enddo
          enddo

        else

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call transformVectorToCurvilinear
     .               (i,j,k,igx,igy,igz
     .               ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .               ,covariant
     .               ,vec(1),vec(2),vec(3))

                arr1(i,j,k) = vec(1)
                arr2(i,j,k) = vec(2)
                arr3(i,j,k) = vec(3)
                
              enddo
            enddo
          enddo

        endif

      end subroutine transformVector
