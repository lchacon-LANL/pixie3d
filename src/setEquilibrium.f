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

      implicit none

c Call variables

      real(8)       :: var(0:nxdp,0:nydp,0:nzdp,neqd)

      character*(5) :: label(neqd)

      integer(4)    :: bcs(6,neqd)

c Local variables

      integer(4)    :: i,j,k,ig,jg,kg,ieq,igx,igy,igz
      real(8)       :: ldaflow,xx,yy,zz,car(3)
      real(8)       :: bx,by,bz

      logical       :: covariant,to_cartesian

      real(8)       :: bz0

c Externals

      real(8)       :: Ax,Ay,Az
      external      :: Ax,Ay,Az

c Begin program

      ldaflow = dlambda/rshear

      var = 0d0

      igx = 1
      igy = 1
      igz = 1

c Label variables
     
      label(IRHO) = 'Rho'
      label(IBX)  = 'Bx'
      label(IBY)  = 'By'
      label(IBZ)  = 'Bz'
      label(IVX)  = 'nVx'
      label(IVY)  = 'nVy'
      label(IVZ)  = 'nVz'
      label(ITMP) = 'Temp'

c Define boundary conditions

      call defineBoundaryConditions(neqd,bcs)

c Set initial guess (in Cartesian coordinates)

      select case (trim(equil))

      case ('kh')

c     Kelvin-Helmholtz with constant magnetic field

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd
              var(i,j,k,IRHO) = 1d0

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = 1d0

              call curl(i,j,k,igx,igy,igz,Ax,Ay,Az,bx,by,bz)
              var(i,j,k,IVX)  = vperflow*bx
              var(i,j,k,IVY)  = vperflow*by
              var(i,j,k,IVZ)  = vperflow*bz

              var(i,j,k,ITMP) = 1d0
            enddo
          enddo
        enddo

      case ('tm')

        bz0 = 1d0

c     Kelvin-Helmholtz/Tearing mode

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd

              call curl(i,j,k,igx,igy,igz,Ax,Ay,Az,bx,by,bz)
              var(i,j,k,IBX)  = bx
              var(i,j,k,IBY)  = by
cc              var(i,j,k,IBX)  = 1d0
cc              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = sqrt(bz0**2 - bx**2)

              var(i,j,k,IVX)  = 0d0
              var(i,j,k,IVY)  = 0d0
              var(i,j,k,IVZ)  = 0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,ITMP) = 1d0

            enddo
          enddo
        enddo

cc      case ('alfwv')
cc
ccc     Alfven wave test (perturb PSI)
cc
cc      do i = 0,nxd+1
cc        do j = 0,nyd+1
cc          var(i,j,PSI) = - yl(j,ngrd)
cc          var(i,j,PHI) = - vperflow*yl(j,ngrd)
cc          var(i,j,IVZ) = 0d0
cc          var(i,j,IBZ) = 0d0
cc        enddo
cc      enddo
cc
cc      case ('flxbd')
cc
ccc     Two flux bundles (no perturbation and no source required)
cc
cc      x1 = xlength/2. - 1.3*dlambda
cc      x2 = xlength/2. + 1.3*dlambda
cc      sigma = dlambda**4
cc
cc      bmax = (108./sigma)**.25*exp(-.75)
cccc      jmax = 4.*exp(-.25)/sqrt(sigma)
cc
cc      do j = 0,nyd+1
cc        do i = 0,nxd+1
cc          var(i,j,PSI) = 1./bmax*
cc     .     (exp( -((xl(i,ngrd)-x1)**2+(yl(j,ngrd)-.5)**2)**2/sigma)
cc     .     +exp( -((xl(i,ngrd)-x2)**2+(yl(j,ngrd)-.5)**2)**2/sigma))
cc          var(i,j,PHI) = - vperflow*yl(j,ngrd)
cc          var(i,j,IVZ) = vparflow*yl(j,ngrd)
cc          var(i,j,IBZ) = 0d0
cc        enddo
cc      enddo
cc
cc      source = .false.

      case default

        write (*,*)
        write (*,*) 'Equilibrium ',trim(equil),' undefined'
        write (*,*) 'Aborting...'
        stop

      end select

c Transform vectors to curvilinear coordinates

      to_cartesian = .false.
      covariant    = .false.

      !Velocity

      call transformVector(igx,igy,igz,0,nxdp,0,nydp,0,nzdp
     .                    ,var(:,:,:,IVX)
     .                    ,var(:,:,:,IVY)
     .                    ,var(:,:,:,IVZ)
     .                    ,covariant,to_cartesian)

      var(:,:,:,IVX) = var(:,:,:,IRHO)*var(:,:,:,IVX)
      var(:,:,:,IVY) = var(:,:,:,IRHO)*var(:,:,:,IVY)
      var(:,:,:,IVZ) = var(:,:,:,IRHO)*var(:,:,:,IVZ)

      !Magnetic field

      call transformVector(igx,igy,igz,0,nxdp,0,nydp,0,nzdp
     .                    ,var(:,:,:,IBX)
     .                    ,var(:,:,:,IBY)
     .                    ,var(:,:,:,IBZ)
     .                    ,covariant,to_cartesian)

c End program

      end subroutine setEquilibrium

c curl
c####################################################################
      subroutine curl(i,j,k,igx,igy,igz,ax,ay,az,bx,by,bz)

c--------------------------------------------------------------------
c     Calculates curl of vector (ax,ay,az) at position i,j,k, and 
c     returns curl components in (bx,by,bz). The input
c     vector components are given by external functions.
c--------------------------------------------------------------------

      use grid

      implicit none

c Call variables

      integer(4) :: i,j,k,igx,igy,igz
      real(8)    :: bx,by,bz
      real(8)    :: ax,ay,az
      external   :: ax,ay,az

c Local variables

      logical    :: sing_point
      integer(4) :: ig,jg,kg
      real(8)    :: dx,dy,dz,jac
      real(8)    :: Ez_ipjp,Ez_imjp,Ez_ipjm,Ez_imjm
     .             ,Ey_ipkp,Ey_imkp,Ey_ipkm,Ey_imkm
     .             ,Ex_jpkp,Ex_jmkp,Ex_jpkm,Ex_jmkm
      real(8)    :: Ez_jp,Ez_jm,Ey_kp,Ey_km,Ex_kp,Ex_km
     .             ,Ez_ip,Ez_im,Ey_ip,Ey_im,Ex_jp,Ex_jm

c Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dx = grid_params%dxh(ig)
      dy = grid_params%dyh(jg)
      dz = grid_params%dzh(kg)

      jac = jacobian(grid_params%xx(ig)
     .              ,grid_params%yy(jg)
     .              ,grid_params%zz(kg))

      sing_point = .false.
      if (jac == 0d0) sing_point = .true.

c Find vector components at edges

cc      if (.not.sing_point) then
        Ex_jpkp = ax(i,j  ,k  ,igx,igy,igz)
        Ex_jmkp = ax(i,j-1,k  ,igx,igy,igz)
        Ex_jpkm = ax(i,j  ,k-1,igx,igy,igz)
        Ex_jmkm = ax(i,j-1,k-1,igx,igy,igz)

        Ey_ipkp = ay(i  ,j,k  ,igx,igy,igz)
        Ey_imkp = ay(i-1,j,k  ,igx,igy,igz)
        Ey_ipkm = ay(i  ,j,k-1,igx,igy,igz)
        Ey_imkm = ay(i-1,j,k-1,igx,igy,igz)

        Ez_ipjp = az(i  ,j  ,k,igx,igy,igz)
        Ez_imjp = az(i-1,j  ,k,igx,igy,igz)
        Ez_ipjm = az(i  ,j-1,k,igx,igy,igz)
        Ez_imjm = az(i-1,j-1,k,igx,igy,igz)
cc      else
cc        Ex_jpkp = ax(i,1,k  ,igx,igy,igz)
cc        Ex_jmkp = ax(i,0,k  ,igx,igy,igz)
cc        Ex_jpkm = ax(i,1,k-1,igx,igy,igz)
cc        Ex_jmkm = ax(i,0,k-1,igx,igy,igz)
cc
cc        Ez_ipjp = az(i,1,k,igx,igy,igz)
cc        Ez_ipjm = az(i,0,k,igx,igy,igz)
cc
cc        Ey_ipkp = ay(i,1,k  ,igx,igy,igz)
cc        Ey_ipkm = ay(i,1,k-1,igx,igy,igz)
cc      endif

c     Bx

cc      if (.not.sing_point) then
        Ez_jp = 0.5*(Ez_ipjp + Ez_imjp)
        Ez_jm = 0.5*(Ez_ipjm + Ez_imjm)

        Ey_kp = 0.5*(Ey_ipkp + Ey_imkp)
        Ey_km = 0.5*(Ey_ipkm + Ey_imkm)

        bx = ( Ez_jp - Ez_jm )/dy - ( Ey_kp - Ey_km )/dz
cc      else
cc
cc        Ez_ip = 0.5*(Ez_ipjp + Ez_ipjm)
cc
cc        Ey_kp = Ey_ipkp
cc        Ey_km = Ey_ipkm
cc
cc        bx = ( Ez_jp - Ez_jm )/dy - ( Ey_kp - Ey_km )/dz
cc      endif

c     By

      Ex_kp = 0.5*(Ex_jpkp + Ex_jmkp)
      Ex_km = 0.5*(Ex_jpkm + Ex_jmkm)

      Ez_ip = 0.5*(Ez_ipjp + Ez_ipjm)
      Ez_im = 0.5*(Ez_imjp + Ez_imjm)

      by = ( Ex_kp - Ex_km )/dz - ( Ez_ip - Ez_im )/dx

c     Bz

      Ey_ip = 0.5*(Ey_ipkp + Ey_ipkm)
      Ey_im = 0.5*(Ey_imkp + Ey_imkm)

      Ex_jp = 0.5*(Ex_jpkp + Ex_jpkm)
      Ex_jm = 0.5*(Ex_jmkp + Ex_jmkm)

      bz = ( Ey_ip - Ey_im )/dx - ( Ex_jp - Ex_jm )/dy

c End program

      end subroutine curl

c Ax
c####################################################################
      real(8) function Ax(i,j,k,igx,igy,igz)

c--------------------------------------------------------------------
c     Gives X-component of vector potential at (j+1/2,k+1/2) edges.
c--------------------------------------------------------------------

      use grid

      use equilibrium

      implicit none

c Call variables

      integer(4) :: i,j,k,igx,igy,igz

c Local variables

      integer(4) :: ig,jg,kg
      real(8)    :: xx,yy,zz,car(3)

c Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      xx = grid_params%xx(ig)
      yy = 0.5*(grid_params%yy(jg) + grid_params%yy(jg+1))
      zz = 0.5*(grid_params%zz(kg) + grid_params%zz(kg+1))

      car = inverse_map(xx,yy,zz)

      xx = car(1)
      yy = car(2)
      zz = car(3)

c Find vector components at edges

      select case (trim(equil))

      case default

        ax = 0d0

      end select

c End program

      end function Ax

c Ay
c####################################################################
      real(8) function Ay(i,j,k,igx,igy,igz)

c--------------------------------------------------------------------
c     Gives X-component of vector potential at (i+1/2,k+1/2) edges.
c--------------------------------------------------------------------

      use grid

      use equilibrium

      implicit none

c Call variables

      integer(4) :: i,j,k,igx,igy,igz

c Local variables

      integer(4) :: ig,jg,kg
      real(8)    :: xx,yy,zz,car(3)

c Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      xx = 0.5*(grid_params%xx(ig) + grid_params%xx(ig+1))
      yy =      grid_params%yy(jg)
      zz = 0.5*(grid_params%zz(kg) + grid_params%zz(kg+1))

      car = inverse_map(xx,yy,zz)

      xx = car(1)
      yy = car(2)
      zz = car(3)

c Find vector components at edges

      select case (trim(equil))

      case default

        ay = 0d0

      end select

c End program

      end function Ay

c Az
c####################################################################
      real(8) function Az(i,j,k,igx,igy,igz)

c--------------------------------------------------------------------
c     Gives X-component of vector potential at (i+1/2,j+1/2) edges.
c--------------------------------------------------------------------

      use grid

      use equilibrium

      implicit none

c Call variables

      integer(4) :: i,j,k,igx,igy,igz

c Local variables

      integer(4) :: ig,jg,kg
      real(8)    :: xx,yy,zz,car(3)

c Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      xx = 0.5*(grid_params%xx(ig) + grid_params%xx(ig+1))
      yy = 0.5*(grid_params%yy(jg) + grid_params%yy(jg+1))
      zz =      grid_params%zz(kg)

cc      car = inverse_map(xx,yy,zz)
cc
cc      xx = car(1)
cc      yy = car(2)
cc      zz = car(3)

c Find vector components at edges

      select case (trim(equil))

      case default

        select case (coords)
        case ('car')

          az =  dlambda*dlog(dcosh((yy-0.5d0)/dlambda)) 

        case ('scl')

          az =  tanh(0.5/dlambda)*(yy**2/ylength-yy) 

        case default

          write (*,*) 'Vector potential not implemented'
          write (*,*) 'Aborting...'
          stop

        end select

      end select

c End program

      end function Az

c transformVector
c######################################################################
      subroutine transformVector(igx,igy,igz
     .                          ,imin,imax,jmin,jmax,kmin,kmax
     .                          ,arr1,arr2,arr3,covariant,to_cartesian)

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
        logical    :: covariant,to_cartesian
        real(8)    :: arr1(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr2(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr3(imin:imax,jmin:jmax,kmin:kmax)

c     Local variables

        integer(4) :: i,j,k,ig,jg,kg
        real(8)    :: vec(3),xx,yy,zz

c     Begin program

        if (to_cartesian) then

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

                xx = grid_params%xx(ig)
                yy = grid_params%yy(jg)
                zz = grid_params%zz(kg)

                call transformVectorToCartesian
     .               (xx,yy,zz
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

                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

                xx = grid_params%xx(ig)
                yy = grid_params%yy(jg)
                zz = grid_params%zz(kg)

                call transformVectorToCurvilinear
     .               (xx,yy,zz
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
