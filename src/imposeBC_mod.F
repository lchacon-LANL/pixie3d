c module BCS
c####################################################################
      module BCS

        use grid

        use grid_aliases

        use auxiliaryVariables

        use operators

        use icond

        use equilibrium

        use transport_params

        use constants

        use variables

        integer(4) :: nnvar,imax,imin,jmax,jmin,kmax,kmin
        integer(4) :: iimax,iimin,jjmax,jjmin,kkmax,kkmin

        real(8),allocatable,dimension(:,:) :: rhs

        real(8),allocatable,dimension(:,:,:,:) :: v_cov,v_cnv,v0

        type :: curv_vector
          real(8),pointer,dimension(:,:,:,:) :: vcomp
          logical :: cov
        end type curv_vector

      contains

c     findLoopLimits
c     ##############################################################
      subroutine findLoopLimits(dim,loc,i1,i2,j1,j2,k1,k2)
c     --------------------------------------------------------------
c     Finds limits for BC loops. On input, it takes (i1,i2),
c     (j1,j2), (k1,k2) and outputs loop limits (imin,imax),
c     (jmin,jmax), (kmin,kmax) (defined in BCS module)
c     --------------------------------------------------------------

c     Call variables

      integer(4),intent(IN) :: dim,loc,i1,i2,j1,j2,k1,k2

c     Begin program

      select case(dim)
      case (1)
        imin = 1 + loc*(nx-1)
        imax = imin

        jmin = j1
        jmax = j2

        kmin = k1
        kmax = k2

      case(2)
        imin = i1
        imax = i2

        jmin = 1 + loc*(ny-1)
        jmax = jmin

        kmin = k1
        kmax = k2

      case(3)
        imin = i1
        imax = i2

        jmin = j1
        jmax = j2

        kmin = 1 + loc*(nz-1)
        kmax = kmin

      end select

c     End program

      end subroutine findLoopLimits

      end module BCS

c module singularBCinterface
c####################################################################
      module singularBCinterface

        use BCS

        INTERFACE singularBC
          module procedure scalarSingularBC,vectorSingularBC
        end INTERFACE

      contains

c     scalarSingularBC
c     #################################################################
      subroutine scalarSingularBC(array,order)
c     -----------------------------------------------------------------
c     Imposes singular point BC. On input:
c        * array: contains variable on which singular BC is imposed
c        * order: order of interpolation towards singular point
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: order
      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,order1
      real(8)    :: avg_q,avg_vol,rho0,vol,x0

c     Begin program

      if (order == 3) then
        order1 = 3
      else
        order1 = order
      endif

      do k=1,nz
        avg_q   = 0d0
        avg_vol = 0d0
        do j=1,ny
          i = 1
          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
          vol = volume(i,j,k,igx,igy,igz)
          x0  = xx(ig-1)

          call IntDriver1d(order1+1,xx(ig),array(i:i+order1+1,j,k)
     .                    ,1,x0,rho0,order)

          avg_q   = avg_q   + vol*rho0
          avg_vol = avg_vol + vol
        enddo
        array(0,:,k) = avg_q/avg_vol
      enddo

c     End program

      end subroutine scalarSingularBC

c     vectorSingularBC
c     #################################################################
      subroutine vectorSingularBC(vec,cov,order)
c     -----------------------------------------------------------------
c     Averages vector components around singular point and calculates
c     curvilinear components at singular point.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: order
      real(8)    :: vec(0:nx+1,0:ny+1,0:nz+1,3)
      logical    :: cov

c     Local variables

      integer(4) :: i,j,k

      integer(4) :: ic,ig,jg,kg,order1
      real(8)    :: x0,avg_q,avg_vol,vol
      real(8)    :: cx1,cy1,cz1
      real(8),allocatable,dimension(:) :: ax0,ay0,az0,cx,cy,cz

c     External

c     Begin program

      allocate(ax0(nz),ay0(nz),az0(nz))

c     Find average cartesian coordinates

      if (order == 3) then
        order1 = 3
      else
        order1 = order
      endif

      do k=1,nz
        ax0(k) = 0d0
        ay0(k) = 0d0
        az0(k) = 0d0
        avg_vol = 0d0
        do j=1,ny
          i = 1
          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
          vol = volume(i,j,k,igx,igy,igz)
          x0  = xx(ig-1)

          allocate(cx(order1+1),cy(order1+1),cz(order1+1))
          do i=1,order1+1
            call transformVectorToCartesian(i,j,k,igx,igy,igz
     .           ,vec(i,j,k,1),vec(i,j,k,2),vec(i,j,k,3),cov
     .           ,cx(i),cy(i),cz(i))
          enddo

          call IntDriver1d(order1+1,xx(ig),cx,1,x0,cx1,order)
          call IntDriver1d(order1+1,xx(ig),cy,1,x0,cy1,order)
          call IntDriver1d(order1+1,xx(ig),cz,1,x0,cz1,order)
          deallocate(cx,cy,cz)

          ax0(k)  = ax0(k)  + vol*cx1
          ay0(k)  = ay0(k)  + vol*cy1
          az0(k)  = az0(k)  + vol*cz1
          avg_vol = avg_vol + vol
        enddo
        ax0(k) = ax0(k)/avg_vol
        ay0(k) = ay0(k)/avg_vol
        az0(k) = az0(k)/avg_vol
      enddo
      
c     Transform to curvilinear components at SP

      i = 0
      do k=1,nz
        do j=1,ny
          call transformVectorToCurvilinear(i,j,k,igx,igy,igz
     .               ,ax0(k),ay0(k),az0(k),cov
     .               ,vec(i,j,k,1),vec(i,j,k,2),vec(i,j,k,3))
        enddo
      enddo

      deallocate(ax0,ay0,az0)

c     End program

      end subroutine vectorSingularBC

      end module singularBCinterface

c module dirichletBCinterface
c####################################################################
      module dirichletBCinterface

        use BCS

        INTERFACE dirichletBC
          module procedure scalarDirichletBC,vectorDirichletBC
        end INTERFACE

      contains

c     scalarDirichletBC
c     #################################################################
      subroutine scalarDirichletBC(array,array0,ieq,dim,loc,order)
c     -----------------------------------------------------------------
c     Imposes dirichlet BC. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c        * order -> order of extrapolation (when used)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Dirichlet BC.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ieq,dim,loc,order
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: icomp
      integer(4) :: i,j,k,ig,jg,kg,nvar,ibc
      real(8)    :: x1,x2,x3,dh(3),diver

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      rhs = 0d0

      call interpolate(array,array0,ibc,order)

c     End program

      end subroutine scalarDirichletBC

c     vectorDirichletBC
c     #################################################################
      subroutine vectorDirichletBC(ivar,array,array0,ieq,dim,loc,order)
c     -----------------------------------------------------------------
c     Imposes dirichlet BC. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c        * order -> order of extrapolation (when used)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Dirichlet BC.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ieq,dim,loc,order,ivar
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1,*)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1,*)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,nvar,ibc
      real(8)    :: x1,x2,x3,dh(3),diver

      real(8)    :: gsuper(3,3),jac0

c     Begin program

      rhs = 0d0

      nvar = nnvar

      ibc = (1+loc)+2*(dim-1)

      select case (ieq)
      case (IVX,IVY,IVZ)

        call interpolate(array(:,:,:,ivar),array0(:,:,:,ivar),ibc,order)

      case (IBX,IBY,IBZ,IJX,IJY,IJZ) !Imposes divergence-free constraint on B and J

        if (ivar /= dim) then

          call interpolate(array (:,:,:,ivar),array0(:,:,:,ivar)
     .                    ,ibc,order)

        else

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
                jac0 = gmetric%grid(igx)%jac(i,j,k)

                dh(1) = 2.*dxh(ig)
                dh(2) = 2.*dyh(jg)
                dh(3) = 2.*dzh(kg)

                select case (ibc)
                case (1)
                  array(i-1,j,k,dim) = array(i+1,j,k,dim)
                  diver = jac0*div(i,j,k,nx,ny,nz,array(:,:,:,1)
     .                            ,array(:,:,:,2),array(:,:,:,3))
                  rhs(j,k) = array(i+1,j,k,dim) + dh(dim)*diver
                case (2)
                  array(i+1,j,k,dim) = array(i-1,j,k,dim)
                  diver = jac0*div(i,j,k,nx,ny,nz,array(:,:,:,1)
     .                             ,array(:,:,:,2),array(:,:,:,3))
                  rhs(j,k) = array(i-1,j,k,dim) - dh(dim)*diver
                case (3)
                  array(i,j-1,k,dim) = array(i,j+1,k,dim)
                  diver = jac0*div(i,j,k,nx,ny,nz,array(:,:,:,1)
     .                             ,array(:,:,:,2),array(:,:,:,3))
                  rhs(i,k) = array(i,j+1,k,dim) + dh(dim)*diver
                case (4)
                  array(i,j+1,k,dim) = array(i,j-1,k,dim)
                  diver = jac0*div(i,j,k,nx,ny,nz,array(:,:,:,1)
     .                             ,array(:,:,:,2),array(:,:,:,3))
                  rhs(i,k) = array(i,j-1,k,dim) - dh(dim)*diver
                case (5)
                  array(i,j,k-1,dim) = array(i,j,k+1,dim)
                  diver = jac0*div(i,j,k,nx,ny,nz,array(:,:,:,1)
     .                             ,array(:,:,:,2),array(:,:,:,3))
                  rhs(i,j) = array(i,j,k+1,dim) + dh(dim)*diver
                case (6)
                  array(i,j,k+1,dim) = array(i,j,k-1,dim)
                  diver = jac0*div(i,j,k,nx,ny,nz,array(:,:,:,1)
     .                             ,array(:,:,:,2),array(:,:,:,3))
                  rhs(i,j) = array(i,j,k-1,dim) - dh(dim)*diver
                end select

              enddo
            enddo
          enddo
        endif

      case default

        write (*,*) 'Error in vectorDirichletBC'
        write (*,*) 'Equation',ieq,'does not exist'
        stop

      end select

c     End program

      end subroutine vectorDirichletBC

c     interpolate
c     #######################################################################
      subroutine interpolate(array,array0,ibc,order)
c     -----------------------------------------------------------------
c     Fills ghost nodes by extrapolation across relevant boundary.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: order,ibc
        real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
     .               ,array0(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer(4) :: i,j,k,ig,jg,kg

c     Begin program

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              select case (ibc)
              case (1)

                rhs(j,k) = 
     .             quad_int(xx(ig-1)+dxh(ig-1),xx(ig),xx(ig+1),xx(ig+2)
     .                     ,array0(i-1,j,k),array(i,j,k)
     .                     ,array (i+1,j,k),array(i+2,j,k)
     .                     ,xx(ig-1),order )
              case (2)

                rhs(j,k) =
     .              quad_int(xx(ig+1)-dxh(ig+1),xx(ig),xx(ig-1),xx(ig-2)
     .                      ,array0(i+1,j,k),array(i,j,k)
     .                      ,array (i-1,j,k),array(i-2,j,k)
     .                      ,xx(ig+1),order )
              case (3)

                rhs(i,k) =
     .              quad_int(yy(jg-1)+dyh(jg-1),yy(jg),yy(jg+1),yy(jg+2)
     .                      ,array0(i,j-1,k),array(i,j,k)
     .                      ,array (i,j+1,k),array(i,j+2,k)
     .                      ,yy(jg-1),order )
              case (4)

                rhs(i,k) =
     .              quad_int(yy(jg+1)-dyh(jg+1),yy(jg),yy(jg-1),yy(jg-2)
     .                      ,array0(i,j+1,k),array(i,j,k)
     .                      ,array (i,j-1,k),array(i,j-2,k)
     .                      ,yy(jg+1),order )
              case (5)

                rhs(i,j) =
     .              quad_int(zz(kg-1)+dzh(kg-1),zz(kg),zz(kg+1),zz(kg+2)
     .                      ,array0(i,j,k-1),array(i,j,k)
     .                      ,array (i,j,k+1),array(i,j,k+2)
     .                      ,zz(kg-1),order )
              case (6)

                rhs(i,j) =
     .              quad_int(zz(kg+1)-dzh(kg+1),zz(kg),zz(kg-1),zz(kg-2)
     .                      ,array0(i,j,k+1),array(i,j,k)
     .                      ,array (i,j,k-1),array(i,j,k-2)
     .                      ,zz(kg+1),order )
              end select
            enddo
          enddo
        enddo

      end subroutine interpolate

c     quad_int
c     #################################################################
      real(8) function quad_int(x0,x1,x2,x3,y0,y1,y2,y3,x,order)
     .        result(y)
c     -----------------------------------------------------------------
c     Interpolation (extrapolation) routine, up to cubic order.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: order
      real(8)    :: x0,x1,x2,x3,y0,y1,y2,y3,x

c     Local variables

c     Begin program

      select case (order)
      case (3)
        y = y0*(x-x1)*(x-x2)*(x-x3)/(x0-x1)/(x0-x2)/(x0-x3)
     .     +y1*(x-x0)*(x-x2)*(x-x3)/(x1-x0)/(x1-x2)/(x1-x3)
     .     +y2*(x-x0)*(x-x1)*(x-x3)/(x2-x0)/(x2-x1)/(x2-x3)
     .     +y3*(x-x0)*(x-x1)*(x-x2)/(x3-x0)/(x3-x1)/(x3-x2)
      case (2)
        y = y0*(x-x1)*(x-x2)/(x0-x1)/(x0-x2)
     .     +y1*(x-x0)*(x-x2)/(x1-x0)/(x1-x2)
     .     +y2*(x-x0)*(x-x1)/(x2-x0)/(x2-x1)
      case (1)
        y = y0*(x-x1)/(x0-x1)
     .     +y1*(x-x0)/(x1-x0)
      case (0)
        y = y0
      case default
        write (*,*) 'Order of interpolation not implemented in quad_int'
        write (*,*) 'Aborting...'
        stop
      end select

c     End program

      end function quad_int

      end module dirichletBCinterface

c module neumannBCinterface
c #####################################################################
      module neumannBCinterface

       use BCS

        INTERFACE neumannBC
          module procedure scalarNeumannBC,vectorNeumannBC
        end INTERFACE

      contains

c     scalarNeumannBC
c     #################################################################
      subroutine scalarNeumannBC(array,ieq,dim,loc)
c     -----------------------------------------------------------------
c     Imposes neumann BC for a scalar. On input:
c        * array -> variable upon which BCs are imposed
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Neumann BC.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ieq,dim,loc
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,ip,im,jp,jm,kp,km,nvar,ibc,icomp
      real(8)    :: x1,x2,x3,dh(3),jac0
      real(8)    :: gsuper(3,3),hessian(3,3,3)

c     Begin program

      nvar = 1

      ibc = (1+loc)+2*(dim-1)

      rhs = 0d0

      select case (ieq)
      case (IRHO)

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              gsuper = gmetric%grid(igx)%gsup(i,j,k,:,:)

              ip = min(i+1,nx)
              im = max(i-1,1)
              jp = min(j+1,ny)
              jm = max(j-1,1)
              kp = min(k+1,nz)
              km = max(k-1,1)

              dh(1) = 2.*dxh(ig)
              if (i == nx) dh(1) = dx(ig-1)
              if (i == 1 ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == ny) dh(2) = dy(jg-1)
              if (j == 1 ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nz) dh(3) = dz(kg-1)
              if (k == 1 ) dh(3) = dz(kg)

              if (dim == 1) then
                rhs(j,k) = -dh(dim)
     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3))
     .              /gsuper(dim,dim)
              elseif (dim == 2) then
                rhs(i,k) = -dh(dim)*
     .              (gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3))
     .              /gsuper(dim,dim)
              elseif (dim == 3) then
                rhs(i,j) = -dh(dim)
     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
     .              +gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2))
     .              /gsuper(dim,dim)
              endif

            enddo
          enddo
        enddo

      case (ITMP)

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              gsuper = gmetric%grid(igx)%gsup(i,j,k,:,:)
              jac0   = gmetric%grid(igx)%jac (i,j,k)

              hessian = gmetric%grid(igx)%Gamma(i,j,k,:,:,:)

              ip = min(i+1,nx)
              im = max(i-1,1)
              jp = min(j+1,ny)
              jm = max(j-1,1)
              kp = min(k+1,nz)
              km = max(k-1,1)

              dh(1) = 2.*dxh(ig)
              if (i == nx) dh(1) = dx(ig-1)
              if (i == 1 ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == ny) dh(2) = dy(jg-1)
              if (j == 1 ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nz) dh(3) = dz(kg-1)
              if (k == 1 ) dh(3) = dz(kg)

              if (dim == 1) then
                if (gamma > 1d0) then
                  rhs(j,k) =  hessian(1,1,1)*vx(i,j,k)*vx(i,j,k)
     .                       +hessian(1,2,2)*vy(i,j,k)*vy(i,j,k)
     .                       +hessian(1,3,3)*vz(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(1,1,2)*vx(i,j,k)*vy(i,j,k)
     .                    +2.*hessian(1,1,3)*vx(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(1,2,3)*vy(i,j,k)*vz(i,j,k)
                endif
                rhs(j,k) = -dh(dim)
     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
     .              -0.5/jac0*rhs(j,k))/gsuper(dim,dim)
              elseif (dim == 2) then
                if (gamma > 1d0) then
                  rhs(i,k) =  hessian(2,1,1)*vx(i,j,k)*vx(i,j,k)
     .                       +hessian(2,2,2)*vy(i,j,k)*vy(i,j,k)
     .                       +hessian(2,3,3)*vz(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(2,1,2)*vx(i,j,k)*vy(i,j,k)
     .                    +2.*hessian(2,1,3)*vx(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(2,2,3)*vy(i,j,k)*vz(i,j,k)
                endif
                rhs(i,k) = -dh(dim)
     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
     .              -0.5/jac0*rhs(i,k))/gsuper(dim,dim)
              elseif (dim == 3) then
                if (gamma > 1d0) then
                  rhs(i,j) =  hessian(3,1,1)*vx(i,j,k)*vx(i,j,k)
     .                       +hessian(3,2,2)*vy(i,j,k)*vy(i,j,k)
     .                       +hessian(3,3,3)*vz(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(3,1,2)*vx(i,j,k)*vy(i,j,k)
     .                    +2.*hessian(3,1,3)*vx(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(3,2,3)*vy(i,j,k)*vz(i,j,k)
                endif
                rhs(i,j) = -dh(dim)
     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
     .              +gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
     .              -0.5/jac0*rhs(i,j))/gsuper(dim,dim)
              endif

            enddo
          enddo
        enddo

      case default

        write (*,*) 'Error in scalarNeumannBC'
        stop

      end select

c     Assign value

      select case (ibc)
      case (1)
        rhs(jmin:jmax,kmin:kmax) = array(1,jmin:jmax,kmin:kmax)
     $                           - rhs(jmin:jmax,kmin:kmax)
      case (2)
        rhs(jmin:jmax,kmin:kmax) = array(nx,jmin:jmax,kmin:kmax)
     $                           + rhs(jmin:jmax,kmin:kmax)
      case (3)
        rhs(imin:imax,kmin:kmax) = array(imin:imax,1,kmin:kmax)
     $                           - rhs(imin:imax,kmin:kmax)
      case (4)
        rhs(imin:imax,kmin:kmax) = array(imin:imax,ny,kmin:kmax)
     $                           + rhs(imin:imax,kmin:kmax)
      case (5)
        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,1)
     $                           - rhs(imin:imax,jmin:jmax)
      case (6)
        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,nz)
     $                           + rhs(imin:imax,jmin:jmax)
      end select

c     End program

      end subroutine scalarNeumannBC

c     vectorNeumannBC
c     #################################################################
      subroutine vectorNeumannBC(ivar,array,ieq,dim,loc)
c     -----------------------------------------------------------------
c     Imposes neumann BC for a scalar. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Neumann BC.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ieq,dim,loc,ivar
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1,*)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,ip,im,jp,jm,kp,km,ibc,icomp
      real(8)    :: x1,x2,x3,dh(3),jac0,jxx,jyy,jzz
      real(8)    :: gsuper(3,3),hessian(3,3,3)

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      rhs = 0d0

      select case (ieq)
      case (IVX,IVY,IVZ) !Velocity components

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              gsuper = gmetric%grid(igx)%gsup(i,j,k,:,:)

              hessian = gmetric%grid(igx)%Gamma(i,j,k,:,:,:)
  
              ip = min(i+1,nx)
              im = max(i-1,1)
              jp = min(j+1,ny)
              jm = max(j-1,1)
              kp = min(k+1,nz)
              km = max(k-1,1)

              dh(1) = 2.*dxh(ig)
              if (i == nx) dh(1) = dx(ig-1)
              if (i == 1 ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == ny) dh(2) = dy(jg-1)
              if (j == 1 ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nz) dh(3) = dz(kg-1)
              if (k == 1 ) dh(3) = dz(kg)


              if (dim == 1) then

                if (ivar /= dim) then

                  rhs(j,k) =
     .                     gsuper(dim,1)
     .                      *(hessian(1,ivar,1)*array(i,j,k,1)
     .                       +hessian(2,ivar,1)*array(i,j,k,2)
     .                       +hessian(3,ivar,1)*array(i,j,k,3))
     .                    +gsuper(dim,2)
     .                      *(hessian(1,ivar,2)*array(i,j,k,1)
     .                       +hessian(2,ivar,2)*array(i,j,k,2)
     .                       +hessian(3,ivar,2)*array(i,j,k,3))
     .                    +gsuper(dim,3)
     .                      *(hessian(1,ivar,3)*array(i,j,k,1)
     .                       +hessian(2,ivar,3)*array(i,j,k,2)
     .                       +hessian(3,ivar,3)*array(i,j,k,3))

                  rhs(j,k) = -dh(dim)
     .                *(gsuper(dim,2)
     .                 *(array(i,jp,k,ivar)-array(i,jm,k,ivar))/dh(2)
     .                +gsuper(dim,3)
     .                 *(array(i,j,kp,ivar)-array(i,j,km,ivar))/dh(3)
     .                +rhs(j,k))/gsuper(dim,dim)

                endif

              elseif (dim == 2) then

                if (ivar /= dim) then

                  rhs(i,k) =
     .                     gsuper(dim,1)
     .                      *(hessian(1,ivar,1)*array(i,j,k,1)
     .                       +hessian(2,ivar,1)*array(i,j,k,2)
     .                       +hessian(3,ivar,1)*array(i,j,k,3))
     .                    +gsuper(dim,2)
     .                      *(hessian(1,ivar,2)*array(i,j,k,1)
     .                       +hessian(2,ivar,2)*array(i,j,k,2)
     .                       +hessian(3,ivar,2)*array(i,j,k,3))
     .                    +gsuper(dim,3)
     .                      *(hessian(1,ivar,3)*array(i,j,k,1)
     .                       +hessian(2,ivar,3)*array(i,j,k,2)
     .                       +hessian(3,ivar,3)*array(i,j,k,3))

                  rhs(i,k) = -dh(dim)
     .                 *(gsuper(dim,1)
     .                   *(array(ip,j,k,ivar)-array(im,j,k,ivar))/dh(1)
     .                 +gsuper(dim,3)
     .                   *(array(i,j,kp,ivar)-array(i,j,km,ivar))/dh(3)
     .                 +rhs(i,k))/gsuper(dim,dim)
                endif

              elseif (dim == 3) then

                if (ivar /= dim) then

                  rhs(i,j) =
     .                     gsuper(dim,1)
     .                      *(hessian(1,ivar,1)*array(i,j,k,1)
     .                       +hessian(2,ivar,1)*array(i,j,k,2)
     .                       +hessian(3,ivar,1)*array(i,j,k,3))
     .                    +gsuper(dim,2)
     .                      *(hessian(1,ivar,2)*array(i,j,k,1)
     .                       +hessian(2,ivar,2)*array(i,j,k,2)
     .                       +hessian(3,ivar,2)*array(i,j,k,3))
     .                    +gsuper(dim,3)
     .                      *(hessian(1,ivar,3)*array(i,j,k,1)
     .                       +hessian(2,ivar,3)*array(i,j,k,2)
     .                       +hessian(3,ivar,3)*array(i,j,k,3))

                  rhs(i,j) = -dh(dim)
     .               *(gsuper(dim,1)
     .                 *(array(ip,j,k,ivar)-array(im,j,k,ivar))/dh(1)
     .               +gsuper(dim,2)
     .                 *(array(i,jp,k,ivar)-array(i,jm,k,ivar))/dh(2)
     .               +rhs(i,j))/gsuper(dim,dim)

                endif

              endif

            enddo
          enddo
        enddo

      case (IBX,IBY,IBZ) 

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              gsuper = gmetric%grid(igx)%gsup(i,j,k,:,:)

              ip = min(i+1,nx)
              im = max(i-1,1)
              jp = min(j+1,ny)
              jm = max(j-1,1)
              kp = min(k+1,nz)
              km = max(k-1,1)

              dh(1) = 2.*dxh(ig)
              if (i == nx) dh(1) = dx(ig-1)
              if (i == 1 ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == ny) dh(2) = dy(jg-1)
              if (j == 1 ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nz) dh(3) = dz(kg-1)
              if (k == 1 ) dh(3) = dz(kg)

              if (dim == 1) then

                jxx = (array(i,jp,k,3)-array(i,jm,k,3))/dh(2)
     .               -(array(i,j,kp,2)-array(i,j,km,2))/dh(3)

                if (ivar == 2) then
                  rhs(j,k) = dh(dim)
     .               *( (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                 +gsuper(dim,3)*jxx/gsuper(dim,dim) )
                elseif (ivar == 3) then
                  rhs(j,k) = dh(dim)
     .               *( (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                 -gsuper(dim,2)*jxx/gsuper(dim,dim) )
                endif

              elseif (dim == 2) then

                jyy = (array(i,j,kp,1)-array(i,j,km,1))/dh(3)
     .               -(array(ip,j,k,3)-array(im,j,k,3))/dh(1)

                if (ivar == 3) then
                  rhs(i,k) = dh(dim)
     .               *( (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                 +gsuper(dim,1)*jyy/gsuper(dim,dim) )
                elseif (ivar == 1) then
                  rhs(i,k) = dh(dim)
     .               *( (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                 -gsuper(dim,3)*jyy/gsuper(dim,dim) )
                endif

              elseif (dim == 3) then

                jzz = (array(ip,j,k,2)-array(im,j,k,2))/dh(1)
     .               -(array(i,jp,k,1)-array(i,jm,k,1))/dh(2)

                if (ivar == 1) then
                  rhs(i,j) = dh(dim)
     .               *( (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                 +gsuper(dim,2)*jzz/gsuper(dim,dim) )
                elseif (ivar == 2) then
                  rhs(i,j) = dh(dim)
     .               *( (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                 -gsuper(dim,1)*jzz/gsuper(dim,dim) )
                endif

              endif

            enddo
          enddo
        enddo

      case default

        write (*,*) 'Error in vectorNeumannBC'
        stop

      end select

c     Assign value

      select case (ibc)
      case (1)
        rhs(jmin:jmax,kmin:kmax) = array(1,jmin:jmax,kmin:kmax,ivar)
     $                           - rhs(jmin:jmax,kmin:kmax)
      case (2)
        rhs(jmin:jmax,kmin:kmax) = array(nx,jmin:jmax,kmin:kmax,ivar)
     $                           + rhs(jmin:jmax,kmin:kmax)
      case (3)
        rhs(imin:imax,kmin:kmax) = array(imin:imax,1,kmin:kmax,ivar)
     $                           - rhs(imin:imax,kmin:kmax)
      case (4)
        rhs(imin:imax,kmin:kmax) = array(imin:imax,ny,kmin:kmax,ivar)
     $                           + rhs(imin:imax,kmin:kmax)
      case (5)
        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,1,ivar)
     $                           - rhs(imin:imax,jmin:jmax)
      case (6)
        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,nz,ivar)
     $                           + rhs(imin:imax,jmin:jmax)
      end select

c     End program

      end subroutine vectorNeumannBC

      end module neumannBCinterface

c module imposeBCinterface
c #####################################################################
      module imposeBCinterface

       use BCS

       use singularBCinterface

       use dirichletBCinterface

       use neumannBCinterface

        INTERFACE setBC
          module procedure imposeBConScalar,imposeBConVector
        end INTERFACE

        integer(4),private :: bctype,dim,loc,ibc,i,j,k

      contains

c     imposeBConScalar
c     #################################################################
      subroutine imposeBConScalar(ieq,array,array0,bcond)
c     -----------------------------------------------------------------
c     Imposes BC on density
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ieq,bcond(6)
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

c     Begin program

c     Impose BCs

      do bctype=1,4             !Enforces a particular order in the BCs (see grid_mod.f)
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
            if (bcond(ibc) == bctype) then
              call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
            endif
          enddo
        enddo
      enddo

c     Singular point boundary condition

      if (bcond(1) == SP) call singularBC(array,2)

c     Synchronize periodic boundaries

cc      bctype=PER
cc
cc      do dim=1,3
cc        do loc=0,1
cc          ibc = (1+loc)+2*(dim-1)
cc          if (bcond(ibc) == bctype) then
cc            call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
cc          endif
cc        enddo
cc      enddo

c     Fill edges

      do dim=1,3
        do loc=0,1
          call fillEdges(array,dim,loc)
        enddo
      enddo

c     Fill corners

      do dim=1,3
        do loc=0,1
          call fillCorners(array,dim,loc)
        enddo
      enddo

c     End program

      end subroutine imposeBConScalar

c     imposeBConVector
c     #################################################################
      subroutine imposeBConVector(fcomp,v_cnv,v_cov,var0,bcond)
c     -----------------------------------------------------------------
c     Imposes BC on velocity field
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: bcond(6,3),fcomp
      real(8)    :: v_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,var0 (0:nx+1,0:ny+1,0:nz+1,3)
     .             ,v_cov(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

      integer(4) :: ivar,ieq,ibc,loc,dim,bctype
      logical    :: cov_to_cnv

c     Begin program

c     Find cov components (we assume cnv components are defined)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .            ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .            ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),.false.)
          enddo
        enddo
      enddo

c     Impose BCs

      cov_to_cnv = .false.

      do bctype=1,4            !Enforces a particular order in the BCs (see grid_mod.f)
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            do ivar = 1,3
              ieq = ivar + fcomp - 1

              if (abs(bcond(ibc,ivar)) == bctype) then

                if (bcond(ibc,ivar) < 0) then
                  call FillGhostNodes(ieq,ivar,3,dim,loc,bctype,v_cov
     .                               ,var0)
                  cov_to_cnv = .true.
                else
                  call FillGhostNodes(ieq,ivar,3,dim,loc,bctype,v_cnv
     .                               ,var0)
                endif

             endif

            enddo

          enddo
        enddo
      enddo

c     Synchronize covariant and contravariant components

      if (cov_to_cnv) call synchronize_cnv(v_cnv,v_cov,bcond)

c     Impose vector singular point BCs

      if (bcond(1,1) == SP) call singularBC(v_cnv,.false.,2)

c     Synchronize periodic boundaries

cc      bctype=PER
cc
cc      do dim=1,3
cc        do loc=0,1
cc          ibc = (1+loc)+2*(dim-1)
cc          do ivar = 1,3
cc            ieq = ivar + fcomp - 1
cc            if (bcond(ibc,ivar) == bctype) then
cc              call FillGhostNodes(ieq,ivar,3,dim,loc,bctype,v_cnv,var0)
cc            endif
cc          enddo
cc        enddo
cc      enddo

c     Fill edges

      do dim=1,3
        do loc=0,1
          do ivar=1,3
            call fillEdges(v_cnv(:,:,:,ivar),dim,loc)
          enddo
        enddo
      enddo

c     Fill corners

      do dim=1,3
        do loc=0,1
          do ivar=1,3
            call fillCorners(v_cnv(:,:,:,ivar),dim,loc)
          enddo
        enddo
      enddo

c     Find covariant components at ALL boundaries

      do dim = 1,3
        do loc = 0,1

          call findLoopLimits(dim,loc,iimin-1,iimax+1
     .                               ,jjmin-1,jjmax+1
     .                               ,kkmin-1,kkmax+1)

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax
                call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .            ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .            ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),.false.)
              enddo
            enddo
          enddo

        enddo
      enddo

c     End program

      end subroutine imposeBConVector

c     synchronize_cnv
c     #################################################################
      subroutine synchronize_cnv(v_cnv,v_cov,bcond)
c     -----------------------------------------------------------------
c     Finds all contravariant components at Neumann boundaries.
c     On input, tangential covariant components and normal contravariant
c     components are known at ghost cells. On output, all contravariant
c     components are known at ghost cells.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: bcond(6,3)
      real(8)    :: v_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,v_cov(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

      integer(4) :: i,j,k,dim,loc,ig,jg,kg,ivar
      real(8)    :: x1,x2,x3,gsuper(3,3),gsub(3,3)

c     Begin program

      do dim = 1,3
        do loc = 0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,3

            if (ivar == dim) then  !Select tangential components
              cycle
            elseif (bcond(ibc,ivar) < 0) then

              call findLoopLimits(dim,loc,iimin,iimax
     .                                   ,jjmin,jjmax
     .                                   ,kkmin,kkmax)

              select case (ibc)
              case (1)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igx)%gsup(i-1,j,k,:,:)

                      v_cov(i-1,j,k,1) = -(gsuper(1,2)*v_cov(i-1,j,k,2)
     .                                    +gsuper(1,3)*v_cov(i-1,j,k,3)
     .                                    -v_cnv(i-1,j,k,1))/gsuper(1,1)

                      call transformFromCurvToCurv(i-1,j,k,igx,igy,igz
     .               ,v_cov(i-1,j,k,1),v_cov(i-1,j,k,2),v_cov(i-1,j,k,3)
     .               ,v_cnv(i-1,j,k,1),v_cnv(i-1,j,k,2),v_cnv(i-1,j,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (2)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igx)%gsup(i+1,j,k,:,:)

                      v_cov(i+1,j,k,1) = -(gsuper(1,2)*v_cov(i+1,j,k,2)
     .                                    +gsuper(1,3)*v_cov(i+1,j,k,3)
     .                                    -v_cnv(i+1,j,k,1))/gsuper(1,1)

                      call transformFromCurvToCurv(i+1,j,k,igx,igy,igz
     .               ,v_cov(i+1,j,k,1),v_cov(i+1,j,k,2),v_cov(i+1,j,k,3)
     .               ,v_cnv(i+1,j,k,1),v_cnv(i+1,j,k,2),v_cnv(i+1,j,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (3)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igx)%gsup(i,j-1,k,:,:)

                      v_cov(i,j-1,k,2) = -(gsuper(2,1)*v_cov(i,j-1,k,1)
     .                                    +gsuper(2,3)*v_cov(i,j-1,k,3)
     .                                    -v_cnv(i,j-1,k,2))/gsuper(2,2)

                      call transformFromCurvToCurv(i,j-1,k,igx,igy,igz
     .               ,v_cov(i,j-1,k,1),v_cov(i,j-1,k,2),v_cov(i,j-1,k,3)
     .               ,v_cnv(i,j-1,k,1),v_cnv(i,j-1,k,2),v_cnv(i,j-1,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (4)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igx)%gsup(i,j+1,k,:,:)

                      v_cov(i,j+1,k,2) = -(gsuper(2,1)*v_cov(i,j+1,k,1)
     .                                    +gsuper(2,3)*v_cov(i,j+1,k,3)
     .                                    -v_cnv(i,j+1,k,2))/gsuper(2,2)

                      call transformFromCurvToCurv(i,j+1,k,igx,igy,igz
     .               ,v_cov(i,j+1,k,1),v_cov(i,j+1,k,2),v_cov(i,j+1,k,3)
     .               ,v_cnv(i,j+1,k,1),v_cnv(i,j+1,k,2),v_cnv(i,j+1,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (5)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igx)%gsup(i,j,k-1,:,:)

                      v_cov(i,j,k-1,3) = -(gsuper(3,1)*v_cov(i,j,k-1,1)
     .                                    +gsuper(3,2)*v_cov(i,j,k-1,2)
     .                                    -v_cnv(i,j,k-1,3))/gsuper(3,3)

                      call transformFromCurvToCurv(i,j,k-1,igx,igy,igz
     .               ,v_cov(i,j,k-1,1),v_cov(i,j,k-1,2),v_cov(i,j,k-1,3)
     .               ,v_cnv(i,j,k-1,1),v_cnv(i,j,k-1,2),v_cnv(i,j,k-1,3)
     .               ,.true.)
                    enddo
                  enddo
                enddo

              case (6)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igx)%gsup(i,j,k+1,:,:)

                      v_cov(i,j,k+1,3) = -(gsuper(3,1)*v_cov(i,j,k+1,1)
     .                                    +gsuper(3,2)*v_cov(i,j,k+1,2)
     .                                    -v_cnv(i,j,k+1,3))/gsuper(3,3)

                      call transformFromCurvToCurv(i,j,k+1,igx,igy,igz
     .               ,v_cov(i,j,k+1,1),v_cov(i,j,k+1,2),v_cov(i,j,k+1,3)
     .               ,v_cnv(i,j,k+1,1),v_cnv(i,j,k+1,2),v_cnv(i,j,k+1,3)
     .               ,.true.)
                    enddo
                  enddo
                enddo

              end select

            endif
          enddo

        enddo
      enddo

c     End program

      end subroutine synchronize_cnv

c     fillGhostNodes
c     ###############################################################
      subroutine fillGhostNodes(ieq,ivar,nvar,dim,loc,bctype
     .                         ,array,array0)

c     ---------------------------------------------------------------
c     Sets adequate boundary conditions on array.
c
c     On input:
c       * ieq    -> equation identifier
c       * ivar   -> vector component
c       * nvar   -> vector dimension
c       * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c       * loc    -> location in dimension (0 -> right, 1 -> left)
c       * bctype -> type of BC (dirichlet, neumann, periodic, etc.)
c       * array  -> real array with ghost-nodes
c       * array0 -> auxiliary real array
c     ---------------------------------------------------------------

      implicit none       !For safe fortran

c     Call variables

      integer(4) :: ieq,dim,loc,bctype,nvar,ivar
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1,nvar)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1,nvar)

c     Local variables

      integer(4) :: neq,ibc,i,j,k,ig,jg,kg

c     Begin program

      nnvar = nvar

c     Determine boundary limits

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax)

c     Allocate rhs

      select case(dim)
      case(1)
        allocate(rhs(0:ny+1,0:nz+1))
cc        allocate(rhs(jmin:jmax,kmin:kmax))
      case(2)
        allocate(rhs(0:nx+1,0:nz+1))
cc        allocate(rhs(imin:imax,kmin:kmax))
      case(3)
        allocate(rhs(0:nx+1,0:ny+1))
cc        allocate(rhs(imin:imax,jmin:jmax))
      end select

c     Find BC update

      ibc = (1+loc)+2*(dim-1)

      if (nvar == 1) then
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,nvar),ibc)
        case(EQU)
          call dirichletBC(array(:,:,:,nvar),array0(:,:,:,nvar)
     .                    ,ieq,dim,loc,0)
        case(DIR)
          call dirichletBC(array(:,:,:,nvar),zeros,ieq,dim,loc,1)
        case(NEU)
          call neumannBC(array(:,:,:,nvar),ieq,dim,loc)
        case default
          write (*,*) 'BC',bctype,' not implemented'
          stop
        end select
      else
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,ivar),ibc)
        case(EQU)
          call dirichletBC(ivar,array,array0,ieq,dim,loc,0)
        case(DIR)
          call dirichletBC(ivar,array,vzeros,ieq,dim,loc,1)
        case(NEU)
          call neumannBC(ivar,array,ieq,dim,loc)
        case default
          write (*,*) 'BC',bctype,' not implemented'
          stop
        end select
      endif

c     Update BC ghost nodes

      select case (ibc)
      case (1)                  !x0
        array(0   ,jmin:jmax,kmin:kmax,ivar) = rhs(jmin:jmax,kmin:kmax)
      case (2)                  !x1
        array(nx+1,jmin:jmax,kmin:kmax,ivar) = rhs(jmin:jmax,kmin:kmax)
      case (3)                  !y0
        array(imin:imax,0   ,kmin:kmax,ivar) = rhs(imin:imax,kmin:kmax)
      case (4)                  !y1
        array(imin:imax,ny+1,kmin:kmax,ivar) = rhs(imin:imax,kmin:kmax)
      case (5)                  !z0
        array(imin:imax,jmin:jmax,0   ,ivar) = rhs(imin:imax,jmin:jmax)
      case (6)                  !z1
        array(imin:imax,jmin:jmax,nz+1,ivar) = rhs(imin:imax,jmin:jmax)
      case default
        write (*,*) 'Boundary',ibc,' non existent'
        stop
      end select

cc      select case (ibc)
cc      case (1)                  !x0
cc        array(0   ,:,:,ivar) = rhs(:,:)
cc      case (2)                  !x1
cc        array(nx+1,:,:,ivar) = rhs(:,:)
cc      case (3)                  !y0
cc        array(:,0   ,:,ivar) = rhs(:,:)
cc      case (4)                  !y1
cc        array(:,ny+1,:,ivar) = rhs(:,:)
cc      case (5)                  !z0
cc        array(:,:,0   ,ivar) = rhs(:,:)
cc      case (6)                  !z1
cc        array(:,:,nz+1,ivar) = rhs(:,:)
cc      case default
cc        write (*,*) 'Boundary',ibc,' non existent'
cc        stop
cc      end select

      deallocate(rhs)

c     End

      end subroutine fillGhostNodes

c     periodicBC
c     #################################################################
      subroutine periodicBC(array,ibc)
c     -----------------------------------------------------------------
c     Imposes periodic BC. On input:
c        * array -> 3D array to impose BC's on.
c        * ibc   -> integer boundary identifier
c                   (0 -> X1, 1->X2, 2->y1, 3->y2, 4->z1, 5->z2)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ibc
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)

c     Local variables

c     Begin program

cc#if !defined(petsc)

cc      select case (ibc)
cc      case (1)
cc        rhs(:,:) = array(nx,:,:)
cc      case (2)
cc        rhs(:,:) = array(1,:,:)
cc      case (3)
cc        rhs(:,:) = array(:,ny,:)
cc      case (4)
cc        rhs(:,:) = array(:,1,:)
cc      case (5)
cc        rhs(:,:) = array(:,:,nz)
cc      case (6)
cc        rhs(:,:) = array(:,:,1)
cc      end select

      select case (ibc)
      case (1)
        rhs(jmin:jmax,kmin:kmax) = array(nx,jmin:jmax,kmin:kmax)
      case (2)
        rhs(jmin:jmax,kmin:kmax) = array(1,jmin:jmax,kmin:kmax)
      case (3)
        rhs(imin:imax,kmin:kmax) = array(imin:imax,ny,kmin:kmax)
      case (4)
        rhs(imin:imax,kmin:kmax) = array(imin:imax,1,kmin:kmax)
      case (5)
        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,nz)
      case (6)
        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,1)
      end select

cc#else !Do nothing

cc      select case (ibc)
cc      case (1)
cc        rhs(jmin:jmax,kmin:kmax) = array(0   ,jmin:jmax,kmin:kmax)
cc      case (2)
cc        rhs(jmin:jmax,kmin:kmax) = array(nx+1,jmin:jmax,kmin:kmax)
cc      case (3)
cc        rhs(imin:imax,kmin:kmax) = array(imin:imax,0   ,kmin:kmax)
cc      case (4)
cc        rhs(imin:imax,kmin:kmax) = array(imin:imax,ny+1,kmin:kmax)
cc      case (5)
cc        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,0   )
cc      case (6)
cc        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,nz+1)
cc      end select

cccc      select case (ibc)
cccc      case (1)
cccc        rhs(:,:) = array(0,:,:)
cccc      case (2)
cccc        rhs(:,:) = array(nx+1,:,:)
cccc      case (3)
cccc        rhs(:,:) = array(:,0,:)
cccc      case (4)
cccc        rhs(:,:) = array(:,ny+1,:)
cccc      case (5)
cccc        rhs(:,:) = array(:,:,0)
cccc      case (6)
cccc        rhs(:,:) = array(:,:,nz+1)
cccc      end select
cc
cc#endif

c     End program

      end subroutine periodicBC

c     fillEdges
c     #################################################################
      subroutine fillEdges(array,dim,loc)
c     -----------------------------------------------------------------
c     Fills ghost cells in edges (12 in 3D). On input:
c        * array -> 3D array to impose BC's on.
c        * ibc   -> integer boundary identifier
c                   (1->X1, 2->X2, 3->y1, 4->y2, 5->z1, 6->z2)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: dim,loc
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,nnodes,ibc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

c     Determine boundary limits

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax)

c     Select edge

      nnodes = 3

      select case (ibc)
      case (1) !Fill edges y=jmin-1,y=jmax+1, vary z

        if (imin == 1) then
          ito   = imin-1
          ifrom = min(imin+nnodes,nx)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          if (jmin == 1) then
            jto   = jmin-1
            jfrom = min(jmin+nnodes,ny)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              call extrapolateToCorner(ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k))
            enddo
          endif

          if (jmax == ny) then
            jto   = jmax+1
            jfrom = max(jmax-nnodes,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              call extrapolateToCorner(ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k))
            enddo
          endif

        endif

      case (2) !Fill edges y=jmin-1,y=jmax+1, vary z

        if (imax == nx) then
          ito   = imax+1
          ifrom = max(imax-nnodes,1)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          if (jmin == 1) then
            jto   = jmin-1
            jfrom = min(jmin+nnodes,ny)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              call extrapolateToCorner(ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k))
            enddo
          endif

          if (jmax == ny) then
            jto   = jmax+1
            jfrom = max(jmax-nnodes,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              call extrapolateToCorner(ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k))
            enddo
          endif

        endif

      case (3) !Fill edges z=kmin-1,z=kmax+1, vary x

        if (jmin == 1) then

          jto   = jmin-1
          jfrom = min(jmin+nnodes,ny)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          if (kmin == 1) then

            kto   = kmin-1
            kfrom = min(kmin+nnodes,nz)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              call extrapolateToCorner(i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx))
            enddo
          endif

          if (kmax == nz) then

            kto   = kmax+1
            kfrom = max(kmax-nnodes,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              call extrapolateToCorner(i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx))
            enddo
          endif

        endif

      case (4) !Fill edges z=kmin-1,z=kmax+1, vary x

        if (jmax == ny) then

          jto   = jmax+1
          jfrom = max(jmax-nnodes,1)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          if (kmin == 1) then

            kto   = kmin-1
            kfrom = min(kmin+nnodes,nz)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              call extrapolateToCorner(i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx))
            enddo
          endif

          if (kmax == nz) then

            kto   = kmax+1
            kfrom = max(kmax-nnodes,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              call extrapolateToCorner(i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx))
            enddo
          endif

        endif

      case (5) !Fill edges x=imin-1,x=imax+1, vary y

        if (kmin == 1) then

          kto   = kmin-1
          kfrom = min(kmin+nnodes,nz)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          if (imin == 1) then

            ito   = imin-1
            ifrom = min(imin+nnodes,nx)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              call extrapolateToCorner(ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx))
            enddo
          endif

          if (imax == nx) then

            ito   = imax+1
            ifrom = max(imax-nnodes,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              call extrapolateToCorner(ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx))
            enddo
          endif

        endif

      case (6) !Fill edges x=imin-1,x=imax+1, vary y

        if (kmax == nz) then

          kto   = kmax+1
          kfrom = max(kmax-nnodes,1)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          if (imin == 1) then

            ito   = imin-1
            ifrom = min(imin+nnodes,nx)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              call extrapolateToCorner(ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx))
            enddo
          endif

          if (imax == nx) then

            ito   = imax+1
            ifrom = max(imax-nnodes,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              call extrapolateToCorner(ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx))
            enddo
          endif

        endif

      end select

c     End program

      end subroutine fillEdges

c     fillCorners
c     #################################################################
      subroutine fillCorners(array,dim,loc)
c     -----------------------------------------------------------------
c     Imposes ghost cells in corners (8 in 3D). On input:
c        * array -> 3D array to impose BC's on.
c        * ibc   -> integer boundary identifier
c                   (0 -> X1, 1->X2, 2->y1, 3->y2, 4->z1, 5->z2)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: dim,loc
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,nnodes,ibc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

c     Determine boundary limits

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax)

c     Select edge

      nnodes = 3

      if     (imin == 1  .and. jmin == 1  .and. kmin == 1) then
        ifrom = min(imin+nnodes,nx)
        jfrom = min(jmin+nnodes,ny)
        kfrom = min(kmin+nnodes,nz)
        ito   = imin-1
        jto   = jmin-1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ito,jto,kto,imn,jmn,kmn
     .                  ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx))
      endif

      if (imin == 1  .and. jmax == ny .and. kmin == 1) then

        ifrom = min(imin+nnodes,nx)
        jfrom = max(jmax-nnodes,1 )
        kfrom = min(kmin+nnodes,nz)
        ito   = imin-1
        jto   = jmax+1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ito,jto,kto,imn,jmn,kmn
     .                  ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx))
      endif

      if (imax == nx .and. jmin == 1  .and. kmin == 1) then

        ifrom = max(imax-nnodes,1 )
        jfrom = min(jmin+nnodes,ny)
        kfrom = min(kmin+nnodes,nz)
        ito   = imax+1
        jto   = jmin-1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ito,jto,kto,imn,jmn,kmn
     .                  ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx))
      endif

      if (imax == nx .and. jmax == ny .and. kmin == 1) then

        ifrom = max(imax-nnodes,1 )
        jfrom = max(jmax-nnodes,1 )
        kfrom = min(kmin+nnodes,nz)
        ito   = imax+1
        jto   = jmax+1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ito,jto,kto,imn,jmn,kmn
     .                  ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx))
      endif

      if (imin == 1  .and. jmin == 1  .and. kmax == nz) then

        ifrom = min(imin+nnodes,nx)
        jfrom = min(jmin+nnodes,ny)
        kfrom = max(kmax-nnodes,1 )
        ito   = imin-1
        jto   = jmin-1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ito,jto,kto,imn,jmn,kmn
     .                  ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx))
      endif

      if (imin == 1  .and. jmax == ny .and. kmax == nz) then

        ifrom = min(imin+nnodes,nx)
        jfrom = max(jmax-nnodes,1 )
        kfrom = max(kmax-nnodes,1 )
        ito   = imin-1
        jto   = jmax+1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ito,jto,kto,imn,jmn,kmn
     .                  ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx))
      endif

      if (imax == nx .and. jmin == 1  .and. kmax == nz) then

        ifrom = max(imax-nnodes,1 )
        jfrom = min(jmin+nnodes,ny)
        kfrom = max(kmax-nnodes,1 )
        ito   = imax+1
        jto   = jmin-1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ito,jto,kto,imn,jmn,kmn
     .                  ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx))
      endif

      if (imax == nx .and. jmax == ny .and. kmax == nz) then

        ifrom = max(imax-nnodes,1 )
        jfrom = max(jmax-nnodes,1 )
        kfrom = max(kmax-nnodes,1 )
        ito   = imax+1
        jto   = jmax+1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ito,jto,kto,imn,jmn,kmn
     .                  ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx))

      endif

c     End program

      end subroutine fillCorners

c     extrapolateToCorner
c     #################################################################
      subroutine extrapolateToCorner(ito,jto,kto,imn,jmn,kmn,imx,jmx,kmx
     .                              ,array)
c     -----------------------------------------------------------------

c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ito,jto,kto,imn,jmn,kmn,imx,jmx,kmx
      real(8)    :: array(imn:imx,jmn:jmx,kmn:kmx)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,order,inx,iny,inz
     .             ,orderx,ordery,orderz
     .             ,imni,imxi,jmni,jmxi,kmni,kmxi,num

      real(8)    :: valx,valy,valz,xto,yto,zto

      logical    :: flgx,flgy,flgz

c     Extrapolation

      real(8),allocatable,dimension(:) :: xx,yy,zz

c     Begin program

      inx = imx-imn
      iny = jmx-jmn
      inz = kmx-kmn

      flgx = .true.
      flgy = .true.
      flgz = .true.
      if (inx == 0) flgx = .false.
      if (iny == 0) flgy = .false.
      if (inz == 0) flgz = .false.

      allocate(xx(inx),yy(iny),zz(inz))

      valx = 0d0
      valy = 0d0
      valz = 0d0

      num = 0

c     X extrapolation

      if (flgx) then
        if (ito == imn) then
          call getMGmap(imn+1,jmn,kmn,igx,igy,igz,ig,jg,kg)
          xx = grid_params%xx(ig:ig+inx-1)
          imni = imn+1
          imxi = imx
        else
          call getMGmap(imx-1,jmn,kmn,igx,igy,igz,ig,jg,kg)
          xx = grid_params%xx(ig-(inx-1):ig)
          imni = imn
          imxi = imx-1
        endif

        call getMGmap(ito,jto,kto,igx,igy,igz,ig,jg,kg)
        xto = grid_params%xx(ig)

        orderx = min(2,inx-1)
        call IntDriver1d(inx,xx,array(imni:imxi,jto,kto)
     .                  ,1,xto,valx,orderx)

        num = num + 1
      endif

c     Y extrapolation

      if (flgy) then

        if (jto == jmn) then
          call getMGmap(imn,jmn+1,kmn,igx,igy,igz,ig,jg,kg)
          yy = grid_params%yy(jg:jg+iny-1)
          jmni = jmn+1
          jmxi = jmx
        else
          call getMGmap(imn,jmx-1,kmn,igx,igy,igz,ig,jg,kg)
          yy = grid_params%yy(jg-(iny-1):jg)
          jmni = jmn
          jmxi = jmx-1
        endif

        call getMGmap(ito,jto,kto,igx,igy,igz,ig,jg,kg)
        yto = grid_params%yy(jg)

        ordery = min(2,iny-1)
        call IntDriver1d(iny,yy,array(ito,jmni:jmxi,kto)
     .                  ,1,yto,valy,ordery)
        
        num = num + 1
      endif

c     Z extrapolation

      if (flgz) then

        if (kto == kmn) then
          call getMGmap(imn,jmn,kmn+1,igx,igy,igz,ig,jg,kg)
          zz = grid_params%zz(kg:kg+inz-1)
          kmni = kmn+1
          kmxi = kmx
        else
          call getMGmap(imn,jmn,kmx-1,igx,igy,igz,ig,jg,kg)
          zz = grid_params%zz(kg-(inz-1):kg)
          kmni = kmn
          kmxi = kmx-1
        endif

        call getMGmap(ito,jto,kto,igx,igy,igz,ig,jg,kg)
        zto = grid_params%zz(kg)

        orderz = min(2,inz-1)
        call IntDriver1d(inz,zz,array(ito,jto,kmni:kmxi)
     .                  ,1,zto,valz,orderz)

        num = num + 1
      endif

c     Average results

cc      write (*,*) ito,jto,kto
cc      write (*,*) ito,jto,kto,valx,valy,valz,num
cc      do k=kmn,kmx
cc        do j=jmn,jmx
cc          do i=imn,imx
cc            write (*,*) i,j,k,array(i,j,k)
cc          enddo
cc        enddo
cc      enddo
      array(ito,jto,kto) = (valx + valy + valz)/num

c     End program

      deallocate(xx,yy,zz)

      end subroutine extrapolateToCorner

      end module imposeBCinterface
