c defineBoundaryConditions
c####################################################################
      subroutine defineBoundaryConditions (neq,bbcs)
c--------------------------------------------------------------------
c     Defines boundary conditions of physical quantities.
c     On input:
c       * neq -> number of equations
c     On output:
c       * bbcs -> real array of size (6,neq) containing BC setup:
c           + bbcs(1) ---> at x0
c           + bbcs(2) ---> at x1
c           + bbcs(3) ---> at y0
c           + bbcs(4) ---> at y1
c           + bbcs(5) ---> at z0
c           + bbcs(6) ---> at z1
c     Definition of BC identifiers is given in "grid_mod.f". In vectors,
c     a negative BC identifier means that BCs are to be imposed on
c     covariant components instead of on (default) contravariant comps.
c--------------------------------------------------------------------

      use icond

      use grid

      use equilibrium

      implicit none

c Call variables

      integer(4) :: neq,bbcs(6,neq)

c Local variables

      integer(4) :: ieq,bcsq(6)

c Begin program

c Default boundary conditions

      bcsq = bbcs(:,IRHO)
      where (bcsq == DEF) bcsq = NEU
      bbcs(:,IRHO) = bcsq

      bcsq = bbcs(:,IVX)
      where (bcsq == DEF) bcsq = DIR
      bbcs(:,IVX) = bcsq

      bcsq = bbcs(:,IVY)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IVY) = bcsq

      bcsq = bbcs(:,IVZ)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IVZ) = bcsq

      bcsq = bbcs(:,IBX)
      where (bcsq == DEF) bcsq = DIR
      bbcs(:,IBX) = bcsq

      bcsq = bbcs(:,IBY)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IBY) = bcsq

      bcsq = bbcs(:,IBZ)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IBZ) = bcsq

      bcsq = bbcs(:,ITMP)
      where (bcsq == DEF) bcsq = NEU !To allow isothermal case
      bbcs(:,ITMP) = bcsq

c Exceptions for specific equilibria

      select case (equil)

      case ('rfp1')

        bbcs(2,IBY) = EQU  !Imposed by equilibrium

        bbcs(2,IBZ) = EQU  !Imposed by equilibrium

      end select

c End

      end subroutine defineBoundaryConditions

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

        real(8),allocatable,dimension(:,:) :: rhs

        real(8),allocatable,dimension(:,:,:,:) :: v_cov,v_cnv,vzeros

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

cc      use grid

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
      logical    :: cartesian

      integer(4) :: ic,ig,jg,kg,order1
      real(8)    :: x0,avg_q,avg_vol,vol
      real(8)    :: cx1,cy1,cz1
      real(8),allocatable,dimension(:) :: ax0,ay0,az0,cx,cy,cz

c     External

      real(8) :: quad_int
      external   quad_int

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

c     Begin program

      rhs = 0d0

      nvar = nnvar

      ibc = (1+loc)+2*(dim-1)

      select case (ieq)
      case (IVX,IVY,IVZ)

        call interpolate(array(:,:,:,ivar),array0(:,:,:,ivar),ibc,order)

      case (IBX,IBY,IBZ,IJX,IJY,IJZ) !Imposes divergence-free constraint on B-field

        if (ivar /= dim) then

          call interpolate(array (:,:,:,ivar),array0(:,:,:,ivar)
     .                    ,ibc,order)

        else

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

                dh(1) = 2.*dxh(ig)
                dh(2) = 2.*dyh(jg)
                dh(3) = 2.*dzh(kg)

                select case (ibc)
                case (1)
                  array(i-1,j,k,dim) = array(i+1,j,k,dim)
                  diver = div(i,j,k,array(:,:,:,1)
     .                             ,array(:,:,:,2),array(:,:,:,3))
                  rhs(j,k) = array(i+1,j,k,dim) + dh(dim)*diver
                case (2)
                  array(i+1,j,k,dim) = array(i-1,j,k,dim)
                  diver = div(i,j,k,array(:,:,:,1)
     .                             ,array(:,:,:,2),array(:,:,:,3))
                  rhs(j,k) = array(i-1,j,k,dim) - dh(dim)*diver
                case (3)
                  array(i,j-1,k,dim) = array(i,j+1,k,dim)
                  diver = div(i,j,k,array(:,:,:,1)
     .                             ,array(:,:,:,2),array(:,:,:,3))
                  rhs(i,k) = array(i,j+1,k,dim) + dh(dim)*diver
                case (4)
                  array(i,j+1,k,dim) = array(i,j-1,k,dim)
                  diver = div(i,j,k,array(:,:,:,1)
     .                             ,array(:,:,:,2),array(:,:,:,3))
                  rhs(i,k) = array(i,j-1,k,dim) - dh(dim)*diver
                case (5)
                  array(i,j,k-1,dim) = array(i,j,k+1,dim)
                  diver = div(i,j,k,array(:,:,:,1)
     .                             ,array(:,:,:,2),array(:,:,:,3))
                  rhs(i,j) = array(i,j,k+1,dim) + dh(dim)*diver
                case (6)
                  array(i,j,k+1,dim) = array(i,j,k-1,dim)
                  diver = div(i,j,k,array(:,:,:,1)
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
      real(8)    :: gsuper(3,3),hessian1(3,3)
     .             ,hessian2(3,3),hessian3(3,3)
      logical    :: cartesian

c     Begin program

      nvar = 1

      ibc = (1+loc)+2*(dim-1)

      rhs = 0d0

      select case (ieq)
      case (IRHO)

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
     .                           ,cartesian)

              gsuper = g_super(x1,x2,x3,cartesian)

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

              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
     .                           ,cartesian)

              gsuper = g_super (x1,x2,x3,cartesian)
              jac0   = jacobian(x1,x2,x3,cartesian)

              hessian1 = hessian(1,x1,x2,x3,cartesian)
              hessian2 = hessian(2,x1,x2,x3,cartesian)
              hessian3 = hessian(3,x1,x2,x3,cartesian)
        
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
                  rhs(j,k) =  hessian1(1,1)*vx(i,j,k)*vx(i,j,k)
     .                            +hessian1(2,2)*vy(i,j,k)*vy(i,j,k)
     .                            +hessian1(3,3)*vz(i,j,k)*vz(i,j,k)
     .                         +2.*hessian1(1,2)*vx(i,j,k)*vy(i,j,k)
     .                         +2.*hessian1(1,3)*vx(i,j,k)*vz(i,j,k)
     .                         +2.*hessian1(2,3)*vy(i,j,k)*vz(i,j,k)
                endif
                rhs(j,k) = -dh(dim)
     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
     .              -0.5/jac0*rhs(j,k))/gsuper(dim,dim)
              elseif (dim == 2) then
                if (gamma > 1d0) then
                  rhs(i,k) =  hessian2(1,1)*vx(i,j,k)*vx(i,j,k)
     .                            +hessian2(2,2)*vy(i,j,k)*vy(i,j,k)
     .                            +hessian2(3,3)*vz(i,j,k)*vz(i,j,k)
     .                         +2.*hessian2(1,2)*vx(i,j,k)*vy(i,j,k)
     .                         +2.*hessian2(1,3)*vx(i,j,k)*vz(i,j,k)
     .                         +2.*hessian2(2,3)*vy(i,j,k)*vz(i,j,k)
                endif
                rhs(i,k) = -dh(dim)
     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
     .              -0.5/jac0*rhs(i,k))/gsuper(dim,dim)
              elseif (dim == 3) then
                if (gamma > 1d0) then
                  rhs(i,j) =  hessian3(1,1)*vx(i,j,k)*vx(i,j,k)
     .                            +hessian3(2,2)*vy(i,j,k)*vy(i,j,k)
     .                            +hessian3(3,3)*vz(i,j,k)*vz(i,j,k)
     .                         +2.*hessian3(1,2)*vx(i,j,k)*vy(i,j,k)
     .                         +2.*hessian3(1,3)*vx(i,j,k)*vz(i,j,k)
     .                         +2.*hessian3(2,3)*vy(i,j,k)*vz(i,j,k)
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
        rhs(:,:) = array(1,:,:)  - rhs(:,:)
      case (2)
        rhs(:,:) = array(nx,:,:) + rhs(:,:)
      case (3)
        rhs(:,:) = array(:,1,:)  - rhs(:,:)
      case (4)
        rhs(:,:) = array(:,ny,:) + rhs(:,:)
      case (5)
        rhs(:,:) = array(:,:,1)  - rhs(:,:)
      case (6)
        rhs(:,:) = array(:,:,nz) + rhs(:,:)
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
      real(8)    :: gsuper(3,3),hessian1(3,3)
     .             ,hessian2(3,3),hessian3(3,3)
      logical    :: cartesian

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      rhs = 0d0

      select case (ieq)
      case (IVX,IVY,IVZ) !Velocity components

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
     .                           ,cartesian)

              gsuper   = g_super(x1,x2,x3,cartesian)

              hessian1 = hessian(1,x1,x2,x3,cartesian)
              hessian2 = hessian(2,x1,x2,x3,cartesian)
              hessian3 = hessian(3,x1,x2,x3,cartesian)

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
     .                      *(hessian1(ivar,1)*array(i,j,k,1)
     .                       +hessian2(ivar,1)*array(i,j,k,2)
     .                       +hessian3(ivar,1)*array(i,j,k,3))
     .                    +gsuper(dim,2)
     .                      *(hessian1(ivar,2)*array(i,j,k,1)
     .                       +hessian2(ivar,2)*array(i,j,k,2)
     .                       +hessian3(ivar,2)*array(i,j,k,3))
     .                    +gsuper(dim,3)
     .                      *(hessian1(ivar,3)*array(i,j,k,1)
     .                       +hessian2(ivar,3)*array(i,j,k,2)
     .                       +hessian3(ivar,3)*array(i,j,k,3))

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
     .                      *(hessian1(ivar,1)*array(i,j,k,1)
     .                       +hessian2(ivar,1)*array(i,j,k,2)
     .                       +hessian3(ivar,1)*array(i,j,k,3))
     .                    +gsuper(dim,2)
     .                      *(hessian1(ivar,2)*array(i,j,k,1)
     .                       +hessian2(ivar,2)*array(i,j,k,2)
     .                       +hessian3(ivar,2)*array(i,j,k,3))
     .                    +gsuper(dim,3)
     .                      *(hessian1(ivar,3)*array(i,j,k,1)
     .                       +hessian2(ivar,3)*array(i,j,k,2)
     .                       +hessian3(ivar,3)*array(i,j,k,3))

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
     .                      *(hessian1(ivar,1)*array(i,j,k,1)
     .                       +hessian2(ivar,1)*array(i,j,k,2)
     .                       +hessian3(ivar,1)*array(i,j,k,3))
     .                    +gsuper(dim,2)
     .                      *(hessian1(ivar,2)*array(i,j,k,1)
     .                       +hessian2(ivar,2)*array(i,j,k,2)
     .                       +hessian3(ivar,2)*array(i,j,k,3))
     .                    +gsuper(dim,3)
     .                      *(hessian1(ivar,3)*array(i,j,k,1)
     .                       +hessian2(ivar,3)*array(i,j,k,2)
     .                       +hessian3(ivar,3)*array(i,j,k,3))

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

              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
     .                           ,cartesian)

              gsuper = g_super(x1,x2,x3,cartesian)

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
        rhs(:,:) = array(1,:,:,ivar)  - rhs(:,:)
      case (2)
        rhs(:,:) = array(nx,:,:,ivar) + rhs(:,:)
      case (3)
        rhs(:,:) = array(:,1,:,ivar)  - rhs(:,:)
      case (4)
        rhs(:,:) = array(:,ny,:,ivar) + rhs(:,:)
      case (5)
        rhs(:,:) = array(:,:,1,ivar)  - rhs(:,:)
      case (6)
        rhs(:,:) = array(:,:,nz,ivar) + rhs(:,:)
      end select

c     End program

      end subroutine vectorNeumannBC

      end module neumannBCinterface

c imposeBoundaryConditions
c####################################################################
      subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use BCS

      use singularBCinterface

      implicit none

c Call variables

      integer(4) :: iigx,iigy,iigz

      type (var_array) :: varray

c Local variables

      integer(4) :: neq,ieq,i,j,k,icomp
      integer(4) :: dim,loc,bctype,ibc
      integer(4) :: bcnd(6,3)

      real(8) :: mag

c Begin program

      igx = iigx
      igy = iigy
      igz = iigz

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      allocate(v_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,v_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,vzeros(0:nx+1,0:ny+1,0:nz+1,3))

c Density BC

      call imposeBConScalar(IRHO,varray%array_var(IRHO)%array
     .                          ,u_0   %array_var(IRHO)%array
     .                          ,varray%array_var(IRHO)%bconds)

c Velocity BC

      vzeros = 0d0

      bcnd(:,1) = varray%array_var(IVX)%bconds
      bcnd(:,2) = varray%array_var(IVY)%bconds
      bcnd(:,3) = varray%array_var(IVZ)%bconds

      where (varray%array_var(IRHO)%array /= 0d0)
        v_cnv(:,:,:,1) = varray%array_var(IVX )%array
     .                  /varray%array_var(IRHO)%array
        v_cnv(:,:,:,2) = varray%array_var(IVY )%array
     .                  /varray%array_var(IRHO)%array
        v_cnv(:,:,:,3) = varray%array_var(IVZ )%array
     .                  /varray%array_var(IRHO)%array
      end where

      call imposeBConVector(IVX,v_cnv,v_cov,vzeros,bcnd)

      vx_cov = v_cov(:,:,:,1)
      vy_cov = v_cov(:,:,:,2)
      vz_cov = v_cov(:,:,:,3)

      vx     = v_cnv(:,:,:,1)
      vy     = v_cnv(:,:,:,2)
      vz     = v_cnv(:,:,:,3)

      varray%array_var(IVX)%array = v_cnv(:,:,:,1)
     .                             *varray%array_var(IRHO)%array
      varray%array_var(IVY)%array = v_cnv(:,:,:,2)
     .                             *varray%array_var(IRHO)%array
      varray%array_var(IVZ)%array = v_cnv(:,:,:,3)
     .                             *varray%array_var(IRHO)%array

c Magnetic field BC

c     BC setup

      bcnd(:,1) = varray%array_var(IBX)%bconds
      bcnd(:,2) = varray%array_var(IBY)%bconds
      bcnd(:,3) = varray%array_var(IBZ)%bconds

      v_cnv(:,:,:,1) = varray%array_var(IBX)%array
      v_cnv(:,:,:,2) = varray%array_var(IBY)%array
      v_cnv(:,:,:,3) = varray%array_var(IBZ)%array

      vzeros(:,:,:,1) = u_0%array_var(IBX)%array
      vzeros(:,:,:,2) = u_0%array_var(IBY)%array
      vzeros(:,:,:,3) = u_0%array_var(IBZ)%array

c     Fill ghost nodes

      call imposeBConVector(IBX,v_cnv,v_cov,vzeros,bcnd)

c     Postprocessing

      bx_cov = v_cov(:,:,:,1)
      by_cov = v_cov(:,:,:,2)
      bz_cov = v_cov(:,:,:,3)

      varray%array_var(IBX)%array = v_cnv(:,:,:,1)
      varray%array_var(IBY)%array = v_cnv(:,:,:,2)
      varray%array_var(IBZ)%array = v_cnv(:,:,:,3)

c Current BC

c     BC setup

      bcnd(:,1) = varray%array_var(IBX)%bconds
      bcnd(:,2) = varray%array_var(IBY)%bconds
      bcnd(:,3) = varray%array_var(IBZ)%bconds
      where (bcnd == -NEU)
        bcnd = -DIR  !Use covariant components for tangential dirichlet
      end where

      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            do icomp=1,3
              v_cnv(i,j,k,icomp)=curl2(i,j,k,bx_cov,by_cov,bz_cov,icomp)
            enddo
          enddo
        enddo
      enddo

      vzeros = v_cnv

c     Fill ghost nodes

      call imposeBConVector(IJX,v_cnv,v_cov,vzeros,bcnd)

c     Postprocessing

      jx_cov = v_cov(:,:,:,1)
      jy_cov = v_cov(:,:,:,2)
      jz_cov = v_cov(:,:,:,3)

      jx = v_cnv(:,:,:,1)
      jy = v_cnv(:,:,:,2)
      jz = v_cnv(:,:,:,3)

c Temperature BCs

      call imposeBConScalar(ITMP,varray%array_var(ITMP)%array
     .                          ,u_0   %array_var(ITMP)%array
     .                          ,varray%array_var(ITMP)%bconds)

c Deallocate variables 

      deallocate(v_cnv,v_cov,vzeros)

c End

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

      bctype=PER

      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)
          if (bcond(ibc) == bctype) then
            call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
          endif
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

c     Preprocess velocity field

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

      if (cov_to_cnv) call synchronize(v_cnv,v_cov,bcond)

c     Impose vector singular point BCs

      if (bcond(1,1) == SP) call singularBC(v_cnv,.false.,2)

c     Synchronize periodic boundaries

      bctype=PER

      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)
          do ivar = 1,3
            ieq = ivar + fcomp - 1
            if (bcond(ibc,ivar) == bctype) then
              call FillGhostNodes(ieq,ivar,3,dim,loc,bctype,v_cnv,var0)
cc              call FillGhostNodes(ieq,ivar,3,dim,loc,bctype,v_cov,var0)
            endif
          enddo
        enddo
      enddo

c     Find covariant components at ALL boundaries

      do dim = 1,3
        do loc = 0,1
          if (dim == 1) then
            imin=0 + loc*(nx+1)
            imax=0 + loc*(nx+1)
            jmin=0
            jmax=ny+1
            kmin=0
            kmax=nz+1
          elseif (dim == 2) then
            imin=0
            imax=nx+1
            jmin=0 + loc*(ny+1)
            jmax=0 + loc*(ny+1)
            kmin=0
            kmax=nz+1
          elseif (dim == 3) then
            imin=0
            imax=nx+1
            jmin=0
            jmax=ny+1
            kmin=0 + loc*(nz+1)
            kmax=0 + loc*(nz+1)
          endif

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

c     synchronize
c     #################################################################
      subroutine synchronize(v_cnv,v_cov,bcond)
c     -----------------------------------------------------------------
c     Finds all contravariant components at Neumann boundaries.
c     On input, tangential covariant components and normal contravariant
c     components are known at ghost cells. On output, all contravarian
c     components are known at ghost cells.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: bcond(6,3)
      real(8)    :: v_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,v_cov(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

      integer(4) :: i,j,k,dim,loc,ig,jg,kg,ivar
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      real(8)    :: x1,x2,x3,gsuper(3,3),gsub(3,3)
      logical    :: cartesian

c     Begin program

      do dim = 1,3
        do loc = 0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,3
            if (ivar == dim) then  !Select tangential components
              cycle
            elseif (bcond(ibc,ivar) < 0) then
              if (dim == 1) then
                imin=1  +    loc *(nx-1)
                imax=nx + (1-loc)*(1-nx)
                jmin=1
                jmax=ny
                kmin=1
                kmax=nz
              elseif (dim == 2) then
                imin=1 
                imax=nx
                jmin=1  +    loc *(ny-1)
                jmax=ny + (1-loc)*(1-ny)
                kmin=1
                kmax=nz
              elseif (dim == 3) then
                imin=1 
                imax=nx
                jmin=1
                jmax=ny
                kmin=1  +    loc *(nz-1)
                kmax=nz + (1-loc)*(1-nz)
              endif

              select case (ibc)
              case (1)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      call getCoordinates(i-1,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,x2,x3,cartesian)

                      gsuper = g_super(x1,x2,x3,cartesian)

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

                      call getCoordinates(i+1,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,x2,x3,cartesian)

                      gsuper = g_super(x1,x2,x3,cartesian)

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

                      call getCoordinates(i,j-1,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,x2,x3,cartesian)

                      gsuper = g_super(x1,x2,x3,cartesian)

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

                      call getCoordinates(i,j+1,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,x2,x3,cartesian)

                      gsuper = g_super(x1,x2,x3,cartesian)

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

                      call getCoordinates(i,j,k-1,igx,igy,igz,ig,jg,kg
     .                                   ,x1,x2,x3,cartesian)

                      gsuper = g_super(x1,x2,x3,cartesian)

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

                      call getCoordinates(i,j,k+1,igx,igy,igz,ig,jg,kg
     .                                   ,x1,x2,x3,cartesian)

                      gsuper = g_super(x1,x2,x3,cartesian)

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

      end subroutine synchronize

      end subroutine imposeBoundaryConditions


c FillGhostNodes
c####################################################################
      subroutine FillGhostNodes(ieq,ivar,nvar,dim,loc,bctype
     .                         ,array,array0)

c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array.
c
c     On input:
c       * ieq    -> equation identifier
c       * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c       * loc    -> location in dimension (0 -> right, 1 -> left)
c       * bctype -> type of BC (dirichlet, neumann, periodic, etc.)
c       * array  -> real array with ghost-nodes
c       * array0 -> auxiliary real array
c--------------------------------------------------------------------

      use BCS

      use singularBCinterface

      use dirichletBCinterface

      use neumannBCinterface

      implicit none       !For safe fortran

c Call variables

      integer(4) :: ieq,dim,loc,bctype,nvar,ivar
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1,nvar)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1,nvar)

c Local variables

      integer(4) :: neq,ibc
      integer(4) :: i,j,k,ig,jg,kg

c Begin program

      nnvar = nvar

cc      write (*,*) 'FillGhostNodes, nvar=',nvar

c Determine boundary limits

      if (dim == 1) then
        imin=1  +    loc *(nx-1)
        imax=nx + (1-loc)*(1-nx)
        jmin=1
        jmax=ny
        kmin=1
        kmax=nz
        allocate(rhs(0:ny+1,0:nz+1))
      elseif (dim == 2) then
        imin=1 
        imax=nx
        jmin=1  +    loc *(ny-1)
        jmax=ny + (1-loc)*(1-ny)
        kmin=1
        kmax=nz
        allocate(rhs(0:nx+1,0:nz+1))
      elseif (dim == 3) then
        imin=1 
        imax=nx
        jmin=1
        jmax=ny
        kmin=1  +    loc *(nz-1)
        kmax=nz + (1-loc)*(1-nz)
        allocate(rhs(0:nx+1,0:ny+1))
      endif

c Find BC update

      ibc = (1+loc)+2*(dim-1)

      if (nvar == 1) then
        ivar = nvar
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,ivar),ibc)
        case(EQU)
          call dirichletBC(array(:,:,:,ivar),array0(:,:,:,ivar)
     .                    ,ieq,dim,loc,0)
        case(DIR)
          call dirichletBC(array(:,:,:,ivar),array0(:,:,:,ivar)
     .                    ,ieq,dim,loc,1)
        case(NEU)
          call neumannBC(array(:,:,:,ivar),ieq,dim,loc)
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
          call dirichletBC(ivar,array,array0,ieq,dim,loc,1)
        case(NEU)
          call neumannBC(ivar,array,ieq,dim,loc)
        case default
          write (*,*) 'BC',bctype,' not implemented'
          stop
        end select
      endif

c Update BC ghost nodes

      select case (ibc)
      case (1)                  !x0
        array(0   ,:,:,ivar) = rhs(:,:)
      case (2)                  !x1
        array(nx+1,:,:,ivar) = rhs(:,:)
      case (3)                  !y0
        array(:,0   ,:,ivar) = rhs(:,:)
      case (4)                  !y1
        array(:,ny+1,:,ivar) = rhs(:,:)
      case (5)                  !z0
        array(:,:,0   ,ivar) = rhs(:,:)
      case (6)                  !z1
        array(:,:,nz+1,ivar) = rhs(:,:)
      case default
        write (*,*) 'Boundary',ibc,' non existent'
        stop
      end select

      deallocate(rhs)

c End

      contains

c     periodicBC
c     #################################################################
      subroutine periodicBC(array,ibc)
c     -----------------------------------------------------------------
c     Imposes singular point BC. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ibc,ivar
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)

c     Local variables

c     Begin program

      select case (ibc)
      case (1)
        rhs(:,:) = array(nx,:,:)
      case (2)
        rhs(:,:) = array(1,:,:)
      case (3)
        rhs(:,:) = array(:,ny,:)
      case (4)
        rhs(:,:) = array(:,1,:)
      case (5)
        rhs(:,:) = array(:,:,nz)
      case (6)
        rhs(:,:) = array(:,:,1)
      end select

c     End program

      end subroutine periodicBC

ccc     neumannBC
ccc     #################################################################
cc      subroutine neumannBC(array,ieq,ibc)
ccc     -----------------------------------------------------------------
ccc     Imposes neumann BC for a scalar. On input:
ccc        * ieq -> equation number (i.e., vector component)
ccc        * dim -> dimension we are imposing BC on (X,Y,Z)
ccc        * loc -> boundary location (0 -> left, 1->right)
ccc     This routine fills up the bi-dimensional array rhs, which 
ccc     contains the right hand side of the Neumann BC.
ccc     -----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer(4) :: ieq,ibc
cc      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc      integer(4) :: i,j,k,ip,im,jp,jm,kp,km,icomp
cc      real(8)    :: x1,x2,x3,dh(3),jac0
cc      real(8)    :: gsuper(3,3),hessian1(3,3)
cc     .             ,hessian2(3,3),hessian3(3,3)
cc      logical    :: cartesian
cc
ccc     Begin program
cc
cc      rhs = 0d0
cc
cc      select case (ieq)
cc      case (IRHO)
cc
cc        do i=imin,imax
cc          do j=jmin,jmax
cc            do k=kmin,kmax
cc
cc              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
cc     .                           ,cartesian)
cc
cc              gsuper = g_super(x1,x2,x3,cartesian)
cc
cc              ip = min(i+1,nx)
cc              im = max(i-1,1)
cc              jp = min(j+1,ny)
cc              jm = max(j-1,1)
cc              kp = min(k+1,nz)
cc              km = max(k-1,1)
cc
cc              dh(1) = 2.*dxh(ig)
cc              if (i == nx) dh(1) = dx(ig-1)
cc              if (i == 1 ) dh(1) = dx(ig)
cc
cc              dh(2) = 2.*dyh(jg)
cc              if (j == ny) dh(2) = dy(jg-1)
cc              if (j == 1 ) dh(2) = dy(jg)
cc
cc              dh(3) = 2.*dzh(kg)
cc              if (k == nz) dh(3) = dz(kg-1)
cc              if (k == 1 ) dh(3) = dz(kg)
cc
cc              if (dim == 1) then
cc                rhs(j,k) = -dh(dim)
cc     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
cc     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3))
cc     .              /gsuper(dim,dim)
cc              elseif (dim == 2) then
cc                rhs(i,k) = -dh(dim)*
cc     .              (gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
cc     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3))
cc     .              /gsuper(dim,dim)
cc              elseif (dim == 3) then
cc                rhs(i,j) = -dh(dim)
cc     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
cc     .              +gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2))
cc     .              /gsuper(dim,dim)
cc              endif
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc      case (ITMP)
cc
cc        do i=imin,imax
cc          do j=jmin,jmax
cc            do k=kmin,kmax
cc
cc              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
cc     .                           ,cartesian)
cc
cc              gsuper = g_super (x1,x2,x3,cartesian)
cc              jac0   = jacobian(x1,x2,x3,cartesian)
cc
cc              hessian1 = hessian(1,x1,x2,x3,cartesian)
cc              hessian2 = hessian(2,x1,x2,x3,cartesian)
cc              hessian3 = hessian(3,x1,x2,x3,cartesian)
cc        
cc              ip = min(i+1,nx)
cc              im = max(i-1,1)
cc              jp = min(j+1,ny)
cc              jm = max(j-1,1)
cc              kp = min(k+1,nz)
cc              km = max(k-1,1)
cc
cc              dh(1) = 2.*dxh(ig)
cc              if (i == nx) dh(1) = dx(ig-1)
cc              if (i == 1 ) dh(1) = dx(ig)
cc
cc              dh(2) = 2.*dyh(jg)
cc              if (j == ny) dh(2) = dy(jg-1)
cc              if (j == 1 ) dh(2) = dy(jg)
cc
cc              dh(3) = 2.*dzh(kg)
cc              if (k == nz) dh(3) = dz(kg-1)
cc              if (k == 1 ) dh(3) = dz(kg)
cc
cc              if (dim == 1) then
cc                if (gamma > 1d0) then
cc                  rhs(j,k) =  hessian1(1,1)*vx(i,j,k)*vx(i,j,k)
cc     .                       +hessian1(2,2)*vy(i,j,k)*vy(i,j,k)
cc     .                       +hessian1(3,3)*vz(i,j,k)*vz(i,j,k)
cc     .                    +2.*hessian1(1,2)*vx(i,j,k)*vy(i,j,k)
cc     .                    +2.*hessian1(1,3)*vx(i,j,k)*vz(i,j,k)
cc     .                    +2.*hessian1(2,3)*vy(i,j,k)*vz(i,j,k)
cc                endif
cc                rhs(j,k) = -dh(dim)
cc     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
cc     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
cc     .              -0.5/jac0*rhs(j,k))/gsuper(dim,dim)
cc              elseif (dim == 2) then
cc                if (gamma > 1d0) then
cc                  rhs(i,k) =  hessian2(1,1)*vx(i,j,k)*vx(i,j,k)
cc     .                       +hessian2(2,2)*vy(i,j,k)*vy(i,j,k)
cc     .                       +hessian2(3,3)*vz(i,j,k)*vz(i,j,k)
cc     .                    +2.*hessian2(1,2)*vx(i,j,k)*vy(i,j,k)
cc     .                    +2.*hessian2(1,3)*vx(i,j,k)*vz(i,j,k)
cc     .                    +2.*hessian2(2,3)*vy(i,j,k)*vz(i,j,k)
cc                endif
cc                rhs(i,k) = -dh(dim)
cc     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
cc     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
cc     .              -0.5/jac0*rhs(i,k))/gsuper(dim,dim)
cc              elseif (dim == 3) then
cc                if (gamma > 1d0) then
cc                  rhs(i,j) =  hessian3(1,1)*vx(i,j,k)*vx(i,j,k)
cc     .                       +hessian3(2,2)*vy(i,j,k)*vy(i,j,k)
cc     .                       +hessian3(3,3)*vz(i,j,k)*vz(i,j,k)
cc     .                    +2.*hessian3(1,2)*vx(i,j,k)*vy(i,j,k)
cc     .                    +2.*hessian3(1,3)*vx(i,j,k)*vz(i,j,k)
cc     .                    +2.*hessian3(2,3)*vy(i,j,k)*vz(i,j,k)
cc                endif
cc                rhs(i,j) = -dh(dim)
cc     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
cc     .              +gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
cc     .              -0.5/jac0*rhs(i,j))/gsuper(dim,dim)
cc              endif
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc      case (IVX,IVY,IVZ) !Velocity components
cc
cc        if (ieq == IVX) icomp = 1
cc        if (ieq == IVY) icomp = 2
cc        if (ieq == IVZ) icomp = 3
cc
cccc        do i=imin,imax
cccc          do j=jmin,jmax
cccc            do k=kmin,kmax
cccc
cccc              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
cccc     .                           ,cartesian)
cccc
cccc              gsuper = g_super(x1,x2,x3,cartesian)
cccc
cccccc              if (loc == 0) then
cccccc                if (dim==1) then
cccccc                  nabla_v = fnabla_v_bc(i-1,j,k,dim)
cccccc                elseif (dim==2) then
cccccc                  nabla_v = fnabla_v_bc(i,j-1,k,dim)
cccccc                else
cccccc                  nabla_v = fnabla_v_bc(i,j,k-1,dim)
cccccc                endif
cccccc              else
cccc                nabla_v = fnabla_v_bc(i,j,k,x1,x2,x3,cartesian,dim)
cccccc              endif
cccc
cccc              dh(1) = 2.*dxh(ig)
cccc              if (i == nx) dh(1) = dx(ig-1)
cccc              if (i == 1 ) dh(1) = dx(ig)
cccc
cccc              dh(2) = 2.*dyh(jg)
cccc              if (j == ny) dh(2) = dy(jg-1)
cccc              if (j == 1 ) dh(2) = dy(jg)
cccc
cccc              dh(3) = 2.*dzh(kg)
cccc              if (k == nz) dh(3) = dz(kg-1)
cccc              if (k == 1 ) dh(3) = dz(kg)
cccc
cccc              if (dim == 1) then
cccc                rhs(j,k) = -dh(dim)
cccc     .               *(gsuper(dim,1)*nabla_v(1,icomp)
cccc     .                +gsuper(dim,2)*nabla_v(2,icomp)
cccc     .                +gsuper(dim,3)*nabla_v(3,icomp))
cccc     .                /gsuper(dim,dim)
cccc              elseif (dim == 2) then
cccc                rhs(i,k) = -dh(dim)
cccc     .               *(gsuper(dim,1)*nabla_v(1,icomp)
cccc     .                +gsuper(dim,2)*nabla_v(2,icomp)
cccc     .                +gsuper(dim,3)*nabla_v(3,icomp))
cccc     .                /gsuper(dim,dim)
cccc              elseif (dim == 3) then
cccc                rhs(i,j) = -dh(dim)
cccc     .               *(gsuper(dim,1)*nabla_v(1,icomp)
cccc     .                +gsuper(dim,2)*nabla_v(2,icomp)
cccc     .                +gsuper(dim,3)*nabla_v(3,icomp))
cccc     .                /gsuper(dim,dim)
cccc              endif
cccc
cccc            enddo
cccc          enddo
cccc        enddo
cc
cc        do i=imin,imax
cc          do j=jmin,jmax
cc            do k=kmin,kmax
cc
cc              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
cc     .                           ,cartesian)
cc
cc              gsuper   = g_super(x1,x2,x3,cartesian)
cc
cc              hessian1 = hessian(1,x1,x2,x3,cartesian)
cc              hessian2 = hessian(2,x1,x2,x3,cartesian)
cc              hessian3 = hessian(3,x1,x2,x3,cartesian)
cc
cc              ip = min(i+1,nx)
cc              im = max(i-1,1)
cc              jp = min(j+1,ny)
cc              jm = max(j-1,1)
cc              kp = min(k+1,nz)
cc              km = max(k-1,1)
cc
cc              dh(1) = 2.*dxh(ig)
cc              if (i == nx) dh(1) = dx(ig-1)
cc              if (i == 1 ) dh(1) = dx(ig)
cc
cc              dh(2) = 2.*dyh(jg)
cc              if (j == ny) dh(2) = dy(jg-1)
cc              if (j == 1 ) dh(2) = dy(jg)
cc
cc              dh(3) = 2.*dzh(kg)
cc              if (k == nz) dh(3) = dz(kg-1)
cc              if (k == 1 ) dh(3) = dz(kg)
cc
cc
cc              if (dim == 1) then
cc
cc                rhs(j,k) = gsuper(dim,1)
cc     .                      *(hessian1(icomp,1)*vx_cov(i,j,k)
cc     .                       +hessian2(icomp,1)*vy_cov(i,j,k)
cc     .                       +hessian3(icomp,1)*vz_cov(i,j,k))
cc     .                    +gsuper(dim,2)
cc     .                      *(hessian1(icomp,2)*vx_cov(i,j,k)
cc     .                       +hessian2(icomp,2)*vy_cov(i,j,k)
cc     .                       +hessian3(icomp,2)*vz_cov(i,j,k))
cc     .                    +gsuper(dim,3)
cc     .                      *(hessian1(icomp,3)*vx_cov(i,j,k)
cc     .                       +hessian2(icomp,3)*vy_cov(i,j,k)
cc     .                       +hessian3(icomp,3)*vz_cov(i,j,k))
cc
cc                if (icomp == 2) then
cc
cc                  rhs(j,k) = -dh(dim)
cc     .             *(gsuper(dim,2)*(vy_cov(i,jp,k)-vy_cov(i,jm,k))/dh(2)
cc     .              +gsuper(dim,3)*(vy_cov(i,j,kp)-vy_cov(i,j,km))/dh(3)
cc     .              +rhs(j,k))/gsuper(dim,dim)
cc
cc                  if (loc == 0) then
cc                    vy_cov(i-1,j,k) = vy_cov(i,j,k) - rhs(j,k)
cc                  else
cc                    vy_cov(i+1,j,k) = vy_cov(i,j,k) + rhs(j,k)
cc                  endif
cc
cc                elseif (icomp == 3) then
cc
cc                  rhs(j,k) = -dh(dim)
cc     .             *(gsuper(dim,2)*(vz_cov(i,jp,k)-vz_cov(i,jm,k))/dh(2)
cc     .              +gsuper(dim,3)*(vz_cov(i,j,kp)-vz_cov(i,j,km))/dh(3)
cc     .              +rhs(j,k))/gsuper(dim,dim)
cc
cc                  if (loc == 0) then
cc                    vz_cov(i-1,j,k) = vz_cov(i,j,k) - rhs(j,k)
cc                  else
cc                    vz_cov(i+1,j,k) = vz_cov(i,j,k) + rhs(j,k)
cc                  endif
cc                endif
cc
cc              elseif (dim == 2) then
cc
cc                rhs(i,k) = gsuper(dim,1)
cc     .                      *(hessian1(icomp,1)*vx_cov(i,j,k)
cc     .                       +hessian2(icomp,1)*vy_cov(i,j,k)
cc     .                       +hessian3(icomp,1)*vz_cov(i,j,k))
cc     .                    +gsuper(dim,2)
cc     .                      *(hessian1(icomp,2)*vx_cov(i,j,k)
cc     .                       +hessian2(icomp,2)*vy_cov(i,j,k)
cc     .                       +hessian3(icomp,2)*vz_cov(i,j,k))
cc     .                    +gsuper(dim,3)
cc     .                      *(hessian1(icomp,3)*vx_cov(i,j,k)
cc     .                       +hessian2(icomp,3)*vy_cov(i,j,k)
cc     .                       +hessian3(icomp,3)*vz_cov(i,j,k))
cc
cc                if (icomp == 3) then
cc
cc                  rhs(i,k) = -dh(dim)
cc     .             *(gsuper(dim,1)*(vz_cov(ip,j,k)-vz_cov(im,j,k))/dh(1)
cc     .              +gsuper(dim,3)*(vz_cov(i,j,kp)-vz_cov(i,j,km))/dh(3)
cc     .              +rhs(i,k))/gsuper(dim,dim)
cc
cc                  if (loc == 0) then
cc                    vz_cov(i,j-1,k) = vz_cov(i,j,k) - rhs(i,k)
cc                  else
cc                    vz_cov(i,j+1,k) = vz_cov(i,j,k) + rhs(i,k)
cc                  endif
cc
cc                elseif (icomp == 1) then
cc
cc                  rhs(i,k) = -dh(dim)
cc     .             *(gsuper(dim,1)*(vx_cov(ip,j,k)-vx_cov(im,j,k))/dh(1)
cc     .              +gsuper(dim,3)*(vx_cov(i,j,kp)-vx_cov(i,j,km))/dh(3)
cc     .              +rhs(i,k))/gsuper(dim,dim)
cc
cc                  if (loc == 0) then
cc                    vx_cov(i,j-1,k) = vx_cov(i,j,k) - rhs(i,k)
cc                  else
cc                    vx_cov(i,j+1,k) = vx_cov(i,j,k) + rhs(i,k)
cc                  endif
cc
cc                endif
cc
cc              elseif (dim == 3) then
cc
cc                rhs(i,j) = gsuper(dim,1)
cc     .                      *(hessian1(icomp,1)*vx_cov(i,j,k)
cc     .                       +hessian2(icomp,1)*vy_cov(i,j,k)
cc     .                       +hessian3(icomp,1)*vz_cov(i,j,k))
cc     .                    +gsuper(dim,2)
cc     .                      *(hessian1(icomp,2)*vx_cov(i,j,k)
cc     .                       +hessian2(icomp,2)*vy_cov(i,j,k)
cc     .                       +hessian3(icomp,2)*vz_cov(i,j,k))
cc     .                    +gsuper(dim,3)
cc     .                      *(hessian1(icomp,3)*vx_cov(i,j,k)
cc     .                       +hessian2(icomp,3)*vy_cov(i,j,k)
cc     .                       +hessian3(icomp,3)*vz_cov(i,j,k))
cc
cc                if (icomp == 1) then
cc
cc                  rhs(i,j) = -dh(dim)
cc     .             *(gsuper(dim,1)*(vx_cov(ip,j,k)-vx_cov(im,j,k))/dh(1)
cc     .              +gsuper(dim,2)*(vx_cov(i,jp,k)-vx_cov(i,jm,k))/dh(2)
cc     .              +rhs(i,j))/gsuper(dim,dim)
cc
cc                  if (loc == 0) then
cc                    vx_cov(i,j,k-1) = vx_cov(i,j,k) - rhs(i,j)
cc                  else
cc                    vx_cov(i,j,k+1) = vx_cov(i,j,k) + rhs(i,j)
cc                  endif
cc
cc                elseif (icomp == 2) then
cc
cc                  rhs(i,j) = -dh(dim)
cc     .             *(gsuper(dim,1)*(vy_cov(ip,j,k)-vy_cov(im,j,k))/dh(1)
cc     .              +gsuper(dim,2)*(vy_cov(i,jp,k)-vy_cov(i,jm,k))/dh(2)
cc     .              +rhs(i,j))/gsuper(dim,dim)
cc
cc                  if (loc == 0) then
cc                    vy_cov(i,j,k-1) = vy_cov(i,j,k) - rhs(i,j)
cc                  else
cc                    vy_cov(i,j,k+1) = vy_cov(i,j,k) + rhs(i,j)
cc                  endif
cc                endif
cc
cc              endif
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc      case (IBX,IBY,IBZ) 
cc
cc        if (ieq == IBX) icomp = 1
cc        if (ieq == IBY) icomp = 2
cc        if (ieq == IBZ) icomp = 3
cc
cc        if (icomp == dim) then
cc          write (*,*) 'Error in B in neumannBC: icomp = dim=',dim
cc          stop
cc        endif
cc
cc        do i=imin,imax
cc          do j=jmin,jmax
cc            do k=kmin,kmax
cc
cc              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
cc     .                           ,cartesian)
cc
cc              gsuper = g_super(x1,x2,x3,cartesian)
cc
cc              ip = min(i+1,nx)
cc              im = max(i-1,1)
cc              jp = min(j+1,ny)
cc              jm = max(j-1,1)
cc              kp = min(k+1,nz)
cc              km = max(k-1,1)
cc
cc              dh(1) = 2.*dxh(ig)
cc              if (i == nx) dh(1) = dx(ig-1)
cc              if (i == 1 ) dh(1) = dx(ig)
cc
cc              dh(2) = 2.*dyh(jg)
cc              if (j == ny) dh(2) = dy(jg-1)
cc              if (j == 1 ) dh(2) = dy(jg)
cc
cc              dh(3) = 2.*dzh(kg)
cc              if (k == nz) dh(3) = dz(kg-1)
cc              if (k == 1 ) dh(3) = dz(kg)
cc
cc              if (dim == 1) then
cc
cc                jx(i,j,k) = (bz_cov(i,jp,k)-bz_cov(i,jm,k))/dh(2)
cc     .                     -(by_cov(i,j,kp)-by_cov(i,j,km))/dh(3)
cc
cc                if (icomp == 2) then
cc                  rhs(j,k) = dh(dim)
cc     .               *( (bx_cov(i,jp,k)-bx_cov(i,jm,k))/dh(2)
cc     .                 +gsuper(dim,3)*jx(i,j,k)/gsuper(dim,dim) )
cc                  if (loc == 0) then
cc                    by_cov(i-1,j,k) = by_cov(i,j,k) - rhs(j,k)
cc                  else
cc                    by_cov(i+1,j,k) = by_cov(i,j,k) + rhs(j,k)
cc                  endif
cc                elseif (icomp == 3) then
cc                  rhs(j,k) = dh(dim)
cc     .               *( (bx_cov(i,j,kp)-bx_cov(i,j,km))/dh(3)
cc     .                 -gsuper(dim,2)*jx(i,j,k)/gsuper(dim,dim) )
cc                  if (loc == 0) then
cc                    bz_cov(i-1,j,k) = bz_cov(i,j,k) - rhs(j,k)
cc                  else
cc                    bz_cov(i+1,j,k) = bz_cov(i,j,k) + rhs(j,k)
cc                  endif
cc                endif
cc
cc              elseif (dim == 2) then
cc
cc                jy(i,j,k) = (bx_cov(i,j,kp)-bx_cov(i,j,km))/dh(3)
cc     .                     -(bz_cov(ip,j,k)-bz_cov(im,j,k))/dh(1)
cc
cc                if (icomp == 3) then
cc                  rhs(i,k) = dh(dim)
cc     .               *( (by_cov(i,j,kp)-by_cov(i,j,km))/dh(3)
cc     .                 +gsuper(dim,1)*jy(i,j,k)/gsuper(dim,dim) )
cc                  if (loc == 0) then
cc                    bz_cov(i,j-1,k) = bz_cov(i,j,k) - rhs(i,k)
cc                  else
cc                    bz_cov(i,j+1,k) = bz_cov(i,j,k) + rhs(i,k)
cc                  endif
cc                elseif (icomp == 1) then
cc                  rhs(i,k) = dh(dim)
cc     .               *( (by_cov(ip,j,k)-by_cov(im,j,k))/dh(1)
cc     .                 -gsuper(dim,3)*jy(i,j,k)/gsuper(dim,dim) )
cc                  if (loc == 0) then
cc                    bx_cov(i,j-1,k) = bx_cov(i,j,k) - rhs(i,k)
cc                  else
cc                    bx_cov(i,j+1,k) = bx_cov(i,j,k) + rhs(i,k)
cc                  endif
cc                endif
cc
cc              elseif (dim == 3) then
cc
cc                jz(i,j,k) = (by_cov(ip,j,k)-by_cov(im,j,k))/dh(1)
cc     .                     -(bx_cov(i,jp,k)-bx_cov(i,jm,k))/dh(2)
cc
cc                if (icomp == 1) then
cc                  rhs(i,j) = dh(dim)
cc     .               *( (bz_cov(ip,j,k)-bz_cov(im,j,k))/dh(1)
cc     .                 +gsuper(dim,2)*jz(i,j,k)/gsuper(dim,dim) )
cc                  if (loc == 0) then
cc                    bx_cov(i,j,k-1) = bx_cov(i,j,k) - rhs(i,j)
cc                  else
cc                    bx_cov(i,j,k+1) = bx_cov(i,j,k) + rhs(i,j)
cc                  endif
cc                elseif (icomp == 2) then
cc                  rhs(i,j) = dh(dim)
cc     .               *( (bz_cov(i,jp,k)-bz_cov(i,jm,k))/dh(2)
cc     .                 -gsuper(dim,1)*jz(i,j,k)/gsuper(dim,dim) )
cc                  if (loc == 0) then
cc                    by_cov(i,j,k-1) = by_cov(i,j,k) - rhs(i,j)
cc                  else
cc                    by_cov(i,j,k+1) = by_cov(i,j,k) + rhs(i,j)
cc                  endif
cc                endif
cc
cc              endif
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc      case (-IBX,-IBY,-IBZ) !Finds current components at boundaries
cc
cc        if (ieq == -IBX) icomp = 1
cc        if (ieq == -IBY) icomp = 2
cc        if (ieq == -IBZ) icomp = 3
cc
cc        if (icomp /= dim) then !Tangential components
cc
cc          do i=imin,imax
cc            do j=jmin,jmax
cc              do k=kmin,kmax
cc
cc              select case (ibc)
cc              case (1)
cc                call getCoordinates(i-1,j,k,igx,igy,igz,ig,jg,kg
cc     .                             ,x1,x2,x3,cartesian)
cc                gsuper = g_super(x1,x2,x3,cartesian)
cc
cc                rhs(j,k) =-gsuper(dim,icomp)/gsuper(dim,dim)*jx(i-1,j,k)
cc     .                    +array(1,j,k)
cc              case (2)
cc                call getCoordinates(i+1,j,k,igx,igy,igz,ig,jg,kg
cc     .                             ,x1,x2,x3,cartesian)
cc                gsuper = g_super(x1,x2,x3,cartesian)
cc
cc                rhs(j,k) = gsuper(dim,icomp)/gsuper(dim,dim)*jx(i+1,j,k)
cc     .                    -array(nx,j,k)
cc              case (3)
cc                call getCoordinates(i,j-1,k,igx,igy,igz,ig,jg,kg
cc     .                             ,x1,x2,x3,cartesian)
cc                gsuper = g_super(x1,x2,x3,cartesian)
cc
cc                rhs(i,k) =-gsuper(dim,icomp)/gsuper(dim,dim)*jy(i,j-1,k)
cc     .                    +array(i,1,k)
cc              case (4)
cc                call getCoordinates(i,j+1,k,igx,igy,igz,ig,jg,kg
cc     .                             ,x1,x2,x3,cartesian)
cc                gsuper = g_super(x1,x2,x3,cartesian)
cc
cc                rhs(i,k) = gsuper(dim,icomp)/gsuper(dim,dim)*jy(i,j+1,k)
cc     .                    -array(i,ny,k)
cc              case (5)
cc                call getCoordinates(i,j,k-1,igx,igy,igz,ig,jg,kg
cc     .                             ,x1,x2,x3,cartesian)
cc                gsuper = g_super(x1,x2,x3,cartesian)
cc
cc                rhs(i,j) =-gsuper(dim,icomp)/gsuper(dim,dim)*jz(i,j,k-1)
cc     .                    +array(i,j,1)
cc              case (6)
cc                call getCoordinates(i,j,k+1,igx,igy,igz,ig,jg,kg
cc     .                             ,x1,x2,x3,cartesian)
cc                gsuper = g_super(x1,x2,x3,cartesian)
cc
cc                rhs(i,j) = gsuper(dim,icomp)/gsuper(dim,dim)*jz(i,j,k+1)
cc     .                    -array(i,j,nz)
cc              end select
cc
cc              enddo
cc            enddo
cc          enddo
cc
cc        endif
cc
cc      case default
cc
cc        write (*,*) 'Error in neumannBC'
cc        stop
cc
cc      end select
cc
ccc     Assign value
cc
cc      select case (ibc)
cc      case (1)
cc        rhs(:,:) = array(1,:,:)  - rhs(:,:)
cc      case (2)
cc        rhs(:,:) = array(nx,:,:) + rhs(:,:)
cc      case (3)
cc        rhs(:,:) = array(:,1,:)  - rhs(:,:)
cc      case (4)
cc        rhs(:,:) = array(:,ny,:) + rhs(:,:)
cc      case (5)
cc        rhs(:,:) = array(:,:,1)  - rhs(:,:)
cc      case (6)
cc        rhs(:,:) = array(:,:,nz) + rhs(:,:)
cc      end select
cc
ccc     End program
cc
cc      end subroutine neumannBC

ccc     dirichletBC
ccc     #################################################################
cc      subroutine dirichletBC(array,array0,ieq,ibc,order)
ccc     -----------------------------------------------------------------
ccc     Imposes dirichlet BC. On input:
ccc        * ieq -> equation number (i.e., vector component)
ccc        * dim -> dimension we are imposing BC on (X,Y,Z)
ccc        * loc -> boundary location (0 -> left, 1->right)
ccc        * order -> order of extrapolation (when used)
ccc     This routine fills up the bi-dimensional array rhs, which 
ccc     contains the right hand side of the Dirichlet BC.
ccc     -----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer(4) :: ieq,ibc,order
cc      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,array0(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc      integer(4) :: i,j,k,ip,im,jp,jm,kp,km,icomp
cc      real(8)    :: x1,x2,x3,dh(3),nabla_v(3,3),jac
cc      logical    :: cartesian
cc
ccc     Begin program
cc
cc      rhs = 0d0
cc
cc      select case (ieq)
cc      case (IRHO,ITMP,IVX,IVY,IVZ)
cc
cc        call interpolate(array,array0,ibc,order)
cc
cc      case (IBX,IBY,IBZ) !Imposes divergence-free constraint on B-field
cc
cc        if (ieq == IBX) icomp = 1
cc        if (ieq == IBY) icomp = 2
cc        if (ieq == IBZ) icomp = 3
cc
cc        if (icomp /= dim) then
cc
cc          call interpolate(array,array0,ibc,order)
cc
cc        else
cc
cc          do i=imin,imax
cc            do j=jmin,jmax
cc              do k=kmin,kmax
cc
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              dh(1) = 2.*dxh(ig)
cc              dh(2) = 2.*dyh(jg)
cc              dh(3) = 2.*dzh(kg)
cc
cc              select case (ibc)
cc              case (1)
cc                array(i-1,j,k) = array(i+1,j,k)
cc                rhs(j,k) = array(i+1,j,k) + dh(1)*div(i,j,k,array,by,bz)
cc              case (2)
cc                array(i+1,j,k) = array(i-1,j,k)
cc                rhs(j,k) = array(i-1,j,k) - dh(1)*div(i,j,k,array,by,bz)
cc              case (3)
cc                array(i,j-1,k) = array(i,j+1,k)
cc                rhs(i,k) = array(i,j+1,k) + dh(2)*div(i,j,k,bx,array,bz)
cc              case (4)
cc                array(i,j+1,k) = array(i,j-1,k)
cc                rhs(i,k) = array(i,j-1,k) - dh(2)*div(i,j,k,bx,array,bz)
cc              case (5)
cc                array(i,j,k-1) = array(i,j,k+1)
cc                rhs(i,j) = array(i,j,k+1) + dh(3)*div(i,j,k,bx,by,array)
cc              case (6)
cc                array(i,j,k+1) = array(i,j,k-1)
cc                rhs(i,j) = array(i,j,k-1) - dh(3)*div(i,j,k,bx,by,array)
cc              end select
cc
cc              enddo
cc            enddo
cc          enddo
cc
cc        endif
cc
cc      case (-IBX,-IBY,-IBZ) !Finds current normal components at boundaries
cc
cc        if (ieq == -IBX) icomp = 1
cc        if (ieq == -IBY) icomp = 2
cc        if (ieq == -IBZ) icomp = 3
cc
cc        if (icomp == dim) then !Normal components (div(J)=0)
cc
cc          do i=imin,imax
cc            do j=jmin,jmax
cc              do k=kmin,kmax
cc
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              dh(1) = 2.*dxh(ig)
cc              dh(2) = 2.*dyh(jg)
cc              dh(3) = 2.*dzh(kg)
cc
cc              select case (ibc)
cc              case (1)
cc                array(i-1,j,k) = array(i+1,j,k)
cc                rhs(j,k) = array(i+1,j,k) + dh(1)*div(i,j,k,array,jy,jz)
cc              case (2)                                                  
cc                array(i+1,j,k) = array(i-1,j,k)                         
cc                rhs(j,k) = array(i-1,j,k) - dh(1)*div(i,j,k,array,jy,jz)
cc              case (3)                                                  
cc                array(i,j-1,k) = array(i,j+1,k)                         
cc                rhs(i,k) = array(i,j+1,k) + dh(2)*div(i,j,k,jx,array,jz)
cc              case (4)                                                  
cc                array(i,j+1,k) = array(i,j-1,k)                         
cc                rhs(i,k) = array(i,j-1,k) - dh(2)*div(i,j,k,jx,array,jz)
cc              case (5)                                                  
cc                array(i,j,k-1) = array(i,j,k+1)                         
cc                rhs(i,j) = array(i,j,k+1) + dh(3)*div(i,j,k,jx,jy,array)
cc              case (6)                                                  
cc                array(i,j,k+1) = array(i,j,k-1)                         
cc                rhs(i,j) = array(i,j,k-1) - dh(3)*div(i,j,k,jx,jy,array)
cc              end select
cc
cc              enddo
cc            enddo
cc          enddo
cc
cc        endif
cc
cc      case default
cc
cc        write (*,*) 'Error in dirichletBC'
cc        stop
cc
cc      end select
cc
ccc     End program
cc
cc      end subroutine dirichletBC

      end subroutine FillGhostNodes

c imposeBConFluxes
c####################################################################
      subroutine imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                            ,flxkp,flxkm,bconds)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on fluxes
c--------------------------------------------------------------------

      use BCS

      implicit none

c Call variables

      integer(4) :: i,j,k,bconds(6)
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

c Local variables

c Begin program

      if (i == 1 .and. bconds(1) == SP) flxim = 0d0

c End

      end subroutine imposeBConfluxes
