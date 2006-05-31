c module BCS_variables
c #####################################################################
#include "../common/BC/BCS_variables_mod.f"
c end module BCS_variables

c module local_BCS_variables
c####################################################################
      module local_BCS_variables

        use BCS_variables

        use auxiliaryVariables

        use operators

        use icond

        use equilibrium

        use timestepping

      end module local_BCS_variables

c module dirichletBCinterface
c####################################################################
      module dirichletBCinterface

        use local_BCS_variables

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
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

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
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,*)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1,*)

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

      case (IBX,IBY,IBZ,IJX,IJY,IJZ,IAX,IAY,IAZ) !Imposes divergence-free constraint on B and J

        if (ivar /= dim) then

          call interpolate(array (:,:,:,ivar),array0(:,:,:,ivar)
     .                    ,ibc,order)

        else

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = 2.*dxh(ig)
                dh(2) = 2.*dyh(jg)
                dh(3) = 2.*dzh(kg)

                select case (ibc)
                case (1)
                  array(i-1,j,k,dim) = array(i+1,j,k,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(j,k) = array(i+1,j,k,dim) + dh(dim)*diver
                case (2)
                  array(i+1,j,k,dim) = array(i-1,j,k,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(j,k) = array(i-1,j,k,dim) - dh(dim)*diver
                case (3)
                  array(i,j-1,k,dim) = array(i,j+1,k,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(i,k) = array(i,j+1,k,dim) + dh(dim)*diver
                case (4)
                  array(i,j+1,k,dim) = array(i,j-1,k,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(i,k) = array(i,j-1,k,dim) - dh(dim)*diver
                case (5)
                  array(i,j,k-1,dim) = array(i,j,k+1,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(i,j) = array(i,j,k+1,dim) + dh(dim)*diver
                case (6)
                  array(i,j,k+1,dim) = array(i,j,k-1,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(i,j) = array(i,j,k-1,dim) - dh(dim)*diver
                end select

              enddo
            enddo
          enddo

c       Update BC ghost nodes

          select case (ibc)
          case (1)              !x0
          array(0     ,jmin:jmax,kmin:kmax,dim)=rhs(jmin:jmax,kmin:kmax)
          case (2)              !x1
          array(nxbc+1,jmin:jmax,kmin:kmax,dim)=rhs(jmin:jmax,kmin:kmax)
          case (3)              !y0
          array(imin:imax,0     ,kmin:kmax,dim)=rhs(imin:imax,kmin:kmax)
          case (4)              !y1
          array(imin:imax,nybc+1,kmin:kmax,dim)=rhs(imin:imax,kmin:kmax)
          case (5)              !z0
          array(imin:imax,jmin:jmax,0     ,dim)=rhs(imin:imax,jmin:jmax)
          case (6)              !z1
          array(imin:imax,jmin:jmax,nzbc+1,dim)=rhs(imin:imax,jmin:jmax)
          case default
            messg = 'Boundary'//int2char(ibc)//' non existent'
            call pstop('vectorDirichletBC',messg)
          end select

        endif

      case default

        messg = 'Equation'//int2char(ieq)//' does not exist'
        call pstop('vectorDirichletBC',messg)

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
        real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .               ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

        integer(4) :: i,j,k,ig,jg,kg

c     Begin program

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

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

c       Update BC ghost nodes

        select case (ibc)
        case (1)                !x0
         array(0     ,jmin:jmax,kmin:kmax) = rhs(jmin:jmax,kmin:kmax)
        case (2)               !x1
         array(nxbc+1,jmin:jmax,kmin:kmax) = rhs(jmin:jmax,kmin:kmax)
        case (3)                !y0
         array(imin:imax,0     ,kmin:kmax) = rhs(imin:imax,kmin:kmax)
        case (4)                !y1
         array(imin:imax,nybc+1,kmin:kmax) = rhs(imin:imax,kmin:kmax)
        case (5)                !z0
         array(imin:imax,jmin:jmax,0     ) = rhs(imin:imax,jmin:jmax)
        case (6)                !z1
         array(imin:imax,jmin:jmax,nzbc+1) = rhs(imin:imax,jmin:jmax)
        case default
         messg = 'Boundary'//int2char(ibc)//' non existent'
         call pstop('interpolate',messg)
        end select

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
        messg = 'Order of interpolation not implemented'
        call pstop('quad_int',messg)
      end select

c     End program

      end function quad_int

      end module dirichletBCinterface

c module neumannBCinterface
c #####################################################################
      module neumannBCinterface

       use local_BCS_variables

       use transport_params

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
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

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

              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

              gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)

              ip = min(i+1,nxbc)
              im = max(i-1,1)
              jp = min(j+1,nybc)
              jm = max(j-1,1)
              kp = min(k+1,nzbc)
              km = max(k-1,1)

              dh(1) = 2.*dxh(ig)
              if (i == nxbc) dh(1) = dx(ig-1)
              if (i == 1   ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == nybc) dh(2) = dy(jg-1)
              if (j == 1   ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nzbc) dh(3) = dz(kg-1)
              if (k == 1   ) dh(3) = dz(kg)

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

              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

              gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)
              jac0   = gmetric%grid(igxbc)%jac (i,j,k)

              hessian = gmetric%grid(igxbc)%Gamma(i,j,k,:,:,:)

              ip = min(i+1,nxbc)
              im = max(i-1,1   )
              jp = min(j+1,nybc)
              jm = max(j-1,1   )
              kp = min(k+1,nzbc)
              km = max(k-1,1   )

              dh(1) = 2.*dxh(ig)
              if (i == nxbc) dh(1) = dx(ig-1)
              if (i == 1   ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == nybc) dh(2) = dy(jg-1)
              if (j == 1   ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nzbc) dh(3) = dz(kg-1)
              if (k == 1   ) dh(3) = dz(kg)

              if (dim == 1) then
                if (gamma > 1d0) then
                  rhs(j,k) =  hessian(1,1,1)*vx(i,j,k)*vx(i,j,k)
     .                       +hessian(1,2,2)*vy(i,j,k)*vy(i,j,k)
     .                       +hessian(1,3,3)*vz(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(1,1,2)*vx(i,j,k)*vy(i,j,k)
     .                    +2.*hessian(1,1,3)*vx(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(1,2,3)*vy(i,j,k)*vz(i,j,k)
                else
                  rhs(j,k) = 0d0
                endif
                rhs(j,k) = -dh(dim)
     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
     .             -0.5/jac0*rhs(j,k))/gsuper(dim,dim)
              elseif (dim == 2) then
                if (gamma > 1d0) then
                  rhs(i,k) =  hessian(2,1,1)*vx(i,j,k)*vx(i,j,k)
     .                       +hessian(2,2,2)*vy(i,j,k)*vy(i,j,k)
     .                       +hessian(2,3,3)*vz(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(2,1,2)*vx(i,j,k)*vy(i,j,k)
     .                    +2.*hessian(2,1,3)*vx(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(2,2,3)*vy(i,j,k)*vz(i,j,k)
                else
                  rhs(i,k) = 0d0
                endif
                rhs(i,k) = -dh(dim)
     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
     .             -0.5/jac0*rhs(i,k))/gsuper(dim,dim)
              elseif (dim == 3) then
                if (gamma > 1d0) then
                  rhs(i,j) =  hessian(3,1,1)*vx(i,j,k)*vx(i,j,k)
     .                       +hessian(3,2,2)*vy(i,j,k)*vy(i,j,k)
     .                       +hessian(3,3,3)*vz(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(3,1,2)*vx(i,j,k)*vy(i,j,k)
     .                    +2.*hessian(3,1,3)*vx(i,j,k)*vz(i,j,k)
     .                    +2.*hessian(3,2,3)*vy(i,j,k)*vz(i,j,k)
                else
                  rhs(i,j) = 0d0
                endif
                rhs(i,j) = -dh(dim)
     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
     .              +gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
     .             -0.5/jac0*rhs(i,j))/gsuper(dim,dim)
              endif

            enddo
          enddo
        enddo

      case default

        call pstop('scalarNeumannBC','Variable not implemented')

      end select

c     Assign value

      select case (ibc)
      case (1)
        array(0     ,jmin:jmax,kmin:kmax)
     $                           = array(1,jmin:jmax,kmin:kmax)
     $                           - rhs(jmin:jmax,kmin:kmax)
      case (2)
        array(nxbc+1,jmin:jmax,kmin:kmax)
     $                           = array(nxbc,jmin:jmax,kmin:kmax)
     $                           + rhs(jmin:jmax,kmin:kmax)
      case (3)
        array(imin:imax,0     ,kmin:kmax)
     $                           = array(imin:imax,1,kmin:kmax)
     $                           - rhs(imin:imax,kmin:kmax)
      case (4)
        array(imin:imax,nybc+1,kmin:kmax)
     $                           = array(imin:imax,nybc,kmin:kmax)
     $                           + rhs(imin:imax,kmin:kmax)
      case (5)
        array(imin:imax,jmin:jmax,0     )
     $                           = array(imin:imax,jmin:jmax,1)
     $                           - rhs(imin:imax,jmin:jmax)
      case (6)
        array(imin:imax,jmin:jmax,nzbc+1)
     $                           = array(imin:imax,jmin:jmax,nzbc)
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
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,*)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,ip,im,jp,jm,kp,km,ibc,icomp,isgn
      real(8)    :: x1,x2,x3,dh(3),jac0,jxx,jyy,jzz,ee,alp,vec(3),cnv(3)
     .             ,cov(3),J0(3),bnorm
      real(8)    :: gsuper(3,3),hessian(3,3,3)
     .             ,gsupip(3,3),gsupim(3,3)
     .             ,gsupjp(3,3),gsupjm(3,3)
     .             ,gsupkp(3,3),gsupkm(3,3)

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      rhs = 0d0

      select case (ieq)
cc      case (IVX,IVY,IVZ,IBX,IBY,IBZ) !Velocity components
      case (IVX,IVY,IVZ,IJX,IJY,IJZ) !Velocity components
cc      case (IVX,IVY,IVZ) !Velocity components

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

              gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)

              hessian = gmetric%grid(igxbc)%Gamma(i,j,k,:,:,:)
  
              ip = min(i+1,nxbc)
              im = max(i-1,1   )
              jp = min(j+1,nybc)
              jm = max(j-1,1   )
              kp = min(k+1,nzbc)
              km = max(k-1,1   )

              dh(1) = 2.*dxh(ig)
              if (i == nxbc) dh(1) = dx(ig-1)
              if (i == 1   ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == nybc) dh(2) = dy(jg-1)
              if (j == 1   ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nzbc) dh(3) = dz(kg-1)
              if (k == 1   ) dh(3) = dz(kg)

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

              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

              gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)
              gsupip = gmetric%grid(igxbc)%gsup(i+1,j,k,:,:)
              gsupim = gmetric%grid(igxbc)%gsup(i-1,j,k,:,:)
              gsupjp = gmetric%grid(igxbc)%gsup(i,j+1,k,:,:)
              gsupjm = gmetric%grid(igxbc)%gsup(i,j-1,k,:,:)
              gsupkp = gmetric%grid(igxbc)%gsup(i,j,k+1,:,:)
              gsupkm = gmetric%grid(igxbc)%gsup(i,j,k-1,:,:)

              ip = min(i+1,nxbc)
              im = max(i-1,1   )
              jp = min(j+1,nybc)
              jm = max(j-1,1   )
              kp = min(k+1,nzbc)
              km = max(k-1,1   )

              dh(1) = 2.*dxh(ig)
              if (i == nxbc) dh(1) = dx(ig-1)
              if (i == 1   ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == nybc) dh(2) = dy(jg-1)
              if (j == 1   ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nzbc) dh(3) = dz(kg-1)
              if (k == 1   ) dh(3) = dz(kg)

              if (eeta(i,j,k) > 0d0) then
                J0 = E0/eeta(i,j,k)
              else
                J0 = 0d0
              endif

              if (equil == 'ppnch' .or. equil == 'ppn3d') then

                !THIS WON'T WORK FOR LINEARIZED EQUATIONS; SHOULD USE ACTUAL B
                !   BUT THIS DOESN'T WORK YET.
                cov = (/array(i,j,k,1),array(i,j,k,2),array(i,j,k,3)/)
cc                cov = (/varray%array(IBX)%array(i,j,k)
cc     .                 ,varray%array(IBY)%array(i,j,k)
cc     .                 ,varray%array(IBZ)%array(i,j,k)/)

                cnv = gsuper(:,1)*cov(1)
     .               +gsuper(:,2)*cov(2)
     .               +gsuper(:,3)*cov(3)

                bnorm = scalarProduct(i,j,k,igxbc,igybc,igzbc
     .                               ,cov(1),cov(2),cov(3)
     .                               ,cnv(1),cnv(2),cnv(3))

                if (dim == 1) then

                  jxx = (array(i,jp,k,3)-array(i,jm,k,3))/dh(2)
     .                 -(array(i,j,kp,2)-array(i,j,km,2))/dh(3)

                  if (ivar == 2) then
                    icomp = 3
                    isgn  = 1

                    rhs(j,k) = dh(dim)
     .                *( (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                    +gsuper(dim,3)*jxx/gsuper(dim,dim))

                  elseif (ivar == 3) then
                    icomp = 2
                    isgn  =-1

                    rhs(j,k) = dh(dim)
     .                *( (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                    -gsuper(dim,2)*jxx/gsuper(dim,dim))
                  endif

                  if (eta > 0d0 .and. bnorm > 0d0) then
                    alp = scalarProduct(i,j,k,igxbc,igybc,igzbc
     .                                 ,J0 (1),J0 (2),J0 (3)
     .                                 ,cnv(1),cnv(2),cnv(3))
     .                   /bnorm
                  else
                    alp = 0d0
                  endif

                  ee = alp*cnv(icomp)

                  rhs(j,k)=(rhs(j,k) + dh(dim)*isgn*ee)

                elseif (dim == 2) then

                  jyy = (array(i,j,kp,1)-array(i,j,km,1))/dh(3)
     .                 -(array(ip,j,k,3)-array(im,j,k,3))/dh(1)

                  if (ivar == 3) then
                    icomp = 1
                    isgn  = 1

                    rhs(i,k) = dh(dim)
     .                 *( (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                    +gsuper(dim,1)*jyy/gsuper(dim,dim))
                  elseif (ivar == 1) then
                    icomp = 3
                    isgn  =-1

                    rhs(i,k) = dh(dim)
     .                 *( (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                    -gsuper(dim,3)*jyy/gsuper(dim,dim))
                  endif

                  if (eta > 0d0 .and. bnorm > 0d0) then
                    alp = scalarProduct(i,j,k,igxbc,igybc,igzbc
     .                                 ,J0 (1),J0 (2),J0 (3)
     .                                 ,cnv(1),cnv(2),cnv(3))
     .                   /bnorm
                  else
                    alp = 0d0
                  endif

                  ee = alp*cnv(icomp)

                  rhs(i,k) = (rhs(i,k) + dh(dim)*isgn*ee)

                elseif (dim == 3) then

                  jzz = (array(ip,j,k,2)-array(im,j,k,2))/dh(1)
     .                 -(array(i,jp,k,1)-array(i,jm,k,1))/dh(2)

                  if (ivar == 1) then
                    icomp = 2
                    isgn  = 1

                    rhs(i,j) = dh(dim)
     .                 *( (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                    +gsuper(dim,2)*jzz/gsuper(dim,dim))
                  elseif (ivar == 2) then
                    icomp = 1
                    isgn  =-1

                    rhs(i,j) = dh(dim)
     .                 *( (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                    -gsuper(dim,1)*jzz/gsuper(dim,dim))
                  endif

                  if (eta > 0d0 .and. bnorm > 0d0) then
                    alp = scalarProduct(i,j,k,igxbc,igybc,igzbc
     .                                 ,J0 (1),J0 (2),J0 (3)
     .                                 ,cnv(1),cnv(2),cnv(3))
     .                   /bnorm
                  else
                    alp = 0d0
                  endif

                  ee = alp*cnv(icomp)

                  rhs(i,j) = (rhs(i,j) + dh(dim)*isgn*ee)

                endif

              else

                if (dim == 1) then

                  jxx = (array(i,jp,k,3)-array(i,jm,k,3))/dh(2)
     .                 -(array(i,j,kp,2)-array(i,j,km,2))/dh(3)

                  if (ivar == 2) then
                    icomp = 3
                    isgn  = 1

                    rhs(j,k) = dh(dim)
     .                 *( (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                    +gsuper(dim,3)*jxx/gsuper(dim,dim))
                  elseif (ivar == 3) then
                    icomp = 2
                    isgn  =-1

                    rhs(j,k) = dh(dim)
     .                 *( (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                    -gsuper(dim,2)*jxx/gsuper(dim,dim))
                  endif

                  if (ibc == 1) then
                    vec = gsuper(icomp,:)
     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
                    vec = 0.5*(vec + gsupim(icomp,:)
     .                 -gsupim(icomp,dim)/gsupim(dim,dim)*gsupim(dim,:))
                    ee = dot_product(vec,J0)

                    rhs(j,k) = rhs(j,k) + dh(dim)*isgn*ee
                  elseif (ibc == 2) then
                    vec = gsuper(icomp,:)
     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
                    vec = 0.5*(vec + gsupip(icomp,:)
     .                 -gsupip(icomp,dim)/gsupip(dim,dim)*gsupip(dim,:))
                    ee = dot_product(vec,J0)

                    rhs(j,k) = rhs(j,k) + dh(dim)*isgn*ee
                  endif

                elseif (dim == 2) then

                  jyy = (array(i,j,kp,1)-array(i,j,km,1))/dh(3)
     .                 -(array(ip,j,k,3)-array(im,j,k,3))/dh(1)

                  if (ivar == 3) then
                    icomp = 1
                    isgn  = 1

                    rhs(i,k) = dh(dim)
     .                 *( (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                    +gsuper(dim,1)*jyy/gsuper(dim,dim))
                  elseif (ivar == 1) then
                    icomp = 3
                    isgn  =-1

                    rhs(i,k) = dh(dim)
     .                 *( (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                    -gsuper(dim,3)*jyy/gsuper(dim,dim))
                  endif

                  if (ibc == 3) then
                    vec = gsuper(icomp,:)
     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
                    vec = 0.5*(vec + gsupjm(icomp,:)
     .                 -gsupjm(icomp,dim)/gsupjm(dim,dim)*gsupjm(dim,:))
                    ee = dot_product(vec,J0)

                    rhs(i,k) = rhs(i,k) + dh(dim)*isgn*ee
                  elseif (ibc == 4) then
                    vec = gsuper(icomp,:)
     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
                    vec = 0.5*(vec + gsupjp(icomp,:)
     .                 -gsupjp(icomp,dim)/gsupjp(dim,dim)*gsupjp(dim,:))
                    ee = dot_product(vec,J0)

                    rhs(i,k) = rhs(i,k) + dh(dim)*isgn*ee
                  endif

                elseif (dim == 3) then

                  jzz = (array(ip,j,k,2)-array(im,j,k,2))/dh(1)
     .                 -(array(i,jp,k,1)-array(i,jm,k,1))/dh(2)

                  if (ivar == 1) then
                    icomp = 2
                    isgn  = 1

                    rhs(i,j) = dh(dim)
     .                 *( (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                    +gsuper(dim,2)*jzz/gsuper(dim,dim))
                  elseif (ivar == 2) then
                    icomp = 1
                    isgn  =-1

                    rhs(i,j) = dh(dim)
     .                 *( (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                    -gsuper(dim,1)*jzz/gsuper(dim,dim))
                  endif

                  if (ibc == 5) then
                    vec = gsuper(icomp,:)
     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
                    vec = 0.5*(vec + gsupkm(icomp,:)
     .                 -gsupkm(icomp,dim)/gsupkm(dim,dim)*gsupkm(dim,:))
                    ee = dot_product(vec,J0)

                    rhs(i,j) = rhs(i,j) + dh(dim)*isgn*ee
                  elseif (ibc == 6) then
                    vec = gsuper(icomp,:)
     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
                    vec = 0.5*(vec + gsupkp(icomp,:)
     .                 -gsupkp(icomp,dim)/gsupkp(dim,dim)*gsupkp(dim,:))
                    ee = dot_product(vec,J0)

                    rhs(i,j) = rhs(i,j) + dh(dim)*isgn*ee
                  endif

                endif

              endif

            enddo
          enddo
        enddo

      case default

        call pstop('vectorNeumanBC','Variable not implemented')

      end select

c     Assign value

      select case (ibc)
      case (1)
        array(0     ,jmin:jmax,kmin:kmax,ivar)
     $                           = array(1,jmin:jmax,kmin:kmax,ivar)
     $                           - rhs(jmin:jmax,kmin:kmax)
      case (2)
        array(nxbc+1,jmin:jmax,kmin:kmax,ivar)
     $                           = array(nxbc,jmin:jmax,kmin:kmax,ivar)
     $                           + rhs(jmin:jmax,kmin:kmax)
      case (3)
        array(imin:imax,0     ,kmin:kmax,ivar)
     $                           = array(imin:imax,1,kmin:kmax,ivar)
     $                           - rhs(imin:imax,kmin:kmax)
      case (4)
        array(imin:imax,nybc+1,kmin:kmax,ivar)
     $                           = array(imin:imax,nybc,kmin:kmax,ivar)
     $                           + rhs(imin:imax,kmin:kmax)
      case (5)
        array(imin:imax,jmin:jmax,0     ,ivar)
     $                           = array(imin:imax,jmin:jmax,1,ivar)
     $                           - rhs(imin:imax,jmin:jmax)
      case (6)
        array(imin:imax,jmin:jmax,nzbc+1,ivar)
     $                           = array(imin:imax,jmin:jmax,nzbc,ivar)
     $                           + rhs(imin:imax,jmin:jmax)
      end select

c     End program

      end subroutine vectorNeumannBC

      end module neumannBCinterface

c module imposeBCinterface
c #####################################################################
cc#include "../common/BC/imposeBCinterface_mod.F"
c end module imposeBCinterface

c module singularBCinterface
c####################################################################
      module singularBCinterface

        use BCS_variables

      contains

c     scalarSingularBC2
c     #################################################################
      subroutine scalarSingularBC2(array,order)
c     -----------------------------------------------------------------
c     Imposes singular point BC. On input:
c        * array: contains variable on which singular BC is imposed
c        * order: order of interpolation towards singular point
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: order
      real(8)    :: array(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,order1
      real(8)    :: lavg_q,lavg_vol,avg_q,avg_vol,rho0(1),vol,x0(1)

c     Begin program

c     Check that we are domain adjacent to singular point

      if (ilog /= 1) then
        return
      else
        i = ilog
      endif

c     Set interpolation order

      order1 = min(order,ihig-ilog)

c     Find scalar value at singular point by averaging neighboring cells

      do k=1,nzbc

c       Perform local sum
        lavg_q   = 0d0
        lavg_vol = 0d0

        do j=1,nybc
          call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
          vol = volume(i,j,k,igxbc,igybc,igzbc)
          x0  = grid_params%xx(ig-1)

          call IntDriver1d(order1+1,grid_params%xx(ig:ig+order1)
     .                    ,array(i:i+order1+1,j,k)
     .                    ,1,x0,rho0,order1)

          lavg_q   = lavg_q   + vol*rho0(1)
          lavg_vol = lavg_vol + vol
        enddo

c       Perform global sum (and broadcast it to all processors)

#if defined(petsc)
        call MPI_Allreduce(lavg_q  ,avg_q  ,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
        call MPI_Allreduce(lavg_vol,avg_vol,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
#else
        avg_q   = lavg_q
        avg_vol = lavg_vol
#endif

c       Perform average

        array(0,:,k) = avg_q/avg_vol

      enddo

c     End program

      end subroutine scalarSingularBC2

c     scalarSingularBC
c     #################################################################
      subroutine scalarSingularBC(array)
c     -----------------------------------------------------------------
c     Imposes singular point BC. On input:
c        * array: contains variable on which singular BC is imposed
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: array(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,jj

c     Begin program

c     Check that we are domain adjacent to singular point

      if (ilog /= 1) then
        return
      else
        i = ilog
      endif

c     Element exchange in poloidal direction

      if (mod(nybc,2) /= 0) then
        call pstop('scalarSingularBC'
     .            ,'Poloidal dimension not a multiple of 2')
      endif

      do k=1,nzbc
        do j=1,nybc
          jj = mod(j+nybc/2,nybc)
          if (jj == 0) jj = nybc
          array(i-1,j,k) = array(i,jj,k)
        enddo
      enddo

c     End program

      end subroutine scalarSingularBC

c     vectorSingularBC2
c     #################################################################
      subroutine vectorSingularBC2(vec,cov,order)
c     -----------------------------------------------------------------
c     Averages vector components around singular point and calculates
c     curvilinear components at singular point.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: order
      real(8)    :: vec(0:nxbc+1,0:nybc+1,0:nzbc+1,3)
      logical    :: cov

c     Local variables

      integer(4) :: i,j,k

      integer(4) :: ic,ig,jg,kg,order1
      real(8)    :: vol,lavg_x,lavg_y,lavg_z,lavg_vol
     .                 , avg_x, avg_y, avg_z, avg_vol
      real(8)    :: daxdr(1),daydr(1),dummy1,dummy2
      real(8),allocatable,dimension(:) :: ax0,ay0,az0,cx,cy,cz,xx

c     External

c     Begin program

c     Check that we are domain adjacent to singular point

      if (ilog /= 1) then
        return
      else
        i = ilog
      endif

c     Set interpolation order

      order1 = min(order,ihig-ilog)

c     Allocate quantities for interpolation

      allocate(ax0(nzbc),ay0(nzbc),az0(nzbc))
      allocate(cx(0:order1+1),cy(0:order1+1),cz(0:order1+1)
     .        ,xx(0:order1+1))

      call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

      xx(0:order1+1) = grid_params%xx(ig-1:ig+order1)

c     Find average cartesian coordinates

      do k=1,nzbc

c       Perform local sum
        lavg_x   = 0d0
        lavg_y   = 0d0
        lavg_z   = 0d0
        lavg_vol = 0d0

        do j=1,nybc
          call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
          vol = volume(i,j,k,igxbc,igybc,igzbc)

          do ic=1,order1+1
            call transformVectorToCartesian(ic,j,k,igxbc,igybc,igzbc
     .           ,vec(ic,j,k,1),vec(ic,j,k,2),vec(ic,j,k,3),cov
     .           ,cx(ic),cy(ic),cz(ic))
          enddo

          call IntDriver1d(order1+1,grid_params%xx(ig:ig+order1)
     .                    ,cx(1:order1+1),1,grid_params%xx(ig-1:ig-1)
     .                    ,cx(0),order1)
          call IntDriver1d(order1+1,grid_params%xx(ig:ig+order1)
     .                    ,cy(1:order1+1),1,grid_params%xx(ig-1:ig-1)
     .                    ,cy(0),order1)
          call IntDriver1d(order1+1,grid_params%xx(ig:ig+order1)
     .                    ,cz(1:order1+1),1,grid_params%xx(ig-1:ig-1)
     .                    ,cz(0),order1)

          lavg_x   = lavg_x   + vol*cx(0)
          lavg_y   = lavg_y   + vol*cy(0)
          lavg_z   = lavg_z   + vol*cz(0)
          lavg_vol = lavg_vol + vol

        enddo

c       Perform global sum (and broadcast it to all processors)

#if defined(petsc)
        call MPI_Allreduce(lavg_x  ,avg_x  ,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
        call MPI_Allreduce(lavg_y  ,avg_y  ,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
        call MPI_Allreduce(lavg_z  ,avg_z  ,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
        call MPI_Allreduce(lavg_vol,avg_vol,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
#else
        avg_x   = lavg_x
        avg_y   = lavg_y
        avg_z   = lavg_z
        avg_vol = lavg_vol
#endif

c       Perform average

        ax0(k) = avg_x/avg_vol
        ay0(k) = avg_y/avg_vol
        az0(k) = avg_z/avg_vol

      enddo

      i = ilog-1
      do k=1,nzbc
        do j=0,nybc+1           !This also fills ghost cells
          call transformVectorToCurvilinear(i,j,k,igxbc,igybc,igzbc
     .                 ,ax0(k),ay0(k),az0(k),cov
     .                 ,vec(i,j,k,1),vec(i,j,k,2),vec(i,j,k,3))
        enddo
      enddo

      deallocate(ax0,ay0,az0)
      deallocate(cx,cy,cz,xx)

c     End program

      end subroutine vectorSingularBC2

c     vectorSingularBC
c     #################################################################
      subroutine vectorSingularBC(ivar,vec,cov)
c     -----------------------------------------------------------------
c     Averages vector components around singular point and calculates
c     curvilinear components at singular point.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ivar
      real(8)    :: vec(0:nxbc+1,0:nybc+1,0:nzbc+1)
      logical    :: cov

c     Local variables

      integer(4) :: i,j,k,jj

c     External

c     Begin program

c     Check that we are domain adjacent to singular point

      if (ilog /= 1) then
        return
      else
        i = ilog
      endif

c     Element exchange in poloidal direction

      if (mod(nybc,2) /= 0) then
        call pstop('vectorSingularBC'
     .            ,'Poloidal dimension not a multiple of 2')
      endif

      if (coords == 'hl2') then
        do k=1,nzbc
          do j=1,nybc
            jj = mod(j+nybc/2,nybc)
            if (jj == 0) jj = nybc
            vec(i-1,j,k) =  vec(i,jj,k)
          enddo
        enddo
      else
        do k=1,nzbc
          do j=1,nybc
            jj = mod(j+nybc/2,nybc)
            if (jj == 0) jj = nybc
            if (cov) then
              select case(ivar)
              case(1)
                vec(i-1,j,k) = -vec(i,jj,k)
              case(2,3)
                vec(i-1,j,k) =  vec(i,jj,k)
              end select
            else
              select case(ivar)
              case(1)
                vec(i-1,j,k) =  vec(i,jj,k)
              case(2,3)
                vec(i-1,j,k) = -vec(i,jj,k)
              end select
            endif
          enddo
        enddo
      endif

c     End program

      end subroutine vectorSingularBC

      end module singularBCinterface

c module imposeBCinterface
c #####################################################################
      module imposeBCinterface

        use BCS_variables

        use singularBCinterface

        use dirichletBCinterface

        use neumannBCinterface

        INTERFACE setBC
          module procedure imposeBConScalar,imposeBConVector
        end INTERFACE

        integer(4) :: iimax,iimin,jjmax,jjmin,kkmax,kkmin

        integer(4),private :: bctype,dim,loc,ibc,i,j,k

        !Define order of interpolation for EDGE and CORNER treatment
        !(linear to preserve nearest-neighbors stencil; this is 
        ! important for colored diagonal formation)
        integer(4),private :: order=1

      contains

c     findLoopLimits
c     ##############################################################
      subroutine findLoopLimits(dim,loc,i1 ,i2 ,j1 ,j2 ,k1 ,k2
     .                                 ,imn,imx,jmn,jmx,kmn,kmx)
c     --------------------------------------------------------------
c     Finds local limits for BC loops. On input, it takes global 
c     limits (i1,i2), (j1,j2), (k1,k2) and outputs local loop limits
c     (imn,imx), (jmn,jmx), (kmn,kmx).
c     --------------------------------------------------------------

c     Call variables

      integer(4),intent(IN) :: dim,loc,i1,i2,j1,j2,k1,k2
      integer(4),intent(OUT):: imn,imx,jmn,jmx,kmn,kmx

c     Begin program

      select case(dim)
      case (1)

        call selectLimits(i1,i2,grid_params%ilo(igxbc)
     .                   ,nxbc,grid_params%nxgl(igxbc),imn,imx)

        jmn = j1 - grid_params%jlo(igybc) + 1
        jmx = j2 - grid_params%jlo(igybc) + 1

        kmn = k1 - grid_params%klo(igzbc) + 1
        kmx = k2 - grid_params%klo(igzbc) + 1

      case(2)
        imn = i1 - grid_params%ilo(igxbc) + 1
        imx = i2 - grid_params%ilo(igxbc) + 1

        call selectLimits(j1,j2,grid_params%jlo(igybc)
     .                   ,nybc,grid_params%nygl(igybc),jmn,jmx)

        kmn = k1 - grid_params%klo(igzbc) + 1
        kmx = k2 - grid_params%klo(igzbc) + 1

      case(3)
        imn = i1 - grid_params%ilo(igxbc) + 1
        imx = i2 - grid_params%ilo(igxbc) + 1

        jmn = j1 - grid_params%jlo(igybc) + 1
        jmx = j2 - grid_params%jlo(igybc) + 1

        call selectLimits(k1,k2,grid_params%klo(igzbc)
     .                   ,nzbc,grid_params%nzgl(igzbc),kmn,kmx)

      end select

c     End program

      contains

c     selectLimits
c     #####################################################################
      subroutine selectLimits(igl,igr,igmin,nxl,nxg,ill,ilr)

c     ---------------------------------------------------------------------
c     Selects LOCAL limits for BC treatment at boundary location determined
c     by (dim,loc).
c     ---------------------------------------------------------------------

        implicit none

        integer(4) :: igl,igr,igmin,nxl,nxg,ill,ilr

c     Begin program

cc        do_nothing = .false.

        select case(loc)
        case(0)
cc          if (igl > 1) do_nothing = .true.
          ill = min(igl - igmin + 1,1)
          ilr = ill
        case(1)
cc          if (igr < nxg) do_nothing = .true.
          ill = max(igr - igmin + 1,nxl)
          ilr = ill
        end select

      end subroutine selectLimits

      end subroutine findLoopLimits

c     do__nothing
c     #####################################################################
      function do__nothing(dim,loc,i1,i2,j1,j2,k1,k2)

c     ---------------------------------------------------------------------
c     Selects whether to fill or not ghost cells of a given local domain
c     based on its position in the global grid.
c     ---------------------------------------------------------------------

        implicit none

        logical    :: do__nothing
        integer(4),intent(IN) :: dim,loc,i1,i2,j1,j2,k1,k2

        integer(4) :: ibc

c     Begin program

        do__nothing = .false.
        if (asm) return      !Return if additive Schwartz method is used (PC)

        select case(dim)
        case(1)
          call checkLimits(i1,i2,grid_params%nxgl(igxbc))
        case(2)
          call checkLimits(j1,j2,grid_params%nygl(igybc))
        case(3)
          call checkLimits(k1,k2,grid_params%nzgl(igzbc))
        end select

        !Do nothing for periodic BCs in parallel (PETSc does it)
        ibc = (1+loc)+2*(dim-1)
cc        if ((.not.do__nothing).and.(bcond(ibc) == PER.and.asm_dir(dim)))
        if((.not.do__nothing).and.(bcond(ibc)==PER).and.(np>1))
     .       do__nothing =.true. 

c     End program

      contains

c     checkLimits
c     #####################################################################
      subroutine checkLimits(igl,igr,nxg)

c     ---------------------------------------------------------------------
c     Selects LOCAL limits for BC treatment at boundary location determined
c     by (dim,loc).
c     ---------------------------------------------------------------------

        implicit none

        integer(4) :: igl,igr,nxg

c     Begin program

        select case(loc)
        case(0)
          if (igl > 1  ) do__nothing = .true.
        case(1)
          if (igr < nxg) do__nothing = .true.
        end select

      end subroutine checkLimits

      end function do__nothing

c     selectFace
c     #####################################################################
      function selectFace(ibc)

c     ---------------------------------------------------------------------
c     Selects LOCAL limits for BC treatment at boundary location determined
c     by (dim,loc).
c     ---------------------------------------------------------------------

      implicit none

      logical    :: selectFace

      integer(4) :: dim,loc,ibc

c     Begin program

      loc = abs(mod(ibc,2)-2)-1
      dim = 1 + (ibc-1-loc)/2

cc#if defined(petsc)
cc      selectFace= .not.((bcond(ibc) == PER).or.
cc     .     do__nothing(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax))
cc#else
      selectFace= 
     .     .not.do__nothing(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax)
cc#endif

cc      selectFace = .true.

      end function selectFace

c     setup_petsc_BC
c     ######################################################################
      subroutine setup_petsc_BC

c     -----------------------------------------------------------------
c     Sets flags to indicate that parallel BCs are NOT to be imposed
c     on the fortran side (are fully provided by PETSc).
c     -----------------------------------------------------------------

      implicit none

      par_bc = .false.  !No parallel BC (provided by PETSc)
      asm    = .false.  !No ASM (no extrapolation to corner, edges)

      end subroutine setup_petsc_BC

c     setASMflag
c     ######################################################################
      subroutine setASMflag(asm_flag)

c     -----------------------------------------------------------------
c     Sets flag to indicate Additive Schwartz method.
c     -----------------------------------------------------------------

      implicit none

      logical :: asm_flag

      asm = asm_flag

      par_bc = .not.asm

      end subroutine setASMflag

c     setASM_PC_BC
c     ######################################################################
      subroutine setASM_PC_BC(igx,igy,igz,bcs,IIBC)

c     -----------------------------------------------------------------
c     Sets boundary conditions for Additive Schwartz method PC.
c     -----------------------------------------------------------------

      implicit none

      integer(4) :: igx,igy,igz,bcs(:,:)

      integer(4),optional :: IIBC

      integer(4) :: dim,loc,ibc,igmin,igmax,jgmin,jgmax,kgmin,kgmax

      logical    :: dn

      if (np == 1) then
        call setASMflag(.false.)
        return
      endif

c     Reset BC grid level

      igxbc = igx
      igybc = igy
      igzbc = igx

c     Find grid local limits

      igmin = grid_params%ilo(igx)
      igmax = grid_params%ihi(igx)
      jgmin = grid_params%jlo(igy)
      jgmax = grid_params%jhi(igy)
      kgmin = grid_params%klo(igz)
      kgmax = grid_params%khi(igz)

c     Determine ASM BCs (requires ASM flag set to false)

      call setASMflag(.false.)

      if (PRESENT(IIBC)) then
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            dn=do__nothing(dim,loc,igmin,igmax,jgmin,jgmax,kgmin,kgmax)

cc            if (dn) bcs(ibc,:) = IIBC
            if (dn) bcs(ibc,:) = IIBC*bcs(ibc,:)/abs(bcs(ibc,:))
          enddo
        enddo
cc      else
cc        do dim=1,3
cc          do loc=0,1
cc            ibc = (1+loc)+2*(dim-1)
cc            if (bcond(ibc) == PER.and.asm_dir(dim)) then !Take care of PER boundaries in ASM
cc              bcs(ibc,:) = EXT
cc            endif
cc          enddo
cc        enddo
      endif

c     Set ASM flag to true to perform ASM in PC

      call setASMflag(.true.)

      end subroutine setASM_PC_BC

c     imposeBConScalar
c     #################################################################
      subroutine imposeBConScalar(ieq,nx,ny,nz,array,array0,bcond
     .                           ,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                           ,iorder)
c     -----------------------------------------------------------------
c     Imposes BC on scalar quantities.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ieq,bcond(6),nx,ny,nz,iigx,iigy,iigz
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1)
      integer(4),optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder

c     Local variables

c     Begin program

c     Initialize interpolation order

      if (PRESENT(iorder)) then
        order = iorder
      else
        order = 2
      endif

c     Initialize BC grid quantities

      igxbc = iigx
      igybc = iigy
      igzbc = iigz

      nxbc = nx
      nybc = ny
      nzbc = nz

c     Initialize global limits

      if (     PRESENT(i1).and.PRESENT(i2)
     .    .and.PRESENT(j1).and.PRESENT(j2)
     .    .and.PRESENT(k1).and.PRESENT(k2)) then
        iimin = i1
        iimax = i2
        jjmin = j1
        jjmax = j2
        kkmin = k1
        kkmax = k2
      else
        iimin = grid_params%ilo(iigx)
        iimax = grid_params%ihi(iigx)
        jjmin = grid_params%jlo(iigy)
        jjmax = grid_params%jhi(iigy)
        kkmin = grid_params%klo(iigz)
        kkmax = grid_params%khi(iigz)
      endif

c     Fill PETSc ghost cells

#if defined(petsc)
      if ( (np > 1) .and. (.not.asm) .and. par_bc) then
        call fillPetscGhostCells(array,igxbc)
      endif
#endif

c     Impose BCs

      do bctype=1,BCLIM             !Enforces a particular order in the BCs (see grid_mod.f)
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
cc            if (bcond(ibc) == bctype) then
cc              call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
cc            endif
            if (abs(bcond(ibc)) == bctype) then
              if (bcond(ibc) < 0 .and. bctype == SYM) then
                symm = .true.
                call FillGhostNodes(ieq,1,1,dim,loc,DIR,array,array0)
                symm = .false.
              elseif (bcond(ibc) > 0 .and. bctype == SYM) then
                symm = .true.
                call FillGhostNodes(ieq,1,1,dim,loc,NEU,array,array0)
                symm = .false.
              else
                call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
              endif
            endif
          enddo
        enddo
      enddo

c     Fill edges

      do dim=1,3
        do loc=0,1
          call fillEdges(array,array0,dim,loc,bcond)
        enddo
      enddo

c     Fill corners

      do dim=1,3
        do loc=0,1
          call fillCorners(array,array0,dim,loc,bcond)
        enddo
      enddo

c     Periodic boundary synchronizing

cc      bctype=PER
cc      do dim=1,3
cc        do loc=0,1
cc          ibc = (1+loc)+2*(dim-1)
cc          if (bcond(ibc) == bctype) then
cc            call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
cc          endif
cc        enddo
cc      enddo

c     End program

      end subroutine imposeBConScalar

c     imposeBConVector
c     #################################################################
      subroutine imposeBConVector(fcomp,ncomp,nx,ny,nz,v_cnv,v_cov,var0
     .                          ,bcond,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                          ,iorder,is_cnv,is_vec)
c     -----------------------------------------------------------------
c     Imposes BC on vector field
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ncomp,bcond(6,ncomp),fcomp,nx,ny,nz,iigx,iigy,iigz
      real(8)    :: v_cnv(0:nx+1,0:ny+1,0:nz+1,ncomp)
     .             ,var0 (0:nx+1,0:ny+1,0:nz+1,ncomp)
     .             ,v_cov(0:nx+1,0:ny+1,0:nz+1,ncomp)
      integer(4),optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder
      logical   ,optional,intent(IN) :: is_cnv,is_vec

c     Local variables

      integer(4) :: ivar,ieq,ibc,loc,dim,bctype
     .             ,ilmin,jlmin,klmin,ilmax,jlmax,klmax
      logical    :: cov_to_cnv,to_cnv,is__vec

c     Begin program

      !By default, we assume input are vector components
      if (PRESENT(is_vec)) then
        is__vec = is_vec
      else
        is__vec = .true.
      endif

      if (is__vec .and. ncomp < 3) then
        call pstop('imposeBConVector'
     .            ,'Cannot impose BC on vector; # components < 3')
      endif

      !By default, we assume input is in cnv representation
      if (PRESENT(is_cnv)) then
        to_cnv = .not.is_cnv
      else
        to_cnv = .false.
      endif

c     Initialize interpolation order

      if (PRESENT(iorder)) then
        order = iorder
      else
        order = 2
      endif

c     Initialize BC grid quantities

      igxbc = iigx
      igybc = iigy
      igzbc = iigz

      nxbc = nx
      nybc = ny
      nzbc = nz

c     Initialize global limits

      if (     PRESENT(i1).and.PRESENT(i2)
     .    .and.PRESENT(j1).and.PRESENT(j2)
     .    .and.PRESENT(k1).and.PRESENT(k2)) then
        iimin = i1
        iimax = i2
        jjmin = j1
        jjmax = j2
        kkmin = k1
        kkmax = k2
      else
        iimin = grid_params%ilo(iigx)
        iimax = grid_params%ihi(iigx)
        jjmin = grid_params%jlo(iigy)
        jjmax = grid_params%jhi(iigy)
        kkmin = grid_params%klo(iigz)
        kkmax = grid_params%khi(iigz)
      endif

c     Fill PETSc ghost cells

#if defined(petsc)
      if ( (np > 1) .and. (.not.asm) .and. par_bc) then
        if (is__vec .and. to_cnv) then
          do ieq=1,ncomp
            call fillPetscGhostCells(v_cov(:,:,:,ieq),igxbc)
          enddo
        else
          do ieq=1,ncomp
            call fillPetscGhostCells(v_cnv(:,:,:,ieq),igxbc)
          enddo
        endif
      endif
#endif

c     Find complete curvilinear representation

      if (is__vec) then

#if defined(petsc)
       call fromGlobalToLocalLimits(iimin,jjmin,kkmin,ilmin,jlmin,klmin
     $                             ,igxbc,igybc,igzbc)
       call fromGlobalToLocalLimits(iimax,jjmax,kkmax,ilmax,jlmax,klmax
     $                             ,igxbc,igybc,igzbc)

       do k=klmin-1,klmax+1
         do j=jlmin-1,jlmax+1
           do i=ilmin-1,ilmax+1
#else
       do k=kkmin-1,kkmax+1
         do j=jjmin-1,jjmax+1
           do i=iimin-1,iimax+1
#endif
              call transformFromCurvToCurv(i,j,k,igxbc,igybc,igzbc
     .            ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .            ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),to_cnv)
           enddo
         enddo
       enddo

      endif

c     Impose BCs

      cov_to_cnv = .false.

      do bctype=1,BCLIM            !Enforces a particular order in the BCs (see grid_mod.f)
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            do ivar = 1,ncomp
              ieq = ivar + fcomp - 1

              if (abs(bcond(ibc,ivar)) == bctype) then

                if (bcond(ibc,ivar) < 0 .and. is__vec) then
                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
     .                               ,v_cov,var0,cnv=.false.)
                  cov_to_cnv = .true.
                else
                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
     .                               ,v_cnv,var0)
                endif

             endif

            enddo

          enddo
        enddo
      enddo

c     Synchronize covariant and contravariant components => contravariant

      if (cov_to_cnv) call synchronize_cnv(v_cnv,v_cov,bcond)

c     Fill edges

      do dim=1,3
        do loc=0,1
          do ivar=1,ncomp
            call fillEdges(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                    ,dim,loc,bcond(:,ivar))
          enddo
        enddo
      enddo

c     Fill corners

      do dim=1,3
        do loc=0,1
          do ivar=1,ncomp
            call fillCorners(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                      ,dim,loc,bcond(:,ivar))
          enddo
        enddo
      enddo

c     Periodic boundary synchronizing

cc      write (*,*) 'DIAG -- imposeBConVector'
cc      bctype=PER
cc      do dim=1,3
cc        do loc=0,1
cc          ibc = (1+loc)+2*(dim-1)
cc
cc          do ivar = 1,ncomp
cc            ieq = ivar + fcomp - 1
cc
cc            if (abs(bcond(ibc,ivar)) == bctype) then
cc              call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype,v_cnv
cc     .                               ,var0)
cc            endif
cc          enddo
cc
cc        enddo
cc      enddo

c     Find covariant components at ALL boundary ghost cells

      if (is__vec) then

        do dim = 1,3
          do loc = 0,1

          call findLoopLimits(dim,loc,iimin-1,iimax+1
     .                               ,jjmin-1,jjmax+1
     .                               ,kkmin-1,kkmax+1
     .                               ,imin,imax,jmin,jmax,kmin,kmax)

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax
                call transformFromCurvToCurv(i,j,k,igxbc,igybc,igzbc
     .            ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .            ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),.false.)
              enddo
            enddo
          enddo

          enddo
        enddo

      endif

c     End program

      end subroutine imposeBConVector

c     synchronize_cnv
c     #################################################################
      subroutine synchronize_cnv(v_cnv,v_cov,bcond)
c     -----------------------------------------------------------------
c     Finds all contravariant components at boundaries with BC in terms
c     of covariant components.
c
c     On input, tangential covariant components and normal contravariant
c     components are known at ghost cells. On output, all contravariant
c     components are known at ghost cells.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: bcond(6,3)
      real(8)    :: v_cnv(0:nxbc+1,0:nybc+1,0:nzbc+1,3)
     .             ,v_cov(0:nxbc+1,0:nybc+1,0:nzbc+1,3)

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

              if (do__nothing(dim,loc
     .                       ,iimin,iimax,jjmin,jjmax,kkmin,kkmax)) exit

              call findLoopLimits(dim,loc,iimin,iimax
     .                                   ,jjmin,jjmax
     .                                   ,kkmin,kkmax
     .                                   ,imin,imax,jmin,jmax,kmin,kmax)

              select case (ibc)
              case (1)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i-1,j,k,:,:)

                      v_cov(i-1,j,k,1) = -(gsuper(1,2)*v_cov(i-1,j,k,2)
     .                                    +gsuper(1,3)*v_cov(i-1,j,k,3)
     .                                    -v_cnv(i-1,j,k,1))/gsuper(1,1)

                      call transformFromCurvToCurv(i-1,j,k
     .               ,igxbc,igybc,igzbc
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

                      gsuper = gmetric%grid(igxbc)%gsup(i+1,j,k,:,:)

                      v_cov(i+1,j,k,1) = -(gsuper(1,2)*v_cov(i+1,j,k,2)
     .                                    +gsuper(1,3)*v_cov(i+1,j,k,3)
     .                                    -v_cnv(i+1,j,k,1))/gsuper(1,1)

                      call transformFromCurvToCurv(i+1,j,k
     .               ,igxbc,igybc,igzbc
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

                      gsuper = gmetric%grid(igxbc)%gsup(i,j-1,k,:,:)

                      v_cov(i,j-1,k,2) = -(gsuper(2,1)*v_cov(i,j-1,k,1)
     .                                    +gsuper(2,3)*v_cov(i,j-1,k,3)
     .                                    -v_cnv(i,j-1,k,2))/gsuper(2,2)

                      call transformFromCurvToCurv(i,j-1,k
     .               ,igxbc,igybc,igzbc
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

                      gsuper = gmetric%grid(igxbc)%gsup(i,j+1,k,:,:)

                      v_cov(i,j+1,k,2) = -(gsuper(2,1)*v_cov(i,j+1,k,1)
     .                                    +gsuper(2,3)*v_cov(i,j+1,k,3)
     .                                    -v_cnv(i,j+1,k,2))/gsuper(2,2)

                      call transformFromCurvToCurv(i,j+1,k
     .               ,igxbc,igybc,igzbc
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

                      gsuper = gmetric%grid(igxbc)%gsup(i,j,k-1,:,:)

                      v_cov(i,j,k-1,3) = -(gsuper(3,1)*v_cov(i,j,k-1,1)
     .                                    +gsuper(3,2)*v_cov(i,j,k-1,2)
     .                                    -v_cnv(i,j,k-1,3))/gsuper(3,3)

                      call transformFromCurvToCurv(i,j,k-1
     .               ,igxbc,igybc,igzbc
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

                      gsuper = gmetric%grid(igxbc)%gsup(i,j,k+1,:,:)

                      v_cov(i,j,k+1,3) = -(gsuper(3,1)*v_cov(i,j,k+1,1)
     .                                    +gsuper(3,2)*v_cov(i,j,k+1,2)
     .                                    -v_cnv(i,j,k+1,3))/gsuper(3,3)

                      call transformFromCurvToCurv(i,j,k+1
     .               ,igxbc,igybc,igzbc
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
     .                         ,array,array0,cnv)

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
c       * cov    -> whether vector components are covariant (optional)
c     ---------------------------------------------------------------

      implicit none       !For safe fortran

c     Call variables

      integer(4) :: ieq,dim,loc,bctype,nvar,ivar
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,nvar)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1,nvar)
      logical,optional :: cnv

c     Local variables

      integer(4) :: neq,ibc,i,j,k,ig,jg,kg
      logical    :: cov

c     Begin program

      nnvar = nvar

      if (PRESENT(cnv)) then
        cov = .not.cnv
      else
        cov = .false.
      endif
      
c     Determine boundary limits

      if(do__nothing(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax))return

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                           ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

c     Allocate rhs

      select case(dim)
      case(1)
        allocate(rhs(0:nybc+1,0:nzbc+1))
      case(2)
        allocate(rhs(0:nxbc+1,0:nzbc+1))
      case(3)
        allocate(rhs(0:nxbc+1,0:nybc+1))
      end select

c     Find BC update

      ibc = (1+loc)+2*(dim-1)

      if (nvar == 1) then
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,nvar),ibc)
        case(EQU)
          call dirichletBC(array(:,:,:,nvar),array0(:,:,:,nvar)
     .                    ,ieq,dim,loc,order)
        case(DIR)
          call dirichletBC(array(:,:,:,nvar),zeros,ieq,dim,loc,1)
        case(NEU)
          call neumannBC(array(:,:,:,nvar),ieq,dim,loc)
cc        case(SYM)
cc          call neumannBC(array(:,:,:,nvar),ieq,dim,loc)
        case(EXT)
          call fillFaces(array(:,:,:,nvar),dim,loc)
        case(SP)
          call scalarSingularBC(array(:,:,:,nvar))
        case(IFC)
          call interpolateToFace(array(:,:,:,nvar),ibc)
        case default
          call pstop('fillGhostNodes'
     .             ,'BC'//int2char(bctype)//' not implemented')
        end select
      else
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,ivar),ibc)
        case(EQU)
          call dirichletBC(ivar,array,array0,ieq,dim,loc,order)
        case(DIR)
          call dirichletBC(ivar,array,vzeros,ieq,dim,loc,1)
        case(NEU)
          call neumannBC(ivar,array,ieq,dim,loc)
cc        case(SYM)
cc          call neumannBC(ivar,array,ieq,dim,loc,symm=.true.)
        case(EXT)
          call fillFaces(array(:,:,:,ivar),dim,loc)
        case(SP)
          call vectorSingularBC(ivar,array(:,:,:,ivar),cov)
        case(IFC)
          call interpolateToFace(array(:,:,:,ivar),ibc)
        case default
          call pstop('fillGhostNodes'
     .             ,'BC'//int2char(bctype)//' not implemented')
        end select
      endif

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
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

c     Begin program

      if (np == 1 .or. asm) then

      select case (ibc)
      case (1)                  !x0
       if (.not.asm_dir(1))
     . array(0     ,jmin:jmax,kmin:kmax)=array(nxbc,jmin:jmax,kmin:kmax)
      case (2)                  !x1                                       
       if (.not.asm_dir(1))
     . array(nxbc+1,jmin:jmax,kmin:kmax)=array(1   ,jmin:jmax,kmin:kmax)
      case (3)                  !y0                                       
       if (.not.asm_dir(2))
     . array(imin:imax,0     ,kmin:kmax)=array(imin:imax,nybc,kmin:kmax)
      case (4)                  !y1                                       
       if (.not.asm_dir(2))
     . array(imin:imax,nybc+1,kmin:kmax)=array(imin:imax,1   ,kmin:kmax)
      case (5)                  !z0                                       
       if (.not.asm_dir(3))
     . array(imin:imax,jmin:jmax,0     )=array(imin:imax,jmin:jmax,nzbc)
      case (6)                  !z1                                       
       if (.not.asm_dir(3))
     . array(imin:imax,jmin:jmax,nzbc+1)=array(imin:imax,jmin:jmax,1   )
      case default
       call pstop('periodicBC'
     .           ,'Boundary '//int2char(ibc)//' non existent')
      end select

cc      write (*,*) 'DIAG -- periodicBC'
cc
cc      select case (ibc)
cc      case (1)                  !x0
cc       array(0     ,:,:)=array(nxbc,:,:)
cc      case (2)                  !x1
cc       array(nxbc+1,:,:)=array(1   ,:,:)
cc      case (3)                  !y0
cc       array(:,0     ,:)=array(:,nybc,:)
cc      case (4)                  !y1
cc       array(:,nybc+1,:)=array(:,1   ,:)
cc      case (5)                  !z0
cc       array(:,:,0     )=array(:,:,nzbc)
cc      case (6)                  !z1
cc       array(:,:,nzbc+1)=array(:,:,1   )
cc      case default
cc       write (*,*) 'Boundary',ibc,' non existent'
cc       stop
cc      end select

      endif

c     End program

      end subroutine periodicBC

c     interpolateToFace
c     #################################################################
      subroutine interpolateToFace(array,ibc)
c     -----------------------------------------------------------------
c     Interpolates values to face for EQU BC. On input:
c        * array -> 3D array to impose BC's on.
c        * ibc   -> integer boundary identifier
c                   (0 -> X1, 1->X2, 2->y1, 3->y2, 4->z1, 5->z2)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ibc
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

c     Begin program

cc      if (np == 1 .or. asm) then

      select case (ibc)
      case (1)                  !x0
cc       if (.not.asm_dir(1))
cc     . array(0     ,jmin:jmax,kmin:kmax)=
       array(0     ,jmin:jmax,kmin:kmax)=
     .      0.5*(array(0,jmin:jmax,kmin:kmax)
     .          +array(1,jmin:jmax,kmin:kmax))
      case (2)                  !x1                                       
cc       if (.not.asm_dir(1))
cc     . array(nxbc+1,jmin:jmax,kmin:kmax)=
       array(nxbc+1,jmin:jmax,kmin:kmax)=
     .      0.5*(array(nxbc+1,jmin:jmax,kmin:kmax)
     .          +array(nxbc  ,jmin:jmax,kmin:kmax))
      case (3)                  !y0                                       
cc       if (.not.asm_dir(2))
cc     . array(imin:imax,0     ,kmin:kmax)=
       array(imin:imax,0     ,kmin:kmax)=
     .      0.5*(array(imin:imax,0,kmin:kmax)
     .          +array(imin:imax,1,kmin:kmax))
      case (4)                  !y1                                       
cc       if (.not.asm_dir(2))
cc     . array(imin:imax,nybc+1,kmin:kmax)=
       array(imin:imax,nybc+1,kmin:kmax)=
     .      0.5*(array(imin:imax,nybc+1,kmin:kmax)
     .          +array(imin:imax,nybc  ,kmin:kmax))
      case (5)                  !z0                                       
cc       if (.not.asm_dir(3))
cc     . array(imin:imax,jmin:jmax,0     )=
       array(imin:imax,jmin:jmax,0     )=
     .      0.5*(array(imin:imax,jmin:jmax,0)
     .          +array(imin:imax,jmin:jmax,1))
      case (6)                  !z1                                       
cc       if (.not.asm_dir(3))
cc     . array(imin:imax,jmin:jmax,nzbc+1)=
       array(imin:imax,jmin:jmax,nzbc+1)=
     .      0.5*(array(imin:imax,jmin:jmax,nzbc+1)
     .          +array(imin:imax,jmin:jmax,nzbc  ))
      case default
       call pstop('interpolateToFace'
     .           ,'Boundary '//int2char(ibc)//' non existent')
      end select

cc      endif

c     End program

      end subroutine interpolateToFace

c     fillFaces
c     #################################################################
      subroutine fillFaces(array,dim,loc)
c     -----------------------------------------------------------------
c     Fills ghost cells in faces by extrapolation (6 in 3D). On input:
c        * array  -> 3D array to impose BC's on.
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     This routine is only used when no good information is available
c     to fill the ghost cells of a particular quantity (for instance,
c     during MG restriction of arrays that contain EQU BCs).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: dim,loc
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer(4) :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc,bcond(6)

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      bcond = 0
      bcond(ibc) = EXT

c     Determine boundary limits

      if(do__nothing(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax))return

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                           ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      if (bc_debug) then
cc      write (*,*) 'Proc',my_rank,iimin,iimax,jjmin,jjmax,kkmin,kkmax
        write(*,*) '  fillFaces: Proc',my_rank
     .            ,'; Boundary: dim =',dim,' loc =',loc
     .            ,'; Local limits:',imin,imax,jmin,jmax,kmin,kmax
      endif

c     Select edge

      select case (ibc)
      case (1) !On face #1, fill face i=0

        if (imin == 1) then

          ito   = imin-1
          ifrom = min(imin+order,nxbc)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          do j=jmin,jmax
            do k=kmin,kmax
              call extrapolateToCorner(ibc,ito,j,k,imn,j,k
     .                  ,imx,j,k,array(imn:imx,j,k),bcond)
            enddo
          enddo

        endif

      case (2) !On face #2, fill face i=nx+1

        if (imax == nxbc) then

          ito   = imax+1
          ifrom = max(imax-order,1)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          do j=jmin,jmax
            do k=kmin,kmax
              call extrapolateToCorner(ibc,ito,j,k,imn,j,k
     .                  ,imx,j,k,array(imn:imx,j,k),bcond)
            enddo
          enddo

        endif

      case (3) !On face #3, fill face j=0

        if (jmin == 1) then

          jto   = jmin-1
          jfrom = min(jmin+order,nybc)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          do k=kmin,kmax
            do i=imin,imax
                call extrapolateToCorner(ibc,i,jto,k,i,jmn,k
     .                  ,i,jmx,k,array(i,jmn:jmx,k),bcond)
            enddo
          enddo

        endif

      case (4) !On face #4, fill face j=ny+1

        if (jmax == nybc) then

          jto   = jmax+1
          jfrom = max(jmax-order,1)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          do k=kmin,kmax
            do i=imin,imax
              call extrapolateToCorner(ibc,i,jto,k,i,jmn,k
     .                  ,i,jmx,k,array(i,jmn:jmx,k),bcond)
            enddo
          enddo

        endif

      case (5) !On face #5, fill face k=0

        if (kmin == 1) then

          kto   = kmin-1
          kfrom = min(kmin+order,nzbc)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          do i=imin,imax
            do j=jmin,jmax
              call extrapolateToCorner(ibc,i,j,kto,i,j,kmn
     .                  ,i,j,kmx,array(i,j,kmn:kmx),bcond)
            enddo
          enddo

        endif

      case (6) !On face #6, fill face k=nz+1

        if (kmax == nzbc) then

          kto   = kmax+1
          kfrom = max(kmax-order,1)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          do i=imin,imax
            do j=jmin,jmax
              call extrapolateToCorner(ibc,i,j,kto,i,j,kmn
     .                  ,i,j,kmx,array(i,j,kmn:kmx),bcond)
            enddo
          enddo

        endif

      end select

c     End program

      end subroutine fillFaces

c     fillEdges
c     #################################################################
      subroutine fillEdges(array,array0,dim,loc,bcond)
c     -----------------------------------------------------------------
c     Fills ghost cells in edges (12 in 3D). On input:
c        * array  -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: dim,loc,bcond(6)
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer(4) :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

c     Determine boundary limits

      if(do__nothing(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax))return

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                           ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      if (bc_debug) then
cc      write (*,*) 'Proc',my_rank,iimin,iimax,jjmin,jjmax,kkmin,kkmax
        write(*,*) '  fillEdges: Proc',my_rank
     .            ,'; Boundary: dim =',dim,' loc =',loc
     .            ,'; Local limits:',imin,imax,jmin,jmax,kmin,kmax
      endif

c     Select edge

      select case (ibc)
      case (1) !On face #1, fill edges y=jmin-1,y=jmax+1, vary z

        if (imin == 1) then
          ito   = imin-1
          ifrom = min(imin+order,nxbc)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          if (jmin == 1.and.(selectFace(ibc).or.selectFace(3))) then  !Check faces #1,#3
            jto   = jmin-1
            jfrom = min(jmin+order,nybc)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(1)==EQU.or.bcond(3)==EQU) .and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

          if (jmax == nybc.and.(selectFace(ibc).or.selectFace(4))) then  !Check faces #1,#4
            jto   = jmax+1
            jfrom = max(jmax-order,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(1)==EQU.or.bcond(4)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

        endif

      case (2) !On face #2, fill edges y=jmin-1,y=jmax+1, vary z

        if (imax == nxbc) then
          ito   = imax+1
          ifrom = max(imax-order,1)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          if (jmin == 1.and.(selectFace(ibc).or.selectFace(3))) then  !Check faces #2,#3
            jto   = jmin-1
            jfrom = min(jmin+order,nybc)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(2)==EQU.or.bcond(3)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

          if (jmax == nybc.and.(selectFace(ibc).or.selectFace(4))) then  !Check faces #2,#4
            jto   = jmax+1
            jfrom = max(jmax-order,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(2)==EQU.or.bcond(4)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

        endif

      case (3) !On face #3, fill edges z=kmin-1,z=kmax+1, vary x

        if (jmin == 1) then

          jto   = jmin-1
          jfrom = min(jmin+order,nybc)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          if (kmin == 1.and.(selectFace(ibc).or.selectFace(5))) then  !Check faces #3,#5

            kto   = kmin-1
            kfrom = min(kmin+order,nzbc)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(3)==EQU.or.bcond(5)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (kmax == nzbc.and.(selectFace(ibc).or.selectFace(6))) then  !Check faces #3,#6

            kto   = kmax+1
            kfrom = max(kmax-order,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(3)==EQU.or.bcond(6)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      case (4) !On face #4, fill edges z=kmin-1,z=kmax+1, vary x

        if (jmax == nybc) then

          jto   = jmax+1
          jfrom = max(jmax-order,1)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          if (kmin == 1.and.(selectFace(ibc).or.selectFace(5))) then  !Check faces #4,#5

            kto   = kmin-1
            kfrom = min(kmin+order,nzbc)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(4)==EQU.or.bcond(5)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (kmax == nzbc.and.(selectFace(ibc).or.selectFace(6))) then  !Check faces #4,#6

            kto   = kmax+1
            kfrom = max(kmax-order,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(4)==EQU.or.bcond(6)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      case (5) !On face #5, fill edges x=imin-1,x=imax+1, vary y

        if (kmin == 1) then

          kto   = kmin-1
          kfrom = min(kmin+order,nzbc)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          if (imin == 1.and.(selectFace(ibc).or.selectFace(1))) then  !Check faces #5,#1

            ito   = imin-1
            ifrom = min(imin+order,nxbc)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(5)==EQU.or.bcond(1)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (imax == nxbc.and.(selectFace(ibc).or.selectFace(2))) then  !Check faces #5,#2

            ito   = imax+1
            ifrom = max(imax-order,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(5)==EQU.or.bcond(2)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      case (6) !On face #6, fill edges x=imin-1,x=imax+1, vary y

        if (kmax == nzbc) then

          kto   = kmax+1
          kfrom = max(kmax-order,1)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          if (imin == 1.and.(selectFace(ibc).or.selectFace(1))) then  !Check faces #6,#1

            ito   = imin-1
            ifrom = min(imin+order,nxbc)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(6)==EQU.or.bcond(1)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (imax == nxbc.and.(selectFace(ibc).or.selectFace(2))) then  !Check faces #6,#2

            ito   = imax+1
            ifrom = max(imax-order,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(6)==EQU.or.bcond(2)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      end select

c     End program

      end subroutine fillEdges

c     fillCorners
c     #################################################################
      subroutine fillCorners(array,array0,dim,loc,bcond)
c     -----------------------------------------------------------------
c     Imposes ghost cells in corners (8 in 3D). On input:
c        * array -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: dim,loc,bcond(6)
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer(4) :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

c     Determine boundary limits

      if(do__nothing(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax))return

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                           ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

c     Select corner

      if (      imin == 1  .and. jmin == 1  .and. kmin == 1
     .    .and.(selectFace(1).or.selectFace(3).or.selectFace(5)) ) then   !Check concurrent faces
        ifrom = min(imin+order,nxbc)
        jfrom = min(jmin+order,nybc)
        kfrom = min(kmin+order,nzbc)
        ito   = imin-1
        jto   = jmin-1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write (*,*) '  fillCorners: Proc',my_rank
     .               ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .          ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imin == 1  .and. jmax == nybc .and. kmin == 1
     .    .and.(selectFace(1).or.selectFace(4).or.selectFace(5))) then

        ifrom = min(imin+order,nxbc)
        jfrom = max(jmax-order,1   )
        kfrom = min(kmin+order,nzbc)
        ito   = imin-1
        jto   = jmax+1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .           ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imax == nxbc .and. jmin == 1  .and. kmin == 1
     .    .and.(selectFace(2).or.selectFace(3).or.selectFace(5))) then

        ifrom = max(imax-order,1   )
        jfrom = min(jmin+order,nybc)
        kfrom = min(kmin+order,nzbc)
        ito   = imax+1
        jto   = jmin-1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imax == nxbc .and. jmax == nybc .and. kmin == 1
     .    .and.(selectFace(2).or.selectFace(4).or.selectFace(5))) then

        ifrom = max(imax-order,1   )
        jfrom = max(jmax-order,1   )
        kfrom = min(kmin+order,nzbc)
        ito   = imax+1
        jto   = jmax+1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imin == 1  .and. jmin == 1  .and. kmax == nzbc
     .    .and.(selectFace(1).or.selectFace(3).or.selectFace(6))) then

        ifrom = min(imin+order,nxbc)
        jfrom = min(jmin+order,nybc)
        kfrom = max(kmax-order,1   )
        ito   = imin-1
        jto   = jmin-1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imin == 1  .and. jmax == nybc .and. kmax == nzbc
     .    .and.(selectFace(1).or.selectFace(4).or.selectFace(6))) then

        ifrom = min(imin+order,nxbc)
        jfrom = max(jmax-order,1   )
        kfrom = max(kmax-order,1   )
        ito   = imin-1
        jto   = jmax+1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imax == nxbc .and. jmin == 1  .and. kmax == nzbc
     .    .and.(selectFace(2).or.selectFace(3).or.selectFace(6))) then

        ifrom = max(imax-order,1   )
        jfrom = min(jmin+order,nybc)
        kfrom = max(kmax-order,1   )
        ito   = imax+1
        jto   = jmin-1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imax == nxbc .and. jmax == nybc .and. kmax == nzbc
     .    .and.(selectFace(2).or.selectFace(4).or.selectFace(6))) then

        ifrom = max(imax-order,1 )
        jfrom = max(jmax-order,1 )
        kfrom = max(kmax-order,1 )
        ito   = imax+1
        jto   = jmax+1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)

      endif

c     End program

      end subroutine fillCorners

c     extrapolateToCorner
c     #################################################################
      subroutine extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .                              ,imx,jmx,kmx,array,bcond)
c     -----------------------------------------------------------------
c     Routine that extrapolates information to edge or vertex corners.
c     It is the workhorse of the fillEdge and fillCorner routines.
c     In the call sequence:
c        * ibc: boundary identifier (integer)
c        * ito,jto,kto: coordinates of ghost cell to be filled (integer)
c        * imn,jmn,kmn: node coordinates identifying start of
c                       interpolation region (integer)
c        * imx,jmx,kmx: node coordinates identifying end of interpolation
c                       region (integer)
c        * array: local array defined on interpolation region, including
c                 ghost cell (ito,jto,kto) (real)
c     The order of extrapolation 'order' is defined in the module header.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ibc,ito,jto,kto,imn,jmn,kmn,imx,jmx,kmx,bcond(6)
      real(8)    :: array(imn:imx,jmn:jmx,kmn:kmx)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,inx,iny,inz
     .             ,order_int(3),offset
     .             ,imni,imxi,jmni,jmxi,kmni,kmxi,num
     .             ,ibcl(3),i_bc,i1,i2

      real(8)    :: valx(1),valy(1),valz(1),xto(1),yto(1),zto(1)
     .             ,smooth(3)

      logical    :: flg(3),select_bc(3)

c     Extrapolation

      real(8),allocatable,dimension(:) :: xx,yy,zz

c     Begin program

c     Identify directions for extrapolation.

      inx = imx-imn
      iny = jmx-jmn
      inz = kmx-kmn

      flg = .true.

      !Order in each direction
      order_int(1) = min(order,max(inx-1,0))
      order_int(2) = min(order,max(iny-1,0))
      order_int(3) = min(order,max(inz-1,0))

      !Smoothness indicator
      smooth = 0d0
      if (order >= 2) then
        if (order_int(1) == order) then
          if (ito == imn) then
            smooth(1) = array(ito+1,jto,kto)-array(ito+2,jto,kto)
          else
            smooth(1) = array(ito-1,jto,kto)-array(ito-2,jto,kto)
          endif
        endif

        if (order_int(2) == order) then
          if (jto == jmn) then
            smooth(2) = array(ito,jto+1,kto)-array(ito,jto+2,kto)
          else                                              
            smooth(2) = array(ito,jto-1,kto)-array(ito,jto-2,kto)
          endif
        endif

        if (order_int(3) == order) then
          if (kto == kmn) then
            smooth(3) = array(ito,jto,kto+1)-array(ito,jto,kto+2)
          else                                                  
            smooth(3) = array(ito,jto,kto-1)-array(ito,jto,kto-2)
          endif
        endif

        smooth = abs(smooth)
      endif

      !Do not interpolate **along** direction if
      !   1) not enough points,
      if (inx == 0) flg(1) = .false.
      if (iny == 0) flg(2) = .false.
      if (inz == 0) flg(3) = .false.

      !Select boundaries involved 
      ![ibc = (1+loc)+2*(dim-1), with dim=1,2,3, and loc=0(left),1(right)]

      ibcl = 0

      if (flg(1) .and. flg(2) .and. flg(3)) then !Corner

        if (ito == 0 .and. jto == 0 .and. kto == 0) then
          ibcl(1) = 1
          ibcl(2) = 3
          ibcl(3) = 5
        elseif (ito == 0 .and. jto == 0 .and. kto /= 0) then
          ibcl(1) = 1
          ibcl(2) = 3
          ibcl(3) = 6
        elseif (ito == 0 .and. jto /= 0 .and. kto == 0) then
          ibcl(1) = 1
          ibcl(2) = 4
          ibcl(3) = 5
        elseif (ito /= 0 .and. jto == 0 .and. kto == 0) then
          ibcl(1) = 2
          ibcl(2) = 3
          ibcl(3) = 5
        elseif (ito /= 0 .and. jto /= 0 .and. kto == 0) then
          ibcl(1) = 2
          ibcl(2) = 4
          ibcl(3) = 5
        elseif (ito /= 0 .and. jto == 0 .and. kto /= 0) then
          ibcl(1) = 2
          ibcl(2) = 3
          ibcl(3) = 6
        elseif (ito == 0 .and. jto /= 0 .and. kto /= 0) then
          ibcl(1) = 1
          ibcl(2) = 4
          ibcl(3) = 6
        elseif (ito /= 0 .and. jto /= 0 .and. kto /= 0) then
          ibcl(1) = 2
          ibcl(2) = 4
          ibcl(3) = 6
        endif

      else !Edge

        select case(ibc)
        case (1,2)
          ibcl(1) = ibc
          if (jto == 0) then
            ibcl(2) = 3
          else
            ibcl(2) = 4
          endif
        case (3,4)
          ibcl(2) = ibc
          if (kto == 0) then
            ibcl(3) = 5
          else
            ibcl(3) = 6
          endif
        case (5,6)
          ibcl(3) = ibc
          if (ito == 0) then
            ibcl(1) = 1
          else
            ibcl(1) = 2
          endif
        end select

      endif
        
      !Treat special boundaries

      !Identify special boundaries to interpolate **along**
      !   1) equilibrium boundary
      !   2) SP boundary
      select_bc = .false.
      do i=1,3
        i_bc = ibcl(i)
        if (i_bc == 0) cycle
        select_bc(i) =    (bcond(i_bc) == EQU)
     .                .or.(bcSP(ibc=i_bc)    )
      enddo

      call select_along

      !Identify special boundaries to interpolate **across**
      !   1) periodic boundary
      !   2) Symmetry boundary
      select_bc = .false.
      do i=1,3
        i_bc = ibcl(i)
        if (i_bc == 0) cycle
        select_bc(i) =    (    bcond(i_bc)  == PER)
     .                .or.(abs(bcond(i_bc)) == SYM)
      enddo

      call select_across

      !Check some direction for interpolation exists

      if ((.not.flg(1)).and.(.not.flg(2)).and.(.not.flg(3))) then
        write (*,*) 'Extrapolation directions',flg
        write (*,*) 'Boundaries involved',ibcl
        write (*,*) 'No extrapolation directions left'
        call pstop('extrapolateToCorner','')
      endif

c     Initialize accumulators

      valx = 0d0
      valy = 0d0
      valz = 0d0

      num = 0

c     X extrapolation

      if (flg(1)) then

        allocate(xx(inx))

        if (ito == imn) then
          call getMGmap(imn+1,jmn,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          xx = grid_params%xx(ig:ig+inx-1)
          imni = imn+1
          imxi = imx
        else
          call getMGmap(imx-1,jmn,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          xx = grid_params%xx(ig-(inx-1):ig)
          imni = imn
          imxi = imx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        xto = grid_params%xx(ig)

        if (abs(bcond(ibcl(1))) == SYM) then  !Exception for symmetry BC
          if (ito == imn) then
            valx(1) = sign(1,bcond(1))*array(ito+2,jto,kto)
          else
            valx(1) = sign(1,bcond(1))*array(ito-2,jto,kto)
          endif
        else
          call IntDriver1d(inx,xx,array(imni:imxi,jto,kto)
     .                  ,1,xto,valx,order_int(1))
        endif

        num = num + 1

        deallocate(xx)
      endif

c     Y extrapolation

      if (flg(2)) then

        allocate(yy(iny))

        if (jto == jmn) then
          call getMGmap(imn,jmn+1,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          yy = grid_params%yy(jg:jg+iny-1)
          jmni = jmn+1
          jmxi = jmx
        else
          call getMGmap(imn,jmx-1,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          yy = grid_params%yy(jg-(iny-1):jg)
          jmni = jmn
          jmxi = jmx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        yto = grid_params%yy(jg)

        if (abs(bcond(ibcl(2))) == SYM) then  !Exception for symmetry BC
          if (jto == jmn) then
            valy(1) = sign(1,bcond(1))*array(ito,jto+2,kto)
          else                                                   
            valy(1) = sign(1,bcond(1))*array(ito,jto-2,kto)
          endif
        else
          call IntDriver1d(iny,yy,array(ito,jmni:jmxi,kto)
     .                  ,1,yto,valy,order_int(2))
        endif

        num = num + 1

        deallocate(yy)

      endif

c     Z extrapolation

      if (flg(3)) then

        allocate(zz(inz))

        if (kto == kmn) then
          call getMGmap(imn,jmn,kmn+1,igxbc,igybc,igzbc,ig,jg,kg)
          zz = grid_params%zz(kg:kg+inz-1)
          kmni = kmn+1
          kmxi = kmx
        else
          call getMGmap(imn,jmn,kmx-1,igxbc,igybc,igzbc,ig,jg,kg)
          zz = grid_params%zz(kg-(inz-1):kg)
          kmni = kmn
          kmxi = kmx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        zto = grid_params%zz(kg)

        if (abs(bcond(ibcl(3))) == SYM) then  !Exception for symmetry BC
          if (kto == kmn) then
            valz(1) = sign(1,bcond(1))*array(ito,jto,kto+2)
          else                                                   
            valz(1) = sign(1,bcond(1))*array(ito,jto,kto-2)
          endif
        else
          call IntDriver1d(inz,zz,array(ito,jto,kmni:kmxi)
     .                  ,1,zto,valz,order_int(3))
        endif

        num = num + 1

        deallocate(zz)
      endif

c     Average extrapolation results and fill ghost cell

      array(ito,jto,kto) = (valx(1) + valy(1) + valz(1))/num

c     End program

      contains

c     select_along
c     ######################################################################
      subroutine select_along
c     -----------------------------------------------------------------------
c     Selects interpolation directions so that interpolation occurs ALONG
c     selected boundary only.
c     -----------------------------------------------------------------------

      do i=1,3

        !Identify other directions
        i1 = mod(i+1,3)
        i2 = mod(i+2,3)
        if (i1 == 0) i1 = 3
        if (i2 == 0) i2 = 3

        !Identify cases in which an edge has two or more selected directions:
        !  can't eliminate both interpolation directions simultaneously;
        !  resolve indetermination by choosing smoothest direction
        if (select_bc(i) .and. (select_bc(i1) .or. select_bc(i2))) then

cc          if (select_bc(i1).and.select_bc(i2)) then
cc            if (order_int(i1) <= order_int(i2)) then
cc              select_bc(i1) = .false.
cc            else
cc              select_bc(i2) = .false.
cc            endif
cc          endif
cc
cc          if (select_bc(i1)) then
cc            if (order_int(i) <= order_int(i1)) then
cc              cycle
cc            else
cc              flg(i) = .not.select_bc(i)
cc            endif
cc          else
cc            if (order_int(i) <= order_int(i2)) then
cc              cycle
cc            else
cc              flg(i) = .not.select_bc(i)
cc            endif
cc          endif

          if (select_bc(i1).and.select_bc(i2)) then
            if (smooth(i1) <= smooth(i2)) then
              select_bc(i1) = .false.
            else
              select_bc(i2) = .false.
            endif
          endif

          if (select_bc(i1)) then
            if (smooth(i) <= smooth(i1)) then
              cycle
            else
              flg(i) = .not.select_bc(i)
            endif
          else
            if (smooth(i) <= smooth(i2)) then
              cycle
            else
              flg(i) = .not.select_bc(i)
            endif
          endif

        elseif (select_bc(i)) then
          flg(i) = .not.select_bc(i)
        endif

      enddo

      end subroutine select_along

c     select_across
c     ######################################################################
      subroutine select_across
c     -----------------------------------------------------------------------
c     Selects interpolation directions so that interpolation occurs ACROSS
c     selected boundary only.
c     -----------------------------------------------------------------------

      do i=1,3

        !Identify other directions
        i1 = mod(i+1,3)
        i2 = mod(i+2,3)
        if (i1 == 0) i1 = 3
        if (i2 == 0) i2 = 3

        if (select_bc(i) .and. flg(i)) then
          flg(i1) = .false.
          flg(i2) = .false.
        endif

      enddo

      end subroutine select_across

      end subroutine extrapolateToCorner

      end module imposeBCinterface
