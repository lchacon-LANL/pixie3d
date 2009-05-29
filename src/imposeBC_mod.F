c module local_BCS_variables
c####################################################################
      module local_BCS_variables

        use BCS_variables

        use imposeBCinterface, ONLY: extrapolateBC

        use auxiliaryVariables

        use operators

        use icond

        use timestepping

        use generalPurposeFunctions

        use equilibrium

        logical :: pinch_flow=.false.

        integer :: bc_order=1

      contains

c     setBCsched
c     ###################################################################
      subroutine setBCsched

c     -------------------------------------------------------------------
c     Sets BC schedule for dependent and auxiliary variables, which
c     are identified by a 2-dim integer array, bc_seq. The first
c     dimension sets order of treatment. The second dimension (with two
c     2 components) has the following convention:
c       * First component identifies whether quantities are dependent (>0)
c         or auxiliary (<0) variables.
c       * Second component identifies scalars (=0) or vectors (=1).
c     -------------------------------------------------------------------

c     Begin program

      if (.not.associated(gv%bc_seq)) then

#if defined(vec_pot)
      gv%nbc_seq = 6
#else
      gv%nbc_seq = 5
#endif

      allocate(gv%bc_seq(gv%nbc_seq,2))

#if defined(vec_pot)
      gv%bc_seq(1,1) = IRHO !Dep. Variable (>0)
      gv%bc_seq(1,2) = 0    !Scalar

      gv%bc_seq(2,1) = IAX  !Dep. variable (>0)
      gv%bc_seq(2,2) = 1    !Vector

      gv%bc_seq(3,1) =-IBX  !Aux. variable (<0)
      gv%bc_seq(3,2) = 1    !Vector

      gv%bc_seq(4,1) =-IJX  !Aux. variable (<0)
      gv%bc_seq(4,2) = 1    !Vector

      gv%bc_seq(5,1) = IVX  !Dep. variable (>0)
      gv%bc_seq(5,2) = 1    !Vector

      gv%bc_seq(6,1) = ITMP !Dep. variable (>0)
      gv%bc_seq(6,2) = 0    !Scalar
#else
      gv%bc_seq(1,1) = IRHO !Dep. Variable (>0)
      gv%bc_seq(1,2) = 0    !Scalar

      gv%bc_seq(2,1) = IBX  !Dep. variable (>0)
      gv%bc_seq(2,2) = 1    !Vector

      gv%bc_seq(3,1) =-IJX  !Aux. variable (<0)
      gv%bc_seq(3,2) = 1    !Vector

      gv%bc_seq(4,1) = IVX  !Dep. variable (>0)
      gv%bc_seq(4,2) = 1    !Vector

      gv%bc_seq(5,1) = ITMP !Dep. variable (>0)
      gv%bc_seq(5,2) = 0    !Scalar
#endif

      endif

c     End program

      end subroutine setBCsched

c     ExB_flow_bc
c     #################################################################
      subroutine ExB_flow_bc(igrid,nx,ny,nz,bcnv,bcov,v0,dbcov)

c     -----------------------------------------------------------------
c     Finds ExB flow boundary condition when Eo is imposed at boundary
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,nx,ny,nz
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3) :: bcnv,bcov,v0
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3),optional :: dbcov

c     Local variables

      integer :: i,j,k
      real(8) :: b2,bbcnv(3),bbcov(3),dum(3)

c     Begin program

      if (.not.pinch_flow) return

      i = nx

      do k=1,nz
        do j=1,ny
          bbcnv = 0.5*(bcnv(i,j,k,:)+bcnv(i+1,j,k,:))
          bbcov = 0.5*(bcov(i,j,k,:)+bcov(i+1,j,k,:))

          b2 = dot_product(bbcnv,bbcov)

          if (b2 > 0d0) then
            if (PRESENT(dbcov)) then  !Linearize ExB flow BC
              dum = 0.5*(dbcov(i,j,k,:)+dbcov(i+1,j,k,:))
              bbcov = dum - 2*dot_product(bbcnv,dum)/b2*bbcov
            endif

            b2 = 2*b2/(gmetric%grid(igrid)%jac(i  ,j,k)
     .                +gmetric%grid(igrid)%jac(i+1,j,k))

            v0(i+1,j,k,1) = -(bbcov(2)*E0(3)-bbcov(3)*E0(2))/b2 !ExB velocity
          endif
        enddo
      enddo

c     End program

      end subroutine ExB_flow_bc

      end module local_BCS_variables

c scalarDirichletBC
c #####################################################################
      subroutine scalarDirichletBC(array,array0,ieq,dim,loc,order)
c ---------------------------------------------------------------------
c     Imposes dirichlet BC. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c        * order -> order of extrapolation (when used)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Dirichlet BC.
c ---------------------------------------------------------------------

      use local_BCS_variables

      implicit none

c     Call variables

      integer    :: ieq,dim,loc,order
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

c     Begin program

cc      call extrapolateBC(array,array0,dim,loc,order)
      call extrapolateBC(array,array0,dim,loc,1)

c     End program

      end subroutine scalarDirichletBC

c vectorDirichletBC
c #####################################################################
      subroutine vectorDirichletBC(ivar,array,array0,ieq,dim,loc,order)
c ---------------------------------------------------------------------
c     Imposes dirichlet BC. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c        * order -> order of extrapolation (when used)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Dirichlet BC.
c ---------------------------------------------------------------------
        
      use local_BCS_variables

      implicit none

c     Call variables

      integer    :: ieq,dim,loc,order,ivar
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,*)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1,*)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: i,j,k,ig,jg,kg,nvar,ibc
      real(8)    :: x1,x2,x3,dh(3),diver

      real(8)    :: gsuper(3,3),jac0

c     Begin program

      rhs = 0d0

      nvar = nnvar  !Where is this needed?

      ibc = (1+loc)+2*(dim-1)

cc      call findBCLoopLimits(dim,loc,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      select case (ieq)
      case (IVX,IVY,IVZ)

        call extrapolateBC(array(:,:,:,ivar),array0(:,:,:,ivar)
cc     .                    ,dim,loc,order)
     .                    ,dim,loc,1)

      case (IBX,IBY,IBZ,IJX,IJY,IJZ) !Imposes divergence-free constraint on B, J

        if (ivar /= dim) then

          call extrapolateBC(array (:,:,:,ivar),array0(:,:,:,ivar)
     .                      ,dim,loc,order)

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

      case (IAX,IAY,IAZ) !Imposes divergence-free constraint on A

        if (ivar /= dim) then

          call extrapolateBC(array (:,:,:,ivar),array0(:,:,:,ivar)
     .                      ,dim,loc,order)

        else

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = dxh(ig)
                dh(2) = dyh(jg)
                dh(3) = dzh(kg)

                select case (ibc)
                case (1)
                  array(i-1,j,k,dim) = array(i+1,j,k,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(j,k) = array(i,j,k,dim) + dh(dim)*diver
                case (2)
                  array(i+1,j,k,dim) = array(i-1,j,k,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(j,k) = array(i,j,k,dim) - dh(dim)*diver
                case (3)
                  array(i,j-1,k,dim) = array(i,j+1,k,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(i,k) = array(i,j,k,dim) + dh(dim)*diver
                case (4)
                  array(i,j+1,k,dim) = array(i,j-1,k,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(i,k) = array(i,j,k,dim) - dh(dim)*diver
                case (5)
                  array(i,j,k-1,dim) = array(i,j,k+1,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(i,j) = array(i,j,k,dim) + dh(dim)*diver
                case (6)
                  array(i,j,k+1,dim) = array(i,j,k-1,dim)
                  diver =jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                            ,array(:,:,:,1)
     .                            ,array(:,:,:,2)
     .                            ,array(:,:,:,3))
                  rhs(i,j) = array(i,j,k,dim) - dh(dim)*diver
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

c scalarNeumannBC
c #####################################################################
      subroutine scalarNeumannBC(array,ieq,dim,loc)
c ---------------------------------------------------------------------
c     Imposes neumann BC for a scalar. On input:
c        * array -> variable upon which BCs are imposed
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Neumann BC.
c ---------------------------------------------------------------------

      use local_BCS_variables

      implicit none

c     Call variables

      integer    :: ieq,dim,loc
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: i,j,k,ig,jg,kg,ip,im,jp,jm,kp,km,nvar,ibc,icomp
      real(8)    :: x1,x2,x3,dh(3),gsuper(3,3),jac0,hessian(3,3,3)

c     Begin program

      rhs = 0d0

      nvar = 1    !Where is this needed?

      ibc = (1+loc)+2*(dim-1)

      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      select case (ieq)
cc      case (IRHO,ITMP)
      case (IRHO)

        rhs = 0d0   !This shuts off advective fluxes for v.n=0

cc        do i=imin,imax
cc          do j=jmin,jmax
cc            do k=kmin,kmax
cc
cc              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
cc
cc              gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)
cc
cc              ip = min(i+1,nxbc)
cc              im = max(i-1,1)
cc              jp = min(j+1,nybc)
cc              jm = max(j-1,1)
cc              kp = min(k+1,nzbc)
cc              km = max(k-1,1)
cc
cc              dh(1) = 2.*dxh(ig)
cc              if (i == nxbc) dh(1) = dx(ig-1)
cc              if (i == 1   ) dh(1) = dx(ig)
cc
cc              dh(2) = 2.*dyh(jg)
cc              if (j == nybc) dh(2) = dy(jg-1)
cc              if (j == 1   ) dh(2) = dy(jg)
cc
cc              dh(3) = 2.*dzh(kg)
cc              if (k == nzbc) dh(3) = dz(kg-1)
cc              if (k == 1   ) dh(3) = dz(kg)
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
                  rhs(j,k) =  hessian(1,1,1)*vcnv(i,j,k,1)*vcnv(i,j,k,1)
     .                       +hessian(1,2,2)*vcnv(i,j,k,2)*vcnv(i,j,k,2)
     .                       +hessian(1,3,3)*vcnv(i,j,k,3)*vcnv(i,j,k,3)
     .                    +2.*hessian(1,1,2)*vcnv(i,j,k,1)*vcnv(i,j,k,2)
     .                    +2.*hessian(1,1,3)*vcnv(i,j,k,1)*vcnv(i,j,k,3)
     .                    +2.*hessian(1,2,3)*vcnv(i,j,k,2)*vcnv(i,j,k,3)
                else
                  rhs(j,k) = 0d0
                endif
                rhs(j,k) = -dh(dim)
     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
     .             -0.5/jac0*rhs(j,k))/gsuper(dim,dim)
              elseif (dim == 2) then
                if (gamma > 1d0) then
                  rhs(i,k) =  hessian(2,1,1)*vcnv(i,j,k,1)*vcnv(i,j,k,1)
     .                       +hessian(2,2,2)*vcnv(i,j,k,2)*vcnv(i,j,k,2)
     .                       +hessian(2,3,3)*vcnv(i,j,k,3)*vcnv(i,j,k,3)
     .                    +2.*hessian(2,1,2)*vcnv(i,j,k,1)*vcnv(i,j,k,2)
     .                    +2.*hessian(2,1,3)*vcnv(i,j,k,1)*vcnv(i,j,k,3)
     .                    +2.*hessian(2,2,3)*vcnv(i,j,k,2)*vcnv(i,j,k,3)
                else
                  rhs(i,k) = 0d0
                endif
                rhs(i,k) = -dh(dim)
     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
     .             -0.5/jac0*rhs(i,k))/gsuper(dim,dim)
              elseif (dim == 3) then
                if (gamma > 1d0) then
                  rhs(i,j) =  hessian(3,1,1)*vcnv(i,j,k,1)*vcnv(i,j,k,1)
     .                       +hessian(3,2,2)*vcnv(i,j,k,2)*vcnv(i,j,k,2)
     .                       +hessian(3,3,3)*vcnv(i,j,k,3)*vcnv(i,j,k,3)
     .                    +2.*hessian(3,1,2)*vcnv(i,j,k,1)*vcnv(i,j,k,2)
     .                    +2.*hessian(3,1,3)*vcnv(i,j,k,1)*vcnv(i,j,k,3)
     .                    +2.*hessian(3,2,3)*vcnv(i,j,k,2)*vcnv(i,j,k,3)
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

c vectorNeumannBC
c #####################################################################
      subroutine vectorNeumannBC(ivar,array,array0,ieq,dim,loc)
c     -----------------------------------------------------------------
c     Imposes neumann BC for a scalar. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Neumann BC.
c     -----------------------------------------------------------------

      use local_BCS_variables, vec_mod=>vec

      implicit none

c     Call variables

      integer :: ieq,dim,loc,ivar
      real(8) :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,3)
     .          ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1,3)

c     Local variables

      integer :: imax,imin,jmax,jmin,kmax,kmin
      integer :: i,j,k,ig,jg,kg,ip,im,jp,jm,kp,km,ibc,icomp,isgn
     .          ,igp,igm,jgp,jgm,kgp,kgm,i1,j1,k1
      real(8) :: x1,x2,x3,dh(3),jac0,ee,alp,vec(3),cnv(3)
     .          ,cov(3),J0(3),bnorm,jn,etal,B00(3)
      real(8) :: gsuper(3,3),hessian(3,3,3)
     .          ,gsupip(3,3),gsupim(3,3)
     .          ,gsupjp(3,3),gsupjm(3,3)
     .          ,gsupkp(3,3),gsupkm(3,3),gsupavg(3,3)
     .          ,cov00(3),covip(3),covim(3)
     .          ,covjp(3),covjm(3)
     .          ,covkp(3),covkm(3)

      logical :: nonlinear

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      rhs = 0d0

      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      select case (ieq)
      case (IVX,IVY,IVZ) !Velocity components

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

                  if (     coords == 'cyl'
     .                .or. coords == 'hel'
     .                .or. coords == 'tor') then

                    rhs(j,k) = dh(dim)/gsuper(dim,dim)*array(i,j,k,ivar)

                  elseif ( coords == 'car') then

                    rhs(j,k) = 0d0

                  else

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

                    rhs(j,k) = dh(dim)/gsuper(dim,dim)
     .                     *(-gsuper(dim,2)*(array(i,jp,k,ivar)
     .                                      -array(i,jm,k,ivar))/dh(2)
     .                       -gsuper(dim,3)*(array(i,j,kp,ivar)
     .                                      -array(i,j,km,ivar))/dh(3)
     .                       +rhs(j,k) )

                  endif

                endif

              elseif (dim == 2) then

                if (ivar /= dim) then

                  if ( coords == 'car') then

                    rhs(i,k) = 0d0

                  else

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

                    rhs(i,k) = dh(dim)/gsuper(dim,dim)
     .                      *(-gsuper(dim,1)*(array(ip,j,k,ivar)
     .                                       -array(im,j,k,ivar))/dh(1)
     .                        -gsuper(dim,3)*(array(i,j,kp,ivar)
     .                                       -array(i,j,km,ivar))/dh(3)
     .                        +rhs(i,k))
                  endif

                endif

              elseif (dim == 3) then

                if (ivar /= dim) then

                  if ( coords == 'car') then

                    rhs(i,j) = 0d0

                  else

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

                    rhs(i,j) = dh(dim)/gsuper(dim,dim)
     .                      *(-gsuper(dim,1)*(array(ip,j,k,ivar)
     .                                       -array(im,j,k,ivar))/dh(1)
     .                        -gsuper(dim,2)*(array(i,jp,k,ivar)
     .                                       -array(i,jm,k,ivar))/dh(2)
     .                        +rhs(i,j))

                  endif

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

              gsuper = gmetric%grid(igxbc)%gsup(i  ,j,k,:,:)
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

              if (dim == 1) then

                if (ibc == 1) then
                  gsupavg = 0.5*(gsuper+gsupim)

cc                  J0 = 0.5*(array0(i-1,j,k,:)
cc     .                     +array0(i  ,j,k,:))
                elseif (ibc == 2) then
                  gsupavg = 0.5*(gsuper+gsupip)

cc                  J0 = 0.5*(array0(i+1,j,k,:)
cc     .                     +array0(i  ,j,k,:))
                endif

                jn = (array(i,jp,k,3)-array(i,jm,k,3))/dh(2)
     .              -(array(i,j,kp,2)-array(i,j,km,2))/dh(3)

                if (ivar == 2) then
                  icomp = 3
                  isgn  = 1

                  rhs(j,k) = 
     .                (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                   +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                elseif (ivar == 3) then
                  icomp = 2
                  isgn  =-1

                  rhs(j,k) =
     .                (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                   +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                endif

                if (pinch_flow) then     !Finite flow at wall

                  if (ibc == 1) then
                    etal =
     .               0.5*(res(i  ,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc)
     .                   +res(i-1,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc))

                    B00 = 0.5*(array0(i-1,j,k,:)
     .                        +array0(i  ,j,k,:))
                  elseif (ibc == 2) then
                    etal =
     .               0.5*(res(i  ,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc)
     .                   +res(i+1,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc))

                    B00 = 0.5*(array0(i+1,j,k,:)
     .                        +array0(i  ,j,k,:))
                  endif

                  J0 = E0/etal   !Applied current at wall

                  nonlinear = (sum(B00**2) == 0d0)

cc                  write (*,*) 'DIAG -- imposeBC',nonlinear

                  if (nonlinear) then
                    cov = array(i,j,k,:)
                    cnv = matmul(gsuper,cov)
                  else
                    cnv = B00
                    cov =matmul(gmetric%grid(igxbc)%gsub(i,j,k,:,:),cnv)
                  endif

                  bnorm = dot_product(cov,cnv)

                  if (bnorm > 0d0) then
                    alp = dot_product(J0,cnv)/bnorm

                    if (.not.nonlinear) then
                      J0 = J0 - 2*alp*cov  !Correct linearizde current (Bnorm effect)

                      cov = array(i,j,k,:)
                      cnv = matmul(gsuper,cov)
                    endif

                    ee = alp*cnv(icomp)

                    if (.not.nonlinear) then
                      ee = ee + dot_product(J0,cnv)/bnorm*B00(icomp)
                    endif
                  else
                    ee = 0d0
                  endif

                else       !Finite current at wall

                  if (ibc == 1) then
                    J0 = 0.5*(array0(i-1,j,k,:)
     .                       +array0(i  ,j,k,:))
                  elseif (ibc == 2) then
                    J0 = 0.5*(array0(i+1,j,k,:)
     .                       +array0(i  ,j,k,:))
                  endif

                  ee =  dot_product(gsupavg(icomp,:),J0)

                endif

                rhs(j,k)=dh(dim)*(rhs(j,k) + isgn*ee)

              elseif (dim == 2) then

                if (ibc == 3) then
                  gsupavg = 0.5*(gsuper+gsupjm)

                  J0 = 0.5*(array0(i,j-1,k,:)
     .                     +array0(i,j  ,k,:))

                elseif (ibc == 4) then
                  gsupavg = 0.5*(gsuper+gsupjp)

                  J0 = 0.5*(array0(i,j+1,k,:)
     .                     +array0(i,j  ,k,:))
                endif

                jn = (array(i,j,kp,1)-array(i,j,km,1))/dh(3)
     .              -(array(ip,j,k,3)-array(im,j,k,3))/dh(1)

                if (ivar == 3) then
                  icomp = 1
                  isgn  = 1

                  rhs(i,k) =
     .                 (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                  +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                elseif (ivar == 1) then
                  icomp = 3
                  isgn  =-1

                  rhs(i,k) = 
     .                (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                  +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                endif

                ee =  dot_product(gsupavg(icomp,:),J0)

                rhs(i,k) = dh(dim)*(rhs(i,k) + isgn*ee)

              elseif (dim == 3) then

                if (ibc == 5) then
                  gsupavg = 0.5*(gsuper+gsupkm)

                  J0 = 0.5*(array0(i,j,k-1,:)
     .                     +array0(i,j,k  ,:))

                elseif (ibc == 6) then
                  gsupavg = 0.5*(gsuper+gsupkp)

                  J0 = 0.5*(array0(i,j,k+1,:)
     .                     +array0(i,j,k  ,:))
                endif

                jn = (array(ip,j,k,2)-array(im,j,k,2))/dh(1)
     .              -(array(i,jp,k,1)-array(i,jm,k,1))/dh(2)

                if (ivar == 1) then
                  icomp = 2
                  isgn  = 1

                  rhs(i,j) =
     .                 (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                   +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                elseif (ivar == 2) then
                  icomp = 1
                  isgn  =-1

                  rhs(i,j) = 
     .                 (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                   +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                endif

                ee =  dot_product(gsupavg(icomp,:),J0)

                rhs(i,j) = dh(dim)*(rhs(i,j) + isgn*ee)

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
     $                         = array(1,jmin:jmax,kmin:kmax,ivar)
     $                         - rhs(jmin:jmax,kmin:kmax)
      case (2)
        array(nxbc+1,jmin:jmax,kmin:kmax,ivar)
     $                         = array(nxbc,jmin:jmax,kmin:kmax,ivar)
     $                         + rhs(jmin:jmax,kmin:kmax)
      case (3)
        array(imin:imax,0     ,kmin:kmax,ivar)
     $                         = array(imin:imax,1,kmin:kmax,ivar)
     $                         - rhs(imin:imax,kmin:kmax)
      case (4)
        array(imin:imax,nybc+1,kmin:kmax,ivar)
     $                         = array(imin:imax,nybc,kmin:kmax,ivar)
     $                         + rhs(imin:imax,kmin:kmax)
      case (5)
        array(imin:imax,jmin:jmax,0     ,ivar)
     $                         = array(imin:imax,jmin:jmax,1,ivar)
     $                         - rhs(imin:imax,jmin:jmax)
      case (6)
        array(imin:imax,jmin:jmax,nzbc+1,ivar)
     $                         = array(imin:imax,jmin:jmax,nzbc,ivar)
     $                         + rhs(imin:imax,jmin:jmax)
      end select

c     End program

      end subroutine vectorNeumannBC
