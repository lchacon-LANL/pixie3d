c module local_BCS_variables
c####################################################################
      module local_BCS_variables

        use BCS_variables

        use imposeBCinterface, ONLY: extrapolateBC

        use auxiliaryVariables

        use operators

        use icond

        use equilibrium

        use timestepping

        use generalPurposeFunctions

        use transport_params

      contains

c     setBCsched
c     ###################################################################
      subroutine setBCsched

c     -------------------------------------------------------------------
c     Sets BC schedule for dependent and auxiliary variables
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

      gv%bc_seq(5,1) = IVX  !Dep. variable (<0)
      gv%bc_seq(5,2) = 1    !Vector

      gv%bc_seq(6,1) = ITMP !Dep. variable (>0)
      gv%bc_seq(6,2) = 0    !Scalar
#else
      gv%bc_seq(1,1) = IRHO !Dep. Variable (>0)
      gv%bc_seq(1,2) = 0    !Scalar

      gv%bc_seq(2,1) = IBX  !Dep. variable (<0)
      gv%bc_seq(2,2) = 1    !Vector

      gv%bc_seq(3,1) =-IJX  !Aux. variable (<0)
      gv%bc_seq(3,2) = 1    !Vector

      gv%bc_seq(4,1) = IVX  !Dep. variable (<0)
      gv%bc_seq(4,2) = 1    !Vector

      gv%bc_seq(5,1) = ITMP !Dep. variable (>0)
      gv%bc_seq(5,2) = 0    !Scalar
#endif

      endif

c     End program

      end subroutine setBCsched

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
        
      use local_BCS_variables, jac_mod=>jac0,gsup_mod=>gsuper

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

      call findBCLoopLimits(dim,loc,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      select case (ieq)
      case (IVX,IVY,IVZ)
cc      case (IVX,IVY,IVZ,IAX,IAY,IAZ)

        call extrapolateBC(array(:,:,:,ivar),array0(:,:,:,ivar)
cc     .                    ,dim,loc,order)
     .                    ,dim,loc,1)

cc      case (IBX,IBY,IBZ,IJX,IJY,IJZ,IAX,IAY,IAZ) !Imposes divergence-free constraint on B, J, and A
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

      use local_BCS_variables, jac_mod=>jac0,gsup_mod=>gsuper

      implicit none

c     Call variables

      integer    :: ieq,dim,loc
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: i,j,k,ig,jg,kg,ip,im,jp,jm,kp,km,nvar,ibc,icomp
      real(8)    :: x1,x2,x3,dh(3),gsuper(3,3),jac0

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
      subroutine vectorNeumannBC(ivar,array,ieq,dim,loc)
c     -----------------------------------------------------------------
c     Imposes neumann BC for a scalar. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Neumann BC.
c     -----------------------------------------------------------------

      use local_BCS_variables, jac_mod=>jac0,gsup_mod=>gsuper
     .                       ,hess_mod=>hessian, vec_mod=>vec

      implicit none

c     Call variables

      integer    :: ieq,dim,loc,ivar
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,*)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: i,j,k,ig,jg,kg,ip,im,jp,jm,kp,km,ibc,icomp,isgn
      real(8)    :: x1,x2,x3,dh(3),jac0,ee,alp,vec(3),cnv(3)
     .             ,cov(3),J0(3),bnorm
      real(8)    :: jxx,jyy,jzz,jxm,jym,jzm,jxp,jyp,jzp
      real(8)    :: gsuper(3,3),hessian(3,3,3)
     .             ,gsupip(3,3),gsupim(3,3)
     .             ,gsupjp(3,3),gsupjm(3,3)
     .             ,gsupkp(3,3),gsupkm(3,3)
     .             ,cov00(3),covip(3),covim(3)
     .             ,covjp(3),covjm(3)
     .             ,covkp(3),covkm(3)

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      rhs = 0d0

      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      select case (ieq)
cc      case (IVX,IVY,IVZ,IBX,IBY,IBZ) !Velocity components
cc      case (IVX,IVY,IVZ,IJX,IJY,IJZ) !Velocity components
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

cc              if (    equil == 'ppnch' .or. equil == 'ppn3d'
cc     .            .or.equil == 'ppnsl' .or. equil == 'ppnst' ) then

                !THIS WON'T WORK FOR LINEARIZED EQUATIONS; SHOULD USE ACTUAL B
                !   BUT THIS DOESN'T WORK YET.
cc                cov00 =(/array(i,j,k,1),array(i,j,k,2),array(i,j,k,3)/)
cc                covip
cc     .           =(/array(i+1,j,k,1),array(i+1,j,k,2),array(i+1,j,k,3)/)
cc                covim
cc     .           =(/array(i-1,j,k,1),array(i-1,j,k,2),array(i-1,j,k,3)/)
cc                covjp
cc     .           =(/array(i,j+1,k,1),array(i,j+1,k,2),array(i,j+1,k,3)/)
cc                covjm
cc     .           =(/array(i,j-1,k,1),array(i,j-1,k,2),array(i,j-1,k,3)/)
cc                covkp
cc     .           =(/array(i,j,k+1,1),array(i,j,k+1,2),array(i,j,k+1,3)/)
cc                covkm
cc     .           =(/array(i,j,k-1,1),array(i,j,k-1,2),array(i,j,k-1,3)/)

                cov =(/array(i,j,k,1),array(i,j,k,2),array(i,j,k,3)/)
                cnv = matmul(gsuper,cov)

                bnorm = scalarProduct(i,j,k,igxbc,igybc,igzbc
     .                               ,cov(1),cov(2),cov(3)
     .                               ,cnv(1),cnv(2),cnv(3))

                if (dim == 1) then

                  jxx = (array(i,jp,k,3)-array(i,jm,k,3))/dh(2)
     .                 -(array(i,j,kp,2)-array(i,j,km,2))/dh(3)

                  if (ivar == 2) then
                    icomp = 3
                    isgn  = 1

                    if (coords /= 'car') then
                      rhs(j,k) = dh(dim)
     .                 *( (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                     +gsuper(dim,3)*jxx/gsuper(dim,dim))
                    endif

                  elseif (ivar == 3) then
                    icomp = 2
                    isgn  =-1

                    if (coords /= 'car') then
                      rhs(j,k) = dh(dim)
     .                 *( (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                     -gsuper(dim,2)*jxx/gsuper(dim,dim))
                    endif

                  endif

cc                  if (ibc == 1) then
cc                    cov = 0.5*(cov00 + covim)
cc                    cnv = 0.5*(matmul(gsuper,cov00)
cc     .                        +matmul(gsupim,covim))
cc                  elseif (ibc == 2) then
cc                    cov = 0.5*(cov00 + covip)
cc                    cnv = 0.5*(matmul(gsuper,cov00)
cc     .                        +matmul(gsupip,covip))
cc                  endif
cc
cc                  bnorm = scalarProduct(i,j,k,igxbc,igybc,igzbc
cc     .                                 ,cov(1),cov(2),cov(3)
cc     .                                 ,cnv(1),cnv(2),cnv(3))

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

                    if (coords /= 'car') then
                      rhs(i,k) = dh(dim)
     .                 *( (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                    +gsuper(dim,1)*jyy/gsuper(dim,dim))
                    endif
                  elseif (ivar == 1) then
                    icomp = 3
                    isgn  =-1

                    if (coords /= 'car') then
                      rhs(i,k) = dh(dim)
     .                 *( (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                    -gsuper(dim,3)*jyy/gsuper(dim,dim))
                    endif
                  endif

cc                  if (ibc == 3) then
cc                    cov = 0.5*(cov00 + covjm)
cc                    cnv = 0.5*(matmul(gsuper,cov00)
cc     .                        +matmul(gsupjm,covjm))
cc                  elseif (ibc == 4) then
cc                    cov = 0.5*(cov00 + covjp)
cc                    cnv = 0.5*(matmul(gsuper,cov00)
cc     .                        +matmul(gsupjp,covjp))
cc                  endif
cc
cc                  bnorm = scalarProduct(i,j,k,igxbc,igybc,igzbc
cc     .                                 ,cov(1),cov(2),cov(3)
cc     .                                 ,cnv(1),cnv(2),cnv(3))

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

                    if (coords /= 'car') then
                      rhs(i,j) = dh(dim)
     .                 *( (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                    +gsuper(dim,2)*jzz/gsuper(dim,dim))
                    endif
                  elseif (ivar == 2) then
                    icomp = 1
                    isgn  =-1

                    if (coords /= 'car') then
                      rhs(i,j) = dh(dim)
     .                 *( (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                    -gsuper(dim,1)*jzz/gsuper(dim,dim))
                    endif
                  endif

cc                  if (ibc == 5) then
cc                    cov = 0.5*(cov00 + covkm)
cc                    cnv = 0.5*(matmul(gsuper,cov00)
cc     .                        +matmul(gsupkm,covkm))
cc                  elseif (ibc == 6) then
cc                    cov = 0.5*(cov00 + covkp)
cc                    cnv = 0.5*(matmul(gsuper,cov00)
cc     .                        +matmul(gsupkp,covkp))
cc                  endif
cc
cc                  bnorm = scalarProduct(i,j,k,igxbc,igybc,igzbc
cc     .                                 ,cov(1),cov(2),cov(3)
cc     .                                 ,cnv(1),cnv(2),cnv(3))

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

cc              else
cc
cc                if (dim == 1) then
cc
cc                  jxx = (array(i,jp,k,3)-array(i,jm,k,3))/dh(2)
cc     .                 -(array(i,j,kp,2)-array(i,j,km,2))/dh(3)
cc
cc                  if (ivar == 2) then
cc                    icomp = 3
cc                    isgn  = 1
cc
cc                    rhs(j,k) = dh(dim)
cc     .                 *( (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
cc     .                    +gsuper(dim,3)*jxx/gsuper(dim,dim))
cc                  elseif (ivar == 3) then
cc                    icomp = 2
cc                    isgn  =-1
cc
cc                    rhs(j,k) = dh(dim)
cc     .                 *( (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
cc     .                    -gsuper(dim,2)*jxx/gsuper(dim,dim))
cc                  endif
cc
cc                  if (ibc == 1) then
cc                    vec = gsuper(icomp,:)
cc     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
cc                    vec = 0.5*(vec + gsupim(icomp,:)
cc     .                 -gsupim(icomp,dim)/gsupim(dim,dim)*gsupim(dim,:))
cc                    ee = dot_product(vec,J0)
cc
cc                    rhs(j,k) = rhs(j,k) + dh(dim)*isgn*ee
cc                  elseif (ibc == 2) then
cc                    vec = gsuper(icomp,:)
cc     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
cc                    vec = 0.5*(vec + gsupip(icomp,:)
cc     .                 -gsupip(icomp,dim)/gsupip(dim,dim)*gsupip(dim,:))
cc                    ee = dot_product(vec,J0)
cc
cc                    rhs(j,k) = rhs(j,k) + dh(dim)*isgn*ee
cc                  endif
cc
cc                elseif (dim == 2) then
cc
cc                  jyy = (array(i,j,kp,1)-array(i,j,km,1))/dh(3)
cc     .                 -(array(ip,j,k,3)-array(im,j,k,3))/dh(1)
cc
cc                  if (ivar == 3) then
cc                    icomp = 1
cc                    isgn  = 1
cc
cc                    rhs(i,k) = dh(dim)
cc     .                 *( (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
cc     .                    +gsuper(dim,1)*jyy/gsuper(dim,dim))
cc                  elseif (ivar == 1) then
cc                    icomp = 3
cc                    isgn  =-1
cc
cc                    rhs(i,k) = dh(dim)
cc     .                 *( (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
cc     .                    -gsuper(dim,3)*jyy/gsuper(dim,dim))
cc                  endif
cc
cc                  if (ibc == 3) then
cc                    vec = gsuper(icomp,:)
cc     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
cc                    vec = 0.5*(vec + gsupjm(icomp,:)
cc     .                 -gsupjm(icomp,dim)/gsupjm(dim,dim)*gsupjm(dim,:))
cc                    ee = dot_product(vec,J0)
cc
cc                    rhs(i,k) = rhs(i,k) + dh(dim)*isgn*ee
cc                  elseif (ibc == 4) then
cc                    vec = gsuper(icomp,:)
cc     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
cc                    vec = 0.5*(vec + gsupjp(icomp,:)
cc     .                 -gsupjp(icomp,dim)/gsupjp(dim,dim)*gsupjp(dim,:))
cc                    ee = dot_product(vec,J0)
cc
cc                    rhs(i,k) = rhs(i,k) + dh(dim)*isgn*ee
cc                  endif
cc
cc                elseif (dim == 3) then
cc
cc                  jzz = (array(ip,j,k,2)-array(im,j,k,2))/dh(1)
cc     .                 -(array(i,jp,k,1)-array(i,jm,k,1))/dh(2)
cc
cc                  if (ivar == 1) then
cc                    icomp = 2
cc                    isgn  = 1
cc
cc                    rhs(i,j) = dh(dim)
cc     .                 *( (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
cc     .                    +gsuper(dim,2)*jzz/gsuper(dim,dim))
cc                  elseif (ivar == 2) then
cc                    icomp = 1
cc                    isgn  =-1
cc
cc                    rhs(i,j) = dh(dim)
cc     .                 *( (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
cc     .                    -gsuper(dim,1)*jzz/gsuper(dim,dim))
cc                  endif
cc
cc                  if (ibc == 5) then
cc                    vec = gsuper(icomp,:)
cc     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
cc                    vec = 0.5*(vec + gsupkm(icomp,:)
cc     .                 -gsupkm(icomp,dim)/gsupkm(dim,dim)*gsupkm(dim,:))
cc                    ee = dot_product(vec,J0)
cc
cc                    rhs(i,j) = rhs(i,j) + dh(dim)*isgn*ee
cc                  elseif (ibc == 6) then
cc                    vec = gsuper(icomp,:)
cc     .                 -gsuper(icomp,dim)/gsuper(dim,dim)*gsuper(dim,:)
cc                    vec = 0.5*(vec + gsupkp(icomp,:)
cc     .                 -gsupkp(icomp,dim)/gsupkp(dim,dim)*gsupkp(dim,:))
cc                    ee = dot_product(vec,J0)
cc
cc                    rhs(i,j) = rhs(i,j) + dh(dim)*isgn*ee
cc                  endif
cc
cc                endif
cc
cc              endif

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
