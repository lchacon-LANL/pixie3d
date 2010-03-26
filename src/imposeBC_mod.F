c module local_BCS_variables
c####################################################################
      module local_BCS_variables

        use BCS_variables

        use imposeBCinterface, ONLY: extrapolateBC

        use auxiliaryVariables

        use operators

        use timestepping

        use generalPurposeFunctions

        use equilibrium

        use nlfunction_setup, ONLY: res,vis,find_E,ion_hall

        logical :: pinch_flow=.false.,lag_pinch_bc=.false.
     .            ,save_Bn=.false.,save_tn=.false.

#if defined(vec_pot)
        integer :: bc_order=2
#else
cc        integer :: bc_order=1
        integer :: bc_order=2
#endif
        integer,parameter :: VECTOR=1,SCALAR=0,AUX=-1

      contains

c     setBCsched
c     ###################################################################
      subroutine setBCsched(itime)

c     -------------------------------------------------------------------
c     Sets BC schedule for dependent and auxiliary variables, which
c     are identified by a 2-dim integer array, bc_seq. The first
c     dimension is of size nbc_seq, and indicates order of treatment
c     lexicographically (i.e., first component of first dimension is
c     to be treated first, second component of second dimension is
c     to be treated second, and so on).
c
c     The second dimension is of size two, and has the following convention:
c       * First component of the second dimension identifies variable,
c         and its sign indicates whether the variable is dependent (>0) or
c         auxiliary (<0).
c       * Second component of the second dimension identifies scalars (=0)
c         or vectors (=1).
c     -------------------------------------------------------------------

c     Call variables

      integer :: itime

c     Local variables

      integer :: ierr,offset=0

c     Begin program

c     Select offset based on whether we are at the beginning of calculation

      if (itime == 0) then
        offset = 2

#if defined(vec_pot)
        gv%nbc_seq = 18
#else
        gv%nbc_seq = 16
#endif

      elseif (itime > 0) then
        deallocate(gv%bc_seq)
        offset = 0

#if defined(vec_pot)
        gv%nbc_seq = 14
#else
        gv%nbc_seq = 12
#endif

      endif

c     Setup schedule

      if (.not.associated(gv%bc_seq)) then

        allocate(gv%bc_seq(gv%nbc_seq,2),STAT=ierr)

        gv%bc_seq = 0

        if (itime == 0) then
          gv%bc_seq(1,1) = IJ0CNV*AUX       !Aux. variable
          gv%bc_seq(1,2) = VECTOR           !Vector

          gv%bc_seq(2,1) = IJ0COV*AUX       !Aux. variable
          gv%bc_seq(2,2) = VECTOR           !Vector
        endif

#if defined(vec_pot)
        gv%bc_seq(1+offset,1) = IRHO        !Dep. Variable
        gv%bc_seq(1+offset,2) = SCALAR      !Scalar

        gv%bc_seq(2+offset,1) = ITMP        !Dep. variable
        gv%bc_seq(2+offset,2) = SCALAR      !Scalar

        gv%bc_seq(3+offset,1) = IETA*AUX    !Aux. variable
        gv%bc_seq(3+offset,2) = SCALAR      !Scalar

        gv%bc_seq(4+offset,1) = IACNV*AUX   !Aux. variable
        gv%bc_seq(4+offset,2) = VECTOR      !Vector

        gv%bc_seq(5+offset,1) = IAX         !Dep. variable
        gv%bc_seq(5+offset,2) = VECTOR      !Vector

        gv%bc_seq(6+offset,1) = IBCNV*AUX   !Aux. variable
        gv%bc_seq(6+offset,2) = VECTOR      !Vector

        gv%bc_seq(7+offset,1) = IBCOV*AUX   !Aux. variable
        gv%bc_seq(7+offset,2) = VECTOR      !Vector

        gv%bc_seq(8+offset,1)= IJCNV*AUX   !Aux. variable
        gv%bc_seq(8+offset,2)= VECTOR      !Vector

        gv%bc_seq(9+offset,1) = IJCOV*AUX   !Aux. variable
        gv%bc_seq(9+offset,2) = VECTOR      !Vector

        gv%bc_seq(10+offset,1)= IVX         !Dep. variable
        gv%bc_seq(10+offset,2)= VECTOR      !Vector

        gv%bc_seq(11+offset,1)= IVCOV*AUX   !Aux. variable
        gv%bc_seq(11+offset,2)= VECTOR      !Vector

        if (itime == 0) then
          gv%bc_seq(12+offset,1) = IVE0CNV*AUX !Aux. variable
          gv%bc_seq(12+offset,2) = VECTOR      !Vector

          gv%bc_seq(13+offset,1) = IVECOV*AUX !Aux. variable
          gv%bc_seq(13+offset,2) = VECTOR      !Vector

          gv%bc_seq(14+offset,1) = IENI0*AUX   !Aux. variable
          gv%bc_seq(14+offset,2) = VECTOR      !Vector

          gv%bc_seq(15+offset,1) = IENI*AUX !Aux. variable
          gv%bc_seq(15+offset,2) = VECTOR      !Vector

          gv%bc_seq(16+offset,1) = IHETA*AUX   !Aux. variable
          gv%bc_seq(16+offset,2) = SCALAR      !Scalar
        else
          gv%bc_seq(12+offset,1) = IVECOV*AUX  !Aux. variable
          gv%bc_seq(12+offset,2) = VECTOR      !Vector

          gv%bc_seq(13+offset,1) = IVECNV*AUX  !Aux. variable
          gv%bc_seq(13+offset,2) = VECTOR      !Vector

          gv%bc_seq(14+offset,1) = IENI*AUX    !Dep. variable
          gv%bc_seq(14+offset,2) = VECTOR      !Vector
        endif

#else
        gv%bc_seq(1+offset,1) = IRHO        !Dep. Variable
        gv%bc_seq(1+offset,2) = SCALAR      !Scalar

        gv%bc_seq(2+offset,1) = ITMP        !Dep. variable
        gv%bc_seq(2+offset,2) = SCALAR      !Scalar

        gv%bc_seq(3+offset,1) = IETA*AUX    !Aux. variable
        gv%bc_seq(3+offset,2) = SCALAR      !Scalar

        gv%bc_seq(4+offset,1) = IBX         !Dep. variable
        gv%bc_seq(4+offset,2) = VECTOR      !Vector

        gv%bc_seq(5+offset,1) = IBCOV*AUX   !Aux. variable
        gv%bc_seq(5+offset,2) = VECTOR      !Vector

        gv%bc_seq(6+offset,1) = IJCNV*AUX   !Aux. variable
        gv%bc_seq(6+offset,2) = VECTOR      !Vector

        gv%bc_seq(7+offset,1) = IJCOV*AUX   !Aux. variable
        gv%bc_seq(7+offset,2) = VECTOR      !Vector

        gv%bc_seq(8+offset,1) = IVX         !Dep. variable
        gv%bc_seq(8+offset,2) = VECTOR      !Vector

        gv%bc_seq(9+offset,1) = IVCOV*AUX !Aux. variable
        gv%bc_seq(9+offset,2) = VECTOR      !Vector

        if (itime == 0) then
          gv%bc_seq(10+offset,1) = IVE0CNV*AUX   !Aux. variable
          gv%bc_seq(10+offset,2) = VECTOR        !Vector

          gv%bc_seq(11+offset,1) = IVECOV*AUX !Aux. variable
          gv%bc_seq(11+offset,2) = VECTOR        !Vector

          gv%bc_seq(12+offset,1) = IENI0*AUX     !Aux. variable
          gv%bc_seq(12+offset,2) = VECTOR        !Vector

          gv%bc_seq(13+offset,1) = IENI*AUX      !Dep. variable
          gv%bc_seq(13+offset,2) = VECTOR        !Vector
          
          gv%bc_seq(14+offset,1) = IHETA*AUX     !Aux. variable
          gv%bc_seq(14+offset,2) = SCALAR        !Scalar
        else
          gv%bc_seq(10+offset,1) = IVECNV*AUX    !Aux. variable
          gv%bc_seq(10+offset,2) = VECTOR        !Vector

          gv%bc_seq(11+offset,1) = IVECOV*AUX !Aux. variable
          gv%bc_seq(11+offset,2) = VECTOR        !Vector

          gv%bc_seq(12+offset,1) = IENI*AUX      !Dep. variable
          gv%bc_seq(12+offset,2) = VECTOR        !Vector
        endif
#endif

      endif

c     End program

      end subroutine setBCsched

c     ExB_flow_bc
c     #################################################################
      subroutine ExB_flow_bc(igrid,nx,ny,nz,bcnv,prs,eta,v0,dbcov)

c     -----------------------------------------------------------------
c     Finds ExB flow boundary condition when Eo is imposed at boundary
c     (returns contravariant components).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,nx,ny,nz
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1)   :: prs,eta
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3) :: bcnv,v0
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3),optional :: dbcov

c     Local variables

      integer :: i,j,k,ig,jg,kg
      real(8) :: b2,bbcnv(3),bbcov(3),dum(3),grad_p(3),etah,gsuper(3,3)

c     Begin program

      if (.not.pinch_flow) return

      if (igrid /= 1) then
        write (*,*) 'ExB_flow_bc not ready for MG treatment'
        stop
      endif

      i = nx

      do k=1,nz
        do j=1,ny
          bbcnv =0.5*(bcnv(i+1,j,k,:)/gmetric%grid(igrid)%jac(i+1,j,k)
     .               +bcnv(i  ,j,k,:)/gmetric%grid(igrid)%jac(i  ,j,k))
          bbcov =0.5*(matmul(gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .                      ,bcnv(i  ,j,k,:))
     .               +matmul(gmetric%grid(igrid)%gsub(i+1,j,k,:,:)
     .                      ,bcnv(i+1,j,k,:)))

          b2 = dot_product(bbcnv,bbcov)

          if (b2 > 0d0) then
            if (PRESENT(dbcov)) then  !Linearize ExB flow BC
              dum = 0.5*(dbcov(i,j,k,:)+dbcov(i+1,j,k,:))
              bbcov = dum - 2*dot_product(bbcnv,dum)/b2*bbcov
            endif

cc            etah = 0.5*(eta(i+1,j,k)+eta(i,j,k))
            etah = 2d0/(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

            grad_p = grad2(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,1,0,0)

            gsuper = 0.5*(gmetric%grid(igrid)%gsup(i+1,j,k,:,:)
     .                   +gmetric%grid(igrid)%gsup(i  ,j,k,:,:))

            grad_p = matmul(gsuper,grad_p)

            !ExB velocity
            v0(i+1,j,k,1)=(E0(2)*bbcov(3)-E0(3)*bbcov(2)
     .                    -etah*grad_p(1))/b2
          endif
        enddo
      enddo

c     End program

      end subroutine ExB_flow_bc

c     J_pinch_bc
c     #################################################################
      subroutine J_pinch_bc(igrid,nx,ny,nz,bcnv,prs,eta,j0,first_order)

c     -----------------------------------------------------------------
c     Finds current boundary condition when Eo is imposed at boundary
c     (returns covariant current).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,nx,ny,nz
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1)   :: prs,eta
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3) :: bcnv,j0

      logical,optional :: first_order

c     Local variables

      integer :: i,j,k
      real(8) :: b2,bbcnv(3),bbcov(3),grad_p(3),alpha,jach
     .          ,gsub(3,3),ietah
      logical :: fo

c     Begin program

      if (.not.pinch_flow) return

      if (igrid /= 1) then
        write (*,*) 'J_pinch_bc not ready for MG treatment'
        stop
      endif

      if (PRESENT(first_order)) then
        fo = first_order
      else
        fo = .false.
cc        fo = .true.
      endif

      i = nx

      if (.not.fo) then
        do k=1,nz
          do j=1,ny
cc            jach  = 0.5*(gmetric%grid(igrid)%jac(i  ,j,k)
cc     .                  +gmetric%grid(igrid)%jac(i+1,j,k))
cc            bbcnv = 0.5*(bcnv(i+1,j,k,:)+bcnv(i,j,k,:))/jach
            bbcnv = 0.5*
     .             (bcnv(i+1,j,k,:)/gmetric%grid(igrid)%jac(i+1,j,k)
     .             +bcnv(i  ,j,k,:)/gmetric%grid(igrid)%jac(i  ,j,k))
            bbcov = 0.5*(matmul(gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .                         ,bcnv(i  ,j,k,:))
     .                  +matmul(gmetric%grid(igrid)%gsub(i+1,j,k,:,:)
     .                         ,bcnv(i+1,j,k,:)))

            b2 = dot_product(bbcov,bbcnv)

            if (b2 > 0d0) then
cc              ietah = 2./(eta(i+1,j,k)+eta(i,j,k))
              ietah =0.5*(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

              grad_p= grad2(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,1,0,0)

              alpha = dot_product(E0,bbcnv)*ietah

              j0(i+1,j,k,:)=(alpha*bbcnv+cross_product(bbcov,grad_p))/b2
            endif
          enddo
        enddo
      else
        do k=1,nz
          do j=1,ny
cc            jach = gmetric%grid(igrid)%jac(i,j,k)
            jach  = 0.5*(gmetric%grid(igrid)%jac(i  ,j,k)
     .                  +gmetric%grid(igrid)%jac(i+1,j,k))
            bbcnv = bcnv(i,j,k,:)/jach

cc            gsub = gmetric%grid(igrid)%gsub(i,j,k,:,:)
            gsub = 0.5*(gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .                 +gmetric%grid(igrid)%gsub(i+1,j,k,:,:))
            bbcov = matmul(gsub,bcnv(i,j,k,:))

            b2 = dot_product(bbcov,bbcnv)
            if (b2 > 0d0) then
cc              ietah = 2./(eta(i+1,j,k)+eta(i,j,k))
              ietah =0.5*(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

              grad_p = grad2(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,1,0,0)

              alpha = dot_product(E0,bbcnv)*ietah

              j0(i+1,j,k,:)=(alpha*bbcnv+cross_product(bbcov,grad_p))/b2
            endif
          enddo
        enddo
      endif

c     End program

      end subroutine J_pinch_bc

      end module local_BCS_variables

c setup_app_BC
c #####################################################################
      subroutine setup_app_BC(it)
c ---------------------------------------------------------------------
c     Sets up application BC
c ---------------------------------------------------------------------

      use local_BCS_variables

      implicit none

c Call variables

      integer :: it

c Local variables

c Begin program

c Set BC schedule

      call setBCsched(it)

c Allocate auxiliary variables

#if defined(samrai)
      call allocAuxVariables
#endif

c End program

      end subroutine setup_app_BC

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

      call extrapolateBC(array,array0,dim,loc,order)
cc      call extrapolateBC(array,array0,dim,loc,1)

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

      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      select case (ieq)
      case (IVX,IVY,IVZ)

        call extrapolateBC(array(:,:,:,ivar),array0(:,:,:,ivar)
     .                    ,dim,loc,order)
cc     .                    ,dim,loc,1)

      case (IBX,IBY,IBZ,IJX,IJY,IJZ) !Imposes divergence-free constraint on B, J

        if (ivar /= dim) then

          call extrapolateBC(array (:,:,:,ivar),array0(:,:,:,ivar)
     .                      ,dim,loc,order)
cc     .                      ,dim,loc,1)
cc     .                      ,dim,loc,2)

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
cc     .                      ,dim,loc,1)

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
      case (IRHO,ITMP)
cc      case (IRHO)

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
cc
cc      case (ITMP)
cc
cc        do i=imin,imax
cc          do j=jmin,jmax
cc            do k=kmin,kmax
cc
cc              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
cc
cc              gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)
cc              jac0   = gmetric%grid(igxbc)%jac (i,j,k)
cc
cc              hessian = gmetric%grid(igxbc)%Gamma(i,j,k,:,:,:)
cc
cc              ip = min(i+1,nxbc)
cc              im = max(i-1,1   )
cc              jp = min(j+1,nybc)
cc              jm = max(j-1,1   )
cc              kp = min(k+1,nzbc)
cc              km = max(k-1,1   )
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
cc                if (gamma > 1d0) then
cc                  rhs(j,k) =  hessian(1,1,1)*vcnv(i,j,k,1)*vcnv(i,j,k,1)
cc     .                       +hessian(1,2,2)*vcnv(i,j,k,2)*vcnv(i,j,k,2)
cc     .                       +hessian(1,3,3)*vcnv(i,j,k,3)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(1,1,2)*vcnv(i,j,k,1)*vcnv(i,j,k,2)
cc     .                    +2.*hessian(1,1,3)*vcnv(i,j,k,1)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(1,2,3)*vcnv(i,j,k,2)*vcnv(i,j,k,3)
cc                else
cc                  rhs(j,k) = 0d0
cc                endif
cc                rhs(j,k) = -dh(dim)
cc     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
cc     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
cc     .             -0.5/jac0*rhs(j,k))/gsuper(dim,dim)
cc              elseif (dim == 2) then
cc                if (gamma > 1d0) then
cc                  rhs(i,k) =  hessian(2,1,1)*vcnv(i,j,k,1)*vcnv(i,j,k,1)
cc     .                       +hessian(2,2,2)*vcnv(i,j,k,2)*vcnv(i,j,k,2)
cc     .                       +hessian(2,3,3)*vcnv(i,j,k,3)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(2,1,2)*vcnv(i,j,k,1)*vcnv(i,j,k,2)
cc     .                    +2.*hessian(2,1,3)*vcnv(i,j,k,1)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(2,2,3)*vcnv(i,j,k,2)*vcnv(i,j,k,3)
cc                else
cc                  rhs(i,k) = 0d0
cc                endif
cc                rhs(i,k) = -dh(dim)
cc     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
cc     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
cc     .             -0.5/jac0*rhs(i,k))/gsuper(dim,dim)
cc              elseif (dim == 3) then
cc                if (gamma > 1d0) then
cc                  rhs(i,j) =  hessian(3,1,1)*vcnv(i,j,k,1)*vcnv(i,j,k,1)
cc     .                       +hessian(3,2,2)*vcnv(i,j,k,2)*vcnv(i,j,k,2)
cc     .                       +hessian(3,3,3)*vcnv(i,j,k,3)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(3,1,2)*vcnv(i,j,k,1)*vcnv(i,j,k,2)
cc     .                    +2.*hessian(3,1,3)*vcnv(i,j,k,1)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(3,2,3)*vcnv(i,j,k,2)*vcnv(i,j,k,3)
cc                else
cc                  rhs(i,j) = 0d0
cc                endif
cc                rhs(i,j) = -dh(dim)
cc     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
cc     .              +gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
cc     .             -0.5/jac0*rhs(i,j))/gsuper(dim,dim)
cc              endif
cc
cc            enddo
cc          enddo
cc        enddo

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

cc      real(8),pointer,dimension(:,:,:) :: eta_bc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      rhs = 0d0  !Default

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
                elseif (ibc == 2) then
                  gsupavg = 0.5*(gsuper+gsupip)
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

                !!!!!!!Finite flow at wall!!!!!!!
                if (pinch_flow) then

                  ee = array0(i+1,j,k,icomp)

                !!!!!!!Finite current at wall!!!!!!!
                else
                  if (ibc == 1) then
                    ee = 0.5*(array0(i-1,j,k,icomp)
     .                       +array0(i  ,j,k,icomp))
                  elseif (ibc == 2) then
                    ee = 0.5*(array0(i+1,j,k,icomp)
     .                       +array0(i  ,j,k,icomp))
                  endif
                endif

                rhs(j,k)=dh(dim)*(rhs(j,k) + isgn*ee)

              elseif (dim == 2) then

                if (ibc == 3) then
                  gsupavg = 0.5*(gsuper+gsupjm)

                  ee = 0.5*(array0(i,j-1,k,icomp)
     .                     +array0(i,j  ,k,icomp))
                elseif (ibc == 4) then
                  gsupavg = 0.5*(gsuper+gsupjp)

                  ee = 0.5*(array0(i,j+1,k,icomp)
     .                     +array0(i,j  ,k,icomp))
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

                rhs(i,k) = dh(dim)*(rhs(i,k) + isgn*ee)

              elseif (dim == 3) then

                if (ibc == 5) then
                  gsupavg = 0.5*(gsuper+gsupkm)

                  ee = 0.5*(array0(i,j,k-1,icomp)
     .                     +array0(i,j,k  ,icomp))
                elseif (ibc == 6) then
                  gsupavg = 0.5*(gsuper+gsupkp)

                  ee = 0.5*(array0(i,j,k+1,icomp)
     .                     +array0(i,j,k  ,icomp))
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

                rhs(i,j) = dh(dim)*(rhs(i,j) + isgn*ee)

              endif

            enddo
          enddo
        enddo

      case default

        write (*,*) 'Variable=',ieq
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
