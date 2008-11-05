
c module precond_variables
c ######################################################################
      module precond_variables

        use timeStepping

        use precond_setup

        use transport_params

        use equilibrium

        use constants

        use auxiliaryVariables

        use operators

        use variables

        use grid

        use mgarraySetup

        use newtongm

        use grid_debug

        integer :: i,j,k,ig,jg,kg,ieq
        real(8) :: x1,y1,z1
        logical :: cartsn,covariant,to_cartsn,to_cnv

        integer :: ntotd2p,si_it,ores=2,oprol=2
        real(8) :: omeg_smth
        character(2) :: smooth
        logical :: gm_smooth,mg_glrkin,si_car,gm_coarse
     .            ,line_relax,inv_B,divCln=.false.,v_solve

cc        real(8), allocatable, dimension(:,:,:,:):: p0

        real(8), allocatable, dimension(:,:)    :: rho_diag,tmp_diag
     .                                            ,b_diag,v_diag,v_diag2
     .                                            ,a_diag,b_diag_hyp

        real(8), allocatable, dimension(:,:),target :: d_pois

        real(8), allocatable, dimension(:,:)    :: mgj0cnv,mgadvdiffV0

        real(8), allocatable, dimension(:,:,:)  :: mgnablaV0

        real(8), allocatable, dimension(:)      :: mgdivV0

        integer, allocatable, dimension(:,:)    :: bcs,bcs_hyp

        type(var_array),pointer :: varray

        type (mg_array ),target :: gp0,gb0,gv0,ga0,grho0,gb0_cov,gv0_cov

#if !defined(vec_pot)
        real(8),pointer,dimension(:,:,:):: rho,px,py,pz,bx,by,bz,tmp
#else
        real(8),pointer,dimension(:,:,:):: rho,px,py,pz,ax,ay,az,tmp

#endif

      contains

c     allocPrecVariables
c     ###################################################################
      subroutine allocPrecVariables

c     -------------------------------------------------------------------
c     Allocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer    :: alloc_stat

c     Begin program

        ntotd2p = 2*ntotd/neqd

        allocate (mgj0cnv    (ntotd2p,3))
        allocate (mgadvdiffV0(ntotd2p,3))
        allocate (mgnablaV0  (ntotd2p,3,3))
        allocate (mgdivV0    (ntotd2p))

        allocate (bcs(6,neqd+3),bcs_hyp(6,6))

        allocate (rho_diag  (1,  ntotd2p)
     .           ,tmp_diag  (1,  ntotd2p)
#if !defined(vec_pot)
     .           ,  b_diag    (3,3*ntotd2p)
     .           ,  b_diag_hyp(6,6*ntotd2p)
#else
     .           ,  a_diag  (3,3*ntotd2p)
#endif
     .           ,  v_diag  (3,3*ntotd2p)
     .           ,  v_diag2 (3,3*ntotd2p),STAT=alloc_stat)

        call allocateMGArray(1,gp0)
        call allocateMGArray(1,grho0)
        call allocateMGArray(3,gv0)
        call allocateMGArray(3,gv0_cov)
        call allocateMGArray(3,gb0)
        call allocateMGArray(3,gb0_cov)
        call allocateMGArray(3,ga0)   !B + di curl(v)

c     End program

      end subroutine allocPrecVariables

c     deallocPrecVariables
c     ###################################################################
      subroutine deallocPrecVariables

c     -------------------------------------------------------------------
c     Deallocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

        deallocate (mgj0cnv)
        deallocate (mgdivV0)
        deallocate (mgadvdiffV0)
        deallocate (mgnablaV0)

        deallocate (bcs,bcs_hyp)

        deallocate (rho_diag,tmp_diag  
#if !defined(vec_pot)
     .             ,b_diag,b_diag_hyp
#else
     .             ,a_diag 
#endif
     .             ,v_diag,v_diag2)

        call deallocateMGArray(gp0)
        call deallocateMGArray(grho0)
        call deallocateMGArray(gb0)
        call deallocateMGArray(gb0_cov)
        call deallocateMGArray(gv0)
        call deallocateMGArray(gv0_cov)
        call deallocateMGArray(ga0)

        call deallocateDerivedType(varray)

c     End program

      end subroutine deallocPrecVariables

c     gather
c     ###################################################################
      subroutine gather(xout,xin)

c     -------------------------------------------------------------------
c     Gathers magnitudes from different arrays to vector.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        real(8),intent(IN),dimension(:,:) :: xin
        real(8),intent(OUT),dimension(:)  :: xout

c     Local variables

        integer    :: i,j,k,ii,iii,ieq
        real(8)    :: dvol

c     Begin program

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              do ieq=1,neqd
                iii = ieq + neqd*(ii-1)
                xout(iii) = xin(ii,ieq)
              enddo
            enddo
          enddo
        enddo

c     End program

      end subroutine gather

c     scatter
c     ###################################################################
      subroutine scatter(xout,xin)

c     -------------------------------------------------------------------
c     Scatters magnitudes in vector to different arrays indexed by the
c     variable identifier ieq.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        real(8),intent(IN)  :: xin(:)
        real(8),intent(OUT) :: xout(:,:)

c     Local variables

        integer    :: i,j,k,ii,iii,ieq

c     Begin program

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              do ieq=1,neqd
                iii = ieq + neqd*(ii-1)
                xout(ii,ieq) = xin(iii)
              enddo

            enddo
          enddo
        enddo

c     End program

      end subroutine scatter

c     XtoCar
c     ###################################################################
      subroutine XtoCar(igrid,x,xout,cov)

c     -------------------------------------------------------------------
c     Transforms vector to Cartesian components
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: igrid
        real(8),intent(INOUT),dimension(:,:) :: x
        real(8),intent(OUT),optional,dimension(:,:) :: xout
        logical, optional :: cov

c     Local variables

        integer    :: i,j,k,ii,nx,ny,nz
        real(8)    :: cnv(3)
        logical    :: covar

c     Begin program

        if (PRESENT(cov)) then
          covar = cov
        else
          covar = .false.
        endif

        if (size(x,2) /= 3) then
          call pstop('XtoCar','Vector has wrong dimensions')
        endif

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        if (PRESENT(xout)) then
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                cnv(:) = x(ii,:)

                call transformVectorToCartesian
     .               (i,j,k,igrid,igrid,igrid
     .               ,cnv(1),cnv(2),cnv(3),covar
     .               ,xout(ii,1),xout(ii,2),xout(ii,3))

              enddo
            enddo
          enddo
        else
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                cnv(:) = x(ii,:)

                call transformVectorToCartesian
     .               (i,j,k,igrid,igrid,igrid
     .               ,cnv(1),cnv(2),cnv(3),covar
     .               ,x(ii,1),x(ii,2),x(ii,3))

              enddo
            enddo
          enddo
        endif

c     End program

      end subroutine XtoCar

c     XtoCurv
c     ###################################################################
      subroutine XtoCurv(igrid,x,xout,cov)

c     -------------------------------------------------------------------
c     Transform vector components from Cartesian to curvilinear
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: igrid
        real(8),intent(INOUT),dimension(:,:) :: x
        real(8),intent(OUT),optional,dimension(:,:) :: xout
        logical, optional :: cov

c     Local variables

        integer    :: i,j,k,ii,nx,ny,nz
        real(8)    :: car(3)
        logical    :: covar

c     Begin program

        if (PRESENT(cov)) then
          covar = cov
        else
          covar = .false.
        endif

        if (size(x,2) /= 3) then
          call pstop('XtoCar','Vector has wrong dimensions')
        endif

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        if (PRESENT(xout)) then
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                call transformVectorToCurvilinear
     .               (i,j,k,igrid,igrid,igrid
     .               ,x(ii,1),x(ii,2),x(ii,3)
     .               ,covar
     .               ,xout(ii,1),xout(ii,2),xout(ii,3))
              enddo
            enddo
          enddo
        else
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                car = x(ii,:)

                call transformVectorToCurvilinear
     .               (i,j,k,igrid,igrid,igrid
     .               ,car(1),car(2),car(3)
     .               ,covar
     .               ,x(ii,1),x(ii,2),x(ii,3))
              enddo
            enddo
          enddo
        endif

c     End program

      end subroutine XtoCurv

c     fillArray
c     ###################################################################
      subroutine fillArray(igrid,icomp,neq,x,arr,bcs,iorder,arr_cov
     .                    ,is_vec,is_cnv)

c     -------------------------------------------------------------------
c     Fills array arr (with BCs) using vector x.
c     -------------------------------------------------------------------

        use setMGBC_interface

        implicit none

c     Call variables

        integer    :: igrid,icomp,neq,iorder,bcs(:,:)
        real(8),dimension(:,:) :: x
        real(8),dimension(:,:,:,:) :: arr
        real(8),optional,dimension(:,:,:,:) :: arr_cov
        logical,optional :: is_vec,is_cnv

c     Local variables

        integer    :: ieq
        logical    :: isvec,iscnv

c     Begin program

        if (PRESENT(is_vec)) then
          isvec = is_vec
        else
          isvec = (neq > 1)
        endif

        if (PRESENT(is_cnv)) then
          iscnv = is_cnv
        else
          iscnv = .true.
        endif

        do ieq=1,neq
          call mapMGVectorToArray(0,1,x(:,ieq),nx,ny,nz
     .                       ,arr(:,:,:,ieq),igrid,.false.)
        enddo

        if (PRESENT(arr_cov)) then
          call setMGBC(0,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=icomp,is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder,arr_cov=arr_cov)
        else
          call setMGBC(0,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=icomp,is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder)
        endif

c     End program

      end subroutine fillArray

c     findPCcoeffs
c     ###################################################################
      subroutine findPCcoeffs

c     -------------------------------------------------------------------
c     Finds coefficients for linearized systems in preconditioner
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer    :: order,nxx,nyy,nzz,igrid,ii,ivar,igr
     .               ,bc_pc(size(bcs,1),size(bcs,2))
     .               ,i,j,k,ip,im,jp,jm,kp,km,nxg,nyg,nzg

        real(8)    :: dvol,veclap(3),cnv(3),cov(3),dxx,dyy,dzz

        real(8), allocatable, dimension(:,:,:,:) :: vector
        real(8), allocatable, dimension(:,:,:,:,:) :: tensor

c     Debug

        real(8),allocatable,dimension(:,:,:) :: debug2
        character(13) :: debugfile

c     Begin program

        nxx = grid_params%nxv(igx)
        nyy = grid_params%nyv(igy)
        nzz = grid_params%nzv(igz)

        nx  = nxx
        ny  = nyy
        nz  = nzz

        igrid = igx

        order = 0

c     Setup BCs

        !Set boundary conditions
        do ieq=1,neqd
          bc_pc(:,ieq) = varray%array_var(ieq)%bconds(:)
        enddo

#if !defined(vec_pot)
        bc_pc(:,IJX:IJZ) = bc_pc(:,IBX:IBZ)
        where (bc_pc(:,IJX:IJZ) == -NEU) bc_pc(:,IJX:IJZ) = -DIR !Use covariant components 
#else
        where (bc_pc(:,IAX:IAZ) == -EQU) bc_pc(:,IAX:IAZ) = -DIR !Use covariant components 

        !Magnetic field conditions
        bc_pc(:,IBX:IBZ) = bc_pc(:,IAX:IAZ)
        where (bc_pc(:,IBX:IBZ) == -DIR) bc_pc(:,IBX:IBZ) = -NEU !Use covariant components 
#endif

        call setASM_PC_BC(igx,igy,igz,bc_pc,IIBC=EXT) !THIS DOES NOT COMMUNICATE COEFF.

c     Store density in all grids (w/o BCs)

        grho0%grid(igrid)%array(:,:,:,1) = rho

        call restrictMGArray(IRHO,1,grho0,bc_pc(:,IRHO),igrid,order)

c     Store magnetic field cnv components in all grids (w/ BCs)

cc        call restrict_B(.true.)

        gb0%grid(igrid)%array = bcnv

        call restrictMGArray(IBX,3,gb0,bc_pc(:,IBX:IBZ),igrid,order)

c     Store magnetic field cov components in all grids (w/ BCs)

        gb0_cov%grid(igrid)%array = bcov

        call restrictMGArray(IBX,3,gb0_cov,bc_pc(:,IBX:IBZ),igrid,order
     .                      ,iscnv=.false.)

c     Store ion velocity components in all grids (w/ BCs)

        gv0%grid(igrid)%array = vcnv

        call restrictMGArray(IVX,3,gv0,bc_pc(:,IVX:IVZ),igrid,order)

c     Store magnetic field cov components in all grids (w/ BCs)

        if (di > 0d0) then
          gv0_cov%grid(igrid)%array = vcov

          call restrictMGArray(IVX,3,gv0_cov,bc_pc(:,IVX:IVZ),igrid
     .                        ,order,iscnv=.false.)
        endif

c     Store (B + di curl(v) in all grids

        if (di > 0d0) then

          bc_pc(:,IJX:IJZ) = bc_pc(:,IBX:IBZ)
          where (bc_pc(:,IJX:IJZ) == -NEU) bc_pc(:,IJX:IJZ) = EXT !Use contravariant components 
cc          bc_pc(:,IJX) = bcond
cc          bc_pc(:,IJY) = bcond
cc          bc_pc(:,IJY) = bcond
cc
cc          where (bc_pc(:,IJX:IJZ) == DEF) bc_pc(:,IJX:IJZ) = EXT
cc
          do k=0,nz+1
            do j=0,ny+1
              do i=0,nx+1

                do ieq=1,3
                  cnv(ieq) = curl2(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                           ,gv0_cov%grid(igrid)%array(:,:,:,1)
     .                           ,gv0_cov%grid(igrid)%array(:,:,:,2)
     .                           ,gv0_cov%grid(igrid)%array(:,:,:,3)
     .                           ,ieq)
                enddo

                ga0%grid(igrid)%array(i,j,k,:)
     .                 = gb0%grid(igrid)%array(i,j,k,:) + di*cnv
              enddo
            enddo
          enddo

          call restrictMGArray(IJX,3,ga0,bc_pc(:,IJX:IJZ),igrid,order)

cc          do igr=1,grid_params%ngrid
cc            nxg = grid_params%nxv(igr)
cc            nyg = grid_params%nyv(igr)
cc            nzg = grid_params%nzv(igr)
cc
cc            do k=0,nzg+1
cc              do j=0,nyg+1
cc                do i=0,nxg+1
cc
cc                  do ieq=1,3
cc                    cnv(ieq) = curl2(i,j,k,nxg,nyg,nzg,igr,igr,igr
cc     .                           ,gv0_cov%grid(igr)%array(:,:,:,1)
cc     .                           ,gv0_cov%grid(igr)%array(:,:,:,2)
cc     .                           ,gv0_cov%grid(igr)%array(:,:,:,3)
cc     .                           ,ieq)
cc                  enddo
cc
cc                  ga0%grid(igr)%array(i,j,k,:)
cc     .                 = gb0%grid(igr)%array(i,j,k,:) + di*cnv
cc                enddo
cc              enddo
cc            enddo
cc
cc          enddo

        else

          call equateMGArray(ga0,gb0)

        endif

c     Find auxiliary quantities and store them in all grids

        !Velocity divergence (w/o BCs)
        allocate(divrgV(0:nxx+1,0:nyy+1,0:nzz+1))

        do k=1,nzz
          do j=1,nyy
            do i=1,nxx
              divrgV(i,j,k) = div(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                           ,vcnv(:,:,:,1)
     .                           ,vcnv(:,:,:,2)
     .                           ,vcnv(:,:,:,3))
            enddo
          enddo
        enddo

        call restrictArrayToMGVector(1,nxx,nyy,nzz,divrgV,mgdivV0,igrid
     .                            ,order,.false.)

        deallocate(divrgV)

        !pressure (w/ BCs)
        gp0%grid(igrid)%array(:,:,:,1) = a_p*rho*tmp

        call restrictMGArray(IRHO,1,gp0,bc_pc(:,IRHO),igrid,order)

        !v0/dt+theta(v0.grad(v0)-mu*veclap(v0)) and nabla_v0
        allocate(vector(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .          ,tensor(0:nxx+1,0:nyy+1,0:nzz+1,3,3))

        do k = 1,nzz
          do j = 1,nyy
            do i = 1,nxx

              ip = i+1
              im = i-1
              jp = j+1
              jm = j-1
              kp = k+1
              km = k-1

              ii  = i + nxx*(j-1) + nxx*nyy*(k-1)

              jac    = gmetric%grid(igrid)%jac(i,j,k)

              nabla_v= fnabla_v(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                         ,vcnv(:,:,:,1)
     .                         ,vcnv(:,:,:,2)
     .                         ,vcnv(:,:,:,3),0)
              tensor(i,j,k,:,:) = nabla_v

              if (nc_eom_v) then
                veclap = veclaplacian(i,j,k,nxx,nyy,nzz
     .                               ,igrid,igrid,igrid
     .                               ,vcnv,alt_eom,nuu,vol=.false.)

                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

                dxx = dxh(ig)
                dyy = dyh(jg)
                dzz = dzh(kg)

                !Grad(p)
                cov(1) = 0.5*( rho(i,j,k)*(tmp(ip,j,k)-tmp(im,j,k))/dxx
     .                        +tmp(i,j,k)*(rho(ip,j,k)-rho(im,j,k))/dxx)

                cov(2) = 0.5*( rho(i,j,k)*(tmp(i,jp,k)-tmp(i,jm,k))/dyy
     .                        +tmp(i,j,k)*(rho(i,jp,k)-rho(i,jm,k))/dyy)
        
                cov(3) = 0.5*( rho(i,j,k)*(tmp(i,j,kp)-tmp(i,j,km))/dzz
     .                        +tmp(i,j,k)*(rho(i,j,kp)-rho(i,j,km))/dzz)

                cov = cov*a_p    !Multiply by alpha_p=1 + Ti/Te

                call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                                      ,cov(1),cov(2),cov(3)
     .                                      ,cnv(1),cnv(2),cnv(3)
     .                                      ,.true.)

                !JxB
                cnv(1) = jcov(i,j,k,2)*bcov(i,j,k,3)
     .                 - jcov(i,j,k,3)*bcov(i,j,k,2) - cnv(1)

                cnv(2) = jcov(i,j,k,3)*bcov(i,j,k,1)
     .                 - jcov(i,j,k,1)*bcov(i,j,k,3) - cnv(2)

                cnv(3) = jcov(i,j,k,1)*bcov(i,j,k,2)
     .                 - jcov(i,j,k,2)*bcov(i,j,k,1) - cnv(3)

                vector(i,j,k,:) = (cnv-veclap)/rho(i,j,k)**2
              else
                do ivar=1,3
                  vector(i,j,k,ivar) = 
     .                          cnp*gv0%grid(igrid)%array(i,j,k,ivar)/dt
     .                        + alpha*vcnv(i,j,k,1)*nabla_v(1,ivar)/jac
     .                        + alpha*vcnv(i,j,k,2)*nabla_v(2,ivar)/jac
     .                        + alpha*vcnv(i,j,k,3)*nabla_v(3,ivar)/jac
                enddo
              endif

            enddo
          enddo
        enddo

        do ivar=1,3
          call restrictArrayToMGVector(1,nxx,nyy,nzz
     .                              ,vector(:,:,:,ivar)
     .                              ,mgadvdiffV0(:,ivar)
     .                              ,igrid,order,.false.)
        enddo

        do i=1,3
          do j=1,3
            call restrictArrayToMGVector(1,nxx,nyy,nzz
     .                              ,tensor(:,:,:,i,j)
     .                              ,mgnablaV0(:,i,j)
     .                              ,igrid,order,.false.)
          enddo
        enddo

        deallocate(vector,tensor)

c diag ****
cc      write (*,*) 'DIAG -- findCoeffs'
cc
cc      debugfile = 'debug_p'//trim(int2char(my_rank))//'.bin'
cc      open(unit=110,file=debugfile,form='unformatted'
cc     .    ,status='replace')
cc
cc      igx = 2
cc      k = 2
cc
cc      nx = grid_params%nxv(igx)
cc      ny = grid_params%nyv(igx)
cc      nz = grid_params%nzv(igx)
cc      allocate(debug2(0:nx+1,0:ny+1,0:nz+1))
cc
cc      do ieq=1,3
cc        debug2 = gb0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,lxmin,lxmax,lymin,lymax,ieq-1,110)
cc      enddo
cc
cc      do ieq=1,3
cc        debug2 = gb0_cov%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,lxmin,lxmax,lymin,lymax,1,110)
cc      enddo
cc
cc      do ieq=1,3
cc        debug2 = gv0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,lxmin,lxmax,lymin,lymax,1,110)
cc      enddo
cc
cc      debug2 = grho0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,1)
cc      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .            ,lxmin,lxmax,lymin,lymax,1,110)
cc
cc      debug2 = gp0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,1)
cc      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .            ,lxmin,lxmax,lymin,lymax,1,110)
cc
cc      deallocate(debug2)
cc      close(110)
cc
cc      call PetscEnd(mpierr)
cc      stop
c diag ****

c     End program

      contains

c     restrict_B
c     ###############################################################
cc      subroutine restrictCurvVec(icomp,ax,ay,az,gvec,order,cartesian)
      subroutine restrict_B(cartesian)

c     ---------------------------------------------------------------
c     This routines restricts B Cartesian components, and then finds
c     corresponding cov, cnv components with BCs.
c     ---------------------------------------------------------------

      use setMGBC_interface

      implicit none

c     Call variables

cc      integer    :: order,icomp
cc      type(mg_array)  :: gvec
      logical    :: cartesian

c     Local variables

      integer    :: igr,bcmod(6,3),nxc,nyc,nzc,neq

c     Begin program

      neq = 3

      igr = 1

      nxc = grid_params%nxv(igr)
      nyc = grid_params%nyv(igr)
      nzc = grid_params%nzv(igr)

c     Find cartesian coordinates

      gb0%grid(1)%array = bcnv
cc      gb0%grid(1)%array(:,:,:,1) = bx
cc      gb0%grid(1)%array(:,:,:,2) = by
cc      gb0%grid(1)%array(:,:,:,3) = bz

      if (cartesian)
     .   call transformVector(igr,igr,igr,0,nxc+1,0,nyc+1,0,nzc+1
     .                       ,gb0%grid(1)%array(:,:,:,1)
     .                       ,gb0%grid(1)%array(:,:,:,2)
     .                       ,gb0%grid(1)%array(:,:,:,3)
     .                       ,'cnv','car')

c     Restrict components (w/o BCs)

      call restrictMGArray(0,neq,gb0,bcs(:,IBX:IBZ),igr,order)
cc        call restrictMGArray(IBX,3,gb0,bcs(:,IBX:IBZ),igrid,order
cc     .                      ,iscnv=.true.,isvec=.not.cartesian
cc     .                      ,result_is_vec=.true.)

c     Impose BCs, and return cov, cnv components

      bcmod = bcs(:,IBX:IBZ)
      where (bcmod == EQU)
        bcmod = EXT
      end where

      do igr=1,grid_params%ngrid
        nxc = grid_params%nxv(igr)
        nyc = grid_params%nyv(igr)
        nzc = grid_params%nzv(igr)

        call setMGBC(0,neq,nxc,nyc,nzc,igr
     .              ,gb0%grid(igr)%array,bcmod       !Input/output cnv array
     .              ,arr_cov=gb0_cov%grid(igr)%array !Input/output cov array
     .              ,icomp=IBX
     .              ,is_cnv=.true.          !This informs that input is in cnv array
     .              ,is_vec=.not.cartesian  !This informs whether input is Cartesian
     .              ,result_is_vec=.true.   !This returns curvilinear components
     .              ,iorder=order)
      enddo

c     End program

      end subroutine restrict_B

      end subroutine findPCcoeffs

c     cSolver
c #   ##################################################################
      subroutine cSolver(neq,ntotp,b,x,bcnd,igrid,out,guess
     .                  ,matvec,dg,ncolors,line_relax,gm_coarse_solve
     .                  ,gm_driver,cvrg_tst,order_res,order_prol,smooth
     .                  ,omega,mvcyc,mu_cycle,galerkin)
c     ---------------------------------------------------------------
c     This subroutine solves a coupled system of neq equations. 
c     In call sequence:
c       * neq: number of coupled equations
c       * ntotp: number of mesh points
c       * b: rhs
c       * x: solution
c       * bcnd: boundary condition defs.
c       * igrid: MG grid level (igrid=1 is finest level)
c       * out: level of output information
c       * guess: whether a non-trivial initial guess is provided
c               (iguess=1) or not (iguess=0)
c       * matvec (external): matrix-vector product definition.
c
c     Optional variables:
c       * dg: matrix neq*neq diagonal block (for stationary its).
c       * ncolors: number of colors in grid (for GS).
c       * line_relax: whether we want line relaxation
c       * gm_smooth: whether we want GMRES as a smoother
c       * gm_driver: whether we want GMRES as a driver
c       * cvrg_tst: whether we want to perform a convergence test
c       * order_res: restriction order
c       * order_prol: restriction order
c       * smooth: determines smoothing (jb,gs,gm; default jb)
c       * omega: damping for JB smoother
c       * mvcyc: maximum number of MG V-cycles (default 1)
c       * mu_cycle: number of mu_cycles in MG (default 1)
c       * galerkin: whether to do Galerkin coarsening of operator
c                  (true) or rediscretization (default: false)
c     ---------------------------------------------------------------

        use mlsolverSetup

        implicit none

c     Call variables

        integer    :: neq,ntotp,igrid,bcnd(6,neq),out,guess
        real(8)    :: x(ntotp,neq),b(ntotp,neq)

        real(8)   ,optional,intent(IN) :: dg(neq,2*neq*ntotp),omega
        integer   ,optional,intent(IN) :: ncolors,order_res,order_prol
     .                                   ,mvcyc,mu_cycle
        logical   ,optional,intent(IN) :: line_relax,gm_coarse_solve
     .                                   ,cvrg_tst,gm_driver,galerkin

        character(2),optional :: smooth

        external   matvec

c     Local variables

        integer    :: ntot,iorder_res,iorder_prol,vcyc,mu_cyc,nclrs
        real(8)    :: xi(ntotp*neq),bi(ntotp*neq),omg
        real(8), target :: ddg(neq,2*neq*ntotp)
        logical    :: gm_crse_solve,cvrg_test,lrelax,gm_drvr,glkin

        character(2) :: smth

c     Begin program

c     Process optional arguments

        if (PRESENT(dg)) ddg = dg

        if (PRESENT(gm_coarse_solve)) then
          gm_crse_solve = gm_coarse_solve
        else
          gm_crse_solve = .false.
cc          write (*,*) 'DIAG -- cSolver'
cc          gm_crse_solve = .true.
        endif

        if (PRESENT(gm_driver)) then
          gm_drvr = gm_driver
        else
          gm_drvr = .false.
        endif

        if (PRESENT(cvrg_tst)) then
          cvrg_test = cvrg_tst
        else
          cvrg_test = .false.
        endif

        if (PRESENT(line_relax)) then
          lrelax = line_relax
        else
          lrelax = .false.
        endif

        if(PRESENT(order_res)) then
          iorder_res = order_res
        else
          iorder_res = 2
        endif

        if(PRESENT(order_prol)) then
          iorder_prol = order_prol
        else
          iorder_prol = 2
        endif

        if(PRESENT(mvcyc)) then
          vcyc = mvcyc
        else
          vcyc = maxvcyc
        endif

        if(PRESENT(mu_cycle)) then
          mu_cyc = mu_cycle
        else
cc          if (bcSP()) then
cc            mu_cyc = 2
cc          else
            mu_cyc = 1
cc          endif
        endif

        if(PRESENT(smooth)) then
          smth = smooth
        else
          smth = 'jb'
        endif

        if(PRESENT(omega)) then
          omg = omega
        else
          omg = 0.7
        endif

        if(PRESENT(galerkin)) then
          glkin = galerkin
        else
          glkin = .false.
        endif

        if (PRESENT(ncolors)) then
          nclrs = ncolors
        else
          nclrs = 8
        endif

c     Gather variables for coupled solve

        if (.not.cvrg_test) then
          do i=1,ntotp
            do ieq=1,neq
              xi(neq*(i-1)+ieq) = x(i,ieq)
              bi(neq*(i-1)+ieq) = b(i,ieq)
            enddo
          enddo
        else
          bi = 0d0
cc          call random_number(xi)

          xi = 1d0
          write (*,*) 'DIAG -- cSolver: initial guess'
cc          open(unit=11000,file='inputvec',status='unknown')
cc          do i=1,ntotp
cc            do ieq=1,neq
cccc              write (11000,*) xi(neq*(i-1)+ieq)
cc              read (11000,*) xi(neq*(i-1)+ieq)
cc            enddo
cc          enddo
cc          close(11000)

          guess = 1
        endif

c diag ****
cc        if (cvrg_test) then
cc          write (*,*) 'DIAG -- cSolver: JB/GS solve'
cc
ccc         Initialize solver
cc
cc          call solverInit
cc
ccc         Setup solver
cc
cc          call solverOptionsInit
cc
cc          solverOptions%iter    = maxvcyc
cc          solverOptions%tol     = mgtol
cc
cc          if (PRESENT(dg)) solverOptions%diag => ddg
cc
cccc          select case (smth)
cccc          case('gs')
cccc            solverOptions%omega   = 1d0
cc            solverOptions%ncolors = ncolors
cc            call assembleSolverHierarchy('gs')
cccc          case('gm')
cccc            solverOptions%stp_test = 1 
cccc            solverOptions%krylov_subspace = solverOptions%iter 
cccc            call assembleSolverHierarchy('gm')
cccc
cccccc            solverOptions%iter    = 1
cccccc            solverOptions%omega   = omg
cccccc            call assembleSolverHierarchy('jb')
cccc            call assembleSolverHierarchy('id')
cccc          case default !JB
cccccc            solverOptions%omega   = omg
cccc            solverOptions%omega   = 1d0
cccc            call assembleSolverHierarchy('jb')
cccc          end select
cc
ccc         Invoke solver
cc
cc          ntot=neq*ntotp
cc          call getSolver(neq,ntot,bi,xi,matvec,igrid,bcnd,guess,out,1)
cc
ccc         Kill solver
cc
cc          call solverKill
cc
ccc         Scatter solution
cc
cc          do i = 1,ntotp
cc            do ieq=1,neq
cc              x(i,ieq) = xi(neq*(i-1)+ieq)
cc            enddo
cc          enddo
cc
cc          return
cc        endif
c diag ****

c     Solve coupled MG

c       Initialize solver

        call solverInit

c       Check whether GMRES is driving solve

        if (gm_drvr.and.vcyc > 1) then
          call solverOptionsInit

          solverOptions%iter    = vcyc
          solverOptions%tol     = mgtol

          solverOptions%stp_test = 1 
          solverOptions%krylov_subspace = solverOptions%iter 
          call assembleSolverHierarchy('gm')
        endif

c       MG solver options

        call solverOptionsInit

        solverOptions%tol      = mgtol
        if (gm_drvr .and. vcyc > 1) then
          solverOptions%vcyc     = 1
        else
          solverOptions%vcyc     = vcyc
        endif
cc        if (np > 1) then
cc          solverOptions%igridmin = 1
cc        else
cc          write (*,*) 'DIAG -- cSolver'
cc          solverOptions%igridmin = 4
          solverOptions%igridmin = 2
cc          solverOptions%igridmin = 6
cc        endif
        solverOptions%orderres = iorder_res
        solverOptions%orderprol= iorder_prol
        solverOptions%mg_mu    = mu_cyc
        solverOptions%vol_res  = vol_wgt
        solverOptions%galerkin = glkin

        if (PRESENT(dg)) solverOptions%diag => ddg

        solverOptions%ncolors  = nclrs

        !Coarse solve
        if (gm_crse_solve) then
        !With GMRES as coarse solver (defined below)
          if (         (gm_drvr  .and. smth/='gm')
     .       .or.((.not.gm_drvr) .and. smth=='gm') ) then
            solverOptions%mg_coarse_solver_depth = 4
          elseif (gm_drvr .and. smth=='gm') then
            solverOptions%mg_coarse_solver_depth = 5
          else
            solverOptions%mg_coarse_solver_depth = 3
          endif
        else
          solverOptions%mg_coarse_solver_depth = 0 !Use smoother as coarse solver
cc          write (*,*) 'DIAG -- cSolver'
cc          solverOptions%mg_coarse_solver_depth = 5
        endif

        !Vertex relaxation
cc        solverOptions%vertex_based_relax = .true.

        !Plane/line relaxation
        if (lrelax) then
          solverOptions%igridmin        = 3
          solverOptions%orderres        = 0
          solverOptions%orderprol       = 2
          solverOptions%mg_mu           = 1    
          solverOptions%mg_line_relax   = .true.
          solverOptions%mg_line_nsweep  = 1
          solverOptions%mg_line_vcyc    = 100
          solverOptions%mg_line_tol     = 1d-1
          solverOptions%mg_line_omega   = 1d0
          solverOptions%mg_line_x       = .true.
          solverOptions%mg_line_y       = .true.
          solverOptions%mg_line_z       = .true.
          solverOptions%mg_line_solve   = "mg"
          solverOptions%mg_line_coarse_solver_depth = 0
        endif

        call assembleSolverHierarchy('mg')

c       MG smoother

        call solverOptionsInit

        solverOptions%iter    = nsweep
        solverOptions%tol     = mgtol

        if (lrelax) then
          solverOptions%iter  = 2
        endif

        select case (smth)
        case('gs')
          solverOptions%omega   = 1d0
          solverOptions%ncolors = nclrs
          call assembleSolverHierarchy('gs')
        case('gm')
          solverOptions%stp_test        = 1 
          solverOptions%krylov_subspace = solverOptions%iter 
          call assembleSolverHierarchy('gm')

          if (gm_smooth) then
            call assembleSolverHierarchy('id')
          else
            call solverOptionsInit
            solverOptions%iter    = 1
            solverOptions%omega   = omg
            call assembleSolverHierarchy('jb')
          endif
        case default !JB
          solverOptions%omega   = omg
          call assembleSolverHierarchy('jb')
        end select

c       MG coarsest grid solve

        call solverOptionsInit

        solverOptions%tol             = 1d-5
        solverOptions%krylov_subspace = 100
        solverOptions%iter            = 100
        solverOptions%stp_test        = 1 
        solverOptions%omega           = 1d0

        call assembleSolverHierarchy('gm')
        call assembleSolverHierarchy('id') !GMRES preconditioner

c       Coarsest grid solve for inner line/plane MG

        call solverOptionsInit

        solverOptions%omega   = 1d0
        solverOptions%iter    = 100
        solverOptions%tol     = 1d-4
        solverOptions%ncolors = nclrs

cc        call assembleSolverHierarchy('jb')
        call assembleSolverHierarchy('gs')

c       Invoke solver

        ntot=neq*ntotp
        call getSolver(neq,ntot,bi,xi,matvec,igrid,bcnd,guess,out,1)

c       Get output data

cc      call getSolverOptions(1)

c       Kill solver

        call solverKill

c     Scatter solution for output

        do i = 1,ntotp
          do ieq=1,neq
            x(i,ieq) = xi(neq*(i-1)+ieq)
          enddo
        enddo

c     End program

      end subroutine cSolver

c     poisson
c     ###############################################################
      subroutine poisson(nn,b,x,iter,tol,bcnd,iout,guess)

        use mlsolverSetup

c     Call variables

        integer    :: nn,iter,bcnd(6,1),iout,guess
        real(8)    :: b(nn),x(nn),tol

c     Local variables

        integer    :: nn2,ncolors

c     Externals

        external    lap_mtvc

c     Begin program

        ncolors = 2

        nn2 = 2*nn

        if(.not.allocated(d_pois)) then
          allocate(d_pois(1,nn2))
          call find_mf_diag_colored(1,ntotdp,lap_mtvc,1
     .                             ,bcnd,d_pois,ncolors)
        endif

c     Initialize solver

        call solverInit

c     Upper_level solver options

        call solverOptionsInit

        solverOptions%tol      = tol
        solverOptions%iter     = iter
        solverOptions%stp_test = 0         !Use rhs for convergence criteria

        call assembleSolverHierarchy('gm')

c     Next level solver

        call solverOptionsInit

        solverOptions%vcyc     = 1
cc        solverOptions%vcyc     = iter
        solverOptions%igridmin = 2
        solverOptions%orderres = 0
        solverOptions%orderprol= 2
        solverOptions%mg_mu    = 1
        solverOptions%ncolors  = ncolors
        solverOptions%diag     => d_pois
        solverOptions%vol_res  = .true.
        solverOptions%mg_coarse_solver_depth = 0

        call assembleSolverHierarchy('mg')

c     Next level solver

        call solverOptionsInit

        solverOptions%omega = 0.8
        solverOptions%iter  = 3
        solverOptions%tol   = 1d-10

        call assembleSolverHierarchy('jb')

cc        solverOptions%omega = 1d0
cc        solverOptions%iter  = 3
cc        solverOptions%tol   = 1d-10
cc        solverOptions%ncolors = ncolors
cc
cc        call assembleSolverHierarchy('gs')

c       Coarsest grid solve for outer MG

        call solverOptionsInit

        solverOptions%tol             = 1d-5
        solverOptions%krylov_subspace = 1000
        solverOptions%iter            = 1000
        solverOptions%stp_test        = 1 
        solverOptions%omega           = 1d0

c diag ****
cc        call assembleSolverHierarchy('jb')
c diag ****
        call assembleSolverHierarchy('gm')
        call assembleSolverHierarchy('id') !GMRES preconditioner

c     Invoke solver

        call getSolver(1,nn,b,x,lap_mtvc,1,bcnd,guess,iout,1)

c     Get output data

        call getSolverOptions(1)

        iter = solverOptions%iter_out

c     Kill solver

        call solverKill

c     End

      end subroutine poisson

c     vecA
c     ###################################################################
      function vecA(i,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,half_elem)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of vector
c     
c       A = curl(bb x vv + di Dv vv)
c
c     at the grid node (i,j,k). One sided derivatives are employed when
c     half_elem=1 (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,nx,ny,nz,half_elem,igx,igy,igz
        real(8)    :: cnv(3)
        real(8)    :: vv    (0:nx+1,0:ny+1,0:nz+1,3)
     .               ,vv_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .               ,bb    (0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

        integer    :: ig,jg,kg
        real(8)    :: vxx,vyy,vzz,minv,idx,idy,idz

c     Begin program

c     Resistive MHD part

        cnv = curl_bxv(i,j,k,nx,ny,nz,igx,igy,igz,vv,bb,half_elem)

c     Hall MHD part

        if (di > 0d0) then

          !Hall term coefficient
cc          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc          idx  = 1./grid_params%dxh(ig)
cc          if (nxd == 1) idx = 0d0
cc          idy  = 1./grid_params%dyh(jg)
cc          if (nyd == 1) idy = 0d0
cc          idz  = 1./grid_params%dzh(kg)
cc          if (nzd == 1) idz = 0d0
cc
cc          vxx = abs(vv(i,j,k,1))
cc          vyy = abs(vv(i,j,k,2))
cc          vzz = abs(vv(i,j,k,3))
cc
cc          minv = cnp/dt
cc     .           +alpha*(scalarProduct(i,j,k,igx,igy,igz
cc     .                              ,idx,idy,idz,vxx,vyy,vzz)
cc     .           +nu*vectorNorm(i,j,k,igx,igy,igz,idx,idy,idz,.true.))

cc          minv = minv*di/alpha
          minv = cnp/dt*di/alpha

          !Hall term correction
          cnv = cnv
     .        + minv*curl(i,j,k,nx,ny,nz,igx,igy,igz
     .                   ,vv_cov(:,:,:,1)
     .                   ,vv_cov(:,:,:,2)
     .                   ,vv_cov(:,:,:,3),he=half_elem)
     .         
        endif

c     End program

      end function vecA

c     adv_diff_v
c     #####################################################################
      function adv_diff_v(i,j,k,nxx,nyy,nzz,igx,igy,igz,rho0,v0,da)
     .         result(psiv)

c     ---------------------------------------------------------------------
c     Find advection/diffusion component of SI operator
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: i,j,k,igx,igy,igz,nxx,nyy,nzz
      real(8)    :: da(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .             ,v0(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .             ,rho0(0:nxx+1,0:nyy+1,0:nzz+1,1),psiv(3)

c     Local Variables

      integer    :: hex,hey,hez
      real(8)    :: upwind,mul,nabla_vv0(3,3),veclap(3)

c     Begin program

      hex = 1
      hey = 1
      hez = 1
      if (v0(i,j,k,1) > 0d0) hex = -1
      if (v0(i,j,k,2) > 0d0) hey = -1
      if (v0(i,j,k,3) > 0d0) hez = -1
cc      hex = 0
cc      hey = 0
cc      hez = 0

      nabla_v = fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                       ,da(:,:,:,1),da(:,:,:,2),da(:,:,:,3)
     .                       ,hex,hey,hez)

      nabla_vv0 =fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                        ,v0(:,:,:,1)
     .                        ,v0(:,:,:,2)
     .                        ,v0(:,:,:,3)
     .                        ,0,0,0)

      mul = vis(i,j,k,nxx,nyy,nzz,igx,igy,igz)

      veclap = veclaplacian(i,j,k,nxx,nyy,nzz,igx,igy,igz,da
     .                     ,alt_eom,vol=.false.)

      do ieq=1,3
        upwind =( v0(i,j,k,1)*nabla_v(1,ieq)
     .           +v0(i,j,k,2)*nabla_v(2,ieq)
     .           +v0(i,j,k,3)*nabla_v(3,ieq))/jac

        upwind = upwind
     .           +( da(i,j,k,1)*nabla_vv0(1,ieq)
     .             +da(i,j,k,2)*nabla_vv0(2,ieq)
     .             +da(i,j,k,3)*nabla_vv0(3,ieq))/jac

cc        if (.not.nc_eom_v) then
cc          hex = floor(sign(1d0,-mgadvdiffV0(ijkg,ieq)))
cccc          hex = 0
cc          upwind = upwind
cc     .          -dt*mgadvdiffV0(ijkg,ieq)*div_upwd(da,hex)/rho0(i,j,k,1)
cc        endif

        psiv(ieq) = cnp*da(i,j,k,ieq)/dt
     .            + alpha*upwind
     .            - alpha*mul*veclap(ieq)/rho0(i,j,k,1)
      enddo

      end function adv_diff_v

c     diagonalScaling
c     #####################################################################
      subroutine diagonalScaling(neq,ntotp,idiag,y,x,igrid)
c     ---------------------------------------------------------------
c     Performs x=idiag*y, where idiag contains the inverse of the
c     diagonal.
c     ---------------------------------------------------------------

cc      use mg_internal

      implicit none

c     Call variables

      integer    :: neq,ntotp,igrid
      real(8)    :: idiag(neq,*),y(ntotp,neq),x(ntotp,neq)

c     Local variables

      integer    :: ii,iii,iig,isig
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

      isig = istart(igrid)

      do ii = 1,ntotp

        iii = neq*(ii-1)
        iig = iii + isig - 1

        x(ii,:) = matmul(transpose(idiag(:,iig+1:iig+neq)),y(ii,:))

      enddo

      call deallocPointers(fpointers)

c     End program

      end subroutine diagonalScaling

      end module precond_variables

c  module matvec
c ###################################################################
      module matvec

        use setMGBC_interface

        use precond_variables

      end module matvec
