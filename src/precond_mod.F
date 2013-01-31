
c module precond_variables
c ######################################################################
      module precond_variables

        use timeStepping

        use precond_setup

        use transport_params

        use equilibrium

        use auxiliaryVariables

        use operators

        use variables

        use grid

#if !defined(samrai)
        use newtongm
#endif
        use mg_solver

        use nlfunction_setup, ONLY: EOM_divPi_ijk
     .                             ,EOM_divPe_ijk
     .                             ,viscous_heat_src
     .                             ,fake_ve
     .                             ,nc_eom_v
     .                             ,nc_eom_jxb
     .                             ,nc_eom_gp
     .                             ,solve_rho
     .                             ,use_p_eom
     .                             ,advect
     .                             ,v_advect
     .                             ,adiabatic

        use local_BCS_variables, ONLY: order_bc,set_vec_symm_bc

        integer :: ntotd2p,nwt_iter
        logical :: cartsn,covariant,to_cartsn,to_cnv

        !PC storage
        real(8), allocatable, dimension(:,:)    :: rho_diag,tmp_diag
     .                                            ,v_diag,v_diag_hyp
     .                                            ,v_diag2
     .                                            ,b_diag,b_diag_hyp
     .                                            ,a_diag,a_diag_hyp

        real(8), allocatable, dimension(:,:) :: lap_diag,vlap_diag

        real(8), allocatable, dimension(:,:),target :: d_pois

        real(8), allocatable, dimension(:,:)    :: mgj0cnv,mgadvdiffV0

        real(8), allocatable, dimension(:,:,:)  :: mgnablaV0

        real(8), allocatable, dimension(:)      :: mgdivV0

        integer, allocatable, dimension(:,:)    :: bcs,bcs_hyp

        type(mg_array ),target  :: gp0,gb0,gv0,ga0,grho0,gb0_cov,gv0_cov
     .                            ,gdb_cov,gvis,geta,gheta,gve0,gbs0,gj0
     .                            ,gj0_cov,gtmp0,gEh0_cov,gvs0

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

#if !defined(vec_pot)
        allocate (bcs(6,neqd+3),bcs_hyp(6,6))
#else
        allocate (bcs(6,neqd+6),bcs_hyp(6,6))
#endif

        allocate (rho_diag(1,  ntotd2p)
     .           ,tmp_diag(1,  ntotd2p)
#if !defined(vec_pot)
     .           ,b_diag    (3,3*ntotd2p)
     .           ,b_diag_hyp(6,6*ntotd2p)
#else
     .           ,a_diag    (3,3*ntotd2p)
     .           ,a_diag_hyp(6,6*ntotd2p)
#endif
     .           ,v_diag    (3,3*ntotd2p)
     .           ,v_diag_hyp(6,6*ntotd2p)
     .           ,v_diag2   (3,3*ntotd2p),STAT=alloc_stat)

        call allocateMGArray(1,gp0)
        call allocateMGArray(1,grho0)
        call allocateMGArray(1,gtmp0)
        call allocateMGArray(3,gv0)
        call allocateMGArray(3,gv0_cov)
        call allocateMGArray(3,gve0)  !v_e = v_i - di*j/rho
        call allocateMGArray(3,gvs0)  !v_s = v_i - di*j/rho/a_p
        call allocateMGArray(3,gv0)
        call allocateMGArray(3,gb0)
        call allocateMGArray(3,gb0_cov)
        call allocateMGArray(3,gj0)
        call allocateMGArray(3,gj0_cov)
        call allocateMGArray(3,ga0)   !B + di curl(v)
        call allocateMGArray(3,gbs0)  !B - de^2/di curl(v_e)
        call allocateMGArray(3,gEh0_cov)  !Hall electric field
        call allocateMGArray(1,gvis)
        call allocateMGArray(1,geta)
        call allocateMGArray(1,gheta)

        call allocateMGArray(3,gdb_cov)

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
     .             ,a_diag,a_diag_hyp
#endif
     .             ,v_diag,v_diag2,v_diag_hyp)

        call deallocateMGArray(gp0)
        call deallocateMGArray(grho0)
        call deallocateMGArray(gtmp0)
        call deallocateMGArray(gb0)
        call deallocateMGArray(gb0_cov)
        call deallocateMGArray(gj0)
        call deallocateMGArray(gj0_cov)
        call deallocateMGArray(gve0)
        call deallocateMGArray(gvs0)
        call deallocateMGArray(gv0)
        call deallocateMGArray(gv0_cov)
        call deallocateMGArray(ga0)
        call deallocateMGArray(gbs0)
        call deallocateMGArray(gEh0_cov)
        call deallocateMGArray(gdb_cov)
        call deallocateMGArray(gvis)
        call deallocateMGArray(geta)
        call deallocateMGArray(gheta)

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

c     Begin program

        do ii=1,size(xin,1)
          do ieq=1,size(xin,2)
            iii = ieq + size(xin,2)*(ii-1)
            xout(iii) = xin(ii,ieq)
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

        do ii=1,size(xout,1)
          do ieq=1,size(xout,2)
            iii = ieq + size(xout,2)*(ii-1)
            xout(ii,ieq) = xin(iii)
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

                xout(ii,:)= XformToCar(i,j,k,igrid,x(ii,:),covar)
              enddo
            enddo
          enddo
        else
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                cnv(:) = x(ii,:)

                x(ii,:)= XformToCar(i,j,k,igrid,cnv,covar)
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
          call pstop('XtoCurv','Vector has wrong dimensions')
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
     .               ,.not.covar
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
     .               ,.not.covar
     .               ,x(ii,1),x(ii,2),x(ii,3))
              enddo
            enddo
          enddo
        endif

c     End program

      end subroutine XtoCurv

c     fillArray
c     ###################################################################
      subroutine fillArray(igrid,icomp,neq,x,arr,bcs,iorder,arr_cov,arr0
     .                    ,is_vec,is_cnv)

c     -------------------------------------------------------------------
c     Fills array arr (with BCs) using vector x.
c     -------------------------------------------------------------------

        use setMGBC_interface

        implicit none

c     Call variables

        integer :: igrid,icomp,neq,iorder,bcs(:,:)
        real(8),dimension(:,:) :: x
        real(8),dimension(:,:,:,:) :: arr
        real(8),optional,dimension(:,:,:,:) :: arr_cov,arr0
        logical,optional :: is_vec,is_cnv

c     Local variables

        integer :: ieq,nx,ny,nz
        logical :: isvec,iscnv

c     Begin program

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

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
     .                           ,arr(:,:,:,ieq),igrid,.false.)
        enddo

cc        if (icomp == IVX .and. alt__eom()) then
cc          arr(:,:,:,2) = arr(:,:,:,2)/gmetric%grid(igrid)%jac
cc        endif
cc
cc        if (icomp == IRHO .or. icomp == ITMP) then
cc          arr(:,:,:,1) = arr(:,:,:,1)/gmetric%grid(igrid)%jac
cc        endif

        call setMGBC(0,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=(/icomp/),is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder,arr_cov=arr_cov,arr0=arr0)

c     End program

      end subroutine fillArray

c     findPCcoeffs
c     ###################################################################
      subroutine findPCcoeffs(igx,igy,igz,varray,vaux)

c     -------------------------------------------------------------------
c     Finds coefficients for linearized systems in preconditioner
c     -------------------------------------------------------------------

        use xdraw_io

        implicit none

c     Call variables

        integer :: igx,igy,igz
        type(var_array),pointer :: varray
        type(aux_array),pointer :: vaux

c     Local variables

        integer    :: order,nx,ny,nz,igrid,ii,ivar,igr
     .               ,bc_pc(size(bcs,1),size(bcs,2)),ig,jg,kg
     .               ,bc_sc(6,1),bc_vec_odd(6,3),bc_vec_even(6,3)
     .               ,bc_cov(6,3)
     .               ,i,j,k,ip,im,jp,jm,kp,km,nxg,nyg,nzg,ieq

        real(8)    :: veclap(3),cnv(3),cov(3),dxx,dyy,dzz
     .               ,nabla_v(3,3),jac

        real(8), allocatable, dimension(:,:,:,:) :: vector
        real(8), allocatable, dimension(:,:,:,:,:) :: tensor

c     Debug

        real(8),allocatable,dimension(:,:,:) :: debug2
        character(13) :: debugfile

c     Begin program

        nx = grid_params%nxv(igx)
        ny = grid_params%nyv(igy)
        nz = grid_params%nzv(igz)

        igrid = igx

        order = order_bc
cc        order = 0

cc        call setASMflag(.false.)

c     Fill LOCAL auxiliary variables

        call fillLocalAuxVars(igrid,varray,vaux)

c     Setup BCs

        !Set boundary conditions
        do ieq=1,neqd
          bc_pc(:,ieq) = varray%array_var(ieq)%bconds(:)
        enddo

#if defined(vec_pot)
        bc_pc(:,IBX:IBZ) = vaux%vec_list(IBCNV)%bconds

        where (bc_pc(:,IAX:IAZ) == -EQU) bc_pc(:,IAX:IAZ) = -DIR !Use covariant components 
#endif

        bc_pc(:,IJX:IJZ) = vaux%vec_list(IJCNV)%bconds

        call setASM_PC_BC(igx,igy,igz,bc_pc,.true.,IIBC=EXT) !THIS EXTRAPOLATES COEFF. AT LOGICAL BDRYS

        bc_sc(:,1) = bcond
        where (bc_sc == DEF) bc_sc = EXT
        where (bc_sc ==FSYM) bc_sc = NEU

        call setASM_PC_BC(igx,igy,igz,bc_sc,.true.,IIBC=EXT) !THIS EXTRAPOLATES COEFF. AT LOGICAL BDRYS

        bc_vec_odd(:,1) = bcond
        bc_vec_odd(:,2) = bcond
        bc_vec_odd(:,3) = bcond
        where (bc_vec_odd == DEF) bc_vec_odd = EXT
        call set_vec_symm_bc(bc_vec_odd,-1)   !Odd normal component

        call setASM_PC_BC(igx,igy,igz,bc_vec_odd,.true.,IIBC=EXT) !THIS EXTRAPOLATES COEFF. AT LOGICAL BDRYS

        bc_vec_even(:,1) = bcond
        bc_vec_even(:,2) = bcond
        bc_vec_even(:,3) = bcond
        where (bc_vec_even == DEF) bc_vec_even = EXT
        call set_vec_symm_bc(bc_vec_even,1)   !Even normal component

        call setASM_PC_BC(igx,igy,igz,bc_vec_even,.true.,IIBC=EXT) !THIS EXTRAPOLATES COEFF. AT LOGICAL BDRYS

        bc_cov(:,1) = bcond
        bc_cov(:,2) = bcond
        bc_cov(:,3) = bcond
        where (bc_cov == DEF) bc_cov = EXT
        call set_vec_symm_bc(bc_cov,-1)        !Odd normal component

        call setASM_PC_BC(igx,igy,igz,bc_cov,.true.,IIBC=EXT) !THIS EXTRAPOLATES COEFF. AT LOGICAL BDRYS

c     Store density in all grids (w/o BCs)

        grho0%grid(igrid)%array(:,:,:,1) = rho

        call restrictMGArray(IRHO,1,grho0,bc_pc(:,IRHO),igrid,order)

c     Store temperature (w/ BCs)

        gtmp0%grid(igrid)%array(:,:,:,1) = tmp

        call restrictMGArray(ITMP,1,gtmp0,bc_pc(:,ITMP),igrid,order)

c     Store pressure (w/ BCs)

        !pressure (w/ BCs)
        gp0%grid(igrid)%array(:,:,:,1) = a_p*rho*tmp

        call restrictMGArray(ITMP,1,gp0,bc_pc(:,ITMP),igrid,order)

cc        do igr=igrid,grid_params%ngrid
cc          gp0%grid(igrid)%array = a_p*grho0%grid(igrid)%array
cc     .                               *gtmp0%grid(igrid)%array
cc        enddo

c     Store magnetic field cnv components in all grids (w/ BCs)

        gb0%grid(igrid)%array = bcnv

        call restrictMGArray(IBX,3,gb0,bc_pc(:,IBX:IBZ),igrid,order)

c     Store magnetic field cov components in all grids (w/ BCs)

        gb0_cov%grid(igrid)%array = bcov

        call restrictMGArray(IBX,3,gb0_cov,bc_pc(:,IBX:IBZ),igrid,order
     .                      ,iscnv=.false.)

c     Store ion velocity components in all grids (w/ BCs)

        gv0%grid(igrid)%array = vcnv

        call restrictMGArray(IVX,3,gv0,bc_pc(:,IVX:IVZ),igrid,order)

c     Store current cnv components in all grids (w/ BCs)

        gj0%grid(igrid)%array = jcnv

        call restrictMGArray(IJX,3,gj0,bc_pc(:,IJX:IJZ),igrid,order)

c     Store current cov components in all grids (w/ BCs)

        gj0_cov%grid(igrid)%array = jcov

        call restrictMGArray(IJX,3,gj0_cov,bc_pc(:,IJX:IJZ),igrid
     .                      ,order,iscnv=.false.)

c     Store viscosity

        if (nu > 0d0) then
          
          gvis%grid(igrid)%array(:,:,:,1) = 1d0/nuu

          call restrictMGArray(IRHO,1,gvis,bc_sc,igrid,order)

          do igr=igrid,grid_params%ngrid
            gvis%grid(igr)%array = 1d0/gvis%grid(igr)%array
          enddo

        else
          do igr=igrid,grid_params%ngrid
            gvis%grid(igr)%array = 0d0
          enddo
        endif

c     Store resistivity

        geta%grid(igrid)%array(:,:,:,1) = eeta

        call restrictMGArray(IRHO,1,geta,bc_sc,igrid,order)

c     Store Hall MHD auxiliary variables

        if (di > 0d0) then

          !Store hyperresisitivity
          gheta%grid(igrid)%array(:,:,:,1) = h_eta

          call restrictMGArray(IRHO,1,gheta,bc_sc,igrid,order)

          !Store electron velocity components in all grids (w/ BCs)
          gve0%grid(igrid)%array = vecnv

          call restrictMGArray(IVX,3,gve0,bc_vec_odd,igrid,order)
cc          call restrictMGArray(IVX,3,gve0,bc_pc(:,IVX:IVZ),igrid,order)

          !Store hybrid velocity components in all grids (w/ BCs)
          gvs0%grid(igrid)%array = vscnv

          call restrictMGArray(IVX,3,gvs0,bc_vec_odd,igrid,order)

          !Store Hall electric field
          gEh0_cov%grid(igrid)%array = E_h
cc          call restrictMGArray(IJX,3,gEh0_cov,bc_pc(:,IJX:IJZ),igrid
          call restrictMGArray(IJX,3,gEh0_cov,bc_pc(:,IVX:IVZ),igrid
cc          call restrictMGArray(IJX,3,gEh0_cov,bc_cov,igrid
     .                        ,order,iscnv=.false.)

          !Store b_star = b - de^2/di curl(ve)
          if (de > 0d0) then
            do k=1,nz
              do j=1,ny
                do i=1,nx
                  cnv = curl(i,j,k,nx,ny,nz,igrid,igrid,igrid,vecov)
cc                  do ieq=1,3
cc                    cnv(ieq) = curl2(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                              ,vecov(:,:,:,1)
cc     .                              ,vecov(:,:,:,2)
cc     .                              ,vecov(:,:,:,3)
cc     .                              ,ieq)
cc                  enddo

                  gbs0%grid(igrid)%array(i,j,k,:)
     .                   = gb0%grid(igrid)%array(i,j,k,:) -de**2/di*cnv
                enddo
              enddo
            enddo

            call restrictMGArray(IBX,3,gbs0,bc_vec_even,igrid,order)
          else
            call equateMGArray(gbs0,gb0)
          endif

          !Auxiliary quantities for new Hall MHD formulation: A = B + di*curl(v) 
          if (pc_ion_hall) then
            do k=0,nz+1
              do j=0,ny+1
                do i=0,nx+1
                  cnv = curl(i,j,k,nx,ny,nz,igrid,igrid,igrid,vcov)
                  ga0%grid(igrid)%array(i,j,k,:)
     .                   = gb0%grid(igrid)%array(i,j,k,:) + di*cnv
                enddo
              enddo
            enddo

            call restrictMGArray(IJX,3,ga0,bc_vec_even,igrid,order)

          !Auxiliary quantities for std Hall MHD formulation: A = B_star
          else
            call equateMGArray(ga0,gbs0)
          endif

        !Resistive MHD
        else
          call equateMGArray(ga0 ,gb0)
          call equateMGArray(gbs0,gb0)
          call equateMGArray(gvs0,gv0)
        endif

c     Find other auxiliary quantities and store them in all grids

        !Velocity divergence (w/o BCs)
        allocate(debug2(0:nx+1,0:ny+1,0:nz+1))

        debug2 = div(igrid,vcnv)

        call restrictArrayToMGVector(1,nx,ny,nz,debug2,mgdivV0,igrid
     .                            ,order,.false.)

        deallocate(debug2)

        !jxB-grad(p)
        allocate(vector(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,tensor(0:nx+1,0:ny+1,0:nz+1,3,3))

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx

              ip = i+1
              im = i-1
              jp = j+1
              jm = j-1
              kp = k+1
              km = k-1

              ii  = i + nx*(j-1) + nx*ny*(k-1)

              jac    = gmetric%grid(igrid)%jac(i,j,k)

              nabla_v= fnabla_v(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                         ,vcnv(:,:,:,1)
     .                         ,vcnv(:,:,:,2)
     .                         ,vcnv(:,:,:,3),0)
              tensor(i,j,k,:,:) = nabla_v

              if (nc_eom_v) then
                veclap =-div_pi(i,j,k,:)

                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

                dxx = grid_params%dxh(ig)
                dyy = grid_params%dyh(jg)
                dzz = grid_params%dzh(kg)

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
          call restrictArrayToMGVector(1,nx,ny,nz
     .                              ,vector(:,:,:,ivar)
     .                              ,mgadvdiffV0(:,ivar)
     .                              ,igrid,order,.false.)
        enddo

cc        do i=1,3
cc          do j=1,3
cc            call restrictArrayToMGVector(1,nx,ny,nz
cc     .                              ,tensor(:,:,:,i,j)
cc     .                              ,mgnablaV0(:,i,j)
cc     .                              ,igrid,order,.false.)
cc          enddo
cc        enddo

        deallocate(vector,tensor)

c diag ****
cc      write (*,*) 'DIAG -- findCoeffs'
cc
cc      debugfile = 'debug_p'//trim(int2char(my_rank))//'.bin'
cc      open(unit=110,file=debugfile,form='unformatted'
cc     .    ,status='replace')
cc
cc      igx = 2
cc      k = 1
cc
cc      nx = grid_params%nxv(igx)
cc      ny = grid_params%nyv(igx)
cc      nz = grid_params%nzv(igx)
cc      allocate(debug2(0:nx+1,0:ny+1,0:nz+1))
cc
cc      do ieq=1,3
cc        debug2 = gb0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,grid_params%gxmin,grid_params%gxmax
cc     .              ,grid_params%gymin,grid_params%gymax,ieq-1,110)
cc      enddo
cc
cc      do ieq=1,3
cc        debug2 = gb0_cov%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,grid_params%gxmin,grid_params%gxmax
cc     .              ,grid_params%gymin,grid_params%gymax,1,110)
cc      enddo
cc
cc      do ieq=1,3
cc        debug2 = gv0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,grid_params%gxmin,grid_params%gxmax
cc     .              ,grid_params%gymin,grid_params%gymax,1,110)
cc      enddo
cc
cc      do ieq=1,3
cc        debug2 = gj0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,grid_params%gxmin,grid_params%gxmax
cc     .              ,grid_params%gymin,grid_params%gymax,1,110)
cc      enddo
cc
cc      debug2 = grho0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,1)
cc      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .            ,grid_params%gxmin,grid_params%gxmax
cc     .            ,grid_params%gymin,grid_params%gymax,1,110)
cc
cc      debug2 = gp0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,1)
cc      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .            ,grid_params%gxmin,grid_params%gxmax
cc     .            ,grid_params%gymin,grid_params%gymax,1,110)
cc
cc      debug2 = geta%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,1)
cc      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .            ,grid_params%gxmin,grid_params%gxmax
cc     .            ,grid_params%gymin,grid_params%gymax,1,110)
cc
cc      debug2 = gvis%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,1)
cc      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .            ,grid_params%gxmin,grid_params%gxmax
cccc     .            ,grid_params%gymin,grid_params%gymax,0,110)
cc     .            ,grid_params%gymin,grid_params%gymax,1,110)
cc
cc      deallocate(debug2)
cc      close(110)
cc
cc      call PetscEnd(mpierr)
cc      stop
c diag ****

c     Deallocate LOCAL auxiliary variables

        call deallocLocalAuxVars

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
     .              ,icomp=(/ IBX /)
     .              ,is_cnv=.true.          !This informs that input is in cnv array
     .              ,is_vec=.not.cartesian  !This informs whether input is Cartesian
     .              ,result_is_vec=.true.   !This returns curvilinear components
     .              ,iorder=order)
      enddo

c     End program

      end subroutine restrict_B

      end subroutine findPCcoeffs

c     poisson
c     ###############################################################
      subroutine poisson(nn,b,x,iter,tol,bcnd,iout,guess)

        use mlsolverSetup

c     Call variables

        integer :: nn,iter,bcnd(6,1),iout,guess
        real(8) :: b(nn),x(nn),tol

c     Local variables

        integer :: nn2,ncolors
        logical :: vw

c     Externals

        external    lap_mtvc

c     Begin program

        vw = .true. !use volume weighting

        ncolors = 2

        nn2 = 2*nn

        if(.not.allocated(d_pois)) then
          allocate(d_pois(1,nn2))
cc          call find_mf_diag_colored(1,ntotdp,lap_mtvc,1
cc     .                             ,bcnd,d_pois,ncolors)
          call find_mf_diag(1,1,ntotdp,lap_mtvc,1
     .                     ,bcnd,d_pois,ncolors,vw,.false.)
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
        solverOptions%mg_gmin  = 2
        solverOptions%orderres = 0
        solverOptions%orderprol= 2
        solverOptions%mg_mu    = 1
        solverOptions%ncolors  = ncolors
        solverOptions%diag     => d_pois
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

c     D_v
c     #####################################################################
      function D_v(i,j,k,nx,ny,nz,igr,rho0,v0_cnv,dv,dp,conserv,stress)
     .         result(psiv)

c     ---------------------------------------------------------------------
c     Find linearized advection/diffusion component of SI operator
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,nx,ny,nz,igr
      real(8) :: psiv(3)
      real(8),target :: dv    (0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,dp    (0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,v0_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,rho0  (0:nx+1,0:ny+1,0:nz+1,1)
      logical :: conserv
      logical, optional :: stress

c     Local variables

      logical :: strss

c     Begin program

      if (PRESENT(stress)) then
        strss = stress
      else
        strss = .true.
      endif

      if (.not.conserv) then
        !TS
        psiv = cnp*dv(i,j,k,:)/dt

        !Advection
        psiv = psiv
     .       + alpha*(tensor_nc_advc(i,j,k,nx,ny,nz,igr,igr,igr
     .                              ,v0_cnv,dv,1)
     .       +        tensor_nc_advc(i,j,k,nx,ny,nz,igr,igr,igr
     .                              ,dv,v0_cnv,v_advect))

        !Diffusion
        if (nu > 0d0 .and. strss) then
          psiv = psiv
     .         +alpha*EOM_divPi_ijk(i,j,k,nx,ny,nz,igr,igr,igr,dv
     .                             ,gvis%grid(igr)%array(:,:,:,1))
     .                     /rho0(i,j,k,1)
        endif
      else
       !TS
        psiv = cnp*dp(i,j,k,:)/dt

        !Advection
        psiv = psiv
     .       + alpha*(tensor_advc(i,j,k,nx,ny,nz,igr,igr,igr,v0_cnv,dp
     .                           ,1)
cc     .                           ,advect)
     .               +tensor_advc(i,j,k,nx,ny,nz,igr,igr,igr,dp,v0_cnv
     .                           ,v_advect))

        !Diffusion
        if (nu > 0d0 .and. strss) then
          psiv = psiv
     .         +alpha*EOM_divPi_ijk(i,j,k,nx,ny,nz,igr,igr,igr,dv
     .                             ,gvis%grid(igr)%array(:,:,:,1))
        endif
      endif

      end function D_v

c     U_Av_ijk
c     ###################################################################
      function U_Av_ijk(i,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,b0,he
     .                 ,stress) result(uav)

c     -------------------------------------------------------------------
c     Finds covariant components (a1,a2,a3) of vector
c     
c       U_Av(dv) = alpha*(-v x b0 + di (dv/dt + div_pi(v)/rho))
c
c     at the grid node (i,j,k). One sided derivatives are employed when
c     half_elem=1 (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,nx,ny,nz,igx,igy,igz
        real(8)    :: uav(3)
        real(8)    :: vv    (0:nx+1,0:ny+1,0:nz+1,3)
     .               ,vv_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .               ,b0    (0:nx+1,0:ny+1,0:nz+1,3)

        integer,optional :: he
        logical,optional :: stress

c     Local variables

        integer    :: ig,jg,kg,half_elem,ip,jp,kp,hx,hy,hz
        real(8)    :: vxx,vyy,vzz,minv,idx,idy,idz,cnv(3),cov(3)
     .               ,veh(3),curlv(3),half=0.5,ijacip
        logical    :: str

c     Begin program

        if (PRESENT(he)) then
          half_elem = he
        else
          half_elem = 0
        endif

        hx = 0 ; hy = 0 ; hz = 0

        select case(half_elem)
        case(1)
          hx = 1
        case(2)
          hy = 1
        case(3)
          hz = 1
        end select

        if (PRESENT(stress)) then
          str = stress
        else
          str = .false.
        endif

        ip = i+1
        jp = j+1
        kp = k+1

c     Resistive MHD part

        select case(half_elem)
        case(1)

          if (isSP2(ip,igx,dim=1)) then
            cov = 0d0
          else
            cov(1) = (vv(i,j,k,2)+vv(ip,j,k,2))
     .              *(b0(i,j,k,3)+b0(ip,j,k,3))
     .              -(vv(i,j,k,3)+vv(ip,j,k,3))
     .              *(b0(i,j,k,2)+b0(ip,j,k,2))

            cov(2) = (vv(i,j,k,3)+vv(ip,j,k,3))
     .              *(b0(i,j,k,1)+b0(ip,j,k,1))
     .              -(vv(i,j,k,1)+vv(ip,j,k,1))
     .              *(b0(i,j,k,3)+b0(ip,j,k,3))

            cov(3) = (vv(i,j,k,1)+vv(ip,j,k,1))
     .              *(b0(i,j,k,2)+b0(ip,j,k,2))
     .              -(vv(i,j,k,2)+vv(ip,j,k,2))
     .              *(b0(i,j,k,1)+b0(ip,j,k,1))

            cov = 0.5*cov/(gmetric%grid(igx)%jac(i ,j,k)
     .                    +gmetric%grid(igx)%jac(ip,j,k))
          endif

        case(2)

          cov(1) = (vv(i,j,k,2)+vv(i,jp,k,2))
     .            *(b0(i,j,k,3)+b0(i,jp,k,3))
     .            -(vv(i,j,k,3)+vv(i,jp,k,3))
     .            *(b0(i,j,k,2)+b0(i,jp,k,2))
                                       
          cov(2) = (vv(i,j,k,3)+vv(i,jp,k,3))
     .            *(b0(i,j,k,1)+b0(i,jp,k,1))
     .            -(vv(i,j,k,1)+vv(i,jp,k,1))
     .            *(b0(i,j,k,3)+b0(i,jp,k,3))
                                       
          cov(3) = (vv(i,j,k,1)+vv(i,jp,k,1))
     .            *(b0(i,j,k,2)+b0(i,jp,k,2))
     .            -(vv(i,j,k,2)+vv(i,jp,k,2))
     .            *(b0(i,j,k,1)+b0(i,jp,k,1))

          cov = 0.5*cov/(gmetric%grid(igx)%jac(i,j ,k)
     .                  +gmetric%grid(igx)%jac(i,jp,k))

        case(3)

          cov(1) = (vv(i,j,k,2)+vv(i,j,kp,2))
     .            *(b0(i,j,k,3)+b0(i,j,kp,3))
     .            -(vv(i,j,k,3)+vv(i,j,kp,3))
     .            *(b0(i,j,k,2)+b0(i,j,kp,2))
                                         
          cov(2) = (vv(i,j,k,3)+vv(i,j,kp,3))
     .            *(b0(i,j,k,1)+b0(i,j,kp,1))
     .            -(vv(i,j,k,1)+vv(i,j,kp,1))
     .            *(b0(i,j,k,3)+b0(i,j,kp,3))
                                         
          cov(3) = (vv(i,j,k,1)+vv(i,j,kp,1))
     .            *(b0(i,j,k,2)+b0(i,j,kp,2))
     .            -(vv(i,j,k,2)+vv(i,j,kp,2))
     .            *(b0(i,j,k,1)+b0(i,j,kp,1))

          cov = 0.5*cov/(gmetric%grid(igx)%jac(i,j,k )
     .                  +gmetric%grid(igx)%jac(i,j,kp))

        case default
          cov(1) = vv(i,j,k,2)*b0(i,j,k,3)
     .           - vv(i,j,k,3)*b0(i,j,k,2)

          cov(2) = vv(i,j,k,3)*b0(i,j,k,1)
     .           - vv(i,j,k,1)*b0(i,j,k,3)

          cov(3) = vv(i,j,k,1)*b0(i,j,k,2)
     .           - vv(i,j,k,2)*b0(i,j,k,1)

          cov = cov/gmetric%grid(igx)%jac(i,j,k)
        end select

        uav = -cov

c     Hall MHD part

        if (di > 0d0) then

          minv = cnp/dt/alpha

          if (pc_ion_hall) then

            select case(half_elem)
            case(1)
              cov = minv*0.5*(vv_cov(i,j,k,:)+vv_cov(ip,j,k,:))
            case(2)
              cov = minv*0.5*(vv_cov(i,j,k,:)+vv_cov(i,jp,k,:))
            case(3)
              cov = minv*0.5*(vv_cov(i,j,k,:)+vv_cov(i,j,kp,:))
            case default
              cov = minv*vv_cov(i,j,k,:)
            end select

            uav = uav + di*cov

            if (str .and. half_elem == 0) then
cc              select case(half_elem)
cc              case(1)
cc                cnv = (EOM_divPi(i ,j,k,nx,ny,nz,igx,igy,igz,vv
cc     .                               ,gvis%grid(igx)%array)
cc     .                +EOM_divPi(ip,j,k,nx,ny,nz,igx,igy,igz,vv
cc     .                               ,gvis%grid(igx)%array))
cc     .                /(grho0%grid(igx)%array(i ,j,k,1)
cc     .                 +grho0%grid(igx)%array(ip,j,k,1))
cc              case(2)
cc                cnv = (EOM_divPi(i,j ,k,nx,ny,nz,igx,igy,igz,vv
cc     .                          ,gvis%grid(igx)%array)
cc     .                +EOM_divPi(i,jp,k,nx,ny,nz,igx,igy,igz,vv
cc     .                          ,gvis%grid(igx)%array))
cc     .                /(grho0%grid(igx)%array(i,j ,k,1)
cc     .                 +grho0%grid(igx)%array(i,jp,k,1))
cc              case(3)
cc                cnv = (EOM_divPi(i,j,k ,nx,ny,nz,igx,igy,igz,vv
cc     .                          ,gvis%grid(igx)%array)
cc     .                +EOM_divPi(i,j,kp,nx,ny,nz,igx,igy,igz,vv
cc     .                          ,gvis%grid(igx)%array))
cc     .                /(grho0%grid(igx)%array(i,j,k ,1)
cc     .                 +grho0%grid(igx)%array(i,j,kp,1))
cc              case default
                if (nu > 0d0) then
                  cnv = EOM_divPi_ijk(i,j,k ,nx,ny,nz,igx,igy,igz,vv
     .                         ,gvis%grid(igx)%array)
     .               /grho0%grid(igx)%array(i,j,k,1)
                else
                  cnv = 0d0
                endif
cc              end select

                call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                                ,cov(1),cov(2),cov(3)
     .                                ,cnv(1),cnv(2),cnv(3)
     .                                ,.false.,half_elem=half_elem)

                uav = uav + di*cov

            endif

          else
            if (de > 0d0) then

              curlv = curl(i,j,k,nx,ny,nz,igx,igy,igz,vv_cov
     .                    ,hex=hx,hey=hy,hez=hz)

              select case(half_elem)
              case(1)
                veh =(gve0%grid(igx)%array(i ,j,k,:)
     .               +gve0%grid(igx)%array(ip,j,k,:))
     .               /(gmetric%grid(igx)%jac(i ,j,k)
     .                +gmetric%grid(igx)%jac(ip,j,k))

                cov(1) = veh(2)*curlv(3) - veh(3)*curlv(2)
                cov(2) = veh(3)*curlv(1) - veh(1)*curlv(3)
                cov(3) = veh(1)*curlv(2) - veh(2)*curlv(1)

                cov = -cov + minv*0.5*(vv_cov(i,j,k,:)+vv_cov(ip,j,k,:))
cc     .              + (EOM_divPe(i ,j,k,nx,ny,nz,igx,igy,igz,vv)
cc     .                +EOM_divPe(ip,j,k,nx,ny,nz,igx,igy,igz,vv))
cc     .                /(grho0%grid(igx)%array(i ,j,k,1)
cc     .                 +grho0%grid(igx)%array(ip,j,k,1))
              case(2)
                veh =(gve0%grid(igx)%array(i,j ,k,:)
     .               +gve0%grid(igx)%array(i,jp,k,:))
     .               /(gmetric%grid(igx)%jac(i,j ,k)
     .                +gmetric%grid(igx)%jac(i,jp,k))

                cov(1) = veh(2)*curlv(3) - veh(3)*curlv(2)
                cov(2) = veh(3)*curlv(1) - veh(1)*curlv(3)
                cov(3) = veh(1)*curlv(2) - veh(2)*curlv(1)

                cov = -cov + minv*0.5*(vv_cov(i,j,k,:)+vv_cov(i,jp,k,:))
cc     .              + (EOM_divPe(i,j ,k,nx,ny,nz,igx,igy,igz,vv)
cc     .                +EOM_divPe(i,jp,k,nx,ny,nz,igx,igy,igz,vv))
cc     .                /(grho0%grid(igx)%array(i,j ,k,1)
cc     .                 +grho0%grid(igx)%array(i,jp,k,1))
              case(3)
                veh =(gve0%grid(igx)%array(i,j,k ,:)
     .               +gve0%grid(igx)%array(i,j,kp,:))
     .               /(gmetric%grid(igx)%jac(i,j,k )
     .                +gmetric%grid(igx)%jac(i,j,kp))

                cov(1) = veh(2)*curlv(3) - veh(3)*curlv(2)
                cov(2) = veh(3)*curlv(1) - veh(1)*curlv(3)
                cov(3) = veh(1)*curlv(2) - veh(2)*curlv(1)

                cov = -cov + minv*0.5*(vv_cov(i,j,k,:)+vv_cov(i,j,kp,:))
cc     .              + (EOM_divPe(i,j,k ,nx,ny,nz,igx,igy,igz,vv)
cc     .                +EOM_divPe(i,j,kp,nx,ny,nz,igx,igy,igz,vv))
cc     .                /(grho0%grid(igx)%array(i,j,k ,1)
cc     .                 +grho0%grid(igx)%array(i,j,kp,1))
              case default
                veh = gve0%grid(igx)%array(i,j,k ,:)
     .               /gmetric%grid(igx)%jac(i,j,k )

                cov(1) = veh(2)*curlv(3) - veh(3)*curlv(2)
                cov(2) = veh(3)*curlv(1) - veh(1)*curlv(3)
                cov(3) = veh(1)*curlv(2) - veh(2)*curlv(1)

                cov = -cov + minv*vv_cov(i,j,k,:)
cc     .              + EOM_divPe(i,j,k ,nx,ny,nz,igx,igy,igz,vv)
cc     .               /grho0%grid(igx)%array(i,j,k ,1)
              end select

              cov = -de**2/di**2*cov

              uav = uav + di*cov

            endif

            if (str .and. half_elem == 0) then
cc              cnv =-EOM_divPe_ijk(i,j,k,nx,ny,nz,igx,igy,igz,vv)
cc     .                      /grho0%grid(igx)%array(i,j,k,1)
cc              cnv =-heta*EOM_divPi_ijk(i,j,k,nx,ny,nz,igx,igy,igz,vv
cc     .                                ,grho0%grid(igx)%array(:,:,:,1))
cc     .                  /grho0%grid(igx)%array(i,j,k,1)
              cnv =-EOM_divPe_ijk(i,j,k,nx,ny,nz,igx,igy,igz,vv
     .                           ,gheta%grid(igx)%array(:,:,:,1))
     .             /grho0%grid(igx)%array(i,j,k,1)

              call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                                    ,cov(1),cov(2),cov(3)
     .                                    ,cnv(1),cnv(2),cnv(3)
     .                                    ,.false.,half_elem=0)

              uav = uav + di*cov
            endif

          endif

        endif

        uav = alpha*uav

c     End program

      end function U_Av_ijk

ccc     U_Bv
ccc     ###################################################################
cc      function U_Bv(i,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,half_elem)
cc     .         result(cnv)
cc
ccc     -------------------------------------------------------------------
ccc     Finds contravariant components (a1,a2,a3) of vector
ccc     
ccc       A = curl(U_Av(dv))
ccc
ccc     at the grid node (i,j,k). One sided derivatives are employed when
ccc     half_elem=1 (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
ccc     -------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,nx,ny,nz,half_elem,igx,igy,igz
cc        real(8)    :: cnv(3)
cc        real(8)    :: vv    (0:nx+1,0:ny+1,0:nz+1,3)
cc     .               ,vv_cov(0:nx+1,0:ny+1,0:nz+1,3)
cc     .               ,bb    (0:nx+1,0:ny+1,0:nz+1,3)
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg
cc        real(8)    :: vxx,vyy,vzz,minv,idx,idy,idz
cc
ccc     Begin program
cc
ccc     Resistive MHD part
cc
cc        cnv = alpha*curl_bxv(i,j,k,nx,ny,nz,igx,igy,igz,vv,bb,half_elem)
cc
ccc     Hall MHD part
cc
cc        if (di > 0d0) then
cc
cc          if (pc_ion_hall) then
cc             !Hall term coefficient
cccc             call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cccc
cccc             idx  = 1./grid_params%dxh(ig)
cccc             if (nxd == 1) idx = 0d0
cccc             idy  = 1./grid_params%dyh(jg)
cccc             if (nyd == 1) idy = 0d0
cccc             idz  = 1./grid_params%dzh(kg)
cccc             if (nzd == 1) idz = 0d0
cccc
cccc             vxx = abs(vv(i,j,k,1))
cccc             vyy = abs(vv(i,j,k,2))
cccc             vzz = abs(vv(i,j,k,3))
cccc
cccc             minv = cnp/dt
cccc     .              +alpha*(scalarProduct(i,j,k,igx,igy,igz
cccc     .                                 ,idx,idy,idz,vxx,vyy,vzz)
cccc     .              +nu*vectorNorm(i,j,k,igx,igy,igz,idx,idy,idz,.true.))
cc
cccc             minv = minv*di/alpha
cc            minv = cnp/dt*di
cc          else
cc            minv =-de**2/di/dt
cc          endif
cc
cc          cnv = cnv + minv*curl(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                         ,vv_cov(:,:,:,1)
cc     .                         ,vv_cov(:,:,:,2)
cc     .                         ,vv_cov(:,:,:,3),he=half_elem)
cc        endif
cc
ccc     End program
cc
cc      end function U_Bv

c     U_Bv_ijk
c     ###################################################################
      function U_Bv_ijk(i,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he,strs)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of vector
c     
c       U_Bv(dv) = curl(U_Av(dv))
c
c     at the grid node (i,j,k). One sided derivatives are employed when
c     half_elem=1 (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,nx,ny,nz,he,igx,igy,igz
        real(8)    :: cnv(3)
        real(8)    :: vv    (0:nx+1,0:ny+1,0:nz+1,3)
     .               ,vv_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .               ,bb    (0:nx+1,0:ny+1,0:nz+1,3)
        logical,optional :: strs

c     Local variables

        integer    :: ig,jg,kg,im,ip,jm,jp,km,kp
        real(8)    :: idhx,idhy,idhz
     .               ,aip(3),aim(3),ajp(3),ajm(3),akp(3),akm(3)

        real(8)    :: one=1d0,half=0.5d0

c     Begin program

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        idhx = half/grid_params%dxh(ig)
        idhy = half/grid_params%dyh(jg)
        idhz = half/grid_params%dzh(kg)

        select case(he)
        case(1)
          idhx = one/grid_params%dx(ig)

          aip = U_Av_ijk(ip,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=0
     .                  ,stress=strs)
          aim = U_Av_ijk(i ,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=0
     .                  ,stress=strs)

          ajp = U_Av_ijk(i,jp,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          ajm = U_Av_ijk(i,jm,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          
          akp = U_Av_ijk(i,j,kp,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          akm = U_Av_ijk(i,j,km,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)

        case(2)
          idhy = one/grid_params%dy(jg)

          aip = U_Av_ijk(ip,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          aim = U_Av_ijk(im,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)

          ajp = U_Av_ijk(i,jp,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=0
     .                  ,stress=strs)
          ajm = U_Av_ijk(i,j ,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=0
     .                  ,stress=strs)
          
          akp = U_Av_ijk(i,j,kp,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          akm = U_Av_ijk(i,j,km,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)

        case(3)
          idhz = one/grid_params%dz(kg)

          aip = U_Av_ijk(ip,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          aim = U_Av_ijk(im,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)

          ajp = U_Av_ijk(i,jp,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          ajm = U_Av_ijk(i,jm,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          
          akp = U_Av_ijk(i,j,kp,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=0
     .                  ,stress=strs)
          akm = U_Av_ijk(i,j,k ,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=0
     .                  ,stress=strs)
        case default
          aip = U_Av_ijk(ip,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          aim = U_Av_ijk(im,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)

          ajp = U_Av_ijk(i,jp,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          ajm = U_Av_ijk(i,jm,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          
          akp = U_Av_ijk(i,j,kp,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
          akm = U_Av_ijk(i,j,k ,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he=he
     .                  ,stress=strs)
        end select

        cnv(1) = (ajp(3)-ajm(3))*idhy
     .          -(akp(2)-akm(2))*idhz

        cnv(2) = (akp(1)-akm(1))*idhz
     .          -(aip(3)-aim(3))*idhx

        cnv(3) = (aip(2)-aim(2))*idhx
     .          -(ajp(1)-ajm(1))*idhy

c     End program

      end function U_Bv_ijk

c     U_rho_v
c     ###################################################################
      subroutine U_rho_v(nx,ny,nz,dv_cnv,crhs,igrid)
c     -------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the density solve,
c     div(rh0*dv). Not weighed by volume or jacobian factor!
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,nx,ny,nz
      real(8) :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz)

c     Local variables

      integer :: i,j,k,ii

c     Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            crhs(ii) = alpha
     .                *flx_advec(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                      ,dv_cnv,grho0%grid(igrid)%array(:,:,:,1)
     .                      ,advect,vol=.false.)
          enddo
        enddo
      enddo

c     End program

      end subroutine U_rho_v

c     U_Tv
c     ####################################################################
      subroutine U_Tv(nx,ny,nz,dv_cnv,crhs,igrid)
c     --------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the temperature
c     solve. Not weighed by volume or jacobian factor.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,nx,ny,nz
      real(8) :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz)

c     Local variables

      integer :: i,j,k,ii
      real(8),pointer,dimension(:,:,:,:) :: tmp0,nu0,v0,rho0,p0
      real(8) :: heat_src

c     Begin program

c     Evaluate rhs correction: gamma*div(dv*T0) - (gamma-1)*dv*grad(dT)

      p0 => gp0%grid(igrid)%array

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            crhs(ii) = 
     .         gamma    *flx_advec(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                    ,dv_cnv,p0(:,:,:,1),advect)
cc     .       +(gamma-1.)*flx_nc_advec(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                   ,-dv_cnv,p0(:,:,:,1),1,vol=.false.)
     .       +(gamma-1.)*flx_nc_advec(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                   ,dv_cnv,p0(:,:,:,1),1,reverse=.true.)
          enddo
        enddo
      enddo

      if (.not.adiabatic) then
        v0   => gv0  %grid(igrid)%array
        rho0 => grho0%grid(igrid)%array

        allocate(nu0(0:nx+1,0:ny+1,0:nz+1,1))

        nu0 = gvis %grid(igrid)%array
     .       +gheta%grid(igrid)%array        !Add ion+electron contributions

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              heat_src = viscous_heat_src(i,j,k,nx,ny,nz
     .                         ,igrid,igrid,igrid,nu0,dv_cnv,v0)
     .                  +viscous_heat_src(i,j,k,nx,ny,nz
     .                         ,igrid,igrid,igrid,nu0,v0,dv_cnv)

              crhs(ii) =crhs(ii)-(gamma-1.)*heat_src
            enddo
          enddo
        enddo

        nullify(rho0,v0)
        deallocate(nu0)
      endif

      crhs = alpha*crhs

      nullify(tmp0,p0)

c     End program

      end subroutine U_Tv

c     L_T
c     #####################################################################
      subroutine L_T(nx,ny,nz,x,crhs,igrid,vol)
c     ---------------------------------------------------------------
c     This subroutine corrects the T-field update with B, rho info.  It
c     returns dimensions of T if vol=.false., residual dimensions
c     (volume weighed) otherwise.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,igrid
      real(8) :: x(nx*ny*nz,neqd),crhs(nx*ny*nz)
      logical :: vol

c     Local variables

      integer :: i,j,k,ii,isig,ivar
      real(8) :: dvol,heat_src

      real(8),allocatable,dimension(:,:,:,:) ::  db_cnv
     .                                          ,db_cov
     .                                          ,da_cov
     .                                          ,dj_cnv
     .                                          ,dj_cov
     .                                          ,drho
     .                                          ,dve_cnv
     .                                          ,vzeros

      real(8),pointer,dimension(:,:,:,:) :: tmp0,rho0,j0,eta0,ve0,p0

c     Begin program

cc      !To recover old TM Hall performance results
cc      return

      isig = grid_params%istartp(igrid)

      allocate (db_cnv (0:nx+1,0:ny+1,0:nz+1,3)
     .         ,db_cov (0:nx+1,0:ny+1,0:nz+1,3)
     .         ,dj_cnv (0:nx+1,0:ny+1,0:nz+1,3)
     .         ,dj_cov (0:nx+1,0:ny+1,0:nz+1,3)
     .         ,dve_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .         ,vzeros (0:nx+1,0:ny+1,0:nz+1,3)
     .         ,drho   (0:nx+1,0:ny+1,0:nz+1,1))

      vzeros = 0d0

c     Find drho

      call fillArray(igrid,IRHO,1,x(:,IRHO:IRHO),drho,bcs(:,IRHO:IRHO)
     .              ,order_bc)

c     Find dB

#if defined (vec_pot)
      !Find dA
      allocate (da_cov(0:nx+1,0:ny+1,0:nz+1,3))

      call fillArray(igrid,IAX,3,x(:,IAX:IAZ),da_cov,bcs(:,IAX:IAZ)
     .              ,order_bc,is_vec=.true.,is_cnv=.false.)

      !Find dB
      db_cnv = curl(igrid,da_cov)

      deallocate(da_cov)

      db_cov = XformVector(igrid,db_cnv,.false.)
cc      db_cov = db_cnv
cc      call transformVector(igrid,igrid,igrid,0,nx+1,0,ny+1,0,nz+1
cc     .                    ,db_cov(:,:,:,1)
cc     .                    ,db_cov(:,:,:,2)
cc     .                    ,db_cov(:,:,:,3)
cc     .                    ,'cnv','cov')
#else
      call fillArray(igrid,IBX,3,x(:,IBX:IBZ),db_cnv,bcs(:,IBX:IBZ)
     .              ,order_bc,arr_cov=db_cov,is_vec=.true.
     .              ,is_cnv=.true.)
#endif

c     Find dj* from dB*

      dj_cnv = curl(igrid,db_cov)
      dj_cov = XformVector(igrid,dj_cnv,.false.)

      !Find covariant components of dj with BCs
      call setBC(IJX,3,nx,ny,nz,dj_cnv,dj_cov,vzeros,bcs(:,IJX:IJZ)
     .          ,igrid,igrid,igrid,iorder=order_bc)

c     Find Joule heating correction

      eta0 => geta %grid(igrid)%array
      j0   => gj0  %grid(igrid)%array
      rho0 => grho0%grid(igrid)%array

      if (.not.adiabatic) then
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx

              ii  = i + nx*(j-1) + nx*ny*(k-1)

              if (vol) then
                dvol = gmetric%grid(igrid)%dvol(i,j,k)
              else
                dvol = 1d0
              endif

              heat_src = 2*eta0(i,j,k,1)
     .                    *dot_product(j0(i,j,k,:),dj_cov(i,j,k,:))
     .                    /gmetric%grid(igrid)%jac(i,j,k)

              crhs(ii)=crhs(ii) + alpha*dvol*(gamma-1.)*heat_src
            enddo
          enddo
        enddo
      endif

c     Hall MHD correction

      if (di > 0d0) then

        !Find dv_star
        dve_cnv(:,:,:,1)=-di*(dj_cnv(:,:,:,1)/rho0(:,:,:,1)
     .                       -j0    (:,:,:,1)*drho(:,:,:,1)
     .                                       /rho0(:,:,:,1)**2)/a_p  !From v_star def
        dve_cnv(:,:,:,2)=-di*(dj_cnv(:,:,:,2)/rho0(:,:,:,1)
     .                       -j0    (:,:,:,2)*drho(:,:,:,1)
     .                                       /rho0(:,:,:,1)**2)/a_p  !From v_star def
        dve_cnv(:,:,:,3)=-di*(dj_cnv(:,:,:,3)/rho0(:,:,:,1)
     .                       -j0    (:,:,:,3)*drho(:,:,:,1)
     .                                       /rho0(:,:,:,1)**2)/a_p  !From v_star def

        !Find correction to tmp RHS
        p0 => gp0%grid(igrid)%array

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx

              ii  = i + nx*(j-1) + nx*ny*(k-1)
              if (vol) then
                dvol = gmetric%grid(igrid)%dvol(i,j,k)
              else
                dvol = 1d0
              endif

              crhs(ii)=crhs(ii)-alpha*dvol
     .            *(   gamma *flx_advec(i,j,k,nx,ny,nz
     .                                 ,igrid,igrid,igrid
     .                                 ,dve_cnv,p0(:,:,:,1),advect
     .                                 ,vol=.false.)
     .            +(1.-gamma)*dot_product(dve_cnv(i,j,k,:)
     .                ,grad(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                     ,p0(:,:,:,1)
     .                     ,hex=int(sign(1d0,dve_cnv(i,j,k,1)))
     .                     ,hey=int(sign(1d0,dve_cnv(i,j,k,2)))
     .                     ,hez=int(sign(1d0,dve_cnv(i,j,k,3)))))
     .                /gmetric%grid(igrid)%jac(i,j,k))
             enddo
          enddo
        enddo

        if (.not.adiabatic) then

          ve0  => gve0%grid(igrid)%array

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx

                ii  = i + nx*(j-1) + nx*ny*(k-1)

                if (vol) then
                  dvol = gmetric%grid(igrid)%dvol(i,j,k)
                else
                  dvol = 1d0
                endif

cc                  if (fake_ve) then
cc                    heat_src = heat_src
cc     .                     +heta*viscous_heat_src(i,j,k,nx,ny,nz
cc     .                                 ,igrid,igrid,igrid,rho0,dvef_cnv,ve0)
cc     .                     +heta*viscous_heat_src(i,j,k,nx,ny,nz
cc     .                                 ,igrid,igrid,igrid,rho0,vefcnv,dve_cnv)
cc                  else
cc                    heat_src = heat_src
cc     .                       +heta*viscous_heat_src(i,j,k,nx,ny,nz
cc     .                              ,igrid,igrid,igrid,rho0,ve0,dve_cnv)
cc     .                       +heta*viscous_heat_src(i,j,k,nx,ny,nz
cc     .                              ,igrid,igrid,igrid,rho0,dve_cnv,ve0)
                heat_src= viscous_heat_src(i,j,k,nx,ny,nz
     .                            ,igrid,igrid,igrid
     .                            ,gheta%grid(igrid)%array(:,:,:,1)
     .                            ,ve0,dve_cnv)
     .                   +viscous_heat_src(i,j,k,nx,ny,nz
     .                            ,igrid,igrid,igrid
     .                            ,gheta%grid(igrid)%array(:,:,:,1)
     .                            ,dve_cnv,ve0)
cc                  endif

                crhs(ii) = crhs(ii) + alpha*dvol*(gamma-1.)*heat_src
              enddo
            enddo
          enddo

        endif

      endif

c     End program

      nullify(tmp0,rho0,j0,eta0,ve0,p0)

      deallocate(db_cnv
     .          ,db_cov
     .          ,dj_cnv
     .          ,dj_cov
     .          ,drho
     .          ,dve_cnv
     .          ,vzeros)

      end subroutine L_T

c     U_Bv
c     ###################################################################
      subroutine U_Bv(nx,ny,nz,dv_cnv,dv_cov,crhs,igrid)
c     -------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the magnetic field
c     solve: div(dv B0 - B0 dv)
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: igrid,nx,ny,nz
      real(8),target :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
      real(8)    :: dv_cov(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz,3)

c     Local variables

      integer    :: i,j,k,ii

c     Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            crhs(ii,:) = U_Bv_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                       ,dv_cnv,dv_cov
     .                       ,ga0%grid(igrid)%array,0)
          enddo
        enddo
      enddo

c     End program

      end subroutine U_Bv

c     U_Av
c     ###################################################################
      subroutine U_Av(nx,ny,nz,dv_cnv,dv_cov,crhs,igrid)
c     --------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the vector potential
c     solve: U_Av(dv)
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: igrid,nx,ny,nz
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,crhs(nx*ny*nz,3)

c     Local variables

      integer    :: i,j,k,ii
      real(8)    :: minv,cnv(3),cov(3)
     .             ,dp_cnv(0:nx+1,0:ny+1,0:nz+1,3)

c     Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            crhs(ii,:) = U_Av_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                       ,dv_cnv,dv_cov,ga0%grid(igrid)%array
     .                       ,stress=.true.)
          enddo
        enddo
      enddo

c     End program

      end subroutine U_Av

c     diagonalScaling
c     #####################################################################
      subroutine diagonalScaling(neq,ntotp,idiag,y,x,igrid)
c     ---------------------------------------------------------------
c     Performs x=idiag*y, where idiag contains the inverse of the
c     diagonal.
c     ---------------------------------------------------------------

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

        use local_BCS_variables, ONLY: order_bc,nc_eom_v

      end module matvec
