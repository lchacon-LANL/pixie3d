
c module PC_setup
c ######################################################################
      module PC_setup

        use grid

        use ts_setup

        use nlfunction_setup, ONLY: ion_hall,electron_hall,sym_st

        !SI solver parameters
        integer :: mg_ores_si
     .            ,mg_oprol_si
     .            ,sm_iter_si
     .            ,iter_si

        logical :: sm_zebra_relax_si

        real(8) :: sm_omega_si

        character(2) :: sm_type_si

        !Configurable setup (in readInput)
        logical :: pc_B_solve = .true.
     .            ,pc_divclean=.false.

        !Hardwired setup
        logical :: pc_v_solve=.false.          !Whether to perform final v corrector step

        !Internal variables
        logical :: pc_ion_hall
     .            ,pc_e_hall
 
        integer :: iDBG

      contains

c     PC_config
c     ##################################################################
      subroutine PC_config

      implicit none

c     ------------------------------------------------------------------
c     Hardwire setup variables for PC and perform consistency checks.
c     ------------------------------------------------------------------

c     Begin program

      pc_ion_hall  = ion_hall       !Whether to use ion Ohm's law in PC
      pc_e_hall    = electron_hall  !Whether to use electron's Ohm's law in PC

cc      mg_gm_coarse = .true.

      !SCIDAC 08 conference data for Hall MHD
cc      pc_B_solve  = .true.
cc      mg_gm_coarse = .true.
cccc      ores   = 2
cccc      pc_B_solve  = .false.
cccc      divCln = .false.

c     SI solver configuration

      sm_type_si  = sm_type
      sm_omega_si = sm_omega
      sm_iter_si  = sm_iter
      sm_zebra_relax_si = sm_zebra_relax
      mg_ores_si  = mg_ores
      mg_oprol_si = mg_oprol

c     Consistency check

      call TS_solver_consistency_chk

      end subroutine PC_config

      end module PC_setup

#if !defined(plotting) && !defined(_equ)
c module PC_var
c ######################################################################
      module PC_var

        use timeStepping

        use PC_setup

        use transport_params

        use equilibrium

        use auxiliaryVariables

        use operators

        use variables

#if !defined(samrai)
        use nk
#endif
        use mg_solver

        use problem_def, ONLY: use_p_eom

        use nlfunction_setup, ONLY: EOM_divPi_ijk
     .                             ,viscous_heat_src
     .                             ,nc_eom_v
     .                             ,nc_eom_jxb
     .                             ,nc_eom_gp
     .                             ,advect
     .                             ,v_advect
     .                             ,form_ve
     .                             ,vlap_etaj

        use local_BCS_variables, ONLY: order_bc,set_vec_symm_bc

        integer :: ntotd2p,nwt_iter
        logical :: cartsn,covariant,to_cartsn,to_cnv
     $            ,pc_vars_allcd=.false.

        !PC storage
        real(8), allocatable, dimension(:,:) :: rho_diag,tmp_diag
     .                                         ,v_diag,v_diag_hyp
     .                                         ,v_diag2
     .                                         ,b_diag,b_diag_hyp
     .                                         ,a_diag,a_diag_hyp

        real(8), allocatable, dimension(:,:) :: lap_diag!,vlap_diag

        real(8), allocatable, dimension(:,:),target :: d_pois

        integer, allocatable, dimension(:,:) :: bcs,bcs_hyp

        type(mg_array),target  :: gp0,gb0,gv0,ga0,grho0,gb0_cov
     .                           ,gdb_cov,gvis,geta,gheta,gve0,gbs0,gj0
     .                           ,gj0_cov,gtmp0,gEh0_cov,gvs0,gbhat
     .                           ,g_advdiffV0

      contains

c     allocPCVariables
c     ###################################################################
      subroutine allocPCVariables

c     -------------------------------------------------------------------
c     Allocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer    :: alloc_stat

c     Begin program

        ntotd2p = 2*ntotd/neqd

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

        call allocateMGArray(gv%gparams,1,gp0)
        call allocateMGArray(gv%gparams,1,grho0)
        call allocateMGArray(gv%gparams,1,gtmp0)
        call allocateMGArray(gv%gparams,3,gv0)
        call allocateMGArray(gv%gparams,3,gve0)  !v_e = v_i - di*j/rho
        call allocateMGArray(gv%gparams,3,gvs0)  !v_s = v_i - di*j/rho/a_p
        call allocateMGArray(gv%gparams,3,gv0)
        call allocateMGArray(gv%gparams,3,gb0)
        call allocateMGArray(gv%gparams,3,gb0_cov)
        call allocateMGArray(gv%gparams,3,gj0)
        call allocateMGArray(gv%gparams,3,gj0_cov)
        call allocateMGArray(gv%gparams,3,ga0)   !B + di curl(v)
        call allocateMGArray(gv%gparams,3,gbs0)  !B - de^2/di curl(v_e)
        call allocateMGArray(gv%gparams,3,gEh0_cov)  !Hall electric field
        call allocateMGArray(gv%gparams,1,gvis)
        call allocateMGArray(gv%gparams,1,geta)
        call allocateMGArray(gv%gparams,1,gheta)
        call allocateMGArray(gv%gparams,3,gbhat)
        call allocateMGArray(gv%gparams,3,gdb_cov)
        call allocateMGArray(gv%gparams,3,g_advdiffV0)

        pc_vars_allcd = .true.

c     End program

      end subroutine allocPCVariables

c     deallocPCVariables
c     ###################################################################
      subroutine deallocPCVariables

c     -------------------------------------------------------------------
c     Deallocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

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
        call deallocateMGArray(ga0)
        call deallocateMGArray(gbs0)
        call deallocateMGArray(gEh0_cov)
        call deallocateMGArray(gdb_cov)
        call deallocateMGArray(gvis)
        call deallocateMGArray(geta)
        call deallocateMGArray(gheta)
        call deallocateMGArray(gbhat)
        call deallocateMGArray(g_advdiffV0)

        pc_vars_allcd = .false.

c     End program

      end subroutine deallocPCVariables

c     XtoCar
c     ###################################################################
      subroutine XtoCar(g_def,igrid,x,xout,cov)

c     -------------------------------------------------------------------
c     Transforms vector to Cartesian components
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

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

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        if (PRESENT(xout)) then
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                xout(ii,:)= XformToCar(g_def,i,j,k,igrid,x(ii,:),covar)
              enddo
            enddo
          enddo
        else
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                cnv(:) = x(ii,:)

                x(ii,:)= XformToCar(g_def,i,j,k,igrid,cnv,covar)
              enddo
            enddo
          enddo
        endif

c     End program

      end subroutine XtoCar

c     XtoCurv
c     ###################################################################
      subroutine XtoCurv(g_def,igrid,x,xout,cov)

c     -------------------------------------------------------------------
c     Transform vector components from Cartesian to curvilinear
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

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

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        if (PRESENT(xout)) then
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                call transformVectorToCurvilinear
     .               (g_def,i,j,k,igrid,igrid,igrid
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
     .               (g_def,i,j,k,igrid,igrid,igrid
     .               ,car(1),car(2),car(3)
     .               ,.not.covar
     .               ,x(ii,1),x(ii,2),x(ii,3))
              enddo
            enddo
          enddo
        endif

c     End program

      end subroutine XtoCurv

c     findPCcoeffs
c     ###################################################################
      subroutine findPCcoeffs(igrid,varray,vaux)

c     -------------------------------------------------------------------
c     Finds coefficients for linearized systems in preconditioner
c     -------------------------------------------------------------------

        use xdraw_io

        implicit none

c     Call variables

        integer :: igrid
        type(var_array),pointer :: varray
        type(aux_array),pointer :: vaux

c     Local variables

        integer :: order,nx,ny,nz,ii,ivar,igr,igx,igy,igz
     .            ,bc_pc(size(bcs,1),size(bcs,2)),ig,jg,kg
     .            ,bc_sc(6,1),bc_vec_odd(6,3),bc_vec_even(6,3)
     .            ,bc_cov(6,3)
     .            ,i,j,k,ip,im,jp,jm,kp,km,nxg,nyg,nzg,ieq

        real(8) :: cnv(3),cov(3),dxx,dyy,dzz,jac,bnorm

        logical :: extrp_cfss

        integer :: EVEN=1,ODD=-1

c     Begin program

        extrp_cfss=.false.

        igx = igrid
        igy = igrid
        igz = igrid

        nx = gv%gparams%nxv(igx)
        ny = gv%gparams%nyv(igy)
        nz = gv%gparams%nzv(igz)

cc        order = order_bc
        order = 0

c     Fill LOCAL auxiliary variables

        call fillLocalAuxVars(igrid,varray,vaux)

c     Setup BCs

        call set_force_corners(.true.) !Force corners owing to div(A)=0 BC

        !Set boundary conditions
        do ieq=1,neqd
          bc_pc(:,ieq) = varray%array_var(ieq)%bconds(:)
        enddo

#if defined(vec_pot)
        bc_pc(:,IBX:IBZ) = vaux%vec_list(IBCNV)%bconds

        where (bc_pc(:,IAX:IAZ) == -EQU) bc_pc(:,IAX:IAZ) = -DIR !Use covariant components 
#endif

        bc_pc(:,IJX:IJZ) = vaux%vec_list(IJCNV)%bconds

        call setASM_PC_BC(gv%gparams,igx,igy,igz,bc_pc,extrp_cfss
     .                   ,IIBC=EXT)         !THIS EXTRAPOLATES COEFF. AT LOGICAL BDRYS

        bc_sc(:,1) = bcond
        where (bc_sc == DEF) bc_sc = EXT
        where (bc_sc ==FSYM) bc_sc = NEU

        call setASM_PC_BC(gv%gparams,igx,igy,igz,bc_sc,extrp_cfss
     .                   ,IIBC=EXT)         !THIS EXTRAPOLATES COEFF. AT LOGICAL BDRYS

        bc_vec_odd(:,1) = bcond
        bc_vec_odd(:,2) = bcond
        bc_vec_odd(:,3) = bcond
        where (bc_vec_odd == DEF) bc_vec_odd = EXT
        call set_vec_symm_bc(bc_vec_odd,ODD)   !Odd normal component

        call setASM_PC_BC(gv%gparams,igx,igy,igz,bc_vec_odd,extrp_cfss
     .                   ,IIBC=EXT)         !THIS EXTRAPOLATES COEFF. AT LOGICAL BDRYS

        bc_vec_even(:,1) = bcond
        bc_vec_even(:,2) = bcond
        bc_vec_even(:,3) = bcond
        where (bc_vec_even == DEF) bc_vec_even = EXT
        call set_vec_symm_bc(bc_vec_even,EVEN)   !Even normal component

        call setASM_PC_BC(gv%gparams,igx,igy,igz,bc_vec_even,extrp_cfss
     .                   ,IIBC=EXT)         !THIS EXTRAPOLATES COEFF. AT LOGICAL BDRYS

        bc_cov(:,1) = bcond
        bc_cov(:,2) = bcond
        bc_cov(:,3) = bcond
        where (bc_cov == DEF) bc_cov = EXT
        call set_vec_symm_bc(bc_cov,ODD)        !Odd normal component

        call setASM_PC_BC(gv%gparams,igx,igy,igz,bc_cov,extrp_cfss
     .                   ,IIBC=EXT)         !THIS EXTRAPOLATES COEFF. AT LOGICAL BDRYS

c     Store density in all grids (w/o BCs)

        grho0%grid(igrid)%array(:,:,:,1) = rho

        call restrictMGArray(IRHO,1,grho0,bc_pc(:,IRHO),igrid,order
     .                      ,order_bc=order_bc)

c     Store total temperature (w/ BCs)

        gtmp0%grid(igrid)%array(:,:,:,1) = tmp

        call restrictMGArray(IPRS,1,gtmp0,bc_pc(:,IPRS),igrid,order
     .                      ,order_bc=order_bc)

c     Store pressure (w/ BCs)

        !total pressure (w/ BCs)
        gp0%grid(igrid)%array(:,:,:,1) = rho*tmp

        call restrictMGArray(IPRS,1,gp0,bc_pc(:,IPRS),igrid,order
     .                      ,order_bc=order_bc)

c     Store magnetic field cnv components in all grids (w/ BCs)

        gb0    %grid(igrid)%array = bcnv
        gb0_cov%grid(igrid)%array = bcov

        call restrictMGArray(IBX,3,gb0,bc_pc(:,IBX:IBZ),igrid,order
     .                      ,mgvcov=gb0_cov,order_bc=order_bc)

c     Store magnetic field cov components in all grids (w/ BCs)

c$$$        gb0_cov%grid(igrid)%array = bcov
c$$$
c$$$        call restrictMGArray(IBX,3,gb0_cov,bc_pc(:,IBX:IBZ),igrid,order
c$$$     .                      ,iscnv=.false.)

c     Store magnetic field cnv components in all grids (w/ BCs)

        if (chi_par /= chi .and. (.not.lagrangian)) then
          STOP 'need to change this to be coarse-MG friendly'
          gbhat%grid(igrid)%array = bhat

          do igr=igrid+1,gv%gparams%ngrid

            do k=0,gv%gparams%nzv(igr)+1
              do j=0,gv%gparams%nyv(igr)+1
                do i=0,gv%gparams%nxv(igr)+1
                  bnorm = vectorNorm(gv%gparams,i,j,k,igr
     $                     ,gb0%grid(igr)%array(i,j,k,:),.false.)

                  gbhat%grid(igr)%array(i,j,k,:)
     $                 = gb0%grid(igr)%array(i,j,k,:)/sqrt(bnorm)
                enddo
              enddo
            enddo

          enddo

        endif

c     Store ion velocity components in all grids (w/ BCs)

        gv0%grid(igrid)%array = vcnv

        call restrictMGArray(IVX,3,gv0,bc_pc(:,IVX:IVZ),igrid,order
     .                      ,order_bc=order_bc)

c     Store current cnv components in all grids (w/ BCs)

        gj0    %grid(igrid)%array = jcnv
        gj0_cov%grid(igrid)%array = jcov

        call restrictMGArray(IJX,3,gj0,bc_pc(:,IJX:IJZ),igrid,order
     .                      ,mgvcov=gj0_cov,order_bc=order_bc)

c     Store current cov components in all grids (w/ BCs)

cc        gj0_cov%grid(igrid)%array = jcov
cc
cc        call restrictMGArray(IJX,3,gj0_cov,bc_pc(:,IJX:IJZ),igrid
cc     .                      ,order,iscnv=.false.)

c     Store viscosity

cc        gvis%grid(igrid)%array(:,:,:,1) = nuu + h_eta
        gvis%grid(igrid)%array(:,:,:,1) = nuu
        call restrictMGArray(IRHO,1,gvis,bc_sc,igrid,order
     .                      ,order_bc=order_bc)

c     Store resistivity

        geta%grid(igrid)%array(:,:,:,1) = eeta

        call restrictMGArray(IRHO,1,geta,bc_sc,igrid,order
     .                      ,order_bc=order_bc)

c     Store Hall MHD auxiliary variables

        if (di > 0d0) then

          !Store hyperresisitivity
          gheta%grid(igrid)%array(:,:,:,1) = h_eta

          call restrictMGArray(IRHO,1,gheta,bc_sc,igrid,order
     .                      ,order_bc=order_bc)

          !Store electron velocity components in all grids (w/ BCs)
          gve0%grid(igrid)%array = vecnv

          call restrictMGArray(IVX,3,gve0,bc_vec_odd,igrid,order
     .                      ,order_bc=order_bc)
cc          call restrictMGArray(IVX,3,gve0,bc_pc(:,IVX:IVZ),igrid,order)

          !Store hybrid velocity components in all grids (w/ BCs)
          gvs0%grid(igrid)%array = vscnv

          call restrictMGArray(IVX,3,gvs0,bc_vec_odd,igrid,order
     .                      ,order_bc=order_bc)

          !Store Hall electric field
          gEh0_cov%grid(igrid)%array = E_h
cc          call restrictMGArray(IJX,3,gEh0_cov,bc_pc(:,IJX:IJZ),igrid
          call restrictMGArray(IJX,3,gEh0_cov,bc_pc(:,IVX:IVZ),igrid
cc          call restrictMGArray(IJX,3,gEh0_cov,bc_cov,igrid
     .                        ,order,iscnv=.false.,order_bc=order_bc)

          !Store b_star = b - de^2/di curl(ve)
          if (de > 0d0) then
            gbs0%grid(igrid)%array = gb0%grid(igrid)%array
     $                            -de**2/di*curl(gv%gparams,igrid,vecov)

            call setMGBC(gv%gparams,0,3,nx,ny,nz,igrid
     .              ,gbs0%grid(igrid)%array,bc_vec_even
     .              ,icomp=(/ IBX /)
     .              ,is_cnv=.true.           !This informs that input is in cnv array
     .              ,is_curv=.true.          !This informs whether input is Cartesian
     .              ,result_is_curv=.true.   !This returns curvilinear components
     .              ,iorder=order)

            call restrictMGArray(IBX,3,gbs0,bc_vec_even,igrid,order
     .                          ,order_bc=order_bc)
          else
            call equateMGArray(gbs0,gb0)
          endif

          !Auxiliary quantities for new Hall MHD formulation: A = B + di*curl(v) 
          if (pc_ion_hall) then
            ga0%grid(igrid)%array = gb0%grid(igrid)%array
     .                            + di*curl(gv%gparams,igrid,vcov)

            call setMGBC(gv%gparams,0,3,nx,ny,nz,igrid
     .              ,ga0%grid(igrid)%array,bc_vec_even
     .              ,icomp=(/ IBX /)
     .              ,is_cnv=.true.           !This informs that input is in cnv array
     .              ,is_curv=.true.          !This informs whether input is Cartesian
     .              ,result_is_curv=.true.   !This returns curvilinear components
     .              ,iorder=order)

            call restrictMGArray(IBX,3,ga0,bc_vec_even,igrid,order
     .                          ,order_bc=order_bc)

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

        if (nc_eom_v) then

          !jxB-grad(p)
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

                jac    = gv%gparams%gmetric%grid(igrid)%jac(i,j,k)

                call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)

                dxx = gv%gparams%dxh(ig)
                dyy = gv%gparams%dyh(jg)
                dzz = gv%gparams%dzh(kg)

                !Grad(p)
                cov(1) =0.5*(rho(i,j,k)*(tmp(ip,j,k)-tmp(im,j,k))/dxx
     .                      +tmp(i,j,k)*(rho(ip,j,k)-rho(im,j,k))/dxx)

                cov(2) =0.5*(rho(i,j,k)*(tmp(i,jp,k)-tmp(i,jm,k))/dyy
     .                      +tmp(i,j,k)*(rho(i,jp,k)-rho(i,jm,k))/dyy)
        
                cov(3) =0.5*(rho(i,j,k)*(tmp(i,j,kp)-tmp(i,j,km))/dzz
     .                      +tmp(i,j,k)*(rho(i,j,kp)-rho(i,j,km))/dzz)

                cnv = XformToCnv_ijk(gv%gparams,i,j,k,igx,cov)

                !JxB
                cnv(1) = jcov(i,j,k,2)*bcov(i,j,k,3)
     .                 - jcov(i,j,k,3)*bcov(i,j,k,2) - cnv(1)

                cnv(2) = jcov(i,j,k,3)*bcov(i,j,k,1)
     .                 - jcov(i,j,k,1)*bcov(i,j,k,3) - cnv(2)

                cnv(3) = jcov(i,j,k,1)*bcov(i,j,k,2)
     .                 - jcov(i,j,k,2)*bcov(i,j,k,1) - cnv(3)

                g_advdiffV0%grid(igrid)%array(i,j,k,:)
     .                 = (cnv+div_pi(i,j,k,:))/rho(i,j,k)**2

              enddo
            enddo
          enddo

          call restrictMGArray(IVX,3,g_advdiffV0,bc_vec_odd,igrid,order
     .                        ,order_bc=order_bc)

        endif

        call set_force_corners(.false.) !Force corners owing to div(A)=0 BC

c diag ****
cc        call dump_pc_coeffs
c diag ****

c     End program

      contains

c     dump_pc_coeffs
c     ###############################################################
      subroutine dump_pc_coeffs

      implicit none

      real(8),allocatable,dimension(:,:,:) :: debug2
      character(13) :: debugfile

      write (*,*) 'DIAG -- findCoeffs'

      debugfile = 'debug_p'//trim(int2char(my_rank))//'.bin'
      open(unit=110,file=debugfile,form='unformatted'
     .    ,status='replace')

      igx = 1
      k = 1

      nx = gv%gparams%nxv(igx)
      ny = gv%gparams%nyv(igx)
      nz = gv%gparams%nzv(igx)
      allocate(debug2(0:nx+1,0:ny+1,0:nz+1))

      do ieq=1,3
        debug2 = gb0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
     .              ,gv%gparams%lxmin,gv%gparams%lxmax
     .              ,gv%gparams%lymin,gv%gparams%lymax,ieq-1,110)
      enddo

      do ieq=1,3
        debug2 = gve0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
     .              ,gv%gparams%lxmin,gv%gparams%lxmax
     .              ,gv%gparams%lymin,gv%gparams%lymax,1,110)
      enddo

c$$$      do ieq=1,3
c$$$        debug2 = gbs0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
c$$$        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
c$$$     .              ,gv%gparams%lxmin,gv%gparams%lxmax
c$$$     .              ,gv%gparams%lymin,gv%gparams%lymax,1,110)
c$$$      enddo

      do ieq=1,3
        debug2 = ga0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
     .              ,gv%gparams%lxmin,gv%gparams%lxmax
     .              ,gv%gparams%lymin,gv%gparams%lymax,1,110)
      enddo

      do ieq=1,3
        debug2 = gj0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,ieq)
        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
     .              ,gv%gparams%lxmin,gv%gparams%lxmax
     .              ,gv%gparams%lymin,gv%gparams%lymax,1,110)
      enddo

      debug2 = grho0%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,1)
      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
     .            ,gv%gparams%lxmin,gv%gparams%lxmax
     .            ,gv%gparams%lymin,gv%gparams%lymax,1,110)

      debug2 = gheta%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,1)
      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
     .            ,gv%gparams%lxmin,gv%gparams%lxmax
     .            ,gv%gparams%lymin,gv%gparams%lymax,1,110)

      debug2 = geta%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,1)
      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
     .            ,gv%gparams%lxmin,gv%gparams%lxmax
     .            ,gv%gparams%lymin,gv%gparams%lymax,1,110)

      debug2 = gvis%grid(igx)%array(0:nx+1,0:ny+1,0:nz+1,1)
      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
     .            ,gv%gparams%lxmin,gv%gparams%lxmax
     .            ,gv%gparams%lymin,gv%gparams%lymax,1,110)

      deallocate(debug2)
      close(110)

#if defined(petsc)
      call PetscEnd(mpierr)
#endif
      stop

      end subroutine dump_pc_coeffs

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

      nxc = gv%gparams%nxv(igr)
      nyc = gv%gparams%nyv(igr)
      nzc = gv%gparams%nzv(igr)

c     Find cartesian coordinates

      gb0%grid(1)%array = bcnv

      if (cartesian)
     .   call transformVector(gv%gparams
     .                       ,igr,igr,igr,0,nxc+1,0,nyc+1,0,nzc+1
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

      do igr=1,gv%gparams%ngrid
        nxc = gv%gparams%nxv(igr)
        nyc = gv%gparams%nyv(igr)
        nzc = gv%gparams%nzv(igr)

        call setMGBC(gv%gparams,0,neq,nxc,nyc,nzc,igr
     .              ,gb0%grid(igr)%array,bcmod       !Input/output cnv array
     .              ,arr_cov=gb0_cov%grid(igr)%array !Input/output cov array
     .              ,icomp=(/ IBX /)
     .              ,is_cnv=.true.          !This informs that input is in cnv array
     .              ,is_curv=.not.cartesian  !This informs whether input is Cartesian
     .              ,result_is_curv=.true.   !This returns curvilinear components
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

        call MG_set_res_vol_wgt(vw)

        ncolors = 2

        nn2 = 2*nn

        if(.not.allocated(d_pois)) then
          allocate(d_pois(1,nn2))
cc          call find_mf_diag_colored(1,ntotdp,lap_mtvc,1
cc     .                             ,bcnd,d_pois,ncolors)
          call find_mf_diag(gv%gparams,1,1,ntotdp,lap_mtvc,1
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
cc        solverOptions%mg_coarse_grid_res  = 2  !Deprecated
        solverOptions%orderres = 0
        solverOptions%orderprol= 2
        solverOptions%mg_mu    = 1
        solverOptions%ncolors  = ncolors
        solverOptions%diag     => d_pois
        solverOptions%mg_coarse_solver_depth = 0

        call assembleSolverHierarchy('mg')

c     Next level solver

        call solverOptionsInit

        solverOptions%omega = 0.8d0
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
      function D_v(g_def,i,j,k,nx,ny,nz,igr,rho0,v0_cnv,dv,dp,conserv
     .            ,stress) result(psiv)

c     ---------------------------------------------------------------------
c     Find linearized advection/diffusion component of SI operator
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igr
      real(8) :: psiv(3)
      real(8),target :: dv    (0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,dp    (0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,v0_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,rho0  (0:nx+1,0:ny+1,0:nz+1,1)
      real(8),pointer,dimension(:,:,:,:) :: nus
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

      call point_to_MG_array(nx,ny,nz,igr,nus,gvis)

      if (.not.conserv) then
        !TS
        psiv = cnp*dv(i,j,k,:)/dt

        !Advection
        psiv = psiv
     .       + alpha*(tensor_nc_advc(g_def,i,j,k,nx,ny,nz,igr,igr,igr
     .                              ,v0_cnv,dv,1)
     .       +        tensor_nc_advc(g_def,i,j,k,nx,ny,nz,igr,igr,igr
     .                              ,dv,v0_cnv,min(v_advect,2)))

        !Diffusion
        if (strss) then
          psiv = psiv
     .        +alpha/rho0(i,j,k,1)
     .        *EOM_divPi_ijk(g_def,i,j,k,nx,ny,nz,igr,dv,nus(:,:,:,1)
     .                      ,sym_st)
        endif
      else
        !TS
        psiv = cnp*dp(i,j,k,:)/dt

        !Advection
        psiv = psiv
     .       + alpha*(tensor_advc(g_def,i,j,k,nx,ny,nz,igr,igr,igr
     .                           ,v0_cnv,dp,1)
     .               +tensor_advc(g_def,i,j,k,nx,ny,nz,igr,igr,igr
     .                           ,dp,v0_cnv,v_advect))

        !Diffusion
        if (strss) then
          psiv = psiv
     .        +alpha
     .        *EOM_divPi_ijk(g_def,i,j,k,nx,ny,nz,igr,dv,nus(:,:,:,1)
     .                      ,sym_st)
        endif
      endif

      end function D_v

c     U_Av_ijk
c     ###################################################################
      function U_Av_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,b0,he
     .                 ,stress) result(uav)

c     -------------------------------------------------------------------
c     Finds covariant components (a1,a2,a3) of vector
c     
c       U_Av(dv) = alpha*(-dv x b0 + di (ddv/dt + div_pi(dv)/rho0)) (ION HALL)
c
c     or
c
c       U_Av(dv) = alpha*(-dv x b0 + de^2/di ve0 x curl(dv)
c                        - de^2/di dv/dt - di (div_pe(dv)/rho0)) (e HALL)
c
c     at the grid node (i,j,k). One sided derivatives are employed when
c     half_elem=1 (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

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
     .               ,veh(3),curlv(3),half=0.5d0,ijacip
     .               ,jac0,jacp,ijac0,ijacp,ijac
        logical    :: str

        real(8),dimension(:,:,:,:),pointer :: rho0,nu0

c     Begin program

        if (PRESENT(he)) then
          half_elem = he
        else
          half_elem = 0
        endif

        if (PRESENT(stress)) then
          str = stress
        else
          str = .false.
        endif

        ip = i+1
        jp = j+1
        kp = k+1

c     Ideal MHD part

        select case(half_elem)
        case(1)

          jacp = g_def%gmetric%grid(igx)%jac(ip,j,k)
          jac0 = g_def%gmetric%grid(igx)%jac(i ,j,k)

          if (isSP2(g_def,ip,igx,dim=1)) then
            ijac0 = 1d0/jac0
            ijacp = 1d0/jacp

            cov(1) = (vv(i,j,k,2)+vv(ip,j,k,2))
     .              *(b0(i,j,k,3)*ijac0+b0(ip,j,k,3)*ijacp)
     .              -(vv(i,j,k,3)*ijac0+vv(ip,j,k,3)*ijacp)
     .              *(b0(i,j,k,2)+b0(ip,j,k,2))

            cov(2) = 0d0

            cov(3) = (vv(i,j,k,1)*ijac0+vv(ip,j,k,1)*ijacp)
     .              *(b0(i,j,k,2)+b0(ip,j,k,2))
     .              -(vv(i,j,k,2)+vv(ip,j,k,2))
     .              *(b0(i,j,k,1)*ijac0+b0(ip,j,k,1)*ijacp)

            uav = -0.25*cov
          else
            ijac = 2d0/(jacp+jac0)

            uav(1) =-0.25*ijac
     $             *((vv(i,j,k,2)+vv(ip,j,k,2))
     .              *(b0(i,j,k,3)+b0(ip,j,k,3))
     .              -(vv(i,j,k,3)+vv(ip,j,k,3))
     .              *(b0(i,j,k,2)+b0(ip,j,k,2)))

            uav(2) =-0.25*ijac
     $             *((vv(i,j,k,3)+vv(ip,j,k,3))
     .              *(b0(i,j,k,1)+b0(ip,j,k,1))
     .              -(vv(i,j,k,1)+vv(ip,j,k,1))
     .              *(b0(i,j,k,3)+b0(ip,j,k,3)))

            uav(3) =-0.25*ijac
     $             *((vv(i,j,k,1)+vv(ip,j,k,1))
     .              *(b0(i,j,k,2)+b0(ip,j,k,2))
     .              -(vv(i,j,k,2)+vv(ip,j,k,2))
     .              *(b0(i,j,k,1)+b0(ip,j,k,1)))

cc            uav = -0.5*cov/(jacp+jac0)
          endif

        case(2)

          jacp = g_def%gmetric%grid(igx)%jac(i,jp,k)
          jac0 = g_def%gmetric%grid(igx)%jac(i,j ,k)

          ijac = 2d0/(jacp+jac0)

          uav(1) =-0.25*ijac
     $           *((vv(i,j,k,2)+vv(i,jp,k,2))
     .            *(b0(i,j,k,3)+b0(i,jp,k,3))
     .            -(vv(i,j,k,3)+vv(i,jp,k,3))
     .            *(b0(i,j,k,2)+b0(i,jp,k,2)))
                                       
          uav(2) =-0.25*ijac
     $           *((vv(i,j,k,3)+vv(i,jp,k,3))
     .            *(b0(i,j,k,1)+b0(i,jp,k,1))
     .            -(vv(i,j,k,1)+vv(i,jp,k,1))
     .            *(b0(i,j,k,3)+b0(i,jp,k,3)))
                                       
          uav(3) =-0.25*ijac
     $           *((vv(i,j,k,1)+vv(i,jp,k,1))
     .            *(b0(i,j,k,2)+b0(i,jp,k,2))
     .            -(vv(i,j,k,2)+vv(i,jp,k,2))
     .            *(b0(i,j,k,1)+b0(i,jp,k,1)))

cc          uav = -0.5*cov/(jacp+jac0)

        case(3)

          jacp = g_def%gmetric%grid(igx)%jac(i,j,kp)
          jac0 = g_def%gmetric%grid(igx)%jac(i,j,k )

          ijac = 2d0/(jacp+jac0)

          uav(1) =-0.25*ijac
     $           *((vv(i,j,k,2)+vv(i,j,kp,2))
     .            *(b0(i,j,k,3)+b0(i,j,kp,3))
     .            -(vv(i,j,k,3)+vv(i,j,kp,3))
     .            *(b0(i,j,k,2)+b0(i,j,kp,2)))
                                         
          uav(2) =-0.25*ijac
     $           *((vv(i,j,k,3)+vv(i,j,kp,3))
     .            *(b0(i,j,k,1)+b0(i,j,kp,1))
     .            -(vv(i,j,k,1)+vv(i,j,kp,1))
     .            *(b0(i,j,k,3)+b0(i,j,kp,3)))
                                         
          uav(3) =-0.25*ijac
     $           *((vv(i,j,k,1)+vv(i,j,kp,1))
     .            *(b0(i,j,k,2)+b0(i,j,kp,2))
     .            -(vv(i,j,k,2)+vv(i,j,kp,2))
     .            *(b0(i,j,k,1)+b0(i,j,kp,1)))

cc          uav = -0.5*cov/(jacp+jac0)

        case default
          ijac0 = 1d0/g_def%gmetric%grid(igx)%jac(i,j,k)

          uav(1) =-(vv(i,j,k,2)*b0(i,j,k,3)
     .             -vv(i,j,k,3)*b0(i,j,k,2))*ijac0

          uav(2) =-(vv(i,j,k,3)*b0(i,j,k,1)
     .             -vv(i,j,k,1)*b0(i,j,k,3))*ijac0

          uav(3) =-(vv(i,j,k,1)*b0(i,j,k,2)
     .             -vv(i,j,k,2)*b0(i,j,k,1))*ijac0
        end select

cc        uav = -cov !-dv x B0

c     Hall MHD part

        if (di > 0d0) then

          hx = 0 ; hy = 0 ; hz = 0

          select case(half_elem)
          case(1)
             hx = 1
          case(2)
             hy = 1
          case(3)
             hz = 1
          end select

          minv = cnp/dt/alpha

          if (pc_ion_hall) then  !Ion hall

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
c$$$            if (str) then
              call point_to_MG_array(nx,ny,nz,igx,rho0,grho0)
              call point_to_MG_array(nx,ny,nz,igx,nu0 ,gvis)

c$$$              select case(half_elem)
c$$$              case(1)
c$$$                cnv =(EOM_divPi_ijk(g_def,i ,j,k,nx,ny,nz,igx,vv,nu0)
c$$$     .               +EOM_divPi_ijk(g_def,ip,j,k,nx,ny,nz,igx,vv,nu0))
c$$$     .               /(rho0(i ,j,k,1)
c$$$     .                +rho0(ip,j,k,1))
c$$$              case(2)
c$$$                cnv =(EOM_divPi_ijk(g_def,i,j ,k,nx,ny,nz,igx,vv,nu0)
c$$$     .               +EOM_divPi_ijk(g_def,i,jp,k,nx,ny,nz,igx,vv,nu0))
c$$$     .               /(rho0(i,j ,k,1)
c$$$     .                +rho0(i,jp,k,1))
c$$$              case(3)
c$$$                cnv =(EOM_divPi_ijk(g_def,i,j,k ,nx,ny,nz,igx,vv,nu0)
c$$$     .               +EOM_divPi_ijk(g_def,i,j,kp,nx,ny,nz,igx,vv,nu0))
c$$$     .               /(rho0(i,j,k ,1)
c$$$     .                +rho0(i,j,kp,1))
c$$$              case default
                cnv = EOM_divPi_ijk(g_def,i,j,k ,nx,ny,nz,igx,vv,nu0
     .                             ,sym_st)
     .               /rho0(i,j,k,1)
c$$$              end select

              call transformFromCurvToCurv(g_def,i,j,k,igx,igy,igz
     .                                    ,cov(1),cov(2),cov(3)
     .                                    ,cnv(1),cnv(2),cnv(3)
     .                                    ,.false.,half_elem=half_elem)

              uav = uav + di*cov

            endif

          else  !Electron hall

cc            !Electron inertia term: -(de/di)^2*(dv/dt - ve0 x curl(dv))
cc            if (de > 0d0) then
cc
cc              curlv = curl(i,j,k,nx,ny,nz,igx,igy,igz,vv_cov
cc     .                    ,hex=hx,hey=hy,hez=hz)
cc
cc              select case(half_elem)
cc              case(1)
cc                veh =(gve0%grid(igx)%array(i ,j,k,:)
cc     .               +gve0%grid(igx)%array(ip,j,k,:))
cc     .               /(g_def%gmetric%grid(igx)%jac(i ,j,k)
cc     .                +g_def%gmetric%grid(igx)%jac(ip,j,k))
cc
cc                cov(1) = veh(2)*curlv(3) - veh(3)*curlv(2)
cc                cov(2) = veh(3)*curlv(1) - veh(1)*curlv(3)
cc                cov(3) = veh(1)*curlv(2) - veh(2)*curlv(1)
cc
cc                cov = -cov + minv*0.5*(vv_cov(i,j,k,:)+vv_cov(ip,j,k,:))
cc              case(2)
cc                veh =(gve0%grid(igx)%array(i,j ,k,:)
cc     .               +gve0%grid(igx)%array(i,jp,k,:))
cc     .               /(g_def%gmetric%grid(igx)%jac(i,j ,k)
cc     .                +g_def%gmetric%grid(igx)%jac(i,jp,k))
cc
cc                cov(1) = veh(2)*curlv(3) - veh(3)*curlv(2)
cc                cov(2) = veh(3)*curlv(1) - veh(1)*curlv(3)
cc                cov(3) = veh(1)*curlv(2) - veh(2)*curlv(1)
cc
cc                cov = -cov + minv*0.5*(vv_cov(i,j,k,:)+vv_cov(i,jp,k,:))
cc              case(3)
cc                veh =(gve0%grid(igx)%array(i,j,k ,:)
cc     .               +gve0%grid(igx)%array(i,j,kp,:))
cc     .               /(g_def%gmetric%grid(igx)%jac(i,j,k )
cc     .                +g_def%gmetric%grid(igx)%jac(i,j,kp))
cc
cc                cov(1) = veh(2)*curlv(3) - veh(3)*curlv(2)
cc                cov(2) = veh(3)*curlv(1) - veh(1)*curlv(3)
cc                cov(3) = veh(1)*curlv(2) - veh(2)*curlv(1)
cc
cc                cov = -cov + minv*0.5*(vv_cov(i,j,k,:)+vv_cov(i,j,kp,:))
cc              case default
cc                veh = gve0%grid(igx)%array(i,j,k ,:)
cc     .               /g_def%gmetric%grid(igx)%jac(i,j,k )
cc
cc                cov(1) = veh(2)*curlv(3) - veh(3)*curlv(2)
cc                cov(2) = veh(3)*curlv(1) - veh(1)*curlv(3)
cc                cov(3) = veh(1)*curlv(2) - veh(2)*curlv(1)
cc
cc                cov = -cov + minv*vv_cov(i,j,k,:)
cc              end select
cc
cc              uav = uav - de**2/di*cov

            endif

cc            !Stress tensor term: -di/rho0*div_pe(dv)
cc            if (str) then
cc              select case(half_elem)
cc              case(1)
cc                cnv = (EOM_divPi_ijk(i ,j,k,nx,ny,nz,igx,igy,igz,vv
cc     .                               ,heta0)
cc     .                +EOM_divPi_ijk(ip,j,k,nx,ny,nz,igx,igy,igz,vv
cc     .                               ,heta0))
cc     .                /(rho0(i ,j,k,1)
cc     .                 +rho0(ip,j,k,1))
cc              case(2)
cc                cnv = (EOM_divPi_ijk(i,j ,k,nx,ny,nz,igx,igy,igz,vv
cc     .                              ,heta0)
cc     .                +EOM_divPi_ijk(i,jp,k,nx,ny,nz,igx,igy,igz,vv
cc     .                              ,heta0))
cc     .                /(rho0(i,j ,k,1)
cc     .                 +rho0(i,jp,k,1))
cc              case(3)
cc                cnv = (EOM_divPi_ijk(i,j,k ,nx,ny,nz,igx,igy,igz,vv
cc     .                              ,heta0)
cc     .                +EOM_divPi_ijk(i,j,kp,nx,ny,nz,igx,igy,igz,vv
cc     .                              ,heta0))
cc     .                /(rho0(i,j,k ,1)
cc     .                 +rho0(i,j,kp,1))
cc              case default
cc                cnv = EOM_divPi_ijk(i,j,k,nx,ny,nz,igx,igy,igz,vv
cc     .                            ,heta0)
cc     .               /rho0(i,j,k,1)
cc
cc              end select
cc
cc              call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                                ,cov(1),cov(2),cov(3)
cc     .                                ,cnv(1),cnv(2),cnv(3)
cc     .                                ,.false.,half_elem=half_elem)
cc              
cc              uav = uav - di*cov
cc            endif
cc
cc          endif

        endif

        uav = alpha*uav

c     End program

      end function U_Av_ijk

c     U_Bv_ijk
c     ###################################################################
      function U_Bv_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,he
     .                 ,strs) result(cnv)

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

        type(grid_mg_def),pointer :: g_def

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

        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

        idhx = half/g_def%dxh(ig)
        idhy = half/g_def%dyh(jg)
        idhz = half/g_def%dzh(kg)

        select case(he)
        case(1)
          idhx = one/g_def%dx(ig)

          aip = U_Av_ijk(g_def,ip,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=0,stress=strs)
          aim = U_Av_ijk(g_def,i ,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=0,stress=strs)

          ajp = U_Av_ijk(g_def,i,jp,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          ajm = U_Av_ijk(g_def,i,jm,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          
          akp = U_Av_ijk(g_def,i,j,kp,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          akm = U_Av_ijk(g_def,i,j,km,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)

        case(2)
          idhy = one/g_def%dy(jg)

          aip = U_Av_ijk(g_def,ip,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          aim = U_Av_ijk(g_def,im,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)

          ajp = U_Av_ijk(g_def,i,jp,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=0,stress=strs)
          ajm = U_Av_ijk(g_def,i,j ,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=0,stress=strs)
          
          akp = U_Av_ijk(g_def,i,j,kp,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          akm = U_Av_ijk(g_def,i,j,km,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)

        case(3)
          idhz = one/g_def%dz(kg)

          aip = U_Av_ijk(g_def,ip,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          aim = U_Av_ijk(g_def,im,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)

          ajp = U_Av_ijk(g_def,i,jp,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          ajm = U_Av_ijk(g_def,i,jm,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          
          akp = U_Av_ijk(g_def,i,j,kp,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=0,stress=strs)
          akm = U_Av_ijk(g_def,i,j,k ,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=0,stress=strs)
        case default
          aip = U_Av_ijk(g_def,ip,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          aim = U_Av_ijk(g_def,im,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)

          ajp = U_Av_ijk(g_def,i,jp,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          ajm = U_Av_ijk(g_def,i,jm,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          
          akp = U_Av_ijk(g_def,i,j,kp,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
          akm = U_Av_ijk(g_def,i,j,km,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb
     .                  ,he=he,stress=strs)
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
      subroutine U_rho_v(g_def,igrid,dv_cnv,crhs)
c     -------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the density solve,
c     div(rh0*dv). Not weighed by volume or jacobian factor!
c     -------------------------------------------------------------------

      implicit none

c     Call variables
      
      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: dv_cnv(0:,0:,0:,:),crhs(:)

c     Local variables

      integer :: i,j,k,ii,nx,ny,nz

      real(8),dimension(:,:,:,:),pointer :: rho0

c     Begin program

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      call point_to_MG_array(nx,ny,nz,igrid,rho0,grho0)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            crhs(ii) = alpha
     .                *flx_advec(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                      ,dv_cnv,rho0(:,:,:,1),advect,vol=.false.)
          enddo
        enddo
      enddo

c     End program

      end subroutine U_rho_v

c     L_rho_v
c     #########################################################################
      subroutine L_rho_v(g_def,igrid,xxx,yyy,crhs,vwgt)
c     --------------------------------------------------------------------
c     This subroutine finds L_rho_v (form depends on whether we
c     are doing conservative EOM or not).
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid,nx,ny,nz
      real(8) :: xxx(:,:),yyy(:,:),crhs(:,:)
      logical :: vwgt

c     Local variables

      integer :: nn,i,j,k,ii,iig,isig,ivar,ip,im
      real(8) :: dvol,cov(3),cnv(3),idx,idy,idz,eom_rhs(size(xxx,1),3)

      real(8),pointer,dimension(:,:,:,:) :: dp_cnv,drho,v0

c     Begin program

      isig = g_def%istartp(igrid)

      nn = size(crhs,1)

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      allocate(dp_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,drho  (0:nx+1,0:ny+1,0:nz+1,1))

c     Find rhs_v

      call fillArray(g_def,igrid,1,xxx(:,IRHO:IRHO),drho
     .              ,bcs(:,IRHO:IRHO),order_bc,icomp=(/IRHO/))

      if (.not.nc_eom_v) then
        call point_to_MG_array(nx,ny,nz,igrid,v0,gv0)

        dp_cnv(:,:,:,1)=drho(:,:,:,1)*v0(:,:,:,1)
        dp_cnv(:,:,:,2)=drho(:,:,:,1)*v0(:,:,:,2)
        dp_cnv(:,:,:,3)=drho(:,:,:,1)*v0(:,:,:,3)

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              iig = ii + isig - 1
              if (vwgt) then
                dvol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
              else
                dvol = 1d0
              endif

              cnv = tensor_advc(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                         ,v0,dp_cnv,v_advect)

              crhs(ii,:) = dvol*(cnp*dp_cnv(i,j,k,:)/dt - alpha*cnv)

              if (alt__eom()) then
                crhs(ii,2) = crhs(ii,2)
     .                      *g_def%gmetric%grid(igrid)%jac(i,j,k)
              endif
            enddo
          enddo
        enddo

      else
        call point_to_MG_array(nx,ny,nz,igrid,v0,g_advdiffV0)

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              if (vwgt) then
                dvol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
              else
                dvol = 1d0
              endif

              crhs(ii,:)=dvol*drho(i,j,k,1)*v0(i,j,k,:)

              if (alt__eom()) then
                crhs(ii,2) = crhs(ii,2)
     .               *g_def%gmetric%grid(igrid)%jac(i,j,k)
              endif
            enddo
          enddo
        enddo

      endif

c     Free memory

      deallocate(dp_cnv,drho)

c     End program

      end subroutine L_rho_v

c     U_Tv
c     ####################################################################
      subroutine U_Tv(g_def,igr,dv_cnv,crhs)
c     --------------------------------------------------------------------
c     This subroutine finds the U_Tv block using dv. Not weighed by
c     volume or jacobian factor.
c     --------------------------------------------------------------------

      implicit none

c     Call variables
      
      type(grid_mg_def),pointer :: g_def

      integer :: igr
      real(8) :: dv_cnv(0:,0:,0:,:),crhs(:)

c     Local variables

      integer :: i,j,k,ii,nx,ny,nz
      real(8),pointer,dimension(:,:,:,:) :: tmp0,nu0,v0,rho0,p0
      real(8) :: heat_src

c     Begin program

      nx = g_def%nxv(igr)
      ny = g_def%nyv(igr)
      nz = g_def%nzv(igr)

c     Evaluate rhs correction: gamma*div(dv*T0) - (gamma-1)*dv*grad(dT)

      call point_to_MG_array(nx,ny,nz,igr,p0,gp0)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            crhs(ii) = 
     .        gamma    *flx_advec(g_def,i,j,k,nx,ny,nz,igr,igr,igr
     .                            ,dv_cnv,p0(:,:,:,1),advect)
     .      +(gamma-1.)*flx_nc_advec(g_def,i,j,k,nx,ny,nz,igr,igr,igr
     .                            ,dv_cnv,p0(:,:,:,1),1,reverse=.true.)
          enddo
        enddo
      enddo

      if (.not.adiabatic) then
        call point_to_MG_array(nx,ny,nz,igr,v0  ,gv0)
        call point_to_MG_array(nx,ny,nz,igr,rho0,grho0)
        call point_to_MG_array(nx,ny,nz,igr,nu0 ,gvis)

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              heat_src =
     .          viscous_heat_src(g_def,i,j,k,igr,nu0(:,:,:,1)
     .                          ,dv_cnv,v0,sym_st)
     .         +viscous_heat_src(g_def,i,j,k,igr,nu0(:,:,:,1)
     .                          ,v0,dv_cnv,sym_st)

              crhs(ii) =crhs(ii)-(gamma-1d0)*heat_src
            enddo
          enddo
        enddo

        nullify(rho0,v0,nu0)
      endif

      crhs = alpha*crhs

      nullify(tmp0,p0)

c     End program

      end subroutine U_Tv

c     L_T
c     #####################################################################
      subroutine L_T(g_def,igrid,x,crhs)
c     ---------------------------------------------------------------------
c     This subroutine computes L_T temperature block with dB, drho. 
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: x(:,:),crhs(:)

c     Local variables

      integer :: nn,nx,ny,nz,i,j,k,ii,isig,ivar
      real(8) :: heat_src

      real(8),allocatable,dimension(:,:,:,:) ::  db_cnv
     .                                          ,db_cov
     .                                          ,da_cov
     .                                          ,dj_cnv
     .                                          ,dj_cov
     .                                          ,drho
     .                                          ,dve_cnv
     .                                          ,vzeros

      real(8),pointer,dimension(:,:,:,:) :: tmp0,rho0,j0,eta0,ve0,p0
     .                                     ,heta0

c     Begin program

      nn = size(crhs)

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      crhs = 0d0
cc      if (adiabatic .or. (di == 0d0) ) return

cc      !To recover old TM Hall performance results
cc      return

      isig = g_def%istartp(igrid)

      allocate (db_cnv (0:nx+1,0:ny+1,0:nz+1,3)
     .         ,db_cov (0:nx+1,0:ny+1,0:nz+1,3)
     .         ,dj_cnv (0:nx+1,0:ny+1,0:nz+1,3)
     .         ,dj_cov (0:nx+1,0:ny+1,0:nz+1,3)
     .         ,dve_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .         ,vzeros (0:nx+1,0:ny+1,0:nz+1,3)
     .         ,drho   (0:nx+1,0:ny+1,0:nz+1,1))

      vzeros = 0d0

c     Find drho

      call fillArray(g_def,igrid,1,x(:,IRHO:IRHO),drho
     .              ,bcs(:,IRHO:IRHO),order_bc,icomp=(/IRHO/))

c     Find dB

#if defined (vec_pot)
      !Find dA
      allocate (da_cov(0:nx+1,0:ny+1,0:nz+1,3))

      call fillArray(g_def,igrid,3,x(:,IAX:IAZ),da_cov,bcs(:,IAX:IAZ)
     .              ,order_bc,is_vec=.true.,is_cnv=.false.
     .              ,icomp=(/IAX/))

      !Find dB
      db_cnv = curl(g_def,igrid,da_cov)

      deallocate(da_cov)

      db_cov = XformVector(g_def,igrid,db_cnv,.false.)
#else
      call fillArray(g_def,igrid,3,x(:,IBX:IBZ),db_cnv,bcs(:,IBX:IBZ)
     .              ,order_bc,arr_cov=db_cov,is_vec=.true.
     .              ,is_cnv=.true.,icomp=(/IBX/))
#endif

c     Find dj* from dB*

      dj_cnv = curl(g_def,igrid,db_cov)
      dj_cov = XformVector(g_def,igrid,dj_cnv,.false.)

      !Find covariant components of dj with BCs
      call setBC(g_def,IJX,3,nx,ny,nz,dj_cnv,dj_cov,vzeros
     .          ,bcs(:,IJX:IJZ),igrid,igrid,igrid,iorder=order_bc)

c     Find Joule heating correction

      call point_to_MG_array(nx,ny,nz,igrid,eta0,geta)
      call point_to_MG_array(nx,ny,nz,igrid,j0,gj0)
      call point_to_MG_array(nx,ny,nz,igrid,rho0,grho0)

      if (.not.adiabatic) then
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx

              ii  = i + nx*(j-1) + nx*ny*(k-1)

              heat_src = 2*eta0(i,j,k,1)
     .                   *scalarProduct_ijk(g_def,i,j,k,igrid
     .                                     ,j0    (i,j,k,:)
     .                                     ,dj_cov(i,j,k,:))

              crhs(ii) =-(gamma-1d0)*heat_src
            enddo
          enddo
        enddo
      endif

c     Hall MHD correction

      if (di > 0d0) then

        if (.not.(e_cons_vstar.and.adiabatic)) then
          dve_cnv = form_ve(vzeros,dj_cnv,rho0(:,:,:,1))
        endif

cc        !dv_star contribution
cc        if (.not.e_cons_vstar) then
cc
cc          !Find correction to tmp RHS
cc          p0 => gp0%grid(igrid)%array
cc
cc          do k = 1,nz
cc            do j = 1,ny
cc              do i = 1,nx
cc
cc                ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc                crhs(ii)=crhs(ii)
cc     .              +(gamma    *flx_advec(i,j,k,nx,ny,nz
cc     .                                   ,igrid,igrid,igrid
cc     .                                   ,dve_cnv,p0(:,:,:,1),advect)
cc     .              +(gamma-1.)*flx_nc_advec(i,j,k,nx,ny,nz
cc     .                                   ,igrid,igrid,igrid
cc     .                                   ,dve_cnv,p0(:,:,:,1),1
cc     .                                   ,reverse=.true.)  )/a_p  !From v_star def
cc               enddo
cc            enddo
cc          enddo
cc
cc        endif

        !Heat source contribution
        if (.not.adiabatic) then

          call point_to_MG_array(nx,ny,nz,igrid,ve0,gve0)
          call point_to_MG_array(nx,ny,nz,igrid,heta0,gheta)

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx

                ii  = i + nx*(j-1) + nx*ny*(k-1)

                heat_src= viscous_heat_src(g_def,i,j,k,igrid
     .                             ,heta0(:,:,:,1),ve0,dve_cnv,.false.)
     .                   +viscous_heat_src(g_def,i,j,k,igrid
     .                             ,heta0(:,:,:,1),dve_cnv,ve0,.false.)

                crhs(ii) = crhs(ii) - (gamma-1d0)*heat_src
              enddo
            enddo
          enddo

        endif

      endif

      crhs = alpha*crhs

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

c     L_A
c     #########################################################################
      subroutine L_A(nn,xxx,rhs_a,igrid)
c     --------------------------------------------------------------------------
c     This subroutine corrects the rhs for the std Hall vector potential Eq.
c     including density and temperature/pressure corrections.
c     --------------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nn,igrid
      real(8) :: xxx(nn,neqd),rhs_a(nn,3)

c     Local variables

      integer :: i,j,k,ii,ivar,ip,im,nx,ny,nz
      real(8) :: cov(3),cnv(3),idx,idy,idz,eom_rhs(nn,3)

      real(8),allocatable,dimension(:,:,:,:),target :: dpe,drho,vdum
      real(8),pointer    ,dimension(:,:,:,:) :: j0,rho0,tmp0,Eh0

c     Begin program

      nx = gv%gparams%nxv(igrid)
      ny = gv%gparams%nyv(igrid)
      nz = gv%gparams%nzv(igrid)

c     Auxiliary variables

      allocate(dpe (0:nx+1,0:ny+1,0:nz+1,1)
     .        ,drho(0:nx+1,0:ny+1,0:nz+1,1)
     .        ,vdum(0:nx+1,0:ny+1,0:nz+1,3 ))

      call point_to_MG_array(nx,ny,nz,igrid,rho0 ,grho0)
      call point_to_MG_array(nx,ny,nz,igrid,j0   ,gj0)
      call point_to_MG_array(nx,ny,nz,igrid,tmp0 ,gtmp0)
      call point_to_MG_array(nx,ny,nz,igrid,Eh0  ,gEh0_cov)

      call fillArray(gv%gparams,igrid,1,xxx(:,IRHO:IRHO),drho
     .              ,bcs(:,IRHO:IRHO),order_bc,icomp=(/IRHO/))

cc      vdum(:,:,:,1) = di*j0(:,:,:,1)/rho0(:,:,:,1)**2*drho(:,:,:,1)
cc      vdum(:,:,:,2) = di*j0(:,:,:,2)/rho0(:,:,:,1)**2*drho(:,:,:,1)
cc      vdum(:,:,:,3) = di*j0(:,:,:,3)/rho0(:,:,:,1)**2*drho(:,:,:,1)

      call fillArray(gv%gparams,igrid,1,xxx(:,IPRS:IPRS),dpe
     .              ,bcs(:,IPRS:IPRS),order_bc,icomp=(/IPRS/))

      dpe  = dpe/a_p

c     Find rhs_a

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            cov = grad(gv%gparams,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                ,dpe(:,:,:,1))

cc            cnv= EOM_divPi_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid,vdum
cc     .                        ,gheta%grid(igrid)%array(:,:,:,1))
cc            cov = cov + matmul(gmetric%grid(igrid)%gsub(i,j,k,:,:),cnv)

            rhs_a(ii,:)=alpha*di
     .                 *(drho(i,j,k,1)*Eh0(i,j,k,:)+cov)/rho0(i,j,k,1)
          enddo
        enddo
      enddo

c     Free memory

      deallocate(dpe,drho,vdum)

c     End program

      end subroutine L_A

c     U_Bv
c     ###################################################################
      subroutine U_Bv(g_def,igrid,dv_cnv,dv_cov,crhs)
c     -------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the magnetic field
c     solve: div(dv B0 - B0 dv)
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer    :: igrid
      real(8),target :: dv_cnv(0:,0:,0:,:)
      real(8)    :: dv_cov(0:,0:,0:,:),crhs(:,:)

c     Local variables

      integer    :: nx,ny,nz,i,j,k,ii
      real(8),pointer,dimension(:,:,:,:) :: a0_cnv
      
c     Begin program

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      call point_to_MG_array(nx,ny,nz,igrid,a0_cnv,ga0)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            crhs(ii,:) = U_Bv_ijk(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                           ,dv_cnv,dv_cov,a0_cnv,0)
          enddo
        enddo
      enddo

c     End program

      end subroutine U_Bv

c     U_Av
c     ###################################################################
      subroutine U_Av(g_def,igrid,dv_cnv,dv_cov,crhs)
c     --------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the vector potential
c     solve: U_Av(dv)
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer    :: igrid
      real(8)    :: dv_cnv(0:,0:,0:,:)
     .             ,dv_cov(0:,0:,0:,:)
     .             ,crhs(:,:)

c     Local variables

      integer    :: nx,ny,nz,i,j,k,ii
      real(8),pointer,dimension(:,:,:,:) :: a0_cnv,b0_cnv

c     Begin program

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      if (pc_ion_hall) then
        call point_to_MG_array(nx,ny,nz,igrid,a0_cnv,ga0)

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              crhs(ii,:) = U_Av_ijk(g_def,i,j,k,nx,ny,nz
     .                             ,igrid,igrid,igrid
     .                             ,dv_cnv,dv_cov,a0_cnv
     .                             ,stress=.true.)
            enddo
          enddo
        enddo
      else
        call point_to_MG_array(nx,ny,nz,igrid,b0_cnv,gbs0)

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              crhs(ii,:) = U_Av_ijk(g_def,i,j,k,nx,ny,nz
     .                             ,igrid,igrid,igrid
     .                             ,dv_cnv,dv_cov,b0_cnv
     .                             ,stress=.false.)
            enddo
          enddo
        enddo
      endif

c     End program

      end subroutine U_Av

c     diagonalScaling
c     #####################################################################
      subroutine diagonalScaling(g_def,igrid,idiag,y,x)
c     ---------------------------------------------------------------
c     Performs x=idiag*y, where idiag contains the inverse of the
c     diagonal.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer    :: igrid
      real(8)    :: idiag(:,:),y(:,:),x(:,:)

c     Local variables

      integer    :: ii,iii,iig,isig,neq,ntotp
      logical    :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c     Begin program

      neq = size(idiag,1)

      ntotp = size(y,1)

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = mg_ctx_l%istart(igrid)

      do ii = 1,ntotp

        iii = neq*(ii-1)
        iig = iii + isig - 1

        x(ii,:) = matmul(transpose(idiag(:,iig+1:iig+neq)),y(ii,:))

      enddo

      call deallocMGPointers(fpointers,mg_ctx_l)

c     End program

      end subroutine diagonalScaling

      end module PC_var

c  module matvec
c ###################################################################
      module matvec

        use setMGBC_interface

        use PC_var

        use local_BCS_variables, ONLY: order_bc,nc_eom_v

      end module matvec
#endif
