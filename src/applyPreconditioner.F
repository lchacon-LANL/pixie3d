#if !defined(samrai)

c apply_app_PC
c###################################################################
      subroutine apply_app_PC(igrid,nn,neq,yyy,xxx,iout)

c-------------------------------------------------------------------
c     This subroutine solves P x = y for the vector x. 
c     The parameter iout controls output level:
c       * iout <= 0 --> No output
c       * iout >  0 --> Level of output design by user.
c-------------------------------------------------------------------

      use parameters

      use grid

      use precond_variables

      use iosetup

      use imposeBCinterface

      use local_BCS_variables, ONLY: ExB_flow_bc

      implicit none

c Call variables

      integer :: igrid,nn,neq,iout
      real(8) :: xxx(nn,neq),yyy(nn,neq)

c Local variables

      integer :: i,j,k,ieq,ii,iii,guess
     .          ,igx,igy,igz,nx,ny,nz

      real(8) :: dvol,car(3),cnv(3),omg
     .          ,r(nn*neq),x(nn*neq),rrr(nn,neq),dxxx(nn,neq)

      logical :: split_EMHD_after=.false.,electron_hall

c Externals

      external   tmp_mtvc,rho_mtvc,b_mtvc,v_mtvc,v_mtvc2,test_mtvc
     .          ,v_mtvc_af,hall_mtvc,a_mtvc,b_mtvc_hyp,a_mtvc_hyp
     .          ,v_mtvc_hyp

c Begin program

      igx = igrid
      igy = igrid
      igz = igrid

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      if (nn /= nx*ny*nz) then
        write (*,*)'Grid sizes do not agree in apply_app_PC'
        write (*,*)'Aborting...'
        stop
      endif

      xxx = 0d0

      electron_hall = (di > 0d0).and.(.not.pc_ion_hall)

c     Create auxiliary arrays

      guess = 0

      select case(pc_type)

c *******************************************************************
c     Semi-implicit preconditioner version I (w/o flow)
c *******************************************************************

      case('s1')

c     Split PC

        call split_PC

c     Postprocessing of velocity -> momentum

        if (pc_v_solve)
     .       call correctVfld(nn,xxx,yyy(:,IVX:IVZ),igrid,iout)

        if (.not.nc_eom_v) call get_dP

c *******************************************************************
c     Semi-implicit preconditioner version II (iterative, with flow)
c *******************************************************************

      case('s2')

c     Perform iteration

        rrr = yyy

        do iter_si=1,pc_iter

          !SI norm calculation
          if (iout >= 0) call norm_test(rrr,iter_si)

          !SI step
          call si_PC(rrr,dxxx,1)

          !Postprocessing of velocity -> momentum
          if (.not.nc_eom_v) call get_dP

          !Map solution xxx to x (jacobi iterate)
          omg = 1d0
          xxx = xxx + omg*dxxx

          !Form residual vector rr=y-Ax
          call gatherVector(x,xxx)
          call matrixFreeMatVec(nn*neq,x,r)
          call scatter(rrr,r)
          rrr = yyy - rrr

        enddo

c     Final SI norm calculation

        if (iout >= 0) call norm_test(rrr,pc_iter)

c *******************************************************************
c     Semi-implicit preconditioner version III (with flow)
c *******************************************************************

      case('s3')

c     Commuted Schur complement PC

        call star_PC

c     Postprocessing of velocity -> momentum

        if (.not.nc_eom_v) call get_dP

c *******************************************************************
c     EMHD preconditioner
c *******************************************************************

      case('s4')

#if !defined(vec_pot)
        call pstop('applyPC','s4 PC not implemented for B')
#endif

        !MHD step
        rrr = yyy

        call invert_M(rrr,xxx,0)

        call invert_P_si(rrr,xxx)

        call find_dx(rrr,xxx,3)

        !EMHD step
#if defined(vec_pot)
        if (pc_B_solve.and.electron_hall.and.split_EMHD_after) then

          if (iout > 0.and.my_rank ==0) then
            write (*,*) '  --> Solving EMHD step for A...'
          endif

          call rhs_A_hall(nn,xxx,rrr(:,IAX:IAZ),igrid)
          rrr(:,IAX:IAZ) = cnp/dt*xxx(:,IAX:IAZ) - rrr(:,IAZ:IAZ)

          if (vol_wgt) call scaleVec(igrid,3,rrr(:,IAX:IAZ),'mul')

          call EMHD_step(rrr(:,IAX:IAZ),xxx(:,IAX:IAZ))
        endif
#endif

c     Postprocessing of velocity -> momentum

        if (pc_v_solve)
     .       call correctVfld(nn,xxx,yyy(:,IVX:IVZ),igrid,iout)

        if (.not.nc_eom_v) call get_dP

      case default

        call pstop('applyPC','Cannot use PC option '//pc_type)

      end select

c     Apply geometric factors to PC solution to match nonlinear function

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            xxx(ii,IRHO) = xxx(ii,IRHO)*gmetric%grid(igrid)%jac(i,j,k)
            xxx(ii,ITMP) = xxx(ii,ITMP)*gmetric%grid(igrid)%jac(i,j,k)
            if (alt__eom()) then
              xxx(ii,IVY) = xxx(ii,IVY)*gmetric%grid(igrid)%jac(i,j,k)
            endif
          enddo
        enddo
      enddo

c     Deallocate variables

c End program

      contains

c     find_dx
c     ##############################################################
      subroutine find_dx(rr,x,flag)

c     --------------------------------------------------------------
c     Finds update dx=(drho,dB,dT), for:
c       * flag = 0: M^-1 ~ dt update (no inv required)
c       * flag = 1: diagonal splitting (requires inv(D_M))
c       * flag = 2: Approximation of U to [M (D_M)^-1 U] (requires inv(D_M))
c       * flag = 3: standard Schur complement (requires inv(M))
c       * flag = 4: standard Schur complement (requires inv(M), but reuses
c                   y^star -- performs better)
c     --------------------------------------------------------------

        implicit none

c     Call variables

        integer :: flag
        real(8) :: x(nn,neq),rr(nn,neq)

c     Local variables

        real(8) :: udv(nn,neq),dx(nn,neq),dum(nn,neq)

c     Begin program

        !Find U.dv --> udv (dimensions of *physical* dependent variables)
        call find_Udv(x,udv)

        !Solve for dx
        select case(flag)
        case(0)

          x = x - dt/cnp*udv

        case(1)

          call correct_dx_rhs(rr,udv)

          call invert_D_M(udv,x)

        case(2)

          call correct_dx_rhs(rr,udv)

          call invert_D_M(udv,dx)

          x = x + dx

        case(3)

          call correct_dx_rhs(rr,udv)

          call invert_M(udv,x,1)

        case(4)

          dum = 0d0
          call correct_dx_rhs(dum,udv)

          call invert_M(udv,dx,0)

          x = x + dx

        case default

          call pstop('find_dx'
     .              ,'Stationary iteration type not implemented')

        end select

c     End program

      end subroutine find_dx

c     correct_dx_rhs
c     ##############################################################
      subroutine correct_dx_rhs(r,y)

c     --------------------------------------------------------------
c     Corrects the rhs of dx=(drho,dB,dT) equation based on residual
c     r=Gx and y=U.dv (input with dimensions of physical dependent
c     variables, i.e., without volume or jacobian weighing), and stores
c     it in y=r-U.dv
c     --------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: y(nn,neq),r(nn,neq)

c     Local variables

        integer :: ii,i,j,k
        real(8) :: dvol,jac

c     Begin program

        !Find corrected residuals (r-U.dv)
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              if (vol_wgt) then !Residuals are volume-weighed by default
                dvol = gmetric%grid(igx)%dlvol(i,j,k)
                jac  = gmetric%grid(igx)%jac  (i,j,k)
              else
                dvol = 1d0
              endif

              if (solve_rho) y(ii,IRHO)=r(ii,IRHO)-dvol*y(ii,IRHO)*jac
#if !defined(vec_pot)
              y(ii,IBX:IBZ)=r(ii,IBX:IBZ)-dvol*y(ii,IBX:IBZ)
#else
              y(ii,IAX:IAZ)=r(ii,IAX:IAZ)-dvol*y(ii,IAX:IAZ)
#endif
              y(ii,ITMP   )=r(ii,ITMP   )-dvol*y(ii,ITMP   )*jac

            enddo
          enddo
        enddo

c     End program

      end subroutine correct_dx_rhs

c     find_Udv
c     ##############################################################
      subroutine find_Udv(x,rhs)

c     --------------------------------------------------------------
c     Finds the operator U(dv)
c     --------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: x(nn,neq),rhs(nn,neq)

c     Local variables

        real(8) :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .            ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .            ,dv0cnv(0:nx+1,0:ny+1,0:nz+1,3)

c     Begin program

        rhs = 0d0

c     Store velocity solution in array format --> dv_cnv

        dv0cnv = 0d0

        call fillArray(igrid,IVX,3,x(:,IVX:IVZ),dv_cnv,bcs(:,IVX:IVZ)
     .                ,order_bc,arr_cov=dv_cov,is_vec=.true.
     .                ,is_cnv=.true.,arr0=dv0cnv)

c     Correct rx to find Deltax (correction for rho, B, T)

        !Find corrections to rho,B,tmp right hand sides (Udv --> rhs)
        if (solve_rho) then
          call U_rho_v(nx,ny,nz,dv_cnv  ,rhs(:,IRHO)   ,igrid)
        endif

#if !defined(vec_pot)
        call U_Bv(nx,ny,nz,dv_cnv,dv_cov,rhs(:,IBX:IBZ),igrid)
#else
        call U_Av(nx,ny,nz,dv_cnv,dv_cov,rhs(:,IAX:IAZ),igrid)
#endif

        call U_Tv(nx,ny,nz,dv_cnv       ,rhs(:,ITMP)   ,igrid)

c     End program

      end subroutine find_Udv

c     invert_M
c     ##############################################################
      subroutine invert_M(rhs,x,guess)

c     --------------------------------------------------------------
c     Inverts r,B,T diagonal blocks. 
c     --------------------------------------------------------------

        implicit none

c     Call variables

        integer :: guess
        real(8) :: x(nn,neq),rhs(nn,neq),dvol

c     Local variables

        real(8) :: xc(nn,6),rhsc(nn,6),dum(nn)

c     Begin program

        if (iout > 0.and.my_rank ==0) then
          write (*,*) '>>>>>>>>>>>>>>>>>>>>>>'
          write (*,*) 'Solving for M block...'
          write (*,*) '>>>>>>>>>>>>>>>>>>>>>>'
        endif

c       Density

        if (solve_rho) then
          if (iout > 0.and.my_rank ==0) then
            write (*,*) '  --> Solving for rho...'
          endif

          call cSolver(1,nn,rhs(:,IRHO),x(:,IRHO),bcs(:,IRHO)
     .                ,igrid,iout,guess,rho_mtvc,vol_wgt
     .                ,tol           = pc_tol
     .                ,dg            = rho_diag
     .                ,mg_debug      = mg_debug
     .                ,mg_galerkin   = mg_glrkin
     .                ,mg_order_res  = mg_ores
     .                ,mg_order_prol = mg_oprol
     .                ,mg_coarse_grid_size = mg_coarse_size
     .                ,mg_gm_coarse_solve  = mg_gm_coarse
     .                ,mg_vcyc       = mg_vcyc
     .                ,mg_smooth     = sm_type
     .                ,sm_it         = sm_iter
     .                ,sm_omega      = sm_omega
     .                ,sm_ncolors    = sm_ncolors
     .                ,sm_line_relax = sm_line_relax
     .                ,sm_zebra_relax= sm_zebra_relax
     .                )  
        endif

#if !defined(vec_pot)
c       Magnetic field

        if (pc_B_solve) then
          if (iout > 0.and.my_rank ==0) then
            write (*,*) '  --> Solving for B...'
          endif

          if (.not.electron_hall) then
            call cSolver(3,nn,rhs(:,IBX:IBZ),x(:,IBX:IBZ)
     .              ,bcs(:,IBX:IBZ),igrid,iout,guess,b_mtvc,vol_wgt
     .              ,tol           = pc_tol
     .              ,dg            = b_diag
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = mg_ores
     .              ,mg_order_prol = mg_oprol
     .              ,mg_coarse_grid_size = mg_coarse_size
     .              ,mg_gm_coarse_solve  = mg_gm_coarse
     .              ,mg_vcyc       = mg_vcyc
     .              ,mg_smooth     = sm_type
     .              ,sm_it         = sm_iter
     .              ,sm_omega      = sm_omega
     .              ,sm_ncolors    = sm_ncolors
     .              ,sm_line_relax = sm_line_relax
     .              ,sm_zebra_relax= sm_zebra_relax
     .              )
          else
            xc(:,1:3) = x(:,IBX:IBZ)
            xc(:,4:6) = 0d0

            rhsc(:,1:3) = rhs(:,IBX:IBZ)
            rhsc(:,4:6) = 0d0

            call cSolver(6,nn,rhsc,xc,bcs_hyp
     .              ,igrid,iout,guess,b_mtvc_hyp,vol_wgt
     .              ,tol           = pc_tol
     .              ,dg            = b_diag_hyp
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = mg_ores
     .              ,mg_order_prol = mg_oprol
     .              ,mg_coarse_grid_size = mg_coarse_size
     .              ,mg_gm_coarse_solve  = mg_gm_coarse
     .              ,mg_vcyc       = mg_vcyc
     .              ,mg_smooth     = sm_type
     .              ,sm_it         = sm_iter
     .              ,sm_omega      = sm_omega
     .              ,sm_ncolors    = sm_ncolors
     .              ,sm_line_relax = sm_line_relax
     .              ,sm_zebra_relax= sm_zebra_relax
     .              )

            x(:,IBX:IBZ) = xc(:,1:3)
          endif

          if (pc_divclean) call divClean(iout,nn,x(:,IBX:IBZ),igrid)
cc          if (pc_divclean) then
cc            call correctBfld(nn,igrid,rhs(:,IBX:IBZ),x,iout)
cc          endif

        else

          x(:,IBX:IBZ) = dt*rhs(:,IBX:IBZ)/cnp
          call scaleVec(igrid,3,x(:,IBX:IBZ),'div')

        endif
#else
c       Vector potential

        if (pc_B_solve.and.electron_hall
     .      .and.(.not.split_EMHD_after)) then

          if (iout > 0.and.my_rank ==0) then
            write (*,*) '  --> Solving EMHD step for A...'
          endif

          call rhs_A_hall(nn,x,rhsc(:,1:3),igrid)
          if (vol_wgt) call scaleVec(igrid,3,rhsc(:,1:3),'mul')
          rhsc(:,1:3) = rhs(:,IAX:IAZ) - rhsc(:,1:3)

          call EMHD_step(rhsc(:,1:3),x(:,IAX:IAZ))

        elseif (pc_B_solve.and.(.not.electron_hall)) then

          if (iout > 0.and.my_rank ==0) then
            write (*,*) '  --> Solving for A...'
          endif

          call cSolver(3,nn,rhs(:,IAX:IAZ),x(:,IAX:IAZ)
     .           ,bcs(:,IAX:IAZ),igrid,iout,guess,a_mtvc,vol_wgt
     .           ,tol           = pc_tol
     .           ,dg            = a_diag
     .           ,mg_debug      = mg_debug
     .           ,mg_galerkin   = mg_glrkin
     .           ,mg_order_res  = mg_ores
     .           ,mg_order_prol = mg_oprol
     .           ,mg_coarse_grid_size = mg_coarse_size
     .           ,mg_gm_coarse_solve  = mg_gm_coarse
     .           ,mg_vcyc       = mg_vcyc
     .           ,mg_smooth     = sm_type
     .           ,sm_it         = sm_iter
     .           ,sm_omega      = sm_omega
     .           ,sm_ncolors    = sm_ncolors
     .           ,sm_line_relax = sm_line_relax
     .           ,sm_zebra_relax= sm_zebra_relax
     .           )

        else

          x(:,IAX:IAZ) = dt*rhs(:,IAX:IAZ)/cnp
          call scaleVec(igrid,3,x(:,IAX:IAZ),'div')

        endif
#endif

c       Temperature

        if (iout > 0.and.my_rank ==0) then
          write (*,*) '  --> Solving for T...'
        endif

        !Correct for electron velocity/source contribution
        call L_T(nn,x,dum,igrid)
        if (vol_wgt) call scaleVec(igrid,1,dum,'mul')
        dum = rhs(:,ITMP) - dum

        call cSolver(1,nn,dum,x(:,ITMP),bcs(:,ITMP)
     .            ,igrid,iout,guess,tmp_mtvc,vol_wgt
     .            ,tol           = pc_tol
     .            ,dg            = tmp_diag
     .            ,mg_debug      = mg_debug
     .            ,mg_galerkin   = mg_glrkin
     .            ,mg_order_res  = mg_ores
     .            ,mg_order_prol = mg_oprol
     .            ,mg_coarse_grid_size = mg_coarse_size
     .            ,mg_gm_coarse_solve  = mg_gm_coarse
     .            ,mg_vcyc       = mg_vcyc
     .            ,mg_smooth     = sm_type
     .            ,sm_it         = sm_iter
     .            ,sm_omega      = sm_omega
     .            ,sm_ncolors    = sm_ncolors
     .            ,sm_line_relax = sm_line_relax
     .            ,sm_zebra_relax= sm_zebra_relax
     .            )

c     End program

      end subroutine invert_M

c     invert_D_M
c     ##############################################################
      subroutine invert_D_M(rhs,x)

c     --------------------------------------------------------------
c     Inverts block diagonal part of M
c     --------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: x(nn,neq),rhs(nn,neq)

c     Local variables

c     Begin program
        if (iout > 0.and.my_rank ==0) then
          write (*,*) '>>>>>>>>>>>>>>>>>>>'
          write (*,*) 'Diagonal scaling...'
          write (*,*) '>>>>>>>>>>>>>>>>>>>'
        endif

        call diagonalScaling(1,nn,rho_diag,rhs(:,IRHO)
     .                      ,x(:,IRHO)   ,igrid)
        call diagonalScaling(1,nn,tmp_diag,rhs(:,ITMP)
     .                      ,x(:,ITMP)   ,igrid)
#if !defined(vec_pot)
        call diagonalScaling(3,nn,b_diag  ,rhs(:,IBX:IBZ)
     .                      ,x(:,IBX:IBZ),igrid)
#else
        call diagonalScaling(3,nn,a_diag  ,rhs(:,IAX:IAZ)
     .                      ,x(:,IAX:IAZ),igrid)
#endif

c     End program

      end subroutine invert_D_M

c     invert_P_si
c     ##############################################################
      subroutine invert_P_si(y,x)

c     --------------------------------------------------------------
c     Inverts Schur complement
c     --------------------------------------------------------------

        implicit none

c     Call variables

        real(8)   :: x(nn,neq),y(nn,neq)

c     Local variables

        real(8)   :: xc(nn,6),rhsc(nn,6),rhssi(nn,3)

c     Begin program

        if (iout > 0.and.my_rank ==0) then
          write (*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
          write (*,*) 'Solving for small-flow SI system...'
          write (*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        endif

        !Form SI rhs -(Gv+Ldx)
        call formSIrhs(nn,x,y(:,IVX:IVZ),rhssi,igrid)

        !Solve Schur-complement SI system
        if (di == 0d0 .or. (.not.pc_ion_hall)) then

          call cSolver(3,nn,rhssi,x(:,IVX:IVZ)
     .            ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,vol_wgt
     .            ,tol           = pc_tol
     .            ,dg            = v_diag
     .            ,cvrg_tst      = mg_cvrg_test
     .            ,mg_debug      = mg_debug
     .            ,mg_galerkin   = mg_glrkin
     .            ,mg_order_res  = mg_ores_si
     .            ,mg_order_prol = mg_oprol_si
     .            ,mg_coarse_grid_size = mg_coarse_size
     .            ,mg_gm_coarse_solve  = mg_gm_coarse
     .            ,mg_smooth     = sm_type_si
     .            ,mg_vcyc       = mg_vcyc
     .            ,sm_it         = sm_iter_si
     .            ,sm_omega      = sm_omega_si
     .            ,sm_ncolors    = sm_ncolors
     .            ,sm_line_relax = sm_line_relax
     .            ,sm_zebra_relax= sm_zebra_relax_si
     .            )

        else

cc          bcs_hnew(:,1:3) = bcs(:,IVX:IVZ)
cc          bcs_hnew(:,4:6) = bcs(:,IJX:IJZ)

          xc(:,1:3) = x(:,IVX:IVZ)
          xc(:,4:6) = 0d0

          rhsc(:,1:3) = rhssi
          rhsc(:,4:6) = 0d0

          call cSolver(6,nn,rhsc,xc,bcs_hyp
     .            ,igrid,iout,guess,v_mtvc_hyp,vol_wgt
     .            ,tol           = pc_tol
     .            ,dg            = v_diag_hyp
     .            ,mg_debug      = mg_debug
     .            ,mg_galerkin   = mg_glrkin
     .            ,mg_order_res  = mg_ores_si
     .            ,mg_order_prol = mg_oprol_si
     .            ,mg_coarse_grid_size = mg_coarse_size
     .            ,mg_gm_coarse_solve  = mg_gm_coarse
     .            ,mg_smooth     = sm_type_si
     .            ,mg_vcyc       = mg_vcyc
     .            ,sm_it         = sm_iter_si
     .            ,sm_omega      = sm_omega_si
     .            ,sm_ncolors    = sm_ncolors
     .            ,sm_line_relax = sm_line_relax
     .            ,sm_zebra_relax= sm_zebra_relax
     .            )

          x(:,IVX:IVZ) = xc(:,1:3)
        endif

        !Transform velocity field for consistency with NLF
cc        if (alt__eom()) then
cc
cc          do k = 1,nz
cc            do j = 1,ny
cc              do i = 1,nx
cc                ii  = i + nx*(j-1) + nx*ny*(k-1)
cc                x(ii,IVY) = x(ii,IVY)*gmetric%grid(igx)%jac(i,j,k)
cc              enddo
cc            enddo
cc          enddo
cc
cc        endif

c     End program

      end subroutine invert_P_si

c     invert_P_star_si
c     ##############################################################
      subroutine invert_P_star_si(y,x)

c     --------------------------------------------------------------
c     Inverts Schur complement
c     --------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: x(nn,neq),y(nn,neq),rhspc(nn,neq)

c     Begin program

        if (iout > 0.and.my_rank ==0) then
          write (*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
          write (*,*) 'Solving for arbitrary-flow SI system...'
          write (*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        endif

        !Form SI rhs (-Gv-Ldx)
        call formSIrhs(nn,x,y(:,IVX:IVZ),rhspc(:,IVX:IVZ),igrid)

        !Solve Schur-complement SI system
        call cSolver(3,nn,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
     .              ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc_af,vol_wgt
     .              ,tol           = pc_tol
     .              ,dg            = v_diag
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = mg_ores
     .              ,mg_order_prol = mg_oprol
     .              ,mg_coarse_grid_size = mg_coarse_size
     .              ,mg_gm_coarse_solve  = mg_gm_coarse
     .              ,mg_smooth     = sm_type_si
     .              ,mg_vcyc       = mg_vcyc
     .              ,sm_it         = sm_iter_si
     .              ,sm_omega      = sm_omega_si
     .              ,sm_ncolors    = sm_ncolors
     .              ,sm_line_relax = sm_line_relax
     .              ,sm_zebra_relax= sm_zebra_relax
     .              )

        !Post-process solution with M_star
        call postProcessSchur(3,nn,x(:,IVX:IVZ),igrid,bcs(:,IVX:IVZ))

c     End program

      end subroutine invert_P_star_si

c     get_dP
c     ##############################################################
      subroutine get_dP

c     --------------------------------------------------------------
c     Gets momentum update
c     --------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        real(8),pointer,dimension(:,:,:,:) :: rho0,vcnv0

c     Begin program

        rho0  => grho0%grid(igrid)%array
        vcnv0 => gv0  %grid(igrid)%array

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
             ii  = i + nx*(j-1) + nx*ny*(k-1)

             xxx(ii,IVX)=xxx(ii,IVX )*rho0 (i,j,k,1)
     .                  +xxx(ii,IRHO)*vcnv0(i,j,k,1)
             xxx(ii,IVY)=xxx(ii,IVY )*rho0 (i,j,k,1)
     .                  +xxx(ii,IRHO)*vcnv0(i,j,k,2)
             xxx(ii,IVZ)=xxx(ii,IVZ )*rho0 (i,j,k,1)
     .                  +xxx(ii,IRHO)*vcnv0(i,j,k,3)

            enddo
          enddo
        enddo

c     End program

      end subroutine get_dP

c     split_PC
c     ##############################################################
      subroutine split_PC

c     --------------------------------------------------------------
c     Gets momentum update
c     --------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer    :: ii,iii,i,j,k
        real(8)    :: dum(nn,neq)

c     Begin program

c     Predictor step

        call invert_M(yyy,xxx,0)

c     SI step: Deltav --> xxx(:,IVX:IVZ)

        call invert_P_si(yyy,xxx)

c     Corrector step

cc        if (di == 0.or.pc_ion_hall) then
          dum = 0d0
          call find_dx(dum,xxx,0)
cc        else
cc          call find_dx(yyy,xxx,4)
cc        endif

c     End program

      end subroutine split_PC

c     si_PC
c     ##############################################################
      subroutine si_PC(y,dx,version)

c     --------------------------------------------------------------
c     Gets momentum update using a stationary iteration to avoid
c     the M^-1 factor in the Schur complement. The parameter
c     version controls the flavor of stationary iteration:
c        - version = 1 --> Simple diagonal split of M
c        - version = 2 --> Approximation of U to [M (D_M)^-1 U]
c        - version = 3 --> Approximation of L to [L (D_M)^-1 M]
c        - version = 4 --> Experimental, solving both inv(M) steps
c                          but using inv(D_M) in Schur complement
c     --------------------------------------------------------------

        implicit none

c     Call variables

        integer :: version
        real(8) :: y(nn,neq),dx(nn,neq)

c     Local variables

c     Begin program

c     Form (D_m^-1)rx vector ---> dx

        if (version == 2.or.version == 4) then
          call invert_M(y,dx,0)
        else
          call invert_D_M(y,dx)
        endif

c     SI step: Deltav --> dx(:,IVX:IVZ)

        call invert_P_si(y,dx)

c     Solve for dx

        call find_dx(y,dx,version)

c     End program

      end subroutine si_PC

c     star_PC
c     ##############################################################
      subroutine star_PC

c     --------------------------------------------------------------
c     Gets momentum update
c     --------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer :: ii,iii,i,j,k
        real(8) :: dum(nn,neq)

c     Begin program

c     Predictor step

        call invert_M(yyy,xxx,1)

c     SI step: Deltav --> xxx(:,IVX:IVZ)

        call invert_P_star_si(yyy,xxx)

c     Solve for dx

        dum = 0d0
        call find_dx(dum,xxx,4)

c     End program

      end subroutine star_PC

c     norm_test
c     ##############################################################
      subroutine norm_test(r,it)

c     --------------------------------------------------------------
c     Calculates SI iteration norm
c     --------------------------------------------------------------

        implicit none

c     Call variables

        integer :: it
        real(8) :: r(nn,neq)

c     Local variables

        real(8)      :: mag,magi(neq),ratio
        real(8),save :: mag0,mag_n

c     Begin program

        !Global residual
        do ieq=1,neq
          magi(ieq) = dot2(r(:,ieq),r(:,ieq))
        enddo
        mag = sqrt(sum(magi))

        if (it == 1) then
          mag0  = mag
          mag_n = -1d0
        endif

        if (my_rank == 0 .and. iout >= 1) write (*,*)

        ratio = max(mag/mag_n,0d0)

        if (my_rank == 0)
     .    write (*,100) '   SI iter:',it-1,';  Res: ',mag
     .                 ,';  Rel.res.:',mag/mag0,';  Ratio:',ratio

        mag_n = mag

        !Individual residual eqs
        if (iout >=2.and.my_rank==0) then
          do ieq=1,neq
            write (*,*) '     Res. eqn',ieq,':',sqrt(magi(ieq))
          enddo
          write (*,*)
        endif

 100    format (a,i3,a,1p,e10.2,a,e10.2,a,e10.2)

      end subroutine norm_test

#if defined(vec_pot)
c     EMHD_step
c     ###################################################################
      subroutine EMHD_step(b,x)
c     --------------------------------------------------------------------
c     This subroutine performs EMHD step for vector potential, and
c     returns updated RHS to feed into MHD step.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: x(nn,3),b(nn,3)

c     Local variables

      real(8) :: xc(nn,6),rhsc(nn,6)

c     Begin program

c     Transfer initial guess and residual

      xc(:,1:3) = x
      xc(:,4:6) = 0d0

      rhsc(:,1:3) = b
      rhsc(:,4:6) = 0d0

c     Solve EMHD system

      guess = 1

      call cSolver(6,nn,rhsc,xc,bcs_hyp
     .            ,igrid,iout,guess,a_mtvc_hyp,vol_wgt
     .            ,dg            = a_diag_hyp
     .            ,mg_debug      = mg_debug
     .            ,mg_galerkin   = mg_glrkin
     .            ,mg_order_res  = mg_ores
     .            ,mg_order_prol = mg_oprol
     .            ,mg_gm_coarse_solve  = mg_gm_coarse
     .            ,mg_vcyc       = mg_vcyc
     .            ,sm_ncolors    = sm_ncolors
     .            ,sm_line_relax = sm_line_relax
     .            ,sm_zebra_relax= sm_zebra_relax
     .            !Production
     .            ,tol           = pc_tol
     .            ,mg_coarse_grid_size = mg_coarse_size
     .            ,mg_smooth     = sm_type
     .            ,sm_omega      = sm_omega
     .            ,sm_it         = sm_iter
cc     .            !Test
cc     .            ,mg_smooth     = 'gs'
cc     .            ,sm_omega      = 1d0
cc     .            ,sm_it         = 5
cc     .            !Diag
cc     .            ,cvrg_tst=.false.
cccc     .            ,mg_smooth     = 'jb'
cccc     .            ,sm_omega      = 0.5d0
cc     .            ,mg_smooth     = 'gs'
cc     .            ,sm_omega      = 0.7d0
cc     .            ,sm_it         = 1000
cc     .            ,mg_coarse_grid_size = 6 !2^n x 2^n mesh
cc     .            ,tol           = 1d-10
cccc     .            ,mg_smooth     = 'gm'
     .            )
cc      stop

c     Return solution and updated residual for second stage of splitting

      x = xc(:,1:3)

c     End program

      end subroutine EMHD_step
#endif

      end subroutine apply_app_PC

c formSIrhs
c #########################################################################
      subroutine formSIrhs(nn,xxx,yyy,rhs_si,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the rhs for the velocity SI solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer :: nn,igrid
      real(8) :: xxx(nn,neqd),yyy(nn,3),rhs_si(nn,3)

c Local variables

      real(8) :: eom_rhs(nn,3)

c Begin program

c Find density lower block (because we are solving for dv instead of dp)

      call L_rho_v(nn,xxx,eom_rhs,eom_rhs,igrid,vol_wgt)

      rhs_si = yyy - eom_rhs

c Find linearized EOM rhs [-jxB+grad(p)] correction

      call L_y_star(nn,xxx,eom_rhs,igrid,'lfz',vol_wgt)

      rhs_si = rhs_si - eom_rhs

c End program

      end subroutine formSIrhs

c L_y_star
c #########################################################################
      subroutine L_y_star(nn,xxx,eom_rhs,igrid,part,vwgt)
c--------------------------------------------------------------------
c     This subroutine finds the EM part of the rhs for the velocity
c     SI solve.
c--------------------------------------------------------------------

      use precond_variables

      use nlfunction_setup, ONLY: EOM_jxb,EOM_gp

      use imposeBCinterface

      implicit none

c Call variables

      integer    :: nn,igrid
      real(8)    :: xxx(nn,neqd),eom_rhs(nn,3)
      character(3) :: part
      logical    :: vwgt

c Local variables

      integer    :: i,j,k,ii,isig,iig,ivar,ip,im,nx,ny,nz
      real(8)    :: cov(3),cnv(3),idx,idy,idz,dvol,nabla_v(3,3)
     .             ,nabla_vv0(3,3),mul,vlap(3),advc(3),jac
      real(8),allocatable,dimension(:,:,:,:),target :: dj_cnv
     .                                                ,dj_cov
     .                                                ,db_cnv
     .                                                ,db_cov
     .                                                ,dv_cnv
     .                                                ,vzeros
     .                                                ,dp_cnv
     .                                                ,dve_cnv
     .                                                ,dve_cov
     .                                                ,dvef_cnv
     .                                                ,dE_cnv
     .                                                ,drho  
     .                                                ,dprs  
#if defined(vec_pot)
     .                                                ,da_cov                           
#endif

      real(8),pointer,dimension(:,:,:,:)   :: rho0,tmp0,bcov0,bcnv0
     .                                       ,jcnv0,jcov0,vcnv0,vcov0
     .                                       ,heta0

c Begin program

      isig = grid_params%istartp(igrid)

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      if (      part /= 'lfz'
     .    .and. part /= 'adv'
     .    .and. part /= 'all') then
        call pstop ('L_y_star','Unknown character specifier '//part)
      endif

      rho0  => grho0  %grid(igrid)%array
      tmp0  => gtmp0  %grid(igrid)%array
      vcnv0 => gv0    %grid(igrid)%array
      vcov0 => gv0_cov%grid(igrid)%array
      bcnv0 => gb0    %grid(igrid)%array
      bcov0 => gb0_cov%grid(igrid)%array
      jcnv0 => gj0    %grid(igrid)%array
      jcov0 => gj0_cov%grid(igrid)%array
      heta0 => gheta  %grid(igrid)%array

c Allocate variables

      allocate(dj_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,db_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,vzeros(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dp_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dve_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dve_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dE_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,drho  (0:nx+1,0:ny+1,0:nz+1,1)
     .        ,dprs  (0:nx+1,0:ny+1,0:nz+1,1)
#if defined(vec_pot)
     .        ,da_cov(0:nx+1,0:ny+1,0:nz+1,3)
#endif
     .        )

      vzeros = 0d0

#if !defined(vec_pot)
c Find dB* (w/ BCs)

      call fillArray(igrid,IBX,3,xxx(:,IBX:IBZ),db_cnv,bcs(:,IBX:IBZ)
     .           ,order_bc,arr_cov=db_cov,is_vec=.true.,is_cnv=.true.)

#else
c Find dA* (w/ BCs)

      call fillArray(igrid,IAX,3,xxx(:,IAX:IAZ),da_cov,bcs(:,IAX:IAZ)
     .              ,order_bc,is_vec=.true.,is_cnv=.false.)

c Find dB* from dA*

      db_cnv = curl(igrid,da_cov)
      db_cov = XformVector(igrid,db_cnv,.false.)

      !Find covariant components of dj with BCs
      call setBC(IBX,3,nx,ny,nz,db_cnv,db_cov,vzeros,bcs(:,IBX:IBZ)
     .          ,igrid,igrid,igrid,iorder=order_bc)
#endif

c Find dv (w/ Bcs)

      call fillArray(igrid,IVX,3,xxx(:,IVX:IVZ),dv_cnv,bcs(:,IVX:IVZ)
     .              ,order_bc,is_vec=.true.,is_cnv=.true.,arr0=vzeros)

c Find dj* from dB*

      !Find contravariant current (without BCs)
      dj_cnv = curl(igrid,db_cov)
      dj_cov = XformVector(igrid,dj_cnv,.false.)

      !Find covariant components of dj with BCs
      call setBC(IJX,3,nx,ny,nz,dj_cnv,dj_cov,vzeros,bcs(:,IJX:IJZ)
     .          ,igrid,igrid,igrid,iorder=order_bc)

c Find drho, dprs arrays with BCs

      call fillArray(igrid,IRHO,1,xxx(:,IRHO:IRHO),drho
     .              ,bcs(:,IRHO:IRHO),order_bc)
      call fillArray(igrid,ITMP,1,xxx(:,ITMP:ITMP),dprs
     .              ,bcs(:,ITMP:ITMP),order_bc)

c Find dp_cnv

      if (.not.nc_eom_v .and.(part == 'adv'.or.part == 'all')) then
        dp_cnv(:,:,:,1)=rho0(:,:,:,1)*dv_cnv(:,:,:,1)
     .                 +drho(:,:,:,1)* vcnv0(:,:,:,1)
        dp_cnv(:,:,:,2)=rho0(:,:,:,1)*dv_cnv(:,:,:,2)
     .                 +drho(:,:,:,1)* vcnv0(:,:,:,2)
        dp_cnv(:,:,:,3)=rho0(:,:,:,1)*dv_cnv(:,:,:,3)
     .                 +drho(:,:,:,1)* vcnv0(:,:,:,3)
      endif

c Find electron velocity

cc      do ieq=1,3
cc        dj_cnv(:,:,:,ieq) = dj_cnv(:,:,:,ieq)
cc     .                    - jcnv0(:,:,:,ieq)*drho(:,:,:,1)
cc     .                                      /rho0(:,:,:,1)
cc      enddo

      dve_cnv = form_ve(dv_cnv,dj_cnv,rho0(:,:,:,1))

c Find rhs_v'

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx

            ii  = i + nx*(j-1) + nx*ny*(k-1)
            iig = ii + isig - 1

            jac = gmetric%grid(igrid)%jac(i,j,k)

            if (vwgt) then
              dvol = gmetric%grid(igrid)%dlvol(i,j,k)
            else
              dvol = 1d0
            endif

            if (part == 'lfz' .or. part == 'all') then

              if (nc_eom_jxb) then
                cnv = EOM_jxb(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                     ,db_cov,jcov0,.false.,nc_eom_jxb
     .                     ,vol=.false.)
     .               +EOM_jxb(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                     ,bcov0,dj_cov,.false.,nc_eom_jxb
     .                     ,vol=.false.)
              else
                cnv = EOM_jxb(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                     ,db_cnv,bcnv0,.true.,nc_eom_jxb
     .                     ,vol=.false.)
     .               +EOM_jxb(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                     ,bcnv0,db_cnv,.true.,nc_eom_jxb
     .                     ,vol=.false.)
              endif

              cnv = cnv
     .               -EOM_gp(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                     ,rho0(:,:,:,1),dprs(:,:,:,1),nc_eom_gp
     .                     ,vol=.false.)

              !Substract div_pe piece corresponding to (-di*j/rho)
              if (part == 'lfz') then
                cnv = cnv
     .               -EOM_divPi_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                             ,dve_cnv,heta0(:,:,:,1))
              endif

              if (nc_eom_v) then
                cnv = -cnv/rho0(i,j,k,1)
              else
                cnv = -cnv
              endif
            else
              cnv = 0d0
            endif

            if (part == 'adv' .or. part == 'all') then

              if (nc_eom_v) then
                cnv = cnv
     .              +(tensor_nc_advc(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                              ,vcnv0,dv_cnv,min(v_advect,2))
     .               +tensor_nc_advc(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                              ,dv_cnv,vcnv0,min(v_advect,2)))
              else
                cnv = cnv
     .              +(tensor_advc(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,vcnv0,dp_cnv,v_advect)
     .               +tensor_advc(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,dp_cnv,vcnv0,v_advect))
              endif

              !Stress tensors
              cnv = cnv
     .            + EOM_divPi_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                         ,dv_cnv,gvis%grid(igrid)%array(:,:,:,1))

              if (nc_eom_v) cnv = cnv/rho0(i,j,k,1)
            endif

            if (alt__eom()) cnv(2) = cnv(2)*jac

            !Correct rhs_v
            eom_rhs(ii,:) = alpha*dvol*cnv(:)

          enddo
        enddo
      enddo

c Free memory

      nullify(rho0,tmp0,bcov0,bcnv0,jcnv0,jcov0,vcnv0,vcov0)

      deallocate(dj_cnv
     .          ,dj_cov
     .          ,db_cnv
     .          ,db_cov
     .          ,dv_cnv
     .          ,vzeros
     .          ,dp_cnv
     .          ,dve_cnv
     .          ,dve_cov
     .          ,dE_cnv
     .          ,drho  
     .          ,dprs  
#if defined(vec_pot)
     .          ,da_cov
#endif
     .          )

c End program

      end subroutine L_y_star

c correctVfld
c #########################################################################
      subroutine correctVfld(nn,x,crhs,igrid,iout)
c--------------------------------------------------------------------
c     This subroutine corrects the flow field with all-new data.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: nn,igrid,iout
      real(8)    :: x(nn,neqd),crhs(nn,3)

c Local variables

      integer    :: i,j,k,ii,iig,isig,ivar,ip,im,guess,nx,ny,nz
      real(8)    :: nrhs(nn,3),dvol,cov(3),cnv(3),idx,idy,idz
     .             ,nabla_vv0(3,3),mul,vlap(3),advc(3)

      external   :: v_mtvc2

c Begin program

      isig = grid_params%istartp(igrid)

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

c Find SI rhs (but with "complete" quantities)

      call formSIrhs(nn,x,crhs,nrhs,igrid)

c Solve for advection-diffusion part

      guess = 1

      call cSolver(3,nn,nrhs,x(:,IVX:IVZ)
     .            ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc2,vol_wgt
     .            ,tol           = pc_tol
     .            ,dg            = v_diag2
     .            ,mg_debug      = mg_debug
     .            ,mg_galerkin   = mg_glrkin
     .            ,mg_order_res  = mg_ores
     .            ,mg_order_prol = mg_oprol
     .            ,mg_coarse_grid_size = mg_coarse_size
     .            ,mg_gm_coarse_solve  = mg_gm_coarse
     .            ,mg_smooth     = sm_type
     .            ,mg_vcyc       = mg_vcyc
     .            ,sm_it         = sm_iter
     .            ,sm_omega      = sm_omega
     .            ,sm_ncolors    = sm_ncolors
     .            ,sm_line_relax = sm_line_relax
     .            ,sm_zebra_relax= sm_zebra_relax
     .            )

      return

c Correction w/o solve

      call L_y_star(nn,x,nrhs,igrid,'all',vol_wgt)

      nrhs = crhs - nrhs

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            if (vol_wgt) then
              dvol = gmetric%grid(igrid)%dlvol(i,j,k)
            else
              dvol = 1d0
            endif
            x(ii,IVX:IVZ) = dt*nrhs(ii,:)/dvol/cnp
          enddo
        enddo
      enddo

c End program

      end subroutine correctVfld

c correctBfld
c #########################################################################
      subroutine correctBfld(nn,igrid,crhs,x,iout)
c--------------------------------------------------------------------
c     This subroutine corrects the B-field update with all-new info,
c     in a divergence-free manner.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use xdraw_io

      implicit none

c Call variables

      integer    :: nn,igrid,iout
      real(8)    :: x(nn,neqd),crhs(nn,3)

c Local variables

      integer    :: i,j,k,ii,isig,ivar,nx,ny,nz,bcnd(6,3)
      real(8)    :: dvol,cov(3),cnv(3)

      real(8),allocatable,dimension(:,:,:,:) ::  db_cnv
     .                                          ,db_cov
     .                                          ,dj_cnv
     .                                          ,dj_cov
     .                                          ,de_cov
     .                                          ,de_cnv
     .                                          ,vzeros

      real(8),pointer,dimension(:,:,:,:) :: v0cnv

c Begin program

      if (my_rank == 0 .and. iout > 0)
     .     write (*,*) 'Correcting B-field...'

      isig = grid_params%istartp(igrid)

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      allocate ( db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,db_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dj_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,de_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,de_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,vzeros(0:nx+1,0:ny+1,0:nz+1,3))

c Find dB

      call fillArray(igrid,IBX,3,x(:,IBX:IBZ),db_cnv,bcs(:,IBX:IBZ)
     .              ,order_bc,arr_cov=db_cov,is_vec=.true.
     .              ,is_cnv=.true.)

c Find dj* from dB*

      dj_cnv = curl(igrid,db_cov)
      dj_cov = XformVector(igrid,dj_cnv,.false.)

      !Find covariant components of dj with BCs
      call setBC(IJX,3,nx,ny,nz,dj_cnv,dj_cov,vzeros,bcs(:,IJX:IJZ)
     .          ,igrid,igrid,igrid,iorder=order_bc)

c Find dE* (electric field update)

      v0cnv => gv0%grid(igrid)%array

      de_cov =-crossProduct(igrid,v0cnv,db_cnv,.true.)

      de_cov(:,:,:,1) = de_cov(:,:,:,1)
     .     + geta%grid(igrid)%array(:,:,:,1)*dj_cov(:,:,:,1)
      de_cov(:,:,:,2) = de_cov(:,:,:,2)
     .     + geta%grid(igrid)%array(:,:,:,1)*dj_cov(:,:,:,2)
      de_cov(:,:,:,3) = de_cov(:,:,:,3)
     .     + geta%grid(igrid)%array(:,:,:,1)*dj_cov(:,:,:,3)

cc      do k = 0,nz+1
cc        do j = 0,ny+1
cc          do i = 0,nx+1
cc            cov(1) = (v0cnv(i,j,k,2)*db_cnv(i,j,k,3)
cc     .               -v0cnv(i,j,k,3)*db_cnv(i,j,k,2))
cc            cov(2) = (v0cnv(i,j,k,3)*db_cnv(i,j,k,1)
cc     .               -v0cnv(i,j,k,1)*db_cnv(i,j,k,3))
cc            cov(3) = (v0cnv(i,j,k,1)*db_cnv(i,j,k,2)
cc     .               -v0cnv(i,j,k,2)*db_cnv(i,j,k,1))
cc
cc            de_cov(i,j,k,:)= -cov/gmetric%grid(igrid)%jac(i,j,k)
cc     .                       +geta%grid(igrid)%array(i,j,k,1)
cc     .                       *dj_cov(i,j,k,:)
cccc            de_cnv(i,j,k,:) = matmul(gmetric%grid(igrid)%gsub(i,j,k,:,:)
cccc     .                              ,de_cov(i,j,k,:))
cc          enddo
cc        enddo
cc      enddo

cc      write (*,*) 'DIAG -- correctBfld'
cc      open(unit=110,file='divclean.bin',form='unformatted'
cc     .    ,status='replace')
cc      call contour(de_cov(:,:,1,1),nx+2,ny+2
cc     .                ,grid_params%gxmin,grid_params%gxmax
cc     .                ,grid_params%gymin,grid_params%gymax
cc     .                ,0,110)
cc      call contour(de_cov(:,:,1,2),nx+2,ny+2
cc     .                ,grid_params%gxmin,grid_params%gxmax
cc     .                ,grid_params%gymin,grid_params%gymax
cc     .                ,1,110)
cc      call contour(de_cov(:,:,1,3),nx+2,ny+2
cc     .                ,grid_params%gxmin,grid_params%gxmax
cc     .                ,grid_params%gymin,grid_params%gymax
cc     .                ,1,110)
cc      close(110)

cc      !Impose topological BCs
cc      bcnd(:,1) = bcond
cc      bcnd(:,2) = bcond
cc      bcnd(:,3) = bcond
cc
cc      call setBC(IJX,3,nx,ny,nz,de_cnv,de_cov,vzeros,bcnd
cc     .          ,igrid,igrid,igrid,iorder=order_bc
cc     .          ,is_cnv=.false.)

c Find correction

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            if (vol_wgt) then
              dvol = gmetric%grid(igrid)%dlvol(i,j,k)
            else
              dvol = 1d0
            endif
            x(ii,IBX:IBZ) = dt/cnp*(crhs(ii,:)/dvol
     .             -alpha*curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,de_cov))
          enddo
        enddo
      enddo

c End program

      deallocate(db_cnv
     .          ,db_cov
     .          ,dj_cnv
     .          ,dj_cov
     .          ,de_cov
     .          ,de_cnv
     .          ,vzeros)

      end subroutine correctBfld

c preprocessBflRHS
c ###################################################################
      subroutine preprocessBflRHS(nn,neq,yyy,igrid)
c--------------------------------------------------------------------
c     This subroutine preprocesses B rhs to add Hall correction.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: igrid,nn,neq

      real(8)    :: yyy(nn,neq)

c Local variables

      integer :: i,j,k,ii,nx,ny,nz,ieq

      real(8) :: cnv(3),dvol

      real(8),allocatable,dimension(:,:,:,:) :: dgv_cnv,dgv_cov

      real(8),pointer,dimension(:,:,:,:) :: rho0

c Begin program

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      rho0 => grho0%grid(igrid)%array

      allocate(dgv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dgv_cov(0:nx+1,0:ny+1,0:nz+1,3))

c Store velocity RHS in array format w/ BCs

      call fillArray(igrid,IVX,3,yyy(:,IVX:IVZ),dgv_cnv,bcs(:,IVX:IVZ)
     .              ,order_bc,arr_cov=dgv_cov,is_vec=.true.
     .              ,is_cnv=.true.)

c Divide by density and volume

      if (vol_wgt) then
        do ieq=1,3
          dgv_cov(:,:,:,ieq) = dgv_cov(:,:,:,ieq)
     .                        /gmetric%grid(igrid)%dlvol
        enddo
      endif

      if (.not.nc_eom_v) then
        do ieq=1,3
          dgv_cov(:,:,:,ieq) = dgv_cov(:,:,:,ieq)/rho0(:,:,:,1)
        enddo
      endif

c Find curl

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            if (vol_wgt) then   !Residuals are volume-weighed by default
              dvol = gmetric%grid(igrid)%dlvol(i,j,k)
            else
              dvol = 1d0
            endif
            yyy(ii,IBX:IBZ) = yyy(ii,IBX:IBZ)
     .                   + di*dvol*curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                                 ,dgv_cov)
          enddo
        enddo
      enddo

      deallocate(dgv_cnv,dgv_cov)

c End program

      end subroutine preprocessBflRHS

c preprocessAflRHS
c ###################################################################
      subroutine preprocessAflRHS(nn,neq,yyy,igrid)
c--------------------------------------------------------------------
c     This subroutine preprocesses B rhs to add Hall correction.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer :: igrid,nn,neq

      real(8) :: yyy(nn,neq)

c Local variables

      integer :: i,j,k,ii,nx,ny,nz
      real(8) :: cov(3),cnv(3)
      real(8),pointer,dimension(:,:,:,:) :: rho0

c Begin program

      rho0 => grho0%grid(igrid)%array

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            cnv = yyy(ii,IVX:IVZ)
            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .           ,cov(1),cov(2),cov(3)
     .           ,cnv(1),cnv(2),cnv(3),.false.)

            if (nc_eom_v) cov = cov/rho0(i,j,k,1)

            yyy(ii,IAX:IAZ) = yyy(ii,IAX:IAZ) + di*cov
          enddo
        enddo
      enddo

c End program

      end subroutine preprocessAflRHS

c divClean
c ###################################################################
      subroutine divClean(iout,nn,xb,igrid)
c--------------------------------------------------------------------
c     This subroutine divergence-cleans B by integrating a vector
c     potential and the differentiating it. It takes B in cartesian
c     or curvilinear and returns Bcnv.
c--------------------------------------------------------------------

      use B_tools

      implicit none

c Call variables

      integer :: iout,nn,igrid
      real(8) :: xb(nn,3)

c Local variables

      integer :: i,j,k,ii,igx,igy,igz,nx,ny,nz,its
      real(8),allocatable,dimension(:,:,:,:) :: db_cnv,db_cov

c Begin program

      if (iout > 0 .and. my_rank == 0)
     .   write (*,*) '  --> Performing divergence cleaning in PC...'

      igx = igrid
      igy = igrid
      igz = igrid

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      allocate(db_cnv (0:nx+1,0:ny+1,0:nz+1,3)
     .        ,db_cov (0:nx+1,0:ny+1,0:nz+1,3))

c Find dB (w/ BCs)

      call fillArray(igrid,IBX,3,xb,db_cnv,bcs(:,IBX:IBZ)
     .              ,order_bc,is_vec=.true.,is_cnv=.true.
     .              ,arr_cov=db_cov)

c Solve dB*=-curl(lap^(-1)(curl(dB)))

      its = 0

c$$$      call curl_div_clean(igx,db_cnv,pc_debug)

      call vlap_div_clean(iout,nx,ny,nz,igx,igy,igz
     .                   ,db_cnv,dbcov=db_cov,pc_call=.true.
     .                   ,diag_plots=pc_debug,gm_it=its)

      it_divcl = it_divcl + its

c Map array to vector

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            xb(ii,:)=db_cnv(i,j,k,:)
          enddo
        enddo
      enddo

      deallocate(db_cnv,db_cov)

c End program

      end subroutine divClean

c postProcessSchur
c ###################################################################
      subroutine postProcessSchur(neq,nn,xxx,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine postprocess SI solution to find final
c     solution to Schur complement system.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use imposeBCinterface

      implicit none

c Call variables

      integer    :: nn,igrid,neq,bcnd(6,neq)
      real(8)    :: xxx(nn,neq)

c Local variables

      integer    :: ieq,ii,ntot
      real(8)    :: xi(nn*neq),bi(nn*neq)

c Begin program

      ntot = nn*neq

c Gather variables

      do ii=1,nn
        do ieq=1,3
          xi(neq*(ii-1)+ieq) = xxx(ii,ieq)
        enddo
      enddo

c Perform matvec with M_star

      call mstar_mtvc(0,3,nn*neq,xi,bi,igrid,bcnd)

c Scatter solution

      do ii=1,nn
        do ieq=1,neq
          xxx(ii,ieq) = bi(neq*(ii-1)+ieq)
        enddo
      enddo

c End program

      end subroutine postProcessSchur

#endif
