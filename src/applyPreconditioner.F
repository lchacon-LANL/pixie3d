c applyPreconditioner
c###################################################################
      subroutine applyPreconditioner(ntot,y,x,iout)

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

      integer    :: ntot,iout
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer    :: ii,iii,igrid,ntotp,guess
     .             ,igx,igy,igz,nx,ny,nz

      real(8)    :: xxx(ntot/neqd,neqd),yyy(ntot/neqd,neqd)
     .             ,rhspc(ntot/neqd,neqd),rr(ntot),dxk(ntot)
     .             ,dvol,car(3),cnv(3),omg(neqd)

      real(8),allocatable,dimension(:,:,:,:) :: dv_cnv,dv_cov

      logical    :: pc_diag

c Externals

      external   tmp_mtvc,rho_mtvc,b_mtvc,v_mtvc,v_mtvc2,test_mtvc
     .          ,v_mtvc3,hall_mtvc,a_mtvc,b_mtvc_hyp,a_mtvc_hyp

c Begin program

      igx = 1
      igy = 1
      igz = 1

      igrid = igx

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      ntotp = ntot/neqd

      if (ntotp /= nx*ny*nz) then
        write (*,*)'Grid sizes do not agree in applyPreconditioner'
        write (*,*)'Aborting...'
        stop
      endif

c *******************************************************************
c     Identity preconditioner
c *******************************************************************

      select case(precon)
      case('id')

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              do ieq=1,neqd
                iii = ieq + neqd*(ii-1)
                if (vol_wgt) then !Residuals are volume-weighed by default
                  dvol = gmetric%grid(igx)%dvol(i,j,k)
                else
                  dvol = 1d0
                endif
                x(iii) = y(iii)*dt/dvol/cnp
              enddo
            enddo
          enddo
        enddo

        return

c *******************************************************************
c     Semi-implicit preconditioner version I (w/o flow)
c *******************************************************************

      case('s1')

        xxx = 0d0

c     Scatter residuals

        call scatter(yyy,y)

        guess = 0

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3))

c     Preprocess B/A rhs

cc        if (di > 0d0 .and. (.not.hall_new)) then
cc#if !defined(vec_pot)
cc          call preprocessBflRHS(ntotp,neqd,yyy,igrid)
cc#else
cc          call preprocessAflRHS(ntotp,neqd,yyy,igrid)
cc#endif
cc        endif

c     Split PC

        call split_PC

c     Postprocessing of velocity -> momentum

        if (v_solve)
     .       call correctVfldRHS(ntotp,xxx,yyy(:,IVX:IVZ),igrid,iout)

        if (.not.nc_eom_v) call get_dP

c     Postprocessing of B

#if !defined(vec_pot)
cc        if (eigen) call divClean(iout,ntotp,xxx(:,IBX:IBZ),igrid)
#endif

c     Gather solution in "x"

        call gather(x,xxx)

c *******************************************************************
c     Semi-implicit preconditioner version II (with flow)
c *******************************************************************

      case('s2')

c     Set up vectors

        x   = 0d0
        xxx = 0d0

        guess = 0   !Equations are in residual form

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3))

c     Preprocess B rhs

        if (di > 0d0) then
          call preprocessBflRHS(ntotp,neqd,y,igrid)
        endif

#if !defined(petsc_c) && !defined(samrai)

        do si_it=1,precpass

c       Form residual vector rr=y-Ax

          if (si_it > 1) then
            call matrixFreeMatVec(ntot,x,rr)
            rr = y - rr
          else
            rr = y
          endif

          if (iout >= 0) call norm_test

#else
          si_it = mod(si_it,precpass)+1

          rr = y
#endif

c       Map residual to individual components (rr --> yyy)

          call scatter(yyy,rr)

c       SI step

          call si_PC(1)

c       Postprocessing of velocity -> momentum

          if (.not.nc_eom_v) call get_dP

c       Postprocessing of B

#if !defined(vec_pot)
          if (eigen) call divClean(iout,ntotp,xxx(:,IBX:IBZ),igrid)
#endif

c       Map solution xxx to x (jacobi iterate)

#if !defined(petsc_c)  && !defined(samrai)
          omg = 1d0
          omg(IVX:IVZ) = 1d0

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                do ieq=1,neqd
                  iii = ieq + neqd*(ii-1)
                  x(iii) = x(iii) + omg(ieq)*xxx(ii,ieq)
                enddo
              enddo
            enddo
          enddo

        enddo

c     Final SI norm calculation

        if (iout >= 0 .and. si_it > precpass) then
          call matrixFreeMatVec(ntot,x,rr)
          rr = y - rr
          call norm_test
        endif

#else

        call gather(x,xxx)

#endif

c *******************************************************************
c     Semi-implicit preconditioner version III (with flow)
c *******************************************************************

      case('s3')

c     Set up vectors

        xxx = 0d0

c     Scatter residuals

        call scatter(yyy,y)

        guess = 0

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3))

c     Preprocess B/A rhs

cc        if (di > 0d0 .and. (.not.hall_new)) then
cc#if !defined(vec_pot)
cc          call preprocessBflRHS(ntotp,neqd,yyy,igrid)
cc#else
cc          call preprocessAflRHS(ntotp,neqd,yyy,igrid)
cc#endif
cc        endif

c     Commuted Schur complement PC

        call star_PC

c     Postprocessing of velocity -> momentum

        if (.not.nc_eom_v) call get_dP

c     Postprocessing of B

#if !defined(vec_pot)
cc        if (divCln) call divClean(iout,ntotp,xxx(:,IBX:IBZ),igrid)
#endif

c     Gather solution in "x"

        call gather(x,xxx)

      case default

        call pstop('applyPC','Cannot use PC option'//precon)

      end select

c     Diagnostics

        if (debug) then
          if (debug_it > 0) then
            if (jit == debug_it) call PCdiagnostics
          else
            if (my_rank == 0) then
              write (*,*) 'Run PC diagnostics (t or f)?'
              read(*,*) pc_diag
            endif

#if defined(petsc) || defined(petsc_c)
            call MPI_Bcast(pc_diag,1,MPI_LOGICAL,0,MPI_COMM_WORLD
     .                    ,mpierr)
#endif
            if (pc_diag) call PCdiagnostics
          endif
        endif

c     Deallocate variables

        deallocate(dv_cnv,dv_cov)

c End program

      contains

c     find_dx
c     ##############################################################
      subroutine find_dx(x,rr,flag)

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

        integer    :: flag
        real(8)    :: x(ntot/neqd,neqd),yyy(ntot/neqd,neqd),rr(ntot)

c     Local variables

        real(8)    :: zz(ntot/neqd,neqd)

c     Begin program

        guess = 1

        !Correct dx rhs

        !Solve for dx
        select case(flag)
        case(0)

          call find_Udv(x,zz)

          if (solve_rho) x(:,IRHO)=x(:,IRHO)-dt*alpha/cnp*zz(:,IRHO)

          x(:,ITMP)   = x(:,ITMP)    - dt*alpha/cnp*zz(:,ITMP)

#if !defined(vec_pot)
          x(:,IBX:IBZ)= x(:,IBX:IBZ) - dt*alpha/cnp*zz(:,IBX:IBZ)
#else
          x(:,IAX:IAZ)= x(:,IAX:IAZ) - dt*alpha/cnp*zz(:,IAX:IAZ)
#endif
        case(1)

          call correct_dx_rhs(x,rr,yyy)
          call invert_D_M(yyy,x)

        case(2)
          zz = 0d0

          call correct_dx_rhs(x,zz,yyy)

          call invert_D_M(yyy,zz)

          if (solve_rho)  x(:,IRHO) = x(:,IRHO) + zz(:,IRHO)
          x(:,ITMP) = x(:,ITMP) + zz(:,ITMP)
#if !defined(vec_pot)
          x(:,IBX:IBZ) = x(:,IBX:IBZ) + zz(:,IBX:IBZ)
#else
          x(:,IAX:IAZ) = x(:,IAX:IAZ) + zz(:,IAX:IAZ)
#endif
        case(3)

          call correct_dx_rhs(x,rr,yyy)
          call invert_M(yyy,x)

        case(4)
          zz = 0d0

          call correct_dx_rhs(x,zz,yyy)

          call invert_M(yyy,zz)

          if (solve_rho) x(:,IRHO) = x(:,IRHO) + zz(:,IRHO)
          x(:,ITMP) = x(:,ITMP) + zz(:,ITMP)
#if !defined(vec_pot)
          x(:,IBX:IBZ) = x(:,IBX:IBZ) + zz(:,IBX:IBZ)
#else
          x(:,IAX:IAZ) = x(:,IAX:IAZ) + zz(:,IAX:IAZ)
#endif

        case default

          call pstop('find_dx'
     .              ,'Stationary iteration type not implemented')

        end select

c     End program

      end subroutine find_dx

c     correct_dx_rhs
c     ##############################################################
      subroutine correct_dx_rhs(x,rx,yyy)

c     --------------------------------------------------------------
c     Finds the operator U(dv) and corrects the rhs of dx=(drho,dB,dT)
c     equation.
c     --------------------------------------------------------------

        implicit none

c     Call variables

        real(8)   :: x(ntot/neqd,neqd),yyy(ntot/neqd,neqd)
     .              ,rhs(ntot/neqd,neqd)
     .              ,rx(ntot)

c     Local variables

c     Begin program

        !Find U(dv)
        call find_Udv(x,rhs)

        !Find corrected residuals (Gx-U.dv)
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              iii = neqd*(ii-1)

              if (vol_wgt) then !Residuals are volume-weighed by default
                dvol = gmetric%grid(igx)%dvol(i,j,k)
              else
                dvol = 1d0
              endif

              yyy(ii,IRHO)   =rx(iii+IRHO)-alpha*dvol*rhs(ii,IRHO)
#if !defined(vec_pot)
              yyy(ii,IBX:IBZ)=rx(iii+IBX:iii+IBZ)
     .                                    -alpha*dvol*rhs(ii,IBX:IBZ)
#else
              yyy(ii,IAX:IAZ)=rx(iii+IAX:iii+IAZ)
     .                                    -alpha*dvol*rhs(ii,IAX:IAZ)
#endif
              yyy(ii,ITMP)   =rx(iii+ITMP)-alpha*dvol*rhs(ii,ITMP)

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

        use local_BCS_variables, ONLY: pinch_flow

        implicit none

c     Call variables

        real(8) :: x(ntot/neqd,neqd),rhs(ntot/neqd,neqd)

c     Local variables

        real(8) :: db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .            ,db_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .            ,dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .            ,dv0cnv(0:nx+1,0:ny+1,0:nz+1,3)

c     Begin program

c     Find dB* (w/ BCs)

        if (pinch_flow) then
          call fillArray(igrid,IBX,3,x(:,IBX:IBZ),db_cnv,bcs(:,IBX:IBZ)
     .                ,2,arr_cov=db_cov,is_vec=.true.,is_cnv=.true.
     .                ,arr0=gb0%grid(igrid)%array)
        endif

c     Store velocity solution in array format --> dv_cnv

        dv0cnv = 0d0

        call ExB_flow_bc(igrid,nx,ny,nz,bcnv,bcov,dv0cnv,dbcov=db_cov)

        call fillArray(igrid,IVX,3,x(:,IVX:IVZ),dv_cnv,bcs(:,IVX:IVZ)
     .                ,2,arr_cov=dv_cov,is_vec=.true.,is_cnv=.true.
     .                ,arr0=dv0cnv)

c     Correct rx to find Deltax (correction for rho, B, T)

        !Find corrections to rho,B,tmp right hand sides (Udv --> rhs)
        call correctRhoRHS(nx,ny,nz,dv_cnv       ,rhs(:,IRHO)   ,igrid)
        call correctTmpRHS(nx,ny,nz,dv_cnv       ,rhs(:,ITMP)   ,igrid)
#if !defined(vec_pot)
        call correctBflRHS(nx,ny,nz,dv_cnv,dv_cov,rhs(:,IBX:IBZ),igrid)
#else
        call correctAflRHS(nx,ny,nz,dv_cnv,dv_cov,rhs(:,IAX:IAZ),igrid)
#endif

c     End program

      end subroutine find_Udv

c     invert_M
c     ##############################################################
      subroutine invert_M(rhs,x)

c     --------------------------------------------------------------
c     Inverts r,B,T diagonal blocks. 
c     --------------------------------------------------------------

        implicit none

c     Call variables

        real(8)    :: x(ntot/neqd,neqd),rhs(ntot/neqd,neqd),dvol

c     Local variables

        real(8)    :: xc(ntot/neqd,6),rhsc(ntot/neqd,6)

c     Begin program

cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc              if (vol_wgt) then !Residuals are volume-weighed by default
cc                dvol = gmetric%grid(igx)%dvol(i,j,k)
cc              else
cc                dvol = 1d0
cc              endif
cc
cc              x(ii,IRHO)    = dt*rhs(ii,IRHO)/cnp/dvol
cc              x(ii,ITMP)    = dt*rhs(ii,ITMP)/cnp/dvol
cc#if !defined(vec_pot)
cc              x(ii,IBX:IBZ) = dt*rhs(ii,IBX:IBZ)/cnp/dvol
cc#else
cc              x(ii,IAX:IAZ) = dt*rhs(ii,IAX:IAZ)/cnp/dvol
cc#endif               
cc            enddo
cc          enddo
cc        enddo
cc        return

        !Temperature
        if (gamma > 1d0) then
          call cSolver(1,ntotp,rhs(:,ITMP),x(:,ITMP),bcs(:,ITMP)
     .              ,igrid,iout,guess,tmp_mtvc,vol_wgt
     .              ,tol           = mgtol
     .              ,dg            = tmp_diag
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = ores
     .              ,mg_order_prol = oprol
     .              ,mg_coarse_grid_size = mg_coarse_size
     .              ,mg_gm_coarse_solve  = gm_coarse
     .              ,mg_smooth     = smooth
     .              ,mg_vcyc       = maxvcyc
     .              ,sm_it         = nsweep
     .              ,sm_omega      = omeg_smth
     .              ,sm_ncolors    = 2  !Star stencil
     .              ,sm_line_relax = line_relax
     .              ,sm_zebra_relax= zebra_relax
     .              )
        endif

        !Density
        if (solve_rho) then
          call cSolver(1,ntotp,rhs(:,IRHO),x(:,IRHO),bcs(:,IRHO)
     .              ,igrid,iout,guess,rho_mtvc,vol_wgt
     .              ,tol           = mgtol
     .              ,dg            = rho_diag
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = ores
     .              ,mg_order_prol = oprol
     .              ,mg_coarse_grid_size = mg_coarse_size
     .              ,mg_gm_coarse_solve  = gm_coarse
     .              ,mg_smooth     = smooth
     .              ,mg_vcyc       = maxvcyc
     .              ,sm_it         = nsweep
     .              ,sm_omega      = omeg_smth
     .              ,sm_ncolors    = 2  !Star stencil
     .              ,sm_line_relax = line_relax
     .              ,sm_zebra_relax= zebra_relax
     .              )
        endif

#if !defined(vec_pot)
        !Magnetic field
        if (inv_B) then

          if (si_car) call XtoCar(igrid,rhs(:,IBX:IBZ))

          if (di*heta == 0d0) then
            call cSolver(3,ntotp,rhs(:,IBX:IBZ),x(:,IBX:IBZ)
     .              ,bcs(:,IBX:IBZ),igrid,iout,guess,b_mtvc,vol_wgt
     .              ,tol           = mgtol
     .              ,dg            = b_diag
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = ores
     .              ,mg_order_prol = oprol
     .              ,mg_coarse_grid_size = mg_coarse_size
     .              ,mg_gm_coarse_solve  = gm_coarse
     .              ,mg_smooth     = smooth
     .              ,mg_vcyc       = maxvcyc
     .              ,sm_it         = nsweep
     .              ,sm_omega      = omeg_smth
     .              ,sm_ncolors    = ncolors
     .              ,sm_line_relax = line_relax
     .              ,sm_zebra_relax= zebra_relax
     .              )
          else
            xc(:,1:3) = x(:,IBX:IBZ)
            xc(:,4:6) = 0d0

            rhsc(:,1:3) = rhs(:,IBX:IBZ)
            rhsc(:,4:6) = 0d0

            call cSolver(6,ntotp,rhsc,xc,bcs_hyp
     .              ,igrid,iout,guess,b_mtvc_hyp,vol_wgt
     .              ,tol           = mgtol
     .              ,dg            = b_diag_hyp
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = ores
     .              ,mg_order_prol = oprol
     .              ,mg_coarse_grid_size = mg_coarse_size
     .              ,mg_gm_coarse_solve  = gm_coarse
     .              ,mg_smooth     = smooth
     .              ,mg_vcyc       = maxvcyc
     .              ,sm_it         = nsweep
     .              ,sm_omega      = omeg_smth
     .              ,sm_ncolors    = ncolors
     .              ,sm_line_relax = line_relax
     .              ,sm_zebra_relax= zebra_relax
     .              )

            x(:,IBX:IBZ) = xc(:,1:3)
          endif

          if (si_car) call XtoCurv(igrid,x(:,IBX:IBZ))

cc          if (divCln) call divClean(iout,ntotp,x(:,IBX:IBZ),igrid)

        else

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                if (vol_wgt) then
                  dvol = gmetric%grid(igx)%dvol(i,j,k)
                else
                  dvol = 1d0
                endif

                x(ii,IBX:IBZ) = dt*rhs(ii,IBX:IBZ)/cnp/dvol
              enddo
            enddo
          enddo

        endif
#else
        !Vector potential
        if (inv_B) then

          if (di == 0d0 .or. hall_new) then
            if (si_car) call XtoCar(igrid,rhs(:,IAX:IAZ))

            call cSolver(3,ntotp,rhs(:,IAX:IAZ),x(:,IAX:IAZ)
     .                ,bcs(:,IAX:IAZ),igrid,iout,guess,a_mtvc,vol_wgt
     .                ,tol           = mgtol
     .                ,dg            = a_diag
     .                ,mg_debug      = mg_debug
     .                ,mg_galerkin   = mg_glrkin
     .                ,mg_order_res  = ores
     .                ,mg_order_prol = oprol
     .                ,mg_coarse_grid_size = mg_coarse_size
     .                ,mg_gm_coarse_solve  = gm_coarse
     .                ,mg_smooth     = smooth
     .                ,mg_vcyc       = maxvcyc
     .                ,sm_it         = nsweep
     .                ,sm_omega      = omeg_smth
     .                ,sm_ncolors    = ncolors
     .                ,sm_line_relax = line_relax
     .                ,sm_zebra_relax= zebra_relax
     .                )
          else
            xc(:,1:3) = x(:,IAX:IAZ)
            xc(:,4:6) = 0d0

            call formArhs_hall(ntotp,x,rhs(:,IAX:IAZ),rhsc(:,1:3),igrid)
cc            rhsc(:,1:3) = rhs(:,IAX:IAZ)
            rhsc(:,4:6) = 0d0

            call cSolver(6,ntotp,rhsc,xc,bcs_hyp
     .              ,igrid,iout,guess,a_mtvc_hyp,vol_wgt
     .              ,tol           = mgtol
     .              ,dg            = a_diag_hyp
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = ores
     .              ,mg_order_prol = oprol
     .              ,mg_coarse_grid_size = mg_coarse_size
     .              ,mg_gm_coarse_solve  = gm_coarse
     .              ,mg_smooth     = smooth
     .              ,mg_vcyc       = maxvcyc
     .              ,sm_it         = nsweep
     .              ,sm_omega      = omeg_smth
     .              ,sm_ncolors    = ncolors
     .              ,sm_line_relax = line_relax
     .              ,sm_zebra_relax= zebra_relax
cc     .              ,cvrg_tst=.true.
cc     .              ,mg_smooth     = 'jb'
cc     .              ,sm_it         = 100
cc     .              ,sm_omega      = 1d0
cc     .              ,mg_smooth     = 'gm'
     .              )

cc           stop
            x(:,IAX:IAZ) = xc(:,1:3)
          endif

          if (si_car) call XtoCurv(igrid,x(:,IAX:IAZ))

        else

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                if (vol_wgt) then
                  dvol = gmetric%grid(igx)%dvol(i,j,k)
                else
                  dvol = 1d0
                endif

                x(ii,IAX:IAZ) = dt*rhs(ii,IAX:IAZ)/cnp/dvol
              enddo
            enddo
          enddo

        endif
#endif

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

        real(8)   :: x(ntot/neqd,neqd),rhs(ntot/neqd,neqd)

c     Local variables

c     Begin program

        call diagonalScaling(1,ntotp,rho_diag,rhs(:,IRHO)
     .                      ,x(:,IRHO)   ,igrid)
        call diagonalScaling(1,ntotp,tmp_diag,rhs(:,ITMP)
     .                      ,x(:,ITMP)   ,igrid)
#if !defined(vec_pot)
        call diagonalScaling(3,ntotp,b_diag  ,rhs(:,IBX:IBZ)
     .                      ,x(:,IBX:IBZ),igrid)
#else
        call diagonalScaling(3,ntotp,a_diag  ,rhs(:,IAX:IAZ)
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

        real(8)   :: x(ntot/neqd,neqd),y(ntot/neqd,neqd)

c     Begin program

        !Form SI rhs -(Gv+Ldx)
        call formSIrhs(ntotp,x,y(:,IVX:IVZ),rhspc(:,IVX:IVZ),igrid)

cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc              if (vol_wgt) then !Residuals are volume-weighed by default
cc                dvol = gmetric%grid(igx)%dvol(i,j,k)
cc              else
cc                dvol = 1d0
cc              endif
cc
cc              x(ii,IVX:IVZ) = dt*rhspc(ii,IVX:IVZ)/cnp/dvol
cccc              x(ii,IVX:IVZ) = dt*y(ii,IVX:IVZ)/cnp/dvol
cc            enddo
cc          enddo
cc        enddo
cc        return

        !Solve Schur-complement SI system
        call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
     .              ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,vol_wgt
     .              ,tol           = mgtol
     .              ,dg            = v_diag
     .              ,cvrg_tst      = cvrg_test
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = ores
     .              ,mg_order_prol = oprol
     .              ,mg_coarse_grid_size = mg_coarse_size
     .              ,mg_gm_coarse_solve  = gm_coarse
     .              ,mg_smooth     = smooth
     .              ,mg_vcyc       = maxvcyc
     .              ,sm_it         = nsweep
     .              ,sm_omega      = omeg_smth
     .              ,sm_ncolors    = ncolors
     .              ,sm_line_relax = line_relax
     .              ,sm_zebra_relax= zebra_relax
     .              )

        if (si_car) call XtoCurv(igrid,x(:,IVX:IVZ))

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

        real(8)   :: x(ntot/neqd,neqd),y(ntot/neqd,neqd)

c     Begin program

        !Form SI rhs (-Gv-Ldx)
        call formSIrhs(ntotp,x,y(:,IVX:IVZ),rhspc(:,IVX:IVZ),igrid)

        !Solve Schur-complement SI system
        call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
     .              ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,vol_wgt
     .              ,tol           = mgtol
     .              ,dg            = v_diag
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = ores
     .              ,mg_order_prol = oprol
     .              ,mg_coarse_grid_size = mg_coarse_size
     .              ,mg_gm_coarse_solve  = gm_coarse
     .              ,mg_smooth     = smooth
     .              ,mg_vcyc       = maxvcyc
     .              ,sm_it         = nsweep
     .              ,sm_omega      = omeg_smth
     .              ,sm_ncolors    = ncolors
     .              ,sm_line_relax = line_relax
     .              ,sm_zebra_relax= zebra_relax
     .              )

        !Post-process solution with M_star

        call postProcessSchur(3,ntotp,x(:,IVX:IVZ),igrid
     .                       ,bcs(:,IVX:IVZ))

        !Transform coordinates
        if (si_car) call XtoCurv(igrid,x(:,IVX:IVZ))

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

c     Begin program

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
             ii  = i + nx*(j-1) + nx*ny*(k-1)
             xxx(ii,IVX)=rho(i,j,k)*xxx(ii,IVX)
     .                  +xxx(ii,IRHO)*vcnv(i,j,k,1)
cc     .                  +xxx(ii,IRHO)*xxx(ii,IVX)
             xxx(ii,IVY)=rho(i,j,k)*xxx(ii,IVY)
     .                  +xxx(ii,IRHO)*vcnv(i,j,k,2)
cc     .                  +xxx(ii,IRHO)*xxx(ii,IVY)
             xxx(ii,IVZ)=rho(i,j,k)*xxx(ii,IVZ)
     .                  +xxx(ii,IRHO)*vcnv(i,j,k,3)
cc     .                  +xxx(ii,IRHO)*xxx(ii,IVZ)
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
        real(8)    :: dum(ntot/neqd,neqd),dvol

c     Begin program

c     Predictor step

        call invert_M(yyy,xxx)

c     SI step: Deltav --> xxx(:,IVX:IVZ)

        call invert_P_si(yyy,xxx)

c     Corrector step

cc        if (di > 0d0 .and. (.not.hall_new)) then
cc          call find_dx(xxx,rr,4)
cc        else
          call find_dx(xxx,rr,0)
cc        endif

c     End program

      end subroutine split_PC

c     si_PC
c     ##############################################################
      subroutine si_PC(version)

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

c     Local variables

c     Begin program

c     Form (D_m^-1)rx vector ---> xxx

        if (version == 2.or.version == 4) then
          call invert_M(yyy,xxx)
        else
          call invert_D_M(yyy,xxx)
        endif

c     SI step: Deltav --> xxx(:,IVX:IVZ)

        call invert_P_si(yyy,xxx)

c     Solve for dx

        call find_dx(xxx,rr,version)

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

        integer    :: ii,iii,i,j,k
        real(8)    :: dum(ntot/neqd,neqd),dvol

c     Begin program

c     Predictor step

        call invert_M(yyy,xxx)

c     SI step: Deltav --> xxx(:,IVX:IVZ)

        call invert_P_star_si(yyy,xxx)

c     Solve for dx

        call find_dx(xxx,rr,4)

c     End program

      end subroutine star_PC

c     norm_test
c     ##############################################################
      subroutine norm_test

c     --------------------------------------------------------------
c     Calculates SI iteration norm
c     --------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        real(8)      :: mag,ratio
        real(8),save :: mag0,mag_n

c     Begin program

        !Global residual
        mag = sqrt(dot(ntot,rr,rr))

        if (si_it == 1) then
          mag0  = mag
          mag_n = -1d0
        endif

        if (my_rank == 0 .and. iout >= 1) write (*,*)

        ratio = max(mag/mag_n,0d0)

        if (my_rank == 0)
     .    write (*,100) '   SI iter:',si_it-1,';  Res: ',mag
     .                 ,';  Rel.res.:',mag/mag0,';  Ratio:',ratio

        mag_n = mag

        !Individual residual eqs
        if (iout >=2) then
          call scatter(yyy,rr)
          do ieq=1,neqd
            mag = sqrt(dot(ntotdp,yyy(:,ieq),yyy(:,ieq)))
cc            mag = sqrt(sum(yyy(:,ieq)**2))
            if (my_rank == 0) write (*,*) '     Res. eqn',ieq,':',mag
          enddo
          if (my_rank == 0) write (*,*)
        endif

 100    format (a,i3,a,1p,e10.2,a,e10.2,a,e10.2)

      end subroutine norm_test

c     PCdiagnostics
c     ##############################################################
      subroutine PCdiagnostics

c     --------------------------------------------------------------
c     Performs diagnostics on PC
c     --------------------------------------------------------------

        use grid_debug,ONLY:C_contour,cont_conf,dbg,createDrawInCfile

        use local_BCS_variables, ONLY: pinch_flow

        implicit none

        character(3) :: xlabel,ylabel

        integer :: iplot,imin,imax,jmin,jmax,kmin,kmax

        real(8) :: mag0,mag,gxmin,gxmax,gymin,gymax

        real(8),allocatable,dimension(:,:,:,:) :: db_cnv,db_cov,dv0cnv

        character(10) :: solf,resf,rhsf,crsf,fldf,divf
        character(4)  :: dbg_desc(neqd)

#if defined(vec_pot)
        integer    :: ivar

        real(8),allocatable,dimension(:,:,:,:) :: da_cov

        allocate(da_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,db_cov(0:nx+1,0:ny+1,0:nz+1,3))
#endif

c     Begin program

c     Read graphics configuration

        tst_flg = .true.

cc        call readGraphicsInput

        !WARNING: This is using default config. Need to use one in pixie3d.in
        xlabel=cont_conf%label(1)
        ylabel=cont_conf%label(2)

c     Create output files

        solf='sol_p'//trim(int2char(my_rank))//'.bin'
        resf='res_p'//trim(int2char(my_rank))//'.bin'
        rhsf='rhs_p'//trim(int2char(my_rank))//'.bin'
        crsf='crs_p'//trim(int2char(my_rank))//'.bin'
        fldf='fld_p'//trim(int2char(my_rank))//'.bin'
        divf='div_p'//trim(int2char(my_rank))//'.bin'

        do iplot=1,neqd
          dbg_desc(iplot) = 'dgb'//trim(int2char(iplot))
        enddo

        call createDrawInCfile(neqd,trim(solf),'Solution'
     .          ,'t',trim(xlabel),trim(ylabel),dbg_desc,'-c -X0 -L57'
     .          ,'drawsol_p'//trim(int2char(my_rank))//'.in')

        call createDrawInCfile(neqd,trim(resf),'Residual'
     .          ,'t',trim(xlabel),trim(ylabel),dbg_desc,'-c -X0 -L57'
     .          ,'drawres_p'//trim(int2char(my_rank))//'.in')

        call createDrawInCfile(neqd,trim(rhsf),'RHS     '
     .          ,'t',trim(xlabel),trim(ylabel),dbg_desc,'-c -X0 -L57'
     .          ,'drawrhs_p'//trim(int2char(my_rank))//'.in')

        call createDrawInCfile(neqd,trim(crsf),'Corr. residual'
     .          ,'t',trim(xlabel),trim(ylabel),dbg_desc,'-c -X0 -L57'
     .          ,'drawcrhs_p'//trim(int2char(my_rank))//'.in')

        call createDrawInCfile(6,trim(fldf),'Vector fields'
     .          ,'t',trim(xlabel),trim(ylabel)
     .          ,(/'V1','V2','V3','B1','B2','B3'/)
     .          ,'-c -X0 -L57'
     .          ,'drawvfld_p'//trim(int2char(my_rank))//'.in')

        call createDrawInCfile(1,trim(divf),'div(B)'
     .          ,'t',trim(xlabel),trim(ylabel),(/'div(B)'/)
     .          ,'-c -X0 -L57'
     .          ,'drawdiv_p'//trim(int2char(my_rank))//'.in')

c     Dump graphics

        allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dv0cnv(0:nx+1,0:ny+1,0:nz+1,3))

        !Solution plot
        call MGplot(neqd,x,igrid,0,solf)

        !Predictor RHS plot
        call MGplot(neqd,y,igrid,0,rhsf)

        !Residual plot
#if !defined(petsc_c) && !defined(samrai)
        call matrixFreeMatVec(ntot,x,rr)
#else
        rr = 0d0
#endif
        mag0 = sqrt(dot(ntot,y,y))
        y = y - rr
        mag  = sqrt(dot(ntot,y,y))
        call MGplot(neqd,y,igrid,0,resf)

        if (my_rank == 0) then
           write (*,*)
           write (*,*) 'LINEAR RESIDUAL'
           write (*,*) 'abs res=',mag
           write (*,*) 'abs res (per node)=',mag/ntot
     $                ,', rel res=',mag/mag0
        endif

        !Individual residual eqs
        call scatter(yyy,y)
        do ieq=1,neqd
          mag = sqrt(dot(ntotdp,yyy(:,ieq),yyy(:,ieq)))
          if (my_rank == 0) write (*,*) 'Res. eqn',ieq,':',mag
        enddo

        !Corrector rhs plot
        if (si_car) call XtoCurv(igrid,rhspc(:,IVX:IVZ))

        call gather(rr,rhspc)
        call MGplot(neqd,rr,igrid,0,crsf)

        !Vector components plot
#if !defined(vec_pot)
        !Find dB* (w/ BCs)
        if (pinch_flow) then
          call fillArray(igrid,IBX,3,xxx(:,IBX:IBZ),db_cnv
     .                ,bcs(:,IBX:IBZ),2,is_vec=.true.,is_cnv=.true.
     .                ,arr0=gb0%grid(igrid)%array)
        else
          call fillArray(igrid,IBX,3,xxx(:,IBX:IBZ),db_cnv
     .                  ,bcs(:,IBX:IBZ),2,is_vec=.true.,is_cnv=.true.)
        endif
#else
        call fillArray(igrid,IAX,3,xxx(:,IAX:IAZ),da_cov,bcs(:,IAX:IAZ)
     .                ,2,is_vec=.true.,is_cnv=.false.)

        db_cnv = 0d0
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              do ivar=1,3
                db_cnv(i,j,k,ivar)=curl2(i,j,k,nx,ny,nz
     .                                ,igrid,igrid,igrid
     .                                ,da_cov(:,:,:,1)
     .                                ,da_cov(:,:,:,2)
     .                                ,da_cov(:,:,:,3)
     .                                ,ivar)
              enddo
            enddo
          enddo
        enddo

        !Find covariant components of dj with BCs
        call setBC(IBX,3,nx,ny,nz,db_cnv,db_cov,vzeros,bcs(:,IBX:IBZ)
     .          ,igrid,igrid,igrid,iorder=2)
#endif

        dv0cnv = 0d0

        call ExB_flow_bc(igrid,nx,ny,nz,bcnv,bcov,dv0cnv,dbcov=db_cov)

        !Find dv (w/ BCs)
        call fillArray(igrid,IVX,3,xxx(:,IVX:IVZ),dv_cnv,bcs(:,IVX:IVZ)
     .                ,2,is_vec=.true.,is_cnv=.true.,arr0=dv0cnv)

        open(unit=110,file=fldf,form='unformatted'
     .     ,status='replace')
        call C_contour(dv_cnv(1:nx,1:ny,1:nz,1)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,0,110)
        call C_contour(dv_cnv(1:nx,1:ny,1:nz,2)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,1,110)
        call C_contour(dv_cnv(1:nx,1:ny,1:nz,3)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,1,110)
cc        call contour(dv_cnv(imin:imax,jmin:jmax,kmin:kmax,1),nx,ny,gxmin
cc     .              ,gxmax,gymin,gymax,0,110)
cc        call contour(dv_cnv(imin:imax,jmin:jmax,kmin:kmax,2),nx,ny,gxmin
cc     .              ,gxmax,gymin,gymax,1,110)
cc        call contour(dv_cnv(imin:imax,jmin:jmax,kmin:kmax,3),nx,ny,gxmin
cc     .              ,gxmax,gymin,gymax,1,110)
#if !defined(vec_pot)
cc        call contour(db_cnv(imin:imax,jmin:jmax,kmin:kmax,1),nx,ny,gxmin
cc     .              ,gxmax,gymin,gymax,1,110)
cc        call contour(db_cnv(imin:imax,jmin:jmax,kmin:kmax,2),nx,ny,gxmin
cc     .              ,gxmax,gymin,gymax,1,110)
cc        call contour(db_cnv(imin:imax,jmin:jmax,kmin:kmax,3),nx,ny,gxmin
cc     .              ,gxmax,gymin,gymax,1,110)
        call C_contour(db_cnv(1:nx,1:ny,1:nz,1)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,1,110)
        call C_contour(db_cnv(1:nx,1:ny,1:nz,2)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,1,110)
        call C_contour(db_cnv(1:nx,1:ny,1:nz,3)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,1,110)
#else
cc        call contour(da_cov(imin:imax,jmin:jmax,kmin:kmax,1),nx,ny,gxmin
cc     .              ,gxmax,gymin,gymax,1,110)
cc        call contour(da_cov(imin:imax,jmin:jmax,kmin:kmax,2),nx,ny,gxmin
cc     .              ,gxmax,gymin,gymax,1,110)
cc        call contour(da_cov(imin:imax,jmin:jmax,kmin:kmax,3),nx,ny,gxmin
cc     .              ,gxmax,gymin,gymax,1,110)
        call C_contour(da_cov(1:nx,1:ny,1:nz,1)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,1,110)
        call C_contour(da_cov(1:nx,1:ny,1:nz,2)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,1,110)
        call C_contour(da_cov(1:nx,1:ny,1:nz,3)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,1,110)
#endif
        close(110)

        !diag B-field divergence
        allocate(dbg(nx,ny,nz,1))
        dbg = 0d0
        do k=1,nz
          do j=1,ny
            do i=1,nx
              dbg(i,j,k,1) = div(i,j,k,nx,ny,nz,igrid,igrid,igrid
     $                          ,db_cnv(:,:,:,1)
     $                          ,db_cnv(:,:,:,2)
     $                          ,db_cnv(:,:,:,3))
            enddo
          enddo
        enddo

        mag = sqrt(dot(ntotdp,dbg,dbg))
        if (my_rank == 0) write (*,*) 'L-2 norm of div(dB)',mag

        open(unit=110,file=divf,form='unformatted',status='replace')
cc        call contour(dbg(imin:imax,jmin:jmax,kmin:kmax,1),nx,ny
cc     .              ,gxmin,gxmax,gymin,gymax,0,110)
        call C_contour(dbg(1:nx,1:ny,1:nz,1)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,0,110)
        close(110)

        deallocate(dbg)

        call pstop('PCdiagnostics','Finished PC diagnostics')

      end subroutine PCdiagnostics

      end subroutine applyPreconditioner

c findEOMrhs
c #########################################################################
      subroutine findEOMrhs(ntotp,xxx,eom_rhs,igrid,part,vwgt)
c--------------------------------------------------------------------
c     This subroutine finds the EM part of the rhs for the velocity
c     SI solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use imposeBCinterface

      use local_BCS_variables, ONLY: ExB_flow_bc,pinch_flow

      implicit none

c Call variables

      integer    :: ntotp,igrid
      real(8)    :: xxx(ntotp,neqd),eom_rhs(ntotp,3)
      character(3) :: part
      logical    :: vwgt

c Local variables

      integer    :: ii,isig,iig,ivar,ip,im,nx,ny,nz
      real(8)    :: cov(3),cnv(3),idx,idy,idz,dvol,nabla_v(3,3)
     .             ,nabla_vv0(3,3),mul,veclap(3),advc(3),jac
      real(8),allocatable,dimension(:,:,:,:),target :: dj_cnv
     .                                                ,dj_cov
     .                                                ,db_cnv
     .                                                ,db_cov
     .                                                ,dv_cnv
     .                                                ,dv0cnv
     .                                                ,dp_cnv
     .                                                ,dve_cnv
     .                                                ,drho  
     .                                                ,dtmp  
#if defined(vec_pot)
     .                                                ,da_cov                           
#endif

c Begin program

      isig = grid_params%istartp(igrid)

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      if (      part /= 'lfz'
     .    .and. part /= 'adv'
     .    .and. part /= 'all') then
        call pstop ('findEOMrhs','Unknown character specifier '//part)
      endif

c Allocate variables

      allocate(dj_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,db_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dv0cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dp_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dve_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,drho  (0:nx+1,0:ny+1,0:nz+1,1)
     .        ,dtmp  (0:nx+1,0:ny+1,0:nz+1,1)
#if defined(vec_pot)
     .        ,da_cov(0:nx+1,0:ny+1,0:nz+1,3)
#endif
     .        )

#if !defined(vec_pot)
c Find dB* (w/ BCs)

      if (pinch_flow) then
        call fillArray(igrid,IBX,3,xxx(:,IBX:IBZ),db_cnv,bcs(:,IBX:IBZ)
     .              ,2,arr_cov=db_cov,is_vec=.true.,is_cnv=.true.
     .              ,arr0=gb0%grid(igrid)%array)
      else
        call fillArray(igrid,IBX,3,xxx(:,IBX:IBZ),db_cnv,bcs(:,IBX:IBZ)
     .              ,2,arr_cov=db_cov,is_vec=.true.,is_cnv=.true.)
      endif

#else
c Find dA* (w/ BCs)

      call fillArray(igrid,IAX,3,xxx(:,IAX:IAZ),da_cov,bcs(:,IAX:IAZ)
     .              ,2,is_vec=.true.,is_cnv=.false.)

c Find dB* from dA*

      db_cnv = 0d0
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            do ivar=1,3
              db_cnv(i,j,k,ivar)=curl2(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                                ,da_cov(:,:,:,1)
     .                                ,da_cov(:,:,:,2)
     .                                ,da_cov(:,:,:,3)
     .                                ,ivar)
            enddo
          enddo
        enddo
      enddo

      !Find covariant components of dj with BCs
      call setBC(IBX,3,nx,ny,nz,db_cnv,db_cov,vzeros,bcs(:,IBX:IBZ)
     .          ,igrid,igrid,igrid,iorder=2)
#endif

c Store db_cov for ExB flow BCs

      if (pinch_flow) then
        gdb_cov%grid(igrid)%array = db_cov

        call restrictMGArray(IBX,3,gdb_cov,bcs(:,IBX:IBZ),igrid
     .                      ,2,iscnv=.false.,equilb=gb0)
      endif

c Find dv (w/ Bcs)

      dv0cnv = 0d0

      call ExB_flow_bc(igrid,nx,ny,nz,bcnv,bcov,dv0cnv,dbcov=db_cov)

      call fillArray(igrid,IVX,3,xxx(:,IVX:IVZ),dv_cnv,bcs(:,IVX:IVZ)
     .              ,2,is_vec=.true.,is_cnv=.true.,arr0=dv0cnv)

c Find dj* from dB*

      !Find contravariant current (without BCs)
      dj_cnv = 0d0
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            do ivar=1,3
              dj_cnv(i,j,k,ivar)=curl2(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                                ,db_cov(:,:,:,1)
     .                                ,db_cov(:,:,:,2)
     .                                ,db_cov(:,:,:,3)
     .                                ,ivar)
            enddo
cc            dj_cnv(i,j,k,:)=curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                              ,db_cov(:,:,:,1)
cc     .                              ,db_cov(:,:,:,2)
cc     .                              ,db_cov(:,:,:,3))
          enddo
        enddo
      enddo

      !Find covariant components of dj with BCs
      call setBC(IJX,3,nx,ny,nz,dj_cnv,dj_cov,vzeros,bcs(:,IJX:IJZ)
     .          ,igrid,igrid,igrid,iorder=2)

c Find drho, dtmp arrays with BCs

      call fillArray(igrid,IRHO,1,xxx(:,IRHO:IRHO),drho
     .              ,bcs(:,IRHO:IRHO),2)
      call fillArray(igrid,ITMP,1,xxx(:,ITMP:ITMP),dtmp
     .              ,bcs(:,ITMP:ITMP),2)

c Find dp_cnv

      if (.not.nc_eom_v .and.(part == 'adv'.or.part == 'all')) then
        dp_cnv(:,:,:,1)=rho*dv_cnv(:,:,:,1)+drho(:,:,:,1)*vcnv(:,:,:,1)
        dp_cnv(:,:,:,2)=rho*dv_cnv(:,:,:,2)+drho(:,:,:,1)*vcnv(:,:,:,2)
        dp_cnv(:,:,:,3)=rho*dv_cnv(:,:,:,3)+drho(:,:,:,1)*vcnv(:,:,:,3)
      endif

c Find electron velocity

      if (di*heta > 0d0) then
        dve_cnv(:,:,:,1)=dv_cnv(:,:,:,1)-di*dj_cnv(:,:,:,1)/rho
        dve_cnv(:,:,:,2)=dv_cnv(:,:,:,2)-di*dj_cnv(:,:,:,2)/rho
        dve_cnv(:,:,:,3)=dv_cnv(:,:,:,3)-di*dj_cnv(:,:,:,3)/rho
      endif

c Find rhs_v'

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx

            ii  = i + nx*(j-1) + nx*ny*(k-1)
            iig = ii + isig - 1

            jac = gmetric%grid(igrid)%jac(i,j,k)

            if (vwgt) then
              dvol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            if (part == 'lfz' .or. part == 'all') then

              if (nc_eom_jxb) then
                if (.not.adiabatic) then
                  cnv = eom_f(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                     ,drho(:,:,:,1),tmp,db_cov,jcov,.false.
     .                     ,nc_eom_jxb,nc_eom_gp,.false.
     .                     ,vol=.false.)
                else
                  cnv = eom_jxb(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                     ,db_cov,jcov,.false.,nc_eom_jxb,.false.
     .                     ,vol=.false.)
                endif

                cnv = cnv
     .               +eom_f(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                     ,rho,dtmp(:,:,:,1),bcov,dj_cov,.false.
     .                     ,nc_eom_jxb,nc_eom_gp,.false.
     .                     ,vol=.false.)
              else
                if (.not.adiabatic) then
                  cnv = eom_f(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                   ,drho(:,:,:,1),tmp,db_cnv,bcnv,.true.
     .                   ,nc_eom_jxb,nc_eom_gp,.false.
     .                   ,vol=.false.)
                else
                  cnv = eom_jxb(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                   ,db_cnv,bcnv,.true.,nc_eom_jxb,.false.
     .                   ,vol=.false.)
                endif

                cnv = cnv
     .               +eom_f(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                   ,rho,dtmp(:,:,:,1),bcnv,db_cnv,.true.
     .                   ,nc_eom_jxb,nc_eom_gp,.false.
     .                   ,vol=.false.)
              endif

              !Substract div_pe piece corresponding to (-di*j/rho)
              if (di*heta > 0d0 .and. part == 'lfz') then
                cnv = cnv
     .              - heta*EOM_divPi(i,j,k,nx,ny,nz
     .                              ,igrid,igrid,igrid,dve_cnv)
              endif

              if (nc_eom_v) then
                cnv = -cnv/rho(i,j,k)
              else
                cnv = -cnv
              endif
            else
              cnv = 0d0
            endif

            if (part == 'adv' .or. part == 'all') then
              
              if (nc_eom_v) then
                cnv = cnv
     .              +(EOM_advc(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,vcnv,dv_cnv,upwind=(advect == 1))
     .               +EOM_advc(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,dv_cnv,vcnv,upwind=(advect == 1)))
              else
                cnv = cnv
     .              +(EOM_advc(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,vcnv,dp_cnv,upwind=(advect == 1))
     .               +EOM_advc(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,dp_cnv,vcnv,upwind=(advect == 1)))
              endif

              !Ion stress tensor
              cnv = cnv
     .            + EOM_divPi(i,j,k,nx,ny,nz,igrid,igrid,igrid,dv_cnv
     .                       ,diff=gvis%grid(igrid)%array(:,:,:,1))

              !Electron stress tensor
              if (heta > 0d0) then
                cnv = cnv
     .            + heta*EOM_divPi(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                            ,dve_cnv)
              endif

              if (nc_eom_v) cnv = cnv/rho(i,j,k)     
            endif

            !Correct rhs_v
            eom_rhs(ii,:) = alpha*dvol*cnv(:)

          enddo
        enddo
      enddo

c Free memory

      deallocate(dj_cnv
     .          ,dj_cov
     .          ,db_cnv
     .          ,db_cov
     .          ,dv_cnv
     .          ,dv0cnv
     .          ,dp_cnv
     .          ,dve_cnv
     .          ,drho  
     .          ,dtmp  
#if defined(vec_pot)
     .          ,da_cov
#endif
     .          )

c End program

      end subroutine findEOMrhs

c formSIrhs
c #########################################################################
      subroutine formSIrhs(ntotp,xxx,yyy,rhs_si,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the rhs for the velocity SI solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: ntotp,igrid,nx,ny,nz
      real(8)    :: xxx(ntotp,neqd),yyy(ntotp,3),rhs_si(ntotp,3)

c Local variables

      integer    :: ii,iig,isig,ivar,ip,im
      real(8)    :: dvol,cov(3),cnv(3),idx,idy,idz,eom_rhs(ntotp,3)

      real(8),allocatable,dimension(:,:,:,:),target :: dp_cnv,drho

c Begin program

      isig = grid_params%istartp(igrid)

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      allocate(dp_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,drho  (0:nx+1,0:ny+1,0:nz+1,1))

c Find rhs_v

      if (.not.nc_eom_v) then

        call fillArray(igrid,IRHO,1,xxx(:,IRHO:IRHO),drho
     .                ,bcs(:,IRHO:IRHO),2)

        dp_cnv(:,:,:,1)=drho(:,:,:,1)*gv0%grid(igrid)%array(:,:,:,1)
        dp_cnv(:,:,:,2)=drho(:,:,:,1)*gv0%grid(igrid)%array(:,:,:,2)
        dp_cnv(:,:,:,3)=drho(:,:,:,1)*gv0%grid(igrid)%array(:,:,:,3)

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              iig = ii + isig - 1
              if (vol_wgt) then
                dvol = gmetric%grid(igrid)%dvol(i,j,k)
              else
                dvol = 1d0
              endif

cc              rhs_si(ii,:) = yyy(ii,:)
cc     .                     - dvol*xxx(ii,IRHO)*mgadvdiffV0(iig,:)

              vec1 => gv0%grid(igrid)%array
              vec2 => dp_cnv
              cnv = div_tensor(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,alt_eom,eom_advc_x,eom_advc_y,eom_advc_z
     .                        ,vol=vol_wgt)
              nullify(vec1,vec2)

              rhs_si(ii,:) = yyy(ii,:)
     .                     - dvol*(cnp*dp_cnv(i,j,k,:)/dt) - alpha*cnv

            enddo
          enddo
        enddo

      else

cc        rhs_si = yyy

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              iig = ii + isig - 1
              if (vol_wgt) then
                dvol = gmetric%grid(igrid)%dvol(i,j,k)
              else
                dvol = 1d0
              endif
              rhs_si(ii,:) = yyy(ii,:)
     .                     - dvol*xxx(ii,IRHO)*mgadvdiffV0(iig,:)
            enddo
          enddo
        enddo

      endif

c Find linearized EOM rhs [jxB-grad(p)] correction

      call findEOMrhs(ntotp,xxx,eom_rhs,igrid,'lfz',vol_wgt)

c Add rhs_v + EOM rhs

      rhs_si = rhs_si - eom_rhs

c Transform to cartesian if applicable

      if (si_car) call XtoCar(igrid,rhs_si)

c Free memory

      deallocate(dp_cnv,drho)

c End program

      end subroutine formSIrhs

c formArhs_hall
c #########################################################################
      subroutine formArhs_hall(ntotp,xxx,yyy,rhs_a,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the rhs for the velocity SI solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: ntotp,igrid
      real(8)    :: xxx(ntotp,neqd),yyy(ntotp,3),rhs_a(ntotp,3)

c Local variables

      integer    :: ii,iig,isig,ivar,ip,im,nx,ny,nz
      real(8)    :: dvol,cov(3),cnv(3),idx,idy,idz,eom_rhs(ntotp,3)

      real(8),allocatable,dimension(:,:,:,:),target :: dpe,drho,vdum
      real(8),pointer    ,dimension(:,:,:,:) :: j0,rho0,tmp0,Eh0

c Begin program

      isig = grid_params%istartp(igrid)

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

c Auxiliary variables

      allocate(dpe (0:nx+1,0:ny+1,0:nz+1,1)
     .        ,drho(0:nx+1,0:ny+1,0:nz+1,1)
     .        ,vdum(0:nx+1,0:ny+1,0:nz+1,3 ))

      j0   => gj0  %grid(igrid)%array
      rho0 => grho0%grid(igrid)%array
      tmp0 => gtmp0%grid(igrid)%array
      Eh0  => gEh0_cov%grid(igrid)%array

      call fillArray(igrid,IRHO,1,xxx(:,IRHO:IRHO),drho
     .              ,bcs(:,IRHO:IRHO),2)

      call fillArray(igrid,ITMP,1,xxx(:,ITMP:ITMP),dpe
     .              ,bcs(:,ITMP:ITMP),2)

      dpe  = rho0*dpe + tmp0*drho

      vdum(:,:,:,1) = di*j0(:,:,:,1)/rho0(:,:,:,1)**2*drho(:,:,:,1)
      vdum(:,:,:,2) = di*j0(:,:,:,2)/rho0(:,:,:,1)**2*drho(:,:,:,1)
      vdum(:,:,:,3) = di*j0(:,:,:,3)/rho0(:,:,:,1)**2*drho(:,:,:,1)

c Find rhs_a

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            iig = ii + isig - 1
            if (vol_wgt) then
              dvol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            call grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,dpe
     .               ,cov(1),cov(2),cov(3))

            cnv=heta*EOM_divPi(i,j,k,nx,ny,nz,igrid,igrid,igrid,vdum)

            cov = cov + matmul(gmetric%grid(igrid)%gsub(i,j,k,:,:),cnv)

            rhs_a(ii,:) = yyy(ii,:) - alpha*dvol
     .                  *(-drho(i,j,k,1)/rho0(i,j,k,1)*Eh0(i,j,k,:)
     .                    -di/rho(i,j,k)*cov)
          enddo
        enddo
      enddo

c Transform to cartesian if applicable

      if (si_car) call XtoCar(igrid,rhs_a)

c Free memory

      deallocate(dpe,drho,vdum)

c End program

      end subroutine formArhs_hall

c correctVfldRHS
c #########################################################################
      subroutine correctVfldRHS(ntotp,x,crhs,igrid,iout)
c--------------------------------------------------------------------
c     This subroutine finds the rhs for the velocity SI solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: ntotp,igrid,iout
      real(8)    :: x(ntotp,neqd),crhs(ntotp,3)

c Local variables

      integer    :: ii,iig,isig,ivar,ip,im,guess,nx,ny,nz
      real(8)    :: nrhs(ntotp,3),dvol,cov(3),cnv(3),idx,idy,idz
     .             ,nabla_vv0(3,3),mul,veclap(3),advc(3)

      external   :: v_mtvc2

c Begin program

      isig = grid_params%istartp(igrid)

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

c Find SI rhs (but with "complete" quantities)

      call formSIrhs(ntotp,x,crhs,nrhs,igrid)

c Solve for advection-diffusion part

      guess = 1

      call cSolver(3,ntotp,nrhs,x(:,IVX:IVZ)
     .            ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc2,vol_wgt
     .            ,tol           = mgtol
     .            ,dg            = v_diag2
     .            ,mg_debug      = mg_debug
     .            ,mg_galerkin   = mg_glrkin
     .            ,mg_order_res  = ores
     .            ,mg_order_prol = oprol
     .            ,mg_coarse_grid_size = mg_coarse_size
     .            ,mg_gm_coarse_solve  = gm_coarse
     .            ,mg_smooth     = smooth
     .            ,mg_vcyc       = maxvcyc
     .            ,sm_it         = nsweep
     .            ,sm_omega      = omeg_smth
     .            ,sm_ncolors    = ncolors
     .            ,sm_line_relax = line_relax
     .            ,sm_zebra_relax= zebra_relax
     .            )

      return

c Correction w/o solve

      call findEOMrhs(ntotp,x,nrhs,igrid,'all',vol_wgt)

      nrhs = crhs - nrhs

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            if (vol_wgt) then
              dvol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              dvol = 1d0
            endif
            x(ii,IVX:IVZ) = dt*nrhs(ii,:)/dvol/cnp
          enddo
        enddo
      enddo

c End program

      end subroutine correctVfldRHS

c correctRhoRHS
c #########################################################################
      subroutine correctRhoRHS(nx,ny,nz,dv_cnv,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the density solve,
c     div(rh0*dv)
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: igrid,nx,ny,nz
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz)

c Local variables

      integer    :: ii

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            crhs(ii) = flx_advec(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                      ,dv_cnv(:,:,:,1)
     .                      ,dv_cnv(:,:,:,2)
     .                      ,dv_cnv(:,:,:,3),rho,advect,vol=.false.)
          enddo
        enddo
      enddo

c End program

      end subroutine correctRhoRHS

c correctTmpRHS
c #########################################################################
      subroutine correctTmpRHS(nx,ny,nz,dv_cnv,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the temperature
c     solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: igrid,nx,ny,nz
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz)

c Local variables

      integer    :: ii
      real(8)    :: jac

c Begin program

c Evaluate rhs correction: dv*grad(T0) + (gamma-1)*T0*div(dv)

      if (adiabatic) then
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              jac = gmetric%grid(igrid)%jac(i,j,k)

              crhs(ii) = flx_advec(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,dv_cnv(:,:,:,1)
     .                        ,dv_cnv(:,:,:,2)
     .                        ,dv_cnv(:,:,:,3),tmp,advect,vol=.false.)
     .                  +(gamma-1)*tmp(i,j,k)
     $                            *div(i,j,k,nx,ny,nz,igrid,igrid,igrid
     $                                ,dv_cnv(:,:,:,1)
     $                                ,dv_cnv(:,:,:,2)
     $                                ,dv_cnv(:,:,:,3))
            enddo
          enddo
        enddo
      else
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              jac = gmetric%grid(igrid)%jac(i,j,k)

              crhs(ii) = flx_advec(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,dv_cnv(:,:,:,1)
     .                        ,dv_cnv(:,:,:,2)
     .                        ,dv_cnv(:,:,:,3),tmp,advect,vol=.false.)
     .                  +(gamma-2)*tmp(i,j,k)
     $                            *div(i,j,k,nx,ny,nz,igrid,igrid,igrid
     $                                ,dv_cnv(:,:,:,1)
     $                                ,dv_cnv(:,:,:,2)
     $                                ,dv_cnv(:,:,:,3))
            enddo
          enddo
        enddo
      endif

c End program

      end subroutine correctTmpRHS

c correctBflRHS
c ###################################################################
      subroutine correctBflRHS(nx,ny,nz,dv_cnv,dv_cov,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the magnetic field
c     solve: div(dv B0 - B0 dv)
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: igrid,nx,ny,nz
      real(8),target :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
      real(8)    :: dv_cov(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz,3)

c Local variables

      integer    :: ii

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            crhs(ii,:) = vecA(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                       ,dv_cnv,dv_cov
     .                       ,ga0%grid(igrid)%array,0)
          enddo
        enddo
      enddo

      return

cc      vec1 => dv_cnv
cc      vec2 => gb0%grid(igrid)%array
cc
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            ii  = i + nx*(j-1) + nx*ny*(k-1)
cccc            crhs(ii,:) = curl_bxv(i,j,k,nx,ny,nz,igrid,igrid,igrid
cccc     .                       ,dv_cnv
cccc     .                       ,gb0%grid(igrid)%array,0)
cc            crhs(ii,:) = div_tensor(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                       ,.false.,btensor_x,btensor_y,btensor_z
cc     .                       ,vol=.false.)
cc          enddo
cc        enddo
cc      enddo
cc
cc      nullify(vec1,vec2)
cc
cc      return

ccc Find Dv.dv --> dummy
cc
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            dummy(i,j,k,:) = diag_v(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                             ,gv0%grid(igrid)%array,dv_cnv)
cc          enddo
cc        enddo
cc      enddo
cc
ccc Find Hall contribution
cc
cc      call setMGBC(0,3,nx,ny,nz,igrid,dummy,bcs(:,IVX:IVZ)
cc     .            ,arr_cov=dum_cov,icomp=IVX,is_vec=.true.
cc     .            ,is_cnv=.true.,iorder=3)
cc
cccc      dum_cov = cnp/dt/alpha*dv_cov
cc
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            ii  = i + nx*(j-1) + nx*ny*(k-1)
cc            crhs(ii,:) = crhs(ii,:)
cc     .                 + di/alpha*curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                                ,dum_cov(:,:,:,1)
cc     .                                ,dum_cov(:,:,:,2)
cc     .                                ,dum_cov(:,:,:,3))
cc          enddo
cc        enddo
cc      enddo

c End program

      end subroutine correctBflRHS

c correctAflRHS
c ###################################################################
      subroutine correctAflRHS(nx,ny,nz,dv_cnv,dv_cov,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the vector potential
c     solve: div(dv B0 - B0 dv)
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: igrid,nx,ny,nz
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,crhs(nx*ny*nz,3)

c Local variables

      integer    :: ii
      real(8)    :: minv,cnv(3),cov(3)
     .             ,dp_cnv(0:nx+1,0:ny+1,0:nz+1,3)

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            crhs(ii,1) = dv_cnv(i,j,k,2)*gbs0%grid(igrid)%array(i,j,k,3)
     .                 - dv_cnv(i,j,k,3)*gbs0%grid(igrid)%array(i,j,k,2)

            crhs(ii,2) = dv_cnv(i,j,k,3)*gbs0%grid(igrid)%array(i,j,k,1)
     .                 - dv_cnv(i,j,k,1)*gbs0%grid(igrid)%array(i,j,k,3)

            crhs(ii,3) = dv_cnv(i,j,k,1)*gbs0%grid(igrid)%array(i,j,k,2)
     .                 - dv_cnv(i,j,k,2)*gbs0%grid(igrid)%array(i,j,k,1)

            crhs(ii,:) = -crhs(ii,:)/gmetric%grid(igrid)%jac(i,j,k)

            !Hall MHD part
            if (di > 0d0) then
              if (hall_new) then
                cnv = findPsiv(i,j,k,nx,ny,nz,igrid
     .                        ,grho0%grid(igrid)%array
     .                        ,gv0  %grid(igrid)%array,dv_cnv)

                call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                      ,cov(1),cov(2),cov(3)
     .                                      ,cnv(1),cnv(2),cnv(3)
     .                                      ,.false.)

                crhs(ii,:) = crhs(ii,:) + di/alpha*cov

cc                minv = di*cnp/dt/alpha
cc                crhs(ii,:) = crhs(ii,:) + minv*dv_cov(i,j,k,:)
              else
                if (heta > 0d0) then
                  cnv = di*heta/grho0%grid(igrid)%array(i,j,k,1)
     .               *EOM_divPi(i,j,k,nx,ny,nz,igrid,igrid,igrid,dv_cnv)

                  call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                      ,cov(1),cov(2),cov(3)
     .                                      ,cnv(1),cnv(2),cnv(3)
     .                                      ,.false.)

                  crhs(ii,:) = crhs(ii,:) + cov
                endif
              endif
            endif

          enddo
        enddo
      enddo

c End program

      end subroutine correctAflRHS

c preprocessBflRHS
c ###################################################################
      subroutine preprocessBflRHS(ntotp,neq,yyy,igrid)
c--------------------------------------------------------------------
c     This subroutine preprocesses B rhs to add Hall correction.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: igrid,ntotp,neq

      real(8)    :: yyy(ntotp,neq)

c Local variables

      integer :: ii,nx,ny,nz

      real(8) :: cnv(3),dvol

      real(8),allocatable,dimension(:,:,:,:) :: dgv_cnv,dgv_cov

c Begin program

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      allocate(dgv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dgv_cov(0:nx+1,0:ny+1,0:nz+1,3))

c Store velocity RHS in array format w/ BCs

      call fillArray(igrid,IVX,3,yyy(:,IVX:IVZ),dgv_cnv,bcs(:,IVX:IVZ)
     .              ,2,arr_cov=dgv_cov,is_vec=.true.,is_cnv=.true.)

c Divide by density and volume

      if (vol_wgt) then
        do ieq=1,3
          dgv_cov(:,:,:,ieq) = dgv_cov(:,:,:,ieq)
     .                        /gmetric%grid(igrid)%dvol
        enddo
      endif

      if (.not.nc_eom_v) then
        do ieq=1,3
          dgv_cov(:,:,:,ieq) = dgv_cov(:,:,:,ieq)/rho
        enddo
      endif

c Find curl

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            if (vol_wgt) then   !Residuals are volume-weighed by default
              dvol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              dvol = 1d0
            endif
            yyy(ii,IBX:IBZ) = yyy(ii,IBX:IBZ)
     .                   + di*dvol*curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                                 ,dgv_cov(:,:,:,1)
     .                                 ,dgv_cov(:,:,:,2)
     .                                 ,dgv_cov(:,:,:,3))
          enddo
        enddo
      enddo

      deallocate(dgv_cnv,dgv_cov)

c End program

      end subroutine preprocessBflRHS

c preprocessAflRHS
c ###################################################################
      subroutine preprocessAflRHS(ntotp,neq,yyy,igrid)
c--------------------------------------------------------------------
c     This subroutine preprocesses B rhs to add Hall correction.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      implicit none

c Call variables

      integer    :: igrid,ntotp,neq

      real(8)    :: yyy(ntotp,neq)

c Local variables

      integer    :: ii,nx,ny,nz
      real(8)    :: cov(3),cnv(3)

c Begin program

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

            if (nc_eom_v) cov = cov/rho(i,j,k)

            yyy(ii,IAX:IAZ) = yyy(ii,IAX:IAZ) + di*cov
          enddo
        enddo
      enddo

c End program

      end subroutine preprocessAflRHS

c divClean
c ###################################################################
      subroutine divClean(iout,ntotp,xxx,igrid)
c--------------------------------------------------------------------
c     This subroutine divergence-cleans B by integrating a vector
c     potential and the differentiating it. It takes B in cartesian
c     or curvilinear and returns Bcnv.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use imposeBCinterface

      use setMGBC_interface

      implicit none

c Call variables

      integer    :: iout,ntotp,igrid
      real(8)    :: xxx(ntotp,3)

c Local variables

      integer    :: ii,ivar,nx,ny,nz
      real(8)    :: a1,a2,a3,etal
      real(8),allocatable,dimension(:,:,:,:) :: db_cnv,da_cov

c Begin program

      if (iout > 0) write (*,*) 'Performing divergence cleaning...'

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,da_cov(0:nx+1,0:ny+1,0:nz+1,3))

c Find dB* (w/ BCs)

      call fillArray(igrid,IBX,3,xxx,db_cnv,bcs(:,IBX:IBZ)
     .              ,2,is_vec=.true.,is_cnv=.true.)

c Find dA* from dB* (w/ ghost cells)

      call curl_inv(nx,ny,nz,igrid,igrid,igrid,db_cnv,da_cov)

c Find clean dB=curl(dA)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            xxx(ii,:)=curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                    ,da_cov(:,:,:,1)
     .                    ,da_cov(:,:,:,2)
     .                    ,da_cov(:,:,:,3))

          enddo
        enddo
      enddo

      deallocate(db_cnv,da_cov)

c End program

      end subroutine divClean

c postProcessSchur
c ###################################################################
      subroutine postProcessSchur(neq,ntotp,xxx,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine postprocess SI solution to find final
c     solution to Schur complement system.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use imposeBCinterface

      implicit none

c Call variables

      integer    :: ntotp,igrid,neq,bcnd(6,neq)
      real(8)    :: xxx(ntotp,neq)

c Local variables

      integer    :: ii,ntot
      real(8)    :: xi(ntotp*neq),bi(ntotp*neq)

c Begin program

      ntot = ntotp*neq

c Gather variables

      do ii=1,ntotp
        do ieq=1,3
          xi(neq*(ii-1)+ieq) = xxx(ii,ieq)
        enddo
      enddo

c Perform matvec with M_star

      call mstar_mtvc(0,3,ntotp*neq,xi,bi,igrid,bcnd)

c Scatter solution

      do ii=1,ntotp
        do ieq=1,neq
          xxx(ii,ieq) = bi(neq*(ii-1)+ieq)
        enddo
      enddo

c End program

      end subroutine postProcessSchur
