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

      use matvec

      use imposeBCinterface

      implicit none

c Call variables

      integer    :: ntot,iout
      real(8)    :: x(ntot),y(ntot)

c Local variables

      real(8)    :: xxx(ntot/neqd,neqd),yyy(ntot/neqd,neqd)
     .             ,rhspc(ntot/neqd,neqd),rr(ntot),dxk(ntot)
     .             ,dvol,car(3),cnv(3),omg(neqd)

      real(8),allocatable,dimension(:,:,:,:) :: dv_cnv,dv_cov

      integer    :: ii,iii,igrid,ntotp,guess

      logical      :: hall

c Debug

      real(8),allocatable,dimension(:,:,:,:) :: dbg

c Externals

      external   tmp_mtvc,rho_mtvc,b_mtvc,v_mtvc,test_mtvc
     .          ,v_mtvc3,hall_mtvc,a_mtvc

c Begin program

cc      call setASMflag(asm_PC.and.(np>1))
cc
cc#if defined(petsc)
cc      if (iout >= 0) then
cc        write (*,*) 'Proc ',my_rank,': asm    in applyPC',asm
cc        write (*,*) 'Proc ',my_rank,': par_bc in applyPC',par_bc
cc      endif
cc#endif

      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      ntotp = ntot/neqd

      if (ntotp /= nx*ny*nz) then
        write (*,*)'Grid sizes do not agree in applyPreconditioner'
        write (*,*)'Aborting...'
        stop
      endif

      igrid = igx

cc      hall = (di > 0d0)
cc      hall = .true.
      hall = .false.  !Controls whether to invert M or not

c *******************************************************************
c     Identity preconditioner
c *******************************************************************

      if (precon.eq.'id') then

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

c *******************************************************************
c     Semi-implicit preconditioner version I (w/o flow)
c *******************************************************************

      elseif (precon == 's1') then

        xxx = 0d0

c     Scatter residuals

        call scatter(yyy,y)

        guess = 0

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3))

c diag****
cc        if (debug) then
cc        write (*,*) 'DIAG -- applyPC'
cccc        call cSolver(3,ntotp,rhspc(:,IVX:IVZ),xxx(:,IVX:IVZ)
cccc     .            ,bcs(:,IVX:IVZ),igrid,iout,guess,hall_mtvc
cccc     .            ,dg=v_diag,ncolors=4,line_relax=line_relax
cccc     .            ,order_res=2,order_prol=2,cvrg_tst=.true.
cccc     .            ,smooth='gm',omega=1d0)
cc        call cSolver(3,ntotp,rhspc(:,IVX:IVZ),xxx(:,IVX:IVZ)
cc     .               ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,dg=v_diag
cc     .               ,ncolors=4,line_relax=line_relax,order_res=0
cc     .               ,order_prol=2,smooth='gs',cvrg_tst=.true.
cc     .               ,omega=1d0,gm_driver=.false.
cc     .               ,gm_coarse_solve=.false.)
cc
cc        call gather(x,xxx)
cc
cc        call PCdiagnostics
cc        stop
cc        endif
c diag****

c     Preprocess B/A rhs

        if (di > 0d0) then
#if !defined(vec_pot)
          call preprocessBflRHS(yyy,igrid)
#else
          call preprocessAflRHS(yyy,igrid)
#endif
cc          call gather(rr,yyy)
cc        else
cc          rr = y
        endif

c     Split PC

        call split_PC

c     Postprocessing of velocity -> momentum

        if (.not.nc_eom_v) call get_dP

c     Postprocessing of B

#if !defined(vec_pot)
        if (solenoidal .and. (.not.hall))
     .       call divClean(ntotp,xxx(:,IBX:IBZ),igrid)
#endif

c     Gather solution in "x"

        call gather(x,xxx)

c     Diagnostics

#if defined(petsc_c)
        si_it = mod(si_it,max(precpass,1))+1
        if (si_it == debug_it .and. debug) call PCdiagnostics
#else
        if (jit == debug_it .and. debug) call PCdiagnostics
#endif

c     Deallocate variables

        deallocate(dv_cnv,dv_cov)

c *******************************************************************
c     Semi-implicit preconditioner version II (with flow)
c *******************************************************************

      elseif (precon == 's2') then

c     Set up vectors

        x   = 0d0
        xxx = 0d0

        guess = 0   !Equations are in residual form

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3))

c     Preprocess B rhs

cc        if (di > 0d0) then
cc          call preprocessBflRHS(y,igrid)
cc        endif

#if !defined(petsc_c)

        do si_it=1,precpass

c       Form residual vector rr=y-Ax

          if (si_it > 1) then
            call matrixFreeMatVec(ntot,x,rr)
            rr = y - rr
          else
            rr = y
          endif

          if (rat_limit) call norm_test

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
cc          if (solenoidal) call divClean(ntotp,xxx(:,IBX:IBZ),igrid)
#endif

c       Map solution xxx to x (jacobi iterate)

#if !defined(petsc_c)
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

        if (      iout >= 0
     .      .and. si_it > precpass
     .      .and. rat_limit       ) then
          call matrixFreeMatVec(ntot,x,rr)
          rr = y - rr
          call norm_test
        endif

#else

        call gather(x,xxx)

#endif

c     Diagnostics

        if (jit == debug_it .and. debug) call PCdiagnostics

c     Deallocate variables

        deallocate(dv_cnv,dv_cov)

c *******************************************************************
c     Semi-implicit preconditioner version III (with flow)
c *******************************************************************

      elseif (precon == 's3') then

c     Set up vectors

        xxx = 0d0

c     Scatter residuals

        call scatter(yyy,y)

        guess = 0

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3))

c     Preprocess B/A rhs

        if (di > 0d0) then
#if !defined(vec_pot)
          call preprocessBflRHS(yyy,igrid)
#else
          call preprocessAflRHS(yyy,igrid)
#endif
cc          call gather(rr,yyy)
cc        else
cc          rr = y
        endif

c     Commuted Schur complement PC

        call star_PC

c     Postprocessing of velocity -> momentum

        if (.not.nc_eom_v) call get_dP

c     Postprocessing of B

#if !defined(vec_pot)
        if (solenoidal .and. (.not.hall))
     .       call divClean(ntotp,xxx(:,IBX:IBZ),igrid)
#endif

c     Gather solution in "x"

        call gather(x,xxx)

c     Diagnostics

#if defined(petsc_c)
        si_it = mod(si_it,max(precpass,1))+1
        if (si_it == debug_it .and. debug) call PCdiagnostics
#else
        if (jit == debug_it .and. debug) call PCdiagnostics
#endif

c     Deallocate variables

        deallocate(dv_cnv,dv_cov)

      endif

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
        if (flag == 0) then

          call find_Udv(x,zz)

          x(:,IRHO)   = x(:,IRHO)    - dt*alpha/cnp*zz(:,IRHO)
          x(:,ITMP)   = x(:,ITMP)    - dt*alpha/cnp*zz(:,ITMP)
#if !defined(vec_pot)
          x(:,IBX:IBZ)= x(:,IBX:IBZ) - dt*alpha/cnp*zz(:,IBX:IBZ)
#else
          x(:,IAX:IAZ)= x(:,IAX:IAZ) - dt*alpha/cnp*zz(:,IAX:IAZ)
#endif
        elseif (flag == 1) then

          call correct_dx_rhs(x,rr,yyy)
          call invert_D_M(yyy,x)

        elseif (flag == 2) then
          zz = 0d0

          call correct_dx_rhs(x,zz,yyy)

          call invert_D_M(yyy,zz)

          x(:,IRHO) = x(:,IRHO) + zz(:,IRHO)
          x(:,ITMP) = x(:,ITMP) + zz(:,ITMP)
#if !defined(vec_pot)
          x(:,IBX:IBZ) = x(:,IBX:IBZ) + zz(:,IBX:IBZ)
#else
          x(:,IAX:IAZ) = x(:,IAX:IAZ) + zz(:,IAX:IAZ)
#endif
        elseif (flag == 3) then

          call correct_dx_rhs(x,rr,yyy)
          call invert_M(yyy,x)

        elseif (flag == 4) then
          zz = 0d0

          call correct_dx_rhs(x,zz,yyy)

          call invert_M(yyy,zz)

          x(:,IRHO) = x(:,IRHO) + zz(:,IRHO)
          x(:,ITMP) = x(:,ITMP) + zz(:,ITMP)
#if !defined(vec_pot)
          x(:,IBX:IBZ) = x(:,IBX:IBZ) + zz(:,IBX:IBZ)
#else
          x(:,IAX:IAZ) = x(:,IAX:IAZ) + zz(:,IAX:IAZ)
#endif

        else
          call pstop('find_dx'
     .              ,'Stationary iteration type not implemented')
        endif

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

        implicit none

c     Call variables

        real(8)   :: x(ntot/neqd,neqd),rhs(ntot/neqd,neqd)

c     Local variables

c     Begin program

c     Store velocity solution in array format --> dv_cnv

        call fillArray(igrid,IVX,3,x(:,IVX:IVZ),dv_cnv,bcs(:,IVX:IVZ)
     .                ,2,arr_cov=dv_cov,is_vec=.true.,is_cnv=.true.)

c     Correct rx to find Deltax (correction for rho, B, T)

        !Find corrections to rho,B,tmp right hand sides (Udv --> rhs)
        call correctRhoRHS(dv_cnv       ,rhs(:,IRHO)   ,igrid)
        call correctTmpRHS(dv_cnv       ,rhs(:,ITMP)   ,igrid)
#if !defined(vec_pot)
        call correctBflRHS(dv_cnv,dv_cov,rhs(:,IBX:IBZ),igrid)
#else
        call correctAflRHS(dv_cnv,dv_cov,rhs(:,IAX:IAZ),igrid)
#endif

c     End program

      end subroutine find_Udv

c     invert_M
c     ##############################################################
      subroutine invert_M(rhs,x)

c     --------------------------------------------------------------
c     Inverts r,B,T diagonal blocks. The variable 'hall' selects whether
c     to use MG (.false.) or a trivial 1/dt inverse (.true.).
c     --------------------------------------------------------------

        implicit none

c     Call variables

        real(8)    :: x(ntot/neqd,neqd),rhs(ntot/neqd,neqd),dvol

c     Local variables

c     Begin program

        if (.not.hall) then
          !Temperature
          if (bcSP()) then
            call cSolver(1,ntotp,rhs(:,ITMP),x(:,ITMP),bcs(:,ITMP)
     .                ,igrid,iout,guess,tmp_mtvc,dg=tmp_diag,ncolors=2
     .                ,line_relax=line_relax,order_res=0,order_prol=2
     .                ,smooth='gs')
          else
            call cSolver(1,ntotp,rhs(:,ITMP),x(:,ITMP),bcs(:,ITMP)
     .                ,igrid,iout,guess,tmp_mtvc,dg=tmp_diag
     .                ,ncolors=2,line_relax=line_relax,order_res=2
     .                ,order_prol=2,smooth='jb',gm_coarse_solve=.true.)
cccc     .                ,line_relax=line_relax,order_res=0,order_prol=2)
cc     .                ,line_relax=line_relax)
          endif

          !Density
          if (bcSP()) then
cc            write (*,*) 'DIAG -- applyPC'
            call cSolver(1,ntotp,rhs(:,IRHO),x(:,IRHO),bcs(:,IRHO)
     .                ,igrid,iout,guess,rho_mtvc,dg=rho_diag,ncolors=2
     .                ,line_relax=line_relax,order_res=0,order_prol=2
cc     .                ,galerkin=mg_glrkin,gm_coarse_solve=.true.
cc     .                ,galerkin=mg_glrkin,gm_coarse_solve=.false.
cc     .                ,line_relax=line_relax,order_res=0,order_prol=0
     .                ,smooth='gs')
          else
            call cSolver(1,ntotp,rhs(:,IRHO),x(:,IRHO),bcs(:,IRHO)
     .                ,igrid,iout,guess,rho_mtvc,dg=rho_diag
     .                ,ncolors=2,line_relax=line_relax,order_res=2
     .                ,order_prol=2,smooth='jb',gm_coarse_solve=.true.)
cc     .                ,line_relax=line_relax,order_res=0,order_prol=2)
cc     .                ,line_relax=line_relax)
          endif

#if !defined(vec_pot)
          !Magnetic field
          if (si_car) call XtoCar(igrid,rhs(:,IBX:IBZ))

          if (gm_smooth) then
            call cSolver(3,ntotp,rhs(:,IBX:IBZ),x(:,IBX:IBZ)
     .                ,bcs(:,IBX:IBZ),igrid,iout,guess,b_mtvc
     .                ,ncolors=4,smooth='gm',line_relax=line_relax)
          else
            if (bcSP()) then
              call cSolver(3,ntotp,rhs(:,IBX:IBZ),x(:,IBX:IBZ)
     .                ,bcs(:,IBX:IBZ),igrid,iout,guess,b_mtvc,dg=b_diag
     .                ,ncolors=4,line_relax=line_relax,order_res=0
     .                ,order_prol=2,smooth='gs')
            else
              call cSolver(3,ntotp,rhs(:,IBX:IBZ),x(:,IBX:IBZ)
     .                ,bcs(:,IBX:IBZ),igrid,iout,guess,b_mtvc,dg=b_diag
     .                ,ncolors=4,line_relax=line_relax,order_res=2
     .                ,order_prol=2,smooth='jb',gm_coarse_solve=.true.)
cc     .                ,ncolors=4,line_relax=line_relax,order_res=0
cc     .                ,order_prol=2)
cc     .                ,ncolors=4,line_relax=line_relax)
            endif
          endif

          if (si_car) call XtoCurv(igrid,x(:,IBX:IBZ))
#else
          if (si_car) call XtoCar(igrid,rhs(:,IAX:IAZ))

          if (gm_smooth) then
            call cSolver(3,ntotp,yyy(:,IAX:IAZ),xxx(:,IAX:IAZ)
     .              ,bcs(:,IAX:IAZ),igrid,iout,guess,a_mtvc
     .              ,ncolors=4,smooth='gm',line_relax=line_relax)
          else
            if (bcSP()) then
              call cSolver(3,ntotp,yyy(:,IAX:IAZ),xxx(:,IAX:IAZ)
     .              ,bcs(:,IAX:IAZ),igrid,iout,guess,a_mtvc,dg=a_diag
     .              ,ncolors=4,line_relax=line_relax,order_res=0
     .              ,order_prol=2,smooth='gs')
            else
              call cSolver(3,ntotp,yyy(:,IAX:IAZ),xxx(:,IAX:IAZ)
     .              ,bcs(:,IAX:IAZ),igrid,iout,guess,a_mtvc,dg=a_diag
     .              ,ncolors=4,line_relax=line_relax)
            endif
          endif

          if (si_car) call XtoCurv(igrid,x(:,IAX:IAZ))
#endif

        else

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                if (vol_wgt) then !Residuals are volume-weighed by default
                  dvol = gmetric%grid(igx)%dvol(i,j,k)
                else
                  dvol = 1d0
                endif

                x(ii,IRHO)    = dt*rhs(ii,IRHO)/cnp/dvol
                x(ii,ITMP)    = dt*rhs(ii,ITMP)/cnp/dvol
#if !defined(vec_pot)
                x(ii,IBX:IBZ) = dt*rhs(ii,IBX:IBZ)/cnp/dvol
#else
                x(ii,IAX:IAZ) = dt*rhs(ii,IAX:IAZ)/cnp/dvol
#endif               
              enddo
            enddo
          enddo

        endif

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

cc        if (.not.(jit == debug_it .and. debug)) then

        !Solve Schur-complement SI system
        if (gm_smooth) then
          call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc
     .                ,ncolors=4,smooth='gm',line_relax=line_relax)
        else
          if (bcSP()) then
cc            if (debug .and. jit == debug_it) then
cc              iout = max(3,iout)
cc              call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
cc     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc
cccc     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc_pr
cc     .                ,dg=v_diag,ncolors=4,order_res=2
cccc     .                ,order_prol=2,gm_driver=.true.,mvcyc=maxvcyc
cc     .                ,order_prol=2,gm_driver=.false.,mvcyc=maxvcyc
cc     .                ,smooth='jb',cvrg_tst=.true.)
cccc     .                ,order_prol=2,mvcyc=1,gm_drvr=.true.)
cccc     .                ,ncolors=4,line_relax=line_relax,cvrg_tst=debug)
cc            elseif (.not.si_car) then
            if (.not.si_car) then
              call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,dg=v_diag
cc     .                ,ncolors=4,line_relax=line_relax,order_res=2
cc     .                ,order_prol=2,gm_driver=.true.,galerkin=mg_glrkin)
     .                ,ncolors=4,line_relax=line_relax,order_res=0
     .                ,order_prol=2,smooth='gs',omega=1d0
     .                ,galerkin=mg_glrkin,gm_coarse_solve=.true.)
cc     .                ,order_prol=2,gm_drvr=.false.,galerkin=mg_glrkin)
cc     .                ,order_prol=2,gm_drvr=.false.,smooth='gs')
cc     .                ,order_prol=1,gm_drvr=.true.,mvcyc=2)
            else
              call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,dg=v_diag
cc     .                ,ncolors=4,line_relax=line_relax,order_res=2
cc     .                ,order_prol=2,gm_driver=.true.,galerkin=mg_glrkin)
     .                ,ncolors=4,line_relax=line_relax,order_res=2
cc     .                ,order_prol=2,smooth='gs',omega=1d0)
     .                ,order_prol=3,smooth='jb',omega=0.7d0)
cc     .                ,order_prol=2,gm_drvr=.false.,galerkin=mg_glrkin)
cc     .                ,order_prol=2,gm_drvr=.false.,smooth='gs')
cc     .                ,order_prol=1,gm_drvr=.true.,mvcyc=2)
            endif
          else
            call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,dg=v_diag
     .                ,ncolors=4,line_relax=line_relax,order_res=2
cc     .                ,ncolors=4,line_relax=line_relax,order_res=2
     .                ,order_prol=2,smooth='jb',gm_coarse_solve=.true.)
cc     .                ,order_prol=2,smooth='jb')
cc     .                ,ncolors=4,line_relax=line_relax,order_res=0
cc     .                ,order_prol=1,smooth='gm')
          endif
        endif

        if (si_car) call XtoCurv(igrid,x(:,IVX:IVZ))

cc        else
cc
cc          write (*,*) 'DIAG -- invert_P_si'
cc
cc          do k = 1,nz
cc            do j = 1,ny
cc              do i = 1,nx
cc                ii  = i + nx*(j-1) + nx*ny*(k-1)
cc                if (vol_wgt) then
cc                  dvol = gmetric%grid(igrid)%dvol(i,j,k)
cc                else
cc                  dvol = 1d0
cc                endif
cc                x(ii,IVX:IVZ) = dt*rhspc(ii,IVX:IVZ)/dvol/cnp
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

        real(8)   :: x(ntot/neqd,neqd),y(ntot/neqd,neqd)

c     Begin program

        !Form SI rhs (-Gv-Ldx)
        call formSIrhs(ntotp,x,y(:,IVX:IVZ),rhspc(:,IVX:IVZ),igrid)

        !Solve Schur-complement SI system
        if (gm_smooth) then
          call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc3
     .                ,ncolors=4,smooth='gm'
     .                ,line_relax=line_relax)
cc                      ,line_relax=line_relax,cvrg_tst=debug)
        else
          call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
     .             ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc3,dg=v_diag
     .             ,ncolors=4,line_relax=line_relax,order_res=2
     .             ,order_prol=2,smooth='jb',gm_coarse_solve=.true.)
cc     .             ,ncolors=4,line_relax=line_relax)
cc     .            ,ncolors=4,line_relax=line_relax,cvrg_tst=debug)
        endif

        !Post-process solution with M_star

        call postProcessSchur(3,ntotp,x(:,IVX:IVZ),igrid
     .                       ,bcs(:,IVX:IVZ))

        !Transform coordinates
        if (si_car) call XtoCurv(igrid,x(:,IVX:IVZ))

cc        !Form SI rhs -(Gv+Ldx)
cc        call formSIrhs(ntotp,x,y(:,IVX:IVZ),rhspc(:,IVX:IVZ),igrid)
cc
cccc        if (.not.(jit == debug_it .and. debug)) then
cc
cc        !Solve Schur-complement SI system
cc        if (gm_smooth) then
cc          call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
cc     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc
cc     .                ,ncolors=4,smooth='gm',line_relax=line_relax)
cc        else
cc          if (bcSP()) then
cccc            if (debug .and. jit == debug_it) then
cccc              iout = max(3,iout)
cccc              call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
cccc     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc
cccccc     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc_pr
cccc     .                ,dg=v_diag,ncolors=4,order_res=2
cccccc     .                ,order_prol=2,gm_driver=.true.,mvcyc=maxvcyc
cccc     .                ,order_prol=2,gm_driver=.false.,mvcyc=maxvcyc
cccc     .                ,smooth='jb',cvrg_tst=.true.)
cccccc     .                ,order_prol=2,mvcyc=1,gm_drvr=.true.)
cccccc     .                ,ncolors=4,line_relax=line_relax,cvrg_tst=debug)
cccc            elseif (.not.si_car) then
cc            if (.not.si_car) then
cc              call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
cc     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,dg=v_diag
cccc     .                ,ncolors=4,line_relax=line_relax,order_res=2
cccc     .                ,order_prol=2,gm_driver=.true.,galerkin=mg_glrkin)
cc     .                ,ncolors=4,line_relax=line_relax,order_res=0
cc     .                ,order_prol=2,smooth='gs',omega=1d0
cc     .                ,galerkin=mg_glrkin,gm_coarse_solve=.true.)
cccc     .                ,order_prol=2,gm_drvr=.false.,galerkin=mg_glrkin)
cccc     .                ,order_prol=2,gm_drvr=.false.,smooth='gs')
cccc     .                ,order_prol=1,gm_drvr=.true.,mvcyc=2)
cc            else
cc              call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
cc     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,dg=v_diag
cccc     .                ,ncolors=4,line_relax=line_relax,order_res=2
cccc     .                ,order_prol=2,gm_driver=.true.,galerkin=mg_glrkin)
cc     .                ,ncolors=4,line_relax=line_relax,order_res=2
cccc     .                ,order_prol=2,smooth='gs',omega=1d0)
cc     .                ,order_prol=3,smooth='jb',omega=0.7d0)
cccc     .                ,order_prol=2,gm_drvr=.false.,galerkin=mg_glrkin)
cccc     .                ,order_prol=2,gm_drvr=.false.,smooth='gs')
cccc     .                ,order_prol=1,gm_drvr=.true.,mvcyc=2)
cc            endif
cc          else
cc            call cSolver(3,ntotp,rhspc(:,IVX:IVZ),x(:,IVX:IVZ)
cc     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,dg=v_diag
cc     .                ,ncolors=4,line_relax=line_relax,order_res=2
cc     .                ,order_prol=2,smooth='jb',gm_coarse_solve=.true.)
cccc     .                ,ncolors=4,line_relax=line_relax,order_res=2
cccc     .                ,order_prol=2,smooth='jb')
cccc     .                ,ncolors=4,line_relax=line_relax,order_res=0
cccc     .                ,order_prol=1,smooth='gm')
cc          endif
cc        endif
cc
cc        if (si_car) call XtoCurv(igrid,x(:,IVX:IVZ))
cc
cccc        else
cc
cccc          write (*,*) 'DIAG -- invert_P_si'
cccc
cccc          do k = 1,nz
cccc            do j = 1,ny
cccc              do i = 1,nx
cccc                ii  = i + nx*(j-1) + nx*ny*(k-1)
cccc                if (vol_wgt) then
cccc                  dvol = volume(i,j,k,igrid,igrid,igrid)
cccc                else
cccc                  dvol = 1d0
cccc                endif
cccc                x(ii,IVX:IVZ) = dt*rhspc(ii,IVX:IVZ)/dvol/cnp
cccc              enddo
cccc            enddo
cccc          enddo
cc
cccc        endif

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
             xxx(ii,IVY)=rho(i,j,k)*xxx(ii,IVY)
     .                  +xxx(ii,IRHO)*vcnv(i,j,k,2)
             xxx(ii,IVZ)=rho(i,j,k)*xxx(ii,IVZ)
     .                  +xxx(ii,IRHO)*vcnv(i,j,k,3)
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

        call find_dx(xxx,rr,0)

cc        if (precon == 's2') then
cccc          call find_Udv(xxx,rhspc)
cccc
cccc          call invert_D_M(rhspc,dum)
cccc
cccc          do k = 1,nz
cccc            do j = 1,ny
cccc              do i = 1,nx
cccc                ii  = i + nx*(j-1) + nx*ny*(k-1)
cccc
cccc                iii = neqd*(ii-1)
cccc
cccc                if (vol_wgt) then !Residuals are volume-weighed by default
cccc                  dvol = volume(i,j,k,igx,igy,igz)
cccc                else
cccc                  dvol = 1d0
cccc                endif
cccc
cccc                xxx(ii,IRHO)   =xxx(iii+IRHO)-dvol*alpha*dum(ii,IRHO)
cccc                xxx(ii,IBX:IBZ)=xxx(iii+IBX:iii+IBZ)
cccc     .                                       -dvol*alpha*dum(ii,IBX:IBZ)
cccc                xxx(ii,ITMP)   =xxx(iii+ITMP)-dvol*alpha*dum(ii,ITMP)
cccc
cccc              enddo
cccc            enddo
cccc          enddo
cc
cccc          xxx(:,IRHO)    = xxx(:,IRHO)    - dum(:,IRHO)
cccc          xxx(:,ITMP)    = xxx(:,ITMP)    - dum(:,ITMP)
cccc          xxx(:,IBX:IBZ) = xxx(:,IBX:IBZ) - dum(:,IBX:IBZ)
cc
cc          call find_dx(xxx,rr,1)
cc        else
cc          call find_Udv(xxx,rhspc)
cc
cc          xxx(:,IRHO)   = xxx(:,IRHO)    - dt*alpha/cnp*rhspc(:,IRHO)
cc          xxx(:,ITMP)   = xxx(:,ITMP)    - dt*alpha/cnp*rhspc(:,ITMP)
cc#if !defined(vec_pot)
cc          xxx(:,IBX:IBZ)= xxx(:,IBX:IBZ) - dt*alpha/cnp*rhspc(:,IBX:IBZ)
cc#else
cc          xxx(:,IAX:IAZ)= xxx(:,IAX:IAZ) - dt*alpha/cnp*rhspc(:,IAX:IAZ)
cc#endif
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
        mag = sqrt(sum(rr*rr))

        if (si_it == 1) then
          mag0  = mag
          mag_n = -1d0
        endif

        if (iout >= 1) write (*,*)

        ratio = max(mag/mag_n,0d0)

        if (iout >= 0) 
     .    write (*,100) '   SI iter:',si_it-1,';  Res: ',mag
     .                 ,';  Rel.res.:',mag/mag0,';  Ratio:',ratio

        mag_n = mag

        !Individual residual eqs
        if (iout >=2) then
          call scatter(yyy,rr)
          do ieq=1,neqd
            mag = sqrt(sum(yyy(:,ieq)**2))
            write (*,*) '     Res. eqn',ieq,':',mag
          enddo
          write (*,*)
        endif

 100    format (a,i3,a,1p,e10.2,a,e10.2,a,e10.2)

      end subroutine norm_test

c     PCdiagnostics
c     ##############################################################
      subroutine PCdiagnostics

c     --------------------------------------------------------------
c     Performs diagnostics on PC
c     --------------------------------------------------------------

        implicit none

        character(1) :: plot

        real(8) :: mag0,mag,lxmin,lxmax,lymin,lymax

        real(8),allocatable,dimension(:,:,:,:) :: db_cnv,db_cov

        character(10) :: solf,resf,rhsf,crsf,fldf,divf

#if defined(vec_pot)
        integer    :: ivar

        real(8),allocatable,dimension(:,:,:,:) :: da_cov

        allocate(da_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,db_cov(0:nx+1,0:ny+1,0:nz+1,3))
#endif

c     Begin program

        lxmin = grid_params%lxmin
        lxmax = grid_params%lxmax

        lymin = grid_params%lymin
        lymax = grid_params%lymax

c     Create output files

        solf='sol_p'//trim(int2char(my_rank))//'.bin'
        resf='res_p'//trim(int2char(my_rank))//'.bin'
        rhsf='rhs_p'//trim(int2char(my_rank))//'.bin'
        crsf='crs_p'//trim(int2char(my_rank))//'.bin'
        fldf='fld_p'//trim(int2char(my_rank))//'.bin'
        divf='div_p'//trim(int2char(my_rank))//'.bin'

c     Dump graphics

        allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3))

        !Solution plot
        call MGplot(neqd,x,igrid,0,solf)

        !Predictor RHS plot
        call MGplot(neqd,y,igrid,0,rhsf)

        !Residual plot
#if !defined(petsc)
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
        call fillArray(igrid,IBX,3,xxx(:,IBX:IBZ),db_cnv,bcs(:,IBX:IBZ)
     .                ,2,is_vec=.true.,is_cnv=.true.)
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

        !Find dv (w/ BCs)
        call fillArray(igrid,IVX,3,xxx(:,IVX:IVZ),dv_cnv,bcs(:,IVX:IVZ)
     .                ,2,is_vec=.true.,is_cnv=.true.)

        open(unit=110,file=fldf,form='unformatted'
     .     ,status='replace')
        call contour(dv_cnv(0:nx+1,1:ny,1,1),nx+2,ny,lxmin
     .              ,lxmax,lymin,lymax,0,110)
        call contour(dv_cnv(0:nx+1,1:ny,1,2),nx+2,ny,lxmin
     .              ,lxmax,lymin,lymax,1,110)
        call contour(dv_cnv(0:nx+1,1:ny,1,3),nx+2,ny,lxmin
     .              ,lxmax,lymin,lymax,1,110)
#if !defined(vec_pot)
        call contour(db_cnv(0:nx+1,1:ny,1,1),nx+2,ny,lxmin
     .              ,lxmax,lymin,lymax,1,110)
        call contour(db_cnv(0:nx+1,1:ny,1,2),nx+2,ny,lxmin
     .              ,lxmax,lymin,lymax,1,110)
        call contour(db_cnv(0:nx+1,1:ny,1,3),nx+2,ny,lxmin
     .              ,lxmax,lymin,lymax,1,110)
#else
        call contour(da_cov(0:nx+1,1:ny,1,1),nx+2,ny,lxmin
     .              ,lxmax,lymin,lymax,1,110)
        call contour(da_cov(0:nx+1,1:ny,1,2),nx+2,ny,lxmin
     .              ,lxmax,lymin,lymax,1,110)
        call contour(da_cov(0:nx+1,1:ny,1,3),nx+2,ny,lxmin
     .              ,lxmax,lymin,lymax,1,110)
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
        call contour(dbg(:,:,1,1),nx,ny,lxmin,lxmax,lymin,lymax
     .              ,0,110)
        close(110)

        deallocate(dbg)

        call pstop('PCdiagnostics','Finished PC diagnostics')

      end subroutine PCdiagnostics

      end subroutine applyPreconditioner

c findEOMrhs
c #########################################################################
      subroutine findEOMrhs(ntotp,xxx,eom_rhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the EM part of the rhs for the velocity
c     SI solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use imposeBCinterface

cc      use mg_internal

      implicit none

c Call variables

      integer    :: ntotp,igrid
      real(8)    :: xxx(ntotp,neqd),eom_rhs(ntotp,3)

c Local variables

      integer    :: ii,isig,ivar,ip,im
      real(8)    :: dvol,cov(3),cnv(3),idx,idy,idz
     .             ,dj_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,db_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,drho  (0:nx+1,0:ny+1,0:nz+1,1)
     .             ,dtmp  (0:nx+1,0:ny+1,0:nz+1,1)
#if defined(vec_pot)
     .             ,da_cov(0:nx+1,0:ny+1,0:nz+1,3)
#endif

c Begin program

      isig = grid_params%istartp(igrid)

#if !defined(vec_pot)
c Find dB* (w/ BCs)

      call fillArray(igrid,IBX,3,xxx(:,IBX:IBZ),db_cnv,bcs(:,IBX:IBZ)
     .              ,2,arr_cov=db_cov,is_vec=.true.,is_cnv=.true.)

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

c Find rhs_v'

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jac = gmetric%grid(igrid)%jac(i,j,k)

            if (vol_wgt) then
              dvol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            !Pressure correction
            ip = i+1
            im = i-1

            idx = 0.5/dxh(ig)
            idy = 0.5/dyh(jg)
            idz = 0.5/dzh(kg)

cc            cov(1) = jy(i,j,k)*db_cnv(i,j,k,3)/jac
cc     .             - jz(i,j,k)*db_cnv(i,j,k,2)/jac
cc     .             - by(i,j,k)*dj_cnv(i,j,k,3)/jac
cc     .             + bz(i,j,k)*dj_cnv(i,j,k,2)/jac
cccc     .             -(dpres(ip,j,k)-dpres(im,j,k))*idx
cc     .             -a_p*(rho(i,j,k)*(dtmp(i+1,j,k)-dtmp(i-1,j,k))*idx
cc     .                  +tmp(i,j,k)*(drho(i+1,j,k)-drho(i-1,j,k))*idx)
cc
cc            cov(2) = jz(i,j,k)*db_cnv(i,j,k,1)/jac
cc     .             - jx(i,j,k)*db_cnv(i,j,k,3)/jac
cc     .             - bz(i,j,k)*dj_cnv(i,j,k,1)/jac
cc     .             + bx(i,j,k)*dj_cnv(i,j,k,3)/jac
cccc     .             -(dpres(i,j+1,k)-dpres(i,j-1,k))*idy
cc     .             -a_p*(rho(i,j,k)*(dtmp(i,j+1,k)-dtmp(i,j-1,k))*idy
cc     .                  +tmp(i,j,k)*(drho(i,j+1,k)-drho(i,j-1,k))*idy)
cc
cc            cov(3) = jx(i,j,k)*db_cnv(i,j,k,2)/jac
cc     .             - jy(i,j,k)*db_cnv(i,j,k,1)/jac
cc     .             - bx(i,j,k)*dj_cnv(i,j,k,2)/jac
cc     .             + by(i,j,k)*dj_cnv(i,j,k,1)/jac
cccc     .             -(dpres(i,j,k+1)-dpres(i,j,k-1))*idz
cc     .             -a_p*(rho(i,j,k)*(dtmp(i,j,k+1)-dtmp(i,j,k-1))*idz
cc     .                  +tmp(i,j,k)*(drho(i,j,k+1)-drho(i,j,k-1))*idz)

            cov(1) =
     .         -a_p*(rho(i,j,k)*(dtmp(i+1,j,k,1)-dtmp(i-1,j,k,1))*idx
     .              +tmp(i,j,k)*(drho(i+1,j,k,1)-drho(i-1,j,k,1))*idx)

            cov(2) =
     .         -a_p*(rho(i,j,k)*(dtmp(i,j+1,k,1)-dtmp(i,j-1,k,1))*idy
     .              +tmp(i,j,k)*(drho(i,j+1,k,1)-drho(i,j-1,k,1))*idy)

            cov(3) =
     .         -a_p*(rho(i,j,k)*(dtmp(i,j,k+1,1)-dtmp(i,j,k-1,1))*idz
     .              +tmp(i,j,k)*(drho(i,j,k+1,1)-drho(i,j,k-1,1))*idz)

            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3),.true.)

            !Add Lorentz force correction
            cnv(1) = cnv(1)
     .             + jcov(i,j,k,2)*db_cov(i,j,k,3)
     .             - jcov(i,j,k,3)*db_cov(i,j,k,2)
     .             - bcov(i,j,k,2)*dj_cov(i,j,k,3)
     .             + bcov(i,j,k,3)*dj_cov(i,j,k,2)

            cnv(2) = cnv(2)
     .             + jcov(i,j,k,3)*db_cov(i,j,k,1)
     .             - jcov(i,j,k,1)*db_cov(i,j,k,3)
     .             - bcov(i,j,k,3)*dj_cov(i,j,k,1)
     .             + bcov(i,j,k,1)*dj_cov(i,j,k,3)

            cnv(3) = cnv(3)
     .             + jcov(i,j,k,1)*db_cov(i,j,k,2)
     .             - jcov(i,j,k,2)*db_cov(i,j,k,1)
     .             - bcov(i,j,k,1)*dj_cov(i,j,k,2)
     .             + bcov(i,j,k,2)*dj_cov(i,j,k,1)

            !Correct rhs_v
            if (nc_eom_v) then
              eom_rhs(ii,:) = alpha*dvol*cnv(:)/rho(i,j,k)
            else
              eom_rhs(ii,:) = alpha*dvol*cnv(:)
            endif
          enddo
        enddo
      enddo

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

cc      use mg_internal

      implicit none

c Call variables

      integer    :: ntotp,igrid
      real(8)    :: xxx(ntotp,neqd),yyy(ntotp,3),rhs_si(ntotp,3)

c Local variables

      integer    :: ii,iig,isig,ivar,ip,im
      real(8)    :: dvol,cov(3),cnv(3),idx,idy,idz,eom_rhs(ntotp,3)

c Begin program

      isig = grid_params%istartp(igrid)

c Find rhs_v

      if (.not.nc_eom_v) then
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
      else
        rhs_si = yyy
      endif

c Find linearized EOM rhs [jxB-grad(p)] correction

      call findEOMrhs(ntotp,xxx,eom_rhs,igrid)

c Add rhs_v + EOM rhs

      rhs_si = rhs_si + eom_rhs

c Transform to cartesian if applicable

      if (si_car) call XtoCar(igrid,rhs_si)

c End program

      end subroutine formSIrhs

c correctVfldRHS
c #########################################################################
      subroutine correctVfldRHS(ntotp,dv,xxx,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the rhs for the velocity SI solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

cc      use mg_internal

      implicit none

c Call variables

      integer    :: ntotp,igrid
      real(8)    :: xxx(ntotp,neqd),crhs(ntotp,3)
     .             ,dv(0:nx+1,0:ny+1,0:nz+1,3)

c Local variables

      integer    :: ii,iig,isig,ivar,ip,im
      real(8)    :: dvol,cov(3),cnv(3),idx,idy,idz,eom_rhs(ntotp,3)
     .             ,nabla_vv0(3,3),mul,veclap(3),advc(3)

c Begin program

      isig = grid_params%istartp(igrid)

c Find linearized EOM rhs [jxB-grad(p)] part

      call findEOMrhs(ntotp,xxx,eom_rhs,igrid)

c Find advection-diffusion part

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            iig = ii + isig - 1

            jac = gmetric%grid(igrid)%jac(i,j,k)

            if (vol_wgt) then
              dvol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            !EOM residual part + EM part
            cnv = -eom_rhs(ii,:)/dvol/alpha

            !First advection contribution
            nabla_v   = fnabla_v(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,dv(:,:,:,1),dv(:,:,:,2),dv(:,:,:,3),0)

            nabla_vv0 = fnabla_v(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                          ,vcnv(:,:,:,1)
     .                          ,vcnv(:,:,:,2)
     .                          ,vcnv(:,:,:,3),0)

            mul = vis(i,j,k,nx,ny,nz,igrid,igrid,igrid)

            veclap = veclaplacian(i,j,k,nx,ny,nz,igrid,igrid,igrid,dv
     .                           ,alt_eom,vol=.false.)

            advc  =( vcnv(i,j,k,1)*nabla_v  (1,:)
     .              +vcnv(i,j,k,2)*nabla_v  (2,:)
     .              +vcnv(i,j,k,3)*nabla_v  (3,:)
     .              +dv  (i,j,k,1)*nabla_vv0(1,:)
     .              +dv  (i,j,k,2)*nabla_vv0(2,:)
     .              +dv  (i,j,k,3)*nabla_vv0(3,:))/jac

            cnv = cnv + rho(i,j,k)*(advc - mul*veclap)

            !Second advection contribution
            veclap = veclaplacian(i,j,k,nx,ny,nz,igrid,igrid,igrid,vcnv
     .                           ,alt_eom,vol=.false.)

            advc  =( vcnv(i,j,k,1)*nabla_vv0(1,:)
     .              +vcnv(i,j,k,2)*nabla_vv0(2,:)
     .              +vcnv(i,j,k,3)*nabla_vv0(3,:))/jac

            cnv = cnv + xxx(ii,IRHO)*(advc - mul*veclap)

            crhs(ii,:) = cnv

          enddo
        enddo
      enddo

c End program

      end subroutine correctVfldRHS

c correctRhoRHS
c #########################################################################
      subroutine correctRhoRHS(dv_cnv,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the density solve,
c     div(rh0*dv)
c--------------------------------------------------------------------

      use precond_variables

      use operators

cc      use mg_internal

      implicit none

c Call variables

      integer    :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz)

c Local variables

      integer    :: ii
cc      real(8)    :: rhodv_cnv(0:nx+1,0:ny+1,0:nz+1,3)

c Begin program

cc      rhodv_cnv(:,:,:,1) = dv_cnv(:,:,:,1)*rho
cc      rhodv_cnv(:,:,:,2) = dv_cnv(:,:,:,2)*rho
cc      rhodv_cnv(:,:,:,3) = dv_cnv(:,:,:,3)*rho

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            crhs(ii) = flx_advec(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                      ,dv_cnv(:,:,:,1)
     .                      ,dv_cnv(:,:,:,2)
     .                      ,dv_cnv(:,:,:,3),rho,4,vol=.false.)

cc            crhs(ii) = div(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                    ,rhodv_cnv(:,:,:,1)
cc     .                    ,rhodv_cnv(:,:,:,2)
cccc     .                    ,rhodv_cnv(:,:,:,3),sp=bcSP())
cc     .                    ,rhodv_cnv(:,:,:,3))

cc            crhs(ii) = ( rho(i,j,k)
cc     $                           *div(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     $                               ,dv_cnv(:,:,:,1)
cc     $                               ,dv_cnv(:,:,:,2)
cc     $                               ,dv_cnv(:,:,:,3))
cc     $                   +dv_cnv(i,j,k,1)/jac
cc     $                       *(rho(i+1,j,k)-rho(i-1,j,k))/dxh(ig)/2.
cc     $                   +dv_cnv(i,j,k,2)/jac
cc     $                       *(rho(i,j+1,k)-rho(i,j-1,k))/dyh(jg)/2.
cc     $                   +dv_cnv(i,j,k,3)/jac
cc     $                       *(rho(i,j,k+1)-rho(i,j,k-1))/dzh(kg)/2. )

cc            crhs(ii) = 0d0
 
          enddo
        enddo
      enddo

c End program

      end subroutine correctRhoRHS

c correctTmpRHS
c #########################################################################
      subroutine correctTmpRHS(dv_cnv,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the temperature
c     solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

cc      use mg_internal

      implicit none

c Call variables

      integer    :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz)

c Local variables

      integer    :: ii

c Begin program

c Evaluate rhs correction: dv*grad(T0) + (gamma-1)*T0*div(dv)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jac = gmetric%grid(igrid)%jac(i,j,k)

            crhs(ii) = ( (gamma-1)*tmp(i,j,k)
     $                           *div(i,j,k,nx,ny,nz,igrid,igrid,igrid
     $                               ,dv_cnv(:,:,:,1)
     $                               ,dv_cnv(:,:,:,2)
     $                               ,dv_cnv(:,:,:,3))
     $                   +dv_cnv(i,j,k,1)/jac
     $                       *(tmp(i+1,j,k)-tmp(i-1,j,k))/dxh(ig)/2.
     $                   +dv_cnv(i,j,k,2)/jac
     $                       *(tmp(i,j+1,k)-tmp(i,j-1,k))/dyh(jg)/2.
     $                   +dv_cnv(i,j,k,3)/jac
     $                       *(tmp(i,j,k+1)-tmp(i,j,k-1))/dzh(kg)/2. )
          enddo
        enddo
      enddo

c End program

      end subroutine correctTmpRHS

c correctBflRHS
c ###################################################################
      subroutine correctBflRHS(dv_cnv,dv_cov,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the magnetic field
c     solve: div(dv B0 - B0 dv)
c--------------------------------------------------------------------

      use precond_variables

      use operators

cc      use mg_internal

      use matvec

      implicit none

c Call variables

      integer    :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz,3)

c Local variables

cc      real(8)    :: dummy  (0:nx+1,0:ny+1,0:nz+1,3)
cc     .             ,dum_cov(0:nx+1,0:ny+1,0:nz+1,3)
      integer    :: ii

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            crhs(ii,:) = vecA(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                       ,dv_cnv,dv_cov
     .                       ,gb0%grid(igrid)%array,0)
          enddo
        enddo
      enddo

      return

cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            ii  = i + nx*(j-1) + nx*ny*(k-1)
cc            crhs(ii,:) = curl_bxv(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                       ,dv_cnv
cc     .                       ,gb0%grid(igrid)%array,0)
cc          enddo
cc        enddo
cc      enddo
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
      subroutine correctAflRHS(dv_cnv,dv_cov,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the vector potential
c     solve: div(dv B0 - B0 dv)
c--------------------------------------------------------------------

      use precond_variables

      use operators

cc      use mg_internal

      implicit none

c Call variables

      integer    :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dv_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,crhs(nx*ny*nz,3)

c Local variables

      integer    :: ii
      real(8)    :: minv

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            jac = gmetric%grid(igrid)%jac(i,j,k)

            crhs(ii,1) = dv_cnv(i,j,k,2)*gb0%grid(igrid)%array(i,j,k,3)
     .                 - dv_cnv(i,j,k,3)*gb0%grid(igrid)%array(i,j,k,2)

            crhs(ii,2) = dv_cnv(i,j,k,3)*gb0%grid(igrid)%array(i,j,k,1)
     .                 - dv_cnv(i,j,k,1)*gb0%grid(igrid)%array(i,j,k,3)

            crhs(ii,3) = dv_cnv(i,j,k,1)*gb0%grid(igrid)%array(i,j,k,2)
     .                 - dv_cnv(i,j,k,2)*gb0%grid(igrid)%array(i,j,k,1)

            crhs(ii,:) = -crhs(ii,:)/jac

            !Hall MHD part
            if (di > 0d0) then
              minv = di*cnp/dt/alpha
              crhs(ii,:) = crhs(ii,:) + minv*dv_cov(i,j,k,:)
            endif

          enddo
        enddo
      enddo


c End program

      end subroutine correctAflRHS

c preprocessBflRHS
c ###################################################################
      subroutine preprocessBflRHS(yyy,igrid)
c--------------------------------------------------------------------
c     This subroutine preprocesses B rhs to add Hall correction.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use matvec

      implicit none

c Call variables

      integer    :: igrid

      real(8)    :: yyy(nx*ny*nz,neqd)

c Local variables

      integer    :: ii

      real(8)    :: dgv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dgv_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,cnv(3)

c Begin program

c Store velocity RHS in array format w/ BCs

      call fillArray(igrid,IVX,3,yyy(:,IVX:IVZ),dgv_cnv,bcs(:,IVX:IVZ)
     .              ,2,arr_cov=dgv_cov,is_vec=.true.,is_cnv=.true.)

c Divide by density

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
            yyy(ii,IBX:IBZ) = yyy(ii,IBX:IBZ)
     .                        + di*curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                                 ,dgv_cov(:,:,:,1)
     .                                 ,dgv_cov(:,:,:,2)
     .                                 ,dgv_cov(:,:,:,3))
          enddo
        enddo
      enddo

c End program

      end subroutine preprocessBflRHS

c preprocessAflRHS
c ###################################################################
      subroutine preprocessAflRHS(yyy,igrid)
c--------------------------------------------------------------------
c     This subroutine preprocesses B rhs to add Hall correction.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use matvec

cc      use mg_internal

      implicit none

c Call variables

      integer    :: igrid

      real(8)    :: yyy(nx*ny*nz,neqd)

c Local variables

      integer    :: ii
      real(8)    :: cov(3),cnv(3)

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            cnv = yyy(ii,IVX:IVZ)
            call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .           ,cov(1),cov(2),cov(3)
     .           ,cnv(1),cnv(2),cnv(3),.false.)

            yyy(ii,IAX:IAZ) = yyy(ii,IAX:IAZ) + di*cov/rho(i,j,k)
          enddo
        enddo
      enddo

c End program

      end subroutine preprocessAflRHS

c divClean
c ###################################################################
      subroutine divClean(ntotp,xxx,igrid)
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

      integer    :: ntotp,igrid
      real(8)    :: xxx(ntotp,3)

c Local variables

      integer    :: ii,ivar,bcnd(6,3)
      real(8)    :: a1,a2,a3,etal
     .             ,db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,da_cov(0:nx+1,0:ny+1,0:nz+1,3)

c Begin program

cc      write (*,*) 'DIAG -- divClean: no divergence cleaning'
cc      return

c Find dB* (w/ BCs)

      call fillArray(igrid,IBX,3,xxx,db_cnv,bcs(:,IBX:IBZ)
     .              ,2,is_vec=.true.,is_cnv=.true.)

c Find dA* from dB* (w/ ghost cells)

      call curl_inv(nx,ny,nz,igrid,igrid,igrid,db_cnv,da_cov)

#if !defined(petsc)
      if (isSP(1,1,1,igrid,igrid,igrid)) then
        !Impose BCs
        bcnd(:,1) = bcond
        bcnd(:,2) = bcond
        bcnd(:,3) = bcond
        where (bcnd /= SP) bcnd = DEF
cc        where (bcnd == DEF) bcnd = EXT

        call setMGBC(0,3,nx,ny,nz,igrid,da_cov,bcnd,icomp=IBX
     .              ,is_vec=.true.,is_cnv=.false.,iorder=2)
cc        call setBC(IBX,3,nx,ny,nz,da_cov,da_cnv,vzeros,bcnd,igrid,igrid,igrid
cc     .            ,is_cnv=.false.,iorder=2)

      endif
#endif

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

cc      use mg_internal

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
