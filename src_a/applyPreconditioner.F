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

      use mg_internal

      implicit none

c Call variables

      integer(4) :: ntot,iout
      real(8)    :: x(ntot),y(ntot)

c Local variables

      real(8)    :: xxx(ntot/neqd,neqd),yyy(ntot/neqd,neqd)
     .             ,rhs(ntot/neqd,neqd),rr(ntot),dxk(ntot)
     .             ,dvol,car(3),cnv(3),omega

      real(8),allocatable,dimension(:,:,:,:) :: dv_cnv,db_cnv,dj_cov

      integer(4) :: ii,iii,igrid,ntotp,guess

      real(8)      :: mag
      real(8),save :: mag0

c Debug

      real(8),allocatable,dimension(:,:,:,:) :: dbg

c Externals

      external   tmp_mtvc,rho_mtvc,a_mtvc,v_mtvc,test_mtvc,v_mtvc2

c Begin program

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
                if (vol_wgt) then
                  dvol = volume(i,j,k,igx,igy,igz)
                else
                  dvol = 1d0
                endif
                x(iii) = y(iii)*dt/dvol
              enddo
            enddo
          enddo
        enddo

c *******************************************************************
c     Semi-implicit preconditioner version I
c *******************************************************************

      elseif (precon == 's1') then

        xxx = 0d0

c     Scatter residuals

        call scatter(yyy,y)

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3))

c     Predictor step

        guess = 0

        !Temperature
        call cSolver(1,ntotp,yyy(:,ITMP),xxx(:,ITMP),bcs(:,ITMP)
     .              ,igrid,iout,guess,tmp_mtvc,dg=tmp_diag,ncolors=2
     .              ,line_relax=line_relax)

        !Density
        call cSolver(1,ntotp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
     .              ,igrid,iout,guess,rho_mtvc,dg=rho_diag,ncolors=2
     .              ,line_relax=line_relax)

        !Vector pot
        if (gm_smooth) then
          call cSolver(3,ntotp,yyy(:,IAX:IAZ),xxx(:,IAX:IAZ)
     .              ,bcs(:,IAX:IAZ),igrid,iout,guess,a_mtvc
     .              ,ncolors=4,gm_smth=gm_smooth,line_relax=line_relax)
        else
          call cSolver(3,ntotp,yyy(:,IAX:IAZ),xxx(:,IAX:IAZ)
     .              ,bcs(:,IAX:IAZ),igrid,iout,guess,a_mtvc,dg=a_diag
     .              ,ncolors=4,line_relax=line_relax)
        endif

cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc              xxx(ii,IRHO)    = dt*yyy(ii,IRHO)
cc     .                           /volume(i,j,k,igx,igy,igz)
cc              xxx(ii,ITMP)    = dt*yyy(ii,ITMP)
cc     .                           /volume(i,j,k,igx,igy,igz)
cc              xxx(ii,IAX:IAZ) = dt*yyy(ii,IAX:IAZ)
cc     .                           /volume(i,j,k,igx,igy,igz)
cc            enddo
cc          enddo
cc        enddo

c     SI step

        !Form SI rhs -(Gv+Ldx)
        call formSIrhs(ntotp,xxx,yyy(:,IVX:IVZ),rhs(:,IVX:IVZ),igrid)

        !Solve Schur-complement SI system
        if (gm_smooth) then
          call cSolver(3,ntotp,rhs(:,IVX:IVZ),xxx(:,IVX:IVZ)
     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc
     .                ,ncolors=4,gm_smth=gm_smooth
     .                ,line_relax=line_relax)
cc     .                ,line_relax=line_relax,cvrg_tst=debug)
        else
          call cSolver(3,ntotp,rhs(:,IVX:IVZ),xxx(:,IVX:IVZ)
     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,dg=v_diag
cc     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc2,dg=v_diag2
     .                ,ncolors=4,line_relax=line_relax)
cc     .                ,ncolors=4,line_relax=line_relax,cvrg_tst=debug)
        endif

cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc              xxx(ii,IVX:IVZ) = dt*rhs(ii,IVX:IVZ)
cc     .                         /volume(i,j,k,igx,igy,igz)
cc            enddo
cc          enddo
cc        enddo

        if (si_car) then
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                car = xxx(ii,IVX:IVZ)
                call transformVectorToCurvilinear
     .               (i,j,k,igrid,igrid,igrid
     .               ,car(1),car(2),car(3)
     .               ,.false.
     .               ,xxx(ii,IVX),xxx(ii,IVY),xxx(ii,IVZ))
              enddo
            enddo
          enddo
        endif

c     Store velocity solution in array format

        !xxx is NOT a MG vector
        do ieq=1,3
          call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
     .                           ,dv_cnv(:,:,:,ieq),igrid,.false.)
        enddo

        call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ)
     .              ,icomp=IVX,is_vec=.true.,is_cnv=.true.,iorder=3)

c     Corrector step

        !Find corrections to rho,B,tmp right hand sides (Udv)
        call correctRhoRHS(dv_cnv,rhs(:,IRHO)   ,igrid)
        call correctTmpRHS(dv_cnv,rhs(:,ITMP)   ,igrid)
        call correctAflRHS(dv_cnv,rhs(:,IAX:IAZ),igrid)

        !Find better initial guesses
        xxx(:,IRHO)    = xxx(:,IRHO)    - dt*alpha*rhs(:,IRHO)
        xxx(:,ITMP)    = xxx(:,ITMP)    - dt*alpha*rhs(:,ITMP)
        xxx(:,IAX:IAZ) = xxx(:,IAX:IAZ) - dt*alpha*rhs(:,IAX:IAZ)

cc        !Find (D_m^-1)(Udv) --> yyy
cc        call diagonalScaling(1,ntotp,rho_diag,rhs(:,IRHO)
cc     .                      ,yyy(:,IRHO)   ,igrid)
cc        call diagonalScaling(3,ntotp,a_diag,rhs(:,IAX:IAZ)
cc     .                      ,yyy(:,IAX:IAZ),igrid)
cc        call diagonalScaling(1,ntotp,tmp_diag,rhs(:,ITMP)
cc     .                      ,yyy(:,ITMP)   ,igrid)
cc
cc        !Find Deltax = xxx - yyy
cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc              dvol = volume(i,j,k,igx,igy,igz)
cc              xxx(ii,IRHO)    = 
cc     .                     xxx(ii,IRHO)    -dvol*alpha*yyy(ii,IRHO)
cc              xxx(ii,ITMP)    =
cc     .                     xxx(ii,ITMP)    -dvol*alpha*yyy(ii,ITMP)
cc              xxx(ii,IAX:IAZ) =
cc     .                     xxx(ii,IAX:IAZ) -dvol*alpha*yyy(ii,IAX:IAZ)
cc            enddo
cc          enddo
cc        enddo

c     Alternate corrector step

cc        guess = 0
cc
cc        !Find corrected residuals
cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc              iii = neqd*(ii-1)
cc
cc              if (vol_wgt)
cc                dvol = volume(i,j,k,igx,igy,igz)
cc              else
cc                dvol = 1d0
cc              endif
cc
cc              yyy(ii,IRHO)   =y(iii+IRHO)       -alpha*dvol*rhs(ii,IRHO)
cc              yyy(ii,IAX:IAZ)=y(iii+IAX:iii+IAZ)-alpha*dvol*rhs(ii
cc     .                                                         ,IAX:IAZ)
cc              yyy(ii,ITMP)   =y(iii+ITMP)       -alpha*dvol*rhs(ii,ITMP)
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc        !Temperature
cc        call cSolver(1,ntotp,yyy(:,ITMP),xxx(:,ITMP),bcs(:,ITMP)
cc     .              ,igrid,iout,guess,tmp_mtvc,dg=tmp_diag,ncolors=2
cccc     .              ,line_relax=.false.)
cc     .              ,line_relax=line_relax)
cc
cc        !Density
cc        call cSolver(1,ntotp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
cc     .              ,igrid,iout,guess,rho_mtvc,dg=rho_diag,ncolors=2
cccc     .              ,line_relax=.false.)
cc     .              ,line_relax=line_relax)
cc
cc        !Magnetic field
cc        call cSolver(3,ntotp,yyy(:,IAX:IAZ),xxx(:,IAX:IAZ)
cc     .              ,bcs(:,IAX:IAZ),igrid,iout,guess,b_mtvc,dg=b_diag
cc     .              ,ncolors=4)

c     Postprocessing of velocity -> momentum

cc        call correctVfldRHS(ntotp,dv_cnv,xxx,rhs(:,IVX:IVZ),igrid)
cc
cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc              if (vol_wgt) then !Residuals are volume-weighed by default
cc                dvol = volume(i,j,k,igx,igy,igz)
cc              else
cc                dvol = 1d0
cc              endif
cc
cc              do ieq = IVX,IVZ
cc                iii = ieq + neqd*(ii-1)
cc                xxx(ii,ieq) = dt*(y(iii)/dvol - alpha*rhs(ii,ieq))
cc              enddo
cc            enddo
cc          enddo
cc        enddo

        if (.not.nc_eom_v) then
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
               ii  = i + nx*(j-1) + nx*ny*(k-1)
               xxx(ii,IVX)=rho(i,j,k)*xxx(ii,IVX)+xxx(ii,IRHO)*vx(i,j,k)
               xxx(ii,IVY)=rho(i,j,k)*xxx(ii,IVY)+xxx(ii,IRHO)*vy(i,j,k)
               xxx(ii,IVZ)=rho(i,j,k)*xxx(ii,IVZ)+xxx(ii,IRHO)*vz(i,j,k)
              enddo
            enddo
          enddo
        endif

c     Gather solution in "x"

        call gather(x,xxx)

c     Diagnostics

        if (jit == 5 .and. debug) call PCdiagnostics

c     Deallocate variables

        deallocate(dv_cnv)

c *******************************************************************
c     Semi-implicit preconditioner version II
c *******************************************************************

      elseif (precon == 's2') then

c     Set up vectors

        x   = 0d0
        xxx = 0d0

        guess = 0   !Equations are in residual form

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3))

        allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3))

#if !defined(petsc)
        do si_it=1,precpass

c       Form residual vector rr=y-Ax

          if (si_it > 1) then
            call matrixFreeMatVec(ntot,x,rr)
            rr = y - rr
          else
            rr = y
          endif

#else
          si_it = mod(si_it,precpass)+1

          rr = y
#endif

          if (iout >= 0) then
            mag = sqrt(sum(rr*rr))
            if (si_it == 1) mag0 = mag
            mag = mag/mag0
            write (*,*) 'PC si iter',si_it-1,'  ;   Rel. res.',mag
            write (*,*)
            !Individual residual eqs
            call scatter(yyy,rr)
            do ieq=1,neqd
              mag = sqrt(sum(yyy(:,ieq)**2))
              write (*,*) 'Res. eqn',ieq,':',mag
            enddo
          endif

c       Map residual to individual components (rr --> yyy)

          call scatter(yyy,rr)

c       Form (D_m^-1)rx vector ---> xxx

cc          if (si_it == 1) then
cc            write (*,*) 'here'
cc            !Temperature
cc            call cSolver(1,ntotp,yyy(:,ITMP),xxx(:,ITMP),bcs(:,ITMP)
cc     .                  ,igrid,iout,guess,tmp_mtvc,dg=tmp_diag,ncolors=2
cc     .                  ,line_relax=line_relax)
cc
cc            !Density
cc            call cSolver(1,ntotp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
cc     .                  ,igrid,iout,guess,rho_mtvc,dg=rho_diag,ncolors=2
cc     .                  ,line_relax=line_relax)
cc
cc            !Magnetic field
cc            if (gm_smooth) then
cc              call cSolver(3,ntotp,yyy(:,IAX:IAZ),xxx(:,IAX:IAZ)
cc     .               ,bcs(:,IAX:IAZ),igrid,iout,guess,b_mtvc
cc     .               ,ncolors=4,gm_smth=gm_smooth,line_relax=line_relax)
cc            else
cc              call cSolver(3,ntotp,yyy(:,IAX:IAZ),xxx(:,IAX:IAZ)
cc     .                ,bcs(:,IAX:IAZ),igrid,iout,guess,b_mtvc,dg=b_diag
cc     .                ,ncolors=4,line_relax=line_relax)
cc            endif
cc          else
            call diagonalScaling(1,ntotp,rho_diag,yyy(:,IRHO)
     .                          ,xxx(:,IRHO)   ,igrid)
            call diagonalScaling(3,ntotp,a_diag  ,yyy(:,IAX:IAZ)
     .                          ,xxx(:,IAX:IAZ),igrid)
            call diagonalScaling(1,ntotp,tmp_diag,yyy(:,ITMP)
     .                          ,xxx(:,ITMP)   ,igrid)
cc          endif

c       SI step: Deltav --> xxx(:,IVX:IVZ)

          !Form SI rhs (-Gv-Ldx)
          call formSIrhs(ntotp,xxx,yyy(:,IVX:IVZ),rhs(:,IVX:IVZ),igrid)

          !Solve Schur-complement SI system
          if (gm_smooth) then
            call cSolver(3,ntotp,rhs(:,IVX:IVZ),xxx(:,IVX:IVZ)
     .                  ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc2
     .                  ,ncolors=4,gm_smth=gm_smooth
     .                  ,line_relax=line_relax)
cc     .                  ,line_relax=line_relax,cvrg_tst=debug)
          else
            call cSolver(3,ntotp,rhs(:,IVX:IVZ),xxx(:,IVX:IVZ)
     .               ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc2,dg=v_diag2
     .               ,ncolors=4,line_relax=line_relax)
cc     .                ,ncolors=4,line_relax=line_relax,cvrg_tst=debug)
          endif

          if (si_car) then
            do k = 1,nz
              do j = 1,ny
                do i = 1,nx
                  ii  = i + nx*(j-1) + nx*ny*(k-1)
                  car = xxx(ii,IVX:IVZ)
                  call transformVectorToCurvilinear
     .               (i,j,k,igrid,igrid,igrid
     .               ,car(1),car(2),car(3)
     .               ,.false.
     .               ,xxx(ii,IVX),xxx(ii,IVY),xxx(ii,IVZ))
                enddo
              enddo
            enddo
          endif

c       Store velocity solution in array format --> dv_cnv

          !xxx is NOT a MG vector
          do ieq=1,3
            call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
     .                             ,dv_cnv(:,:,:,ieq),igrid,.false.)
          enddo

          call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ)
     .                ,icomp=IVX,is_vec=.true.,is_cnv=.true.,iorder=3)

c       Correct rx to find Deltax (correction for rho, B, T)

          !Find corrections to rho,B,tmp right hand sides (Udv --> rhs)
          call correctRhoRHS(dv_cnv,rhs(:,IRHO)   ,igrid)
          call correctTmpRHS(dv_cnv,rhs(:,ITMP)   ,igrid)
          call correctAflRHS(dv_cnv,rhs(:,IAX:IAZ),igrid)

          !Find (D_m^-1)(Udv) --> yyy
          call diagonalScaling(1,ntotp,rho_diag,rhs(:,IRHO)
     .                        ,yyy(:,IRHO)   ,igrid)
          call diagonalScaling(3,ntotp,a_diag,rhs(:,IAX:IAZ)
     .                        ,yyy(:,IAX:IAZ),igrid)
          call diagonalScaling(1,ntotp,tmp_diag,rhs(:,ITMP)
     .                        ,yyy(:,ITMP)   ,igrid)

          !Find Deltax = xxx - yyy
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                dvol = volume(i,j,k,igx,igy,igz)
                xxx(ii,IRHO)    = 
     .                      xxx(ii,IRHO)   -dvol*alpha*yyy(ii,IRHO)
                xxx(ii,ITMP)    =
     .                      xxx(ii,ITMP)   -dvol*alpha*yyy(ii,ITMP)
                xxx(ii,IAX:IAZ) =
     .                      xxx(ii,IAX:IAZ)-dvol*alpha*yyy(ii,IAX:IAZ)
              enddo
            enddo
          enddo

c       Postprocessing of velocity -> momentum

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
               ii  = i + nx*(j-1) + nx*ny*(k-1)
               xxx(ii,IVX)=rho(i,j,k)*xxx(ii,IVX)+xxx(ii,IRHO)*vx(i,j,k)
               xxx(ii,IVY)=rho(i,j,k)*xxx(ii,IVY)+xxx(ii,IRHO)*vy(i,j,k)
               xxx(ii,IVZ)=rho(i,j,k)*xxx(ii,IVZ)+xxx(ii,IRHO)*vz(i,j,k)
              enddo
            enddo
          enddo

c      Clean up of rho,T,B

cc          if (si_it == precpass) then
cc
cc            guess = 1
cc
cc            !Find corrected residuals
cc            do k = 1,nz
cc              do j = 1,ny
cc                do i = 1,nx
cc                  ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc                  iii = neqd*(ii-1)
cc
cc                  dvol = volume(i,j,k,igx,igy,igz)
cc
cc                  yyy(ii,IRHO)   =y(iii+IRHO)-alpha*dvol*rhs(ii,IRHO)
cc                  yyy(ii,IAX:IAZ)=y(iii+IAX:iii+IAZ)
cc     .                                       -alpha*dvol*rhs(ii,IAX:IAZ)
cc                  yyy(ii,ITMP)   =y(iii+ITMP)-alpha*dvol*rhs(ii,ITMP)
cc
cc
cc                  if (.not.vol_wgt) yyy(ii,:) = yyy(ii,:)/dvol
cc
cc                enddo
cc              enddo
cc            enddo
cc
cc            !Temperature
cc            call cSolver(1,ntotp,yyy(:,ITMP),xxx(:,ITMP),bcs(:,ITMP)
cc     .                  ,igrid,iout,guess,tmp_mtvc,dg=tmp_diag,ncolors=2
cccc     .                  ,line_relax=.false.)
cc     .                  ,line_relax=line_relax)
cc
cc            !Density
cc            call cSolver(1,ntotp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
cc     .                  ,igrid,iout,guess,rho_mtvc,dg=rho_diag,ncolors=2
cccc     .                  ,line_relax=.false.)
cc     .                  ,line_relax=line_relax)
cc
cc            !Magnetic field
cc            call cSolver(3,ntotp,yyy(:,IAX:IAZ),xxx(:,IAX:IAZ)
cc     .                 ,bcs(:,IAX:IAZ),igrid,iout,guess,b_mtvc,dg=b_diag
cc     .                 ,ncolors=4)
cc
cc          endif

c     Postprocessing of magnetic field: divergence cleaning

cc        if (solenoidal.and.(.not.gm_smooth)) then
cc          call findDivfreeRHS(ntotp,dv_cnv,xxx(:,IAX:IAZ),rhs(:,IAX:IAZ)
cc     .                       ,igrid)
cc
cc        !Divergence-free correction for magnetic field
cc          do k = 1,nz
cc            do j = 1,ny
cc              do i = 1,nx
cc                ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc                if (vol_wgt) then !Residuals are volume-weighed by default
cc                  dvol = volume(i,j,k,igx,igy,igz)
cc                else
cc                  dvol = 1d0
cc                endif
cc
cc                do ieq = IAX,IAZ
cc                  iii = ieq + neqd*(ii-1)
cc                  xxx(ii,ieq) = dt*(y(iii)/dvol - alpha*rhs(ii,ieq))
cc                enddo
cc              enddo
cc            enddo
cc          enddo
cc
cc        endif

c       Map solution xxx to x (jacobi iterate)

#if !defined(petsc)
          omega = 1d0

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                do ieq=1,neqd
                  iii = ieq + neqd*(ii-1)
                  x(iii) = x(iii) + omega*xxx(ii,ieq)
                enddo
              enddo
            enddo
          enddo

        enddo

        if (iout >= 0) then
          call matrixFreeMatVec(ntot,x,rr)
          rr = y - rr
          mag = sqrt(sum(rr*rr))/mag0
          write (*,*) 'PC si iter',si_it-1,'  ;   Rel. res.',mag
          write (*,*)
          !Individual residual eqs
          call scatter(yyy,rr)
          do ieq=1,neqd
            mag = sqrt(sum(yyy(:,ieq)**2))
            write (*,*) 'Res. eqn',ieq,':',mag
          enddo
        endif
#else

        call gather(x,xxx)

#endif

c     Diagnostics

        if (jit == 1 .and. debug) call PCdiagnostics

c     Deallocate variables

        deallocate(dv_cnv,db_cnv,dj_cov)

      endif

c End program

      contains

c     PCdiagnostics
c     ##############################################################
      subroutine PCdiagnostics

c     --------------------------------------------------------------
c     Performs diagnostics on PC
c     --------------------------------------------------------------

        implicit none

        real(8),allocatable,dimension(:,:,:,:) :: da_cov

c     Begin program

#if !defined(petsc)

        allocate(da_cov(0:nx+1,0:ny+1,0:nz+1,3))

        !Solution plot
        call MGplot(neqd,x,igrid,0,'sol.bin')

        !Predictor RHS plot
        call MGplot(neqd,y,igrid,0,'rhs.bin')

        !Residual plot
        call matrixFreeMatVec(ntot,x,rr)
        mag0 = sqrt(sum(y*y))
        y = y - rr
        mag = sqrt(sum(y*y))
        write (*,*) 'abs res=',mag/ntot,', rel res=',mag/mag0
        call MGplot(neqd,y,igrid,0,'res.bin')

        !Individual residual eqs
        call scatter(yyy,y)
        do ieq=1,neqd
          mag = sqrt(sum(yyy(:,ieq)**2))/nx/ny/nz
          write (*,*) 'Res. eqn',ieq,':',mag
        enddo

        !Corrector rhs plot
        if (si_car) then
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                car = rhs(ii,IVX:IVZ)
                call transformVectorToCurvilinear
     .             (i,j,k,igrid,igrid,igrid
     .             ,car(1)     ,car(2)     ,car(3)
     .             ,.false.
     .             ,rhs(ii,IVX),rhs(ii,IVY),rhs(ii,IVZ))
              enddo
            enddo
          enddo
        endif

        call MGplot(1,rhs(:,IRHO),igrid,0,'cor_rhs.bin')
        call MGplot(1,rhs(:,IVX) ,igrid,1,'cor_rhs.bin')
        call MGplot(1,rhs(:,IVY) ,igrid,1,'cor_rhs.bin')
        call MGplot(1,rhs(:,IVZ) ,igrid,1,'cor_rhs.bin')
        call MGplot(1,rhs(:,IAX) ,igrid,1,'cor_rhs.bin')
        call MGplot(1,rhs(:,IAY) ,igrid,1,'cor_rhs.bin')
        call MGplot(1,rhs(:,IAZ) ,igrid,1,'cor_rhs.bin')
        call MGplot(1,rhs(:,ITMP),igrid,1,'cor_rhs.bin')

        !Vector components plot
        !Find dA* (w/ BCs)
        do ieq=1,3
          call mapMGVectorToArray(0,1,xxx(:,IAX+ieq-1),nx,ny,nz
cc          call mapMGVectorToArray(0,1,rhs(:,IAX+ieq-1),nx,ny,nz
     .                       ,da_cov(:,:,:,ieq),igrid,.false.)
        enddo

        call setMGBC(0,3,nx,ny,nz,igrid,da_cov,bcs(:,IAX:IAZ)
     .              ,icomp=IAX,is_vec=.true.,is_cnv=.false.,iorder=2)

        !Find dv (w/ BCs)
        do ieq=1,3
          call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
     .                         ,dv_cnv(:,:,:,ieq),igrid,.false.)
        enddo

        call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ)
     .            ,icomp=IVX,is_vec=.true.,is_cnv=.true.,iorder=2)

        open(unit=110,file='debug4.bin',form='unformatted'
     .     ,status='replace')
        call contour(dv_cnv(0:nx+1,1:ny,1,1),nx+2,ny,0d0
     .              ,xmax,0d0,ymax,0,110)
        call contour(dv_cnv(0:nx+1,1:ny,1,2),nx+2,ny,0d0
     .              ,xmax,0d0,ymax,1,110)
        call contour(dv_cnv(0:nx+1,1:ny,1,3),nx+2,ny,0d0
     .              ,xmax,0d0,ymax,1,110)
        call contour(da_cov(0:nx+1,1:ny,1,1),nx+2,ny,0d0
     .              ,xmax,0d0,ymax,1,110)
        call contour(da_cov(0:nx+1,1:ny,1,2),nx+2,ny,0d0
     .              ,xmax,0d0,ymax,1,110)
        call contour(da_cov(0:nx+1,1:ny,1,3),nx+2,ny,0d0
     .              ,xmax,0d0,ymax,1,110)
        close(110)

        deallocate(da_cov)
        stop
#endif

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

      use mg_internal

      implicit none

c Call variables

      integer(4) :: ntotp,igrid
      real(8)    :: xxx(ntotp,neqd),eom_rhs(ntotp,3)

c Local variables

      integer(4) :: ii,isig,ivar,ip,im
      real(8)    :: dvol,cov(3),cnv(3),idx,idy,idz
     .             ,da_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,da_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dj_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,db_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dpres (0:nx+1,0:ny+1,0:nz+1)
     .             ,drho  (0:nx+1,0:ny+1,0:nz+1)
     .             ,dtmp  (0:nx+1,0:ny+1,0:nz+1)

c Begin program

      isig = grid_params%istartp(igrid)

c Find dA* (w/ BCs)

      da_cov = 0d0
      !xxx is NOT a MG vector
      do ieq=1,3
        call mapMGVectorToArray(0,1,xxx(:,IAX+ieq-1),nx,ny,nz
     .                         ,da_cov(:,:,:,ieq),igrid,.false.)
      enddo

      !Find covariant components of db with BCs
      call setBC(IAX,3,nx,ny,nz,da_cnv,da_cov,vzeros,bcs(:,IAX:IAZ)
     .          ,igrid,igrid,igrid,is_cnv=.false.)

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
     .          ,igrid,igrid,igrid)

c Find dj* from dB*

      !Find contravariant current (without BCs)
      dj_cnv = 0d0
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            do ivar=1,3
              dj_cnv(i,j,k,ivar)=curl2(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                                ,db_cov(:,:,:,1)
     .                                ,db_cov(:,:,:,2)
     .                                ,db_cov(:,:,:,3)
     .                                ,ivar)
              call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .         ,dj_cov(i,j,k,1),dj_cov(i,j,k,2),dj_cov(i,j,k,3)
     .         ,dj_cnv(i,j,k,1),dj_cnv(i,j,k,2),dj_cnv(i,j,k,3),.false.)
            enddo
          enddo
        enddo
      enddo

      !Find covariant components of dj with BCs
cc      call setBC(IJX,3,nx,ny,nz,dj_cnv,dj_cov,vzeros,bcs(:,IJX:IJZ)
cc     .          ,igrid,igrid,igrid)

c Find dpres array with BCs

cc      dpres = 0d0
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            ii = i + nx*(j-1) + nx*ny*(k-1)
cc            dpres(i,j,k)= 2.*(rho(i,j,k)*xxx(ii,ITMP)
cc     .                       +tmp(i,j,k)*xxx(ii,IRHO))
cc          enddo
cc        enddo
cc      enddo
cc
cc      !Find BCs
cc      call setBC(IRHO,nx,ny,nz,dpres,zeros,bcs(:,IRHO)
cc     .          ,igrid,igrid,igrid)

      drho = 0d0
      dtmp = 0d0
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii = i + nx*(j-1) + nx*ny*(k-1)
            drho(i,j,k) = xxx(ii,IRHO)
            dtmp(i,j,k) = xxx(ii,ITMP)
          enddo
        enddo
      enddo

      !Find BCs
      call setBC(IRHO,nx,ny,nz,drho,zeros,bcs(:,IRHO)
     .          ,igrid,igrid,igrid)

      call setBC(ITMP,nx,ny,nz,dtmp,zeros,bcs(:,ITMP)
     .          ,igrid,igrid,igrid)

c Find rhs_v'

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jac = gmetric%grid(igrid)%jac(i,j,k)

            if (vol_wgt) then
              dvol = volume(i,j,k,igrid,igrid,igrid)
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
cc     .             -2*(rho(i,j,k)*(dtmp(i+1,j,k)-dtmp(i-1,j,k))*idx
cc     .                +tmp(i,j,k)*(drho(i+1,j,k)-drho(i-1,j,k))*idx)
cc
cc            cov(2) = jz(i,j,k)*db_cnv(i,j,k,1)/jac
cc     .             - jx(i,j,k)*db_cnv(i,j,k,3)/jac
cc     .             - bz(i,j,k)*dj_cnv(i,j,k,1)/jac
cc     .             + bx(i,j,k)*dj_cnv(i,j,k,3)/jac
cccc     .             -(dpres(i,j+1,k)-dpres(i,j-1,k))*idy
cc     .             -2*(rho(i,j,k)*(dtmp(i,j+1,k)-dtmp(i,j-1,k))*idy
cc     .                +tmp(i,j,k)*(drho(i,j+1,k)-drho(i,j-1,k))*idy)
cc
cc            cov(3) = jx(i,j,k)*db_cnv(i,j,k,2)/jac
cc     .             - jy(i,j,k)*db_cnv(i,j,k,1)/jac
cc     .             - bx(i,j,k)*dj_cnv(i,j,k,2)/jac
cc     .             + by(i,j,k)*dj_cnv(i,j,k,1)/jac
cccc     .             -(dpres(i,j,k+1)-dpres(i,j,k-1))*idz
cc     .             -2*(rho(i,j,k)*(dtmp(i,j,k+1)-dtmp(i,j,k-1))*idz
cc     .                +tmp(i,j,k)*(drho(i,j,k+1)-drho(i,j,k-1))*idz)

cc            cov(1) = -(dpres(ip,j,k)-dpres(im,j,k))*idx
            cov(1) = -2*(rho(i,j,k)*(dtmp(i+1,j,k)-dtmp(i-1,j,k))*idx
     .                  +tmp(i,j,k)*(drho(i+1,j,k)-drho(i-1,j,k))*idx)

cc            cov(2) = -(dpres(i,j+1,k)-dpres(i,j-1,k))*idy
            cov(2) = -2*(rho(i,j,k)*(dtmp(i,j+1,k)-dtmp(i,j-1,k))*idy
     .                  +tmp(i,j,k)*(drho(i,j+1,k)-drho(i,j-1,k))*idy)

cc            cov(3) = -(dpres(i,j,k+1)-dpres(i,j,k-1))*idz
            cov(3) = -2*(rho(i,j,k)*(dtmp(i,j,k+1)-dtmp(i,j,k-1))*idz
     .                  +tmp(i,j,k)*(drho(i,j,k+1)-drho(i,j,k-1))*idz)

            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3),.true.)

            !Add Lorentz force correction
cc            vec1 => gb0%grid(igrid)%array
cc            vec2 => db_cnv
cc            cnv = cnv
cc     .          + div_tensor(i,j,k,nx,ny,nz,igrid,igrid,igrid,.false.
cc     .                      ,lf_x,lf_y,lf_z,vol=.false.)
cc            nullify(vec1,vec2)

            cnv(1) = cnv(1)
     .             + jy_cov(i,j,k)*db_cov(i,j,k,3)
     .             - jz_cov(i,j,k)*db_cov(i,j,k,2)
     .             - by_cov(i,j,k)*dj_cov(i,j,k,3)
     .             + bz_cov(i,j,k)*dj_cov(i,j,k,2)

            cnv(2) = cnv(2)
     .             + jz_cov(i,j,k)*db_cov(i,j,k,1)
     .             - jx_cov(i,j,k)*db_cov(i,j,k,3)
     .             - bz_cov(i,j,k)*dj_cov(i,j,k,1)
     .             + bx_cov(i,j,k)*dj_cov(i,j,k,3)

            cnv(3) = cnv(3)
     .             + jx_cov(i,j,k)*db_cov(i,j,k,2)
     .             - jy_cov(i,j,k)*db_cov(i,j,k,1)
     .             - bx_cov(i,j,k)*dj_cov(i,j,k,2)
     .             + by_cov(i,j,k)*dj_cov(i,j,k,1)

            !Correct rhs_v
            eom_rhs(ii,:) = alpha*dvol*cnv(:)

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

      use mg_internal

      implicit none

c Call variables

      integer(4) :: ntotp,igrid
      real(8)    :: xxx(ntotp,neqd),yyy(ntotp,3),rhs_si(ntotp,3)

c Local variables

      integer(4) :: ii,iig,isig,ivar,ip,im
      real(8)    :: dvol,cov(3),cnv(3),idx,idy,idz,eom_rhs(ntotp,3)

c Begin program

      isig = grid_params%istartp(igrid)

c Find rhs_v

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            iig = ii + isig - 1
            if (vol_wgt) then
              dvol = volume(i,j,k,igrid,igrid,igrid)
            else
              dvol = 1d0
            endif
            rhs_si(ii,:) = yyy(ii,:)
     .                   - dvol*xxx(ii,IRHO)*mgadvdiffV0(iig,:)
          enddo
        enddo
      enddo

c Find linearized EOM rhs [jxB-grad(p)] correction

      call findEOMrhs(ntotp,xxx,eom_rhs,igrid)

c Add rhs_v + EOM rhs

      rhs_si = rhs_si + eom_rhs

c Transform to cartesian if applicable

      if (si_car) then  

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              cnv(:) = rhs_si(ii,:)

              call transformVectorToCartesian
     .               (i,j,k,igrid,igrid,igrid
     .               ,cnv(1),cnv(2),cnv(3),.false.
     .               ,rhs_si(ii,1),rhs_si(ii,2),rhs_si(ii,3))


            enddo
          enddo
        enddo

      endif

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

      use mg_internal

      implicit none

c Call variables

      integer(4) :: ntotp,igrid
      real(8)    :: xxx(ntotp,neqd),crhs(ntotp,3)
     .             ,dv(0:nx+1,0:ny+1,0:nz+1,3)

c Local variables

      integer(4) :: ii,iig,isig,ivar,ip,im
      real(8)    :: dvol,cov(3),cnv(3),idx,idy,idz,eom_rhs(ntotp,3)
     .             ,nabla_vv0(3,3),mul,veclap(3),advec(3)

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
              dvol = volume(i,j,k,igrid,igrid,igrid)
            else
              dvol = 1d0
            endif

            !EOM residual part + EM part
            cnv = -eom_rhs(ii,:)/dvol/alpha

            !First advection contribution
            nabla_v   = fnabla_v(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,dv(:,:,:,1),dv(:,:,:,2),dv(:,:,:,3),0)

            nabla_vv0 = fnabla_v(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,vx,vy,vz,0)

            mul = vis(i,j,k,nx,ny,nz,igrid,igrid,igrid)

            veclap = veclaplacian(i,j,k,nx,ny,nz,igrid,igrid,igrid,dv
     .                           ,alt_eom,vol=.false.)

            advec =( vx(i,j,k)  *nabla_v  (1,:)
     .              +vy(i,j,k)  *nabla_v  (2,:)
     .              +vz(i,j,k)  *nabla_v  (3,:)
     .              +dv(i,j,k,1)*nabla_vv0(1,:)
     .              +dv(i,j,k,2)*nabla_vv0(2,:)
     .              +dv(i,j,k,3)*nabla_vv0(3,:))/jac

            cnv = cnv + rho(i,j,k)*(advec - mul*veclap)

            !Second advection contribution
            veclap = veclaplacian(i,j,k,nx,ny,nz,igrid,igrid,igrid,vcnv
     .                           ,alt_eom,vol=.false.)

            advec =( vx(i,j,k)*nabla_vv0(1,:)
     .              +vy(i,j,k)*nabla_vv0(2,:)
     .              +vz(i,j,k)*nabla_vv0(3,:))/jac

            cnv = cnv + xxx(ii,IRHO)*(advec - mul*veclap)

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
c     This subroutine finds the corrected rhs for the density solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use mg_internal

      implicit none

c Call variables

      integer(4) :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz)

c Local variables

      integer(4) :: ii,iig,ip,im,jp,jm,kp,km
      real(8)    :: dh(3),mag,tvol,vol
      real(8)    :: rhodv_cnv(0:nx+1,0:ny+1,0:nz+1,3)

c Begin program

c Find rho0*dv

      rhodv_cnv(:,:,:,1) = dv_cnv(:,:,:,1)*rho
      rhodv_cnv(:,:,:,2) = dv_cnv(:,:,:,2)*rho
      rhodv_cnv(:,:,:,3) = dv_cnv(:,:,:,3)*rho

c Evaluate rhs correction: div(rh0*dv)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            crhs(ii) = div(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                    ,rhodv_cnv(:,:,:,1)
     .                    ,rhodv_cnv(:,:,:,2)
     .                    ,rhodv_cnv(:,:,:,3))

cc            ip = i+1
cc            im = i-1
cc            jp = j+1
cc            jm = j-1
cc            kp = k+1
cc            km = k-1
cc
cc            jac = gmetric%grid(igrid)%jac(i,j,k)
cc
cc            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
cc
cc            dh(1) = 2*dxh(ig)
cc            dh(2) = 2*dyh(jg)
cc            dh(3) = 2*dzh(kg)
cc
cc            crhs(ii) = rho(i,j,k)
cc     .                *div(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                    ,dv_cnv(:,:,:,1)
cc     .                    ,dv_cnv(:,:,:,2)
cc     .                    ,dv_cnv(:,:,:,3))
cc     .              +dv_cnv(i,j,k,1)*(rho(ip,j,k)-rho(im,j,k))/dh(1)/jac
cc     .              +dv_cnv(i,j,k,2)*(rho(i,jp,k)-rho(i,jm,k))/dh(2)/jac
cc     .              +dv_cnv(i,j,k,3)*(rho(i,j,kp)-rho(i,j,km))/dh(3)/jac
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

      use mg_internal

      implicit none

c Call variables

      integer(4) :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz)

c Local variables

      integer(4) :: ii

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

c correctAflRHS
c ###################################################################
      subroutine correctAflRHS(dv_cnv,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the vector potential
c     solve: div(dv B0 - B0 dv)
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use mg_internal

      implicit none

c Call variables

      integer(4) :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz,3)

c Local variables

      integer(4) :: ii

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
          enddo
        enddo
      enddo

c End program

      end subroutine correctAflRHS

c diagonalScaling
c #########################################################################
      subroutine diagonalScaling(neq,ntotp,idiag,y,x,igrid)
c--------------------------------------------------------------------
c     Performs x=idiag*y, where idiag contains the inverse of the
c     diagonal.
c--------------------------------------------------------------------

      use mg_internal

      implicit none

c Call variables

      integer(4) :: neq,ntotp,igrid
      real(8)    :: idiag(neq,*),y(ntotp,neq),x(ntotp,neq)

c Local variables

      integer(4) :: ii,iii,iig,isig
      logical    :: fpointers

c Begin program

      call allocPointers(neq,fpointers)

      isig = istart(igrid)

      do ii = 1,ntotp

        iii = neq*(ii-1)
        iig = iii + isig - 1

        x(ii,:) = matmul(idiag(:,iig+1:iig+neq),y(ii,:))

      enddo

      call deallocPointers(fpointers)

c End program

      end subroutine diagonalScaling
