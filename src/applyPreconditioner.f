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

      use newton_gmres

      use mg_internal

      implicit none

c Call variables

      integer(4) :: ntot,iout
      real(8)    :: x(ntot),y(ntot)

c Local variables

      real(8)    :: xxx(ntot/neqd,neqd),yyy(ntot/neqd,neqd)
     .             ,rhs(ntot/neqd,neqd),rr(ntot),dxk(ntot)
     .             ,dvol,car(3),cnv(3)

      real(8),allocatable,dimension(:,:,:,:) :: dv_cnv,db_cnv,dj_cov

      integer(4) :: ii,iii,igrid,ntotp,guess,si_it

c Debug

      real(8)    :: mag,mag2,mag0
      real(8),allocatable,dimension(:,:,:,:) :: dbg

c Externals

      external   tmp_mtvc,rho_mtvc,b_mtvc,v_mtvc,test_mtvc

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
                x(iii) = y(iii)*dt/volume(i,j,k,igx,igy,igz)
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

        yyy = y     !Overloaded operation

        guess = 0

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3))

        allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3))

c     Predictor step

        !Temperature
        call cSolver(1,ntotp,yyy(:,ITMP),xxx(:,ITMP),bcs(:,ITMP)
     .              ,igrid,iout,guess,tmp_mtvc,dg=tmp_diag,ncolors=2
     .              ,line_relax=line_relax)

        !Density
        call cSolver(1,ntotp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
     .              ,igrid,iout,guess,rho_mtvc,dg=rho_diag,ncolors=2
     .              ,line_relax=line_relax)

        call cSolver(3,ntotp,yyy(:,IBX:IBZ),xxx(:,IBX:IBZ)
     .              ,bcs(:,IBX:IBZ),igrid,iout,guess,b_mtvc,dg=b_diag
cc     .              ,ncolors=4)
     .              ,ncolors=4,line_relax=line_relax)

cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc              xxx(ii,IBX:IBZ) = xxx(ii,IBX:IBZ) + dt*yyy(ii,IBX:IBZ)
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
     .                ,ncolors=4,line_relax=line_relax)
cc     .                ,ncolors=4,line_relax=line_relax,cvrg_tst=debug)
        endif

cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc              xxx(ii,IVX:IVZ) = xxx(ii,IVX:IVZ) + dt*rhs(ii,IVX:IVZ)
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
     .              ,icomp=IVX,is_vec=.true.,is_cnv=.true.,iorder=2)

c     Corrector step

        !Find corrections to rho,B,tmp right hand sides (Udv)
        call correctRhoRHS(dv_cnv,rhs(:,IRHO)   ,igrid)
        call correctTmpRHS(dv_cnv,rhs(:,ITMP)   ,igrid)
        call correctBflRHS(dv_cnv,rhs(:,IBX:IBZ),igrid)

        !Find better initial guesses
        xxx(:,IRHO)    = xxx(:,IRHO)    - dt*alpha*rhs(:,IRHO)
        xxx(:,ITMP)    = xxx(:,ITMP)    - dt*alpha*rhs(:,ITMP)
        xxx(:,IBX:IBZ) = xxx(:,IBX:IBZ) - dt*alpha*rhs(:,IBX:IBZ)

c     Postprocessing of magnetic field: divergence cleaning

c$$$        if (.not.zip_eom_b) then
c$$$          call findDivfreeRHS(dv_cnv,xxx(:,IBX:IBZ),rhs(:,IBX:IBZ)
c$$$     .                       ,igrid)
c$$$
c$$$        !Divergence-free correction for magnetic field
c$$$          do k = 1,nz
c$$$            do j = 1,ny
c$$$              do i = 1,nx
c$$$                ii  = i + nx*(j-1) + nx*ny*(k-1)
c$$$                dvol = 1d0/volume(i,j,k,igx,igy,igz)
c$$$                do ieq = IBX,IBZ
c$$$                  iii = ieq + neqd*(ii-1)
c$$$                  xxx(ii,ieq) = dt*(y(iii)*dvol - alpha*rhs(ii,ieq))
c$$$                enddo
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$
c$$$        endif

ccc     Alternate corrector step
cc
cc        guess = 1
cc
cc        !Find corrected residuals
cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc              iii = neqd*(ii-1)
cc
cc              dvol = volume(i,j,k,igx,igy,igz)
cc
cc              yyy(ii,IRHO)   =y(iii+IRHO)       -alpha*dvol*rhs(ii,IRHO)
cc              yyy(ii,IBX:IBZ)=y(iii+IBX:iii+IBZ)-alpha*dvol*rhs(ii,IBX:IBZ)
cc              yyy(ii,ITMP)   =y(iii+ITMP)       -alpha*dvol*rhs(ii,ITMP)
cc
cc
cc              if (.not.vol_wgt) yyy(ii,:) = yyy(ii,:)/dvol
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
cc        call cSolver(3,ntotp,yyy(:,IBX:IBZ),xxx(:,IBX:IBZ)
cc     .              ,bcs(:,IBX:IBZ),igrid,iout,guess,b_mtvc,dg=b_diag
cc     .              ,ncolors=4)

c     Postprocessing of velocity -> momentum

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              xxx(ii,IVX) =rho(i,j,k)*xxx(ii,IVX)+xxx(ii,IRHO)*vx(i,j,k)
              xxx(ii,IVY) =rho(i,j,k)*xxx(ii,IVY)+xxx(ii,IRHO)*vy(i,j,k)
              xxx(ii,IVZ) =rho(i,j,k)*xxx(ii,IVZ)+xxx(ii,IRHO)*vz(i,j,k)
            enddo
          enddo
        enddo

c     Gather solution in "x"

        x = xxx         !Overloaded operation

c     Diagnostics

c diag ****
        if (jit == 1.and.debug) then

          vol_wgt = .false.

          !Solution plot
          call MGplot(neqd,x,igrid,0,'debug.bin')

          !Predictor RHS plot
          call MGplot(neqd,y,igrid,0,'debug3.bin')

          !Residual plot
          call matrixFreeMatVec(ntot,x,rr)
          mag0 = sqrt(sum(y*y))
          y = y - rr
          mag = sqrt(sum(y*y))
          write (*,*) 'abs res=',mag/ntot,', rel res=',mag/mag0

          call MGplot(neqd,y,igrid,0,'debug1.bin')
c$$$          yyy = y
c$$$          call MGplot(1,yyy(:,IRHO),igrid,0,'debug1.bin')
c$$$          call MGplot(1,yyy(:,IVX) ,igrid,1,'debug1.bin')
c$$$          call MGplot(1,yyy(:,IVY) ,igrid,1,'debug1.bin')
c$$$          call MGplot(1,yyy(:,IVZ) ,igrid,1,'debug1.bin')
c$$$          call MGplot(1,yyy(:,IBX) ,igrid,1,'debug1.bin')
c$$$          call MGplot(1,yyy(:,IBY) ,igrid,1,'debug1.bin')
c$$$          call MGplot(1,yyy(:,IBZ) ,igrid,1,'debug1.bin')
c$$$          call MGplot(1,yyy(:,ITMP),igrid,1,'debug1.bin')

          !Corrector rhs plot
          if (si_car) then
            do k = 1,nz
              do j = 1,ny
                do i = 1,nx
                  ii  = i + nx*(j-1) + nx*ny*(k-1)
                  car = rhs(ii,IVX:IVZ)
                  call transformVectorToCurvilinear
     .               (i,j,k,igrid,igrid,igrid
     .               ,car(1)     ,car(2)     ,car(3)
     .               ,.false.
     .               ,rhs(ii,IVX),rhs(ii,IVY),rhs(ii,IVZ))
                enddo
              enddo
            enddo
          endif

          call MGplot(1,rhs(:,IRHO),igrid,0,'debug2.bin')
          call MGplot(1,rhs(:,IVX) ,igrid,1,'debug2.bin')
          call MGplot(1,rhs(:,IVY) ,igrid,1,'debug2.bin')
          call MGplot(1,rhs(:,IVZ) ,igrid,1,'debug2.bin')
          call MGplot(1,rhs(:,IBX) ,igrid,1,'debug2.bin')
          call MGplot(1,rhs(:,IBY) ,igrid,1,'debug2.bin')
          call MGplot(1,rhs(:,IBZ) ,igrid,1,'debug2.bin')
          call MGplot(1,rhs(:,ITMP),igrid,1,'debug2.bin')

          !Vector components plot
          !Find dB* (w/ BCs)
          do ieq=1,3
            call mapMGVectorToArray(0,1,xxx(:,IBX+ieq-1),nx,ny,nz
     .                         ,db_cnv(:,:,:,ieq),igrid,.false.)
          enddo

          call setMGBC(0,3,nx,ny,nz,igrid,db_cnv,bcs(:,IBX:IBZ)
     .                ,icomp=IBX,is_vec=.true.,is_cnv=.true.,iorder=2)

          !Find dv (w/ BCs)
          do ieq=1,3
            call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
     .                           ,dv_cnv(:,:,:,ieq),igrid,.false.)
          enddo

          call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ)
     .              ,icomp=IVX,is_vec=.true.,is_cnv=.true.,iorder=2)

          open(unit=110,file='debug4.bin',form='unformatted'
     .       ,status='replace')
          call contour(dv_cnv(0:nx+1,1:ny,1,1),nx+2,ny,0d0
     .                ,xmax,0d0,ymax,0,110)
          call contour(dv_cnv(0:nx+1,1:ny,1,2),nx+2,ny,0d0
     .                ,xmax,0d0,ymax,1,110)
          call contour(dv_cnv(0:nx+1,1:ny,1,3),nx+2,ny,0d0
     .                ,xmax,0d0,ymax,1,110)
          call contour(db_cnv(0:nx+1,1:ny,1,1),nx+2,ny,0d0
     .                ,xmax,0d0,ymax,1,110)
          call contour(db_cnv(0:nx+1,1:ny,1,2),nx+2,ny,0d0
     .                ,xmax,0d0,ymax,1,110)
          call contour(db_cnv(0:nx+1,1:ny,1,3),nx+2,ny,0d0
     .                ,xmax,0d0,ymax,1,110)
          close(110)

         !diag B-field divergence
          allocate(dbg(0:nx+1,0:ny+1,0:nz+1,1))
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

          call setMGBC(0,1,nx,ny,nz,igrid,dbg,bcond,icomp=IRHO,iorder=2)

          do k=1,nz
            dbg(1,:,k) = dbg(0,:,k)
          enddo

          mag = sqrt(sum(dbg(1:nx,1:ny,1:nz)**2))
          write (*,*) 'L-2 norm of div(dB)',mag
          deallocate(dbg)

cc        write (*,*) 'plot div(b)?'
cc        read (*,'(a)') plot
cc        if (plot == 'y') then
cc          open(unit=110,file='debug.bin',form='unformatted'
cc     .       ,status='replace')
cc          call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,0,110)
cc          close(110)
cc          stop
cc        endif

          stop
        endif
c diag ******

c     Deallocate variables

        deallocate(dv_cnv,db_cnv,dj_cov)

c *******************************************************************
c     Semi-implicit preconditioner version II
c *******************************************************************

      elseif (precon == 's2') then

c     Set up vectors

        x   = 0d0
        xxx = 0d0
        rr  = y

        guess = 0   !Equations are in residual form

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3))

        allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3))

        do si_it=1,precpass

c       Form residual vector rr=y-Ax

          if (si_it > 1) then
            call matrixFreeMatVec(ntot,x,rr)
            rr = y - rr

            mag = sqrt(sum(rr*rr))
            mag0=sqrt(sum(y*y))
            mag = mag/mag0
            write (*,*) si_it-1,mag
          endif

c       Map residual to individual components (rr --> yyy)

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                dvol = volume(i,j,k,igx,igy,igz)
                do ieq=1,neqd
                  iii = ieq + neqd*(ii-1)
                  yyy(ii,ieq) = rr(iii)*dvol !Volume weighing of residuals
                enddo
              enddo
            enddo
          enddo

c       Form (D_m^-1)rx vector ---> xxx

          call diagonalScaling(1,ntotp,rho_diag,yyy(:,IRHO)
     .                        ,xxx(:,IRHO)   ,igrid)
          call diagonalScaling(3,ntotp,b_diag  ,yyy(:,IBX:IBZ)
     .                        ,xxx(:,IBX:IBZ),igrid)
          call diagonalScaling(1,ntotp,tmp_diag,yyy(:,ITMP)
     .                        ,xxx(:,ITMP)   ,igrid)
          
c       SI step: Deltav --> xxx(:,IVX:IVZ)

          !Form SI rhs (rv-L.(Dm^-1)rx) --> rhs(:,IVX:IVZ) (also form arrays db_cnv,dj_cov)
          call formSIrhs(ntotp,xxx,yyy(:,IVX:IVZ),rhs(:,IVX:IVZ)
     .                  ,db_cnv,dj_cov,igrid)

          !Solve Schur-complement SI system ---> xxx(:,IVX:IVZ)
          call cSolver(3,ntotp,rhs(:,IVX:IVZ),xxx(:,IVX:IVZ)
     .           ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,dg=v_diag
     .           ,ncolors=4,line_relax=line_relax)

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

          !Set grid=1 because xxx is NOT a MG vector
          do ieq=1,3
            call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
     .                             ,dv_cnv(:,:,:,ieq),igrid,.false.)
          enddo

          call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ)
     .                ,icomp=IVX,is_vec=.true.,is_cnv=.true.,iorder=2)

c       Correct rx to find Deltax (correction for rho, B, T)

          !Find corrections to rho,B,tmp right hand sides (Udv --> rhs)
          call correctRhoRHS(dv_cnv,rhs(:,IRHO)   ,igrid)
          call correctTmpRHS(dv_cnv,rhs(:,ITMP)   ,igrid)
          call correctBflRHS(dv_cnv,rhs(:,IBX:IBZ),igrid)

          !Find (D_m^-1)(Udv) --> yyy
          call diagonalScaling(1,ntotp,rho_diag,rhs(:,IRHO)
     .                        ,yyy(:,IRHO)   ,igrid)
          call diagonalScaling(3,ntotp,b_diag,rhs(:,IBX:IBZ)
     .                        ,yyy(:,IBX:IBZ),igrid)
          call diagonalScaling(1,ntotp,tmp_diag,rhs(:,ITMP)
     .                        ,yyy(:,ITMP)   ,igrid)

c diag ******
cc          call MGplot(1,rhs(:,IBX),igrid,0,'debug.bin')
cc          call MGplot(1,rhs(:,IBY),igrid,1,'debug.bin')
cc          call MGplot(1,rhs(:,IBZ),igrid,1,'debug.bin')
cc          call MGplot(1,yyy(:,IBX),igrid,1,'debug.bin')
cc          call MGplot(1,yyy(:,IBY),igrid,1,'debug.bin')
cc          call MGplot(1,yyy(:,IBZ),igrid,1,'debug.bin')
cc          stop
c diag ******

          !Find Deltax = xxx - yyy
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                dvol = volume(i,j,k,igx,igy,igz)
                xxx(ii,IRHO)    = 
     .                       xxx(ii,IRHO)    -dvol*alpha*yyy(ii,IRHO)
                xxx(ii,ITMP)    =
     .                       xxx(ii,ITMP)    -dvol*alpha*yyy(ii,ITMP)
                xxx(ii,IBX:IBZ) =
     .                       xxx(ii,IBX:IBZ) -dvol*alpha*yyy(ii,IBX:IBZ)
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

c     Postprocessing of magnetic field: divergence cleaning

c$$$          call findDivfreeRHS(dv_cnv,db_cnv,dj_cov,rhs(:,IBX:IBZ),igrid)
c$$$
c$$$          !Divergence-free correction for magnetic field (needs changes in correctBflRHS)
c$$$          do k = 1,nz
c$$$            do j = 1,ny
c$$$              do i = 1,nx
c$$$                ii  = i + nx*(j-1) + nx*ny*(k-1)
c$$$                do ieq = IBX,IBZ
c$$$                  iii = ieq + neqd*(ii-1)
c$$$                  xxx(ii,ieq) = dt*(y(iii) - alpha*rhs(ii,ieq))
c$$$                enddo
c$$$              enddo
c$$$            enddo
c$$$          enddo

c       Map solution xxx to x (jacobi iterate)

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                do ieq=1,neqd
                  iii = ieq + neqd*(ii-1)
                  x(iii) = x(iii) + xxx(ii,ieq)
                enddo
              enddo
            enddo
          enddo

        enddo

c diag ******
        call matrixFreeMatVec(ntot,x,rr)

        rr = y - rr

        mag = sqrt(sum(rr*rr))
        mag = mag/mag0
        write (*,*) si_it-1,mag
c diag ******

c plot ************
cc        if (mag > 1d0) then
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                do ieq=1,neqd
                  iii = ieq + neqd*(ii-1)
crhs                  xxx(ii,ieq) = y(iii)
                  xxx(ii,ieq) = x(iii)
                enddo
              enddo
            enddo
          enddo
          
          call MGplot(1,xxx(:,IRHO),igrid,0,'debug.bin')
          call MGplot(1,xxx(:,IVX) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,IVY) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,IVZ) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,IBX) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,IBY) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,IBZ) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,ITMP),igrid,1,'debug.bin')

          stop
cc        endif
c plot *************

ccc     Store velocity solution in array format
cc
cc        !Set grid=1 because xxx is NOT a MG vector
cc        do ieq=1,3
cc          call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
cc     .                           ,dv_cnv(:,:,:,ieq),1)
cc        enddo
cc
cc        call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ),icomp=IVX)
cc
ccc     Corrector step
cc
cc        !Find corrections to rho,B,tmp right hand sides (Udv)
cc        call correctRhoRHS(dv_cnv,rhs(:,IRHO)   ,igrid)
cc        call correctTmpRHS(dv_cnv,rhs(:,ITMP)   ,igrid)
cc        call correctBflRHS(dv_cnv,rhs(:,IBX:IBZ),igrid)
cc
cc        !Find better initial guesses
cc        xxx(:,IRHO)    = xxx(:,IRHO)    - dt*alpha*rhs(:,IRHO)
cc        xxx(:,ITMP)    = xxx(:,ITMP)    - dt*alpha*rhs(:,ITMP)
cc        xxx(:,IBX:IBZ) = xxx(:,IBX:IBZ) - dt*alpha*rhs(:,IBX:IBZ)
cc
cc        guess = 1
cc
cc        !Find corrected residuals
cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc              vol = volume(i,j,k,igx,igy,igz)
cc
cc              yyy(ii,IRHO)   =yyy(ii,IRHO)   -vol*alpha*rhs(ii,IRHO)
cc              yyy(ii,IBX:IBZ)=yyy(ii,IBX:IBZ)-vol*alpha*rhs(ii,IBX:IBZ)
cc              yyy(ii,ITMP)   =yyy(ii,ITMP)   -vol*alpha*rhs(ii,ITMP)
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc        !Temperature
cc        call cSolver(1,ntotp,yyy(:,ITMP),xxx(:,ITMP),bcs(:,ITMP)
cc     .              ,igrid,iout,guess,tmp_mtvc,tmp_diag,2,.false.)
cc
cc        !Density
cc        call cSolver(1,ntotp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
cc     .              ,igrid,iout,guess,rho_mtvc,rho_diag,2,.false.)
cc
cc        !Magnetic field
cc        call cSolver(3,ntotp,yyy(:,IBX:IBZ),xxx(:,IBX:IBZ)
cc     .           ,bcs(:,IBX:IBZ),igrid,iout,guess,b_mtvc,b_diag,2
cc     .           ,.false.)

c diag B-field divergence

cc        do ieq=1,3
cc          call mapMGVectorToArray(0,1,xxx(:,IBX+ieq-1),nx,ny,nz
cc     .                           ,db_cnv(:,:,:,ieq),1)
cc        enddo
cc
cc        call setMGBC(0,3,nx,ny,nz,igrid,db_cnv,bcs(:,IBX:IBZ),icomp=IBX)
cc
cc        do k=1,nz
cc          do j=1,ny
cc            do i=1,nx
cc              debug(i,j,k) = div(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     $                          ,db_cnv(:,:,:,1)
cc     $                          ,db_cnv(:,:,:,2)
cc     $                          ,db_cnv(:,:,:,3))
cc            enddo
cc          enddo
cc        enddo
cc
cc        mag = sqrt(sum(debug(1:nx,1:ny,1:nz)**2))
cc        write (*,*) mag

cc        write (*,*) 'plot div(b)?'
cc        read (*,'(a)') plot
cc        if (plot == 'y') then
cc          open(unit=110,file='debug.bin',form='unformatted'
cc     .       ,status='replace')
cc          call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,0,110)
cc          close(110)
cc          stop
cc        endif

c     Deallocate variables

        deallocate(dv_cnv,db_cnv,dj_cov)

      endif

c End program

      end subroutine applyPreconditioner

c averageScalarUpdates
c #########################################################################
      subroutine averageScalarUpdates(ntotp,xxx,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the rhs for the velocity SI solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use imposeBCinterface

      use mg_internal

      implicit none

c Call variables

      integer(4) :: ntotp,igrid
      real(8)    :: xxx(ntotp)

c Local variables

      integer(4) :: ii
      real(8)    :: mag,vol,tvol

c Average around singular point

      if (bcond(1) == SP) then
        i = 1
        do k=1,nz
          mag  = 0d0
          tvol = 0d0
          do j=1,ny
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            vol = volume(i,j,k,igrid,igrid,igrid)
            mag = mag + xxx(ii)*vol
            tvol = tvol + vol
          enddo
          mag = mag/tvol
          do j=1,ny
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            xxx(ii) = mag
          enddo
        enddo
      endif

c End program

      end subroutine averageScalarUpdates

c formSIrhs
c #########################################################################
      subroutine formSIrhs(ntotp,xxx,yyy,rhs_si,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the rhs for the velocity SI solve.
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use imposeBCinterface

      use mg_internal

      implicit none

c Call variables

      integer(4) :: ntotp,igrid
      real(8)    :: xxx(ntotp,neqd),yyy(ntotp,3),rhs_si(ntotp,3)

c Local variables

      integer(4) :: ii,iig,isig,ivar,ip,im
      real(8)    :: dvol,cov(3),cnv(3),idx,idy,idz
     .             ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,db_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dj_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dpres (0:nx+1,0:ny+1,0:nz+1)
      real(8),target :: db_cnv(0:nx+1,0:ny+1,0:nz+1,3)

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

c Find dB* (w/ BCs)

      !xxx is NOT a MG vector
      do ieq=1,3
        call mapMGVectorToArray(0,1,xxx(:,IBX+ieq-1),nx,ny,nz
     .                         ,db_cnv(:,:,:,ieq),igrid,.false.)
      enddo

      !Find covariant components of db with BCs
      call setBC(IBX,nx,ny,nz,db_cnv,db_cov,vzeros,bcs(:,IBX:IBZ)
     .          ,igrid,igrid,igrid,iorder=1)

c Find dj* from dB*

      !Find contravariant current (without BCs)
cc      dj_cnv = 0d0
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            do ivar=1,3
cc              dj_cnv(i,j,k,ivar)=curl2(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                                ,db_cov(:,:,:,1)
cc     .                                ,db_cov(:,:,:,2)
cc     .                                ,db_cov(:,:,:,3)
cc     .                                ,ivar)
cccc              dj_cnv(i,j,k,:)=curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
cccc     .                                ,db_cov(:,:,:,1)
cccc     .                                ,db_cov(:,:,:,2)
cccc     .                                ,db_cov(:,:,:,3))
cc            enddo
cc          enddo
cc        enddo
cc      enddo
cc
cc      !Find covariant components of dj with BCs
cc      call setBC(IJX,nx,ny,nz,dj_cnv,dj_cov,vzeros,bcs(:,IJX:IJZ)
cc     .          ,igrid,igrid,igrid,iorder=2)

c Find dpres array with BCs

      dpres = 0d0
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            dpres(i,j,k)= 2.*(rho(i,j,k)*xxx(ii,ITMP)
     .                       +tmp(i,j,k)*xxx(ii,IRHO))
          enddo
        enddo
      enddo

      !Find BCs
      call setBC(IRHO,nx,ny,nz,dpres,zeros,bcs(:,IRHO)
     .          ,igrid,igrid,igrid,iorder=1)

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

cc            if (isSP(i,j,k,igrid,igrid,igrid)) then
cc              idx = 1d0/dx(ig)
cc              im  = i
cc            endif
cc
cc            cov(1) = jy(i,j,k)*db_cnv(i,j,k,3)/jac
cc     .             - jz(i,j,k)*db_cnv(i,j,k,2)/jac
cc     .             - by(i,j,k)*dj_cnv(i,j,k,3)/jac
cc     .             + bz(i,j,k)*dj_cnv(i,j,k,2)/jac
cc     .             -(dpres(ip,j,k)-dpres(im,j,k))*idx
cc
cc            cov(2) = jz(i,j,k)*db_cnv(i,j,k,1)/jac
cc     .             - jx(i,j,k)*db_cnv(i,j,k,3)/jac
cc     .             - bz(i,j,k)*dj_cnv(i,j,k,1)/jac
cc     .             + bx(i,j,k)*dj_cnv(i,j,k,3)/jac
cc     .             -(dpres(i,j+1,k)-dpres(i,j-1,k))*idy
cc
cc            cov(3) = jx(i,j,k)*db_cnv(i,j,k,2)/jac
cc     .             - jy(i,j,k)*db_cnv(i,j,k,1)/jac
cc     .             - bx(i,j,k)*dj_cnv(i,j,k,2)/jac
cc     .             + by(i,j,k)*dj_cnv(i,j,k,1)/jac
cc     .             -(dpres(i,j,k+1)-dpres(i,j,k-1))*idz

            if (isSP(i,j,k,igrid,igrid,igrid)) then
              cov(1) = -(dpres(ip ,j,k)+dpres(i,j,k)
     .                  -2*dpres(im ,j,k))*idx
            else
              cov(1) = -(dpres(ip ,j,k)-dpres(im ,j,k))*idx
            endif

            cov(2) = -(dpres(i,j+1,k)-dpres(i,j-1,k))*idy

            cov(3) = -(dpres(i,j,k+1)-dpres(i,j,k-1))*idz

            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3),.true.)

            !Add Lorentz force correction
            vec1 => gb0%grid(igrid)%array
            vec2 => db_cnv
            cnv = cnv
     .          + div_tensor(i,j,k,nx,ny,nz,igrid,igrid,igrid,.false.
     .                      ,lf_x,lf_y,lf_z,vol=.false.)
            nullify(vec1,vec2)

            !Correct rhs_v
            rhs_si(ii,:) = rhs_si(ii,:)+alpha*dvol*cnv(:)

            !Transform to cartesian
            if (si_car) then  

              cnv(:) = rhs_si(ii,:)

              call transformVectorToCartesian
     .               (i,j,k,igrid,igrid,igrid
     .               ,cnv(1),cnv(2),cnv(3)
     .               ,.false.
     .               ,rhs_si(ii,1),rhs_si(ii,2),rhs_si(ii,3))

            endif

          enddo
        enddo
      enddo

c End program

      end subroutine formSIrhs

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
cc      real(8)    :: rhodv_cnv(0:nx+1,0:ny+1,0:nz+1,3)

c Begin program

c Find rho0*dv

cc      rhodv_cnv(:,:,:,1) = dv_cnv(:,:,:,1)*rho
cc      rhodv_cnv(:,:,:,2) = dv_cnv(:,:,:,2)*rho
cc      rhodv_cnv(:,:,:,3) = dv_cnv(:,:,:,3)*rho

c Evaluate rhs correction: div(rh0*dv)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

cc            crhs(ii) = div(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                    ,rhodv_cnv(:,:,:,1)
cc     .                    ,rhodv_cnv(:,:,:,2)
cc     .                    ,rhodv_cnv(:,:,:,3))

            ip = i+1
            im = i-1
            jp = j+1
            jm = j-1
            kp = k+1
            km = k-1

            jac = gmetric%grid(igrid)%jac(i,j,k)

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            dh(1) = 2*dxh(ig)
            dh(2) = 2*dyh(jg)
            dh(3) = 2*dzh(kg)

            if (isSP(i,j,k,igrid,igrid,igrid)) then
              dh(1) = dx(ig)
              im = i
            endif

            crhs(ii) = rho(i,j,k)
     .                *div(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                    ,dv_cnv(:,:,:,1)
     .                    ,dv_cnv(:,:,:,2)
     .                    ,dv_cnv(:,:,:,3))
     .              +dv_cnv(i,j,k,1)*(rho(ip,j,k)-rho(im,j,k))/dh(1)/jac
     .              +dv_cnv(i,j,k,2)*(rho(i,jp,k)-rho(i,jm,k))/dh(2)/jac
     .              +dv_cnv(i,j,k,3)*(rho(i,j,kp)-rho(i,j,km))/dh(3)/jac
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

c correctBflRHS
c ###################################################################
      subroutine correctBflRHS(dv_cnv,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds the corrected rhs for the magnetic field
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
            call find_curl_vxb(i,j,k,nx,ny,nz
     .                        ,dv_cnv
     .                        ,gb0%grid(igrid)%array
     $                        ,crhs(ii,1),crhs(ii,2),crhs(ii,3),0,igrid)
          enddo
        enddo
      enddo

c End program

      end subroutine correctBflRHS

c findDivfreeRHS
c ###################################################################
      subroutine findDivfreeRHS(dv_cnv,xxx,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds a div-free rhs correction for the magnetic
c     field:
c        rhs = div(dv B0 - B0 dv) + div(V0 db - db V0) + curl(eta dj)
c
c     The div-free magnetic field is found from:
c        db=dt*(-Gb - theta*rhs)
c
c     Note: div(dv B0 - B0 dv) is already calculated and stored in
c           crhs on input
c--------------------------------------------------------------------

      use precond_variables

      use operators

      use imposeBCinterface

      use mg_internal

      implicit none

c Call variables

      integer(4) :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,crhs(nx*ny*nz,3),xxx(nx*ny*nz,3)

c Local variables

      integer(4) :: ii,ivar
      real(8)    :: a1,a2,a3,etal
     .             ,db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,db_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dj_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3)

c Begin program

c Find dB* (w/ BCs)

      !xxx is NOT a MG vector
      do ieq=1,3
        call mapMGVectorToArray(0,1,xxx(:,IBX+ieq-1),nx,ny,nz
     .                         ,db_cnv(:,:,:,ieq),igrid,.false.)
      enddo

      !Find covariant components of db with BCs
      call setBC(IBX,nx,ny,nz,db_cnv,db_cov,vzeros,bcs(:,IBX:IBZ)
     .          ,igrid,igrid,igrid,iorder=2)

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
cc              dj_cnv(i,j,k,:)=curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                                ,db_cov(:,:,:,1)
cc     .                                ,db_cov(:,:,:,2)
cc     .                                ,db_cov(:,:,:,3))
            enddo
          enddo
        enddo
      enddo

      !Find covariant components of dj with BCs
      call setBC(IJX,nx,ny,nz,dj_cnv,dj_cov,vzeros,bcs(:,IJX:IJZ)
     .          ,igrid,igrid,igrid,iorder=2)

c Evaluate eta*dj

      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            etal = res(i,j,k,nx,ny,nz,igrid,igrid,igrid)
            dj_cov(i,j,k,:) = etal*dj_cov(i,j,k,:)
          enddo
        enddo
      enddo

c Evaluate div-free rhs correction: div(V0 db - db V0) + curl(eta dj)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            !div(dv B0 - B0 dv) ==> Already stored in crhs on input
c$$$            call find_curl_vxb(i,j,k,nx,ny,nz
c$$$     .                        ,dv_cnv
c$$$     .                        ,gb0%grid(igrid)%array
c$$$     $                        ,a1,a2,a3,0,igrid)
c$$$            crhs(ii,1) = a1
c$$$            crhs(ii,2) = a2
c$$$            crhs(ii,3) = a3

            !div(dv B0 - B0 dv) (<= crhs on input)  + div(V0 db - db V0)
            call find_curl_vxb(i,j,k,nx,ny,nz
     .                        ,gv0%grid(igrid)%array
     .                        ,db_cnv
     $                        ,a1,a2,a3,0,igrid)
            crhs(ii,1) = crhs(ii,1)+a1
            crhs(ii,2) = crhs(ii,2)+a2
            crhs(ii,3) = crhs(ii,3)+a3

            !curl(eta dj)
            a1= curl2(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .              ,dj_cov(:,:,:,1)
     .              ,dj_cov(:,:,:,2)
     .              ,dj_cov(:,:,:,3),1)
            a2= curl2(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .              ,dj_cov(:,:,:,1)
     .              ,dj_cov(:,:,:,2)
     .              ,dj_cov(:,:,:,3),2)
            a3= curl2(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .              ,dj_cov(:,:,:,1)
     .              ,dj_cov(:,:,:,2)
     .              ,dj_cov(:,:,:,3),3)
            crhs(ii,1) = crhs(ii,1)+a1
            crhs(ii,2) = crhs(ii,2)+a2
            crhs(ii,3) = crhs(ii,3)+a3

          enddo
        enddo
      enddo

c End program

      end subroutine findDivfreeRHS

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
