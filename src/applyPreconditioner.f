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
     .             ,rhs(ntot/neqd,neqd),rr(ntot),dxk(ntot),vol

      real(8),allocatable,dimension(:,:,:,:) :: dv_cnv,db_cnv,dj_cov

      integer(4) :: ii,iii,igrid,ntotp,guess,si_it

c Debug

      real(8)    :: debug(0:nxd+1,0:nyd+1,0:nzd+1),mag,mag2

c Externals

      external   :: tmp_mtvc,rho_mtvc,b_mtvc,v_mtvc,test_mtvc

c Interface

      INTERFACE
        subroutine cSolver(neq,ntotp,b,x,bcnd,igrid,out,guess
     $                    ,icmp,matvec,dg,ncolors)
         integer(4) :: neq,ntotp,igrid,bcnd(6,neq),out,icmp,guess
     $                ,ncolors
         real(8)    :: x(ntotp,neq),b(ntotp,neq)
         real(8), target :: dg(neq,2*neq*ntotp)
         external   :: matvec
        end subroutine cSolver
      END INTERFACE

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

        x = y*dt   !Diagonal scaling

c *******************************************************************
c     Semi-implicit preconditioner version I
c *******************************************************************

      elseif (precon == 's1') then

        xxx = 0d0

c     Set up size of operator split vectors

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              do ieq=1,neqd
                iii = ieq + neqd*(ii-1)

                yyy(ii,ieq) = y(iii)*volume(i,j,k,igx,igy,igz)  !Volume weighing of residuals
              enddo
            enddo
          enddo
        enddo

        guess = 0

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3))

        allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3))

c     Predictor step

        !Temperature
        call cSolver(1,ntotp,yyy(:,ITMP),xxx(:,ITMP),bcs(:,ITMP)
     .              ,igrid,iout,guess,ITMP,tmp_mtvc,tmp_diag,2)

        !Density
        call cSolver(1,ntotp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
     .              ,igrid,iout,guess,IRHO,rho_mtvc,rho_diag,2)

        !Magnetic field
        call cSolver(3,ntotp,yyy(:,IBX:IBZ),xxx(:,IBX:IBZ)
     .           ,bcs(:,IBX:IBZ),igrid,iout,guess,IBX,b_mtvc,b_diag,2)

c     SI step

        !Form SI rhs -(Gv+Ldx) (and also returns arrays db_cnv,dj_cov)
        call formSIrhs(ntotp,xxx,yyy(:,IVX:IVZ),rhs(:,IVX:IVZ)
     .                ,db_cnv,dj_cov,igrid)

        !Solve Schur-complement SI system
        call cSolver(3,ntotp,rhs(:,IVX:IVZ),xxx(:,IVX:IVZ)
     .           ,bcs(:,IVX:IVZ),igrid,iout,guess,IVX,v_mtvc,v_diag,4)

c     Store velocity solution in array format

        !xxx is NOT a MG vector
        do ieq=1,3
          call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
     .                           ,dv_cnv(:,:,:,ieq),igrid,.false.)
        enddo

        icomp = IVX
        call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ))

c     Corrector step

        !Find corrections to rho,B,tmp right hand sides (Udv)
        call correctRhoRHS(dv_cnv,rhs(:,IRHO)   ,igrid)
        call correctTmpRHS(dv_cnv,rhs(:,ITMP)   ,igrid)
        call correctBflRHS(dv_cnv,rhs(:,IBX:IBZ),igrid)

        !Find better initial guesses
        xxx(:,IRHO)    = xxx(:,IRHO)    - dt*alpha*rhs(:,IRHO)
        xxx(:,ITMP)    = xxx(:,ITMP)    - dt*alpha*rhs(:,ITMP)
        xxx(:,IBX:IBZ) = xxx(:,IBX:IBZ) - dt*alpha*rhs(:,IBX:IBZ)

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
cc     .              ,igrid,iout,guess,ITMP,tmp_mtvc,tmp_diag,2)
cc
cc        !Density
cc        call cSolver(1,ntotp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
cc     .              ,igrid,iout,guess,IRHO,rho_mtvc,rho_diag,2)
cc
cc        !Magnetic field
cc        call cSolver(3,ntotp,yyy(:,IBX:IBZ),xxx(:,IBX:IBZ)
cc     .           ,bcs(:,IBX:IBZ),igrid,iout,guess,IBX,b_mtvc,b_diag,2)

c     Postprocessing of magnetic field: divergence cleaning

cc        call findDivfreeRHS(dv_cnv,db_cnv,dj_cov,rhs(:,IBX:IBZ),igrid)
cc
cc        !Divergence-free correction for magnetic field (needs changes in correctBflRHS)
cc        do k = 1,nz
cc          do j = 1,ny
cc            do i = 1,nx
cc              ii  = i + nx*(j-1) + nx*ny*(k-1)
cc              do ieq = IBX,IBZ
cc                iii = ieq + neqd*(ii-1)
cc                xxx(ii,ieq) = dt*(y(iii) - alpha*rhs(ii,ieq))
cc              enddo
cc            enddo
cc          enddo
cc        enddo

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

cc        open(unit=110,file='debug.bin',form='unformatted'
cc     .       ,status='replace')
cc        call mapMGVectorToArray(0,1,xxx(:,IRHO),nx,ny,nz,debug,igrid,.false.)
cc        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,0,110)
cc        call mapMGVectorToArray(0,1,xxx(:,IVX),nx,ny,nz,debug,igrid,.false.)
cc        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc        call mapMGVectorToArray(0,1,xxx(:,IVY),nx,ny,nz,debug,igrid,.false.)
cc        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc        call mapMGVectorToArray(0,1,xxx(:,IVZ),nx,ny,nz,debug,igrid,.false.)
cc        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc        call mapMGVectorToArray(0,1,xxx(:,IBX),nx,ny,nz,debug,igrid,.false.)
cc        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc        call mapMGVectorToArray(0,1,xxx(:,IBY),nx,ny,nz,debug,igrid,.false.)
cc        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc        call mapMGVectorToArray(0,1,xxx(:,IBZ),nx,ny,nz,debug,igrid,.false.)
cc        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc        call mapMGVectorToArray(0,1,xxx(:,ITMP),nx,ny,nz,debug,igrid,.false.)
cc        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc        close(110)
cc        stop

c diag B-field divergence

        do ieq=1,3
          call mapMGVectorToArray(0,1,xxx(:,IBX+ieq-1),nx,ny,nz
     .                           ,db_cnv(:,:,:,ieq),igrid,.false.)
        enddo

        icomp = IBX
        call setMGBC(0,3,nx,ny,nz,igrid,db_cnv,bcs(:,IBX:IBZ))

        do k=1,nz
          do j=1,ny
            do i=1,nx
              debug(i,j,k) = div(i,j,k,nx,ny,nz,db_cnv(:,:,:,1)
     $                                         ,db_cnv(:,:,:,2)
     $                                         ,db_cnv(:,:,:,3))
            enddo
          enddo
        enddo

        mag = sqrt(sum(debug(1:nx,1:ny,1:nz)**2))
        write (*,*) mag

cc        write (*,*) 'plot div(b)?'
cc        read (*,'(a)') plot
cc        if (plot == 'y') then
cc          open(unit=110,file='debug.bin',form='unformatted'
cc     .       ,status='replace')
cc          call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,0,110)
cc          close(110)
cc          stop
cc        endif

        deallocate(dv_cnv,db_cnv,dj_cov)

c     Store solution in "x"

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              do ieq=1,neqd
                iii = ieq + neqd*(ii-1)
                x(iii) = xxx(ii,ieq)
              enddo
            enddo
          enddo
        enddo

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
            call matrixFreeMatVec(ntot,xk,x,rr)
            rr = y - rr

            mag = sqrt(sum(rr*rr))
            mag = mag/sqrt(sum(y*y))
            write (*,*) si_it-1,mag
          endif

c       Map residual to individual components (rr --> yyy)

          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                vol = volume(i,j,k,igx,igy,igz)
                do ieq=1,neqd
                  iii = ieq + neqd*(ii-1)
                  yyy(ii,ieq) = rr(iii)*vol !Volume weighing of residuals
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

          !Laplacian test
cc          bcs(1,IRHO) = DIR
cc          bcs(2,IRHO) = DIR
cc          bcs(3,IRHO) = DIR
cc          bcs(4,IRHO) = DIR
          call cSolver(1,ntotp,rhs(:,IRHO),xxx(:,IRHO)
     .          ,bcs(:,IRHO),igrid,iout,guess,IRHO,test_mtvc,rho_diag,2)
          call MGplot(1,xxx(:,IRHO),igrid,0,'debug.bin')
          stop

          !Solve Schur-complement SI system ---> xxx(:,IVX:IVZ)
          call cSolver(3,ntotp,rhs(:,IVX:IVZ),xxx(:,IVX:IVZ)
     .           ,bcs(:,IVX:IVZ),igrid,iout,guess,IVX,v_mtvc,v_diag,4)


        call MGplot(1,xxx(:,IVX),igrid,0,'debug.bin')
        call MGplot(1,xxx(:,IVY),igrid,1,'debug.bin')
        call MGplot(1,xxx(:,IVZ),igrid,1,'debug.bin')

        stop

c       Store velocity solution in array format --> dv_cnv

          !Set grid=1 because xxx is NOT a MG vector
          do ieq=1,3
            call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
     .                             ,dv_cnv(:,:,:,ieq),igrid,.false.)
          enddo

          icomp = IVX
          call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ))

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

          !Find Deltax = xxx - yyy
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                vol = volume(i,j,k,igx,igy,igz)
                xxx(ii,IRHO)    = 
     .                       xxx(ii,IRHO)    - vol*alpha*yyy(ii,IRHO)
                xxx(ii,ITMP)    =
     .                       xxx(ii,ITMP)    - vol*alpha*yyy(ii,ITMP)
                xxx(ii,IBX:IBZ) =
     .                       xxx(ii,IBX:IBZ) - vol*alpha*yyy(ii,IBX:IBZ)
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
        call matrixFreeMatVec(ntot,xk,x,rr)

        rr = y - rr

        mag = sqrt(sum(rr*rr))
        mag = mag/sqrt(sum(y*y))
        write (*,*) si_it-1,mag
c diag ******

c plot ************
c$$$        if (mag > 1d0) then
c$$$        do k = 1,nz
c$$$          do j = 1,ny
c$$$            do i = 1,nx
c$$$              ii  = i + nx*(j-1) + nx*ny*(k-1)
c$$$              do ieq=1,neqd
c$$$                iii = ieq + neqd*(ii-1)
c$$$                xxx(ii,ieq) = y(iii)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$
c$$$        open(unit=110,file='debug.bin',form='unformatted'
c$$$     .       ,status='replace')
c$$$        call mapMGVectorToArray(0,1,xxx(:,IRHO),nx,ny,nz,debug,1)
c$$$        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,0,110)
c$$$        call mapMGVectorToArray(0,1,xxx(:,IVX),nx,ny,nz,debug,1)
c$$$        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
c$$$        call mapMGVectorToArray(0,1,xxx(:,IVY),nx,ny,nz,debug,1)
c$$$        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
c$$$        call mapMGVectorToArray(0,1,xxx(:,IVZ),nx,ny,nz,debug,1)
c$$$        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
c$$$        call mapMGVectorToArray(0,1,xxx(:,IBX),nx,ny,nz,debug,1)
c$$$        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
c$$$        call mapMGVectorToArray(0,1,xxx(:,IBY),nx,ny,nz,debug,1)
c$$$        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
c$$$        call mapMGVectorToArray(0,1,xxx(:,IBZ),nx,ny,nz,debug,1)
c$$$        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
c$$$        call mapMGVectorToArray(0,1,xxx(:,ITMP),nx,ny,nz,debug,1)
c$$$        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
c$$$        close(110)
c$$$        stop
c$$$        endif
c plot *************

ccc     Store velocity solution in array format
cc
cc        !Set grid=1 because xxx is NOT a MG vector
cc        do ieq=1,3
cc          call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
cc     .                           ,dv_cnv(:,:,:,ieq),1)
cc        enddo
cc
cc        icomp = IVX
cc        call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ))
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
cc     .              ,igrid,iout,guess,ITMP,tmp_mtvc,tmp_diag,2)
cc
cc        !Density
cc        call cSolver(1,ntotp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
cc     .              ,igrid,iout,guess,IRHO,rho_mtvc,rho_diag,2)
cc
cc        !Magnetic field
cc        call cSolver(3,ntotp,yyy(:,IBX:IBZ),xxx(:,IBX:IBZ)
cc     .           ,bcs(:,IBX:IBZ),igrid,iout,guess,IBX,b_mtvc,b_diag,2)

c diag B-field divergence

cc        do ieq=1,3
cc          call mapMGVectorToArray(0,1,xxx(:,IBX+ieq-1),nx,ny,nz
cc     .                           ,db_cnv(:,:,:,ieq),1)
cc        enddo
cc
cc        icomp = IBX
cc        call setMGBC(0,3,nx,ny,nz,igrid,db_cnv,bcs(:,IBX:IBZ))
cc
cc        do k=1,nz
cc          do j=1,ny
cc            do i=1,nx
cc              debug(i,j,k) = div(i,j,k,nx,ny,nz,db_cnv(:,:,:,1)
cc     $                                         ,db_cnv(:,:,:,2)
cc     $                                         ,db_cnv(:,:,:,3))
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

        deallocate(dv_cnv,db_cnv,dj_cov)

      endif

c End program

      end subroutine applyPreconditioner

c cSolver
c #########################################################################
      subroutine cSolver(neq,ntotp,b,x,bcnd,igrid,out,guess
     $                  ,icmp,matvec,dg,ncolors)
c--------------------------------------------------------------------
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
c       * icmp: defines magnitudes that are being solved. Relevant
c               for boundary condition subroutine.
c       * matvec (external): matrix-vector product definition.
c       * dg: matrix neq*neq diagonal block (for stationary its).
c       * ncolors: number of colors in grid (for GS).
c--------------------------------------------------------------------

      use precond_variables

      use mlsolverSetup

      implicit none

c Call variables

      integer(4) :: neq,ntotp,igrid,bcnd(6,neq),out,icmp,guess,ncolors
      real(8)    :: x(ntotp,neq),b(ntotp,neq)
      real(8), target :: dg(neq,2*neq*ntotp)

      external   :: matvec

c Local variables

      integer(4) :: ntot
      real(8)    :: xi(ntotp*neq),bi(ntotp*neq)

      real(8)    :: drand
      external   :: drand

c Begin program

      icomp = icmp  !Define icomp for setMGBC

c Interlace variables for coupled solve

cc      do i=1,ntotp
cc        do ieq=1,neq
cc          xi(neq*(i-1)+ieq) = x(i,ieq)
cc          bi(neq*(i-1)+ieq) = b(i,ieq)
cc        enddo
cc      enddo

c diag *** Convergence test
      bi = 0d0
      call random_number(xi)
      guess = 1
c diag ***

c Solve coupled MG

c     Initialize solver

      call solverInit

c     Upper_level solver options

c$$$      call solverOptionsInit
c$$$
c$$$      solverOptions%tol      = mgtol
c$$$
c$$$      call assembleSolverHierarchy('gm')

      call solverOptionsInit

      solverOptions%tol      = mgtol
      solverOptions%vcyc     = maxvcyc
      solverOptions%igridmin = 3
      solverOptions%orderres = 0
      solverOptions%orderprol= 2
      solverOptions%mg_mu    = 1
      solverOptions%vol_res  = .true.
cc      solverOptions%diag     => dg
      solverOptions%mg_coarse_solver_depth = 3
cc      solverOptions%vertex_based_relax = .true.
      solverOptions%mg_line_relax = .true.

      call assembleSolverHierarchy('mg')

c     Next level solver (smoother)

cc      call solverOptionsInit
cc
cc      solverOptions%tol             = mgtol
cc      solverOptions%stp_test        = 1 
cc      solverOptions%krylov_subspace = nsweep
cc      solverOptions%iter            = nsweep
cccc      solverOptions%krylov_subspace = 100
cccc      solverOptions%iter            = 100
cc
cc      call assembleSolverHierarchy('gm')

c     Next level solver (preconditioner for GMRES)

      call solverOptionsInit

      solverOptions%iter    = nsweep
      solverOptions%omega   = 1d0
      solverOptions%tol     = mgtol
cc      solverOptions%ncolors = ncolors
cc      solverOptions%diag    => dg
cc      solverOptions%vertex_based_relax = .false.

      call assembleSolverHierarchy('gs')

cc      solverOptions%omega   = 0.7
cc      solverOptions%iter    = 2*nsweep/solverOptions%omega
cc      solverOptions%tol     = mgtol
cccc      solverOptions%diag    => dg
cc
cc      call assembleSolverHierarchy('jb')

c     Coarsest grid solve for outer MG

      call solverOptionsInit

      solverOptions%tol             = 1d-10
      solverOptions%krylov_subspace = 1000
      solverOptions%iter            = 1000
      solverOptions%stp_test        = 1 

      call assembleSolverHierarchy('gm')
      call assembleSolverHierarchy('id') !GMRES preconditioner

c     Coarsest grid solve for line MG

      call solverOptionsInit

      solverOptions%iter    = 1000
      solverOptions%omega   = 1d0
      solverOptions%tol     = 1d-4
      solverOptions%ncolors = 1
cc      solverOptions%ncolors = ncolors
cc      solverOptions%diag    => dg
cc      solverOptions%vertex_based_relax = .false.

      call assembleSolverHierarchy('jb')

c     Invoke solver

      ntot=neq*ntotp
      call getSolver(neq,ntot,bi,xi,matvec,igrid,bcnd,guess,out,1)

c     Get output data

cc      call getSolverOptions(1)

c     Kill solver

      call solverKill

c Unravel solution for output

      do i = 1,ntotp
        do ieq=1,neq
          x(i,ieq) = xi(neq*(i-1)+ieq)
        enddo
      enddo

c End program

      end subroutine cSolver

c formSIrhs
c #########################################################################
      subroutine formSIrhs(ntotp,xxx,yyy,rhs_si,db_cnv,dj_cov,igrid)
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
      real(8)    :: db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3)

c Local variables

      integer(4) :: ii,iig,isig
      real(8)    :: dvol,cov(3),cnv(3)
     .             ,db_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dj_cnv(0:nx+1,0:ny+1,0:nz+1,3)
      real(8),allocatable,dimension(:,:,:)   :: dpres

c Begin program

      isig = grid_params%istartp(igrid)

c Find rhs_v

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            iig = ii + isig - 1
            dvol   = volume(i,j,k,igx,igy,igz)
            rhs_si(ii,:) = yyy(ii,:)
     .                   - dvol*xxx(ii,IRHO)*mgadvdiffV0(iig,:)

          enddo
        enddo
      enddo

c Find dj* from dB*

      !xxx is NOT a MG vector
      do ieq=1,3
        call mapMGVectorToArray(0,1,xxx(:,IBX+ieq-1),nx,ny,nz
     .                         ,db_cnv(:,:,:,ieq),igrid,.false.)
      enddo

      !Find covariant components of db with BCs
      iimin = 1
      iimax = nx
      jjmin = 1
      jjmax = ny
      kkmin = 1
      kkmax = nz
      call setBC(IBX,db_cnv,db_cov,vzeros,bcs(:,IBX:IBZ))

      !Find contravariant current (without BCs)
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            do icomp=1,3
              dj_cnv(i,j,k,icomp)=curl2(i,j,k,nx,ny,nz,db_cov(:,:,:,1)
     .                                                ,db_cov(:,:,:,2)
     .                                                ,db_cov(:,:,:,3)
     .                                 ,icomp)
            enddo
          enddo
        enddo
      enddo

      !Find covariant components of dj with BCs
      iimin = 1
      iimax = nx
      jjmin = 1
      jjmax = ny
      kkmin = 1
      kkmax = nz
      call setBC(IJX,dj_cnv,dj_cov,vzeros,bcs(:,IJX:IJZ))

c Find dpres array with BCs

      allocate(dpres(0:nx+1,0:ny+1,0:nz+1))

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
      iimin = 1
      iimax = nx
      jjmin = 1
      jjmax = ny
      kkmin = 1
      kkmax = nz
      call setBC(IRHO,dpres,zeros,bcs(:,IRHO))

c Find rhs_v'

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            jac    = gmetric%grid(igx)%jac(i,j,k)
            dvol   = volume(i,j,k,igx,igy,igz)

            cov(1) = jy(i,j,k)*db_cnv(i,j,k,3)/jac
     .             - jz(i,j,k)*db_cnv(i,j,k,2)/jac
     .             - by(i,j,k)*dj_cnv(i,j,k,3)/jac
     .             + bz(i,j,k)*dj_cnv(i,j,k,2)/jac
     .             -(dpres(i+1,j,k)-dpres(i-1,j,k))/dxh(ig)*0.5

            cov(2) = jz(i,j,k)*db_cnv(i,j,k,1)/jac
     .             - jx(i,j,k)*db_cnv(i,j,k,3)/jac
     .             - bz(i,j,k)*dj_cnv(i,j,k,1)/jac
     .             + bx(i,j,k)*dj_cnv(i,j,k,3)/jac
     .             -(dpres(i,j+1,k)-dpres(i,j-1,k))/dyh(jg)*0.5

            cov(3) = jx(i,j,k)*db_cnv(i,j,k,2)/jac
     .             - jy(i,j,k)*db_cnv(i,j,k,1)/jac
     .             - bx(i,j,k)*dj_cnv(i,j,k,2)/jac
     .             + by(i,j,k)*dj_cnv(i,j,k,1)/jac
     .             -(dpres(i,j,k+1)-dpres(i,j,k-1))/dzh(kg)*0.5

            !Transform to cov to cnv
            call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3),.true.)

            !Correct rhs_v
            rhs_si(ii,:) = rhs_si(ii,:)+alpha*dvol*cnv(:)

          enddo
        enddo
      enddo

      deallocate(dpres)

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

      use imposeBCinterface

      use mg_internal

      implicit none

c Call variables

      integer(4) :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz)

c Local variables

      integer(4) :: ii,iig
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
            crhs(ii) = div(i,j,k,nx,ny,nz,rhodv_cnv(:,:,:,1)
     $                                   ,rhodv_cnv(:,:,:,2)
     $                                   ,rhodv_cnv(:,:,:,3))
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

      use imposeBCinterface

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
     $                           *div(i,j,k,nx,ny,nz,dv_cnv(:,:,:,1)
     $                                              ,dv_cnv(:,:,:,2)
     $                                              ,dv_cnv(:,:,:,3))
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

      use imposeBCinterface

      use mg_internal

      implicit none

c Call variables

      integer(4) :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz,3)

c Local variables

      integer(4) :: ii
      real(8)    :: a1,a2,a3,etal

c Externals

      real(8)    :: res
      external   :: res

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            call find_curl_vxb(i,j,k,nx,ny,nz
     .                        ,dv_cnv
     .                        ,gb0%grid(igrid)%array
     $                        ,a1,a2,a3,0,igrid)
            crhs(ii,1) = a1
            crhs(ii,2) = a2
            crhs(ii,3) = a3
          enddo
        enddo
      enddo

c End program

      end subroutine correctBflRHS

c findDivfreeRHS
c ###################################################################
      subroutine findDivfreeRHS(dv_cnv,db_cnv,dj_cov,crhs,igrid)
c--------------------------------------------------------------------
c     This subroutine finds a div-free rhs correction for the magnetic
c     field:
c        rhs = div(dv B0 - B0 dv) + div(V0 db - db V0) + curl(eta dj)
c
c     The div-free magnetic field is found from:
c        db=dt*(-Gb + theta*rhs)
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
     .             ,db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,crhs(nx*ny*nz,3)

c Local variables

      integer(4) :: ii
      real(8)    :: a1,a2,a3,etal

c Externals

      real(8)    :: res
      external   :: res

c Begin program

c Evaluate eta*dj

      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            do icomp=1,3
              etal = res(i,j,k,nx,ny,nz,igx,igy,igz)
              dj_cov(i,j,k,:) = etal*dj_cov(i,j,k,:)
            enddo
          enddo
        enddo
      enddo

c Evaluate div-free rhs correction: div(V0 db - db V0) + curl(eta dj)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
cc            call find_curl_vxb(i,j,k,nx,ny,nz
cc     .                        ,dv_cnv
cc     .                        ,gb0%grid(igrid)%array
cc     $                        ,a1,a2,a3,0,igrid)
cc            crhs(ii,1) = a1
cc            crhs(ii,2) = a2
cc            crhs(ii,3) = a3

            !div(V0 db - db V0)
            call find_curl_vxb(i,j,k,nx,ny,nz
     .                        ,gv0%grid(igrid)%array
     .                        ,db_cnv
     $                        ,a1,a2,a3,0,igrid)
            crhs(ii,1) = crhs(ii,1)+a1
            crhs(ii,2) = crhs(ii,2)+a2
            crhs(ii,3) = crhs(ii,3)+a3

            !curl(eta dj)
            a1= curl(i,j,k,nx,ny,nz,dj_cov(:,:,:,1)
     .                             ,dj_cov(:,:,:,2)
     .                             ,dj_cov(:,:,:,3),1)
            a2= curl(i,j,k,nx,ny,nz,dj_cov(:,:,:,1)
     .                             ,dj_cov(:,:,:,2)
     .                             ,dj_cov(:,:,:,3),2)
            a3= curl(i,j,k,nx,ny,nz,dj_cov(:,:,:,1)
     .                             ,dj_cov(:,:,:,2)
     .                             ,dj_cov(:,:,:,3),3)
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
      real(8)    :: dummy(neq)
      logical    :: fpointers

c Begin program

      call allocPointers(neq,grid_params,fpointers)

      isig = istart(igrid)

      do ii = 1,ntotp

        iii = neq*(ii-1)
        iig = iii + isig - 1

        dummy = y(ii,:)
        x(ii,:) = matmul(idiag(:,iig+1:iig+neq),dummy)

      enddo

      call deallocPointers(fpointers)

c End program

      end subroutine diagonalScaling
