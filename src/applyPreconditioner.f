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
     .             ,rhs(ntot/neqd,neqd),rr(ntot),dxk(ntot),dvol

      real(8),allocatable,dimension(:,:,:,:) :: dv_cnv,db_cnv,dj_cov

      integer(4) :: ii,iii,igrid,ntotp,guess,si_it

c Debug

      real(8)    :: mag,mag2,mag0

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

        x = y*dt   !Diagonal scaling

c *******************************************************************
c     Semi-implicit preconditioner version I
c *******************************************************************

      elseif (precon == 's1') then

        xxx = 0d0

c     Scatter residuals

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              if (vol_wgt) then  !Volume weighing of residuals

                dvol = volume(i,j,k,igx,igy,igz)
                do ieq=1,neqd
                  iii = ieq + neqd*(ii-1)
                  yyy(ii,ieq) = y(iii)*dvol
                enddo

cc                !Scalars
cc                dvol = volume(i,j,k,igx,igy,igz)
cc
cc                ieq=IRHO
cc                iii = ieq + neqd*(ii-1)
cc                yyy(ii,ieq) = y(iii)*dvol
cc
cc                ieq=ITMP
cc                iii = ieq + neqd*(ii-1)
cc                yyy(ii,ieq) = y(iii)*dvol
cc
cc                !Velocity
cc                dvol = dvol/gmetric%grid(igrid)%jac(i,j,k)
cc
cc                ieq=IVX
cc                iii = ieq + neqd*(ii-1)
cc                yyy(ii,ieq) = y(iii)*dvol
cc
cc                ieq=IVY
cc                iii = ieq + neqd*(ii-1)
cc                if (alt_eom) then
cc                  yyy(ii,ieq)=y(iii)*dvol*gmetric%grid(igrid)%jac(i,j,k)
cc                else
cc                  yyy(ii,ieq)=y(iii)*dvol
cc                endif
cc
cc                ieq=IVZ
cc                iii = ieq + neqd*(ii-1)
cc                yyy(ii,ieq) = y(iii)*dvol
cc
cc                !Magnetic field
cc                do ieq=IBX,IBZ
cc                  iii = ieq + neqd*(ii-1)
cc                  yyy(ii,ieq) = y(iii)*dvol
cc                enddo

              else

                do ieq=1,neqd
                  iii = ieq + neqd*(ii-1)
                  yyy(ii,ieq) = y(iii)
                enddo

              endif
            enddo
          enddo
        enddo

        guess = 0

c diag ****
cc        if (jit == 4) then
cc        !Solution plot
cc          call MGplot(1,yyy(:,IRHO),igrid,0,'debug.bin')
cc          call MGplot(1,yyy(:,IVX) ,igrid,1,'debug.bin')
cc          call MGplot(1,yyy(:,IVY) ,igrid,1,'debug.bin')
cc          call MGplot(1,yyy(:,IVZ) ,igrid,1,'debug.bin')
cc          call MGplot(1,yyy(:,IBX) ,igrid,1,'debug.bin')
cc          call MGplot(1,yyy(:,IBY) ,igrid,1,'debug.bin')
cc          call MGplot(1,yyy(:,IBZ) ,igrid,1,'debug.bin')
cc          call MGplot(1,yyy(:,ITMP),igrid,1,'debug.bin')
cc          stop
cc        endif
c diag ****

c     Create auxiliary arrays

        allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3))

        allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3))

c     Predictor step

        !Temperature
        call cSolver(1,ntotp,yyy(:,ITMP),xxx(:,ITMP),bcs(:,ITMP)
     .              ,igrid,iout,guess,tmp_mtvc,dg=tmp_diag,ncolors=2)

        !Density
        call cSolver(1,ntotp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
     .              ,igrid,iout,guess,rho_mtvc,dg=rho_diag,ncolors=2)

        !Magnetic field
        call cSolver(3,ntotp,yyy(:,IBX:IBZ),xxx(:,IBX:IBZ)
     .              ,bcs(:,IBX:IBZ),igrid,iout,guess,b_mtvc,dg=b_diag
     .              ,ncolors=4)

c     SI step

        !Form SI rhs -(Gv+Ldx) (and also returns arrays db_cnv,dj_cov)
        call formSIrhs(ntotp,xxx,yyy(:,IVX:IVZ),rhs(:,IVX:IVZ)
     .                ,db_cnv,dj_cov,igrid)

        !Solve Schur-complement SI system
        if (gm_smooth) then
          call cSolver(3,ntotp,rhs(:,IVX:IVZ),xxx(:,IVX:IVZ)
     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc
     .                ,ncolors=4,gm_smth=gm_smooth,line_relax=.false.)
        else
          call cSolver(3,ntotp,rhs(:,IVX:IVZ),xxx(:,IVX:IVZ)
     .                ,bcs(:,IVX:IVZ),igrid,iout,guess,v_mtvc,dg=v_diag
     .                ,ncolors=4,line_relax=.false.)
        endif

c diag ****
cc        !Solution plot
c         if (debug) then
cc          call MGplot(1,xxx(:,IRHO),igrid,0,'debug.bin')
cc          call MGplot(1,xxx(:,IVX) ,igrid,1,'debug.bin')
cc          call MGplot(1,xxx(:,IVY) ,igrid,1,'debug.bin')
cc          call MGplot(1,xxx(:,IVZ) ,igrid,1,'debug.bin')
cc          call MGplot(1,xxx(:,IBX) ,igrid,1,'debug.bin')
cc          call MGplot(1,xxx(:,IBY) ,igrid,1,'debug.bin')
cc          call MGplot(1,xxx(:,IBZ) ,igrid,1,'debug.bin')
cc          call MGplot(1,xxx(:,ITMP),igrid,1,'debug.bin')
cc          stop
cc        endif
c diag ****

c     Store velocity solution in array format

        !xxx is NOT a MG vector
        do ieq=1,3
          call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
     .                           ,dv_cnv(:,:,:,ieq),igrid,.false.)
        enddo

        call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ),icomp=IVX)

c     Corrector step

        !Find corrections to rho,B,tmp right hand sides (Udv)
        call correctRhoRHS(dv_cnv,rhs(:,IRHO)   ,igrid)
        call correctTmpRHS(dv_cnv,rhs(:,ITMP)   ,igrid)
        call correctBflRHS(dv_cnv,rhs(:,IBX:IBZ),igrid)

        !Find better initial guesses
        xxx(:,IRHO)    = xxx(:,IRHO)    - dt*alpha*rhs(:,IRHO)
        xxx(:,ITMP)    = xxx(:,ITMP)    - dt*alpha*rhs(:,ITMP)
        xxx(:,IBX:IBZ) = xxx(:,IBX:IBZ) - dt*alpha*rhs(:,IBX:IBZ)

c diag ****
cc        if (jit == 4.and.debug) then
cc          call MGplot(1,rhs(:,IRHO),igrid,0,'debug.bin')
cc          call MGplot(1,rhs(:,IBX) ,igrid,1,'debug.bin')
cc          call MGplot(1,rhs(:,IBY) ,igrid,1,'debug.bin')
cc          call MGplot(1,rhs(:,IBZ) ,igrid,1,'debug.bin')
cc          stop
cc        endif
c diag ****

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

c     Postprocessing of magnetic field: divergence cleaning

cc        call findDivfreeRHS(dv_cnv,db_cnv,dj_cov,rhs(:,IBX:IBZ),igrid)
cc
cc        !Divergence-free correction for magnetic field
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
cc              xxx(ii,IVX)=(rho(i,j,k)+xxx(ii,IRHO))
cc     .                   *(xxx(ii,IVX)+vx(i,j,k)) - rho(i,j,k)*vx(i,j,k)
cc              xxx(ii,IVY)=(rho(i,j,k)+xxx(ii,IRHO))
cc     .                   *(xxx(ii,IVY)+vy(i,j,k)) - rho(i,j,k)*vy(i,j,k)
cc              xxx(ii,IVZ)=(rho(i,j,k)+xxx(ii,IRHO))
cc     .                   *(xxx(ii,IVZ)+vz(i,j,k)) - rho(i,j,k)*vz(i,j,k)
            enddo
          enddo
        enddo

c     Diagnostics

c diag ****
        if (jit == 1.and.debug) then
        !Solution plot
          call MGplot(1,xxx(:,IRHO),igrid,0,'debug.bin')
          call MGplot(1,xxx(:,IVX) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,IVY) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,IVZ) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,IBX) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,IBY) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,IBZ) ,igrid,1,'debug.bin')
cc          call MGplot(1,rhs(:,IBX) ,igrid,1,'debug.bin')
cc          call MGplot(1,rhs(:,IBY) ,igrid,1,'debug.bin')
cc          call MGplot(1,rhs(:,IBZ) ,igrid,1,'debug.bin')
          call MGplot(1,xxx(:,ITMP),igrid,1,'debug.bin')
          stop
        endif
c diag ****

        !diag B-field divergence
c diag ****
cc        do ieq=1,3
cc          call mapMGVectorToArray(0,1,xxx(:,IBX+ieq-1),nx,ny,nz
cc     .                           ,db_cnv(:,:,:,ieq),igrid,.false.)
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
cc
cc        write (*,*) 'plot div(b)?'
cc        read (*,'(a)') plot
cc        if (plot == 'y') then
cc          open(unit=110,file='debug.bin',form='unformatted'
cc     .       ,status='replace')
cc          call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,0,110)
cc          close(110)
cc          stop
cc        endif
c diag ****

c     Deallocate variables

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
     .           ,ncolors=4,line_relax=.true.)

c       Store velocity solution in array format --> dv_cnv

          !Set grid=1 because xxx is NOT a MG vector
          do ieq=1,3
            call mapMGVectorToArray(0,1,xxx(:,IVX+ieq-1),nx,ny,nz
     .                             ,dv_cnv(:,:,:,ieq),igrid,.false.)
          enddo

          call setMGBC(0,3,nx,ny,nz,igrid,dv_cnv,bcs(:,IVX:IVZ)
     .                ,icomp=IVX)

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
        call matrixFreeMatVec(ntot,xk,x,rr)

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

      integer(4) :: ii,iig,isig,ivar,ip,im
      real(8)    :: dvol(3),cov(3),cnv(3),idx,idy,idz
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
            if (vol_wgt) then
              dvol = volume(i,j,k,igrid,igrid,igrid)
cc              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
cc              dvol = dxh(ig)*dyh(jg)*dzh(kg)
cc              if (alt_eom) dvol(2) = volume(i,j,k,igrid,igrid,igrid)
            else
              dvol = 1d0
            endif
            rhs_si(ii,:) = yyy(ii,:)
     .                   - dvol(:)*xxx(ii,IRHO)*mgadvdiffV0(iig,:)

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
      call setBC(IBX,nx,ny,nz,db_cnv,db_cov,vzeros,bcs(:,IBX:IBZ)
     .          ,igrid,igrid,igrid)

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
          enddo
        enddo
      enddo

      !Find covariant components of dj with BCs
      call setBC(IJX,nx,ny,nz,dj_cnv,dj_cov,vzeros,bcs(:,IJX:IJZ)
     .          ,igrid,igrid,igrid)

c Find dpres array with BCs

      allocate(dpres(0:nx+1,0:ny+1,0:nz+1))

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
cc              dvol = dxh(ig)*dyh(jg)*dzh(kg)
cc              if (alt_eom) dvol(2) = volume(i,j,k,igrid,igrid,igrid)
            else
              dvol = 1d0
            endif

            ip = i+1
            im = i-1

            idx = 0.5/dxh(ig)
            idy = 0.5/dyh(jg)
            idz = 0.5/dzh(kg)

            if (isSP(i,j,k,igrid,igrid,igrid)) then
              idx=1./dx(ig)
              im = i
            endif

            cov(1) = jy(i,j,k)*db_cnv(i,j,k,3)/jac
     .             - jz(i,j,k)*db_cnv(i,j,k,2)/jac
     .             - by(i,j,k)*dj_cnv(i,j,k,3)/jac
     .             + bz(i,j,k)*dj_cnv(i,j,k,2)/jac
     .             -(dpres(ip,j,k)-dpres(im,j,k))*idx

            cov(2) = jz(i,j,k)*db_cnv(i,j,k,1)/jac
     .             - jx(i,j,k)*db_cnv(i,j,k,3)/jac
     .             - bz(i,j,k)*dj_cnv(i,j,k,1)/jac
     .             + bx(i,j,k)*dj_cnv(i,j,k,3)/jac
     .             -(dpres(i,j+1,k)-dpres(i,j-1,k))*idy

            cov(3) = jx(i,j,k)*db_cnv(i,j,k,2)/jac
     .             - jy(i,j,k)*db_cnv(i,j,k,1)/jac
     .             - bx(i,j,k)*dj_cnv(i,j,k,2)/jac
     .             + by(i,j,k)*dj_cnv(i,j,k,1)/jac
     .             -(dpres(i,j,k+1)-dpres(i,j,k-1))*idz

            !Transform to cov to cnv
            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3),.true.)

            !Correct rhs_v
            rhs_si(ii,:) = rhs_si(ii,:)+alpha*dvol(:)*cnv(:)

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

      use mg_internal

      implicit none

c Call variables

      integer(4) :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3),crhs(nx*ny*nz)

c Local variables

      integer(4) :: ii,iig,ip,im,jp,jm,kp,km
      real(8)    :: rhodv_cnv(0:nx+1,0:ny+1,0:nz+1,3),dh(3)

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
     $                    ,rhodv_cnv(:,:,:,2)
     $                    ,rhodv_cnv(:,:,:,3))

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
cc            dh(1) = dxh(ig)
cc            dh(2) = dyh(jg)
cc            dh(3) = dzh(kg)
cc
cc            if (isSP(i,j,k,igx,igy,igz)) then
cc              dh(1) = dx(ig)
cc              im = i
cc            endif
cc
cc            crhs(ii) = rho(i,j,k)
cc     .                *div(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                    ,dv_cnv(:,:,:,1)
cc     $                    ,dv_cnv(:,:,:,2)
cc     $                    ,dv_cnv(:,:,:,3))
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
      real(8)    :: a1,a2,a3,etal

c Externals

      real(8)    :: res
      external      res

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

      use mg_internal

      implicit none

c Call variables

      integer(4) :: igrid
      real(8)    :: dv_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,dj_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,crhs(nx*ny*nz,3)

c Local variables

      integer(4) :: ii,ivar
      real(8)    :: a1,a2,a3,etal

c Externals

      real(8)    :: res
      external      res

c Begin program

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
            a1= curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .              ,dj_cov(:,:,:,1)
     .              ,dj_cov(:,:,:,2)
     .              ,dj_cov(:,:,:,3),1)
            a2= curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .              ,dj_cov(:,:,:,1)
     .              ,dj_cov(:,:,:,2)
     .              ,dj_cov(:,:,:,3),2)
            a3= curl(i,j,k,nx,ny,nz,igrid,igrid,igrid
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
