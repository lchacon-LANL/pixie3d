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

      implicit none

c Call variables

      integer(4) :: ntot,iout
      real(8)    :: x(ntot),y(ntot)

c Local variables

      real(8)    :: xxx(ntotdp,neqd)
      real(8)    :: yyy(ntotdp,neqd)

      integer(4) :: ii,iii,igrid

c Debug

      real(8)    :: debug(0:nx+1,0:ny+1,0:nz+1)

c Externals

      external   :: tmp_mtvc,rho_mtvc,b_mtvc

c Begin program

!This is known from setupNonlinearFunction
cc      nx = grid_params%nxv(igx)
cc      ny = grid_params%nyv(igy)
cc      nz = grid_params%nzv(igz)

      igrid = igx

c *******************************************************************
c     Identity preconditioner
c *******************************************************************

      if (precon.eq.'id') then

        x = y*dt   !Diagonal scaling

c *******************************************************************
c     Semi-implicit preconditioner
c *******************************************************************

      else

c     Set up size of operator split vectors

        do ieq=1,neqd
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                iii = ieq + neqd*(i-1 + nx*(j-1) + nx*ny*(k-1))

                call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                             ,cartsn)

                yyy(ii,ieq) = y(iii)*volume(i,j,k,igx,igy,igz)  !Volume weighing of residuals
                
cc                yyy(ii,ieq) = sin(2*pi*x1/xmax)*cos(2*pi*y1/ymax)
cc     .                       *dxh(ig)*dyh(jg)*dzh(kg)
cc                write (*,*) yyy(ii,ieq)

                xxx(ii,ieq) = 0d0
              enddo
            enddo
          enddo
        enddo

c     Predictor step

        !Temperature
        call cSolver(1,ntotdp,yyy(:,ITMP),xxx(:,ITMP),bcs(:,ITMP)
     .                ,igrid,iout,ITMP,tmp_mtvc,tmp_diag)

        !Density
        call cSolver(1,ntotdp,yyy(:,IRHO),xxx(:,IRHO),bcs(:,IRHO)
     .                ,igrid,iout,IRHO,rho_mtvc,rho_diag)

        !Magnetic field
        call cSolver(3,ntotdp,yyy(:,IBX:IBZ),xxx(:,IBX:IBZ)
     .                  ,bcs(:,IBX:IBZ),igrid,iout,IBX,b_mtvc,b_diag)

        open(unit=110,file='debug.bin',form='unformatted'
     .       ,status='replace')
        call mapMGVectorToArray(0,1,1,yyy(:,IBX),nx,ny,nz,debug,1)
        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,0,110)
        call mapMGVectorToArray(0,1,1,xxx(:,IBX),nx,ny,nz,debug,1)
        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
        call mapMGVectorToArray(0,1,1,yyy(:,IBY),nx,ny,nz,debug,1)
        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
        call mapMGVectorToArray(0,1,1,xxx(:,IBY),nx,ny,nz,debug,1)
        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
        call mapMGVectorToArray(0,1,1,yyy(:,IBZ),nx,ny,nz,debug,1)
        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
        call mapMGVectorToArray(0,1,1,xxx(:,IBZ),nx,ny,nz,debug,1)
        call contour(debug(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
        close(110)
        stop

c     SI step


c     Corrector step



c     Store solution in "x"

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              do ieq=1,neqd
                ii  = i + nx*(j-1) + nx*ny*(k-1)
                iii = ieq + neqd*(i-1 + nx*(j-1) + nx*ny*(k-1))
                x(iii) = xxx(ii,ieq)
              enddo
            enddo
          enddo
        enddo

      endif

c End program

      end subroutine applyPreconditioner

c scalarSolver
c #########################################################################
      subroutine scalarSolver(ntot,b,x,bcond,igrid,out,matvec)
c--------------------------------------------------------------------
c     This subroutine test the matrix-light solver, solving Ax=b,
c     where the action of A on a vector is given in routine test_mtvc.
c--------------------------------------------------------------------

      use precond_setup

      use mlsolverSetup

      implicit none

c Call variables

      integer(4) :: ntot,igrid,bcond(6,*),out
      real(8)    :: x(ntot),b(ntot)
      external   :: matvec

c Local variables

      integer(4) :: guess

c Begin program

c     Initialize solver

      call solverInit

c     Upper_level solver options

      call solverOptionsInit

      solverOptions%tol      = mgtol
      solverOptions%vcyc     = maxvcyc
      solverOptions%igridmin = 3
      solverOptions%orderres = 2
      solverOptions%orderprol= 2
      solverOptions%vol_res  = .true.
      call assembleSolverHierarchy('mg')

c     Next level solver

      call solverOptionsInit

      solverOptions%iter = nsweep
      solverOptions%omega= 0.7

      call assembleSolverHierarchy('jb')

c     Invoke solver

      call getSolver(1,ntot,b,x,matvec,igrid,bcond,guess,out,1)

c     Get output data

cc      call getSolverOptions(1)

c     Kill solver

      call solverKill

c End program

      end subroutine scalarSolver

c cSolver
c #########################################################################
      subroutine cSolver(neq,ntotp,b,x,bcond,igrid,out,icmp,matvec,diag)
c--------------------------------------------------------------------
c     This subroutine test the matrix-light solver, solving Ax=b,
c     where the action of A on a vector is given in routine test_mtvc.
c--------------------------------------------------------------------

      use precond_variables

      use mlsolverSetup

      implicit none

c Call variables

      integer(4) :: neq,ntotp,igrid,bcond(6,neq),out,icmp
      real(8)    :: x(ntotp,neq),b(ntotp,neq)
      real(8), target :: diag(neq,2*neq*ntotp)

      external   :: matvec

c Local variables

      integer(4) :: guess,ntot
      real(8)    :: xi(ntotp*neq),bi(ntotp*neq)

c Begin program

      icomp = icmp  !Define icomp for setMGBC

c Interlace variables for coupled solve

      do i=1,ntotp
        do ieq=1,neq
          xi(neq*(i-1)+ieq) = x(i,ieq)
          bi(neq*(i-1)+ieq) = b(i,ieq)
        enddo
      enddo

c Solve coupled MG

c     Initialize solver

      call solverInit

c     Upper_level solver options

      call solverOptionsInit

      solverOptions%tol      = mgtol
      solverOptions%vcyc     = maxvcyc
      solverOptions%igridmin = 3
      solverOptions%orderres = 2
      solverOptions%orderprol= 2
      solverOptions%vol_res  = .true.
cc      solverOptions%diag     => diag

      call assembleSolverHierarchy('mg')

c     Next level solver

      call solverOptionsInit

      solverOptions%iter = nsweep
      solverOptions%omega= 0.7

      call assembleSolverHierarchy('jb')

c     Invoke solver

      ntot=neq*ntotp
      call getSolver(neq,ntot,bi,xi,matvec,igrid,bcond,guess,out,1)

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

c     contour
c     #####################################################################
      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)
      implicit none               !For safe fortran
c     ---------------------------------------------------------------------
c     Contours arr in xdraw format
c     Notes:
c      put the next 2 lines in main
c      open(unit=nunit,file='contour.bin',form='unformatted') before
c      close(unit=nunit) after
c     ---------------------------------------------------------------------

c     Call variables

        integer*4      nx,ny,iopt,nunit
        real*8         arr(nx,ny),xmin,xmax,ymin,ymax

c     Local variables

        integer*4      i,j

c     Begin program

        if(iopt.eq.0) then
          write(nunit) nx-1,ny-1,0
          write(nunit) real(xmin,4),real(xmax,4)
     .                ,real(ymin,4),real(ymax,4) 
        endif
        write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)

      end subroutine contour
