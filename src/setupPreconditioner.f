c setupPreconditioner
c####################################################################
      subroutine setupPreconditioner(ntot,x)

c--------------------------------------------------------------------
c     Finds velocity and magnetic field components in all grids
c     for matrix-free preconditioning.
c--------------------------------------------------------------------

      use grid

      use precond_variables

      use matvec

      use operators

      use newton_gmres

      implicit none

c Call variables

      integer(4) :: ntot
      real(8)    :: x(ntot)

c Local variables

      integer(4) :: ncolors

c Debug

      real(8)    :: mag,debug(0:nx+1,0:ny+1,0:nz+1)
      integer(4) :: nt,ii
      real(8),allocatable,dimension(:,:) :: v_mat,debug2

c Externals

      external   :: v_mtvc,b_mtvc,rho_mtvc,tmp_mtvc

c Begin program

c diag ******* NUMERICAL GRID QUANTITIES
cc      open(unit=110,file='debug.bin',form='unformatted'
cc     .    ,status='replace')
cc      igx = 1
cc      do i=1,3
cc        do j=1,3
cc          nx = grid_params%nxv(igx)
cc          ny = grid_params%nyv(igx)
cc          nz = grid_params%nzv(igx)
cc          debug = gmetric%grid(igx)%cnv(0:nx+1,0:ny+1,0:nz+1,i,j)
cccc          debug = gmetric%grid(igx)%Gamma(0:nx+1,0:ny+1,0:nz+1,1,i,j)
cccc     .           *gmetric%grid(igx)%jac
cc
cc          call contour(debug(0:nx+1,0:ny+1,1),nx+2,ny+2,0d0
cc     .                ,xmax,0d0,ymax,i+j-2,110)
cc        enddo
cc      enddo
cc      close(110)
cc      stop
c diag *******

      if (precon == 'id') return

      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Unpack vector x

      varray = x  !This allocates varray; overloaded assignment

      call imposeBoundaryConditions(varray,igx,igy,igz)

c Extract arrays and BC's

      rho => varray%array_var(IRHO)%array
      rvx => varray%array_var(IVX )%array
      rvy => varray%array_var(IVY )%array
      rvz => varray%array_var(IVZ )%array
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
      tmp => varray%array_var(ITMP)%array

      do ieq=1,neqd
        bcs(:,ieq) = varray%array_var(ieq)%bconds(:)
      enddo

      !Current boundary conditions
      bcs(:,IJX:IJZ) = bcs(:,IBX:IBZ)
      where (bcs(:,IJX:IJZ) == -NEU)
        bcs(:,IJX:IJZ) = -DIR  !Use covariant components for tangential dirichlet
      end where

c Find coefficients for linearized systems

      if (jit == 1) call findCoeffs
cc      call findCoeffs

c Find required diagonals in all grids

      ncolors = 4

      if (jit == 1) then
        call find_mf_diag_colored(1,ntotdp,rho_mtvc,1,bcs(:,IRHO)
     .                           ,rho_diag,ncolors)

        call find_mf_diag_colored(1,ntotdp,tmp_mtvc,1,bcs(:,ITMP)
     .                           ,tmp_diag,ncolors)

        call find_mf_diag_colored(3,3*ntotdp,b_mtvc,1,bcs(:,IBX:IBZ)
     .                           ,b_diag,ncolors)

        if (.not.gm_smooth) then
          call find_mf_diag_colored(3,3*ntotdp,v_mtvc,1,bcs(:,IVX:IVZ)
     .                             ,v_diag,ncolors)
        endif

        form_diag = .false.
      endif

c diag ***** DIAGONAL GENERATION FOR XDRAW PLOTTING
cc      do k=1,nzz
cc        do j=1,nyy
cc          do i=1,nxx
cc            ii = i + nxx*(j-1)+nxx*nyy*(k-1)
cc            debug(i,j,k) = 1./sqrt(sum(v_diag(3,3*ii-2:3*ii)**2)
cc          enddo
cc        enddo
cc      enddo
cc
cc      open(unit=110,file='debug.bin',form='unformatted'
cc     .    ,status='replace')
cc      call contour(debug(1:nxx,1:nyy,1),nxx,nyy,0d0,xmax,0d0,ymax,0,110)
cc      close(110)
cc      stop
c diag *****

c diag ***** MATRIX GENERATION FOR MATLAB PROCESSING
cc      igrid = 1
cc
cc      nxx = grid_params%nxv(igrid)
cc      nyy = grid_params%nyv(igrid)
cc      nzz = grid_params%nzv(igrid)
cc
cc      nt = 3*nxx*nyy*nzz
cc
cc      allocate(v_mat(nt,nt))
cc
cc      icomp = IVX
cc      call find_mf_mat(3,nt,v_mtvc,igrid,bcs(:,IVX:IVZ),v_mat)
cc
cccc      open(unit=110,file='debug.bin',form='unformatted'
cccc     .       ,status='replace')
cccc      call contour(v_mat,nt,nt,0d0,1d0,0d0,1d0,0,110)
cccc      close(110)
cc
cc      open(unit=110,file='debug.mat',status='replace')
cc      do i=1,nt
cc        write (110,*) (v_mat(i,j),j=1,nt)
cc      enddo
cc      close(110)
cc      stop
c diag *****

c End program

      end subroutine setupPreconditioner

c killPreconditioner
c###################################################################
      subroutine killPreconditioner

      use precond_variables

c     Deallocate variables

      if (precon == 'id') return

cc      call deallocPrecVariables

      end subroutine killPreconditioner
