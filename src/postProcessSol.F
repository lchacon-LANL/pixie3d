c postProcessSol
c####################################################################
      subroutine postProcessSol(v_arr)

c--------------------------------------------------------------------
c     Postprocess solution (e.g. filtering). In this instance, we
c     filter with a few (1-2) passes of JB on the system:
c                 lap(psi*) = psi ; lap(bz*) = bz
c     where star quantities are filtered ones. We implement this
c     in terms of increments with respect to the equilibrium.
c
c     In parallel, the routine should perform communication to
c     provide adequate BCs to v_arr for subsequent postprocessing.
c--------------------------------------------------------------------

cc      use matvec
cc
cc      use nlfunction_setup
cc
cc      use imposeBCinterface
cc
cc      use local_BCS_variables, ONLY:order_bc
cc
cc      implicit none
cc
ccc Call variables
cc
cc      type(var_array) :: v_arr
cc
ccc Local variables
cc
cc      integer(4) :: ii,nn,guess,iout,bcpsi(6),bcbz(6),ncolors
cccc     .             ,filter_it
cc      real(8)    :: dvol,opt_om
cc      real(8),allocatable,dimension(:,:) :: xpsi,xbz
cc     .                                     ,xrhs_bz,xrhs_psi
cc
cc      real(8),pointer,dimension(:,:,:) :: dpsi,dbz,dpsio,dbzo
cc
ccc Externals
cc
cc      external    hz_mtvc
cc
ccc Begin program
cc
cc      if (di == 0d0) return
cc
cc      !Setup grid info
cc      igx = 1
cc      igy = 1
cc      igz = 1
cc
cc      nx = grid_params%nxv(igx)
cc      ny = grid_params%nyv(igy)
cc      nz = grid_params%nzv(igz)
cc
cc      nn = nx*ny*nz
cc
ccc Setup JB filter 
cc
cc      !Number of colors on grid
cc      ncolors = 2
cc
cc      !Optimal JB damping parameter (to remove Nyquist mode on grid)
cc      opt_om = dt*di*(((xmax-xmin)/nxd)**(-2)
cc     .               +((ymax-ymin)/nyd)**(-2)
cc     .               +((zmax-zmin)/nzd)**(-2))
cc
cc      opt_om = (1+2*opt_om)/(1+4*opt_om)
cc
ccc Enable communication for parallel BC
cc
cc      call setASMflag(.false.)
cc
ccc Copy solution state info
cc
cc      dpsi  => v_arr%array_var(IPSI)%array
cc      dpsio => u_0  %array_var(IPSI)%array
cc      bcpsi =  v_arr%array_var(IPSI)%bconds
cc
cc      dbz  => v_arr%array_var(IBZ)%array
cc      dbzo => u_0  %array_var(IBZ)%array
cc      bcbz =  v_arr%array_var(IBZ)%bconds
cc
ccc Form required diagonals (only once)
cc
cc      if(.not.allocated(d_hz_psi)) then
cc        allocate(d_hz_psi(1,2*nn),d_hz_bz(1,2*nn))
cc        call find_mf_diag_colored(1,nn,hz_mtvc,1,bcpsi,d_hz_psi,ncolors)
cc        call find_mf_diag_colored(1,nn,hz_mtvc,1,bcbz ,d_hz_bz ,ncolors)
cc      endif
cc
ccc Define solution vectors and rhs
ccc (defined in terms of increments to have homogeneous Dirichlet BCs)
cc
cc      allocate(xpsi(nn,1),xrhs_psi(nn,1))
cc      allocate(xbz (nn,1),xrhs_bz (nn,1))
cc
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc            if (vol_wgt) then
cc              dvol = volume(i,j,k,igx,igy,igz)
cc            else
cc              dvol = 1d0
cc            endif
cc
cc            xpsi    (ii,1) = dpsi(i,j,k)-dpsio(i,j,k)
cc            xbz     (ii,1) = dbz (i,j,k)-dbzo (i,j,k)
cc
cc            xrhs_psi(ii,1) = xpsi(ii,1)*dvol
cc            xrhs_bz (ii,1) = xbz (ii,1)*dvol
cc          enddo
cc        enddo
cc      enddo
cc
ccc Postprocess relevant quantities
cc
cc      iout = ilevel - 1
cc      guess = 1
cc
cc      if (ilevel > 0 .and. my_rank == 0) then
cccc        if (iout > 0 .and. my_rank == 0) write (*,*)
cc        write (*,*) '>>>>> Psi smoothing...'
cc        if (iout > 0 .and. my_rank == 0) write (*,*)
cc      endif
cc
cc      !Psi
cc      call cSolver(1,nn,xrhs_psi,xpsi,bcpsi
cc     .            ,igx,iout,guess,hz_mtvc,dg=d_hz_psi
cc     .            ,ncolors=ncolors,sm_it=filter_it
cc     .            ,igmin=grid_params%ngrid+1
cc     .            ,smoother='jb',omega=opt_om)
cc
cc      !Bz
cc      if (ilevel > 0 .and. my_rank == 0) then
cc        if (iout > 0 .and. my_rank == 0) write (*,*)
cc        write (*,*) '>>>>> Bz  smoothing...'
cc        if (iout > 0 .and. my_rank == 0) write (*,*)
cc      endif
cc      call cSolver(1,nn,xrhs_bz,xbz,bcbz
cc     .            ,igx,iout,guess,hz_mtvc,dg=d_hz_bz
cc     .            ,ncolors=ncolors,sm_it=filter_it
cc     .            ,igmin=grid_params%ngrid+1
cc     .            ,smoother='jb',omega=opt_om)
cc
cc      if (ilevel > 0 .and. my_rank == 0) then
cc        if (iout > 0 .and. my_rank == 0) write (*,*)
cc        write (*,*) '>>>>> Done!'
cc      endif
cc
ccc Gather relevant quantities
cc
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            ii  = i + nx*(j-1) + nx*ny*(k-1)
cc            dpsi(i,j,k) = xpsi(ii,1)+dpsio(i,j,k)
cc            dbz (i,j,k) = xbz (ii,1)+dbzo (i,j,k)
cc          enddo
cc        enddo
cc      enddo
cc
ccc Set BCs on postprocessed arrays
ccc (communication required to update PETSc ghost cells)
cc
cc      if (np > 1) then
cc        ieq = IPSI
cc        call setBC(ieq,nx,ny,nz
cc     .                 ,v_arr%array_var(ieq)%array
cc     .                 ,u_0  %array_var(ieq)%array
cc     .                 ,v_arr%array_var(ieq)%bconds
cc     .                 ,igx,igy,igz,iorder=order_bc)
cc
cc        ieq = IBZ
cc        call setBC(ieq,nx,ny,nz
cc     .               ,v_arr%array_var(ieq)%array
cc     .               ,u_0  %array_var(ieq)%array
cc     .               ,v_arr%array_var(ieq)%bconds
cc     .               ,igx,igy,igz,iorder=order_bc)
cc      endif
cc
ccc End program
cc
cc      deallocate(xpsi,xrhs_psi)
cc      deallocate(xbz ,xrhs_bz )

      end subroutine postProcessSol
