c module precond_setup
c ######################################################################
      module precond_setup

        integer(4)    ::  precpass,nsweep,maxvcyc,ndiagdp
        parameter (ndiagdp=3)

        character*(10)::  precon

      end module precond_setup

c  module matvec
c ###################################################################
      module matvec

        use grid

        use mg_internal

      contains

c     limits
c     ###############################################################
      subroutine limits(elem,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax)
      implicit none
c     ---------------------------------------------------------------
c     Finds limits on loops for matvec routines. Used in finding 
c     diagonal from matvec.
c     ---------------------------------------------------------------

c     Call variables

      integer(4) :: elem,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax

c     Local variables

      integer(4) :: el1

c     Begin program

      if (elem.eq.0) then
        imin = 1
        imax = nx
        jmin = 1
        jmax = ny
        kmin = 1
        kmax = nz
      else
        el1  = mod(elem,nx*ny)
        if (el1 == 0) el1 = nx*ny

        imin = mod(el1 ,nx)
        if (imin == 0) imin = nx
        imax = imin

        jmin = 1 + (el1 - imin)/nx
        jmax = jmin

        kmin = 1 + (elem - imin - nx*(jmin-1))/nx*ny
        kmax = kmin
      endif

c     End program

      end subroutine limits

c     restrictArray
c     #################################################################
      subroutine restrictArray(ieq,neq,nx,ny,nz,array,mgvector,igr0
     .                        ,order,volf)
c     -----------------------------------------------------------------
c     Restricts array to mgvector in all grids (without ghost nodes),
c     starting at grid igr0.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer(4) :: ieq,neq,nx,ny,nz,order
      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1)
      real(8)    :: mgvector(*)
      logical    :: volf

c     Local variables

      integer(4) :: nxf,nyf,nzf,nxc,nyc,nzc,igridc,igridf,igr0
     .             ,isigf,isigc

c     Begin program

c     Consistency check

      nxf = grid_params%nxv(igr0)
      nyf = grid_params%nyv(igr0)
      nzf = grid_params%nzv(igr0)

      if (nxf /= nx .or. nyf /= ny .or. nzf /= nz) then
        write (*,*) 'Grid mismatch in restrictArray:'
        write (*,*) 'Aborting...'
        stop
      endif

c     Map array in initial grid onto MG vector

      call mapArrayToMGVector(ieq,neq,nx,ny,nz,array,mgvector,igr0)

c     Restrict array to coarser grids

      do igridc = igr0+1,grid_params%ngrid

        igridf = igridc-1

c       Characterize coarse and fine grids

        nxf = grid_params%nxv(igridf)
        nyf = grid_params%nyv(igridf)
        nzf = grid_params%nzv(igridf)

        nxc = grid_params%nxv(igridc)
        nyc = grid_params%nyv(igridc)
        nzc = grid_params%nzv(igridc)

        isigc = grid_params%istartp(igridc)
        isigf = grid_params%istartp(igridf)

c       Restrict MG vector

        call restrict(mgvector(isigc),nxc,nyc,nzc
     .               ,mgvector(isigf),nxf,nyf,nzf
     .               ,order,igridf,volf)

      enddo

      end subroutine restrictArray

c     symm_test
c     ###############################################################
      subroutine symm_test(neq,igrid,matvec,bcnd)
c     ---------------------------------------------------------------
c     Performs symmetry test of matvec on grid "igrid".
c     ---------------------------------------------------------------

      implicit none     !For safe fortran

c     Call variables

      integer(4) :: neq,igrid,bcnd(6,neq)

      external      matvec

c     Local variables

      real(8),allocatable,dimension(:)::x1,dummy,dummy2

      real(8)    :: dd1,dd2,error
      integer(4) :: nx,ny,nz,nn,ii,jj,i1,j1,i2,j2,ix1,iy1,ix2,iy2,ieq

c     Begin program

c     Initialize variables

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      nn = neq*nx*ny*nz

      write (*,*) 'Performing symmetry test of system matrix ',
     .            'on grid:',nx,'x',ny,'x',nz,'...'

      allocate(x1(nn),dummy(nn),dummy2(nn))

      error = 0d0

      do ii = 1,nn
        x1    (ii) = 0d0
        dummy (ii) = 0d0
        dummy2(ii) = 0d0
      enddo

c     Check symmetry

      do ii = 1,nn/neq

        do ieq =1,neq

c       Find column vector ii

          call findBaseVector(ii,ieq,neq,nn,x1,1d0)

          call matvec(0,nn,x1,dummy,igrid,bcnd)

          call findBaseVector(ii,ieq,neq,nn,x1,0d0)

c       Compare column vector ii with corresponding row vector (intersect in
c       diagonal)

          do jj = ii,nn/neq

            call findBaseVector(jj,ieq,neq,nn,x1,1d0)

            call matvec(ii,nn,x1,dummy2,igrid,bcnd)

            call findBaseVector(jj,ieq,neq,nn,x1,0d0)

            dd1 = abs(dummy(jj) - dummy2(ii))
            if(abs(dummy(jj)).gt.1d-15.or.abs(dummy2(ii)).gt.1d-15) then
              write(*,15) jj,ii,dummy(jj),ii,jj,dummy2(ii),dd1
     .             ,100*dd1/max(abs(dummy(jj)),abs(dummy2(ii)))
              error = error + dd1
            endif

          enddo

        enddo
      enddo

      write (*,20) error

      stop

c     End program

      deallocate (x1,dummy,dummy2)

 15   format ('(',i3,',',i3,'):',1pe10.2,'; (',i3,',',i3,'):',e10.2,
     .        '  Error:',e10.2,'  %error:',0pf7.2)
 20   format (/,'Total relative error:',1pe10.3)

      end subroutine symm_test

ccc shiftIndices
ccc####################################################################
cc      subroutine shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,isig,bcnd)
cc      implicit none
ccc--------------------------------------------------------------------
ccc     Shifts indices in a 5pt stencil for both arrays (isig = 0) 
ccc     and MG vectors (isig = pointer to grid), accounting for boundary 
ccc     conditions as defined in integer array bcnd.
ccc--------------------------------------------------------------------
cc
ccc Call variables
cc
cc      integer :: i,j,ii,jj,im,ip,jm,jp,nx,ny,isig,bcnd(6)
cc
ccc Local variables
cc
cc      integer :: nxx,row
cc
ccc Functions
cc
cc      row(i,j) = i + nx*(j-1) + isig - 1
cc
ccc Begin program
cc
cc      if (bcnd(3) == 0 .or. bcnd(4) == 0) then
cc        call    periodicBC(i,nx,ii,im,ip)
cc      else
cc        call nonperiodicBC(i,nx,ii,im,ip,bcnd(3),bcnd(4))
cc      endif
cc
cc      if (bcnd(1) == 0 .or. bcnd(2) == 0) then
cc        call    periodicBC(j,ny,jj,jm,jp)
cc      else
cc        call nonperiodicBC(j,ny,jj,jm,jp,bcnd(1),bcnd(2))
cc      endif
cc
cc      if (isig.gt.0) then    !Transform to MG vector coordinates
cc        ip = min(row(ip,jj),row(nx,ny))
cc        im = max(row(im,jj),row(1 , 1))
cc        jp = min(row(ii,jp),row(nx,ny))
cc        jm = max(row(ii,jm),row(1 , 1))
cc        ii = row(ii,jj)
cc        jj = ii
cc      endif
cc
ccc End program
cc
cc      contains
cc
cc      subroutine periodicBC(i,nx,ii,im,ip)
cc
cc      implicit none
cc
cc      integer :: i,ii,im,ip,nx
cc
cc      if (i.lt.1) then
cc        ii = nx + i - 1 
cc      elseif (i.gt.nx) then
cc        ii = i - nx + 1 
cc      else
cc        ii = i 
cc      endif
cc
cc      if (i.eq.nx) then
cc        ip = 2 
cc      else
cc        ip = ii + 1
cc      endif
cc
cc      if (i.eq.1 ) then
cc        im = nx - 1 
cc      else
cc        im = ii - 1
cc      endif
cc
cc      end subroutine
cc
cc      subroutine nonperiodicBC(i,nx,ii,im,ip,bcs1,bcs2)
cc
cc      implicit none
cc
cc      integer :: i,ii,im,ip,nx,bcs1,bcs2
cc
cc      ii = i
cc      if (i == 1 .and. bcs1 == 2) then
cccc        im = i
cc        im = i+1
cc      else 
cc        im = i-1
cc      endif
cc
cc      if (i == nx .and. bcs2 == 2) then
cccc        ip = i
cc        ip = i-1
cc      else
cc        ip = i+1
cc      endif
cc
cc      end subroutine
cc
cc      end subroutine

      end module matvec

c module precond_variables
c ######################################################################
      module precond_variables

        use timeStepping

        use precond_setup

        use transport_params

        use equilibrium

        use parameters

        use constants

        use grid_aliases

        use auxiliaryVariables

        integer(4) :: icomp

        integer(4) :: i,j,k,ig,jg,kg,ieq
        real(8)    :: x1,y1,z1,jac
        logical    :: cartsn,covariant,to_cartsn,to_cnv

        real(8), allocatable, dimension(:,:) :: bcnv,vcnv

        real(8), allocatable, dimension(:) :: divV

cc        real(8), allocatable, dimension(:,:,:) :: pp

        integer(4), allocatable, dimension(:,:) :: bcs

cc        real(8), allocatable, dimension(:):: diag_mu
cc        real(8),allocatable,dimension(:,:):: d_sc,d_bc,d_pc,d_bh,d_alf
cc
cc        logical :: formdiag
cc
cc        integer :: iiout
cc
cc        integer :: bc_sc(4,1),bc_cpld(4,2)

      end module precond_variables
