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
     .                        ,order,bcond,volf)
c     -----------------------------------------------------------------
c     Restricts array to mgvector in all grids (without ghost nodes),
c     starting at grid igr0.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer(4) :: ieq,neq,nx,ny,nz,order,bcond(6)
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
     .               ,order,igridf,bcond,volf)

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

c test_mtvc
c####################################################################
      subroutine test_mtvc(gpos,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the Alfven SI system.
c--------------------------------------------------------------------

      use grid

      use grid_aliases

      use matvec

      use precond_variables

      use constants

      implicit none

c Call variables

      integer*4    ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: i,j,k,ig,jg,kg,isig,ieq,ijg,neq
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: xarr

      real(8)    :: upwind,nu2

c External

      real*8       lap_vec,vis
      external     lap_vec,vis

c Begin program

      isig = grid_params%istartp(igrid)

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      neq = ntot/nx/ny/nz

c Find limits for loops

      call limits(gpos,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(xarr(0:nx+1,0:ny+1,0:nz+1,neq))

      !Set grid=1 because x is NOT a MG vector
      do ieq=1,neq
        call mapMGVectorToArray(ieq,neq,x,nx,ny,nz,xarr(:,:,:,ieq),1)
      enddo

      call setMGBC(neq,nx,ny,nz,igrid,xarr,bcnd)

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(i,j,k
     .            ,min(igrid,ngrdx),min(igrid,ngrdy),min(igrid,ngrdz)
     .            ,ig,jg,kg)

            ijg   = i + nx*(j-1) + nx*ny*(k-1) + isig - 1

cc            upwind = .5*(vxx(ijg)+abs(vxx(ijg)))
cc     .                 *( array(i  ,j,k,1) - array(i-1,j,k,1) )/dx(ig-1)
cc     .              +.5*(vxx(ijg)-abs(vxx(ijg)))
cc     .                 *( array(i+1,j,k,1) - array(i  ,j,k,1) )/dx(ig)
cc     .              +.5*(vyy(ijg)+abs(vyy(ijg)))
cc     .                 *( array(i,j  ,k,1) - array(i,j-1,k,1) )/dy(jg-1)
cc     .              +.5*(vyy(ijg)-abs(vyy(ijg)))
cc     .                 *( array(i,j+1,k,1) - array(i,j  ,k,1) )/dy(jg)
cc     .              +.5*(vzz(ijg)+abs(vzz(ijg)))
cc     .                 *( array(i,j,k  ,1) - array(i,j,k-1,1) )/dz(kg-1)
cc     .              +.5*(vzz(ijg)-abs(vzz(ijg)))
cc     .                 *( array(i,j,k+1,1) - array(i,j,k  ,1) )/dz(kg)
cc
cc            
cc          nu2 = vis(i,j,nx,ny,igrid)
cc
cc          y(ij) = dx1*dy1*( x(ij)/dt + alpha*upwind 
cc     .                    -(1.-cnfactor*cnf_d)
cc     .                   *nu2*lap_vec(i,j,nx,ny,dx1,dy1,x,zeros,igrid) )

          enddo
        enddo
      enddo

c End program

      deallocate(xarr)

      end subroutine test_mtvc

ccc bihar_mtvc_cpl
ccc####################################################################
cc      subroutine bihar_mtvc_cpl(elem,ntot,x,y,igrid,bcond)
ccc--------------------------------------------------------------------
ccc     This subroutine calculates, for given x, y = A(psi)x  matrix-free
ccc     for the advection-hyperdiffusion coupled system.
ccc--------------------------------------------------------------------
cc
cc      use precond_variables
cc
cc      use mg_setup
cc
cc      use constants
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer*4    ntot,igrid,elem,bcond(4,2)
cc      real*8       x(ntot),y(ntot)
cc
ccc Local variables
cc
cc      integer*4    i,j,isig,nx,ny,neq,dum
cc      integer(4)   imin,imax,jmin,jmax
cc      integer(4)   iimin,iimax,jjmin,jjmax
cc
cc      integer*4    ii,im,ip,jj,jm,jp
cc      integer*4    ijg,ij
cc
cc      real(8)       dx1,dy1,bgrad2,upwind,lap,chyp
cc      real(8),allocatable,dimension(:,:) :: array1,array2
cc
ccc External
cc
cc      real*8       lap_vec,bgrad_vec,laplace,curr_vec,res
cc      external     lap_vec,bgrad_vec,laplace,curr_vec,res
cc
ccc Begin program
cc
cc      isig = istartp(igrid)
cc
cc      dx1 = dx(igrid)
cc      dy1 = dy(igrid)
cc
cc      nx = nxvp(igrid)
cc      ny = nyvp(igrid)
cc
cc      neq = ntot/nx/ny
cc
cc      allocate (array1(0:nx+1,0:ny+1))
cc      allocate (array2(0:nx+1,0:ny+1))
cc
cc      chyp = (1.-cnfactor*cnf_d)*ddhyper
cc
ccc Find limits for loops
cc
cc      call limits(elem,nx,ny,imin,imax,jmin,jmax)
cc
ccc Unravel quantities in separate vectors (x1 --> Bz, x2 --> xi)
cc
cc      jjmin = max(jmin-1,1)
cc      jjmax = min(jmax+1,ny)
cc      iimin = max(imin-1,1)
cc      iimax = min(imax+1,nx)
cc
cc      if (elem.ne.0) then
cc        do j=jjmin,jjmax
cc          do i=iimin,iimax
cc            ij = i + nx*(j-1)
cc            array1(i,j) = x(neq*(ij-1)+1)
cc            array2(i,j) = x(neq*(ij-1)+2)
cc          enddo
cc        enddo
cc      else
cc        do j=1,ny
cc          do i=1,nx
cc            ij = i + nx*(j-1)
cc            array1(i,j) = x(neq*(ij-1)+1)
cc            array2(i,j) = x(neq*(ij-1)+2)
cc          enddo
cc        enddo
cc      endif
cc
ccc Set boundary conditions in ghost nodes
cc
cc      call setBoundaryConditions(array1,iimin,iimax,jjmin,jjmax,nx,ny
cc     .                          ,bcond(:,1))
cc      call setBoundaryConditions(array2,iimin,iimax,jjmin,jjmax,nx,ny
cc     .                          ,bcond(:,2))
cc
ccc Calculate matrix-vector product
cc
cc      do j = jmin,jmax
cc        do i = imin,imax
cc
cc          call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond(1,1))
cc
cc          ij   = ii + nx*(jj-1)
cc          ijg  = ij  + isig - 1
cc
cc          upwind = .5*(vxx(ijg)+abs(vxx(ijg)))
cc     .                 *( array1(ii,jj)-array1(im,jj) )/dx1
cc     .            +.5*(vxx(ijg)-abs(vxx(ijg)))
cc     .                 *( array1(ip,jj)-array1(ii,jj) )/dx1
cc     .            +.5*(vyy(ijg)+abs(vyy(ijg)))
cc     .                 *( array1(ii,jj)-array1(ii,jm) )/dy1
cc     .            +.5*(vyy(ijg)-abs(vyy(ijg)))
cc     .                 *( array1(ii,jp)-array1(ii,jj) )/dy1
cc
cc          lap = laplace(i,j,nx,ny,dx1,dy1,array1)
cc
cc          y(neq*(ij-1)+1) = dx1*dy1*( array1(ii,jj)/dt + alpha*upwind 
cc     .                    -(1.-cnfactor*cnf_d)*res(i,j,nx,ny,igrid)*lap
cc     .                    +laplace(i,j,nx,ny,dx1,dy1,array2) )
cc
cc          y(neq*(ij-1)+2) = dx1*dy1*( array2(ii,jj) - chyp*lap )
cc        enddo
cc      enddo
cc
ccc End program
cc
cc      deallocate(array1,array2)
cc
cc      return
cc      end
