      subroutine order_tst

c##########################################################################
c     Tests spatial convergence
c##########################################################################

      use parameters

      use variables

      use grid

      use iosetup

      use equilibrium

      implicit none

c Local variables

      integer(4)     :: i,j,k,ig,jg,kg,itime
      integer(4)     :: ierr,system,nplot,igx,igy,igz
      character*(40) :: command

      real(8)        :: mag,time,dt
      real(8),allocatable,dimension(:,:,:,:) :: pos 
      real(8),pointer,dimension(:,:,:) :: array

      type (var_array),target :: v256,vcur

c Interpolation

      real(8)    :: xp,yp,zp,ff

      real(8),allocatable,dimension(:) :: xx,yy,zz

      integer(4) ::  kx,ky,kz,nx,ny,nz,dim,flg,order
      real(8), dimension(:),allocatable:: tx,ty,tz,work,q
      real(8), dimension(:,:,:),allocatable:: bcoef

      real(8)    :: db3val
      external      db3val

c Begin program

      write (*,*) 'Start order test'

      igx = 1
      igy = 1
      igz = 1

      call deallocateGridStructure(grid_params)

c Read reference grid info

      open(unit=urecord,file='record-256x256.bin'
     .    ,form='unformatted',status='old')

      read (urecord) nxd
      read (urecord) nyd
      read (urecord) nzd

      !Create grid
      call createGrid(nxd,nyd,nzd)

      !Read solutions
      call allocateDerivedType(v256)

      write (*,*) ' Reading reference file...'
      write (*,*) ' Grid: ',nxd,'x',nyd,'x',nzd

      do
        call readRecord(urecord,itime,time,dt,v256,ierr)

        if (ierr /= 0) exit
      enddo

      close (urecord)

c Spline reference solution

      order = 2

      allocate(xx(nxd+2),yy(nxd+2),zz(nxd+2))

      call getMGmap(1,1,1,igx,igy,igz,ig,jg,kg)

      xx(1:nxd+2) = grid_params%xx(ig-1:ig+nxd)
      yy(1:nyd+2) = grid_params%yy(jg-1:jg+nyd)
      zz(1:nzd+2) = grid_params%zz(kg-1:kg+nzd)

      flg = 0
      nx = nxd+2
      ny = nyd+2
      nz = nzd+2
      kx = min(order+1,nx-1)
      ky = min(order+1,ny-1)
      kz = min(order+1,nz-1)

      dim = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))

      allocate(tx(nx+kx))
      allocate(ty(ny+ky))
      allocate(tz(nz+kz))
      allocate(work(dim))
      allocate(bcoef(nx,ny,nz))

      array => v256%array_var(IVX)%array

      call db3ink(xx,nx,yy,ny,zz,nz,array(0:nxd+1,0:nyd+1,0:nzd+1)
     .           ,nx,ny,kx,ky,kz,tx,ty,tz,bcoef,work,flg)

      call deallocateGridStructure(grid_params)

c Initialize current grid info and read current grid solution

      open(urecord,file=recordfile,form='unformatted',status='unknown')

      read (urecord) nxd
      read (urecord) nyd
      read (urecord) nzd

      call setVectorDimensions

      call createGrid(nxd,nyd,nzd)

      call allocateDerivedType(vcur)

      write (*,*) ' Reading solution file...'
      write (*,*) ' Grid: ',nxd,'x',nyd,'x',nzd

      do
        call readRecord(urecord,itime,time,dt,vcur,ierr)

        if (ierr /= 0) exit
      enddo

      close (urecord)

      array => vcur%array_var(IVX)%array

c Calculate difference

      mag = 0d0

      do k = 1,nzd
        do j = 1,nyd
          do i = 1,nxd
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            xp = grid_params%xx(ig)
            yp = grid_params%yy(jg)
            zp = grid_params%zz(kg)

            mag = mag +
     .            (array(i,j,k)
     .            -db3val(xp,yp,zp,0,0,0,tx,ty,tz,nx,ny,nz
     .                   ,kx,ky,kz,bcoef,work))**2

          enddo
        enddo
      enddo

      deallocate(tx,ty,tz,work,bcoef,xx,yy,zz)

      mag = sqrt(mag/nxd/nyd/nzd)

      write (*,*) ' L2-norm error:',mag

c End program

      end subroutine order_tst
