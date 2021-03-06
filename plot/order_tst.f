      subroutine order_tst

c##########################################################################
c     Tests spatial convergence
c##########################################################################

      use parameters

      use variables

      use grid

      use var_io

      use ts_setup

      use equilibrium

      use slatec_splines

      implicit none

c Local variables

      integer        :: i,j,k,ig,jg,kg,nx,ny,nz
      integer        :: ierr,nplot,igx,igy,igz

      real(8)        :: mag,mag0,jac
      real(8),allocatable,dimension(:,:,:,:) :: pos 
      real(8),pointer,dimension(:,:,:) :: array

      type(var_array),pointer :: vref,vsol

c Interpolation

      real(8)    :: xp,yp,zp,interp
      real(8),allocatable,dimension(:) :: xx,yy,zz

cc      integer    :: kx,ky,kz,nx,ny,nz,dim,flg,order
cc
cc      real(8), dimension(:),allocatable:: tx,ty,tz,work
cc      real(8), dimension(:,:,:),allocatable:: bcoef
cc
cc      real(8)    :: db2val,db3val
cc      external      db2val,db3val

c Begin program

      write (*,*)
      write (*,*) ' Starting convergence test...'

      igx = 1
      igy = 1
      igz = 1

      call deallocateGridStructure(gv%gparams)
      call deallocateGridMetric(gv%gparams%gmetric)

c Read reference grid info

      open(unit=urecord,file='record-256x256.bin'
     .    ,form='unformatted',status='old')

      read (urecord) nxd,ilog,ihig
      read (urecord) nyd,jlog,jhig
      read (urecord) nzd,klog,khig

      write (*,*)
      write (*,*) ' Reading reference file...'
      write (*,*) ' Grid: ',nxd,'x',nyd,'x',nzd

      call setVectorDimensions

      !Create grid
cc      call createGrid(1,nxd,1,nyd,1,nzd,nxd,nyd,nzd,g_pack)

      !Read solutions
      call allocateDerivedType(vref)


      do
        ierr = readRecordFile(urecord,itime,time,dt,gammat,vref)

        if (ierr /= 0) exit
      enddo

      close (urecord)

      array => vref%array_var(IVX)%array

c Spline reference solution

      order = 2

      nx = nxd+2
      ny = nyd+2
      nz = nzd+2

      allocate(xx(nx),yy(ny),zz(nz))

      call getMGmap(gv%gparams,1,1,1,igx,igy,igz,ig,jg,kg)

      xx(1:nx) = gv%gparams%xx(ig-1:ig+nxd)
      yy(1:ny) = gv%gparams%yy(jg-1:jg+nyd)
      zz(1:nz) = gv%gparams%zz(kg-1:kg+nzd)

      flg = 0

      kx = min(order+1,nx-1)
      ky = min(order+1,ny-1)
      kz = min(order+1,nz-1)

      dim = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))

      allocate(tx(nx+kx))
      allocate(ty(ny+ky))
      allocate(tz(nz+kz))
      allocate(work(dim))
      allocate(bcoef(nx,ny,nz))

      call db3ink(xx,nx,yy,ny,zz,nz,array(0:nxd+1,0:nyd+1,0:nzd+1)
     .           ,nx,ny,kx,ky,kz,tx,ty,tz,bcoef,work,flg)

      call deallocateGridMetric(gv%gparams%gmetric)
      call deallocateGridStructure(gv%gparams)

c Initialize current grid info and read current grid solution

      igx = 1
      igy = 1
      igz = 1

      open(urecord,file=recordfile,form='unformatted',status='unknown')

      read (urecord) nxd,ilog,ihig
      read (urecord) nyd,jlog,jhig
      read (urecord) nzd,klog,khig

      write (*,*)
      write (*,*) ' Reading solution file...'
      write (*,*) ' Grid: ',nxd,'x',nyd,'x',nzd

      call setVectorDimensions

cc      call createGrid(1,nxd,1,nyd,1,nzd,nxd,nyd,nzd,g_pack)

      call allocateDerivedType(vsol)

      do
        ierr = readRecordFile(urecord,itime,time,dt,gammat,vsol)

        if (ierr /= 0) exit
      enddo

      close (urecord)

      array => vsol%array_var(IVX)%array

c Calculate difference

      mag0= 0d0
      mag = 0d0

      do k = 1,nzd
        do j = 1,nyd
          do i = 1,nxd
            call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)

            xp = gv%gparams%xx(ig)
            yp = gv%gparams%yy(jg)
            zp = gv%gparams%zz(kg)

            jac = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

            interp = db3val(xp,yp,zp,0,0,0,tx,ty,tz,nx,ny,nz
     .                     ,kx,ky,kz,bcoef,work)

cc            mag0 = mag0 + jac* array(i,j,k)**2
            mag0 = mag0 + jac* interp**2
            mag  = mag  + jac*(array(i,j,k)-interp)**2

          enddo
        enddo
      enddo

      mag = sqrt(mag/mag0)

      write (*,*)
      write (*,*) ' L2-norm relative error:',mag

c End program

      deallocate(tx,ty,tz,work,bcoef,xx,yy,zz)

      end subroutine order_tst
