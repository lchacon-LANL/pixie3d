c setMGBC
c####################################################################
      subroutine setMGBC(gpos,neq,nnx,nny,nnz,iig,array,bcnd)
c--------------------------------------------------------------------
c     Interfaces BC routines with MG code, for preconditioning.
c--------------------------------------------------------------------

      use imposeBCinterface

      use precond_variables

      use constants

      use mg_internal

      implicit none

c Call variables

      integer(4) :: nnx,nny,nnz,neq,bcnd(6,neq),iig,gpos

      real(8)    :: array(0:nnx+1,0:nny+1,0:nnz+1,neq)

c Local variables

      integer(4) :: offset

c Begin program

      igx = iig
      igy = iig
      igz = iig
      
      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      if (nx /= nnx .or. ny /= nny .or. nz /= nnz) then
        write (*,*) 'Grid sizes do not agree in setMGBC'
        write (*,*) 'Aborting...'
        stop
      endif

c Find loop limits for BC (passed to BC routines via imposeBCinterface module)

      call limits(gpos,nx,ny,nz,iimin,iimax,jjmin,jjmax,kkmin,kkmax)

      offset = 1

      iimin = max(iimin-offset,1)
      iimax = min(iimax+offset,nx)
      jjmin = max(jjmin-offset,1)
      jjmax = min(jjmax+offset,ny)
      kkmin = max(kkmin-offset,1)
      kkmax = min(kkmax+offset,nz)

c Select operation

      select case (neq)
      case(1)

        call setBC(IRHO,array(:,:,:,neq),zeros,bcnd(:,1))

cc        call setBC(icomp,array(:,:,:,neq),zeros,bcnd(:,1)) !Velocities for T BC unknown

      case(3)

        allocate(v_cov (0:nx+1,0:ny+1,0:nz+1,neq))

        select case (icomp)  !'icomp' is passed by the module precond_variables
        case(IVX,IVY,IVZ)
          call setBC(IVX,array,v_cov,vzeros,bcnd)
        case(IBX,IBY,IBZ)
          call setBC(IBX,array,v_cov,vzeros,bcnd)
        case(IJX,IJY,IJZ)
          call setBC(IJX,array,v_cov,vzeros,bcnd)
        case default
          write (*,*) 'icomp undefined in setMGBC'
          write (*,*) 'Aborting...'
          stop
        end select

        deallocate(v_cov)

      case default
        write (*,*) 'Number of equations not implemented in setMGBC'
        write (*,*) 'Aborting...'
        stop
      end select

c End

      end subroutine setMGBC

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

      use operators

      use mlsolverSetup

      implicit none

c Call variables

      integer*4    ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,neq
      integer(4) ::  imin, imax, jmin, jmax, kmin, kmax
      integer(4) :: iimin,iimax,jjmin,jjmax,kkmin,kkmax

      real(8),allocatable,dimension(:,:,:,:) :: xarr

      real(8)    :: upwind,nu2

c External

      real*8       lap_vec,vis
      external     lap_vec,vis

c Begin program

      isig = grid_params%istartp(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      neq = ntot/nx/ny/nz

c Find limits for loops

      call limits(gpos,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(xarr(0:nx+1,0:ny+1,0:nz+1,neq))

      !Set grid=1 because x is NOT a MG vector
      do ieq=1,neq
        call mapMGVectorToArray(gpos,ieq,neq,x,nx,ny,nz,xarr(:,:,:,ieq)
     .                         ,1)
      enddo

      call setMGBC(gpos,neq,nx,ny,nz,igrid,xarr,bcnd)

cc      write (*,*) bcnd(:,1)
cc      stop

cc      if (igrid == 4) then
cc        open(unit=110,file='debug.bin',form='unformatted'
cc     $       ,status='replace')
cc        call contour(xarr(0:nx+1,1:ny+1,1,1),nx+2,ny+1
cc     .       ,0d0,xmax,0d0,ymax,0,110)
cc        close(110)
cc        stop
cc      endif

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartsn)
            jac = jacobian(x1,y1,z1,cartsn)

            ijk    = i + nx*(j-1) + nx*ny*(k-1)
            ijkg   = ijk + isig - 1

            upwind = .5*(vcnv(ijkg,1)+abs(vcnv(ijkg,1)))
     .                 *( xarr(i  ,j,k,1) - xarr(i-1,j,k,1) )/dx(ig-1)
     .              +.5*(vcnv(ijkg,1)-abs(vcnv(ijkg,1)))
     .                 *( xarr(i+1,j,k,1) - xarr(i  ,j,k,1) )/dx(ig)
     .              +.5*(vcnv(ijkg,2)+abs(vcnv(ijkg,2)))
     .                 *( xarr(i,j  ,k,1) - xarr(i,j-1,k,1) )/dy(jg-1)
     .              +.5*(vcnv(ijkg,2)-abs(vcnv(ijkg,2)))
     .                 *( xarr(i,j+1,k,1) - xarr(i,j  ,k,1) )/dy(jg)
     .              +.5*(vcnv(ijkg,3)+abs(vcnv(ijkg,3)))
     .                 *( xarr(i,j,k  ,1) - xarr(i,j,k-1,1) )/dz(kg-1)
     .              +.5*(vcnv(ijkg,3)-abs(vcnv(ijkg,3)))
     .                 *( xarr(i,j,k+1,1) - xarr(i,j,k  ,1) )/dz(kg)

            upwind = upwind/jac
cc            nu2 = 0.1

            y(ijk) = (1./dt + alpha*(gamma-1.)*divV(ijkg))*x(ijk)
     .              + alpha*upwind 
cc     .              - alpha*nu2*laplacian(i,j,k,xarr)

            y(ijk) = y(ijk)*volume(i,j,k,igx,igy,igz)

          enddo
        enddo
      enddo

c End program

      deallocate(xarr)

      end subroutine test_mtvc

c tmp_mtvc
c####################################################################
      subroutine tmp_mtvc(gpos,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the Alfven SI system.
c--------------------------------------------------------------------

      use grid

      use grid_aliases

      use matvec

      use precond_variables

      use constants

      use operators

      use mlsolverSetup

      implicit none

c Call variables

      integer*4    ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,neq
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: dtmp

      real(8)    :: upwind,nu2

c External

      real*8       lap_vec,vis
      external     lap_vec,vis

c Begin program

      isig = grid_params%istartp(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      neq = ntot/nx/ny/nz

c Find limits for loops

      call limits(gpos,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(dtmp(0:nx+1,0:ny+1,0:nz+1,neq))

      !Set grid=1 because x is NOT a MG vector
      do ieq=1,neq
        call mapMGVectorToArray(gpos,ieq,neq,x,nx,ny,nz,dtmp(:,:,:,ieq)
     .                         ,1)
      enddo

      call setMGBC(gpos,neq,nx,ny,nz,igrid,dtmp,bcnd)

cc      if (gpos == 0 .and. igrid == 1) then
cc        open(unit=110,file='debug.bin',form='unformatted'
cc     .       ,status='replace')
cc        call contour(dtmp(1:nx,1:ny,1,1),nx,ny,0d0,xmax,0d0,ymax,0,110)
cc        close(110)
cc        stop
cc      endif

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartsn)
            jac = jacobian(x1,y1,z1,cartsn)

            ijk    = i + nx*(j-1) + nx*ny*(k-1)
            ijkg   = ijk + isig - 1

            upwind = .5*(vcnv(ijkg,1)+abs(vcnv(ijkg,1)))
     .                 *( dtmp(i  ,j,k,1) - dtmp(i-1,j,k,1) )/dx(ig-1)
     .              +.5*(vcnv(ijkg,1)-abs(vcnv(ijkg,1)))
     .                 *( dtmp(i+1,j,k,1) - dtmp(i  ,j,k,1) )/dx(ig)
     .              +.5*(vcnv(ijkg,2)+abs(vcnv(ijkg,2)))
     .                 *( dtmp(i,j  ,k,1) - dtmp(i,j-1,k,1) )/dy(jg-1)
     .              +.5*(vcnv(ijkg,2)-abs(vcnv(ijkg,2)))
     .                 *( dtmp(i,j+1,k,1) - dtmp(i,j  ,k,1) )/dy(jg)
     .              +.5*(vcnv(ijkg,3)+abs(vcnv(ijkg,3)))
     .                 *( dtmp(i,j,k  ,1) - dtmp(i,j,k-1,1) )/dz(kg-1)
     .              +.5*(vcnv(ijkg,3)-abs(vcnv(ijkg,3)))
     .                 *( dtmp(i,j,k+1,1) - dtmp(i,j,k  ,1) )/dz(kg)

            upwind = upwind/jac

            y(ijk) = ((1./dt + alpha*(gamma-1.)*divV(ijkg))*x(ijk)
     .               + alpha*upwind )*volume(i,j,k,igx,igy,igz)
     .               - alpha*chi*laplacian(i,j,k,dtmp)

          enddo
        enddo
      enddo

c End program

      deallocate(dtmp)

      end subroutine tmp_mtvc

c rho_mtvc
c####################################################################
      subroutine rho_mtvc(gpos,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the Alfven SI system.
c--------------------------------------------------------------------

      use grid

      use grid_aliases

      use matvec

      use precond_variables

      use constants

      use operators

      use mlsolverSetup

      implicit none

c Call variables

      integer*4    ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,neq
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: drho

      real(8)    :: upwind,nu2

c External

      real(8)    :: lap_vec,vis
      external   :: lap_vec,vis

c Begin program

      isig = grid_params%istartp(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      neq = ntot/nx/ny/nz

c Find limits for loops

      call limits(gpos,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(drho(0:nx+1,0:ny+1,0:nz+1,neq))

      !Set grid=1 because x is NOT a MG vector
      do ieq=1,neq
        call mapMGVectorToArray(gpos,ieq,neq,x,nx,ny,nz,drho(:,:,:,ieq)
     .                         ,1)
      enddo

      call setMGBC(gpos,neq,nx,ny,nz,igrid,drho,bcnd)

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartsn)
            jac = jacobian(x1,y1,z1,cartsn)

            ijk    = i + nx*(j-1) + nx*ny*(k-1)
            ijkg   = ijk + isig - 1

            upwind = .5*(vcnv(ijkg,1)+abs(vcnv(ijkg,1)))
     .                 *( drho(i  ,j,k,1) - drho(i-1,j,k,1) )/dx(ig-1)
     .              +.5*(vcnv(ijkg,1)-abs(vcnv(ijkg,1)))
     .                 *( drho(i+1,j,k,1) - drho(i  ,j,k,1) )/dx(ig)
     .              +.5*(vcnv(ijkg,2)+abs(vcnv(ijkg,2)))
     .                 *( drho(i,j  ,k,1) - drho(i,j-1,k,1) )/dy(jg-1)
     .              +.5*(vcnv(ijkg,2)-abs(vcnv(ijkg,2)))
     .                 *( drho(i,j+1,k,1) - drho(i,j  ,k,1) )/dy(jg)
     .              +.5*(vcnv(ijkg,3)+abs(vcnv(ijkg,3)))
     .                 *( drho(i,j,k  ,1) - drho(i,j,k-1,1) )/dz(kg-1)
     .              +.5*(vcnv(ijkg,3)-abs(vcnv(ijkg,3)))
     .                 *( drho(i,j,k+1,1) - drho(i,j,k  ,1) )/dz(kg)

            upwind = upwind/jac

            y(ijk) = ( (1./dt + alpha*divV(ijkg))*x(ijk)
     .               + alpha*upwind )*volume(i,j,k,igx,igy,igz)
     .               - alpha*dd*laplacian(i,j,k,drho)

          enddo
        enddo
      enddo

c End program

      deallocate(drho)

      end subroutine rho_mtvc

c b_mtvc
c####################################################################
      subroutine b_mtvc(gpos,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the Alfven SI system.
c--------------------------------------------------------------------

      use grid

      use grid_aliases

      use matvec

      use precond_variables

      use constants

      use operators

      use mlsolverSetup

      implicit none

c Call variables

      integer*4    ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,neq,ip,im,jp,jm,kp,km
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

      real(8),allocatable,dimension(:,:,:,:) :: bb
      real(8),allocatable,dimension(:,:,:)   :: one

      real(8)    :: upwind,etal,flxip,flxim,flxjp,flxjm,flxkp,flxkm

      real(8)    :: xim,yim,zim,xip,yip,zip
     .             ,xjm,yjm,zjm,xjp,yjp,zjp
     .             ,xkm,ykm,zkm,xkp,ykp,zkp

      real(8)    :: x0,y0,z0,xh,yh,zh

      real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .             ,jacp,jacm,jach,jac0

c External

      real(8)    :: res
      external   :: res

c Begin program

      isig = grid_params%istartp(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      neq = ntot/nx/ny/nz

c Find limits for loops

      call limits(gpos,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(bb(0:nx+1,0:ny+1,0:nz+1,neq))
      allocate(one(0:nx+1,0:ny+1,0:nz+1))
      one=1d0

      !Set grid=1 because x is NOT a MG vector
      do ieq=1,neq
        call mapMGVectorToArray(gpos,ieq,neq,x,nx,ny,nz,bb(:,:,:,ieq)
     .                         ,1)
      enddo

      call setMGBC(gpos,neq,nx,ny,nz,igrid,bb,bcnd)

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            ip = min(i+1,nx)
            im = max(i-1,1)
            jp = min(j+1,ny)
            jm = max(j-1,1)
            kp = min(k+1,nz)
            km = max(k-1,1)

            call getCoordinates(im,j,k,igx,igy,igz,ig,jg,kg,xim,yim,zim
     .                         ,cartsn)
            call getCoordinates(ip,j,k,igx,igy,igz,ig,jg,kg,xip,yip,zip
     .                         ,cartsn)
            call getCoordinates(i,jm,k,igx,igy,igz,ig,jg,kg,xjm,yjm,zjm
     .                         ,cartsn)
            call getCoordinates(i,jp,k,igx,igy,igz,ig,jg,kg,xjp,yjp,zjp
     .                         ,cartsn)
            call getCoordinates(i,j,km,igx,igy,igz,ig,jg,kg,xkm,ykm,zkm
     .                         ,cartsn)
            call getCoordinates(i,j,kp,igx,igy,igz,ig,jg,kg,xkp,ykp,zkp
     .                         ,cartsn)
            call getCoordinates(i,j,k ,igx,igy,igz,ig,jg,kg,x0,y0,z0
     .                         ,cartsn)

            jacip  = jacobian(xip,yip,zip,cartsn)
            jacim  = jacobian(xim,yim,zim,cartsn)
            jacjp  = jacobian(xjp,yjp,zjp,cartsn)
            jacjm  = jacobian(xjm,yjm,zjm,cartsn)
            jackp  = jacobian(xkp,ykp,zkp,cartsn)
            jackm  = jacobian(xkm,ykm,zkm,cartsn)
            jac    = jacobian(x0 ,y0 ,z0 ,cartsn)

            ijk    = i + nx*(j-1) + nx*ny*(k-1)

            ijkg    = i  + nx*(j -1) + nx*ny*(k -1) + isig - 1
            ipjkg   = ip + nx*(j -1) + nx*ny*(k -1) + isig - 1
            imjkg   = im + nx*(j -1) + nx*ny*(k -1) + isig - 1
            ijpkg   = i  + nx*(jp-1) + nx*ny*(k -1) + isig - 1
            ijmkg   = i  + nx*(jm-1) + nx*ny*(k -1) + isig - 1
            ijkpg   = i  + nx*(j -1) + nx*ny*(kp-1) + isig - 1
            ijkmg   = i  + nx*(j -1) + nx*ny*(km-1) + isig - 1

            etal = res(i,j,k,nx,ny,nz,igx,igy,igz)
cc            etal = eta

            flxjp = 0.5/(jac+jacjp)*(
     .           (    (vcnv(ijkg,2)+vcnv(ijpkg,2))
     .            +abs(vcnv(ijkg,2)+vcnv(ijpkg,2)) )*bb(i,j ,k,1)
     .          +(    (vcnv(ijkg,2)+vcnv(ijpkg,2))         
     .            -abs(vcnv(ijkg,2)+vcnv(ijpkg,2)) )*bb(i,jp,k,1) )
     .             -0.5/(jac+jacjp)*(
     .           (    (vcnv(ijkg,1)+vcnv(ijpkg,1))
     .            +abs(vcnv(ijkg,1)+vcnv(ijpkg,1)) )*bb(i,j ,k,2)
     .          +(    (vcnv(ijkg,1)+vcnv(ijpkg,1))          
     .            -abs(vcnv(ijkg,1)+vcnv(ijpkg,1)) )*bb(i,jp,k,2) )
            flxjm = 0.5/(jac+jacjm)*(
     .           (    (vcnv(ijkg,2)+vcnv(ijmkg,2))
     .            +abs(vcnv(ijkg,2)+vcnv(ijmkg,2)) )*bb(i,j ,k,1)
     .          +(    (vcnv(ijkg,2)+vcnv(ijmkg,2))          
     .            -abs(vcnv(ijkg,2)+vcnv(ijmkg,2)) )*bb(i,jm,k,1) )
     .             -0.5/(jac+jacjm)*(
     .           (    (vcnv(ijkg,1)+vcnv(ijmkg,1))
     .            +abs(vcnv(ijkg,1)+vcnv(ijmkg,1)) )*bb(i,j ,k,2)
     .          +(    (vcnv(ijkg,1)+vcnv(ijmkg,1))          
     .            -abs(vcnv(ijkg,1)+vcnv(ijmkg,1)) )*bb(i,jm,k,2) )

            flxkp = 0.5/(jac+jackp)*(
     .           (    (vcnv(ijkg,3)+vcnv(ijkpg,3))
     .            +abs(vcnv(ijkg,3)+vcnv(ijkpg,3)) )*bb(i,j,k ,1)
     .          +(    (vcnv(ijkg,3)+vcnv(ijkpg,3))          
     .            -abs(vcnv(ijkg,3)+vcnv(ijkpg,3)) )*bb(i,j,kp,1) )
     .             -0.5/(jac+jackp)*(
     .           (    (vcnv(ijkg,1)+vcnv(ijkpg,1))
     .            +abs(vcnv(ijkg,1)+vcnv(ijkpg,1)) )*bb(i,j,k ,3)
     .          +(    (vcnv(ijkg,1)+vcnv(ijkpg,1))          
     .            -abs(vcnv(ijkg,1)+vcnv(ijkpg,1)) )*bb(i,j,kp,3) )
            flxkm = 0.5/(jac+jackm)*(
     .           (    (vcnv(ijkg,3)+vcnv(ijkmg,3))
     .            +abs(vcnv(ijkg,3)+vcnv(ijkmg,3)) )*bb(i,j,k ,1)
     .          +(    (vcnv(ijkg,3)+vcnv(ijkmg,3))          
     .            -abs(vcnv(ijkg,3)+vcnv(ijkmg,3)) )*bb(i,j,km,1) )
     .             -0.5/(jac+jackm)*(
     .           (    (vcnv(ijkg,1)+vcnv(ijkmg,1))
     .            +abs(vcnv(ijkg,1)+vcnv(ijkmg,1)) )*bb(i,j,k ,3)
     .          +(    (vcnv(ijkg,1)+vcnv(ijkmg,1))          
     .            -abs(vcnv(ijkg,1)+vcnv(ijkmg,1)) )*bb(i,j,km,3) )

            upwind =  (flxjp-flxjm)/dyh(jg)
     .               +(flxkp-flxkm)/dzh(kg)

            y(neq*(ijk-1)+1) = ( bb(i,j,k,1)/dt
     .               + alpha*upwind )*volume(i,j,k,igx,igy,igz)
     .               - alpha*etal
cc     .                      *laplacian(i,j,k,bb(:,:,:,1))
     .                      *veclaplacian(i,j,k,bb(:,:,:,1)
     .                                         ,bb(:,:,:,2)
     .                                         ,bb(:,:,:,3)
     .                                         ,ones,.false.,1)

            flxip = 0.5/(jac+jacip)*(
     .           (    (vcnv(ijkg,1)+vcnv(ipjkg,1))
     .            +abs(vcnv(ijkg,1)+vcnv(ipjkg,1)) )*bb(i ,j,k,2)
     .          +(    (vcnv(ijkg,1)+vcnv(ipjkg,1))          
     .            -abs(vcnv(ijkg,1)+vcnv(ipjkg,1)) )*bb(ip,j,k,2) )
     .             -0.5/(jac+jacip)*(
     .           (    (vcnv(ijkg,2)+vcnv(ipjkg,2))
     .            +abs(vcnv(ijkg,2)+vcnv(ipjkg,2)) )*bb(i ,j,k,1)
     .          +(    (vcnv(ijkg,2)+vcnv(ipjkg,2))          
     .            -abs(vcnv(ijkg,2)+vcnv(ipjkg,2)) )*bb(ip,j,k,1) )
            flxim = 0.5/(jac+jacim)*(
     .           (    (vcnv(ijkg,1)+vcnv(imjkg,1))
     .            +abs(vcnv(ijkg,1)+vcnv(imjkg,1)) )*bb(i ,j,k,2)
     .          +(    (vcnv(ijkg,1)+vcnv(imjkg,1))          
     .            -abs(vcnv(ijkg,1)+vcnv(imjkg,1)) )*bb(im,j,k,2) )
     .             -0.5/(jac+jacim)*(
     .           (    (vcnv(ijkg,2)+vcnv(imjkg,2))
     .            +abs(vcnv(ijkg,2)+vcnv(imjkg,2)) )*bb(i ,j,k,1)
     .          +(    (vcnv(ijkg,2)+vcnv(imjkg,2))          
     .            -abs(vcnv(ijkg,2)+vcnv(imjkg,2)) )*bb(im,j,k,1) )

            flxkp = 0.5/(jac+jackp)*(
     .           (    (vcnv(ijkg,3)+vcnv(ijkpg,3))
     .            +abs(vcnv(ijkg,3)+vcnv(ijkpg,3)) )*bb(i,j,k ,2)
     .          +(    (vcnv(ijkg,3)+vcnv(ijkpg,3))          
     .            -abs(vcnv(ijkg,3)+vcnv(ijkpg,3)) )*bb(i,j,kp,2) )
     .             -0.5/(jac+jackp)*(
     .           (    (vcnv(ijkg,2)+vcnv(ijkpg,2))
     .            +abs(vcnv(ijkg,2)+vcnv(ijkpg,2)) )*bb(i,j,k ,3)
     .          +(    (vcnv(ijkg,2)+vcnv(ijkpg,2))          
     .            -abs(vcnv(ijkg,2)+vcnv(ijkpg,2)) )*bb(i,j,kp,3) )
            flxkm = 0.5/(jac+jackm)*(
     .           (    (vcnv(ijkg,3)+vcnv(ijkmg,3))
     .            +abs(vcnv(ijkg,3)+vcnv(ijkmg,3)) )*bb(i,j,k ,2)
     .          +(    (vcnv(ijkg,3)+vcnv(ijkmg,3))          
     .            -abs(vcnv(ijkg,3)+vcnv(ijkmg,3)) )*bb(i,j,km,2) )
     .             -0.5/(jac+jackm)*(
     .           (    (vcnv(ijkg,2)+vcnv(ijkmg,2))
     .            +abs(vcnv(ijkg,2)+vcnv(ijkmg,2)) )*bb(i,j,k ,3)
     .          +(    (vcnv(ijkg,2)+vcnv(ijkmg,2))          
     .            -abs(vcnv(ijkg,2)+vcnv(ijkmg,2)) )*bb(i,j,km,3) )

            upwind =  (flxip-flxim)/dxh(ig)
     .               +(flxkp-flxkm)/dzh(kg)

            y(neq*(ijk-1)+2) = ( bb(i,j,k,2)/dt
     .               + alpha*upwind )*volume(i,j,k,igx,igy,igz)
     .               - alpha*etal
cc     .                      *laplacian(i,j,k,bb(:,:,:,2))
     .                      *veclaplacian(i,j,k,bb(:,:,:,1)
     .                                         ,bb(:,:,:,2)
     .                                         ,bb(:,:,:,3)
     .                                         ,ones,.false.,2)

            flxip = 0.5/(jac+jacip)*(
     .           (    (vcnv(ipjkg,1)+vcnv(ijkg,1))
     .            +abs(vcnv(ipjkg,1)+vcnv(ijkg,1)) )*bb(i ,j,k,3)
     .          +(    (vcnv(ipjkg,1)+vcnv(ijkg,1))          
     .            -abs(vcnv(ipjkg,1)+vcnv(ijkg,1)) )*bb(ip,j,k,3) )
     .             -0.5/(jac+jacip)*(
     .           (    (vcnv(ipjkg,3)+vcnv(ijkg,3))
     .            +abs(vcnv(ipjkg,3)+vcnv(ijkg,3)) )*bb(i ,j,k,1)
     .          +(    (vcnv(ipjkg,3)+vcnv(ijkg,3))          
     .            -abs(vcnv(ipjkg,3)+vcnv(ijkg,3)) )*bb(ip,j,k,1) )
            flxim = 0.5/(jac+jacim)*(
     .           (    (vcnv(imjkg,1)+vcnv(ijkg,1))
     .            +abs(vcnv(imjkg,1)+vcnv(ijkg,1)) )*bb(i ,j,k,3)
     .          +(    (vcnv(imjkg,1)+vcnv(ijkg,1))          
     .            -abs(vcnv(imjkg,1)+vcnv(ijkg,1)) )*bb(im,j,k,3) )
     .             -0.5/(jac+jacim)*(
     .           (    (vcnv(imjkg,3)+vcnv(ijkg,3))
     .            +abs(vcnv(imjkg,3)+vcnv(ijkg,3)) )*bb(i ,j,k,1)
     .          +(    (vcnv(imjkg,3)+vcnv(ijkg,3))          
     .            -abs(vcnv(imjkg,3)+vcnv(ijkg,3)) )*bb(im,j,k,1) )

            flxjp = 0.5/(jac+jacjp)*(
     .           (    (vcnv(ijkg,2)+vcnv(ijpkg,2))
     .            +abs(vcnv(ijkg,2)+vcnv(ijpkg,2)) )*bb(i,j ,k,3)
     .          +(    (vcnv(ijkg,2)+vcnv(ijpkg,2))          
     .            -abs(vcnv(ijkg,2)+vcnv(ijpkg,2)) )*bb(i,jp,k,3) )
     .             -0.5/(jac+jacjp)*(
     .           (    (vcnv(ijkg,3)+vcnv(ijpkg,3))
     .            +abs(vcnv(ijkg,3)+vcnv(ijpkg,3)) )*bb(i,j ,k,2)
     .          +(    (vcnv(ijkg,3)+vcnv(ijpkg,3))          
     .            -abs(vcnv(ijkg,3)+vcnv(ijpkg,3)) )*bb(i,jp,k,2) )
            flxjm = 0.5/(jac+jacjm)*(
     .           (    (vcnv(ijkg,2)+vcnv(ijmkg,2))
     .            +abs(vcnv(ijkg,2)+vcnv(ijmkg,2)) )*bb(i,j ,k,3)
     .          +(    (vcnv(ijkg,2)+vcnv(ijmkg,2))          
     .            -abs(vcnv(ijkg,2)+vcnv(ijmkg,2)) )*bb(i,jm,k,3) )
     .             -0.5/(jac+jacjm)*(
     .           (    (vcnv(ijkg,3)+vcnv(ijmkg,3))
     .            +abs(vcnv(ijkg,3)+vcnv(ijmkg,3)) )*bb(i,j ,k,2)
     .          +(    (vcnv(ijkg,3)+vcnv(ijmkg,3))          
     .            -abs(vcnv(ijkg,3)+vcnv(ijmkg,3)) )*bb(i,jm,k,2) )

            upwind =  (flxip-flxim)/dxh(ig)
     .               +(flxjp-flxjm)/dyh(jg)

            y(neq*(ijk-1)+3) = ( bb(i,j,k,3)/dt
     .               + alpha*upwind )*volume(i,j,k,igx,igy,igz)
     .               - alpha*etal
cc     .                      *laplacian(i,j,k,bb(:,:,:,3))
     .                      *veclaplacian(i,j,k,bb(:,:,:,1)
     .                                         ,bb(:,:,:,2)
     .                                         ,bb(:,:,:,3)
     .                                         ,one,.false.,3)
          enddo
        enddo
      enddo

c End program

      deallocate(bb)
      deallocate(one)

      end subroutine b_mtvc

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
