c setMGBC
c####################################################################
      subroutine setMGBC(gpos,neq,nnx,nny,nnz,iig,array,bcnd)
c--------------------------------------------------------------------
c     Interfaces BC routines with MG code, for preconditioning.
c     NOTE: we assume that array, when representing a vector quantity,
c           is in the contravariant representation.
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

c Define variables for BC routines

      !igx,igy,igz,nx,ny,nz are needed for internal BC routines
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

c$$$      imin = 1
c$$$      imax = nnx
c$$$      jmin = 1
c$$$      jmax = nny
c$$$      kmin = 1
c$$$      kmax = nnz

      call limits(gpos,nnx,nny,nnz,imin,imax,jmin,jmax,kmin,kmax)

c Check limits (return if not close to boundaries)

      if (     (imin > 1 .and. imax < nnx)
     .    .and.(jmin > 1 .and. jmax < nny)
     .    .and.(kmin > 1 .and. kmax < nnz)) return

c Offset limits to cover stencil

      offset = 1

      iimin = max(imin-offset,1)
      iimax = min(imax+offset,nnx)
      jjmin = max(jmin-offset,1)  
      jjmax = min(jmax+offset,nny)
      kkmin = max(kmin-offset,1)  
      kkmax = min(kmax+offset,nnz)

c Select operation

      select case (neq)
      case(1)

        call setBC(IRHO,array(:,:,:,neq),zeros,bcnd(:,1))

cc        call setBC(icomp,array(:,:,:,neq),zeros,bcnd(:,1)) !Velocities for T BC unknown

      case(3)

        allocate(v_cov (0:nnx+1,0:nny+1,0:nnz+1,neq))

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

cc        if (gpos == 0) icomp = -1  !Set icomp to -1 for safety, to capture undefined instances

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
c     This subroutine is a matvec test, y = A.x.
c     In call:
c      * gpos: vector index of position on the numerical grid
c             gpos = i + nx*(j-1) + ny*nx*(k-1)
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c--------------------------------------------------------------------

      use matvec

      use precond_variables

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
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,nu2

c External

      real*8       vis
      external     vis

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
      call mapMGVectorToArray(gpos,neq,x,nx,ny,nz,xarr,1)

      icomp = IRHO
      call setMGBC(gpos,neq,nx,ny,nz,igrid,xarr,bcnd)

c Map velocity components

      v0_cnv => gv0%grid(igrid)%array

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartsn)
            jac = jacobian(x1,y1,z1,cartsn)

            ijk    = i + nx*(j-1) + nx*ny*(k-1)
            ijkg   = ijk + isig - 1

            upwind = .5*(v0_cnv(i,j,k,1)+abs(v0_cnv(i,j,k,1)))
     .                 *( xarr(i  ,j,k,1) - xarr(i-1,j,k,1) )/dx(ig-1)
     .              +.5*(v0_cnv(i,j,k,1)-abs(v0_cnv(i,j,k,1)))
     .                 *( xarr(i+1,j,k,1) - xarr(i  ,j,k,1) )/dx(ig)
     .              +.5*(v0_cnv(i,j,k,2)+abs(v0_cnv(i,j,k,2)))
     .                 *( xarr(i,j  ,k,1) - xarr(i,j-1,k,1) )/dy(jg-1)
     .              +.5*(v0_cnv(i,j,k,2)-abs(v0_cnv(i,j,k,2)))
     .                 *( xarr(i,j+1,k,1) - xarr(i,j  ,k,1) )/dy(jg)
     .              +.5*(v0_cnv(i,j,k,3)+abs(v0_cnv(i,j,k,3)))
     .                 *( xarr(i,j,k  ,1) - xarr(i,j,k-1,1) )/dz(kg-1)
     .              +.5*(v0_cnv(i,j,k,3)-abs(v0_cnv(i,j,k,3)))
     .                 *( xarr(i,j,k+1,1) - xarr(i,j,k  ,1) )/dz(kg)

            upwind = upwind/jac
cc            nu2 = 0.1

            y(ijk) = (1./dt + alpha*(gamma-1.)*mgdivV0(ijkg))*x(ijk)
     .              + alpha*upwind 
cc     .              - alpha*nu2*laplacian(i,j,k,nx,ny,nz,xarr)

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
c     for the energy equation.
c     In call:
c      * gpos: vector index of position on the numerical grid
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D
c              , 27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c--------------------------------------------------------------------

      use matvec

      use precond_variables

      use mlsolverSetup

      implicit none

c Call variables

      integer*4    ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,neq
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: dtmp
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,nu2

c External

      real*8       vis
      external     vis

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

      call limits(abs(gpos),nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(dtmp(0:nx+1,0:ny+1,0:nz+1,neq))

      dtmp = 0d0

      !Set grid=1 because x is NOT a MG vector
      call mapMGVectorToArray(max(0,gpos),neq,x,nx,ny,nz,dtmp,1)

      icomp = ITMP
      call setMGBC(max(0,gpos),neq,nx,ny,nz,igrid,dtmp,bcnd)

c Map velocity components

      v0_cnv => gv0%grid(igrid)%array

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartsn)
            jac = jacobian(x1,y1,z1,cartsn)

            ijk    = i + nx*(j-1) + nx*ny*(k-1)
            ijkg   = ijk + isig - 1

            upwind = .5*(v0_cnv(i,j,k,1)+abs(v0_cnv(i,j,k,1)))
     .                 *( dtmp(i  ,j,k,1) - dtmp(i-1,j,k,1) )/dx(ig-1)
     .              +.5*(v0_cnv(i,j,k,1)-abs(v0_cnv(i,j,k,1)))
     .                 *( dtmp(i+1,j,k,1) - dtmp(i  ,j,k,1) )/dx(ig)
     .              +.5*(v0_cnv(i,j,k,2)+abs(v0_cnv(i,j,k,2)))
     .                 *( dtmp(i,j  ,k,1) - dtmp(i,j-1,k,1) )/dy(jg-1)
     .              +.5*(v0_cnv(i,j,k,2)-abs(v0_cnv(i,j,k,2)))
     .                 *( dtmp(i,j+1,k,1) - dtmp(i,j  ,k,1) )/dy(jg)
     .              +.5*(v0_cnv(i,j,k,3)+abs(v0_cnv(i,j,k,3)))
     .                 *( dtmp(i,j,k  ,1) - dtmp(i,j,k-1,1) )/dz(kg-1)
     .              +.5*(v0_cnv(i,j,k,3)-abs(v0_cnv(i,j,k,3)))
     .                 *( dtmp(i,j,k+1,1) - dtmp(i,j,k  ,1) )/dz(kg)

            upwind = upwind/jac

            y(ijk) = ((1./dt + 0*alpha*(gamma-1.)*mgdivV0(ijkg))*x(ijk)
     .               + alpha*upwind )*volume(i,j,k,igx,igy,igz)
     .               - alpha*chi*laplacian(i,j,k,nx,ny,nz,dtmp(:,:,:,1))

          enddo
        enddo
      enddo

c End program

      deallocate(dtmp)
      nullify(v0_cnv)

      end subroutine tmp_mtvc

c rho_mtvc
c####################################################################
      subroutine rho_mtvc(gpos,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A.x  matrix-free
c     for the continuity equation.
c     In call:
c      * gpos: vector index of position on the numerical grid
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D
c              , 27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c--------------------------------------------------------------------

      use matvec

      use precond_variables

      use mlsolverSetup

      implicit none

c Call variables

      integer*4    ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,neq
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: drho
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,nu2

c External

      real(8)    :: vis
      external   :: vis

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

      call limits(abs(gpos),nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(drho(0:nx+1,0:ny+1,0:nz+1,neq))

      drho = 0d0

      !Set grid=1 because x is NOT a MG vector
      !For GS, gpos < 0 so that the whole vector x is mapped
      !For finding the diagonal, gpos > 0
      call mapMGVectorToArray(max(0,gpos),neq,x,nx,ny,nz,drho,1)

      icomp = IRHO
      call setMGBC(max(0,gpos),1,nx,ny,nz,igrid,drho,bcnd)

c Map velocity components

      v0_cnv => gv0%grid(igrid)%array

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartsn)
            jac = jacobian(x1,y1,z1,cartsn)

            ijk    = i + nx*(j-1) + nx*ny*(k-1)
            ijkg   = ijk + isig - 1

            upwind = .5*(v0_cnv(i,j,k,1)+abs(v0_cnv(i,j,k,1)))
     .                 *( drho(i  ,j,k,1) - drho(i-1,j,k,1) )/dx(ig-1)
     .              +.5*(v0_cnv(i,j,k,1)-abs(v0_cnv(i,j,k,1)))
     .                 *( drho(i+1,j,k,1) - drho(i  ,j,k,1) )/dx(ig)
     .              +.5*(v0_cnv(i,j,k,2)+abs(v0_cnv(i,j,k,2)))
     .                 *( drho(i,j  ,k,1) - drho(i,j-1,k,1) )/dy(jg-1)
     .              +.5*(v0_cnv(i,j,k,2)-abs(v0_cnv(i,j,k,2)))
     .                 *( drho(i,j+1,k,1) - drho(i,j  ,k,1) )/dy(jg)
     .              +.5*(v0_cnv(i,j,k,3)+abs(v0_cnv(i,j,k,3)))
     .                 *( drho(i,j,k  ,1) - drho(i,j,k-1,1) )/dz(kg-1)
     .              +.5*(v0_cnv(i,j,k,3)-abs(v0_cnv(i,j,k,3)))
     .                 *( drho(i,j,k+1,1) - drho(i,j,k  ,1) )/dz(kg)

            upwind = upwind/jac

            y(ijk) = ( (1./dt + alpha*mgdivV0(ijkg))*x(ijk)
     .               + alpha*upwind )*volume(i,j,k,igx,igy,igz)
     .               - alpha*dd*laplacian(i,j,k,nx,ny,nz,drho(:,:,:,1))

          enddo
        enddo
      enddo

c End program

      deallocate(drho)
      nullify(v0_cnv)

      end subroutine rho_mtvc

c b_mtvc
c####################################################################
      subroutine b_mtvc(gpos,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for Faraday's law.
c     In call:
c      * gpos: vector index of position on the numerical grid
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D
c              , 27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c--------------------------------------------------------------------

      use matvec

      use precond_variables

      use mlsolverSetup

      implicit none

c Call variables

      integer*4    ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,neq,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

      real(8)    :: upwind,etal,flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol

      real(8),allocatable,dimension(:,:,:,:) :: db
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

c External

      real(8)    :: res
      external   :: res

c Begin program

      isig = grid_params%istartp(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

      nxx = grid_params%nxv(igx)
      nyy = grid_params%nyv(igy)
      nzz = grid_params%nzv(igz)

      neq = ntot/nxx/nyy/nzz

c Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(db(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      db = 0d0

      !Set grid=1 because x is NOT a MG vector
      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,db,1)

      icomp = IBX
      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,db,bcnd)

c Velocity field (including BCs)

      v0_cnv => gv0%grid(igrid)%array

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            ip = i+1
            im = i-1
            jp = j+1
            jm = j-1
            kp = k+1
            km = k-1

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

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            vol = volume(i,j,k,igx,igy,igz)

            etal = res(i,j,k,nxx,nyy,nzz,igx,igy,igz)

            flxjp = 0.5/(jac+jacjp)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,j ,k,1)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))         
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,jp,k,1))
     .             -0.5/(jac+jacjp)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1)) )*db(i,j ,k,2)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1)) )*db(i,jp,k,2))
            flxjm = 0.5/(jac+jacjm)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,j ,k,1)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,jm,k,1))
     .             -0.5/(jac+jacjm)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1)) )*db(i,j ,k,2)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1)) )*db(i,jm,k,2))

            flxkp = 0.5/(jac+jackp)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,k ,1)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,kp,1))
     .             -0.5/(jac+jackp)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1)) )*db(i,j,k ,3)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1)) )*db(i,j,kp,3))
            flxkm = 0.5/(jac+jackm)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,k ,1)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,km,1))
     .             -0.5/(jac+jackm)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1)) )*db(i,j,k ,3)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1)) )*db(i,j,km,3))

            upwind =  (flxjp-flxjm)/dyh(jg)
     .               +(flxkp-flxkm)/dzh(kg)

            y(neq*(ijk-1)+1) = ( db(i,j,k,1)/dt
     .               + alpha*upwind )*vol
     .               - alpha*etal
     .                      *veclaplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,1)
     .                                                  ,db(:,:,:,2)
     .                                                  ,db(:,:,:,3)
     .                                                  ,ones,.false.,1)
cc     .                      *laplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,1))

            flxip = 0.5/(jac+jacip)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(i ,j,k,2)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(ip,j,k,2))
     .             -0.5/(jac+jacip)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2)) )*db(i ,j,k,1)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2)) )*db(ip,j,k,1))
            flxim = 0.5/(jac+jacim)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(i ,j,k,2)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(im,j,k,2))
     .             -0.5/(jac+jacim)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )*db(i ,j,k,1)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )*db(im,j,k,1))

            flxkp = 0.5/(jac+jackp)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,k ,2)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,kp,2))
     .             -0.5/(jac+jackp)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2)) )*db(i,j,k ,3)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2)) )*db(i,j,kp,3))
            flxkm = 0.5/(jac+jackm)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,k ,2)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,km,2))
     .             -0.5/(jac+jackm)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2)) )*db(i,j,k ,3)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2)) )*db(i,j,km,3))

            upwind =  (flxip-flxim)/dxh(ig)
     .               +(flxkp-flxkm)/dzh(kg)

            y(neq*(ijk-1)+2) = ( db(i,j,k,2)/dt
     .               + alpha*upwind )*vol
     .               - alpha*etal
     .                      *veclaplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,1)
     .                                                  ,db(:,:,:,2)
     .                                                  ,db(:,:,:,3)
     .                                                  ,ones,.false.,2)
cc     .                      *laplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,2))

            flxip = 0.5/(jac+jacip)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(i ,j,k,3)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(ip,j,k,3))
     .             -0.5/(jac+jacip)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3)) )*db(i ,j,k,1)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3)) )*db(ip,j,k,1))
            flxim = 0.5/(jac+jacim)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(i ,j,k,3)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(im,j,k,3))
     .             -0.5/(jac+jacim)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3)) )*db(i ,j,k,1)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3)) )*db(im,j,k,1))

            flxjp = 0.5/(jac+jacjp)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,j ,k,3)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,jp,k,3))
     .             -0.5/(jac+jacjp)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,jp,k,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,jp,k,3)) )*db(i,j ,k,2)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,jp,k,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,jp,k,3)) )*db(i,jp,k,2))
            flxjm = 0.5/(jac+jacjm)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,j ,k,3)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,jm,k,3))
     .             -0.5/(jac+jacjm)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3)) )*db(i,j ,k,2)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3)) )*db(i,jm,k,2))

            upwind =  (flxip-flxim)/dxh(ig)
     .               +(flxjp-flxjm)/dyh(jg)

            y(neq*(ijk-1)+3) = ( db(i,j,k,3)/dt
     .               + alpha*upwind )*vol
     .               - alpha*etal
     .                      *veclaplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,1)
     .                                                     ,db(:,:,:,2)
     .                                                     ,db(:,:,:,3)
     .                                                 ,ones,.false.,3)
cc     .                      *laplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,3))

          enddo
        enddo
      enddo

c End program

      deallocate(db)
      nullify(v0_cnv)

      end subroutine b_mtvc

c v_mtvc
c####################################################################
      subroutine v_mtvc(gpos,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the velocity SI system.
c     In call:
c      * gpos: vector index of position on the numerical grid:
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D
c              , 27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c--------------------------------------------------------------------

      use matvec

      use precond_variables

      use mlsolverSetup

      use nlfunction_setup

      implicit none

c Call variables

      integer*4    ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,neq,ip,im,jp,jm,kp,km,hex,hey,hez,icompsave
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax,nxx,nyy,nzz,ijk,ijkg

      real(8),allocatable,dimension(:,:,:,:) :: dv
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv,pp0,b0_cnv,rho0

      real(8)    :: upwind,mul,vol,nabla_vv0(3,3)
     $             ,flxip,flxim,flxjp,flxjm,flxkp,flxkm
     $             ,psiv(3),psib(3),psit(3),cov(3),cnv(3)
     $             ,divip,divim,divjp,divjm,divkp,divkm

      real(8)    :: xiph,ximh,xjph,xjmh,xkph,xkmh
     .             ,yiph,yimh,yjph,yjmh,ykph,ykmh
     .             ,ziph,zimh,zjph,zjmh,zkph,zkmh

      real(8)    :: a10,a20,a30
     .             ,a1covip,a2covip,a3covip,a1covim,a2covim,a3covim
     .             ,a1covjp,a2covjp,a3covjp,a1covjm,a2covjm,a3covjm
     .             ,a1covkp,a2covkp,a3covkp,a1covkm,a2covkm,a3covkm
     .             ,a1cnvip,a2cnvip,a3cnvip,a1cnvim,a2cnvim,a3cnvim
     .             ,a1cnvjp,a2cnvjp,a3cnvjp,a1cnvjm,a2cnvjm,a3cnvjm
     .             ,a1cnvkp,a2cnvkp,a3cnvkp,a1cnvkm,a2cnvkm,a3cnvkm

c External

      real(8)    :: vis
      external   :: vis

c Begin program

      isig = grid_params%istartp(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

      nxx = grid_params%nxv(igx)
      nyy = grid_params%nyv(igy)
      nzz = grid_params%nzv(igz)

      neq = ntot/nxx/nyy/nzz

c Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(dv(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      dv = 0d0

      !Set grid=1 because x is NOT a MG vector
      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,1)

      icomp = IVX
      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd)

c Define pointers to MG arrays

      rho0   => grho0%grid(igrid)%array

      pp0    => gp0%grid(igrid)%array

      v0_cnv => gv0%grid(igrid)%array

      b0_cnv => gb0%grid(igrid)%array

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            !Preparations
            ip = i+1
            im = i-1
            jp = j+1
            jm = j-1
            kp = k+1
            km = k-1

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

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = i  + nxx*(j -1) + nxx*nyy*(k -1) + isig - 1

            vol = volume(i,j,k,igx,igy,igz)

            mul = vis(i,j,k,nxx,nyy,nzz,igx,igy,igz)

            !P_si^v  ******************************

            hex = 1
            hey = 1
            hez = 1
            if (v0_cnv(i,j,k,1) > 0d0) hex = -1
            if (v0_cnv(i,j,k,2) > 0d0) hey = -1
            if (v0_cnv(i,j,k,3) > 0d0) hez = -1

            nabla_v = fnabla_v_upwd(i,j,k,nxx,nyy,nzz,dv(:,:,:,1)
     .                                               ,dv(:,:,:,2)
     .                                               ,dv(:,:,:,3)
     .                                               ,hex,hey,hez)

            nabla_vv0 = fnabla_v_upwd(i,j,k,nxx,nyy,nzz,v0_cnv(:,:,:,1)
     .                                                 ,v0_cnv(:,:,:,2)
     .                                                 ,v0_cnv(:,:,:,3)
     .                                                 ,0,0,0)

            do ieq=1,3
              upwind =( v0_cnv(i,j,k,1)*nabla_v(1,ieq)
     .                 +v0_cnv(i,j,k,2)*nabla_v(2,ieq)
     .                 +v0_cnv(i,j,k,3)*nabla_v(3,ieq))/jac

              upwind = upwind
     .               +( dv(i,j,k,1)*nabla_vv0(1,ieq)
     .                 +dv(i,j,k,2)*nabla_vv0(2,ieq)
     .                 +dv(i,j,k,3)*nabla_vv0(3,ieq))/jac

              hex = floor(sign(1d0,-mgadvdiffV0(ijkg,ieq)))
              upwind = upwind
     .            -dt*mgadvdiffV0(ijkg,ieq)*div_upwd(hex)/rho0(i,j,k,1)

              psiv(ieq) = alpha*upwind*vol
     .                  - alpha*mul
     .                      *veclaplacian(i,j,k,nxx,nyy,nzz,dv(:,:,:,1)
     .                                                     ,dv(:,:,:,2)
     .                                                     ,dv(:,:,:,3)
     .                                   ,ones,alt_eom,ieq)
            enddo

            psiv = rho0(i,j,k,1)*psiv

            !P_si^T  ******************************

            !Fluxes at faces for calculation of grad(dv.grad(p0))
            flxip =( (dv(i,j,k,1)+dv(ip,j,k,1))
     .                 *(pp0(ip,j,k,1)-pp0(i,j,k,1))/dx(ig)
     .              +(dv(i,j,k,2)+dv(ip,j,k,2))
     .                 *(pp0(ip,jp,k,1)-pp0(ip,jm,k,1)
     .                  +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
     .              +(dv(i,j,k,3)+dv(ip,j,k,3))
     .                 *(pp0(ip,j,kp,1)-pp0(ip,j,km,1)
     .                  +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4. )
     $              /(jac+jacip)
            flxim =( (dv(i,j,k,1)+dv(im,j,k,1))
     .                  *(pp0(i ,j,k,1)-pp0(im,j,k,1))/dx(ig-1)
     .              +(dv(i,j,k,2)+dv(im,j,k,2))
     .                  *(pp0(im,jp,k,1)-pp0(im,jm,k,1)
     .                   +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
     .              +(dv(i,j,k,3)+dv(im,j,k,3))
     .                  *(pp0(im,j,kp,1)-pp0(im,j,km,1)
     .                   +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4.)
     $              /(jac+jacim)

            flxjp =( (dv(i,j,k,1)+dv(i,jp,k,1))
     .                  *(pp0(ip,jp,k,1)-pp0(im,jp,k,1)
     .                   +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/dxh(ig)/4.
     .              +(dv(i,j,k,2)+dv(i,jp,k,2))
     .                  *(pp0(i,jp,k,1)-pp0(i,j,k,1))/dy(jg)
     .              +(dv(i,j,k,3)+dv(i,jp,k,3))
     .                  *(pp0(i,jp,kp,1)-pp0(i,jp,km,1)
     .                   +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/dzh(kg)/4.)
     $              /(jac+jacjp)
            flxjm =( (dv(i,j,k,1)+dv(i,jm,k,1))
     .                  *(pp0(ip,jm,k,1)-pp0(im,jm,k,1)
     .                   +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/dxh(ig)/4.
     .              +(dv(i,j,k,2)+dv(i,jm,k,2))
     .                  *(pp0(i,j ,k,1)-pp0(i,jm,k,1))/dy(jg-1)
     .              +(dv(i,j,k,3)+dv(i,jm,k,3))
     .                  *(pp0(i,jm,kp,1)-pp0(i,jm,km,1)
     .                   +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/dzh(kg)/4.)
     $              /(jac+jacjm)

            flxkp =( (dv(i,j,k,1)+dv(i,j,kp,1))
     .                  *(pp0(ip,j,kp,1)-pp0(im,j,kp,1)
     .                   +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/dxh(ig)/4.
     .              +(dv(i,j,k,2)+dv(i,j,kp,2))
     .                  *(pp0(i,jp,kp,1)-pp0(i,jm,kp,1)
     .                   +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/dyh(jg)/4.
     .              +(dv(i,j,k,3)+dv(i,j,kp,3))
     .                  *(pp0(i,j,kp,1)-pp0(i,j,k,1))/dz(kg) )
     $              /(jac+jackp)
            flxkm =( (dv(i,j,k,1)+dv(i,j,km,1))
     .                  *(pp0(ip,j,km,1)-pp0(im,j,km,1)
     .                   +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/dxh(ig)/4.
     .              +(dv(i,j,k,2)+dv(i,j,km,2))
     .                  *(pp0(i,jp,km,1)-pp0(i,jm,km,1)
     .                   +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/dyh(jg)/4.
     .              +(dv(i,j,k,3)+dv(i,j,km,3))
     .                  *(pp0(i,j,k ,1)-pp0(i,j,km,1))/dz(kg-1) )
     $              /(jac+jackm)

c$$$            flxip = 0d0
c$$$            flxim = 0d0
c$$$            flxjp = 0d0
c$$$            flxjm = 0d0
c$$$            flxkp = 0d0
c$$$            flxkm = 0d0

            !Fluxes at faces for calculation of grad(gamma*p0*div(dv))
cc
cc            ip = min(i+1,nx)
cc            im = max(i-1,1)
cc            jp = min(j+1,ny)
cc            jm = max(j-1,1)
cc            kp = min(k+1,nz)
cc            km = max(k-1,1)

            !!Divergence at faces i+-1/2, etc.
            divip = (dv(ip,j ,k,1)-dv(i ,j ,k,1))/dx(ig)
     .             +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .              +dv(ip,jp,k,2)-dv(ip,jm,k,2))/dyh(jg)/4.
     .             +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .              +dv(ip,j,kp,3)-dv(ip,j,km,3))/dzh(kg)/4.
            divim = (dv(i ,j ,k,1)-dv(im,j ,k,1))/dx(ig-1)
     .             +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .              +dv(im,jp,k,2)-dv(im,jm,k,2))/dyh(jg)/4.
     .             +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .              +dv(im,j,kp,3)-dv(im,j,km,3))/dzh(kg)/4.

            divjp = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .              +dv(ip,jp,k,1)-dv(im,jp,k,1))/dxh(ig)/4.
     .             +(dv(i ,jp,k,2)-dv(i ,j ,k,2))/dy(jg)
     .             +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .              +dv(i,jp,kp,3)-dv(i,jp,km,3))/dzh(kg)/4.
            divjm = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .              +dv(ip,jm,k,1)-dv(im,jm,k,1))/dxh(ig)/4.
     .             +(dv(i ,j ,k,2)-dv(i ,jm,k,2))/dy(jg-1)
     .             +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .              +dv(i,jm,kp,3)-dv(i,jm,km,3))/dzh(kg)/4.

            divkp = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .              +dv(ip,j,kp,1)-dv(im,j,kp,1))/dxh(ig)/4.
     .             +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .              +dv(i,jp,kp,2)-dv(i,jm,kp,2))/dyh(jg)/4.
     .             +(dv(i,j ,kp,3)-dv(i,j ,k ,3))/dz(kg)
            divkm = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .              +dv(ip,j,km,1)-dv(im,j,km,1))/dxh(ig)/4.
     .             +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .              +dv(i,jp,km,2)-dv(i,jm,km,2))/dyh(jg)/4.
     .             +(dv(i,j ,k ,3)-dv(i,j ,km,3))/dz(kg-1)

            flxip = flxip
     .            + gamma*(pp0(i,j,k,1)+pp0(ip,j,k,1))*divip/(jac+jacip)
            flxim = flxim
     .            + gamma*(pp0(i,j,k,1)+pp0(im,j,k,1))*divim/(jac+jacim)

            flxjp = flxjp
     .            + gamma*(pp0(i,j,k,1)+pp0(i,jp,k,1))*divjp/(jac+jacjp)
            flxjm = flxjm
     .            + gamma*(pp0(i,j,k,1)+pp0(i,jm,k,1))*divjm/(jac+jacjm)

            flxkp = flxkp
     .            + gamma*(pp0(i,j,k,1)+pp0(i,j,kp,1))*divkp/(jac+jackp)
            flxkm = flxkm
     .            + gamma*(pp0(i,j,k,1)+pp0(i,j,km,1))*divkm/(jac+jackm)

            cov(1) = (flxip - flxim)/dxh(ig)
            cov(2) = (flxjp - flxjm)/dyh(jg)
            cov(3) = (flxkp - flxkm)/dzh(kg)

            !Transform cov to cnv
            call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3),.true.)

            psit = -dt*alpha**2*cnv*vol

            !P_si^B  ******************************

c$$$            ip = i+1
c$$$            im = i-1
c$$$            jp = j+1
c$$$            jm = j-1
c$$$            kp = k+1
c$$$            km = k-1

            !(j0 x a) part; a = curl(dv x B0); UPWIND for diagonal dominance
cc            hex = 1
cc            hey = 1
cc            hez = 1
cc            if (mgj0cnv(ijkg,1) > 0d0) hey = -1
cc            if (mgj0cnv(ijkg,2) > 0d0) hez = -1
cc            if (mgj0cnv(ijkg,3) > 0d0) hex = -1
cc            call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
cc     .                        ,a10,a20,a30,hex,hey,hez,igrid)
cc
cc            cov(1) = -mgj0cnv(ijkg,2)*a30 
cc            cov(2) = -mgj0cnv(ijkg,3)*a10 
cc            cov(3) = -mgj0cnv(ijkg,1)*a20 
cc
cc            hex = 1
cc            hey = 1
cc            hez = 1
cc            if (mgj0cnv(ijkg,1) > 0d0) hez = -1
cc            if (mgj0cnv(ijkg,2) > 0d0) hex = -1
cc            if (mgj0cnv(ijkg,3) > 0d0) hey = -1
cc            hex = 0
cc            hey = 0
cc            hez = 0
cc
cc            call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
cc     .                        ,a10,a20,a30,hex,hey,hez,igrid)
cc
cc            cov(1) = cov(1) + mgj0cnv(ijkg,3)*a20
cc            cov(2) = cov(2) + mgj0cnv(ijkg,1)*a30
cc            cov(3) = cov(3) + mgj0cnv(ijkg,2)*a10

            call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                        ,a10,a20,a30,0,igrid)

            cov(1) = mgj0cnv(ijkg,2)*a30 - mgj0cnv(ijkg,3)*a20
            cov(2) = mgj0cnv(ijkg,3)*a10 - mgj0cnv(ijkg,1)*a30
            cov(3) = mgj0cnv(ijkg,1)*a20 - mgj0cnv(ijkg,2)*a10

            !B0 x curl(a) part

            !!Find contravariant components of curl(dv x B0) at faces

            call find_curl_vxb(i ,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                         ,a1cnvip,a2cnvip,a3cnvip,1,igrid)
            call find_curl_vxb(im,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                         ,a1cnvim,a2cnvim,a3cnvim,1,igrid)

            call find_curl_vxb(i,j ,k,nxx,nyy,nzz,dv,b0_cnv
     .                         ,a1cnvjp,a2cnvjp,a3cnvjp,2,igrid)
            call find_curl_vxb(i,jm,k,nxx,nyy,nzz,dv,b0_cnv
     .                         ,a1cnvjm,a2cnvjm,a3cnvjm,2,igrid)

            call find_curl_vxb(i,j,k ,nxx,nyy,nzz,dv,b0_cnv
     .                         ,a1cnvkp,a2cnvkp,a3cnvkp,3,igrid)
            call find_curl_vxb(i,j,km,nxx,nyy,nzz,dv,b0_cnv
     .                         ,a1cnvkm,a2cnvkm,a3cnvkm,3,igrid)

            !!Find covariant components of curl(dv x B0) at faces
            xiph = (xip+x0)*0.5
            ximh = (xim+x0)*0.5
            xjph = (xjp+x0)*0.5
            xjmh = (xjm+x0)*0.5
            xkph = (xkp+x0)*0.5
            xkmh = (xkm+x0)*0.5

            yiph = (yip+y0)*0.5
            yimh = (yim+y0)*0.5
            yjph = (yjp+y0)*0.5
            yjmh = (yjm+y0)*0.5
            ykph = (ykp+y0)*0.5
            ykmh = (ykm+y0)*0.5

            ziph = (zip+z0)*0.5
            zimh = (zim+z0)*0.5
            zjph = (zjp+z0)*0.5
            zjmh = (zjm+z0)*0.5
            zkph = (zkp+z0)*0.5
            zkmh = (zkm+z0)*0.5

            call transformCurvToCurv(xiph,yiph,ziph,cartsn
     .                              ,a1covip,a2covip,a3covip
     .                              ,a1cnvip,a2cnvip,a3cnvip,.false.)
            call transformCurvToCurv(ximh,yimh,zimh,cartsn
     .                              ,a1covim,a2covim,a3covim
     .                              ,a1cnvim,a2cnvim,a3cnvim,.false.)

            call transformCurvToCurv(xjph,yjph,zjph,cartsn
     .                              ,a1covjp,a2covjp,a3covjp
     .                              ,a1cnvjp,a2cnvjp,a3cnvjp,.false.)
            call transformCurvToCurv(xjmh,yjmh,zjmh,cartsn
     .                              ,a1covjm,a2covjm,a3covjm
     .                              ,a1cnvjm,a2cnvjm,a3cnvjm,.false.)

            call transformCurvToCurv(xkph,ykph,zkph,cartsn
     .                              ,a1covkp,a2covkp,a3covkp
     .                              ,a1cnvkp,a2cnvkp,a3cnvkp,.false.)
            call transformCurvToCurv(xkmh,ykmh,zkmh,cartsn
     .                              ,a1covkm,a2covkm,a3covkm
     .                              ,a1cnvkm,a2cnvkm,a3cnvkm,.false.)

            !!Assemble B0 x curl(a) and j0 x a
            cov(1) =-cov(1)+b0_cnv(i,j,k,2)*( (a2covip-a2covim)/dxh(ig)
     .                                       -(a1covjp-a1covjm)/dyh(jg))
     .                     +b0_cnv(i,j,k,3)*( (a3covip-a3covim)/dxh(ig)
     .                                       -(a1covkp-a1covkm)/dzh(kg))

            cov(2) =-cov(2)+b0_cnv(i,j,k,1)*( (a1covjp-a1covjm)/dyh(jg)
     .                                       -(a2covip-a2covim)/dxh(ig))
     .                     +b0_cnv(i,j,k,3)*( (a3covjp-a3covjm)/dyh(jg)
     .                                       -(a2covkp-a2covkm)/dzh(kg))

            cov(3) =-cov(3)+b0_cnv(i,j,k,1)*( (a1covkp-a1covkm)/dzh(kg)
     .                                       -(a3covip-a3covim)/dxh(ig))
     .                     +b0_cnv(i,j,k,2)*( (a2covkp-a2covkm)/dzh(kg)
     .                                       -(a3covjp-a3covjm)/dyh(jg))

            !Introduce jacobian factor
            cov = cov/jac

            !Transform cov to cnv

            call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3),.true.)

            psib = -dt*alpha**2*cnv*vol

            !Add all contributions to form matvec ***********************

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = dv(i,j,k,ieq)/dt*vol*rho0(i,j,k,1)
     .                         + psiv(ieq) + psit(ieq) + psib(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(dv)

      nullify(pp0,rho0)
      nullify(b0_cnv)
      nullify(v0_cnv)

      contains

c     div_upwd
c     #####################################################################
      real(8) function div_upwd(half_elem)

        integer(4) :: half_elem

        integer(4) :: ip,im,jp,jm,kp,km
        real(8)    :: dxx,dyy,dzz,axp,axm,ayp,aym,azp,azm

c     Begin program

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        if (half_elem > 0) then
          dxx = dx(ig-1)
          dyy = dy(jg-1)
          dzz = dz(kg-1)
          ip = i
          jp = j
          kp = k
        elseif (half_elem < 0) then
          dxx = dx(ig)
          dyy = dy(jg)
          dzz = dz(kg)
          im = i
          jm = j
          km = k
        else
          div_upwd = 0d0
          return
        endif

        axp = dv(ip,j,k,1)*rho0(ip,j,k,1)
        axm = dv(im,j,k,1)*rho0(im,j,k,1)
        ayp = dv(i,jp,k,2)*rho0(i,jp,k,1)
        aym = dv(i,jm,k,2)*rho0(i,jm,k,1)
        azp = dv(i,j,kp,3)*rho0(i,j,kp,1)
        azm = dv(i,j,km,3)*rho0(i,j,km,1)

        div_upwd = ( (axp-axm)/dxx
     .              +(ayp-aym)/dyy
     .              +(azp-azm)/dzz )/jac
      
      end function div_upwd

      end subroutine v_mtvc
