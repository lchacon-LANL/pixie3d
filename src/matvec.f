c setMGBC
c####################################################################
      subroutine setMGBC(gpos,neq,nnx,nny,nnz,iig,array,bcnd,icomp)
cc     .                  ,array0)
c--------------------------------------------------------------------
c     Interfaces BC routines with MG code, for preconditioning.
c     NOTE: we assume that array, when representing a vector quantity,
c           is in the contravariant representation.
c     WARNING: if call sequence is modified, the INTERFACE block 
c           in module mg_internal needs to be updated accordingly.
c--------------------------------------------------------------------

      use grid

      use mg_internal

      use imposeBCinterface

      implicit none

c Call variables

      integer(4) :: nnx,nny,nnz,neq,bcnd(6,neq),iig,gpos

      real(8)    :: array(0:nnx+1,0:nny+1,0:nnz+1,neq)

      integer(4),optional,intent(IN) :: icomp

cc      real(8),optional,intent(IN) :: array0(0:nnx+1,0:nny+1,0:nnz+1,neq)

c Local variables

      integer(4) :: stencil_width,imn,imx,jmn,jmx,kmn,kmx

      integer(4),save :: ivar  !Needed to "remember' ivar from previous calls
                               !when icomp is not provided.

cc      real(8)    :: dummy(0:nnx+1,0:nny+1,0:nnz+1,neq)

c Begin program

c Consistency check

      if (PRESENT(icomp)) ivar = icomp

cc      if (PRESENT(array0)) then
cc        dummy = array0
cc      else
cc        dummy = 0d0
cc      endif

      if (     grid_params%nxv(iig) /= nnx
     .    .or. grid_params%nyv(iig) /= nny
     .    .or. grid_params%nzv(iig) /= nnz) then
        write (*,*) 'Grid sizes do not agree in setMGBC'
        write (*,*) 'Aborting...'
        stop
      endif

c Find grid node LOCAL position (if gpos > 0, else return local domain limits)

      call limits(gpos,nnx,nny,nnz,iig,imn,imx,jmn,jmx,kmn,kmx)

c Check LOCAL limits (return if not close to local domain boundaries)

      if (     (imn > 1 .and. imx < nnx)
     .    .and.(jmn > 1 .and. jmx < nny)
     .    .and.(kmn > 1 .and. kmx < nnz)) return

c Find GLOBAL limits (required for BC routine setBC)

      stencil_width = 1

      imn = max(imn-stencil_width,1)   + grid_params%ilo(iig)-1
      imx = min(imx+stencil_width,nnx) + grid_params%ilo(iig)-1
      jmn = max(jmn-stencil_width,1)   + grid_params%jlo(iig)-1
      jmx = min(jmx+stencil_width,nny) + grid_params%jlo(iig)-1
      kmn = max(kmn-stencil_width,1)   + grid_params%klo(iig)-1
      kmx = min(kmx+stencil_width,nnz) + grid_params%klo(iig)-1

c Select operation

      select case (neq)
      case(1)

cc        call setBC(IRHO,nnx,nny,nnz,array(:,:,:,neq),dummy(:,:,:,neq)
        call setBC(IRHO,nnx,nny,nnz,array(:,:,:,neq),zeros
     .            ,bcnd(:,neq),iig,iig,iig
     .            ,i1=imn,i2=imx,j1=jmn,j2=jmx,k1=kmn,k2=kmx)

cc        !Warning: Velocities for T BC unknown
cc        call setBC(ivar,nnx,nny,nnz,array(:,:,:,neq),dummy(:,:,:,neq)
cc     .            ,bcnd(:,1),iig,iig,iig
cc     .            ,i1=imn,i2=imx,j1=jmn,j2=jmx,k1=kmn,k2=kmx)

      case(3)

        allocate(v_cov (0:nnx+1,0:nny+1,0:nnz+1,neq))

        select case (ivar)
        case(IVX,IVY,IVZ)
cc          call setBC(IVX,nnx,nny,nnz,array,v_cov,dummy,bcnd
          call setBC(IVX,nnx,nny,nnz,array,v_cov,vzeros,bcnd
     .              ,iig,iig,iig
     .              ,i1=imn,i2=imx,j1=jmn,j2=jmx,k1=kmn,k2=kmx)
        case(IBX,IBY,IBZ)
cc          call setBC(IBX,nnx,nny,nnz,array,v_cov,dummy,bcnd
          call setBC(IBX,nnx,nny,nnz,array,v_cov,vzeros,bcnd
     .              ,iig,iig,iig
     .              ,i1=imn,i2=imx,j1=jmn,j2=jmx,k1=kmn,k2=kmx)
        case(IJX,IJY,IJZ)
cc          call setBC(IJX,nnx,nny,nnz,array,v_cov,dummy,bcnd
          call setBC(IJX,nnx,nny,nnz,array,v_cov,vzeros,bcnd
     .              ,iig,iig,iig
     .              ,i1=imn,i2=imx,j1=jmn,j2=jmx,k1=kmn,k2=kmx)
        case default
          write (*,*) 'Undefined variable in setMGBC'
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
      subroutine test_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine is a matvec test, y = A.x.
c     In call:
c      * gpos: vector index of position on the numerical grid
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D
c              , 27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c--------------------------------------------------------------------

      use matvec

      use precond_variables

      use mg_internal

      implicit none

c Call variables

      integer*4    neq,ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,nxx,nyy,nzz
      integer(4) ::  imin, imax, jmin, jmax, kmin, kmax
      integer(4) :: iimin,iimax,jjmin,jjmax,kkmin,kkmax

      real(8),allocatable,dimension(:,:,:,:) :: xarr
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,nu2

      logical    :: fpointers

c External

      real*8       vis
      external     vis

c Begin program

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

cc      nxx = MGgrid%nxv(igx)
cc      nyy = MGgrid%nyv(igy)
cc      nzz = MGgrid%nzv(igz)
      nxx = grid_params%nxv(igrid)
      nyy = grid_params%nyv(igrid)
      nzz = grid_params%nzv(igrid)

c Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz,igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(xarr(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      xarr = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,xarr,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,xarr,bcnd
     .            ,icomp=IRHO)

c Map velocity components

      v0_cnv => gv0%grid(igrid)%array

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jac = gmetric%grid(igrid)%jac(i,j,k)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)
            ijkg   = ijk + isig - 1

cc            upwind = .5*(v0_cnv(i,j,k,1)+abs(v0_cnv(i,j,k,1)))
cc     .                 *( xarr(i  ,j,k,1) - xarr(i-1,j,k,1) )/dx(ig-1)
cc     .              +.5*(v0_cnv(i,j,k,1)-abs(v0_cnv(i,j,k,1)))
cc     .                 *( xarr(i+1,j,k,1) - xarr(i  ,j,k,1) )/dx(ig)
cc     .              +.5*(v0_cnv(i,j,k,2)+abs(v0_cnv(i,j,k,2)))
cc     .                 *( xarr(i,j  ,k,1) - xarr(i,j-1,k,1) )/dy(jg-1)
cc     .              +.5*(v0_cnv(i,j,k,2)-abs(v0_cnv(i,j,k,2)))
cc     .                 *( xarr(i,j+1,k,1) - xarr(i,j  ,k,1) )/dy(jg)
cc     .              +.5*(v0_cnv(i,j,k,3)+abs(v0_cnv(i,j,k,3)))
cc     .                 *( xarr(i,j,k  ,1) - xarr(i,j,k-1,1) )/dz(kg-1)
cc     .              +.5*(v0_cnv(i,j,k,3)-abs(v0_cnv(i,j,k,3)))
cc     .                 *( xarr(i,j,k+1,1) - xarr(i,j,k  ,1) )/dz(kg)
cc
cc            upwind = upwind/jac
cc            nu2 = 0.1

cc            y(ijk) = (1./dt + alpha*(gamma-1.)*mgdivV0(ijkg))*x(ijk)
cc     .              + alpha*upwind 
cccc     .              - alpha*nu2*laplacian(i,j,k,nxx,nyy,nzz,xarr)

cc            y(ijk) = y(ijk)*volume(i,j,k,igrid,igrid,igrid)

            y(ijk) = laplacian(i,j,k,nxx,nyy,nzz,xarr) !Laplacian includes volume factor

          enddo
        enddo
      enddo

c End program

      deallocate(xarr)

      call deallocPointers(fpointers)

      end subroutine test_mtvc

c tmp_mtvc
c####################################################################
      subroutine tmp_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
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
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c--------------------------------------------------------------------

      use matvec

      use precond_variables

      use mg_internal

      implicit none

c Call variables

      integer*4    neq,ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,nnx,nny,nnz
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: dtmp
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,nu2

      logical    :: fpointers

c External

      real*8       vis
      external     vis

c Begin program

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

cc      nx = MGgrid%nxv(igx)
cc      ny = MGgrid%nyv(igy)
cc      nz = MGgrid%nzv(igz)
      nnx = grid_params%nxv(igrid)
      nny = grid_params%nyv(igrid)
      nnz = grid_params%nzv(igrid)

c Find limits for loops

      call limits(abs(gpos),nnx,nny,nnz,igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(dtmp(0:nnx+1,0:nny+1,0:nnz+1,neq))

      dtmp = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nnx,nny,nnz,dtmp,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nnx,nny,nnz,igrid,dtmp,bcnd
     .            ,icomp=ITMP)

c Map velocity components

      v0_cnv => gv0%grid(igrid)%array

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jac = gmetric%grid(igrid)%jac(i,j,k)

            ijk    = i + nnx*(j-1) + nnx*nny*(k-1)
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

            y(ijk) = ((1./dt + alpha*(gamma-1.)*mgdivV0(ijkg))*x(ijk)
     .            + alpha*upwind )*volume(i,j,k,igrid,igrid,igrid)
     .            - alpha*chi*laplacian(i,j,k,nnx,nny,nnz,dtmp(:,:,:,1))

          enddo
        enddo
      enddo

c End program

      deallocate(dtmp)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine tmp_mtvc

c rho_mtvc
c####################################################################
      subroutine rho_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
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
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c--------------------------------------------------------------------

      use matvec

      use precond_variables

      use mg_internal

      implicit none

c Call variables

      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,*)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,nnx,nny,nnz
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: drho
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,nu2

      logical    :: fpointers

c External

      real(8)    :: vis
      external   :: vis

c Begin program

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

cc      nx = MGgrid%nxv(igx)
cc      ny = MGgrid%nyv(igy)
cc      nz = MGgrid%nzv(igz)
      nnx = grid_params%nxv(igrid)
      nny = grid_params%nyv(igrid)
      nnz = grid_params%nzv(igrid)

c Find limits for loops

      call limits(abs(gpos),nnx,nny,nnz,igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(drho(0:nnx+1,0:nny+1,0:nnz+1,neq))

      drho = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped
      !For finding the diagonal, gpos > 0
      call mapMGVectorToArray(max(0,gpos),neq,x,nnx,nny,nnz,drho,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nnx,nny,nnz,igrid,drho,bcnd
     .            ,icomp=IRHO)

c Map velocity components

      v0_cnv => gv0%grid(igrid)%array

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jac = gmetric%grid(igrid)%jac(i,j,k)

            ijk    = i + nnx*(j-1) + nnx*nny*(k-1)
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
     .            + alpha*upwind )*volume(i,j,k,igrid,igrid,igrid)
     .            - alpha*dd*laplacian(i,j,k,nnx,nny,nnz,drho(:,:,:,1))

          enddo
        enddo
      enddo

c End program

      deallocate(drho)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine rho_mtvc

c b_mtvc
c####################################################################
      subroutine b_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
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
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c--------------------------------------------------------------------

      use matvec

      use precond_variables

      use mg_internal

      implicit none

c Call variables

      integer*4    neq,ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

      real(8)    :: upwind,etal,flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol

      real(8),allocatable,dimension(:,:,:,:) :: db
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      logical    :: fpointers

c External

      real(8)    :: res
      external   :: res

c Begin program

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

cc      nxx = MGgrid%nxv(igx)
cc      nyy = MGgrid%nyv(igy)
cc      nzz = MGgrid%nzv(igz)
      nxx = grid_params%nxv(igrid)
      nyy = grid_params%nyv(igrid)
      nzz = grid_params%nzv(igrid)

c Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz
     .           ,igrid,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(db(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      db = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,db,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,db,bcnd
     .            ,icomp=IBX)

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

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jacip  = gmetric%grid(igrid)%jac(ip,j,k)
            jacim  = gmetric%grid(igrid)%jac(im,j,k)
            jacjp  = gmetric%grid(igrid)%jac(i,jp,k)
            jacjm  = gmetric%grid(igrid)%jac(i,jm,k)
            jackp  = gmetric%grid(igrid)%jac(i,j,kp)
            jackm  = gmetric%grid(igrid)%jac(i,j,km)
            jac    = gmetric%grid(igrid)%jac(i,j,k)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            vol    = volume(i,j,k,igrid,igrid,igrid)

            etal   = res(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid)

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
cc     .                      *veclaplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,1)
cc     .                                                     ,db(:,:,:,2)
cc     .                                                     ,db(:,:,:,3)
cc     .                                                 ,ones,.false.,1)
     .                      *laplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,1))

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
cc     .                      *veclaplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,1)
cc     .                                                     ,db(:,:,:,2)
cc     .                                                     ,db(:,:,:,3)
cc     .                                                 ,ones,.false.,2)
     .                      *laplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,2))

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
cc     .                      *veclaplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,1)
cc     .                                                     ,db(:,:,:,2)
cc     .                                                     ,db(:,:,:,3)
cc     .                                                 ,ones,.false.,3)
     .                      *laplacian(i,j,k,nxx,nyy,nzz,db(:,:,:,3))

          enddo
        enddo
      enddo

c End program

      deallocate(db)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine b_mtvc

c v_mtvc
c####################################################################
      subroutine v_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
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
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c--------------------------------------------------------------------

      use matvec

      use precond_variables

      use mg_internal

      use nlfunction_setup

      implicit none

c Call variables

      integer*4    neq,ntot,igrid,gpos,bcnd(6,*)
      real*8       x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ip,im,jp,jm,kp,km,hex,hey,hez
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

      logical    :: fpointers

c External

      real(8)    :: vis
      external   :: vis

c Begin program

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

cc      nxx = MGgrid%nxv(igrid)
cc      nyy = MGgrid%nyv(igrid)
cc      nzz = MGgrid%nzv(igrid)
      nxx = grid_params%nxv(igrid)
      nyy = grid_params%nyv(igrid)
      nzz = grid_params%nzv(igrid)

c Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz,igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(dv(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      dv = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
     .            ,icomp=IVX)

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

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jacip  = gmetric%grid(igrid)%jac(ip,j,k)
            jacim  = gmetric%grid(igrid)%jac(im,j,k)
            jacjp  = gmetric%grid(igrid)%jac(i,jp,k)
            jacjm  = gmetric%grid(igrid)%jac(i,jm,k)
            jackp  = gmetric%grid(igrid)%jac(i,j,kp)
            jackm  = gmetric%grid(igrid)%jac(i,j,km)
            jac    = gmetric%grid(igrid)%jac(i,j,k)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            vol    = volume(i,j,k,igrid,igrid,igrid)

            mul    = vis(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid)

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
cc     .                      *veclaplacian(i,j,k,nxx,nyy,nzz,dv(:,:,:,1)
cc     .                                                     ,dv(:,:,:,2)
cc     .                                                     ,dv(:,:,:,3)
cc     .                                   ,ones,alt_eom,ieq)
     .                      *laplacian(i,j,k,nxx,nyy,nzz,dv(:,:,:,ieq))
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

            !Fluxes at faces for calculation of grad(gamma*p0*div(dv))

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

            psit = -dt*alpha**2*cov*vol

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

            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                              ,a1covip,a2covip,a3covip
     .                              ,a1cnvip,a2cnvip,a3cnvip,.false.)
            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                              ,a1covim,a2covim,a3covim
     .                              ,a1cnvim,a2cnvim,a3cnvim,.false.)

            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                              ,a1covjp,a2covjp,a3covjp
     .                              ,a1cnvjp,a2cnvjp,a3cnvjp,.false.)
            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                              ,a1covjm,a2covjm,a3covjm
     .                              ,a1cnvjm,a2cnvjm,a3cnvjm,.false.)

            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                              ,a1covkp,a2covkp,a3covkp
     .                              ,a1cnvkp,a2cnvkp,a3cnvkp,.false.)
            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
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

            psib = -dt*alpha**2*cov*vol

            !Add all contributions to form matvec ***********************

            !Transform cov to cnv

            cov = psib + psit

            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3),.true.)

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = dv(i,j,k,ieq)/dt*vol*rho0(i,j,k,1)
     .                            + psiv(ieq) + cnv(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(dv)

      nullify(pp0,rho0)
      nullify(b0_cnv)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

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
