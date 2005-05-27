c setMGBC
c####################################################################
      subroutine setMGBC(gpos,neq,nnx,nny,nnz,iig,array,bcnd,arr_cov
     .                  ,arr0,icomp,is_cnv,is_vec,result_is_vec,iorder)
c--------------------------------------------------------------------
c     Interfaces BC routines with MG code, for preconditioning.
c     On input:
c       *gpos: node number (whole grid if gpos=0)
c       *neq,nnx,nny,nnz: number of degrees of freedom.
c       *iig: grid level
c       *array: array containing values of magnitudes for all grid
c               points. The magnitudes may be scalars or vector
c               components (the latter assumed contravariant by default)
c       *bcnd: array containing BC information
c
c     Optional parameters:
c       *arr_cov: covariant components of vector components.
c       *arr0: array containing dirichlet BCs (zero by default)
c       *icomp: integer specifying magnitude of interest (if not
c               provided, use last used value).
c       *is_cnv: logical variable indicating whether array
c                contains covariant or contravariant components
c                (if not provided, use last used value).
c       *is_vec: logical variable specifying whether array
c          represents a curvilinear vector or a set of scalars
c          (if not provided, use last used value):
c             -If is_vec=.true., array is taken as curvilinear vector
c              components (covariant or contravariant, according
c              to variable is_cnv).
c             -If is_vec=.false., array is taken as coupled scalars
c              (also Cartesian vector components).
c       *result_is_vec: logical variable specifying whether output
c                       in array should be curvilinear (.true.) or
c                       not (default=is_vec)
c       *iorder: order of inter/extrapolation in BC routines (linear
c                by default).
c
c     On output, array is returned in the same representation.
c
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

      real(8),optional,intent(INOUT) ::
     .                         arr_cov(0:nnx+1,0:nny+1,0:nnz+1,neq)
      real(8),optional,intent(IN) ::
     .                         arr0   (0:nnx+1,0:nny+1,0:nnz+1,neq)

      integer(4),optional,intent(IN) :: icomp,iorder

      logical   ,optional,intent(IN) :: is_cnv,is_vec,result_is_vec

c Local variables

      integer(4) :: i,j,k,stencil_width,order
     .             ,imng,imxg,jmng,jmxg,kmng,kmxg
     .             ,imnl,imxl,jmnl,jmxl,kmnl,kmxl

      integer(4),save :: ivar  !Needed to "remember' ivar from previous calls
                               !when icomp is not provided.

      logical,save    :: iscnv,isvec !Needed to "remember' these from previous calls
                                     !when they are not provided.

      logical         :: res_is_vec

      real(8)    :: car(3),curv(3)

c Begin program

c Optional arguments

      if (PRESENT(icomp))  ivar  = icomp

      if (PRESENT(is_cnv)) iscnv = is_cnv

      if (PRESENT(is_vec)) isvec = is_vec

      if (PRESENT(result_is_vec)) then
        res_is_vec = result_is_vec
      else
        res_is_vec = isvec
      endif

      if (PRESENT(iorder)) then
        order = iorder
      else
        order = 2
      endif

c Consistency check

      if (     grid_params%nxv(iig) /= nnx
     .    .or. grid_params%nyv(iig) /= nny
     .    .or. grid_params%nzv(iig) /= nnz) then
        write (*,*) 'Grid sizes do not agree in setMGBC'
        write (*,*) 'Aborting...'
        stop
      endif

      if (.not.isvec) iscnv=.true.

c Find grid node LOCAL position (if gpos > 0, else return local domain limits)

      call limits(gpos,nnx,nny,nnz,iig,imnl,imxl,jmnl,jmxl,kmnl,kmxl)

c Check LOCAL limits (return if not close to local domain boundaries)

      if (     (imnl > 1 .and. imxl < nnx)
     .    .and.(jmnl > 1 .and. jmxl < nny)
     .    .and.(kmnl > 1 .and. kmxl < nnz)) return

c Find GLOBAL limits (required for BC routine setBC)

      stencil_width = 1

      imng = max(imnl-stencil_width,1)   + grid_params%ilo(iig)-1
      imxg = min(imxl+stencil_width,nnx) + grid_params%ilo(iig)-1
      jmng = max(jmnl-stencil_width,1)   + grid_params%jlo(iig)-1
      jmxg = min(jmxl+stencil_width,nny) + grid_params%jlo(iig)-1
      kmng = max(kmnl-stencil_width,1)   + grid_params%klo(iig)-1
      kmxg = min(kmxl+stencil_width,nnz) + grid_params%klo(iig)-1

c Select operation

      select case (neq)
      case(1)

        allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,1))

        if (PRESENT(arr0)) then
          v0 = arr0
        else
          v0 = 0d0
        endif

        call setBC(IRHO,nnx,nny,nnz,array(:,:,:,1),v0(:,:,:,1)
     .            ,bcnd(:,neq),iig,iig,iig
     .            ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .            ,iorder=order)

cc        !Warning: Velocities for T BC unknown
cc        call setBC(ivar,nnx,nny,nnz,array(:,:,:,neq),v0(:,:,:,neq)
cc     .            ,bcnd(:,1),iig,iig,iig
cc     .            ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg)
cc
        deallocate(v0)

      case(3)

        allocate(v_cnv (0:nnx+1,0:nny+1,0:nnz+1,neq)
     .          ,v_cov (0:nnx+1,0:nny+1,0:nnz+1,neq)
     .          ,v0    (0:nnx+1,0:nny+1,0:nnz+1,neq))

        v_cnv = 0d0
        v_cov = 0d0

        if (PRESENT(arr0)) then
          v0 = arr0
        else
          v0 = 0d0
        endif

        if (isvec) then

          if (PRESENT(arr_cov)) then
            v_cnv = array
            v_cov = arr_cov
          else
            if (iscnv) then
              v_cnv = array
            else
              v_cov = array
            endif
          endif

        else

          do k = kmnl-1,kmxl+1
            do j = jmnl-1,jmxl+1
              do i = imnl-1,imxl+1
                car = array(i,j,k,:)
                call transformVectorToCurvilinear
     .               (i,j,k,iig,iig,iig
     .               ,car(1),car(2),car(3)
     .               ,.not.iscnv
     .               ,array(i,j,k,1),array(i,j,k,2),array(i,j,k,3))
              enddo
            enddo
          enddo

          if (iscnv) then
            v_cnv = array
          else
            v_cov = array
          endif

        endif

        call setBC(ivar,neq,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd
     .            ,iig,iig,iig
     .            ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .            ,is_cnv=iscnv,is_vec=.true.,iorder=order)

        if (isvec .or. res_is_vec) then

          if (PRESENT(arr_cov)) then
            array   = v_cnv
            arr_cov = v_cov
          else
            if (iscnv) then
              array = v_cnv
            else
              array = v_cov
            endif
          endif

        else

          do k = kmnl-1,kmxl+1
            do j = jmnl-1,jmxl+1
              do i = imnl-1,imxl+1
                curv = v_cnv(i,j,k,:)
                call transformVectorToCartesian
     .               (i,j,k,iig,iig,iig
     .               ,curv(1),curv(2),curv(3)
     .               ,.false.
     .               ,array(i,j,k,1),array(i,j,k,2),array(i,j,k,3))
              enddo
            enddo
          enddo

        endif

        deallocate(v_cov,v_cnv,v0)

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

      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,nxx,nyy,nzz
      integer(4) ::  imin, imax, jmin, jmax, kmin, kmax
      integer(4) :: iimin,iimax,jjmin,jjmax,kkmin,kkmax

      real(8),allocatable,dimension(:,:,:,:) :: xarr
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,nu2

      logical    :: fpointers

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
cccc     .              - alpha*nu2*laplacian(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,xarr)

cc            y(ijk) = y(ijk)*volume(i,j,k,igrid,igrid,igrid)

            y(ijk) = laplacian(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,xarr
     .                        ,vol=vol_wgt)

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

      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,nnx,nny,nnz
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: dtmp
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,nu2,dvol

      logical    :: fpointers

c Begin program

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

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

            if (vol_wgt) then
              dvol = volume(i,j,k,igrid,igrid,igrid)
            else
              dvol = 1d0
            endif

            y(ijk) = ((1./dt + alpha*(gamma-1.)*mgdivV0(ijkg))*x(ijk)
     .            + alpha*upwind )*dvol
     .            - alpha*chi
     .                 *laplacian(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                           ,dtmp(:,:,:,1),vol=vol_wgt)

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

      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ijk,ijkg,nnx,nny,nnz,ip,im,jp,jm,kp,km
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
 
      real(8),allocatable,dimension(:,:,:,:) :: drho
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,lap,dvol,flxip,flxim,flxjp,flxjm,flxkp,flxkm

      logical    :: fpointers

c Begin program

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

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

            ijk = i + nnx*(j-1) + nnx*nny*(k-1)

            if (vol_wgt) then
              dvol = volume(i,j,k,igrid,igrid,igrid)
            else
              dvol = 1d0
            endif

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
     .            + alpha*upwind )*dvol
     .            - alpha*dd
     .                *laplacian(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                          ,drho(:,:,:,1),vol=vol_wgt)

cc            ip = i+1
cc            im = i-1
cc            jp = j+1
cc            jm = j-1
cc            kp = k+1
cc            km = k-1
cc
cc            jacip  = gmetric%grid(igrid)%jac(ip,j,k)
cc            jacim  = gmetric%grid(igrid)%jac(im,j,k)
cc            jacjp  = gmetric%grid(igrid)%jac(i,jp,k)
cc            jacjm  = gmetric%grid(igrid)%jac(i,jm,k)
cc            jackp  = gmetric%grid(igrid)%jac(i,j,kp)
cc            jackm  = gmetric%grid(igrid)%jac(i,j,km)
cc
cc            if (isSP(i,j,k,igrid,igrid,igrid)) then
cc              flxip = 0.125*(jac+jacip)*(
cc     .         (    (v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip)
cc     .          +abs(v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip) )
cc     .                                                  *drho(i ,j,k,1)
cc     .        +(    (v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip)          
cc     .          -abs(v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip) )
cc     .                                                  *drho(ip,j,k,1))
cc            else
cc              flxip = 0.25*(
cc     .         (    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))
cc     .          +abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*drho(i ,j,k,1)
cc     .        +(    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))          
cc     .          -abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*drho(ip,j,k,1))
cc
cc            endif
cc            if (isSP(im,j,k,igrid,igrid,igrid)) then
cc              flxim = 0.125*(jac+jacim)*(
cc     .         (    (v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim)
cc     .          +abs(v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim) )
cc     .                                                  *drho(i ,j,k,1)
cc     .        +(    (v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim)          
cc     .          -abs(v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim) )
cc     .                                                  *drho(im,j,k,1))
cc            elseif (isSP(i,j,k,igrid,igrid,igrid)) then
cc              flxim = 0.5*(
cc     .         (    (v0_cnv(im,j,k,1))
cc     .          +abs(v0_cnv(im,j,k,1)) )*drho(i ,j,k,1)
cc     .        +(    (v0_cnv(im,j,k,1))          
cc     .          -abs(v0_cnv(im,j,k,1)) )*drho(im,j,k,1))
cc            else
cc              flxim = 0.25*(
cc     .         (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
cc     .          +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*drho(i ,j,k,1)
cc     .        +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
cc     .          -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*drho(im,j,k,1))
cc            endif
cc
cc            flxjp = 0.25*(
cc     .         (    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))
cc     .          +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*drho(i,j ,k,1)
cc     .        +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))          
cc     .          -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*drho(i,jp,k,1))
cc            flxjm = 0.25*(
cc     .         (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
cc     .          +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*drho(i,j ,k,1)
cc     .        +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
cc     .          -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*drho(i,jm,k,1))
cc
cc            flxkp = 0.25*(
cc     .         (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))
cc     .          +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*drho(i,j,k ,1)
cc     .        +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))              
cc     .          -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*drho(i,j,kp,1))
cc            flxkm = 0.25*(
cc     .         (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
cc     .          +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*drho(i,j,k ,1)
cc     .        +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))              
cc     .          -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*drho(i,j,km,1))
cc
cc            upwind =( (flxip-flxim)/dxh(ig)
cc     .               +(flxjp-flxjm)/dyh(jg)
cc     .               +(flxkp-flxkm)/dzh(kg) )/jac
cc
cc            lap    = laplacian(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
cc     .                        ,drho(:,:,:,1),vol=vol_wgt)
cc
cc            y(ijk) = (drho(i,j,k,1)/dt+alpha*upwind)*dvol - alpha*dd*lap

          enddo
        enddo
      enddo

c End program

      deallocate(drho)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine rho_mtvc

ccc b_mtvc
ccc####################################################################
cc      subroutine b_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
ccc--------------------------------------------------------------------
ccc     This subroutine calculates, for given x, y = A(psi)x  matrix-free
ccc     for Faraday's law.
ccc     In call:
ccc      * gpos: vector index of position on the numerical grid
ccc            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
ccc              surrounding stencil is filled (9-pt stencil in 2D
ccc              , 27-pt stencil in 3D).
ccc            + If gpos = 0, all the grid is considered.
ccc            + If gpos < 0, all grid is mapped, but operations are 
ccc              restricted to stencil of abs(gpos) (useful for
ccc              matrix-light GS)
ccc      * neq: number of coupled equations
ccc      * ntot: total number of unknowns: neq*nx*ny*nz
ccc      * x(ntot): input vector
ccc      * y(ntot): output vector
ccc      * igrid: grid level
ccc      * bcnf: boundary conditions on x vector.
ccc--------------------------------------------------------------------
cc
cc      use matvec
cc
cc      use precond_variables
cc
cc      use mg_internal
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq)
cc      real(8)    :: x(ntot),y(ntot)
cc
ccc Local variables
cc
cc      integer(4) :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
cc      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
cc      integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
cc
cc      real(8)    :: upwind,etal,flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol
cc     .             ,veclap(3),cnv(3),cov(3),car(3)
cc
cc      real(8),allocatable,dimension(:,:,:,:) :: db
cc      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv
cc
cc      logical    :: fpointers,alt_eom_b,is_cnv
cc
ccc Begin program
cc
cc      is_cnv = .true.
cc
cc      alt_eom_b = .false.
cc
cc      call allocPointers(neq,fpointers)
cc
cc      isig = MGgrid%istartp(igrid)
cc
cc      nxx = grid_params%nxv(igrid)
cc      nyy = grid_params%nyv(igrid)
cc      nzz = grid_params%nzv(igrid)
cc
cc      igx = igrid
cc      igy = igrid
cc      igz = igrid
cc
ccc Find limits for loops
cc
cc      call limits(abs(gpos),nxx,nyy,nzz
cc     .           ,igrid,imin,imax,jmin,jmax,kmin,kmax)
cc
ccc Map vector x to array for processing
cc
cc      allocate(db(0:nxx+1,0:nyy+1,0:nzz+1,neq))
cc
cc      db = 0d0
cc
cc      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
cc      !For finding the diagonal, gpos > 0
cc      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,db,igrid
cc     .                       ,.false.)
cc
cc      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,db,bcnd
cc     .            ,icomp=IBX,is_cnv=is_cnv,is_vec=.true.)
cccc     .            ,icomp=IBX,is_cnv=is_cnv,is_vec=.true.,iorder=2)
cc
ccc Velocity field (including BCs)
cc
cc      v0_cnv => gv0%grid(igrid)%array
cc
ccc Calculate matrix-vector product
cc
cc      do k = kmin,kmax
cc        do j = jmin,jmax
cc          do i = imin,imax
cc
cc            ip = i+1
cc            im = i-1
cc            jp = j+1
cc            jm = j-1
cc            kp = k+1
cc            km = k-1
cc
cccc            if (isSP(i,j,k,igrid,igrid,igrid)) im = i
cc
cc            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
cc
cc            jacip  = gmetric%grid(igrid)%jac(ip,j,k)
cc            jacim  = gmetric%grid(igrid)%jac(im,j,k)
cc            jacjp  = gmetric%grid(igrid)%jac(i,jp,k)
cc            jacjm  = gmetric%grid(igrid)%jac(i,jm,k)
cc            jackp  = gmetric%grid(igrid)%jac(i,j,kp)
cc            jackm  = gmetric%grid(igrid)%jac(i,j,km)
cc            jac    = gmetric%grid(igrid)%jac(i,j,k)
cc
cc            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + 1d-10
cc
cc            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)
cc
cc            if (vol_wgt) then
cc              vol = volume(i,j,k,igrid,igrid,igrid)
cc            else
cc              vol = 1d0
cc            endif
cc
cc            veclap = veclaplacian(i,j,k,nxx,nyy,nzz
cc     .                           ,igrid,igrid,igrid,db
cc     .                           ,alt_eom_b,vol=.false.)
cc
cc            etal   = res(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid)
cc
cccc            call find_curl_vxb(i,j,k,nxx,nyy,nzz,v0_cnv,db
cccc     .                      ,cnv(1),cnv(2),cnv(3),0,igrid)
cccc
cccc            cnv = (db(i,j,k,:)/dt + alpha*cnv - alpha*etal*veclap)*vol
cccc
cccc            do ieq=1,3
cccc              y(neq*(ijk-1)+ieq) = cnv(ieq)
cccc            enddo
cccc            cycle
cc
cc            flxjp = 0.5/(jac+jacjp)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,j ,k,1)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))         
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,jp,k,1))
cc     .             -0.5/(jac+jacjp)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1)) )*db(i,j ,k,2)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1)) )*db(i,jp,k,2))
cccc            flxjm = 0.5/(jac+jacjm)*(
cccc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
cccc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,j ,k,1)
cccc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
cccc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,jm,k,1))
cccc     .             -0.5/(jac+jacjm)*(
cccc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1))
cccc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1)) )*db(i,j ,k,2)
cccc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1))          
cccc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1)) )*db(i,jm,k,2))
cc            flxjm = 0.5/(jac+jacjm)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,jm,k,1)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,j ,k,1))
cc     .             -0.5/(jac+jacjm)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1)) )*db(i,jm,k,2)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1)) )*db(i,j ,k,2))
cc
cc            flxkp = 0.5/(jac+jackp)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,k ,1)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,kp,1))
cc     .             -0.5/(jac+jackp)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1)) )*db(i,j,k ,3)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1)) )*db(i,j,kp,3))
cccc            flxkm = 0.5/(jac+jackm)*(
cccc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
cccc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,k ,1)
cccc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))          
cccc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,km,1))
cccc     .             -0.5/(jac+jackm)*(
cccc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1))
cccc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1)) )*db(i,j,k ,3)
cccc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1))          
cccc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1)) )*db(i,j,km,3))
cc            flxkm = 0.5/(jac+jackm)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,km,1)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,k ,1))
cc     .             -0.5/(jac+jackm)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1)) )*db(i,j,km,3)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1)) )*db(i,j,k ,3))
cc
cc            upwind =  (flxjp-flxjm)/dyh(jg)
cc     .               +(flxkp-flxkm)/dzh(kg)
cc
cc            cnv(1) = ( db(i,j,k,1)/dt + alpha*upwind 
cc     .                - alpha*etal*veclap(1))*vol
cc
cccc            if (isSP(i,j,k,igrid,igrid,igrid)) then
cccc              flxip = 0.25*(
cccc     .           (    (v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip)
cccc     .            +abs(v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip) )
cccc     .                                              *db(i ,j,k,2)
cccc     .          +(    (v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip)          
cccc     .            -abs(v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip) )
cccc     .                                              *db(ip,j,k,2))
cccc     .               -0.25*(
cccc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2))
cccc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2)) )
cccc     .                                              *db(i ,j,k,1)/jac
cccc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2))          
cccc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2)) )
cccc     .                                              *db(ip,j,k,1)/jacip)
cccc            else
cc              flxip = 0.5/(jac+jacip)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(i ,j,k,2)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(ip,j,k,2))
cc     .             -0.5/(jac+jacip)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2)) )*db(i ,j,k,1)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2)) )*db(ip,j,k,1))
cccc            endif
cc
cccc            if (isSP(im,j,k,igrid,igrid,igrid)) then
cccccc              flxim = 0.25*(
cccccc     .          (    (v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim)
cccccc     .           +abs(v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim) )
cccccc     .                                             *db(i ,j,k,2)
cccccc     .         +(    (v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim)          
cccccc     .           -abs(v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim) )
cccccc     .                                             *db(im,j,k,2)      )
cccccc     .              -0.25*(
cccccc     .          (    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))
cccccc     .           +abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )
cccccc     .                                             *db(i ,j,k,1)/jac
cccccc     .         +(    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))          
cccccc     .           -abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )
cccccc     .                                             *db(im,j,k,1)/jacim)
cccc              flxim = 0.25*(
cccc     .          (    (v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim)
cccc     .           +abs(v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim) )
cccc     .                                             *db(im,j,k,2)
cccc     .         +(    (v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim)          
cccc     .           -abs(v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim) )
cccc     .                                             *db(i ,j,k,2)      )
cccc     .              -0.25*(
cccc     .          (    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))
cccc     .           +abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )
cccc     .                                             *db(im,j,k,1)/jac
cccc     .         +(    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))          
cccc     .           -abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )
cccc     .                                             *db(i ,j,k,1)/jacim)
cccc            elseif (isSP(i,j,k,igrid,igrid,igrid)) then
cccccc              flxim = 0.5/jacim*(
cccccc     .           (    (v0_cnv(im,j,k,1))
cccccc     .            +abs(v0_cnv(im,j,k,1)) )*db(i ,j,k,2)
cccccc     .          +(    (v0_cnv(im,j,k,1))          
cccccc     .            -abs(v0_cnv(im,j,k,1)) )*db(im,j,k,2))
cccccc     .               -0.5*(
cccccc     .           (    (v0_cnv(im,j,k,2))
cccccc     .            +abs(v0_cnv(im,j,k,2)) )*db(i ,j,k,1)/jac
cccccc     .          +(    (v0_cnv(im,j,k,2))          
cccccc     .            -abs(v0_cnv(im,j,k,2)) )*db(im,j,k,1)/jacim)
cccc              flxim = 0.5/jacim*(
cccc     .           (    (v0_cnv(im,j,k,1))
cccc     .            +abs(v0_cnv(im,j,k,1)) )*db(im,j,k,2)
cccc     .          +(    (v0_cnv(im,j,k,1))          
cccc     .            -abs(v0_cnv(im,j,k,1)) )*db(i ,j,k,2))
cccc     .               -0.5*(
cccc     .           (    (v0_cnv(im,j,k,2))
cccc     .            +abs(v0_cnv(im,j,k,2)) )*db(im,j,k,1)/jac
cccc     .          +(    (v0_cnv(im,j,k,2))          
cccc     .            -abs(v0_cnv(im,j,k,2)) )*db(i ,j,k,1)/jacim)
cccc            else
cccc              flxim = 0.5/(jac+jacim)*(
cccc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
cccc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(i ,j,k,2)
cccc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
cccc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(im,j,k,2))
cccc     .               -0.5/(jac+jacim)*(
cccc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))
cccc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )*db(i ,j,k,1)
cccc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))          
cccc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )*db(im,j,k,1))
cc              flxim = 0.5/(jac+jacim)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(im,j,k,2)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(i ,j,k,2))
cc     .               -0.5/(jac+jacim)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )*db(im,j,k,1)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )*db(i ,j,k,1))
cccc            endif
cc
cc            flxkp = 0.5/(jac+jackp)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,k ,2)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,kp,2))
cc     .             -0.5/(jac+jackp)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2)) )*db(i,j,k ,3)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2)) )*db(i,j,kp,3))
cccc            flxkm = 0.5/(jac+jackm)*(
cccc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
cccc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,k ,2)
cccc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))          
cccc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,km,2))
cccc     .             -0.5/(jac+jackm)*(
cccc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2))
cccc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2)) )*db(i,j,k ,3)
cccc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2))          
cccc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2)) )*db(i,j,km,3))
cc            flxkm = 0.5/(jac+jackm)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,km,2)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,k ,2))
cc     .             -0.5/(jac+jackm)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2)) )*db(i,j,km,3)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2)) )*db(i,j,k ,3))
cc
cc            upwind =  (flxip-flxim)/dxh(ig)
cc     .               +(flxkp-flxkm)/dzh(kg)
cc
cc            cnv(2) = ( db(i,j,k,2)/dt + alpha*upwind 
cc     .                - alpha*etal*veclap(2) )*vol
cc
cccc            if (isSP(i,j,k,igrid,igrid,igrid)) then
cccc              flxip = 0.125*(jac+jacip)*(
cccc     .           (    (v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip)
cccc     .            +abs(v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip) )
cccc     .                                              *db(i ,j,k,3)/jac
cccc     .          +(    (v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip)          
cccc     .            -abs(v0_cnv(i,j,k,1)/jac+v0_cnv(ip,j,k,1)/jacip) )
cccc     .                                              *db(ip,j,k,3)/jacip)
cccc     .               -0.125*(jac+jacip)*(
cccc     .           (    (v0_cnv(i,j,k,3)/jac+v0_cnv(ip,j,k,3)/jacip)
cccc     .            +abs(v0_cnv(i,j,k,3)/jac+v0_cnv(ip,j,k,3)/jacip) )
cccc     .                                              *db(i ,j,k,1)/jac
cccc     .          +(    (v0_cnv(i,j,k,3)/jac+v0_cnv(ip,j,k,3)/jacip)          
cccc     .            -abs(v0_cnv(i,j,k,3)/jac+v0_cnv(ip,j,k,3)/jacip) )
cccc     .                                              *db(ip,j,k,1)/jacip)
cccc            else
cc              flxip = 0.5/(jac+jacip)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(i ,j,k,3)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(ip,j,k,3))
cc     .               -0.5/(jac+jacip)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3)) )*db(i ,j,k,1)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3)) )*db(ip,j,k,1))
cccc            endif
cc
cccc            if (isSP(im,j,k,igrid,igrid,igrid)) then
cccccc              flxim = 0.125*(jac+jacim)*(
cccccc     .           (    (v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim)
cccccc     .            +abs(v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim) )
cccccc     .                                              *db(i ,j,k,3)/jac
cccccc     .          +(    (v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim)          
cccccc     .            -abs(v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim) )
cccccc     .                                              *db(im,j,k,3)/jacim)
cccccc     .               -0.125*(jac+jacim)*(
cccccc     .           (    (v0_cnv(i,j,k,3)/jac+v0_cnv(im,j,k,3)/jacim)
cccccc     .            +abs(v0_cnv(i,j,k,3)/jac+v0_cnv(im,j,k,3)/jacim) )
cccccc     .                                              *db(i ,j,k,1)/jac
cccccc     .          +(    (v0_cnv(i,j,k,3)/jac+v0_cnv(im,j,k,3)/jacim)          
cccccc     .            -abs(v0_cnv(i,j,k,3)/jac+v0_cnv(im,j,k,3)/jacim) )
cccccc     .                                              *db(im,j,k,1)/jacim)
cccc              flxim = 0.125*(jac+jacim)*(
cccc     .           (    (v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim)
cccc     .            +abs(v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim) )
cccc     .                                              *db(im,j,k,3)/jac
cccc     .          +(    (v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim)          
cccc     .            -abs(v0_cnv(i,j,k,1)/jac+v0_cnv(im,j,k,1)/jacim) )
cccc     .                                              *db(i ,j,k,3)/jacim)
cccc     .               -0.125*(jac+jacim)*(
cccc     .           (    (v0_cnv(i,j,k,3)/jac+v0_cnv(im,j,k,3)/jacim)
cccc     .            +abs(v0_cnv(i,j,k,3)/jac+v0_cnv(im,j,k,3)/jacim) )
cccc     .                                              *db(im,j,k,1)/jac
cccc     .          +(    (v0_cnv(i,j,k,3)/jac+v0_cnv(im,j,k,3)/jacim)          
cccc     .            -abs(v0_cnv(i,j,k,3)/jac+v0_cnv(im,j,k,3)/jacim) )
cccc     .                                              *db(i ,j,k,1)/jacim)
cccc            elseif (isSP(i,j,k,igrid,igrid,igrid)) then
cccccc              flxim = 0.5/jacim*(
cccccc     .           (    (v0_cnv(im,j,k,1))
cccccc     .            +abs(v0_cnv(im,j,k,1)) )*db(i ,j,k,3)
cccccc     .          +(    (v0_cnv(im,j,k,1))          
cccccc     .            -abs(v0_cnv(im,j,k,1)) )*db(im,j,k,3))
cccccc     .               -0.5/jacim*(
cccccc     .           (    (v0_cnv(im,j,k,3))
cccccc     .            +abs(v0_cnv(im,j,k,3)) )*db(i ,j,k,1)
cccccc     .          +(    (v0_cnv(im,j,k,3))          
cccccc     .            -abs(v0_cnv(im,j,k,3)) )*db(im,j,k,1))
cccc              flxim = 0.5/jacim*(
cccc     .           (    (v0_cnv(im,j,k,1))
cccc     .            +abs(v0_cnv(im,j,k,1)) )*db(im,j,k,3)
cccc     .          +(    (v0_cnv(im,j,k,1))          
cccc     .            -abs(v0_cnv(im,j,k,1)) )*db(i ,j,k,3))
cccc     .               -0.5/jacim*(
cccc     .           (    (v0_cnv(im,j,k,3))
cccc     .            +abs(v0_cnv(im,j,k,3)) )*db(im,j,k,1)
cccc     .          +(    (v0_cnv(im,j,k,3))          
cccc     .            -abs(v0_cnv(im,j,k,3)) )*db(i ,j,k,1))
cccc            else
cccc              flxim = 0.5/(jac+jacim)*(
cccc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
cccc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(i ,j,k,3)
cccc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
cccc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(im,j,k,3))
cccc     .               -0.5/(jac+jacim)*(
cccc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3))
cccc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3)) )*db(i ,j,k,1)
cccc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3))          
cccc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3)) )*db(im,j,k,1))
cc              flxim = 0.5/(jac+jacim)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(im,j,k,3)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(i ,j,k,3))
cc     .               -0.5/(jac+jacim)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3)) )*db(im,j,k,1)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3)) )*db(i ,j,k,1))
cccc            endif
cc
cc            flxjp = 0.5/(jac+jacjp)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,j ,k,3)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,jp,k,3))
cc     .             -0.5/(jac+jacjp)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,jp,k,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,jp,k,3)) )*db(i,j ,k,2)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,jp,k,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,jp,k,3)) )*db(i,jp,k,2))
cccc            flxjm = 0.5/(jac+jacjm)*(
cccc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
cccc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,j ,k,3)
cccc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
cccc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,jm,k,3))
cccc     .             -0.5/(jac+jacjm)*(
cccc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3))
cccc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3)) )*db(i,j ,k,2)
cccc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3))          
cccc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3)) )*db(i,jm,k,2))
cc            flxjm = 0.5/(jac+jacjm)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,jm,k,3)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,j ,k,3))
cc     .             -0.5/(jac+jacjm)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3)) )*db(i,jm,k,2)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3)) )*db(i,j ,k,2))
cc
cc            upwind =  (flxip-flxim)/dxh(ig)
cc     .               +(flxjp-flxjm)/dyh(jg)
cc
cc            cnv(3) = ( db(i,j,k,3)/dt + alpha*upwind
cc     .                - alpha*etal*veclap(3) )*vol
cc
cc            do ieq=1,3
cc              y(neq*(ijk-1)+ieq) = cnv(ieq)
cc            enddo
cc
cc          enddo
cc        enddo
cc      enddo
cc
ccc End program
cc
cc      deallocate(db)
cc      nullify(v0_cnv)
cc
cc      call deallocPointers(fpointers)
cc
cc      end subroutine b_mtvc

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

      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

      real(8)    :: upwind,etal,flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol
     .             ,veclap(3),cnv(3),cov(3),car(3)

      real(8),allocatable,dimension(:,:,:,:) :: db
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      logical    :: fpointers,alt_eom_b,is_cnv

c Begin program

      is_cnv = .true.

      alt_eom_b = .false.

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

      nxx = grid_params%nxv(igrid)
      nyy = grid_params%nyv(igrid)
      nzz = grid_params%nzv(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

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
     .            ,icomp=IBX,is_cnv=is_cnv,is_vec=.true.)
cc     .            ,icomp=IBX,is_cnv=is_cnv,is_vec=.true.,iorder=2)

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

            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + 1d-10

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            if (vol_wgt) then
              vol = volume(i,j,k,igrid,igrid,igrid)
            else
              vol = 1d0
            endif

            veclap = veclaplacian(i,j,k,nxx,nyy,nzz
     .                           ,igrid,igrid,igrid,db
     .                           ,alt_eom_b,vol=.false.)

            etal   = res(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid)

            if (gm_smooth) then

              call find_curl_vxb(i,j,k,nxx,nyy,nzz,v0_cnv,db
     .                          ,cnv(1),cnv(2),cnv(3),0,igrid)

            else


            flxjp = 0.5/(jac+jacjp)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,j ,k,1)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))         
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,jp,k,1))
cc     .             -0.5/(jac+jacjp)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1)) )*db(i,j ,k,2)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,jp,k,1)) )*db(i,jp,k,2))
     .             -0.5*(v0_cnv(i,j ,k,1)*db(i,j ,k,2)/jac
     .                  +v0_cnv(i,jp,k,1)*db(i,jp,k,2)/jacjp)
            flxjm = 0.5/(jac+jacjm)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,jm,k,1)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,j ,k,1))
cc     .             -0.5/(jac+jacjm)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1)) )*db(i,j ,k,2)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,jm,k,1)) )*db(i,jm,k,2))
     .             -0.5*(v0_cnv(i,j ,k,1)*db(i,j ,k,2)/jac
     .                  +v0_cnv(i,jm,k,1)*db(i,jm,k,2)/jacjm)

            flxkp = 0.5/(jac+jackp)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,k ,1)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,kp,1))
cc     .             -0.5/(jac+jackp)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1)) )*db(i,j,k ,3)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,kp,1)) )*db(i,j,kp,3))
     .             -0.5*(v0_cnv(i,j,k ,1)*db(i,j,k ,3)/jac
     .                  +v0_cnv(i,j,kp,1)*db(i,j,kp,3)/jackp)
            flxkm = 0.5/(jac+jackm)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,km,1)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,k ,1))
cc     .             -0.5/(jac+jackm)*(
cc     .           (    (v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1))
cc     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1)) )*db(i,j,k ,3)
cc     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1))          
cc     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(i,j,km,1)) )*db(i,j,km,3))
     .             -0.5*(v0_cnv(i,j,k ,1)*db(i,j,k ,3)/jac
     .                  +v0_cnv(i,j,km,1)*db(i,j,km,3)/jackm)

            cnv(1) =  (flxjp-flxjm)/dyh(jg)
     .               +(flxkp-flxkm)/dzh(kg)

            flxip = 0.5/(jac+jacip)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(i ,j,k,2)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(ip,j,k,2))
cc     .             -0.5/(jac+jacip)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2)) )*db(i ,j,k,1)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(ip,j,k,2)) )*db(ip,j,k,1))
     .             -0.5*(v0_cnv(i ,j,k,2)*db(i ,j,k,1)/jac
     .                  +v0_cnv(ip,j,k,2)*db(ip,j,k,1)/jacip)
            flxim = 0.5/(jac+jacim)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(im,j,k,2)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(i ,j,k,2))
cc     .             -0.5/(jac+jacim)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )*db(i ,j,k,1)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(im,j,k,2)) )*db(im,j,k,1))
     .             -0.5*(v0_cnv(i ,j,k,2)*db(i ,j,k,1)/jac
     .                  +v0_cnv(im,j,k,2)*db(im,j,k,1)/jacim)

            flxkp = 0.5/(jac+jackp)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,k ,2)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,kp,2))
cc     .             -0.5/(jac+jackp)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2)) )*db(i,j,k ,3)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,kp,2)) )*db(i,j,kp,3))
     .             -0.5*(v0_cnv(i,j,k ,2)*db(i,j,k ,3)/jac
     .                  +v0_cnv(i,j,kp,2)*db(i,j,kp,3)/jackp)
            flxkm = 0.5/(jac+jackm)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,km,2)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,k ,2))
cc     .             -0.5/(jac+jackm)*(
cc     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2))
cc     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2)) )*db(i,j,k ,3)
cc     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2))          
cc     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,j,km,2)) )*db(i,j,km,3))
     .             -0.5*(v0_cnv(i,j,k ,2)*db(i,j,k ,3)/jac
     .                  +v0_cnv(i,j,km,2)*db(i,j,km,3)/jackm)

            cnv(2) =  (flxip-flxim)/dxh(ig)
     .               +(flxkp-flxkm)/dzh(kg)

            flxip = 0.5/(jac+jacip)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(i ,j,k,3)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(ip,j,k,3))
cc     .             -0.5/(jac+jacip)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3)) )*db(i ,j,k,1)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(ip,j,k,3)) )*db(ip,j,k,1))
     .             -0.5*(v0_cnv(i ,j,k,3)*db(i ,j,k,1)/jac
     .                  +v0_cnv(ip,j,k,3)*db(ip,j,k,1)/jacip)
            flxim = 0.5/(jac+jacim)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(im,j,k,3)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(i ,j,k,3))
cc     .             -0.5/(jac+jacim)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3)) )*db(i ,j,k,1)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(im,j,k,3)) )*db(im,j,k,1))
     .             -0.5*(v0_cnv(i ,j,k,3)*db(i ,j,k,1)/jac
     .                  +v0_cnv(im,j,k,3)*db(im,j,k,1)/jacim)

            flxjp = 0.5/(jac+jacjp)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,j ,k,3)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,jp,k,3))
cc     .             -0.5/(jac+jacjm)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3)) )*db(i,j ,k,2)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3)) )*db(i,jm,k,2))
     .             -0.5*(v0_cnv(i,j ,k,3)*db(i,j ,k,2)/jac
     .                  +v0_cnv(i,jp,k,3)*db(i,jp,k,2)/jacjp)
            flxjm = 0.5/(jac+jacjm)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,jm,k,3)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,j ,k,3))
cc     .             -0.5/(jac+jacjm)*(
cc     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3))
cc     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3)) )*db(i,j ,k,2)
cc     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3))          
cc     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,jm,k,3)) )*db(i,jm,k,2))
     .             -0.5*(v0_cnv(i,j ,k,3)*db(i,j ,k,2)/jac
     .                  +v0_cnv(i,jm,k,3)*db(i,jm,k,2)/jacjm)

            cnv(3) =  (flxip-flxim)/dxh(ig)
     .               +(flxjp-flxjm)/dyh(jg)

            endif

            cnv = (db(i,j,k,:)/dt + alpha*cnv - alpha*etal*veclap)*vol

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = cnv(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(db)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine b_mtvc

c module v_mtvc_mod
c ##################################################################
      module v_mtvc_mod

      use precond_variables

      real(8),allocatable,dimension(:,:,:,:) :: dv

      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv,b0_cnv,b0_cov
     .                                         ,rho0,pp0

      integer(4) :: isig,iisig,ip,im,jp,jm,kp,km,hex,hey,hez
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
     .             ,nxx,nyy,nzz,ijk,ijkg,isp,igr

      real(8),allocatable,dimension(:,:)     :: a_sp,a_vol

      real(8)    :: vol,psiv(3),psib(3),psit(3),cov(3),cnv(3),car(3)
     .             ,mag(3),veclap(3)

      logical    :: fpointers,is_cnv,is_car,postproc

      contains

c     findPsiv
c     #####################################################################
      subroutine findPsiv

      implicit none

      real(8)    :: upwind,mul,nabla_vv0(3,3)

      hex = 1
      hey = 1
      hez = 1
      if (v0_cnv(i,j,k,1) > 0d0) hex = -1
      if (v0_cnv(i,j,k,2) > 0d0) hey = -1
      if (v0_cnv(i,j,k,3) > 0d0) hez = -1
      hex = 0
      hey = 0
      hez = 0

      nabla_v = fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igr,igr,igr
     .                       ,dv(:,:,:,1),dv(:,:,:,2),dv(:,:,:,3)
     .                       ,hex,hey,hez)

cc      nabla_vv0 =fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igr,igr,igr
cc     .                        ,v0_cnv(:,:,:,1)
cc     .                        ,v0_cnv(:,:,:,2)
cc     .                        ,v0_cnv(:,:,:,3)
cc     .                        ,0,0,0)
      nabla_vv0 = mgnablaV0(ijkg,:,:)

      mul = vis(i,j,k,nxx,nyy,nzz,igr,igr,igr)

      veclap = veclaplacian(i,j,k,nxx,nyy,nzz,igr,igr,igr,dv
     .                     ,alt_eom,vol=.false.)

      do ieq=1,3
        upwind =( v0_cnv(i,j,k,1)*nabla_v(1,ieq)
     .           +v0_cnv(i,j,k,2)*nabla_v(2,ieq)
     .           +v0_cnv(i,j,k,3)*nabla_v(3,ieq))/jac

        upwind = upwind
     .           +( dv(i,j,k,1)*nabla_vv0(1,ieq)
     .             +dv(i,j,k,2)*nabla_vv0(2,ieq)
     .             +dv(i,j,k,3)*nabla_vv0(3,ieq))/jac

        hex = floor(sign(1d0,-mgadvdiffV0(ijkg,ieq)))
cc        hex = 0
        upwind = upwind
     .            -dt*mgadvdiffV0(ijkg,ieq)*div_upwd(hex)/rho0(i,j,k,1)

        cnv(ieq) = dv(i,j,k,ieq)/dt
     .           + alpha*upwind
     .           - alpha*mul*veclap(ieq)
      enddo

      if (si_car) then
        call transformVectorToCartesian
     .              (i,j,k,igr,igr,igr
     .               ,cnv(1),cnv(2),cnv(3)
     .               ,.false.
     .               ,car(1),car(2),car(3))
        psiv = rho0(i,j,k,1)*car*vol
      elseif (is_cnv) then
        psiv = rho0(i,j,k,1)*cnv*vol
      else
        call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
        psiv = rho0(i,j,k,1)*cov*vol
      endif

      end subroutine findPsiv

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
        
        dxx = 2*dxh(ig)
        dyy = 2*dyh(jg)
        dzz = 2*dzh(kg)

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
cc        else
cc          div_upwd = 0d0
cc          return
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

c     findPsit
c     #####################################################################
      subroutine findPsit

      implicit none

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
     $             ,divip,divim,divjp,divjm,divkp,divkm

      !Fluxes at faces for calculation of grad[div(dv p0)]
      flxip =( (dv(ip,j ,k ,1)*pp0(ip,j ,k ,1)
     .         -dv(i ,j ,k ,1)*pp0(i ,j ,k ,1))/dx(ig)
     .        +(dv(ip,jp,k ,2)*pp0(ip,jp,k ,1)
     .         -dv(ip,jm,k ,2)*pp0(ip,jm,k ,1)
     .         +dv(i ,jp,k ,2)*pp0(i ,jp,k ,1)
     .         -dv(i ,jm,k ,2)*pp0(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv(ip,j ,kp,3)*pp0(ip,j ,kp,1)
     .         -dv(ip,j ,km,3)*pp0(ip,j ,km,1)
     .         +dv(i ,j ,kp,3)*pp0(i ,j ,kp,1)
     .         -dv(i ,j ,km,3)*pp0(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacip)*2
      flxim =( (dv(i ,j ,k ,1)*pp0(i ,j ,k ,1)
     .         -dv(im,j ,k ,1)*pp0(im,j ,k ,1))/dx(ig-1)
     .        +(dv(im,jp,k ,2)*pp0(im,jp,k ,1)
     .         -dv(im,jm,k ,2)*pp0(im,jm,k ,1)
     .         +dv(i ,jp,k ,2)*pp0(i ,jp,k ,1)
     .         -dv(i ,jm,k ,2)*pp0(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv(im,j ,kp,3)*pp0(im,j ,kp,1)
     .         -dv(im,j ,km,3)*pp0(im,j ,km,1)
     .         +dv(i ,j ,kp,3)*pp0(i ,j ,kp,1)
     .         -dv(i ,j ,km,3)*pp0(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacim)*2

      flxjp =( (dv(ip,jp,k ,1)*pp0(ip,jp,k ,1)
     .         -dv(im,jp,k ,1)*pp0(im,jp,k ,1)
     .         +dv(ip,j ,k ,1)*pp0(ip,j ,k ,1)
     .         -dv(im,j ,k ,1)*pp0(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv(i ,jp,k ,2)*pp0(i ,jp,k ,1)
     .         -dv(i ,j ,k ,2)*pp0(i ,j ,k ,1))/dy(jg)
     .        +(dv(i ,jp,kp,3)*pp0(i ,jp,kp,1)
     .         -dv(i ,jp,km,3)*pp0(i ,jp,km,1)
     .         +dv(i ,j ,kp,3)*pp0(i ,j ,kp,1)
     .         -dv(i ,j ,km,3)*pp0(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacjp)*2
      flxjm =( (dv(ip,jm,k ,1)*pp0(ip,jm,k ,1)
     .         -dv(im,jm,k ,1)*pp0(im,jm,k ,1)
     .         +dv(ip,j ,k ,1)*pp0(ip,j ,k ,1)
     .         -dv(im,j ,k ,1)*pp0(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv(i ,j ,k ,2)*pp0(i ,j ,k ,1)
     .         -dv(i ,jm,k ,2)*pp0(i ,jm,k ,1))/dy(jg-1)
     .        +(dv(i ,jm,kp,3)*pp0(i ,jm,kp,1)
     .         -dv(i ,jm,km,3)*pp0(i ,jm,km,1)
     .         +dv(i ,j ,kp,3)*pp0(i ,j ,kp,1)
     .         -dv(i ,j ,km,3)*pp0(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacjm)*2

      flxkp =( (dv(ip,j ,kp,1)*pp0(ip,j ,kp,1)
     .         -dv(im,j ,kp,1)*pp0(im,j ,kp,1)
     .         +dv(ip,j ,k ,1)*pp0(ip,j ,k ,1)
     .         -dv(im,j ,k ,1)*pp0(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv(i ,jp,kp,2)*pp0(i ,jp,kp,1)
     .         -dv(i ,jm,kp,2)*pp0(i ,jm,kp,1)
     .         +dv(i ,jp,k ,2)*pp0(i ,jp,k ,1)
     .         -dv(i ,jm,k ,2)*pp0(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv(i ,j ,kp,3)*pp0(i ,j ,kp,1)
     .         -dv(i ,j ,k ,3)*pp0(i ,j ,k ,1))/dz(kg) )
     $        /(jac+jackp)*2
      flxkm =( (dv(ip,j ,km,1)*pp0(ip,j ,km,1)
     .         -dv(im,j ,km,1)*pp0(im,j ,km,1)
     .         +dv(ip,j ,k ,1)*pp0(ip,j ,k ,1)
     .         -dv(im,j ,k ,1)*pp0(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv(i ,jp,km,2)*pp0(i ,jp,km,1)
     .         -dv(i ,jm,km,2)*pp0(i ,jm,km,1)
     .         +dv(i ,jp,k ,2)*pp0(i ,jp,k ,1)
     .         -dv(i ,jm,k ,2)*pp0(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv(i ,j ,k ,3)*pp0(i ,j ,k ,1)
     .         -dv(i ,j ,km,3)*pp0(i ,j ,km,1))/dz(kg-1) )
     $        /(jac+jackm)*2

      !Fluxes at faces for calculation of grad[(gamma-1)*p0*div(dv)]

      !!Divergence at faces i+-1/2, etc.
      divip = (dv(ip,j ,k,1)-dv(i ,j ,k,1))/dx(ig)
     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .        +dv(ip,jp,k,2)-dv(ip,jm,k,2))/dyh(jg)/4.
     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .        +dv(ip,j,kp,3)-dv(ip,j,km,3))/dzh(kg)/4.
      divim = (dv(i ,j ,k,1)-dv(im,j ,k,1))/dx(ig-1)
     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .        +dv(im,jp,k,2)-dv(im,jm,k,2))/dyh(jg)/4.
     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .        +dv(im,j,kp,3)-dv(im,j,km,3))/dzh(kg)/4.

      divjp = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .        +dv(ip,jp,k,1)-dv(im,jp,k,1))/dxh(ig)/4.
     .       +(dv(i ,jp,k,2)-dv(i ,j ,k,2))/dy(jg)
     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .        +dv(i,jp,kp,3)-dv(i,jp,km,3))/dzh(kg)/4.
      divjm = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .        +dv(ip,jm,k,1)-dv(im,jm,k,1))/dxh(ig)/4.
     .       +(dv(i ,j ,k,2)-dv(i ,jm,k,2))/dy(jg-1)
     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .        +dv(i,jm,kp,3)-dv(i,jm,km,3))/dzh(kg)/4.

      divkp = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .        +dv(ip,j,kp,1)-dv(im,j,kp,1))/dxh(ig)/4.
     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .        +dv(i,jp,kp,2)-dv(i,jm,kp,2))/dyh(jg)/4.
     .       +(dv(i,j ,kp,3)-dv(i,j ,k ,3))/dz(kg)
      divkm = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .        +dv(ip,j,km,1)-dv(im,j,km,1))/dxh(ig)/4.
     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .        +dv(i,jp,km,2)-dv(i,jm,km,2))/dyh(jg)/4.
     .       +(dv(i,j ,k ,3)-dv(i,j ,km,3))/dz(kg-1)

      flxip = flxip
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(ip,j,k,1))*divip/(jac+jacip)
      flxim = flxim
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(im,j,k,1))*divim/(jac+jacim)

      flxjp = flxjp
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(i,jp,k,1))*divjp/(jac+jacjp)
      flxjm = flxjm
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(i,jm,k,1))*divjm/(jac+jacjm)

      flxkp = flxkp
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(i,j,kp,1))*divkp/(jac+jackp)
      flxkm = flxkm
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(i,j,km,1))*divkm/(jac+jackm)

      cov(1) = (flxip - flxim)/dxh(ig)
      cov(2) = (flxjp - flxjm)/dyh(jg)
      cov(3) = (flxkp - flxkm)/dzh(kg)

      if (si_car) then
        call transformVectorToCartesian(i,j,k,igr,igr,igr
     .                                 ,cov(1),cov(2),cov(3),.true.
     .                                 ,car(1),car(2),car(3))
        psit = -alpha**2*car*vol
      elseif (is_cnv) then
        call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
        psit = -alpha**2*cnv*vol
      else
        psit = -alpha**2*cov*vol
      endif

      end subroutine findPsit

c     findPsit_diag
c     #####################################################################
      subroutine findPsit_diag

      implicit none

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
     $             ,divip,divim,divjp,divjm,divkp,divkm

cc      integer(4) :: ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

      !Fluxes at faces for calculation of grad[div(dv p0)]
      flxip =( (dv(ip,j ,k ,1)*pp0(ip,j ,k ,1)
     .         -dv(i ,j ,k ,1)*pp0(i ,j ,k ,1))/dx(ig)
     .        +(dv(ip,jp,k ,2)*pp0(ip,jp,k ,1)
     .         -dv(ip,jm,k ,2)*pp0(ip,jm,k ,1)
     .         +dv(i ,jp,k ,2)*pp0(i ,jp,k ,1)
     .         -dv(i ,jm,k ,2)*pp0(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv(ip,j ,kp,3)*pp0(ip,j ,kp,1)
     .         -dv(ip,j ,km,3)*pp0(ip,j ,km,1)
     .         +dv(i ,j ,kp,3)*pp0(i ,j ,kp,1)
     .         -dv(i ,j ,km,3)*pp0(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacip)*2
      flxim =( (dv(i ,j ,k ,1)*pp0(i ,j ,k ,1)
     .         -dv(im,j ,k ,1)*pp0(im,j ,k ,1))/dx(ig-1)
     .        +(dv(im,jp,k ,2)*pp0(im,jp,k ,1)
     .         -dv(im,jm,k ,2)*pp0(im,jm,k ,1)
     .         +dv(i ,jp,k ,2)*pp0(i ,jp,k ,1)
     .         -dv(i ,jm,k ,2)*pp0(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv(im,j ,kp,3)*pp0(im,j ,kp,1)
     .         -dv(im,j ,km,3)*pp0(im,j ,km,1)
     .         +dv(i ,j ,kp,3)*pp0(i ,j ,kp,1)
     .         -dv(i ,j ,km,3)*pp0(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacim)*2

      flxjp =( (dv(ip,jp,k ,1)*pp0(ip,jp,k ,1)
     .         -dv(im,jp,k ,1)*pp0(im,jp,k ,1)
     .         +dv(ip,j ,k ,1)*pp0(ip,j ,k ,1)
     .         -dv(im,j ,k ,1)*pp0(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv(i ,jp,k ,2)*pp0(i ,jp,k ,1)
     .         -dv(i ,j ,k ,2)*pp0(i ,j ,k ,1))/dy(jg)
     .        +(dv(i ,jp,kp,3)*pp0(i ,jp,kp,1)
     .         -dv(i ,jp,km,3)*pp0(i ,jp,km,1)
     .         +dv(i ,j ,kp,3)*pp0(i ,j ,kp,1)
     .         -dv(i ,j ,km,3)*pp0(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacjp)*2
      flxjm =( (dv(ip,jm,k ,1)*pp0(ip,jm,k ,1)
     .         -dv(im,jm,k ,1)*pp0(im,jm,k ,1)
     .         +dv(ip,j ,k ,1)*pp0(ip,j ,k ,1)
     .         -dv(im,j ,k ,1)*pp0(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv(i ,j ,k ,2)*pp0(i ,j ,k ,1)
     .         -dv(i ,jm,k ,2)*pp0(i ,jm,k ,1))/dy(jg-1)
     .        +(dv(i ,jm,kp,3)*pp0(i ,jm,kp,1)
     .         -dv(i ,jm,km,3)*pp0(i ,jm,km,1)
     .         +dv(i ,j ,kp,3)*pp0(i ,j ,kp,1)
     .         -dv(i ,j ,km,3)*pp0(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacjm)*2

      flxkp =( (dv(ip,j ,kp,1)*pp0(ip,j ,kp,1)
     .         -dv(im,j ,kp,1)*pp0(im,j ,kp,1)
     .         +dv(ip,j ,k ,1)*pp0(ip,j ,k ,1)
     .         -dv(im,j ,k ,1)*pp0(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv(i ,jp,kp,2)*pp0(i ,jp,kp,1)
     .         -dv(i ,jm,kp,2)*pp0(i ,jm,kp,1)
     .         +dv(i ,jp,k ,2)*pp0(i ,jp,k ,1)
     .         -dv(i ,jm,k ,2)*pp0(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv(i ,j ,kp,3)*pp0(i ,j ,kp,1)
     .         -dv(i ,j ,k ,3)*pp0(i ,j ,k ,1))/dz(kg) )
     $        /(jac+jackp)*2
      flxkm =( (dv(ip,j ,km,1)*pp0(ip,j ,km,1)
     .         -dv(im,j ,km,1)*pp0(im,j ,km,1)
     .         +dv(ip,j ,k ,1)*pp0(ip,j ,k ,1)
     .         -dv(im,j ,k ,1)*pp0(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv(i ,jp,km,2)*pp0(i ,jp,km,1)
     .         -dv(i ,jm,km,2)*pp0(i ,jm,km,1)
     .         +dv(i ,jp,k ,2)*pp0(i ,jp,k ,1)
     .         -dv(i ,jm,k ,2)*pp0(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv(i ,j ,k ,3)*pp0(i ,j ,k ,1)
     .         -dv(i ,j ,km,3)*pp0(i ,j ,km,1))/dz(kg-1) )
     $        /(jac+jackm)*2

      !Fluxes at faces for calculation of grad[(gamma-1)*p0*div(dv)]

      !!Divergence at faces i+-1/2, etc.
      divip = (dv(ip,j ,k,1)-dv(i ,j ,k,1))/dx(ig)
     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .        +dv(ip,jp,k,2)-dv(ip,jm,k,2))/dyh(jg)/4.
     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .        +dv(ip,j,kp,3)-dv(ip,j,km,3))/dzh(kg)/4.
      divim = (dv(i ,j ,k,1)-dv(im,j ,k,1))/dx(ig-1)
     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .        +dv(im,jp,k,2)-dv(im,jm,k,2))/dyh(jg)/4.
     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .        +dv(im,j,kp,3)-dv(im,j,km,3))/dzh(kg)/4.

      divjp = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .        +dv(ip,jp,k,1)-dv(im,jp,k,1))/dxh(ig)/4.
     .       +(dv(i ,jp,k,2)-dv(i ,j ,k,2))/dy(jg)
     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .        +dv(i,jp,kp,3)-dv(i,jp,km,3))/dzh(kg)/4.
      divjm = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .        +dv(ip,jm,k,1)-dv(im,jm,k,1))/dxh(ig)/4.
     .       +(dv(i ,j ,k,2)-dv(i ,jm,k,2))/dy(jg-1)
     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .        +dv(i,jm,kp,3)-dv(i,jm,km,3))/dzh(kg)/4.

      divkp = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .        +dv(ip,j,kp,1)-dv(im,j,kp,1))/dxh(ig)/4.
     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .        +dv(i,jp,kp,2)-dv(i,jm,kp,2))/dyh(jg)/4.
     .       +(dv(i,j ,kp,3)-dv(i,j ,k ,3))/dz(kg)
      divkm = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .        +dv(ip,j,km,1)-dv(im,j,km,1))/dxh(ig)/4.
     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .        +dv(i,jp,km,2)-dv(i,jm,km,2))/dyh(jg)/4.
     .       +(dv(i,j ,k ,3)-dv(i,j ,km,3))/dz(kg-1)

      flxip = flxip
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(ip,j,k,1))*divip/(jac+jacip)
      flxim = flxim
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(im,j,k,1))*divim/(jac+jacim)

      flxjp = flxjp
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(i,jp,k,1))*divjp/(jac+jacjp)
      flxjm = flxjm
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(i,jm,k,1))*divjm/(jac+jacjm)

      flxkp = flxkp
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(i,j,kp,1))*divkp/(jac+jackp)
      flxkm = flxkm
     .      + (gamma-1.)*(pp0(i,j,k,1)+pp0(i,j,km,1))*divkm/(jac+jackm)


      cov(1) = (flxip - flxim)/dxh(ig)
      cov(2) = (flxjp - flxjm)/dyh(jg)
      cov(3) = (flxkp - flxkm)/dzh(kg)

      cov = cov*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol

      if (si_car) then
        call transformVectorToCartesian(i,j,k,igr,igr,igr
     .                                 ,cov(1),cov(2),cov(3),.true.
     .                                 ,car(1),car(2),car(3))
        psit = -alpha**2*car*vol
      elseif (is_cnv) then
        call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
        psit = -alpha**2*cnv*vol
      else
        psit = -alpha**2*cov*vol
      endif

      end subroutine findPsit_diag

ccc     findPsit
ccc     #####################################################################
cc      subroutine findPsit
cc
cc      implicit none
cc
cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc     $             ,divip,divim,divjp,divjm,divkp,divkm
cc
cc      !Fluxes at faces for calculation of grad(dv.grad(p0))
cc      flxip =( (dv(i,j,k,1)+dv(ip,j,k,1))
cc     .           *(pp0(ip,j,k,1)-pp0(i,j,k,1))/dx(ig)
cc     .        +(dv(i,j,k,2)+dv(ip,j,k,2))
cc     .           *(pp0(ip,jp,k,1)-pp0(ip,jm,k,1)
cc     .            +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
cc     .        +(dv(i,j,k,3)+dv(ip,j,k,3))
cc     .           *(pp0(ip,j,kp,1)-pp0(ip,j,km,1)
cc     .            +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4. )
cc     $        /(jac+jacip)
cc      flxim =( (dv(i,j,k,1)+dv(im,j,k,1))
cc     .            *(pp0(i ,j,k,1)-pp0(im,j,k,1))/dx(ig-1)
cc     .        +(dv(i,j,k,2)+dv(im,j,k,2))
cc     .            *(pp0(im,jp,k,1)-pp0(im,jm,k,1)
cc     .             +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
cc     .        +(dv(i,j,k,3)+dv(im,j,k,3))
cc     .            *(pp0(im,j,kp,1)-pp0(im,j,km,1)
cc     .             +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4.)
cc     $        /(jac+jacim)
cc
cc      flxjp =( (dv(i,j,k,1)+dv(i,jp,k,1))
cc     .            *(pp0(ip,jp,k,1)-pp0(im,jp,k,1)
cc     .             +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/dxh(ig)/4.
cc     .        +(dv(i,j,k,2)+dv(i,jp,k,2))
cc     .            *(pp0(i,jp,k,1)-pp0(i,j,k,1))/dy(jg)
cc     .        +(dv(i,j,k,3)+dv(i,jp,k,3))
cc     .            *(pp0(i,jp,kp,1)-pp0(i,jp,km,1)
cc     .             +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/dzh(kg)/4.)
cc     $        /(jac+jacjp)
cc      flxjm =( (dv(i,j,k,1)+dv(i,jm,k,1))
cc     .            *(pp0(ip,jm,k,1)-pp0(im,jm,k,1)
cc     .             +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/dxh(ig)/4.
cc     .        +(dv(i,j,k,2)+dv(i,jm,k,2))
cc     .            *(pp0(i,j ,k,1)-pp0(i,jm,k,1))/dy(jg-1)
cc     .        +(dv(i,j,k,3)+dv(i,jm,k,3))
cc     .            *(pp0(i,jm,kp,1)-pp0(i,jm,km,1)
cc     .             +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/dzh(kg)/4.)
cc     $        /(jac+jacjm)
cc
cc      flxkp =( (dv(i,j,k,1)+dv(i,j,kp,1))
cc     .            *(pp0(ip,j,kp,1)-pp0(im,j,kp,1)
cc     .             +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/dxh(ig)/4.
cc     .        +(dv(i,j,k,2)+dv(i,j,kp,2))
cc     .            *(pp0(i,jp,kp,1)-pp0(i,jm,kp,1)
cc     .             +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/dyh(jg)/4.
cc     .        +(dv(i,j,k,3)+dv(i,j,kp,3))
cc     .            *(pp0(i,j,kp,1)-pp0(i,j,k,1))/dz(kg) )
cc     $        /(jac+jackp)
cc      flxkm =( (dv(i,j,k,1)+dv(i,j,km,1))
cc     .            *(pp0(ip,j,km,1)-pp0(im,j,km,1)
cc     .             +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/dxh(ig)/4.
cc     .        +(dv(i,j,k,2)+dv(i,j,km,2))
cc     .            *(pp0(i,jp,km,1)-pp0(i,jm,km,1)
cc     .             +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/dyh(jg)/4.
cc     .        +(dv(i,j,k,3)+dv(i,j,km,3))
cc     .            *(pp0(i,j,k ,1)-pp0(i,j,km,1))/dz(kg-1) )
cc     $        /(jac+jackm)
cc
cc      !Fluxes at faces for calculation of grad(gamma*p0*div(dv))
cc
cc      !!Divergence at faces i+-1/2, etc.
cc      divip = (dv(ip,j ,k,1)-dv(i ,j ,k,1))/dx(ig)
cc     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
cc     .        +dv(ip,jp,k,2)-dv(ip,jm,k,2))/dyh(jg)/4.
cc     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
cc     .        +dv(ip,j,kp,3)-dv(ip,j,km,3))/dzh(kg)/4.
cc      divim = (dv(i ,j ,k,1)-dv(im,j ,k,1))/dx(ig-1)
cc     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
cc     .        +dv(im,jp,k,2)-dv(im,jm,k,2))/dyh(jg)/4.
cc     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
cc     .        +dv(im,j,kp,3)-dv(im,j,km,3))/dzh(kg)/4.
cc
cc      divjp = (dv(ip,j ,k,1)-dv(im,j ,k,1)
cc     .        +dv(ip,jp,k,1)-dv(im,jp,k,1))/dxh(ig)/4.
cc     .       +(dv(i ,jp,k,2)-dv(i ,j ,k,2))/dy(jg)
cc     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
cc     .        +dv(i,jp,kp,3)-dv(i,jp,km,3))/dzh(kg)/4.
cc      divjm = (dv(ip,j ,k,1)-dv(im,j ,k,1)
cc     .        +dv(ip,jm,k,1)-dv(im,jm,k,1))/dxh(ig)/4.
cc     .       +(dv(i ,j ,k,2)-dv(i ,jm,k,2))/dy(jg-1)
cc     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
cc     .        +dv(i,jm,kp,3)-dv(i,jm,km,3))/dzh(kg)/4.
cc
cc      divkp = (dv(ip,j,k ,1)-dv(im,j,k ,1)
cc     .        +dv(ip,j,kp,1)-dv(im,j,kp,1))/dxh(ig)/4.
cc     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
cc     .        +dv(i,jp,kp,2)-dv(i,jm,kp,2))/dyh(jg)/4.
cc     .       +(dv(i,j ,kp,3)-dv(i,j ,k ,3))/dz(kg)
cc      divkm = (dv(ip,j,k ,1)-dv(im,j,k ,1)
cc     .        +dv(ip,j,km,1)-dv(im,j,km,1))/dxh(ig)/4.
cc     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
cc     .        +dv(i,jp,km,2)-dv(i,jm,km,2))/dyh(jg)/4.
cc     .       +(dv(i,j ,k ,3)-dv(i,j ,km,3))/dz(kg-1)
cc
cc      flxip = flxip
cc     .      + gamma*(pp0(i,j,k,1)+pp0(ip,j,k,1))*divip/(jac+jacip)
cc      flxim = flxim
cc     .      + gamma*(pp0(i,j,k,1)+pp0(im,j,k,1))*divim/(jac+jacim)
cc
cc      flxjp = flxjp
cc     .      + gamma*(pp0(i,j,k,1)+pp0(i,jp,k,1))*divjp/(jac+jacjp)
cc      flxjm = flxjm
cc     .      + gamma*(pp0(i,j,k,1)+pp0(i,jm,k,1))*divjm/(jac+jacjm)
cc
cc      flxkp = flxkp
cc     .      + gamma*(pp0(i,j,k,1)+pp0(i,j,kp,1))*divkp/(jac+jackp)
cc      flxkm = flxkm
cc     .      + gamma*(pp0(i,j,k,1)+pp0(i,j,km,1))*divkm/(jac+jackm)
cc
cc      if (isSP(i,j,k,igx,igy,igz)) flxim = 0d0
cc
cc      cov(1) = (flxip - flxim)/dxh(ig)
cc      cov(2) = (flxjp - flxjm)/dyh(jg)
cc      cov(3) = (flxkp - flxkm)/dzh(kg)
cc
cc      if (si_car) then
cc        call transformVectorToCartesian(i,j,k,igr,igr,igr
cc     .                                 ,cov(1),cov(2),cov(3),.true.
cc     .                                 ,car(1),car(2),car(3))
cc        psit = -alpha**2*car*vol
cc      elseif (is_cnv) then
cc        call transformFromCurvToCurv(i,j,k,igr,igr,igr
cc     .                              ,cov(1),cov(2),cov(3)
cc     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
cc        psit = -alpha**2*cnv*vol
cc      else
cc        psit = -alpha**2*cov*vol
cc      endif
cc
cc      end subroutine findPsit

c     findPsib
c     #####################################################################
      subroutine findPsib

      implicit none

      cnv = div_tensor(i,j,k,nxx,nyy,nzz,igr,igr,igr,alt_eom
     .                ,tensor_x,tensor_y,tensor_z,vol=.false.)
cc      cnv = div_tensor(i,j,k,nxx,nyy,nzz,igr,igr,igr,.false.
cc     .                ,tensor_x,tensor_y,tensor_z,vol=.false.)

      if (si_car) then
        call transformVectorToCartesian(i,j,k,igr,igr,igr
     .                                 ,cnv(1),cnv(2),cnv(3),.false.
     .                                 ,car(1),car(2),car(3))
        psib = alpha**2*car*vol
      elseif (is_cnv) then
        psib = alpha**2*cnv*vol
      else
        call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
        psib = alpha**2*cov*vol
      endif

      end subroutine findPsib

c     tensor_x
c     #############################################################
      subroutine tensor_x(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer(4) :: ip,igrid
        real(8)    :: jac,jac0,jacp,gsuper(3,3)
     .               ,acnv(3),acnvp(3),b0cov(3),b0cnv(3),scalar_prod

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

        if (isSP(i+1,j,k,igx,igy,igz)) jac = 1d-10

        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = gmetric%grid(igrid)%jac(ip,j,k)
          jac0 = gmetric%grid(igrid)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        if (flag /= 0) then
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),1,igrid)
        else
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),0,igrid)
        endif

cc        b0cnv = 0.5*(b0_cnv(ip,j,k,:)+b0_cnv(i,j,k,:))
cc        b0cov = 0.5*(b0_cov(ip,j,k,:)+b0_cov(i,j,k,:))

cc        if (spoint .and. flag /= 0) then


        b0cnv(1) = 0.5*(b0_cnv(ip,j,k,1)/jacp+b0_cnv(i,j,k,1)/jac0)*jac
        b0cnv(2) = 0.5*(b0_cnv(ip,j,k,2)     +b0_cnv(i,j,k,2))
        b0cnv(3) = 0.5*(b0_cnv(ip,j,k,3)/jacp+b0_cnv(i,j,k,3)/jac0)*jac

cc        b0cov(1) = 0.5*(b0_cov(ip,j,k,1)     +b0_cov(i,j,k,1)     )
cc        b0cov(2) = 0.5*(b0_cov(ip,j,k,2)/jacp+b0_cov(i,j,k,2)/jac0)*jac
cc        b0cov(3) = 0.5*(b0_cov(ip,j,k,3)     +b0_cov(i,j,k,3)     )

        if (flag /= 0) then
          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=1)
        else
          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=0)
        endif

        scalar_prod = dot_product(acnv,b0cov)

        t11 =( acnv(1)*b0cnv(1)
     .        +acnv(1)*b0cnv(1)
     .        -gsuper(1,1)*scalar_prod )

        t12 =( acnv(1)*b0cnv(2)
     .        +acnv(2)*b0cnv(1)
     .        -gsuper(1,2)*scalar_prod )

        t13 =( acnv(1)*b0cnv(3)
     .        +acnv(3)*b0cnv(1)
     .        -gsuper(1,3)*scalar_prod )

        if (flag /= 0) then
          t11 = t11/jac
          if (.not.alt_eom) t12 = t12/jac
          t13 = t13/jac
        endif

c     End program

      end subroutine tensor_x

c     tensor_y
c     #############################################################
      subroutine tensor_y(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EM SI operator
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer(4) :: jp,igrid
        real(8)    :: jac,gsuper(3,3),acnv(3),b0cov(3),b0cnv(3)
     .                ,scalar_prod

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
     .               +gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),2,igrid)
        else
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),0,igrid)
        endif

        b0cnv = 0.5*(b0_cnv(i,jp,k,:)+b0_cnv(i,j,k,:))
        b0cov = 0.5*(b0_cov(i,jp,k,:)+b0_cov(i,j,k,:))

cc        if (flag /= 0) then
cc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,b0cov(1),b0cov(2),b0cov(3)
cc     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
cc     .                        ,.false.,half_elem=2)
cc        else
cc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,b0cov(1),b0cov(2),b0cov(3)
cc     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
cc     .                        ,.false.,half_elem=0)
cc        endif

        scalar_prod = dot_product(acnv,b0cov)

        t21 =( acnv(2)*b0cnv(1)
     .        +acnv(1)*b0cnv(2)
     .        -gsuper(2,1)*scalar_prod )

        t22 =( acnv(2)*b0cnv(2)
     .        +acnv(2)*b0cnv(2)
     .        -gsuper(2,2)*scalar_prod )

        t23 =( acnv(2)*b0cnv(3)
     .        +acnv(3)*b0cnv(2)
     .        -gsuper(2,3)*scalar_prod )

        if (flag /= 0) then
          t21 = t21/jac
          if (.not.alt_eom) t22 = t22/jac
          t23 = t23/jac
        endif


c     End program

      end subroutine tensor_y

c     tensor_z
c     #############################################################
      subroutine tensor_z(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EM SI operator
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer(4) :: kp,igrid
        real(8)    :: jac,gsuper(3,3),acnv(3),b0cov(3),b0cnv(3)
     .               ,scalar_prod

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
     .               +gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),3,igrid)
        else
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),0,igrid)
        endif

        b0cnv = 0.5*(b0_cnv(i,j,kp,:)+b0_cnv(i,j,k,:))
        b0cov = 0.5*(b0_cov(i,j,kp,:)+b0_cov(i,j,k,:))

cc        if (flag /= 0) then
cc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,b0cov(1),b0cov(2),b0cov(3)
cc     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
cc     .                        ,.false.,half_elem=3)
cc        else
cc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,b0cov(1),b0cov(2),b0cov(3)
cc     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
cc     .                        ,.false.,half_elem=0)
cc        endif

        scalar_prod = dot_product(acnv,b0cov)

        t31 =( acnv(3)*b0cnv(1)
     .        +acnv(1)*b0cnv(3)
     .        -gsuper(3,1)*scalar_prod )

        t32 =( acnv(3)*b0cnv(2)
     .        +acnv(2)*b0cnv(3)
     .        -gsuper(3,2)*scalar_prod )

        t33 =( acnv(3)*b0cnv(3)
     .        +acnv(3)*b0cnv(3)
     .        -gsuper(3,3)*scalar_prod )

        if (flag /= 0) then
          t31 = t31/jac
          if (.not.alt_eom) t32 = t32/jac
          t33 = t33/jac
        endif

c     End program

      end subroutine tensor_z

c     findPsib_diag
c     #####################################################################
      subroutine findPsib_diag

      implicit none

      integer(4) :: iig

      cnv = div_tensor(i,j,k,nxx,nyy,nzz,igr,igr,igr,alt_eom
     .                ,tensor_x,tensor_y,tensor_z,vol=.false.)
      iig = 3*(ijk-1) + iisig - 1
      cnv = matmul(b_diag(:,iig+1:iig+3),cnv)

cc      cnv = div_tensor(i,j,k,nxx,nyy,nzz,igr,igr,igr,alt_eom
cc     .                ,tensor_x_diag,tensor_y_diag,tensor_z_diag
cc     .                ,vol=.false.)
cc      cnv = div_tensor(i,j,k,nxx,nyy,nzz,igr,igr,igr,.false.
cc     .                ,tensor_x_diag,tensor_y_diag,tensor_z_diag
cc                      ,vol=.false.)


      if (si_car) then
        call transformVectorToCartesian(i,j,k,igr,igr,igr
     .                                 ,cnv(1),cnv(2),cnv(3),.false.
     .                                 ,car(1),car(2),car(3))
        psib = alpha**2*car*vol**2
      elseif (is_cnv) then
        psib = alpha**2*cnv*vol**2
      else
        call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
        psib = alpha**2*cov*vol**2
      endif

      end subroutine findPsib_diag

c     tensor_x_diag
c     #############################################################
      subroutine tensor_x_diag(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer(4) :: ip,igrid,iig
        real(8)    :: jac,jac0,jacp,gsuper(3,3)
     .               ,acnv(3),acnvp(3),b0cov(3),b0cnv(3),scalar_prod

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

        if (isSP(i+1,j,k,igx,igy,igz)) jac = 1d-10

        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = gmetric%grid(igrid)%jac(ip,j,k)
          jac0 = gmetric%grid(igrid)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        if (flag /= 0) then
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),1,igrid)
          iig = 3*(ijk-1) + iisig - 1
          acnv = matmul(b_diag(:,iig+1:iig+3),acnv)
        else !This part does not contribute due to antisymmetry
cc          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
cc     .                      ,acnv(1),acnv(2),acnv(3),0,igrid)
          t11 = 0d0
          t12 = 0d0
          t13 = 0d0
          return
        endif

        b0cnv(1) = 0.5*(b0_cnv(ip,j,k,1)/jacp+b0_cnv(i,j,k,1)/jac0)*jac
        b0cnv(2) = 0.5*(b0_cnv(ip,j,k,2)     +b0_cnv(i,j,k,2))
        b0cnv(3) = 0.5*(b0_cnv(ip,j,k,3)/jacp+b0_cnv(i,j,k,3)/jac0)*jac

cc        b0cov(1) = 0.5*(b0_cov(ip,j,k,1)     +b0_cov(i,j,k,1)     )
cc        b0cov(2) = 0.5*(b0_cov(ip,j,k,2)/jacp+b0_cov(i,j,k,2)/jac0)*jac
cc        b0cov(3) = 0.5*(b0_cov(ip,j,k,3)     +b0_cov(i,j,k,3)     )

        if (flag /= 0) then
          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=1)
        else
          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=0)
        endif

        scalar_prod = dot_product(acnv,b0cov)

        t11 =( acnv(1)*b0cnv(1)
     .        +acnv(1)*b0cnv(1)
     .        -gsuper(1,1)*scalar_prod )

        t12 =( acnv(1)*b0cnv(2)
     .        +acnv(2)*b0cnv(1)
     .        -gsuper(1,2)*scalar_prod )

        t13 =( acnv(1)*b0cnv(3)
     .        +acnv(3)*b0cnv(1)
     .        -gsuper(1,3)*scalar_prod )

        if (flag /= 0) then
          t11 = t11/jac
          if (.not.alt_eom) t12 = t12/jac
          t13 = t13/jac
        endif

c     End program

      end subroutine tensor_x_diag

c     tensor_y_diag
c     #############################################################
      subroutine tensor_y_diag(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EM SI operator
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer(4) :: jp,igrid,iig
        real(8)    :: jac,gsuper(3,3),acnv(3),b0cov(3),b0cnv(3)
     .                ,scalar_prod

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
     .               +gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),2,igrid)
          iig = 3*(ijk-1) + iisig - 1
          acnv = matmul(b_diag(:,iig+1:iig+3),acnv)
        else  !This part does not contribute due to antisymmetry
cc          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
cc     .                      ,acnv(1),acnv(2),acnv(3),0,igrid)
          t21 = 0d0
          t22 = 0d0
          t23 = 0d0
          return
        endif

        b0cnv = 0.5*(b0_cnv(i,jp,k,:)+b0_cnv(i,j,k,:))
        b0cov = 0.5*(b0_cov(i,jp,k,:)+b0_cov(i,j,k,:))

        scalar_prod = dot_product(acnv,b0cov)

        t21 =( acnv(2)*b0cnv(1)
     .        +acnv(1)*b0cnv(2)
     .        -gsuper(2,1)*scalar_prod )

        t22 =( acnv(2)*b0cnv(2)
     .        +acnv(2)*b0cnv(2)
     .        -gsuper(2,2)*scalar_prod )

        t23 =( acnv(2)*b0cnv(3)
     .        +acnv(3)*b0cnv(2)
     .        -gsuper(2,3)*scalar_prod )

        if (flag /= 0) then
          t21 = t21/jac
          if (.not.alt_eom) t22 = t22/jac
          t23 = t23/jac
        endif

c     End program

      end subroutine tensor_y_diag

c     tensor_z_diag
c     #############################################################
      subroutine tensor_z_diag(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EM SI operator
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer(4) :: kp,igrid,iig
        real(8)    :: jac,gsuper(3,3),acnv(3),b0cov(3),b0cnv(3)
     .               ,scalar_prod

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
     .               +gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),3,igrid)
          iig = 3*(ijk-1) + iisig - 1
          acnv = matmul(b_diag(:,iig+1:iig+3),acnv)
        else  !This part does not contribute due to antisymmetry
cc          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
cc     .                      ,acnv(1),acnv(2),acnv(3),0,igrid)
          t31 = 0d0
          t32 = 0d0
          t33 = 0d0
          return
        endif

        b0cnv = 0.5*(b0_cnv(i,j,kp,:)+b0_cnv(i,j,k,:))
        b0cov = 0.5*(b0_cov(i,j,kp,:)+b0_cov(i,j,k,:))

        scalar_prod = dot_product(acnv,b0cov)

        t31 =( acnv(3)*b0cnv(1)
     .        +acnv(1)*b0cnv(3)
     .        -gsuper(3,1)*scalar_prod )

        t32 =( acnv(3)*b0cnv(2)
     .        +acnv(2)*b0cnv(3)
     .        -gsuper(3,2)*scalar_prod )

        t33 =( acnv(3)*b0cnv(3)
     .        +acnv(3)*b0cnv(3)
     .        -gsuper(3,3)*scalar_prod )

        if (flag /= 0) then
          t31 = t31/jac
          if (.not.alt_eom) t32 = t32/jac
          t33 = t33/jac
        endif

c     End program

      end subroutine tensor_z_diag

      end module v_mtvc_mod

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

      use nlfunction_setup

      use mg_internal

      use v_mtvc_mod

      implicit none

c Call variables

      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz

c Local variables

c Begin program

      igr = igrid

      is_cnv = .true.

      call allocPointers(neq,fpointers)

      isig  = MGgrid%istartp(igrid)

      nxx = grid_params%nxv(igrid)
      nyy = grid_params%nyv(igrid)
      nzz = grid_params%nzv(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

c Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz,igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Define pointers to MG arrays

      rho0   => grho0  %grid(igrid)%array

      pp0    => gp0    %grid(igrid)%array

      v0_cnv => gv0    %grid(igrid)%array

      b0_cnv => gb0    %grid(igrid)%array

      b0_cov => gb0_cov%grid(igrid)%array

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dv(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
     .            ,icomp=IVX,is_cnv=is_cnv,is_vec=.not.si_car
     .            ,result_is_vec=.true.,iorder=3)

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
            jac    = gmetric%grid(igrid)%jac(i,j,k )

            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + 1d-10

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (vol_wgt) then
              vol = volume(i,j,k,igrid,igrid,igrid)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************

            call findPsiv

            !P_si^T  ******************************

            call findPsit

            !P_si^B  ******************************

            call findPsib

            !Find Minv ****************************
cc            !!Cov vector k
cc            idx  = 1./grid_params%dxh(ig)
cc            if (nxd == 1) idx = 0d0
cc            idy  = 1./grid_params%dyh(jg)
cc            if (nyd == 1) idy = 0d0
cc            idz  = 1./grid_params%dzh(kg)
cc            if (nzd == 1) idz = 0d0
cc
cc            !!Cnv vector k
cc            call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                  ,idx,idy,idz,idxcnv,idycnv,idzcnv,.true.)
cc
cc            !!Cnv vector v
cc            vxx = abs(v0_cnv(i,j,k,1))
cc            vyy = abs(v0_cnv(i,j,k,2))
cc            vzz = abs(v0_cnv(i,j,k,3))
cc
cc            minv = 1./dt
cc     .            + alpha*(scalarProduct(i,j,k,igx,igy,igz
cc     .                               ,idx,idy,idz,vxx,vyy,vzz)
cc     .            + 0*max(eta,dd,chi)*scalarProduct(i,j,k,igx,igy,igz
cc     .                               ,idx,idy,idz,idxcnv,idycnv,idzcnv))
cc            minv = 1./minv
cc            write (*,*) minv/dt
            minv = dt

c diag ****
cc            psib = 0d0
cc            psit = 0d0
c diag ****
            do ieq=1,3
              y(neq*(ijk-1)+ieq) = psiv(ieq)
     .                           + minv*(psit(ieq) + psib(ieq))
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(dv)

      nullify(pp0,rho0)
      nullify(b0_cnv,b0_cov)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine v_mtvc

c v_mtvc2
c####################################################################
      subroutine v_mtvc2(gpos,neq,ntot,x,y,igrid,bcnd)
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

      use nlfunction_setup

      use mg_internal

      use v_mtvc_mod

      implicit none

c Call variables

      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz

c Local variables

c Begin program

      igr = igrid

      is_cnv = .true.

      call allocPointers(neq,fpointers)

      isig  = MGgrid%istartp(igrid)
      iisig = istart (igrid)

      nxx = grid_params%nxv(igrid)
      nyy = grid_params%nyv(igrid)
      nzz = grid_params%nzv(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

c Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz,igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Define pointers to MG arrays

      rho0   => grho0  %grid(igrid)%array

      pp0    => gp0    %grid(igrid)%array

      v0_cnv => gv0    %grid(igrid)%array

      b0_cnv => gb0    %grid(igrid)%array

      b0_cov => gb0_cov%grid(igrid)%array

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dv(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
     .            ,icomp=IVX,is_cnv=is_cnv,is_vec=.not.si_car
     .            ,result_is_vec=.true.,iorder=3)

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
            jac    = gmetric%grid(igrid)%jac(i,j,k )

            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + 1d-10

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (vol_wgt) then
              vol = volume(i,j,k,igrid,igrid,igrid)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************

            call findPsiv

            !P_si^T  ******************************

            call findPsit_diag

            !P_si^B  ******************************

            call findPsib_diag

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = psiv(ieq) + psit(ieq) + psib(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(dv)

      nullify(pp0,rho0)
      nullify(b0_cnv,b0_cov)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine v_mtvc2
