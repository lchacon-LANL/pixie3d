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

      integer(4) :: i,j,k,stencil_width,order,ieq
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

        if (.not.isvec) then

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

        endif

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

c lap_mtvc
c####################################################################
      subroutine lap_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
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

      logical    :: fpointers

c Begin program

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

      nxx = grid_params%nxv(igrid)
      nyy = grid_params%nyv(igrid)
      nzz = grid_params%nzv(igrid)

c Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz,igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(xarr(0:nxx+1,0:nyy+1,0:nzz+1,neq))

cc      xarr = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,xarr,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,xarr,bcnd,iorder=2)

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            y(ijk) = laplacian(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                        ,xarr(:,:,:,1),vol=vol_wgt)

          enddo
        enddo
      enddo

c End program

      deallocate(xarr)

      call deallocPointers(fpointers)

      end subroutine lap_mtvc

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

            y(ijk) = laplacian(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                        ,xarr(:,:,:,1),vol=vol_wgt)

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

      real(8)    :: upwind,nu2,dvol,lap

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
              dvol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            if (chi > 0d0) then
              lap = laplacian(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                          ,dtmp(:,:,:,1),vol=vol_wgt)
            else
              lap = 0d0
            endif

            y(ijk) = ((cnp/dt + alpha*(gamma-1.)*mgdivV0(ijkg))*x(ijk)
     .            + alpha*upwind )*dvol - alpha*chi*lap

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
              dvol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            !NC advection

cc            upwind = .5*(v0_cnv(i,j,k,1)+abs(v0_cnv(i,j,k,1)))
cc     .                 *( drho(i  ,j,k,1) - drho(i-1,j,k,1) )/dx(ig-1)
cc     .              +.5*(v0_cnv(i,j,k,1)-abs(v0_cnv(i,j,k,1)))
cc     .                 *( drho(i+1,j,k,1) - drho(i  ,j,k,1) )/dx(ig)
cc     .              +.5*(v0_cnv(i,j,k,2)+abs(v0_cnv(i,j,k,2)))
cc     .                 *( drho(i,j  ,k,1) - drho(i,j-1,k,1) )/dy(jg-1)
cc     .              +.5*(v0_cnv(i,j,k,2)-abs(v0_cnv(i,j,k,2)))
cc     .                 *( drho(i,j+1,k,1) - drho(i,j  ,k,1) )/dy(jg)
cc     .              +.5*(v0_cnv(i,j,k,3)+abs(v0_cnv(i,j,k,3)))
cc     .                 *( drho(i,j,k  ,1) - drho(i,j,k-1,1) )/dz(kg-1)
cc     .              +.5*(v0_cnv(i,j,k,3)-abs(v0_cnv(i,j,k,3)))
cc     .                 *( drho(i,j,k+1,1) - drho(i,j,k  ,1) )/dz(kg)
cc
cc            upwind = upwind/jac
cc
cc            y(ijk) = ( (cnp/dt + alpha*mgdivV0(ijkg))*x(ijk)
cc     .            + alpha*upwind )*dvol
cc     .            - alpha*dd
cc     .                *laplacian(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
cc     .                          ,drho(:,:,:,1),vol=vol_wgt)

            !Conservative advection
            upwind = c_advec(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                    ,v0_cnv(:,:,:,1)
     .                    ,v0_cnv(:,:,:,2)
     .                    ,v0_cnv(:,:,:,3)
     .                    ,drho  (:,:,:,1)
     .                    ,vol=vol_wgt,upwind=.true.,sp=bcSP())

            ! Diffusion
            if (dd > 0d0) then
              lap = laplacian(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                          ,drho(:,:,:,1),vol=vol_wgt)
            else
              lap = 0d0
            endif

            y(ijk) = cnp*drho(i,j,k,1)/dt*dvol + alpha*(upwind - dd*lap)

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
cccc            call curl_bxv(i,j,k,nxx,nyy,nzz,v0_cnv,db
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
cc     .            ,icomp=IBX,is_cnv=is_cnv,is_vec=.true.,iorder=2)
     .            ,icomp=IBX,is_cnv=is_cnv,is_vec=.not.si_car
     .            ,result_is_vec=.true.,iorder=2)

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

cc            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + SP_flsv
cc            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + 1d-3

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            if (vol_wgt) then
              vol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              vol = 1d0
            endif

            veclap = veclaplacian(i,j,k,nxx,nyy,nzz
     .                           ,igrid,igrid,igrid,db
     .                           ,alt_eom_b,vol=.false.)

            etal   = res(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid)

            if (gm_smooth) then

              cnv = curl_bxv(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                      ,v0_cnv,db,0)

            else

              !First component
              flxjp = 0.5/(jac+jacjp)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,j ,k,1)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))         
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,jp,k,1))
     .               -0.5*(v0_cnv(i,j ,k,1)*db(i,j ,k,2)/jac
     .                    +v0_cnv(i,jp,k,1)*db(i,jp,k,2)/jacjp)
              flxjm = 0.5/(jac+jacjm)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,jm,k,1)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,j ,k,1))
     .               -0.5*(v0_cnv(i,j ,k,1)*db(i,j ,k,2)/jac
     .                    +v0_cnv(i,jm,k,1)*db(i,jm,k,2)/jacjm)

              flxkp = 0.5/(jac+jackp)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,k ,1)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,kp,1))
     .               -0.5*(v0_cnv(i,j,k ,1)*db(i,j,k ,3)/jac
     .                    +v0_cnv(i,j,kp,1)*db(i,j,kp,3)/jackp)
              flxkm = 0.5/(jac+jackm)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,km,1)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,k ,1))
     .               -0.5*(v0_cnv(i,j,k ,1)*db(i,j,k ,3)/jac
     .                    +v0_cnv(i,j,km,1)*db(i,j,km,3)/jackm)

              cnv(1) =  (flxjp-flxjm)/dyh(jg)
     .                 +(flxkp-flxkm)/dzh(kg)

              !Second component
              flxip = 0.5/(jac+jacip)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(i ,j,k,2)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(ip,j,k,2))
     .               -0.5*(v0_cnv(i ,j,k,2)*db(i ,j,k,1)/jac
     .                    +v0_cnv(ip,j,k,2)*db(ip,j,k,1)/jacip)

              if (isSP(i,j,k,igrid,igrid,igrid)) then
                flxim = 0d0
              else
                flxim = 0.5/(jac+jacim+ SP_flsv)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(im,j,k,2)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(i ,j,k,2))
     .               -0.5*(v0_cnv(i ,j,k,2)*db(i ,j,k,1)/jac
     .                    +v0_cnv(im,j,k,2)*db(im,j,k,1)/jacim)
              endif

              flxkp = 0.5/(jac+jackp)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,k ,2)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,kp,3)) )*db(i,j,kp,2))
     .               -0.5*(v0_cnv(i,j,k ,2)*db(i,j,k ,3)/jac
     .                    +v0_cnv(i,j,kp,2)*db(i,j,kp,3)/jackp)
              flxkm = 0.5/(jac+jackm)*(
     .           (    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))
     .            +abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,km,2)
     .          +(    (v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3))          
     .            -abs(v0_cnv(i,j,k,3)+v0_cnv(i,j,km,3)) )*db(i,j,k ,2))
     .               -0.5*(v0_cnv(i,j,k ,2)*db(i,j,k ,3)/jac
     .                    +v0_cnv(i,j,km,2)*db(i,j,km,3)/jackm)

              cnv(2) =  (flxip-flxim)/dxh(ig)
     .                 +(flxkp-flxkm)/dzh(kg)

              !Third component
              flxip = 0.5/(jac+jacip)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(i ,j,k,3)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(ip,j,k,1)) )*db(ip,j,k,3))
     .               -0.5*(v0_cnv(i ,j,k,3)*db(i ,j,k,1)/jac
     .                    +v0_cnv(ip,j,k,3)*db(ip,j,k,1)/jacip)

              if (isSP(i,j,k,igrid,igrid,igrid)) then
                flxim = 0d0
              else
                flxim = 0.5/(jac+jacim)*(
     .           (    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))
     .            +abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(im,j,k,3)
     .          +(    (v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1))          
     .            -abs(v0_cnv(i,j,k,1)+v0_cnv(im,j,k,1)) )*db(i ,j,k,3))
     .               -0.5*(v0_cnv(i ,j,k,3)*db(i ,j,k,1)/jac
     .                    +v0_cnv(im,j,k,3)*db(im,j,k,1)/jacim)
              endif

              flxjp = 0.5/(jac+jacjp)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,j ,k,3)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jp,k,2)) )*db(i,jp,k,3))
     .               -0.5*(v0_cnv(i,j ,k,3)*db(i,j ,k,2)/jac
     .                    +v0_cnv(i,jp,k,3)*db(i,jp,k,2)/jacjp)
              flxjm = 0.5/(jac+jacjm)*(
     .           (    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))
     .            +abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,jm,k,3)
     .          +(    (v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2))          
     .            -abs(v0_cnv(i,j,k,2)+v0_cnv(i,jm,k,2)) )*db(i,j ,k,3))
     .               -0.5*(v0_cnv(i,j ,k,3)*db(i,j ,k,2)/jac
     .                    +v0_cnv(i,jm,k,3)*db(i,jm,k,2)/jacjm)

              cnv(3) =  (flxip-flxim)/dxh(ig)
     .                 +(flxjp-flxjm)/dyh(jg)

            endif

            cnv=(cnp*db(i,j,k,:)/dt + alpha*cnv - alpha*etal*veclap)*vol

            if (si_car) then
              call transformVectorToCartesian
     .              (i,j,k,igrid,igrid,igrid
     .               ,cnv(1),cnv(2),cnv(3)
     .               ,.false.
     .               ,car(1),car(2),car(3))

              do ieq=1,3
                y(neq*(ijk-1)+ieq) = car(ieq)
              enddo
            else
              do ieq=1,3
                y(neq*(ijk-1)+ieq) = cnv(ieq)
              enddo
            endif

          enddo
        enddo
      enddo

c End program

      deallocate(db)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine b_mtvc

c a_mtvc
c####################################################################
      subroutine a_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the Vector Potential..
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

      real(8),allocatable,dimension(:,:,:,:) :: da,da_cnv
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      logical    :: fpointers,is_cnv

c Begin program

      is_cnv = .false.

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

      allocate(da    (0:nxx+1,0:nyy+1,0:nzz+1,neq))
cc      allocate(da_cnv(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      da = 0d0
cc      da_cnv = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,da,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,da,bcnd
cc     .            ,icomp=IAX,is_cnv=is_cnv,is_vec=.true.,iorder=2)
     .            ,icomp=IAX,is_cnv=is_cnv,is_vec=.not.si_car
     .            ,result_is_vec=.true.,iorder=2)

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

            jac    = gmetric%grid(igrid)%jac(i,j,k)

cc            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + 1d-10

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            if (vol_wgt) then
              vol = volume(i,j,k,igrid,igrid,igrid)
            else
              vol = 1d0
            endif

            !Resistive part
            veclap = veclap_cov(i,j,k,nxx,nyy,nzz
     .                         ,igrid,igrid,igrid,da
     .                         ,vol=.false.)

cc            cnv = veclaplacian(i,j,k,nxx,nyy,nzz
cc     .                           ,igrid,igrid,igrid,da_cnv
cc     .                           ,.false.,vol=.false.)
cc
cc            call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,veclap(1),veclap(2),veclap(3)
cc     .                        ,cnv   (1),cnv   (2),cnv   (3)
cc     .                        ,.false.)

            etal   = res(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid)

            !Convective part

cc            if (i == imax) write (*,*) 'ip',da(i,j,k,3)+da(ip,j,k,3)
cc            if (i == imin) write (*,*) 'im',da(i,j,k,3)+da(im,j,k,3)

            flxip = -v0_cnv(i,j,k,2)*0.5*(da(ip,j,k,2)+da(i,j,k,2))
     .              -v0_cnv(i,j,k,3)*0.5*(da(ip,j,k,3)+da(i,j,k,3))

            flxim = -v0_cnv(i,j,k,2)*0.5*(da(im,j,k,2)+da(i,j,k,2))
     .              -v0_cnv(i,j,k,3)*0.5*(da(im,j,k,3)+da(i,j,k,3))

cc            flxip =-0.5*(    v0_cnv(i,j,k,2)
cc     .                  +abs(v0_cnv(i,j,k,2)))*da(i ,j,k,2)
cc     .             -0.5*(    v0_cnv(i,j,k,2)
cc     .                  -abs(v0_cnv(i,j,k,2)))*da(ip,j,k,2)
cc     .             -0.5*(    v0_cnv(i,j,k,3)
cc     .                  +abs(v0_cnv(i,j,k,3)))*da(i ,j,k,3)
cc     .             -0.5*(    v0_cnv(i,j,k,3)
cc     .                  -abs(v0_cnv(i,j,k,3)))*da(ip,j,k,3)
cc            flxim =-0.5*(    v0_cnv(i,j,k,2)
cc     .                  +abs(v0_cnv(i,j,k,2)))*da(im,j,k,2)
cc     .             -0.5*(    v0_cnv(i,j,k,2)
cc     .                  -abs(v0_cnv(i,j,k,2)))*da(i ,j,k,2)
cc     .             -0.5*(    v0_cnv(i,j,k,3)
cc     .                  +abs(v0_cnv(i,j,k,3)))*da(im,j,k,3)
cc     .             -0.5*(    v0_cnv(i,j,k,3)
cc     .                  -abs(v0_cnv(i,j,k,3)))*da(i ,j,k,3)

            flxjp = 0.5*(    v0_cnv(i,j,k,2)
     .                  +abs(v0_cnv(i,j,k,2)))*da(i,j ,k,1)
     .             +0.5*(    v0_cnv(i,j,k,2)
     .                  -abs(v0_cnv(i,j,k,2)))*da(i,jp,k,1)
            flxjm = 0.5*(    v0_cnv(i,j,k,2)
     .                  +abs(v0_cnv(i,j,k,2)))*da(i,jm,k,1)
     .             +0.5*(    v0_cnv(i,j,k,2)
     .                  -abs(v0_cnv(i,j,k,2)))*da(i,j ,k,1)

            flxkp = 0.5*(    v0_cnv(i,j,k,3)
     .                  +abs(v0_cnv(i,j,k,3)))*da(i,j,k ,1)
     .             +0.5*(    v0_cnv(i,j,k,3)            
     .                  -abs(v0_cnv(i,j,k,3)))*da(i,j,kp,1)
            flxkm = 0.5*(    v0_cnv(i,j,k,3)            
     .                  +abs(v0_cnv(i,j,k,3)))*da(i,j,km,1)
     .             +0.5*(    v0_cnv(i,j,k,3)            
     .                  -abs(v0_cnv(i,j,k,3)))*da(i,j,k ,1)

            cov(1) =  (flxip-flxim)/dxh(ig)
     .               +(flxjp-flxjm)/dyh(jg)
     .               +(flxkp-flxkm)/dzh(kg)


            flxip = 0.5*(    v0_cnv(i,j,k,1)
     .                  +abs(v0_cnv(i,j,k,1)))*da(i ,j,k,2)
     .             +0.5*(    v0_cnv(i,j,k,1)        
     .                  -abs(v0_cnv(i,j,k,1)))*da(ip,j,k,2)
            flxim = 0.5*(    v0_cnv(i,j,k,1)        
     .                  +abs(v0_cnv(i,j,k,1)))*da(im,j,k,2)
     .             +0.5*(    v0_cnv(i,j,k,1)        
     .                  -abs(v0_cnv(i,j,k,1)))*da(i ,j,k,2)

            flxjp = -v0_cnv(i,j,k,1)*0.5*(da(i,jp,k,1)+da(i,j,k,1))
     .              -v0_cnv(i,j,k,3)*0.5*(da(i,jp,k,3)+da(i,j,k,3))

            flxjm = -v0_cnv(i,j,k,1)*0.5*(da(i,jm,k,1)+da(i,j,k,1))
     .              -v0_cnv(i,j,k,3)*0.5*(da(i,jm,k,3)+da(i,j,k,3))

            flxkp = 0.5*(    v0_cnv(i,j,k,3)
     .                  +abs(v0_cnv(i,j,k,3)))*da(i,j,k ,2)
     .             +0.5*(    v0_cnv(i,j,k,3)            
     .                  -abs(v0_cnv(i,j,k,3)))*da(i,j,kp,2)
            flxkm = 0.5*(    v0_cnv(i,j,k,3)            
     .                  +abs(v0_cnv(i,j,k,3)))*da(i,j,km,2)
     .             +0.5*(    v0_cnv(i,j,k,3)            
     .                  -abs(v0_cnv(i,j,k,3)))*da(i,j,k ,2)

            cov(2) =  (flxip-flxim)/dxh(ig)
     .               +(flxjp-flxjm)/dyh(jg)
     .               +(flxkp-flxkm)/dzh(kg)


            flxip = 0.5*(    v0_cnv(i,j,k,1)
     .                  +abs(v0_cnv(i,j,k,1)))*da(i ,j,k,3)
     .             +0.5*(    v0_cnv(i,j,k,1)        
     .                  -abs(v0_cnv(i,j,k,1)))*da(ip,j,k,3)
            flxim = 0.5*(    v0_cnv(i,j,k,1)        
     .                  +abs(v0_cnv(i,j,k,1)))*da(im,j,k,3)
     .             +0.5*(    v0_cnv(i,j,k,1)        
     .                  -abs(v0_cnv(i,j,k,1)))*da(i ,j,k,3)

            flxjp = 0.5*(    v0_cnv(i,j,k,2)
     .                  +abs(v0_cnv(i,j,k,2)))*da(i,j ,k,3)
     .             +0.5*(    v0_cnv(i,j,k,2)
     .                  -abs(v0_cnv(i,j,k,2)))*da(i,jp,k,3)
            flxjm = 0.5*(    v0_cnv(i,j,k,2)
     .                  +abs(v0_cnv(i,j,k,2)))*da(i,jm,k,3)
     .             +0.5*(    v0_cnv(i,j,k,2)
     .                  -abs(v0_cnv(i,j,k,2)))*da(i,j ,k,3)

            flxkp = -v0_cnv(i,j,k,1)*0.5*(da(i,j,kp,1)+da(i,j,k,1))
     .              -v0_cnv(i,j,k,2)*0.5*(da(i,j,kp,2)+da(i,j,k,2))
                                                   
            flxkm = -v0_cnv(i,j,k,1)*0.5*(da(i,j,km,1)+da(i,j,k,1))
     .              -v0_cnv(i,j,k,2)*0.5*(da(i,j,km,2)+da(i,j,k,2))

            cov(3) =  (flxip-flxim)/dxh(ig)
     .               +(flxjp-flxjm)/dyh(jg)
     .               +(flxkp-flxkm)/dzh(kg)

            cov=(cnp*da(i,j,k,:)/dt + alpha*(cov/jac - etal*veclap))*vol

            if (si_car) then
              call transformVectorToCartesian
     .              (i,j,k,igrid,igrid,igrid
     .               ,cov(1),cov(2),cov(3)
     .               ,.true.
     .               ,car(1),car(2),car(3))

              do ieq=1,3
                y(neq*(ijk-1)+ieq) = car(ieq)
              enddo
            else
              do ieq=1,3
                y(neq*(ijk-1)+ieq) = cov(ieq)
              enddo
            endif

          enddo
        enddo
      enddo

c End program

      deallocate(da)
cc      deallocate(da_cnv)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine a_mtvc

c module v_mtvc_mod
c ##################################################################
      module v_mtvc_mod

      use precond_variables

      real(8),allocatable,dimension(:,:,:,:) :: dv,dv_cov

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
      subroutine findPsiv(da)

c     ---------------------------------------------------------------------
c     Find advection/diffusion component of SI operator
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: da(0:nxx+1,0:nyy+1,0:nzz+1,3)

c     Begin program

        psiv = diag_v(i,j,k,nxx,nyy,nzz,igr,igr,igr,v0_cnv,da)

      end subroutine findPsiv

c     div_upwd
c     #####################################################################
      function div_upwd(da,half_elem)

        integer(4) :: half_elem
        real(8)    :: da(0:nxx+1,0:nyy+1,0:nzz+1,*),div_upwd

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

        axp = da(ip,j,k,1)*rho0(ip,j,k,1)
        axm = da(im,j,k,1)*rho0(im,j,k,1)
        ayp = da(i,jp,k,2)*rho0(i,jp,k,1)
        aym = da(i,jm,k,2)*rho0(i,jm,k,1)
        azp = da(i,j,kp,3)*rho0(i,j,kp,1)
        azm = da(i,j,km,3)*rho0(i,j,km,1)

        div_upwd = ( (axp-axm)/dxx
     .              +(ayp-aym)/dyy
     .              +(azp-azm)/dzz )/jac
      
      end function div_upwd

c     findPsit
c     #####################################################################
      subroutine findPsit

      implicit none

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
     .             ,divip,divim,divjp,divjm,divkp,divkm
     .             ,coeff

c     Fluxes at faces for calculation of grad(dv.grad(p0))

      flxip =( (dv(i,j,k,1)/jac+dv(ip,j,k,1)/jacip)
     .           *(pp0(ip,j,k,1)-pp0(i,j,k,1))/dx(ig)
     .        +(dv(i,j,k,2)/jac+dv(ip,j,k,2)/jacip)
     .           *(pp0(ip,jp,k,1)-pp0(ip,jm,k,1)
     .            +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(ip,j,k,3)/jacip)
     .           *(pp0(ip,j,kp,1)-pp0(ip,j,km,1)
     .            +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4. )
     .        *0.5

      flxim =( (dv(i,j,k,1)/jac+dv(im,j,k,1)/jacim)
     .            *(pp0(i ,j,k,1)-pp0(im,j,k,1))/dx(ig-1)
     .        +(dv(i,j,k,2)/jac+dv(im,j,k,2)/jacim)
     .            *(pp0(im,jp,k,1)-pp0(im,jm,k,1)
     .             +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(im,j,k,3)/jacim)
     .            *(pp0(im,j,kp,1)-pp0(im,j,km,1)
     .             +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4.)
     .        *0.5

      flxjp =( (dv(i,j,k,1)/jac+dv(i,jp,k,1)/jacjp)
     .            *(pp0(ip,jp,k,1)-pp0(im,jp,k,1)
     .             +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,jp,k,2)/jacjp)
     .            *(pp0(i,jp,k,1)-pp0(i,j,k,1))/dy(jg)
     .        +(dv(i,j,k,3)/jac+dv(i,jp,k,3)/jacjp)
     .            *(pp0(i,jp,kp,1)-pp0(i,jp,km,1)
     .             +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/dzh(kg)/4.)
     .        *0.5
      flxjm =( (dv(i,j,k,1)/jac+dv(i,jm,k,1)/jacjm)
     .            *(pp0(ip,jm,k,1)-pp0(im,jm,k,1)
     .             +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,jm,k,2)/jacjm)
     .            *(pp0(i,j ,k,1)-pp0(i,jm,k,1))/dy(jg-1)
     .        +(dv(i,j,k,3)/jac+dv(i,jm,k,3)/jacjm)
     .            *(pp0(i,jm,kp,1)-pp0(i,jm,km,1)
     .             +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/dzh(kg)/4.)
     .        *0.5

      flxkp =( (dv(i,j,k,1)/jac+dv(i,j,kp,1)/jackp)
     .            *(pp0(ip,j,kp,1)-pp0(im,j,kp,1)
     .             +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,j,kp,2)/jackp)
     .            *(pp0(i,jp,kp,1)-pp0(i,jm,kp,1)
     .             +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(i,j,kp,3)/jackp)
     .            *(pp0(i,j,kp,1)-pp0(i,j,k,1))/dz(kg) )
     .        *0.5
      flxkm =( (dv(i,j,k,1)/jac+dv(i,j,km,1)/jackm)
     .            *(pp0(ip,j,km,1)-pp0(im,j,km,1)
     .             +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,j,km,2)/jackm)
     .            *(pp0(i,jp,km,1)-pp0(i,jm,km,1)
     .             +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(i,j,km,3)/jackm)
     .            *(pp0(i,j,k ,1)-pp0(i,j,km,1))/dz(kg-1) )
     .        *0.5

ccc     Fluxes at faces for calculation of grad(gamma*p0*div(dv))
cc
cc      !!Divergence at faces i+-1/2, etc.
cc      divip =  (dv(ip,j ,k,1)-dv(i ,j ,k,1))/(jac+jacip)*2/dx(ig)
cc     .       +((dv(i ,jp,k,2)-dv(i ,jm,k,2))/jac
cc     .        +(dv(ip,jp,k,2)-dv(ip,jm,k,2))/jacip)/dyh(jg)/4.
cc     .       +((dv(i ,j,kp,3)-dv(i ,j,km,3))/jac
cc     .        +(dv(ip,j,kp,3)-dv(ip,j,km,3))/jacip)/dzh(kg)/4.
cc      divim = (dv(i ,j ,k,1)-dv(im,j ,k,1))/(jac+jacim)*2/dx(ig-1)
cc     .       +((dv(i ,jp,k,2)-dv(i ,jm,k,2))/jac
cc     .        +(dv(im,jp,k,2)-dv(im,jm,k,2))/jacim)/dyh(jg)/4.
cc     .       +((dv(i ,j,kp,3)-dv(i ,j,km,3))/jac
cc     .        +(dv(im,j,kp,3)-dv(im,j,km,3))/jacim)/dzh(kg)/4.
cc
cc      divjp = ((dv(ip,j ,k,1)-dv(im,j ,k,1))/jac
cc     .        +(dv(ip,jp,k,1)-dv(im,jp,k,1))/jacjp)/dxh(ig)/4.
cc     .       +(dv(i ,jp,k,2)-dv(i ,j ,k,2))/(jac+jacjp)*2/dy(jg)
cc     .       +((dv(i,j ,kp,3)-dv(i,j ,km,3))/jac
cc     .        +(dv(i,jp,kp,3)-dv(i,jp,km,3))/jacjp)/dzh(kg)/4.
cc      divjm = ((dv(ip,j ,k,1)-dv(im,j ,k,1))/jac
cc     .        +(dv(ip,jm,k,1)-dv(im,jm,k,1))/jacjm)/dxh(ig)/4.
cc     .       +(dv(i ,j ,k,2)-dv(i ,jm,k,2))/(jac+jacjm)*2/dy(jg-1)
cc     .       +((dv(i,j ,kp,3)-dv(i,j ,km,3))/jac
cc     .        +(dv(i,jm,kp,3)-dv(i,jm,km,3))/jacjm)/dzh(kg)/4.
cc
cc      divkp = ((dv(ip,j,k ,1)-dv(im,j,k ,1))/jac
cc     .        +(dv(ip,j,kp,1)-dv(im,j,kp,1))/jackp)/dxh(ig)/4.
cc     .       +((dv(i,jp,k ,2)-dv(i,jm,k ,2))/jac
cc     .        +(dv(i,jp,kp,2)-dv(i,jm,kp,2))/jackp)/dyh(jg)/4.
cc     .       +(dv(i,j ,kp,3)-dv(i,j ,k ,3))/(jac+jackp)*2/dz(kg)
cc      divkm = ((dv(ip,j,k ,1)-dv(im,j,k ,1))/jac
cc     .        +(dv(ip,j,km,1)-dv(im,j,km,1))/jackm)/dxh(ig)/4.
cc     .       +((dv(i,jp,k ,2)-dv(i,jm,k ,2))/jac
cc     .        +(dv(i,jp,km,2)-dv(i,jm,km,2))/jackm)/dyh(jg)/4.
cc     .       +(dv(i,j ,k ,3)-dv(i,j ,km,3))/(jac+jackm)*2/dz(kg-1)
cc
cc      flxip = flxip
cc     .      + gamma*(pp0(i,j,k,1)+pp0(ip,j,k,1))*divip*0.5
cc      if (.not.isSP(i,j,k,igx,igy,igz)) then
cccc        flxim = 0d0
cccc      else
cc        flxim = flxim
cc     .      + gamma*(pp0(i,j,k,1)+pp0(im,j,k,1))*divim*0.5
cc      endif
cc
cc      flxjp = flxjp
cc     .      + gamma*(pp0(i,j,k,1)+pp0(i,jp,k,1))*divjp*0.5
cc      flxjm = flxjm
cc     .      + gamma*(pp0(i,j,k,1)+pp0(i,jm,k,1))*divjm*0.5
cc
cc      flxkp = flxkp
cc     .      + gamma*(pp0(i,j,k,1)+pp0(i,j,kp,1))*divkp*0.5
cc      flxkm = flxkm
cc     .      + gamma*(pp0(i,j,k,1)+pp0(i,j,km,1))*divkm*0.5

cc      !Fluxes at faces for calculation of grad(dv.grad(p0))
cc      flxip =( (dv(i,j,k,1)+dv(ip,j,k,1))
cc     .           *(pp0(ip,j,k,1)-pp0(i,j,k,1))/dx(ig)
cc     .        +(dv(i,j,k,2)+dv(ip,j,k,2))
cc     .           *(pp0(ip,jp,k,1)-pp0(ip,jm,k,1)
cc     .            +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
cc     .        +(dv(i,j,k,3)+dv(ip,j,k,3))
cc     .           *(pp0(ip,j,kp,1)-pp0(ip,j,km,1)
cc     .            +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4. )
cc     .        /(jac+jacip)
cc
cc      if (.not.isSP(i,j,k,igx,igy,igz)) then
cc        flxim =( (dv(i,j,k,1)+dv(im,j,k,1))
cc     .              *(pp0(i ,j,k,1)-pp0(im,j,k,1))/dx(ig-1)
cc     .          +(dv(i,j,k,2)+dv(im,j,k,2))
cc     .              *(pp0(im,jp,k,1)-pp0(im,jm,k,1)
cc     .               +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
cc     .          +(dv(i,j,k,3)+dv(im,j,k,3))
cc     .              *(pp0(im,j,kp,1)-pp0(im,j,km,1)
cc     .               +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4.)
cc     .          /(jac+jacim)
cc      else
cc        flxim =( (dv(i,j,k,1)/jac+dv(im,j,k,1)/jacim)
cc     .              *(pp0(i ,j,k,1)-pp0(im,j,k,1))/dx(ig-1)
cc     .          +(dv(i,j,k,2)/jac+dv(im,j,k,2)/jacim)
cc     .              *(pp0(im,jp,k,1)-pp0(im,jm,k,1)
cc     .               +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
cc     .          +(dv(i,j,k,3)/jac+dv(im,j,k,3)/jacim)
cc     .              *(pp0(im,j,kp,1)-pp0(im,j,km,1)
cc     .               +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4.)
cc     .          *0.5
cc      endif
cc
cc      flxjp =( (dv(i,j,k,1)+dv(i,jp,k,1))
cc     .            *(pp0(ip,jp,k,1)-pp0(im,jp,k,1)
cc     .             +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/dxh(ig)/4.
cc     .        +(dv(i,j,k,2)+dv(i,jp,k,2))
cc     .            *(pp0(i,jp,k,1)-pp0(i,j,k,1))/dy(jg)
cc     .        +(dv(i,j,k,3)+dv(i,jp,k,3))
cc     .            *(pp0(i,jp,kp,1)-pp0(i,jp,km,1)
cc     .             +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/dzh(kg)/4.)
cc     .        /(jac+jacjp)
cc      flxjm =( (dv(i,j,k,1)+dv(i,jm,k,1))
cc     .            *(pp0(ip,jm,k,1)-pp0(im,jm,k,1)
cc     .             +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/dxh(ig)/4.
cc     .        +(dv(i,j,k,2)+dv(i,jm,k,2))
cc     .            *(pp0(i,j ,k,1)-pp0(i,jm,k,1))/dy(jg-1)
cc     .        +(dv(i,j,k,3)+dv(i,jm,k,3))
cc     .            *(pp0(i,jm,kp,1)-pp0(i,jm,km,1)
cc     .             +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/dzh(kg)/4.)
cc     .        /(jac+jacjm)
cc
cc      flxkp =( (dv(i,j,k,1)+dv(i,j,kp,1))
cc     .            *(pp0(ip,j,kp,1)-pp0(im,j,kp,1)
cc     .             +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/dxh(ig)/4.
cc     .        +(dv(i,j,k,2)+dv(i,j,kp,2))
cc     .            *(pp0(i,jp,kp,1)-pp0(i,jm,kp,1)
cc     .             +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/dyh(jg)/4.
cc     .        +(dv(i,j,k,3)+dv(i,j,kp,3))
cc     .            *(pp0(i,j,kp,1)-pp0(i,j,k,1))/dz(kg) )
cc     .        /(jac+jackp)
cc      flxkm =( (dv(i,j,k,1)+dv(i,j,km,1))
cc     .            *(pp0(ip,j,km,1)-pp0(im,j,km,1)
cc     .             +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/dxh(ig)/4.
cc     .        +(dv(i,j,k,2)+dv(i,j,km,2))
cc     .            *(pp0(i,jp,km,1)-pp0(i,jm,km,1)
cc     .             +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/dyh(jg)/4.
cc     .        +(dv(i,j,k,3)+dv(i,j,km,3))
cc     .            *(pp0(i,j,k ,1)-pp0(i,j,km,1))/dz(kg-1) )
cc     .        /(jac+jackm)

c      Fluxes at faces for calculation of grad(gamma*p0*div(dv))

      !!Divergence at faces i+-1/2, etc.
      divip = (dv(ip,j ,k,1)-dv(i ,j ,k,1))/dx(ig)
     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .        +dv(ip,jp,k,2)-dv(ip,jm,k,2))/dyh(jg)*0.25
     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .        +dv(ip,j,kp,3)-dv(ip,j,km,3))/dzh(kg)*0.25
      divim = (dv(i ,j ,k,1)-dv(im,j ,k,1))/dx(ig-1)
     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .        +dv(im,jp,k,2)-dv(im,jm,k,2))/dyh(jg)*0.25
     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .        +dv(im,j,kp,3)-dv(im,j,km,3))/dzh(kg)*0.25

      divjp = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .        +dv(ip,jp,k,1)-dv(im,jp,k,1))/dxh(ig)*0.25
     .       +(dv(i ,jp,k,2)-dv(i ,j ,k,2))/dy(jg)
     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .        +dv(i,jp,kp,3)-dv(i,jp,km,3))/dzh(kg)*0.25
      divjm = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .        +dv(ip,jm,k,1)-dv(im,jm,k,1))/dxh(ig)*0.25
     .       +(dv(i ,j ,k,2)-dv(i ,jm,k,2))/dy(jg-1)
     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .        +dv(i,jm,kp,3)-dv(i,jm,km,3))/dzh(kg)*0.25

      divkp = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .        +dv(ip,j,kp,1)-dv(im,j,kp,1))/dxh(ig)*0.25
     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .        +dv(i,jp,kp,2)-dv(i,jm,kp,2))/dyh(jg)*0.25
     .       +(dv(i,j ,kp,3)-dv(i,j ,k ,3))/dz(kg)
      divkm = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .        +dv(ip,j,km,1)-dv(im,j,km,1))/dxh(ig)*0.25
     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .        +dv(i,jp,km,2)-dv(i,jm,km,2))/dyh(jg)*0.25
     .       +(dv(i,j ,k ,3)-dv(i,j ,km,3))/dz(kg-1)

      flxip = flxip
     .      + gamma*(pp0(i,j,k,1)+pp0(ip,j,k,1))*divip/(jac+jacip)
      if (.not.isSP(i,j,k,igx,igy,igz)) then
        flxim = flxim
     .      + gamma*(pp0(i,j,k,1)+pp0(im,j,k,1))*divim/(jac+jacim)
      endif

      flxjp = flxjp
     .      + gamma*(pp0(i,j,k,1)+pp0(i,jp,k,1))*divjp/(jac+jacjp)
      flxjm = flxjm
     .      + gamma*(pp0(i,j,k,1)+pp0(i,jm,k,1))*divjm/(jac+jacjm)

      flxkp = flxkp
     .      + gamma*(pp0(i,j,k,1)+pp0(i,j,kp,1))*divkp/(jac+jackp)
      flxkm = flxkm
     .      + gamma*(pp0(i,j,k,1)+pp0(i,j,km,1))*divkm/(jac+jackm)

cc      if (isSP(i,j,k,igx,igy,igz)) flxim = 0d0

c     Transform to contravariant

      cov(1) = (flxip - flxim)/dxh(ig)
      cov(2) = (flxjp - flxjm)/dyh(jg)
      cov(3) = (flxkp - flxkm)/dzh(kg)

      call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                            ,cov(1),cov(2),cov(3)
     .                            ,cnv(1),cnv(2),cnv(3),is_cnv)

      coeff = (alpha**2/cnp + k_si*alpha)

      psit = -coeff*cnv

      end subroutine findPsit

c     findPsib
c     #############################################################
      subroutine findPsib

      implicit none

      real(8) :: coeff

      cnv = div_tensor(i,j,k,nxx,nyy,nzz,igr,igr,igr,alt_eom
cc      cnv = div_tensor(i,j,k,nxx,nyy,nzz,igr,igr,igr,.false.
     .                ,tensor_x,tensor_y,tensor_z,vol=.false.)

      coeff = (alpha**2/cnp + k_si*alpha)

      psib = coeff*cnv

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
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,gsuper(3,3)
     .               ,acnv(3),acnvp(3),b0cov(3),b0cnv(3),scalar_prod

c     Begin program

cc        write (*,*) igx,igy,igz,nxx,nyy,nzz

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

        if (isSP(i+1,j,k,igx,igy,igz).and.flag /= 0) jac = SP_flsv

        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = gmetric%grid(igrid)%jac(ip,j,k)
          jac0 = gmetric%grid(igrid)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

        if (flag /= 0) then
          acnv =
     .      vecA(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov,b0_cnv,1)
        else
          acnv =
     .      vecA(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov,b0_cnv,0)
        endif

cc        b0cnv = 0.5*(b0_cnv(ip,j,k,:)+b0_cnv(i,j,k,:))
cc        b0cov = 0.5*(b0_cov(ip,j,k,:)+b0_cov(i,j,k,:))

        b0cnv(1)=0.5*(b0_cnv(ip,j,k,1)*ijacp+b0_cnv(i,j,k,1)*ijac0)*jac
        b0cnv(2)=0.5*(b0_cnv(ip,j,k,2)      +b0_cnv(i,j,k,2))
        b0cnv(3)=0.5*(b0_cnv(ip,j,k,3)*ijacp+b0_cnv(i,j,k,3)*ijac0)*jac

cc        b0cov(1)=0.5*(b0_cov(ip,j,k,1)      +b0_cov(i,j,k,1)      )
cc        b0cov(2)=0.5*(b0_cov(ip,j,k,2)*ijacp+b0_cov(i,j,k,2)*ijac0)*jac
cc        b0cov(3)=0.5*(b0_cov(ip,j,k,3)      +b0_cov(i,j,k,3)      )

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

cc        if (isSP(i+1,j,k,igx,igy,igz).and.flag /= 0) jac = SP_flsv

        if (flag /= 0) then
          if (isSP(i+1,j,k,igx,igy,igz)) then
            t11 = 0d0
            if (.not.alt_eom) t12 = 0d0
            t13 = 0d0
          else
            ijac  = 1d0/jac
            t11 = t11*ijac
            if (.not.alt_eom) t12 = t12*ijac
            t13 = t13*ijac
          endif
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
        real(8)    :: jac,ijac,gsuper(3,3),acnv(3),b0cov(3),b0cnv(3)
     .                ,scalar_prod

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
     .               +gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        ijac  = 1d0/jac

        if (flag /= 0) then
          acnv =
     .      vecA(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov,b0_cnv,2)
        else
          acnv =
     .      vecA(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov,b0_cnv,0)
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
          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
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
        real(8)    :: jac,ijac,gsuper(3,3),acnv(3),b0cov(3),b0cnv(3)
     .               ,scalar_prod

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
     .               +gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        ijac  = 1d0/jac

        if (flag /= 0) then
          acnv =
     .      vecA(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov,b0_cnv,3)
        else
          acnv =
     .      vecA(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov,b0_cnv,0)
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
          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine tensor_z

c     gradDiv
c     #####################################################################
      function gradDiv(da) result(cnv)

      implicit none

      real(8)    :: cnv(3)
      real(8)    :: da(0:nxx+1,0:nyy+1,0:nzz+1,3)

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
     .             ,divip,divim,divjp,divjm,divkp,divkm
     .             ,coeff

cc      !Fluxes at faces for calculation of grad[div(da p0)]
cc      flxip =( (da(ip,j ,k ,1)*rho0(ip,j ,k ,1)
cc     .         -da(i ,j ,k ,1)*rho0(i ,j ,k ,1))/dx(ig)
cc     .        +(da(ip,jp,k ,2)*rho0(ip,jp,k ,1)
cc     .         -da(ip,jm,k ,2)*rho0(ip,jm,k ,1)
cc     .         +da(i ,jp,k ,2)*rho0(i ,jp,k ,1)
cc     .         -da(i ,jm,k ,2)*rho0(i ,jm,k ,1))/dyh(jg)/4.
cc     .        +(da(ip,j ,kp,3)*rho0(ip,j ,kp,1)
cc     .         -da(ip,j ,km,3)*rho0(ip,j ,km,1)
cc     .         +da(i ,j ,kp,3)*rho0(i ,j ,kp,1)
cc     .         -da(i ,j ,km,3)*rho0(i ,j ,km,1))/dzh(kg)/4.)
cc     $        /(jac+jacip)*2
cc      flxim =( (da(i ,j ,k ,1)*rho0(i ,j ,k ,1)
cc     .         -da(im,j ,k ,1)*rho0(im,j ,k ,1))/dx(ig-1)
cc     .        +(da(im,jp,k ,2)*rho0(im,jp,k ,1)
cc     .         -da(im,jm,k ,2)*rho0(im,jm,k ,1)
cc     .         +da(i ,jp,k ,2)*rho0(i ,jp,k ,1)
cc     .         -da(i ,jm,k ,2)*rho0(i ,jm,k ,1))/dyh(jg)/4.
cc     .        +(da(im,j ,kp,3)*rho0(im,j ,kp,1)
cc     .         -da(im,j ,km,3)*rho0(im,j ,km,1)
cc     .         +da(i ,j ,kp,3)*rho0(i ,j ,kp,1)
cc     .         -da(i ,j ,km,3)*rho0(i ,j ,km,1))/dzh(kg)/4.)
cc     $        /(jac+jacim)*2
cc
cc      flxjp =( (da(ip,jp,k ,1)*rho0(ip,jp,k ,1)
cc     .         -da(im,jp,k ,1)*rho0(im,jp,k ,1)
cc     .         +da(ip,j ,k ,1)*rho0(ip,j ,k ,1)
cc     .         -da(im,j ,k ,1)*rho0(im,j ,k ,1))/dxh(ig)/4.
cc     .        +(da(i ,jp,k ,2)*rho0(i ,jp,k ,1)
cc     .         -da(i ,j ,k ,2)*rho0(i ,j ,k ,1))/dy(jg)
cc     .        +(da(i ,jp,kp,3)*rho0(i ,jp,kp,1)
cc     .         -da(i ,jp,km,3)*rho0(i ,jp,km,1)
cc     .         +da(i ,j ,kp,3)*rho0(i ,j ,kp,1)
cc     .         -da(i ,j ,km,3)*rho0(i ,j ,km,1))/dzh(kg)/4.)
cc     $        /(jac+jacjp)*2
cc      flxjm =( (da(ip,jm,k ,1)*rho0(ip,jm,k ,1)
cc     .         -da(im,jm,k ,1)*rho0(im,jm,k ,1)
cc     .         +da(ip,j ,k ,1)*rho0(ip,j ,k ,1)
cc     .         -da(im,j ,k ,1)*rho0(im,j ,k ,1))/dxh(ig)/4.
cc     .        +(da(i ,j ,k ,2)*rho0(i ,j ,k ,1)
cc     .         -da(i ,jm,k ,2)*rho0(i ,jm,k ,1))/dy(jg-1)
cc     .        +(da(i ,jm,kp,3)*rho0(i ,jm,kp,1)
cc     .         -da(i ,jm,km,3)*rho0(i ,jm,km,1)
cc     .         +da(i ,j ,kp,3)*rho0(i ,j ,kp,1)
cc     .         -da(i ,j ,km,3)*rho0(i ,j ,km,1))/dzh(kg)/4.)
cc     $        /(jac+jacjm)*2
cc
cc      flxkp =( (da(ip,j ,kp,1)*rho0(ip,j ,kp,1)
cc     .         -da(im,j ,kp,1)*rho0(im,j ,kp,1)
cc     .         +da(ip,j ,k ,1)*rho0(ip,j ,k ,1)
cc     .         -da(im,j ,k ,1)*rho0(im,j ,k ,1))/dxh(ig)/4.
cc     .        +(da(i ,jp,kp,2)*rho0(i ,jp,kp,1)
cc     .         -da(i ,jm,kp,2)*rho0(i ,jm,kp,1)
cc     .         +da(i ,jp,k ,2)*rho0(i ,jp,k ,1)
cc     .         -da(i ,jm,k ,2)*rho0(i ,jm,k ,1))/dyh(jg)/4.
cc     .        +(da(i ,j ,kp,3)*rho0(i ,j ,kp,1)
cc     .         -da(i ,j ,k ,3)*rho0(i ,j ,k ,1))/dz(kg) )
cc     $        /(jac+jackp)*2
cc      flxkm =( (da(ip,j ,km,1)*rho0(ip,j ,km,1)
cc     .         -da(im,j ,km,1)*rho0(im,j ,km,1)
cc     .         +da(ip,j ,k ,1)*rho0(ip,j ,k ,1)
cc     .         -da(im,j ,k ,1)*rho0(im,j ,k ,1))/dxh(ig)/4.
cc     .        +(da(i ,jp,km,2)*rho0(i ,jp,km,1)
cc     .         -da(i ,jm,km,2)*rho0(i ,jm,km,1)
cc     .         +da(i ,jp,k ,2)*rho0(i ,jp,k ,1)
cc     .         -da(i ,jm,k ,2)*rho0(i ,jm,k ,1))/dyh(jg)/4.
cc     .        +(da(i ,j ,k ,3)*rho0(i ,j ,k ,1)
cc     .         -da(i ,j ,km,3)*rho0(i ,j ,km,1))/dz(kg-1) )
cc     $        /(jac+jackm)*2
cc
cc      cov(1) = (flxip - flxim)/dxh(ig)
cc      cov(2) = (flxjp - flxjm)/dyh(jg)
cc      cov(3) = (flxkp - flxkm)/dzh(kg)

      !Fluxes at faces for calculation of grad[div(da)]

      !!Divergence at faces i+-1/2, etc.
      divip = (da(ip,j ,k,1)-da(i ,j ,k,1))/dx(ig)
     .       +(da(i ,jp,k,2)-da(i ,jm,k,2)
     .        +da(ip,jp,k,2)-da(ip,jm,k,2))/dyh(jg)/4.
     .       +(da(i ,j,kp,3)-da(i ,j,km,3)
     .        +da(ip,j,kp,3)-da(ip,j,km,3))/dzh(kg)/4.
      divim = (da(i ,j ,k,1)-da(im,j ,k,1))/dx(ig-1)
     .       +(da(i ,jp,k,2)-da(i ,jm,k,2)
     .        +da(im,jp,k,2)-da(im,jm,k,2))/dyh(jg)/4.
     .       +(da(i ,j,kp,3)-da(i ,j,km,3)
     .        +da(im,j,kp,3)-da(im,j,km,3))/dzh(kg)/4.

      divjp = (da(ip,j ,k,1)-da(im,j ,k,1)
     .        +da(ip,jp,k,1)-da(im,jp,k,1))/dxh(ig)/4.
     .       +(da(i ,jp,k,2)-da(i ,j ,k,2))/dy(jg)
     .       +(da(i,j ,kp,3)-da(i,j ,km,3)
     .        +da(i,jp,kp,3)-da(i,jp,km,3))/dzh(kg)/4.
      divjm = (da(ip,j ,k,1)-da(im,j ,k,1)
     .        +da(ip,jm,k,1)-da(im,jm,k,1))/dxh(ig)/4.
     .       +(da(i ,j ,k,2)-da(i ,jm,k,2))/dy(jg-1)
     .       +(da(i,j ,kp,3)-da(i,j ,km,3)
     .        +da(i,jm,kp,3)-da(i,jm,km,3))/dzh(kg)/4.

      divkp = (da(ip,j,k ,1)-da(im,j,k ,1)
     .        +da(ip,j,kp,1)-da(im,j,kp,1))/dxh(ig)/4.
     .       +(da(i,jp,k ,2)-da(i,jm,k ,2)
     .        +da(i,jp,kp,2)-da(i,jm,kp,2))/dyh(jg)/4.
     .       +(da(i,j ,kp,3)-da(i,j ,k ,3))/dz(kg)
      divkm = (da(ip,j,k ,1)-da(im,j,k ,1)
     .        +da(ip,j,km,1)-da(im,j,km,1))/dxh(ig)/4.
     .       +(da(i,jp,k ,2)-da(i,jm,k ,2)
     .        +da(i,jp,km,2)-da(i,jm,km,2))/dyh(jg)/4.
     .       +(da(i,j ,k ,3)-da(i,j ,km,3))/dz(kg-1)

      cov(1) = (divip - divim)/dxh(ig)
      cov(2) = (divjp - divjm)/dyh(jg)
      cov(3) = (divkp - divkm)/dzh(kg)

      call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                            ,cov(1),cov(2),cov(3)
     .                            ,cnv(1),cnv(2),cnv(3),is_cnv)

      end function gradDiv

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

      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq),iig
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

      allocate(dv    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
     .                       ,.false.)

      !Higher-order interpolation than 1 is JB-unstable for cyl.
      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
     .            ,icomp=IVX,is_cnv=is_cnv,is_vec=.not.si_car
     .            ,result_is_vec=.true.,iorder=1,arr_cov=dv_cov)

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

cc            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + SP_flsv

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (vol_wgt) then
              vol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************

            call findPsiv(dv)

            !P_si^T  ******************************

            call findPsit

            !P_si^B  ******************************

            call findPsib

            !Diagonal scaling
             
            if (precon == 's2') then

cc              !Find Minv ****************************
cc              !!Cov vector k
cc              idx  = 1./grid_params%dxh(ig)
cc              if (nxd == 1) idx = 0d0
cc              idy  = 1./grid_params%dyh(jg)
cc              if (nyd == 1) idy = 0d0
cc              idz  = 1./grid_params%dzh(kg)
cc              if (nzd == 1) idz = 0d0
cc
cc              !!Cnv vector v
cc              vxx = abs(v0_cnv(i,j,k,1))
cc              vyy = abs(v0_cnv(i,j,k,2))
cc              vzz = abs(v0_cnv(i,j,k,3))
cc
cc              minv = cnp/dt
cc     .              + alpha*(scalarProduct(i,j,k,igx,igy,igz
cc     .                                 ,idx,idy,idz,vxx,vyy,vzz)
cc     .              + max(eta,dd,chi)*vectorNorm(i,j,k,igx,igy,igz
cc     .                               ,idx,idy,idz,.true.))
cc
cc              minv = 1./minv
cc
cc              cnv = cnv*minv

cc              --------------------------

cc              iig = 3*(ijk-1) + iisig - 1
cc
cc              cnv = matmul(transpose(mstar_diag(:,iig+1:iig+3)),cnv)
cccc              cnv = matmul(transpose(b_diag(:,iig+1:iig+3)),cnv)*vol

cc              --------------------------

cccc              cnv = cnv*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol
              psit = psit*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol
              psib = psib*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol

              cnv = (psit + psib)

cc              --------------------------

cc              psib = dt*psib
cc              psit = dt*psit

cc              cnv = (psit + psib)

cc              --------------------------

cc              psit = psit*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol
cc
cc              iig = 3*(ijk-1) + iisig - 1
cc              psib = matmul(transpose(b_diag(:,iig+1:iig+3)),psib)*vol

cc              cnv = (psit + psib)

            else

              cnv = dt*(psit + psib)

            endif

            !Form total operator

            if (nc_eom_v) cnv = cnv/rho0(i,j,k,1)

            cnv = (psiv + cnv)*vol

            if (si_car) then
              call transformVectorToCartesian
     .              (i,j,k,igrid,igrid,igrid
     .               ,cnv(1),cnv(2),cnv(3),.false.
     .               ,car(1),car(2),car(3))

              do ieq=1,3
                y(neq*(ijk-1)+ieq) = car(ieq)
              enddo
            else
              do ieq=1,3
                y(neq*(ijk-1)+ieq) = cnv(ieq)
              enddo
            endif

          enddo
        enddo
      enddo

c End program

      deallocate(dv,dv_cov)

      nullify(pp0,rho0)
      nullify(b0_cnv,b0_cov)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine v_mtvc

ccc v_mtvc_pr
ccc####################################################################
cc      subroutine v_mtvc_pr(gpos,neq,ntot,x,y,igrid,bcnd)
ccc--------------------------------------------------------------------
ccc     This subroutine calculates, for given x, y = A(psi)x  matrix-free
ccc     for the velocity SI system.
ccc     In call:
ccc      * gpos: vector index of position on the numerical grid:
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
cc      use nlfunction_setup
cc
cc      use mg_internal
cc
cc      use mgarraySetup
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq)
cc      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
cc     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz
cc
ccc Local variables
cc
cc      integer(4) :: igr,nxx,nyy,nzz,ntotf,isig
cc      real(8),allocatable,dimension(:) :: MGx,MGy
cc
cc      logical    :: fpointers
cc
ccc Begin program
cc
cc      call allocPointers(neq,fpointers)
cc
cc      if (igrid > 1) then
cc
cc        isig  = istart(igrid)
cc
ccc     Allocate MG vector
cc
cc        igr = 1
cc
cc        nxx = grid_params%nxv(igr)
cc        nyy = grid_params%nyv(igr)
cc        nzz = grid_params%nzv(igr)
cc
cc        ntotf = neq*nxx*nyy*nzz
cc
cc        allocate(MGx(2*ntotf),MGy(2*ntotf))
cc
ccc     Prolong vector to finest grid
cc
cc        MGx(isig:isig+ntot-1) = x
cc
cc        call prolongMGVector(neq,MGx,igrid,igr,0,bcnd)
cc
ccc     Perform matvec in finest grid
cc
cc        call v_mtvc(0,neq,ntotf,MGx(1:ntotf),MGy(1:ntotf),igr,bcnd)
cc
ccc     Restrict residual
cc
cc        call restrictMGVector(neq,MGy,igr,igrid,0,vol_wgt)
cc
ccc     Map coarse grid residual
cc
cc        y = MGy(isig:isig+ntot-1)
cc
ccc     Deallocate MG vectors
cc
cc        deallocate(MGx,MGy)
cc
cc      else
cc
cc        call v_mtvc(0,neq,ntot,x,y,igrid,bcnd)
cc
cc      endif
cc
cc      call deallocPointers(fpointers)
cc
cc      end subroutine v_mtvc_pr

ccc v_mtvc2
ccc####################################################################
cc      subroutine v_mtvc2(gpos,neq,ntot,x,y,igrid,bcnd)
ccc--------------------------------------------------------------------
ccc     This subroutine calculates, for given x, y = A(psi)x  matrix-free
ccc     for the velocity SI system.
ccc     In call:
ccc      * gpos: vector index of position on the numerical grid:
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
cc      use nlfunction_setup
cc
cc      use mg_internal
cc
cc      use v_mtvc_mod
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq)
cc      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
cc     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz
cc
ccc Local variables
cc
ccc Begin program
cc
cc      igr = igrid
cc
cc      is_cnv = .true.
cc
cc      call allocPointers(neq,fpointers)
cc
cc      isig  = MGgrid%istartp(igrid)
cc      iisig = istart (igrid)
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
cc      call limits(abs(gpos),nxx,nyy,nzz,igrid
cc     .           ,imin,imax,jmin,jmax,kmin,kmax)
cc
ccc Define pointers to MG arrays
cc
cc      rho0   => grho0  %grid(igrid)%array
cc
cc      pp0    => gp0    %grid(igrid)%array
cc
cc      v0_cnv => gv0    %grid(igrid)%array
cc
cc      b0_cnv => gb0    %grid(igrid)%array
cc
cc      b0_cov => gb0_cov%grid(igrid)%array
cc
ccc Map vector x to array for processing (return dv in curvilinear comp.)
cc
cc      allocate(dv    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
cc     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))
cc
cc      dv = 0d0
cc
cc      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
cc     .                       ,.false.)
cc
cc      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
cc     .            ,icomp=IVX,is_cnv=is_cnv,is_vec=.not.si_car
cc     .            ,result_is_vec=.true.,iorder=1,arr_cov=dv_cov)
cccc     .            ,result_is_vec=.true.)
cc
ccc Calculate matrix-vector product
cc
cc      do k = kmin,kmax
cc        do j = jmin,jmax
cc          do i = imin,imax
cc
cc            !Preparations
cc            ip = i+1
cc            im = i-1
cc            jp = j+1
cc            jm = j-1
cc            kp = k+1
cc            km = k-1
cc
cc            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
cc
cc            jacip  = gmetric%grid(igrid)%jac(ip,j,k)
cc            jacim  = gmetric%grid(igrid)%jac(im,j,k)
cc            jacjp  = gmetric%grid(igrid)%jac(i,jp,k)
cc            jacjm  = gmetric%grid(igrid)%jac(i,jm,k)
cc            jackp  = gmetric%grid(igrid)%jac(i,j,kp)
cc            jackm  = gmetric%grid(igrid)%jac(i,j,km)
cc            jac    = gmetric%grid(igrid)%jac(i,j,k )
cc
cccc            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + SP_flsv
cc
cc            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)
cc
cc            ijkg   = ijk + isig - 1
cc
cc            if (vol_wgt) then
cc              vol = gmetric%grid(igrid)%dvol(i,j,k)
cc            else
cc              vol = 1d0
cc            endif
cc
cc            !P_si^v  ******************************
cc
cc            call findPsiv(dv)
cc
cc            !P_si^T  ******************************
cc
cccc            call findPsit_diag
cc            call findPsit
cc
cc            !P_si^B  ******************************
cc
cccc            call findPsib_diag
cc            call findPsib
cc
cc            cnv = (psiv + psit + psib)*vol
cc
cc            if (nc_eom_v) cnv = cnv/rho0(i,j,k,1)
cc
cc            if (si_car) then
cc              call transformVectorToCartesian
cc     .              (i,j,k,igrid,igrid,igrid
cc     .               ,cnv(1),cnv(2),cnv(3),.false.
cc     .               ,car(1),car(2),car(3))
cc
cc              do ieq=1,3
cc                y(neq*(ijk-1)+ieq) = car(ieq)
cc              enddo
cc            else
cc              do ieq=1,3
cc                y(neq*(ijk-1)+ieq) = cnv(ieq)
cc              enddo
cc            endif
cc
cc          enddo
cc        enddo
cc      enddo
cc
ccc End program
cc
cc      deallocate(dv,dv_cov)
cc
cc      nullify(pp0,rho0)
cc      nullify(b0_cnv,b0_cov)
cc      nullify(v0_cnv)
cc
cc      call deallocPointers(fpointers)
cc
cc      end subroutine v_mtvc2

c v_mtvc3
c####################################################################
      subroutine v_mtvc3(gpos,neq,ntot,x,y,igrid,bcnd)
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

      real(8)    :: z(ntot)

      real(8),allocatable,dimension(:,:,:,:) :: dzz

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

      allocate(dv (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dzz(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      !Velocity array
      dv = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
     .            ,icomp=IVX,is_cnv=is_cnv,is_vec=.not.si_car
     .            ,result_is_vec=.true.,iorder=1)
cc     .            ,result_is_vec=.true.)

      !M_star*dv array
      dzz = 0d0

      call mstar_mtvc(gpos,neq,ntot,x,z,igrid,bcnd)

      call mapMGVectorToArray(max(0,gpos),neq,z,nxx,nyy,nzz,dzz,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dzz,bcnd
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

            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + SP_flsv

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (vol_wgt) then
              vol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************

            call findPsiv(dzz)
cc            call findPsiv(dv)

            !P_si^T  ******************************

            call findPsit

            !P_si^B  ******************************

            call findPsib

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = psiv(ieq) + psit(ieq) + psib(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(dv,dzz)

      nullify(pp0,rho0)
      nullify(b0_cnv,b0_cov)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine v_mtvc3

c mstar_mtvc
c####################################################################
      subroutine mstar_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the M_star operator.
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
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz
     .             ,idxcnv,idycnv,idzcnv,vxx,vyy,vzz
     .             ,diff,upwind

c Local variables

      real(8),allocatable,dimension(:,:,:,:) :: da

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

      v0_cnv => gv0    %grid(igrid)%array

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(da(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      da = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,da,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,da,bcnd
     .            ,icomp=IVX,is_cnv=is_cnv,is_vec=.not.si_car
     .            ,result_is_vec=.true.,iorder=2)
cc     .            ,result_is_vec=.true.)

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

            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + SP_flsv

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            diff = max(eta,chi,dd)

            cnv = gradDiv(da)

            do ieq=1,neq
              hex = floor(sign(1d0,v0_cnv(i,j,k,ieq)))
cc              hex = 0

              upwind = v0_cnv(i,j,k,ieq)*div_upwd(da,hex)/rho0(i,j,k,1)

              cnv(ieq) = cnp*da(i,j,k,ieq)/dt
     .                 + alpha*upwind
     .                 - alpha*diff*cnv(ieq)/rho0(i,j,k,1)

              !No volume weighing in M_star
              y(neq*(ijk-1)+ieq) = cnv(ieq)

            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(da)

      nullify(v0_cnv,rho0)

      call deallocPointers(fpointers)

      end subroutine mstar_mtvc

c hall_mtvc
c####################################################################
      subroutine hall_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
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

      b0_cnv => gb0    %grid(igrid)%array

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dv    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
     .                       ,.false.)

cc      where (bcnd == DIR) bcnd = NEU
      !Higher-order than 1 is JB-unstable for cyl.
      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
     .            ,icomp=IVX,is_cnv=is_cnv,is_vec=.not.si_car
     .            ,result_is_vec=.true.,iorder=1,arr_cov=dv_cov)

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

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (vol_wgt) then
              vol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              vol = 1d0
            endif

            cnv(1) = dv(i,j,k,1)/dt
cc     .          - di*b0_cnv(i,j,k,3)
     .          - di
     .           *((dv(ip,j,k,2)-2*dv(i,j,k,2)+dv(im,j,k,2))/dxh(ig)**2
     .            -(dv(ip,jp,k,1)-dv(ip,jm,k,1)
     .             -dv(im,jp,k,1)+dv(im,jm,k,1))/4./dxh(ig)/dyh(jg))

            cnv(2) = dv(i,j,k,2)/dt
cc     .          + di*b0_cnv(i,j,k,3)
     .          + di
     .           *((dv(i,jp,k,1)-2*dv(i,j,k,1)+dv(i,jm,k,1))/dyh(jg)**2
     .            -(dv(ip,jp,k,2)-dv(ip,jm,k,2)
     .             -dv(im,jp,k,2)+dv(im,jm,k,2))/4./dxh(ig)/dyh(jg))

            cnv(3) = dv(i,j,k,3)/dt

            cnv = cnv*vol

            if (si_car) then
              call transformVectorToCartesian
     .              (i,j,k,igrid,igrid,igrid
     .               ,cnv(1),cnv(2),cnv(3),.false.
     .               ,car(1),car(2),car(3))

              do ieq=1,3
                y(neq*(ijk-1)+ieq) = car(ieq)
              enddo
            else
              do ieq=1,3
                y(neq*(ijk-1)+ieq) = cnv(ieq)
              enddo
            endif

          enddo
        enddo
      enddo

c End program

      deallocate(dv,dv_cov)

      nullify(b0_cnv)

      call deallocPointers(fpointers)

      end subroutine hall_mtvc
