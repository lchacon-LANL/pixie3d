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
        order = 1
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

        !This below changes results. Compiler bug?????
cc        allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,1))
        allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,2))

        if (PRESENT(arr0)) then
          v0(:,:,:,1) = arr0(:,:,:,1)
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

        call setBC(ivar,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd
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

      integer(4) :: isig,ijk,ijkg,nnx,nny,nnz
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: drho
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

            if (vol_wgt) then
              dvol = volume(i,j,k,igrid,igrid,igrid)
            else
              dvol = 1d0
            endif

            y(ijk) = ( (1./dt + alpha*mgdivV0(ijkg))*x(ijk)
     .            + alpha*upwind )*dvol
     .            - alpha*dd
     .                *laplacian(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                          ,drho(:,:,:,1),vol=vol_wgt)

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

      integer(4) :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

      real(8)    :: upwind,etal,flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol
     .             ,veclap(3)

      real(8),allocatable,dimension(:,:,:,:) :: db
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      logical    :: fpointers,alt_eom_b

c Begin program

      alt_eom_b = .false.

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

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
     .            ,icomp=IBX,is_cnv=.true.,is_vec=.true.)

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

            if (vol_wgt) then
              vol = volume(i,j,k,igrid,igrid,igrid)
            else
              vol = 1d0
            endif

            veclap = veclaplacian(i,j,k,nxx,nyy,nzz
     .                           ,igrid,igrid,igrid,db
     .                           ,alt_eom_b,vol=.false.)

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

            y(neq*(ijk-1)+1) = ( db(i,j,k,1)/dt + alpha*upwind 
     .                         - alpha*etal*veclap(1))*vol

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

            y(neq*(ijk-1)+2) = ( db(i,j,k,2)/dt + alpha*upwind 
     .                         - alpha*etal*veclap(2) )*vol

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

            y(neq*(ijk-1)+3) = ( db(i,j,k,3)/dt + alpha*upwind
     .                         - alpha*etal*veclap(3) )*vol
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

      contains

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
        logical    :: spoint

c     Begin program

        spoint = isSP(i,j,k,igx,igy,igz)

        igrid = igx

        ip = i+1
        if (flag == 0 .or. isSP(ip,j,k,igx,igy,igz) ) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

cc        if (i + grid_params%ilo(igrid)-1 < grid_params%nxgl(igrid)
cc     .      .and. bcond(1) == SP
cc     .      .and. flag /= 0) then
cc          jacp = gmetric%grid(igrid)%jac(ip,j,k)
cc          jac0 = gmetric%grid(igrid)%jac(i ,j,k)
cc        else
cc          jacp = jac
cc          jac0 = jac
cc        endif
cc
cc        b0cnv(1) = 0.5*(b0_cnv(ip,j,k,1)/jacp+b0_cnv(i,j,k,1)/jac0)*jac
cc        b0cnv(2) = 0.5*(b0_cnv(ip,j,k,2)     +b0_cnv(i,j,k,2))
cc        b0cnv(3) = 0.5*(b0_cnv(ip,j,k,3)/jacp+b0_cnv(i,j,k,3)/jac0)*jac
cc
cc        if (flag /= 0) then
cc          call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
cc     .                        ,b0cov(1),b0cov(2),b0cov(3)
cc     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
cc     .                        ,.false.,half_elem=1)
cc        else
cc          call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
cc     .                        ,b0cov(1),b0cov(2),b0cov(3)
cc     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
cc     .                        ,.false.,half_elem=0)
cc        endif

        if (flag /= 0) then
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),1,igrid)
        else
          call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                      ,acnv(1),acnv(2),acnv(3),0,igrid)
        endif

c diag ****
cc        acnv = 0.5*(dv(i,j,k,:)+dv(ip,j,k,:))
cc        acnv = 0d0
cc        scalar_prod = 0d0
c diag ****

        b0cnv = 0.5*(b0_cnv(ip,j,k,:)+b0_cnv(i,j,k,:))
        b0cov = 0.5*(b0_cov(ip,j,k,:)+b0_cov(i,j,k,:))

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

cc        write (*,*) t11,t12,t13
cc        write (*,*) b0cnv

cc        call find_curl_vxb(i ,j,k,nxx,nyy,nzz,dv,b0_cnv
cc     .                    ,acnv (1),acnv (2),acnv (3),0,igrid)
cc        call find_curl_vxb(ip,j,k,nxx,nyy,nzz,dv,b0_cnv
cc     .                    ,acnvp(1),acnvp(2),acnvp(3),0,igrid)
cc
cc        scalar_prod = 0.5*(dot_product(acnv ,b0_cov(ip,j,k,:))
cc     .                    +dot_product(acnvp,b0_cov(i ,j,k,:)))
cc
cc        t11 =0.5*(acnv (1)*b0_cnv(ip,j,k,1)/jacp/jac0
cc     .           +acnvp(1)*b0_cnv(i ,j,k,1)/jacp/jac0)*jac**2
cc     .      -gsuper(1,1)*scalar_prod
cc
cc        t12 =0.5*(acnv (1)*b0_cnv(ip,j,k,2)/jac0
cc     .           +acnvp(2)*b0_cnv(i ,j,k,1)/jacp)*jac
cc     .      -gsuper(1,2)*scalar_prod
cc
cc        t13 =0.5*(acnv (1)*b0_cnv(ip,j,k,3)/jacp/jac0
cc     .           +acnvp(3)*b0_cnv(i ,j,k,1)/jacp/jac0)*jac**2
cc     .      -gsuper(1,2)*scalar_prod

cc        t11 = dv(i,j,k,1)
cc        t12 = 0d0
cc        t13 = 0d0

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

c diag ****
cc        acnv = 0.5*(dv(i,j,k,:)+dv(i,jp,k,:))
cc        acnv = 0d0
cc        scalar_prod = 0d0
c diag ****

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

cc        t21 = 0d0
cc        t22 = dv(i,j,k,2)
cc        t23 = 0d0

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

c diag ****
cc        acnv = 0.5*(dv(i,j,k,:)+dv(i,j,kp,:))
cc        acnv = 0d0
cc        scalar_prod = 0d0
c diag ****

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

cc        t31 = 0d0
cc        t32 = 0d0
cc        t33 = dv(i,j,k,3)

        if (flag /= 0) then
          t31 = t31/jac
          if (.not.alt_eom) t32 = t32/jac
          t33 = t33/jac
        endif

c     End program

      end subroutine tensor_z

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
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer(4) :: isig,ip,im,jp,jm,kp,km,hex,hey,hez
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax,iming
     .             ,nxx,nyy,nzz,ijk,ijkg,isp

      real(8),allocatable,dimension(:,:)     :: a_sp,a_vol

      real(8)    :: vol,psiv(3),psib(3),psit(3),cov(3),cnv(3),car(3)
     .             ,mag(3),veclap(3)

      logical    :: fpointers,spoint,is_cnv,is_car,postproc

c Begin program

      is_cnv = .true.

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

      nxx = grid_params%nxv(igrid)
      nyy = grid_params%nyv(igrid)
      nzz = grid_params%nzv(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

c Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz,igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dv(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
     .            ,icomp=IVX,is_cnv=is_cnv,is_vec=.not.si_car
     .            ,result_is_vec=.true.)

c Define pointers to MG arrays

      rho0   => grho0  %grid(igrid)%array

      pp0    => gp0    %grid(igrid)%array

      v0_cnv => gv0    %grid(igrid)%array

      b0_cnv => gb0    %grid(igrid)%array

      b0_cov => gb0_cov%grid(igrid)%array

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

            spoint = isSP(i,j,k,igrid,igrid,igrid)

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jacip  = gmetric%grid(igrid)%jac(ip,j,k)
            jacim  = gmetric%grid(igrid)%jac(im,j,k)
            jacjp  = gmetric%grid(igrid)%jac(i,jp,k)
            jacjm  = gmetric%grid(igrid)%jac(i,jm,k)
            jackp  = gmetric%grid(igrid)%jac(i,j,kp)
            jackm  = gmetric%grid(igrid)%jac(i,j,km)
            jac    = gmetric%grid(igrid)%jac(i,j,k )

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
cc            call findPsib_old

c diag ****
cc      psib = 0d0
cc      psit = 0d0
cc      write (*,*) psiv
c diag ****

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

      nabla_v = fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                       ,dv(:,:,:,1),dv(:,:,:,2),dv(:,:,:,3)
     .                       ,hex,hey,hez)

      nabla_vv0 =fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                        ,v0_cnv(:,:,:,1)
     .                        ,v0_cnv(:,:,:,2)
     .                        ,v0_cnv(:,:,:,3)
     .                        ,0,0,0)

      mul = vis(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid)

      veclap = veclaplacian(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv
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
        upwind = upwind
     .            -dt*mgadvdiffV0(ijkg,ieq)*div_upwd(hex)/rho0(i,j,k,1)

        cnv(ieq) = dv(i,j,k,ieq)/dt
     .           + alpha*upwind
     .           - alpha*mul*veclap(ieq)
      enddo

      if (si_car) then
        call transformVectorToCartesian
     .              (i,j,k,igrid,igrid,igrid
     .               ,cnv(1),cnv(2),cnv(3)
     .               ,.false.
     .               ,car(1),car(2),car(3))
        psiv = rho0(i,j,k,1)*car*vol
      elseif (is_cnv) then
        psiv = rho0(i,j,k,1)*cnv*vol
      else
        call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
        psiv = rho0(i,j,k,1)*cov*vol
      endif

      end subroutine findPsiv

c     findPsit
c     #####################################################################
      subroutine findPsit

      implicit none

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
     $             ,divip,divim,divjp,divjm,divkp,divkm

      !Fluxes at faces for calculation of grad(dv.grad(p0))
      flxip =( (dv(i,j,k,1)+dv(ip,j,k,1))
     .           *(pp0(ip,j,k,1)-pp0(i,j,k,1))/dx(ig)
     .        +(dv(i,j,k,2)+dv(ip,j,k,2))
     .           *(pp0(ip,jp,k,1)-pp0(ip,jm,k,1)
     .            +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
     .        +(dv(i,j,k,3)+dv(ip,j,k,3))
     .           *(pp0(ip,j,kp,1)-pp0(ip,j,km,1)
     .            +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4. )
     $        /(jac+jacip)
      flxim =( (dv(i,j,k,1)+dv(im,j,k,1))
     .            *(pp0(i ,j,k,1)-pp0(im,j,k,1))/dx(ig-1)
     .        +(dv(i,j,k,2)+dv(im,j,k,2))
     .            *(pp0(im,jp,k,1)-pp0(im,jm,k,1)
     .             +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/dyh(jg)/4.
     .        +(dv(i,j,k,3)+dv(im,j,k,3))
     .            *(pp0(im,j,kp,1)-pp0(im,j,km,1)
     .             +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/dzh(kg)/4.)
     $        /(jac+jacim)

      flxjp =( (dv(i,j,k,1)+dv(i,jp,k,1))
     .            *(pp0(ip,jp,k,1)-pp0(im,jp,k,1)
     .             +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/dxh(ig)/4.
     .        +(dv(i,j,k,2)+dv(i,jp,k,2))
     .            *(pp0(i,jp,k,1)-pp0(i,j,k,1))/dy(jg)
     .        +(dv(i,j,k,3)+dv(i,jp,k,3))
     .            *(pp0(i,jp,kp,1)-pp0(i,jp,km,1)
     .             +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacjp)
      flxjm =( (dv(i,j,k,1)+dv(i,jm,k,1))
     .            *(pp0(ip,jm,k,1)-pp0(im,jm,k,1)
     .             +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/dxh(ig)/4.
     .        +(dv(i,j,k,2)+dv(i,jm,k,2))
     .            *(pp0(i,j ,k,1)-pp0(i,jm,k,1))/dy(jg-1)
     .        +(dv(i,j,k,3)+dv(i,jm,k,3))
     .            *(pp0(i,jm,kp,1)-pp0(i,jm,km,1)
     .             +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacjm)

      flxkp =( (dv(i,j,k,1)+dv(i,j,kp,1))
     .            *(pp0(ip,j,kp,1)-pp0(im,j,kp,1)
     .             +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/dxh(ig)/4.
     .        +(dv(i,j,k,2)+dv(i,j,kp,2))
     .            *(pp0(i,jp,kp,1)-pp0(i,jm,kp,1)
     .             +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/dyh(jg)/4.
     .        +(dv(i,j,k,3)+dv(i,j,kp,3))
     .            *(pp0(i,j,kp,1)-pp0(i,j,k,1))/dz(kg) )
     $        /(jac+jackp)
      flxkm =( (dv(i,j,k,1)+dv(i,j,km,1))
     .            *(pp0(ip,j,km,1)-pp0(im,j,km,1)
     .             +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/dxh(ig)/4.
     .        +(dv(i,j,k,2)+dv(i,j,km,2))
     .            *(pp0(i,jp,km,1)-pp0(i,jm,km,1)
     .             +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/dyh(jg)/4.
     .        +(dv(i,j,k,3)+dv(i,j,km,3))
     .            *(pp0(i,j,k ,1)-pp0(i,j,km,1))/dz(kg-1) )
     $        /(jac+jackm)

      !Fluxes at faces for calculation of grad(gamma*p0*div(dv))

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
     .      + gamma*(pp0(i,j,k,1)+pp0(ip,j,k,1))*divip/(jac+jacip)
      flxim = flxim
     .      + gamma*(pp0(i,j,k,1)+pp0(im,j,k,1))*divim/(jac+jacim)

      flxjp = flxjp
     .      + gamma*(pp0(i,j,k,1)+pp0(i,jp,k,1))*divjp/(jac+jacjp)
      flxjm = flxjm
     .      + gamma*(pp0(i,j,k,1)+pp0(i,jm,k,1))*divjm/(jac+jacjm)

      flxkp = flxkp
     .      + gamma*(pp0(i,j,k,1)+pp0(i,j,kp,1))*divkp/(jac+jackp)
      flxkm = flxkm
     .      + gamma*(pp0(i,j,k,1)+pp0(i,j,km,1))*divkm/(jac+jackm)

      cov(1) = (flxip - flxim)/dxh(ig)
      cov(2) = (flxjp - flxjm)/dyh(jg)
      cov(3) = (flxkp - flxkm)/dzh(kg)

      if (si_car) then
        call transformVectorToCartesian
     .              (i,j,k,igrid,igrid,igrid
     .               ,cov(1),cov(2),cov(3)
     .               ,.true.
     .               ,car(1),car(2),car(3))
        psit = -dt*alpha**2*car*vol
      elseif (is_cnv) then
        call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
        psit = -dt*alpha**2*cnv*vol
      else
        psit = -dt*alpha**2*cov*vol
      endif

      end subroutine findPsit

c     findPsib_old
c     #####################################################################
      subroutine findPsib_old

      implicit none

      real(8)    :: idx,idy,idz,coeff

      real(8)    :: a10,a20,a30
     .             ,a1covip,a2covip,a3covip,a1covim,a2covim,a3covim
     .             ,a1covjp,a2covjp,a3covjp,a1covjm,a2covjm,a3covjm
     .             ,a1covkp,a2covkp,a3covkp,a1covkm,a2covkm,a3covkm
     .             ,a1cnvip,a2cnvip,a3cnvip,a1cnvim,a2cnvim,a3cnvim
     .             ,a1cnvjp,a2cnvjp,a3cnvjp,a1cnvjm,a2cnvjm,a3cnvjm
     .             ,a1cnvkp,a2cnvkp,a3cnvkp,a1cnvkm,a2cnvkm,a3cnvkm

cc      real(8) :: j0_cnv(3)

      !(j0 x a) part; a = curl(dv x B0); UPWIND for diagonal dominance
cc      hex = 1
cc      hey = 1
cc      hez = 1
cc      if (mgj0cnv(ijkg,1) > 0d0) hey = -1
cc      if (mgj0cnv(ijkg,2) > 0d0) hez = -1
cc      if (mgj0cnv(ijkg,3) > 0d0) hex = -1
cc      call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
cc     .                  ,a10,a20,a30,hex,hey,hez,igrid)
cc
cc      cov(1) = -mgj0cnv(ijkg,2)*a30 
cc      cov(2) = -mgj0cnv(ijkg,3)*a10 
cc      cov(3) = -mgj0cnv(ijkg,1)*a20 
cc
cc      hex = 1
cc      hey = 1
cc      hez = 1
cc      if (mgj0cnv(ijkg,1) > 0d0) hez = -1
cc      if (mgj0cnv(ijkg,2) > 0d0) hex = -1
cc      if (mgj0cnv(ijkg,3) > 0d0) hey = -1
cc      hex = 0
cc      hey = 0
cc      hez = 0
cc
cc      call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
cc     .                  ,a10,a20,a30,hex,hey,hez,igrid)
cc
cc      cov(1) = cov(1) + mgj0cnv(ijkg,3)*a20
cc      cov(2) = cov(2) + mgj0cnv(ijkg,1)*a30
cc      cov(3) = cov(3) + mgj0cnv(ijkg,2)*a10

      call find_curl_vxb(i,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                  ,a10,a20,a30,0,igrid)

      cov(1) = mgj0cnv(ijkg,2)*a30 - mgj0cnv(ijkg,3)*a20
      cov(2) = mgj0cnv(ijkg,3)*a10 - mgj0cnv(ijkg,1)*a30
      cov(3) = mgj0cnv(ijkg,1)*a20 - mgj0cnv(ijkg,2)*a10

c diag ****
cc      cov = 0d0
c diag ****

      !B0 x curl(a) part

      im = i-1
      jm = j-1
      km = k-1

      !!Find contravariant components of curl(dv x B0) at faces

      call find_curl_vxb(i ,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                  ,a1cnvip,a2cnvip,a3cnvip,1,igrid)

cc      if (spoint .and. (.not.form_diag) ) then
cccc        write (*,*) 'here'
cc        call transformVectorToCurvilinear(im,j,k,igrid,igrid,igrid
cc     .               ,a_sp(k,1),a_sp(k,2),a_sp(k,3),.true.
cc     .               ,a1covim,a2covim,a3covim)
cc        coeff = 1d0
cc      elseif (spoint .and. form_diag) then

cc      if (spoint) then
cc        call find_curl_vxb(i ,j,k,nxx,nyy,nzz,dv,b0_cnv
cc     .                    ,a1cnvim,a2cnvim,a3cnvim,0,igrid)
cc        coeff = 2d0
cc      else
        call find_curl_vxb(im,j,k,nxx,nyy,nzz,dv,b0_cnv
     .                    ,a1cnvim,a2cnvim,a3cnvim,1,igrid)
        coeff = 1d0
cc      endif

      if (spoint) then
        a1cnvim = 0d0
        a3cnvim = 0d0
      endif

      call find_curl_vxb(i,j ,k,nxx,nyy,nzz,dv,b0_cnv
     .                   ,a1cnvjp,a2cnvjp,a3cnvjp,2,igrid)
      call find_curl_vxb(i,jm,k,nxx,nyy,nzz,dv,b0_cnv
     .                   ,a1cnvjm,a2cnvjm,a3cnvjm,2,igrid)

      call find_curl_vxb(i,j,k ,nxx,nyy,nzz,dv,b0_cnv
     .                   ,a1cnvkp,a2cnvkp,a3cnvkp,3,igrid)
      call find_curl_vxb(i,j,km,nxx,nyy,nzz,dv,b0_cnv
     .                   ,a1cnvkm,a2cnvkm,a3cnvkm,3,igrid)

      !!Find covariant components of curl(dv x B0) at faces

      call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                        ,a1covip,a2covip,a3covip
     .                        ,a1cnvip,a2cnvip,a3cnvip
     .                        ,.false.,half_elem=1)

cc      if (spoint) then
cc        call transformFromCurvToCurv(i ,j,k,igrid,igrid,igrid
cc     .                        ,a1covim,a2covim,a3covim
cc     .                        ,a1cnvim,a2cnvim,a3cnvim
cc     .                        ,.false.,half_elem=0)
cc      else
        call transformFromCurvToCurv(im,j,k,igrid,igrid,igrid
     .                        ,a1covim,a2covim,a3covim
     .                        ,a1cnvim,a2cnvim,a3cnvim
     .                        ,.false.,half_elem=1)
cc      endif

      call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                        ,a1covjp,a2covjp,a3covjp
     .                        ,a1cnvjp,a2cnvjp,a3cnvjp
     .                        ,.false.,half_elem=2)
      call transformFromCurvToCurv(i,jm,k,igrid,igrid,igrid
     .                        ,a1covjm,a2covjm,a3covjm
     .                        ,a1cnvjm,a2cnvjm,a3cnvjm
     .                        ,.false.,half_elem=2)

      call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                        ,a1covkp,a2covkp,a3covkp
     .                        ,a1cnvkp,a2cnvkp,a3cnvkp
     .                        ,.false.,half_elem=3)
      call transformFromCurvToCurv(i,j,km,igrid,igrid,igrid
     .                        ,a1covkm,a2covkm,a3covkm
     .                        ,a1cnvkm,a2cnvkm,a3cnvkm
     .                        ,.false.,half_elem=3)

      if (spoint) then
        a2covim = 0d0
      endif

      !!Assemble B0 x curl(a) and j0 x a
      idx = coeff/dxh(ig)
      idy = 1d0  /dyh(jg)
      idz = 1d0  /dzh(kg)

      cov(1) =-cov(1)
     .        +b0_cnv(i,j,k,2)*( (a2covip-a2covim)*idx
     .                          -(a1covjp-a1covjm)*idy)
     .        +b0_cnv(i,j,k,3)*( (a3covip-a3covim)*idx
     .                          -(a1covkp-a1covkm)*idz)

      cov(2) =-cov(2)
     .        +b0_cnv(i,j,k,1)*( (a1covjp-a1covjm)*idy
     .                          -(a2covip-a2covim)*idx)
     .        +b0_cnv(i,j,k,3)*( (a3covjp-a3covjm)*idy
     .                          -(a2covkp-a2covkm)*idz)

      cov(3) =-cov(3)
     .        +b0_cnv(i,j,k,1)*( (a1covkp-a1covkm)*idz
     .                          -(a3covip-a3covim)*idx)
     .        +b0_cnv(i,j,k,2)*( (a2covkp-a2covkm)*idz
     .                          -(a3covjp-a3covjm)*idy)
c diag *****
cc      cov = -cov
c diag *****

      !Introduce jacobian factor
      cov = cov/jac

      if (si_car) then
        call transformVectorToCartesian(i,j,k,igrid,igrid,igrid
     .                                 ,cov(1),cov(2),cov(3)
     .                                 ,.true.
     .                                 ,car(1),car(2),car(3))
        psib = -dt*alpha**2*car*vol
      elseif (is_cnv) then
        call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
        psib = -dt*alpha**2*cnv*vol
      else
        psib = -dt*alpha**2*cov*vol
      endif

      end subroutine findPsib_old

c     findPsib
c     #####################################################################
      subroutine findPsib

      implicit none

cc      cnv = divtensor(i,j,k,alt_eom)

cc      cnv = dv(i,j,k,:)

      cnv = div_tensor(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,alt_eom
     .                ,tensor_x,tensor_y,tensor_z,vol=.false.)

      if (si_car) then
        call transformVectorToCartesian(i,j,k,igrid,igrid,igrid
     .                                 ,cnv(1),cnv(2),cnv(3),.false.
     .                                 ,car(1),car(2),car(3))
        psib = dt*alpha**2*car*vol
      elseif (is_cnv) then
        psib = dt*alpha**2*cnv*vol
      else
        call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
        psib = dt*alpha**2*cov*vol
      endif

      end subroutine findPsib

ccc     divtensor
ccc     ###############################################################
cc      function divtensor(i,j,k,alt_eom) result (divt)
cc
ccc     ---------------------------------------------------------------
ccc     Calculates dvol*lap(vector) at cell centers in general non-orthog.
ccc     coordinates.
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      real(8)    :: divt(3)
cc
cc      integer(4) :: i,j,k
cc
cc      logical    :: alt_eom
cc
ccc     Local variables
cc
cc      integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km
cc
cc      real(8)    :: jac,dvol,dS1,dS2,dS3,dxx,dyy,dzz
cc
cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc
cc      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
cc     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
cc     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o
cc
cc      real(8)    :: hess(3,3,3),msource
cc
ccc     Begin program
cc
cc      call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
cc
cc      jac = gmetric%grid(igrid)%jac(i,j,k)
cc
cc      dxx = dxh(ig)
cc      dyy = dyh(jg)
cc      dzz = dzh(kg)
cc
cc      dS1 = dyy*dzz
cc      dS2 = dxx*dzz
cc      dS3 = dxx*dyy
cc
cc      dvol = dxx*dyy*dzz
cc
cc      ip = i+1
cc      im = i-1
cc      jp = j+1
cc      jm = j-1
cc      kp = k+1
cc      km = k-1
cc
cc      if (coords /= 'car') then
cc        hess = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc      endif
cc
cc      call tensor_x(i ,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .             ,t11p,t12p,t13p,1)
cc      call tensor_x(im,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .             ,t11m,t12m,t13m,-1)
cc      if (coords /= 'car') call tensor_x(i,j,k,nxx,nyy,nzz,igx,igy,igz
cc     .                                  ,alt_eom,t11o,t12o,t13o,0)
cc
cc      call tensor_y(i,j ,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .             ,t21p,t22p,t23p,1)
cc      call tensor_y(i,jm,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .             ,t21m,t22m,t23m,-1)
cc      if (coords /= 'car') call tensor_y(i,j,k,nxx,nyy,nzz,igx,igy,igz
cc     .                                  ,alt_eom,t21o,t22o,t23o,0)
cc
cc      call tensor_z(i,j,k ,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .             ,t31p,t32p,t33p,1)
cc      call tensor_z(i,j,km,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .             ,t31m,t32m,t33m,-1)
cc      if (coords /= 'car') call tensor_z(i,j,k,nxx,nyy,nzz,igx,igy,igz
cc     .                                  ,alt_eom,t31o,t32o,t33o,0)
cc
cc      msource = 0d0
cc
cc      !Component 1
cc      flxip = t11p
cc      flxim = t11m
cc
cc      flxjp = t21p
cc      flxjm = t21m
cc
cc      flxkp = t31p
cc      flxkm = t31m
cc
cc      if (spoint) flxim = 0d0
cc
cc      if (coords /= 'car') then
cc        msource =  (t11o*hess(1,1,1)+t12o*hess(1,1,2)+t13o*hess(1,1,3)
cc     .             +t21o*hess(1,2,1)+t22o*hess(1,2,2)+t23o*hess(1,2,3)
cc     .             +t31o*hess(1,3,1)+t32o*hess(1,3,2)+t33o*hess(1,3,3))
cc      endif
cc
cc      divt(1) =  (flxip - flxim)/dxx
cc     .          +(flxjp - flxjm)/dyy
cc     .          +(flxkp - flxkm)/dzz + msource/jac
cc
cc      !Component 2
cc      flxip = t12p
cc      flxim = t12m
cc
cc      flxjp = t22p
cc      flxjm = t22m
cc
cc      flxkp = t32p
cc      flxkm = t32m
cc
cc      if (coords /= 'car') then
cc        if (alt_eom) then
cc
cc          if (spoint) flxim = 0d0
cc
cc          msource=  (t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
cc     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
cc     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3)
cc     .              -t12o*(hess(1,1,1)+hess(2,1,2)+hess(3,1,3))
cc     .              -t22o*(hess(1,2,1)+hess(2,2,2)+hess(3,2,3))
cc     .              -t32o*(hess(1,3,1)+hess(2,3,2)+hess(3,3,3)))
cc
cc          flxip = flxip/jac
cc          flxim = flxim/jac
cc                         
cc          flxjp = flxjp/jac
cc          flxjm = flxjm/jac
cc                         
cc          flxkp = flxkp/jac
cc          flxkm = flxkm/jac
cc
cc        else
cc          msource=  (t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
cc     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
cc     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3))
cc        endif
cc      endif
cc
cc      divt(2) =  (flxip - flxim)/dxx
cc     .          +(flxjp - flxjm)/dyy
cc     .          +(flxkp - flxkm)/dzz  + msource/jac
cc
cc      !Component 3
cc      flxip = t13p
cc      flxim = t13m
cc
cc      flxjp = t23p
cc      flxjm = t23m
cc
cc      flxkp = t33p
cc      flxkm = t33m
cc
cc      if (spoint) flxim = 0d0
cc
cc      if (coords /= 'car') then
cc        msource = (t11o*hess(3,1,1)+t12o*hess(3,1,2)+t13o*hess(3,1,3)
cc     .            +t21o*hess(3,2,1)+t22o*hess(3,2,2)+t23o*hess(3,2,3)
cc     .            +t31o*hess(3,3,1)+t32o*hess(3,3,2)+t33o*hess(3,3,3))
cc      endif
cc
cc      divt(3) = ( (flxip - flxim)/dxx
cc     .           +(flxjp - flxjm)/dyy
cc     .           +(flxkp - flxkm)/dzz ) + msource/jac
cc
ccc     End program
cc
cc      end function divtensor

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
