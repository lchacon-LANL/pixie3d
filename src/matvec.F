
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

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer    :: i,j,k,ig,jg,kg,isig,ijk,ijkg,nxx,nyy,nzz
      integer    :: imin,imax,jmin,jmax,kmin,kmax,order

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

      order = bc_order
      if (    abs(bcnd(1,1)) == NEU
     .    .or.abs(bcnd(2,1)) == NEU
     .    .or.abs(bcnd(3,1)) == NEU
     .    .or.abs(bcnd(4,1)) == NEU
     .    .or.abs(bcnd(5,1)) == NEU
     .    .or.abs(bcnd(6,1)) == NEU) order = 0

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,xarr,bcnd
     .            ,iorder=order,icomp=(/IRHO/))

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

c lap_mtvc2
c####################################################################
      subroutine lap_mtvc2(gpos,neq,ntot,x,y,igrid,bcnd)
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

      use equilibrium

      implicit none

c Call variables

      integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8) :: x(ntot),y(ntot)

c Local variables

      integer :: i,j,k,ig,jg,kg,isig,ijk,ijkg,nxx,nyy,nzz
      integer :: imin,imax,jmin,jmax,kmin,kmax

      real(8),pointer,dimension(:,:,:,:) :: xarr,bz0

      logical :: fpointers

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

      bz0 => gbz0%grid(igrid)%array

      allocate(xarr(0:nxx+1,0:nyy+1,0:nzz+1,neq))

cc      xarr = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,xarr,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,xarr,bcnd
     .            ,iorder=bc_order,arr0=bz0)

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

      end subroutine lap_mtvc2

c del_hel
c####################################################################
      subroutine del_hel(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine is a matvec for the helical Grad-Shafranov
c     operator, del^H. In call:
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

      use equilibrium

      implicit none

c Call variables

      integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8) :: x(ntot),y(ntot)

c Local variables

      integer :: ip,im,jp,jm,kp,km,isig,ijk,ijkg,nxx,nyy,nzz
      integer :: imin,imax,jmin,jmax,kmin,kmax,i,j,k,ig,jg,kg
      real(8) :: xhp,xhm,kk,x1,y1,z1

      real(8),pointer,dimension(:,:,:,:) :: xarr,psi0

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

      psi0 => gpsi0%grid(igrid)%array

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,xarr,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,xarr,bcnd
     .            ,iorder=bc_order,arr0=psi0)

c Calculate matrix-vector product

      kk = grid_params%params(2)

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

cc            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
            call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
     .                                    ,ig,jg,kg,x1,y1,z1)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ip = i+1
            im = i-1
            jp = j+1
            jm = j-1
            kp = k+1
            km = k-1

            xhp = 0.5*(grid_params%xx(ig)+grid_params%xx(ig+1))
            xhm = 0.5*(grid_params%xx(ig)+grid_params%xx(ig-1))
      
            y(ijk) = 1d0/x1/grid_params%dxh(ig)
     .                *( xhp/(1+(kk*xhp)**2)
     .                      *(xarr(ip,j,k,1)-xarr(i,j,k,1))
     .                      /grid_params%dx(ig)
     .                  +xhm/(1+(kk*xhm)**2)
     .                      *(xarr(im,j,k,1)-xarr(i,j,k,1))
     .                      /grid_params%dx(ig-1))
     .              +1d0/x1**2/grid_params%dyh(jg)
     .                *( (xarr(i,jp,k,1)-xarr(i,j,k,1))
     .                   /grid_params%dy(jg)
     .                  +(xarr(i,jm,k,1)-xarr(i,j,k,1))
     .                   /grid_params%dy(jg-1) )

            if (vol_wgt) y(ijk)=y(ijk)*gmetric%grid(igrid)%dvol(i,j,k)
          enddo
        enddo
      enddo

c End program

      deallocate(xarr)

      call deallocPointers(fpointers)

      end subroutine del_hel

c del_star
c####################################################################
      subroutine del_star(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine is a matvec for the toroidal Grad-Shafranov
c     operator, del^*. In call:
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

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer    :: ip,im,jp,jm,kp,km,isig,ijk,ijkg,nxx,nyy,nzz
      integer    :: imin,imax,jmin,jmax,kmin,kmax,i,j,k,ig,jg,kg

      real(8),allocatable,dimension(:,:,:,:) :: xarr
cc      real(8),pointer    ,dimension(:,:,:)   :: ljac
      real(8),allocatable,dimension(:,:,:)   :: iRR

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

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,xarr,bcnd
     .            ,iorder=bc_order)

c Calculate matrix-vector product

      !Initalize toroidal radius (R=R0+r*cos(theta)=(x^2+y^2))
      allocate(iRR(0:nxx+1,0:nyy+1,0:nzz+1))

      do k = kmin-1,kmax+1
        do j = jmin-1,jmax+1
          do i = imin-1,imax+1
cc            call getCartesianCoordinates(i,j,k,igrid,igrid,igrid
cc     .                                  ,ig,jg,kg,x1,y1,z1)
cc
cc            iRR(i,j,k) = 1d0/(x1**2+y1**2)
            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            iRR(i,j,k) =  1d0/(grid_params%params(1)
     .                       + grid_params%xx(ig)
     .                        *cos(grid_params%yy(jg)))**2
          enddo
        enddo
      enddo

      !Main loop
      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            y(ijk) = laplacian(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                        ,xarr(:,:,:,1),dff=iRR,vol=vol_wgt)
          enddo
        enddo
      enddo

c End program

      deallocate(xarr,iRR)

      call deallocPointers(fpointers)

      end subroutine del_star

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

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer    :: isig,ijk,ijkg,nxx,nyy,nzz,i,j,k,ig,jg,kg
      integer    ::  imin, imax, jmin, jmax, kmin, kmax

      real(8),allocatable,dimension(:,:,:,:) :: xarr
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,nu2,jac

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
     .            ,icomp=(/IRHO/))

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

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer    :: isig,ijk,ijkg,nnx,nny,nnz,i,j,k,ig,jg,kg
      integer    :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: dtmp
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv,rho0

      real(8)    :: upwind,nu2,dvol,lap,gm,jac,sgn

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
     .            ,icomp=(/ITMP/),iorder=bc_order)

c Map velocity components

      v0_cnv => gvs0 %grid(igrid)%array
      rho0   => grho0%grid(igrid)%array

c Calculate matrix-vector product

      if (solve_prs) then
        gm = gamma
        sgn =-1d0
      else
        gm = gamma - 1d0
        sgn = 1d0
      endif

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jac = gmetric%grid(igrid)%jac(i,j,k)

            ijk    = i + nnx*(j-1) + nnx*nny*(k-1)
            ijkg   = ijk + isig - 1

cc            upwind = gm*c_advec(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
cc     .                    ,v0_cnv,dtmp  (:,:,:,1)
cc     .                    ,vol=vol_wgt,upwind=.true.,sp=bcSP())
cc     .              +(1.-gm)*dot_product(v0_cnv(i,j,k,:)
cc     .                  ,grad(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
cc     .                       ,dtmp(:,:,:,1)
cc     .                       ,hex=int(sign(1d0,-sgn*v0_cnv(i,j,k,1)))
cc     .                       ,hey=int(sign(1d0,-sgn*v0_cnv(i,j,k,2)))
cc     .                       ,hez=int(sign(1d0,-sgn*v0_cnv(i,j,k,3)))))
cc     .                  /jac
            upwind = gm*flx_advec(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                           ,v0_cnv,dtmp(:,:,:,1),1
     .                           ,vol=vol_wgt,zip_vel=.false.)!Do not use ZIP, since dT is not pos.def.
     .             -sgn*(gm-1.)
     .                 *flx_nc_advec(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                              ,sgn*v0_cnv,dtmp(:,:,:,1),1
     .                              ,vol=vol_wgt)

            if (vol_wgt) then
              dvol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            if (chi > 0d0) then
              lap = laplacian(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                          ,dtmp(:,:,:,1),vol=.false.)
              lap = (gamma-1d0)*chi*lap
            else
              lap = 0d0
            endif

            if (.not.solve_prs) lap=lap/rho0(i,j,k,1)

            y(ijk) =dvol*( cnp/dt*x(ijk) + alpha*(upwind - lap) )

          enddo
        enddo
      enddo

c End program

      deallocate(dtmp)
      nullify(v0_cnv,rho0)

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

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer    :: isig,ijk,ijkg,nnx,nny,nnz,ip,im,jp,jm,kp,km
      integer    :: imin,imax,jmin,jmax,kmin,kmax,i,j,k,ig,jg,kg
 
      real(8),allocatable,dimension(:,:,:,:) :: drho
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,lap,dvol,flxip,flxim,flxjp,flxjm,flxkp,flxkm
     .             ,jac

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
     .            ,icomp=(/IRHO/),iorder=bc_order)

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

            !Conservative advection
            upwind = flx_advec(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                        ,v0_cnv,drho(:,:,:,1),1,vol=vol_wgt
     .                        ,zip_vel=.false.) !Do not use ZIP form because drho is not pos. def.

            ! Diffusion
            if (dd > 0d0) then
              lap = laplacian(i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                       ,drho(:,:,:,1),vol=vol_wgt)
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
cc      implicit none
cc
ccc Call variables
cc
cc      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
cc      real(8)    :: x(ntot),y(ntot)
cc
ccc Local variables
cc
cc      integer    :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
cc      integer    :: imin,imax,jmin,jmax,kmin,kmax
cc      integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
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

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer    :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer    :: imin,imax,jmin,jmax,kmin,kmax,igx,igy,igz
      integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
     .             ,i,j,k,ig,jg,kg,ieq

      real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm,jac
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol
     .             ,veclap(3),cnv(3),cov(3),car(3)
      real(8)    :: upwind,etal

      real(8),allocatable,dimension(:,:,:,:) :: db
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv,b0_cnv

      logical    :: fpointers,is_cnv

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
     .            ,icomp=(/IBX/),is_cnv=is_cnv,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=bc_order)

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
              vol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              vol = 1d0
            endif

            etal   = geta%grid(igrid)%array(i,j,k,1)

            veclap = etal*veclaplacian(i,j,k,nxx,nyy,nzz
     .                           ,igrid,igrid,igrid,db
     .                           ,vol=.false.)
cc            veclap = -etal*curlcurl(i,j,k,nxx,nyy,nzz
cc     .                           ,igrid,igrid,igrid,db
cc     .                           ,vol=.false.)

            if (sm_type == 'gm') then

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

              cnv(1) =  (flxjp-flxjm)/grid_params%dyh(jg)
     .                 +(flxkp-flxkm)/grid_params%dzh(kg)

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
                flxim = 0.5/(jac+jacim)*(
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

              cnv(2) =  (flxip-flxim)/grid_params%dxh(ig)
     .                 +(flxkp-flxkm)/grid_params%dzh(kg)

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

              cnv(3) =  (flxip-flxim)/grid_params%dxh(ig)
     .                 +(flxjp-flxjm)/grid_params%dyh(jg)

            endif

            cnv=(cnp*db(i,j,k,:)/dt+alpha*(cnv-veclap))*vol

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

c b_mtvc_hyp
c####################################################################
      subroutine b_mtvc_hyp(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for Faraday's law, with a hyper-resistive term.
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

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer    :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer    :: imin,imax,jmin,jmax,kmin,kmax,igx,igy,igz
      integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
     .             ,i,j,k,ig,jg,kg,ieq

      real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm,jac
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol
     .             ,veclap_b(3),veclap_b2(3),cnv(3),cov(3),car(3)
      real(8)    :: upwind,etal

      real(8),allocatable,dimension(:,:,:,:) :: db
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      logical    :: fpointers,is_cnv

c Begin program

      if (neq /= 6) call pstop('b_mtvc_hyp','neq /= 6; aborting...')
      
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

      call limits(abs(gpos),nxx,nyy,nzz
     .           ,igrid,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(db(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      db = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,db,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,db
     .            ,bcnd,icomp=(/IBX,IJX/),is_cnv=is_cnv  !Check icomp!!!
cc     .            ,is_vec=.not.si_car,result_is_vec=.true.
     .            ,is_vec=.true.,result_is_vec=.true.
     .            ,iorder=bc_order)

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

            veclap_b = veclaplacian(i,j,k,nxx,nyy,nzz
     .                             ,igrid,igrid,igrid,db(:,:,:,1:3)
     .                             ,vol=.false.)

            veclap_b2= veclaplacian(i,j,k,nxx,nyy,nzz
     .                             ,igrid,igrid,igrid,db(:,:,:,4:6)
     .                             ,vol=.false.)

            etal   = geta%grid(igrid)%array(i,j,k,1)

            if (sm_type == 'gm') then

              cnv = curl_bxv(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                      ,v0_cnv,db(:,:,:,1:3),0)

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

              cnv(1) =  (flxjp-flxjm)/grid_params%dyh(jg)
     .                 +(flxkp-flxkm)/grid_params%dzh(kg)

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

              cnv(2) =  (flxip-flxim)/grid_params%dxh(ig)
     .                 +(flxkp-flxkm)/grid_params%dzh(kg)

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

              cnv(3) =  (flxip-flxim)/grid_params%dxh(ig)
     .                 +(flxjp-flxjm)/grid_params%dyh(jg)

            endif

cc            if (si_car) then
cc              cnv=(cnp*db(i,j,k,1:3)/dt + alpha*cnv
cc     .            - alpha*etal*veclap_b + alpha*di*heta*veclap_b2)*vol
cc
cc              call transformVectorToCartesian
cc     .              (i,j,k,igrid,igrid,igrid
cc     .               ,cnv(1),cnv(2),cnv(3)
cc     .               ,.false.
cc     .               ,car(1),car(2),car(3))
cc
cc              do ieq=1,3
cc                y(neq*(ijk-1)+ieq) = car(ieq)
cc              enddo
cc
cc              cnv = (db(i,j,k,4:6) + veclap_b)*vol
cc
cc              call transformVectorToCartesian
cc     .              (i,j,k,igrid,igrid,igrid
cc     .               ,cnv(1),cnv(2),cnv(3)
cc     .               ,.false.
cc     .               ,car(1),car(2),car(3))
cc
cc              do ieq=4,6
cc                y(neq*(ijk-1)+ieq) = car(ieq)
cc              enddo
cc
cc            else
              do ieq=1,3
                y(neq*(ijk-1)+ieq) = (cnp*db(i,j,k,ieq)/dt
     .                               + alpha*cnv(ieq)
     .                               - alpha*etal*veclap_b(ieq)
     .                               + alpha*di*heta*veclap_b2(ieq))*vol
              enddo

              do ieq=4,6
                y(neq*(ijk-1)+ieq) = (db(i,j,k,ieq)-veclap_b(ieq-3))*vol
              enddo

cc            endif

          enddo
        enddo
      enddo

c End program

      deallocate(db)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine b_mtvc_hyp

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

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer    :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer    :: imin,imax,jmin,jmax,kmin,kmax,igx,igy,igz
      integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
     .             ,i,j,k,ig,jg,kg,ieq

      real(8)    :: upwind,etal,flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol
     .             ,veclap(3),cnv(3),cov(3),car(3),jac

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
     .            ,icomp=(/IAX/),is_cnv=is_cnv,is_vec=.true.
cc     .            ,icomp=(/IAX/),is_cnv=is_cnv,is_vec=.not.si_car
     .            ,result_is_vec=.true.,iorder=bc_order)

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
cc              vol = volume(i,j,k,igrid,igrid,igrid)
              vol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              vol = 1d0
            endif

            !Resistive part
            veclap = veclap_cov(i,j,k,nxx,nyy,nzz
     .                         ,igrid,igrid,igrid,da
     .                         ,vol=.false.)

cc            cnv = veclaplacian(i,j,k,nxx,nyy,nzz
cc     .                           ,igrid,igrid,igrid,da_cnv
cc     .                           ,vol=.false.)
cc
cc            call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,veclap(1),veclap(2),veclap(3)
cc     .                        ,cnv   (1),cnv   (2),cnv   (3)
cc     .                        ,.false.)

            etal   = geta%grid(igrid)%array(i,j,k,1)

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

            cov(1) =  (flxip-flxim)/grid_params%dxh(ig)
     .               +(flxjp-flxjm)/grid_params%dyh(jg)
     .               +(flxkp-flxkm)/grid_params%dzh(kg)


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

            cov(2) =  (flxip-flxim)/grid_params%dxh(ig)
     .               +(flxjp-flxjm)/grid_params%dyh(jg)
     .               +(flxkp-flxkm)/grid_params%dzh(kg)


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

            cov(3) =  (flxip-flxim)/grid_params%dxh(ig)
     .               +(flxjp-flxjm)/grid_params%dyh(jg)
     .               +(flxkp-flxkm)/grid_params%dzh(kg)

            cov=(cnp*da(i,j,k,:)/dt + alpha*(cov/jac - etal*veclap))*vol

cc            if (si_car) then
cc              call transformVectorToCartesian
cc     .              (i,j,k,igrid,igrid,igrid
cc     .               ,cov(1),cov(2),cov(3)
cc     .               ,.true.
cc     .               ,car(1),car(2),car(3))
cc
cc              do ieq=1,3
cc                y(neq*(ijk-1)+ieq) = car(ieq)
cc              enddo
cc            else
              do ieq=1,3
                y(neq*(ijk-1)+ieq) = cov(ieq)
              enddo
cc            endif

          enddo
        enddo
      enddo

c End program

      deallocate(da)
cc      deallocate(da_cnv)
      nullify(v0_cnv)

      call deallocPointers(fpointers)

      end subroutine a_mtvc

c a_mtvc_hyp
c####################################################################
      subroutine a_mtvc_hyp(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for Faraday's law, with a hyper-resistive term.
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

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c Local variables

      integer    :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer    :: imin,imax,jmin,jmax,kmin,kmax,igx,igy,igz
      integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
     .             ,i,j,k,ig,jg,kg,ieq

      real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm,jac
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol
     .             ,a_res(3),curlcurla(3),cnv(3),cov(3),car(3)
      real(8)    :: upwind,etal

      real(8),pointer,dimension(:,:,:,:) :: du,du_cov
      real(8),allocatable,dimension(:,:,:,:) :: das,db_cov
      real(8),pointer    ,dimension(:,:,:,:) :: ve0_cnv,da_cov,dj_cnv
     .                                        ,rho0,bs0_cnv,da_cnv,heta0

      logical    :: fpointers,is_cnv

c Begin program

      if (neq /= 6) call pstop('a_mtvc_hyp','neq /= 6; aborting...')
      
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

      allocate(das   (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,db_cov(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,da_cov(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,da_cnv(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dj_cnv(0:nxx+1,0:nyy+1,0:nzz+1,3))

      allocate(du    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,du_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      du = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,du,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),6,nxx,nyy,nzz,igrid,du,bcnd
     .            ,icomp=(/IAX,IJX/),is_cnv=.false.,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=bc_order,arr_cov=du_cov)

cc      da_cov => du(:,:,:,1:3)
cc      dj_cnv => du(:,:,:,4:6)
      da_cov = du(:,:,:,1:3)
      dj_cnv = du(:,:,:,4:6)

cc      call setMGBC(max(0,gpos),3,nxx,nyy,nzz,igrid,da_cnv,bcnd(:,1:3)
cccc     .            ,icomp=(/IAX/),is_cnv=.false.,is_vec=.not.si_car
cc     .            ,icomp=(/IAX/),is_cnv=.false.,is_vec=.true.
cc     .            ,result_is_vec=.true.,iorder=bc_order,arr_cov=da_cov)
cc
cc      call setMGBC(max(0,gpos),3,nxx,nyy,nzz,igrid,dj_cnv,bcnd(:,4:6)
cccc     .            ,icomp=(/IJX/),is_cnv=.true.,is_vec=.not.si_car
cc     .            ,icomp=(/IJX/),is_cnv=.true.,is_vec=.true.
cc     .            ,result_is_vec=.true.,iorder=bc_order,arr_cov=das)
cc
cc      da_cov = du(:,:,:,1:3)
cc      dj_cnv = du(:,:,:,4:6)

c Auxiliary quantities

      ve0_cnv => gve0%grid(igrid)%array
      bs0_cnv => gbs0%grid(igrid)%array
      rho0    => grho0%grid(igrid)%array
      heta0   => gheta%grid(igrid)%array

      dj_cnv(:,:,:,1) = di*dj_cnv(:,:,:,1)/rho0(:,:,:,1)        !dj => di*dj/rho (contravariant)
      dj_cnv(:,:,:,2) = di*dj_cnv(:,:,:,2)/rho0(:,:,:,1)        !dj => di*dj/rho (contravariant)
      dj_cnv(:,:,:,3) = di*dj_cnv(:,:,:,3)/rho0(:,:,:,1)        !dj => di*dj/rho (contravariant)

cc      das(:,:,:,1) = da_cov(:,:,:,1) + de**2*das(:,:,:,1)/rho0(:,:,:,1)  !da_star (covariant)
cc      das(:,:,:,2) = da_cov(:,:,:,2) + de**2*das(:,:,:,2)/rho0(:,:,:,1)  !da_star (covariant)
cc      das(:,:,:,3) = da_cov(:,:,:,3) + de**2*das(:,:,:,3)/rho0(:,:,:,1)  !da_star (covariant)

      das(:,:,:,1) = du_cov(:,:,:,1)+de**2*du_cov(:,:,:,4)/rho0(:,:,:,1)  !da_star (covariant)
      das(:,:,:,2) = du_cov(:,:,:,2)+de**2*du_cov(:,:,:,5)/rho0(:,:,:,1)  !da_star (covariant)
      das(:,:,:,3) = du_cov(:,:,:,3)+de**2*du_cov(:,:,:,6)/rho0(:,:,:,1)  !da_star (covariant)

      !Magnetic field db_cov
      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax
            cnv = curl(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,da_cov)
            db_cov(i,j,k,:) =
     .           matmul(gmetric%grid(igrid)%gsub(i,j,k,:,:),cnv)
          enddo
        enddo
      enddo

      call setMGBC(max(0,gpos),3,nxx,nyy,nzz,igrid,db_cov,bcs(:,IBX:IBZ)
cc     .            ,icomp=(/IBX/),is_cnv=.false.,is_vec=.not.si_car
     .            ,icomp=(/IBX/),is_cnv=.false.,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=bc_order)

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

            jac = gmetric%grid(igrid)%jac(i,j,k)

            ijk = i + nxx*(j-1) + nxx*nyy*(k-1)

            if (vol_wgt) then
              vol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              vol = 1d0
            endif

            !Temporal piece
            a_res = cnp*das(i,j,k,:)/dt

            !Resistive piece
            cov = veclap_cov(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,da_cov
     .                      ,vol=.false.)

            etal= geta%grid(igrid)%array(i,j,k,1)

            a_res = a_res - alpha*etal*cov

            !Curlcurl piece: curl(curl) = grad(div)-veclap
cc            curlcurla = curl(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
cc     .                     ,db_cov(:,:,:,1)
cc     .                     ,db_cov(:,:,:,2)
cc     .                     ,db_cov(:,:,:,3))

            cov =-cov
     .          +gradDiv(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,da_cnv
     .                  ,vol=.false.)
cc     .          +gradDiv_vrtx(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,da_cnv
cc     .                  ,vol=.false.)

            curlcurla = matmul(gmetric%grid(igrid)%gsup(i,j,k,:,:),cov)

            !Convective part (-ve x curl(da_star)): component 1
            flxip = -ve0_cnv(i,j,k,2)*0.5*(das(ip,j,k,2)+das(i,j,k,2))
     .              -ve0_cnv(i,j,k,3)*0.5*(das(ip,j,k,3)+das(i,j,k,3))

            flxim = -ve0_cnv(i,j,k,2)*0.5*(das(im,j,k,2)+das(i,j,k,2))
     .              -ve0_cnv(i,j,k,3)*0.5*(das(im,j,k,3)+das(i,j,k,3))

            flxjp = 0.5*(    ve0_cnv(i,j,k,2)
     .                  +abs(ve0_cnv(i,j,k,2)))*das(i,j ,k,1)
     .             +0.5*(    ve0_cnv(i,j,k,2)
     .                  -abs(ve0_cnv(i,j,k,2)))*das(i,jp,k,1)
            flxjm = 0.5*(    ve0_cnv(i,j,k,2)
     .                  +abs(ve0_cnv(i,j,k,2)))*das(i,jm,k,1)
     .             +0.5*(    ve0_cnv(i,j,k,2)
     .                  -abs(ve0_cnv(i,j,k,2)))*das(i,j ,k,1)

            flxkp = 0.5*(    ve0_cnv(i,j,k,3)
     .                  +abs(ve0_cnv(i,j,k,3)))*das(i,j,k ,1)
     .             +0.5*(    ve0_cnv(i,j,k,3)            
     .                  -abs(ve0_cnv(i,j,k,3)))*das(i,j,kp,1)
            flxkm = 0.5*(    ve0_cnv(i,j,k,3)            
     .                  +abs(ve0_cnv(i,j,k,3)))*das(i,j,km,1)
     .             +0.5*(    ve0_cnv(i,j,k,3)            
     .                  -abs(ve0_cnv(i,j,k,3)))*das(i,j,k ,1)

            cov(1) =  (flxip-flxim)/grid_params%dxh(ig)
     .               +(flxjp-flxjm)/grid_params%dyh(jg)
     .               +(flxkp-flxkm)/grid_params%dzh(kg)

            !Convective part (-ve x curl(da_star)): component 2
            flxip = 0.5*(    ve0_cnv(i,j,k,1)
     .                  +abs(ve0_cnv(i,j,k,1)))*das(i ,j,k,2)
     .             +0.5*(    ve0_cnv(i,j,k,1)        
     .                  -abs(ve0_cnv(i,j,k,1)))*das(ip,j,k,2)
            flxim = 0.5*(    ve0_cnv(i,j,k,1)        
     .                  +abs(ve0_cnv(i,j,k,1)))*das(im,j,k,2)
     .             +0.5*(    ve0_cnv(i,j,k,1)        
     .                  -abs(ve0_cnv(i,j,k,1)))*das(i ,j,k,2)

            flxjp = -ve0_cnv(i,j,k,1)*0.5*(das(i,jp,k,1)+das(i,j,k,1))
     .              -ve0_cnv(i,j,k,3)*0.5*(das(i,jp,k,3)+das(i,j,k,3))

            flxjm = -ve0_cnv(i,j,k,1)*0.5*(das(i,jm,k,1)+das(i,j,k,1))
     .              -ve0_cnv(i,j,k,3)*0.5*(das(i,jm,k,3)+das(i,j,k,3))

            flxkp = 0.5*(    ve0_cnv(i,j,k,3)
     .                  +abs(ve0_cnv(i,j,k,3)))*das(i,j,k ,2)
     .             +0.5*(    ve0_cnv(i,j,k,3)            
     .                  -abs(ve0_cnv(i,j,k,3)))*das(i,j,kp,2)
            flxkm = 0.5*(    ve0_cnv(i,j,k,3)            
     .                  +abs(ve0_cnv(i,j,k,3)))*das(i,j,km,2)
     .             +0.5*(    ve0_cnv(i,j,k,3)            
     .                  -abs(ve0_cnv(i,j,k,3)))*das(i,j,k ,2)

            cov(2) =  (flxip-flxim)/grid_params%dxh(ig)
     .               +(flxjp-flxjm)/grid_params%dyh(jg)
     .               +(flxkp-flxkm)/grid_params%dzh(kg)

            !Convective part (-ve x curl(da_star)): component 3
            flxip = 0.5*(    ve0_cnv(i,j,k,1)
     .                  +abs(ve0_cnv(i,j,k,1)))*das(i ,j,k,3)
     .             +0.5*(    ve0_cnv(i,j,k,1)        
     .                  -abs(ve0_cnv(i,j,k,1)))*das(ip,j,k,3)
            flxim = 0.5*(    ve0_cnv(i,j,k,1)        
     .                  +abs(ve0_cnv(i,j,k,1)))*das(im,j,k,3)
     .             +0.5*(    ve0_cnv(i,j,k,1)        
     .                  -abs(ve0_cnv(i,j,k,1)))*das(i ,j,k,3)

            flxjp = 0.5*(    ve0_cnv(i,j,k,2)
     .                  +abs(ve0_cnv(i,j,k,2)))*das(i,j ,k,3)
     .             +0.5*(    ve0_cnv(i,j,k,2)
     .                  -abs(ve0_cnv(i,j,k,2)))*das(i,jp,k,3)
            flxjm = 0.5*(    ve0_cnv(i,j,k,2)
     .                  +abs(ve0_cnv(i,j,k,2)))*das(i,jm,k,3)
     .             +0.5*(    ve0_cnv(i,j,k,2)
     .                  -abs(ve0_cnv(i,j,k,2)))*das(i,j ,k,3)

            flxkp = -ve0_cnv(i,j,k,1)*0.5*(das(i,j,kp,1)+das(i,j,k,1))
     .              -ve0_cnv(i,j,k,2)*0.5*(das(i,j,kp,2)+das(i,j,k,2))
                                                   
            flxkm = -ve0_cnv(i,j,k,1)*0.5*(das(i,j,km,1)+das(i,j,k,1))
     .              -ve0_cnv(i,j,k,2)*0.5*(das(i,j,km,2)+das(i,j,k,2))

            cov(3) =( (flxip-flxim)/grid_params%dxh(ig)
     .               +(flxjp-flxjm)/grid_params%dyh(jg)
     .               +(flxkp-flxkm)/grid_params%dzh(kg))/jac

            a_res = a_res + alpha*cov

            !Electron Hall piece: dj x B_star
cc            cnv = di*curlcurla/rho0(i,j,k,1)
            cnv = dj_cnv(i,j,k,:)

            cov(1) = (cnv(2)*bs0_cnv(i,j,k,3)
     .               -cnv(3)*bs0_cnv(i,j,k,2))/jac
            cov(2) = (cnv(3)*bs0_cnv(i,j,k,1)
     .               -cnv(1)*bs0_cnv(i,j,k,3))/jac
            cov(3) = (cnv(1)*bs0_cnv(i,j,k,2)
     .               -cnv(2)*bs0_cnv(i,j,k,1))/jac

            a_res = a_res + alpha*cov

            !Electron pressure tensor piece: di/rho*div(pe(dj))
cc            cnv= EOM_divPe_ijk(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
cc     .                        ,dj_cnv)
cc            cnv= heta*EOM_divPi_ijk(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
cc     .                        ,dj_cnv,rho0)
            cnv= EOM_divPe_ijk(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                        ,dj_cnv,heta0)

            a_res = a_res + alpha*di/rho0(i,j,k,1)
     .             *matmul(gmetric%grid(igrid)%gsub(i,j,k,:,:),cnv)

            !Assemble residual
cc            if (si_car) then
cc              call pstop('a_mtvc_hyp','si_car option not implemented')
cc            else

              do ieq=1,3
                y(neq*(ijk-1)+ieq) = a_res(ieq)*vol

cc                y(neq*(ijk-1)+ieq+3) = vol*(dj_cnv(i,j,k,ieq)
cc     .                                 -di*curlcurla(ieq)/rho0(i,j,k,1))
                y(neq*(ijk-1)+ieq+3) = vol*(x(neq*(ijk-1)+ieq+3)
     .                                     -curlcurla(ieq))
              enddo

cc            endif

          enddo
        enddo
      enddo

c End program

cc      deallocate(du,das,db_cov)
cc      nullify(ve0_cnv,da_cov,dj_cnv,rho0,bs0_cnv)
      deallocate(du,du_cov,das
     .          ,db_cov
     .          ,dj_cnv
     .          ,da_cov
     .          ,da_cnv)
      nullify(ve0_cnv,rho0,bs0_cnv)

      call deallocPointers(fpointers)

      end subroutine a_mtvc_hyp

c module v_mtvc_mod
c ##################################################################
      module v_mtvc_mod

      use matvec

      real(8),allocatable,dimension(:,:,:,:) :: dv,dv_cov

      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv,b0_cnv,b0_cov
     .                                         ,rho0,pp0,a0_cnv,j0_cov
cc     .                                         ,db_cov

      integer    :: isig,iisig,ip,im,jp,jm,kp,km,hex,hey,hez
      integer    :: imin,imax,jmin,jmax,kmin,kmax
     .             ,nxx,nyy,nzz,ijk,ijkg,isp,igr,igx,igy,igz
     .             ,i,j,k,ig,jg,kg,ieq

      real(8),allocatable,dimension(:,:)     :: a_sp,a_vol

      real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm,jac

      real(8)    :: vol,psiv(3),psib(3),psit(3),cov(3),cnv(3),car(3)
     .             ,mag(3),veclap(3)

      logical    :: fpointers,is_cnv,is_car,postproc

      contains

c     defineLocalPointers
c     #####################################################################
      subroutine defineLocalPointers(igrid)

      integer :: igrid

c     Begin program

      rho0   => grho0  %grid(igrid)%array

      pp0    => gp0    %grid(igrid)%array

      v0_cnv => gv0    %grid(igrid)%array

      b0_cnv => gb0    %grid(igrid)%array

      b0_cov => gb0_cov%grid(igrid)%array

      a0_cnv => ga0    %grid(igrid)%array   !B + di curl(v)

cc      db_cov => gdb_cov%grid(igrid)%array

      end subroutine defineLocalPointers

c     nullLocalPointers
c     #####################################################################
      subroutine nullLocalPointers

c     Begin program

      nullify(pp0,rho0)
      nullify(b0_cnv,b0_cov)
      nullify(v0_cnv,a0_cnv)

      end subroutine nullLocalPointers

c     div_upwd
c     #####################################################################
      function div_upwd(da,half_elem)

        integer    :: half_elem
        real(8)    :: da(0:nxx+1,0:nyy+1,0:nzz+1,*),div_upwd

        integer    :: ip,im,jp,jm,kp,km
        real(8)    :: dxx,dyy,dzz,axp,axm,ayp,aym,azp,azm

c     Begin program

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1
        
        dxx = 2*grid_params%dxh(ig)
        dyy = 2*grid_params%dyh(jg)
        dzz = 2*grid_params%dzh(kg)

        if (half_elem > 0) then
          dxx = grid_params%dx(ig-1)
          dyy = grid_params%dy(jg-1)
          dzz = grid_params%dz(kg-1)
          ip = i
          jp = j
          kp = k
        elseif (half_elem < 0) then
          dxx = grid_params%dx(ig)
          dyy = grid_params%dy(jg)
          dzz = grid_params%dz(kg)
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
      function findPsit() result(psit)

      implicit none

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
     .             ,divip,divim,divjp,divjm,divkp,divkm
     .             ,coeff,gamm,psit(3)

c     Fluxes at faces for calculation of grad(dv.grad(p0))

      flxip =( (dv(i,j,k,1)/jac+dv(ip,j,k,1)/jacip)
     .          *(pp0(ip,j,k,1)-pp0(i,j,k,1))/grid_params%dx(ig)
     .        +(dv(i,j,k,2)/jac+dv(ip,j,k,2)/jacip)
     .          *(pp0(ip,jp,k,1)-pp0(ip,jm,k,1)
     .           +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/grid_params%dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(ip,j,k,3)/jacip)
     .          *(pp0(ip,j,kp,1)-pp0(ip,j,km,1)
     .           +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/grid_params%dzh(kg)/4.)
     .        *0.5

      flxim =( (dv(i,j,k,1)/jac+dv(im,j,k,1)/jacim)
     .          *(pp0(i ,j,k,1)-pp0(im,j,k,1))/grid_params%dx(ig-1)
     .        +(dv(i,j,k,2)/jac+dv(im,j,k,2)/jacim)
     .          *(pp0(im,jp,k,1)-pp0(im,jm,k,1)
     .           +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/grid_params%dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(im,j,k,3)/jacim)
     .          *(pp0(im,j,kp,1)-pp0(im,j,km,1)
     .           +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/grid_params%dzh(kg)/4.)
     .        *0.5

      flxjp =( (dv(i,j,k,1)/jac+dv(i,jp,k,1)/jacjp)
     .          *(pp0(ip,jp,k,1)-pp0(im,jp,k,1)
     .           +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/grid_params%dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,jp,k,2)/jacjp)
     .          *(pp0(i,jp,k,1)-pp0(i,j,k,1))/grid_params%dy(jg)
     .        +(dv(i,j,k,3)/jac+dv(i,jp,k,3)/jacjp)
     .          *(pp0(i,jp,kp,1)-pp0(i,jp,km,1)
     .           +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/grid_params%dzh(kg)/4.)
     .        *0.5
      flxjm =( (dv(i,j,k,1)/jac+dv(i,jm,k,1)/jacjm)
     .          *(pp0(ip,jm,k,1)-pp0(im,jm,k,1)
     .           +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/grid_params%dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,jm,k,2)/jacjm)
     .          *(pp0(i,j ,k,1)-pp0(i,jm,k,1))/grid_params%dy(jg-1)
     .        +(dv(i,j,k,3)/jac+dv(i,jm,k,3)/jacjm)
     .          *(pp0(i,jm,kp,1)-pp0(i,jm,km,1)
     .           +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/grid_params%dzh(kg)/4.)
     .        *0.5

      flxkp =( (dv(i,j,k,1)/jac+dv(i,j,kp,1)/jackp)
     .          *(pp0(ip,j,kp,1)-pp0(im,j,kp,1)
     .           +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/grid_params%dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,j,kp,2)/jackp)
     .          *(pp0(i,jp,kp,1)-pp0(i,jm,kp,1)
     .           +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/grid_params%dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(i,j,kp,3)/jackp)
     .          *(pp0(i,j,kp,1)-pp0(i,j,k,1))/grid_params%dz(kg) )
     .        *0.5
      flxkm =( (dv(i,j,k,1)/jac+dv(i,j,km,1)/jackm)
     .          *(pp0(ip,j,km,1)-pp0(im,j,km,1)
     .           +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/grid_params%dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,j,km,2)/jackm)
     .          *(pp0(i,jp,km,1)-pp0(i,jm,km,1)
     .           +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/grid_params%dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(i,j,km,3)/jackm)
     .          *(pp0(i,j,k ,1)-pp0(i,j,km,1))/grid_params%dz(kg-1) )
     .        *0.5

c      Fluxes at faces for calculation of grad(gamma*p0*div(dv))

      !!Divergence at faces i+-1/2, etc.
      divip = (dv(ip,j ,k,1)-dv(i ,j ,k,1))/grid_params%dx(ig)
     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .        +dv(ip,jp,k,2)-dv(ip,jm,k,2))/grid_params%dyh(jg)*0.25
     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .        +dv(ip,j,kp,3)-dv(ip,j,km,3))/grid_params%dzh(kg)*0.25
      divim = (dv(i ,j ,k,1)-dv(im,j ,k,1))/grid_params%dx(ig-1)
     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .        +dv(im,jp,k,2)-dv(im,jm,k,2))/grid_params%dyh(jg)*0.25
     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .        +dv(im,j,kp,3)-dv(im,j,km,3))/grid_params%dzh(kg)*0.25

      divjp = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .        +dv(ip,jp,k,1)-dv(im,jp,k,1))/grid_params%dxh(ig)*0.25
     .       +(dv(i ,jp,k,2)-dv(i ,j ,k,2))/grid_params%dy(jg)
     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .        +dv(i,jp,kp,3)-dv(i,jp,km,3))/grid_params%dzh(kg)*0.25
      divjm = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .        +dv(ip,jm,k,1)-dv(im,jm,k,1))/grid_params%dxh(ig)*0.25
     .       +(dv(i ,j ,k,2)-dv(i ,jm,k,2))/grid_params%dy(jg-1)
     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .        +dv(i,jm,kp,3)-dv(i,jm,km,3))/grid_params%dzh(kg)*0.25

      divkp = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .        +dv(ip,j,kp,1)-dv(im,j,kp,1))/grid_params%dxh(ig)*0.25
     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .        +dv(i,jp,kp,2)-dv(i,jm,kp,2))/grid_params%dyh(jg)*0.25
     .       +(dv(i,j ,kp,3)-dv(i,j ,k ,3))/grid_params%dz(kg)
      divkm = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .        +dv(ip,j,km,1)-dv(im,j,km,1))/grid_params%dxh(ig)*0.25
     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .        +dv(i,jp,km,2)-dv(i,jm,km,2))/grid_params%dyh(jg)*0.25
     .       +(dv(i,j ,k ,3)-dv(i,j ,km,3))/grid_params%dz(kg-1)

      if (solve_rho.or.solve_prs) then
        gamm = gamma
      else
        gamm = gamma - 1
      endif

      flxip = flxip
     .      + gamm*(pp0(i,j,k,1)+pp0(ip,j,k,1))*divip/(jac+jacip)
      if (.not.isSP(i,j,k,igx,igy,igz)) then
        flxim = flxim
     .      + gamm*(pp0(i,j,k,1)+pp0(im,j,k,1))*divim/(jac+jacim)
      endif

      flxjp = flxjp
     .      + gamm*(pp0(i,j,k,1)+pp0(i,jp,k,1))*divjp/(jac+jacjp)
      flxjm = flxjm
     .      + gamm*(pp0(i,j,k,1)+pp0(i,jm,k,1))*divjm/(jac+jacjm)

      flxkp = flxkp
     .      + gamm*(pp0(i,j,k,1)+pp0(i,j,kp,1))*divkp/(jac+jackp)
      flxkm = flxkm
     .      + gamm*(pp0(i,j,k,1)+pp0(i,j,km,1))*divkm/(jac+jackm)

c     Transform to contravariant

      cov(1) = (flxip - flxim)/grid_params%dxh(ig)
      cov(2) = (flxjp - flxjm)/grid_params%dyh(jg)
      cov(3) = (flxkp - flxkm)/grid_params%dzh(kg)

      call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                            ,cov(1),cov(2),cov(3)
     .                            ,cnv(1),cnv(2),cnv(3),is_cnv)

      psit = -alpha**2*cnv

      end function findPsit

c     findPsib
c     #############################################################
      function findPsib() result(psib)

      implicit none

      integer :: ig,jg,kg
      real(8) :: coeff,jac,dxx,dyy,dzz,ijac,psib(3)
     .          ,j0_cnv(3),acnv(3)
     .          ,acnvip(3),acnvim(3)
     .          ,acnvjp(3),acnvjm(3)
     .          ,acnvkp(3),acnvkm(3)
     .          ,acovip(3),acovim(3)
     .          ,acovjp(3),acovjm(3)
     .          ,acovkp(3),acovkm(3)

      if (nc_eom_jxb.and.(.not.bcSP()) ) then  !Force conservative form in SP geometries
cc      if (.false.) then
cc      if (.true.) then
        call getMGmap(i,j,k,igr,igr,igr,ig,jg,kg)

        jac = gmetric%grid(igr)%jac(i,j,k)
        ijac = 1d0/jac

        dxx = grid_params%dxh(ig)
        dyy = grid_params%dyh(jg)
        dzz = grid_params%dzh(kg)

        j0_cnv = curl(i,j,k,nxx,nyy,nzz,igr,igr,igr,b0_cov)

        acnv =
     .      U_Bv_ijk(i  ,j,k,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov,a0_cnv,0)

        acnvip =
     .      U_Bv_ijk(i  ,j,k,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov,a0_cnv,1)
        acnvim =
     .      U_Bv_ijk(i-1,j,k,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov,a0_cnv,1)

        acnvjp =
     .      U_Bv_ijk(i,j  ,k,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov,a0_cnv,2)
        acnvjm =
     .      U_Bv_ijk(i,j-1,k,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov,a0_cnv,2)

        acnvkp =
     .      U_Bv_ijk(i,j,k  ,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov,a0_cnv,3)
        acnvkm =
     .      U_Bv_ijk(i,j,k-1,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov,a0_cnv,3)

        call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                        ,acovip(1),acovip(2),acovip(3)
     .                        ,acnvip(1),acnvip(2),acnvip(3)
     .                        ,.false.,half_elem=1)

        call transformFromCurvToCurv(i-1,j,k,igr,igr,igr
     .                        ,acovim(1),acovim(2),acovim(3)
     .                        ,acnvim(1),acnvim(2),acnvim(3)
     .                        ,.false.,half_elem=1)

        call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                        ,acovjp(1),acovjp(2),acovjp(3)
     .                        ,acnvjp(1),acnvjp(2),acnvjp(3)
     .                        ,.false.,half_elem=2)

        call transformFromCurvToCurv(i,j-1,k,igr,igr,igr
     .                        ,acovjm(1),acovjm(2),acovjm(3)
     .                        ,acnvjm(1),acnvjm(2),acnvjm(3)
     .                        ,.false.,half_elem=2)

        call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                        ,acovkp(1),acovkp(2),acovkp(3)
     .                        ,acnvkp(1),acnvkp(2),acnvkp(3)
     .                        ,.false.,half_elem=3)

        call transformFromCurvToCurv(i,j,k-1,igr,igr,igr
     .                        ,acovkm(1),acovkm(2),acovkm(3)
     .                        ,acnvkm(1),acnvkm(2),acnvkm(3)
     .                        ,.false.,half_elem=3)

        cov(1) =-b0_cnv(i,j,k,2)*( (acovip(2)-acovim(2))/dxx
     .                            -(acovjp(1)-acovjm(1))/dyy)
     .          -b0_cnv(i,j,k,3)*( (acovip(3)-acovim(3))/dxx
     .                            -(acovkp(1)-acovkm(1))/dzz)
     .          +(j0_cnv(2)*acnv(3)-j0_cnv(3)*acnv(2))

        cov(2) =-b0_cnv(i,j,k,1)*( (acovjp(1)-acovjm(1))/dyy
     .                            -(acovip(2)-acovim(2))/dxx)
     .          -b0_cnv(i,j,k,3)*( (acovjp(3)-acovjm(3))/dyy
     .                            -(acovkp(2)-acovkm(2))/dzz)
     .          +(j0_cnv(3)*acnv(1)-j0_cnv(1)*acnv(3))

        cov(3) =-b0_cnv(i,j,k,1)*( (acovkp(1)-acovkm(1))/dzz
     .                            -(acovip(3)-acovim(3))/dxx)
     .          -b0_cnv(i,j,k,2)*( (acovkp(2)-acovkm(2))/dzz
     .                            -(acovjp(3)-acovjm(3))/dyy)
     .          +(j0_cnv(1)*acnv(2)-j0_cnv(2)*acnv(1))

        cov = cov*ijac

        call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)


      else
        cnv = div_tensor(i,j,k,nxx,nyy,nzz,igr,igr,igr,alt_eom
     .                  ,tensor_x,tensor_y,tensor_z,vol=.false.)
      endif

      psib = alpha*cnv

      end function findPsib

c     tensor_x
c     #############################################################
      subroutine tensor_x(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ip,igrid,iisig,iig,ijk
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,gsuper(3,3),vol
     .               ,acnv(3),acnvp(3),b0cov(3),b0cnv(3),scalar_prod

c     Begin program

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
          acnv = U_Bv_ijk(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov
     .                   ,a0_cnv,1)
        else
          acnv = U_Bv_ijk(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov
     .                   ,a0_cnv,0)
        endif

cc        if (precon == 's2') then
cc          iisig = istart (igrid)
cc          ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)
cc          iig = 3*(ijk-1) + iisig - 1
cc
cc          if (vol_wgt) then
cc            vol = gmetric%grid(igrid)%dvol(i,j,k)
cc          else
cc            vol = 1d0
cc          endif
cc
cc          acnv = matmul(transpose(b_diag(:,iig+1:iig+3)),acnv)*vol
cc        endif

        b0cnv = 0.5*(b0_cnv(ip,j,k,:)+b0_cnv(i,j,k,:))
cc        b0cov = 0.5*(b0_cov(ip,j,k,:)+b0_cov(i,j,k,:))

cc        b0cnv(1)=0.5*(b0_cnv(ip,j,k,1)*ijacp+b0_cnv(i,j,k,1)*ijac0)*jac
cc        b0cnv(2)=0.5*(b0_cnv(ip,j,k,2)      +b0_cnv(i,j,k,2))
cc        b0cnv(3)=0.5*(b0_cnv(ip,j,k,3)*ijacp+b0_cnv(i,j,k,3)*ijac0)*jac

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

        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: jp,igrid,iisig,iig,ijk
        real(8)    :: jac,ijac,gsuper(3,3),acnv(3),b0cov(3),b0cnv(3)
     .               ,scalar_prod,vol

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
          acnv = U_Bv_ijk(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov
     .                   ,a0_cnv,2)
        else
          acnv = U_Bv_ijk(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov
     .                   ,a0_cnv,0)
        endif

cc        if (precon == 's2') then
cc          iisig = istart (igrid)
cc          ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)
cc          iig = 3*(ijk-1) + iisig - 1
cc
cc          if (vol_wgt) then
cc            vol = gmetric%grid(igrid)%dvol(i,j,k)
cc          else
cc            vol = 1d0
cc          endif
cc
cc          acnv = matmul(transpose(b_diag(:,iig+1:iig+3)),acnv)*vol
cc        endif

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

        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: kp,igrid,iisig,iig,ijk
        real(8)    :: jac,ijac,gsuper(3,3),acnv(3),b0cov(3),b0cnv(3)
     .               ,scalar_prod,vol

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
          acnv = U_Bv_ijk(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov
     .                   ,a0_cnv,3)
        else
          acnv = U_Bv_ijk(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv,dv_cov
     .                   ,a0_cnv,0)
        endif

cc        if (precon == 's2') then
cc          iisig = istart (igrid)
cc          ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)
cc          iig = 3*(ijk-1) + iisig - 1
cc
cc          if (vol_wgt) then
cc            vol = gmetric%grid(igrid)%dvol(i,j,k)
cc          else
cc            vol = 1d0
cc          endif
cc
cc          acnv = matmul(transpose(b_diag(:,iig+1:iig+3)),acnv)*vol
cc        endif

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

      end module v_mtvc_mod

c v_mtvc
c####################################################################
      subroutine v_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the SMALL-FLOW velocity SI system.
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

      use nlfunction_setup

      use v_mtvc_mod

      use local_BCS_variables, ONLY: ExB_flow_bc!,lag_pinch_bc

      implicit none

c Call variables

      integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8) :: x(ntot),y(ntot),minv,idx,idy,idz
     .          ,idxcnv,idycnv,idzcnv,vxx,vyy,vzz

c Local variables

      integer :: iig
      real(8) :: b2,bbcnv(3),bbcov(3),etal,E00(3)

      real(8),allocatable,dimension(:,:,:,:) :: v0,dpp

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

      call defineLocalPointers(igrid)

c ExB flow BC (if applicable)

      allocate(v0(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      v0 = 0d0

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dpp   (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
     .                       ,.false.)

      if (di > 0) then
        call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=1,arr_cov=dv_cov
     .            ,arr0=v0)
      else
        call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=bc_order,arr_cov=dv_cov
     .            ,arr0=v0)
      endif

      dpp(:,:,:,1)=rho0(:,:,:,1)*dv(:,:,:,1)
      dpp(:,:,:,2)=rho0(:,:,:,1)*dv(:,:,:,2)
      dpp(:,:,:,3)=rho0(:,:,:,1)*dv(:,:,:,3)

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

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (vol_wgt) then
              vol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************

            psiv = D_v(i,j,k,nxx,nyy,nzz,igr,rho0,v0_cnv,dv,dpp
     .                ,.not.nc_eom_v)

            !P_si^T  ******************************

            psit = findPsit()

            !P_si^B  ******************************

            psib = findPsib()

            !Diagonal scaling
             
            if (pc_type == 's2') then
cc              cnv = cnv*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol
              psit = psit*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol
              psib = psib*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol

              cnv = (psit + psib)
            else
              cnv = dt/cnp*(psit + psib)
            endif

            !Form total operator

            if (nc_eom_v) then
              cnv = (psiv + cnv/rho0(i,j,k,1))*vol
            else
              cnv = (psiv + cnv)*vol
            endif

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = cnv(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(dv,dv_cov,v0,dpp)

      call nullLocalPointers

      call deallocPointers(fpointers)

      end subroutine v_mtvc

c v_mtvc_hyp
c####################################################################
      subroutine v_mtvc_hyp(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the Hall-MHD velocity SI system.
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

      use nlfunction_setup

      use v_mtvc_mod

      use local_BCS_variables, ONLY: ExB_flow_bc!,lag_pinch_bc

      implicit none

c Call variables

      integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8) :: x(ntot),y(ntot),minv,idx,idy,idz
     .          ,idxcnv,idycnv,idzcnv,vxx,vyy,vzz

c Local variables

      integer :: iig
      real(8) :: b2,bbcnv(3),bbcov(3),etal,E00(3),v_res(3),divpe(3)
     .          ,dj(3),aip(3),aim(3),ajp(3),ajm(3),akp(3),akm(3)

      real(8),pointer,dimension(:,:,:,:) :: v0,dE_cnv,du,du_cov,dpp
     .                                     ,dve_cnv,divpi,dB_cnv,dB_cov
     .                                     ,heta0

c Begin program

      if (neq /= 6) call pstop('v_mtvc_hyp','neq /= 6; aborting...')

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

      call defineLocalPointers(igrid)

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(du    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,du_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      du = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,du,igrid
     .                       ,.false.)

      allocate(dv     (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dv_cov (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dve_cnv(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dpp    (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,divpi  (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dE_cnv (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dB_cnv (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dB_cov (0:nxx+1,0:nyy+1,0:nzz+1,3))

c BCs

      !ExB flow BC (if applicable)
      allocate(v0(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      v0 = 0d0

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,du,bcnd
     .            ,icomp=(/IVX,IJX/),is_cnv=is_cnv,is_vec=.true.
cc     .            ,result_is_vec=.true.,iorder=bc_order,arr_cov=du_cov
     .            ,result_is_vec=.true.,iorder=1,arr_cov=du_cov
     .            ,arr0=v0)

      dv      = du    (:,:,:,1:3)
      dv_cov  = du_cov(:,:,:,1:3)

      dve_cnv(:,:,:,1) = dv(:,:,:,1) - di*du(:,:,:,4)/rho0(:,:,:,1)
      dve_cnv(:,:,:,2) = dv(:,:,:,2) - di*du(:,:,:,5)/rho0(:,:,:,1)
      dve_cnv(:,:,:,3) = dv(:,:,:,3) - di*du(:,:,:,6)/rho0(:,:,:,1)

c Find div(Pi(dv))

      if (nu > 0d0) then
        divpi=EOM_divPi(nxx,nyy,nzz,igrid,igrid,igrid,dv
     .                 ,gvis%grid(igrid)%array(:,:,:,1))
      else
        divpi = 0d0
      endif

c Find dE=U_Av(dv) (with BCs)

      if (fake_ve) then

cc        do k = kmin,kmax
cc          do j = jmin,jmax
cc            do i = imin,imax
cc              dE_cnv(i,j,k,:)=matmul(gmetric%grid(igrid)%gsup(i,j,k,:,:)
cc     .                              ,U_Av_ijk(i,j,k,nxx,nyy,nzz
cc     .                                        ,igrid,igrid,igrid,dv
cc     .                                        ,dv_cov,a0_cnv,he=0
cc     .                                        ,stress=.false.))
cc            enddo
cc          enddo
cc        enddo

        dE_cnv = di*dv/dt

        dE_cnv(:,:,:,1)=di*(dv(:,:,:,1)/dt+divpi(:,:,:,1)/rho0(:,:,:,1))
        dE_cnv(:,:,:,2)=di*(dv(:,:,:,2)/dt+divpi(:,:,:,2)/rho0(:,:,:,1))
        dE_cnv(:,:,:,3)=di*(dv(:,:,:,3)/dt+divpi(:,:,:,3)/rho0(:,:,:,1))

        call setMGBC(max(0,gpos),3,nxx,nyy,nzz,igrid,dE_cnv
     .            ,gv%aux%vec_list(IENI)%bconds
     .            ,icomp=(/IVX/),is_cnv=.true.,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=1)

cc      else
cc
cc        do k = kmin,kmax
cc          do j = jmin,jmax
cc            do i = imin,imax
cc              dB_cnv(i,j,k,:) = -dt/cnp*U_Bv_ijk(i,j,k,nxx,nyy,nzz
cc     .                                 ,igrid,igrid,igrid,dv
cc     .                                 ,dv_cov,a0_cnv,he=0
cc     .                                 ,strs=.false.)
cc            enddo
cc          enddo
cc        enddo
cc
cc        call setMGBC(max(0,gpos),3,nxx,nyy,nzz,igrid,dB_cnv
cc     .            ,gv%aux%vec_list(IBCNV)%bconds,arr_cov=dB_cov
cc     .            ,icomp=(/IBX/),is_cnv=.true.,is_vec=.true.
cc     .            ,result_is_vec=.true.,iorder=1)
      endif

c Calculate matrix-vector product

      dpp(:,:,:,1)=rho0(:,:,:,1)*dv(:,:,:,1)
      dpp(:,:,:,2)=rho0(:,:,:,1)*dv(:,:,:,2)
      dpp(:,:,:,3)=rho0(:,:,:,1)*dv(:,:,:,3)

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

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (vol_wgt) then
              vol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************

            psiv = D_v(i,j,k,nxx,nyy,nzz,igr,rho0,v0_cnv,dv,dpp
     .                ,.not.nc_eom_v,stress=.false.)

            !Ion viscosity piece ******************

            if (nc_eom_v) then
              psiv = psiv + alpha*divpi(i,j,k,:)/rho0(i,j,k,1)
            else
              psiv = psiv + alpha*divpi(i,j,k,:)
            endif

            !Electron viscosity piece *************

cc            divpe = heta*EOM_divPe_ijk(i,j,k,nxx,nyy,nzz,igr,igr,igr
cc     .                                ,dve_cnv,rho0)
            divpe = EOM_divPe_ijk(i,j,k,nxx,nyy,nzz,igr,igr,igr,dve_cnv
     .                           ,heta0)

            if (nc_eom_v) then
              psiv = psiv + alpha*divpe/rho0(i,j,k,1)
            else
              psiv = psiv + alpha*divpe
            endif

            !P_si^T  ******************************

            psit = findPsit()

            !P_si^B  ******************************

            psib = findPsib()

            !Diagonal scaling
            cnv = psit + psib 

            if (pc_type == 's2') then
              cnv = max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol*cnv
            else
              cnv = dt/cnp*cnv
            endif

            !Form total operator
            if (nc_eom_v) then
              v_res = (psiv + cnv/rho0(i,j,k,1))
            else
              v_res = (psiv + cnv)
            endif

            !Find current update, dj=-dt*curl(curl(dE))
            if (fake_ve) then

              dj = -dt*curlcurl(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                         ,dE_cnv,vol=.false.)
cc              dj = dt*veclaplacian(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
cc     .                          ,dE_cnv,vol=.false.)

            else

cc              !Find dj=curl(db)
cc              dj=curl(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dB_cov)

              !Find dj=-dt*curl(U_Bv_ijk.dv)
              cnv = U_Bv_ijk(i ,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv
     .                  ,dv_cov,a0_cnv,1)
              call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                    ,aip(1),aip(2),aip(3)
     .                                    ,cnv(1),cnv(2),cnv(3)
     .                                    ,.false.,half_elem=1)

              cnv = U_Bv_ijk(im,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv
     .                  ,dv_cov,a0_cnv,1)
              call transformFromCurvToCurv(im,j,k,igrid,igrid,igrid
     .                                    ,aim(1),aim(2),aim(3)
     .                                    ,cnv(1),cnv(2),cnv(3)
     .                                    ,.false.,half_elem=1)
                               
              cnv = U_Bv_ijk(i,j ,k,nxx,nyy,nzz,igrid,igrid,igrid,dv
     .                  ,dv_cov,a0_cnv,2)
              call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                    ,ajp(1),ajp(2),ajp(3)
     .                                    ,cnv(1),cnv(2),cnv(3)
     .                                    ,.false.,half_elem=2)
                               
              cnv = U_Bv_ijk(i,jm,k,nxx,nyy,nzz,igrid,igrid,igrid,dv
     .                  ,dv_cov,a0_cnv,2)
              call transformFromCurvToCurv(i,jm,k,igrid,igrid,igrid
     .                                    ,ajm(1),ajm(2),ajm(3)
     .                                    ,cnv(1),cnv(2),cnv(3)
     .                                    ,.false.,half_elem=2)
                               
              cnv = U_Bv_ijk(i,j,k ,nxx,nyy,nzz,igrid,igrid,igrid,dv
     .                  ,dv_cov,a0_cnv,3)
              call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
     .                                    ,akp(1),akp(2),akp(3)
     .                                    ,cnv(1),cnv(2),cnv(3)
     .                                    ,.false.,half_elem=3)
                               
              cnv = U_Bv_ijk(i,j,km,nxx,nyy,nzz,igrid,igrid,igrid,dv
     .                  ,dv_cov,a0_cnv,3)
              call transformFromCurvToCurv(i,j,km,igrid,igrid,igrid
     .                                    ,akm(1),akm(2),akm(3)
     .                                    ,cnv(1),cnv(2),cnv(3)
     .                                    ,.false.,half_elem=3)

              cnv(1) = (ajp(3)-ajm(3))/grid_params%dyh(jg)
     .                -(akp(2)-akm(2))/grid_params%dzh(kg)

              cnv(2) = (akp(1)-akm(1))/grid_params%dzh(kg)
     .                -(aip(3)-aim(3))/grid_params%dxh(ig)

              cnv(3) = (aip(2)-aim(2))/grid_params%dxh(ig)
     .                -(ajp(1)-ajm(1))/grid_params%dyh(jg)

              dj = -dt*cnv/cnp
            endif

            !Assemble matvec
            do ieq=1,3
              y(neq*(ijk-1)+ieq)   = vol*v_res(ieq)
              y(neq*(ijk-1)+ieq+3) = vol*(du(i,j,k,ieq+3)-dj(ieq))
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(du,du_cov,dpp,dv,dv_cov,v0,dE_cnv,dve_cnv,divpi
     .          ,dB_cnv,dB_cov)

      call nullLocalPointers

      call deallocPointers(fpointers)

      end subroutine v_mtvc_hyp

c v_mtvc2
c####################################################################
      subroutine v_mtvc2(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the ?? velocity SI system.
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

      use nlfunction_setup

      use v_mtvc_mod

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz

c Local variables

      integer    :: iig
      real(8),allocatable,dimension(:,:,:,:) :: dpp

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

      a0_cnv => ga0    %grid(igrid)%array   !B + di curl(v)

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dpp   (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
     .                       ,.false.)

      if (di > 0) then
        call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
cc     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.not.si_car
     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=1,arr_cov=dv_cov)
      else
        call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
cc     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.not.si_car
     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=bc_order,arr_cov=dv_cov)
      endif

      dpp(:,:,:,1)=rho0(:,:,:,1)*dv(:,:,:,1)
      dpp(:,:,:,2)=rho0(:,:,:,1)*dv(:,:,:,2)
      dpp(:,:,:,3)=rho0(:,:,:,1)*dv(:,:,:,3)

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

            psiv = D_v(i,j,k,nxx,nyy,nzz,igr,rho0,v0_cnv,dv,dpp
     .                ,.not.nc_eom_v)

            cnv = psiv*vol

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
              do ieq=1,3
                y(neq*(ijk-1)+ieq) = cnv(ieq)
              enddo
cc            endif

          enddo
        enddo
      enddo

c End program

      deallocate(dpp,dv,dv_cov)

      nullify(pp0,rho0)
      nullify(b0_cnv,b0_cov)
      nullify(v0_cnv,a0_cnv)

      call deallocPointers(fpointers)

      end subroutine v_mtvc2

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
cc      use mgarraySetup
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
cc      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
cc     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz
cc
ccc Local variables
cc
cc      integer    :: igr,nxx,nyy,nzz,ntotf,isig
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
cc      use nlfunction_setup
cc
cc      use v_mtvc_mod
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
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
cc            call D_v(dv)
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

c v_mtvc_af
c####################################################################
      subroutine v_mtvc_af(gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for the ARBITRARY-FLOW (AF) velocity SI system.
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

      use nlfunction_setup

      use v_mtvc_mod

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz

c Local variables

      real(8)    :: z(ntot)

      real(8),allocatable,dimension(:,:,:,:) :: dzz,dzp

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

      a0_cnv => ga0    %grid(igrid)%array   !B + di curl(v)


c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dv    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dzz   (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dzp   (0:nxx+1,0:nyy+1,0:nzz+1,neq))

      !Velocity array
      dv = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dv,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
cc     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.not.si_car
     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=bc_order,arr_cov=dv_cov)

      !M_star*dv array
      dzz = 0d0

      call mstar_mtvc(gpos,neq,ntot,x,z,igrid,bcnd)

      call mapMGVectorToArray(max(0,gpos),neq,z,nxx,nyy,nzz,dzz,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dzz,bcnd
cc     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.not.si_car
     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=bc_order)

      dzp(:,:,:,1)=rho0(:,:,:,1)*dzz(:,:,:,1)
      dzp(:,:,:,2)=rho0(:,:,:,1)*dzz(:,:,:,2)
      dzp(:,:,:,3)=rho0(:,:,:,1)*dzz(:,:,:,3)

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

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (vol_wgt) then
              vol = gmetric%grid(igrid)%dvol(i,j,k)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************

            psiv = D_v(i,j,k,nxx,nyy,nzz,igr,rho0,v0_cnv,dzz,dzp
     .                     ,.not.nc_eom_v)

            !P_si^T  ******************************

            psit = findPsit()

            !P_si^B  ******************************

            psib = findPsib()

            !Form total operator

            cnv = psit + psib

            if (nc_eom_v) then
              cnv = (psiv + cnv/rho0(i,j,k,1))*vol
            else
              cnv = (psiv*rho0(i,j,k,1) + cnv)*vol
            endif

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = cnv(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(dv,dv_cov,dzz,dzp)

      nullify(pp0,rho0)
      nullify(b0_cnv,b0_cov)
      nullify(v0_cnv,a0_cnv)

      call deallocPointers(fpointers)

      end subroutine v_mtvc_af

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

      use nlfunction_setup

      use v_mtvc_mod

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz
     .             ,idxcnv,idycnv,idzcnv,vxx,vyy,vzz
     .             ,diff,upwind

c Local variables

      real(8),allocatable,dimension(:,:,:,:) :: da

      logical :: flg_rho,flg_v

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

      flg_rho = .false.      !To prevent interference with v_mtvc_af
      if (.not.associated(rho0  )) then
        flg_rho = .true.
        rho0   => grho0%grid(igrid)%array
      endif

      flg_v = .false.        !To prevent interference with v_mtvc_af
      if (.not.associated(v0_cnv)) then
        flg_v = .true.
        v0_cnv => gv0  %grid(igrid)%array
      endif

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(da(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      da = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,da,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,da,bcnd
cc     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.not.si_car
     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=bc_order)

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

cc            diff = max(eta,chi,dd)
cc            diff = min(eta,chi,dd)

cc            cnv = gradDiv(i,j,k,nxx,nyy,nzz,igx,igy,igz,da,vol=.false.)
cc            cnv = matmul(gmetric%grid(igrid)%gsup(i,j,k,:,:),cnv)
cc            if (diff /= 0d0) cnv = gradDiv2(da)

            do ieq=1,neq
cc              hex = floor(sign(1d0,v0_cnv(i,j,k,ieq)))
              hex = 0

              upwind = v0_cnv(i,j,k,ieq)*div_upwd(da,hex)

              cnv(ieq) = cnp*da(i,j,k,ieq)/dt
     .                 + alpha*upwind
cc     .                 - alpha*diff*cnv(ieq)

              !No volume weighing in M_star
              y(neq*(ijk-1)+ieq) = cnv(ieq)

            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(da)

      if (flg_rho) nullify(rho0)
      if (flg_v  ) nullify(v0_cnv)

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

      use nlfunction_setup

      use v_mtvc_mod

      implicit none

c Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
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
cc     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.not.si_car
     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=bc_order,arr_cov=dv_cov)

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
     .           *((dv(ip,j,k,2)-2*dv(i,j,k,2)+dv(im,j,k,2))
     .             /grid_params%dxh(ig)**2
     .            -(dv(ip,jp,k,1)-dv(ip,jm,k,1)
     .             -dv(im,jp,k,1)+dv(im,jm,k,1))
     .             /4./grid_params%dxh(ig)/grid_params%dyh(jg))

            cnv(2) = dv(i,j,k,2)/dt
cc     .          + di*b0_cnv(i,j,k,3)
     .          + di
     .           *((dv(i,jp,k,1)-2*dv(i,j,k,1)+dv(i,jm,k,1))
     .             /grid_params%dyh(jg)**2
     .            -(dv(ip,jp,k,2)-dv(ip,jm,k,2)
     .             -dv(im,jp,k,2)+dv(im,jm,k,2))
     .             /4./grid_params%dxh(ig)/grid_params%dyh(jg))

            cnv(3) = dv(i,j,k,3)/dt

            cnv = cnv*vol

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
              do ieq=1,3
                y(neq*(ijk-1)+ieq) = cnv(ieq)
              enddo
cc            endif

          enddo
        enddo
      enddo

c End program

      deallocate(dv,dv_cov)

      nullify(b0_cnv)

      call deallocPointers(fpointers)

      end subroutine hall_mtvc
