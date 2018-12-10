c lap_mtvc
c####################################################################
      subroutine lap_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer :: i,j,k,ig,jg,kg,isig,ijk,ijkg,nxx,nyy,nzz
      integer :: imin,imax,jmin,jmax,kmin,kmax,order

      real(8),allocatable,dimension(:,:,:,:) :: xarr

      logical :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(xarr(0:nxx+1,0:nyy+1,0:nzz+1,neq))

cc      xarr = 0d0

      order = order_bc
      if (    abs(bcnd(1,1)) == NEU
     .    .or.abs(bcnd(2,1)) == NEU
     .    .or.abs(bcnd(3,1)) == NEU
     .    .or.abs(bcnd(4,1)) == NEU
     .    .or.abs(bcnd(5,1)) == NEU
     .    .or.abs(bcnd(6,1)) == NEU) order = 0

      call fillArray(g_def,igrid,neq,x,xarr,bcnd,order,gpos=gpos
     .              ,icomp=(/IRHO/))

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            y(ijk) = lap(g_def,i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                  ,xarr(:,:,:,1),vol=MG_get_res_vol_wgt())

          enddo
        enddo
      enddo

c End program

      deallocate(xarr)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine lap_mtvc

c lap_mtvc2
c####################################################################
      subroutine lap_mtvc2(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer :: i,j,k,ig,jg,kg,isig,ijk,ijkg,nxx,nyy,nzz
      integer :: imin,imax,jmin,jmax,kmin,kmax

      real(8),pointer,dimension(:,:,:,:) :: xarr,bz0

      logical :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      call point_to_MG_array(nxx,nyy,nzz,igrid,bz0,gbz0)
 
      allocate(xarr(0:nxx+1,0:nyy+1,0:nzz+1,neq))

cc      xarr = 0d0

      call fillArray(g_def,igrid,neq,x,xarr,bcnd,order_bc,gpos=gpos
     .              ,arr0=bz0)

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            y(ijk) = lap(g_def,i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                  ,xarr(:,:,:,1),vol=MG_get_res_vol_wgt())

          enddo
        enddo
      enddo

c End program

      deallocate(xarr)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine lap_mtvc2

c del_hel
c####################################################################
      subroutine del_hel(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer :: ip,im,jp,jm,kp,km,isig,ijk,ijkg,nxx,nyy,nzz
      integer :: imin,imax,jmin,jmax,kmin,kmax,i,j,k,ig,jg,kg
      real(8) :: xhp,xhm,kk,x1,y1,z1

      real(8),pointer,dimension(:,:,:,:) :: xarr,psi0

      logical    :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(xarr(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      call point_to_MG_array(nxx,nyy,nzz,igrid,psi0,gpsi0)
 
      call fillArray(g_def,igrid,neq,x,xarr,bcnd,order_bc,gpos=gpos
     .              ,arr0=psi0)

c Calculate matrix-vector product

      kk = g_def%params(2)

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

cc            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)
            call getCurvilinearCoordinates(g_def,i,j,k,igrid,igrid,igrid
     .                                    ,ig,jg,kg,x1,y1,z1)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ip = i+1
            im = i-1
            jp = j+1
            jm = j-1
            kp = k+1
            km = k-1

            xhp = 0.5*(g_def%xx(ig)+g_def%xx(ig+1))
            xhm = 0.5*(g_def%xx(ig)+g_def%xx(ig-1))
      
            y(ijk) = 1d0/x1/g_def%dxh(ig)
     .                *( xhp/(1+(kk*xhp)**2)
     .                      *(xarr(ip,j,k,1)-xarr(i,j,k,1))
     .                      /g_def%dx(ig)
     .                  +xhm/(1+(kk*xhm)**2)
     .                      *(xarr(im,j,k,1)-xarr(i,j,k,1))
     .                      /g_def%dx(ig-1))
     .              +1d0/x1**2/g_def%dyh(jg)
     .                *( (xarr(i,jp,k,1)-xarr(i,j,k,1))
     .                   /g_def%dy(jg)
     .                  +(xarr(i,jm,k,1)-xarr(i,j,k,1))
     .                   /g_def%dy(jg-1) )

            if (MG_get_res_vol_wgt())
     .           y(ijk)=y(ijk)*g_def%gmetric%grid(igrid)%dvol(i,j,k)
          enddo
        enddo
      enddo

c End program

      deallocate(xarr)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine del_hel

c del_star
c####################################################################
      subroutine del_star(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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

      integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8) :: x(ntot),y(ntot)
      type(grid_mg_def),pointer :: g_def

c Local variables

      real(8) :: x1,y1,z1
      integer :: ip,im,jp,jm,kp,km,isig,ijk,ijkg,nxx,nyy,nzz
      integer :: imin,imax,jmin,jmax,kmin,kmax,i,j,k,ig,jg,kg

      real(8),allocatable,dimension(:,:,:,:) :: xarr
cc      real(8),pointer    ,dimension(:,:,:)   :: ljac
      real(8),allocatable,dimension(:,:,:)   :: iRR

      logical    :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(xarr(0:nxx+1,0:nyy+1,0:nzz+1,neq))

cc      xarr = 0d0
      call fillArray(g_def,igrid,neq,x,xarr,bcnd,order_bc,gpos=gpos)

c Calculate matrix-vector product

      !Initalize toroidal radius (R=R0+r*cos(theta)=(x^2+y^2))
      allocate(iRR(0:nxx+1,0:nyy+1,0:nzz+1))

      do k = kmin-1,kmax+1
        do j = jmin-1,jmax+1
          do i = imin-1,imax+1
            call getCartesianCoordinates(g_def,i,j,k,igrid,igrid,igrid
     .                                  ,ig,jg,kg,x1,y1,z1)

            iRR(i,j,k) = 1d0/(x1**2+y1**2)
          enddo
        enddo
      enddo

      !Main loop
      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            y(ijk) = lap(g_def,i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                  ,xarr(:,:,:,1),dff=iRR,vol=MG_get_res_vol_wgt())
          enddo
        enddo
      enddo

c End program

      deallocate(xarr,iRR)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine del_star

c test_mtvc
c####################################################################
      subroutine test_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer    :: isig,ijk,ijkg,nxx,nyy,nzz,i,j,k,ig,jg,kg
      integer    ::  imin, imax, jmin, jmax, kmin, kmax

      real(8),allocatable,dimension(:,:,:,:) :: xarr
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8)    :: upwind,nu2,jac

      logical    :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(xarr(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      xarr = 0d0

      call fillArray(g_def,igrid,neq,x,xarr,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IRHO/))

c Map velocity components

      call point_to_MG_array(nxx,nyy,nzz,igrid,v0_cnv,gv0)
 
c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jac = g_def%gmetric%grid(igrid)%jac(i,j,k)

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
cccc     .              - alpha*nu2*lap(g_def,i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,xarr)

cc            y(ijk) = y(ijk)*volume(i,j,k,igrid,igrid,igrid)

            y(ijk) = lap(g_def,i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                   ,xarr(:,:,:,1),vol=MG_get_res_vol_wgt())

          enddo
        enddo
      enddo

c End program

      deallocate(xarr)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine test_mtvc

c tmp_mtvc
c####################################################################
      subroutine tmp_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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

      integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8) :: x(ntot),y(ntot)
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer :: isig,ijk,ijkg,nnx,nny,nnz,i,j,k,ig,jg,kg
      integer :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:) :: dtmp,dprs
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv,rho0,bh0

      real(8),pointer    ,dimension(:,:,:) :: jac

      real(8) :: upwind,nu2,dvol,slap,plap,gm,sgn

      logical :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nnx = g_def%nxv(igrid)
      nny = g_def%nyv(igrid)
      nnz = g_def%nzv(igrid)

      jac => g_def%gmetric%grid(igrid)%jac

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(dtmp(0:nnx+1,0:nny+1,0:nnz+1,neq)
     .        ,dprs(0:nnx+1,0:nny+1,0:nnz+1,neq))

      dtmp = 0d0
      dprs = 0d0
 
      call fillArray(g_def,igrid,neq,x,dprs,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IPRS/))

c Map velocity components

      call point_to_MG_array(nnx,nny,nnz,igrid,v0_cnv,gvs0)
      call point_to_MG_array(nnx,nny,nnz,igrid,rho0  ,grho0)
      call point_to_MG_array(nnx,nny,nnz,igrid,bh0   ,gbhat)
 
      dtmp = dprs/rho0

c Calculate matrix-vector product

      gm = gamma
      sgn =-1d0

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk    = i + nnx*(j-1) + nnx*nny*(k-1)
            ijkg   = ijk + isig - 1

            if (use_p_eom) then !Need this for strong stability
              upwind = flx_advec(g_def,i,j,k,nnx,nny,nnz
     .                          ,igrid,igrid,igrid,v0_cnv,dprs(:,:,:,1)
     .                          ,1,zip_vel=.false.) !Do not use ZIP, since dT is not pos.def.
     .               +(gm-1.)*dprs(i,j,k,1)
     .                       *div(g_def,i,j,k,nnx,nny,nnz
     .                           ,igrid,igrid,igrid,v0_cnv)
            else
              upwind = gm*flx_advec(g_def,i,j,k,nnx,nny,nnz
     .                           ,igrid,igrid,igrid,v0_cnv,dprs(:,:,:,1)
     .                           ,1,zip_vel=.false.) !Do not use ZIP, since dT is not pos.def.
     .               +(gm-1.)
     .                 *flx_nc_advec(g_def,i,j,k,nnx,nny,nnz
     .                              ,igrid,igrid,igrid,v0_cnv
     .                              ,dprs(:,:,:,1),1,reverse=.true.)
            endif
            
            if (MG_get_res_vol_wgt()) then
              dvol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              dvol = 1d0
            endif

            if (chi > 0d0) then
              if (use_p_diff) then
                slap = lap(g_def,i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                         ,dprs(:,:,:,1))
              else
                slap = lap(g_def,i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                         ,dtmp(:,:,:,1))
              endif
              slap = (gamma-1d0)*chi*slap
            else
              slap = 0d0
            endif

            if (chi_par /= chi .and.(.not.lagrangian)) then
              plap = par_diff(g_def,i,j,k,igrid,dtmp(:,:,:,1),bh0
     .                       ,2,"symm",vol=.false.)
              plap = (gamma-1d0)*(chi_par-chi)/a_p*plap !Factor of a_p because operator acts on Te
            else
              plap = 0d0
            endif

cc            y(ijk)=dvol*( cnp/dt*dtmp(i,j,k,1) + alpha*(upwind - slap) )
            y(ijk)=dvol*( cnp/dt*dprs(i,j,k,1)
     $                  + alpha*(upwind - slap - plap) )*jac(i,j,k)

          enddo
        enddo
      enddo

c End program

      deallocate(dtmp,dprs)
      nullify(v0_cnv,rho0,bh0)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine tmp_mtvc

c rho_mtvc
c####################################################################
      subroutine rho_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer    :: isig,ijk,ijkg,nnx,nny,nnz,ip,im,jp,jm,kp,km
      integer    :: imin,imax,jmin,jmax,kmin,kmax,i,j,k,ig,jg,kg
 
      real(8),allocatable,dimension(:,:,:,:) :: drho
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8),pointer    ,dimension(:,:,:) :: jac

      real(8)    :: upwind,slap,dvol,flxip,flxim,flxjp,flxjm,flxkp,flxkm

      logical    :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nnx = g_def%nxv(igrid)
      nny = g_def%nyv(igrid)
      nnz = g_def%nzv(igrid)

      jac => g_def%gmetric%grid(igrid)%jac

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(drho(0:nnx+1,0:nny+1,0:nnz+1,neq))

      drho = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped
      !For finding the diagonal, gpos > 0
      call fillArray(g_def,igrid,neq,x,drho,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IRHO/))

c Map velocity components

      call point_to_MG_array(nnx,nny,nnz,igrid,v0_cnv,gv0)
 
c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk = i + nnx*(j-1) + nnx*nny*(k-1)

            if (MG_get_res_vol_wgt()) then
              dvol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              dvol = 1d0
            endif

            !Conservative advection (Do not use ZIP form because drho is not pos. def.)
            upwind = flx_advec(g_def,i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                        ,v0_cnv,drho(:,:,:,1),1,zip_vel=.false.)

            ! Diffusion
            if (dd > 0d0) then
              slap = dd*lap(g_def,i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                          ,drho(:,:,:,1))
            else
              slap = 0d0
            endif

            y(ijk)=dvol*(cnp*drho(i,j,k,1)/dt + alpha*(upwind - slap))
     .                 *jac(i,j,k)
          enddo
        enddo
      enddo

c End program

      deallocate(drho)
      nullify(v0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine rho_mtvc

c b_mtvc
c####################################################################
      subroutine b_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer    :: isig,ip,im,jp,jm,kp,km,nx,ny,nz
      integer    :: imin,imax,jmin,jmax,kmin,kmax,igx,igy,igz
      integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
     .             ,i,j,k,ig,jg,kg,ieq

      real(8)    :: vol,vlap(3),cnv(3),cov(3),car(3)
      real(8)    :: upwind

      real(8),allocatable,dimension(:,:,:,:) :: db,db_cov,curl_eta_dj
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv,b0_cnv,eta0

      logical    :: fpointers,is_cnv

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      is_cnv = .true.

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos)
     .           ,igrid,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(db         (0:nx+1,0:ny+1,0:nz+1,neq)
     .        ,db_cov     (0:nx+1,0:ny+1,0:nz+1,neq)
     .        ,curl_eta_dj(0:nx+1,0:ny+1,0:nz+1,neq))

      db = 0d0 ; db_cov = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call fillArray(g_def,igrid,neq,x,db,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IBX/),is_cnv=is_cnv,is_vec=.true.
     .              ,result_is_curv=.true.,arr_cov=db_cov)

c Velocity field (including BCs)

      call point_to_MG_array(nx,ny,nz,igrid,v0_cnv,gv0)
      call point_to_MG_array(nx,ny,nz,igrid,eta0,geta)

c Eta*J term

      if (.not.MG_get_find_diag()) then
        curl_eta_dj = curl(g_def,igrid,db_cov)                    !dj_cnv

        curl_eta_dj(:,:,:,1) = eta0(:,:,:,1)*curl_eta_dj(:,:,:,1) !eta*dj_cnv
        curl_eta_dj(:,:,:,2) = eta0(:,:,:,1)*curl_eta_dj(:,:,:,2) !eta*dj_cnv
        curl_eta_dj(:,:,:,3) = eta0(:,:,:,1)*curl_eta_dj(:,:,:,3) !eta*dj_cnv

        curl_eta_dj = XformToCov(g_def,igrid,curl_eta_dj)         !eta*dj_cov
        curl_eta_dj = curl(g_def,igrid,curl_eta_dj)               !curl(eta*dj)
      else
        curl_eta_dj =-veclap(g_def,igrid,db)
        write (*,*) curl_eta_dj(:,:,:,1)
        curl_eta_dj(:,:,:,1) = eta0(:,:,:,1)*curl_eta_dj(:,:,:,1)
        curl_eta_dj(:,:,:,2) = eta0(:,:,:,1)*curl_eta_dj(:,:,:,2)
        curl_eta_dj(:,:,:,3) = eta0(:,:,:,1)*curl_eta_dj(:,:,:,3)
      endif

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk    = i + nx*(j-1) + nx*ny*(k-1)

            if (MG_get_res_vol_wgt()) then
              vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              vol = 1d0
            endif

cc            if (.not.MG_get_find_diag()) then
cc              cnv = curl_bxv     (g_def,i,j,k,igrid,v0_cnv,db,0)
cc            else
              cnv = curl_bxv_upwd(g_def,i,j,k,igrid,v0_cnv,db)
cc            endif

            cnv=(cnp*db(i,j,k,:)/dt+alpha*(cnv+curl_eta_dj(i,j,k,:)))

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = cnv(ieq)*vol
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(db,db_cov,curl_eta_dj)
      nullify(v0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine b_mtvc

c b_mtvc_hyp
c####################################################################
      subroutine b_mtvc_hyp(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer    :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer    :: imin,imax,jmin,jmax,kmin,kmax
      integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
     .             ,i,j,k,ig,jg,kg,ieq

      real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm,jac
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol
     .             ,veclap_b(3),veclap_b2(3),cnv(3),cov(3),car(3)
      real(8)    :: upwind,etal

      real(8),allocatable,dimension(:,:,:,:) :: db
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv,eta0

      logical    :: fpointers,is_cnv

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      if (neq /= 6) call pstop('b_mtvc_hyp','neq /= 6; aborting...')
      
      is_cnv = .true.

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos)
     .           ,igrid,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(db(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      db = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call fillArray(g_def,igrid,neq,x,db,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IBX,IJX/),is_cnv=is_cnv,is_vec=.true.  !Check icomp!!!
     .              ,result_is_curv=.true.)

c Velocity field (including BCs)

      call point_to_MG_array(nxx,nyy,nzz,igrid,v0_cnv,gv0)
      call point_to_MG_array(nxx,nyy,nzz,igrid,eta0,geta)

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jacip  = g_def%gmetric%grid(igrid)%jac(ip,j,k)
            jacim  = g_def%gmetric%grid(igrid)%jac(im,j,k)
            jacjp  = g_def%gmetric%grid(igrid)%jac(i,jp,k)
            jacjm  = g_def%gmetric%grid(igrid)%jac(i,jm,k)
            jackp  = g_def%gmetric%grid(igrid)%jac(i,j,kp)
            jackm  = g_def%gmetric%grid(igrid)%jac(i,j,km)
            jac    = g_def%gmetric%grid(igrid)%jac(i,j,k)

cc            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + SP_flsv
cc            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + 1d-3

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            if (MG_get_res_vol_wgt()) then
              vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              vol = 1d0
            endif

            veclap_b = veclap(g_def,i,j,k,nxx,nyy,nzz
     .                       ,igrid,igrid,igrid,db(:,:,:,1:3)
     .                       ,vol=.false.)

            veclap_b2= veclap(g_def,i,j,k,nxx,nyy,nzz
     .                       ,igrid,igrid,igrid,db(:,:,:,4:6)
     .                       ,vol=.false.)

            etal   = eta0(i,j,k,1)

            if (sm_type == 'gm') then
              cnv =curl_bxv(g_def,i,j,k,igrid,v0_cnv,db(:,:,:,1:3),0)
            else
              cnv =curl_bxv_upwd(g_def,i,j,k,igrid,v0_cnv,db(:,:,:,1:3))
            endif

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = (cnp*db(i,j,k,ieq)/dt
     .                             + alpha*cnv(ieq)
     .                             - alpha*etal*veclap_b(ieq)
     .                             + alpha*di*heta*veclap_b2(ieq))*vol
            enddo

            do ieq=4,6
              y(neq*(ijk-1)+ieq) = (db(i,j,k,ieq)-veclap_b(ieq-3))*vol
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(db)
      nullify(v0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine b_mtvc_hyp

c a_mtvc
c####################################################################
      subroutine a_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer    :: isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz,ierr
      integer    :: imin,imax,jmin,jmax,kmin,kmax
      integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
     .             ,i,j,k,ig,jg,kg,ieq,meth,bcnd2(6,3)

      real(8)    :: upwind,vol,vlap(3),vlap_cnv(3)
     $             ,cnv(3),cov(3),car(3)

      real(8),allocatable,dimension(:,:,:,:) :: da,da_cnv,db_cnv,db_cov
     .                                         ,dj_cnv,dj_cov
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv,eta0

      logical    :: fpointers,is_cnv

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      is_cnv = .false.

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(da    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dj_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      da = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call fillArray(g_def,igrid,neq,x,da,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IAX/),is_cnv=is_cnv,is_vec=.true.
     .              ,result_is_curv=.true.)

c$$$      allocate(da     (0:nxx+1,0:nyy+1,0:nzz+1,neq)
c$$$     .        ,da_cnv (0:nxx+1,0:nyy+1,0:nzz+1,neq)
c$$$     .        ,db_cnv (0:nxx+1,0:nyy+1,0:nzz+1,neq)
c$$$     .        ,db_cov (0:nxx+1,0:nyy+1,0:nzz+1,neq)
c$$$     .        ,dj_cnv (0:nxx+1,0:nyy+1,0:nzz+1,neq)
c$$$     .        ,dj_cov (0:nxx+1,0:nyy+1,0:nzz+1,neq))
c$$$
c$$$      da = 0d0
c$$$
c$$$      call set_force_corners(.true.) !Force corners owing to div(A)=0 BC
c$$$      call fillArray(g_def,igrid,neq,x,da_cnv,bcnd,order_bc,gpos=gpos
c$$$     .              ,icomp=(/IAX/),is_cnv=is_cnv,is_vec=.true.
c$$$     .              ,result_is_curv=.true.,arr_cov=da)
c$$$      call set_force_corners(.false.) !Force corners owing to div(A)=0 BC
c$$$
c$$$      !Currents
c$$$      if (MG_get_find_diag().or.vlap_etaj) then
c$$$c$$$        dj_cov =-veclap_vrtx(g_def,igrid,da,ret_cnv=.false.)
c$$$        !!!!!
c$$$c$$$        dj_cnv = curlcurl_vrtx(g_def,igrid,da)
c$$$c$$$        dj_cov = XformToCov(g_def,igrid,dj_cnv)
c$$$        !!!!!!
c$$$        dj_cnv =-veclap    (g_def,igrid,da_cnv)
c$$$        dj_cov = XformToCov(g_def,igrid,dj_cnv)
c$$$        !!!!!!!
c$$$c$$$        dj_cnv =-veclap_cnv(g_def,igrid,da_cnv)
c$$$c$$$        dj_cov = XformToCov(g_def,igrid,dj_cnv)
c$$$        !!!!!!!
c$$$c$$$        dj_cov =-veclap_cov(g_def,igrid,da)
c$$$      else
c$$$        db_cnv = curl(g_def,igrid,da)
c$$$        db_cov = XformVector(g_def,igrid,db_cnv,.false.)
c$$$c$$$        call setMGBC(g_def,max(0,gpos),3,nxx,nyy,nzz,igrid,db_cnv
c$$$c$$$     .              ,bcs(:,IBX:IBZ),icomp=(/IBX/),is_cnv=.true.
c$$$c$$$     .              ,is_vec=.true.,result_is_curv=.true.
c$$$c$$$     .              ,iorder=order_bc,arr_cov=db_cov)
c$$$
c$$$        dj_cnv = curl(g_def,igrid,db_cov)
c$$$        dj_cov = XformVector(g_def,igrid,dj_cnv,.false.)
c$$$      endif

c Velocity field (including BCs)

      call point_to_MG_array(nxx,nyy,nzz,igrid,v0_cnv,gv0)
      call point_to_MG_array(nxx,nyy,nzz,igrid,eta0,geta)

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk = i + nxx*(j-1) + nxx*nyy*(k-1)

            if (MG_get_res_vol_wgt()) then
              vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              vol = 1d0
            endif

            !Resistive part
cc            vlap = -eta0(i,j,k,1)*dj_cov(i,j,k,:)
            vlap = eta0(i,j,k,1)
     .            *veclap_cov(g_def,i,j,k,nxx,nyy,nzz
     .                       ,igrid,igrid,igrid,da
     .                       ,vol=.false.)
     .            
            !Convective part
            cov = curla_x_v_upwd(g_def,i,j,k,igrid,v0_cnv(i,j,k,:),da)
            
            !Assemble residual
            cov=(cnp*da(i,j,k,:)/dt+alpha*(cov-vlap))*vol

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = cov(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(da,da_cnv,db_cnv,db_cov,dj_cnv,dj_cov,STAT=ierr)
      nullify(v0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine a_mtvc

c a_mtvc_hyp
c####################################################################
      subroutine a_mtvc_hyp(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
c--------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A(psi)x  matrix-free
c     for EMHD Ohm's law in vector potential form, with a hyper-resistive 
c     term. In call:
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

      integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8) :: x(ntot),y(ntot)
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer :: isig,nxx,nyy,nzz
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
     .          ,i,j,k,ig,jg,kg,ieq

      real(8) :: jac,vol,upwind,etal
     .          ,a_res(3),j_res(3),cc_a_cov(3),cnv(3),cov(3),car(3)
     .          ,lap_a_cov(3)

      real(8),pointer,dimension(:,:,:,:) :: du_cnv,du_cov,das
     .                                     ,db_cnv,db_cov,dve_cnv
     .                                     ,ve0_cnv,da_cov,dj_cnv,dj_cov
     .                                     ,rho0,bs0_cnv,da_cnv,heta0
     .                                     ,eta0,djf_cov,djf_cnv,vzr

      logical    :: fpointers,is_cnv

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      if (neq /= 6) call pstop('a_mtvc_hyp','neq /= 6; aborting...')
      
      is_cnv = .false.

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos)
     .           ,igrid,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(du_cnv(0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,du_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      du_cov = 0d0
      du_cnv = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call set_force_corners(.true.) !Force corners owing to div(A)=0 BC
      call fillArray(g_def,igrid,neq,x,du_cnv,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IAX,IAX/),is_cnv=is_cnv,is_vec=.true.
     .              ,result_is_curv=.true.,arr_cov=du_cov)
      call set_force_corners(.false.) !Force corners owing to div(A)=0 BC

c Auxiliary quantities

      allocate(das    (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,db_cnv (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,db_cov (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dj_cnv (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dj_cov (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dve_cnv(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,vzr    (0:nxx+1,0:nyy+1,0:nzz+1,3))

      vzr = 0d0

      !Variables
cc#if defined(pgf90)
      allocate(da_cov (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,da_cnv (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,djf_cov(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,djf_cnv(0:nxx+1,0:nyy+1,0:nzz+1,3))

      da_cov (0:,0:,0:,1:) = du_cov(0:,0:,0:,1:3)
      da_cnv (0:,0:,0:,1:) = du_cnv(0:,0:,0:,1:3)
      djf_cov(0:,0:,0:,1:) = du_cov(0:,0:,0:,4:6)
      djf_cnv(0:,0:,0:,1:) = du_cnv(0:,0:,0:,4:6)
cc#else
c$$$      da_cov (0:,0:,0:,1:) => du_cov(0:,0:,0:,1:3)
c$$$      da_cnv (0:,0:,0:,1:) => du_cnv(0:,0:,0:,1:3)
c$$$      djf_cov(0:,0:,0:,1:) => du_cov(0:,0:,0:,4:6)
c$$$      djf_cnv(0:,0:,0:,1:) => du_cnv(0:,0:,0:,4:6)
cc#endif

      !Coefficients
      call point_to_MG_array(nxx,nyy,nzz,igrid,ve0_cnv,gve0)
      call point_to_MG_array(nxx,nyy,nzz,igrid,bs0_cnv,gbs0)
      call point_to_MG_array(nxx,nyy,nzz,igrid,rho0   ,grho0)
      call point_to_MG_array(nxx,nyy,nzz,igrid,heta0  ,gheta)
      call point_to_MG_array(nxx,nyy,nzz,igrid,eta0   ,geta)

      !Magnetic fields
      db_cnv = curl(g_def,igrid,da_cov)
      db_cov = XformVector(g_def,igrid,db_cnv,.false.)
cc      call setMGBC(g_def,max(0,gpos),3,nxx,nyy,nzz,igrid,db_cnv
cc     .            ,bcs(:,IBX:IBZ)
cc     .            ,icomp=(/IBX/),is_cnv=.true.,is_vec=.true.
cc     .            ,result_is_curv=.true.,iorder=order_bc,arr_cov=db_cov)

      !Currents
      dj_cnv = curl(g_def,igrid,db_cov)
      dj_cov = XformVector(g_def,igrid,dj_cnv,.false.)
cc      call setMGBC(g_def,max(0,gpos),3,nxx,nyy,nzz,igrid,dj_cnv
cc     .            ,bcs(:,IJX:IJZ)
cc     .            ,icomp=(/IJX/),is_cnv=.true.,is_vec=.true.
cc     .            ,result_is_curv=.true.,iorder=order_bc,arr_cov=dj_cov)

      !Electron velocity (only keep current piece)
      dve_cnv = form_ve(vzr,djf_cnv,rho0(:,:,:,1))

      !da_star (covariant)
      das = da_cov - de**2*form_ve(vzr,djf_cov,rho0(:,:,:,1))

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jac = g_def%gmetric%grid(igrid)%jac(i,j,k)

            ijk = i + nxx*(j-1) + nxx*nyy*(k-1)

            if (MG_get_res_vol_wgt()) then
              vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              vol = 1d0
            endif

            !Temporal piece
            a_res = cnp*das(i,j,k,:)/dt
 
            !Resistive part
            etal = eta0(i,j,k,1)

            lap_a_cov = veclap_cov(g_def,i,j,k,nxx,nyy,nzz
     .                            ,igrid,igrid,igrid,da_cov)

            if (vlap_etaj) then
              a_res = a_res - alpha*etal*lap_a_cov
            else
              a_res = a_res + alpha*etal*djf_cov(i,j,k,:)
            endif

            !Convective part (curl(da_star) x ve)
            a_res = a_res
     .            + alpha*curla_x_v_upwd(g_def,i,j,k,igrid
     .                                  ,ve0_cnv(i,j,k,:),das)

            !Electron Hall piece: -dve x B_star
            cov = -crossProduct_ijk(g_def,i,j,k,igrid
     .                             ,dve_cnv(i,j,k,:)
     .                             ,bs0_cnv(i,j,k,:),.true.)

            a_res = a_res + alpha*cov

            !Electron pressure tensor piece:-di/rho*div(Pi(dve))
            cnv = EOM_divPi_ijk(g_def,i,j,k,nxx,nyy,nzz,igrid
     .                         ,dve_cnv,heta0,.false.)
            cov = matmul(g_def%gmetric%grid(igrid)%gsub(i,j,k,:,:),cnv)

            a_res = a_res - alpha*di*cov/rho0(i,j,k,1)

            !Assemble residual
            if (MG_get_find_diag()) then
c$$$            if (MG_get_find_diag()
c$$$     .          .or.nk_current_Newton_iteration()>1) then
              cc_a_cov =-lap_a_cov
cc              cc_a_cov = dj_cov(i,j,k,:)  !Unstable
cc              cc_a_cov = curlcurl_cov_ijk(i,j,k,nxx,nyy,nzz
cc     .                                   ,igx,igy,igz,da_cov)
            else
              cc_a_cov = dj_cov(i,j,k,:)
            endif

            j_res = djf_cov(i,j,k,:) - cc_a_cov

            do ieq=1,3
              y(neq*(ijk-1)+ieq  ) = vol*a_res(ieq)
              y(neq*(ijk-1)+ieq+3) = vol*j_res(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

cc#if defined(pgf90)
      deallocate(da_cov,da_cnv,djf_cov,djf_cnv)
cc#else
cc      nullify   (da_cov,da_cnv,djf_cov,djf_cnv)
cc#endif

      deallocate(du_cnv
     .          ,du_cov
     .          ,vzr
     .          ,das
     .          ,db_cnv
     .          ,db_cov
     .          ,dj_cov
     .          ,dj_cnv
     .          ,dve_cnv)

      nullify(ve0_cnv,rho0,bs0_cnv,heta0,eta0)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine a_mtvc_hyp

c module v_mtvc_mod
c ##################################################################
      module v_mtvc_mod

      use matvec

      real(8),pointer,dimension(:,:,:,:) :: dv,dv_cov,v0_cnv,b0_cnv
     .                                     ,b0_cov,rho0,pp0,a0_cnv
     .                                     ,j0_cov,vis0

      integer :: isig,iisig,ip,im,jp,jm,kp,km,hex,hey,hez

      real(8),allocatable,dimension(:,:) :: a_sp,a_vol

      real(8) :: jacip,jacim,jacjp,jacjm,jackp,jackm,jac

      real(8) :: vol,psiv(3),psib(3),psit(3),cov(3),cnv(3),car(3)
     .          ,mag(3),vlap(3)

      logical :: fpointers,is_cnv,is_car

      contains

c     defineLocalPointers
c     #####################################################################
      subroutine defineLocalPointers(igrid,nx,ny,nz)

      implicit none

      integer :: igrid,nx,ny,nz

c     Begin program

      call point_to_MG_array(nx,ny,nz,igrid,rho0  ,grho0)
      call point_to_MG_array(nx,ny,nz,igrid,pp0   ,gp0)
      call point_to_MG_array(nx,ny,nz,igrid,v0_cnv,gv0)
      call point_to_MG_array(nx,ny,nz,igrid,b0_cnv,gb0)
      call point_to_MG_array(nx,ny,nz,igrid,b0_cov,gb0_cov)
      call point_to_MG_array(nx,ny,nz,igrid,a0_cnv,ga0)
      call point_to_MG_array(nx,ny,nz,igrid,vis0  ,gvis)

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
      function div_upwd(g_def,i,j,k,igr,da,half_elem)

        type(grid_mg_def),pointer :: g_def

c     Call variables

        integer    :: i,j,k,igr,half_elem
        real(8)    :: da(0:,0:,0:,:),div_upwd

c     Local variables

        integer    :: ip,im,jp,jm,kp,km,ig,jg,kg
        real(8)    :: dxx,dyy,dzz,axp,axm,ayp,aym,azp,azm

c     Begin program

        call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1
        
        dxx = 2*g_def%dxh(ig)
        dyy = 2*g_def%dyh(jg)
        dzz = 2*g_def%dzh(kg)

        if (half_elem > 0) then
          dxx = g_def%dx(ig-1)
          dyy = g_def%dy(jg-1)
          dzz = g_def%dz(kg-1)
          ip = i
          jp = j
          kp = k
        elseif (half_elem < 0) then
          dxx = g_def%dx(ig)
          dyy = g_def%dy(jg)
          dzz = g_def%dz(kg)
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
      function findPsit(g_def,i,j,k,igr) result(psit)

      implicit none

c     Call variables

      integer :: i,j,k,igr

      type(grid_mg_def),pointer :: g_def

c     Local variables

      integer :: ig,jg,kg

      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
     .          ,divip,divim,divjp,divjm,divkp,divkm
     .          ,coeff,gamm,psit(3)

c     Fluxes at faces for calculation of grad(dv.grad(p0))

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      flxip =( (dv(i,j,k,1)/jac+dv(ip,j,k,1)/jacip)
     .          *(pp0(ip,j,k,1)-pp0(i,j,k,1))/g_def%dx(ig)
     .        +(dv(i,j,k,2)/jac+dv(ip,j,k,2)/jacip)
     .          *(pp0(ip,jp,k,1)-pp0(ip,jm,k,1)
     .           +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/g_def%dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(ip,j,k,3)/jacip)
     .          *(pp0(ip,j,kp,1)-pp0(ip,j,km,1)
     .           +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/g_def%dzh(kg)/4.)
     .        *0.5

      flxim =( (dv(i,j,k,1)/jac+dv(im,j,k,1)/jacim)
     .          *(pp0(i ,j,k,1)-pp0(im,j,k,1))/g_def%dx(ig-1)
     .        +(dv(i,j,k,2)/jac+dv(im,j,k,2)/jacim)
     .          *(pp0(im,jp,k,1)-pp0(im,jm,k,1)
     .           +pp0(i ,jp,k,1)-pp0(i ,jm,k,1))/g_def%dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(im,j,k,3)/jacim)
     .          *(pp0(im,j,kp,1)-pp0(im,j,km,1)
     .           +pp0(i ,j,kp,1)-pp0(i ,j,km,1))/g_def%dzh(kg)/4.)
     .        *0.5

      flxjp =( (dv(i,j,k,1)/jac+dv(i,jp,k,1)/jacjp)
     .          *(pp0(ip,jp,k,1)-pp0(im,jp,k,1)
     .           +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/g_def%dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,jp,k,2)/jacjp)
     .          *(pp0(i,jp,k,1)-pp0(i,j,k,1))/g_def%dy(jg)
     .        +(dv(i,j,k,3)/jac+dv(i,jp,k,3)/jacjp)
     .          *(pp0(i,jp,kp,1)-pp0(i,jp,km,1)
     .           +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/g_def%dzh(kg)/4.)
     .        *0.5
      flxjm =( (dv(i,j,k,1)/jac+dv(i,jm,k,1)/jacjm)
     .          *(pp0(ip,jm,k,1)-pp0(im,jm,k,1)
     .           +pp0(ip,j ,k,1)-pp0(im,j ,k,1))/g_def%dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,jm,k,2)/jacjm)
     .          *(pp0(i,j ,k,1)-pp0(i,jm,k,1))/g_def%dy(jg-1)
     .        +(dv(i,j,k,3)/jac+dv(i,jm,k,3)/jacjm)
     .          *(pp0(i,jm,kp,1)-pp0(i,jm,km,1)
     .           +pp0(i,j ,kp,1)-pp0(i,j ,km,1))/g_def%dzh(kg)/4.)
     .        *0.5

      flxkp =( (dv(i,j,k,1)/jac+dv(i,j,kp,1)/jackp)
     .          *(pp0(ip,j,kp,1)-pp0(im,j,kp,1)
     .           +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/g_def%dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,j,kp,2)/jackp)
     .          *(pp0(i,jp,kp,1)-pp0(i,jm,kp,1)
     .           +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/g_def%dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(i,j,kp,3)/jackp)
     .          *(pp0(i,j,kp,1)-pp0(i,j,k,1))/g_def%dz(kg) )
     .        *0.5
      flxkm =( (dv(i,j,k,1)/jac+dv(i,j,km,1)/jackm)
     .          *(pp0(ip,j,km,1)-pp0(im,j,km,1)
     .           +pp0(ip,j,k ,1)-pp0(im,j,k ,1))/g_def%dxh(ig)/4.
     .        +(dv(i,j,k,2)/jac+dv(i,j,km,2)/jackm)
     .          *(pp0(i,jp,km,1)-pp0(i,jm,km,1)
     .           +pp0(i,jp,k ,1)-pp0(i,jm,k ,1))/g_def%dyh(jg)/4.
     .        +(dv(i,j,k,3)/jac+dv(i,j,km,3)/jackm)
     .          *(pp0(i,j,k ,1)-pp0(i,j,km,1))/g_def%dz(kg-1) )
     .        *0.5

c      Fluxes at faces for calculation of grad(gamma*p0*div(dv))

      !!Divergence at faces i+-1/2, etc.
      divip = (dv(ip,j ,k,1)-dv(i ,j ,k,1))/g_def%dx(ig)
     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .        +dv(ip,jp,k,2)-dv(ip,jm,k,2))/g_def%dyh(jg)*0.25
     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .        +dv(ip,j,kp,3)-dv(ip,j,km,3))/g_def%dzh(kg)*0.25
      divim = (dv(i ,j ,k,1)-dv(im,j ,k,1))/g_def%dx(ig-1)
     .       +(dv(i ,jp,k,2)-dv(i ,jm,k,2)
     .        +dv(im,jp,k,2)-dv(im,jm,k,2))/g_def%dyh(jg)*0.25
     .       +(dv(i ,j,kp,3)-dv(i ,j,km,3)
     .        +dv(im,j,kp,3)-dv(im,j,km,3))/g_def%dzh(kg)*0.25

      divjp = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .        +dv(ip,jp,k,1)-dv(im,jp,k,1))/g_def%dxh(ig)*0.25
     .       +(dv(i ,jp,k,2)-dv(i ,j ,k,2))/g_def%dy(jg)
     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .        +dv(i,jp,kp,3)-dv(i,jp,km,3))/g_def%dzh(kg)*0.25
      divjm = (dv(ip,j ,k,1)-dv(im,j ,k,1)
     .        +dv(ip,jm,k,1)-dv(im,jm,k,1))/g_def%dxh(ig)*0.25
     .       +(dv(i ,j ,k,2)-dv(i ,jm,k,2))/g_def%dy(jg-1)
     .       +(dv(i,j ,kp,3)-dv(i,j ,km,3)
     .        +dv(i,jm,kp,3)-dv(i,jm,km,3))/g_def%dzh(kg)*0.25

      divkp = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .        +dv(ip,j,kp,1)-dv(im,j,kp,1))/g_def%dxh(ig)*0.25
     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .        +dv(i,jp,kp,2)-dv(i,jm,kp,2))/g_def%dyh(jg)*0.25
     .       +(dv(i,j ,kp,3)-dv(i,j ,k ,3))/g_def%dz(kg)
      divkm = (dv(ip,j,k ,1)-dv(im,j,k ,1)
     .        +dv(ip,j,km,1)-dv(im,j,km,1))/g_def%dxh(ig)*0.25
     .       +(dv(i,jp,k ,2)-dv(i,jm,k ,2)
     .        +dv(i,jp,km,2)-dv(i,jm,km,2))/g_def%dyh(jg)*0.25
     .       +(dv(i,j ,k ,3)-dv(i,j ,km,3))/g_def%dz(kg-1)

      gamm = gamma

      flxip = flxip
     .      + gamm*(pp0(i,j,k,1)+pp0(ip,j,k,1))*divip/(jac+jacip)
      if (.not.isSP(g_def,i,j,k,igr,igr,igr)) then
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

      cov(1) = (flxip - flxim)/g_def%dxh(ig)
      cov(2) = (flxjp - flxjm)/g_def%dyh(jg)
      cov(3) = (flxkp - flxkm)/g_def%dzh(kg)

      call transformFromCurvToCurv(g_def,i,j,k,igr,igr,igr
     .                            ,cov(1),cov(2),cov(3)
     .                            ,cnv(1),cnv(2),cnv(3),is_cnv)

      psit = -alpha**2*cnv

      end function findPsit

c     findPsib
c     #############################################################
      function findPsib(g_def,i,j,k,igr) result(psib)

      implicit none

c     Call variables

      integer :: i,j,k,igr

      type(grid_mg_def),pointer :: g_def

c     Local variables

      integer :: ig,jg,kg,nxx,nyy,nzz

      real(8) :: coeff,jac,dxx,dyy,dzz,ijac,psib(3)
     .          ,j0_cnv(3),acnv(3)
     .          ,acnvip(3),acnvim(3)
     .          ,acnvjp(3),acnvjm(3)
     .          ,acnvkp(3),acnvkm(3)
     .          ,acovip(3),acovim(3)
     .          ,acovjp(3),acovjm(3)
     .          ,acovkp(3),acovkm(3)

c     Begin program

      nxx = g_def%nxv(igr)
      nyy = g_def%nyv(igr)
      nzz = g_def%nzv(igr)

      if (nc_eom_jxb.and.(.not.bcSP()) ) then  !Force conservative form in SP geometries

        call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

        jac = g_def%gmetric%grid(igr)%jac(i,j,k)
        ijac = 1d0/jac

        dxx = g_def%dxh(ig)
        dyy = g_def%dyh(jg)
        dzz = g_def%dzh(kg)

        j0_cnv = curl(g_def,i,j,k,nxx,nyy,nzz,igr,igr,igr,b0_cov)

        acnv =
     .      U_Bv_ijk(g_def,i  ,j,k,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov
     .              ,a0_cnv,0)

        acnvip =
     .      U_Bv_ijk(g_def,i  ,j,k,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov
     .              ,a0_cnv,1)
        acnvim =
     .      U_Bv_ijk(g_def,i-1,j,k,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov
     .              ,a0_cnv,1)

        acnvjp =
     .      U_Bv_ijk(g_def,i,j  ,k,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov
     .              ,a0_cnv,2)
        acnvjm =
     .      U_Bv_ijk(g_def,i,j-1,k,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov
     .              ,a0_cnv,2)

        acnvkp =
     .      U_Bv_ijk(g_def,i,j,k  ,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov
     .              ,a0_cnv,3)
        acnvkm =
     .      U_Bv_ijk(g_def,i,j,k-1,nxx,nyy,nzz,igr,igr,igr,dv,dv_cov
     .              ,a0_cnv,3)

        call transformFromCurvToCurv(g_def,i,j,k,igr,igr,igr
     .                        ,acovip(1),acovip(2),acovip(3)
     .                        ,acnvip(1),acnvip(2),acnvip(3)
     .                        ,.false.,half_elem=1)

        call transformFromCurvToCurv(g_def,i-1,j,k,igr,igr,igr
     .                        ,acovim(1),acovim(2),acovim(3)
     .                        ,acnvim(1),acnvim(2),acnvim(3)
     .                        ,.false.,half_elem=1)

        call transformFromCurvToCurv(g_def,i,j,k,igr,igr,igr
     .                        ,acovjp(1),acovjp(2),acovjp(3)
     .                        ,acnvjp(1),acnvjp(2),acnvjp(3)
     .                        ,.false.,half_elem=2)

        call transformFromCurvToCurv(g_def,i,j-1,k,igr,igr,igr
     .                        ,acovjm(1),acovjm(2),acovjm(3)
     .                        ,acnvjm(1),acnvjm(2),acnvjm(3)
     .                        ,.false.,half_elem=2)

        call transformFromCurvToCurv(g_def,i,j,k,igr,igr,igr
     .                        ,acovkp(1),acovkp(2),acovkp(3)
     .                        ,acnvkp(1),acnvkp(2),acnvkp(3)
     .                        ,.false.,half_elem=3)

        call transformFromCurvToCurv(g_def,i,j,k-1,igr,igr,igr
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

        call transformFromCurvToCurv(g_def,i,j,k,igr,igr,igr
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),is_cnv)
      else
        cnv = div_tensor(g_def,i,j,k,nxx,nyy,nzz,igr,igr,igr,alt__eom()
     .                  ,tensor_x,tensor_y,tensor_z,vol=.false.)
      endif

      psib = alpha*cnv

      end function findPsib

c     tensor_x
c     #############################################################
      subroutine tensor_x(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ip,igrid
        real(8)    :: jac,ijac,gsuper(3,3),vol
     .               ,acnv(3),acnvp(3),b0cov(3),b0cnv(3),scalar_prod

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(g_def%gmetric%grid(igrid)%jac (ip,j,k)
     .               +g_def%gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +g_def%gmetric%grid(igrid)%gsup(i ,j,k,:,:))

        if (isSP(g_def,i+1,j,k,igx,igy,igz).and.flag /= 0) jac = SP_flsv

        if (flag /= 0) then
          acnv = U_Bv_ijk(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,dv
     .                   ,dv_cov,a0_cnv,1)
        else
          acnv = U_Bv_ijk(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,dv
     .                   ,dv_cov,a0_cnv,0)
        endif

        b0cnv = 0.5*(b0_cnv(ip,j,k,:)+b0_cnv(i,j,k,:))
c$$$        b0cov = 0.5*(b0_cov(ip,j,k,:)+b0_cov(i,j,k,:))

        if (flag /= 0) then
          call transformFromCurvToCurv(g_def,i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=1)
        else
          call transformFromCurvToCurv(g_def,i,j,k,igx,igy,igz
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
          if (isSP(g_def,i+1,j,k,igx,igy,igz)) then
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
      subroutine tensor_y(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EM SI operator
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

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

        jac    = 0.5*(g_def%gmetric%grid(igrid)%jac (i,jp,k)
     .               +g_def%gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +g_def%gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          acnv = U_Bv_ijk(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,dv
     .                   ,dv_cov,a0_cnv,2)
        else
          acnv = U_Bv_ijk(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,dv
     .                   ,dv_cov,a0_cnv,0)
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
          ijac  = 1d0/jac
          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine tensor_y

c     tensor_z
c     #############################################################
      subroutine tensor_z(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EM SI operator
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

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

        jac    = 0.5*(g_def%gmetric%grid(igrid)%jac (i,j,kp)
     .               +g_def%gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +g_def%gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          acnv = U_Bv_ijk(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,dv
     .                   ,dv_cov,a0_cnv,3)
        else
          acnv = U_Bv_ijk(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,dv
     .                   ,dv_cov,a0_cnv,0)
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
          ijac  = 1d0/jac
          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine tensor_z

      end module v_mtvc_mod

c v_mtvc
c####################################################################
      subroutine v_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer :: iig,ieq,ijk,ijkg,i,j,k,ig,jg,kg
     .          ,imin,imax,jmin,jmax,kmin,kmax
     .          ,nxx,nyy,nzz

      real(8) :: b2,bbcnv(3),bbcov(3),etal,E00(3),diag

      real(8),allocatable,dimension(:,:,:,:) :: v0,dpp

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      is_cnv = .true.

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig  = g_def%istartp(igrid)
      iisig = mg_ctx_l%istart(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Define pointers to MG arrays

      call defineLocalPointers(igrid,nxx,nyy,nzz)

c ExB flow BC (if applicable)

      allocate(v0(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      v0 = 0d0

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dpp   (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call fillArray(g_def,igrid,neq,x,dv,bcnd,1,gpos=gpos
     .              ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .              ,result_is_curv=.true.,arr_cov=dv_cov
     .              ,arr0=v0)

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jacip  = g_def%gmetric%grid(igrid)%jac(ip,j,k)
            jacim  = g_def%gmetric%grid(igrid)%jac(im,j,k)
            jacjp  = g_def%gmetric%grid(igrid)%jac(i,jp,k)
            jacjm  = g_def%gmetric%grid(igrid)%jac(i,jm,k)
            jackp  = g_def%gmetric%grid(igrid)%jac(i,j,kp)
            jackm  = g_def%gmetric%grid(igrid)%jac(i,j,km)
            jac    = g_def%gmetric%grid(igrid)%jac(i,j,k )

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (MG_get_res_vol_wgt()) then
              vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************
            psiv = D_v(g_def,i,j,k,nxx,nyy,nzz,igrid,rho0,v0_cnv,dv,dpp
     .                ,.not.nc_eom_v)

            !P_si^T  ******************************
            psit = findPsit(g_def,i,j,k,igrid)

            !P_si^B  ******************************
            psib = findPsib(g_def,i,j,k,igrid)

            !Diagonal scaling
            if (pc_type == 's2') then
              diag = cnp/dt
     .             + alpha*(abs(v0_cnv(i,j,k,1))/g_def%dx(ig)
     .                     +abs(v0_cnv(i,j,k,2))/g_def%dy(jg)
     .                     +abs(v0_cnv(i,j,k,3))/g_def%dz(kg))
     .                     /g_def%gmetric%grid(igrid)%jac(i,j,k)
c$$$     .             + alpha*max(dd,chi)
c$$$     .                   *((g_def%gmetric%grid(igrid)%gsup(i,j ,k ,1,1)
c$$$     .                     +g_def%gmetric%grid(igrid)%gsup(i,jm,k ,1,1))
c$$$     .                     /g_def%dx(ig)**2
c$$$     .                    +(g_def%gmetric%grid(igrid)%gsup(i ,j,k ,2,2)
c$$$     .                     +g_def%gmetric%grid(igrid)%gsup(im,j,k ,2,2))
c$$$     .                     /g_def%dy(jg)**2
c$$$     .                    +(g_def%gmetric%grid(igrid)%gsup(i ,j,k ,3,3)
c$$$     .                     +g_def%gmetric%grid(igrid)%gsup(im,j,k ,3,3))
c$$$     .                     /g_def%dz(kg)**2
c$$$     .                     )/g_def%gmetric%grid(igrid)%jac(i,j,k)
c$$$              cnv = cnv*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol
c$$$              psit = psit*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol
c$$$              psib = psib*max(rho_diag(1,ijkg),tmp_diag(1,ijkg))*vol
            else
              diag = cnp/dt
            endif

            cnv = (psit + psib)/diag

            !Form total operator

            if (nc_eom_v) then
              cnv = (psiv + cnv/rho0(i,j,k,1))*vol
            else
              cnv = (psiv + cnv)*vol
            endif

            if (alt__eom()) cnv(2) = cnv(2)*jac

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = cnv(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(dv,dv_cov,v0,dpp)

      call nullLocalPointers

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine v_mtvc

c v_mtvc_hyp
c####################################################################
      subroutine v_mtvc_hyp(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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

      type(grid_mg_def),pointer :: g_def

      integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8) :: x(ntot),y(ntot),minv,idx,idy,idz
     .          ,idxcnv,idycnv,idzcnv,vxx,vyy,vzz

c Local variables

      integer :: iig,ieq,ijk,ijkg,i,j,k,ig,jg,kg
     .          ,imin,imax,jmin,jmax,kmin,kmax
     .          ,nxx,nyy,nzz

      real(8) :: b2,bbcnv(3),bbcov(3),etal,E00(3),v_res(3),divpe(3)
     .          ,dj(3),aip(3),aim(3),ajp(3),ajm(3),akp(3),akm(3),diag

      real(8),pointer,dimension(:,:,:,:) :: v0,dE_cnv,du,du_cov,dpp
     .                                     ,dve_cnv,dB_cnv,dB_cov
     .                                     ,heta0,dj_cnv

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      if (neq /= 6) call pstop('v_mtvc_hyp','neq /= 6; aborting...')

      is_cnv = .true.

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig  = g_def%istartp(igrid)
      iisig = mg_ctx_l%istart(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Define pointers to MG arrays

      call defineLocalPointers(igrid,nxx,nyy,nzz)

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(du    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,du_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      allocate(dve_cnv(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dpp    (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dE_cnv (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dB_cnv (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dB_cov (0:nxx+1,0:nyy+1,0:nzz+1,3))

c BCs

      !ExB flow BC (if applicable)
      allocate(v0(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      v0 = 0d0

      du = 0d0

      call set_force_corners(.true.) !Force corners owing to div(A)=0 BC

      call fillArray(g_def,igrid,neq,x,du,bcnd,1,gpos=gpos
cc     .              ,icomp=(/IVX,IJX/),is_cnv=is_cnv,is_vec=.true.
     .              ,icomp=(/IVX,IAX/),is_cnv=is_cnv,is_vec=.true.
     .              ,result_is_curv=.true.,arr_cov=du_cov,arr0=v0)

      call set_force_corners(.false.) !Force corners owing to div(A)=0 BC

cc#if defined(pgf90)
      allocate(dv    (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,dj_cnv(0:nxx+1,0:nyy+1,0:nzz+1,3))

      dv    (0:,0:,0:,1:3)  = du    (:,:,:,1:3)   !Initialize module variable
      dv_cov(0:,0:,0:,1:3)  = du_cov(:,:,:,1:3)   !Initialize module variable
      dj_cnv(0:,0:,0:,1:3)  = du    (:,:,:,4:6)
c$$$#else
c$$$      dv    (0:,0:,0:,1:)  => du    (:,:,:,1:3)
c$$$      dv_cov(0:,0:,0:,1:)  => du_cov(:,:,:,1:3)
c$$$      dj_cnv(0:,0:,0:,1:)  => du    (:,:,:,4:6)
c$$$#endif

      dve_cnv = form_ve(dv,dj_cnv,rho0(:,:,:,1))

c Calculate matrix-vector product

      call point_to_MG_array(nxx,nyy,nzz,igrid,heta0,gheta)

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jacip  = g_def%gmetric%grid(igrid)%jac(ip,j,k)
            jacim  = g_def%gmetric%grid(igrid)%jac(im,j,k)
            jacjp  = g_def%gmetric%grid(igrid)%jac(i,jp,k)
            jacjm  = g_def%gmetric%grid(igrid)%jac(i,jm,k)
            jackp  = g_def%gmetric%grid(igrid)%jac(i,j,kp)
            jackm  = g_def%gmetric%grid(igrid)%jac(i,j,km)
            jac    = g_def%gmetric%grid(igrid)%jac(i,j,k )

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (MG_get_res_vol_wgt()) then
              vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************

            psiv = D_v(g_def,i,j,k,nxx,nyy,nzz,igrid,rho0,v0_cnv,dv,dpp
     .                ,.not.nc_eom_v)

            !Electron viscosity piece *************
            divpe = EOM_divPi_ijk(g_def,i,j,k,nxx,nyy,nzz,igrid
     .                           ,dve_cnv,heta0,.false.)

            if (nc_eom_v) then
              psiv = psiv + alpha*divpe/rho0(i,j,k,1)
            else
              psiv = psiv + alpha*divpe
            endif

            !P_si^T  ******************************

            psit = findPsit(g_def,i,j,k,igrid)

            !P_si^B  ******************************

            psib = findPsib(g_def,i,j,k,igrid)

            !Diagonal scaling
            cnv = psit + psib 

c$$$            if (pc_type == 's2') then
              diag = cnp/dt
     .             + alpha*(abs(v0_cnv(i,j,k,1))/g_def%dx(ig)
     .                     +abs(v0_cnv(i,j,k,2))/g_def%dy(jg)
     .                     +abs(v0_cnv(i,j,k,3))/g_def%dz(kg))
     .                     /g_def%gmetric%grid(igrid)%jac(i,j,k)
c$$$     .             + alpha*max(dd,chi)
c$$$     .                   *((g_def%gmetric%grid(igrid)%gsup(i,j ,k ,1,1)
c$$$     .                     +g_def%gmetric%grid(igrid)%gsup(i,jm,k ,1,1))
c$$$     .                     /g_def%dx(ig)**2
c$$$     .                    +(g_def%gmetric%grid(igrid)%gsup(i ,j,k ,2,2)
c$$$     .                     +g_def%gmetric%grid(igrid)%gsup(im,j,k ,2,2))
c$$$     .                     /g_def%dy(jg)**2
c$$$     .                    +(g_def%gmetric%grid(igrid)%gsup(i ,j,k ,3,3)
c$$$     .                     +g_def%gmetric%grid(igrid)%gsup(i ,j,km,3,3))
c$$$     .                     /g_def%dz(kg)**2
c$$$     .                     )
c$$$     .                     /g_def%gmetric%grid(igrid)%jac(i,j,k)
c$$$            else
c$$$              diag = cnp/dt
c$$$            endif

            cnv = cnv/diag

            !Form total operator
            if (nc_eom_v) then
              v_res = (psiv + cnv/rho0(i,j,k,1))
            else
              v_res = (psiv + cnv)
            endif

            !Find current update, dj=-dt*curl(curl(dE))
cc            !Find dj=curl(db)
cc            dj=curl(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dB_cov)

            !Find dj=-dt*curl(U_Bv_ijk.dv)
            cnv = U_Bv_ijk(g_def,i ,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                    ,dv,dv_cov,a0_cnv,1)
            call transformFromCurvToCurv(g_def,i,j,k
     .                                  ,igrid,igrid,igrid
     .                                  ,aip(1),aip(2),aip(3)
     .                                  ,cnv(1),cnv(2),cnv(3)
     .                                  ,.false.,half_elem=1)

            cnv = U_Bv_ijk(g_def,im,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                    ,dv,dv_cov,a0_cnv,1)
            call transformFromCurvToCurv(g_def,im,j,k
     .                                  ,igrid,igrid,igrid
     .                                  ,aim(1),aim(2),aim(3)
     .                                  ,cnv(1),cnv(2),cnv(3)
     .                                  ,.false.,half_elem=1)
                             
            cnv = U_Bv_ijk(g_def,i,j ,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                    ,dv,dv_cov,a0_cnv,2)
            call transformFromCurvToCurv(g_def,i,j,k,igrid,igrid,igrid
     .                                  ,ajp(1),ajp(2),ajp(3)
     .                                  ,cnv(1),cnv(2),cnv(3)
     .                                  ,.false.,half_elem=2)
                             
            cnv = U_Bv_ijk(g_def,i,jm,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                    ,dv,dv_cov,a0_cnv,2)
            call transformFromCurvToCurv(g_def,i,jm,k
     .                                  ,igrid,igrid,igrid
     .                                  ,ajm(1),ajm(2),ajm(3)
     .                                  ,cnv(1),cnv(2),cnv(3)
     .                                  ,.false.,half_elem=2)
                             
            cnv = U_Bv_ijk(g_def,i,j,k ,nxx,nyy,nzz,igrid,igrid,igrid
     .                    ,dv,dv_cov,a0_cnv,3)
            call transformFromCurvToCurv(g_def,i,j,k,igrid,igrid,igrid
     .                                  ,akp(1),akp(2),akp(3)
     .                                  ,cnv(1),cnv(2),cnv(3)
     .                                  ,.false.,half_elem=3)
                             
            cnv = U_Bv_ijk(g_def,i,j,km,nxx,nyy,nzz,igrid,igrid,igrid
     .                    ,dv,dv_cov,a0_cnv,3)
            call transformFromCurvToCurv(g_def,i,j,km
     .                                  ,igrid,igrid,igrid
     .                                  ,akm(1),akm(2),akm(3)
     .                                  ,cnv(1),cnv(2),cnv(3)
     .                                  ,.false.,half_elem=3)

            cnv(1) = (ajp(3)-ajm(3))/g_def%dyh(jg)
     .              -(akp(2)-akm(2))/g_def%dzh(kg)

            cnv(2) = (akp(1)-akm(1))/g_def%dzh(kg)
     .              -(aip(3)-aim(3))/g_def%dxh(ig)

            cnv(3) = (aip(2)-aim(2))/g_def%dxh(ig)
     .              -(ajp(1)-ajm(1))/g_def%dyh(jg)

            dj = -dt*cnv/cnp

            if (alt__eom()) v_res(2) = v_res(2)*jac

            !Assemble matvec
            do ieq=1,3
              y(neq*(ijk-1)+ieq)   = vol*v_res(ieq)
cc              y(neq*(ijk-1)+ieq+3) = vol*(du(i,j,k,ieq+3)-dj(ieq))
              y(neq*(ijk-1)+ieq+3) = vol*(dj_cnv(i,j,k,ieq)-dj(ieq))
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(du,du_cov,dpp,v0,dE_cnv,dve_cnv
     .          ,dB_cnv,dB_cov)

cc#if defined(pgf90)
      deallocate(dv,dv_cov,dj_cnv)
c$$$#else
c$$$      nullify   (dv,dv_cov,dj_cnv)
c$$$#endif

      call nullLocalPointers

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine v_mtvc_hyp

c v_mtvc2
c####################################################################
      subroutine v_mtvc2(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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
      type(grid_mg_def),pointer :: g_def

c Local variables

      integer :: iig,ieq,ijk,ijkg,i,j,k,ig,jg,kg
     .          ,imin,imax,jmin,jmax,kmin,kmax
     .          ,nxx,nyy,nzz

      real(8),allocatable,dimension(:,:,:,:) :: dpp

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      is_cnv = .true.

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig  = g_def%istartp(igrid)
      iisig = mg_ctx_l%istart(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Define pointers to MG arrays

      call point_to_MG_array(nxx,nyy,nzz,igrid,rho0  ,grho0)
      call point_to_MG_array(nxx,nyy,nzz,igrid,pp0   ,gp0)
      call point_to_MG_array(nxx,nyy,nzz,igrid,v0_cnv,gv0)
      call point_to_MG_array(nxx,nyy,nzz,igrid,b0_cnv,gb0)
      call point_to_MG_array(nxx,nyy,nzz,igrid,b0_cov,gb0_cov)
      call point_to_MG_array(nxx,nyy,nzz,igrid,a0_cnv,ga0)
 
c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dpp   (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call fillArray(g_def,igrid,neq,x,dv,bcnd,1,gpos=gpos
     .              ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .              ,result_is_curv=.true.,arr_cov=dv_cov)

cc      call mapMGVectorToArray(mg_ctx_temp,max(0,gpos),x,dv,igrid
cc     .                       ,.false.)
cc
cccc      if (alt__eom()) dv(:,:,:,2) = dv(:,:,:,2)/g_def%gmetric%grid(igrid)%jac
cc
cc      if (di > 0) then
cc        call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
cccc     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.not.si_car
cc     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
cc     .            ,result_is_curv=.true.,iorder=1,arr_cov=dv_cov)
cc      else
cc        call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dv,bcnd
cccc     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.not.si_car
cc     .            ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
cc     .            ,result_is_curv=.true.,iorder=order_bc,arr_cov=dv_cov)
cc      endif

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jacip  = g_def%gmetric%grid(igrid)%jac(ip,j,k)
            jacim  = g_def%gmetric%grid(igrid)%jac(im,j,k)
            jacjp  = g_def%gmetric%grid(igrid)%jac(i,jp,k)
            jacjm  = g_def%gmetric%grid(igrid)%jac(i,jm,k)
            jackp  = g_def%gmetric%grid(igrid)%jac(i,j,kp)
            jackm  = g_def%gmetric%grid(igrid)%jac(i,j,km)
            jac    = g_def%gmetric%grid(igrid)%jac(i,j,k )

cc            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + SP_flsv

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (MG_get_res_vol_wgt()) then
              vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************

            psiv = D_v(g_def,i,j,k,nxx,nyy,nzz,igrid,rho0,v0_cnv,dv
     .                ,dpp,.not.nc_eom_v)

            cnv = psiv*vol

            !Form total operator
            if (alt__eom()) cnv(2) = cnv(2)*jac

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = cnv(ieq)
            enddo

          enddo
        enddo
      enddo

c End program

      deallocate(dpp,dv,dv_cov)

      nullify(pp0,rho0)
      nullify(b0_cnv,b0_cov)
      nullify(v0_cnv,a0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

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
cc      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)
cc
cc      if (igrid > 1) then
cc
cc        isig  = istart(igrid)
cc
ccc     Allocate MG vector
cc
cc        igr = 1
cc
cc        nxx = g_def%nxv(igr)
cc        nyy = g_def%nyv(igr)
cc        nzz = g_def%nzv(igr)
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
cc        call restrictMGVector(neq,MGy,igr,igrid,0,MG_get_res_vol_wgt())
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
cc      call deallocMGPointers(fpointers,mg_ctx_l)
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
cc      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)
cc
cc      isig  = g_def%istartp(igrid)
cc      iisig = mg_ctx_l%istart(igrid)
cc
cc      nxx = g_def%nxv(igrid)
cc      nyy = g_def%nyv(igrid)
cc      nzz = g_def%nzv(igrid)
cc
cc      igx = igrid
cc      igy = igrid
cc      igz = igrid
cc
ccc Find limits for loops
cc
cc      call limits(mg_ctx,abs(gpos),nxx,nyy,nzz,igrid
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
cc            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)
cc
cc            jacip  = g_def%gmetric%grid(igrid)%jac(ip,j,k)
cc            jacim  = g_def%gmetric%grid(igrid)%jac(im,j,k)
cc            jacjp  = g_def%gmetric%grid(igrid)%jac(i,jp,k)
cc            jacjm  = g_def%gmetric%grid(igrid)%jac(i,jm,k)
cc            jackp  = g_def%gmetric%grid(igrid)%jac(i,j,kp)
cc            jackm  = g_def%gmetric%grid(igrid)%jac(i,j,km)
cc            jac    = g_def%gmetric%grid(igrid)%jac(i,j,k )
cc
cccc            if (isSP(i,j,k,igrid,igrid,igrid)) jacim = jacim + SP_flsv
cc
cc            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)
cc
cc            ijkg   = ijk + isig - 1
cc
cc            if (MG_get_res_vol_wgt()) then
cc              vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
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
cc      call deallocMGPointers(fpointers,mg_ctx_l)
cc
cc      end subroutine v_mtvc2

c v_mtvc_af
c####################################################################
      subroutine v_mtvc_af(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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

      type(grid_mg_def),pointer :: g_def

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz

c Local variables

      integer :: iig,ieq,ijk,ijkg,i,j,k,ig,jg,kg
     .          ,imin,imax,jmin,jmax,kmin,kmax
     .          ,nxx,nyy,nzz

      real(8) :: z(ntot)

      real(8),allocatable,dimension(:,:,:,:) :: dzz,dzp

      type(mg_ctx),pointer :: mg_ctx_l => null()

      INTERFACE
        subroutine mstar_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
          use grid
          type(grid_mg_def),pointer :: g_def
          integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
          real(8) :: x(ntot),y(ntot)
        end subroutine mstar_mtvc
      END INTERFACE

c Begin program

      is_cnv = .true.

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig  = g_def%istartp(igrid)
      iisig = mg_ctx_l%istart(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Define pointers to MG arrays

      call defineLocalPointers(igrid,nxx,nyy,nzz)

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dv    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dzz   (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dzp   (0:nxx+1,0:nyy+1,0:nzz+1,neq))

      !Velocity array
      dv = 0d0

      call fillArray(g_def,igrid,neq,x,dv,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .              ,result_is_curv=.true.,arr_cov=dv_cov)

      !M_star*dv array
      dzz = 0d0

      call mstar_mtvc(g_def,gpos,neq,ntot,x,z,igrid,bcnd)

      call fillArray(g_def,igrid,neq,z,dzz,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .              ,result_is_curv=.true.)

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jacip  = g_def%gmetric%grid(igrid)%jac(ip,j,k)
            jacim  = g_def%gmetric%grid(igrid)%jac(im,j,k)
            jacjp  = g_def%gmetric%grid(igrid)%jac(i,jp,k)
            jacjm  = g_def%gmetric%grid(igrid)%jac(i,jm,k)
            jackp  = g_def%gmetric%grid(igrid)%jac(i,j,kp)
            jackm  = g_def%gmetric%grid(igrid)%jac(i,j,km)
            jac    = g_def%gmetric%grid(igrid)%jac(i,j,k )

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (MG_get_res_vol_wgt()) then
              vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              vol = 1d0
            endif

            !P_si^v  ******************************

            psiv = D_v(g_def,i,j,k,nxx,nyy,nzz,igrid,rho0,v0_cnv,dzz,dzp
     .                ,.not.nc_eom_v)

            !P_si^T  ******************************

            psit = findPsit(g_def,i,j,k,igrid)

            !P_si^B  ******************************

            psib = findPsib(g_def,i,j,k,igrid)

            !Form total operator

            cnv = psit + psib

            if (nc_eom_v) then
              cnv = (psiv + cnv/rho0(i,j,k,1))*vol
            else
              cnv = (psiv*rho0(i,j,k,1) + cnv)*vol
            endif

            if (alt__eom()) cnv(2) = cnv(2)*jac

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

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine v_mtvc_af

c mstar_mtvc
c####################################################################
      subroutine mstar_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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

      type(grid_mg_def),pointer :: g_def

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz
     .             ,idxcnv,idycnv,idzcnv,vxx,vyy,vzz
     .             ,diff,upwind

c Local variables

      integer :: iig,ieq,ijk,ijkg,i,j,k,ig,jg,kg
     .          ,imin,imax,jmin,jmax,kmin,kmax
     .          ,nxx,nyy,nzz

      real(8),allocatable,dimension(:,:,:,:) :: da

      logical :: flg_rho,flg_v

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program

      is_cnv = .true.

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig  = g_def%istartp(igrid)
      iisig = mg_ctx_l%istart(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Define pointers to MG arrays

      flg_rho = .false.      !To prevent interference with v_mtvc_af
      if (.not.associated(rho0  )) then
        flg_rho = .true.
        call point_to_MG_array(nxx,nyy,nzz,igrid,rho0  ,grho0)
      endif

      flg_v = .false.        !To prevent interference with v_mtvc_af
      if (.not.associated(v0_cnv)) then
        flg_v = .true.
        call point_to_MG_array(nxx,nyy,nzz,igrid,v0_cnv,gv0)
      endif

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(da(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      da = 0d0

      call fillArray(g_def,igrid,neq,x,da,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .              ,result_is_curv=.true.)

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            jacip  = g_def%gmetric%grid(igrid)%jac(ip,j,k)
            jacim  = g_def%gmetric%grid(igrid)%jac(im,j,k)
            jacjp  = g_def%gmetric%grid(igrid)%jac(i,jp,k)
            jacjm  = g_def%gmetric%grid(igrid)%jac(i,jm,k)
            jackp  = g_def%gmetric%grid(igrid)%jac(i,j,kp)
            jackm  = g_def%gmetric%grid(igrid)%jac(i,j,km)
            jac    = g_def%gmetric%grid(igrid)%jac(i,j,k )

            if (isSP(g_def,i,j,k,igrid,igrid,igrid))
     .           jacim = jacim + SP_flsv

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

cc            diff = max(eta,chi,dd)
cc            diff = min(eta,chi,dd)

cc            cnv = gradDiv(i,j,k,nxx,nyy,nzz,igx,igy,igz,da,vol=.false.)
cc            cnv = matmul(g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:),cnv)
cc            if (diff /= 0d0) cnv = gradDiv2(da)

            do ieq=1,neq
cc              hex = floor(sign(1d0,v0_cnv(i,j,k,ieq)))
              hex = 0

              upwind = v0_cnv(i,j,k,ieq)
     .                *div_upwd(g_def,i,j,k,igrid,da,hex)

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

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine mstar_mtvc

c hall_mtvc
c####################################################################
      subroutine hall_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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

      type(grid_mg_def),pointer :: g_def

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Local variables

      integer :: iig,ieq,ijk,ijkg,i,j,k,ig,jg,kg
     .          ,imin,imax,jmin,jmax,kmin,kmax
     .          ,nxx,nyy,nzz

c Begin program

      is_cnv = .true.

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig  = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Define pointers to MG arrays

      call point_to_MG_array(nxx,nyy,nzz,igrid,b0_cnv,gb0)

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dv    (0:nxx+1,0:nyy+1,0:nzz+1,neq)
     .        ,dv_cov(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call fillArray(g_def,igrid,neq,x,dv,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IVX/),is_cnv=is_cnv,is_vec=.true.
     .              ,result_is_curv=.true.,arr_cov=dv_cov)

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            if (MG_get_res_vol_wgt()) then
              vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              vol = 1d0
            endif

            cnv(1) = dv(i,j,k,1)/dt
cc     .          - di*b0_cnv(i,j,k,3)
     .          - di
     .           *((dv(ip,j,k,2)-2*dv(i,j,k,2)+dv(im,j,k,2))
     .             /g_def%dxh(ig)**2
     .            -(dv(ip,jp,k,1)-dv(ip,jm,k,1)
     .             -dv(im,jp,k,1)+dv(im,jm,k,1))
     .             /4./g_def%dxh(ig)/g_def%dyh(jg))

            cnv(2) = dv(i,j,k,2)/dt
cc     .          + di*b0_cnv(i,j,k,3)
     .          + di
     .           *((dv(i,jp,k,1)-2*dv(i,j,k,1)+dv(i,jm,k,1))
     .             /g_def%dyh(jg)**2
     .            -(dv(ip,jp,k,2)-dv(ip,jm,k,2)
     .             -dv(im,jp,k,2)+dv(im,jm,k,2))
     .             /4./g_def%dxh(ig)/g_def%dyh(jg))

            cnv(3) = dv(i,j,k,3)/dt

            cnv = cnv*vol

            if (alt__eom()) cnv(2) = cnv(2)*jac

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = cnv(ieq)
            enddo
          enddo
        enddo
      enddo

c End program

      deallocate(dv,dv_cov)

      nullify(b0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine hall_mtvc

c sm_mtvc_B
c####################################################################
      subroutine sm_mtvc_B(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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

      type(grid_mg_def),pointer :: g_def

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Local variables

      integer :: iig,ieq,ijk,ijkg,i,j,k,ig,jg,kg
     .          ,imin,imax,jmin,jmax,kmin,kmax
     .          ,nxx,nyy,nzz

      real(8) :: coeff

c Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig  = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Define pointers to MG arrays

cc      call point_to_MG_array(nxx,nyy,nzz,igrid,b0_cnv,gb0)

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dv(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

#if defined(vec_pot)
      call fillArray(g_def,igrid,neq,x,dv,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IAX/),is_cnv=.false.,is_vec=.true.
     .              ,result_is_curv=.true.)
#else
      call fillArray(g_def,igrid,neq,x,dv,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IBX/),is_cnv=.true.,is_vec=.true.
     .              ,result_is_curv=.true.)
#endif

c Calculate matrix-vector product

      coeff = 1d0
c$$$      do k = kmin,kmax
c$$$        do j = jmin,jmax
c$$$          do i = imin,imax
c$$$            coeff= max(coeff
c$$$     .     ,sqrt(vectorNorm(g_def,i,j,k,igrid,b0_cnv(i,j,k,:),.false.)))
c$$$          enddo
c$$$        enddo
c$$$      enddo

      coeff = di*coeff  !v_A*di

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)

            cnv(:) = (dv(ip,j,k,:)-2*dv(i,j,k,:)+dv(im,j,k,:))
     .              +(dv(i,jp,k,:)-2*dv(i,j,k,:)+dv(i,jm,k,:))
     .              +(dv(i,j,kp,:)-2*dv(i,j,k,:)+dv(i,j,km,:))

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = dv(i,j,k,ieq)*vol-dt*coeff*cnv(ieq)
            enddo
          enddo
        enddo
      enddo

c End program

      deallocate(dv)

cc      nullify(b0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine sm_mtvc_B
      
c sm_mtvc_V
c####################################################################
      subroutine sm_mtvc_V(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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

      type(grid_mg_def),pointer :: g_def

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot),minv,idx,idy,idz,
     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Local variables

      integer :: iig,ieq,ijk,ijkg,i,j,k,ig,jg,kg
     .          ,imin,imax,jmin,jmax,kmin,kmax
     .          ,nxx,nyy,nzz

      real(8) :: coeff

c Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig  = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Define pointers to MG arrays

cc      call point_to_MG_array(nxx,nyy,nzz,igrid,b0_cnv,gb0)

c Map vector x to array for processing (return dv in curvilinear comp.)

      allocate(dv(0:nxx+1,0:nyy+1,0:nzz+1,neq))
      
      dv = 0d0

      call fillArray(g_def,igrid,neq,x,dv,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IVX/),is_cnv=.true.,is_vec=.true.
     .              ,result_is_curv=.true.)

c Calculate matrix-vector product

      coeff = 1d0
c$$$      do k = kmin,kmax
c$$$        do j = jmin,jmax
c$$$          do i = imin,imax
c$$$            coeff= max(coeff
c$$$     .     ,sqrt(vectorNorm(g_def,i,j,k,igrid,b0_cnv(i,j,k,:),.false.)))
c$$$          enddo
c$$$        enddo
c$$$      enddo

      coeff = di*coeff  !v_A*di

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

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg   = ijk + isig - 1

            vol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)

            cnv(:) = (dv(ip,j,k,:)-2*dv(i,j,k,:)+dv(im,j,k,:))
     .              +(dv(i,jp,k,:)-2*dv(i,j,k,:)+dv(i,jm,k,:))
     .              +(dv(i,j,kp,:)-2*dv(i,j,k,:)+dv(i,j,km,:))

            do ieq=1,3
              y(neq*(ijk-1)+ieq) = dv(i,j,k,ieq)*vol-dt*coeff*cnv(ieq)
            enddo
          enddo
        enddo
      enddo

c End program

      deallocate(dv)

cc      nullify(b0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine sm_mtvc_V
