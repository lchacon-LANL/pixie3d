c setupPreconditioner
c####################################################################
      subroutine setupPreconditioner(ntot,x)

c--------------------------------------------------------------------
c     Finds velocity and magnetic field components in all grids
c     for matrix-free preconditioning.
c--------------------------------------------------------------------

      use grid

      use precond_variables

      use matvec

      use operators

      use newton_gmres

      implicit none

c Call variables

      integer(4) :: ntot
      real(8)    :: x(ntot)

c Local variables

      integer(4) :: ntotd2p,order,nxx,nyy,nzz,igrid,alloc_stat
      real(8)    :: dvol

      real(8), allocatable, dimension(:,:,:,:) :: vector

c Debug

      real(8)    :: mag,debug(0:nx+1,0:ny+1,0:nz+1)
      integer(4) :: nt,ii
      real(8),allocatable,dimension(:,:) :: v_mat,debug2

c Externals

      external   :: v_mtvc,b_mtvc,rho_mtvc,tmp_mtvc

c Begin program

      if (precon == 'id') return

      igx = 1
      igy = 1
      igz = 1

      nxx = grid_params%nxv(igx)
      nyy = grid_params%nyv(igy)
      nzz = grid_params%nzv(igz)

      nx = nxx
      ny = nyy
      nz = nzz

      ngrid = grid_params%ngrid

      igrid = igx

      ntotd2p = 2*ntot/neqd

      order=0

c Allocate variables

      allocate (mgj0cnv(ntotd2p,3))
      allocate (mgadvdiffV0(ntotd2p,3))
      allocate (mgdivV0(ntotd2p))
      allocate (bcs(6,neqd+3))

      if (jit == 1) then  !Only allocate at the beginning of Newton iteration
        deallocate (rho_diag,tmp_diag,b_diag,v_diag,STAT=alloc_stat)

        allocate (rho_diag(1,ntotd2p)
     .           ,tmp_diag(1,ntotd2p)
     .           ,  b_diag(3,3*ntotd2p)
     .           ,  v_diag(3,3*ntotd2p),STAT=alloc_stat)
      endif

      call allocateMGArray(1,gp0)
      call allocateMGArray(1,grho0)
      call allocateMGArray(3,gv0)
      call allocateMGArray(3,gb0)

c Unpack vector x, taking BCs from t=n solution

      varray = x  !This allocates varray; overloaded assignment

      call imposeBoundaryConditions(varray,igx,igy,igz)

c Extract arrays and BC's

      rho => varray%array_var(IRHO)%array
      rvx => varray%array_var(IVX )%array
      rvy => varray%array_var(IVY )%array
      rvz => varray%array_var(IVZ )%array
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
      tmp => varray%array_var(ITMP)%array

      do ieq=1,neqd
        bcs(:,ieq) = varray%array_var(ieq)%bconds(:)
      enddo
      !Current boundary conditions
      bcs(:,IJX:IJZ) = bcs(:,IBX:IBZ)
      where (bcs(:,IJX:IJZ) == -NEU)
        bcs(:,IJX:IJZ) = -DIR  !Use covariant components for tangential dirichlet
      end where

c Store density in all grids (w/o BCs)

cc      call restrictArrayToMGVector(1,nxx,nyy,nzz,rho,mgrho0,igrid,order
cc     .                            ,.false.)

      grho0%grid(igrid)%array(:,:,:,1) = rho

      call restrictMGArray(IRHO,1,grho0,bcs(:,IRHO),igrid,order)

c Store magnetic field components in all grids (w/ BCs)

      gb0%grid(igrid)%array(:,:,:,1) = bx
      gb0%grid(igrid)%array(:,:,:,2) = by
      gb0%grid(igrid)%array(:,:,:,3) = bz

      call restrictMGArray(IBX,3,gb0,bcs(:,IBX:IBZ),igrid,order)

c Store ion velocity components in all grids (w/ BCs)

      gv0%grid(igrid)%array(:,:,:,1) = vx
      gv0%grid(igrid)%array(:,:,:,2) = vy
      gv0%grid(igrid)%array(:,:,:,3) = vz

      call restrictMGArray(IVX,3,gv0,bcs(:,IVX:IVZ),igrid,order)

c Store current components in all grids (w/o BCs)

      call restrictArrayToMGVector(1,nxx,nyy,nzz,jx,mgj0cnv(:,1),igrid
     .                  ,order,.false.)
      call restrictArrayToMGVector(1,nxx,nyy,nzz,jy,mgj0cnv(:,2),igrid
     .                  ,order,.false.)
      call restrictArrayToMGVector(1,nxx,nyy,nzz,jz,mgj0cnv(:,3),igrid
     .                  ,order,.false.)

c Find auxiliary quantities and store them in all grids

      !Velocity divergence (w/o BCs)
      allocate(divrgV(0:nxx+1,0:nyy+1,0:nzz+1))

      do k=1,nzz
        do j=1,nyy
          do i=1,nxx
            divrgV(i,j,k) = div(i,j,k,nxx,nyy,nzz,vx,vy,vz)
          enddo
        enddo
      enddo

      call restrictArrayToMGVector(1,nxx,nyy,nzz,divrgV,mgdivV0,igrid
     .                            ,order,.false.)

      deallocate(divrgV)

      !pressure (w/ BCs)
      gp0%grid(igrid)%array(:,:,:,1) = 2.*rho*tmp

      call restrictMGArray(IRHO,1,gp0,bcs(:,IRHO),igrid,order)

      !v0/dt+theta(v0.grad(v0)-mu*veclap(v0))

      allocate(vector(0:nxx+1,0:nyy+1,0:nzz+1,3))

      do k = 1,nzz
        do j = 1,nyy
          do i = 1,nxx
            ii  = i + nxx*(j-1) + nxx*nyy*(k-1)

            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                         ,x1,y1,z1,cartsn)
            jac    = jacobian(x1,y1,z1,cartsn)
            nabla_v= fnabla_v(i,j,k,nxx,nyy,nzz,x1,y1,z1
     .                       ,vx,vy,vz,cartsn,0)
            dvol   = volume(i,j,k,igx,igy,igz)

            !Eqn 1
            do icomp=1,3
              vector(i,j,k,icomp) = 
     .                            gv0%grid(igrid)%array(i,j,k,icomp)/dt
     .                          + alpha*vx(i,j,k)*nabla_v(1,icomp)/jac
     .                          + alpha*vy(i,j,k)*nabla_v(2,icomp)/jac
     .                          + alpha*vz(i,j,k)*nabla_v(3,icomp)/jac
     .                          - alpha*veclaplacian(i,j,k,nxx,nyy,nzz
     .                                              ,vx,vy,vz,nuu
     .                                              ,alt_eom,icomp)/dvol
            enddo

          enddo
        enddo
      enddo

      do icomp=1,3
        call restrictArrayToMGVector(1,nxx,nyy,nzz
     .                              ,vector(:,:,:,icomp)
     .                              ,mgadvdiffV0(:,icomp)
     .                              ,igrid,order,.false.)
      enddo

      deallocate(vector)

c Find required diagonals in all grids

      if (jit == 1) then
        icomp = IRHO
        call find_mf_diag(1,ntotdp,rho_mtvc,1,bcs(:,IRHO),rho_diag)

        icomp = ITMP
        call find_mf_diag(1,ntotdp,tmp_mtvc,1,bcs(:,ITMP),tmp_diag)

        icomp = IBX
        call find_mf_diag(3,3*ntotdp,b_mtvc,1,bcs(:,IBX:IBZ),b_diag)

        icomp = IVX
        call find_mf_diag(3,3*ntotdp,v_mtvc,1,bcs(:,IVX:IVZ),v_diag)
      endif

c diag
c DIAGONAL GENERATION FOR XDRAW PLOTTING
cc      do k=1,nzz
cc        do j=1,nyy
cc          do i=1,nxx
cc            ii = i + nxx*(j-1)+nxx*nyy*(k-1)
cc            debug(i,j,k) = 1./sqrt(sum(v_diag(3,3*ii-2:3*ii)**2)
cc          enddo
cc        enddo
cc      enddo
cc
cc      open(unit=110,file='debug.bin',form='unformatted'
cc     .    ,status='replace')
cc      call contour(debug(1:nxx,1:nyy,1),nxx,nyy,0d0,xmax,0d0,ymax,0,110)
cc      close(110)
cc      stop

c diag ***************
c MATRIX GENERATION FOR MATLAB PROCESSING
cc      igrid = 1
cc
cc      nxx = grid_params%nxv(igrid)
cc      nyy = grid_params%nyv(igrid)
cc      nzz = grid_params%nzv(igrid)
cc
cc      nt = 3*nxx*nyy*nzz
cc
cc      allocate(v_mat(nt,nt))
cc
cc      icomp = IVX
cc      call find_mf_mat(3,nt,v_mtvc,igrid,bcs(:,IVX:IVZ),v_mat)
cc
cccc      open(unit=110,file='debug.bin',form='unformatted'
cccc     .       ,status='replace')
cccc      call contour(v_mat,nt,nt,0d0,1d0,0d0,1d0,0,110)
cccc      close(110)
cc
cc      open(unit=110,file='debug.mat',status='replace')
cc      do i=1,nt
cc        write (110,*) (v_mat(i,j),j=1,nt)
cc      enddo
cc      close(110)
cc      stop

c End program

      end subroutine setupPreconditioner
