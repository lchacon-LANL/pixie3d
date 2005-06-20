c defineBoundaryConditions
c####################################################################
      subroutine defineBoundaryConditions (neq,bbcs)
c--------------------------------------------------------------------
c     Defines boundary conditions of physical quantities.
c     On input:
c       * neq -> number of equations
c     On output:
c       * bbcs -> real array of size (6,neq) containing BC setup:
c           + bbcs(1) ---> at x0
c           + bbcs(2) ---> at x1
c           + bbcs(3) ---> at y0
c           + bbcs(4) ---> at y1
c           + bbcs(5) ---> at z0
c           + bbcs(6) ---> at z1
c     Definition of BC identifiers is given in "grid_mod.f". In vectors,
c     a negative BC identifier means that BCs are to be imposed on
c     covariant components instead of on (default) contravariant comps.
c--------------------------------------------------------------------

      use icond

      use grid

      use equilibrium

      implicit none

c Call variables

      integer(4) :: neq,bbcs(6,neq)

c Local variables

      integer(4) :: ieq,bcsq(6)

c Begin program

c Default boundary conditions

      bcsq = bbcs(:,IRHO)
      where (bcsq == DEF) bcsq = NEU
      bbcs(:,IRHO) = bcsq

      bcsq = bbcs(:,IVX)
      where (bcsq == DEF) bcsq = DIR
      bbcs(:,IVX) = bcsq

      bcsq = bbcs(:,IVY)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IVY) = bcsq

      bcsq = bbcs(:,IVZ)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IVZ) = bcsq

      bcsq = bbcs(:,IAX)
      where (bcsq == DEF) bcsq = DIR   !On contravariant components
cc      where (bcsq == SP ) bcsq = -SP
      bbcs(:,IAX) = bcsq

      bcsq = bbcs(:,IAY)
      where (bcsq == DEF) bcsq = -EQU  !On covariant components
cc      where (bcsq == SP ) bcsq = -SP
      bbcs(:,IAY) = bcsq

      bcsq = bbcs(:,IAZ)
      where (bcsq == DEF) bcsq = -EQU  !On covariant components
cc      where (bcsq == SP ) bcsq = -SP
      bbcs(:,IAZ) = bcsq

cc      bcsq = bbcs(:,IBX)
cc      where (bcsq == DEF) bcsq = DIR
cc      bbcs(:,IBX) = bcsq
cc
cc      bcsq = bbcs(:,IBY)
cc      where (bcsq == DEF) bcsq = -NEU  !On covariant components
cc      bbcs(:,IBY) = bcsq
cc
cc      bcsq = bbcs(:,IBZ)
cc      where (bcsq == DEF) bcsq = -NEU  !On covariant components
cc      bbcs(:,IBZ) = bcsq

      bcsq = bbcs(:,ITMP)
      where (bcsq == DEF) bcsq = NEU !To allow isothermal case
      bbcs(:,ITMP) = bcsq

c Exceptions for specific equilibria

      select case (equil)

      case ('rfp1')

cc        bbcs(2,IBY) = EQU  !Imposed by equilibrium
cc
cc        bbcs(2,IBZ) = EQU  !Imposed by equilibrium

      case ('tmcyz','3dkh')

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = -NEU
        where (bbcs(:,IVY) ==-NEU) bbcs(:,IVY) =  DIR

        where (bbcs(:,IBX) == DIR) bbcs(:,IBX) = -NEU
        where (bbcs(:,IBY) ==-NEU) bbcs(:,IBY) =  DIR

      end select

c End

      end subroutine defineBoundaryConditions

c imposeBoundaryConditions
c####################################################################
      subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use imposeBCinterface

      use auxiliaryVariables

      implicit none

c Call variables

      integer(4) :: iigx,iigy,iigz

      type (var_array) :: varray

c Local variables

      integer(4) :: i,j,k,icomp,bcnd(6,3),nnx,nny,nnz,order

c Begin program

      order = 2   !Order of extrapolation in BCs

      call setASMflag(.false.)  !Not additive Schwartz method

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

c Allocate auxiliary variables in local domain

      allocate(v_cnv(0:nnx+1,0:nny+1,0:nnz+1,3)
     .        ,v_cov(0:nnx+1,0:nny+1,0:nnz+1,3)
     .        ,v0   (0:nnx+1,0:nny+1,0:nnz+1,3))

c Density BC

      call setBC(IRHO,nnx,nny,nnz
     .               ,varray%array_var(IRHO)%array
     .               ,u_0   %array_var(IRHO)%array
     .               ,varray%array_var(IRHO)%bconds
     .               ,iigx,iigy,iigz,iorder=order)

c Velocity BC

c     BC setup

      bcnd(:,1) = varray%array_var(IVX)%bconds
      bcnd(:,2) = varray%array_var(IVY)%bconds
      bcnd(:,3) = varray%array_var(IVZ)%bconds

      if (.not.nc_eom_v) then
        where (varray%array_var(IRHO)%array /= 0d0)
          v_cnv(:,:,:,1) = varray%array_var(IVX )%array
     .                    /varray%array_var(IRHO)%array
          v_cnv(:,:,:,2) = varray%array_var(IVY )%array
     .                    /varray%array_var(IRHO)%array
          v_cnv(:,:,:,3) = varray%array_var(IVZ )%array
     .                    /varray%array_var(IRHO)%array
        end where
      else
        v_cnv(:,:,:,1) = varray%array_var(IVX)%array
        v_cnv(:,:,:,2) = varray%array_var(IVY)%array
        v_cnv(:,:,:,3) = varray%array_var(IVZ)%array
      endif

      v0 = v_cnv

c     Fill ghost nodes

      call setBC(IVX,3,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c     Postprocessing

      vx_cov = v_cov(:,:,:,1)
      vy_cov = v_cov(:,:,:,2)
      vz_cov = v_cov(:,:,:,3)

      vx     = v_cnv(:,:,:,1)
      vy     = v_cnv(:,:,:,2)
      vz     = v_cnv(:,:,:,3)

      if (.not.nc_eom_v) then
        varray%array_var(IVX)%array = v_cnv(:,:,:,1)
     .                               *varray%array_var(IRHO)%array
        varray%array_var(IVY)%array = v_cnv(:,:,:,2)
     .                               *varray%array_var(IRHO)%array
        varray%array_var(IVZ)%array = v_cnv(:,:,:,3)
     .                               *varray%array_var(IRHO)%array
      else
        varray%array_var(IVX)%array = v_cnv(:,:,:,1)
        varray%array_var(IVY)%array = v_cnv(:,:,:,2)
        varray%array_var(IVZ)%array = v_cnv(:,:,:,3)
      endif

c Vector Pot. BC

c     BC setup

      bcnd(:,1) = varray%array_var(IAX)%bconds
      bcnd(:,2) = varray%array_var(IAY)%bconds
      bcnd(:,3) = varray%array_var(IAZ)%bconds

      v_cov(:,:,:,1) = varray%array_var(IAX)%array
      v_cov(:,:,:,2) = varray%array_var(IAY)%array
      v_cov(:,:,:,3) = varray%array_var(IAZ)%array

      v0(:,:,:,1) = u_0%array_var(IAX)%array
      v0(:,:,:,2) = u_0%array_var(IAY)%array
      v0(:,:,:,3) = u_0%array_var(IAZ)%array

c     Fill ghost nodes

      call setBC(IAX,3,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd,iigx,iigy,iigz
     .          ,is_cnv=.false.,iorder=order)

c     Postprocessing

      varray%array_var(IAX)%array = v_cov(:,:,:,1)
      varray%array_var(IAY)%array = v_cov(:,:,:,2)
      varray%array_var(IAZ)%array = v_cov(:,:,:,3)

      acnv = v_cnv

c Magnetic field BC

c     BC setup  (only topological constraints)

cc      bcnd(:,1) = bcond
cc      bcnd(:,2) = bcond
cc      bcnd(:,3) = bcond
      bcnd(:,1) = varray%array_var(IAX)%bconds
      bcnd(:,2) = varray%array_var(IAY)%bconds
      bcnd(:,3) = varray%array_var(IAZ)%bconds
      where (bcnd == -EQU)
cc        bcnd = -NEU !Use covariant components for tangential Neumann
        bcnd = EXT !Use covariant components for tangential Neumann
      end where

c     B=curl(A)

      do k = 0,nnz+1
        do j = 0,nny+1
          do i = 0,nnx+1
            do icomp=1,3
              !Here, v_cov=a per previous step
              v_cnv(i,j,k,icomp)=curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                                ,v_cov(:,:,:,1)
     .                                ,v_cov(:,:,:,2)
     .                                ,v_cov(:,:,:,3),icomp)
            enddo
          enddo
        enddo
      enddo

#if defined(petsc) 
      do icomp=1,3
        call fillPetscGhostCells(v_cnv(:,:,:,icomp))
      enddo
#endif

      v0 = v_cnv

c     Fill ghost nodes

      call setBC(IBX,3,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c     Postprocessing

cc      do k = 0,nnz+1
cc        do j = 0,nny+1
cc          do i = 0,nnx+1
cc            call transformFromCurvToCurv(i,j,k,iigx,iigy,iigz
cc     .           ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
cc     .           ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),.false.)
cc
cc          enddo
cc        enddo
cc      enddo

      bx     = v_cnv(:,:,:,1)
      by     = v_cnv(:,:,:,2)
      bz     = v_cnv(:,:,:,3)

      bx_cov = v_cov(:,:,:,1)
      by_cov = v_cov(:,:,:,2)
      bz_cov = v_cov(:,:,:,3)

c Current

c     BC setup (same as B)

cc      bcnd(:,1) = bcond
cc      bcnd(:,2) = bcond
cc      bcnd(:,3) = bcond

cc      bcnd(:,1) = varray%array_var(IAX)%bconds
cc      bcnd(:,2) = varray%array_var(IAY)%bconds
cc      bcnd(:,3) = varray%array_var(IAZ)%bconds
cc      where (bcnd == -EQU)
cccc        bcnd = -NEU !Use covariant components for tangential Neumann
cc        bcnd = EXT !Use covariant components for tangential Neumann
cc      end where

c     J=curl(B)

      do k = 0,nnz+1
        do j = 0,nny+1
          do i = 0,nnx+1
            do icomp=1,3
              !Here, v_cov=a per previous step
              v_cnv(i,j,k,icomp)=curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                                ,v_cov(:,:,:,1)
     .                                ,v_cov(:,:,:,2)
     .                                ,v_cov(:,:,:,3),icomp)
            enddo
          enddo
        enddo
      enddo

#if defined(petsc) 
      do icomp=1,3
        call fillPetscGhostCells(v_cnv(:,:,:,icomp))
      enddo
#endif

      v0 = v_cnv

c     Fill ghost nodes

      call setBC(IBX,3,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c     Postprocessing

cc      do k = 0,nnz+1
cc        do j = 0,nny+1
cc          do i = 0,nnx+1
cc            call transformFromCurvToCurv(i,j,k,iigx,iigy,iigz
cc     .           ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
cc     .           ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),.false.)
cc
cc          enddo
cc        enddo
cc      enddo

      jx     = v_cnv(:,:,:,1)
      jy     = v_cnv(:,:,:,2)
      jz     = v_cnv(:,:,:,3)

      jx_cov = v_cov(:,:,:,1)
      jy_cov = v_cov(:,:,:,2)
      jz_cov = v_cov(:,:,:,3)

ccc Current BC
cc
ccc     BC setup
cc
cc      bcnd(:,1) = varray%array_var(IBX)%bconds
cc      bcnd(:,2) = varray%array_var(IBY)%bconds
cc      bcnd(:,3) = varray%array_var(IBZ)%bconds
cc      where (bcnd == -NEU)
cccc        bcnd = DIR  !Use contravariant components for tangential dirichlet
cc        bcnd = -DIR  !Use covariant components for tangential dirichlet
cc      end where
cc
cc      do k = 0,nnz+1
cc        do j = 0,nny+1
cc          do i = 0,nnx+1
cc            do icomp=1,3
cc              v_cnv(i,j,k,icomp)=curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
cc     .                                ,bx_cov,by_cov,bz_cov,icomp)
cc            enddo
cc          enddo
cc        enddo
cc      enddo
cc
cc#if defined(petsc) 
cc      do icomp=1,3
cc        call fillPetscGhostCells(v_cnv(:,:,:,icomp))
cc      enddo
cc#endif
cc
cc      v0 = v_cnv
cc
ccc     Fill ghost nodes
cc
cc      call setBC(IJX,3,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd,iigx,iigy,iigz
cc     .          ,iorder=order)
cc
ccc     Postprocessing
cc
cc      jx_cov = v_cov(:,:,:,1)
cc      jy_cov = v_cov(:,:,:,2)
cc      jz_cov = v_cov(:,:,:,3)
cc
cc      jx = v_cnv(:,:,:,1)
cc      jy = v_cnv(:,:,:,2)
cc      jz = v_cnv(:,:,:,3)

c Temperature BCs

      call setBC(ITMP,nnx,nny,nnz
     .               ,varray%array_var(ITMP)%array
     .               ,u_0   %array_var(ITMP)%array
     .               ,varray%array_var(ITMP)%bconds
     .               ,iigx,iigy,iigz,iorder=order)

c Deallocate variables 

      deallocate(v_cnv,v_cov,v0)

c End

      end subroutine imposeBoundaryConditions

c imposeBConFluxes
c####################################################################
      subroutine imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                            ,flxkp,flxkm,bconds)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on fluxes
c--------------------------------------------------------------------

      use singularBCinterface

      implicit none

c Call variables

      integer(4) :: i,j,k,bconds(6)
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

c Local variables

c Begin program

c End

      end subroutine imposeBConfluxes
