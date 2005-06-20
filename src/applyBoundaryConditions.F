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

      bcsq = bbcs(:,IBX)
      where (bcsq == DEF) bcsq = DIR
      bbcs(:,IBX) = bcsq

      bcsq = bbcs(:,IBY)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IBY) = bcsq

      bcsq = bbcs(:,IBZ)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IBZ) = bcsq

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

c Magnetic field BC

c     BC setup

      bcnd(:,1) = varray%array_var(IBX)%bconds
      bcnd(:,2) = varray%array_var(IBY)%bconds
      bcnd(:,3) = varray%array_var(IBZ)%bconds

      v_cnv(:,:,:,1) = varray%array_var(IBX)%array
      v_cnv(:,:,:,2) = varray%array_var(IBY)%array
      v_cnv(:,:,:,3) = varray%array_var(IBZ)%array

      v0(:,:,:,1) = u_0%array_var(IBX)%array
      v0(:,:,:,2) = u_0%array_var(IBY)%array
      v0(:,:,:,3) = u_0%array_var(IBZ)%array

c     Fill ghost nodes

      call setBC(IBX,3,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c     Postprocessing

      bx_cov = v_cov(:,:,:,1)
      by_cov = v_cov(:,:,:,2)
      bz_cov = v_cov(:,:,:,3)

      varray%array_var(IBX)%array = v_cnv(:,:,:,1)
      varray%array_var(IBY)%array = v_cnv(:,:,:,2)
      varray%array_var(IBZ)%array = v_cnv(:,:,:,3)

c Current BC

c     BC setup

      bcnd(:,1) = varray%array_var(IBX)%bconds
      bcnd(:,2) = varray%array_var(IBY)%bconds
      bcnd(:,3) = varray%array_var(IBZ)%bconds
      where (bcnd == -NEU)
cc        bcnd = DIR  !Use contravariant components for tangential dirichlet
        bcnd = -DIR  !Use covariant components for tangential dirichlet
      end where

      do k = 0,nnz+1
        do j = 0,nny+1
          do i = 0,nnx+1
            do icomp=1,3
              v_cnv(i,j,k,icomp)=curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                                ,bx_cov,by_cov,bz_cov,icomp)
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

      call setBC(IJX,3,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c     Postprocessing

      jx_cov = v_cov(:,:,:,1)
      jy_cov = v_cov(:,:,:,2)
      jz_cov = v_cov(:,:,:,3)

      jx = v_cnv(:,:,:,1)
      jy = v_cnv(:,:,:,2)
      jz = v_cnv(:,:,:,3)

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
