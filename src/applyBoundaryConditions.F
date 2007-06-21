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

      integer    :: neq,bbcs(6,neq)

c Local variables

      integer    :: ieq,bcsq(6)

c Begin program

c Default boundary conditions

      bcsq = bbcs(:,IRHO)
      where (bcsq == DEF) bcsq = NEU
      bbcs(:,IRHO) = bcsq

      bcsq = bbcs(:,IVX)
      where (bcsq == DEF) bcsq = DIR
      bbcs(:,IVX) = bcsq

      bcsq = bbcs(:,IVY)
      if (     coords == 'cyl'
     .    .or. coords == 'hel'
     .    .or. coords == 'tor'
     .    .or. coords == 'car') then
        where (bcsq == DEF) bcsq =  NEU !On contravariant components
      else
        where (bcsq == DEF) bcsq = -NEU !On covariant components
      endif
      bbcs(:,IVY) = bcsq

      bcsq = bbcs(:,IVZ)
      if (     coords == 'cyl'
     .    .or. coords == 'hel'
     .    .or. coords == 'tor'
     .    .or. coords == 'car') then
        where (bcsq == DEF) bcsq =  NEU !On contravariant components
      else
        where (bcsq == DEF) bcsq = -NEU !On covariant components
      endif
      bbcs(:,IVZ) = bcsq

#if !defined(vec_pot)
      bcsq = bbcs(:,IBX)
      where (bcsq == DEF) bcsq = DIR
      bbcs(:,IBX) = bcsq

      bcsq = bbcs(:,IBY)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IBY) = bcsq

      bcsq = bbcs(:,IBZ)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IBZ) = bcsq
#else
      bcsq = bbcs(:,IAX)
      where (bcsq == DEF) bcsq = DIR   !On contravariant components
cc      where (bcsq == DEF) bcsq = EXT   !On contravariant components
cc      where (bcsq == DEF) bcsq = -EQU  !On contravariant components
      bbcs(:,IAX) = bcsq

      bcsq = bbcs(:,IAY)
      where (bcsq == DEF) bcsq = -EQU  !On covariant components
      bbcs(:,IAY) = bcsq

      bcsq = bbcs(:,IAZ)
      where (bcsq == DEF) bcsq = -EQU  !On covariant components
      bbcs(:,IAZ) = bcsq
#endif

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
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) =  DIR

        where (bbcs(:,IBX) == DIR) bbcs(:,IBX) = -NEU
        where (bbcs(:,IBY) == NEU) bbcs(:,IBY) =  DIR

      case ('ppnch','ppn3d','ppnst')   !No stress BCs

cc        where (bbcs(:,IRHO) == NEU) bbcs(:,IRHO) = EQU

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX ) = EQU

      case ('ppnsl')                   !No slip BCs

cc        where (bbcs(:,IRHO) == NEU) bbcs(:,IRHO) = EQU

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX ) = EQU

        !Dirichlet boundary conditions for tangential components
        where (bbcs(:,IVY ) == NEU) bbcs(:,IVY ) = DIR
        where (bbcs(:,IVZ ) == NEU) bbcs(:,IVZ ) = DIR

c$$$      case ('ic')
c$$$
c$$$        bbcs(2,:) = EQU

      end select

c End

      end subroutine defineBoundaryConditions

c imposeBoundaryConditions
c####################################################################
      subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variables

      implicit none

c Call variables

      integer    :: iigx,iigy,iigz

      type (var_array) :: varray

c Local variables

      integer    :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz,order

      real(8)    :: bbcov(3),bbcnv(3),b2,J0(3),E00(3)

      real(8),allocatable,dimension(:,:,:,:) :: v0

c Begin program

      order = 2   !Order of extrapolation in BCs

cc      call setup_petsc_BC !Setup parallel BC flags

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

#if !defined(vec_pot)

c Magnetic field BC

c     Resistivity (needed for BCs; may depend on T)

      do k=0,nnz+1
        do j=0,nny+1
          do i=0,nnx+1
            eeta (i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
          enddo
        enddo
      enddo

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

      bcov = v_cov
      bcnv = v_cnv

      varray%array_var(IBX)%array = v_cnv(:,:,:,1)
      varray%array_var(IBY)%array = v_cnv(:,:,:,2)
      varray%array_var(IBZ)%array = v_cnv(:,:,:,3)

#else

c Vector Pot. BC

c     Resistivity (needed for BCs; may depend on T)

      do k=0,nnz+1
        do j=0,nny+1
          do i=0,nnx+1
            eeta (i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
          enddo
        enddo
      enddo

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
      acov = v_cov

c Magnetic field BC

c     BC setup  (only topological constraints)

      bcsq = bcond
      where (bcsq == DEF) bcsq = DIR
      bcnd(:,1) = bcsq

cc      if (equil == 'spnch' .or. equil == 'rfp2') then
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cc        bcnd(:,2) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cc        bcnd(:,3) = bcsq
cc      else
        bcsq = bcond
        where (bcsq == DEF) bcsq = -NEU !On covariant components
        bcnd(:,2) = bcsq

        bcsq = bcond
        where (bcsq == DEF) bcsq = -NEU !On covariant components
        bcnd(:,3) = bcsq
cc      endif

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

      v0 = v_cnv

c     Fill ghost nodes

      call setBC(IBX,3,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c     Postprocessing

      bcnv = v_cnv
      bcov = v_cov

#endif

c Current BC

c     BC setup

#if !defined(vec_pot)
      bcnd(:,1) = varray%array_var(IBX)%bconds
      bcnd(:,2) = varray%array_var(IBY)%bconds
      bcnd(:,3) = varray%array_var(IBZ)%bconds
      where (bcnd == -NEU)
cc        bcnd = -DIR  !Use covariant components for tangential BCs
        bcnd = -EQU  !Use covariant components for tangential BCs
      end where
#else
cc      if (equil == 'spnch' .or. equil == 'rfp2') then
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = DIR
cc        bcnd(:,1) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = EXT
cccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cc        bcnd(:,2) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = EXT
cccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cc        bcnd(:,3) = bcsq
cc      else
        bcnd(:,1) = varray%array_var(IAX)%bconds
        bcnd(:,2) = varray%array_var(IAY)%bconds
        bcnd(:,3) = varray%array_var(IAZ)%bconds
cc      endif
#endif

      do k = 0,nnz+1
        do j = 0,nny+1
          do i = 0,nnx+1
            do icomp=1,3
              v_cnv(i,j,k,icomp)=curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                                ,bcov(:,:,:,1)
     .                                ,bcov(:,:,:,2)
     .                                ,bcov(:,:,:,3),icomp)
            enddo
          enddo
        enddo
      enddo

      if (    equil == 'ppnch' .or. equil == 'ppn3d'
     .    .or.equil == 'ppnsl' .or. equil == 'ppnst' ) then       !Jt = alpha*Bt in cov representation
cc      if (sum(abs(E0)) /= 0d0) then  !Jt = alpha*Bt in cov representation

        i = nnx

        do k=1,nnz
          do j=1,nny
#if !defined(vec_pot)
            bbcnv = 0.5*(/ varray%array_var(IBX)%array(i  ,j,k)
     .                    ,varray%array_var(IBY)%array(i  ,j,k)
     .                    ,varray%array_var(IBZ)%array(i  ,j,k) /)
     .             +0.5*(/ varray%array_var(IBX)%array(i+1,j,k)
     .                    ,varray%array_var(IBY)%array(i+1,j,k)
     .                    ,varray%array_var(IBZ)%array(i+1,j,k) /)
#else
            bbcnv = 0.5*(/ bcnv(i  ,j,k,1)
     .                    ,bcnv(i  ,j,k,2)
     .                    ,bcnv(i  ,j,k,3) /)
     .             +0.5*(/ bcnv(i+1,j,k,1)
     .                    ,bcnv(i+1,j,k,2)
     .                    ,bcnv(i+1,j,k,3) /)
#endif
            bbcov = matmul(0.5*(gmetric%grid(iigx)%gsub(i  ,j,k,:,:)
     .                         +gmetric%grid(iigx)%gsub(i+1,j,k,:,:))
     .                    ,bbcnv)

            b2   = 2*dot_product(bbcnv,bbcov)
     .            /(gmetric%grid(iigx)%jac(i  ,j,k)
     .             +gmetric%grid(iigx)%jac(i+1,j,k))

            J0 = E0*2./(eeta(i,j,k)+eeta(i+1,j,k))

            v0(i+1,j,k,:) = scalarProduct(i,j,k,iigx,iigy,iigz
     .                                 ,J0   (1),J0   (2),J0   (3)
     .                                 ,bbcnv(1),bbcnv(2),bbcnv(3))/b2
     .                     *bbcov(:)
          enddo
        enddo

      else

        v0 = 0d0

      endif

c     Fill ghost nodes

      call setBC(IJX,3,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c     Postprocessing

      jcov = v_cov
      jcnv = v_cnv

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

      if (    equil == 'ppnch' .or. equil == 'ppn3d'
     .    .or.equil == 'ppnsl' .or. equil == 'ppnst' ) then

        i = nnx

        do k=1,nnz
          do j=1,nny
#if !defined(vec_pot)
            bbcnv = 0.5*(/ varray%array_var(IBX)%array(i  ,j,k)
     .                    ,varray%array_var(IBY)%array(i  ,j,k)
     .                    ,varray%array_var(IBZ)%array(i  ,j,k) /)
     .             +0.5*(/ varray%array_var(IBX)%array(i+1,j,k)
     .                    ,varray%array_var(IBY)%array(i+1,j,k)
     .                    ,varray%array_var(IBZ)%array(i+1,j,k) /)
#else
            bbcnv = 0.5*(/ bcnv(i  ,j,k,1)
     .                    ,bcnv(i  ,j,k,2)
     .                    ,bcnv(i  ,j,k,3) /)
     .             +0.5*(/ bcnv(i+1,j,k,1)
     .                    ,bcnv(i+1,j,k,2)
     .                    ,bcnv(i+1,j,k,3) /)
cc            bbcnv = 0.5*(/ bx(i  ,j,k)
cc     .                    ,by(i  ,j,k)
cc     .                    ,bz(i  ,j,k) /)
cc     .             +0.5*(/ bx(i+1,j,k)
cc     .                    ,by(i+1,j,k)
cc     .                    ,bz(i+1,j,k) /)
#endif
            bbcov = matmul(0.5*(gmetric%grid(iigx)%gsub(i  ,j,k,:,:)
     .                         +gmetric%grid(iigx)%gsub(i+1,j,k,:,:))
     .                    ,bbcnv)

            b2   = 2*dot_product(bbcnv,bbcov)
     .            /(gmetric%grid(iigx)%jac(i  ,j,k)
     .             +gmetric%grid(iigx)%jac(i+1,j,k))

            if (b2 > 0d0) then
cc              E00(1) = E0(1) - 0.5*(eeta(i  ,j,k)*jcov(i  ,j,k,1)
cc     .                             +eeta(i+1,j,k)*jcov(i+1,j,k,1))
              E00(2) = E0(2) - 0.5*(eeta(i  ,j,k)*jcov(i  ,j,k,2)
     .                             +eeta(i+1,j,k)*jcov(i+1,j,k,2))
              E00(3) = E0(3) - 0.5*(eeta(i  ,j,k)*jcov(i  ,j,k,3)
     .                             +eeta(i+1,j,k)*jcov(i+1,j,k,3))

              v0(i+1,j,k,1) = -(bbcov(2)*E00(3)-bbcov(3)*E00(2))/b2 !ExB velocity
            else
              v0(i+1,j,k,1) = 0d0
            endif
          enddo
        enddo

      else

        v0(:,:,:,1) = u_0%array_var(IVX)%array

      endif

      v0(:,:,:,2) = u_0%array_var(IVY)%array
      v0(:,:,:,3) = u_0%array_var(IVZ)%array

c     Fill ghost nodes

      call setBC(IVX,3,nnx,nny,nnz,v_cnv,v_cov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c     Postprocessing

      vcov = v_cov
      vcnv = v_cnv

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

      integer    :: i,j,k,bconds(6)
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

c Local variables

c Begin program

c End

      end subroutine imposeBConfluxes
