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

      use variables

      use auxiliaryVariables

      use local_BCS_variables

      implicit none

c Call variables

      integer    :: neq,bbcs(6,neq)

c Local variables

      integer    :: ieq,bcsq(6)

c Begin program

c Reset BCs

      do ieq=1,neq
        bbcs(:,ieq) = bcond
      enddo

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
      case ('tmcyz','3dkh')

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = -NEU
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) =  DIR

        where (bbcs(:,IBX) == DIR) bbcs(:,IBX) = -NEU
        where (bbcs(:,IBY) == NEU) bbcs(:,IBY) =  DIR

      case ('ppnch','ppn3d','ppnst','2dhel')   !No stress BCs, pinch flow

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX ) = EQU

      case ('ppnsl','p3nsl','ohtor')           !No slip BCs, pinch flow

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX ) = EQU

        !Dirichlet boundary conditions for tangential components
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY ) = DIR
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ ) = DIR

cc      case ('vmec','tor','gem')
      case ('vmec','tor')

        !Dirichlet boundary conditions for tangential components
        where (abs(bbcs(:,IVY)) == NEU) bbcs(:,IVY) = DIR
        where (abs(bbcs(:,IVZ)) == NEU) bbcs(:,IVZ) = DIR
cc        if (solve_rho) 
cc     .       where (abs(bbcs(:,IRHO)) == NEU) bbcs(:,IRHO) = DIR
        if (gamma > 1d0)
     .       where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = DIR

cc      case ('ntm2d')                   !No slip BCs
cc
cc        !Dirichlet boundary conditions for tangential components
cc        where (abs(bbcs(:,IVY)) == NEU) bbcs(:,IVY) = DIR
cc        where (abs(bbcs(:,IVZ)) == NEU) bbcs(:,IVZ) = DIR

      end select

      pinch_flow = (    equil == 'ppnch' .or. equil == 'ppn3d'
     .              .or.equil == 'ppnsl' .or. equil == 'ppnst'
     .              .or.equil == '2dhel' .or. equil == 'p3nsl'
     .              .or.equil == 'ohtor')

c BCs for auxiliary variables

      call defineAuxBCs(neq,bbcs)

c End subroutine

      end subroutine defineBoundaryConditions

c     defineAuxBCs
c     ###################################################################
      subroutine defineAuxBCs(neq,bcs)

      use variables

      use equil_def

c     -------------------------------------------------------------------
c     Defines BCs for auxiliary quantities
c     -------------------------------------------------------------------

c     Call variables

      integer    :: bcs(6,neq)

c     Local variables

      integer    :: ieq,bcsq(6)

c     Begin program

#if defined(vec_pot)
      !Vector potential
      gv%aux%vec_list(IACNV)%bconds(:,1) = bcs(:,IAX)
      gv%aux%vec_list(IACNV)%bconds(:,2) = bcs(:,IAY)
      gv%aux%vec_list(IACNV)%bconds(:,3) = bcs(:,IAZ)

      gv%aux%vec_list(IACOV)%bconds = gv%aux%vec_list(IACNV)%bconds

      !Magnetic field
      bcsq = bcond
      where (bcsq == DEF) bcsq = DIR
      gv%aux%vec_list(IBCNV)%bconds(:,1) = bcsq

      bcsq = bcond
      where (bcsq == DEF) bcsq = -NEU !On covariant components
      gv%aux%vec_list(IBCNV)%bconds(:,2) = bcsq

      bcsq = bcond
      where (bcsq == DEF) bcsq = -NEU !On covariant components
      gv%aux%vec_list(IBCNV)%bconds(:,3) = bcsq

      gv%aux%vec_list(IBCOV)%bconds = gv%aux%vec_list(IBCNV)%bconds

      !Current
      gv%aux%vec_list(IJCNV)%bconds = gv%aux%vec_list(IACNV)%bconds
      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds
#else
      !Magnetic field
      gv%aux%vec_list(IBCNV)%bconds(:,1) = bcs(:,IBX)
      gv%aux%vec_list(IBCNV)%bconds(:,2) = bcs(:,IBY)
      gv%aux%vec_list(IBCNV)%bconds(:,3) = bcs(:,IBZ)

      gv%aux%vec_list(IBCOV)%bconds = gv%aux%vec_list(IBCNV)%bconds


      !Current
      gv%aux%vec_list(IJCNV)%bconds(:,1) = bcs(:,IBX)
      gv%aux%vec_list(IJCNV)%bconds(:,2) = bcs(:,IBY)
      gv%aux%vec_list(IJCNV)%bconds(:,3) = bcs(:,IBZ)

      where (gv%aux%vec_list(IJCNV)%bconds == -NEU)
        gv%aux%vec_list(IJCNV)%bconds = -EQU  !Use covariant components for tangential BCs
      end where

      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds
#endif

      !Equilibrium Current (extrapolate covariant components)
#if defined(vec_pot)
      gv%aux%vec_list(IJ0CNV)%bconds = gv%aux%vec_list(IACNV)%bconds
      where (gv%aux%vec_list(IJ0CNV)%bconds == -EQU)
        gv%aux%vec_list(IJ0CNV)%bconds = -EXT
      end where
#else
      gv%aux%vec_list(IJ0CNV)%bconds = gv%aux%vec_list(IBCNV)%bconds
      where (gv%aux%vec_list(IJ0CNV)%bconds == -NEU)
        gv%aux%vec_list(IJ0CNV)%bconds = -EXT
      end where
#endif
      gv%aux%vec_list(IJ0COV)%bconds = gv%aux%vec_list(IJ0CNV)%bconds

      !Flow
      gv%aux%vec_list(IVCNV)%bconds(:,1) = bcs(:,IVX)
      gv%aux%vec_list(IVCNV)%bconds(:,2) = bcs(:,IVY)
      gv%aux%vec_list(IVCNV)%bconds(:,3) = bcs(:,IVZ)

      gv%aux%vec_list(IVCOV)%bconds = gv%aux%vec_list(IVCNV)%bconds

c     End program

      end subroutine defineAuxBCs

cc# if 0
cc
ccc imposeBoundaryConditions
ccc####################################################################
cc      subroutine imposeBoundaryConditions(varray,iigx,iigy,iigz,obc)
ccc--------------------------------------------------------------------
ccc     Sets adequate boundary conditions on array structure varray.
ccc--------------------------------------------------------------------
cc
cc      use local_BCS_variables
cc
cc      use imposeBCinterface
cc
cc      use auxiliaryVariables
cc
cc      use variable_setup
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer    :: iigx,iigy,iigz
cc      integer,optional :: obc
cc
cc      type(var_array),pointer :: varray
cc
ccc Local variables
cc
cc      integer    :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz,order
cc
cc      real(8)    :: bbcov(3),bbcnv(3),b2,J0(3),E00(3)
cc
cc      real(8),allocatable,dimension(:,:,:,:) :: v0
cc
ccc Begin program
cc
cc      if (PRESENT(obc)) then
cc        order = obc
cc      else
cc        order = bc_order     !Order of extrapolation in BCs (defined in imposeBC_mod.F)
cc      endif
cc
ccc Local grid sizes
cc
cc      nnx = grid_params%nxv(iigx) 
cc      nny = grid_params%nyv(iigy)
cc      nnz = grid_params%nzv(iigz)
cc
ccc Allocate auxiliary variables in local domain
cc
cc#if defined(samrai)
cc      call allocAuxVariables
cc#endif
cc
cc      allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,3))
cc
ccc Density BC
cc
cc      call setBC(IRHO,nnx,nny,nnz
cc     .               ,varray%array_var(IRHO)%array
cc     .               ,u_0   %array_var(IRHO)%array
cc     .               ,varray%array_var(IRHO)%bconds
cc     .               ,iigx,iigy,iigz,iorder=order)
cc
cc#if !defined(vec_pot)
cc
ccc Magnetic field BC
cc
ccc     Resistivity (needed for BCs; may depend on T)
cc
cc      do k=0,nnz+1
cc        do j=0,nny+1
cc          do i=0,nnx+1
cc            eeta (i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
cc          enddo
cc        enddo
cc      enddo
cc
ccc     BC setup
cc
cc      bcnv(:,:,:,1) = varray%array_var(IBX)%array
cc      bcnv(:,:,:,2) = varray%array_var(IBY)%array
cc      bcnv(:,:,:,3) = varray%array_var(IBZ)%array
cc
cc      if (pinch_flow) then
cc        v0 = 0d0
cc      else
cc        v0 = jcov_0             !To fix current at boundary
cc      endif
cc
ccc     BC (Neumann)
cc
cc      bcnd(:,1) = varray%array_var(IBX)%bconds
cc      bcnd(:,2) = varray%array_var(IBY)%bconds
cc      bcnd(:,3) = varray%array_var(IBZ)%bconds
cc
cc      call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd,iigx,iigy,iigz
cc     .          ,iorder=order)
cc
ccc     Postprocessing
cc
cc      varray%array_var(IBX)%array = bcnv(:,:,:,1)
cc      varray%array_var(IBY)%array = bcnv(:,:,:,2)
cc      varray%array_var(IBZ)%array = bcnv(:,:,:,3)
cc
cc#else
cc
ccc Vector Pot. BC
cc
ccc     Resistivity (needed for BCs; may depend on T)
cc
cc      do k=0,nnz+1
cc        do j=0,nny+1
cc          do i=0,nnx+1
cc            eeta (i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
cc          enddo
cc        enddo
cc      enddo
cc
ccc     BC setup
cc
cc      bcnd(:,1) = varray%array_var(IAX)%bconds
cc      bcnd(:,2) = varray%array_var(IAY)%bconds
cc      bcnd(:,3) = varray%array_var(IAZ)%bconds
cc
cc      acov(:,:,:,1) = varray%array_var(IAX)%array
cc      acov(:,:,:,2) = varray%array_var(IAY)%array
cc      acov(:,:,:,3) = varray%array_var(IAZ)%array
cc
cc      v0(:,:,:,1) = u_0%array_var(IAX)%array
cc      v0(:,:,:,2) = u_0%array_var(IAY)%array
cc      v0(:,:,:,3) = u_0%array_var(IAZ)%array
cc
ccc     Fill ghost nodes
cc
cc      call setBC(IAX,3,nnx,nny,nnz,acnv,acov,v0,bcnd,iigx,iigy,iigz
cc     .          ,is_cnv=.false.,iorder=order)
cc
ccc     Postprocessing
cc
cc      varray%array_var(IAX)%array = acov(:,:,:,1)
cc      varray%array_var(IAY)%array = acov(:,:,:,2)
cc      varray%array_var(IAZ)%array = acov(:,:,:,3)
cc
ccc Magnetic field BC
cc
ccc     BC setup  (only topological constraints)
cc
cc      bcsq = bcond
cc      where (bcsq == DEF) bcsq = DIR
cc      bcnd(:,1) = bcsq
cc
cc      bcsq = bcond
cc      where (bcsq == DEF) bcsq = -NEU !On covariant components
cc      bcnd(:,2) = bcsq
cc
cc      bcsq = bcond
cc      where (bcsq == DEF) bcsq = -NEU !On covariant components
cc      bcnd(:,3) = bcsq
cc
ccc     B=curl(A)
cc
cc      do k = 0,nnz+1
cc        do j = 0,nny+1
cc          do i = 0,nnx+1
cc            do icomp=1,3
cc              !Here, v_cov=a per previous step
cc              bcnv(i,j,k,icomp)=B0(icomp)
cc     .                         +curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
cc     .                                ,acov(:,:,:,1)
cc     .                                ,acov(:,:,:,2)
cc     .                                ,acov(:,:,:,3),icomp)
cc            enddo
cc          enddo
cc        enddo
cc      enddo
cc
cc      v0 = jcov_0   !To fix current at boundary
cc
ccc     Fill ghost nodes
cc
cc      call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd,iigx,iigy,iigz
cc     .          ,iorder=order)
cc
cc#endif
cc
ccc Current BC
cc
ccc     BC setup
cc
cc#if !defined(vec_pot)
cc      bcnd(:,1) = varray%array_var(IBX)%bconds
cc      bcnd(:,2) = varray%array_var(IBY)%bconds
cc      bcnd(:,3) = varray%array_var(IBZ)%bconds
cc      where (bcnd == -NEU)
cc        bcnd = -EQU             !Use covariant components for tangential BCs
cc      end where
cc#else
cc      bcnd(:,1) = varray%array_var(IAX)%bconds
cc      bcnd(:,2) = varray%array_var(IAY)%bconds
cc      bcnd(:,3) = varray%array_var(IAZ)%bconds
cc#endif
cc
cc#if defined(vec_pot)
cc      do k = 0,nnz+1
cc        do j = 0,nny+1
cc          do i = 0,nnx+1
cc            do icomp=1,3
cc              jcnv(i,j,k,icomp)=curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
cc     .                                ,bcov(:,:,:,1)
cc     .                                ,bcov(:,:,:,2)
cc     .                                ,bcov(:,:,:,3),icomp)
cc            enddo
cc          enddo
cc        enddo
cc      enddo
cccc      do k = 1,nnz
cccc        do j = 1,nny
cccc          do i = 1,nnx
cccc              jcnv(i,j,k,:)= matmul(gmetric%grid(iigx)%gsup(i,j,k,:,:),
cccc     .             -veclap_cov(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,acov
cccc     .                        ,vol=.false.)
cccc     .                +graddiv(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,acnv
cccc     .                        ,vol=.false.))
cccc          enddo
cccc        enddo
cccc      enddo
cc#else
cc      do k = 0,nnz+1
cc        do j = 0,nny+1
cc          do i = 0,nnx+1
cc            do icomp=1,3
cc              jcnv(i,j,k,icomp)=curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
cc     .                                ,bcov(:,:,:,1)
cc     .                                ,bcov(:,:,:,2)
cc     .                                ,bcov(:,:,:,3),icomp)
cc            enddo
cc          enddo
cc        enddo
cc      enddo
cc#endif
cc
cc      if (pinch_flow) then  !Jt = alpha*Bt in cov representation
cc
cc        i = nnx
cc
cc        do k=1,nnz
cc          do j=1,nny
cc            bbcnv = 0.5*(bcnv(i+1,j,k,:)+bcnv(i,j,k,:))
cc            bbcov = 0.5*(bcov(i+1,j,k,:)+bcov(i,j,k,:))
cc
cc            J0 = 0.5*(jcov_0(i+1,j,k,:)+jcov_0(i,j,k,:))
cc
cc            v0(i+1,j,k,:) = dot_product(J0   ,bbcnv)
cc     .                     /dot_product(bbcov,bbcnv)*bbcov
cc          enddo
cc        enddo
cc
cc      else
cc
cc        v0(nnx+1,:,:,:) = 0.5*(jcov_0(nnx  ,:,:,:)
cc     .                        +jcov_0(nnx+1,:,:,:))
cc        v0(0    ,:,:,:) = 0.5*(jcov_0(1    ,:,:,:)
cc     .                        +jcov_0(0    ,:,:,:))
cc
cc        v0(:,nny+1,:,:) = 0.5*(jcov_0(:,nny  ,:,:)
cc     .                        +jcov_0(:,nny+1,:,:))
cc        v0(:,0    ,:,:) = 0.5*(jcov_0(:,1    ,:,:)
cc     .                        +jcov_0(:,0    ,:,:))
cc
cc        v0(:,:,nnz+1,:) = 0.5*(jcov_0(:,:,nnz  ,:)
cc     .                        +jcov_0(:,:,nnz+1,:))
cc        v0(:,:,0    ,:) = 0.5*(jcov_0(:,:,1    ,:)
cc     .                        +jcov_0(:,:,0    ,:))
cc
cc      endif
cc
ccc     Fill ghost nodes
cc
cc      call setBC(IJX,3,nnx,nny,nnz,jcnv,jcov,v0,bcnd,iigx,iigy,iigz
cc     .          ,iorder=order)
cc
ccc Velocity BC
cc
ccc     BC setup
cc
cc      bcnd(:,1) = varray%array_var(IVX)%bconds
cc      bcnd(:,2) = varray%array_var(IVY)%bconds
cc      bcnd(:,3) = varray%array_var(IVZ)%bconds
cc
cc      if (.not.nc_eom_v) then
cc        where (varray%array_var(IRHO)%array /= 0d0)
cc          vcnv(:,:,:,1) = varray%array_var(IVX )%array
cc     .                   /varray%array_var(IRHO)%array
cc          vcnv(:,:,:,2) = varray%array_var(IVY )%array
cc     .                   /varray%array_var(IRHO)%array
cc          vcnv(:,:,:,3) = varray%array_var(IVZ )%array
cc     .                   /varray%array_var(IRHO)%array
cc        end where
cc      else
cc        vcnv(:,:,:,1) = varray%array_var(IVX)%array
cc        vcnv(:,:,:,2) = varray%array_var(IVY)%array
cc        vcnv(:,:,:,3) = varray%array_var(IVZ)%array
cc      endif
cc
cc      v0(:,:,:,1) = u_0%array_var(IVX)%array
cc      v0(:,:,:,2) = u_0%array_var(IVY)%array
cc      v0(:,:,:,3) = u_0%array_var(IVZ)%array
cc
cc      call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv,bcov,v0)
cc
ccc     Fill ghost nodes
cc
cc      call setBC(IVX,3,nnx,nny,nnz,vcnv,vcov,v0,bcnd,iigx,iigy,iigz
cc     .          ,iorder=order)
cc
ccc     Postprocessing
cc
cc      if (.not.nc_eom_v) then
cc        varray%array_var(IVX)%array = vcnv(:,:,:,1)
cc     .                               *varray%array_var(IRHO)%array
cc        varray%array_var(IVY)%array = vcnv(:,:,:,2)
cc     .                               *varray%array_var(IRHO)%array
cc        varray%array_var(IVZ)%array = vcnv(:,:,:,3)
cc     .                               *varray%array_var(IRHO)%array
cc      else
cc        varray%array_var(IVX)%array = vcnv(:,:,:,1)
cc        varray%array_var(IVY)%array = vcnv(:,:,:,2)
cc        varray%array_var(IVZ)%array = vcnv(:,:,:,3)
cc      endif
cc
ccc Temperature BCs
cc
cc      call setBC(ITMP,nnx,nny,nnz
cc     .               ,varray%array_var(ITMP)%array
cc     .               ,u_0   %array_var(ITMP)%array
cc     .               ,varray%array_var(ITMP)%bconds
cc     .               ,iigx,iigy,iigz,iorder=order)
cc
ccc Deallocate variables 
cc
cc      deallocate(v0)
cc
ccc End
cc
cc      end subroutine imposeBoundaryConditions
cc
cc#else

#if defined(samrai)

c applybc
c######################################################################
      subroutine applybc(i,patch_var)

c----------------------------------------------------------------------
c     Applies the boundary conditions for SAMRAI
c----------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer :: i,obc
      type(patch), TARGET :: patch_var

c Local variables

      integer :: order

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

c Interfaces

      INTERFACE
        subroutine applyDepVarBC(bc_inf,varray,vaux,iigx,iigy,iigz
     .                          ,order)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine applyDepVarBC
      END INTERFACE

      INTERFACE
        subroutine applyAuxVarBC(bc_inf,varray,vaux,iigx,iigy,iigz
     .                          ,order)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine applyAuxVarBC
      END INTERFACE

c Begin program

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      varray => gv%u_n
      vaux   => gv%aux
      u_0    => gv%u_0

c Order of extrapolation in BCs (bc_order defined in imposeBC_mod.F)

      order = bc_order     

      call allocAuxVariables

c Set boundary conditions

      if (gv%bc_seq(i+1,1)*AUX < 0) then     !Process dep. variable
        call applyDepVarBC(gv%bc_seq(i+1,:),varray,vaux,1,1,1
     .                    ,order=order)
      elseif (gv%bc_seq(i+1,1)*AUX > 0) then !Process aux. variable
        call applyAuxVarBC(gv%bc_seq(i+1,:),varray,vaux,1,1,1
     .                    ,order=order)
      endif

c End program

      end subroutine applybc

#else

c imposeBoundaryConditions
c####################################################################
      subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz,obc)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer    :: iigx,iigy,iigz
      integer,optional :: obc

      type(var_array),pointer :: varray

c Local variables

      integer :: i,order

c Interfaces

      INTERFACE
        subroutine applyDepVarBC(bc_inf,varray,vaux,iigx,iigy,iigz
     .                          ,order)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine applyDepVarBC
      END INTERFACE

      INTERFACE
        subroutine applyAuxVarBC(bc_inf,varray,vaux,iigx,iigy,iigz
     .                          ,order)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine applyAuxVarBC
      END INTERFACE

c Begin program

      if (PRESENT(obc)) then
        order = obc
      else
        order = bc_order     !Order of extrapolation in BCs (defined in imposeBC_mod.F)
      endif

c Set application BCs (see imposeBC_mod.F)

      call setup_app_BC(itime)

c Impose BCs

      do i = 1,gv%nbc_seq
        if (gv%bc_seq(i,1)*AUX < 0) then      !Process dep. variable
          call applyDepVarBC(gv%bc_seq(i,:),varray,gv%aux
     .                      ,iigx,iigy,iigz,order=order)
        elseif (gv%bc_seq(i,1)*AUX > 0) then  !Process aux. variable
          call applyAuxVarBC(gv%bc_seq(i,:),varray,gv%aux
     .                      ,iigx,iigy,iigz,order=order)
        endif
      enddo

c End program

      end subroutine imposeBoundaryConditions
#endif

c applyDepVarBC
c####################################################################
      subroutine applyDepVarBC(bc_inf,varray,vaux,iigx,iigy,iigz,order)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer    :: iigx,iigy,iigz,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

      integer,optional :: order

c Local variables

      integer    :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz

cc      real(8)    :: bbcov(3),bbcnv(3),b2,J0(3),E00(3)

      real(8),allocatable,dimension(:,:,:,:) :: v0

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

c Allocate auxiliary variables in local domain

      allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,3))

c Select variable

      select case(abs(bc_inf(1)))
      case(IRHO)

c     Density BCs

        if (bc_inf(2) /= SCALAR) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be scalar (0)')
        endif

        call setBC(IRHO,nnx,nny,nnz
     .               ,varray%array_var(IRHO)%array
     .               ,u_0   %array_var(IRHO)%array
     .               ,varray%array_var(IRHO)%bconds
     .               ,iigx,iigy,iigz,iorder=order)

#if defined(vec_pot)
      case(IAX)

c     Vector Pot. BC

        if (bc_inf(2) /= VECTOR) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be vector (1)')
        endif

c       BC setup

        bcnd(:,1) = varray%array_var(IAX)%bconds
        bcnd(:,2) = varray%array_var(IAY)%bconds
        bcnd(:,3) = varray%array_var(IAZ)%bconds

        acov(:,:,:,1) = varray%array_var(IAX)%array
        acov(:,:,:,2) = varray%array_var(IAY)%array
        acov(:,:,:,3) = varray%array_var(IAZ)%array

        v0(:,:,:,1) = u_0%array_var(IAX)%array
        v0(:,:,:,2) = u_0%array_var(IAY)%array
        v0(:,:,:,3) = u_0%array_var(IAZ)%array

c       Fill ghost nodes

        call setBC(IAX,3,nnx,nny,nnz,acnv,acov,v0,bcnd,iigx,iigy,iigz
     .            ,is_cnv=.false.,iorder=order)

        varray%array_var(IAX)%array = acov(:,:,:,1)
        varray%array_var(IAY)%array = acov(:,:,:,2)
        varray%array_var(IAZ)%array = acov(:,:,:,3)

c       Postprocessing: B=curl(A)

        do k = 0,nnz+1
          do j = 0,nny+1
            do i = 0,nnx+1
              do icomp=1,3
                bcnv(i,j,k,icomp) = B0(icomp)
     .                 + curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                        ,acov(:,:,:,1)
     .                        ,acov(:,:,:,2)
     .                        ,acov(:,:,:,3),icomp)
              enddo
cc                call transformFromCurvToCurv(i,j,k,iigx,iigy,iigz
cc     .            ,b_cov(i,j,k,1),b_cov(i,j,k,2),b_cov(i,j,k,3)
cc     .            ,b_cnv(i,j,k,1),b_cnv(i,j,k,2),b_cnv(i,j,k,3)
cc     .            ,.false.)

            enddo
          enddo
        enddo

        if (save_Bn) bcnv_n = bcnv
      endif

#else
      case(IBX)

c     Magnetic field BC

        if (bc_inf(2) /= VECTOR) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be vector (1)')
        endif

c       BC setup

        bcnd(:,1) = varray%array_var(IBX)%bconds
        bcnd(:,2) = varray%array_var(IBY)%bconds
        bcnd(:,3) = varray%array_var(IBZ)%bconds

        bcnv(:,:,:,1) = varray%array_var(IBX)%array
        bcnv(:,:,:,2) = varray%array_var(IBY)%array
        bcnv(:,:,:,3) = varray%array_var(IBZ)%array

        if (pinch_flow) then
          v0 = 0d0
        else
          v0 = jcov_0             !To fix current at boundary
        endif

        if (save_Bn) bcnv_n = bcnv  !Save without BCs

c       Fill ghost nodes

        call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd,iigx,iigy,iigz
     .            ,iorder=order)

        varray%array_var(IBX)%array = bcnv(:,:,:,1)
        varray%array_var(IBY)%array = bcnv(:,:,:,2)
        varray%array_var(IBZ)%array = bcnv(:,:,:,3)

        if (save_Bn) bcnv_n = bcnv  !Save with updated BCs

c       Postprocessing: J=curl(B)

        do k = 0,nnz+1
          do j = 0,nny+1
            do i = 0,nnx+1
              do icomp=1,3
                jcnv(i,j,k,icomp)=
     .                     curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                          ,bcov(:,:,:,1)
     .                          ,bcov(:,:,:,2)
     .                          ,bcov(:,:,:,3),icomp)
              enddo
cc                call transformFromCurvToCurv(i,j,k,iigx,iigy,iigz
cc     .            ,j_cov(i,j,k,1),j_cov(i,j,k,2),j_cov(i,j,k,3)
cc     .            ,j_cnv(i,j,k,1),j_cnv(i,j,k,2),j_cnv(i,j,k,3)
cc     .            ,.false.)
            enddo
          enddo
        enddo
#endif

      case(IVX)

c     Velocity BC

        if (bc_inf(2) /= VECTOR) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be vector (1)')
        endif

c       BC setup

        bcnd(:,1) = varray%array_var(IVX)%bconds
        bcnd(:,2) = varray%array_var(IVY)%bconds
        bcnd(:,3) = varray%array_var(IVZ)%bconds

        if (.not.nc_eom_v) then
          where (varray%array_var(IRHO)%array /= 0d0)
            vcnv(:,:,:,1) = varray%array_var(IVX )%array
     .                     /varray%array_var(IRHO)%array
            vcnv(:,:,:,2) = varray%array_var(IVY )%array
     .                     /varray%array_var(IRHO)%array
            vcnv(:,:,:,3) = varray%array_var(IVZ )%array
     .                     /varray%array_var(IRHO)%array
          end where
        else
          vcnv(:,:,:,1) = varray%array_var(IVX)%array
          vcnv(:,:,:,2) = varray%array_var(IVY)%array
          vcnv(:,:,:,3) = varray%array_var(IVZ)%array
        endif

        v0(:,:,:,1) = u_0%array_var(IVX)%array
        v0(:,:,:,2) = u_0%array_var(IVY)%array
        v0(:,:,:,3) = u_0%array_var(IVZ)%array

        if (rfx_pinch_bc) then
          call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv_n,v0)
        else
          call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv,v0)
        endif

c       Fill ghost nodes

        call setBC(IVX,3,nnx,nny,nnz,vcnv,vcov,v0,bcnd,iigx,iigy,iigz
     .            ,iorder=order)

c       Postprocessing

        if (.not.nc_eom_v) then
          varray%array_var(IVX)%array = vcnv(:,:,:,1)
     .                                 *varray%array_var(IRHO)%array
          varray%array_var(IVY)%array = vcnv(:,:,:,2)
     .                                 *varray%array_var(IRHO)%array
          varray%array_var(IVZ)%array = vcnv(:,:,:,3)
     .                                 *varray%array_var(IRHO)%array
        else
          varray%array_var(IVX)%array = vcnv(:,:,:,1)
          varray%array_var(IVY)%array = vcnv(:,:,:,2)
          varray%array_var(IVZ)%array = vcnv(:,:,:,3)
        endif

      case(ITMP)

c     Temperature BCs

        if (bc_inf(2) /= SCALAR) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be scalar (0)')
        endif

        call setBC(ITMP,nnx,nny,nnz
     .               ,varray%array_var(ITMP)%array
     .               ,u_0   %array_var(ITMP)%array
     .               ,varray%array_var(ITMP)%bconds
     .               ,iigx,iigy,iigz,iorder=order)

      end select

c Deallocate variables 

      deallocate(v0)

c End

      end subroutine applyDepVarBC

c applyAuxVarBC
c####################################################################
      subroutine applyAuxVarBC(bc_inf,varray,vaux,iigx,iigy,iigz,order)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

c diag ****
      use app_iosetup
c diag ****

      implicit none

c Call variables

      integer    :: iigx,iigy,iigz,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

      integer,optional :: order

c Local variables

      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),etal

      real(8),allocatable,dimension(:,:,:,:) :: v0

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

        select case(abs(bc_inf(1)))
        case(IETA)

c       Resistivity (needed for BCs; may depend on T)

          do k=0,nnz+1
            do j=0,nny+1
              do i=0,nnx+1
                eeta(i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
              enddo
            enddo
          enddo

        case default
          call pstop('applyDepVarBC'
     .              ,'Scalar auxiliary variable does not exist')
        end select

      case(VECTOR)

c       Allocate auxiliary variables in local domain

        allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,3))

        select case(abs(bc_inf(1)))
        case(IJX,IJCNV)

c       Current BC

          if (bc_inf(2) /= 1) then
            call pstop('applyDepVarBC'
     .                ,'Inconsistent qty def: should be vector (1)')
          endif

          bcnd = vaux%vec_list(IJCNV)%bconds

          if (pinch_flow) then  !Jt = alpha*Bt in cov representation

            i = nnx

            do k=1,nnz
              do j=1,nny
                if (rfx_pinch_bc) then
                  bbcnv = 0.5*(bcnv_n(i+1,j,k,:)+bcnv_n(i,j,k,:))
                  bbcov =
     .               0.5*(matmul(gmetric%grid(iigx)%gsub(i  ,j,k,:,:)
     .                          ,bcnv_n(i  ,j,k,:))
     .                   +matmul(gmetric%grid(iigx)%gsub(i+1,j,k,:,:)
     .                          ,bcnv_n(i+1,j,k,:)))

                else
                  bbcnv = 0.5*(bcnv(i+1,j,k,:)+bcnv(i,j,k,:))
                  bbcov = 0.5*(bcov(i+1,j,k,:)+bcov(i,j,k,:))
                endif

                J0 = 0.5*(jcov_0(i+1,j,k,:)+jcov_0(i,j,k,:))

                v0(i+1,j,k,:) = dot_product(J0   ,bbcnv)
     .                         /dot_product(bbcov,bbcnv)*bbcov
              enddo
            enddo

          else

            v0(nnx+1,:,:,:) = 0.5*(jcov_0(nnx  ,:,:,:)
     .                            +jcov_0(nnx+1,:,:,:))
            v0(0    ,:,:,:) = 0.5*(jcov_0(1    ,:,:,:)
     .                            +jcov_0(0    ,:,:,:))

            v0(:,nny+1,:,:) = 0.5*(jcov_0(:,nny  ,:,:)
     .                            +jcov_0(:,nny+1,:,:))
            v0(:,0    ,:,:) = 0.5*(jcov_0(:,1    ,:,:)
     .                            +jcov_0(:,0    ,:,:))

            v0(:,:,nnz+1,:) = 0.5*(jcov_0(:,:,nnz  ,:)
     .                            +jcov_0(:,:,nnz+1,:))
            v0(:,:,0    ,:) = 0.5*(jcov_0(:,:,1    ,:)
     .                            +jcov_0(:,:,0    ,:))

          endif

c         Fill ghost nodes

          call setBC(IJX,3,nnx,nny,nnz,jcnv,jcov,v0,bcnd,iigx,iigy,iigz
     .              ,iorder=order)

        case(IJ0CNV)

          v0 = 0d0

c       Current BC

          if (bc_inf(2) /= 1) then
            call pstop('applyDepVarBC'
     .                ,'Inconsistent qty def: should be vector (1)')
          endif

c diag ****
          if (test) then

            jcov_0 = 0d0
            jcnv_0 = 0d0
c diag ****

          elseif (sum(E0) /= 0d0) then

            !Enforce imposed E field
            do k=0,nnz+1
              do j=0,nny+1
                do i=0,nnx+1
                  etal = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
cc                  etal = eeta(i,j,k)
                  jcov_0(i,j,k,:) = E0/etal
                  jcnv_0(i,j,k,:) =
     .                 matmul(gmetric%grid(iigx)%gsup(i,j,k,:,:)
     .                       ,jcov_0(i,j,k,:))
                enddo
              enddo
            enddo

          else

            bcnd = vaux%vec_list(IJ0CNV)%bconds

c           Find B0

#if defined(vec_pot)
            acov(:,:,:,1) = varray%array_var(IAX)%array
            acov(:,:,:,2) = varray%array_var(IAY)%array
            acov(:,:,:,3) = varray%array_var(IAZ)%array

            call setBC(IAX,3,nnx,nny,nnz,acnv,acov,v0,bcnd
     .                ,iigx,iigy,iigz,is_cnv=.false.,iorder=2)

            do k = 0,nnz+1
              do j = 0,nny+1
                do i = 0,nnx+1
                  do icomp=1,3
                    bcnv(i,j,k,icomp)=B0(icomp)
     .                               +curl2(i,j,k,nnx,nny,nnz
     .                                      ,iigx,iigy,iigz
     .                                      ,acov(:,:,:,1)
     .                                      ,acov(:,:,:,2)
     .                                      ,acov(:,:,:,3),icomp)
                  enddo
                enddo
              enddo
            enddo
#else
            bcnv(:,:,:,1) = varray%array_var(IBX)%array
            bcnv(:,:,:,2) = varray%array_var(IBY)%array
            bcnv(:,:,:,3) = varray%array_var(IBZ)%array
#endif

            call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd
     .              ,iigx,iigy,iigz,is_cnv=.true.,iorder=2)

c           Find J0

            do k = 0,nnz+1
              do j = 0,nny+1
                do i = 0,nnx+1
                  do icomp=1,3
                    jcnv_0(i,j,k,icomp)=curl2(i,j,k,nnx,nny,nnz
     .                                      ,iigx,iigy,iigz
     .                                      ,bcov(:,:,:,1)
     .                                      ,bcov(:,:,:,2)
     .                                      ,bcov(:,:,:,3),icomp)
                  enddo
                enddo
              enddo
            enddo

            call setBC(IJX,3,nnx,nny,nnz,jcnv_0,jcov_0,v0,bcnd
     .                ,iigx,iigy,iigz,is_cnv=.true.,iorder=2)

cc            call setBC(IJX,3,nnx,nny,nnz,jcnv_0,jcov_0,v0,bcnd
cc     .                ,iigx,iigy,iigz,iorder=order)
          endif

#if defined(vec_pot)
        case(IBX)

c       Magnetic field BC (with fixed current at boundary)

          if (bc_inf(2) /= 1) then
            call pstop('applyDepVarBC'
     .                ,'Inconsistent qty def: should be vector (1)')
          endif

          bcnd = vaux%vec_list(IBCNV)%bconds

          call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,jcov_0,bcnd
     .              ,iigx,iigy,iigz,iorder=order)

c         Postprocessing: J=curl(B)=curl(curl(A))

          do k = 0,nnz+1
            do j = 0,nny+1
              do i = 0,nnx+1
                do icomp=1,3
                  jcnv(i,j,k,icomp)=
     .                       curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                            ,bcov(:,:,:,1)
     .                            ,bcov(:,:,:,2)
     .                            ,bcov(:,:,:,3),icomp)
                enddo
cc              jcnv(i,j,k,:)= matmul(gmetric%grid(iigx)%gsup(i,j,k,:,:),
cc     .             -veclap_cov(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,acov
cc     .                        ,vol=.false.)
cc     .                +graddiv(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,acnv
cc     .                        ,vol=.false.))

cc              call transformFromCurvToCurv(i,j,k,iigx,iigy,iigz
cc     .            ,j_cov(i,j,k,1),j_cov(i,j,k,2),j_cov(i,j,k,3)
cc     .            ,j_cnv(i,j,k,1),j_cnv(i,j,k,2),j_cnv(i,j,k,3)
cc     .            ,.false.)

              enddo
            enddo
          enddo

#endif
        case default
          call pstop('applyDepVarBC'
     .              ,'Vector auxiliary variable does not exist')

        end select

c       Deallocate variables 

        deallocate(v0)

      end select

c End

      end subroutine applyAuxVarBC

cc#endif
