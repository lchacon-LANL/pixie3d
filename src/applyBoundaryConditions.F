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
      bbcs(:,IAX) = bcsq

      bcsq = bbcs(:,IAY)
      where (bcsq == DEF) bcsq = -EQU  !On covariant components
      bbcs(:,IAY) = bcsq

      bcsq = bbcs(:,IAZ)
      where (bcsq == DEF) bcsq = -EQU  !On covariant components
      bbcs(:,IAZ) = bcsq
#endif

      bcsq = bbcs(:,ITMP)
      where (bcsq == DEF) bcsq = NEU
cc      if (gamma /= 1d0 .and. (.not.adiabatic)) then
cc        where (bcsq == DEF) bcsq = DIR
cc      else
cc        where (bcsq == DEF) bcsq = NEU !To allow isothermal case
cc      endif
      bbcs(:,ITMP) = bcsq

c Exceptions for specific equilibria

      select case (equil)
      case ('tmcyz','3dkh')

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = -NEU
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) =  DIR

        where (bbcs(:,IBX) == DIR) bbcs(:,IBX) = -NEU
        where (bbcs(:,IBY) ==-NEU) bbcs(:,IBY) =  DIR

      case ('ppnch','ppn3d','ppnst','2dhel')   !No stress BCs, pinch flow

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX ) = EQU

      case ('ppnsl','p3nsl','ohtor')           !No slip BCs, pinch flow

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX ) = EQU

        !Dirichlet boundary conditions for tangential components
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY ) = DIR
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ ) = DIR

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
cc        !Dirichlet boundary conditions for tangential velocity components
cc        where (abs(bbcs(:,IVY)) == NEU) bbcs(:,IVY) = DIR
cc        where (abs(bbcs(:,IVZ)) == NEU) bbcs(:,IVZ) = DIR

      case('tmcar','khcar')
cc
cc        where (abs(bbcs(:,IVY )) == NEU) bbcs(:,IVY ) = DIR
cc        where (abs(bbcs(:,IVZ )) == NEU) bbcs(:,IVZ ) = DIR

        where (abs(bbcs(:,ITMP)) == DIR) bbcs(:,ITMP) = NEU

      case('gem')

        !Dirichlet boundary conditions for tangential velocity components
cc        where (abs(bbcs(:,IVY )) == NEU) bbcs(:,IVY ) = DIR
cc        where (abs(bbcs(:,IVZ )) == NEU) bbcs(:,IVZ ) = DIR

        where (abs(bbcs(:,ITMP)) == DIR) bbcs(:,ITMP) = NEU
cc        where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = DIR

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

      use problem_def

      use local_BCS_variables

c     -------------------------------------------------------------------
c     Defines BCs for auxiliary quantities
c     -------------------------------------------------------------------

c     Call variables

      integer    :: bcs(6,neq)

c     Local variables

      integer    :: ieq,bcsq(6)

c     Begin program

      !Flow
      gv%aux%vec_list(IVCNV)%bconds(:,1) = bcs(:,IVX)
      gv%aux%vec_list(IVCNV)%bconds(:,2) = bcs(:,IVY)
      gv%aux%vec_list(IVCNV)%bconds(:,3) = bcs(:,IVZ)

      gv%aux%vec_list(IVCOV)%bconds = gv%aux%vec_list(IVCNV)%bconds

#if defined(vec_pot)
      !Vector potential
      gv%aux%vec_list(IACNV)%bconds(:,1) = bcs(:,IAX)
      gv%aux%vec_list(IACNV)%bconds(:,2) = bcs(:,IAY)
      gv%aux%vec_list(IACNV)%bconds(:,3) = bcs(:,IAZ)

      gv%aux%vec_list(IACOV)%bconds = gv%aux%vec_list(IACNV)%bconds

      !Magnetic field
      gv%aux%vec_list(IBCNV)%bconds = gv%aux%vec_list(IACNV)%bconds
      where (gv%aux%vec_list(IBCNV)%bconds == -EQU)
        gv%aux%vec_list(IBCNV)%bconds = -NEU
      end where

      gv%aux%vec_list(IBCOV)%bconds = gv%aux%vec_list(IBCNV)%bconds
#else
      !Magnetic field
      gv%aux%vec_list(IBCNV)%bconds(:,1) = bcs(:,IBX)
      gv%aux%vec_list(IBCNV)%bconds(:,2) = bcs(:,IBY)
      gv%aux%vec_list(IBCNV)%bconds(:,3) = bcs(:,IBZ)

      gv%aux%vec_list(IBCOV)%bconds = gv%aux%vec_list(IBCNV)%bconds
#endif

      !Current
      gv%aux%vec_list(IJCNV)%bconds = gv%aux%vec_list(IBCNV)%bconds
      where (gv%aux%vec_list(IJCNV)%bconds == -NEU)
        gv%aux%vec_list(IJCNV)%bconds = EQU !Use contravariant components for tangential BCs
      end where

      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds

      !Equilibrium Current (extrapolate covariant components)
      gv%aux%vec_list(IJ0CNV)%bconds = gv%aux%vec_list(IJCNV)%bconds
      where (gv%aux%vec_list(IJ0CNV)%bconds == EQU)
        gv%aux%vec_list(IJ0CNV)%bconds = EXT
      end where
      gv%aux%vec_list(IJ0COV)%bconds = gv%aux%vec_list(IJ0CNV)%bconds

      !Electric field
cc      gv%aux%vec_list(IENI)%bconds = gv%aux%vec_list(IVCNV)%bconds
      gv%aux%vec_list(IENI)%bconds = gv%aux%vec_list(IBCNV)%bconds
      where (gv%aux%vec_list(IENI)%bconds == -NEU)
cc        gv%aux%vec_list(IENI)%bconds =-EQU
cc        gv%aux%vec_list(IENI)%bconds =-EXT
        gv%aux%vec_list(IENI)%bconds =-DIR
      end where

      !Equilibrium electric field
      gv%aux%vec_list(IENI0)%bconds = gv%aux%vec_list(IBCNV)%bconds
      where (gv%aux%vec_list(IENI0)%bconds == -NEU)
        gv%aux%vec_list(IENI0)%bconds =-EXT
      end where

c     End program

      end subroutine defineAuxBCs

#if defined(samrai)

c applybc
c######################################################################
      subroutine applybc(i,patch_var, it)

c----------------------------------------------------------------------
c     Applies the boundary conditions for SAMRAI
c----------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer :: i,obc, it
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

      INTERFACE
         subroutine setup_app_BC(it)
         use local_BCS_variables
         implicit none
         integer :: it
         end subroutine setup_app_BC
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

c Set application BCs (see imposeBC_mod.F)

      call setup_app_BC(it)

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

      integer :: i,order,igrp

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

      INTERFACE
        subroutine preproc_AuxVar(bc_inf,varray,vaux,iigx,iigy,iigz
     .                          ,order)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine preproc_AuxVar
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

cc      do i = 1,gv%nbc_seq
cc        if (gv%bc_seq(i,1)*AUX < 0) then      !Process dep. variable
cc          call applyDepVarBC(gv%bc_seq(i,:),varray,gv%aux
cc     .                      ,iigx,iigy,iigz,order=order)
cc        elseif (gv%bc_seq(i,1)*AUX > 0) then  !Process aux. variable
cc          call applyAuxVarBC(gv%bc_seq(i,:),varray,gv%aux
cc     .                      ,iigx,iigy,iigz,order=order)
cc        endif
cc      enddo

      do igrp = 1,gv%nbc_grp
        !Fill auxiliary quantities
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
            call preproc_AuxVar(gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                         ,gv%aux,iigx,iigy,iigz,order=order)
          endif
        enddo

        !Impose BCs on group
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX < 0) then     !Process dep. variable
            call applyDepVarBC(gv%bc_grp(igrp)%bc_seq(i,:),varray,gv%aux
     .                      ,iigx,iigy,iigz,order=order)
          elseif (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then !Process aux. variable
            call applyAuxVarBC(gv%bc_grp(igrp)%bc_seq(i,:),varray,gv%aux
     .                      ,iigx,iigy,iigz,order=order)
          endif
        enddo
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

      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),grad_p(3)

      real(8),allocatable,dimension(:,:,:,:) :: v0

      real(8),pointer,dimension(:,:,:) :: rho,tmp,prs

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

      if (lag_pinch_bc) then
        rho => u_n%array_var(IRHO)%array
        tmp => u_n%array_var(ITMP)%array
      else
        rho => varray%array_var(IRHO)%array
        tmp => varray%array_var(ITMP)%array
      endif

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

#else
      case(IBX)

        allocate(prs(0:nnx+1,0:nny+1,0:nnz+1))

        if (adiabatic) then
          prs = a_p*tmp
        else
          prs = a_p*rho*tmp
        endif

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

        if (save_Bn) bcnv_n = bcnv  !Save previous time step

        if (pinch_flow) then  !Jt = alpha*Bt in cnv representation
          if (lag_pinch_bc) then
            call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0)
          else
            call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
     .                     ,first_order=.true.)
          endif
        else
          v0 = jcnv_0         !To fix current at boundary
        endif

c       Fill ghost nodes

        call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd,iigx,iigy,iigz
     .            ,iorder=order)

        varray%array_var(IBX)%array = bcnv(:,:,:,1)
        varray%array_var(IBY)%array = bcnv(:,:,:,2)
        varray%array_var(IBZ)%array = bcnv(:,:,:,3)

cc        if (save_Bn) bcnv_n = bcnv  !Save previous time step

        deallocate(prs)
#endif

      case(IVX)

c     Velocity BC

        if (bc_inf(2) /= VECTOR) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be vector (1)')
        endif

        allocate(prs(0:nnx+1,0:nny+1,0:nnz+1))

        if (adiabatic) then
          prs = a_p*tmp
        else
          prs = a_p*rho*tmp
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

        if (pinch_flow) then
          if (lag_pinch_bc) then
            call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0)
          else
            call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0)
          endif
        else
          v0(:,:,:,1) = u_0%array_var(IVX)%array
          v0(:,:,:,2) = u_0%array_var(IVY)%array
          v0(:,:,:,3) = u_0%array_var(IVZ)%array
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

        deallocate(prs)

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

ccc applyAuxVarBC
ccc####################################################################
cc      subroutine applyAuxVarBC(bc_inf,varray,vaux,iigx,iigy,iigz,order)
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
ccc diag ****
cc      use app_iosetup
ccc diag ****
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer    :: iigx,iigy,iigz,bc_inf(2)
cc
cc      type(var_array),pointer :: varray
cc      type(aux_array),pointer :: vaux
cc
cc      integer,optional :: order
cc
ccc Local variables
cc
cc      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz
cc     .          ,ig,jg,kg
cc
cc      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),grad_p(3),etal
cc     .          ,idx,idy,idz,k2,kk,k_par,bnorm,lheta,vol,lvol
cc
cc      real(8),allocatable,dimension(:,:,:,:) :: v0,vdum
cc
cc      real(8),pointer,dimension(:,:,:) :: rho,tmp,prs,bx,by,bz
cc
ccc Begin program
cc
ccc Local grid sizes
cc
cc      nnx = grid_params%nxv(iigx) 
cc      nny = grid_params%nyv(iigy)
cc      nnz = grid_params%nzv(iigz)
cc
ccc Select variable
cc
cc      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
cc      case(SCALAR)
cc
cc        select case(abs(bc_inf(1)))
cc        case(IETA,INU)
cc
cc          ones  = 1d0
cc          zeros = 0d0
cc          vzeros= 0d0
cc
ccc       Resistivity and viscosity (needed for BCs; may depend on T)
cc 
cc          do k=0,nnz+1
cc            do j=0,nny+1
cc              do i=0,nnx+1
cc                eeta(i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,psi)
cc                nuu (i,j,k) = vis(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
cc              enddo
cc            enddo
cc          enddo
cc
cc        case(IHETA)
cc
ccc       Find hyperresistivity (e viscosity)
cc
cc          if (di /= 0d0 .and. heta == 0d0 .and. (.not.test)) then
cc            do k=1,nnz
cc              do j=1,nny
cc                do i=1,nnx
cc                  call getMGmap(i,j,k,iigx,iigy,iigz,ig,jg,kg)
cc
cc                  !kk
cc                  idx  = pi/grid_params%dx(ig)
cc                  if (nnx == 1) idx = 0d0
cc                  idy  = pi/grid_params%dy(jg)
cc                  if (nny == 1) idy = 0d0
cc                  idz  = pi/grid_params%dz(kg)
cc                  if (nnz == 1) idz = 0d0
cc
cc                  k2=vectorNorm(i,j,k,iigx,iigy,iigz,idx,idy,idz,.true.)
cc                  kk=sqrt(k2)
cc
cc                  !k_par
cc                  bnorm = vectorNorm(i,j,k,iigx,iigy,iigz
cc     .                              ,bcnv(i,j,k,1)
cc     .                              ,bcnv(i,j,k,2)
cc     .                              ,bcnv(i,j,k,3)
cc     .                              ,.false.)
cccc     .                          ,bx(i,j,k),by(i,j,k),bz(i,j,k),.false.)
cc
cc                  if (bnorm > 0d0) then
cc                    k_par=scalarProduct(i,j,k,iigx,iigy,iigz,idx,idy,idz
cc     .                                 ,bcnv(i,j,k,1)
cc     .                                 ,bcnv(i,j,k,2)
cc     .                                 ,bcnv(i,j,k,3))**2/bnorm
cc                    k_par = sqrt(k_par)
cccc     .                          ,bx(i,j,k),by(i,j,k),bz(i,j,k))/bnorm
cc                  else
cc                    k_par=0d0
cc                  endif
cc
cc                  heta = heta + 0.1*di/kk/k_par
cc     .                             *gmetric%grid(iigx)%dvol(i,j,k)
cc                  vol  = vol + gmetric%grid(iigx)%dvol(i,j,k)
cc
cc                enddo
cc              enddo
cc            enddo
cc
cc#if defined(petsc) && !defined(samrai)
cc            lheta = heta
cc            call MPI_Allreduce(lheta,heta,1,MPI_DOUBLE_PRECISION
cc     .                        ,MPI_SUM,MPI_COMM_WORLD,mpierr)
cc            lvol = vol
cc            call MPI_Allreduce(lvol,vol,1,MPI_DOUBLE_PRECISION
cc     .                        ,MPI_SUM,MPI_COMM_WORLD,mpierr)
cc#else
cc#if defined(samrai)
cc            call pstop('applyAuxVarBC'
cc     .           ,'Need to fix hyperresistivity calculation in SAMRAI')
cc#endif
cc#endif
cc
cc            heta = heta/vol/di**2
cc
cc            if (my_rank == 0) write (*,*) 'Hyperresistivity=',heta
cc
cc          endif
cc
cc        case default
cc          call pstop('applyDepVarBC'
cc     .              ,'Scalar auxiliary variable does not exist')
cc        end select
cc
cc      case(VECTOR)
cc
cc        if (lag_pinch_bc) then
cc          rho => u_n%array_var(IRHO)%array
cc          tmp => u_n%array_var(ITMP)%array
cc        else
cc          rho => varray%array_var(IRHO)%array
cc          tmp => varray%array_var(ITMP)%array
cc        endif
cc
cc        allocate(prs(0:nnx+1,0:nny+1,0:nnz+1))
cc
cc        if (adiabatic) then
cc          prs = a_p*tmp
cc        else
cc          prs = a_p*rho*tmp
cc        endif
cc
ccc       Allocate auxiliary variables in local domain
cc
cc        allocate(v0  (0:nnx+1,0:nny+1,0:nnz+1,3)
cc     .          ,vdum(0:nnx+1,0:nny+1,0:nnz+1,3))
cc
cc        select case(abs(bc_inf(1)))
cc        case(IJCOV)
cc
cc          do k = 0,nnz+1
cc            do j = 0,nny+1
cc              do i = 0,nnx+1
cc                jcov(i,j,k,:)=matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
cc     .                              ,jcnv(i,j,k,:))
cc              enddo
cc            enddo
cc          enddo
cc
cc        case(IJCNV)
cc
cc          bcnd = vaux%vec_list(IJCNV)%bconds
cc
cc          jcnv=curl(nnx,nny,nnz,iigx,iigy,iigz,bcov)
cc          
cc          if (pinch_flow) then  !Jt = alpha*Bt in cov representation
cc            if (lag_pinch_bc) then
cc              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0)
cc            else
cc              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
cc     .                       ,first_order=.true.)
cc            endif
cc          else
cc            v0(nnx+1,:,:,:) = 0.5*(jcnv_0(nnx  ,:,:,:)
cc     .                            +jcnv_0(nnx+1,:,:,:))
cc            v0(0    ,:,:,:) = 0.5*(jcnv_0(1    ,:,:,:)
cc     .                            +jcnv_0(0    ,:,:,:))
cc
cc            v0(:,nny+1,:,:) = 0.5*(jcnv_0(:,nny  ,:,:)
cc     .                            +jcnv_0(:,nny+1,:,:))
cc            v0(:,0    ,:,:) = 0.5*(jcnv_0(:,1    ,:,:)
cc     .                            +jcnv_0(:,0    ,:,:))
cc
cc            v0(:,:,nnz+1,:) = 0.5*(jcnv_0(:,:,nnz  ,:)
cc     .                            +jcnv_0(:,:,nnz+1,:))
cc            v0(:,:,0    ,:) = 0.5*(jcnv_0(:,:,1    ,:)
cc     .                            +jcnv_0(:,:,0    ,:))
cc          endif
cc
ccc         Fill ghost nodes
cc
cc          call setBC(IJX,3,nnx,nny,nnz,jcnv,jcov,v0,bcnd,iigx,iigy,iigz
cc     .              ,iorder=1)
cc
cc        case(IJ0COV)
cc
cc          do k = 0,nnz+1
cc            do j = 0,nny+1
cc              do i = 0,nnx+1
cc                jcov_0(i,j,k,:)=
cc     .               matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
cc     .                     ,jcnv_0(i,j,k,:))
cc              enddo
cc            enddo
cc          enddo
cc
cc        case(IJ0CNV)
cc
cc          v0 = 0d0
cc
ccc       Current BC
cc
ccc diag ****
cc          if (test .and. equil=='tmcar') then 
cc
cc            jcnv_0 = 0d0
cc
ccc diag ****
cc
cc          elseif (sum(E0) /= 0d0) then
cc
cc            !Enforce imposed E field
cc            do k=0,nnz+1
cc              do j=0,nny+1
cc                do i=0,nnx+1
cc                  etal = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,psi)
cc                  jcov_0(i,j,k,:) = E0/etal
cc                  jcnv_0(i,j,k,:) =
cc     .                 matmul(gmetric%grid(iigx)%gsup(i,j,k,:,:)
cc     .                       ,jcov_0(i,j,k,:))
cc                enddo
cc              enddo
cc            enddo
cc
cc          elseif (sum(E0) == 0d0) then
cc
cc            bcnd = vaux%vec_list(IJ0CNV)%bconds
cc
ccc           Find B0
cc
cc#if defined(vec_pot)
cc            acov(:,:,:,1) = varray%array_var(IAX)%array
cc            acov(:,:,:,2) = varray%array_var(IAY)%array
cc            acov(:,:,:,3) = varray%array_var(IAZ)%array
cc
cc            call setBC(IAX,3,nnx,nny,nnz,acnv,acov,v0,bcnd
cc     .                ,iigx,iigy,iigz,is_cnv=.false.,iorder=2)
cc
cc            do k = 0,nnz+1
cc              do j = 0,nny+1
cc                do i = 0,nnx+1
cc                  bcnv(i,j,k,:) = B0
cc     .                + curl(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,acov)
cc                enddo
cc              enddo
cc            enddo
cc#else
cc            bcnv(:,:,:,1) = varray%array_var(IBX)%array
cc            bcnv(:,:,:,2) = varray%array_var(IBY)%array
cc            bcnv(:,:,:,3) = varray%array_var(IBZ)%array
cc#endif
cc
cc            call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd
cc     .              ,iigx,iigy,iigz,is_cnv=.true.,iorder=2)
cc
ccc           Find J0
cc
cc            do k = 0,nnz+1
cc              do j = 0,nny+1
cc                do i = 0,nnx+1
cc                  jcnv_0(i,j,k,:) = curl(i,j,k,nnx,nny,nnz
cc     .                                  ,iigx,iigy,iigz,bcov)
cc                enddo
cc              enddo
cc            enddo
cc
cc            call setBC(IJX,3,nnx,nny,nnz,jcnv_0,jcov_0,v0,bcnd
cc     .                ,iigx,iigy,iigz,is_cnv=.true.,iorder=1)
cc          endif
cc
cc        case(IVCOV)
cc
ccc$$$          if (.not.nc_eom_v) then
ccc$$$            where (varray%array_var(IRHO)%array /= 0d0)
ccc$$$              vcnv(:,:,:,1) = varray%array_var(IVX )%array
ccc$$$     .                       /varray%array_var(IRHO)%array
ccc$$$              vcnv(:,:,:,2) = varray%array_var(IVY )%array
ccc$$$     .                       /varray%array_var(IRHO)%array
ccc$$$              vcnv(:,:,:,3) = varray%array_var(IVZ )%array
ccc$$$     .                       /varray%array_var(IRHO)%array
ccc$$$            end where
ccc$$$          else
ccc$$$            vcnv(:,:,:,1) = varray%array_var(IVX)%array
ccc$$$            vcnv(:,:,:,2) = varray%array_var(IVY)%array
ccc$$$            vcnv(:,:,:,3) = varray%array_var(IVZ)%array
ccc$$$          endif
ccc$$$
cc          do k = 0,nnz+1
cc            do j = 0,nny+1
cc              do i = 0,nnx+1
cc                vcov(i,j,k,:)=matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
cc     .                              ,vcnv(i,j,k,:))
cc              enddo
cc            enddo
cc          enddo
cc         
cc          if (save_tn .and. ion_hall) vcov_n = vcov
cc
cc        case(IBCOV)
cc
cc          do k = 0,nnz+1
cc            do j = 0,nny+1
cc              do i = 0,nnx+1
cc                bcov(i,j,k,:)=matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
cc     .                              ,bcnv(i,j,k,:))
cc              enddo
cc            enddo
cc          enddo
cc
cc#if defined(vec_pot)
cc       case(IACNV)
cc
cc          do k = 0,nnz+1
cc            do j = 0,nny+1
cc              do i = 0,nnx+1
cc                acnv(i,j,k,:)=matmul(gmetric%grid(iigx)%gsup(i,j,k,:,:)
cc     .                              ,acov(i,j,k,:))
cc              enddo
cc            enddo
cc          enddo
cc
cc        case(IBCNV)
cc
ccc       Magnetic field BC (with fixed current at boundary)
cc
cc          do k = 0,nnz+1
cc            do j = 0,nny+1
cc              do i = 0,nnx+1
cc                bcnv(i,j,k,:) = B0
cc     .                 + curl(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,acov)
cc              enddo
cc            enddo
cc          enddo
cc         
cc          bcnd = vaux%vec_list(IBCNV)%bconds
cc
cc          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation
cc            if (lag_pinch_bc) then
cc             call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0)
cc            else
cc             call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
cc     .                      ,first_order=.true.)
cc            endif
cc          else
cc            v0 = jcnv_0         !To fix current at boundary (average performed in BC routine)
cc          endif
cc
cc          call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd
cc     .              ,iigx,iigy,iigz,iorder=order)
cc#endif
cc
cc        case(IVECOV)
cc
cc          if (di > 0d0) then
cc            vecov(:,:,:,1) = vcov(:,:,:,1)-di*jcov(:,:,:,1)/rho
cc            vecov(:,:,:,2) = vcov(:,:,:,2)-di*jcov(:,:,:,2)/rho
cc            vecov(:,:,:,3) = vcov(:,:,:,3)-di*jcov(:,:,:,3)/rho
cc          else
cc            vecov = vcov
cc          endif
cc
cc          if (save_tn .and. de > 0d0) vecov_n = vecov
cc
cc        case(IVE0CNV,IVECNV)
cc
ccc       Electron velocity BC
cc
cc          !Electron velocity piece
cc          if (di > 0d0) then
ccc           ----------------------------------------------------
cc            !Uncomment this to have j~curl(curl(-di*v)) for div_pe(ve)
cc            do k=1,nnz
cc              do j=1,nny
cc                do i=1,nnx
cc                  vecnv(i,j,k,:) = -di*curlcurl(i,j,k,nnx,nny,nnz
cc     .                                  ,iigx,iigy,iigz,vcnv
cc     .                                  ,vol=.false.)
cccc                  vecnv(i,j,k,:) = di*veclaplacian(i,j,k,nnx,nny,nnz
cccc     .                                  ,iigx,iigy,iigz,vcnv
cccc     .                                  ,vol=.false.)
cc                enddo
cc              enddo
cc            enddo
cc
cc            v0 = vecnv_0
cc
cc            if (abs(bc_inf(1)) == IVE0CNV) then
cc              bcnd = vaux%vec_list(IJ0CNV)%bconds !J BCs (extrapolation)
cc            else
cc              bcnd = vaux%vec_list(IJCNV )%bconds !J BCs (equilibrium)
cc            endif
cc
cc            call setBC(IJX,3,nnx,nny,nnz,vecnv,vdum,v0,bcnd
cc     .                ,iigx,iigy,iigz,iorder=order)
cc
cccc            v0 = 0d0
cccc
cccc            bcnd = vaux%vec_list(IBCNV )%bconds !B BCs (equilibrium)
cccc
cccc            call setBC(IBX,3,nnx,nny,nnz,vecnv,vecov,v0,bcnd
cccc     .                ,iigx,iigy,iigz,iorder=order)
cc
ccc           ----------------------------------------------------
cc            !Uncomment this to have ve=v-di*j/rho for div_pe(ve)
cccc            vecnv = jcnv
ccc           ----------------------------------------------------
cc
cc            vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*vecnv(:,:,:,1)/rho
cc            vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*vecnv(:,:,:,2)/rho
cc            vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*vecnv(:,:,:,3)/rho
cc          else
cc            vecnv = vcnv
cc          endif
cc
cc          if (abs(bc_inf(1)) == IVE0CNV) vecnv_0 = vecnv
cc
cc        case(IENI0)
cc
ccc       Electric field BC (time=0)
cc
cc          E_ni_0 = 0d0
cc
cc          call find_E(nnx,nny,nnz,iigx,iigy,iigz,rho,tmp,E_ni_0)
cc
cc#if !defined(vec_pot)
cc          if (di > 0d0) then
cc            bcnd = vaux%vec_list(IENI0)%bconds !A BCs with extrapolation
cc
cc            call setBC(IAX,3,nnx,nny,nnz,v0,E_ni_0,vdum,bcnd
cccc     .              ,iigx,iigy,iigz,iorder=2,is_cnv=.false.)
cc     .              ,iigx,iigy,iigz,iorder=order,is_cnv=.false.)
cc          endif
cc#endif
cccc          E_ni = E_ni_0
cc
cccc#if !defined(vec_pot)
cccc          if (di > 0d0) then
cccc            bcnd(:,1) = bcond
cccc            bcnd(:,2) = bcond
cccc            bcnd(:,3) = bcond
cccc            where (bcnd == DEF) bcnd = -IFC !Average to face
cccc
cccc            call setBC(IAX,3,nnx,nny,nnz,v0,E_ni_0,E_ni,bcnd
cccc     .              ,iigx,iigy,iigz,iorder=order,is_cnv=.false.)
cccc          endif
cccc#endif
cc        case(IENI)
cc
ccc       Electric field BC
cc
cc          call find_E(nnx,nny,nnz,iigx,iigy,iigz,rho,tmp,E_ni)
cc
cc#if !defined(vec_pot)
cc          if (di > 0d0) then
cc            v0(nnx+1,:,:,:) = 0.5*(E_ni_0(nnx  ,:,:,:)
cc     .                            +E_ni_0(nnx+1,:,:,:))
cc            v0(0    ,:,:,:) = 0.5*(E_ni_0(1    ,:,:,:)
cc     .                            +E_ni_0(0    ,:,:,:))
cc
cc            v0(:,nny+1,:,:) = 0.5*(E_ni_0(:,nny  ,:,:)
cc     .                            +E_ni_0(:,nny+1,:,:))
cc            v0(:,0    ,:,:) = 0.5*(E_ni_0(:,1    ,:,:)
cc     .                            +E_ni_0(:,0    ,:,:))
cc
cc            v0(:,:,nnz+1,:) = 0.5*(E_ni_0(:,:,nnz  ,:)
cc     .                            +E_ni_0(:,:,nnz+1,:))
cc            v0(:,:,0    ,:) = 0.5*(E_ni_0(:,:,1    ,:)
cc     .                            +E_ni_0(:,:,0    ,:))
cc
cc            bcnd = vaux%vec_list(IENI)%bconds
cc
cc            call setBC(IAX,3,nnx,nny,nnz,vdum,E_ni,v0,bcnd
cccc            call setBC(IVX,3,nnx,nny,nnz,vdum,E_ni,v0,bcnd
cc     .              ,iigx,iigy,iigz,iorder=order,is_cnv=.false.)
cccc     .              ,iigx,iigy,iigz,iorder=2,is_cnv=.false.)
cc          endif
cc#endif
cc
cc        case default
cc          call pstop('applyDepVarBC'
cc     .              ,'Vector auxiliary variable does not exist')
cc
cc        end select
cc
ccc       Deallocate variables 
cc
cc        deallocate(v0,vdum,prs)
cc
cc      end select
cc
ccc End
cc
cc      end subroutine applyAuxVarBC

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
     .          ,ig,jg,kg

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),grad_p(3),etal
     .          ,idx,idy,idz,k2,kk,k_par,bnorm,lheta,vol,lvol

      real(8),allocatable,dimension(:,:,:,:) :: v0,vdum

      real(8),pointer,dimension(:,:,:) :: rho,tmp,prs,bx,by,bz

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

        select case(abs(bc_inf(1)))
        case(IETA,INU)

        case(IHETA)

        case default
          call pstop('applyAuxVarBC'
     .              ,'Scalar auxiliary variable does not exist')
        end select

      case(VECTOR)

        if (lag_pinch_bc) then
          rho => u_n%array_var(IRHO)%array
          tmp => u_n%array_var(ITMP)%array
        else
          rho => varray%array_var(IRHO)%array
          tmp => varray%array_var(ITMP)%array
        endif

        allocate(prs(0:nnx+1,0:nny+1,0:nnz+1))

        if (adiabatic) then
          prs = a_p*tmp
        else
          prs = a_p*rho*tmp
        endif

c       Allocate auxiliary variables in local domain

        allocate(v0  (0:nnx+1,0:nny+1,0:nnz+1,3)
     .          ,vdum(0:nnx+1,0:nny+1,0:nnz+1,3))

        select case(abs(bc_inf(1)))
        case(IJCOV)

        case(IJCNV)

          bcnd = vaux%vec_list(IJCNV)%bconds
          
          if (pinch_flow) then  !Jt = alpha*Bt in cov representation
            if (lag_pinch_bc) then
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0)
            else
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
     .                       ,first_order=.true.)
            endif
          else
            v0(nnx+1,:,:,:) = 0.5*(jcnv_0(nnx  ,:,:,:)
     .                            +jcnv_0(nnx+1,:,:,:))
            v0(0    ,:,:,:) = 0.5*(jcnv_0(1    ,:,:,:)
     .                            +jcnv_0(0    ,:,:,:))

            v0(:,nny+1,:,:) = 0.5*(jcnv_0(:,nny  ,:,:)
     .                            +jcnv_0(:,nny+1,:,:))
            v0(:,0    ,:,:) = 0.5*(jcnv_0(:,1    ,:,:)
     .                            +jcnv_0(:,0    ,:,:))

            v0(:,:,nnz+1,:) = 0.5*(jcnv_0(:,:,nnz  ,:)
     .                            +jcnv_0(:,:,nnz+1,:))
            v0(:,:,0    ,:) = 0.5*(jcnv_0(:,:,1    ,:)
     .                            +jcnv_0(:,:,0    ,:))
          endif

c         Fill ghost nodes

          call setBC(IJX,3,nnx,nny,nnz,jcnv,jcov,v0,bcnd,iigx,iigy,iigz
     .              ,iorder=1)

        case(IJ0COV)

        case(IJ0CNV)

c       Current BC

          if (.not.(test .and. equil=='tmcar')
     .        .and.(sum(E0) == 0d0)) then

            bcnd = vaux%vec_list(IJ0CNV)%bconds

            v0 = 0d0

            call setBC(IJX,3,nnx,nny,nnz,jcnv_0,jcov_0,v0,bcnd
     .                ,iigx,iigy,iigz,is_cnv=.true.,iorder=1)
          endif

        case(IVCOV)

        case(IVCOV_N)

          if (save_tn .and. ion_hall) vcov_n = vcov

        case(IBCOV)

#if defined(vec_pot)
        case(IACNV)

        case(IBCNV)
         
          bcnd = vaux%vec_list(IBCNV)%bconds

          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation
            if (lag_pinch_bc) then
             call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0)
            else
             call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
     .                      ,first_order=.true.)
            endif
          else
            v0 = jcnv_0         !To fix current at boundary (average performed in BC routine)
          endif

          call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd
     .              ,iigx,iigy,iigz,iorder=order)
#endif

        case(IVE0CNV)

c       Electron velocity BC

          if (di > 0d0) then
c           ----------------------------------------------------
            !Uncomment this to have j~curl(curl(-di*v)) for div_pe(ve)
            bcnd = vaux%vec_list(IJ0CNV)%bconds !J BCs (extrapolation)

            call setBC(IJX,3,nnx,nny,nnz,vecnv_0,vdum,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order)
c           ----------------------------------------------------
          endif

        case(IVECOV_N)

          if (save_tn .and. de > 0d0) vecov_n = vecov

        case(IVECOV)

        case(IVECNV)

c       Electron velocity BC

          if (di > 0d0) then
c           ----------------------------------------------------
            !Uncomment this to have j~curl(curl(-di*v)) for div_pe(ve)
            v0 = vecnv_0

            bcnd = vaux%vec_list(IJCNV )%bconds !J BCs (equilibrium)

            call setBC(IJX,3,nnx,nny,nnz,vecnv,vdum,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order)
c           ----------------------------------------------------
          endif

        case(IENI0)

c       Electric field BC (time=0)

#if !defined(vec_pot)
          if (di > 0d0) then
            bcnd = vaux%vec_list(IENI0)%bconds !A BCs with extrapolation

            call setBC(IAX,3,nnx,nny,nnz,v0,E_ni_0,vdum,bcnd
     .                ,iigx,iigy,iigz,iorder=order,is_cnv=.false.)
          endif
#endif

        case(IENI)

c       Electric field BC

#if !defined(vec_pot)
          if (di > 0d0) then
            v0(nnx+1,:,:,:) = 0.5*(E_ni_0(nnx  ,:,:,:)
     .                            +E_ni_0(nnx+1,:,:,:))
            v0(0    ,:,:,:) = 0.5*(E_ni_0(1    ,:,:,:)
     .                            +E_ni_0(0    ,:,:,:))

            v0(:,nny+1,:,:) = 0.5*(E_ni_0(:,nny  ,:,:)
     .                            +E_ni_0(:,nny+1,:,:))
            v0(:,0    ,:,:) = 0.5*(E_ni_0(:,1    ,:,:)
     .                            +E_ni_0(:,0    ,:,:))

            v0(:,:,nnz+1,:) = 0.5*(E_ni_0(:,:,nnz  ,:)
     .                            +E_ni_0(:,:,nnz+1,:))
            v0(:,:,0    ,:) = 0.5*(E_ni_0(:,:,1    ,:)
     .                            +E_ni_0(:,:,0    ,:))

            bcnd = vaux%vec_list(IENI)%bconds

            call setBC(IAX,3,nnx,nny,nnz,vdum,E_ni,v0,bcnd
cc            call setBC(IVX,3,nnx,nny,nnz,vdum,E_ni,v0,bcnd
     .              ,iigx,iigy,iigz,iorder=order,is_cnv=.false.)
cc     .              ,iigx,iigy,iigz,iorder=2,is_cnv=.false.)
          endif
#endif

        case default
          write (*,*) abs(bc_inf(1))

          call pstop('applyAuxVarBC'
     .              ,'Vector auxiliary variable does not exist')

        end select

c       Deallocate variables 

        deallocate(v0,vdum,prs)

      end select

c End

      end subroutine applyAuxVarBC

c preproc_AuxVar
c####################################################################
      subroutine preproc_AuxVar(bc_inf,varray,vaux,iigx,iigy,iigz,order)
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
     .          ,ig,jg,kg

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),grad_p(3),etal
     .          ,idx,idy,idz,k2,kk,k_par,bnorm,lheta,vol,lvol

      real(8),allocatable,dimension(:,:,:,:) :: v0

      real(8),pointer,dimension(:,:,:) :: rho,tmp,prs,bx,by,bz

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

        select case(abs(bc_inf(1)))
        case(IETA,INU)

          ones  = 1d0
          zeros = 0d0
          vzeros= 0d0

c       Resistivity and viscosity (needed for BCs; may depend on T)
 
          do k=0,nnz+1
            do j=0,nny+1
              do i=0,nnx+1
                eeta(i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,psi)
                nuu (i,j,k) = vis(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
              enddo
            enddo
          enddo

        case(IHETA)

c       Find hyperresistivity (e viscosity)

          if (di /= 0d0 .and. heta == 0d0 .and. (.not.test)) then
            do k=1,nnz
              do j=1,nny
                do i=1,nnx
                  call getMGmap(i,j,k,iigx,iigy,iigz,ig,jg,kg)

                  !kk
                  idx  = pi/grid_params%dx(ig)
                  if (nnx == 1) idx = 0d0
                  idy  = pi/grid_params%dy(jg)
                  if (nny == 1) idy = 0d0
                  idz  = pi/grid_params%dz(kg)
                  if (nnz == 1) idz = 0d0

                  k2=vectorNorm(i,j,k,iigx,iigy,iigz,idx,idy,idz,.true.)
                  kk=sqrt(k2)

                  !k_par
                  bnorm = vectorNorm(i,j,k,iigx,iigy,iigz
     .                              ,bcnv(i,j,k,1)
     .                              ,bcnv(i,j,k,2)
     .                              ,bcnv(i,j,k,3)
     .                              ,.false.)
cc     .                          ,bx(i,j,k),by(i,j,k),bz(i,j,k),.false.)

                  if (bnorm > 0d0) then
                    k_par=scalarProduct(i,j,k,iigx,iigy,iigz,idx,idy,idz
     .                                 ,bcnv(i,j,k,1)
     .                                 ,bcnv(i,j,k,2)
     .                                 ,bcnv(i,j,k,3))**2/bnorm
                    k_par = sqrt(k_par)
cc     .                          ,bx(i,j,k),by(i,j,k),bz(i,j,k))/bnorm
                  else
                    k_par=0d0
                  endif

                  heta = heta + 0.1*di/kk/k_par
     .                             *gmetric%grid(iigx)%dvol(i,j,k)
                  vol  = vol + gmetric%grid(iigx)%dvol(i,j,k)

                enddo
              enddo
            enddo

#if defined(petsc) && !defined(samrai)
            lheta = heta
            call MPI_Allreduce(lheta,heta,1,MPI_DOUBLE_PRECISION
     .                        ,MPI_SUM,MPI_COMM_WORLD,mpierr)
            lvol = vol
            call MPI_Allreduce(lvol,vol,1,MPI_DOUBLE_PRECISION
     .                        ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#else
#if defined(samrai)
            call pstop('applyAuxVarBC'
     .           ,'Need to fix hyperresistivity calculation in SAMRAI')
#endif
#endif

            heta = heta/vol/di**2

            if (my_rank == 0) write (*,*) 'Hyperresistivity=',heta

          endif

        case default
          call pstop('applyDepVarBC'
     .              ,'Scalar auxiliary variable does not exist')
        end select

      case(VECTOR)

        if (lag_pinch_bc) then
          rho => u_n%array_var(IRHO)%array
          tmp => u_n%array_var(ITMP)%array
        else
          rho => varray%array_var(IRHO)%array
          tmp => varray%array_var(ITMP)%array
        endif

c       Allocate auxiliary variables in local domain

        allocate(v0  (0:nnx+1,0:nny+1,0:nnz+1,3))

        select case(abs(bc_inf(1)))
        case(IJCOV)

          do k = 1,nnz
            do j = 1,nny
              do i = 1,nnx
                jcov(i,j,k,:)=matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
     .                              ,jcnv(i,j,k,:))
              enddo
            enddo
          enddo

        case(IJCNV)

          jcnv=curl(nnx,nny,nnz,iigx,iigy,iigz,bcov)

        case(IJ0COV)

          do k = 1,nnz
            do j = 1,nny
              do i = 1,nnx
                jcov_0(i,j,k,:)=
     .               matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
     .                     ,jcnv_0(i,j,k,:))
              enddo
            enddo
          enddo

        case(IJ0CNV)

          v0 = 0d0

c       Current BC

c diag ****
          if (test .and. equil=='tmcar') then 

            jcnv_0 = 0d0

c diag ****

          elseif (sum(E0) /= 0d0) then

            !Enforce imposed E field
            do k=0,nnz+1
              do j=0,nny+1
                do i=0,nnx+1
                  etal = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,psi)
                  jcov_0(i,j,k,:) = E0/etal
                  jcnv_0(i,j,k,:) =
     .                 matmul(gmetric%grid(iigx)%gsup(i,j,k,:,:)
     .                       ,jcov_0(i,j,k,:))
                enddo
              enddo
            enddo

          elseif (sum(E0) == 0d0) then

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
                  bcnv(i,j,k,:) = B0
     .                + curl(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,acov)
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

            do k = 1,nnz
              do j = 1,nny
                do i = 1,nnx
                  jcnv_0(i,j,k,:) = curl(i,j,k,nnx,nny,nnz
     .                                  ,iigx,iigy,iigz,bcov)
                enddo
              enddo
            enddo

          endif

        case(IVCOV_N)

        case(IVCOV)

          do k = 1,nnz
            do j = 1,nny
              do i = 1,nnx
                vcov(i,j,k,:)=matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
     .                              ,vcnv(i,j,k,:))
              enddo
            enddo
          enddo

        case(IBCOV)

          do k = 1,nnz
            do j = 1,nny
              do i = 1,nnx
                bcov(i,j,k,:)=matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
     .                              ,bcnv(i,j,k,:))
              enddo
            enddo
          enddo

#if defined(vec_pot)
       case(IACNV)

          do k = 1,nnz
            do j = 1,nny
              do i = 1,nnx
                acnv(i,j,k,:)=matmul(gmetric%grid(iigx)%gsup(i,j,k,:,:)
     .                              ,acov(i,j,k,:))
              enddo
            enddo
          enddo

        case(IBCNV)

c       Magnetic field BC (with fixed current at boundary)

          do k = 1,nnz
            do j = 1,nny
              do i = 1,nnx
                bcnv(i,j,k,:) = B0
     .                 + curl(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,acov)
              enddo
            enddo
          enddo
#endif

        case(IVECOV_N)

        case(IVECOV)

          if (di > 0d0) then
            vecov(:,:,:,1) = vcov(:,:,:,1)-di*jcov(:,:,:,1)/rho
            vecov(:,:,:,2) = vcov(:,:,:,2)-di*jcov(:,:,:,2)/rho
            vecov(:,:,:,3) = vcov(:,:,:,3)-di*jcov(:,:,:,3)/rho
          else
            vecov = vcov
          endif

        case(IVE0CNV,IVECNV)

c       Electron velocity BC

          !Electron velocity piece
          if (di > 0d0) then
c           ----------------------------------------------------
            !Uncomment this to have j~curl(curl(-di*v)) for div_pe(ve)
            do k=1,nnz
              do j=1,nny
                do i=1,nnx
                  vecnv(i,j,k,:) = -di*curlcurl(i,j,k,nnx,nny,nnz
     .                                  ,iigx,iigy,iigz,vcnv
     .                                  ,vol=.false.)
cc                  vecnv(i,j,k,:) = di*veclaplacian(i,j,k,nnx,nny,nnz
cc     .                                  ,iigx,iigy,iigz,vcnv
cc     .                                  ,vol=.false.)
                enddo
              enddo
            enddo

c           ----------------------------------------------------
            !Uncomment this to have ve=v-di*j/rho for div_pe(ve)
cc            vecnv = jcnv
c           ----------------------------------------------------

            vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*vecnv(:,:,:,1)/rho
            vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*vecnv(:,:,:,2)/rho
            vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*vecnv(:,:,:,3)/rho

          else
            vecnv = vcnv
          endif

          if (abs(bc_inf(1)) == IVE0CNV) vecnv_0 = vecnv

        case(IENI0)

c       Electric field BC (time=0)

          E_ni_0 = 0d0

          call find_E(nnx,nny,nnz,iigx,iigy,iigz,rho,tmp,E_ni_0)

        case(IENI)

c       Electric field BC

          call find_E(nnx,nny,nnz,iigx,iigy,iigz,rho,tmp,E_ni)

        case default
          call pstop('preproc_AuxVar'
     .              ,'Vector auxiliary variable does not exist')

        end select

c       Deallocate variables 

        deallocate(v0)

      end select

c End

      end subroutine preproc_AuxVar
