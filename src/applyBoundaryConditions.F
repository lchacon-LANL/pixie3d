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
      where (bcsq == DEF .or. bcsq == FSYM) bcsq = NEU
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

      call vector_symm_bc(bbcs(:,IVX:IVZ),-1)

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

      call vector_symm_bc(bbcs(:,IBX:IBZ),1)
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

      call vector_symm_bc(bbcs(:,IAX:IAZ),-1)
#endif

      bcsq = bbcs(:,ITMP)
      where (bcsq == DEF .or. bcsq == FSYM) bcsq = NEU
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

cc        if (di > 0d0) then
cc          where (abs(bbcs(:,IVY )) == NEU) bbcs(:,IVY ) = DIR
cc          where (abs(bbcs(:,IVZ )) == NEU) bbcs(:,IVZ ) = DIR
cc        endif

        where (abs(bbcs(:,ITMP)) == DIR) bbcs(:,ITMP) = NEU

      case('gem')

        !Dirichlet boundary conditions for tangential velocity components
        where (abs(bbcs(:,IVY )) == NEU) bbcs(:,IVY ) = DIR
        where (abs(bbcs(:,IVZ )) == NEU) bbcs(:,IVZ ) = DIR

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

      implicit none

c     Call variables

      integer    :: neq,bcs(6,neq)

c     Local variables

      integer    :: ieq,bcsq(6)

c     Begin program

      !Equilibrium density and temperature
      gv%aux%var_list(IRHO_0)%bconds = bcs(:,IRHO)
      gv%aux%var_list(ITMP_0)%bconds = bcs(:,ITMP)

      !Flow
      gv%aux%vec_list(IVCNV)%bconds(:,1) = bcs(:,IVX)
      gv%aux%vec_list(IVCNV)%bconds(:,2) = bcs(:,IVY)
      gv%aux%vec_list(IVCNV)%bconds(:,3) = bcs(:,IVZ)

      gv%aux%vec_list(IVCOV)%bconds = gv%aux%vec_list(IVCNV)%bconds

      gv%aux%vec_list(IVCNV_0)%bconds = gv%aux%vec_list(IVCNV)%bconds

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

      where (bcond == FSYM)
        gv%aux%vec_list(IJCNV)%bconds(:,1) = FSYM !Recover symmetry BCs
        gv%aux%vec_list(IJCNV)%bconds(:,2) = FSYM !Recover symmetry BCs
        gv%aux%vec_list(IJCNV)%bconds(:,3) = FSYM !Recover symmetry BCs
      end where

      where (gv%aux%vec_list(IJCNV)%bconds == -NEU)
        gv%aux%vec_list(IJCNV)%bconds = EQU !Use contravariant components for tangential BCs
      end where

      call vector_symm_bc(gv%aux%vec_list(IJCNV)%bconds,-1)  !Same as velocity

      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds

      !Equilibrium Current (extrapolate covariant components)
      gv%aux%vec_list(IJCNV_0)%bconds = gv%aux%vec_list(IJCNV)%bconds
      where (gv%aux%vec_list(IJCNV_0)%bconds == EQU)
        gv%aux%vec_list(IJCNV_0)%bconds =-EXT
      end where
      gv%aux%vec_list(IJCOV_0)%bconds = gv%aux%vec_list(IJCNV_0)%bconds

      gv%aux%vec_list(IBCNV_0)%bconds = gv%aux%vec_list(IJCNV_0)%bconds

#if defined(vec_pot)
      gv%aux%vec_list(IACOV_0)%bconds = gv%aux%vec_list(IJCNV_0)%bconds
#endif

      !Pressure tensors (extrapolation)
      gv%aux%vec_list(IDIVPI)%bconds(:,1) = bcond
      gv%aux%vec_list(IDIVPI)%bconds(:,2) = bcond
      gv%aux%vec_list(IDIVPI)%bconds(:,3) = bcond
      where (gv%aux%vec_list(IDIVPI)%bconds == DEF)
        gv%aux%vec_list(IDIVPI)%bconds = EXT
      end where

      call vector_symm_bc(gv%aux%vec_list(IDIVPI)%bconds,-1)  !Same as velocity

      gv%aux%vec_list(IDIVPE)%bconds = gv%aux%vec_list(IDIVPI)%bconds

      !Electric field (non-ideal piece)
      gv%aux%vec_list(IENI)%bconds = gv%aux%vec_list(IVCOV)%bconds

cc      gv%aux%vec_list(IENI)%bconds(:,1) = bcond
cc      gv%aux%vec_list(IENI)%bconds(:,2) = bcond
cc      gv%aux%vec_list(IENI)%bconds(:,3) = bcond
cc      where (gv%aux%vec_list(IENI)%bconds == DEF)
cccc        gv%aux%vec_list(IENI)%bconds =-EXT
cc        gv%aux%vec_list(IENI)%bconds =-DIR
cc      end where
cc
cc      call vector_symm_bc(gv%aux%vec_list(IENI)%bconds,-1)  !Same as current

      !Equilibrium electric field (non-ideal piece)
      gv%aux%vec_list(IENI_0)%bconds = gv%aux%vec_list(IENI)%bconds
cc      gv%aux%vec_list(IENI_0)%bconds(:,1) = bcond
cc      gv%aux%vec_list(IENI_0)%bconds(:,2) = bcond
cc      gv%aux%vec_list(IENI_0)%bconds(:,3) = bcond
cc      where (gv%aux%vec_list(IENI_0)%bconds == DEF)
cc        gv%aux%vec_list(IENI_0)%bconds =-EXT
cc      end where
cc
cc      call vector_symm_bc(gv%aux%vec_list(IENI_0)%bconds,-1)  !Same as current

      !Electron velocity
cc      gv%aux%vec_list(IVECNV)%bconds = gv%aux%vec_list(IDIVPE)%bconds
      gv%aux%vec_list(IVECNV)%bconds = gv%aux%vec_list(IJCNV )%bconds
cc      gv%aux%vec_list(IVECNV)%bconds = gv%aux%vec_list(IVCNV )%bconds
      gv%aux%vec_list(IVECOV)%bconds = gv%aux%vec_list(IVECNV)%bconds
cc      gv%aux%vec_list(IVECNV_0)%bconds = gv%aux%vec_list(IVECNV)%bconds
      gv%aux%vec_list(IVECNV_0)%bconds = gv%aux%vec_list(IJCNV_0)%bconds

c     End program

      end subroutine defineAuxBCs

#if defined(samrai)

c applybc
c######################################################################
      subroutine applybc(igrp,patch_var,it)

c----------------------------------------------------------------------
c     Applies the boundary conditions for SAMRAI
c----------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer :: igrp,obc, it
      type(patch), TARGET :: patch_var

c Local variables

      integer :: i
      type(var_array),pointer :: varray => null()

c Interfaces

      INTERFACE
        subroutine applyDepVarBC(bc_inf,varray,iigx,iigy,iigz,order)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
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

cc      INTERFACE
cc         subroutine setup_app_BC(it)
cc         use local_BCS_variables
cc         implicit none
cc         integer :: it
cc         end subroutine setup_app_BC
cc      END INTERFACE

c Begin program

c Setup pointers

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      varray => gv%u_n
      u_0    => gv%u_0

c Set application BCs

cc      10/20/2010: I don't think this call is needed, since it seems to be 
cc      called stand-alone by the SAMRAI driver before this routine
cc      (and initializeAuxVar below, which also needs it).
cc
      call setupVariableInitializationSequence(gv,it)

c Allocate auxiliary variables

      call allocAuxVariables(gv%aux)

c Set boundary conditions

      !Impose BCs on group
      do i = 1,gv%bc_grp(igrp)%nbc_seq
        if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX < 0) then !Process dep. variable
          call applyDepVarBC(gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                      ,1,1,1,order=bc_order)
        elseif (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then !Process aux. variable
          call applyAuxVarBC(gv%bc_grp(igrp)%bc_seq(i,:),varray,gv%aux
     .                      ,1,1,1,order=bc_order)
        endif
      enddo

c Nullify aux pointers

      call deallocAuxVariables(gv%aux)

c End program

      end subroutine applybc

c initializeAuxVar
c######################################################################
      subroutine initializeAuxVar(patch_var,igrp)

c----------------------------------------------------------------------
c     Fills auxiliary variable arrays
c----------------------------------------------------------------------

      use local_BCS_variables
      use parameters
      use grid
      use variable_setup

      implicit none

c Call variables

      type(patch), TARGET :: patch_var
      integer:: igrp,i

c Local variables

      type(var_array),pointer :: varray => null()

      INTERFACE
      subroutine preproc_AuxVar(bc_inf,varray,vaux,iigx,iigy,iigz,order)
        use variable_setup
        integer    :: iigx,iigy,iigz,bc_inf(2)             
        type(var_array),pointer :: varray
        type(aux_array),pointer :: vaux               
        integer,optional :: order
       end subroutine preproc_AuxVar
      END INTERFACE

c Begin program

c Setup pointers

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      varray => gv%u_n
      u_0    => gv%u_0

      call allocAuxVariables(gv%aux)

c Fill auxiliary variable arrays

      do i = 1,gv%bc_grp(igrp)%nbc_seq
        if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
          call preproc_AuxVar(gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                       ,gv%aux,1,1,1,order=bc_order)
        endif
      enddo

c Nullify pointers

      call deallocAuxVariables(gv%aux)

      end subroutine initializeAuxVar

c setupVariableInitializationSequence
c #####################################################################
      subroutine setupVariableInitializationSequence(patch_var,it)
c ---------------------------------------------------------------------
c     Sets up application BC for SAMRAI
c ---------------------------------------------------------------------

      use parameters
      use grid
      use variable_setup
      use local_BCS_variables

      implicit none

c Call variables

      integer :: it

      type(patch), TARGET :: patch_var

c Local variables

c Begin program

      gv => patch_var

c Set BC schedule (see imposeBC_mod.F)

      call setBCsched(it)

c End program

      end subroutine setupVariableInitializationSequence

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
        subroutine applyDepVarBC(bc_inf,varray,iigx,iigy,iigz,order)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
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

      call setBCsched(itime)

c Allocate auxiliary variables

      call allocAuxVariables(gv%aux)

c Impose BCs

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
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX < 0) then     !Process dep. variables
            call applyDepVarBC(gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                        ,iigx,iigy,iigz,order=order)
          elseif (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then !Process aux. variables
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
      subroutine applyDepVarBC(bc_inf,varray,iigx,iigy,iigz,order)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer :: iigx,iigy,iigz,bc_inf(2)

      type(var_array),pointer :: varray

      integer,optional :: order

c Local variables

      integer :: nnx,nny,nnz

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

c Select variable

      select case(abs(bc_inf(1)))
      case(IRHO)

c     Density BCs

        if (bc_inf(2) /= SCALAR) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be scalar (0)')
        endif

        call setBC(IRHO,nnx,nny,nnz
     .            ,varray%array_var(IRHO)%array
     .            ,u_0   %array_var(IRHO)%array
     .            ,varray%array_var(IRHO)%bconds
     .            ,iigx,iigy,iigz,iorder=order)

#if defined(vec_pot)
      case(IAX,IAY,IAZ)

        varray%array_var(IAX)%array = acov(:,:,:,1)
        varray%array_var(IAY)%array = acov(:,:,:,2)
        varray%array_var(IAZ)%array = acov(:,:,:,3)

#else
      case(IBX,IBY,IBZ)

        varray%array_var(IBX)%array = bcnv(:,:,:,1)
        varray%array_var(IBY)%array = bcnv(:,:,:,2)
        varray%array_var(IBZ)%array = bcnv(:,:,:,3)
#endif

      case(IVX,IVY,IVZ)

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
     .            ,varray%array_var(ITMP)%array
     .            ,u_0   %array_var(ITMP)%array
     .            ,varray%array_var(ITMP)%bconds
     .            ,iigx,iigy,iigz,iorder=order)

      end select

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
        case(IETA,INU,IETA_0)

        case(IHETA)

        case(IRHO_0)

c       Density BCs

          if (bc_inf(2) /= SCALAR) then
            call pstop('applyAuxVarBC'
     .                ,'Inconsistent qty def: should be scalar (0)')
          endif

          call setBC(IRHO,nnx,nny,nnz,rho_0
     .              ,u_0%array_var(IRHO)%array
     .              ,vaux%var_list(IRHO_0)%bconds
     .              ,iigx,iigy,iigz,iorder=order)

        case(ITMP_0)

c       Temperature BCs

          if (bc_inf(2) /= SCALAR) then
            call pstop('applyAuxVarBC'
     .                ,'Inconsistent qty def: should be scalar (0)')
          endif

          call setBC(ITMP,nnx,nny,nnz,tmp_0
     .              ,u_0%array_var(ITMP)%array
     .              ,vaux%var_list(ITMP_0)%bconds
     .              ,iigx,iigy,iigz,iorder=order)

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
          
          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation
            if (lag_pinch_bc) then
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0
     .                       ,first_order=.true.)
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

#if defined(vec_pot)
        case(IACOV_0)

          bcnd = vaux%vec_list(IACOV_0)%bconds

          v0 = 0d0

          call setBC(IAX,3,nnx,nny,nnz,acnv,acov_0,v0,bcnd
     .              ,iigx,iigy,iigz
cc     .              ,is_cnv=.false.
     .              ,is_cnv=vaux%vec_list(IACOV_0)%cnv
     .              ,iorder=2)
#endif

        case(IBCOV_0)

        case(IBCNV_0)

          bcnd = vaux%vec_list(IBCNV_0)%bconds

          v0 = 0d0

          call setBC(IBX,3,nnx,nny,nnz,bcnv_0,bcov_0,v0,bcnd
     .              ,iigx,iigy,iigz
cc     .              ,is_cnv=.true.
     .              ,is_cnv=vaux%vec_list(IBCNV_0)%cnv
     .              ,iorder=2)

        case(IJCOV_0)

        case(IJCNV_0)

c       Current BC

          bcnd = vaux%vec_list(IJCNV_0)%bconds

          v0 = 0d0

          call setBC(IJX,3,nnx,nny,nnz,jcnv_0,jcov_0,v0,bcnd
     .              ,iigx,iigy,iigz
cc     .              ,is_cnv=.true.
     .              ,is_cnv=vaux%vec_list(IJCNV_0)%cnv
     .              ,iorder=1)

        case(IJCNV_N)

          if (save_tn .and. ion_hall) jcnv_n = jcnv

        case(IVCOV)

        case(IVCOV_N)

          if (save_tn .and. ion_hall) vcov_n = vcov

        case(IVCNV_N)

          if (save_tn .and. ion_hall) vcnv_n = vcnv

        case(IVCNV)

c       Velocity BC

          if (bc_inf(2) /= VECTOR) then
            call pstop('applyAuxVarBC'
     .                ,'Inconsistent qty def: should be vector (1)')
          endif

c         BC setup

          bcnd = vaux%vec_list(IVCNV)%bconds

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

c         Fill ghost nodes

          call setBC(IVX,3,nnx,nny,nnz,vcnv,vcov,v0,bcnd,iigx,iigy,iigz
     .              ,iorder=order)

        case(IVCNV_0)

c       Velocity BC

          if (bc_inf(2) /= VECTOR) then
            call pstop('applyAuxVarBC'
     .                ,'Inconsistent qty def: should be vector (1)')
          endif

c         BC setup

          bcnd = vaux%vec_list(IVCNV_0)%bconds

          if (adiabatic) then
            prs = a_p*tmp_0
          else
            prs = a_p*rho_0*tmp_0
          endif

          if (pinch_flow) then
            call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv_0,prs,eeta_0,v0)
          else
            v0(:,:,:,1) = u_0%array_var(IVX)%array
            v0(:,:,:,2) = u_0%array_var(IVY)%array
            v0(:,:,:,3) = u_0%array_var(IVZ)%array
          endif

c         Fill ghost nodes

          call setBC(IVX,3,nnx,nny,nnz,vcnv_0,vcov,v0,bcnd
     .              ,iigx,iigy,iigz,iorder=order)

#if defined(vec_pot)

        case(IACOV)

c       Vector Pot. BC

          if (bc_inf(2) /= VECTOR) then
            call pstop('applyAuxVarBC'
     .                ,'Inconsistent qty def: should be vector (1)')
          endif

c         BC setup

          bcnd = vaux%vec_list(IACOV)%bconds

cc          bcnd(:,1) = varray%array_var(IAX)%bconds
cc          bcnd(:,2) = varray%array_var(IAY)%bconds
cc          bcnd(:,3) = varray%array_var(IAZ)%bconds

          v0(:,:,:,1) = u_0%array_var(IAX)%array
          v0(:,:,:,2) = u_0%array_var(IAY)%array
          v0(:,:,:,3) = u_0%array_var(IAZ)%array

c         Fill ghost nodes

          call setBC(IAX,3,nnx,nny,nnz,acnv,acov,v0,bcnd,iigx,iigy,iigz
     .              ,iorder=order,is_cnv=vaux%vec_list(IACOV)%cnv)

        case(IACNV)

cc        case(IBCNV)
cc         
cc          bcnd = vaux%vec_list(IBCNV)%bconds
cc
cc          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation
cc            if (lag_pinch_bc) then
cc             call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0
cc     .                      ,first_order=.true.)
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
cc#else
#endif
        case(IBCNV)

c       Magnetic field BC

          if (bc_inf(2) /= VECTOR) then
            call pstop('applyAuxVarBC'
     .              ,'Inconsistent qty def: should be vector (1)')
          endif

          bcnd = vaux%vec_list(IBCNV)%bconds

          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation

            if (lag_pinch_bc) then
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0
     .                       ,first_order=.true.)
            else
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
     .                       ,first_order=.true.)
            endif

          else
            v0 = vaux%vec_list(IJCNV_0)%vec         !To fix current at boundary
          endif

c       Fill ghost nodes

          call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd,iigx,iigy,iigz
     .              ,iorder=order,is_cnv=vaux%vec_list(IBCNV)%cnv)
cc#endif
        case(IBCOV)

        case(IBCNV_N)

          if (save_tn) bcnv_n = bcnv !Save previous time step

        case(IVECNV_0)

ccc       Electron velocity BC
cc
cc          if (di > 0d0 .and. fake_ve) then
cc            bcnd = vaux%vec_list(IVECNV_0)%bconds
cc
cc            call setBC(IVX,3,nnx,nny,nnz,vecnv_0,vdum,v0,bcnd
cc     .                ,iigx,iigy,iigz,iorder=order)
cc          endif

        case(IVECOV_N)

          if (save_tn .and. de > 0d0) vecov_n = vecov

        case(IVECOV)

        case(IVECNV)

c       Electron velocity BC

          if (fake_ve) then
            !ve=v-di*j/rho with j=-dt*curl(curl(E)) for div_pe(ve)
            v0 = 0d0

            bcnd = vaux%vec_list(IVECNV)%bconds
            call setBC(IJX,3,nnx,nny,nnz,vecnv,vdum,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order
     .                ,is_cnv=vaux%vec_list(IVECNV)%cnv)

            vecnv = vecnv + jcnv_n

            vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*vecnv(:,:,:,1)/rho
            vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*vecnv(:,:,:,2)/rho
            vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*vecnv(:,:,:,3)/rho
          endif

#if !defined(vec_pot)
        case(IDIVPI)

          if (di > 0d0.and.ion_hall) then
            bcnd = vaux%vec_list(IDIVPI)%bconds !A BCs with extrapolation

            call setBC(IVX,3,nnx,nny,nnz,div_pi,vdum,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order
     .                ,is_cnv=vaux%vec_list(IDIVPI)%cnv)
          endif
#endif

#if !defined(vec_pot)
        case(IDIVPE)

          if (di > 0d0.and.(.not.ion_hall)) then
            bcnd = vaux%vec_list(IDIVPE)%bconds !A BCs with extrapolation

            call setBC(IVX,3,nnx,nny,nnz,div_pe,vdum,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order
     .                ,is_cnv=vaux%vec_list(IDIVPE)%cnv)
          endif
#endif

        case(IENI_0)

        case(IENI)

c       Electric field BC

          if (di > 0d0) then
#if !defined(vec_pot)
            bcnd = vaux%vec_list(IENI)%bconds
            v0   = 0d0
            vdum = 0d0
            call setBC(IVX,3,nnx,nny,nnz,vdum,E_h,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order,is_cnv=.false.)
#endif
            E_ni = E_ni + di*E_h
          endif

          E_ni = E_ni - E_ni_0  !Substract imposed electric field

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

        zeros = 0d0
        vzeros = 0d0
        ones = 1d0

        select case(abs(bc_inf(1)))
        case(IETA,INU)

c       Resistivity and viscosity (needed for BCs; may depend on T)
 
          nuu = vis(nnx,nny,nnz,iigx,iigy,iigz)
          eeta= res(nnx,nny,nnz,iigx,iigy,iigz,psi)

        case(IETA_0)

c       Resistivity and viscosity (needed for BCs; may depend on T)
 
          eeta_0=res(nnx,nny,nnz,iigx,iigy,iigz,psi)

        case(IHETA)

c       Find hyperresistivity (electron viscosity)

          if (di /= 0d0 .and. heta == 0d0 .and. (.not.test)) then
            heta=hres(nnx,nny,nnz,iigx,iigy,iigz,bcnv_0,rho_0)
          endif

        case(IRHO_0)

          rho_0 = u_0%array_var(IRHO)%array

        case(ITMP_0)

          tmp_0 = u_0%array_var(ITMP)%array

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

c       Allocate auxiliary variables in local domain

        allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,3))

        select case(abs(bc_inf(1)))
        case(IJCOV)

          do k=0,nnz+1
            do j=0,nny+1
              do i=0,nnx+1
                jcov(i,j,k,:)=matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
     .                              ,jcnv(i,j,k,:))
              enddo
            enddo
          enddo

        case(IJCNV)

          jcnv=curl(nnx,nny,nnz,iigx,iigy,iigz,bcov)

        case(IJCNV_N)

          if (save_tn) jcnv_n = jcnv !Save previous time step

#if defined(vec_pot)
        case(IACOV_0)

          acov_0(:,:,:,1) = u_0%array_var(IAX)%array
          acov_0(:,:,:,2) = u_0%array_var(IAY)%array
          acov_0(:,:,:,3) = u_0%array_var(IAZ)%array
#endif

        case(IBCNV_0)

#if defined(vec_pot)
          do k = 0,nnz+1
            do j = 0,nny+1
              do i = 0,nnx+1
                bcnv_0(i,j,k,:) = B0
     .                + curl(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,acov_0)
              enddo
            enddo
          enddo
#else
          bcnv_0(:,:,:,1) = u_0%array_var(IBX)%array
          bcnv_0(:,:,:,2) = u_0%array_var(IBY)%array
          bcnv_0(:,:,:,3) = u_0%array_var(IBZ)%array
#endif

        case(IBCOV_0)

          do k=0,nnz+1
            do j=0,nny+1
              do i=0,nnx+1
                bcov_0(i,j,k,:)=
     .               matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
     .                     ,bcnv_0(i,j,k,:))
              enddo
            enddo
          enddo

        case(IJCOV_0)

          do k=0,nnz+1
            do j=0,nny+1
              do i=0,nnx+1
                jcov_0(i,j,k,:)=
     .               matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
     .                     ,jcnv_0(i,j,k,:))
              enddo
            enddo
          enddo

        case(IJCNV_0)

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
                  jcnv_0(i,j,k,:) =
     .                 matmul(gmetric%grid(iigx)%gsup(i,j,k,:,:),E0)
     .                /eeta_0(i,j,k)
                enddo
              enddo
            enddo

          elseif (sum(E0) == 0d0) then

cc            bcnd = vaux%vec_list(IJCNV_0)%bconds
cc
ccc           Find B0
cc
cc#if defined(vec_pot)
cc            acov(:,:,:,1) = u_0%array_var(IAX)%array
cc            acov(:,:,:,2) = u_0%array_var(IAY)%array
cc            acov(:,:,:,3) = u_0%array_var(IAZ)%array
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
cc            bcnv(:,:,:,1) = u_0%array_var(IBX)%array
cc            bcnv(:,:,:,2) = u_0%array_var(IBY)%array
cc            bcnv(:,:,:,3) = u_0%array_var(IBZ)%array
cc#endif
cc
cc            call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd
cc     .              ,iigx,iigy,iigz,is_cnv=.true.,iorder=2)

c           Find J0

            jcnv_0 = curl(nnx,nny,nnz,iigx,iigy,iigz,bcov_0)

          endif

        case(IVCOV_N,IVCNV_N)

        case(IVCOV)

          do k = 0,nnz+1
            do j = 0,nny+1
              do i = 0,nnx+1
                vcov(i,j,k,:)=matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
     .                              ,vcnv(i,j,k,:))
              enddo
            enddo
          enddo

        case(IVCNV)

          if (.not.nc_eom_v) then
            where (varray%array_var(IRHO)%array /= 0d0)
              vcnv(:,:,:,1) = varray%array_var(IVX )%array
     .                       /varray%array_var(IRHO)%array
              vcnv(:,:,:,2) = varray%array_var(IVY )%array
     .                       /varray%array_var(IRHO)%array
              vcnv(:,:,:,3) = varray%array_var(IVZ )%array
     .                       /varray%array_var(IRHO)%array
            end where
          else
            vcnv(:,:,:,1) = varray%array_var(IVX)%array
            vcnv(:,:,:,2) = varray%array_var(IVY)%array
            vcnv(:,:,:,3) = varray%array_var(IVZ)%array
          endif

        case(IVCNV_0)

          if (.not.nc_eom_v) then
            where (rho_0 /= 0d0)
              vcnv_0(:,:,:,1) = u_0%array_var(IVX)%array/rho_0
              vcnv_0(:,:,:,2) = u_0%array_var(IVY)%array/rho_0
              vcnv_0(:,:,:,3) = u_0%array_var(IVZ)%array/rho_0
            end where
          else
            vcnv_0(:,:,:,1) = u_0%array_var(IVX)%array
            vcnv_0(:,:,:,2) = u_0%array_var(IVY)%array
            vcnv_0(:,:,:,3) = u_0%array_var(IVZ)%array
          endif

#if defined(vec_pot)

        case(IACOV)

          acov(:,:,:,1) = varray%array_var(IAX)%array
          acov(:,:,:,2) = varray%array_var(IAY)%array
          acov(:,:,:,3) = varray%array_var(IAZ)%array

        case(IACNV)

          do k = 0,nnz+1
            do j = 0,nny+1
              do i = 0,nnx+1
                acnv(i,j,k,:)=matmul(gmetric%grid(iigx)%gsup(i,j,k,:,:)
     .                              ,acov(i,j,k,:))
              enddo
            enddo
          enddo

        case(IBCNV)

c       Magnetic field BC (with fixed current at boundary)

          do k = 0,nnz+1
            do j = 0,nny+1
              do i = 0,nnx+1
                bcnv(i,j,k,:) = B0
     .                 + curl(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,acov)
              enddo
            enddo
          enddo

#else
        case(IBCNV)

          bcnv(:,:,:,1) = varray%array_var(IBX)%array
          bcnv(:,:,:,2) = varray%array_var(IBY)%array
          bcnv(:,:,:,3) = varray%array_var(IBZ)%array

#endif
        case(IBCNV_N)

          if (save_tn) bcnv_n = bcnv !Save previous time step

        case(IBCOV)

          do k = 0,nnz+1
            do j = 0,nny+1
              do i = 0,nnx+1
                bcov(i,j,k,:)=matmul(gmetric%grid(iigx)%gsub(i,j,k,:,:)
     .                              ,bcnv(i,j,k,:))
              enddo
            enddo
          enddo
          
        case(IVECOV_N)

        case(IVECOV)

          if (di > 0d0) then
            vecov(:,:,:,1) = vcov(:,:,:,1)-di*jcov(:,:,:,1)/rho
            vecov(:,:,:,2) = vcov(:,:,:,2)-di*jcov(:,:,:,2)/rho
            vecov(:,:,:,3) = vcov(:,:,:,3)-di*jcov(:,:,:,3)/rho
          else
            vecov = vcov
          endif

        case(IVECNV_0)

c       Electron velocity

          vecnv_0(:,:,:,1) = vcnv_0(:,:,:,1)-di*jcnv_0(:,:,:,1)/rho_0
          vecnv_0(:,:,:,2) = vcnv_0(:,:,:,2)-di*jcnv_0(:,:,:,2)/rho_0
          vecnv_0(:,:,:,3) = vcnv_0(:,:,:,3)-di*jcnv_0(:,:,:,3)/rho_0

        case(IVECNV)

c       Electron velocity

          if (fake_ve) then

            !ve=v-di*j/rho with j=jn-dt*curl(curl(E)) for div_pe(ve)
            if (save_tn) then

              vecnv = 0d0

            else
cc              !Equivalent to fake_ve=.false.
cc              vecnv = jcnv - jcnv_n

cc              !Use full electric field
cc              v0 =-dt*E_ni  !First-order approximation to find j^(n+1)
cc              call transformVector(iigx,iigy,iigz
cc     .                    ,0,nnx+1,0,nny+1,0,nnz+1
cc     .                    ,v0(:,:,:,1)
cc     .                    ,v0(:,:,:,2)
cc     .                    ,v0(:,:,:,3)
cc     .                    ,'cov','cnv')

cc              !Pick certain pieces of electric field: IDEAL
cc              v0(:,:,:,1) = (vcov(:,:,:,2)*bcov(:,:,:,3)
cc     .                      -vcov(:,:,:,3)*bcov(:,:,:,2))
cc              v0(:,:,:,2) = (vcov(:,:,:,3)*bcov(:,:,:,1)
cc     .                      -vcov(:,:,:,1)*bcov(:,:,:,3))
cc              v0(:,:,:,3) = (vcov(:,:,:,1)*bcov(:,:,:,2)
cc     .                      -vcov(:,:,:,2)*bcov(:,:,:,1))

cc              !Resistive
cc              v0(:,:,:,1) = v0(:,:,:,1) + eeta*jcnv(:,:,:,1)
cc              v0(:,:,:,2) = v0(:,:,:,2) + eeta*jcnv(:,:,:,2)
cc              v0(:,:,:,3) = v0(:,:,:,3) + eeta*jcnv(:,:,:,3)

              !Ion inertia
cc              v0 = -v0 + di*(vcnv-vcnv_n)/dt
              v0 =  di*(vcnv-vcnv_n)/dt

              !Ion pressure tensor
              v0(:,:,:,1) = v0(:,:,:,1) + di*div_pi(:,:,:,1)/rho
              v0(:,:,:,2) = v0(:,:,:,2) + di*div_pi(:,:,:,2)/rho
              v0(:,:,:,3) = v0(:,:,:,3) + di*div_pi(:,:,:,3)/rho

              do k=1,nnz
                do j=1,nny
                  do i=1,nnx
cc                    vecnv(i,j,k,:) = -dt*curlcurl(i,j,k,nnx,nny,nnz
cc     .                                           ,iigx,iigy,iigz,v0
cc     .                                           ,vol=.false.)
                    vecnv(i,j,k,:) = dt*veclaplacian(i,j,k,nnx,nny,nnz
     .                                           ,iigx,iigy,iigz,v0
     .                                           ,vol=.false.)
                  enddo
                enddo
              enddo

            endif

          else

           !ve=v-di*j/rho with j=curl(B)
            vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv(:,:,:,1)/rho
            vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv(:,:,:,2)/rho
            vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv(:,:,:,3)/rho

          endif

        case(IDIVPI)

          !Divergence of ion stress tensor
          if (nu > 0d0) then
            div_pi = EOM_divPi(nnx,nny,nnz,iigx,iigy,iigz,vcnv,nuu)
          else
            div_pi = 0d0
          endif

        case(IDIVPE)

          div_pe = EOM_divPe(nnx,nny,nnz,iigx,iigy,iigz,vecnv-vecnv_0)

        case(IENI_0)

c       Electric field (time=0, initial EOM equ => only ideal+resistive contributions)

#if defined(vec_pot)
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,acov_0)
#else
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,jcov_0)
#endif
cc          if (di > 0d0) then
cc            E_h  = E_fld_Hall(nnx,nny,nnz,iigx,iigy,iigz,rho_0,tmp_0)
cc          else
cc            E_h  = 0d0
cc          endif

        case(IENI)

c       Electric field

#if defined(vec_pot)
          E_ni = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,acov)
#else
          E_ni = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,jcov)
#endif
          if (di > 0d0) then
            E_h  = E_fld_Hall(nnx,nny,nnz,iigx,iigy,iigz,rho,tmp)
          else
            E_h  = 0d0
          endif

        case default

          write (*,*) abs(bc_inf(1))

          call pstop('preproc_AuxVar'
     .              ,'Vector auxiliary variable does not exist')

        end select

c       Deallocate variables 

        deallocate(v0)

      end select

c End

      end subroutine preproc_AuxVar
