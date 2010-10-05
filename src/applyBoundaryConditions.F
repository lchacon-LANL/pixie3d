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
      subroutine applybc(igrp,patch_var, it)

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

      integer :: order
      integer :: i
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

c$$$      do igrp = 1,gv%nbc_grp
c$$$        !Fill auxiliary quantities
c$$$        do i = 1,gv%bc_grp(igrp)%nbc_seq
c$$$          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
c$$$            call preproc_AuxVar(gv%bc_grp(igrp)%bc_seq(i,:),varray
c$$$     .                         ,gv%aux,iigx,iigy,iigz,order=order)
c$$$          endif
c$$$        enddo
c$$$
        !Impose BCs on group
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX < 0) then     !Process dep. variable
             call applyDepVarBC(gv%bc_grp(igrp)%bc_seq(i,:),varray,vaux, 
     .       1,1,1,order=order)
          elseif (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then !Process aux. variable
             call applyAuxVarBC(gv%bc_grp(igrp)%bc_seq(i,:),varray,vaux,
     .       1,1,1,order=order)
          endif
        enddo

c End program

      end subroutine applybc

c####################################################
c     fills aux var array
      subroutine initializeAuxVar(patch_var, igrp)

      use local_BCS_variables
cc      use imposeBCinterface
cc      use auxiliaryVariables
      use parameters
      use grid
      use variable_setup
c diag ****
cc      use app_iosetup

      implicit none

      type(patch), TARGET :: patch_var
      integer:: igrp,i
      integer :: order

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

      INTERFACE
      subroutine preproc_AuxVar(bc_inf,varray,vaux,iigx,iigy,iigz,order)
               use local_BCS_variables
               use imposeBCinterface
               use auxiliaryVariables
               use variable_setup
c diag ****
               use app_iosetup
c diag ****
               implicit none
               
               integer    :: iigx,iigy,iigz,bc_inf(2)             
               type(var_array),pointer :: varray
               type(aux_array),pointer :: vaux               
               integer,optional :: order

       end subroutine preproc_AuxVar
      END INTERFACE

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      varray => gv%u_n
      vaux   => gv%aux
      u_0    => gv%u_0

c Order of extrapolation in BCs (bc_order defined in imposeBC_mod.F)

#if defined(samrai)
      call allocAuxVariables
#endif

      order = bc_order

      do i = 1,gv%bc_grp(igrp)%nbc_seq
         if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
            call preproc_AuxVar(gv%bc_grp(igrp)%bc_seq(i,:),varray
     .           ,vaux,1,1,1,order=order)
         endif
      enddo
      
      end subroutine initializeAuxVar

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
            call applyDepVarBC(gv%bc_grp(igrp)%bc_seq(i,:),varray,gv%aux
     .                      ,iigx,iigy,iigz,order=order)
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

      integer :: iigx,iigy,iigz,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

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
      case(IAX)

        varray%array_var(IAX)%array = acov(:,:,:,1)
        varray%array_var(IAY)%array = acov(:,:,:,2)
        varray%array_var(IAZ)%array = acov(:,:,:,3)

#else
      case(IBX)

        varray%array_var(IBX)%array = bcnv(:,:,:,1)
        varray%array_var(IBY)%array = bcnv(:,:,:,2)
        varray%array_var(IBZ)%array = bcnv(:,:,:,3)
#endif

      case(IVX)

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

        case(IJCOV_0)

        case(IJCNV_0)

c       Current BC

          bcnd = vaux%vec_list(IJCNV_0)%bconds

          v0 = 0d0

          call setBC(IJX,3,nnx,nny,nnz,jcnv_0,jcov_0,v0,bcnd
     .              ,iigx,iigy,iigz,is_cnv=.true.,iorder=1)
        case(IJCNV_N)

          if (save_tn .and. ion_hall) jcnv_n = jcnv

        case(IVCOV)

        case(IVCOV_N)

          if (save_tn .and. ion_hall) vcov_n = vcov
cc          if (save_tn .and. ion_hall) then
cc            vcov_n  = vcov
cc            vcov_nm = vcov_n
cc          endif

        case(IVCNV_N)

          if (save_tn .and. ion_hall) vcnv_n = vcnv

        case(IVCNV)

c       Velocity BC

          if (bc_inf(2) /= VECTOR) then
            call pstop('applyDepVarBC'
     .                ,'Inconsistent qty def: should be vector (1)')
          endif

c         BC setup

          bcnd(:,1) = varray%array_var(IVX)%bconds
          bcnd(:,2) = varray%array_var(IVY)%bconds
          bcnd(:,3) = varray%array_var(IVZ)%bconds

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

#if defined(vec_pot)

        case(IACOV)

c       Vector Pot. BC

          if (bc_inf(2) /= VECTOR) then
            call pstop('applyDepVarBC'
     .                ,'Inconsistent qty def: should be vector (1)')
          endif

c         BC setup

          bcnd(:,1) = varray%array_var(IAX)%bconds
          bcnd(:,2) = varray%array_var(IAY)%bconds
          bcnd(:,3) = varray%array_var(IAZ)%bconds

          v0(:,:,:,1) = u_0%array_var(IAX)%array
          v0(:,:,:,2) = u_0%array_var(IAY)%array
          v0(:,:,:,3) = u_0%array_var(IAZ)%array

c         Fill ghost nodes

          call setBC(IAX,3,nnx,nny,nnz,acnv,acov,v0,bcnd,iigx,iigy,iigz
     .              ,is_cnv=.false.,iorder=order)

        case(IACNV)

        case(IBCNV)
         
          bcnd = vaux%vec_list(IBCNV)%bconds

          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation
            if (lag_pinch_bc) then
             call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0
     .                      ,first_order=.true.)
            else
             call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
     .                      ,first_order=.true.)
            endif
          else
            v0 = jcnv_0         !To fix current at boundary (average performed in BC routine)
          endif

          call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd
     .              ,iigx,iigy,iigz,iorder=order)
#else
        case(IBCNV)

c       Magnetic field BC

          if (bc_inf(2) /= VECTOR) then
            call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be vector (1)')
          endif

          bcnd(:,1) = varray%array_var(IBX)%bconds
          bcnd(:,2) = varray%array_var(IBY)%bconds
          bcnd(:,3) = varray%array_var(IBZ)%bconds

          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation

            if (lag_pinch_bc) then
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0
     .                       ,first_order=.true.)
            else
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
     .                       ,first_order=.true.)
            endif

          else
            v0 = jcnv_0         !To fix current at boundary
          endif

c       Fill ghost nodes

        call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd,iigx,iigy,iigz
     .            ,iorder=order)
#endif
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
     .                ,iigx,iigy,iigz,iorder=order,is_cnv=.true.)

            vecnv = vecnv + jcnv_n

            vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*vecnv(:,:,:,1)/rho
            vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*vecnv(:,:,:,2)/rho
            vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*vecnv(:,:,:,3)/rho
          endif

        case(IDIVPI)

#if !defined(vec_pot)
          if (di > 0d0.and.ion_hall) then
            bcnd = vaux%vec_list(IDIVPI)%bconds !A BCs with extrapolation

            call setBC(IVX,3,nnx,nny,nnz,div_pi,vdum,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order,is_cnv=.true.)
          endif
#endif

        case(IDIVPE)

#if !defined(vec_pot)
          if (di > 0d0.and.(.not.ion_hall)) then
            bcnd = vaux%vec_list(IDIVPE)%bconds !A BCs with extrapolation

            call setBC(IVX,3,nnx,nny,nnz,div_pe,vdum,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order,is_cnv=.true.)
          endif
#endif

        case(IENI_0)

c       Electric field BC (time=0)

          if (di > 0d0) then
#if !defined(vec_pot)
            bcnd = vaux%vec_list(IENI_0)%bconds !A BCs with extrapolation
            v0   = 0d0
            vdum = 0d0
            call setBC(IVX,3,nnx,nny,nnz,vdum,E_h,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order,is_cnv=.false.)
#endif
            E_ni_0 = E_ni_0 + di*E_h
          endif

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
            rho => varray%array_var(IRHO)%array
            heta=hres(nnx,nny,nnz,iigx,iigy,iigz,bcnv,rho)
            nullify(rho)
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
                  etal = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz,psi)
                  jcov_0(i,j,k,:) = E0/etal
                  jcnv_0(i,j,k,:) =
     .                 matmul(gmetric%grid(iigx)%gsup(i,j,k,:,:)
     .                       ,jcov_0(i,j,k,:))
                enddo
              enddo
            enddo

          elseif (sum(E0) == 0d0) then

            bcnd = vaux%vec_list(IJCNV_0)%bconds

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

c          do k = 1,nnz
c            do j = 1,nny
c              do i = 1,nnx
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

          vecnv_0(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv_0(:,:,:,1)/rho
          vecnv_0(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv_0(:,:,:,2)/rho
          vecnv_0(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv_0(:,:,:,3)/rho

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

c       Electric field (time=0)

#if defined(vec_pot)
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,acov)
#else
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,jcov)
#endif
          if (di > 0d0) then
            E_h  = E_fld_Hall(nnx,nny,nnz,iigx,iigy,iigz,rho,tmp)
          else
            E_h  = 0d0
          endif

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
