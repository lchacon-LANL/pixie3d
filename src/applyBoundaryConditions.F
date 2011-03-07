
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

c Begin program

c Setup pointers

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      varray => gv%u_n
      u_0    => gv%u_0

c Set application BC sequence

      call setupVarInitSeq(gv,it)

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

cc      call deallocAuxVariables(gv%aux)
      call deallocLocalAuxVars

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

cc      call deallocAuxVariables(gv%aux)
      call deallocLocalAuxVars

      end subroutine initializeAuxVar

c setupVarInitSeq
c #####################################################################
      subroutine setupVarInitSeq(patch_var,it)
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

c Define boundary conditions and dependencies

      call defineDepBCs(gv%u_n)

      call defineDepBC_Deps(gv%u_n)

      call defineDepDom_Deps(gv%u_n)

c Allocate auxiliary variables

      call allocAuxVariables(gv%aux)

c Set BC schedule (see imposeBC_mod.F)

      call setBCsched(gv%u_n,it)

c End program

      end subroutine setupVarInitSeq

#else

c applyBC
c####################################################################
      subroutine applyBC(varray,iigx,iigy,iigz,obc)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer :: iigx,iigy,iigz
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

      INTERFACE
        subroutine postproc_AuxVar(bc_inf,varray,vaux,iigx,iigy,iigz)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine postproc_AuxVar
      END INTERFACE

c Begin program

      if (PRESENT(obc)) then
        order = obc
      else
        order = bc_order     !Order of extrapolation in BCs (defined in imposeBC_mod.F)
      endif

c Allocate auxiliary variables

      call allocAuxVariables(gv%aux)

c Set application BCs (see imposeBC_mod.F)

      call setBCsched(varray,itime)

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
     .                        ,iigx,iigy,iigz,order=order)
          endif
        enddo

        !Postprocess auxiliary quantities
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
            call postproc_AuxVar(gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                          ,gv%aux,iigx,iigy,iigz)
          endif
        enddo

      enddo

c End program

      end subroutine applyBC

c applyBCplot
c####################################################################
      subroutine applyBCplot(varray,iigx,iigy,iigz)
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

      type(var_array),pointer :: varray

c Local variables

      integer :: i,order,igrp

c Interfaces

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

      INTERFACE
        subroutine postproc_AuxVar(bc_inf,varray,vaux,iigx,iigy,iigz)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine postproc_AuxVar
      END INTERFACE

c Begin program

      order = bc_order          !Order of extrapolation in BCs (defined in imposeBC_mod.F)

c Allocate auxiliary variables

      call allocAuxVariables(gv%aux)

c Set application BCs (see imposeBC_mod.F)

      call setBCsched(varray,itime)

c Impose BCs

      do igrp = 1,gv%nbc_grp
        !Fill auxiliary quantities
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
            call preproc_AuxVar(gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                         ,gv%aux,iigx,iigy,iigz,order=order)
          endif
        enddo

        !Postprocess auxiliary quantities
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
            call postproc_AuxVar(gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                          ,gv%aux,iigx,iigy,iigz)
          endif
        enddo

      enddo

c End program

      end subroutine applyBCplot
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
      real(8),pointer,dimension(:,:,:) :: rho

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

c Select variable

      select case(abs(bc_inf(1)))
      case(IRHO)

c     Density BCs

        call setBC(IRHO,nnx,nny,nnz
     .            ,varray%array_var(IRHO)%array
     .            ,u_0   %array_var(IRHO)%array
     .            ,varray%array_var(IRHO)%bconds
     .            ,iigx,iigy,iigz,iorder=order)

#if defined(vec_pot)
      case(IAX)

        call XferBCs(acov(:,:,:,1),varray%array_var(IAX)%array,.false.)

      case(IAY)

        call XferBCs(acov(:,:,:,2),varray%array_var(IAY)%array,.false.)

      case(IAZ)

        call XferBCs(acov(:,:,:,3),varray%array_var(IAZ)%array,.false.)
#else
      case(IBX)

        call XferBCs(bcnv(:,:,:,1),varray%array_var(IBX)%array,.false.)

      case(IBY)

        call XferBCs(bcnv(:,:,:,2),varray%array_var(IBY)%array,.false.)

      case(IBZ)

        call XferBCs(bcnv(:,:,:,3),varray%array_var(IBZ)%array,.false.)
#endif

      case(IVX)

        if (.not.nc_eom_v) then

          rho => varray%array_var(IRHO)%array

          varray%array_var(IVX)%array(0    ,:,:) = vcnv(0    ,:,:,1)
     .                                            *rho (0    ,:,:)
          varray%array_var(IVX)%array(nnx+1,:,:) = vcnv(nnx+1,:,:,1)
     .                                            *rho (nnx+1,:,:)
          varray%array_var(IVX)%array(:,0    ,:) = vcnv(:,0    ,:,1)
     .                                            *rho (:,0    ,:)
          varray%array_var(IVX)%array(:,nny+1,:) = vcnv(:,nny+1,:,1)
     .                                            *rho (:,nny+1,:)
          varray%array_var(IVX)%array(:,:,0    ) = vcnv(:,:,0    ,1)
     .                                            *rho (:,:,0    )
          varray%array_var(IVX)%array(:,:,nnz+1) = vcnv(:,:,nnz+1,1)
     .                                            *rho (:,:,nnz+1)

          nullify(rho)

        else

          call XferBCs(vcnv(:,:,:,1),varray%array_var(IVX)%array
     .                ,.false.)

        endif

      case(IVY)

        if (.not.nc_eom_v) then

          rho => varray%array_var(IRHO)%array

          varray%array_var(IVY)%array(0    ,:,:) = vcnv(0    ,:,:,2)
     .                                            *rho (0    ,:,:)
          varray%array_var(IVY)%array(nnx+1,:,:) = vcnv(nnx+1,:,:,2)
     .                                            *rho (nnx+1,:,:)
          varray%array_var(IVY)%array(:,0    ,:) = vcnv(:,0    ,:,2)
     .                                            *rho (:,0    ,:)
          varray%array_var(IVY)%array(:,nny+1,:) = vcnv(:,nny+1,:,2)
     .                                            *rho (:,nny+1,:)
          varray%array_var(IVY)%array(:,:,0    ) = vcnv(:,:,0    ,2)
     .                                            *rho (:,:,0    )
          varray%array_var(IVY)%array(:,:,nnz+1) = vcnv(:,:,nnz+1,2)
     .                                            *rho (:,:,nnz+1)

          nullify(rho)

        else

          call XferBCs(vcnv(:,:,:,2),varray%array_var(IVY)%array
     .                ,.false.)

        endif

       case(IVZ)

        if (.not.nc_eom_v) then

          rho => varray%array_var(IRHO)%array

          varray%array_var(IVZ)%array(0    ,:,:) = vcnv(0    ,:,:,3)
     .                                            *rho (0    ,:,:)
          varray%array_var(IVZ)%array(nnx+1,:,:) = vcnv(nnx+1,:,:,3)
     .                                            *rho (nnx+1,:,:)
          varray%array_var(IVZ)%array(:,0    ,:) = vcnv(:,0    ,:,3)
     .                                            *rho (:,0    ,:)
          varray%array_var(IVZ)%array(:,nny+1,:) = vcnv(:,nny+1,:,3)
     .                                            *rho (:,nny+1,:)
          varray%array_var(IVZ)%array(:,:,0    ) = vcnv(:,:,0    ,3)
     .                                            *rho (:,:,0    )
          varray%array_var(IVZ)%array(:,:,nnz+1) = vcnv(:,:,nnz+1,3)
     .                                            *rho (:,:,nnz+1)

          nullify(rho)

        else

          call XferBCs(vcnv(:,:,:,3),varray%array_var(IVZ)%array
     .                ,.false.)

        endif

      case(ITMP)

c     Temperature BCs

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

      integer :: iigx,iigy,iigz,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

      integer,optional :: order

c Local variables

      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz
     .          ,ig,jg,kg

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),grad_p(3),etal

      real(8),allocatable,dimension(:,:,:,:) :: v0,vdum

      real(8),pointer,dimension(:,:,:) :: rho,tmp,prs

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

      rho => varray%array_var(IRHO)%array
      tmp => varray%array_var(ITMP)%array

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

        select case(abs(bc_inf(1)))
        case(IETA)

#if defined(RFX)
          eeta= res(nnx,nny,nnz,iigx,psi)
#else
          eeta= res(nnx,nny,nnz,iigx,tmp)
#endif

        case(INU)

#if defined(RFX)
          nuu = vis(nnx,nny,nnz,iigx,tmp)
#else
          if (slava) then
            nuu = rho*vis(nnx,nny,nnz,iigx,tmp)
          else
            nuu = vis(nnx,nny,nnz,iigx,tmp)
          endif
#endif

        case(IETA_0)

#if defined(RFX)
          eeta_0=res(nnx,nny,nnz,iigx,psi)
#else
          eeta_0=res(nnx,nny,nnz,iigx,tmp_0)
#endif

        case(IHETA)

          h_eta=hres(nnx,nny,nnz,iigx,iigy,iigz,bcnv_0,rho_0)

        case(IRHO_0)

c       Density BCs

          call setBC(IRHO,nnx,nny,nnz,rho_0
     .              ,u_0%array_var(IRHO)%array
     .              ,vaux%var_list(IRHO_0)%bconds
     .              ,iigx,iigy,iigz,iorder=order)

        case(ITMP_0)

c       Temperature BCs

          call setBC(ITMP,nnx,nny,nnz,tmp_0
     .              ,u_0%array_var(ITMP)%array
     .              ,vaux%var_list(ITMP_0)%bconds
     .              ,iigx,iigy,iigz,iorder=order)

        case default
          call pstop('applyAuxVarBC'
     .              ,'Scalar auxiliary variable does not exist')
        end select

      case(VECTOR)

        allocate(prs(0:nnx+1,0:nny+1,0:nnz+1))

c       Allocate auxiliary variables in local domain

        allocate(v0  (0:nnx+1,0:nny+1,0:nnz+1,3)
     .          ,vdum(0:nnx+1,0:nny+1,0:nnz+1,3))

        v0   = 0d0
        vdum = 0d0

        select case(abs(bc_inf(1)))
#if defined(vec_pot)
        case(IACOV_0)

          vdum = XformVector(iigx,acov_0,.true.)
          call setBC(IAX,3,nnx,nny,nnz,vdum,acov_0,v0
     .              ,vaux%vec_list(IACOV_0)%bconds,iigx,iigy,iigz
     .              ,is_cnv=vaux%vec_list(IACOV_0)%cnv
     .              ,iorder=2)
#endif

        case(IBCNV_0)

          call setBC(IBX,3,nnx,nny,nnz,bcnv_0,bcov_0,v0
     .              ,vaux%vec_list(IBCNV_0)%bconds,iigx,iigy,iigz
     .              ,is_cnv=vaux%vec_list(IBCNV_0)%cnv
     .              ,iorder=2)
cc     .              ,iorder=order)

        case(IBCOV_0)   !Done with bcnv_0

cc          call XformVector_BC(iigx,bcnv_0,bcov_0,.false.)

        case(IVCNV_0)

c       Velocity BC

          if (adiabatic) then
            prs = a_p*tmp_0
          else
            prs = a_p*rho_0*tmp_0
          endif

          if (pinch_flow) then
            call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv_0,prs,eeta_0,v0)
          else
            !Average has already been performed in u_0
            call XferBCs(u_0%array_var(IVX)%array,v0(:,:,:,1),.false.)
            call XferBCs(u_0%array_var(IVY)%array,v0(:,:,:,2),.false.)
            call XferBCs(u_0%array_var(IVZ)%array,v0(:,:,:,3),.false.)
          endif

c         Fill ghost nodes

          vdum = XformVector(iigx,vcnv_0,.false.)
          call setBC(IVX,3,nnx,nny,nnz,vcnv_0,vdum,v0
     .              ,vaux%vec_list(IVCNV_0)%bconds
     .              ,iigx,iigy,iigz,iorder=order)

        case(IJCNV_0)

c       Current BC

          call setBC(IJX,3,nnx,nny,nnz,jcnv_0,jcov_0,v0
     .              ,vaux%vec_list(IJCNV_0)%bconds,iigx,iigy,iigz
     .              ,is_cnv=vaux%vec_list(IJCNV_0)%cnv
     .              ,iorder=order)

        case(IJCOV_0)  !Done when computing jcov_0

cc          call XformVector_BC(iigx,jcnv_0,jcov_0,.false.)

        case(IENI_0)

#if defined(vec_pot)
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,acov_0,eeta_0)
#else
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,jcov_0,eeta_0)
#endif

        case(IVECOV_N)

cc          if (save_tn.and.(.not.ion_hall)) then
cc            call XferBCs(vecov,vecov_n,.false.)
cc          endif

        case(IVCOV_N)

cc          if (save_tn.and.ion_hall) call XferBCs(vcov,vcov_n,.false.)

        case(IVCNV_N)

cc          if (save_tn.and.ion_hall) call XferBCs(vcnv,vcnv_n,.false.)

        case(IVCNV)

c       Velocity BC

          if (pinch_flow) then

cc           if (lag_pinch_bc) then
cc             rho => u_n%array_var(IRHO)%array
cc             tmp => u_n%array_var(ITMP)%array
cc           else
             rho => varray%array_var(IRHO)%array
             tmp => varray%array_var(ITMP)%array
cc           endif

            if (adiabatic) then
              prs = a_p*tmp
            else
              prs = a_p*rho*tmp
            endif

cc            if (lag_pinch_bc) then
cc              call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0)
cc            else
              call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0)
cc            endif
          else
            call XferBCs(vcnv_0,v0,.true.)
          endif

          call setBC(IVX,3,nnx,nny,nnz,vcnv,vcov,v0
     .              ,vaux%vec_list(IVCNV)%bconds,iigx,iigy,iigz
     .              ,iorder=order)

        case(IVCOV)   !Done when computing vcnv

cc          call XformVector_BC(iigx,vcnv,vcov,.false.)

#if defined(vec_pot)

        case(IACOV)

c       Vector Pot. BC

          call setBC(IAX,3,nnx,nny,nnz,acnv,acov,acov_0
     .              ,vaux%vec_list(IACOV)%bconds,iigx,iigy,iigz
     .              ,iorder=order,is_cnv=vaux%vec_list(IACOV)%cnv)

        case(IACNV)  !Done when computing acov

cc          call XformVector_BC(iigx,acov,acnv,.true.)

#endif
        case(IBCNV)

c       Magnetic field BC

          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation

cc           if (lag_pinch_bc) then
cc             rho => u_n%array_var(IRHO)%array
cc             tmp => u_n%array_var(ITMP)%array
cc           else
             rho => varray%array_var(IRHO)%array
             tmp => varray%array_var(ITMP)%array
cc           endif

            if (adiabatic) then
              prs = a_p*tmp
            else
              prs = a_p*rho*tmp
            endif

cc            if (lag_pinch_bc) then
cc              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0
cc     .                       ,first_order=.true.)
cc            else
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
     .                       ,first_order=.true.)
cc            endif
          else
            call XferBCs(jcnv_0,v0,.true.)
cc            call XferBCs(jcov_0,v0,.true.)
cc            v0 = 0d0
          endif

          call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0
     .              ,vaux%vec_list(IBCNV)%bconds,iigx,iigy,iigz
     .              ,iorder=order,is_cnv=vaux%vec_list(IBCNV)%cnv)

        case(IBCOV) !Done when computing bcnv

cc          call XformVector_BC(iigx,bcnv,bcov,.false.)

        case(IJCNV)

c       Current BC

          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation

cc           if (lag_pinch_bc) then
cc             rho => u_n%array_var(IRHO)%array
cc             tmp => u_n%array_var(ITMP)%array
cc           else
             rho => varray%array_var(IRHO)%array
             tmp => varray%array_var(ITMP)%array
cc           endif

            if (adiabatic) then
              prs = a_p*tmp
            else
              prs = a_p*rho*tmp
            endif

cc            if (lag_pinch_bc) then
cc              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0
cc     .                       ,first_order=.true.)
cc            else
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
     .                       ,first_order=.true.)
cc            endif
          else
            call XferBCs(jcnv_0,v0,.true.)
cc            call XferBCs(jcov_0,v0,.true.)
cc            v0 = 0d0
          endif

          call setBC(IJX,3,nnx,nny,nnz,jcnv,jcov,v0
     .              ,vaux%vec_list(IJCNV)%bconds,iigx,iigy,iigz
     .              ,iorder=order)

        case(IJCOV) !Done when computing jcnv

cc          call XformVector_BC(iigx,jcnv,jcov,.false.)

        case(IVECNV)

c       Electron velocity BC: ve=v-di*j/rho with j=curl(B)

          vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv(:,:,:,1)/rho
          vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv(:,:,:,2)/rho
          vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv(:,:,:,3)/rho

        case(IVECOV)

cc          call XformVector_BC(iigx,vecnv,vecov,.false.)

        case(IVEFCNV)

c       Fake electron velocity BC: ve=v-di*j*/rho with j*=-dt*curl(curl(E_fake))

          bcnd = vaux%vec_list(IVEFCNV)%bconds
          vdum = XformVector(iigx,vefcnv,.false.)
          call setBC(IJX,3,nnx,nny,nnz,vefcnv,vdum,v0,bcnd
     .              ,iigx,iigy,iigz,iorder=order
     .              ,is_cnv=vaux%vec_list(IVEFCNV)%cnv)

          !Comment out for old TM Hall performance results
          vefcnv = vefcnv + jcnv_0
cc          vefcnv = vefcnv + jcnv_n

          vefcnv(:,:,:,1) = vcnv(:,:,:,1)-di*vefcnv(:,:,:,1)/rho
          vefcnv(:,:,:,2) = vcnv(:,:,:,2)-di*vefcnv(:,:,:,2)/rho
          vefcnv(:,:,:,3) = vcnv(:,:,:,3)-di*vefcnv(:,:,:,3)/rho

        case(IDIVPI)

#if !defined(vec_pot)
          if (di > 0d0.and.ion_hall) then
            bcnd = vaux%vec_list(IDIVPI)%bconds
            vdum = XformVector(iigx,div_pi,.false.)
            call setBC(IVX,3,nnx,nny,nnz,div_pi,vdum,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order
     .                ,is_cnv=vaux%vec_list(IDIVPI)%cnv)
          endif
#endif

        case(IDIVPE)

#if !defined(vec_pot)
          if (di > 0d0.and.(.not.ion_hall)) then
            bcnd = vaux%vec_list(IDIVPE)%bconds
            vdum = XformVector(iigx,div_pe,.false.)
            call setBC(IVX,3,nnx,nny,nnz,div_pe,vdum,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order
     .                ,is_cnv=vaux%vec_list(IDIVPE)%cnv)
          endif
#endif

        case(IEH)

c       Hall electric field BC

#if !defined(vec_pot)
          if (di > 0d0) then
            bcnd = vaux%vec_list(IEH)%bconds
            vdum = XformVector(iigx,E_h,.true.)
            call setBC(IVX,3,nnx,nny,nnz,vdum,E_h,v0,bcnd
     .                ,iigx,iigy,iigz,iorder=order,is_cnv=.false.)
          endif
#endif

        case(IENI)

c       Total electric field BC

#if defined(vec_pot)
          E_ni=E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,acov,eeta)
#else
          E_ni=E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,jcov,eeta)
#endif

cc          E_ni = E_ni - E_ni_0  !Substract imposed electric field

          if (di > 0d0) E_ni = E_ni + di*E_h !Add Hall electric field

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

cc      use imposeBCinterface

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

      real(8),pointer,dimension(:,:,:) :: rho,tmp

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

      rho => varray%array_var(IRHO)%array
      tmp => varray%array_var(ITMP)%array

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

        select case(abs(bc_inf(1)))
        case(IRHO_0)

          rho_0 = u_0%array_var(IRHO)%array

        case(ITMP_0)

          tmp_0 = u_0%array_var(ITMP)%array

        case(IETA_0)

c       Resistivity at initial time
 
#if defined(RFX)
          eeta_0=res(nnx,nny,nnz,iigx,psi)
#else
          eeta_0=res(nnx,nny,nnz,iigx,tmp_0)
#endif
        case(INU)

c       Viscosity (needed for BCs; may depend on T)
 
#if defined(RFX)
          nuu = vis(nnx,nny,nnz,iigx,tmp)
#else
          if (slava) then
            nuu = rho*vis(nnx,nny,nnz,iigx,tmp)
          else
            nuu = vis(nnx,nny,nnz,iigx,tmp)
          endif
#endif
        case(IETA)

c       Resistivity (needed for BCs; may depend on T)
 
#if defined(RFX)
          eeta= res(nnx,nny,nnz,iigx,psi)
#else
          eeta= res(nnx,nny,nnz,iigx,tmp)
#endif
        case(IHETA)

c       Find hyperresistivity (electron viscosity)

          h_eta=hres(nnx,nny,nnz,iigx,iigy,iigz,bcnv_0,rho_0)

        case default
          call pstop('applyAuxVarBC'
     .              ,'Scalar auxiliary variable does not exist')
        end select

      case(VECTOR)

cc        if (lag_pinch_bc) then
cc          rho => u_n%array_var(IRHO)%array
cc          tmp => u_n%array_var(ITMP)%array
cc        else
cc          rho => varray%array_var(IRHO)%array
cc          tmp => varray%array_var(ITMP)%array
cc        endif

c       Allocate auxiliary variables in local domain

        allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,3))

        v0 = 0d0

        select case(abs(bc_inf(1)))
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

          bcov_0 = XformVector(iigx,bcnv_0,.false.)

        case(IJCOV_0)

          jcov_0 = XformVector(iigx,jcnv_0,.false.)

        case(IJCNV_0)

c       Equilibrium current

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

          else

            jcnv_0 = curl(iigx,bcov_0)

          endif

        case(IENI_0)

c       Electric field (time=0, initial EOM equ => only ideal+resistive contributions)

#if defined(vec_pot)
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,acov_0,eeta_0)
#else
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,jcov_0,eeta_0)
#endif

        case(IVCOV_N)

          if (save_tn .and. ion_hall) vcov_n = vcov

        case(IVCNV_N)

          if (save_tn .and. ion_hall) vcnv_n = vcnv

        case(IVECOV_N)

          if (save_tn .and. (.not.ion_hall)) vecov_n = vecov

cc        case(IBCNV_N)
cc
cc          if (save_tn) bcnv_n = bcnv !Save previous time step
cc
cc        case(IJCNV_N)
cc
cc          if (save_tn .and. ion_hall) jcnv_n = jcnv !Save previous time step

        case(IVCOV)

          vcov = XformVector(iigx,vcnv,.false.)

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

          acnv = XformVector(iigx,acov,.true.)

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
        case(IBCOV)

          bcov = XformVector(iigx,bcnv,.false.)

        case(IJCOV)

          jcov = XformVector(iigx,jcnv,.false.)

        case(IJCNV)

          jcnv = curl(iigx,bcov)

        case(IVECNV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          where (rho /= 0d0)
            vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv(:,:,:,1)/rho
            vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv(:,:,:,2)/rho
            vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv(:,:,:,3)/rho
          end where

        case(IVECOV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          vecov = XformVector(iigx,vecnv,.false.)

        case(IVEFCNV)

c       Fake electron velocity: find -dt*curl(curl(E_fake))

          if (save_tn) then !Set dt=0 at time t=t_n

            vefcnv = 0d0

          else
cc            !Equivalent to fake_ve=.false.
cc            vefcnv = jcnv - jcnv_n

cc            !Use full electric field
cc            v0 =-dt*E_ni  !First-order approximation to find j^(n+1)
cc            call transformVector(iigx,iigy,iigz
cc     .                  ,0,nnx+1,0,nny+1,0,nnz+1
cc     .                  ,v0(:,:,:,1)
cc     .                  ,v0(:,:,:,2)
cc     .                  ,v0(:,:,:,3)
cc     .                  ,'cov','cnv')
cc
cc            !Pick certain pieces of electric field: IDEAL
cc            v0 =-crossProduct(iigx,vcov,bcov,.false.)
cc
cc            !Resistive
cc            v0(:,:,:,1) = v0(:,:,:,1) + eeta*jcnv(:,:,:,1)
cc            v0(:,:,:,2) = v0(:,:,:,2) + eeta*jcnv(:,:,:,2)
cc            v0(:,:,:,3) = v0(:,:,:,3) + eeta*jcnv(:,:,:,3)

            !Ion inertia
cc            v0 = v0 + di*(vcnv-vcnv_n)/dt
            v0 =  di*(vcnv-vcnv_n)/dt
cc            v0 =  di*(vcnv-vcnv_0)/dt

            !Ion pressure tensor
            v0(:,:,:,1) = v0(:,:,:,1) + di*div_pi(:,:,:,1)/rho
            v0(:,:,:,2) = v0(:,:,:,2) + di*div_pi(:,:,:,2)/rho
            v0(:,:,:,3) = v0(:,:,:,3) + di*div_pi(:,:,:,3)/rho

            do k=1,nnz
              do j=1,nny
                do i=1,nnx
cc                  vefcnv(i,j,k,:) = -dt*curlcurl(i,j,k,nnx,nny,nnz
cc     .                                          ,iigx,iigy,iigz,v0
cc     .                                          ,vol=.false.)
                  vefcnv(i,j,k,:) = dt*veclaplacian(i,j,k,nnx,nny,nnz
     .                                             ,iigx,iigy,iigz,v0
     .                                             ,vol=.false.)
                enddo
              enddo
            enddo

          endif

c old TM Hall MHD performance results
cc          do k=1,nnz
cc            do j=1,nny
cc              do i=1,nnx
cccc                vefcnv(i,j,k,:) = -dt*curlcurl(i,j,k,nnx,nny,nnz
cccc     .                                       ,iigx,iigy,iigz,vcnv
cccc     .                                       ,vol=.false.)
cc                vefcnv(i,j,k,:) = di*veclaplacian(i,j,k,nnx,nny,nnz
cc     .                                       ,iigx,iigy,iigz,vcnv
cc     .                                       ,vol=.false.)
cc              enddo
cc            enddo
cc          enddo
c old TM Hall performance results

        case(IDIVPI)

c       Divergence of ion stress tensor

          if (nu > 0d0) then
            div_pi = EOM_divPi(nnx,nny,nnz,iigx,iigy,iigz,vcnv,nuu)
          else
            div_pi = 0d0
          endif

        case(IDIVPE)

c       Divergence of electron stress tensor

          if (di > 0d0) then
            if (fake_ve) then
              div_pe=EOM_divPe(nnx,nny,nnz,iigx,iigy,iigz,vefcnv,h_eta)
            else
              div_pe=EOM_divPe(nnx,nny,nnz,iigx,iigy,iigz,vecnv ,h_eta)
            endif
          else
            div_pe = 0d0
          endif

        case(IEH)

c       Hall electric field

          if (di > 0d0) then
            E_h  = E_fld_Hall(nnx,nny,nnz,iigx,iigy,iigz,rho,tmp)
          else
            E_h  = 0d0
          endif

        case(IENI)

c       Ideal+resistive E field

#if defined(vec_pot)
          E_ni=E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,acov,eeta)
#else
          E_ni=E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,jcov,eeta)
#endif

        case default

          if (my_rank == 0) then
            write (*,*) 'Variable identifier=',abs(bc_inf(1))
          endif

          call pstop('preproc_AuxVar'
     .              ,'Vector auxiliary variable does not exist')

        end select

c       Deallocate variables 

        deallocate(v0)

      end select

c End

      end subroutine preproc_AuxVar

c postproc_AuxVar
c####################################################################
      subroutine postproc_AuxVar(bc_inf,varray,vaux,iigx,iigy,iigz)
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

c Local variables

      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz
     .          ,ig,jg,kg

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),grad_p(3),etal

      real(8),pointer,dimension(:,:,:) :: rho,tmp

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

      rho => varray%array_var(IRHO)%array
      tmp => varray%array_var(ITMP)%array

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

        select case(abs(bc_inf(1)))
        case(IRHO_0)

          u_0%array_var(IRHO)%array = rho_0

        case(ITMP_0)

          u_0%array_var(ITMP)%array = tmp_0

        case(IETA_0)

c       Resistivity at initial time
 
#if defined(RFX)
          eeta_0=res(nnx,nny,nnz,iigx,psi)
#else
          eeta_0=res(nnx,nny,nnz,iigx,tmp_0)
#endif
        case(INU)

c       Viscosity (needed for BCs; may depend on T)
 
#if defined(RFX)
          nuu = vis(nnx,nny,nnz,iigx,tmp)
#else
          if (slava) then
            nuu = rho*vis(nnx,nny,nnz,iigx,tmp)
          else
            nuu = vis(nnx,nny,nnz,iigx,tmp)
          endif
#endif
        case(IETA)

c       Resistivity (needed for BCs; may depend on T)
 
#if defined(RFX)
          eeta= res(nnx,nny,nnz,iigx,psi)
#else
          eeta= res(nnx,nny,nnz,iigx,tmp)
#endif
        case(IHETA)

c       Find hyperresistivity (electron viscosity)

          h_eta=hres(nnx,nny,nnz,iigx,iigy,iigz,bcnv_0,rho_0)

        case default
          call pstop('applyAuxVarBC'
     .              ,'Scalar auxiliary variable does not exist')
        end select

      case(VECTOR)

c       Allocate auxiliary variables in local domain

        select case(abs(bc_inf(1)))
        case(IVCNV_0)

          if (.not.nc_eom_v) then
            where (rho_0 /= 0d0)
              u_0%array_var(IVX)%array = vcnv_0(:,:,:,1)*rho_0
              u_0%array_var(IVY)%array = vcnv_0(:,:,:,2)*rho_0
              u_0%array_var(IVZ)%array = vcnv_0(:,:,:,3)*rho_0
            end where
          else
            u_0%array_var(IVX)%array = vcnv_0(:,:,:,1)
            u_0%array_var(IVY)%array = vcnv_0(:,:,:,2)
            u_0%array_var(IVZ)%array = vcnv_0(:,:,:,3)
          endif

#if defined(vec_pot)
        case(IACOV_0)
#endif

        case(IBCNV_0)

          u_0%array_var(IBX)%array = bcnv_0(:,:,:,1)
          u_0%array_var(IBY)%array = bcnv_0(:,:,:,2)
          u_0%array_var(IBZ)%array = bcnv_0(:,:,:,3)

        case(IBCOV_0)

cc          call XformVector_BC(iigx,bcnv_0,bcov_0,.false.)

        case(IJCOV_0)

cc          call XformVector_BC(iigx,jcnv_0,jcov_0,.false.)

        case(IJCNV_0)

        case(IENI_0)

c       Electric field (time=0, initial EOM equ => only ideal+resistive contributions)

#if defined(vec_pot)
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,acov_0,eeta_0)
#else
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,jcov_0,eeta_0)
#endif

        case(IVCOV_N)

          if (save_tn.and.ion_hall) call XferBCs(vcov,vcov_n,.false.)

        case(IVCNV_N)

          if (save_tn.and.ion_hall) call XferBCs(vcnv,vcnv_n,.false.)

        case(IVECOV_N)

          if (save_tn.and.(.not.ion_hall)) then
            call XferBCs(vecov,vecov_n,.false.)
          endif

        case(IVCOV)

cc          call XformVector_BC(iigx,vcnv,vcov,.false.)

        case(IVCNV)

#if defined(vec_pot)

        case(IACOV)

        case(IACNV)

cc          call XformVector_BC(iigx,acov,acnv,.true.)

#endif

        case(IBCNV)

        case(IBCOV)

cc          call XformVector_BC(iigx,bcnv,bcov,.false.)

        case(IJCOV)

        case(IJCNV)

        case(IVECNV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv(:,:,:,1)/rho
          vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv(:,:,:,2)/rho
          vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv(:,:,:,3)/rho

        case(IVECOV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          call XformVector_BC(iigx,vecnv,vecov,.false.)

        case(IVEFCNV)

c       Fake electron velocity: 

          !Comment out for old TM Hall performance results
          vefcnv = vefcnv + jcnv_0
cc          vefcnv = vefcnv + jcnv_n

          vefcnv(:,:,:,1) = vcnv(:,:,:,1)-di*vefcnv(:,:,:,1)/rho
          vefcnv(:,:,:,2) = vcnv(:,:,:,2)-di*vefcnv(:,:,:,2)/rho
          vefcnv(:,:,:,3) = vcnv(:,:,:,3)-di*vefcnv(:,:,:,3)/rho

        case(IDIVPI)

        case(IDIVPE)

        case(IEH)

        case(IENI)

c       Ideal+resistive E field

#if defined(vec_pot)
          E_ni=E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,acov,eeta)
#else
          E_ni=E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,jcov,eeta)
#endif

cc          E_ni = E_ni - E_ni_0  !Substract imposed electric field

          if (di > 0d0) E_ni = E_ni + di*E_h !Add Hall electric field

        case default

          if (my_rank == 0) then
            write (*,*) 'Variable identifier=',abs(bc_inf(1))
          endif

          call pstop('postproc_AuxVar'
     .              ,'Vector auxiliary variable does not exist')

        end select

      end select

c End

      end subroutine postproc_AuxVar

c scalarDirichletBC
c #####################################################################
      subroutine scalarDirichletBC(array,array0,ieq,dim,loc,order)
c ---------------------------------------------------------------------
c     Imposes dirichlet BC. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c        * order -> order of extrapolation (when used)
c ---------------------------------------------------------------------

      use local_BCS_variables

      implicit none

c     Call variables

      integer    :: ieq,dim,loc,order
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

c     Begin program

      call extrapolateBC(array,array0,dim,loc,order)
cc      call extrapolateBC(array,array0,dim,loc,1)

c     End program

      end subroutine scalarDirichletBC

c vectorDirichletBC
c #####################################################################
      subroutine vectorDirichletBC(ivar,array,array0,ieq,dim,loc,order)
c ---------------------------------------------------------------------
c     Imposes dirichlet BC. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c        * order -> order of extrapolation (when used)
c ---------------------------------------------------------------------
        
      use local_BCS_variables

      implicit none

c     Call variables

      integer    :: ieq,dim,loc,order,ivar
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,*)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1,*)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: i,j,k,ig,jg,kg,ibc,hex,hey,hez
      real(8)    :: x1,x2,x3,dh(3),diver

      real(8)    :: gsuper(3,3),jac0

      real(8),allocatable,dimension(:,:) :: rhs

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      select case (ieq)
      case (IVX,IVY,IVZ)

        call extrapolateBC(array (:,:,:,ivar)
     .                    ,array0(:,:,:,ivar)
cc     .                    ,dim,loc,order)
     .                    ,dim,loc,1)  !Required for mass conservation

      case (IBX,IBY,IBZ,IJX,IJY,IJZ)

        if (ivar /= dim) then

          call extrapolateBC(array (:,:,:,ivar)
     .                      ,array0(:,:,:,ivar)
     .                      ,dim,loc,order)

        else  !Imposes divergence-free constraint on B, J

          call solen_bc(array,dim,loc,.false.)

        endif

      case (IAX,IAY,IAZ)

        if (ivar /= dim) then

          call extrapolateBC(array (:,:,:,ivar)
     .                      ,array0(:,:,:,ivar)
     .                      ,dim,loc,order)

        else  !Imposes divergence-free constraint on B, J

          call solen_bc(array,dim,loc,.true.)

        endif

      case default

        messg = 'Equation'//int2char(ieq)//' does not exist'
        call pstop('vectorDirichletBC',messg)

      end select

c     End program

      end subroutine vectorDirichletBC

c scalarNeumannBC
c #####################################################################
      subroutine scalarNeumannBC(array,ieq,dim,loc)
c ---------------------------------------------------------------------
c     Imposes neumann BC for a scalar. On input:
c        * array -> variable upon which BCs are imposed
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c ---------------------------------------------------------------------

      use local_BCS_variables

      implicit none

c     Call variables

      integer    :: ieq,dim,loc
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin,offset
      integer    :: i,j,k,ig,jg,kg,ip,im,jp,jm,kp,km,ibc,icomp
      real(8)    :: x1,x2,x3,dh(3),gsuper(3,3),jac0,hessian(3,3,3)

c     Begin program

      offset = 0
      if (symm) offset = 1

      ibc = (1+loc)+2*(dim-1)

      select case (ibc)
      case (1)
        array(0     ,:,:) = array(1   +offset,:,:)
      case (2)
        array(nxbc+1,:,:) = array(nxbc-offset,:,:)
      case (3)
        array(:,0     ,:) = array(:,1   +offset,:)
      case (4)
        array(:,nybc+1,:) = array(:,nybc-offset,:)
      case (5)
        array(:,:,0     ) = array(:,:,1   +offset)
      case (6)
        array(:,:,nzbc+1) = array(:,:,nzbc-offset)
      end select

cc      select case(dim)
cc      case(1)
cc        allocate(rhs(0:nybc+1,0:nzbc+1))
cc      case(2)
cc        allocate(rhs(0:nxbc+1,0:nzbc+1))
cc      case(3)
cc        allocate(rhs(0:nxbc+1,0:nybc+1))
cc      end select
cc
cc      rhs = 0d0
cc
cc      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
cc     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax )
cc
cc      select case (ieq)
cc      case (IRHO)
cc
cc        do i=imin,imax
cc          do j=jmin,jmax
cc            do k=kmin,kmax
cc
cc              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
cc
cc              gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)
cc
cc              ip = i+1
cc              im = i-1
cc              dh(1) = 2.*grid_params%dxh(ig)
cc              if(isBdry(i,igxbc,1).and.(.not.isBC(i,igxbc,1,PER)))then
cc                dh(1) = grid_params%dx(ig)
cc                im = i
cc              endif
cc              if(isBdry(i,igxbc,2).and.(.not.isBC(i,igxbc,2,PER)))then
cc                dh(1) = grid_params%dx(ig-1)
cc                ip = i
cc              endif
cc
cc              jp = j+1
cc              jm = j-1
cc              dh(2) = 2.*grid_params%dyh(jg)
cc              if(isBdry(j,igybc,3).and.(.not.isBC(j,igybc,3,PER)))then
cc                dh(2) = grid_params%dy(jg)
cc                jm = j
cc              endif
cc              if(isBdry(j,igybc,4).and.(.not.isBC(j,igybc,4,PER)))then
cc                dh(2) = grid_params%dy(jg-1)
cc                jp = j
cc              endif
cc
cc              kp = k+1
cc              km = k-1
cc              dh(3) = 2.*grid_params%dzh(kg)
cc              if(isBdry(k,igzbc,5).and.(.not.isBC(k,igzbc,5,PER)))then
cc                dh(3) = grid_params%dz(kg)
cc                km = k
cc              endif
cc              if(isBdry(k,igzbc,6).and.(.not.isBC(k,igzbc,6,PER)))then
cc                dh(3) = grid_params%dz(kg-1)
cc                kp = k
cc              endif
cccc              ip = min(i+1,nxbc)
cccc              im = max(i-1,1)
cccc              jp = min(j+1,nybc)
cccc              jm = max(j-1,1)
cccc              kp = min(k+1,nzbc)
cccc              km = max(k-1,1)
cccc
cccc              dh(1) = 2.*grid_params%dxh(ig)
cccc              if (i == nxbc) dh(1) = grid_params%dx(ig-1)
cccc              if (i == 1   ) dh(1) = grid_params%dx(ig)
cccc
cccc              dh(2) = 2.*grid_params%dyh(jg)
cccc              if (j == nybc) dh(2) = grid_params%dy(jg-1)
cccc              if (j == 1   ) dh(2) = grid_params%dy(jg)
cccc
cccc              dh(3) = 2.*grid_params%dzh(kg)
cccc              if (k == nzbc) dh(3) = grid_params%dz(kg-1)
cccc              if (k == 1   ) dh(3) = grid_params%dz(kg)
cc
cc              if (dim == 1) then
cc                rhs(j,k) = -dh(dim)
cc     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
cc     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3))
cc     .              /gsuper(dim,dim)
cc              elseif (dim == 2) then
cc                rhs(i,k) = -dh(dim)*
cc     .              (gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
cc     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3))
cc     .              /gsuper(dim,dim)
cc              elseif (dim == 3) then
cc                rhs(i,j) = -dh(dim)
cc     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
cc     .              +gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2))
cc     .              /gsuper(dim,dim)
cc              endif
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc      case (ITMP)
cc
cc        do i=imin,imax
cc          do j=jmin,jmax
cc            do k=kmin,kmax
cc
cc              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
cc
cc              gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)
cc              jac0   = gmetric%grid(igxbc)%jac (i,j,k)
cc
cc              hessian = gmetric%grid(igxbc)%Gamma(i,j,k,:,:,:)
cc
cc              ip = i+1
cc              im = i-1
cc              dh(1) = 2.*grid_params%dxh(ig)
cc              if(isBdry(i,igxbc,1).and.(.not.isBC(i,igxbc,1,PER)))then
cc                dh(1) = grid_params%dx(ig)
cc                im = i
cc              endif
cc              if(isBdry(i,igxbc,2).and.(.not.isBC(i,igxbc,2,PER)))then
cc                dh(1) = grid_params%dx(ig-1)
cc                ip = i
cc              endif
cc
cc              jp = j+1
cc              jm = j-1
cc              dh(2) = 2.*grid_params%dyh(jg)
cc              if(isBdry(j,igybc,3).and.(.not.isBC(j,igybc,3,PER)))then
cc                dh(2) = grid_params%dy(jg)
cc                jm = j
cc              endif
cc              if(isBdry(j,igybc,4).and.(.not.isBC(j,igybc,4,PER)))then
cc                dh(2) = grid_params%dy(jg-1)
cc                jp = j
cc              endif
cc
cc              kp = k+1
cc              km = k-1
cc              dh(3) = 2.*grid_params%dzh(kg)
cc              if(isBdry(k,igzbc,5).and.(.not.isBC(k,igzbc,5,PER)))then
cc                dh(3) = grid_params%dz(kg)
cc                km = k
cc              endif
cc              if(isBdry(k,igzbc,6).and.(.not.isBC(k,igzbc,6,PER)))then
cc                dh(3) = grid_params%dz(kg-1)
cc                kp = k
cc              endif
cccc              ip = min(i+1,nxbc)
cccc              im = max(i-1,1   )
cccc              jp = min(j+1,nybc)
cccc              jm = max(j-1,1   )
cccc              kp = min(k+1,nzbc)
cccc              km = max(k-1,1   )
cccc
cccc              dh(1) = 2.*grid_params%dxh(ig)
cccc              if (i == nxbc) dh(1) = grid_params%dx(ig-1)
cccc              if (i == 1   ) dh(1) = grid_params%dx(ig)
cccc
cccc              dh(2) = 2.*grid_params%dyh(jg)
cccc              if (j == nybc) dh(2) = grid_params%dy(jg-1)
cccc              if (j == 1   ) dh(2) = grid_params%dy(jg)
cccc
cccc              dh(3) = 2.*grid_params%dzh(kg)
cccc              if (k == nzbc) dh(3) = grid_params%dz(kg-1)
cccc              if (k == 1   ) dh(3) = grid_params%dz(kg)
cc
cc              if (dim == 1) then
cc                if (gamma > 1d0) then
cc                  rhs(j,k) =  hessian(1,1,1)*vcnv(i,j,k,1)*vcnv(i,j,k,1)
cc     .                       +hessian(1,2,2)*vcnv(i,j,k,2)*vcnv(i,j,k,2)
cc     .                       +hessian(1,3,3)*vcnv(i,j,k,3)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(1,1,2)*vcnv(i,j,k,1)*vcnv(i,j,k,2)
cc     .                    +2.*hessian(1,1,3)*vcnv(i,j,k,1)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(1,2,3)*vcnv(i,j,k,2)*vcnv(i,j,k,3)
cc                else
cc                  rhs(j,k) = 0d0
cc                endif
cc                rhs(j,k) = -dh(dim)
cc     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
cc     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
cc     .             -0.5/jac0*rhs(j,k))/gsuper(dim,dim)
cc              elseif (dim == 2) then
cc                if (gamma > 1d0) then
cc                  rhs(i,k) =  hessian(2,1,1)*vcnv(i,j,k,1)*vcnv(i,j,k,1)
cc     .                       +hessian(2,2,2)*vcnv(i,j,k,2)*vcnv(i,j,k,2)
cc     .                       +hessian(2,3,3)*vcnv(i,j,k,3)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(2,1,2)*vcnv(i,j,k,1)*vcnv(i,j,k,2)
cc     .                    +2.*hessian(2,1,3)*vcnv(i,j,k,1)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(2,2,3)*vcnv(i,j,k,2)*vcnv(i,j,k,3)
cc                else
cc                  rhs(i,k) = 0d0
cc                endif
cc                rhs(i,k) = -dh(dim)
cc     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
cc     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3)
cc     .             -0.5/jac0*rhs(i,k))/gsuper(dim,dim)
cc              elseif (dim == 3) then
cc                if (gamma > 1d0) then
cc                  rhs(i,j) =  hessian(3,1,1)*vcnv(i,j,k,1)*vcnv(i,j,k,1)
cc     .                       +hessian(3,2,2)*vcnv(i,j,k,2)*vcnv(i,j,k,2)
cc     .                       +hessian(3,3,3)*vcnv(i,j,k,3)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(3,1,2)*vcnv(i,j,k,1)*vcnv(i,j,k,2)
cc     .                    +2.*hessian(3,1,3)*vcnv(i,j,k,1)*vcnv(i,j,k,3)
cc     .                    +2.*hessian(3,2,3)*vcnv(i,j,k,2)*vcnv(i,j,k,3)
cc                else
cc                  rhs(i,j) = 0d0
cc                endif
cc                rhs(i,j) = -dh(dim)
cc     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
cc     .              +gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
cc     .             -0.5/jac0*rhs(i,j))/gsuper(dim,dim)
cc              endif
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc      case default
cc
cc        call pstop('scalarNeumannBC','Variable not implemented')
cc
cc      end select
cc
ccc     Assign value
cc
cc      select case (ibc)
cc      case (1)
cc        array(0     ,jmin:jmax,kmin:kmax)
cc     $                    = array(1+offset,jmin:jmax,kmin:kmax)
cc     $                    - rhs(jmin:jmax,kmin:kmax)
cc      case (2)
cc        array(nxbc+1,jmin:jmax,kmin:kmax)
cc     $                    = array(nxbc-offset,jmin:jmax,kmin:kmax)
cc     $                    + rhs(jmin:jmax,kmin:kmax)
cc      case (3)
cc        array(imin:imax,0     ,kmin:kmax)
cc     $                    = array(imin:imax,1+offset,kmin:kmax)
cc     $                    - rhs(imin:imax,kmin:kmax)
cc      case (4)
cc        array(imin:imax,nybc+1,kmin:kmax)
cc     $                    = array(imin:imax,nybc-offset,kmin:kmax)
cc     $                    + rhs(imin:imax,kmin:kmax)
cc      case (5)
cc        array(imin:imax,jmin:jmax,0     )
cc     $                    = array(imin:imax,jmin:jmax,1+offset)
cc     $                    - rhs(imin:imax,jmin:jmax)
cc      case (6)
cc        array(imin:imax,jmin:jmax,nzbc+1)
cc     $                    = array(imin:imax,jmin:jmax,nzbc-offset)
cc     $                    + rhs(imin:imax,jmin:jmax)
cc      end select
cc
cc      deallocate(rhs)

c     End program

      end subroutine scalarNeumannBC

c vectorNeumannBC
c #####################################################################
      subroutine vectorNeumannBC(ivar,array,array0,ieq,dim,loc)
c     -----------------------------------------------------------------
c     Imposes neumann BC for a scalar. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     -----------------------------------------------------------------

      use local_BCS_variables, vec_mod=>vec

      implicit none

c     Call variables

      integer :: ieq,dim,loc,ivar
      real(8) :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,3)
     .          ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1,3)

c     Local variables

      integer :: imax,imin,jmax,jmin,kmax,kmin,offset
      integer :: i,j,k,ig,jg,kg,ip,im,jp,jm,kp,km,ibc,icomp,isgn
     .          ,igp,igm,jgp,jgm,kgp,kgm,i1,j1,k1,hex,hey,hez
      real(8) :: x1,x2,x3,dh(3),jac0,ee,alp,vec(3),cnv(3)
     .          ,cov(3),J0(3),bnorm,jn,etal,B00(3),diver
      real(8) :: gsuper(3,3),hessian(3,3,3)
     .          ,gsupip(3,3),gsupim(3,3)
     .          ,gsupjp(3,3),gsupjm(3,3)
     .          ,gsupkp(3,3),gsupkm(3,3),gsupavg(3,3)
     .          ,cov00(3),covip(3),covim(3)
     .          ,covjp(3),covjm(3)
     .          ,covkp(3),covkm(3)

      logical :: nonlinear

      real(8),allocatable,dimension(:,:) :: rhs

c     Begin program

      offset = 0
      if (symm) offset = 1

      ibc = (1+loc)+2*(dim-1)

      select case(dim)
      case(1)
        allocate(rhs(0:nybc+1,0:nzbc+1))
      case(2)
        allocate(rhs(0:nxbc+1,0:nzbc+1))
      case(3)
        allocate(rhs(0:nxbc+1,0:nybc+1))
      end select

      rhs = 0d0  !Default

      call findBCLoopLimits(dim,loc,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      select case (ieq)
      case (IVX,IVY,IVZ) !Velocity components

        if (coords/='car'.and.(.not.bcSYM(ibc))) then

          call findBCLoopLimits(dim,loc
     .                         ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                         ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

                gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)

                hessian = gmetric%grid(igxbc)%Gamma(i,j,k,:,:,:)

                ip = i+1
                im = i-1
                dh(1) = 2.*grid_params%dxh(ig)
                if(isBdry(i,igxbc,1).and.(.not.isBC(i,igxbc,1,PER)))then
                  dh(1) = grid_params%dx(ig)
                  im = i
                endif
                if(isBdry(i,igxbc,2).and.(.not.isBC(i,igxbc,2,PER)))then
                  dh(1) = grid_params%dx(ig-1)
                  ip = i
                endif

                jp = j+1
                jm = j-1
                dh(2) = 2.*grid_params%dyh(jg)
                if(isBdry(j,igybc,3).and.(.not.isBC(j,igybc,3,PER)))then
                  dh(2) = grid_params%dy(jg)
                  jm = j
                endif
                if(isBdry(j,igybc,4).and.(.not.isBC(j,igybc,4,PER)))then
                  dh(2) = grid_params%dy(jg-1)
                  jp = j
                endif

                kp = k+1
                km = k-1
                dh(3) = 2.*grid_params%dzh(kg)
                if(isBdry(k,igzbc,5).and.(.not.isBC(k,igzbc,5,PER)))then
                  dh(3) = grid_params%dz(kg)
                  km = k
                endif
                if(isBdry(k,igzbc,6).and.(.not.isBC(k,igzbc,6,PER)))then
                  dh(3) = grid_params%dz(kg-1)
                  kp = k
                endif

                if (dim == 1) then

                  if (ivar /= dim) then

                    if (     coords == 'cyl'
     .                  .or. coords == 'hel'
     .                  .or. coords == 'tor') then

                      rhs(j,k)=dh(dim)/gsuper(dim,dim)*array(i,j,k,ivar)

                    else

                      rhs(j,k) =
     .                       gsuper(dim,1)
     .                        *(hessian(1,ivar,1)*array(i,j,k,1)
     .                         +hessian(2,ivar,1)*array(i,j,k,2)
     .                         +hessian(3,ivar,1)*array(i,j,k,3))
     .                      +gsuper(dim,2)
     .                        *(hessian(1,ivar,2)*array(i,j,k,1)
     .                         +hessian(2,ivar,2)*array(i,j,k,2)
     .                         +hessian(3,ivar,2)*array(i,j,k,3))
     .                      +gsuper(dim,3)
     .                        *(hessian(1,ivar,3)*array(i,j,k,1)
     .                         +hessian(2,ivar,3)*array(i,j,k,2)
     .                         +hessian(3,ivar,3)*array(i,j,k,3))

                      rhs(j,k) = dh(dim)/gsuper(dim,dim)
     .                       *(-gsuper(dim,2)*(array(i,jp,k,ivar)
     .                                        -array(i,jm,k,ivar))/dh(2)
     .                         -gsuper(dim,3)*(array(i,j,kp,ivar)
     .                                        -array(i,j,km,ivar))/dh(3)
     .                         +rhs(j,k) )

                    endif

                  endif

                elseif (dim == 2) then

                  if (ivar /= dim) then
                    rhs(i,k) =
     .                     gsuper(dim,1)
     .                      *(hessian(1,ivar,1)*array(i,j,k,1)
     .                       +hessian(2,ivar,1)*array(i,j,k,2)
     .                       +hessian(3,ivar,1)*array(i,j,k,3))
     .                    +gsuper(dim,2)
     .                      *(hessian(1,ivar,2)*array(i,j,k,1)
     .                       +hessian(2,ivar,2)*array(i,j,k,2)
     .                       +hessian(3,ivar,2)*array(i,j,k,3))
     .                    +gsuper(dim,3)
     .                      *(hessian(1,ivar,3)*array(i,j,k,1)
     .                       +hessian(2,ivar,3)*array(i,j,k,2)
     .                       +hessian(3,ivar,3)*array(i,j,k,3))

                    rhs(i,k) = dh(dim)/gsuper(dim,dim)
     .                      *(-gsuper(dim,1)*(array(ip,j,k,ivar)
     .                                       -array(im,j,k,ivar))/dh(1)
     .                        -gsuper(dim,3)*(array(i,j,kp,ivar)
     .                                       -array(i,j,km,ivar))/dh(3)
     .                        +rhs(i,k))
                  endif

                elseif (dim == 3) then

                  if (ivar /= dim) then
                    rhs(i,j) =
     .                     gsuper(dim,1)
     .                      *(hessian(1,ivar,1)*array(i,j,k,1)
     .                       +hessian(2,ivar,1)*array(i,j,k,2)
     .                       +hessian(3,ivar,1)*array(i,j,k,3))
     .                    +gsuper(dim,2)
     .                      *(hessian(1,ivar,2)*array(i,j,k,1)
     .                       +hessian(2,ivar,2)*array(i,j,k,2)
     .                       +hessian(3,ivar,2)*array(i,j,k,3))
     .                    +gsuper(dim,3)
     .                      *(hessian(1,ivar,3)*array(i,j,k,1)
     .                       +hessian(2,ivar,3)*array(i,j,k,2)
     .                       +hessian(3,ivar,3)*array(i,j,k,3))

                    rhs(i,j) = dh(dim)/gsuper(dim,dim)
     .                      *(-gsuper(dim,1)*(array(ip,j,k,ivar)
     .                                       -array(im,j,k,ivar))/dh(1)
     .                        -gsuper(dim,2)*(array(i,jp,k,ivar)
     .                                       -array(i,jm,k,ivar))/dh(2)
     .                        +rhs(i,j))
                  endif

                endif

              enddo
            enddo
          enddo

        endif

        select case (ibc)
        case (1)
          array(0     ,:,:,ivar) = array(1+offset   ,:,:,ivar) - rhs
        case (2)
          array(nxbc+1,:,:,ivar) = array(nxbc-offset,:,:,ivar) + rhs
        case (3)
          array(:,0     ,:,ivar) = array(:,1+offset   ,:,ivar) - rhs
        case (4)
          array(:,nybc+1,:,ivar) = array(:,nybc-offset,:,ivar) + rhs
        case (5)
          array(:,:,0     ,ivar) = array(:,:,1+offset   ,ivar) - rhs
        case (6)
          array(:,:,nzbc+1,ivar) = array(:,:,nzbc-offset,ivar) + rhs
        end select

      case (IJX,IJY,IJZ) !J components (for symmetry BCs)

        if (ivar /= dim) then

          select case (ibc)
          case (1)
            array(0     ,:,:,ivar) = array(1+offset   ,:,:,ivar) - rhs
          case (2)                                                    
            array(nxbc+1,:,:,ivar) = array(nxbc-offset,:,:,ivar) + rhs
          case (3)                                                    
            array(:,0     ,:,ivar) = array(:,1+offset   ,:,ivar) - rhs
          case (4)                                                    
            array(:,nybc+1,:,ivar) = array(:,nybc-offset,:,ivar) + rhs
          case (5)                                                    
            array(:,:,0     ,ivar) = array(:,:,1+offset   ,ivar) - rhs
          case (6)                                                    
            array(:,:,nzbc+1,ivar) = array(:,:,nzbc-offset,ivar) + rhs
          end select

        else  !Imposes divergence-free constraint on J

          call solen_BC(array,dim,loc,.false.)

        endif

      case (IAX,IAY,IAZ) !A components (for symmetry BCs)

        if (ivar /= dim) then

          select case (ibc)
          case (1)
            array(0     ,:,:,ivar) = array(1+offset   ,:,:,ivar) - rhs
          case (2)                                                    
            array(nxbc+1,:,:,ivar) = array(nxbc-offset,:,:,ivar) + rhs
          case (3)                                                    
            array(:,0     ,:,ivar) = array(:,1+offset   ,:,ivar) - rhs
          case (4)                                                    
            array(:,nybc+1,:,ivar) = array(:,nybc-offset,:,ivar) + rhs
          case (5)                                                    
            array(:,:,0     ,ivar) = array(:,:,1+offset   ,ivar) - rhs
          case (6)                                                    
            array(:,:,nzbc+1,ivar) = array(:,:,nzbc-offset,ivar) + rhs
          end select

        else  !Imposes divergence-free constraint on A

          call solen_BC(array,dim,loc,.true.)

        endif

      case (IBX,IBY,IBZ) 

        if (ivar /= dim) then

          call findBCLoopLimits(dim,loc
     .                         ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                         ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

                gsuper = gmetric%grid(igxbc)%gsup(i  ,j,k,:,:)
                gsupip = gmetric%grid(igxbc)%gsup(i+1,j,k,:,:)
                gsupim = gmetric%grid(igxbc)%gsup(i-1,j,k,:,:)
                gsupjp = gmetric%grid(igxbc)%gsup(i,j+1,k,:,:)
                gsupjm = gmetric%grid(igxbc)%gsup(i,j-1,k,:,:)
                gsupkp = gmetric%grid(igxbc)%gsup(i,j,k+1,:,:)
                gsupkm = gmetric%grid(igxbc)%gsup(i,j,k-1,:,:)

                ip = i+1
                im = i-1
                dh(1) = 2.*grid_params%dxh(ig)
cc                if(isBdry(i,igxbc,1).and.(.not.isBC(i,igxbc,1,PER)))then  !Gives bdry noise in Soloveev case
                if(isBdry(i,igxbc,1)) then
                  dh(1) = grid_params%dx(ig)
                  im = i
                endif
cc                if(isBdry(i,igxbc,2).and.(.not.isBC(i,igxbc,2,PER)))then  !Gives bdry noise in Soloveev case
                if(isBdry(i,igxbc,2)) then
                  dh(1) = grid_params%dx(ig-1)
                  ip = i
                endif

                jp = j+1
                jm = j-1
                dh(2) = 2.*grid_params%dyh(jg)
cc                if(isBdry(j,igybc,3).and.(.not.isBC(j,igybc,3,PER)))then  !Gives bdry noise in Soloveev case
                if(isBdry(j,igybc,3)) then
                  dh(2) = grid_params%dy(jg)
                  jm = j
                endif
cc                if(isBdry(j,igybc,4).and.(.not.isBC(j,igybc,4,PER)))then  !Gives bdry noise in Soloveev case
                if(isBdry(j,igybc,4)) then
                  dh(2) = grid_params%dy(jg-1)
                  jp = j
                endif

                kp = k+1
                km = k-1
                dh(3) = 2.*grid_params%dzh(kg)
cc                if(isBdry(k,igzbc,5).and.(.not.isBC(k,igzbc,5,PER)))then  !Gives bdry noise in Soloveev case
                if(isBdry(k,igzbc,5)) then
                  dh(3) = grid_params%dz(kg)
                  km = k
                endif
cc                if(isBdry(k,igzbc,6).and.(.not.isBC(k,igzbc,6,PER)))then  !Gives bdry noise in Soloveev case
                if(isBdry(k,igzbc,6)) then
                  dh(3) = grid_params%dz(kg-1)
                  kp = k
                endif

                if (dim == 1) then

                  if (ibc == 1) then
                    gsupavg = 0.5*(gsuper+gsupim)
                  elseif (ibc == 2) then
                    gsupavg = 0.5*(gsuper+gsupip)
                  endif

                  jn = (array(i,jp,k,3)-array(i,jm,k,3))/dh(2)
     .                -(array(i,j,kp,2)-array(i,j,km,2))/dh(3)

                  if (ivar == 2) then
                    icomp = 3
                    isgn  = 1

                    rhs(j,k) = 
     .                  (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                     +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                  elseif (ivar == 3) then
                    icomp = 2
                    isgn  =-1

                    rhs(j,k) =
     .                  (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                     +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                  endif

                  if (ibc == 1) then
                    ee = array0(i-1,j,k,icomp)
                  elseif (ibc == 2) then
                    ee = array0(i+1,j,k,icomp)
                  endif

                  rhs(j,k)=dh(dim)*(rhs(j,k) + isgn*ee)

                elseif (dim == 2) then

                  if (ibc == 3) then
                    gsupavg = 0.5*(gsuper+gsupjm)
                  elseif (ibc == 4) then
                    gsupavg = 0.5*(gsuper+gsupjp)
                  endif

                  jn = (array(i,j,kp,1)-array(i,j,km,1))/dh(3)
     .                -(array(ip,j,k,3)-array(im,j,k,3))/dh(1)

                  if (ivar == 3) then
                    icomp = 1
                    isgn  = 1

                    rhs(i,k) =
     .                   (array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                    +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                  elseif (ivar == 1) then
                    icomp = 3
                    isgn  =-1

                    rhs(i,k) = 
     .                  (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                    +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                  endif

                  if (ibc == 3) then
                    ee = array0(i,j-1,k,icomp)
                  elseif (ibc == 4) then
                    ee = array0(i,j+1,k,icomp)
                  endif
 
                  rhs(i,k) = dh(dim)*(rhs(i,k) + isgn*ee)

                elseif (dim == 3) then

                  if (ibc == 5) then
                    gsupavg = 0.5*(gsuper+gsupkm)
                  elseif (ibc == 6) then
                    gsupavg = 0.5*(gsuper+gsupkp)
                  endif

                  jn = (array(ip,j,k,2)-array(im,j,k,2))/dh(1)
     .                -(array(i,jp,k,1)-array(i,jm,k,1))/dh(2)

                  if (ivar == 1) then
                    icomp = 2
                    isgn  = 1

                    rhs(i,j) =
     .                   (array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                     +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                  elseif (ivar == 2) then
                    icomp = 1
                    isgn  =-1

                    rhs(i,j) = 
     .                   (array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                     +isgn*gsupavg(dim,icomp)*jn/gsupavg(dim,dim)
                  endif

                  if (ibc == 5) then
                    ee = array0(i,j,k-1,icomp)
                  elseif (ibc == 6) then
                    ee = array0(i,j,k+1,icomp)
                  endif

                  rhs(i,j) = dh(dim)*(rhs(i,j) + isgn*ee)

                endif

              enddo
            enddo
          enddo

          select case (ibc)
          case (1)
            array(0     ,jmin:jmax,kmin:kmax,ivar)
     $                     = array(1+offset,jmin:jmax,kmin:kmax,ivar)
     $                     - rhs(jmin:jmax,kmin:kmax)
          case (2)
            array(nxbc+1,jmin:jmax,kmin:kmax,ivar)
     $                     = array(nxbc-offset,jmin:jmax,kmin:kmax,ivar)
     $                     + rhs(jmin:jmax,kmin:kmax)
          case (3)
            array(imin:imax,0     ,kmin:kmax,ivar)
     $                     = array(imin:imax,1+offset,kmin:kmax,ivar)
     $                     - rhs(imin:imax,kmin:kmax)
          case (4)
            array(imin:imax,nybc+1,kmin:kmax,ivar)
     $                     = array(imin:imax,nybc-offset,kmin:kmax,ivar)
     $                     + rhs(imin:imax,kmin:kmax)
          case (5)
            array(imin:imax,jmin:jmax,0     ,ivar)
     $                     = array(imin:imax,jmin:jmax,1+offset,ivar)
     $                     - rhs(imin:imax,jmin:jmax)
          case (6)
            array(imin:imax,jmin:jmax,nzbc+1,ivar)
     $                     = array(imin:imax,jmin:jmax,nzbc-offset,ivar)
     $                     + rhs(imin:imax,jmin:jmax)
          end select

        else  !Imposes divergence-free constraint on B

          call solen_BC(array,dim,loc,.false.)

        endif

      case default

        write (*,*) 'Variable=',ieq
        call pstop('vectorNeumanBC','Variable not implemented')

      end select

c     Assign value

      deallocate(rhs)

c     End program

      end subroutine vectorNeumannBC

c solen_BC
c #####################################################################
      subroutine solen_BC(array,dim,loc,face)
c ---------------------------------------------------------------------
c     Imposes dirichlet BC. On input:
c        * array -> contains vector
c        * dim   -> dimension we are imposing BC on (X,Y,Z)
c        * loc   -> boundary location (0 -> left, 1->right)
c ---------------------------------------------------------------------
        
      use local_BCS_variables

      implicit none

c     Call variables

      integer :: dim,loc
      real(8) :: array(0:nxbc+1,0:nybc+1,0:nzbc+1,3)
      logical :: face

c     Local variables

      integer :: imax,imin,jmax,jmin,kmax,kmin,offset
      integer :: i,j,k,ig,jg,kg,ibc,hex,hey,hez
      real(8) :: x1,x2,x3,dh(3),diver,coeff,gsuper(3,3),jac0

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      if (face) then
        offset= 0
        coeff = 1d0
      else
        offset= 1
        coeff = 2d0
      endif

      call findBCLoopLimits(dim,loc
cc     .                     ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                     ,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                     ,imin ,imax ,jmin ,jmax ,kmin ,kmax )

      select case (ibc)
      case (1)
 
        i = 1

        array(i-1,jmin:jmax,kmin:kmax,dim)
     .       = array(i+1,jmin:jmax,kmin:kmax,dim) !To cancel normal derivative

        do k=kmin,kmax
          do j=jmin,jmax

            call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
            jac0 = gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*grid_params%dxh(ig)
            dh(2) = coeff*grid_params%dyh(jg)
            dh(3) = coeff*grid_params%dzh(kg)

            diver=jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)

            array(i-1,j,k,dim) = array(i+offset,j,k,dim) + dh(dim)*diver

          enddo
        enddo
 
      case (2)

        i = nxbc

        array(i+1,jmin:jmax,kmin:kmax,dim)
     .       = array(i-1,jmin:jmax,kmin:kmax,dim) !To cancel normal derivative

        do k=kmin,kmax
          do j=jmin,jmax

            call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

            jac0 = gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*grid_params%dxh(ig)
            dh(2) = coeff*grid_params%dyh(jg)
            dh(3) = coeff*grid_params%dzh(kg)

            diver=jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)

            array(i+1,j,k,dim) = array(i-offset,j,k,dim) - dh(dim)*diver

          enddo
        enddo

      case (3)

        j = 1

        array(imin:imax,j-1,kmin:kmax,dim)
     .       = array(imin:imax,j+1,kmin:kmax,dim) !To cancel normal derivative

        do k=kmin,kmax
          do i=imin,imax

            call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
            jac0 = gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*grid_params%dxh(ig)
            dh(2) = coeff*grid_params%dyh(jg)
            dh(3) = coeff*grid_params%dzh(kg)

            diver=jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)

            array(i,j-1,k,dim) = array(i,j+offset,k,dim) + dh(dim)*diver
              
          enddo
        enddo

      case (4)

        j = nybc

        array(imin:imax,j+1,kmin:kmax,dim)
     .       = array(imin:imax,j-1,kmin:kmax,dim) !To cancel normal derivative

        do k=kmin,kmax
          do i=imin,imax

            call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
            jac0 = gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*grid_params%dxh(ig)
            dh(2) = coeff*grid_params%dyh(jg)
            dh(3) = coeff*grid_params%dzh(kg)

            diver=jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)

            array(i,j+1,k,dim) = array(i,j-offset,k,dim) - dh(dim)*diver
              
          enddo
        enddo

      case (5)

        k = 1

        array(imin:imax,jmin:jmax,k-1,dim)
     .       = array(imin:imax,jmin:jmax,k+1,dim) !To cancel normal derivative

        do j=jmin,jmax
          do i=imin,imax

            call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
            jac0 = gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*grid_params%dxh(ig)
            dh(2) = coeff*grid_params%dyh(jg)
            dh(3) = coeff*grid_params%dzh(kg)

            diver=jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)

            array(i,j,k-1,dim) = array(i,j,k+offset,dim) + dh(dim)*diver

          enddo
        enddo

      case (6)

        k = nzbc

        array(imin:imax,jmin:jmax,k+1,dim)
     .       = array(imin:imax,jmin:jmax,k-1,dim) !To cancel normal derivative

        do j=jmin,jmax
          do i=imin,imax

            call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
            jac0 = gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*grid_params%dxh(ig)
            dh(2) = coeff*grid_params%dyh(jg)
            dh(3) = coeff*grid_params%dzh(kg)

            diver=jac0*div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)

            array(i,j,k+1,dim) = array(i,j,k-offset,dim) - dh(dim)*diver

          enddo
        enddo

      case default

        messg = 'Boundary'//int2char(ibc)//' non existent'
        call pstop('solen_BC',messg)

      end select

      end subroutine solen_BC

c create_bc_schedule
c######################################################################
      subroutine create_bc_schedule(u,aux,nbc_grp,bc_grp,
     &   inherit,interior_only,print_used,reverse)

      use variable_setup, only: var_array,aux_array,bcgrp

c     Input variables
      type(var_array) :: u                  ! We need nvar, bc_dep_list, dom_dep_list, descr for each variable
      type(aux_array) :: aux                ! We need nvar, nvec, bc_dep_list, dom_dep_list, descr for each variable
      logical :: inherit                    ! Do we want to inherit the dependencies for the groups (set to false)
      logical :: interior_only              ! Do we want to allow some variables to fill interiors only (potential performance increase if true)
      logical :: print_used                 ! Do we want to print how the variables are used
      logical :: reverse                    ! Do we want to reverse the default order within a group (useful for debugging, dependencies will always be preserved)

c     Output variables
      integer :: nbc_grp                    ! Number of boundary condition groups
      type(bcgrp), pointer :: bc_grp(:)     ! Number of boundary condition groups

c     Internal variables
      integer :: N, i, j, k, m, nn, N_groups
      integer :: error, n_var(4)
      logical :: test
      INTEGER, POINTER :: dep_list(:,:)
      integer,pointer :: bc_dep(:,:) ,dom_dep(:,:), tot_dep(:,:)    ! Pointers used to store the dependency tables
      integer,pointer :: n_seq(:), bc_groups(:,:), used(:)          ! Pointers used to store the internal group information
      integer,pointer :: tmp1(:), tmp2(:,:)                         ! Pointers used for temporary memory

      n_var(1) = u%nvar     ! Number of scalar dependent variables
      n_var(2) = 0          ! Number of vector dependent variables
      n_var(3) = aux%nvar   ! Number of scalar auxiliary variables
      n_var(4) = aux%nvec   ! Number of vector auxiliary variables
      N = n_var(1) + n_var(2) + n_var(3) + n_var(4)

c     Search for invalid or duplicate entries (this is bad practice and may indicate an error)
      error = 0;
      do i = 1,n_var(1)
         if ( size(u%array_var(i)%bc_dep_list,1)/=N .OR. 
     &        size(u%array_var(i)%dom_dep_list,1)/=N ) then
            call pstop('create_bc_schedule',
     &         'Unexpected size for dependency list')
         endif
         error = check_dep_list(N,u%array_var(i)%bc_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in dependent variable (bc_dep_list)')
         end if
         error = check_dep_list(N,u%array_var(i)%dom_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in dependent variable (dom_dep_list)')
         end if
      end do
      do i = 1,n_var(2)
         !!!! Place holder for dependent vectors
      end do
      do i = 1,n_var(3)
         if ( size(aux%var_list(i)%bc_dep_list,1)/=N .OR. 
     &        size(aux%var_list(i)%dom_dep_list,1)/=N ) then
            call pstop('create_bc_schedule',
     &         'Unexpected size for dependency list')
         endif
         error = check_dep_list(N,aux%var_list(i)%bc_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in scalar auxiliary '//
     &      'variable (bc_dep_list)')
         endif
         error = check_dep_list(N,aux%var_list(i)%dom_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in scalar auxiliary '//
     &      'variable (dom_dep_list)')
         endif
      enddo
      do i = 1,n_var(4)
         if ( size(aux%vec_list(i)%bc_dep_list,1)/=N .OR. 
     &        size(aux%vec_list(i)%dom_dep_list,1)/=N ) then
            call pstop('create_bc_schedule',
     &         'Unexpected size for dependency list')
         endif
         error = check_dep_list(N,aux%vec_list(i)%bc_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in vector auxiliary '//
     &      'variable (bc_dep_list)')
         endif
         error = check_dep_list(N,aux%vec_list(i)%dom_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in vector auxiliary '//
     &      'variable (dom_dep_list)')
         endif
      enddo

c     Create complete dependency tables
      allocate ( bc_dep(N,N) )
      allocate ( dom_dep(N,N) )
      allocate ( tot_dep(N,N) )
      bc_dep = 0
      dom_dep = 0
      tot_dep = 0
c     Add the dependencies for the scalar dependent variables
      do i = 1,n_var(1)
         j = 1
         dep_list => u%array_var(i)%bc_dep_list
         do while ( dep_list(j,1)/=0 )
            k = get_dep_index( dep_list(j,:), n_var )
            bc_dep(i,k) = dep_list(j,3)+1
            j = j+1
         enddo
         j = 1
         dep_list => u%array_var(i)%dom_dep_list
         do while ( dep_list(j,1)/=0 )
            k = get_dep_index( dep_list(j,:), n_var )
            dom_dep(i,k) = dep_list(j,3)+1
            j = j+1
         enddo
      enddo
c     Add the dependencies for the vector dependent variables
      do i = 1,n_var(2)
         !!!! Place holder for dependent vectors
      enddo
c     Add the dependencies for the scalar auxiliary variables
      do i = 1,n_var(3)
         j = 1
         dep_list => aux%var_list(i)%bc_dep_list
         do while ( dep_list(j,1)/=0 )
            k = get_dep_index( dep_list(j,:), n_var )
            bc_dep(i+n_var(1)+n_var(2),k) = dep_list(j,3)+1
            j = j+1
         enddo
         j = 1
         dep_list => aux%var_list(i)%dom_dep_list
        do  while ( dep_list(j,1)/=0 )
            k = get_dep_index( dep_list(j,:), n_var )
            dom_dep(i+n_var(1)+n_var(2),k) = dep_list(j,3)+1
            j = j+1
         enddo
      enddo
c     Add the dependencies for the vector auxiliary variables
      do i = 1,n_var(4)
         j = 1
         dep_list => aux%vec_list(i)%bc_dep_list
         do while ( dep_list(j,1)/=0 ) 
            k = get_dep_index( dep_list(j,:), n_var )
            bc_dep(i+n_var(1)+n_var(2)+n_var(3),k) = dep_list(j,3)+1
            j = j+1
         enddo
         j = 1
         dep_list => aux%vec_list(i)%dom_dep_list
         do while ( dep_list(j,1)/=0 )
            k = get_dep_index( dep_list(j,:), n_var )
            dom_dep(i+n_var(1)+n_var(2)+n_var(3),k) = dep_list(j,3)+1
            j = j+1
         enddo
      enddo

c     Check which variables are used from the dependency list and fill the total dependency table
c     (0: Variable is not needed, 1: Only interior needs to be filled, 2: Interior and BC are needed)
      allocate ( used(N) )
      used = 0
      tot_dep = 0
c     All dependent variables need to have BC applied 
c     Note: we do not need to respect the interior dependencies for determining the groups and group order)
      nn = n_var(1)+n_var(2)
      used(1:nn) = 2
      tot_dep(1:nn,:) = bc_dep(1:nn,:)
c     Set the dependencies needed for the dependent variables (note: we need to respect the domain and bc dependencies)
      do i = nn+1,N
         used(i) = MAX(MAXVAL(dom_dep(1:nn,i)),MAXVAL(bc_dep(1:nn,i)))
         if ( (.NOT.interior_only) .AND. used(i)==1 ) then
c           We are requiring that we will fill bc on all used variables
            used(i) = 2
         endif
         if ( used(i)==1 ) then
c           We only need to fill the interior domain
            tot_dep(i,:) = dom_dep(i,:)
         else if ( used(i)==2 ) then
c           We only need to fill the interior domain
            do j = 1,N
               tot_dep(i,j) = MAX(dom_dep(i,j),bc_dep(i,j))
            enddo
         endif
      enddo
c     Set the dependencies needed for all the used variables
      test = .TRUE.
      do while ( test )
         test = .FALSE.
         do i = 1,N
            if ( used(i)>0 ) then
               do j = 1,N
                  if ( tot_dep(i,j) > used(j) ) then
                     used(j) = tot_dep(i,j)
                     if ( (.NOT.interior_only) .AND. used(j)==1 ) then
c                       We are requiring that we will fill bc on all used variables
                        used(j) = 2
                     endif
                     if ( used(j)==1 ) then
c                       We only need to fill the interior domain
                        tot_dep(j,:) = dom_dep(j,:)
                     else if ( used(j)==2 ) then
c                       We only need to fill the interior domain
                        do k = 1,N
                           tot_dep(j,k) = MAX(dom_dep(j,k),bc_dep(j,k))
                        enddo
                     endif
                     test = .TRUE.
                  endif
               enddo
            endif
         enddo
      enddo
c     Remove self-dependencies (this can create issues when we do the sorting)
      do i = 1,N
         tot_dep(i,i) = 0
      enddo

c     List the variables that are used
      if ( print_used ) then
         allocate ( tmp1(4) )
c        Determine the number of variables with bc, interior, unused, and dep
         tmp1 = 0
         tmp1(1) = n_var(1)+n_var(2)
         do i = n_var(1)+n_var(2)+1,N
            if ( used(i)==2 ) then
               tmp1(2) = tmp1(2)+1
            else if ( used(i)==1 ) then
               tmp1(3) = tmp1(3)+1
            else if ( used(i)==0 ) then
               tmp1(4) = tmp1(4)+1
            endif
         enddo
         write (*,100) N, tmp1(1), tmp1(2), tmp1(3), tmp1(4)
100      format ( 1x, 1I2, ' total variables: ',
     &                1I2, ' dependent variables, ',
     &                1I2, ' bc applied, ',
     &                1I2, ' interiors only, ',
     &                1I2, ' unused' )
c        Write the auxiliary variables and how they are used
         write (*,*) 'Dependent variables: '
         do i = 1,SUM(n_var(1:2))
             call print_var_name(u,aux,i)
         enddo
         write (*,*) 'Variables that need bc applied: '
         do i = SUM(n_var(1:2)),N
            if ( used(i)==2 ) then
               call print_var_name(u,aux,i)
            endif
         enddo
         write (*,*) 'Variables that only need interiors filled: '
         do i = SUM(n_var(1:2)),N
c LC 1/25/11, small bug:            if ( used(i)==2 ) then
            if ( used(i)==1 ) then
               call print_var_name(u,aux,i)
            endif
         enddo
         write (*,*) 'Variables that are unused: '
         do i = SUM(n_var(1:2)),N
            if ( used(i)==0 ) then
               call print_var_name(u,aux,i)
            endif
         enddo
         deallocate ( tmp1 )
      endif

c     Check if any variables have no dependencies and are used (this is an error)
c     Note that if a variable only depends on its self, it is OK (it might depend on the equilibrium data)
      do i = n_var(1)+n_var(2)+1,N      ! Only auxiliary variables need dependencies
         if ( MAXVAL(tot_dep(i,:))==0 .AND. used(i)>0 .AND. 
     &        bc_dep(i,i)==0 .AND. tot_dep(i,i)==0 ) then
            call pstop('create_bc_schedule',
     &        'Variables are used but do not contain any dependencies')
         endif
      enddo

c     Create a first group containing all used variables
      allocate ( n_seq(N+1) )
      allocate ( bc_groups(N+1,N) )
      n_seq = 0
      bc_groups = 0
      N_groups = 1
      do i = 1,N
         if ( used(i)>0 ) then
            bc_groups(1,n_seq(1)+1) = i
            n_seq(1) = n_seq(1)+1
         endif
      enddo

c     While there are any dependencies, create new groups consisting of those dependencies
      i = 1
      do
c        Identify variables that can be in the same group and variables that must be in an earlier group
         do j = n_seq(i),1,-1
            test = .FALSE.
            do k = 1,n_seq(i)
               if ( tot_dep(bc_groups(i,k),bc_groups(i,j))==2 ) then
                  test = .TRUE.
               endif
            end do
            if ( test ) then
               n_seq(i+1) = n_seq(i+1)+1
               bc_groups(i+1,n_seq(i+1)) = bc_groups(i,j)
               N_groups = i+1
               do k = j+1,n_seq(i)
                  bc_groups(i,k-1) = bc_groups(i,k)
               enddo
               bc_groups(i,n_seq(i)) = 0
               n_seq(i) = n_seq(i)-1
            endif
         enddo
c        Move any variables that are needed by the variables in the earlier group to the earlier group
         test = .TRUE.
         do while ( test ) 
            test = .FALSE.
            j = 1
            do while (j <= n_seq(i+1))
               k = 1
               do while (k <= n_seq(i))
                  if ( tot_dep(bc_groups(i+1,j),bc_groups(i,k))>0 ) then
                     n_seq(i+1) = n_seq(i+1)+1
                     bc_groups(i+1,n_seq(i+1)) = bc_groups(i,k)                     
                     do m = k+1,n_seq(i)
                        bc_groups(i,m-1) = bc_groups(i,m)
                     enddo
                     bc_groups(i,n_seq(i)) = 0
                     n_seq(i) = n_seq(i)-1
                     test = .TRUE.
                  endif
                  k = k + 1
               enddo
               j = j + 1
            enddo 
         enddo
         if ( n_seq(i)==0 ) then
c           All variables are in the earlier group, there must be a circular dependence
            write (*,*) 'Problem separating groups, ', 
     &         'circular dependence detected between: '
            do j = 1,n_seq(i+1)
               call print_var_name(u,aux,bc_groups(i+1,j))
            enddo
            call pstop('create_bc_schedule',
     &        'error creating schedule, circular dependence detected')
         endif 
         if ( n_seq(i+1)==0 ) then
c           We are done creating new groups
            exit
         endif
         if ( i==N ) then
c           We have had too many iterations (this should never occur)
            call pstop('create_bc_schedule','Infinite loop')
         endif
         i = i+1
      enddo

c     Change the group order so it reflects the order in which the variables need to be processed
      allocate ( tmp1(N+1) )
      allocate ( tmp2(N+1,N) )
      tmp1 = n_seq
      tmp2 = bc_groups
      do i = 1,N_groups
         bc_groups(i,:) = tmp2(N_groups-i+1,:)
         n_seq(i) = tmp1(N_groups-i+1)
      end do
      deallocate ( tmp1 )
      deallocate ( tmp2 )

c     Sort the variables in each group in ascending order (to create a consistent default order)
      allocate ( tmp1(N) )
      do i = 1,N_groups
         tmp1 = 0
         do j = 1,n_seq(i)
            tmp1(bc_groups(i,j)) = 1
         enddo
         j = 1
         do k = 1,N
            if ( tmp1(k)==1 ) then
               bc_groups(i,j) = k
               j = j+1
            endif
         enddo
      end do
      deallocate ( tmp1 )

c     Reverse the default order within a group (useful for debugging, dependencies will always be preserved)
      if ( reverse ) then
         do i = 1,N_groups
            allocate ( tmp1(n_seq(i)) )
            tmp1 = bc_groups(i,1:n_seq(i))
            do j = 1,n_seq(i)
               bc_groups(i,j) = tmp1(n_seq(i)-j+1)
            enddo
            deallocate ( tmp1 )
         end do
      endif

c     Sort the schedules within groups to reflect the internal dependencies
      do i = 1,N_groups
         do j = 1,n_seq(i)
c           Get the dependencies between the remaining members of the group
            nn = n_seq(i)-j+1
            allocate ( tmp1(nn) )
            allocate ( tmp2(nn,nn) )
            tmp1 = bc_groups(i,j:n_seq(i))
            do k = 1,nn
               do m = 1,nn
                  tmp2(k,m) = tot_dep(tmp1(k),tmp1(m))
               end do
            end do
c           Check that no members of the group depend critically on any other members
            if ( MAXVAL(tmp2)==2 ) then
               call pstop('create_bc_schedule',
     &            'Error in dependency detected')
            end if
c           Find the first variable which does not depend on the others
            k = -1
            do m = nn,1,-1
               if ( MAXVAL(tmp2(m,:))==0 ) then
                  k = m
               endif
            enddo
            if ( k==-1 ) then 
c              The variables all depend weakly on each other, we only need to
c              obey the bc dependencies within a group and re-search
               do k = 1,nn
                  do m = 1,nn
                     tmp2(k,m) = bc_dep(tmp1(k),tmp1(m))
                  end do
                  tmp2(k,k) = 0
               end do
               do m = nn,1,-1
                  if ( MAXVAL(tmp2(m,:))==0 ) then
                     k = m
                  end if
               end do
            end if
            if ( k==-1) then
c              Unable to find a variable that does not depend on the others, a circular dependence is present
               write (*,*) 'Problem ordering groups, ', 
     &            'circular dependence detected between: '
               do m = 1,nn
                  call print_var_name(u,aux,tmp1(m))
               enddo
               call pstop('create_bc_schedule',
     &            'Erroor sorting group, circular dependence detected')
            end if
            bc_groups(i,j) = tmp1(k)
            bc_groups(i,j+1:j+k-1) = tmp1(1:k-1)
            bc_groups(i,j+k:n_seq(i)) = tmp1(k+1:nn)            
            deallocate ( tmp1 )
            deallocate ( tmp2 )
         end do
      end do
      
c     If we need inherite the dependencies from previous groups, add them
c     This is needed for temporary patches when the src and dst ids are not the same.  
c     In this case we assume all of the interior data is compied sucessfully, 
c     and we only need to worry about bc values on the temporary patch.
      if ( inherit ) then
         allocate ( tmp1(N) )
         do i = 2,N_groups
c           Get an array that indicates if each variable is needed by any member of the current group
            tmp1(:) = 0
            do j = 1,n_seq(i)
               k = bc_groups(i,j)
               if ( used(k) == 2 ) then
                  do m = 1,N
                     if ( bc_dep(k,m)>0 ) then
                        tmp1(m) = 1
                     endif
                  end do
               endif
            enddo
c           We need any dependencies that are necessary to fill the bc
c           for each dependency (interior dependencies should not be needed)
            test = .TRUE.
            do while ( test ) 
               test = .FALSE.
               do j = 1,N
                  if ( tmp1(j)>0 ) then
                     do k = 1,N
                        if ( bc_dep(j,k)>0 .AND. tmp1(k)==0 
     &                        .AND. used(k)==2 ) then
                           tmp1(k) = 1
                           test = .TRUE.
                        endif
                     enddo
                  endif
               enddo
            enddo
c           Remove any dependency variables from the current group (these are met by being members of the group)
            do k = 1,n_seq(i)
               m = bc_groups(i,k)
               tmp1(m) = 0
            enddo
c           Make space for the new variables we will add to the group
            nn = SUM(tmp1)       ! This is the number of variables we will be adding to the group
            do k = n_seq(i),1,-1
               bc_groups(i,k+nn) = bc_groups(i,k)
            enddo
            bc_groups(i,1:nn) = 0
            n_seq(i) = n_seq(i)+nn
c           Add the variables to the group preserving their orders
            nn = 1
            do j = 1,i-1
               do k = 1,n_seq(j)
                  m = bc_groups(j,k)
                  if ( tmp1(m)>0 ) then
                     bc_groups(i,nn) = m
                     tmp1(m) = 0
                     nn = nn+1
                  endif
               enddo
            enddo
c           Check the group
            if ( MINVAL(bc_groups(i,1:n_seq(i)))<1 ) then
               call pstop('create_bc_schedule',
     &            'Error inheriting dependencies')
            endif
         enddo
         deallocate ( tmp1 )
      endif

c     Deallocate existing 
      if ( associated(bc_grp) ) then
         do j = 1,size(bc_grp)
            if ( associated(bc_grp(j)%bc_seq) ) then
               deallocate ( bc_grp(j)%bc_seq )
            endif
         enddo
         deallocate ( bc_grp )
      endif

c     Create the output groups
      nbc_grp = N_groups
      allocate ( bc_grp(nbc_grp) )
      do i = 1,nbc_grp
         bc_grp(i)%nbc_seq = n_seq(i)
         allocate ( bc_grp(i)%bc_seq(n_seq(i),3) )
         do j = 1,n_seq(i)
            if ( bc_groups(i,j)<=n_var(1) ) then
               bc_grp(i)%bc_seq(j,1) = bc_groups(i,j)
               bc_grp(i)%bc_seq(j,2) = 0
            else if ( bc_groups(i,j)<=SUM(n_var(1:2)) ) then
               bc_grp(i)%bc_seq(j,1) = bc_groups(i,j)-n_var(1)
               bc_grp(i)%bc_seq(j,2) = 1
            else if ( bc_groups(i,j)<=SUM(n_var(1:3)) ) then
               bc_grp(i)%bc_seq(j,1) = -(bc_groups(i,j)-SUM(n_var(1:2)))
               bc_grp(i)%bc_seq(j,2) = 0
            else if ( bc_groups(i,j)<=SUM(n_var(1:4)) ) then
               bc_grp(i)%bc_seq(j,1) = -(bc_groups(i,j)-SUM(n_var(1:3)))
               bc_grp(i)%bc_seq(j,2) = 1
            else
               call pstop('create_bc_schedule','This should not occur')
            endif
            if ( used(bc_groups(i,j))==1 ) then
               bc_grp(i)%bc_seq(j,3) = 0
            else if ( used(bc_groups(i,j))==2 ) then
               bc_grp(i)%bc_seq(j,3) = 1
            else
               call pstop('create_bc_schedule','This should not occur')
            endif
         enddo
      enddo

      if (print_used) then
        write (*,*) 'BC sequence'
        do i = 1,nbc_grp
          write (*,*) 'Group=',i
          do j = 1,n_seq(i)
            call print_var_name(u,aux,bc_groups(i,j))
            if (bc_grp(i)%bc_seq(j,2) == 0) then
              write (*,*) '    Scalar=',bc_grp(i)%bc_seq(j,1)
            else
              write (*,*) '    Vector=',bc_grp(i)%bc_seq(j,1)
            endif
          enddo
          write (*,*)
        enddo
      endif

c     Free temporary memory
      deallocate ( used )
      deallocate ( bc_dep )
      deallocate ( dom_dep )
      deallocate ( tot_dep )
      deallocate ( n_seq )
      deallocate ( bc_groups )

      contains

c     check_dep_list
c     ################################################################
      integer function check_dep_list(N,dep_list)

      implicit none

      INTEGER :: dep_list(N,3)
      INTEGER i, j, N

      check_dep_list = 0
      i = 1
      do while ( dep_list(i,1).ne.0 ) 
         if ( dep_list(i,1)>0 ) then
            if ( dep_list(i,2).ne.0 ) then
               check_dep_list = 1       ! Invalid entry found
            endif
         else if ( dep_list(i,1)<0 ) then
            if ( dep_list(i,2)<0 .or. dep_list(i,2)<0 ) then
               check_dep_list = 1       ! Invalid entry found
            endif
         endif
         do j = 1,i-1
            if ( (dep_list(i,1)==dep_list(j,1)) .AND. 
     &           (dep_list(i,2)==dep_list(j,2)) ) then
               check_dep_list = 2       ! Duplicate entry found
            endif
         enddo
         i = i+1
      enddo

      end function check_dep_list

c     get_dep_index
c     ################################################################
      function get_dep_index(dep_var,n_var)

      implicit none

      INTEGER :: get_dep_index
      INTEGER, INTENT(IN) :: dep_var(3)
      INTEGER, INTENT(IN) :: n_var(4)
      get_dep_index = -1
      if ( dep_var(1)>0 .AND. dep_var(2)==0 ) then
c        Dependency is a scalar dependent variable
         get_dep_index = dep_var(1)
      else if ( dep_var(1)>0 .AND. dep_var(2)==1 ) then
c        Dependency is a vector dependent variable
         get_dep_index = n_var(1) + dep_var(1)
      else if ( dep_var(1)<0 .AND. dep_var(2)==0 ) then
c        Dependency is a vector dependent variable
         get_dep_index = n_var(1) + n_var(2) - dep_var(1)
      else if ( dep_var(1)<0 .AND. dep_var(2)==1 ) then
c        Dependency is a vector dependent variable
         get_dep_index = n_var(1) + n_var(2) + n_var(3) - dep_var(1)
      endif
      if ( get_dep_index<=0 ) then
c        Bad case
         call pstop('get_dep_index','Bad index')
      endif
      end function get_dep_index

c     print_var_name
c     ################################################################
      subroutine print_var_name(u,aux,k)

      use variable_setup, only: var_array, aux_array

      implicit none
      type(var_array) :: u
      type(aux_array) :: aux
      INTEGER :: k, n_var(4), N
      n_var(1) = u%nvar     ! Number of scalar dependent variables
      n_var(2) = 0          ! Number of vector dependent variables
      n_var(3) = aux%nvar   ! Number of scalar auxiliary variables
      n_var(4) = aux%nvec   ! Number of vector auxiliary variables
      N = n_var(1) + n_var(2) + n_var(3) + n_var(4)
      if ( k <= n_var(1) ) then
         write (*,*) '   ', u%array_var(k)%descr
      else if ( k <= SUM(n_var(1:2)) ) then
         !!!! Place holder for dependent vector
      else if ( k <= SUM(n_var(1:3)) ) then
         write (*,*) '   ', aux%var_list(k-SUM(n_var(1:2)))%descr
      else if ( k <= SUM(n_var(1:4)) ) then
         write (*,*) '   ', aux%vec_list(k-SUM(n_var(1:3)))%descr
      endif
      end subroutine print_var_name

c     pstop
c     ################################################################
      subroutine pstop(routine,message)

c     ---------------------------------------------------------------
c     Stops program at "routine" with "message"
c     ---------------------------------------------------------------

        implicit none

        character(*)  :: routine, message

c     Begin program

          write (*,*)
          write (*,*) trim(message)
          write (*,*) 'Program stopped at routine ',trim(routine)

        call abort

        stop

      end subroutine pstop

      end subroutine create_bc_schedule
