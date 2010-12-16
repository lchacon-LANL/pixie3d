
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

cc      10/20/2010: I don't think this call is needed, since it seems to be 
cc      called stand-alone by the SAMRAI driver before this routine
cc      (and initializeAuxVar below, which also needs it).

      call setupVarInitSeq(gv,it)

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

c Set BC schedule (see imposeBC_mod.F)

      call setBCsched(it)

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

      end subroutine applyBC
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

cc        varray%array_var(IAX)%array = acov(:,:,:,1)
cc        varray%array_var(IAY)%array = acov(:,:,:,2)
cc        varray%array_var(IAZ)%array = acov(:,:,:,3)

        call XferBCs(acov(:,:,:,1),varray%array_var(IAX)%array)
        call XferBCs(acov(:,:,:,2),varray%array_var(IAY)%array)
        call XferBCs(acov(:,:,:,3),varray%array_var(IAZ)%array)

#else
      case(IBX,IBY,IBZ)

cc        varray%array_var(IBX)%array = bcnv(:,:,:,1)
cc        varray%array_var(IBY)%array = bcnv(:,:,:,2)
cc        varray%array_var(IBZ)%array = bcnv(:,:,:,3)

        call XferBCs(bcnv(:,:,:,1),varray%array_var(IBX)%array)
        call XferBCs(bcnv(:,:,:,2),varray%array_var(IBY)%array)
        call XferBCs(bcnv(:,:,:,3),varray%array_var(IBZ)%array)

#endif

      case(IVX,IVY,IVZ)

        if (.not.nc_eom_v) then
cc          varray%array_var(IVX)%array = vcnv(:,:,:,1)
cc     .                                 *varray%array_var(IRHO)%array
cc          varray%array_var(IVY)%array = vcnv(:,:,:,2)
cc     .                                 *varray%array_var(IRHO)%array
cc          varray%array_var(IVZ)%array = vcnv(:,:,:,3)
cc     .                                 *varray%array_var(IRHO)%array

          rho => varray%array_var(IRHO)%array

          varray%array_var(IVX)%array(0    ,:,:) = vcnv(0    ,:,:,1)
     .                                            *rho (0    ,:,:)
          varray%array_var(IVY)%array(0    ,:,:) = vcnv(0    ,:,:,2)
     .                                            *rho (0    ,:,:)
          varray%array_var(IVZ)%array(0    ,:,:) = vcnv(0    ,:,:,3)
     .                                            *rho (0    ,:,:)

          varray%array_var(IVX)%array(nnx+1,:,:) = vcnv(nnx+1,:,:,1)
     .                                            *rho (nnx+1,:,:)
          varray%array_var(IVY)%array(nnx+1,:,:) = vcnv(nnx+1,:,:,2)
     .                                            *rho (nnx+1,:,:)
          varray%array_var(IVZ)%array(nnx+1,:,:) = vcnv(nnx+1,:,:,3)
     .                                            *rho (nnx+1,:,:)

          varray%array_var(IVX)%array(:,0    ,:) = vcnv(:,0    ,:,1)
     .                                            *rho (:,0    ,:)
          varray%array_var(IVY)%array(:,0    ,:) = vcnv(:,0    ,:,2)
     .                                            *rho (:,0    ,:)
          varray%array_var(IVZ)%array(:,0    ,:) = vcnv(:,0    ,:,3)
     .                                            *rho (:,0    ,:)

          varray%array_var(IVX)%array(:,nny+1,:) = vcnv(:,nny+1,:,1)
     .                                            *rho (:,nny+1,:)
          varray%array_var(IVY)%array(:,nny+1,:) = vcnv(:,nny+1,:,2)
     .                                            *rho (:,nny+1,:)
          varray%array_var(IVZ)%array(:,nny+1,:) = vcnv(:,nny+1,:,3)
     .                                            *rho (:,nny+1,:)

          varray%array_var(IVX)%array(:,:,0    ) = vcnv(:,:,0    ,1)
     .                                            *rho (:,:,0    )
          varray%array_var(IVY)%array(:,:,0    ) = vcnv(:,:,0    ,2)
     .                                            *rho (:,:,0    )
          varray%array_var(IVZ)%array(:,:,0    ) = vcnv(:,:,0    ,3)
     .                                            *rho (:,:,0    )

          varray%array_var(IVX)%array(:,:,nnz+1) = vcnv(:,:,nnz+1,1)
     .                                            *rho (:,:,nnz+1)
          varray%array_var(IVY)%array(:,:,nnz+1) = vcnv(:,:,nnz+1,2)
     .                                            *rho (:,:,nnz+1)
          varray%array_var(IVZ)%array(:,:,nnz+1) = vcnv(:,:,nnz+1,3)
     .                                            *rho (:,:,nnz+1)

          nullify(rho)
        else
cc          varray%array_var(IVX)%array = vcnv(:,:,:,1)
cc          varray%array_var(IVY)%array = vcnv(:,:,:,2)
cc          varray%array_var(IVZ)%array = vcnv(:,:,:,3)

          call XferBCs(vcnv(:,:,:,1),varray%array_var(IVX)%array)
          call XferBCs(vcnv(:,:,:,2),varray%array_var(IVY)%array)
          call XferBCs(vcnv(:,:,:,3),varray%array_var(IVZ)%array)
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

      integer :: iigx,iigy,iigz,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

      integer,optional :: order

c Local variables

      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz
     .          ,ig,jg,kg

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),grad_p(3),etal

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

cc        if (lag_pinch_bc) then
cc          rho => u_n%array_var(IRHO)%array
cc          tmp => u_n%array_var(ITMP)%array
cc        else
          rho => varray%array_var(IRHO)%array
          tmp => varray%array_var(ITMP)%array
cc        endif

        allocate(prs(0:nnx+1,0:nny+1,0:nnz+1))

        if (adiabatic) then
          prs = a_p*tmp
        else
          prs = a_p*rho*tmp
        endif

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
            call XferBCs(u_0%array_var(IVX)%array,v0(:,:,:,1))
            call XferBCs(u_0%array_var(IVY)%array,v0(:,:,:,2))
            call XferBCs(u_0%array_var(IVZ)%array,v0(:,:,:,3))
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
     .              ,iorder=1)

        case(IJCOV_0)  !Done when computing jcov_0

cc          call XformVector_BC(iigx,jcnv_0,jcov_0,.false.)

        case(IENI_0)

#if defined(vec_pot)
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,acov_0)
#else
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,jcov_0)
#endif

        case(IVECOV_N)

          if (save_tn.and.(.not.ion_hall)) call XferBCs(vecov,vecov_n)

        case(IVCOV_N)

          if (save_tn.and.ion_hall) call XferBCs(vcov,vcov_n)

        case(IVCNV_N)

          if (save_tn.and.ion_hall) call XferBCs(vcnv,vcnv_n)

        case(IVCNV)

c       Velocity BC

          if (pinch_flow) then
cc            if (lag_pinch_bc) then
cc              call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0)
cc            else
              call ExB_flow_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0)
cc            endif
          else
            call XferBCs(vcnv_0,v0)
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
cc            if (lag_pinch_bc) then
cc              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0
cc     .                       ,first_order=.true.)
cc            else
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
     .                       ,first_order=.true.)
cc            endif
          else
            call XferBCs(jcnv_0,v0)
          endif

          call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0
     .              ,vaux%vec_list(IBCNV)%bconds,iigx,iigy,iigz
     .              ,iorder=order,is_cnv=vaux%vec_list(IBCNV)%cnv)

        case(IBCOV) !Done when computing bcnv

cc          call XformVector_BC(iigx,bcnv,bcov,.false.)

        case(IJCNV)

c       Current BC

          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation
cc            if (lag_pinch_bc) then
cc              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv_n,prs,eeta,v0
cc     .                       ,first_order=.true.)
cc            else
              call J_pinch_bc(iigx,nnx,nny,nnz,bcnv  ,prs,eeta,v0
     .                       ,first_order=.true.)
cc            endif
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

          call setBC(IJX,3,nnx,nny,nnz,jcnv,jcov,v0
     .              ,vaux%vec_list(IJCNV)%bconds,iigx,iigy,iigz
     .              ,iorder=1)

        case(IJCOV) !Done when computing jcnv

cc          call XformVector_BC(iigx,jcnv,jcov,.false.)

        case(IVECNV)

c       Electron velocity BC: ve=v-di*j/rho with j=curl(B)

          vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv(:,:,:,1)/rho
          vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv(:,:,:,2)/rho
          vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv(:,:,:,3)/rho

        case(IVECOV)

          call XformVector_BC(iigx,vecnv,vecov,.false.)

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

cc#if defined(vec_pot)
cc          E_ni = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,acov)
cc#else
cc          E_ni = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,jcov)
cc#endif

          E_ni = E_ni - E_ni_0 + di*E_h  !Substract imposed electric field

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
          eeta_0=res(nnx,nny,nnz,iigx,iigy,iigz,psi)
#else
          eeta_0=res(nnx,nny,nnz,iigx,iigy,iigz,tmp_0)
#endif
        case(INU)

c       Viscosity (needed for BCs; may depend on T)
 
#if defined(RFX)
          nuu =     vis(nnx,nny,nnz,iigx,iigy,iigz)
#else
          nuu = rho*vis(nnx,nny,nnz,iigx,iigy,iigz)
#endif
        case(IETA)

c       Resistivity (needed for BCs; may depend on T)
 
#if defined(RFX)
          eeta= res(nnx,nny,nnz,iigx,iigy,iigz,psi)
#else
          eeta= res(nnx,nny,nnz,iigx,iigy,iigz,tmp)
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
     .                      ,vcnv_0,bcnv_0,acov_0)
#else
          E_ni_0 = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz
     .                      ,vcnv_0,bcnv_0,jcov_0)
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

          vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv(:,:,:,1)/rho
          vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv(:,:,:,2)/rho
          vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv(:,:,:,3)/rho

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
          E_ni = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,acov)
#else
          E_ni = E_fld_res(nnx,nny,nnz,iigx,iigy,iigz,vcnv,bcnv,jcov)
#endif

cc          E_ni = E_ni  - E_ni_0 + di*E_h  !Done in applyBC

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
