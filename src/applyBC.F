c setAppBCs
c ########################################################################
      subroutine setAppBCs(varray,vaux)
c ------------------------------------------------------------------------
c     Defines boundary conditions of all application variables, both
c     dependent and auxiliary. Convention:
c           + bbcs(1) ---> at x0
c           + bbcs(2) ---> at x1
c           + bbcs(3) ---> at y0
c           + bbcs(4) ---> at y1
c           + bbcs(5) ---> at z0
c           + bbcs(6) ---> at z1
c     Definition of BC identifiers is given in "grid_mod.f". In vectors,
c     a negative BC identifier means that BCs are to be imposed on
c     covariant components instead of on (default) contravariant comps.
c ------------------------------------------------------------------------

      use local_BCS_variables

      implicit none

c Call variables

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

c Local variables

c Begin program

c Setup dependent variables

      call defineDepVariables(varray)

c Setup auxiliary variables

      call defineAuxVariables(vaux)

c End program

      end subroutine setAppBCs

c destroyAppBCs
c ########################################################################
      subroutine destroyAppBCs
c ------------------------------------------------------------------------
c     Destroys variables for application BC treatment.
c ------------------------------------------------------------------------

      use local_BCS_variables

      implicit none

c Call variables

c Local variables

c Begin program

      call deallocLocalAuxVars

c End program

      end subroutine destroyAppBCs

c applyDepVarBC
c####################################################################
      subroutine applyDepVarBC(igr,bc_inf,varray,order)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use var_setup

      implicit none

c Call variables

      integer :: igr,bc_inf(2)

      type(var_array),pointer :: varray

      integer,optional :: order

c Local variables

      integer :: nnx,nny,nnz
      real(8),pointer,dimension(:,:,:) :: jac

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(igr) 
      nny = grid_params%nyv(igr)
      nnz = grid_params%nzv(igr)

      jac => gmetric%grid(igr)%jac

c Select variable

      select case(abs(bc_inf(1)))
      case(IRHO)

c     Density BCs

cc        call setBC(IRHO,nnx,nny,nnz
cc     .            ,varray%array_var(IRHO)%array
cc     .            ,u_0   %array_var(IRHO)%array
cc     .            ,varray%array_var(IRHO)%bconds
cc     .            ,igr,igr,igr,iorder=order)

cc        call XferBCs(rho*jac,varray%array_var(IRHO)%array,.false.)
cc
cc        varray%array_var(IRHO)%array(0    ,:,:) = jac(0    ,:,:)
cc     .                                           *rho(0    ,:,:)
cc        varray%array_var(IRHO)%array(nnx+1,:,:) = jac(nnx+1,:,:)
cc     .                                           *rho(nnx+1,:,:)
cc        varray%array_var(IRHO)%array(:,0    ,:) = jac(:,0    ,:)
cc     .                                           *rho(:,0    ,:)
cc        varray%array_var(IRHO)%array(:,nny+1,:) = jac(:,nny+1,:)
cc     .                                           *rho(:,nny+1,:)
cc        varray%array_var(IRHO)%array(:,:,0    ) = jac(:,:,0    )
cc     .                                           *rho(:,:,0    )
cc        varray%array_var(IRHO)%array(:,:,nnz+1) = jac(:,:,nnz+1)
cc     .                                           *rho(:,:,nnz+1)

#if defined(vec_pot)
      case(IAX)

cc        call XferBCs(acov(:,:,:,1),varray%array_var(IAX)%array,.false.)

      case(IAY)

cc        call XferBCs(acov(:,:,:,2),varray%array_var(IAY)%array,.false.)

      case(IAZ)

cc        call XferBCs(acov(:,:,:,3),varray%array_var(IAZ)%array,.false.)
#else
      case(IBX)

cc        call XferBCs(bcnv(:,:,:,1),varray%array_var(IBX)%array,.false.)

      case(IBY)

cc        call XferBCs(bcnv(:,:,:,2),varray%array_var(IBY)%array,.false.)

      case(IBZ)

cc        call XferBCs(bcnv(:,:,:,3),varray%array_var(IBZ)%array,.false.)
#endif

      case(IVX)

cc        if (.not.nc_eom_v) then
cc
cc          varray%array_var(IVX)%array(0    ,:,:) = vcnv(0    ,:,:,1)
cc     .                                            *rho (0    ,:,:)
cc          varray%array_var(IVX)%array(nnx+1,:,:) = vcnv(nnx+1,:,:,1)
cc     .                                            *rho (nnx+1,:,:)
cc          varray%array_var(IVX)%array(:,0    ,:) = vcnv(:,0    ,:,1)
cc     .                                            *rho (:,0    ,:)
cc          varray%array_var(IVX)%array(:,nny+1,:) = vcnv(:,nny+1,:,1)
cc     .                                            *rho (:,nny+1,:)
cc          varray%array_var(IVX)%array(:,:,0    ) = vcnv(:,:,0    ,1)
cc     .                                            *rho (:,:,0    )
cc          varray%array_var(IVX)%array(:,:,nnz+1) = vcnv(:,:,nnz+1,1)
cc     .                                            *rho (:,:,nnz+1)
cc
cc        Else
cc
cc          call XferBCs(vcnv(:,:,:,1),varray%array_var(IVX)%array
cc     .                ,.false.)
cc
cc        endif

      case(IVY)

cc        if (.not.nc_eom_v) then
cc
cc          varray%array_var(IVY)%array(0    ,:,:) = vcnv(0    ,:,:,2)
cc     .                                            *rho (0    ,:,:)
cc          varray%array_var(IVY)%array(nnx+1,:,:) = vcnv(nnx+1,:,:,2)
cc     .                                            *rho (nnx+1,:,:)
cc          varray%array_var(IVY)%array(:,0    ,:) = vcnv(:,0    ,:,2)
cc     .                                            *rho (:,0    ,:)
cc          varray%array_var(IVY)%array(:,nny+1,:) = vcnv(:,nny+1,:,2)
cc     .                                            *rho (:,nny+1,:)
cc          varray%array_var(IVY)%array(:,:,0    ) = vcnv(:,:,0    ,2)
cc     .                                            *rho (:,:,0    )
cc          varray%array_var(IVY)%array(:,:,nnz+1) = vcnv(:,:,nnz+1,2)
cc     .                                            *rho (:,:,nnz+1)
cc
cc        else
cc
cc          call XferBCs(vcnv(:,:,:,2),varray%array_var(IVY)%array
cc     .                ,.false.)
cc
cc        endif
cc
cc        if (alt__eom()) then
cc          varray%array_var(IVY)%array(0    ,:,:) =
cc     .         varray%array_var(IVY)%array(0    ,:,:)*jac (0    ,:,:)
cc          varray%array_var(IVY)%array(nnx+1,:,:) =
cc     .         varray%array_var(IVY)%array(nnx+1,:,:)*jac (nnx+1,:,:)
cc          varray%array_var(IVY)%array(:,0    ,:) =
cc     .         varray%array_var(IVY)%array(:,0    ,:)*jac (:,0    ,:)
cc          varray%array_var(IVY)%array(:,nny+1,:) =
cc     .         varray%array_var(IVY)%array(:,nny+1,:)*jac (:,nny+1,:)
cc          varray%array_var(IVY)%array(:,:,0    ) =
cc     .         varray%array_var(IVY)%array(:,:,0    )*jac (:,:,0    )
cc          varray%array_var(IVY)%array(:,:,nnz+1) =
cc     .         varray%array_var(IVY)%array(:,:,nnz+1)*jac (:,:,nnz+1)
cc        endif

       case(IVZ)

cc        if (.not.nc_eom_v) then
cc
cc          varray%array_var(IVZ)%array(0    ,:,:) = vcnv(0    ,:,:,3)
cc     .                                            *rho (0    ,:,:)
cc          varray%array_var(IVZ)%array(nnx+1,:,:) = vcnv(nnx+1,:,:,3)
cc     .                                            *rho (nnx+1,:,:)
cc          varray%array_var(IVZ)%array(:,0    ,:) = vcnv(:,0    ,:,3)
cc     .                                            *rho (:,0    ,:)
cc          varray%array_var(IVZ)%array(:,nny+1,:) = vcnv(:,nny+1,:,3)
cc     .                                            *rho (:,nny+1,:)
cc          varray%array_var(IVZ)%array(:,:,0    ) = vcnv(:,:,0    ,3)
cc     .                                            *rho (:,:,0    )
cc          varray%array_var(IVZ)%array(:,:,nnz+1) = vcnv(:,:,nnz+1,3)
cc     .                                            *rho (:,:,nnz+1)
cc
cc        else
cc
cc          call XferBCs(vcnv(:,:,:,3),varray%array_var(IVZ)%array
cc     .                ,.false.)
cc
cc        endif

      case(ITMP)

c     Temperature BCs

cc        call setBC(ITMP,nnx,nny,nnz
cc     .            ,varray%array_var(ITMP)%array
cc     .            ,u_0   %array_var(ITMP)%array
cc     .            ,varray%array_var(ITMP)%bconds
cc     .            ,igr,igr,igr,iorder=order)

cc        call XferBCs(rho*tmp*jac,varray%array_var(ITMP)%array,.false.)

cc        varray%array_var(ITMP)%array(0    ,:,:) = jac (0    ,:,:)
cc     .                                           *tmpe(0    ,:,:)*a_p
cc        varray%array_var(ITMP)%array(nnx+1,:,:) = jac (nnx+1,:,:)
cc     .                                           *tmpe(nnx+1,:,:)*a_p
cc        varray%array_var(ITMP)%array(:,0    ,:) = jac (:,0    ,:)
cc     .                                           *tmpe(:,0    ,:)*a_p
cc        varray%array_var(ITMP)%array(:,nny+1,:) = jac (:,nny+1,:)
cc     .                                           *tmpe(:,nny+1,:)*a_p
cc        varray%array_var(ITMP)%array(:,:,0    ) = jac (:,:,0    )
cc     .                                           *tmpe(:,:,0    )*a_p
cc        varray%array_var(ITMP)%array(:,:,nnz+1) = jac (:,:,nnz+1)
cc     .                                           *tmpe(:,:,nnz+1)*a_p

      end select

c End

      end subroutine applyDepVarBC

c applyAuxVarBC
c####################################################################
      subroutine applyAuxVarBC(igr,bc_inf,varray,vaux,order)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use var_setup

c diag ****
      use app_iosetup
c diag ****

      implicit none

c Call variables

      integer :: igr,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

      integer,optional :: order

c Local variables

      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz
     .          ,ig,jg,kg

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),grad_p(3),etal

      real(8),allocatable,dimension(:,:,:,:) :: v0,vdum

      real(8),pointer,dimension(:,:,:) :: s0,prs0,prs1

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(igr) 
      nny = grid_params%nyv(igr)
      nnz = grid_params%nzv(igr)

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

        allocate(s0(0:nnx+1,0:nny+1,0:nnz+1))

        select case(abs(bc_inf(1)))
        case(IETA)

          eeta= res(igr,tmpe)

        case(ICHI)

          cchi= chi_perp(igr,tmpe)

        case(INU)

#if defined(RFX)
          nuu = vis(igr,tmp)
#else
          if (slava) then
            nuu = rho*vis(igr,tmp)
          else
            nuu = vis(igr,tmp)
          endif
#endif

        case(IETA_0)

          eeta_0=res(igr,tmpe_0)

        case(IHETA)

cc          h_eta=hres(nnx,nny,nnz,igr,igr,igr,bcnv_0,rho_0)
          h_eta=hres(igr,bcnv_0,rho)

        case(IRHO_0)

c       Density BCs

          call setBC(IRHO,nnx,nny,nnz,rho_0
     .              ,u_0%array_var(IRHO)%array
     .              ,vaux%var_list(IRHO_0)%bconds
     .              ,igr,igr,igr,iorder=order)

        case(ITMP_0)

c       Temperature BCs

          call setBC(ITMP,nnx,nny,nnz,tmp_0
     .              ,u_0%array_var(ITMP)%array
     .              ,vaux%var_list(ITMP_0)%bconds
     .              ,igr,igr,igr,iorder=order)

        case(ITMPE_0)

c       Temperature BCs

          call setBC(ITMP,nnx,nny,nnz,tmpe_0
     .              ,u_0%array_var(ITMP)%array
     .              ,vaux%var_list(ITMPE_0)%bconds
     .              ,igr,igr,igr,iorder=order)

        case(IRHO_AUX)

c       Density BCs

          call XferBCs(rho_0,s0,.true.)

          call setBC(IRHO,nnx,nny,nnz,rho,s0
     .              ,vaux%var_list(IRHO_AUX)%bconds
     .              ,igr,igr,igr,iorder=order)

        case(ITMP_AUX)

c       Temperature BCs

          call XferBCs(tmp_0,s0,.true.)

          call setBC(ITMP,nnx,nny,nnz,tmp,s0
     .              ,vaux%var_list(ITMP_AUX)%bconds
     .              ,igr,igr,igr,iorder=order)

        case(ITMPE_AUX)

c       Temperature BCs

          call XferBCs(tmpe_0,s0,.true.)

          call setBC(ITMP,nnx,nny,nnz,tmpe,s0
     .              ,vaux%var_list(ITMPE_AUX)%bconds
     .              ,igr,igr,igr,iorder=order)

        case default
          call pstop('applyAuxVarBC'
     .              ,'Scalar auxiliary variable does not exist')
        end select

        deallocate(s0)

      case(VECTOR)

        allocate(prs0(0:nnx+1,0:nny+1,0:nnz+1)
     .          ,prs1(0:nnx+1,0:nny+1,0:nnz+1))

c       Allocate auxiliary variables in local domain

        allocate(v0  (0:nnx+1,0:nny+1,0:nnz+1,3)
     .          ,vdum(0:nnx+1,0:nny+1,0:nnz+1,3))

        v0   = 0d0
        vdum = 0d0

        select case(abs(bc_inf(1)))
#if defined(vec_pot)
        case(IACOV_0)

          vdum = XformVector(igr,acov_0,.true.)
          call setBC(IAX,3,nnx,nny,nnz,vdum,acov_0,v0
     .              ,vaux%vec_list(IACOV_0)%bconds,igr,igr,igr
     .              ,is_cnv=vaux%vec_list(IACOV_0)%cnv
     .              ,iorder=2)
#endif

        case(IBCNV_0)

          call setBC(IBX,3,nnx,nny,nnz,bcnv_0,bcov_0,v0
     .              ,vaux%vec_list(IBCNV_0)%bconds,igr,igr,igr
     .              ,is_cnv=vaux%vec_list(IBCNV_0)%cnv
     .              ,iorder=1)

        case(IBCOV_0)   !Done with bcnv_0

cc          call XformVector_BC(igr,bcnv_0,bcov_0,.false.)

        case(IVCNV_0)

c       Velocity BC

          if (pinch_flow) then
            prs0 = rho_0*tmp_0

            call ExB_flow_bc(igr,nnx,nny,nnz,bcnv_0,prs0,eeta_0,v0)

            !Fill ghost nodes
            vdum = XformVector(igr,vcnv_0,.false.)
            call setBC(IVX,3,nnx,nny,nnz,vcnv_0,vdum,v0
     .                ,vaux%vec_list(IVCNV_0)%bconds
     .                ,igr,igr,igr,iorder=order)
          else
            !Extrapolation
            vdum = XformVector(igr,vcnv_0,.false.)
            call setBC(IVX,3,nnx,nny,nnz,vcnv_0,vdum,v0
     .                ,vaux%vec_list(IVCNV_0)%bconds
     .                ,igr,igr,igr,iorder=order)

          endif

        case(IJCNV_0)

c       Current BC

          call setBC(IJX,3,nnx,nny,nnz,jcnv_0,jcov_0,v0
     .              ,vaux%vec_list(IJCNV_0)%bconds,igr,igr,igr
     .              ,is_cnv=vaux%vec_list(IJCNV_0)%cnv
     .              ,iorder=1)

        case(IJCOV_0)  !Done when computing jcov_0

cc          call XformVector_BC(igr,jcnv_0,jcov_0,.false.)

        case(IENI_0)

#if defined(vec_pot)
          E_ni_0 = E_fld_res(igr,vcnv_0,bcnv_0,acov_0,eeta_0)
#else
          E_ni_0 = E_fld_res(igr,vcnv_0,bcnv_0,jcov_0,eeta_0)
#endif

        case(IVECOV_N)

cc          if (save_tn.and.(.not.ion_hall)) then
cc            call XferBCs(vecov,vecov_n,.false.)
cc          endif

cc        case(IVECOV_NM)

cc        case(IVCNV_NM)

cc          if (save_tn.and.ion_hall) call XferBCs(vcov,vcov_n,.false.)

        case(IVCNV_N)

cc          if (save_tn.and.ion_hall) call XferBCs(vcnv,vcnv_n,.false.)

        case(IJCNV_N)

cc          if (save_tn.and.ion_hall) call XferBCs(jcnv,jcnv_n,.false.)

        case(IVCNV)

c       Velocity BC

          if (pinch_flow) then
            prs1 = rho*tmp
            call ExB_flow_bc(igr,nnx,nny,nnz,bcnv,prs1,eeta,v0)
          else
            call XferBCs(vcnv_0,v0,.true.)
          endif

          call setBC(IVX,3,nnx,nny,nnz,vcnv,vcov,v0
     .              ,vaux%vec_list(IVCNV)%bconds,igr,igr,igr
     .              ,iorder=order)

        case(IVCOV)   !Done when computing vcnv

cc          call XformVector_BC(igr,vcnv,vcov,.false.)

#if defined(vec_pot)

        case(IACOV)

c       Vector Pot. BC

          call XferBCs(acov_0,v0,.true.)

          call setBC(IAX,3,nnx,nny,nnz,acnv,acov,v0
     .              ,vaux%vec_list(IACOV)%bconds,igr,igr,igr
     .              ,iorder=order,is_cnv=vaux%vec_list(IACOV)%cnv)

        case(IACNV)  !Done when computing acov

cc          call XformVector_BC(igr,acov,acnv,.true.)

#endif
        case(IBCNV)

c       Magnetic field BC

          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation
            prs1 = rho*tmp
            call J_pinch_bc(igr,nnx,nny,nnz,bcnv,prs1,eeta,v0)
          else
            call XferBCs(jcnv_0,v0,.true.)
          endif

          call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0
     .              ,vaux%vec_list(IBCNV)%bconds,igr,igr,igr
     .              ,iorder=order,is_cnv=vaux%vec_list(IBCNV)%cnv)

        case(IBCOV) !Done when computing bcnv

cc          call XformVector_BC(igr,bcnv,bcov,.false.)

        case(IJCNV)

c       Current BC

          if (pinch_flow) then  !Jt = alpha*Bt in cnv representation
            prs1 = rho*tmp
            call J_pinch_bc(igr,nnx,nny,nnz,bcnv,prs1,eeta,v0)
          else
            call XferBCs(jcnv_0,v0,.true.)
          endif

          call setBC(IJX,3,nnx,nny,nnz,jcnv,jcov,v0
     .              ,vaux%vec_list(IJCNV)%bconds,igr,igr,igr
     .              ,iorder=order)

        case(IJCOV) !Done when computing jcnv

cc          call XformVector_BC(igr,jcnv,jcov,.false.)

        case(IVECNV)

c       Electron velocity BC: ve=v-di*j/rho with j=curl(B)

          vecnv = form_ve(vcnv,jcnv,rho)

        case(IVECOV)

cc          call XformVector_BC(igr,vecnv,vecov,.false.)

        case(IVEFCNV)

c       Current BC on fake current

          if (fake_ve) then
cc            if (ion_hall) then
cc              v0 = 0d0
cc
cc              bcnd = vaux%vec_list(IVEFCNV)%bconds
cc              vdum = XformVector(igr,vefcnv,.false.)
cc              call setBC(IAX,3,nnx,nny,nnz,vefcnv,vdum,v0,bcnd
cc     .                ,igr,igr,igr,iorder=order
cc     .                ,is_cnv=vaux%vec_list(IVEFCNV)%cnv)
cc
cc            else
              !BCs on j=curl(curl(A))
              if (pinch_flow) then  !Jt = alpha*Bt in cnv representation
                prs1 = rho*tmp
                call J_pinch_bc(igr,nnx,nny,nnz,bcnv,prs1,eeta,v0)
              else
                call XferBCs(jcnv_0,v0,.true.)
              endif

              bcnd = vaux%vec_list(IVEFCNV)%bconds
              vdum = XformVector(igr,vefcnv,.false.)
cc            call setBC(IJX,3,nnx,nny,nnz,vefcnv,vdum,v0,bcnd
cc            call setBC(IVX,3,nnx,nny,nnz,vefcnv,vdum,v0,bcnd
              call setBC(IAX,3,nnx,nny,nnz,vefcnv,vdum,v0,bcnd
     .                ,igr,igr,igr,iorder=order
     .                ,is_cnv=vaux%vec_list(IVEFCNV)%cnv)
cc            endif
          else
            vefcnv = jcnv
          endif

        case(IDIVPI)

#if !defined(vec_pot)
          if (di > 0d0.and.ion_hall) then
            bcnd = vaux%vec_list(IDIVPI)%bconds
            vdum = XformVector(igr,div_pi,.false.)
            call setBC(IVX,3,nnx,nny,nnz,div_pi,vdum,v0,bcnd
     .                ,igr,igr,igr,iorder=order
     .                ,is_cnv=vaux%vec_list(IDIVPI)%cnv)
          endif
#endif

        case(IDIVPE)

#if !defined(vec_pot)
          if (di > 0d0.and.(.not.ion_hall)) then
            bcnd = vaux%vec_list(IDIVPE)%bconds
            vdum = XformVector(igr,div_pe,.false.)
            call setBC(IVX,3,nnx,nny,nnz,div_pe,vdum,v0,bcnd
     .                ,igr,igr,igr,iorder=order
     .                ,is_cnv=vaux%vec_list(IDIVPE)%cnv)
          endif
#endif

        case(IEH)

c       Hall electric field BC

#if !defined(vec_pot)
          if (di > 0d0) then
            bcnd = vaux%vec_list(IEH)%bconds
            vdum = XformVector(igr,E_h,.true.)
            call setBC(IVX,3,nnx,nny,nnz,vdum,E_h,v0,bcnd
     .                ,igr,igr,igr,iorder=order,is_cnv=.false.)
          endif
#endif

        case(IENI)

c       Total electric field BC

#if defined(vec_pot)
          E_ni=E_fld_res(igr,vcnv,bcnv,acov,eeta)
#else
          E_ni=E_fld_res(igr,vcnv,bcnv,jcov,eeta)
#endif

cc          E_ni = E_ni - E_ni_0  !Substract imposed electric field

          if (di > 0d0) E_ni = E_ni + di*E_h !Add Hall electric field

        case default
          write (*,*) abs(bc_inf(1))

          call pstop('applyAuxVarBC'
     .              ,'Vector auxiliary variable does not exist')

        end select

c       Deallocate variables 

        deallocate(v0,vdum,prs0,prs1)

      end select

c End

      end subroutine applyAuxVarBC

c preproc_AuxVar
c####################################################################
      subroutine preproc_AuxVar(igr,bc_inf,varray,vaux,order)
c--------------------------------------------------------------------
c     Preprocesses auxiliary variables
c--------------------------------------------------------------------

      use local_BCS_variables

      use auxiliaryVariables

      use var_setup

c diag ****
      use app_iosetup
c diag ****

      implicit none

c Call variables

      integer :: igr,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

      integer,optional :: order

c Local variables

      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz
     .          ,ig,jg,kg

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),grad_p(3),etal

      real(8),allocatable,dimension(:,:,:,:) :: v0

      real(8),pointer,dimension(:,:,:) :: jac

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(igr) 
      nny = grid_params%nyv(igr)
      nnz = grid_params%nzv(igr)

      jac => gmetric%grid(igr)%jac

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

        select case(abs(bc_inf(1)))
        case(IRHO_0)

          rho_0 = u_0%array_var(IRHO)%array/jac

        case(ITMP_0)

          tmp_0 = tmpe_0*a_p

        case(ITMPE_0)

          where (u_0%array_var(IRHO)%array /= 0d0)
            tmpe_0 = u_0%array_var(ITMP)%array
     .              /u_0%array_var(IRHO)%array/a_p
          end where

        case(IETA_0)

c       Resistivity at initial time
 
          eeta_0=res(igr,tmpe_0)

        case(IRHO_AUX)

          if ((.not.solve_rho).and.adiabatic.and.gamma == 1) then
            rho = varray%array_var(ITMP)%array/(tmp_0*jac)
          else
            rho = varray%array_var(IRHO)%array/jac
          endif

        case(ITMP_AUX)

          tmp = tmpe*a_p
cc          where (varray%array_var(IRHO)%array /= 0d0)
cc            tmp = varray%array_var(ITMP)%array
cc     .           /varray%array_var(IRHO)%array
cc          end where

        case(ITMPE_AUX)

          where (varray%array_var(IRHO)%array /= 0d0)
            tmpe = varray%array_var(ITMP)%array
     .            /varray%array_var(IRHO)%array/a_p
          end where

        case(INU)

c       Viscosity (needed for BCs; may depend on T)
 
#if defined(RFX)
          nuu = vis(igr,tmp)
#else
          if (slava) then
            nuu = rho*vis(igr,tmp)
          else
            nuu = vis(igr,tmp)
          endif
#endif
        case(IETA)

c       Resistivity (needed for BCs; may depend on T)
 
          eeta= res(igr,tmpe)

        case(ICHI)

          cchi= chi_perp(igr,tmpe)

        case(IHETA)

c       Find hyperresistivity (electron viscosity)

cc          h_eta=hres(nnx,nny,nnz,igr,igr,igr,bcnv_0,rho_0)
          h_eta=hres(igr,bcnv_0,rho)

        case default
          call pstop('preproc_AuxVar'
     .              ,'Scalar auxiliary variable does not exist')
        end select

      case(VECTOR)

c       Allocate auxiliary variables in local domain

        allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,3))

        v0 = 0d0

        select case(abs(bc_inf(1)))
        case(IVCNV_0)

          vcnv_0(:,:,:,1) = u_0%array_var(IVX)%array
          vcnv_0(:,:,:,2) = u_0%array_var(IVY)%array
          vcnv_0(:,:,:,3) = u_0%array_var(IVZ)%array

          if (.not.nc_eom_v) then
            where (rho_0 /= 0d0)
              vcnv_0(:,:,:,1) = vcnv_0(:,:,:,1)/rho_0
              vcnv_0(:,:,:,2) = vcnv_0(:,:,:,2)/rho_0
              vcnv_0(:,:,:,3) = vcnv_0(:,:,:,3)/rho_0
            end where
          endif

          if (alt__eom()) vcnv_0(:,:,:,2) = vcnv_0(:,:,:,2)/jac

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
     .                + curl(i,j,k,nnx,nny,nnz,igr,igr,igr,acov_0)
              enddo
            enddo
          enddo
#else
          bcnv_0(:,:,:,1) = u_0%array_var(IBX)%array
          bcnv_0(:,:,:,2) = u_0%array_var(IBY)%array
          bcnv_0(:,:,:,3) = u_0%array_var(IBZ)%array
#endif

        case(IBCOV_0)

          bcov_0 = XformVector(igr,bcnv_0,.false.)

        case(IJCOV_0)

          jcov_0 = XformVector(igr,jcnv_0,.false.)

        case(IJCNV_0)

c       Equilibrium current

c diag ****
          if (test .and. equil=='tmcar') then 

            jcnv_0 = 0d0

          elseif (sum(E0) /= 0d0) then
c diag ****
cc          if (sum(E0) /= 0d0) then

            !Enforce imposed E field
            do k=0,nnz+1
              do j=0,nny+1
                do i=0,nnx+1
                  jcnv_0(i,j,k,:) =
     .                 matmul(gmetric%grid(igr)%gsup(i,j,k,:,:),E0)
     .                /eeta_0(i,j,k)
                enddo
              enddo
            enddo

          else

            jcnv_0 = curl(igr,bcov_0)

          endif

        case(IENI_0)

c       Electric field (time=0, initial EOM equ => only ideal+resistive contributions)

#if defined(vec_pot)
          E_ni_0 = E_fld_res(igr,vcnv_0,bcnv_0,acov_0,eeta_0)
#else
          E_ni_0 = E_fld_res(igr,vcnv_0,bcnv_0,jcov_0,eeta_0)
#endif

cc        case(IVCNV_NM)
cc
cc          if (save_tn.and.ion_hall) vcnv_nm = vcnv_n

        case(IVCNV_N)

          if (save_tn.and.ion_hall) vcnv_n = vcnv  !Initializes for BDF2 as well

cc        case(IVECOV_NM)
cc
cc          if (save_tn.and.(.not.ion_hall)) vecov_nm= vecov_n

        case(IVECOV_N)

          if (save_tn.and.(.not.ion_hall)) vecov_n= vecov  !Initializes for BDF2 as well

        case(IJCNV_N)

cc          if (save_tn.and.ion_hall) jcnv_n = jcnv  !Initializes for BDF2 as well

cc        case(IBCNV_N)
cc
cc          if (save_tn) bcnv_n = bcnv !Save previous time step

        case(IVCOV)

          vcov = XformVector(igr,vcnv,.false.)

        case(IVCNV)

          vcnv(:,:,:,1) = varray%array_var(IVX)%array
          vcnv(:,:,:,2) = varray%array_var(IVY)%array
          vcnv(:,:,:,3) = varray%array_var(IVZ)%array

          if (.not.nc_eom_v) then
            where (rho /= 0d0)
              vcnv(:,:,:,1) = vcnv(:,:,:,1)/rho
              vcnv(:,:,:,2) = vcnv(:,:,:,2)/rho
              vcnv(:,:,:,3) = vcnv(:,:,:,3)/rho
            end where
          endif

          if (alt__eom()) vcnv(:,:,:,2) = vcnv(:,:,:,2)/jac

#if defined(vec_pot)

        case(IACOV)

          acov(:,:,:,1) = varray%array_var(IAX)%array
          acov(:,:,:,2) = varray%array_var(IAY)%array
          acov(:,:,:,3) = varray%array_var(IAZ)%array

        case(IACNV)

          acnv = XformVector(igr,acov,.true.)

        case(IBCNV)

c       Magnetic field BC (with fixed current at boundary)

          do k = 0,nnz+1
            do j = 0,nny+1
              do i = 0,nnx+1
                bcnv(i,j,k,:) = B0
     .                 + curl(i,j,k,nnx,nny,nnz,igr,igr,igr,acov)
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

          bcov = XformVector(igr,bcnv,.false.)

        case(IJCOV)

          jcov = XformVector(igr,jcnv,.false.)

        case(IJCNV)

          jcnv = curl(igr,bcov)

        case(IVECNV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          vecnv = form_ve(vcnv,jcnv,rho)

        case(IVECOV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          vecov = XformVector(igr,vecnv,.false.)

        case(IVEFCNV)

c       Fake current

          if (fake_ve) then
            if (ion_hall) then
              !Ion inertia
              if (save_tn) then !Set dt=0 at time t=t_n
                v0 = 0d0
              else
cc                v0 =  cnp*vcnv+cn*vcnv_n+cnm*vcnv_nm
                v0 =  cnp*(vcnv-vcnv_n)  !Works for BDF2
              endif

              !j=jn + dj_fake
              vefcnv = jcnv_n + form_jfake(igr,v0,div_pi,rho)
            else
#if defined(vec_pot)
              !j=curl(curl(A))
              v0 = curlcurl_cov(igr,acov)
              vefcnv = XformVector(igr,v0,.true.)
#else
              call pstop('applyBC'
     .                  ,'No fake Ve available for B-based E-hall')
#endif
            endif

          else

            vefcnv = jcnv

          endif

c old TM Hall MHD performance results
cc          do k=1,nnz
cc            do j=1,nny
cc              do i=1,nnx
cc                vefcnv(i,j,k,:) = di*veclap(i,j,k,nnx,nny,nnz
cc     .                                       ,igr,igr,igr,vcnv
cc     .                                       ,vol=.false.)
cc              enddo
cc            enddo
cc          enddo
c old TM Hall performance results

        case(IDIVPI)

c       Divergence of ion stress tensor

          if (nu > 0d0) then
            div_pi = EOM_divPi(igr,vcnv,nuu)
          else
            div_pi = 0d0
          endif

        case(IDIVPE)

c       Divergence of electron stress tensor

          if (di > 0d0) then
            div_pe = EOM_divPi(igr,vefcnv,h_eta)
          else
            div_pe = 0d0
          endif

        case(IEH)

c       Hall electric field

          if (di > 0d0) then
            E_h  = E_fld_Hall(nnx,nny,nnz,igr,igr,igr)
          else
            E_h  = 0d0
          endif

        case(IENI)

c       Ideal+resistive E field

#if defined(vec_pot)
          E_ni=E_fld_res(igr,vcnv,bcnv,acov,eeta)
#else
          E_ni=E_fld_res(igr,vcnv,bcnv,jcov,eeta)
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
      subroutine postproc_AuxVar(igr,bc_inf,varray,vaux)
c--------------------------------------------------------------------
c     Postprocesses auxiliary variables after apply BC routine.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use var_setup

      implicit none

c Call variables

      integer :: igr,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

c Local variables

      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz
     .          ,ig,jg,kg

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),grad_p(3),etal

      real(8),pointer,dimension(:,:,:) :: jac

c Begin program

c Local grid sizes

      nnx = grid_params%nxv(igr) 
      nny = grid_params%nyv(igr)
      nnz = grid_params%nzv(igr)

      jac => gmetric%grid(igr)%jac

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

        select case(abs(bc_inf(1)))
        case(IRHO_0)

          u_0%array_var(IRHO)%array = rho_0*jac

        case(ITMP_0)

          u_0%array_var(ITMP)%array = rho_0*tmp_0*jac

        case(ITMPE_0)

        case(IETA_0)

c       Resistivity at initial time
 
          eeta_0=res(igr,tmpe_0)

        case(IRHO_AUX)

          if ((.not.solve_rho).and.adiabatic.and.gamma > 1d0) then
            where (tmp_0 /= 0d0)
              rho = rho_0*(tmp/tmp_0)**(1./(gamma-1.))
            end where
          endif

          varray%array_var(IRHO)%array = rho*jac

        case(ITMP_AUX)

          if ((.not.solve_rho).and.adiabatic.and.gamma == 1d0) then
            varray%array_var(ITMP)%array = rho*tmp_0*jac
          else
            varray%array_var(ITMP)%array = rho*tmp*jac
          endif

        case(ITMPE_AUX)

        case(INU)

c       Viscosity (needed for BCs; may depend on T)
 
#if defined(RFX)
          nuu = vis(igr,tmp)
#else
          if (slava) then
            nuu = rho*vis(igr,tmp)
          else
            nuu = vis(igr,tmp)
          endif
#endif
        case(IETA)

c       Resistivity (needed for BCs; may depend on T)
 
          eeta= res(igr,tmpe)

        case(ICHI)

c       Heat diffusivity (may depend on T)

          cchi= chi_perp(igr,tmpe)

        case(IHETA)

c       Find hyperresistivity (electron viscosity)

cc          h_eta=hres(igr,bcnv_0,rho_0)
          h_eta=hres(igr,bcnv_0,rho)

        case default
          call pstop('postproc_AuxVar'
     .              ,'Scalar auxiliary variable does not exist')
        end select

      case(VECTOR)

c       Allocate auxiliary variables in local domain

        select case(abs(bc_inf(1)))
        case(IVCNV_0)

          u_0%array_var(IVX)%array = vcnv_0(:,:,:,1)
          u_0%array_var(IVY)%array = vcnv_0(:,:,:,2)
          u_0%array_var(IVZ)%array = vcnv_0(:,:,:,3)

          if (.not.nc_eom_v) then
            u_0%array_var(IVX)%array = u_0%array_var(IVX)%array*rho_0
            u_0%array_var(IVY)%array = u_0%array_var(IVY)%array*rho_0
            u_0%array_var(IVZ)%array = u_0%array_var(IVZ)%array*rho_0
          endif

          if (alt__eom()) u_0%array_var(IVY)%array =
     .                    u_0%array_var(IVY)%array*jac

#if defined(vec_pot)
        case(IACOV_0)

          u_0%array_var(IAX)%array = acov_0(:,:,:,1)
          u_0%array_var(IAY)%array = acov_0(:,:,:,2)
          u_0%array_var(IAZ)%array = acov_0(:,:,:,3)
#endif
        case(IBCNV_0)

#if !defined(vec_pot)
          u_0%array_var(IBX)%array = bcnv_0(:,:,:,1)
          u_0%array_var(IBY)%array = bcnv_0(:,:,:,2)
          u_0%array_var(IBZ)%array = bcnv_0(:,:,:,3)
#endif
        case(IBCOV_0)

cc          call XformVector_BC(igr,bcnv_0,bcov_0,.false.)

        case(IJCOV_0)

cc          call XformVector_BC(igr,jcnv_0,jcov_0,.false.)

        case(IJCNV_0)

        case(IENI_0)

c       Electric field (time=0, initial EOM equ => only ideal+resistive contributions)

#if defined(vec_pot)
          E_ni_0 = E_fld_res(igr,vcnv_0,bcnv_0,acov_0,eeta_0)
#else
          E_ni_0 = E_fld_res(igr,vcnv_0,bcnv_0,jcov_0,eeta_0)
#endif

cc        case(IVCNV_NM)
cc
cc          if (save_tn.and.ion_hall) vcnv_nm = vcnv_n

        case(IVCNV_N)

          if (save_tn.and.ion_hall) vcnv_n  = vcnv  !Initializes for BDF2 as well

cc        case(IVECOV_NM)
cc
cc          if (save_tn.and.(.not.ion_hall)) vecov_nm= vecov_n

        case(IVECOV_N)

          if (save_tn.and.(.not.ion_hall)) vecov_n= vecov  !Initializes for BDF2 as well

        case(IJCNV_N)

          if (save_tn.and.ion_hall) jcnv_n = jcnv  !Initializes for BDF2 as well

        case(IVCOV)

cc          call XformVector_BC(igr,vcnv,vcov,.false.)

        case(IVCNV)

          varray%array_var(IVX)%array = vcnv(:,:,:,1)
          varray%array_var(IVY)%array = vcnv(:,:,:,2)
          varray%array_var(IVZ)%array = vcnv(:,:,:,3)

          if (.not.nc_eom_v) then
            varray%array_var(IVX)%array =varray%array_var(IVX)%array*rho
            varray%array_var(IVY)%array =varray%array_var(IVY)%array*rho
            varray%array_var(IVZ)%array =varray%array_var(IVZ)%array*rho
          endif

          if (alt__eom()) varray%array_var(IVY)%array =
     .                    varray%array_var(IVY)%array*jac

#if defined(vec_pot)

        case(IACOV)

          varray%array_var(IAX)%array = acov(:,:,:,1)
          varray%array_var(IAY)%array = acov(:,:,:,2)
          varray%array_var(IAZ)%array = acov(:,:,:,3)

        case(IACNV)

cc          call XformVector_BC(igr,acov,acnv,.true.)

        case(IBCNV)

#else

        case(IBCNV)

          varray%array_var(IBX)%array = bcnv(:,:,:,1)
          varray%array_var(IBY)%array = bcnv(:,:,:,2)
          varray%array_var(IBZ)%array = bcnv(:,:,:,3)
#endif

        case(IBCOV)

cc          call XformVector_BC(igr,bcnv,bcov,.false.)

        case(IJCOV)

        case(IJCNV)

        case(IVECNV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          vecnv = form_ve(vcnv,jcnv,rho)

        case(IVECOV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          call XformVector_BC(igr,vecnv,vecov,.false.)

        case(IVEFCNV)

c       Fake electron velocity: vef = v -di*jfake/rho

          vefcnv = form_ve(vcnv,vefcnv,rho)

        case(IDIVPI)

        case(IDIVPE)

        case(IEH)

        case(IENI)

c       Total E field

#if defined(vec_pot)
          E_ni=E_fld_res(igr,vcnv,bcnv,acov,eeta)
#else
          E_ni=E_fld_res(igr,vcnv,bcnv,jcov,eeta)
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
     .                    ,dim,loc,order)
cc     .                    ,dim,loc,1)  !Required for mass conservation

      case (IBX,IBY,IBZ,IJX,IJY,IJZ)

        if (ivar /= dim) then

          call extrapolateBC(array (:,:,:,ivar)
     .                      ,array0(:,:,:,ivar)
     .                      ,dim,loc,order)

        else  !Imposes divergence-free constraint on B, J

          call solen_bc(array,dim,loc,.false.)

        endif

      case (IAX,IAY,IAZ)

        if (ivar /= dim .or. bcond(ibc) == FSYM ) then

          call extrapolateBC(array (:,:,:,ivar)
     .                      ,array0(:,:,:,ivar)
     .                      ,dim,loc,order)

        else  !Imposes divergence-free constraint on A @ physical boundaries

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

      ibc = (1+loc)+2*(dim-1)

      offset = 0
      if (bcSYM(ibc)) offset = 1

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

      ibc = (1+loc)+2*(dim-1)

      offset = 0
      if (bcSYM(ibc)) offset = 1

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

        if ((.not.no_map).and.(.not.bcSYM(ibc))) then

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
c        * dim   -> dimension we are imposing BC on (X,Y,Z) (1,2,3)
c        * loc   -> boundary location (0 -> left, 1->right)
c        * face  -> whether we are approximating normal derivative at
c                   boundary face or at inner boundary cell center.
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

