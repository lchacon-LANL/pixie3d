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

      integer :: nnx,nny,nnz,bcsq(6)
      real(8),pointer,dimension(:,:,:) :: jac

      real(8),pointer,dimension(:,:,:) :: s0

c Begin program

c Local grid sizes

      nnx = gv%gparams%nxv(igr) 
      nny = gv%gparams%nyv(igr)
      nnz = gv%gparams%nzv(igr)

      jac => gv%gparams%gmetric%grid(igr)%jac

c Select variable

      select case(abs(bc_inf(1)))
      case(IRHO)

c     Density BCs

cc        varray%array_var(IRHO)%array = varray%array_var(IRHO)%array/jac
cc
cc        call setBC(gv%gparams,IRHO,nnx,nny,nnz
cc     .            ,varray%array_var(IRHO)%array
cc     .            ,u_0   %array_var(IRHO)%array
cc     .            ,varray%array_var(IRHO)%bconds
cc     .            ,igr,igr,igr,iorder=order)
cc
cc        varray%array_var(IRHO)%array = varray%array_var(IRHO)%array*jac
        
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

      case(IPRS)

c     Temperature/pressure BCs

c$$$cc        if (adiabatic.and.(.not.adiabatic_rho)) then
c$$$        if (di>0d0) then
c$$$          !Impose BCs on perturbation only
c$$$          allocate(s0(0:nnx+1,0:nny+1,0:nnz+1))
c$$$          s0 = 0d0
c$$$
c$$$          varray%array_var(IPRS)%array=(varray%array_var(IPRS)%array
c$$$     .                                 -u_0  %array_var(IPRS)%array)/jac
c$$$
c$$$          bcsq = varray%array_var(IPRS)%bconds
c$$$          where (bcsq == EQU) bcsq = NEU
c$$$
c$$$          call setBC(gv%gparams,IPRS,nnx,nny,nnz
c$$$     .              ,varray%array_var(IPRS)%array
c$$$     .              ,s0,bcsq,igr,igr,igr,iorder=order)
c$$$
c$$$          varray%array_var(IPRS)%array =varray%array_var(IPRS)%array*jac
c$$$     .                                 +u_0   %array_var(IPRS)%array
c$$$          deallocate(s0)
c$$$        else
          varray%array_var(IPRS)%array =varray%array_var(IPRS)%array/jac

          allocate(s0(0:nnx+1,0:nny+1,0:nnz+1))
          call XferBCs(u_0%array_var(IPRS)%array/jac,s0,.true.)

          call setBC(gv%gparams,IPRS,nnx,nny,nnz
     .              ,varray%array_var(IPRS)%array
     .              ,s0
     .              ,varray%array_var(IPRS)%bconds
     .              ,igr,igr,igr,iorder=order)

          varray%array_var(IPRS)%array =varray%array_var(IPRS)%array*jac
          deallocate(s0)
c$$$        endif

      end select

c End

      end subroutine applyDepVarBC

c applyAuxVarBC
c####################################################################
      subroutine applyAuxVarBC(igr,bc_inf,varray,vaux,order)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use rw_bc

      implicit none

c Call variables

      integer :: igr,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

      integer,optional :: order

c Local variables

      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz
     .          ,ig,jg,kg

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),etal

      real(8),allocatable,dimension(:,:,:,:) :: v0,vdum,acnv_n

      real(8),pointer,dimension(:,:,:) :: s0,prs0,prs1,jac

      logical :: pinch_bc
      
c Begin program

#if defined(plotting)
      save_tn = .true.
#else
      save_tn = old_time_eval .or. source_eval .or. init_time_eval
#endif
      
c Local grid sizes

      nnx = gv%gparams%nxv(igr) 
      nny = gv%gparams%nyv(igr)
      nnz = gv%gparams%nzv(igr)

      jac => gv%gparams%gmetric%grid(igr)%jac

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

#if !defined(samrai)
        no_eq_bc = vaux%var_list(abs(bc_inf(1)))%time_0.and.(time/=0d0)
        no_tn_bc = vaux%var_list(abs(bc_inf(1)))%time_n
     $            .and.(.not.save_tn)
        if ((no_eq_bc.or.no_tn_bc).and.(.not.ts_restart)) return
#endif

        allocate(s0(0:nnx+1,0:nny+1,0:nnz+1))

        select case(abs(bc_inf(1)))
        case(IETA)

c$$$          if (.not.spitzer) then
c$$$            call XferBCs(eeta_0,s0,.true.)
c$$$
c$$$            call setBC(gv%gparams,IETA,nnx,nny,nnz,eeta,s0
c$$$     .                ,vaux%var_list(IETA)%bconds
c$$$     .                ,igr,igr,igr,iorder=0)
c$$$          else
            if (fixed_Te) then
              eeta= res(gv%gparams,igr,tmpe_0)
            elseif (lagged_Te) then
              if (solve_prs) then
                s0 = u_n%array_var(IPRS)%array
     .              /u_n%array_var(IRHO)%array/a_p
              else
                s0 = u_n%array_var(ITMPE)%array/jac
              endif
              eeta= res(gv%gparams,igr,s0)
            else
              eeta= res(gv%gparams,igr,tmpe)
            endif
c$$$          endif
          
        case(IDFF)

          dff     = diff_mesh(gv%gparams,igr,tmpe,rho,bcnv)

        case(ICHI)

          cchi    = chi_perp(gv%gparams,igr,tmpe,rho,bcnv)

        case(ICHIPAR)

          cchi_par= chi_para(gv%gparams,igr,tmpe,cchi)

        case(INU)

          nuu     = vis(gv%gparams,igr,tmp,rho)

        case(IETA_0)

c$$$          if (.not.spitzer) then
c$$$            s0 = 0d0
c$$$            call setBC(gv%gparams,IETA,nnx,nny,nnz,eeta,s0
c$$$     .                ,vaux%var_list(IETA_0)%bconds
c$$$     .                ,igr,igr,igr,iorder=0)
c$$$          else
            eeta_0  = res(gv%gparams,igr,tmpe_0)
c$$$          endif
          
        case(IHETA)

cc          h_eta=hres(gv%gparams,igr,bcnv_0,rho_0)
          h_eta   = hres(gv%gparams,igr,bcnv_0,rho)

        case(IRHO_0)

c       Density BCs
          
          call setBC(gv%gparams,IRHO,nnx,nny,nnz,rho_0
     .              ,u_0%array_var(IRHO)%array
     .              ,vaux%var_list(IRHO_0)%bconds
     .              ,igr,igr,igr,iorder=order)

        case(ITMP_0)

c       Temperature BCs
          
c$$$          call setBC(gv%gparams,IPRS,nnx,nny,nnz,tmp_0
c$$$     .              ,u_0%array_var(IPRS)%array
c$$$     .              ,vaux%var_list(ITMP_0)%bconds
c$$$     .              ,igr,igr,igr,iorder=order)

        case(ITMPE_0)

c       Temperature BCs

          call setBC(gv%gparams,ITMPE,nnx,nny,nnz,tmpe_0
     .              ,u_0%array_var(ITMPE)%array
     .              ,vaux%var_list(ITMPE_0)%bconds
     .              ,igr,igr,igr,iorder=order)

        case(IRHO_AUX)

c       Density BCs

c$$$          if (solve_rho) then
          call XferBCs(rho_0,s0,.true.)

          call setBC(gv%gparams,IRHO,nnx,nny,nnz,rho,s0
     .              ,vaux%var_list(IRHO_AUX)%bconds
     .              ,igr,igr,igr,iorder=order)
c$$$          endif

        case(ITMP_AUX)

c$$$c       Temperature BCs
c$$$
c$$$          if (adiabatic.and.(.not.adiabatic_rho)) then
c$$$          !Impose BCs on perturbation only
c$$$            s0 = 0d0
c$$$            tmp = tmp - tmp_0
c$$$            bcsq = vaux%var_list(ITMP_AUX)%bconds
c$$$            where (bcsq == EQU) bcsq = NEU
c$$$            call setBC(gv%gparams,IPRS,nnx,nny,nnz,tmp,s0
c$$$     .                ,bcsq,igr,igr,igr,iorder=order)
c$$$            tmp = tmp + tmp_0
c$$$          else
c$$$            call XferBCs(tmp_0,s0,.true.)
c$$$            call setBC(gv%gparams,IPRS,nnx,nny,nnz,tmp,s0
c$$$     .                ,vaux%var_list(ITMP_AUX)%bconds
c$$$     .                ,igr,igr,igr,iorder=order)
c$$$          endif
c$$$          
c$$$          tmp = tmpe*a_p

        case(ITMPE_AUX)

c       Temperature BCs

c$$$          if (adiabatic.and.(.not.adiabatic_rho)) then
c$$$          !Impose BCs on perturbation only
c$$$            s0 = 0d0
c$$$            tmpe = tmpe - tmpe_0
c$$$            bcsq = vaux%var_list(ITMPE_AUX)%bconds
c$$$            where (bcsq == EQU) bcsq = NEU
c$$$            call setBC(gv%gparams,IPRS,nnx,nny,nnz,tmpe,s0
c$$$     .                ,bcsq,igr,igr,igr,iorder=order)
c$$$            tmpe = tmpe + tmpe_0
c$$$          else
c$$$            call XferBCs(tmpe_0,s0,.true.)
c$$$            call setBC(gv%gparams,IPRS,nnx,nny,nnz,tmpe,s0
c$$$     .                ,vaux%var_list(ITMPE_AUX)%bconds
c$$$     .                ,igr,igr,igr,iorder=order)
c$$$          endif

        case default
          call pstop('applyAuxVarBC'
     .              ,'Scalar auxiliary variable does not exist')
        end select

        deallocate(s0)

      case(VECTOR)

#if !defined(samrai)
        no_eq_bc = vaux%vec_list(abs(bc_inf(1)))%time_0.and.(time/=0d0)
        no_tn_bc = vaux%vec_list(abs(bc_inf(1)))%time_n
     $            .and.(.not.save_tn)
        if ((no_eq_bc.or.no_tn_bc).and.(.not.ts_restart)) return
#endif

c       Allocate auxiliary variables in local domain

        allocate(prs0(0:nnx+1,0:nny+1,0:nnz+1)
     .          ,prs1(0:nnx+1,0:nny+1,0:nnz+1))

        allocate(v0  (0:nnx+1,0:nny+1,0:nnz+1,3)
     .          ,vdum(0:nnx+1,0:nny+1,0:nnz+1,3))

        v0   = 0d0
        vdum = 0d0

        pinch_bc = pinch_flow.or.res_wall_bc

        select case(abs(bc_inf(1)))
#if defined(vec_pot)
        case(IACOV_0)

          if (fill_A0_bc) then
            vdum = XformVector(gv%gparams,igr,acov_0,.true.)
            call setBC(gv%gparams,IAX,3,nnx,nny,nnz,vdum,acov_0,v0
     .              ,vaux%vec_list(IACOV_0)%bconds,igr,igr,igr
     .              ,is_cnv=vaux%vec_list(IACOV_0)%cnv
     .              ,iorder=2)
          else
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !Provide these in equilibrium routine!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          endif
          
        case(IBCNV_0)

          call setBC(gv%gparams,IBX,3,nnx,nny,nnz,bcnv_0,bcov_0,v0
     .              ,vaux%vec_list(IBCNV_0)%bconds,igr,igr,igr
     .              ,is_cnv=vaux%vec_list(IBCNV_0)%cnv
     .              ,iorder=1)

        case(IBCOV_0)   !Done with bcnv_0
#else
        case(IBCNV_0)

c$$$          if (pinch_bc) then    !Jt = alpha*Bt in cnv representation
c$$$            call setBC(gv%gparams,IBX,3,nnx,nny,nnz,bcnv_0,bcov_0,v0
c$$$     .                ,vaux%vec_list(IBCNV_0)%bconds,igr,igr,igr
c$$$     .                ,is_cnv=vaux%vec_list(IBCNV_0)%cnv
c$$$     .                ,iorder=1)
c$$$          else
c$$$          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c$$$          !Provide these in equilibrium routine!
c$$$          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c$$$          endif

        case(IBCOV_0)

          call XformVector_BC(gv%gparams,igr,bcnv_0,bcov_0,.false.)
#endif

        case(IVCNV_0)

c       Velocity BC

          if (pinch_bc) then
            prs0 = rho_0*tmp_0
            E_bc = 0d0
            call ExB_flow_bc(gv%gparams,igr,bcnv_0,prs0,eeta_0,E_bc,v0)
          else
            v0 = 0d0
          endif

          !Fill ghost nodes
          vdum = XformVector(gv%gparams,igr,vcnv_0,.false.)

          call setBC(gv%gparams,IVX,3,nnx,nny,nnz,vcnv_0,vdum,v0
     .              ,vaux%vec_list(IVCNV_0)%bconds
     .              ,igr,igr,igr,iorder=order)

        case(IJCNV_0)

c       Current BC

          v0 = 0d0

          if (pinch_bc) then  !Jt = alpha*Bt in cnv representation
            prs0 = rho_0*tmp_0
            E_bc = 0d0
            call J_pinch_bc(gv%gparams,igr,bcnv_0,prs0,eeta_0,E_bc,v0)

            call setBC(gv%gparams,IJX,3,nnx,nny,nnz,jcnv_0,jcov_0,v0
     .                ,vaux%vec_list(IJCNV)%bconds,igr,igr,igr
     .                ,iorder=order)

          else
            call setBC(gv%gparams,IJX,3,nnx,nny,nnz,jcnv_0,jcov_0,v0
     .                ,vaux%vec_list(IJCNV_0)%bconds,igr,igr,igr
     .                ,is_cnv=vaux%vec_list(IJCNV_0)%cnv
     .                ,iorder=1)
          endif
          
        case(IJCOV_0)  !Done when computing jcov_0

cc          call XformVector_BC(gv%gparams,igr,jcnv_0,jcov_0,.false.)

        case(IENI_0)

c       Electric field (time=0, initial EOM equ => only ideal+resistive contributions)
          
          if (subtract_E0) then
c$$$#if defined(vec_pot)
c$$$          E_ni_0=E_fld_res_A(gv%gparams,igr,vcnv_0,acov_0,jcov_0,eeta_0)
c$$$#else
            E_ni_0=E_fld_res_B(gv%gparams,igr,vcnv_0,bcnv_0,jcov_0
     .                        ,eeta_0)
c$$$#endif
          endif

        case(IVECNV_0)

c       Electron velocity BC: ve=v-di*j/rho with j=curl(B)

          vecnv_0 = form_ve(vcnv_0,jcnv_0,rho_0)

        case(IVECOV_N)

cc          if (electron_hall) call XferBCs(vecov,vecov_n,.false.)

        case(IVCNV_N)

#if !defined(plotting)
c$$$          call XferBCs(vcnv,vcnv_n,.false.)
#endif

        case(IVCNV)

c       Velocity BC

          if (pinch_bc) then
            prs1 = rho*tmp
            call ExB_flow_bc(gv%gparams,igr,bcnv,prs1,eeta,E_bc,v0)
          else
            call XferBCs(vcnv_0,v0,.true.)
          endif

          call setBC(gv%gparams,IVX,3,nnx,nny,nnz,vcnv,vcov,v0
     .              ,vaux%vec_list(IVCNV)%bconds,igr,igr,igr
     .              ,iorder=order)
 
        case(IVCOV)   !Done when computing vcnv

cc          call XformVector_BC(gv%gparams,igr,vcnv,vcov,.false.)

#if defined(vec_pot)

        case(IACOV)

c       Vector Pot. BC

          if (regularize_A) then
            allocate(acnv_n(0:nnx+1,0:nny+1,0:nnz+1,3))

            !Subtract previous time step (to regularize delta only)
c$$$            acov_n(:,:,:,1) = u_n%array_var(IAX)%array
c$$$            acov_n(:,:,:,2) = u_n%array_var(IAY)%array
c$$$            acov_n(:,:,:,3) = u_n%array_var(IAZ)%array
            acnv_n = XformVector(gv%gparams,igr,acov_n,.true.)

            acov = acov-acov_n
            acnv = acnv-acnv_n

c$$$            acnv_n = XformVector(gv%gparams,igr,acov_0,.true.)
c$$$
c$$$            acov = acov-acov_0
c$$$            acnv = acnv-acnv_n

            !Impose BC
            v0 = 0d0

            call setBC(gv%gparams,IAX,3,nnx,nny,nnz,acnv,acov,v0
     .                ,vaux%vec_list(IACOV)%bconds,igr,igr,igr
     .                ,iorder=order,is_cnv=vaux%vec_list(IACOV)%cnv
     .                ,sp_reg=.true.)

            !Add previous time step
            acov = acov + acov_n
cc            acov = acov + acov_0
            acnv = acnv + acnv_n

            deallocate(acnv_n)
          else
            if (res_wall_bc.and.(.not.source_eval).and.(time/=0d0)) then
#if defined(plotting)
              if (save_tn) call rw_Ebc(gv%gparams,igr,bcnv_n,E_bc)
              call XferBCs(acov,v0,.true.)  !Reuse existing BC for plotting
              call setBC(gv%gparams,IAX,3,nnx,nny,nnz,acnv,acov,v0
     .                  ,vaux%vec_list(IACOV)%bconds,igr,igr,igr
     .                  ,iorder=order,is_cnv=vaux%vec_list(IACOV)%cnv)
#else
              if (save_tn) then
                call rw_Ebc(gv%gparams,igr,bcnv_n,E_bc)
              elseif (rw_implicit_bc) then
                call find_B(gv%gparams,igr,B0,acov,bcov,.true.
     $                     ,bcond=vaux%vec_list(IACOV)%bconds)
                call rw_Ebc(gv%gparams,igr,bcov,E_bc)
              endif

              call XferBCs(acov_n,v0,.true.)
              v0(nnx+1,:,:,:) = v0(nnx+1,:,:,:)-dt*E_bc(nnx+1,:,:,:)

              if (.not.init_time_eval) then
                call setBC(gv%gparams,IAX,3,nnx,nny,nnz,acnv,acov,v0
     .                    ,vaux%vec_list(IACOV)%bconds,igr,igr,igr
     .                    ,iorder=order,is_cnv=vaux%vec_list(IACOV)%cnv)
              endif
#endif
            else
              call XferBCs(acov_0,v0,.true.)
              call setBC(gv%gparams,IAX,3,nnx,nny,nnz,acnv,acov,v0
     .                  ,vaux%vec_list(IACOV)%bconds,igr,igr,igr
     .                  ,iorder=order,is_cnv=vaux%vec_list(IACOV)%cnv)
            endif
          endif

        case(IACNV)  !Done when computing acov

cc          call XformVector_BC(gv%gparams,igr,acov,acnv,.true.)

        case(IACOV_N)
#endif
        case(IBCNV)

c       Magnetic field BC

          if (pinch_bc) then  !Jt = alpha*Bt in cnv representation
#if !defined(vec_pot)
            if (res_wall_bc.and.(.not.source_eval).and.(time/=0d0)) then
              if (save_tn) then
                call rw_Ebc(gv%gparams,igr,bcnv_n,E_bc)
              elseif (rw_implicit_bc) then
                call rw_Ebc(gv%gparams,igr,bcov  ,E_bc,extrp_b=.true.)
              endif
            endif
c$$$#else
c$$$            if (res_wall_bc.and.(.not.source_eval).and.(time/=0d0)) then
c$$$              if ((.not.save_tn).and.rw_implicit_bc) then
c$$$                call rw_Ebc(gv%gparams,igr,bcov  ,E_bc,implicit=.true.)
c$$$              endif
c$$$            endif
#endif
            v0 = 0d0 ; prs1 = rho*tmp
            call J_pinch_bc(gv%gparams,igr,bcnv,prs1,eeta,E_bc,v0)
          else
            call XferBCs(jcnv_0,v0,.true.)
          endif

          call setBC(gv%gparams,IBX,3,nnx,nny,nnz,bcnv,bcov,v0
     .              ,vaux%vec_list(IBCNV)%bconds,igr,igr,igr
     .              ,iorder=order,is_cnv=vaux%vec_list(IBCNV)%cnv)

        case(IBCOV) !Done when computing bcnv

cc          call XformVector_BC(gv%gparams,igr,bcnv,bcov,.false.)

        case(IBCNV_N)

        case(IJCNV)

c       Current BC

          if (pinch_bc) then  !Jt = alpha*Bt in cnv representation
            v0 = 0d0 ; prs1 = rho*tmp
            call J_pinch_bc(gv%gparams,igr,bcnv,prs1,eeta,E_bc,v0)
          else
            call XferBCs(jcnv_0,v0,.true.)
          endif

          call setBC(gv%gparams,IJX,3,nnx,nny,nnz,jcnv,jcov,v0
     .              ,vaux%vec_list(IJCNV)%bconds,igr,igr,igr
     .              ,iorder=order)

        case(IJCOV) !Done when computing jcnv

cc          call XformVector_BC(gv%gparams,igr,jcnv,jcov,.false.)

        case(IVECNV)

c       Electron velocity BC: ve=v-di*j/rho with j=curl(B)

          vecnv = form_ve(vcnv,jcnv,rho)

        case(IVECOV)

cc          call XformVector_BC(gv%gparams,igr,vecnv,vecov,.false.)

        case(IDIVPI)

#if !defined(vec_pot) || defined(plotting)
          if (ion_hall) then
            bcnd = vaux%vec_list(IDIVPI)%bconds
            v0 = 0d0
            vdum = XformVector(gv%gparams,igr,div_pi,.false.)
            call setBC(gv%gparams,IVX,3,nnx,nny,nnz,div_pi,vdum,v0,bcnd
     .                ,igr,igr,igr,iorder=order
     .                ,is_cnv=vaux%vec_list(IDIVPI)%cnv)
          endif
#endif

        case(IDIVPE)

#if !defined(vec_pot) || defined(plotting)
          if (electron_hall) then
            bcnd = vaux%vec_list(IDIVPE)%bconds
            v0 = 0d0
            vdum = XformVector(gv%gparams,igr,div_pe,.false.)
            call setBC(gv%gparams,IVX,3,nnx,nny,nnz,div_pe,vdum,v0,bcnd
     .                ,igr,igr,igr,iorder=order
     .                ,is_cnv=vaux%vec_list(IDIVPE)%cnv)
          endif
#endif

        case(IEH)

c       Hall electric field BC

          if (di > 0d0) then
            prs1 = rho*tmpe*temp_ratio  !Ion pressure
            call XferBCs(grad(gv%gparams,igr,prs1),v0,.true.)

            bcnd = vaux%vec_list(IEH)%bconds
            vdum = XformVector(gv%gparams,igr,E_h,.true.)
            call setBC(gv%gparams,IVX,3,nnx,nny,nnz,vdum,E_h,v0,bcnd
     .                ,igr,igr,igr,iorder=order
     .                ,is_cnv=vaux%vec_list(IEH)%cnv)
          endif

        case(IENI)

c       Total electric field BC

c$$$          bcnd = vaux%vec_list(IENI)%bconds
c$$$          vdum = XformVector(gv%gparams,igr,E_ni,.true.)
c$$$          call setBC(gv%gparams,IVX,3,nnx,nny,nnz,vdum,E_ni,v0,bcnd
c$$$     .              ,igr,igr,igr,iorder=order
c$$$     .              ,is_cnv=vaux%vec_list(IENI)%cnv)
#if defined(vec_pot)
          E_ni=E_fld_res_A(gv%gparams,igr,vcnv,acov,jcov,eeta)
#else
          E_ni=E_fld_res_B(gv%gparams,igr,vcnv,bcnv,jcov,eeta)
#endif

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

      implicit none

c Call variables

      integer :: igr,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

      integer,optional :: order

c Local variables

      integer :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz
     .          ,ig,jg,kg

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),etal
     $          ,lref

      real(8),allocatable,dimension(:,:,:,:) :: v0

      real(8),pointer,dimension(:,:,:) :: jac,s0

c Begin program

#if defined(plotting)
      save_tn = .true.
#else
      save_tn = old_time_eval .or. source_eval .or. init_time_eval
#endif
      
c Local grid sizes

      nnx = gv%gparams%nxv(igr) 
      nny = gv%gparams%nyv(igr)
      nnz = gv%gparams%nzv(igr)

      jac => gv%gparams%gmetric%grid(igr)%jac

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

#if !defined(samrai)
        no_eq_bc = vaux%var_list(abs(bc_inf(1)))%time_0.and.(time/=0d0)
        no_tn_bc = vaux%var_list(abs(bc_inf(1)))%time_n
     $            .and.(.not.save_tn)
        if ((no_eq_bc.or.no_tn_bc).and.(.not.ts_restart)) return
#endif

        select case(abs(bc_inf(1)))
        case(IRHO_0)

          rho_0 = u_0%array_var(IRHO)%array/jac

        case(ITMP_0)

          tmp_0 = tmpe_0*a_p

        case(ITMPE_0)

          if (solve_prs) then
            where (u_0%array_var(IRHO)%array /= 0d0)
              tmpe_0 = u_0%array_var(IPRS)%array
     .               /(u_0%array_var(IRHO)%array*a_p)
            elsewhere
              tmpe_0 = 0d0
            end where
          else
            where (u_0%array_var(IRHO)%array /= 0d0)
              tmpe_0 = u_0%array_var(ITMPE)%array/jac
            elsewhere
              tmpe_0 = 0d0
            end where
          endif
        
          !Reference temperature for transport parameters
          if (spitzer.or.chi_bgnskii) then
            te0_ref = pmax(maxval(tmpe_0))
          endif

        case(IETA_0)

c       Resistivity at initial time
 
          eeta_0=res(gv%gparams,igr,tmpe_0)

        case(IRHO_AUX)

          if (adiabatic_rho) then
            if (solve_prs) then
              rho = varray%array_var(IPRS)%array
     .             /u_0   %array_var(IPRS)%array
              rho = rho_0*(abs(rho))**(1d0/gamma)
            else
              rho = varray%array_var(ITMPE)%array
     .             /u_0   %array_var(ITMPE)%array
              rho = rho_0*(abs(rho))**(1d0/(gamma-1d0))
            endif
          else
            rho = varray%array_var(IRHO)%array/jac
          endif

        case(ITMP_AUX)  !Total temperature

          tmp = tmpe*a_p

        case(ITMPE_AUX) !Electron temperature

          if (gamma == 1d0) then
            tmpe = tmpe_0  !Isothermal
c$$$          elseif((adiabatic.and.solve_rho).and.(.not.adiabatic_rho))then
c$$$            tmpe = varray%array_var(IRHO)%array
c$$$     .            /u_0   %array_var(IRHO)%array
c$$$            tmpe = tmpe_0*(abs(tmpe))**(gamma-1d0)
          else
            if (solve_prs) then
              where (rho /= 0d0)
                tmpe = varray%array_var(IPRS )%array/(rho*jac*a_p)
              end where
            else
              tmpe = varray%array_var(ITMPE)%array/jac
            endif
          endif

        case(INU)

c       Viscosity (needed for BCs; may depend on T)
 
          nuu = vis(gv%gparams,igr,tmp,rho)

        case(IETA)

c       Resistivity (needed for BCs; may depend on T)
 
          if (fixed_Te) then
            eeta= res(gv%gparams,igr,tmpe_0)
          elseif (lagged_Te) then
            allocate(s0(0:nnx+1,0:nny+1,0:nnz+1))
            if (solve_prs) then
              s0 = u_n%array_var(IPRS)%array
     .           /(u_n%array_var(IRHO)%array*a_p)
            else
              s0 = u_n%array_var(ITMPE)%array/jac
            endif
            eeta= res(gv%gparams,igr,s0)
            deallocate(s0)
          else
            eeta= res(gv%gparams,igr,tmpe)
          endif

        case(IDFF)

          dff     = diff_mesh(gv%gparams,igr,tmpe,rho,bcnv)

        case(ICHI)

          cchi    = chi_perp(gv%gparams,igr,tmpe,rho,bcnv)

        case(ICHIPAR)

          cchi_par= chi_para(gv%gparams,igr,tmpe,cchi)

        case(IHETA)

c       Find hyperresistivity (electron viscosity)

cc          h_eta=hres(gv%gparams,igr,bcnv_0,rho_0)
          h_eta=hres(gv%gparams,igr,bcnv_0,rho)

        case default
          call pstop('preproc_AuxVar'
     .              ,'Scalar auxiliary variable does not exist')
        end select

      case(VECTOR)

#if !defined(samrai)
        no_eq_bc = vaux%vec_list(abs(bc_inf(1)))%time_0.and.(time/=0d0)
        no_tn_bc = vaux%vec_list(abs(bc_inf(1)))%time_n
     $            .and.(.not.save_tn)
        if ((no_eq_bc.or.no_tn_bc).and.(.not.ts_restart)) return
#endif

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
          call find_B(gv%gparams,igr,B0,acov_0,bcnv_0,.false.)
#else
          bcnv_0(:,:,:,1) = u_0%array_var(IBX)%array
          bcnv_0(:,:,:,2) = u_0%array_var(IBY)%array
          bcnv_0(:,:,:,3) = u_0%array_var(IBZ)%array
#endif

        case(IBCOV_0)

          bcov_0 = XformVector(gv%gparams,igr,bcnv_0,.false.)

        case(IJCOV_0)

          jcov_0 = XformVector(gv%gparams,igr,jcnv_0,.false.)

        case(IJCNV_0)

c       Equilibrium current

c diag ****
          if (test .and. equil=='tmcar') then 

            jcnv_0 = 0d0

c diag ****
c$$$          elseif (equil=='efit') then
c$$$
c$$$            jcnv_0 = 0d0

          else

            jcnv_0 = curl(gv%gparams,igr,bcov_0)

          endif

        case(IENI_0)

c       Electric field (time=0, initial EOM equ => only ideal+resistive contributions)

          if (subtract_E0) then
c$$$#if defined(vec_pot)
c$$$          E_ni_0=E_fld_res_A(gv%gparams,igr,vcnv_0,acov_0,jcov_0,eeta_0)
c$$$#else
            E_ni_0=E_fld_res_B(gv%gparams,igr,vcnv_0,bcnv_0,jcov_0
     .                        ,eeta_0)
c$$$#endif
          endif
          
        case(IVECNV_0)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          vecnv_0 = form_ve(vcnv_0,jcnv_0,rho_0)

        case(IVCNV_N)

#if defined(plotting)
          if (.not.nc_eom_v) then
              vcnv_n(:,:,:,1) = u_n%array_var(IVX )%array
     .                         /u_n%array_var(IRHO)%array*jac
              vcnv_n(:,:,:,2) = u_n%array_var(IVY )%array
     .                         /u_n%array_var(IRHO)%array*jac
              vcnv_n(:,:,:,3) = u_n%array_var(IVZ )%array
     .                         /u_n%array_var(IRHO)%array*jac
          else
            vcnv_n(:,:,:,1) = u_n%array_var(IVX)%array
            vcnv_n(:,:,:,2) = u_n%array_var(IVY)%array
            vcnv_n(:,:,:,3) = u_n%array_var(IVZ)%array
          endif

          if (alt__eom()) vcnv_n(:,:,:,2) = vcnv_n(:,:,:,2)/jac
#else
          if (.not.nc_eom_v) then
              if (cnm == 0d0) then
                vcnv_n(:,:,:,1) =(u_n %array_var(IVX )%array
     .                           /u_n %array_var(IRHO)%array)*jac
                vcnv_n(:,:,:,2) =(u_n %array_var(IVY )%array
     .                           /u_n %array_var(IRHO)%array)*jac
                vcnv_n(:,:,:,3) =(u_n %array_var(IVZ )%array
     .                           /u_n %array_var(IRHO)%array)*jac
              else
                vcnv_n(:,:,:,1) =(-cn/cnp*u_n %array_var(IVX )%array
     .                                   /u_n %array_var(IRHO)%array
     .                           -cnm/cnp*u_nm%array_var(IVX )%array
     .                                   /u_nm%array_var(IRHO)%array)
     .                           *jac
                vcnv_n(:,:,:,2) =(-cn/cnp*u_n %array_var(IVY )%array
     .                                   /u_n %array_var(IRHO)%array
     .                           -cnm/cnp*u_nm%array_var(IVY )%array
     .                                   /u_nm%array_var(IRHO)%array)
     .                           *jac
                vcnv_n(:,:,:,3) =(-cn/cnp*u_n %array_var(IVZ )%array
     .                                   /u_n %array_var(IRHO)%array
     .                           -cnm/cnp*u_nm%array_var(IVZ )%array
     .                                   /u_nm%array_var(IRHO)%array)
     .                           *jac
             endif
          else
            vcnv_n(:,:,:,1) = -cn/cnp*u_n %array_var(IVX)%array
     .                       -cnm/cnp*u_nm%array_var(IVX)%array
            vcnv_n(:,:,:,2) = -cn/cnp*u_n %array_var(IVY)%array
     .                       -cnm/cnp*u_nm%array_var(IVY)%array
            vcnv_n(:,:,:,3) = -cn/cnp*u_n %array_var(IVZ)%array
     .                       -cnm/cnp*u_nm%array_var(IVZ)%array
          endif

          if (alt__eom()) vcnv_n(:,:,:,2) = vcnv_n(:,:,:,2)/jac

c$$$          vcnv_n = vcnv  !Initializes for BDF2 as well
#endif

        case(IVECOV_N)

          if (electron_hall) vecov_n = vecov     !Does NOT work for BDF2
 
        case(IBCNV_N)

#if defined(plotting)
          if (itime==0) then
#else
          if (source_eval) then
#endif
            bcnv_n = bcnv_0     !Set it to initial cond.
          elseif (save_tn) then
#if defined(vec_pot)
            acov_n(:,:,:,1) = varray%array_var(IAX)%array
            acov_n(:,:,:,2) = varray%array_var(IAY)%array
            acov_n(:,:,:,3) = varray%array_var(IAZ)%array
            call find_B(gv%gparams,igr,B0,acov_n,bcnv_n,.false.)
#else
            bcnv_n(:,:,:,1) = varray%array_var(IBX)%array
            bcnv_n(:,:,:,2) = varray%array_var(IBY)%array
            bcnv_n(:,:,:,3) = varray%array_var(IBZ)%array
#endif
          endif
          
        case(IVCOV)

          vcov = XformVector(gv%gparams,igr,vcnv,.false.)

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
        case(IACOV_N)

#if defined(plotting)
          if (itime==0) then
            acov_n = acov_0     !Set it to initial cond.
          else
            acov_n = acov       !Save previous time step
          endif
#else
          if (source_eval) then
            acov_n = acov_0     !Set it to initial cond.
          else  !This determines bcnv_n also
            acov_n(:,:,:,1) = u_n%array_var(IAX)%array
            acov_n(:,:,:,2) = u_n%array_var(IAY)%array
            acov_n(:,:,:,3) = u_n%array_var(IAZ)%array
          endif
#endif

        case(IACOV)

          acov(:,:,:,1) = varray%array_var(IAX)%array
          acov(:,:,:,2) = varray%array_var(IAY)%array
          acov(:,:,:,3) = varray%array_var(IAZ)%array

        case(IACNV)

          acnv = XformVector(gv%gparams,igr,acov,.true.)

        case(IBCNV)

c       Magnetic field BC (with fixed current at boundary)

          call find_B(gv%gparams,igr,B0,acov,bcnv,.false.)
#else
        case(IBCNV)

          bcnv(:,:,:,1) = varray%array_var(IBX)%array
          bcnv(:,:,:,2) = varray%array_var(IBY)%array
          bcnv(:,:,:,3) = varray%array_var(IBZ)%array
#endif
        case(IBCOV)

          bcov = XformVector(gv%gparams,igr,bcnv,.false.)

        case(IJCOV)

          jcov = XformVector(gv%gparams,igr,jcnv,.false.)

        case(IJCNV)

          jcnv = curl(gv%gparams,igr,bcov)

        case(IVECNV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          vecnv = form_ve(vcnv,jcnv,rho)

        case(IVECOV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          vecov = XformVector(gv%gparams,igr,vecnv,.false.)

        case(IDIVPI)

c       Divergence of ion stress tensor

          if (nu > 0d0) then
            div_pi = EOM_divPi(gv%gparams,igr,vcnv,nuu,sym_st)
          else
            div_pi = 0d0
          endif

        case(IDIVPE)

c       Divergence of electron stress tensor

          if (di > 0d0) then
            if (subtract_E0) then
              div_pe=EOM_divPi(gv%gparams,igr,vecnv-vecnv_0,h_eta
     .                        ,.false.)
            else
              div_pe=EOM_divPi(gv%gparams,igr,vecnv,h_eta,.false.)
            endif
          else
            div_pe=0d0
          endif

        case(IEH)

c       Hall electric field

          if (electron_hall) then
            E_h  = E_fld_eHall(gv%gparams,igr,rho,tmpe,bcnv,jcnv
     .                        ,div_pe,vecov,vecov_n)
          elseif(ion_hall) then
            E_h  = E_fld_iHall(gv%gparams,igr,rho,tmpe
     .                        ,div_pi,vcnv ,vcnv_n)
          endif

        case(IENI)

c       Ideal+resistive E field

#if defined(vec_pot)
          E_ni=E_fld_res_A(gv%gparams,igr,vcnv,acov,jcov,eeta)
#else
          E_ni=E_fld_res_B(gv%gparams,igr,vcnv,bcnv,jcov,eeta)
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

      real(8) :: bbcov(3),bbcnv(3),b2,J0(3),E00(3),etal

      real(8),pointer,dimension(:,:,:) :: jac,s0

c Begin program

#if defined(plotting)
      save_tn = .true.
#else
      save_tn = old_time_eval .or. source_eval .or. init_time_eval
#endif
      
c Local grid sizes

      nnx = gv%gparams%nxv(igr) 
      nny = gv%gparams%nyv(igr)
      nnz = gv%gparams%nzv(igr)

      jac => gv%gparams%gmetric%grid(igr)%jac

c Select variable

      select case(bc_inf(2)) !Whether we are dealing with scalars or vectors
      case(SCALAR)

#if !defined(samrai)
        no_eq_bc = vaux%var_list(abs(bc_inf(1)))%time_0.and.(time/=0d0)
        no_tn_bc = vaux%var_list(abs(bc_inf(1)))%time_n
     $            .and.(.not.save_tn)
        if ((no_eq_bc.or.no_tn_bc).and.(.not.ts_restart)) return
#endif

        select case(abs(bc_inf(1)))
        case(IRHO_0)

          u_0%array_var(IRHO)%array = rho_0*jac

        case(ITMP_0)

          tmp_0 = tmpe_0*a_p

c$$$          u_0%array_var(IPRS)%array = rho_0*tmp_0*jac

        case(ITMPE_0)

c$$$          if (solve_prs) then
c$$$            tmpe_0 = u_0%array_var(IPRS)%array
c$$$     .             /(u_0%array_var(IRHO)%array*a_p)
c$$$          else
c$$$            tmpe_0 = u_0%array_var(ITMPE)%array
c$$$          endif

        case(IETA_0)

c       Resistivity at initial time
 
c$$$          eeta_0=res(gv%gparams,igr,tmpe_0)

        case(IRHO_AUX)

          varray%array_var(IRHO)%array = rho*jac

        case(ITMP_AUX)  !Recover total temperature

c$$$          if (.not.adiabatic_rho) then
c$$$            varray%array_var(IPRS)%array = rho*tmp*jac
c$$$          endif
          
          tmp = tmpe*a_p

        case(ITMPE_AUX)

          if (gamma == 1d0) then
            tmpe = tmpe_0  !Isothermal
c$$$          elseif((adiabatic.and.solve_rho).and.(.not.adiabatic_rho))then
c$$$            tmpe = varray%array_var(IRHO)%array
c$$$     .            /u_0   %array_var(IRHO)%array
c$$$            tmpe = tmpe_0*(abs(tmpe))**(gamma-1d0)
          else
            if (solve_prs) then
              where (rho /= 0d0)
                tmpe = varray%array_var(IPRS)%array/(rho*jac*a_p)
              end where
            else
              tmpe = varray%array_var(ITMPE)%array/jac
            endif
          endif

        case(INU)

c       Viscosity (needed for BCs; may depend on T)
 
          nuu = vis(gv%gparams,igr,tmp,rho)

        case(IETA)

c       Resistivity (needed for BCs; may depend on T)

c$$$          if (fixed_Te) then
c$$$            eeta= res(gv%gparams,igr,tmpe_0)
c$$$          elseif (lagged_Te) then
c$$$            allocate(s0(1:nnx+1,1:nny+1,1:nnz+1))
c$$$            if (solve_prs) then
c$$$              s0 = u_n%array_var(IPRS)%array
c$$$     .           /(u_n%array_var(IRHO)%array*a_p)
c$$$            else
c$$$              s0 = u_n%array_var(ITMPE)%array/jac
c$$$            endif
c$$$            eeta= res(gv%gparams,igr,s0)
c$$$            deallocate(s0)
c$$$          else
c$$$            eeta= res(gv%gparams,igr,tmpe)
c$$$          endif
          
        case(IDFF)

c       Particle diffusivity

          dff = diff_mesh(gv%gparams,igr,tmpe,rho,bcnv)

        case(ICHI)

c       Perpendicular heat diffusivity

          cchi= chi_perp(gv%gparams,igr,tmpe,rho,bcnv)

        case(ICHIPAR)

c       Parallel heat diffusivity

          cchi_par= chi_para(gv%gparams,igr,tmpe,cchi)

        case(IHETA)

c       Find hyperresistivity (electron viscosity)

cc          h_eta=hres(gv%gparams,igr,bcnv_0,rho_0)
          h_eta=hres(gv%gparams,igr,bcnv_0,rho)

        case default
          call pstop('postproc_AuxVar'
     .              ,'Scalar auxiliary variable does not exist')
        end select

      case(VECTOR)

#if !defined(samrai)
        no_eq_bc = vaux%vec_list(abs(bc_inf(1)))%time_0.and.(time/=0d0)
        no_tn_bc = vaux%vec_list(abs(bc_inf(1)))%time_n
     $            .and.(.not.save_tn)
        if ((no_eq_bc.or.no_tn_bc).and.(.not.ts_restart)) return
#endif

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

cc          call XformVector_BC(gv%gparams,igr,bcnv_0,bcov_0,.false.)

        case(IJCOV_0)

cc          call XformVector_BC(gv%gparams,igr,jcnv_0,jcov_0,.false.)

        case(IJCNV_0)

        case(IENI_0)

        case(IVECNV_0)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          vecnv_0 = form_ve(vcnv_0,jcnv_0,rho_0)

        case(IBCNV_N)

        case(IVCNV_N)

cc          vcnv_n  = vcnv  !Initializes for BDF2 as well

        case(IVECOV_N)

cc          if (electron_hall) vecov_n= vecov  !Does NOT work for BDF2

        case(IVCOV)

cc          call XformVector_BC(gv%gparams,igr,vcnv,vcov,.false.)

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
        case(IACOV_N)

        case(IACOV)

          varray%array_var(IAX)%array = acov(:,:,:,1)
          varray%array_var(IAY)%array = acov(:,:,:,2)
          varray%array_var(IAZ)%array = acov(:,:,:,3)

        case(IACNV)

cc          acnv = XformVector(gv%gparams,igr,acov,.true.)

        case(IBCNV)

#else

        case(IBCNV)

          varray%array_var(IBX)%array = bcnv(:,:,:,1)
          varray%array_var(IBY)%array = bcnv(:,:,:,2)
          varray%array_var(IBZ)%array = bcnv(:,:,:,3)
#endif

        case(IBCOV)

cc          call XformVector_BC(gv%gparams,igr,bcnv,bcov,.false.)

        case(IJCOV)

        case(IJCNV)

        case(IVECNV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          vecnv = form_ve(vcnv,jcnv,rho)

        case(IVECOV)

c       Electron velocity ve=v-di*j/rho with j=curl(B)

          call XformVector_BC(gv%gparams,igr,vecnv,vecov,.false.)

        case(IDIVPI)

        case(IDIVPE)

        case(IEH)

        case(IENI)

c       Total E field

          if (subtract_E0) E_ni = E_ni - E_ni_0  !Substract imposed electric field

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
      subroutine scalarDirichletBC(g_def,array,array0,ieq,dim,loc,order)
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

      type(grid_mg_def),pointer :: g_def

      integer    :: ieq,dim,loc,order
      real(8)    :: array (0:,0:,0:)
     .             ,array0(0:,0:,0:)

c     Local variables

c     Begin program

      call extrapolateBC(g_def,array,array0,dim,loc,order)
cc      call extrapolateBC(g_def,array,array0,dim,loc,1)

c     End program

      end subroutine scalarDirichletBC

c vectorDirichletBC
c #####################################################################
      subroutine vectorDirichletBC(g_def,ivar,array,array0,ieq,dim,loc
     .                            ,order)
c ---------------------------------------------------------------------
c     Imposes dirichlet BC. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c        * order -> order of extrapolation (when used)
c ---------------------------------------------------------------------
        
      use rw_bc

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer    :: ieq,dim,loc,order,ivar
      real(8)    :: array (0:,0:,0:,:)
     .             ,array0(0:,0:,0:,:)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: i,j,k,ig,jg,kg,ibc,hex,hey,hez
      real(8)    :: x1,x2,x3,dh(3),diver

      real(8)    :: gsuper(3,3),jac0

      real(8),allocatable,dimension(:,:) :: rhs

      INTERFACE
        subroutine solen_bc(g_def,array,dim,loc,face,safe)
          use grid_def_st
          type(grid_mg_def),pointer :: g_def
          integer :: dim,loc
          logical :: face,safe
          real(8) :: array(0:,0:,0:,:)
        end subroutine solen_bc
      END INTERFACE

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      select case (ieq)
      case (IVX,IVY,IVZ)

        call extrapolateBC(g_def,array (:,:,:,ivar)
     .                    ,array0(:,:,:,ivar)
     .                    ,dim,loc,order)
cc     .                    ,dim,loc,1)  !Required for mass conservation

      case (IBX,IBY,IBZ,IJX,IJY,IJZ)

        if (ivar /= dim) then

          call extrapolateBC(g_def,array (:,:,:,ivar)
     .                      ,array0(:,:,:,ivar)
     .                      ,dim,loc,order)

        else  !Imposes divergence-free constraint on B, J

          call solen_bc(g_def,array,dim,loc,.false.,.false.)

        endif

      case (IAX,IAY,IAZ)

c$$$        if (ivar/=dim.or.bcond(ibc)==FSYM.or.res_wall_bc) then
        if (ivar/=dim.or.bcond(ibc)==FSYM) then

          call extrapolateBC(g_def,array (:,:,:,ivar)
     .                      ,array0(:,:,:,ivar)
     .                      ,dim,loc,order)

        else  !Imposes divergence-free constraint on A @ physical boundaries

          call solen_bc(g_def,array,dim,loc,.true.,.false.)

        endif
        
      case default

        messg = 'Variable '//int2char(ieq)//' does not exist'
        call pstop('vectorDirichletBC',messg)

      end select

c     End program

      end subroutine vectorDirichletBC

c scalarNeumannBC
c #####################################################################
      subroutine scalarNeumannBC(g_def,array,array0,ieq,dim,loc)
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

      type(grid_mg_def),pointer :: g_def

      integer    :: ieq,dim,loc
      real(8)    :: array (0:,0:,0:)
      real(8)    :: array0(0:,0:,0:)

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
cc              gsuper = gv%gparams%gmetric%grid(igxbc)%gsup(i,j,k,:,:)
cc
cc              ip = i+1
cc              im = i-1
cc              dh(1) = 2.*gv%gparams%dxh(ig)
cc              if(isBdry(i,igxbc,1).and.(.not.isBC(i,igxbc,1,PER)))then
cc                dh(1) = gv%gparams%dx(ig)
cc                im = i
cc              endif
cc              if(isBdry(i,igxbc,2).and.(.not.isBC(i,igxbc,2,PER)))then
cc                dh(1) = gv%gparams%dx(ig-1)
cc                ip = i
cc              endif
cc
cc              jp = j+1
cc              jm = j-1
cc              dh(2) = 2.*gv%gparams%dyh(jg)
cc              if(isBdry(j,igybc,3).and.(.not.isBC(j,igybc,3,PER)))then
cc                dh(2) = gv%gparams%dy(jg)
cc                jm = j
cc              endif
cc              if(isBdry(j,igybc,4).and.(.not.isBC(j,igybc,4,PER)))then
cc                dh(2) = gv%gparams%dy(jg-1)
cc                jp = j
cc              endif
cc
cc              kp = k+1
cc              km = k-1
cc              dh(3) = 2.*gv%gparams%dzh(kg)
cc              if(isBdry(k,igzbc,5).and.(.not.isBC(k,igzbc,5,PER)))then
cc                dh(3) = gv%gparams%dz(kg)
cc                km = k
cc              endif
cc              if(isBdry(k,igzbc,6).and.(.not.isBC(k,igzbc,6,PER)))then
cc                dh(3) = gv%gparams%dz(kg-1)
cc                kp = k
cc              endif
cccc              ip = min(i+1,nxbc)
cccc              im = max(i-1,1)
cccc              jp = min(j+1,nybc)
cccc              jm = max(j-1,1)
cccc              kp = min(k+1,nzbc)
cccc              km = max(k-1,1)
cccc
cccc              dh(1) = 2.*gv%gparams%dxh(ig)
cccc              if (i == nxbc) dh(1) = gv%gparams%dx(ig-1)
cccc              if (i == 1   ) dh(1) = gv%gparams%dx(ig)
cccc
cccc              dh(2) = 2.*gv%gparams%dyh(jg)
cccc              if (j == nybc) dh(2) = gv%gparams%dy(jg-1)
cccc              if (j == 1   ) dh(2) = gv%gparams%dy(jg)
cccc
cccc              dh(3) = 2.*gv%gparams%dzh(kg)
cccc              if (k == nzbc) dh(3) = gv%gparams%dz(kg-1)
cccc              if (k == 1   ) dh(3) = gv%gparams%dz(kg)
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
cc      case (IPRS)
cc
cc        do i=imin,imax
cc          do j=jmin,jmax
cc            do k=kmin,kmax
cc
cc              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
cc
cc              gsuper = gv%gparams%gmetric%grid(igxbc)%gsup(i,j,k,:,:)
cc              jac0   = gv%gparams%gmetric%grid(igxbc)%jac (i,j,k)
cc
cc              hessian = gv%gparams%gmetric%grid(igxbc)%Gamma(i,j,k,:,:,:)
cc
cc              ip = i+1
cc              im = i-1
cc              dh(1) = 2.*gv%gparams%dxh(ig)
cc              if(isBdry(i,igxbc,1).and.(.not.isBC(i,igxbc,1,PER)))then
cc                dh(1) = gv%gparams%dx(ig)
cc                im = i
cc              endif
cc              if(isBdry(i,igxbc,2).and.(.not.isBC(i,igxbc,2,PER)))then
cc                dh(1) = gv%gparams%dx(ig-1)
cc                ip = i
cc              endif
cc
cc              jp = j+1
cc              jm = j-1
cc              dh(2) = 2.*gv%gparams%dyh(jg)
cc              if(isBdry(j,igybc,3).and.(.not.isBC(j,igybc,3,PER)))then
cc                dh(2) = gv%gparams%dy(jg)
cc                jm = j
cc              endif
cc              if(isBdry(j,igybc,4).and.(.not.isBC(j,igybc,4,PER)))then
cc                dh(2) = gv%gparams%dy(jg-1)
cc                jp = j
cc              endif
cc
cc              kp = k+1
cc              km = k-1
cc              dh(3) = 2.*gv%gparams%dzh(kg)
cc              if(isBdry(k,igzbc,5).and.(.not.isBC(k,igzbc,5,PER)))then
cc                dh(3) = gv%gparams%dz(kg)
cc                km = k
cc              endif
cc              if(isBdry(k,igzbc,6).and.(.not.isBC(k,igzbc,6,PER)))then
cc                dh(3) = gv%gparams%dz(kg-1)
cc                kp = k
cc              endif
cccc              ip = min(i+1,nxbc)
cccc              im = max(i-1,1   )
cccc              jp = min(j+1,nybc)
cccc              jm = max(j-1,1   )
cccc              kp = min(k+1,nzbc)
cccc              km = max(k-1,1   )
cccc
cccc              dh(1) = 2.*gv%gparams%dxh(ig)
cccc              if (i == nxbc) dh(1) = gv%gparams%dx(ig-1)
cccc              if (i == 1   ) dh(1) = gv%gparams%dx(ig)
cccc
cccc              dh(2) = 2.*gv%gparams%dyh(jg)
cccc              if (j == nybc) dh(2) = gv%gparams%dy(jg-1)
cccc              if (j == 1   ) dh(2) = gv%gparams%dy(jg)
cccc
cccc              dh(3) = 2.*gv%gparams%dzh(kg)
cccc              if (k == nzbc) dh(3) = gv%gparams%dz(kg-1)
cccc              if (k == 1   ) dh(3) = gv%gparams%dz(kg)
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
      subroutine vectorNeumannBC(g_def,ivar,array,array0,ieq,dim,loc)
c     -----------------------------------------------------------------
c     Imposes neumann BC for a scalar. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     -----------------------------------------------------------------

      use local_BCS_variables, vec_mod=>vec

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: ieq,dim,loc,ivar
      real(8) :: array (0:,0:,0:,:)
     .          ,array0(0:,0:,0:,:)

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

      logical :: nonlinear,safe

      real(8),allocatable,dimension(:,:) :: rhs

      INTERFACE
        subroutine solen_bc(g_def,array,dim,loc,face,safe)
          use grid_def_st
          type(grid_mg_def),pointer :: g_def
          integer :: dim,loc
          logical :: face,safe
          real(8) :: array(0:,0:,0:,:)
        end subroutine solen_bc

        subroutine curl_BC(g_def,array,array0,ivar,dim,loc,face,safe)
          use grid_def_st
          type(grid_mg_def),pointer :: g_def
          integer :: dim,loc,ivar
          logical :: face,safe
          real(8) :: array(0:,0:,0:,:),array0(0:,0:,0:,:)
        end subroutine curl_BC
      END INTERFACE

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

      call findBCLoopLimits(g_def
     .                     ,dim,loc,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      select case (ieq)
      case (IVX,IVY,IVZ) !Velocity components

        if ((.not.no_map).and.(.not.bcSYM(ibc))) then

          call findBCLoopLimits(g_def,dim,loc
     .                         ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                         ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call getMGmap(g_def,i,j,k,igxbc,igybc,igzbc
     .                       ,ig,jg,kg)

                gsuper = g_def%gmetric%grid(igxbc)%gsup(i,j,k,:,:)

                hessian= g_def%gmetric%grid(igxbc)
     .                     %Gamma(i,j,k,:,:,:)

                ip = i+1
                im = i-1
                dh(1) = 2.*g_def%dxh(ig)
                if(         isBdry(g_def,i,igxbc,1)
     .             .and.(.not.isBC(g_def,i,igxbc,1,PER)))then
                  dh(1) = g_def%dx(ig)
                  im = i
                endif
                if(         isBdry(g_def,i,igxbc,2)
     .             .and.(.not.isBC(g_def,i,igxbc,2,PER)))then
                  dh(1) = g_def%dx(ig-1)
                  ip = i
                endif

                jp = j+1
                jm = j-1
                dh(2) = 2.*g_def%dyh(jg)
                if(         isBdry(g_def,j,igybc,3)
     .             .and.(.not.isBC(g_def,j,igybc,3,PER)))then
                  dh(2) = g_def%dy(jg)
                  jm = j
                endif
                if(         isBdry(g_def,j,igybc,4)
     .             .and.(.not.isBC(g_def,j,igybc,4,PER)))then
                  dh(2) = g_def%dy(jg-1)
                  jp = j
                endif

                kp = k+1
                km = k-1
                dh(3) = 2.*g_def%dzh(kg)
                if(         isBdry(g_def,k,igzbc,5)
     .             .and.(.not.isBC(g_def,k,igzbc,5,PER)))then
                  dh(3) = g_def%dz(kg)
                  km = k
                endif
                if(         isBdry(g_def,k,igzbc,6)
     .             .and.(.not.isBC(g_def,k,igzbc,6,PER)))then
                  dh(3) = g_def%dz(kg-1)
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

          safe = (dim == 1 .and. bcFSYM(2) .and. bcFSYM(3))
     .       .or.(dim == 2 .and. bcFSYM(1) .and. bcFSYM(3))
     .       .or.(dim == 3 .and. bcFSYM(1) .and. bcFSYM(2))

          call solen_BC(g_def,array,dim,loc,.false.,safe)

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

          call solen_BC(g_def,array,dim,loc,.true.,.true.)

        endif

      case (IBX,IBY,IBZ) 

        if (ivar /= dim) then

          safe = (dim == 1 .and. bcFSYM(2) .and. bcFSYM(3))
     .       .or.(dim == 2 .and. bcFSYM(1) .and. bcFSYM(3))
     .       .or.(dim == 3 .and. bcFSYM(1) .and. bcFSYM(2))
          
          call curl_BC(g_def,array,array0,ivar,dim,loc,.true.,safe)

        else  !Imposes divergence-free constraint on B

          safe = (dim == 1 .and. bcFSYM(2) .and. bcFSYM(3))
     .       .or.(dim == 2 .and. bcFSYM(1) .and. bcFSYM(3))
     .       .or.(dim == 3 .and. bcFSYM(1) .and. bcFSYM(2))

          call solen_BC(g_def,array,dim,loc,.false.,safe)

        endif

      case default

        select case (ibc)
        case (1)
          array(0     ,:,:,ivar) = array(1   +offset,:,:,ivar)
        case (2)                                             
          array(nxbc+1,:,:,ivar) = array(nxbc-offset,:,:,ivar)
        case (3)                                             
          array(:,0     ,:,ivar) = array(:,1   +offset,:,ivar)
        case (4)                                             
          array(:,nybc+1,:,ivar) = array(:,nybc-offset,:,ivar)
        case (5)                                             
          array(:,:,0     ,ivar) = array(:,:,1   +offset,ivar)
        case (6)                                             
          array(:,:,nzbc+1,ivar) = array(:,:,nzbc-offset,ivar)
        end select

      end select

c     Assign value

      deallocate(rhs)

c     End program

      end subroutine vectorNeumannBC

c solen_BC
c #####################################################################
      subroutine solen_BC(g_def,array,dim,loc,face,safe)
c ---------------------------------------------------------------------
c     Imposes solenoidal BC. On input:
c        * array -> contains vector
c        * dim   -> dimension we are imposing BC on (X,Y,Z) (1,2,3)
c        * loc   -> boundary location (0 -> left, 1->right)
c        * face  -> whether we are approximating normal derivative at
c                   boundary face or at inner boundary cell center.
c ---------------------------------------------------------------------
        
      use local_BCS_variables

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: dim,loc
      real(8) :: array(0:,0:,0:,:)
      logical :: face,safe

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

      call findBCLoopLimits(g_def,dim,loc
cc     .                     ,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                     ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                     ,imin ,imax ,jmin ,jmax ,kmin ,kmax )

      hex = 0 ; hey = 0 ; hez = 0

      select case (ibc)
      case (1)
 
        i = 1

        array(i-1,jmin:jmax,kmin:kmax,dim)
     .       = array(i+1,jmin:jmax,kmin:kmax,dim) !To cancel normal derivative

        do k=kmin,kmax
          do j=jmin,jmax

            if (safe) then
              hex = 0
cc              if(isBdry(i,igxbc,1)) hex = 1
cc              if(isBdry(i,igxbc,2)) hex =-1

              hey = 0
              if(isBdry(g_def,j,igybc,3)) hey = 1
              if(isBdry(g_def,j,igybc,4)) hey =-1

              hez = 0
              if(isBdry(g_def,k,igzbc,5)) hez = 1
              if(isBdry(g_def,k,igzbc,6)) hez =-1
            endif

            call getMGmap(g_def,i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
            jac0 = g_def%gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*g_def%dxh(ig)
            dh(2) = coeff*g_def%dyh(jg)
            dh(3) = coeff*g_def%dzh(kg)

            diver=jac0*div(g_def,i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                    ,array,hex=hex,hey=hey,hez=hez)

            array(i-1,j,k,dim) = array(i+offset,j,k,dim) + dh(dim)*diver

          enddo
        enddo
 
      case (2)

        i = nxbc

        array(i+1,jmin:jmax,kmin:kmax,dim)
     .       = array(i-1,jmin:jmax,kmin:kmax,dim) !To cancel normal derivative

        do k=kmin,kmax
          do j=jmin,jmax

            if (safe) then
              hex = 0
cc              if(isBdry(i,igxbc,1)) hex = 1
cc              if(isBdry(i,igxbc,2)) hex =-1

              hey = 0
              if(isBdry(g_def,j,igybc,3)) hey = 1
              if(isBdry(g_def,j,igybc,4)) hey =-1

              hez = 0
              if(isBdry(g_def,k,igzbc,5)) hez = 1
              if(isBdry(g_def,k,igzbc,6)) hez =-1
            endif

            call getMGmap(g_def,i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

            jac0 = g_def%gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*g_def%dxh(ig)
            dh(2) = coeff*g_def%dyh(jg)
            dh(3) = coeff*g_def%dzh(kg)

            diver=jac0*div(g_def,i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                    ,array,hex=hex,hey=hey,hez=hez)

            array(i+1,j,k,dim) = array(i-offset,j,k,dim) - dh(dim)*diver

          enddo
        enddo

      case (3)

        j = 1

        array(imin:imax,j-1,kmin:kmax,dim)
     .       = array(imin:imax,j+1,kmin:kmax,dim) !To cancel normal derivative

        do k=kmin,kmax
          do i=imin,imax

            if (safe) then
              hex = 0
              if(isBdry(g_def,i,igxbc,1)) hex = 1
              if(isBdry(g_def,i,igxbc,2)) hex =-1

              hey = 0
cc              if(isBdry(j,igybc,3)) hey = 1
cc              if(isBdry(j,igybc,4)) hey =-1

              hez = 0
              if(isBdry(g_def,k,igzbc,5)) hez = 1
              if(isBdry(g_def,k,igzbc,6)) hez =-1
            endif

            call getMGmap(g_def,i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
            jac0 = g_def%gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*g_def%dxh(ig)
            dh(2) = coeff*g_def%dyh(jg)
            dh(3) = coeff*g_def%dzh(kg)

            diver=jac0*div(g_def,i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                    ,array,hex=hex,hey=hey,hez=hez)

            array(i,j-1,k,dim) = array(i,j+offset,k,dim) + dh(dim)*diver
              
          enddo
        enddo

      case (4)

        j = nybc

        array(imin:imax,j+1,kmin:kmax,dim)
     .       = array(imin:imax,j-1,kmin:kmax,dim) !To cancel normal derivative

        do k=kmin,kmax
          do i=imin,imax

            if (safe) then
              hex = 0
              if(isBdry(g_def,i,igxbc,1)) hex = 1
              if(isBdry(g_def,i,igxbc,2)) hex =-1

              hey = 0
cc              if(isBdry(j,igybc,3)) hey = 1
cc              if(isBdry(j,igybc,4)) hey =-1

              hez = 0
              if(isBdry(g_def,k,igzbc,5)) hez = 1
              if(isBdry(g_def,k,igzbc,6)) hez =-1
            endif

            call getMGmap(g_def,i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
            jac0 = g_def%gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*g_def%dxh(ig)
            dh(2) = coeff*g_def%dyh(jg)
            dh(3) = coeff*g_def%dzh(kg)

            diver=jac0*div(g_def,i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                    ,array,hex=hex,hey=hey,hez=hez)

            array(i,j+1,k,dim) = array(i,j-offset,k,dim) - dh(dim)*diver
              
          enddo
        enddo

      case (5)

        k = 1

        array(imin:imax,jmin:jmax,k-1,dim)
     .       = array(imin:imax,jmin:jmax,k+1,dim) !To cancel normal derivative

        do j=jmin,jmax
          do i=imin,imax

            if (safe) then
              hex = 0
              if(isBdry(g_def,i,igxbc,1)) hex = 1
              if(isBdry(g_def,i,igxbc,2)) hex =-1

              hey = 0
              if(isBdry(g_def,j,igybc,3)) hey = 1
              if(isBdry(g_def,j,igybc,4)) hey =-1

              hez = 0
cc              if(isBdry(k,igzbc,5)) hez = 1
cc              if(isBdry(k,igzbc,6)) hez =-1
            endif

            call getMGmap(g_def,i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
            jac0 = g_def%gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*g_def%dxh(ig)
            dh(2) = coeff*g_def%dyh(jg)
            dh(3) = coeff*g_def%dzh(kg)

            diver=jac0*div(g_def,i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                    ,array,hex=hex,hey=hey,hez=hez)

            array(i,j,k-1,dim) = array(i,j,k+offset,dim) + dh(dim)*diver

          enddo
        enddo

      case (6)

        k = nzbc

        array(imin:imax,jmin:jmax,k+1,dim)
     .       = array(imin:imax,jmin:jmax,k-1,dim) !To cancel normal derivative

        do j=jmin,jmax
          do i=imin,imax

            if (safe) then
              hex = 0
              if(isBdry(g_def,i,igxbc,1)) hex = 1
              if(isBdry(g_def,i,igxbc,2)) hex =-1

              hey = 0
              if(isBdry(g_def,j,igybc,3)) hey = 1
              if(isBdry(g_def,j,igybc,4)) hey =-1

              hez = 0
cc              if(isBdry(k,igzbc,5)) hez = 1
cc              if(isBdry(k,igzbc,6)) hez =-1
            endif

            call getMGmap(g_def,i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
            jac0 = g_def%gmetric%grid(igxbc)%jac(i,j,k)

            dh(1) = coeff*g_def%dxh(ig)
            dh(2) = coeff*g_def%dyh(jg)
            dh(3) = coeff*g_def%dzh(kg)

            diver=jac0*div(g_def,i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .                    ,array,hex=hex,hey=hey,hez=hez)

            array(i,j,k+1,dim) = array(i,j,k-offset,dim) - dh(dim)*diver

          enddo
        enddo

      case default

        messg = 'Boundary'//int2char(ibc)//' non existent'
        call pstop('solen_BC',messg)

      end select

      end subroutine solen_BC

c curl_BC
c #####################################################################
      subroutine curl_BC(g_def,array,array0,ivar,dim,loc,face,safe)
c ---------------------------------------------------------------------
c     Imposes curl BC. On input:
c        * array,array0 -> contain vectors
c        * ivar  -> defines component to apply BC on
c        * dim   -> dimension we are imposing BC on (X,Y,Z) (1,2,3)
c        * loc   -> boundary location (0 -> left, 1->right)
c        * face  -> whether we are approximating normal derivative at
c                   boundary face or at inner boundary cell center.
c ---------------------------------------------------------------------
        
      use local_BCS_variables

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: dim,loc,ivar
      real(8) :: array(0:,0:,0:,:),array0(0:,0:,0:,:)
      logical :: face,safe

c     Local variables

      integer :: imax,imin,jmax,jmin,kmax,kmin,offset
     .          ,im,jm,km,ip,jp,kp
      integer :: i,j,k,ig,jg,kg,ibc,hex,hey,hez,icomp,isgn
      real(8) :: dh(3),coeff,crl(3)
      real(8) :: gsuper(3,3),gspavg(3,3)
     .          ,gsupip(3,3),gsupim(3,3)
     .          ,gsupjp(3,3),gsupjm(3,3)
     .          ,gsupkp(3,3),gsupkm(3,3)

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      if (face) then
        offset= 0
        coeff = 1d0
      else
        offset= 1
        coeff = 2d0
      endif

cc      call findBCLoopLimits(g_def,dim,loc
cc     .                     ,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
cc     .                     ,imin ,imax ,jmin ,jmax ,kmin ,kmax )
      call findBCLoopLimits(g_def,dim,loc
     .                     ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                     ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      hex = 0 ; hey = 0 ; hez = 0

      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax

            ip = i+1
            im = i-1
            jp = j+1
            jm = j-1
            kp = k+1
            km = k-1

            if (safe) then
              hex = 0
              if(isBdry(g_def,i,igxbc,1).and.(.not.bcPER(1))) then
                hex = 1 ; im = i
              endif
              if(isBdry(g_def,i,igxbc,2).and.(.not.bcPER(1))) then
                hex =-1 ; ip = i
              endif

              hey = 0
              if(isBdry(g_def,j,igybc,3).and.(.not.bcPER(2))) then
                hey = 1 ; jm = j
              endif
              if(isBdry(g_def,j,igybc,4).and.(.not.bcPER(2))) then
                hey =-1 ; jp = j
              endif

              hez = 0
              if(isBdry(g_def,k,igzbc,5).and.(.not.bcPER(3))) then
                hez = 1 ; km = k
              endif
              if(isBdry(g_def,k,igzbc,6).and.(.not.bcPER(3))) then
                hez =-1 ; kp = k
              endif
            endif

            call getMGmap(g_def,i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

            gsuper= g_def%gmetric%grid(igxbc)%gsup(i ,j,k,:,:)
            gsupip= g_def%gmetric%grid(igxbc)%gsup(ip,j,k,:,:)
            gsupim= g_def%gmetric%grid(igxbc)%gsup(im,j,k,:,:)
            gsupjp= g_def%gmetric%grid(igxbc)%gsup(i,jp,k,:,:)
            gsupjm= g_def%gmetric%grid(igxbc)%gsup(i,jm,k,:,:)
            gsupkp= g_def%gmetric%grid(igxbc)%gsup(i,j,kp,:,:)
            gsupkm= g_def%gmetric%grid(igxbc)%gsup(i,j,km,:,:)

            crl=curl(g_def,i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc
     .              ,array,hex=hex,hey=hey,hez=hez)

            dh(1) = (ip-im)*g_def%dxh(ig)
            dh(2) = (jp-jm)*g_def%dyh(jg)
            dh(3) = (kp-km)*g_def%dzh(kg)

            select case (dim)
            case (1)

              dh(1) = coeff*g_def%dxh(ig)

              if (face) then
                if (ibc == 1) then
                  gspavg = 0.5*(gsuper+gsupim)
                elseif (ibc == 2) then
                  gspavg = 0.5*(gsuper+gsupip)
                endif
              else
                gspavg = gsuper
              endif
                
              if (ivar == 2) then
                icomp = 3
                isgn  = 1

                crl(dim) =(array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                  +isgn*gspavg(dim,icomp)*crl(dim)/gspavg(dim,dim)
              elseif (ivar == 3) then
                icomp = 2
                isgn  =-1

                crl(dim) =(array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                  +isgn*gspavg(dim,icomp)*crl(dim)/gspavg(dim,dim)
              endif

            case(2)

              dh(2) = coeff*g_def%dyh(jg)

              if (face) then
                if (ibc == 3) then
                  gspavg = 0.5*(gsuper+gsupjm)
                elseif (ibc == 4) then
                  gspavg = 0.5*(gsuper+gsupjp)
                endif
              else
                gspavg = gsuper
              endif
              
              if (ivar == 3) then
                icomp = 1
                isgn  = 1

                crl(dim) =(array(i,j,kp,dim)-array(i,j,km,dim))/dh(ivar)
     .                  +isgn*gspavg(dim,icomp)*crl(dim)/gspavg(dim,dim)
              elseif (ivar == 1) then
                icomp = 3
                isgn  =-1

                crl(dim) =(array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                  +isgn*gspavg(dim,icomp)*crl(dim)/gspavg(dim,dim)
              endif

            case(3)

              dh(3) = coeff*g_def%dzh(kg)

              if (face) then
                if (ibc == 5) then
                  gspavg = 0.5*(gsuper+gsupkm)
                elseif (ibc == 6) then
                  gspavg = 0.5*(gsuper+gsupkp)
                endif
              else
                gspavg = gsuper
              endif

              if (ivar == 1) then
                icomp = 2
                isgn  = 1

                crl(dim) =(array(ip,j,k,dim)-array(im,j,k,dim))/dh(ivar)
     .                  +isgn*gspavg(dim,icomp)*crl(dim)/gspavg(dim,dim)
              elseif (ivar == 2) then
                icomp = 1
                isgn  =-1

                crl(dim) =(array(i,jp,k,dim)-array(i,jm,k,dim))/dh(ivar)
     .                  +isgn*gspavg(dim,icomp)*crl(dim)/gspavg(dim,dim)
              endif

            end select

            select case (ibc)
            case (1)
              crl(dim) = dh(dim)*(crl(dim) + isgn*array0(i-1,j,k,icomp))

              array(0     ,j,k,ivar) = array(1+offset,j,k,ivar)
     .                               - crl(dim)
            case (2)
              crl(dim) = dh(dim)*(crl(dim) + isgn*array0(i+1,j,k,icomp))

              array(nxbc+1,j,k,ivar) = array(nxbc-offset,j,k,ivar)
     $                               + crl(dim)
            case (3)
              crl(dim) = dh(dim)*(crl(dim) + isgn*array0(i,j-1,k,icomp))

              array(i,0     ,k,ivar) = array(i,1+offset,k,ivar)
     $                               - crl(dim)
            case (4)
              crl(dim) = dh(dim)*(crl(dim) + isgn*array0(i,j+1,k,icomp))

              array(i,nybc+1,k,ivar) = array(i,nybc-offset,k,ivar)
     $                               + crl(dim)
            case (5)
              crl(dim) = dh(dim)*(crl(dim) + isgn*array0(i,j,k-1,icomp))

              array(i,j,0     ,ivar) = array(i,j,1+offset,ivar)
     $                               - crl(dim)
            case (6)
              crl(dim) = dh(dim)*(crl(dim) + isgn*array0(i,j,k+1,icomp))

              array(i,j,nzbc+1,ivar) = array(i,j,nzbc-offset,ivar)
     $                               + crl(dim)
            end select

          enddo
        enddo
      enddo

      end subroutine curl_BC

