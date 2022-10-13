c module local_BCS_variables
c####################################################################
      module local_BCS_variables

        use variables

        use imposeBCinterface, ONLY: extrapolateBC

        use auxiliaryVariables

        use operators

        use ts_setup

        use equilibrium

        use nlfunction_setup, ONLY: res,chi_perp,chi_para,vis,ion_hall
     .                             ,E_fld_res_A,E_fld_res_B,E_fld_ihall
     .                             ,E_fld_ehall,EOM_divPi,hres,diff_mesh
     .                             ,nc_eom_v,nc_eom_jxb,find_B
     .                             ,pinch_flow,sym_st,subtract_E0
     .                             ,form_ve,solenoidal,binom_alph

        use xdraw_io

        use app_iosetup

        integer :: order_bc=1

        logical :: regularize_A = .false.  !Whether we regularize A @ SP
     $            ,fill_A0_bc   = .false.
     $            ,res_wall_bc  = .false.
     $            ,rw_is_vacuum = .false.

        logical :: no_eq_bc,save_tn,no_tn_bc
        
      contains

c     defineDepVariables
c     ####################################################################
      subroutine defineDepVariables(varray)
c     --------------------------------------------------------------------
c     Defines boundary conditions of dependent variables. Convention:
c           + bbcs(1) ---> at x0
c           + bbcs(2) ---> at x1
c           + bbcs(3) ---> at y0
c           + bbcs(4) ---> at y1
c           + bbcs(5) ---> at z0
c           + bbcs(6) ---> at z1
c     Definition of BC identifiers is given in "grid_mod.f". In vectors,
c     a negative BC identifier means that BCs are to be imposed on
c     covariant components instead of on (default) contravariant comps.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray

c     Local variables

c     Begin program

      call defineDepBCs(varray)

      call defineDepBC_Deps(varray)

      call defineDepDom_Deps(varray)

c     End program

      end subroutine defineDepVariables

c     defineDepBCs
c     ####################################################################
      subroutine defineDepBCs(varray)
c     --------------------------------------------------------------------
c     Defines boundary conditions of dependent variables. Convention:
c           + bbcs(1) ---> at x0
c           + bbcs(2) ---> at x1
c           + bbcs(3) ---> at y0
c           + bbcs(4) ---> at y1
c           + bbcs(5) ---> at z0
c           + bbcs(6) ---> at z1
c     Definition of BC identifiers is given in "grid_mod.f". In vectors,
c     a negative BC identifier means that BCs are to be imposed on
c     covariant components instead of on (default) contravariant comps.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray

c     Local variables

      integer :: ieq,bbcs(6,neqd),bcsq(6)

c     Begin program

c     Assign BCs

      do ieq=1,neqd
        bbcs(:,ieq) = varray%array_var(ieq)%bconds
      enddo

c     Default RHO boundary conditions 

      call default_SCL_BCs(bbcs(:,IRHO))  !Neumann at walls

c     Default VEL boundary conditions (no-stress)

      call default_V_BCs(bbcs(:,IVX:IVZ))

c     Default B-A boundary conditions

#if !defined(vec_pot)
      call default_B_BCs(bbcs(:,IBX:IBZ))
#else
      call default_A_BCs(bbcs(:,IAX:IAZ))
#endif

c     Default PRS boundary conditions

      call default_SCL_BCs(bbcs(:,IPRS))  !Neumann at walls

c     Exceptions for specific equilibria

      select case (equil)
      case ('tmcyz','3dkh')

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = -NEU
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) =  DIR

        where (bbcs(:,IBX) == DIR) bbcs(:,IBX) = -NEU
        where (bbcs(:,IBY) ==-NEU) bbcs(:,IBY) =  DIR

      case ('ppnch','ppn3d','ppnst','2dhel','rwhel')   !No stress BCs, pinch flow

        pinch_flow = .true.

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX) = EQU  !Normal

        if (gamma > 1d0)
     .       where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU

      case ('ppnsl','p3nsl','ohtor')        !No slip BCs, pinch flow

        pinch_flow = .true.

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = EQU  !Normal

        !Dirichlet boundary conditions for tangential components
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) = EQU  !Tangential
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ) = EQU  !Tangential

        if (gamma > 1d0)
     .       where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU !Do not change EQU

      case ('tbbr1','tbbr2')        !No slip BCs, no pinch flow

        pinch_flow = .false.

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = EQU

        !Dirichlet boundary conditions for tangential components in resistive MHD
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) = EQU
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ) = EQU

        if (gamma > 1d0.and.(.not.adiabatic).and.(di == 0d0)) then
          where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU !Do not change EQU
        endif

      case ('rsw11') !No slip BCs, no pinch flow

        pinch_flow = .false.

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = EQU

        !Dirichlet boundary conditions for tangential components in resistive MHD
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) = EQU
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ) = EQU

        if (gamma > 1d0.and.(.not.adiabatic).and.(di == 0d0)) then
          where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU !Do not change EQU
        endif

      case ('bpnch','bpchT')        !No slip BCs, no pinch flow

        pinch_flow = .false.

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = EQU

        !Dirichlet boundary conditions for tangential components in resistive MHD
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) = EQU
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ) = EQU

        if (gamma > 1d0.and.(.not.adiabatic).and.(di > 0d0)) then
          where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU !Do not change EQU
        endif

      case ('vmec','tor','stell','gstor','liu20')

cc        where (abs(bbcs(:,IBX:IBZ)) == NEU) bbcs(:,IBX:IBZ) = DEF !for testing

        !Dirichlet boundary conditions for tangential components
        where (abs(bbcs(:,IVY)) == NEU) bbcs(:,IVY) = DIR
        where (abs(bbcs(:,IVZ)) == NEU) bbcs(:,IVZ) = DIR

cc        where (abs(bbcs(:,IRHO)) == NEU) bbcs(:,IRHO) = EQU
        where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU

      case ('efit')

cc        where (abs(bbcs(:,IBX:IBZ)) == NEU) bbcs(:,IBX:IBZ) = DEF !for testing

        !Dirichlet boundary conditions for tangential components
        where (abs(bbcs(:,IVY)) == NEU) bbcs(:,IVY) = DIR
        where (abs(bbcs(:,IVZ)) == NEU) bbcs(:,IVZ) = DIR

cc        where (abs(bbcs(:,IRHO)) == NEU) bbcs(:,IRHO) = EQU
        where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU

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

        where (abs(bbcs(:,IPRS)) == DIR) bbcs(:,IPRS) = NEU

      case('gem','diamh')

cc        !Dirichlet boundary conditions for tangential velocity components
cc        where (abs(bbcs(:,IVY )) == NEU) bbcs(:,IVY ) = DIR
cc        where (abs(bbcs(:,IVZ )) == NEU) bbcs(:,IVZ ) = DIR

        where (abs(bbcs(:,IPRS)) == DIR) bbcs(:,IPRS) = NEU

      case('vela2')

c$$$        if (gamma > 1d0)
c$$$     .       where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU

      case('adv2')

        bbcs(1,IRHO) = EQU
        bbcs(2,IRHO) = EXT
        bbcs(3,IRHO) = EQU
        bbcs(4,IRHO) = EXT
        
        bbcs(1,IPRS) = EQU
        bbcs(2,IPRS) = EXT
        bbcs(3,IPRS) = EQU
        bbcs(4,IPRS) = EXT

        where (bcond == DEF) bbcs(:,IVX) = NEU
        where (bcond == DEF) bbcs(:,IVY) = NEU
 
cc      case('sod')
cc
cc        bbcs(1,IRHO) = NEU
cc        bbcs(2,IRHO) = NEU
cc        
cc        bbcs(1,IPRS) = NEU
cc        bbcs(2,IPRS) = NEU
cc
cc        bbcs(1,IVX)  = DIR
cc        bbcs(2,IVX)  = DIR

      case('khihk')

#if defined(vec_pot)
        if (di > 0d0) then !Force D_t(v x n) = 0 to prevent boundary issues
          where (abs(bbcs(:,IVY)) == NEU) bbcs(:,IVY) =-EQU
          where (abs(bbcs(:,IVZ)) == NEU) bbcs(:,IVZ) =-EQU

c$$$          where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU
c$$$          where (abs(bbcs(:,IRHO)) == NEU) bbcs(:,IRHO) = EQU
        endif
#endif

      end select

c     Reassign BCs back

      do ieq=1,neqd
        varray%array_var(ieq)%bconds = bbcs(:,ieq)
      enddo

c     End subroutine

      end subroutine defineDepBCs

c     defineDepBC_Deps
c     ####################################################################
      subroutine defineDepBC_Deps(varray)
c     --------------------------------------------------------------------
c     Defines BC dependencies for dependent variables. Dependencies
c     could be with other dependent variables, or with auxiliary
c     variables. If the latter, multiply variable identifier by
c     "AUX". Dependencies can be "trivial" (i.e., algebraic, set to '0')
c     or "non-trivial" (i.e., differential, set to '1').
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray

c     Local variables

      integer :: ieq

c     Begin program

      !Rho
      varray%array_var(IRHO)%bc_dep_list(1,1) = IRHO_AUX*AUX
      varray%array_var(IRHO)%bc_dep_list(1,2) = SCALAR
      varray%array_var(IRHO)%bc_dep_list(1,3) = 0 !trivial

      !Tmp
cc      varray%array_var(IPRS)%bc_dep_list(1,1) = ITMP_AUX*AUX
      varray%array_var(IPRS)%bc_dep_list(1,1) = IPRS
      varray%array_var(IPRS)%bc_dep_list(1,2) = SCALAR
      varray%array_var(IPRS)%bc_dep_list(1,3) = 0 !trivial

      !V
      if (nc_eom_v) then
        do ieq=IVX,IVZ
          varray%array_var(ieq)%bc_dep_list(1,1) = IVCNV*AUX
          varray%array_var(ieq)%bc_dep_list(1,2) = VECTOR
          varray%array_var(ieq)%bc_dep_list(1,3) = 0 !trivial
        enddo
      else
        do ieq=IVX,IVZ
          varray%array_var(ieq)%bc_dep_list(1,1) = IRHO_AUX*AUX
          varray%array_var(ieq)%bc_dep_list(1,2) = SCALAR
          varray%array_var(ieq)%bc_dep_list(1,3) = 0 !trivial

          varray%array_var(ieq)%bc_dep_list(2,1) = IVCOV*AUX
          varray%array_var(ieq)%bc_dep_list(2,2) = VECTOR
          varray%array_var(ieq)%bc_dep_list(2,3) = 0 !trivial
        enddo
      endif

#if !defined(vec_pot)
      !B
      do ieq=IBX,IBZ
        varray%array_var(ieq)%bc_dep_list(1,1) = IBCOV*AUX
        varray%array_var(ieq)%bc_dep_list(1,2) = VECTOR
        varray%array_var(ieq)%bc_dep_list(1,3) = 0 !trivial
      enddo
#else
      !A
      do ieq=IAX,IAZ
        varray%array_var(ieq)%bc_dep_list(1,1) = IACNV*AUX
        varray%array_var(ieq)%bc_dep_list(1,2) = VECTOR
        varray%array_var(ieq)%bc_dep_list(1,3) = 0 !trivial
      enddo
#endif

c     End subroutine

      end subroutine defineDepBC_Deps

c     defineDepDom_Deps
c     ####################################################################
      subroutine defineDepDom_Deps(varray)
c     --------------------------------------------------------------------
c     Defines domain dependencies for dependent variables. Dependencies
c     could be with other dependent variables, or with auxiliary variables.
c     If the latter, multiply variable identifier by "AUX"
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray

c     Local variables

      integer :: offset,ieq

c     Begin program

      !Rho
      varray%array_var(IRHO)%dom_dep_list(1,1) = IRHO_AUX*AUX
      varray%array_var(IRHO)%dom_dep_list(1,2) = SCALAR
      varray%array_var(IRHO)%dom_dep_list(1,3) = 1 !non-trivial

      varray%array_var(IRHO)%dom_dep_list(2,1) = IVCNV*AUX
      varray%array_var(IRHO)%dom_dep_list(2,2) = VECTOR
      varray%array_var(IRHO)%dom_dep_list(2,3) = 1 !non-trivial

      varray%array_var(IRHO)%dom_dep_list(3,1) = IDFF*AUX
      varray%array_var(IRHO)%dom_dep_list(3,2) = SCALAR
      varray%array_var(IRHO)%dom_dep_list(3,3) = 1 !non-trivial

      !Tmp
      varray%array_var(IPRS)%dom_dep_list(1,1) = ITMP_AUX*AUX
      varray%array_var(IPRS)%dom_dep_list(1,2) = SCALAR
      varray%array_var(IPRS)%dom_dep_list(1,3) = 1 !non-trivial

      varray%array_var(IPRS)%dom_dep_list(2,1) = IVCNV*AUX
      varray%array_var(IPRS)%dom_dep_list(2,2) = VECTOR
      varray%array_var(IPRS)%dom_dep_list(2,3) = 1 !non-trivial

      varray%array_var(IPRS)%dom_dep_list(3,1) = ICHI*AUX
      varray%array_var(IPRS)%dom_dep_list(3,2) = SCALAR
      varray%array_var(IPRS)%dom_dep_list(3,3) = 1 !non-trivial

      varray%array_var(IPRS)%dom_dep_list(4,1) = ICHIPAR*AUX
      varray%array_var(IPRS)%dom_dep_list(4,2) = SCALAR
      varray%array_var(IPRS)%dom_dep_list(4,3) = 1 !non-trivial

      !V
      do ieq=IVX,IVZ
        varray%array_var(ieq)%dom_dep_list(1,1) = IRHO_AUX*AUX
        varray%array_var(ieq)%dom_dep_list(1,2) = SCALAR
        varray%array_var(ieq)%dom_dep_list(1,3) = 1 !non-trivial

        varray%array_var(ieq)%dom_dep_list(2,1) = ITMP_AUX*AUX
        varray%array_var(ieq)%dom_dep_list(2,2) = SCALAR
        varray%array_var(ieq)%dom_dep_list(2,3) = 1 !non-trivial

        varray%array_var(ieq)%dom_dep_list(3,1) = IVCNV*AUX
        varray%array_var(ieq)%dom_dep_list(3,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(3,3) = 1 !non-trivial

        varray%array_var(ieq)%dom_dep_list(4,1) = IDIVPI*AUX
        varray%array_var(ieq)%dom_dep_list(4,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(4,3) = 0 !trivial

        if (di > 0d0) then
          varray%array_var(ieq)%dom_dep_list(5,1) = IDIVPE*AUX
          varray%array_var(ieq)%dom_dep_list(5,2) = VECTOR
          varray%array_var(ieq)%dom_dep_list(5,3) = 0 !trivial
          offset = 1
        else
          offset = 0
        endif

        if (nc_eom_jxb) then
          varray%array_var(ieq)%dom_dep_list(5+offset,1) = IBCNV*AUX
          varray%array_var(ieq)%dom_dep_list(5+offset,2) = VECTOR
          varray%array_var(ieq)%dom_dep_list(5+offset,3) = 0 !trivial

          varray%array_var(ieq)%dom_dep_list(6+offset,1) = IJCNV*AUX
          varray%array_var(ieq)%dom_dep_list(6+offset,2) = VECTOR
          varray%array_var(ieq)%dom_dep_list(6+offset,3) = 0 !trivial
        else
          varray%array_var(ieq)%dom_dep_list(5+offset,1) = IBCNV*AUX
          varray%array_var(ieq)%dom_dep_list(5+offset,2) = VECTOR
          varray%array_var(ieq)%dom_dep_list(5+offset,3) = 1 !non-trivial
        endif
      enddo

#if !defined(vec_pot)
      !B
      do ieq=IBX,IBZ
        varray%array_var(ieq)%dom_dep_list(1,1) = IENI*AUX
        varray%array_var(ieq)%dom_dep_list(1,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(1,3) = 1 !non-trivial

        varray%array_var(ieq)%dom_dep_list(2,1) = IENI_0*AUX
        varray%array_var(ieq)%dom_dep_list(2,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(2,3) = 1 !non-trivial
      enddo
#else
      !A
      do ieq=IAX,IAZ
        varray%array_var(ieq)%dom_dep_list(1,1) = IENI*AUX
        varray%array_var(ieq)%dom_dep_list(1,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(1,3) = 0 !non-trivial

        varray%array_var(ieq)%dom_dep_list(2,1) = IENI_0*AUX
        varray%array_var(ieq)%dom_dep_list(2,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(2,3) = 0 !non-trivial
      enddo
#endif

c     End subroutine

      end subroutine defineDepDom_Deps

c     defineAuxVariables
c     ###################################################################
      subroutine defineAuxVariables(vaux)

c     -------------------------------------------------------------------
c     Defines aux variables, creates dependencies for BC scheduler.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(aux_array),pointer :: vaux

c     Local variables

        integer :: nx,ny,nz

c     Begin program

        call createAuxVariables(vaux)

c     Define BCs

        call defineAuxBCs(vaux)

c     Setup variable dependencies at boundaries

        call defineAuxBCDeps(vaux)

c     Setup variable dependencies in inner domain

        call defineAuxDomainDeps(vaux)

c     End program

      end subroutine defineAuxVariables

c     defineAuxBCs
c     ###################################################################
      subroutine defineAuxBCs(vaux)

c     -------------------------------------------------------------------
c     Defines BCs for auxiliary quantities
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      type(aux_array),pointer :: vaux

c     Local variables

      integer :: ieq,bcs(6,neqd)

c     Begin program

c     Read dependent variable BC setup

      do ieq=1,neqd
        bcs(:,ieq) = gv%u_0%array_var(ieq)%bconds
      enddo

c     Set auxiliary variable BC

      !Density and temperature
      vaux%var_list(IRHO_0)%bconds = bcs(:,IRHO)
      where (vaux%var_list(IRHO_0)%bconds == EQU)
        vaux%var_list(IRHO_0)%bconds = DEF  !Do nothing
      end where

      vaux%var_list(ITMP_0)%bconds = bcs(:,IPRS)
      where (vaux%var_list(ITMP_0)%bconds == EQU)
        vaux%var_list(ITMP_0)%bconds = DEF  !Do nothing
      end where

      vaux%var_list(ITMPE_0)%bconds = bcs(:,IPRS)
      where (vaux%var_list(ITMPE_0)%bconds == EQU)
        vaux%var_list(ITMPE_0)%bconds = DEF !Do nothing
      end where

      vaux%var_list(IRHO_AUX)%bconds = bcs(:,IRHO)
      vaux%var_list(ITMP_AUX)%bconds = bcs(:,IPRS)
      vaux%var_list(ITMPE_AUX)%bconds = bcs(:,IPRS)
      
      if (.not.spitzer) then
        vaux%var_list(IETA_0)%bconds = bcond
        where (vaux%var_list(IETA_0)%bconds == DEF) 
          vaux%var_list(IETA_0)%bconds = EXT
        end where
      else
        vaux%var_list(IETA_0)%bconds = bcs(:,IPRS)
        where (vaux%var_list(IETA_0)%bconds == DIR)
          vaux%var_list(IETA_0)%bconds = EQU
        end where
      endif

      if (.not.spitzer) then
        vaux%var_list(IETA)%bconds = bcond
        where (vaux%var_list(IETA)%bconds == DEF) 
          vaux%var_list(IETA)%bconds = EXT
        end where
      else
        vaux%var_list(IETA)%bconds = bcs(:,IPRS)
        where (vaux%var_list(IETA)%bconds == DIR)
          vaux%var_list(IETA)%bconds = EQU
        end where
      endif

      !Flow
      vaux%vec_list(IVCNV)%bconds = bcs(:,IVX:IVZ)

      vaux%vec_list(IVCOV)%bconds = vaux%vec_list(IVCNV)%bconds

      vaux%vec_list(IVCNV_0)%bconds = vaux%vec_list(IVCNV)%bconds

      if (.not.pinch_flow) then
        where (vaux%vec_list(IVCNV_0)%bconds == EQU)
          vaux%vec_list(IVCNV_0)%bconds = DEF
        end where
      endif

#if defined(vec_pot)
      !Vector potential
      vaux%vec_list(IACNV)%bconds = bcs(:,IAX:IAZ)

      vaux%vec_list(IACOV)%bconds = vaux%vec_list(IACNV)%bconds

      !Magnetic field
      vaux%vec_list(IBCNV)%bconds = vaux%vec_list(IACNV)%bconds

      where (bcond == FSYM)  !Needed because B has different symmetries than A
        vaux%vec_list(IBCNV)%bconds(:,1) = FSYM !Recover symmetry BCs
        vaux%vec_list(IBCNV)%bconds(:,2) = FSYM !Recover symmetry BCs
        vaux%vec_list(IBCNV)%bconds(:,3) = FSYM !Recover symmetry BCs
      end where

      where (vaux%vec_list(IBCNV)%bconds == -EQU)
        vaux%vec_list(IBCNV)%bconds = -NEU
cc        vaux%vec_list(IBCNV)%bconds =-EXT !For debugging B-field BCs
      end where

      call set_vec_symm_bc(vaux%vec_list(IBCNV)%bconds,1)

      vaux%vec_list(IBCOV)%bconds = vaux%vec_list(IBCNV)%bconds

#else
      !Magnetic field
      vaux%vec_list(IBCNV)%bconds = bcs(:,IBX:IBZ)

cc      !For debugging B-field BCs
cc      where (vaux%vec_list(IBCNV)%bconds == -NEU)
cc        vaux%vec_list(IBCNV)%bconds =-EXT
cc      end where

      vaux%vec_list(IBCOV)%bconds = vaux%vec_list(IBCNV)%bconds
#endif

      !Equilibrium magnetic field
      vaux%vec_list(IBCNV_0)%bconds = vaux%vec_list(IBCNV)%bconds

      where (vaux%vec_list(IBCNV_0)%bconds == -NEU)
        vaux%vec_list(IBCNV_0)%bconds =-EXT   !Extrapolate covariant components
      end where

      vaux%vec_list(IBCOV_0)%bconds = vaux%vec_list(IBCNV_0)%bconds

      !Current
      vaux%vec_list(IJCNV)%bconds = vaux%vec_list(IBCNV)%bconds

      where (bcond == FSYM)  !Needed because J has different symmetries than B
        vaux%vec_list(IJCNV)%bconds(:,1) = FSYM !Recover symmetry BCs
        vaux%vec_list(IJCNV)%bconds(:,2) = FSYM !Recover symmetry BCs
        vaux%vec_list(IJCNV)%bconds(:,3) = FSYM !Recover symmetry BCs
      end where

      where (vaux%vec_list(IJCNV)%bconds == -NEU)
        vaux%vec_list(IJCNV)%bconds = EQU !Use cnv representation in applyBC routine!
cc          vaux%vec_list(IJCNV)%bconds =-EQU !Use cov representation in applyBC routine!
cc          vaux%vec_list(IJCNV)%bconds =-EXT !Use contravariant components for extrapolation BCs
      end where
      
      where (vaux%vec_list(IJCNV)%bconds == -EQU)
cc        vaux%vec_list(IJCNV)%bconds = EQU !Use cnv representation in applyBC routine!
cc        vaux%vec_list(IJCNV)%bconds =-EQU !Use cov representation in applyBC routine!
        vaux%vec_list(IJCNV)%bconds = -EXT !Use cov components for extrapolation BCs
      end where

      call set_vec_symm_bc(vaux%vec_list(IJCNV)%bconds,-1)  !Same as velocity

      vaux%vec_list(IJCOV)%bconds = vaux%vec_list(IJCNV)%bconds

      !Equilibrium Current
      vaux%vec_list(IJCNV_0)%bconds = vaux%vec_list(IJCNV)%bconds

      where (abs(vaux%vec_list(IJCNV_0)%bconds) == EQU)
        vaux%vec_list(IJCNV_0)%bconds = DEF  !Do nothing (curl takes care of it)
cc        vaux%vec_list(IJCNV_0)%bconds =-EXT  !Extrapolate covariant components
cc        vaux%vec_list(IJCNV_0)%bconds = EXT  !Extrapolate contravariant components
      end where

      vaux%vec_list(IJCOV_0)%bconds = vaux%vec_list(IJCNV_0)%bconds

#if defined(vec_pot)
      vaux%vec_list(IACOV_0)%bconds = vaux%vec_list(IJCNV_0)%bconds

      select case(equil)
      case('ohtor','gstor')
        where(abs(vaux%vec_list(IACOV_0)%bconds) == EXT)
          vaux%vec_list(IACOV_0)%bconds = DEF !Do nothing
        end where
      end select
#endif

      !Pressure tensors
cc      vaux%vec_list(IDIVPI)%bconds(:,1) = bcond
cc      vaux%vec_list(IDIVPI)%bconds(:,2) = bcond
cc      vaux%vec_list(IDIVPI)%bconds(:,3) = bcond
cc      where (vaux%vec_list(IDIVPI)%bconds == DEF)
cc        vaux%vec_list(IDIVPI)%bconds = EXT
cc      end where
cc
cc      call set_vec_symm_bc(vaux%vec_list(IDIVPI)%bconds,-1)  !Same as velocity

      vaux%vec_list(IDIVPI)%bconds = vaux%vec_list(IVCNV )%bconds

      vaux%vec_list(IDIVPE)%bconds = vaux%vec_list(IDIVPI)%bconds

      !Electric field (resistive piece)
cc      vaux%vec_list(IENI)%bconds = vaux%vec_list(IACOV)%bconds

      !Electric field (Hall piece)
c$$$      call default_A_BCs(vaux%vec_list(IEH)%bconds)
c$$$      if (ion_hall) then
        vaux%vec_list(IEH )%bconds = vaux%vec_list(IVCOV)%bconds
c$$$      else
c$$$        vaux%vec_list(IEH)%bconds(:,1) = bcond
c$$$        vaux%vec_list(IEH)%bconds(:,2) = bcond
c$$$        vaux%vec_list(IEH)%bconds(:,3) = bcond
c$$$        where (vaux%vec_list(IEH)%bconds == DEF)
c$$$          vaux%vec_list(IEH)%bconds = -NEU  !On covariant components
c$$$        end where
c$$$
c$$$        call set_vec_symm_bc(vaux%vec_list(IEH)%bconds,-1) !Same as velocity
c$$$      endif

cc      !Equilibrium electric field (non-ideal piece)
cc      vaux%vec_list(IENI_0)%bconds = vaux%vec_list(IEH)%bconds

c     End program

      end subroutine defineAuxBCs

c     defineAuxBCDeps
c     ###################################################################
      subroutine defineAuxBCDeps(vaux)

c     -------------------------------------------------------------------
c     Defines intervariable dependencies for auxiliary quantities.
c     When setting dependencies on vector components, use secondary
c     component representation (i.e., the one whose only dependency
c     is the other --primary-- representation). This ensures that they
c     both belong to the same BC group. One example:
c        - bcnv is considered primary (depends on other variables)
c        - bcov is secondary (only depends on bcnv)
c     Trivial (0) and nontrivial (1) dependencies can be set. A dependency
c     is trivial if it involves no derivatives, just algebraic relations.
c     It is nontrivial otherwise.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      type(aux_array),pointer :: vaux

c     Local variables

      logical :: pinch_bc
      
c     Begin program

      pinch_bc = pinch_flow.or.res_wall_bc
     
c     Equilibrium quantities

      !Rho0
      vaux%var_list(IRHO_0)%bc_dep_list(1,1) = IRHO_0*AUX
      vaux%var_list(IRHO_0)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(IRHO_0)%bc_dep_list(1,3) = 0 !trivial

      !Temp0
      vaux%var_list(ITMP_0)%bc_dep_list(1,1) = ITMP_0*AUX
      vaux%var_list(ITMP_0)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(ITMP_0)%bc_dep_list(1,3) = 0 !trivial

      vaux%var_list(ITMPE_0)%bc_dep_list(1,1) = ITMPE_0*AUX
      vaux%var_list(ITMPE_0)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(ITMPE_0)%bc_dep_list(1,3) = 0 !trivial

      !Eta0
      vaux%var_list(IETA_0)%bc_dep_list(1,1) = ITMPE_0*AUX
      vaux%var_list(IETA_0)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(IETA_0)%bc_dep_list(1,3) = 0 !trivial

#if defined(vec_pot) 
      !A0 cov
      vaux%vec_list(IACOV_0)%bc_dep_list(1,1) = IACOV_0*AUX
      vaux%vec_list(IACOV_0)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IACOV_0)%bc_dep_list(1,3) = 0 !trivial

      !B0 cnv
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,1) = IACOV_0*AUX
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,3) = 1 !non-trivial
#else
      !B0 cnv
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,1) = IBCNV_0*AUX
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,3) = 0 !trivial
#endif

      !B0 cov
      vaux%vec_list(IBCOV_0)%bc_dep_list(1,1) = IBCNV_0*AUX
      vaux%vec_list(IBCOV_0)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCOV_0)%bc_dep_list(1,3) = 0 !trivial

      !V0 cnv
      if (pinch_bc) then
        vaux%vec_list(IVCNV_0)%bc_dep_list(1,1) = IRHO_0*AUX
        vaux%vec_list(IVCNV_0)%bc_dep_list(1,2) = SCALAR
        vaux%vec_list(IVCNV_0)%bc_dep_list(1,3) = 1 !non-trivial

        vaux%vec_list(IVCNV_0)%bc_dep_list(2,1) = ITMP_0*AUX
        vaux%vec_list(IVCNV_0)%bc_dep_list(2,2) = SCALAR
        vaux%vec_list(IVCNV_0)%bc_dep_list(2,3) = 1 !non-trivial

        vaux%vec_list(IVCNV_0)%bc_dep_list(3,1) = IETA_0*AUX
        vaux%vec_list(IVCNV_0)%bc_dep_list(3,2) = SCALAR
        vaux%vec_list(IVCNV_0)%bc_dep_list(3,3) = 0 !trivial

        vaux%vec_list(IVCNV_0)%bc_dep_list(4,1) = IBCOV_0*AUX
        vaux%vec_list(IVCNV_0)%bc_dep_list(4,2) = VECTOR
        vaux%vec_list(IVCNV_0)%bc_dep_list(4,3) = 0 !trivial

        vaux%vec_list(IVCNV_0)%bc_dep_list(5,1) = IVCNV_0*AUX
        vaux%vec_list(IVCNV_0)%bc_dep_list(5,2) = VECTOR
        vaux%vec_list(IVCNV_0)%bc_dep_list(5,3) = 0 !trivial
      else
        vaux%vec_list(IVCNV_0)%bc_dep_list(1,1) = IVCNV_0*AUX
        vaux%vec_list(IVCNV_0)%bc_dep_list(1,2) = VECTOR
        vaux%vec_list(IVCNV_0)%bc_dep_list(1,3) = 0 !trivial
      endif
      
      !J0 cnv
      if (pinch_bc) then
        vaux%vec_list(IVCNV_0)%bc_dep_list(1,1) = IRHO_0*AUX
        vaux%vec_list(IVCNV_0)%bc_dep_list(1,2) = SCALAR
        vaux%vec_list(IVCNV_0)%bc_dep_list(1,3) = 1 !non-trivial

        vaux%vec_list(IVCNV_0)%bc_dep_list(2,1) = ITMP_0*AUX
        vaux%vec_list(IVCNV_0)%bc_dep_list(2,2) = SCALAR
        vaux%vec_list(IVCNV_0)%bc_dep_list(2,3) = 1 !non-trivial

        vaux%vec_list(IJCNV_0)%bc_dep_list(3,1) = IETA_0*AUX
        vaux%vec_list(IJCNV_0)%bc_dep_list(3,2) = SCALAR
        vaux%vec_list(IJCNV_0)%bc_dep_list(3,3) = 0 !trivial

        vaux%vec_list(IJCNV_0)%bc_dep_list(4,1) = IBCOV_0*AUX
        vaux%vec_list(IJCNV_0)%bc_dep_list(4,2) = VECTOR
        vaux%vec_list(IJCNV_0)%bc_dep_list(4,3) = 1 !non-trivial

        vaux%vec_list(IJCNV_0)%bc_dep_list(5,1) = IJCNV_0*AUX
        vaux%vec_list(IJCNV_0)%bc_dep_list(5,2) = VECTOR
        vaux%vec_list(IJCNV_0)%bc_dep_list(5,3) = 0 !trivial
      else
        vaux%vec_list(IJCNV_0)%bc_dep_list(1,1) = IJCNV_0*AUX
        vaux%vec_list(IJCNV_0)%bc_dep_list(1,2) = VECTOR
        vaux%vec_list(IJCNV_0)%bc_dep_list(1,3) = 0 !trivial
      endif
      
      !J0 cov
      vaux%vec_list(IJCOV_0)%bc_dep_list(1,1) = IJCNV_0*AUX
      vaux%vec_list(IJCOV_0)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCOV_0)%bc_dep_list(1,3) = 0 !trivial

      !V_n cnv
      vaux%vec_list(IVCNV_N)%bc_dep_list(1,1) = IVCNV*AUX
      vaux%vec_list(IVCNV_N)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IVCNV_N)%bc_dep_list(1,3) = 0 !trivial

      !Ve_n cnv
      vaux%vec_list(IVECOV_N)%bc_dep_list(1,1) = IVECOV*AUX
      vaux%vec_list(IVECOV_N)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IVECOV_N)%bc_dep_list(1,3) = 0 !trivial

      !B_n cov
      vaux%vec_list(IBCNV_N)%bc_dep_list(1,1) = IBCNV_N*AUX
      vaux%vec_list(IBCNV_N)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_N)%bc_dep_list(1,3) = 0 !trivial

      !E-field cov
      vaux%vec_list(IENI_0)%bc_dep_list(1,1) = IETA_0*AUX
      vaux%vec_list(IENI_0)%bc_dep_list(1,2) = SCALAR
      vaux%vec_list(IENI_0)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IENI_0)%bc_dep_list(2,1) = IJCOV_0*AUX
      vaux%vec_list(IENI_0)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IENI_0)%bc_dep_list(2,3) = 0 !non-trivial

      vaux%vec_list(IENI_0)%bc_dep_list(3,1) = IVCNV_0*AUX
      vaux%vec_list(IENI_0)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IENI_0)%bc_dep_list(3,3) = 0 !trivial

c$$$#if defined(vec_pot)
c$$$      vaux%vec_list(IENI_0)%bc_dep_list(4,1) = IACOV_0*AUX
c$$$      vaux%vec_list(IENI_0)%bc_dep_list(4,2) = VECTOR
c$$$      vaux%vec_list(IENI_0)%bc_dep_list(4,3) = 1 !non-trivial
c$$$
c$$$#else
      vaux%vec_list(IENI_0)%bc_dep_list(4,1) = IBCNV_0*AUX
      vaux%vec_list(IENI_0)%bc_dep_list(4,2) = VECTOR
      vaux%vec_list(IENI_0)%bc_dep_list(4,3) = 0 !trivial
c$$$#endif

      !Ve cnv = v - di*j/rho
      vaux%vec_list(IVECNV_0)%bc_dep_list(1,1) = IRHO_0*AUX
      vaux%vec_list(IVECNV_0)%bc_dep_list(1,2) = SCALAR
      vaux%vec_list(IVECNV_0)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IVECNV_0)%bc_dep_list(2,1) = IVCNV_0*AUX
      vaux%vec_list(IVECNV_0)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IVECNV_0)%bc_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IVECNV_0)%bc_dep_list(3,1) = IJCNV_0*AUX
      vaux%vec_list(IVECNV_0)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IVECNV_0)%bc_dep_list(3,3) = 0 !trivial

c     Current variables

      !Rho
      vaux%var_list(IRHO_AUX)%bc_dep_list(1,1) = IRHO_AUX*AUX
      vaux%var_list(IRHO_AUX)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(IRHO_AUX)%bc_dep_list(1,3) = 0 !trivial

      !Temp
      vaux%var_list(ITMP_AUX)%bc_dep_list(1,1) = ITMPE_AUX*AUX
      vaux%var_list(ITMP_AUX)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(ITMP_AUX)%bc_dep_list(1,3) = 0 !trivial

      vaux%var_list(ITMP_AUX)%bc_dep_list(2,1) = ITMP_0*AUX
      vaux%var_list(ITMP_AUX)%bc_dep_list(2,2) = SCALAR
      vaux%var_list(ITMP_AUX)%bc_dep_list(2,3) = 0 !trivial

      !Temp_e
      vaux%var_list(ITMPE_AUX)%bc_dep_list(1,1) = ITMPE_0*AUX
      vaux%var_list(ITMPE_AUX)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(ITMPE_AUX)%bc_dep_list(1,3) = 0 !trivial

      vaux%var_list(ITMPE_AUX)%bc_dep_list(2,1) = IRHO_AUX*AUX
      vaux%var_list(ITMPE_AUX)%bc_dep_list(2,2) = SCALAR
      vaux%var_list(ITMPE_AUX)%bc_dep_list(2,3) = 0 !trivial

      vaux%var_list(ITMPE_AUX)%bc_dep_list(3,1) = IPRS
      vaux%var_list(ITMPE_AUX)%bc_dep_list(3,2) = SCALAR
      vaux%var_list(ITMPE_AUX)%bc_dep_list(3,3) = 0 !trivial

      !Eta
      vaux%var_list(IETA)%bc_dep_list(1,1) = IETA_0*AUX
      vaux%var_list(IETA)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(IETA)%bc_dep_list(1,3) = 0 !trivial

      vaux%var_list(IETA)%bc_dep_list(2,1) = ITMPE_AUX*AUX
      vaux%var_list(IETA)%bc_dep_list(2,2) = SCALAR
      vaux%var_list(IETA)%bc_dep_list(2,3) = 0 !trivial

      !Diff
      vaux%var_list(IDFF)%bc_dep_list(1,1) = ITMPE_0*AUX
      vaux%var_list(IDFF)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(IDFF)%bc_dep_list(1,3) = 0 !trivial

      vaux%var_list(IDFF)%bc_dep_list(2,1) = ITMPE_AUX*AUX
      vaux%var_list(IDFF)%bc_dep_list(2,2) = SCALAR
      vaux%var_list(IDFF)%bc_dep_list(2,3) = 0 !trivial

      vaux%var_list(IDFF)%bc_dep_list(3,1) = IRHO_AUX*AUX
      vaux%var_list(IDFF)%bc_dep_list(3,2) = SCALAR
      vaux%var_list(IDFF)%bc_dep_list(3,3) = 0 !trivial

      vaux%var_list(IDFF)%bc_dep_list(4,1) = IBCNV*AUX
      vaux%var_list(IDFF)%bc_dep_list(4,2) = VECTOR
      vaux%var_list(IDFF)%bc_dep_list(4,3) = 0 !trivial

      !Chi
      vaux%var_list(ICHI)%bc_dep_list(1,1) = ITMPE_0*AUX
      vaux%var_list(ICHI)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(ICHI)%bc_dep_list(1,3) = 0 !trivial

      vaux%var_list(ICHI)%bc_dep_list(2,1) = ITMPE_AUX*AUX
      vaux%var_list(ICHI)%bc_dep_list(2,2) = SCALAR
      vaux%var_list(ICHI)%bc_dep_list(2,3) = 0 !trivial

      vaux%var_list(ICHI)%bc_dep_list(3,1) = IRHO_AUX*AUX
      vaux%var_list(ICHI)%bc_dep_list(3,2) = SCALAR
      vaux%var_list(ICHI)%bc_dep_list(3,3) = 0 !trivial

      vaux%var_list(ICHI)%bc_dep_list(4,1) = IBCNV*AUX
      vaux%var_list(ICHI)%bc_dep_list(4,2) = VECTOR
      vaux%var_list(ICHI)%bc_dep_list(4,3) = 0 !trivial

      !Chipar
      vaux%var_list(ICHIPAR)%bc_dep_list(1,1) = ICHI*AUX
      vaux%var_list(ICHIPAR)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(ICHIPAR)%bc_dep_list(1,3) = 0 !trivial

#if defined(vec_pot)
      !A cov
      vaux%vec_list(IACOV)%bc_dep_list(1,1) = IBCOV_0*AUX
      vaux%vec_list(IACOV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IACOV)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IACOV)%bc_dep_list(2,1) = IACOV*AUX
      vaux%vec_list(IACOV)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IACOV)%bc_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IACOV)%bc_dep_list(3,1) = IBCNV_N*AUX
      vaux%vec_list(IACOV)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IACOV)%bc_dep_list(3,3) = 0 !trivial

cc      vaux%vec_list(IACOV)%bc_dep_list(4,1) = IENI*AUX
cc      vaux%vec_list(IACOV)%bc_dep_list(4,2) = VECTOR
cc      vaux%vec_list(IACOV)%bc_dep_list(4,3) = 0 !trivial

      !A cnv
      vaux%vec_list(IACNV)%bc_dep_list(1,1) = IACOV*AUX
      vaux%vec_list(IACNV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IACNV)%bc_dep_list(1,3) = 0 !trivial

      !A_n cov
      vaux%vec_list(IACOV_N)%bc_dep_list(1,1) = IACOV_N*AUX
      vaux%vec_list(IACOV_N)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IACOV_N)%bc_dep_list(1,3) = 0 !trivial

      !B cnv
      vaux%vec_list(IBCNV)%bc_dep_list(1,1) = IACOV*AUX
      vaux%vec_list(IBCNV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV)%bc_dep_list(1,3) = 1 !non-trivial
#else
      !B cnv
      vaux%vec_list(IBCNV)%bc_dep_list(1,1) = IBCNV*AUX
      vaux%vec_list(IBCNV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV)%bc_dep_list(1,3) = 0 !trivial
#endif

      if (.not.pinch_bc) then
        vaux%vec_list(IBCNV)%bc_dep_list(2,1) = IJCOV_0*AUX
        vaux%vec_list(IBCNV)%bc_dep_list(2,2) = VECTOR
        vaux%vec_list(IBCNV)%bc_dep_list(2,3) = 0 !trivial
      else
        vaux%vec_list(IBCNV)%bc_dep_list(2,1) = IRHO_AUX*AUX
        vaux%vec_list(IBCNV)%bc_dep_list(2,2) = SCALAR
        vaux%vec_list(IBCNV)%bc_dep_list(2,3) = 1 !non-trivial

        vaux%vec_list(IBCNV)%bc_dep_list(3,1) = ITMP_AUX*AUX
        vaux%vec_list(IBCNV)%bc_dep_list(3,2) = SCALAR
        vaux%vec_list(IBCNV)%bc_dep_list(3,3) = 1 !non-trivial

        vaux%vec_list(IBCNV)%bc_dep_list(4,1) = IETA*AUX
        vaux%vec_list(IBCNV)%bc_dep_list(4,2) = SCALAR
        vaux%vec_list(IBCNV)%bc_dep_list(4,3) = 0 !trivial

        vaux%vec_list(IBCNV)%bc_dep_list(5,1) = IBCNV_N*AUX
        vaux%vec_list(IBCNV)%bc_dep_list(5,2) = VECTOR
        vaux%vec_list(IBCNV)%bc_dep_list(5,3) = 0 !trivial
      endif

      !B cov
      vaux%vec_list(IBCOV)%bc_dep_list(1,1) = IBCNV*AUX
      vaux%vec_list(IBCOV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCOV)%bc_dep_list(1,3) = 0 !trivial

      !Nu=rho*nu_i(T)
      vaux%var_list(INU)%bc_dep_list(1,1) = IRHO_AUX*AUX
      vaux%var_list(INU)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(INU)%bc_dep_list(1,3) = 0 !trivial

      vaux%var_list(INU)%bc_dep_list(2,1) = ITMP_AUX*AUX
      vaux%var_list(INU)%bc_dep_list(2,2) = SCALAR
      vaux%var_list(INU)%bc_dep_list(2,3) = 0 !trivial

      !Heta=rho*nu_e(T)
      vaux%var_list(IHETA)%bc_dep_list(1,1) = IRHO*AUX
      vaux%var_list(IHETA)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(IHETA)%bc_dep_list(1,3) = 0 !trivial

      vaux%var_list(IHETA)%bc_dep_list(2,1) = IBCOV_0*AUX
      vaux%var_list(IHETA)%bc_dep_list(2,2) = VECTOR
cc      vaux%var_list(IHETA)%bc_dep_list(2,3) = 0 !trivial
      vaux%var_list(IHETA)%bc_dep_list(2,3) = 1 !non-trivial

      !V cnv
      vaux%vec_list(IVCNV)%bc_dep_list(1,1) = IVCNV*AUX
      vaux%vec_list(IVCNV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IVCNV)%bc_dep_list(1,3) = 0 !trivial

      if (pinch_bc) then
        vaux%vec_list(IVCNV)%bc_dep_list(2,1) = IRHO_AUX*AUX
        vaux%vec_list(IVCNV)%bc_dep_list(2,2) = SCALAR
        vaux%vec_list(IVCNV)%bc_dep_list(2,3) = 1 !non-trivial
                                          
        vaux%vec_list(IVCNV)%bc_dep_list(3,1) = ITMP_AUX*AUX
        vaux%vec_list(IVCNV)%bc_dep_list(3,2) = SCALAR
        vaux%vec_list(IVCNV)%bc_dep_list(3,3) = 1 !non-trivial

        vaux%vec_list(IVCNV)%bc_dep_list(4,1) = IETA*AUX
        vaux%vec_list(IVCNV)%bc_dep_list(4,2) = SCALAR
        vaux%vec_list(IVCNV)%bc_dep_list(4,3) = 0 !trivial

        vaux%vec_list(IVCNV)%bc_dep_list(5,1) = IBCOV*AUX
        vaux%vec_list(IVCNV)%bc_dep_list(5,2) = VECTOR
        vaux%vec_list(IVCNV)%bc_dep_list(5,3) = 0 !trivial

        vaux%vec_list(IVCNV)%bc_dep_list(6,1) = IBCNV_N*AUX
        vaux%vec_list(IVCNV)%bc_dep_list(6,2) = VECTOR
        vaux%vec_list(IVCNV)%bc_dep_list(6,3) = 0 !trivial
      else
        vaux%vec_list(IVCNV)%bc_dep_list(2,1) = IVCNV_0*AUX
        vaux%vec_list(IVCNV)%bc_dep_list(2,2) = VECTOR
        vaux%vec_list(IVCNV)%bc_dep_list(2,3) = 0 !trivial
      endif

      !V cov
      vaux%vec_list(IVCOV)%bc_dep_list(1,1) = IVCOV*AUX
      vaux%vec_list(IVCOV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IVCOV)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IVCOV)%bc_dep_list(2,1) = IVCNV*AUX
      vaux%vec_list(IVCOV)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IVCOV)%bc_dep_list(2,3) = 0 !trivial

      !J cnv
      vaux%vec_list(IJCNV)%bc_dep_list(1,1) = IJCNV*AUX
      vaux%vec_list(IJCNV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCNV)%bc_dep_list(1,3) = 0 !trivial

      if (pinch_bc) then
        vaux%vec_list(IJCNV)%bc_dep_list(2,1) = IRHO_AUX*AUX
        vaux%vec_list(IJCNV)%bc_dep_list(2,2) = SCALAR
        vaux%vec_list(IJCNV)%bc_dep_list(2,3) = 1 !non-trivial
                                          
        vaux%vec_list(IJCNV)%bc_dep_list(3,1) = ITMP_AUX*AUX
        vaux%vec_list(IJCNV)%bc_dep_list(3,2) = SCALAR
        vaux%vec_list(IJCNV)%bc_dep_list(3,3) = 1 !non-trivial
                                          
        vaux%vec_list(IJCNV)%bc_dep_list(4,1) = IETA*AUX
        vaux%vec_list(IJCNV)%bc_dep_list(4,2) = SCALAR
        vaux%vec_list(IJCNV)%bc_dep_list(4,3) = 0 !trivial

        vaux%vec_list(IJCNV)%bc_dep_list(5,1) = IBCOV*AUX
        vaux%vec_list(IJCNV)%bc_dep_list(5,2) = VECTOR
        vaux%vec_list(IJCNV)%bc_dep_list(5,3) = 1 !non-trivial

        vaux%vec_list(IJCNV)%bc_dep_list(6,1) = IBCNV_N*AUX
        vaux%vec_list(IJCNV)%bc_dep_list(6,2) = VECTOR
        vaux%vec_list(IJCNV)%bc_dep_list(6,3) = 1 !non-trivial
      else
        vaux%vec_list(IJCNV)%bc_dep_list(2,1) = IJCNV_0*AUX
        vaux%vec_list(IJCNV)%bc_dep_list(2,2) = VECTOR
        vaux%vec_list(IJCNV)%bc_dep_list(2,3) = 0 !trivial
      endif

      !J cov
      vaux%vec_list(IJCOV)%bc_dep_list(1,1) = IJCNV*AUX
      vaux%vec_list(IJCOV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCOV)%bc_dep_list(1,3) = 0 !trivial

      !Ve cnv = v - di*j/rho
      vaux%vec_list(IVECNV)%bc_dep_list(1,1) = IRHO_AUX*AUX
      vaux%vec_list(IVECNV)%bc_dep_list(1,2) = SCALAR
      vaux%vec_list(IVECNV)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IVECNV)%bc_dep_list(2,1) = IVCOV*AUX
      vaux%vec_list(IVECNV)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IVECNV)%bc_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IVECNV)%bc_dep_list(3,1) = IJCOV*AUX
      vaux%vec_list(IVECNV)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IVECNV)%bc_dep_list(3,3) = 0 !trivial

      !Ve cov
      vaux%vec_list(IVECOV)%bc_dep_list(1,1) = IVECNV*AUX
      vaux%vec_list(IVECOV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IVECOV)%bc_dep_list(1,3) = 0 !trivial

      !Div(pi)
      vaux%vec_list(IDIVPI)%bc_dep_list(1,1) = INU*AUX
      vaux%vec_list(IDIVPI)%bc_dep_list(1,2) = SCALAR
      vaux%vec_list(IDIVPI)%bc_dep_list(1,3) = 1 !non-trivial

      vaux%vec_list(IDIVPI)%bc_dep_list(2,1) = IVCOV*AUX
      vaux%vec_list(IDIVPI)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IDIVPI)%bc_dep_list(2,3) = 1 !non-trivial

      !Div(pe)
      vaux%vec_list(IDIVPE)%bc_dep_list(1,1) = IDIVPE*AUX
      vaux%vec_list(IDIVPE)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IDIVPE)%bc_dep_list(1,3) = 0 !trivial

      !E Hall
      if (di > 0d0) then
        vaux%vec_list(IEH)%bc_dep_list(1,1) = IEH*AUX
        vaux%vec_list(IEH)%bc_dep_list(1,2) = VECTOR
        vaux%vec_list(IEH)%bc_dep_list(1,3) = 0 !trivial

        vaux%vec_list(IEH)%bc_dep_list(2,1) = IPRS
        vaux%vec_list(IEH)%bc_dep_list(2,2) = SCALAR
        vaux%vec_list(IEH)%bc_dep_list(2,3) = 1 !non-trivial

c$$$        vaux%vec_list(IEH)%bc_dep_list(3,1) = ITMP_AUX*AUX
c$$$        vaux%vec_list(IEH)%bc_dep_list(3,2) = SCALAR
c$$$        vaux%vec_list(IEH)%bc_dep_list(3,3) = 1 !non-trivial
      endif

      !Total Electric field
      vaux%vec_list(IENI)%bc_dep_list(1,1) = IENI_0*AUX
      vaux%vec_list(IENI)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IENI)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IENI)%bc_dep_list(2,1) = IETA*AUX
      vaux%vec_list(IENI)%bc_dep_list(2,2) = SCALAR
      vaux%vec_list(IENI)%bc_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IENI)%bc_dep_list(3,1) = IJCOV*AUX
      vaux%vec_list(IENI)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IENI)%bc_dep_list(3,3) = 0 !trivial

      vaux%vec_list(IENI)%bc_dep_list(4,1) = IVCOV*AUX
      vaux%vec_list(IENI)%bc_dep_list(4,2) = VECTOR
      if (binom_alph > 0d0) then
        vaux%vec_list(IENI)%bc_dep_list(4,3) = 1 !non-trivial
      else
        vaux%vec_list(IENI)%bc_dep_list(4,3) = 0 !trivial
      endif
      
#if defined(vec_pot)
      vaux%vec_list(IENI)%bc_dep_list(5,1) = IACOV*AUX
      vaux%vec_list(IENI)%bc_dep_list(5,2) = VECTOR
      vaux%vec_list(IENI)%bc_dep_list(5,3) = 1 !non-trivial
#else
      vaux%vec_list(IENI)%bc_dep_list(5,1) = IBCNV*AUX
      vaux%vec_list(IENI)%bc_dep_list(5,2) = VECTOR
      vaux%vec_list(IENI)%bc_dep_list(5,3) = 0 !trivial
#endif

      if (di > 0d0) then
        vaux%vec_list(IENI)%bc_dep_list(6,1) = IEH*AUX
        vaux%vec_list(IENI)%bc_dep_list(6,2) = VECTOR
        vaux%vec_list(IENI)%bc_dep_list(6,3) = 0 !trivial
      endif

c     End program

      end subroutine defineAuxBCDeps

c     defineAuxDomainDeps
c     ###################################################################
      subroutine defineAuxDomainDeps(vaux)

c     -------------------------------------------------------------------
c     Defines intervariable dependencies for constitutive equations of
c     auxiliary quantities.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      type(aux_array),pointer :: vaux

c     Begin program

c     Equilibrium quantities

      !Rho0
      vaux%var_list(IRHO_0)%dom_dep_list(1,1) = IRHO_0*AUX
      vaux%var_list(IRHO_0)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(IRHO_0)%dom_dep_list(1,3) = 0 !trivial

      !Temp0
      vaux%var_list(ITMPE_0)%dom_dep_list(1,1) = ITMPE_0*AUX
      vaux%var_list(ITMPE_0)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(ITMPE_0)%dom_dep_list(1,3) = 0 !trivial

      vaux%var_list(ITMPE_0)%dom_dep_list(2,1) = IRHO_0*AUX
      vaux%var_list(ITMPE_0)%dom_dep_list(2,2) = SCALAR
      vaux%var_list(ITMPE_0)%dom_dep_list(2,3) = 0 !trivial

      vaux%var_list(ITMP_0)%dom_dep_list(1,1) = ITMPE_0*AUX
      vaux%var_list(ITMP_0)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(ITMP_0)%dom_dep_list(1,3) = 0 !trivial

      !Eta0
      vaux%var_list(IETA_0)%dom_dep_list(1,1) = ITMPE_0*AUX
      vaux%var_list(IETA_0)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(IETA_0)%dom_dep_list(1,3) = 0 !trivial

#if defined(vec_pot)
      !A0 cov
      vaux%vec_list(IACOV_0)%dom_dep_list(1,1) = IACOV_0*AUX
      vaux%vec_list(IACOV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IACOV_0)%dom_dep_list(1,3) = 0 !trivial

      !B0 cnv
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,1) = IACOV_0*AUX
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,3) = 1 !non-trivial

#else
      !B0 cnv
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,1) = IBCNV_0*AUX
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,3) = 0 !trivial
#endif
      !B0 cov
      vaux%vec_list(IBCOV_0)%dom_dep_list(1,1) = IBCNV_0*AUX
      vaux%vec_list(IBCOV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCOV_0)%dom_dep_list(1,3) = 0 !trivial

      !V0 cnv
      vaux%vec_list(IVCNV_0)%dom_dep_list(1,1) = IVCNV_0*AUX
      vaux%vec_list(IVCNV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IVCNV_0)%dom_dep_list(1,3) = 0 !trivial

      !J0 cnv
      vaux%vec_list(IJCNV_0)%dom_dep_list(1,1) = IBCOV_0*AUX
      vaux%vec_list(IJCNV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCNV_0)%dom_dep_list(1,3) = 1 !non-trivial

      !J0 cov
      vaux%vec_list(IJCOV_0)%dom_dep_list(1,1) = IJCNV_0*AUX
      vaux%vec_list(IJCOV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCOV_0)%dom_dep_list(1,3) = 0 !trivial

      !V_n cnv
      vaux%vec_list(IVCNV_N)%dom_dep_list(1,1) = IVCNV*AUX
      vaux%vec_list(IVCNV_N)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IVCNV_N)%dom_dep_list(1,3) = 0 !trivial

      !Ve_n cnv
      vaux%vec_list(IVECOV_N)%dom_dep_list(1,1) = IVECOV*AUX
      vaux%vec_list(IVECOV_N)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IVECOV_N)%dom_dep_list(1,3) = 0 !trivial

#if defined(vec_pot)
      !A_n cov
      vaux%vec_list(IACOV_N)%dom_dep_list(1,1) = IACOV_N*AUX
      vaux%vec_list(IACOV_N)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IACOV_N)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IACOV_N)%dom_dep_list(2,1) = IACOV_0*AUX
      vaux%vec_list(IACOV_N)%dom_dep_list(2,2) = VECTOR
      vaux%vec_list(IACOV_N)%dom_dep_list(2,3) = 0 !trivial

      !B_n cov
      vaux%vec_list(IBCNV_N)%dom_dep_list(1,1) = IACOV_N*AUX
      vaux%vec_list(IBCNV_N)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_N)%dom_dep_list(1,3) = 1 !non-trivial     
#else
      !B_n cov
      vaux%vec_list(IBCNV_N)%dom_dep_list(1,1) = IBCNV_N*AUX
      vaux%vec_list(IBCNV_N)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_N)%dom_dep_list(1,3) = 0 !trivial     
#endif
      
      !E-field cov
      vaux%vec_list(IENI_0)%dom_dep_list(1,1) = IETA_0*AUX
      vaux%vec_list(IENI_0)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IENI_0)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IENI_0)%dom_dep_list(2,1) = IJCOV_0*AUX
      vaux%vec_list(IENI_0)%dom_dep_list(2,2) = VECTOR
      vaux%vec_list(IENI_0)%dom_dep_list(2,3) = 0 !non-trivial

      vaux%vec_list(IENI_0)%dom_dep_list(3,1) = IVCNV_0*AUX
      vaux%vec_list(IENI_0)%dom_dep_list(3,2) = VECTOR
      vaux%vec_list(IENI_0)%dom_dep_list(3,3) = 0 !trivial

c$$$#if defined(vec_pot)
c$$$      vaux%vec_list(IENI_0)%dom_dep_list(4,1) = IACOV_0*AUX
c$$$      vaux%vec_list(IENI_0)%dom_dep_list(4,2) = VECTOR
c$$$      vaux%vec_list(IENI_0)%dom_dep_list(4,3) = 1 !non-trivial
c$$$#else
      vaux%vec_list(IENI_0)%dom_dep_list(4,1) = IBCNV_0*AUX
      vaux%vec_list(IENI_0)%dom_dep_list(4,2) = VECTOR
      vaux%vec_list(IENI_0)%dom_dep_list(4,3) = 0 !trivial
c$$$#endif

      !Ve cnv = v - di*j/rho
      vaux%vec_list(IVECNV_0)%dom_dep_list(1,1) = IRHO_0*AUX
      vaux%vec_list(IVECNV_0)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IVECNV_0)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IVECNV_0)%dom_dep_list(2,1) = IVCNV_0*AUX
      vaux%vec_list(IVECNV_0)%dom_dep_list(2,2) = VECTOR
      vaux%vec_list(IVECNV_0)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IVECNV_0)%dom_dep_list(3,1) = IJCNV_0*AUX
      vaux%vec_list(IVECNV_0)%dom_dep_list(3,2) = VECTOR
      vaux%vec_list(IVECNV_0)%dom_dep_list(3,3) = 0 !trivial

c     Current variables

      !Rho
      if (adiabatic_rho) then
        vaux%var_list(IRHO_AUX)%dom_dep_list(1,1) = IRHO_0*AUX
        vaux%var_list(IRHO_AUX)%dom_dep_list(1,2) = SCALAR
        vaux%var_list(IRHO_AUX)%dom_dep_list(1,3) = 0 !trivial

        vaux%var_list(IRHO_AUX)%dom_dep_list(2,1) = ITMP_0*AUX
        vaux%var_list(IRHO_AUX)%dom_dep_list(2,2) = SCALAR
        vaux%var_list(IRHO_AUX)%dom_dep_list(2,3) = 0 !trivial

        vaux%var_list(IRHO_AUX)%dom_dep_list(3,1) = IPRS
        vaux%var_list(IRHO_AUX)%dom_dep_list(3,2) = SCALAR
        vaux%var_list(IRHO_AUX)%dom_dep_list(3,3) = 0 !trivial
      else
        vaux%var_list(IRHO_AUX)%dom_dep_list(1,1) = IRHO
        vaux%var_list(IRHO_AUX)%dom_dep_list(1,2) = SCALAR
        vaux%var_list(IRHO_AUX)%dom_dep_list(1,3) = 0 !trivial
      endif

      !Temp
      vaux%var_list(ITMP_AUX)%dom_dep_list(1,1) = ITMPE_AUX*AUX
      vaux%var_list(ITMP_AUX)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(ITMP_AUX)%dom_dep_list(1,3) = 0 !trivial

      !Electron temp
      vaux%var_list(ITMPE_AUX)%dom_dep_list(1,1) = IPRS
      vaux%var_list(ITMPE_AUX)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(ITMPE_AUX)%dom_dep_list(1,3) = 0 !trivial

      vaux%var_list(ITMPE_AUX)%dom_dep_list(2,1) = IRHO_AUX*AUX
      vaux%var_list(ITMPE_AUX)%dom_dep_list(2,2) = SCALAR
      vaux%var_list(ITMPE_AUX)%dom_dep_list(2,3) = 0 !trivial

      vaux%var_list(ITMPE_AUX)%dom_dep_list(3,1) = ITMPE_0*AUX
      vaux%var_list(ITMPE_AUX)%dom_dep_list(3,2) = SCALAR
      vaux%var_list(ITMPE_AUX)%dom_dep_list(3,3) = 0 !trivial

      !Eta
      vaux%var_list(IETA)%dom_dep_list(1,1) = ITMPE_AUX*AUX
      vaux%var_list(IETA)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(IETA)%dom_dep_list(1,3) = 0 !trivial

      !Diff
      vaux%var_list(IDFF)%dom_dep_list(1,1) = ITMPE_0*AUX
      vaux%var_list(IDFF)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(IDFF)%dom_dep_list(1,3) = 0 !trivial

      vaux%var_list(IDFF)%dom_dep_list(2,1) = ITMPE_AUX*AUX
      vaux%var_list(IDFF)%dom_dep_list(2,2) = SCALAR
      vaux%var_list(IDFF)%dom_dep_list(2,3) = 0 !trivial

      vaux%var_list(IDFF)%dom_dep_list(3,1) = IRHO_AUX*AUX
      vaux%var_list(IDFF)%dom_dep_list(3,2) = SCALAR
      vaux%var_list(IDFF)%dom_dep_list(3,3) = 0 !trivial

      vaux%var_list(IDFF)%dom_dep_list(4,1) = IBCNV*AUX
      vaux%var_list(IDFF)%dom_dep_list(4,2) = VECTOR
      vaux%var_list(IDFF)%dom_dep_list(4,3) = 0 !trivial

      !Chi
      vaux%var_list(ICHI)%dom_dep_list(1,1) = ITMPE_0*AUX
      vaux%var_list(ICHI)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(ICHI)%dom_dep_list(1,3) = 0 !trivial

      vaux%var_list(ICHI)%dom_dep_list(2,1) = ITMPE_AUX*AUX
      vaux%var_list(ICHI)%dom_dep_list(2,2) = SCALAR
      vaux%var_list(ICHI)%dom_dep_list(2,3) = 0 !trivial

      vaux%var_list(ICHI)%dom_dep_list(3,1) = IRHO_AUX*AUX
      vaux%var_list(ICHI)%dom_dep_list(3,2) = SCALAR
      vaux%var_list(ICHI)%dom_dep_list(3,3) = 0 !trivial

      vaux%var_list(ICHI)%dom_dep_list(4,1) = IBCNV*AUX
      vaux%var_list(ICHI)%dom_dep_list(4,2) = VECTOR
      vaux%var_list(ICHI)%dom_dep_list(4,3) = 0 !trivial

      !Chi
      vaux%var_list(ICHIPAR)%dom_dep_list(1,1) = ICHI*AUX
      vaux%var_list(ICHIPAR)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(ICHIPAR)%dom_dep_list(1,3) = 0 !trivial

#if defined(vec_pot)
      !A cov
      vaux%vec_list(IACOV)%dom_dep_list(1,1) = IAX
      vaux%vec_list(IACOV)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IACOV)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IACOV)%dom_dep_list(2,1) = IAY
      vaux%vec_list(IACOV)%dom_dep_list(2,2) = SCALAR
      vaux%vec_list(IACOV)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IACOV)%dom_dep_list(3,1) = IAZ
      vaux%vec_list(IACOV)%dom_dep_list(3,2) = SCALAR
      vaux%vec_list(IACOV)%dom_dep_list(3,3) = 0 !trivial

      !A cnv
      vaux%vec_list(IACNV)%dom_dep_list(1,1) = IACOV*AUX
      vaux%vec_list(IACNV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IACNV)%dom_dep_list(1,3) = 0 !trivial

      !B cnv
      vaux%vec_list(IBCNV)%dom_dep_list(1,1) = IACOV*AUX
      vaux%vec_list(IBCNV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV)%dom_dep_list(1,3) = 1 !non-trivial
#else
      !B cnv
      vaux%vec_list(IBCNV)%dom_dep_list(1,1) = IBX
      vaux%vec_list(IBCNV)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IBCNV)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IBCNV)%dom_dep_list(2,1) = IBY
      vaux%vec_list(IBCNV)%dom_dep_list(2,2) = SCALAR
      vaux%vec_list(IBCNV)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IBCNV)%dom_dep_list(3,1) = IBZ
      vaux%vec_list(IBCNV)%dom_dep_list(3,2) = SCALAR
      vaux%vec_list(IBCNV)%dom_dep_list(3,3) = 0 !trivial
#endif
      !B cov
      vaux%vec_list(IBCOV)%dom_dep_list(1,1) = IBCNV*AUX
      vaux%vec_list(IBCOV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCOV)%dom_dep_list(1,3) = 0 !trivial

      !Nu=rho*nu_i(T)
      vaux%var_list(INU)%dom_dep_list(1,1) = IRHO_AUX*AUX
      vaux%var_list(INU)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(INU)%dom_dep_list(1,3) = 0 !trivial

      vaux%var_list(INU)%dom_dep_list(2,1) = ITMP_AUX*AUX
      vaux%var_list(INU)%dom_dep_list(2,2) = SCALAR
      vaux%var_list(INU)%dom_dep_list(2,3) = 0 !trivial

      !Heta=rho*nu_e(T)
      vaux%var_list(IHETA)%dom_dep_list(1,1) = IRHO_0*AUX
      vaux%var_list(IHETA)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(IHETA)%dom_dep_list(1,3) = 0 !trivial

      vaux%var_list(IHETA)%dom_dep_list(2,1) = IBCNV_0*AUX
      vaux%var_list(IHETA)%dom_dep_list(2,2) = VECTOR
      vaux%var_list(IHETA)%dom_dep_list(2,3) = 0 !trivial

      !V cnv
      vaux%vec_list(IVCNV)%dom_dep_list(1,1) = IVX
      vaux%vec_list(IVCNV)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IVCNV)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IVCNV)%dom_dep_list(2,1) = IVY
      vaux%vec_list(IVCNV)%dom_dep_list(2,2) = SCALAR
      vaux%vec_list(IVCNV)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IVCNV)%dom_dep_list(3,1) = IVZ
      vaux%vec_list(IVCNV)%dom_dep_list(3,2) = SCALAR
      vaux%vec_list(IVCNV)%dom_dep_list(3,3) = 0 !trivial

      !V cov
      vaux%vec_list(IVCOV)%dom_dep_list(1,1) = IVCNV*AUX
      vaux%vec_list(IVCOV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IVCOV)%dom_dep_list(1,3) = 0 !trivial

      !J cnv
      vaux%vec_list(IJCNV)%dom_dep_list(1,1) = IBCOV*AUX
      vaux%vec_list(IJCNV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCNV)%dom_dep_list(1,3) = 1 !non-trivial

      !J cov
      vaux%vec_list(IJCOV)%dom_dep_list(1,1) = IJCNV*AUX
      vaux%vec_list(IJCOV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCOV)%dom_dep_list(1,3) = 0 !trivial

      !Ve cnv = v - di*j/rho
      vaux%vec_list(IVECNV)%dom_dep_list(1,1) = IRHO_AUX*AUX
      vaux%vec_list(IVECNV)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IVECNV)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IVECNV)%dom_dep_list(2,1) = IVCOV*AUX
      vaux%vec_list(IVECNV)%dom_dep_list(2,2) = VECTOR
      vaux%vec_list(IVECNV)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IVECNV)%dom_dep_list(3,1) = IJCOV*AUX
      vaux%vec_list(IVECNV)%dom_dep_list(3,2) = VECTOR
      vaux%vec_list(IVECNV)%dom_dep_list(3,3) = 0 !trivial

      !Ve cov
      vaux%vec_list(IVECOV)%dom_dep_list(1,1) = IVECNV*AUX
      vaux%vec_list(IVECOV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IVECOV)%dom_dep_list(1,3) = 0 !trivial

      !Div(pi)
      vaux%vec_list(IDIVPI)%dom_dep_list(1,1) = INU*AUX
      vaux%vec_list(IDIVPI)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IDIVPI)%dom_dep_list(1,3) = 1 !non-trivial

      vaux%vec_list(IDIVPI)%dom_dep_list(2,1) = IVCNV*AUX
      vaux%vec_list(IDIVPI)%dom_dep_list(2,2) = VECTOR
      vaux%vec_list(IDIVPI)%dom_dep_list(2,3) = 1 !non-trivial

      !Div(pe)
      vaux%vec_list(IDIVPE)%dom_dep_list(1,1) = IHETA*AUX
      vaux%vec_list(IDIVPE)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IDIVPE)%dom_dep_list(1,3) = 1 !non-trivial

      vaux%vec_list(IDIVPE)%dom_dep_list(2,1) = IVECNV*AUX
      vaux%vec_list(IDIVPE)%dom_dep_list(2,2) = VECTOR
      vaux%vec_list(IDIVPE)%dom_dep_list(2,3) = 1 !non-trivial

      vaux%vec_list(IDIVPE)%dom_dep_list(3,1) = IVECNV_0*AUX
      vaux%vec_list(IDIVPE)%dom_dep_list(3,2) = VECTOR
      vaux%vec_list(IDIVPE)%dom_dep_list(3,3) = 1 !non-trivial

      !E Hall
      if (di > 0d0) then
        if (ion_hall) then
          vaux%vec_list(IEH)%dom_dep_list(1,1) = IVCNV_N*AUX
          vaux%vec_list(IEH)%dom_dep_list(1,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(1,3) = 0 !trivial

          vaux%vec_list(IEH)%dom_dep_list(2,1) = IVCNV*AUX
          vaux%vec_list(IEH)%dom_dep_list(2,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(2,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(3,1) = IRHO_AUX*AUX
          vaux%vec_list(IEH)%dom_dep_list(3,2) = SCALAR
          vaux%vec_list(IEH)%dom_dep_list(3,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(4,1) = ITMPE_AUX*AUX
          vaux%vec_list(IEH)%dom_dep_list(4,2) = SCALAR
          vaux%vec_list(IEH)%dom_dep_list(4,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(5,1) = IDIVPI*AUX
          vaux%vec_list(IEH)%dom_dep_list(5,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(5,3) = 0 !trivial
        else
          vaux%vec_list(IEH)%dom_dep_list(1,1) = IRHO_AUX*AUX
          vaux%vec_list(IEH)%dom_dep_list(1,2) = SCALAR
          vaux%vec_list(IEH)%dom_dep_list(1,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(2,1) = ITMPE_AUX*AUX
          vaux%vec_list(IEH)%dom_dep_list(2,2) = SCALAR
          vaux%vec_list(IEH)%dom_dep_list(2,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(3,1) = IBCOV*AUX
          vaux%vec_list(IEH)%dom_dep_list(3,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(3,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(4,1) = IJCOV*AUX
          vaux%vec_list(IEH)%dom_dep_list(4,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(4,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(5,1) = IDIVPE*AUX
          vaux%vec_list(IEH)%dom_dep_list(5,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(5,3) = 1 !non-trivial

          if (de > 0d0) then
            vaux%vec_list(IEH)%dom_dep_list(6,1) = IVECOV*AUX
            vaux%vec_list(IEH)%dom_dep_list(6,2) = VECTOR
            vaux%vec_list(IEH)%dom_dep_list(6,3) = 1 !non-trivial

            vaux%vec_list(IEH)%dom_dep_list(7,1) = IVECOV_N*AUX
            vaux%vec_list(IEH)%dom_dep_list(7,2) = VECTOR
            vaux%vec_list(IEH)%dom_dep_list(7,3) = 1 !non-trivial
          endif
        endif
      endif

      !E total
      vaux%vec_list(IENI)%dom_dep_list(1,1) = IENI_0*AUX
      vaux%vec_list(IENI)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IENI)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IENI)%dom_dep_list(2,1) = IETA*AUX
      vaux%vec_list(IENI)%dom_dep_list(2,2) = SCALAR
      vaux%vec_list(IENI)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IENI)%dom_dep_list(3,1) = IJCOV*AUX
      vaux%vec_list(IENI)%dom_dep_list(3,2) = VECTOR
      vaux%vec_list(IENI)%dom_dep_list(3,3) = 0 !trivial

      vaux%vec_list(IENI)%dom_dep_list(4,1) = IVCOV*AUX
      vaux%vec_list(IENI)%dom_dep_list(4,2) = VECTOR
      if (binom_alph > 0d0) then
        vaux%vec_list(IENI)%dom_dep_list(4,3) = 1 !non-trivial
      else
        vaux%vec_list(IENI)%dom_dep_list(4,3) = 0 !trivial
      endif
      
#if defined(vec_pot)
      vaux%vec_list(IENI)%dom_dep_list(5,1) = IACOV*AUX
      vaux%vec_list(IENI)%dom_dep_list(5,2) = VECTOR
      vaux%vec_list(IENI)%dom_dep_list(5,3) = 1 !non-trivial
#else
      vaux%vec_list(IENI)%dom_dep_list(5,1) = IBCNV*AUX
      vaux%vec_list(IENI)%dom_dep_list(5,2) = VECTOR
      vaux%vec_list(IENI)%dom_dep_list(5,3) = 0 !trivial
#endif

      if (di > 0d0) then
        vaux%vec_list(IENI)%dom_dep_list(6,1) = IEH*AUX
        vaux%vec_list(IENI)%dom_dep_list(6,2) = VECTOR
        vaux%vec_list(IENI)%dom_dep_list(6,3) = 0 !trivial
      endif

c     End program

      end subroutine defineAuxDomainDeps

c     ExB_flow_bc
c     #################################################################
      subroutine ExB_flow_bc(g_def,igrid,bcnv,prs,eta,E_0,v0,dbcov)

c     -----------------------------------------------------------------
c     Finds ExB flow boundary condition when Eo is imposed at boundary
c     (returns contravariant components).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8),dimension(0:,0:,0:)   :: prs,eta
      real(8),dimension(0:,0:,0:,:) :: bcnv,v0,E_0
      real(8),dimension(0:,0:,0:,:),optional :: dbcov

c     Local variables

      integer :: i,j,k,ig,jg,kg,ord,nx,ny,nz,ibc

      real(8) :: b2,bbcnv(3),bbcov(3),dum(3),grad_p(3),etah,gsup(3,3)
     .          ,gsub(3,3),jac,E00(3)

c     Begin program

      nx = size(bcnv,1)-2
      ny = size(bcnv,2)-2
      nz = size(bcnv,3)-2

      i = nx ; ibc = 2

      if (.not.isBdry(g_def,i,igrid,ibc)) return

      if (igrid /= 1) then
        write (*,*) 'ExB_flow_bc not ready for MG treatment'
        stop
      endif

      !Compute value at face
      do k=0,nz+1
        do j=0,ny+1
          jac = 0.5*(g_def%gmetric%grid(igrid)%jac(i  ,j,k)
     .              +g_def%gmetric%grid(igrid)%jac(i+1,j,k))

          bbcnv = 0.5*(bcnv(i,j,k,:)+bcnv(i+1,j,k,:))
          if (.not.rw_is_vacuum) bbcnv(1) = 0d0 !Zero-out normal B-field
          bbcov = XformToCov_ijk_face(g_def,i,j,k,igrid,bbcnv,1)
          
          b2 = sum(bbcnv*bbcov)/jac

          if (b2 > 0d0) then
            if (PRESENT(dbcov)) then  !Linearize ExB flow BC
              dum = 0.5*(dbcov(i,j,k,:)+dbcov(i+1,j,k,:))
              bbcov = dum - 2*dot_product(bbcnv,dum)/(jac*b2)*bbcov
            endif

            grad_p = grad(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid,prs
     .                   ,hex=1)
            grad_p = XformToCnv_ijk_face(g_def,i,j,k,igrid,grad_p,1)

            etah = 0.5*(eta(i,j,k)+eta(i+1,j,k))

            !Tangential eta*j at wall
            E00 = E0 + E_0(i+1,j,k,:) !E_0(nx+1) contains face value

            !ExB velocity
            v0(i+1,j,k,1)=(E00(2)*bbcov(3)-E00(3)*bbcov(2)
     .                   -etah*grad_p(1))/b2
            v0(i+1,j,k,2)=(E00(3)*bbcov(1)-E00(1)*bbcov(3)
     .                   -etah*grad_p(2))/b2
            v0(i+1,j,k,3)=(E00(1)*bbcov(2)-E00(2)*bbcov(1)
     .                   -etah*grad_p(3))/b2
            
            if (.not.rw_is_vacuum) v0(i+1,j,k,2:3) = 0d0 !Zero-out tangential flow
          endif
        enddo
      enddo

c     End program

      end subroutine ExB_flow_bc

c     J_pinch_bc
c     #################################################################
      subroutine J_pinch_bc(g_def,igrid,bcnv,prs,eta,E_0,j0,dbcnv)

c     -----------------------------------------------------------------
c     Finds current boundary condition when Eo is imposed at boundary
c     (returns covariant current).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8),dimension(0:,0:,0:)   :: prs,eta
      real(8),dimension(0:,0:,0:,:) :: bcnv,j0,E_0
      real(8),dimension(0:,0:,0:,:),optional :: dbcnv

c     Local variables

      integer :: i,j,k,ig,jg,kg,ord,nx,ny,nz,ibc
      real(8) :: Jb2,b2,bbcnv(3),bbcov(3),grad_p(3),alpha
     .          ,gsub(3,3),ietah,jach,E00(3),xh(1)

c     Begin program

      nx = size(bcnv,1)-2
      ny = size(bcnv,2)-2
      nz = size(bcnv,3)-2

      if (igrid /= 1) then
        write (*,*) 'J_pinch_bc not ready for MG treatment'
        stop
      endif

      i = nx ; ord = 1 ; ibc = 2

      if (.not.isBdry(g_def,i,igrid,ibc)) return

      call getMGmap(g_def,i,1,1,igrid,igrid,igrid,ig,jg,kg)

      do k=1,nz
        do j=1,ny
          !First-order extrapolation of B to boundary, arbitrary mesh (critical)
c$$$          bbcnv = bcnv(i-1,j,k,:) + (bcnv(i,j,k,:)-bcnv(i-1,j,k,:))
c$$$     $                             /(g_def%xx(ig) -g_def%xx(ig-1))
c$$$     $                               *(0.5*g_def%xx(ig+1)
c$$$     $                                +0.5*g_def%xx(ig)
c$$$     $                                -    g_def%xx(ig-1))
          xh = 0.5*(g_def%xx(ig+1)+g_def%xx(ig))
          call IntDriver1d(4,g_def%xx(ig-3:ig),bcnv(i-3:i,j,k,1)
     .                    ,1,xh,bbcnv(1),ord,0)
          call IntDriver1d(4,g_def%xx(ig-3:ig),bcnv(i-3:i,j,k,2)
     .                    ,1,xh,bbcnv(2),ord,0)
          call IntDriver1d(4,g_def%xx(ig-3:ig),bcnv(i-3:i,j,k,3)
     .                    ,1,xh,bbcnv(3),ord,0)

          if (.not.rw_is_vacuum) bbcnv(1) = 0d0 !Zero-out normal B-field
          bbcov = XformToCov_ijk_face(g_def,i,j,k,igrid,bbcnv,1)
c$$$          gsub = 0.5*(g_def%gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
c$$$     .               +g_def%gmetric%grid(igrid)%gsub(i+1,j,k,:,:))
c$$$          bbcov = matmul(gsub,bbcnv)

          jach = 0.5*(g_def%gmetric%grid(igrid)%jac(i  ,j,k)
     .               +g_def%gmetric%grid(igrid)%jac(i+1,j,k))

          Jb2 = dot_product(bbcov,bbcnv)
          b2  = Jb2/jach

          if (b2 > 0d0) then
c$$$            ietah = 2./(eta(i+1,j,k)+eta(i,j,k))
            ietah =0.5*(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

            grad_p
     .         = grad(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,hex=1)

            E00 = E0 + E_0(i+1,j,k,:) !Tangential current at wall
            alpha = dot_product(E00,bbcnv)*ietah/Jb2

            if (PRESENT(dbcnv)) bbcnv = dbcnv(i,j,k,:)

            j0(i+1,j,k,:) = alpha*bbcnv + cross_product(bbcov,grad_p)/b2
          endif
        enddo
      enddo

c     End program

      end subroutine J_pinch_bc

c     init_BCs
c     ####################################################################
      subroutine init_BCs(bcsq)
c     --------------------------------------------------------------------
c     Defines default BCs for SCALARS (rho, T)
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: bcsq(:,:)

c     Local variables

      integer :: i

c     Begin program

      do i=1,size(bcsq,2)
        bcsq(:,i) = bcond
      enddo

c     End program

      end subroutine init_BCs

c     default_SCL_BCs
c     ####################################################################
      subroutine default_SCL_BCs(bcsq)
c     --------------------------------------------------------------------
c     Defines default BCs for SCALARS (rho, T)
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: bcsq(6)

c     Local variables

c     Begin program

      where (bcsq == DEF .or. bcsq == FSYM) bcsq = NEU

c     End program

      end subroutine default_SCL_BCs

c     default_B_BCs
c     ####################################################################
      subroutine default_B_BCs(bbcs)
c     --------------------------------------------------------------------
c     Defines default BCs for magnetic field
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: bbcs(6,3)

c     Local variables

      integer :: bcsq(6)

c     Begin program

      call init_BCs(bbcs)

      bcsq = bbcs(:,1)
      where (bcsq(1:2) == DEF) bcsq(1:2) =  DIR
      where (bcsq(3:6) == DEF) bcsq(3:6) = -NEU !On covariant components
c$$$      where (bcsq(3:6) == DEF.and.(.not.rw_bcs(3:6))) bcsq(3:6) = -NEU !On covariant components
c$$$      where (bcsq(3:6) == DEF.and.(     rw_bcs(3:6))) bcsq(3:6) =  EQU !On contravariant components
      bbcs(:,1) = bcsq

      bcsq = bbcs(:,2)
      where (bcsq(3:4) == DEF) bcsq(3:4) =  DIR
      where (bcsq(1:2) == DEF) bcsq(1:2) = -NEU  !On covariant components
      where (bcsq(5:6) == DEF) bcsq(5:6) = -NEU  !On covariant components
c$$$      where (bcsq(1:2) == DEF.and.(.not.rw_bcs(1:2))) bcsq(1:2) = -NEU  !On covariant components
c$$$      where (bcsq(1:2) == DEF.and.(     rw_bcs(1:2))) bcsq(1:2) =  EQU  !On contravariant components
c$$$      where (bcsq(5:6) == DEF.and.(.not.rw_bcs(5:6))) bcsq(5:6) = -NEU  !On covariant components
c$$$      where (bcsq(5:6) == DEF.and.(     rw_bcs(5:6))) bcsq(5:6) =  EQU  !On contravariant components
      bbcs(:,2) = bcsq

      bcsq = bbcs(:,3)
      where (bcsq(1:4) == DEF) bcsq(1:4) = -NEU !On covariant components
c$$$      where (bcsq(1:4) == DEF.and.(.not.rw_bcs(1:4))) bcsq(1:4) = -NEU !On covariant components
c$$$      where (bcsq(1:4) == DEF.and.(     rw_bcs(1:4))) bcsq(1:4) =  EQU !On contravariant components
      where (bcsq(5:6) == DEF) bcsq(5:6) =  DIR
      bbcs(:,3) = bcsq

      call set_vec_symm_bc(bbcs,1)

c     End program

      end subroutine default_B_BCs

c     default_A_BCs
c     ####################################################################
      subroutine default_A_BCs(bbcs)
c     --------------------------------------------------------------------
c     Defines default BCs for vector potential
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: bbcs(6,3)

c     Local variables

      integer :: bcsq(6)

c     Begin program

      call init_BCs(bbcs)

      bcsq = bbcs(:,1)
      where (bcsq(1:2) == DEF) bcsq(1:2) =  DIR !On contravariant components
      where (bcsq(3:6) == DEF) bcsq(3:6) = -EQU !On covariant components
      bbcs(:,1) = bcsq

      bcsq = bbcs(:,2)
      where (bcsq(3:4) == DEF) bcsq(3:4) =  DIR
      where (bcsq(1:2) == DEF) bcsq(1:2) = -EQU !On covariant components
      where (bcsq(5:6) == DEF) bcsq(5:6) = -EQU !On covariant components
      bbcs(:,2) = bcsq

      bcsq = bbcs(:,3)
      where (bcsq(1:4) == DEF) bcsq(1:4) = -EQU !On covariant components
      where (bcsq(5:6) == DEF) bcsq(5:6) =  DIR
      bbcs(:,3) = bcsq

      call set_vec_symm_bc(bbcs,-1)

c     End program

      end subroutine default_A_BCs

c     default_V_BCs
c     ####################################################################
      subroutine default_V_BCs(bbcs)
c     --------------------------------------------------------------------
c     Defines default BCs for velocity
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: bbcs(6,3)

c     Local variables

      integer :: bcsq(6)

c     Begin program

      call init_BCs(bbcs)

c$$$      if (di > 0d0) then        !Force D_t(v x n) = 0 for energy conservation
c$$$        bcsq = bbcs(:,1)
c$$$        where (bcsq(1:2) == DEF) bcsq(1:2) = DIR
c$$$        where (bcsq(3:6) == DEF) bcsq(3:6) =-EQU !On covariant components
c$$$        bbcs(:,1) = bcsq
c$$$
c$$$        bcsq = bbcs(:,2)
c$$$        if (     coords == 'cyl'
c$$$     .      .or. coords == 'hel'
c$$$     .      .or. coords == 'tor'
c$$$     .      .or. coords == 'car') then
c$$$          where (bcsq(1:2) == DEF) bcsq(1:2) =  EQU !On contravariant components
c$$$        else
c$$$          where (bcsq(1:2) == DEF) bcsq(1:2) = -EQU !On covariant components
c$$$        endif
c$$$        where (bcsq(3:4) == DEF) bcsq(3:4) = DIR
c$$$        where (bcsq(5:6) == DEF) bcsq(5:6) =-EQU
c$$$        bbcs(:,2) = bcsq
c$$$
c$$$        bcsq = bbcs(:,3)
c$$$        if (     coords == 'cyl'
c$$$     .      .or. coords == 'hel'
c$$$     .      .or. coords == 'tor'
c$$$     .      .or. coords == 'car') then
c$$$          where (bcsq(1:2) == DEF) bcsq(1:2) =  EQU !On contravariant components
c$$$        else
c$$$          where (bcsq(1:2) == DEF) bcsq(1:2) = -EQU !On covariant components
c$$$        endif
c$$$        where (bcsq(3:4) == DEF) bcsq(3:4) =-EQU
c$$$        where (bcsq(5:6) == DEF) bcsq(5:6) = DIR
c$$$        bbcs(:,3) = bcsq
c$$$      else

        bcsq = bbcs(:,1)
        where (bcsq(1:2) == DEF) bcsq(1:2) = DIR !Normal
        where (bcsq(3:6) == DEF) bcsq(3:6) =-NEU !Tangential, on covariant components
        bbcs(:,1) = bcsq

        bcsq = bbcs(:,2)
        if (     coords == 'cyl'
     .      .or. coords == 'hel'
     .      .or. coords == 'tor'
     .      .or. coords == 'car') then
          where (bcsq(1:2) == DEF) bcsq(1:2) =  NEU !Tangential, on contravariant components
        else
          where (bcsq(1:2) == DEF) bcsq(1:2) = -NEU !Tangential, on covariant components
        endif
        where (bcsq(3:4) == DEF) bcsq(3:4) = DIR !Normal
        where (bcsq(5:6) == DEF) bcsq(5:6) =-NEU !Tangential, on covariant components
        bbcs(:,2) = bcsq

        bcsq = bbcs(:,3)
        if (     coords == 'cyl'
     .      .or. coords == 'hel'
     .      .or. coords == 'tor'
     .      .or. coords == 'car') then
          where (bcsq(1:2) == DEF) bcsq(1:2) =  NEU !Tangential, on contravariant components
        else
          where (bcsq(1:2) == DEF) bcsq(1:2) = -NEU !Tangential, on covariant components
        endif
        where (bcsq(3:4) == DEF) bcsq(3:4) =-NEU !Tangential, on covariant components
        where (bcsq(5:6) == DEF) bcsq(5:6) = DIR !Normal
        bbcs(:,3) = bcsq
c$$$      endif

      call set_vec_symm_bc(bbcs,-1)

c     End program

      end subroutine default_V_BCs

      end module local_BCS_variables

