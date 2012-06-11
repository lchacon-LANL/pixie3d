c module local_BCS_variables
c####################################################################
      module local_BCS_variables

        use variables

        use imposeBCinterface, ONLY: extrapolateBC

        use auxiliaryVariables

        use operators

        use timestepping

        use equilibrium

        use nlfunction_setup, ONLY: res,vis,E_fld_res,E_fld_hall
     .                             ,ion_hall,fake_ve
     .                             ,EOM_divPi,EOM_divPe,hres
     .                             ,nc_eom_v,nc_eom_jxb
     .                             ,save_tn,slava,pinch_flow
     .                             ,res_ijk,vis_ijk

c       Hardwired BC configuration parameters

        !Order of BC interpolation
#if defined(vec_pot)
        integer :: bc_order=1
#else
        integer :: bc_order=1
#endif
      contains

c     setAppBCs
c     ####################################################################
      subroutine setAppBCs(varray,vaux,it)
c     --------------------------------------------------------------------
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
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: it
      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

c     Local variables

c     Begin program

c     Setup dependent variables

      call defineDepVariables(varray)

c     Setup auxiliary variables

      call defineAuxVariables(vaux)

c     Set application BCs

      call setBCsched(varray,vaux,it)

c     End program

      end subroutine setAppBCs

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

      pinch_flow = (    equil == 'ppnch' .or. equil == 'ppn3d'
     .              .or.equil == 'ppnsl' .or. equil == 'ppnst'
     .              .or.equil == '2dhel' .or. equil == 'p3nsl'
     .              .or.equil == 'ohtor')

c     Assign BCs

      do ieq=1,neqd
        bbcs(:,ieq) = varray%array_var(ieq)%bconds
      enddo

c     Default boundary conditions

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

      call set_vec_symm_bc(bbcs(:,IVX:IVZ),-1)

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

      call set_vec_symm_bc(bbcs(:,IBX:IBZ),1)
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

      call set_vec_symm_bc(bbcs(:,IAX:IAZ),-1)
#endif

      bcsq = bbcs(:,ITMP)
      where (bcsq == DEF .or. bcsq == FSYM) bcsq = NEU
cc      if (gamma /= 1d0 .and. (.not.solve_prs)) then
cc        where (bcsq == DEF) bcsq = DIR
cc      else
cc        where (bcsq == DEF) bcsq = NEU !To allow isothermal case
cc      endif
      bbcs(:,ITMP) = bcsq

c     Exceptions for specific equilibria

      select case (equil)
      case ('tmcyz','3dkh')

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = -NEU
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) =  DIR

        where (bbcs(:,IBX) == DIR) bbcs(:,IBX) = -NEU
        where (bbcs(:,IBY) ==-NEU) bbcs(:,IBY) =  DIR

      case ('ppnch','ppn3d','ppnst','2dhel')   !No stress BCs, pinch flow

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX ) = EQU

cc        if (gamma > 1d0)
cc     .       where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = EQU

      case ('ppnsl','p3nsl','ohtor')           !No slip BCs, pinch flow

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX ) = EQU

        !Dirichlet boundary conditions for tangential components
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY ) = DIR
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ ) = DIR

        if (gamma > 1d0)
     .       where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = DIR

      case ('vmec','tor')

        !Dirichlet boundary conditions for tangential components
        where (abs(bbcs(:,IVY)) == NEU) bbcs(:,IVY) = DIR
        where (abs(bbcs(:,IVZ)) == NEU) bbcs(:,IVZ) = DIR

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

cc        !Dirichlet boundary conditions for tangential velocity components
cc        where (abs(bbcs(:,IVY )) == NEU) bbcs(:,IVY ) = DIR
cc        where (abs(bbcs(:,IVZ )) == NEU) bbcs(:,IVZ ) = DIR

        where (abs(bbcs(:,ITMP)) == DIR) bbcs(:,ITMP) = NEU

      case('adv2')

        bbcs(1,IRHO) = EQU
        bbcs(2,IRHO) = EXT
        bbcs(3,IRHO) = EQU
        bbcs(4,IRHO) = EXT
        
        bbcs(1,ITMP) = EQU
        bbcs(2,ITMP) = EXT
        bbcs(3,ITMP) = EQU
        bbcs(4,ITMP) = EXT

        where (bcond == DEF) bbcs(:,IVX) = NEU
        where (bcond == DEF) bbcs(:,IVY) = NEU
 
cc      case('sod')
cc
cc        bbcs(1,IRHO) = NEU
cc        bbcs(2,IRHO) = NEU
cc        
cc        bbcs(1,ITMP) = NEU
cc        bbcs(2,ITMP) = NEU
cc
cc        bbcs(1,IVX)  = DIR
cc        bbcs(2,IVX)  = DIR

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
      varray%array_var(ITMP)%bc_dep_list(1,1) = ITMP_AUX*AUX
      varray%array_var(ITMP)%bc_dep_list(1,2) = SCALAR
      varray%array_var(ITMP)%bc_dep_list(1,3) = 0 !trivial

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

      !Tmp
      varray%array_var(ITMP)%dom_dep_list(1,1) = ITMP_AUX*AUX
      varray%array_var(ITMP)%dom_dep_list(1,2) = SCALAR
      varray%array_var(ITMP)%dom_dep_list(1,3) = 1 !non-trivial

      varray%array_var(ITMP)%dom_dep_list(2,1) = IVCNV*AUX
      varray%array_var(ITMP)%dom_dep_list(2,2) = VECTOR
      varray%array_var(ITMP)%dom_dep_list(2,3) = 1 !non-trivial

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
      !B
      do ieq=IAX,IAZ
        varray%array_var(ieq)%dom_dep_list(1,1) = IENI*AUX
        varray%array_var(ieq)%dom_dep_list(1,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(1,3) = 1 !non-trivial

        varray%array_var(ieq)%dom_dep_list(2,1) = IENI_0*AUX
        varray%array_var(ieq)%dom_dep_list(2,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(2,3) = 1 !non-trivial
      enddo
#endif

c     End subroutine

      end subroutine defineDepDom_Deps

cc#if defined(petsc)
#if 1
c     setBCsched
c     ###################################################################
      subroutine setBCsched(varray,vaux,itime)

c     -------------------------------------------------------------------
c     Sets BC schedule for dependent and auxiliary variables, which
c     are identified by a 2-dim integer array, bc_seq. The first
c     dimension is of size nbc_seq, and indicates order of treatment
c     lexicographically (i.e., first component of first dimension is
c     to be treated first, second component of first dimension is
c     to be treated second, and so on).
c
c     The second dimension is of size two, and has the following convention:
c       * First component of the second dimension identifies variable,
c         and its sign indicates whether the variable is dependent (>0) or
c         auxiliary (<0).
c       * Second component of the second dimension identifies whether
c         variable is a scalar (=0) or a vector (=1).
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: itime

      type(var_array),pointer :: varray

      type(aux_array),pointer :: vaux

c     Local variables

      integer :: ierr,offset=0,igr

c     Begin program

c     Select offset based on whether we are at the beginning of calculation

      if (itime == 0) then
#if defined(vec_pot)
        offset = 8
#else
        offset = 7
#endif
      elseif (itime > 0) then

        offset = 0

      endif

c     Set number of BC groups

#if defined(vec_pot)
      gv%nbc_grp = 9+offset
#else
      if (fake_ve) then
        gv%nbc_grp = 9+offset
      else
        gv%nbc_grp = 8+offset
      endif
#endif

c     Check BC groups

      if (     associated(gv%bc_grp)
     .    .and.size(gv%bc_grp) /= gv%nbc_grp) then

        do igr=1,size(gv%bc_grp)
          deallocate(gv%bc_grp(igr)%bc_seq)
        enddo
        deallocate(gv%bc_grp)

      endif

c     Setup schedule

      if (.not.associated(gv%bc_grp)) then

        allocate(gv%bc_grp(gv%nbc_grp),STAT=ierr)

        !Define group sizes
        if (itime == 0) then
          gv%bc_grp(1)%nbc_seq = 1
          gv%bc_grp(2)%nbc_seq = 1
          gv%bc_grp(3)%nbc_seq = 1

          offset = 0
#if defined(vec_pot)
          gv%bc_grp(4)%nbc_seq = 1
          offset = 1
#endif
          gv%bc_grp(4+offset)%nbc_seq = 2
          gv%bc_grp(5+offset)%nbc_seq = 3
          gv%bc_grp(6+offset)%nbc_seq = 1
          gv%bc_grp(7+offset)%nbc_seq = 1
          offset = offset + 7
        endif

#if defined(vec_pot)
        gv%bc_grp(1+offset)%nbc_seq = 1
        gv%bc_grp(2+offset)%nbc_seq = 1
        gv%bc_grp(3+offset)%nbc_seq = 2
        gv%bc_grp(4+offset)%nbc_seq = 5
        gv%bc_grp(5+offset)%nbc_seq = 2
        gv%bc_grp(6+offset)%nbc_seq = 2
        gv%bc_grp(7+offset)%nbc_seq = 7
        gv%bc_grp(8+offset)%nbc_seq = 4
        gv%bc_grp(9+offset)%nbc_seq = 4
#else
        gv%bc_grp(1+offset)%nbc_seq = 1
        gv%bc_grp(2+offset)%nbc_seq = 1
        gv%bc_grp(3+offset)%nbc_seq = 2
        gv%bc_grp(4+offset)%nbc_seq = 5
        gv%bc_grp(5+offset)%nbc_seq = 2
        gv%bc_grp(6+offset)%nbc_seq = 7
        if (fake_ve) then
          gv%bc_grp(7+offset)%nbc_seq = 3
          gv%bc_grp(8+offset)%nbc_seq = 4
          gv%bc_grp(9+offset)%nbc_seq = 1
        else
          gv%bc_grp(7+offset)%nbc_seq = 3
          gv%bc_grp(8+offset)%nbc_seq = 4
        endif
#endif

        !Allocate groups
        do igr=1,gv%nbc_grp
          allocate(gv%bc_grp(igr)
     .               %bc_seq(gv%bc_grp(igr)%nbc_seq,2),STAT=ierr)

          gv%bc_grp(igr)%bc_seq = 0
        enddo

        !Set groups
        if (itime == 0) then
          !Group 1
          gv%bc_grp(1)%bc_seq(1,1) = IRHO_0*AUX          !Aux. Variable
          gv%bc_grp(1)%bc_seq(1,2) = SCALAR              !Scalar

          !Group 2   
          gv%bc_grp(2)%bc_seq(1,1) = ITMP_0*AUX          !Aux. variable
          gv%bc_grp(2)%bc_seq(1,2) = SCALAR              !Scalar

          !Group 3
          gv%bc_grp(3)%bc_seq(1,1) = IETA_0*AUX          !Aux. variable
          gv%bc_grp(3)%bc_seq(1,2) = SCALAR              !Scalar

          !Group 4
          offset = 0
#if defined(vec_pot)
          gv%bc_grp(4)%bc_seq(1,1) = IACOV_0*AUX         !Aux. variable
          gv%bc_grp(4)%bc_seq(1,2) = VECTOR              !Vector
          offset = 1
#endif
          !Group 4+
          gv%bc_grp(4+offset)%bc_seq(1,1) = IBCNV_0*AUX  !Aux. variable
          gv%bc_grp(4+offset)%bc_seq(1,2) = VECTOR       !Vector

          gv%bc_grp(4+offset)%bc_seq(2,1) = IBCOV_0*AUX  !Aux. variable
          gv%bc_grp(4+offset)%bc_seq(2,2) = VECTOR       !Vector

          !Group 5+
          gv%bc_grp(5+offset)%bc_seq(1,1) = IHETA*AUX    !Aux. variable
          gv%bc_grp(5+offset)%bc_seq(1,2) = SCALAR       !Scalar

          gv%bc_grp(5+offset)%bc_seq(2,1) = IJCNV_0*AUX  !Aux. variable
          gv%bc_grp(5+offset)%bc_seq(2,2) = VECTOR       !Vector

          gv%bc_grp(5+offset)%bc_seq(3,1) = IJCOV_0*AUX  !Aux. variable
          gv%bc_grp(5+offset)%bc_seq(3,2) = VECTOR       !Vector

          !Group 6+
          gv%bc_grp(6+offset)%bc_seq(1,1) = IVCNV_0*AUX  !Aux. variable
          gv%bc_grp(6+offset)%bc_seq(1,2) = VECTOR       !Vector

          !Group 7+
cc          gv%bc_grp(7+offset)%bc_seq(1,1) = IVECNV_0*AUX !Aux. variable
cc          gv%bc_grp(7+offset)%bc_seq(1,2) = VECTOR       !Vector

          gv%bc_grp(7+offset)%bc_seq(1,1) = IENI_0*AUX   !Aux. variable
          gv%bc_grp(7+offset)%bc_seq(1,2) = VECTOR       !Vector

          offset = offset + 7
        endif

#if defined(vec_pot)
        !Group 1
        gv%bc_grp(1+offset)%bc_seq(1,1) = IRHO_AUX*AUX!Aux. Variable
        gv%bc_grp(1+offset)%bc_seq(1,2) = SCALAR      !Scalar

        !Group 2    
        gv%bc_grp(2+offset)%bc_seq(1,1) = ITMP_AUX*AUX!Aux. variable
        gv%bc_grp(2+offset)%bc_seq(1,2) = SCALAR      !Scalar

        !Group 3
        gv%bc_grp(3+offset)%bc_seq(1,1) = IETA*AUX    !Aux. variable
        gv%bc_grp(3+offset)%bc_seq(1,2) = SCALAR      !Scalar

        gv%bc_grp(3+offset)%bc_seq(2,1) = INU*AUX    !Aux. variable
        gv%bc_grp(3+offset)%bc_seq(2,2) = SCALAR      !Scalar

        !Group 4
        gv%bc_grp(4+offset)%bc_seq(1,1) = IACOV*AUX   !Aux. variable
        gv%bc_grp(4+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(4+offset)%bc_seq(2,1) = IACNV*AUX   !Aux. variable
        gv%bc_grp(4+offset)%bc_seq(2,2) = VECTOR      !Vector

        gv%bc_grp(4+offset)%bc_seq(3,1) = IAX         !Dep. variable
        gv%bc_grp(4+offset)%bc_seq(3,2) = SCALAR      !Vector

        gv%bc_grp(4+offset)%bc_seq(4,1) = IAY         !Dep. variable
        gv%bc_grp(4+offset)%bc_seq(4,2) = SCALAR      !Vector

        gv%bc_grp(4+offset)%bc_seq(5,1) = IAZ         !Dep. variable
        gv%bc_grp(4+offset)%bc_seq(5,2) = SCALAR      !Vector

        !Group 5
cc        gv%bc_grp(5+offset)%bc_seq(1,1) = IBCNV_N*AUX !Aux. variable
cc        gv%bc_grp(5+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(5+offset)%bc_seq(1,1) = IBCNV*AUX   !Aux. variable
        gv%bc_grp(5+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(5+offset)%bc_seq(2,1) = IBCOV*AUX   !Aux. variable
        gv%bc_grp(5+offset)%bc_seq(2,2) = VECTOR      !Vector

        !Group 6
        gv%bc_grp(6+offset)%bc_seq(1,1) = IJCNV*AUX   !Aux. variable
        gv%bc_grp(6+offset)%bc_seq(1,2) = VECTOR      !Vector

cc        gv%bc_grp(6+offset)%bc_seq(2,1) = IJCNV_N*AUX !Aux. variable
cc        gv%bc_grp(6+offset)%bc_seq(2,2) = VECTOR      !Vector

        gv%bc_grp(6+offset)%bc_seq(2,1) = IJCOV*AUX   !Aux. variable
        gv%bc_grp(6+offset)%bc_seq(2,2) = VECTOR      !Vector

        !Group 7
        gv%bc_grp(7+offset)%bc_seq(1,1) = IVCNV*AUX   !Aux. variable
        gv%bc_grp(7+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(7+offset)%bc_seq(2,1) = IVCNV_N*AUX !Aux. variable
        gv%bc_grp(7+offset)%bc_seq(2,2) = VECTOR      !Vector

        gv%bc_grp(7+offset)%bc_seq(3,1) = IVCOV*AUX   !Aux. variable
        gv%bc_grp(7+offset)%bc_seq(3,2) = VECTOR      !Vector

        gv%bc_grp(7+offset)%bc_seq(4,1) = IVCOV_N*AUX !Aux. variable
        gv%bc_grp(7+offset)%bc_seq(4,2) = VECTOR      !Vector

        gv%bc_grp(7+offset)%bc_seq(5,1) = IVX         !Dep. variable
        gv%bc_grp(7+offset)%bc_seq(5,2) = SCALAR      !Scalar

        gv%bc_grp(7+offset)%bc_seq(6,1) = IVY         !Dep. variable
        gv%bc_grp(7+offset)%bc_seq(6,2) = SCALAR      !Scalar

        gv%bc_grp(7+offset)%bc_seq(7,1) = IVZ         !Dep. variable
        gv%bc_grp(7+offset)%bc_seq(7,2) = SCALAR      !Scalar

        !Group 8
        gv%bc_grp(8+offset)%bc_seq(1,1) = IVECNV*AUX  !Aux. variable
        gv%bc_grp(8+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(8+offset)%bc_seq(2,1) = IVECOV*AUX  !Aux. variable
        gv%bc_grp(8+offset)%bc_seq(2,2) = VECTOR      !Vector

        gv%bc_grp(8+offset)%bc_seq(3,1) = IVECOV_N*AUX!Aux. variable
        gv%bc_grp(8+offset)%bc_seq(3,2) = VECTOR      !Vector

        gv%bc_grp(8+offset)%bc_seq(4,1) = IVEFCNV*AUX !Aux. variable
        gv%bc_grp(8+offset)%bc_seq(4,2) = VECTOR      !Vector

        !Group 9
        gv%bc_grp(9+offset)%bc_seq(1,1) = IDIVPI*AUX  !Aux. variable
        gv%bc_grp(9+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(9+offset)%bc_seq(2,1) = IDIVPE*AUX  !Aux. variable
        gv%bc_grp(9+offset)%bc_seq(2,2) = VECTOR      !Vector

        gv%bc_grp(9+offset)%bc_seq(3,1) = IEH*AUX     !Aux. variable
        gv%bc_grp(9+offset)%bc_seq(3,2) = VECTOR      !Vector

        gv%bc_grp(9+offset)%bc_seq(4,1) = IENI*AUX    !Aux. variable
        gv%bc_grp(9+offset)%bc_seq(4,2) = VECTOR      !Vector
 
#else
        !Group 1    
        gv%bc_grp(1+offset)%bc_seq(1,1) = IRHO_AUX*AUX!Aux. variable
        gv%bc_grp(1+offset)%bc_seq(1,2) = SCALAR      !Scalar

        !Group 2
        gv%bc_grp(2+offset)%bc_seq(1,1) = ITMP_AUX*AUX!Aux. Variable
        gv%bc_grp(2+offset)%bc_seq(1,2) = SCALAR      !Scalar

        !Group 3
        gv%bc_grp(3+offset)%bc_seq(1,1) = IETA*AUX    !Aux. variable
        gv%bc_grp(3+offset)%bc_seq(1,2) = SCALAR      !Scalar

        gv%bc_grp(3+offset)%bc_seq(2,1) = INU*AUX     !Aux. variable
        gv%bc_grp(3+offset)%bc_seq(2,2) = SCALAR      !Scalar

        !Group 4
cc        gv%bc_grp(4+offset)%bc_seq(1,1) = IBCNV_N*AUX !Aux. variable
cc        gv%bc_grp(4+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(4+offset)%bc_seq(1,1) = IBCNV*AUX   !Aux. variable
        gv%bc_grp(4+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(4+offset)%bc_seq(2,1) = IBCOV*AUX   !Aux. variable
        gv%bc_grp(4+offset)%bc_seq(2,2) = VECTOR      !Vector

        gv%bc_grp(4+offset)%bc_seq(3,1) = IBX         !Dep. variable
        gv%bc_grp(4+offset)%bc_seq(3,2) = SCALAR      !Scalar

        gv%bc_grp(4+offset)%bc_seq(4,1) = IBY         !Dep. variable
        gv%bc_grp(4+offset)%bc_seq(4,2) = SCALAR      !Scalar

        gv%bc_grp(4+offset)%bc_seq(5,1) = IBZ         !Dep. variable
        gv%bc_grp(4+offset)%bc_seq(5,2) = SCALAR      !Scalar

        !Group 5
        gv%bc_grp(5+offset)%bc_seq(1,1) = IJCNV*AUX   !Aux. variable
        gv%bc_grp(5+offset)%bc_seq(1,2) = VECTOR      !Vector

cc        gv%bc_grp(5+offset)%bc_seq(2,1) = IJCNV_N*AUX !Aux. variable
cc        gv%bc_grp(5+offset)%bc_seq(2,2) = VECTOR      !Vector

        gv%bc_grp(5+offset)%bc_seq(2,1) = IJCOV*AUX   !Aux. variable
        gv%bc_grp(5+offset)%bc_seq(2,2) = VECTOR      !Vector

        !Group 6
        gv%bc_grp(6+offset)%bc_seq(1,1) = IVCNV*AUX   !Aux. variable
        gv%bc_grp(6+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(6+offset)%bc_seq(2,1) = IVCNV_N*AUX !Aux. variable
        gv%bc_grp(6+offset)%bc_seq(2,2) = VECTOR      !Vector

        gv%bc_grp(6+offset)%bc_seq(3,1) = IVCOV*AUX   !Aux. variable
        gv%bc_grp(6+offset)%bc_seq(3,2) = VECTOR      !Vector

        gv%bc_grp(6+offset)%bc_seq(4,1) = IVCOV_N*AUX !Aux. variable
        gv%bc_grp(6+offset)%bc_seq(4,2) = VECTOR      !Vector

        gv%bc_grp(6+offset)%bc_seq(5,1) = IVX         !Dep. variable
        gv%bc_grp(6+offset)%bc_seq(5,2) = SCALAR      !Scalar

        gv%bc_grp(6+offset)%bc_seq(6,1) = IVY         !Dep. variable
        gv%bc_grp(6+offset)%bc_seq(6,2) = SCALAR      !Scalar

        gv%bc_grp(6+offset)%bc_seq(7,1) = IVZ         !Dep. variable
        gv%bc_grp(6+offset)%bc_seq(7,2) = SCALAR      !Scalar

        if (fake_ve) then
          !Group 7
          gv%bc_grp(7+offset)%bc_seq(1,1) = IDIVPI*AUX  !Aux. variable
          gv%bc_grp(7+offset)%bc_seq(1,2) = VECTOR      !Vector

          gv%bc_grp(7+offset)%bc_seq(2,1) = IEH*AUX    !Aux. variable
          gv%bc_grp(7+offset)%bc_seq(2,2) = VECTOR      !Vector

          gv%bc_grp(7+offset)%bc_seq(3,1) = IENI*AUX    !Aux. variable
          gv%bc_grp(7+offset)%bc_seq(3,2) = VECTOR      !Vector

          !Group 8
          gv%bc_grp(8+offset)%bc_seq(1,1) = IVECNV*AUX  !Aux. variable
          gv%bc_grp(8+offset)%bc_seq(1,2) = VECTOR      !Vector

          gv%bc_grp(8+offset)%bc_seq(2,1) = IVECOV*AUX  !Aux. variable
          gv%bc_grp(8+offset)%bc_seq(2,2) = VECTOR      !Vector

          gv%bc_grp(8+offset)%bc_seq(3,1) = IVECOV_N*AUX!Aux. variable
          gv%bc_grp(8+offset)%bc_seq(3,2) = VECTOR      !Vector

          gv%bc_grp(8+offset)%bc_seq(4,1) = IVEFCNV*AUX !Aux. variable
          gv%bc_grp(8+offset)%bc_seq(4,2) = VECTOR      !Vector

          !Group 9
          gv%bc_grp(9+offset)%bc_seq(1,1) = IDIVPE*AUX  !Aux. variable
          gv%bc_grp(9+offset)%bc_seq(1,2) = VECTOR      !Vector
        else
          !Group 7
          gv%bc_grp(7+offset)%bc_seq(1,1) = IVECNV*AUX  !Aux. variable
          gv%bc_grp(7+offset)%bc_seq(1,2) = VECTOR      !Vector

          gv%bc_grp(7+offset)%bc_seq(2,1) = IVECOV*AUX  !Aux. variable
          gv%bc_grp(7+offset)%bc_seq(2,2) = VECTOR      !Vector
                     
          gv%bc_grp(7+offset)%bc_seq(3,1) = IVECOV_N*AUX!Aux. variable
          gv%bc_grp(7+offset)%bc_seq(3,2) = VECTOR      !Vector

          !Group 8
          gv%bc_grp(8+offset)%bc_seq(1,1) = IDIVPI*AUX  !Aux. variable
          gv%bc_grp(8+offset)%bc_seq(1,2) = VECTOR      !Vector

          gv%bc_grp(8+offset)%bc_seq(2,1) = IDIVPE*AUX  !Aux. variable
          gv%bc_grp(8+offset)%bc_seq(2,2) = VECTOR      !Vector

          gv%bc_grp(8+offset)%bc_seq(3,1) = IEH*AUX    !Aux. variable
          gv%bc_grp(8+offset)%bc_seq(3,2) = VECTOR      !Vector

          gv%bc_grp(8+offset)%bc_seq(4,1) = IENI*AUX    !Aux. variable
          gv%bc_grp(8+offset)%bc_seq(4,2) = VECTOR      !Vector
        endif
#endif

      endif

c     End program

      end subroutine setBCsched

#endif

c     ExB_flow_bc
c     #################################################################
      subroutine ExB_flow_bc(igrid,nx,ny,nz,bcnv,prs,eta,v0,dbcov)

c     -----------------------------------------------------------------
c     Finds ExB flow boundary condition when Eo is imposed at boundary
c     (returns contravariant components).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,nx,ny,nz
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1)   :: prs,eta
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3) :: bcnv,v0
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3),optional :: dbcov

c     Local variables

      integer :: i,j,k,ig,jg,kg
      real(8) :: b2,bbcnv(3),bbcov(3),dum(3),grad_p(3),etah,gsuper(3,3)
     .          ,EE(3)

c     Begin program

      if (.not.pinch_flow) return

      if (igrid /= 1) then
        write (*,*) 'ExB_flow_bc not ready for MG treatment'
        stop
      endif

      i = nx

      do k=1,nz
        do j=1,ny
          bbcnv =0.5*(bcnv(i+1,j,k,:)/gmetric%grid(igrid)%jac(i+1,j,k)
     .               +bcnv(i  ,j,k,:)/gmetric%grid(igrid)%jac(i  ,j,k))
          bbcov =0.5*(matmul(gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .                      ,bcnv(i  ,j,k,:))
     .               +matmul(gmetric%grid(igrid)%gsub(i+1,j,k,:,:)
     .                      ,bcnv(i+1,j,k,:)))

          b2 = dot_product(bbcnv,bbcov)

          if (b2 > 0d0) then
            if (PRESENT(dbcov)) then  !Linearize ExB flow BC
              dum = 0.5*(dbcov(i,j,k,:)+dbcov(i+1,j,k,:))
              bbcov = dum - 2*dot_product(bbcnv,dum)/b2*bbcov
            endif

cc            etah = 0.5*(eta(i+1,j,k)+eta(i,j,k))
            etah = 2d0/(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

            grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,hex=1)

            gsuper = 0.5*(gmetric%grid(igrid)%gsup(i+1,j,k,:,:)
     .                   +gmetric%grid(igrid)%gsup(i  ,j,k,:,:))

            grad_p = matmul(gsuper,grad_p)

            EE = E0 + perturb_E(i,j,k,igrid)

            !ExB velocity
            v0(i+1,j,k,1)=(EE(2)*bbcov(3)-EE(3)*bbcov(2)
     .                    -etah*grad_p(1))/b2
            v0(i+1,j,k,2)=(EE(3)*bbcov(1)-EE(1)*bbcov(3)
     .                    -etah*grad_p(2))/b2
            v0(i+1,j,k,3)=(EE(1)*bbcov(2)-EE(2)*bbcov(1)
     .                    -etah*grad_p(3))/b2
          endif
        enddo
      enddo

c     End program

      end subroutine ExB_flow_bc

c     J_pinch_bc
c     #################################################################
      subroutine J_pinch_bc(igrid,nx,ny,nz,bcnv,prs,eta,j0,first_order)

c     -----------------------------------------------------------------
c     Finds current boundary condition when Eo is imposed at boundary
c     (returns covariant current).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,nx,ny,nz
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1)   :: prs,eta
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3) :: bcnv,j0

      logical,optional :: first_order

c     Local variables

      integer :: i,j,k,ig,jg,kg
      real(8) :: b2,bbcnv(3),bbcov(3),grad_p(3),alpha
     .          ,gsub(3,3),ietah,EE(3)
      logical :: fo

c     Begin program

      if (.not.pinch_flow) return

      if (igrid /= 1) then
        write (*,*) 'J_pinch_bc not ready for MG treatment'
        stop
      endif

      if (PRESENT(first_order)) then
        fo = first_order
      else
        fo = .false.
      endif

      i = nx

      if (.not.fo) then
        do k=1,nz
          do j=1,ny
            bbcnv = 0.5*(bcnv(i+1,j,k,:)+bcnv(i  ,j,k,:))
            bbcov = 0.5*(matmul(gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .                         ,bcnv(i  ,j,k,:))
     .                  +matmul(gmetric%grid(igrid)%gsub(i+1,j,k,:,:)
     .                         ,bcnv(i+1,j,k,:)))

            b2 = dot_product(bbcov,bbcnv)

            if (b2 > 0d0) then
cc              ietah = 2./(eta(i+1,j,k)+eta(i,j,k))
              ietah =0.5*(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

              grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,hex=1)

cc              if (have_jparB) then
cc                alpha = jpar_B(i,j,k)
cc              else
                EE = E0 + perturb_E(i,j,k,igrid)
                alpha = dot_product(EE,bbcnv)*ietah/b2
cc              endif

              j0(i+1,j,k,:)=alpha*bbcnv + cross_product(bbcov,grad_p)/b2
            endif
          enddo
        enddo
      else
        do k=1,nz
          do j=1,ny
            bbcnv = bcnv(i,j,k,:)

cc            gsub = gmetric%grid(igrid)%gsub(i,j,k,:,:)
            gsub = 0.5*(gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .                 +gmetric%grid(igrid)%gsub(i+1,j,k,:,:))
            bbcov = matmul(gsub,bcnv(i,j,k,:))

            b2 = dot_product(bbcov,bbcnv)
            if (b2 > 0d0) then
cc              ietah = 2./(eta(i+1,j,k)+eta(i,j,k))
              ietah =0.5*(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

              grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,hex=1)

cc              if (have_jparB) then
cc                alpha = jpar_B(i,j,k)
cc              else
                EE = E0 + perturb_E(i,j,k,igrid)
                alpha = dot_product(EE,bbcnv)*ietah/b2
cc              endif

              j0(i+1,j,k,:)=alpha*bbcnv + cross_product(bbcov,grad_p)/b2

            endif
          enddo
        enddo
      endif

c     End program

      end subroutine J_pinch_bc

c     perturb_E
c     ##############################################################
      function perturb_E(i,j,k,igrid) result(dE)

      implicit none

c     Call variables

      integer :: i,j,k,igrid
      real(8) :: dE(3)

c     Local variables

      integer :: ig,jg,kg
      real(8) :: mm,kk,RR,x1,y1,z1,etah

c     Begin program

      dE=0d0

cc      !Tangential E-field to produce radial B
cc      etah =2d0/(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))
cc
cc      select case(coords)
cc      case('hel','tor','cyl')
cc        mm = grid_params%params(1)
cc        kk = grid_params%params(2)
cc        RR = grid_params%params(3)
cc      case default
cc        kk = 0d0
cc        mm = 0d0
cc        RR = 0d0
cc      end select
cc
cc      call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
cc     .                              ,ig,jg,kg,x1,y1,z1)
cc
cc      !Cylindrical components
cc      dE(1) = 0d0
cc      dE(2) = kk*br_pert_bc*cos(mm*y1+kk*z1)
cc      dE(3) =-mm*br_pert_bc*cos(mm*y1+kk*z1)/x1
cc
cc      !Transform to covariant
cc      select case(coords)
cc      case('hel')
cc        dE(2) = x1*dE(2)
cc        dE(3) = dE(3)-kk/mm*dE(2)
cc      case('cyl')
cc        dE(2) = x1*dE(2)
cc      case('tor')
cc        dE(2) = x1*dE(2)
cc        dE(3) = RR*dE(3)
cc      end select
cc
cc      dE = etah*dE

      end function perturb_E

ccc     res_bc
ccc     #############################################################
cc      function res_bc(nx,ny,nz,igrid,field) result(rsvty)
ccc     -------------------------------------------------------------
ccc     This function computes the resistivity at the boundary on
ccc     the grid level igrid.
ccc     -------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer :: nx,ny,nz,igrid
cc      real(8) :: field(0:nx+1,0:ny+1,0:nz+1)
cc     .          ,rsvty(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc      integer :: i,j,k,imin,imax,jmin,jmax,kmin,kmax,dim,loc
cc
ccc     Begin program
cc
cc      do dim=1,3
cc        do loc=0,1
cc
cc          call findBCLoopLimits(dim,loc
cc     .                     ,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
cc     .                     ,imin ,imax ,jmin ,jmax ,kmin ,kmax )
cc
cc          do k=kmin,kmax
cc            do j=jmin,jmax
cc              do i=imin,imax
cc                rsvty(i,j,k)=res_ijk(i,j,k,igrid,field(i,j,k))
cc              enddo
cc            enddo
cc          enddo
cc
cc        enddo
cc      enddo
cc
ccc     End program
cc
cc      end function res_bc
cc
ccc     vis_bc
ccc     #############################################################
cc      function vis_bc(nx,ny,nz,igrid,field) result(vscty)
ccc     -------------------------------------------------------------
ccc     This function computes the viscosity at the boundary
ccc     on the grid level igrid.
ccc     -------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer :: nx,ny,nz,igrid
cc      real(8) :: field(0:nx+1,0:ny+1,0:nz+1)
cc     .          ,vscty(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc      integer :: i,j,k,imin,imax,jmin,jmax,kmin,kmax,dim,loc
cc
ccc     Begin program
cc
cc      do dim=1,3
cc        do loc=0,1
cc
cc          call findBCLoopLimits(dim,loc
cc     .                     ,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
cc     .                     ,imin ,imax ,jmin ,jmax ,kmin ,kmax )
cc
cc          do k=kmin,kmax
cc            do j=jmin,jmax
cc              do i=imin,imax
cc                vscty(i,j,k)=vis_ijk(i,j,k,igrid,field(i,j,k))
cc              enddo
cc            enddo
cc          enddo
cc
cc        enddo
cc      enddo
cc
ccc     End program
cc
cc      end function vis_bc

      end module local_BCS_variables

