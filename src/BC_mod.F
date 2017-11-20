c module local_BCS_variables
c####################################################################
      module local_BCS_variables

        use variables

        use imposeBCinterface, ONLY: extrapolateBC

        use auxiliaryVariables

        use operators

        use ts_setup

        use equilibrium

        use nlfunction_setup, ONLY: res,chi_perp,vis,ion_hall
     .                             ,E_fld_res_A,E_fld_res_B,E_fld_ihall
     .                             ,E_fld_ehall,EOM_divPi,hres
     .                             ,nc_eom_v,nc_eom_jxb
     .                             ,pinch_flow,sym_st,subtract_E0
     .                             ,res_ijk,vis_ijk,form_ve,solenoidal

        integer :: order_bc=1

        logical :: no_eq_bc=.false.,no_tn_bc=.false.

        logical :: save_tn

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

c     Default VEL boundary conditions (no-stress)

      bcsq = bbcs(:,IRHO)
      where (bcsq == DEF .or. bcsq == FSYM) bcsq = NEU
      bbcs(:,IRHO) = bcsq

c$$$      if (di > 0d0) then        !Force D_t(v x n) = 0 for energy conservation
c$$$        bcsq = bbcs(:,IVX)
c$$$        where (bcsq(1:2) == DEF) bcsq(1:2) = DIR
c$$$        where (bcsq(3:6) == DEF) bcsq(3:6) =-EQU !On covariant components
c$$$        bbcs(:,IVX) = bcsq
c$$$
c$$$        bcsq = bbcs(:,IVY)
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
c$$$        bbcs(:,IVY) = bcsq
c$$$
c$$$        bcsq = bbcs(:,IVZ)
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
c$$$        bbcs(:,IVZ) = bcsq
c$$$      else
        bcsq = bbcs(:,IVX)
        where (bcsq(1:2) == DEF) bcsq(1:2) = DIR
        where (bcsq(3:6) == DEF) bcsq(3:6) =-NEU !On covariant components
        bbcs(:,IVX) = bcsq

        bcsq = bbcs(:,IVY)
        if (     coords == 'cyl'
     .      .or. coords == 'hel'
     .      .or. coords == 'tor'
     .      .or. coords == 'car') then
          where (bcsq(1:2) == DEF) bcsq(1:2) =  NEU !On contravariant components
        else
          where (bcsq(1:2) == DEF) bcsq(1:2) = -NEU !On covariant components
        endif
        where (bcsq(3:4) == DEF) bcsq(3:4) = DIR
        where (bcsq(5:6) == DEF) bcsq(5:6) =-NEU
        bbcs(:,IVY) = bcsq

        bcsq = bbcs(:,IVZ)
        if (     coords == 'cyl'
     .      .or. coords == 'hel'
     .      .or. coords == 'tor'
     .      .or. coords == 'car') then
          where (bcsq(1:2) == DEF) bcsq(1:2) =  NEU !On contravariant components
        else
          where (bcsq(1:2) == DEF) bcsq(1:2) = -NEU !On covariant components
        endif
        where (bcsq(3:4) == DEF) bcsq(3:4) =-NEU
        where (bcsq(5:6) == DEF) bcsq(5:6) = DIR
        bbcs(:,IVZ) = bcsq
c$$$      endif

      call set_vec_symm_bc(bbcs(:,IVX:IVZ),-1)

c     Default B-A boundary conditions

#if !defined(vec_pot)
      bcsq = bbcs(:,IBX)
      where (bcsq(1:2) == DEF) bcsq(1:2) =  DIR
      where (bcsq(3:6) == DEF) bcsq(3:6) = -NEU !On covariant components
      bbcs(:,IBX) = bcsq

      bcsq = bbcs(:,IBY)
      where (bcsq(3:4) == DEF) bcsq(3:4) =  DIR
      where (bcsq(1:2) == DEF) bcsq(1:2) = -NEU  !On covariant components
      where (bcsq(5:6) == DEF) bcsq(5:6) = -NEU  !On covariant components
      bbcs(:,IBY) = bcsq

      bcsq = bbcs(:,IBZ)
      where (bcsq(1:4) == DEF) bcsq(1:4) = -NEU !On covariant components
      where (bcsq(5:6) == DEF) bcsq(5:6) =  DIR
      bbcs(:,IBZ) = bcsq

      call set_vec_symm_bc(bbcs(:,IBX:IBZ),1)
#else
      bcsq = bbcs(:,IAX)
      where (bcsq(1:2) == DEF) bcsq(1:2) =  DIR !On contravariant components
      where (bcsq(3:6) == DEF) bcsq(3:6) = -EQU !On covariant components
      bbcs(:,IAX) = bcsq

      bcsq = bbcs(:,IAY)
      where (bcsq(3:4) == DEF) bcsq(3:4) =  DIR
      where (bcsq(1:2) == DEF) bcsq(1:2) = -EQU !On covariant components
      where (bcsq(5:6) == DEF) bcsq(5:6) = -EQU !On covariant components
      bbcs(:,IAY) = bcsq

      bcsq = bbcs(:,IAZ)
      where (bcsq(1:4) == DEF) bcsq(1:4) = -EQU !On covariant components
      where (bcsq(5:6) == DEF) bcsq(5:6) =  DIR
      bbcs(:,IAZ) = bcsq

      call set_vec_symm_bc(bbcs(:,IAX:IAZ),-1)
#endif

c     Default TMP boundary conditions

      bcsq = bbcs(:,ITMP)
      where (bcsq == DEF .or. bcsq == FSYM) bcsq = NEU
cc      if (gamma /= 1d0 .and. (.not.use_p_eom)) then
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

        pinch_flow = .true.

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX ) = EQU

        if (gamma > 1d0)
     .       where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = EQU

      case ('ppnsl','p3nsl','ohtor')           !No slip BCs, pinch flow

        pinch_flow = .true.

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX ) = EQU

        !Dirichlet boundary conditions for tangential components
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY ) = EQU
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ ) = EQU

        if (gamma > 1d0)
     .       where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = EQU

      case ('vmec','tor','stell')

cc        where (abs(bbcs(:,IBX:IBZ)) == NEU) bbcs(:,IBX:IBZ) = DEF !for testing

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

      case('gem','diamh')

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

      case('khihk')

#if defined(vec_pot)
        if (di > 0d0) then !Force D_t(v x n) = 0 to prevent boundary issues
          where (abs(bbcs(:,IVY)) == NEU) bbcs(:,IVY) =-EQU
          where (abs(bbcs(:,IVZ)) == NEU) bbcs(:,IVZ) =-EQU

c$$$          where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = EQU
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

      varray%array_var(ITMP)%dom_dep_list(3,1) = ICHI*AUX
      varray%array_var(ITMP)%dom_dep_list(3,2) = SCALAR
      varray%array_var(ITMP)%dom_dep_list(3,3) = 1 !non-trivial

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
        varray%array_var(ieq)%dom_dep_list(1,3) = 1 !non-trivial

        varray%array_var(ieq)%dom_dep_list(2,1) = IENI_0*AUX
        varray%array_var(ieq)%dom_dep_list(2,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(2,3) = 1 !non-trivial
      enddo
#endif

c     End subroutine

      end subroutine defineDepDom_Deps

c     ExB_flow_bc
c     #################################################################
      subroutine ExB_flow_bc(g_def,igrid,bcnv,prs,eta,v0,dbcov)

c     -----------------------------------------------------------------
c     Finds ExB flow boundary condition when Eo is imposed at boundary
c     (returns contravariant components).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8),dimension(0:,0:,0:)   :: prs,eta
      real(8),dimension(0:,0:,0:,:) :: bcnv,v0
      real(8),dimension(0:,0:,0:,:),optional :: dbcov

c     Local variables

      integer :: i,j,k,ig,jg,kg,ord,nx,ny,nz

      real(8) :: b2,bbcnv(3),bbcov(3),dum(3),grad_p(3),etah,gsup(3,3)
     .          ,gsub(3,3),jac

c     Begin program

      nx = size(bcnv,1)-2
      ny = size(bcnv,2)-2
      nz = size(bcnv,3)-2

      if (igrid /= 1) then
        write (*,*) 'ExB_flow_bc not ready for MG treatment'
        stop
      endif

      i = nx + 1

      do k=0,nz+1
        do j=0,ny+1
          jac = g_def%gmetric%grid(igrid)%jac(i,j,k)

          bbcnv = bcnv(i,j,k,:)
          bbcov = XformToCov_ijk(g_def,i,j,k,igrid,bbcnv)

          b2 = vectorNorm_ijk(g_def,i,j,k,igrid,bbcnv,.false.)

          if (b2 > 0d0) then
            if (PRESENT(dbcov)) then  !Linearize ExB flow BC
              dum = dbcov(i,j,k,:)
              bbcov = dum - 2*dot_product(bbcnv,dum)/(jac*b2)*bbcov
            endif

            grad_p = grad(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid,prs)
            grad_p = XformToCnv_ijk(g_def,i,j,k,igrid,grad_p)

            etah = eta(i,j,k)

            !ExB velocity
            v0(i,j,k,1)=(E0(2)*bbcov(3)-E0(3)*bbcov(2)
     .                 -etah*grad_p(1))/b2
            v0(i,j,k,2)=(E0(3)*bbcov(1)-E0(1)*bbcov(3)
     .                 -etah*grad_p(2))/b2
            v0(i,j,k,3)=(E0(1)*bbcov(2)-E0(2)*bbcov(1)
     .                 -etah*grad_p(3))/b2

          endif
        enddo
      enddo

c     End program

      end subroutine ExB_flow_bc

c     J_pinch_bc
c     #################################################################
      subroutine J_pinch_bc(g_def,igrid,bcnv,prs,eta,j0)

c     -----------------------------------------------------------------
c     Finds current boundary condition when Eo is imposed at boundary
c     (returns covariant current).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8),dimension(0:,0:,0:)   :: prs,eta
      real(8),dimension(0:,0:,0:,:) :: bcnv,j0

c     Local variables

      integer :: i,j,k,ig,jg,kg,ord,nx,ny,nz
      real(8) :: b2,bbcnv(3),bbcov(3),grad_p(3),alpha
     .          ,gsub(3,3),ietah,jach

c$$$      real(8) :: x(5),vec(4),xt(1)

c     Begin program

      nx = size(bcnv,1)-2
      ny = size(bcnv,2)-2
      nz = size(bcnv,3)-2

      if (igrid /= 1) then
        write (*,*) 'J_pinch_bc not ready for MG treatment'
        stop
      endif

      i = nx

      call getMGmap(g_def,i,1,1,igrid,igrid,igrid,ig,jg,kg)

c$$$      ord = 3

c$$$      x = (/g_def%xx(ig-3)
c$$$     $     ,g_def%xx(ig-2)
c$$$     $     ,g_def%xx(ig-1)
c$$$     $     ,g_def%xx(ig  )
c$$$     $     ,g_def%xx(ig+1)/)
c$$$
c$$$      xt = 0.5*(g_def%xx(ig+1)
c$$$     $         +g_def%xx(ig  ))

      do k=1,nz
        do j=1,ny
c ====================
cc          !Zeroth-order extrapolation of B to boundary
cc          bbcnv = bcnv(i  ,j,k,:)
c ====================
cc          !First-order extrapolation of B to boundary, uniform mesh (critical)
cc          bbcnv = 1.5*bcnv(i  ,j,k,:)
cc     $           -0.5*bcnv(i-1,j,k,:)
c ====================
          !First-order extrapolation of B to boundary, arbitrary mesh (critical)
          bbcnv = bcnv(i-1,j,k,:) + (bcnv(i,j,k,:)-bcnv(i-1,j,k,:))
     $                               /(g_def%xx(ig)
     $                                -g_def%xx(ig-1))
     $                               *(0.5*g_def%xx(ig+1)
     $                                +0.5*g_def%xx(ig)
     $                                -    g_def%xx(ig-1))
c ====================
c$$$          vec = (/ bcnv(i-3,j,k,1)
c$$$     $            ,bcnv(i-2,j,k,1)
c$$$     $            ,bcnv(i-1,j,k,1)
c$$$     $            ,bcnv(i  ,j,k,1) /)
c$$$          call IntDriver1d (4,x(1:4),vec,1,xt,bbcnv(1:1),ord)
c$$$          vec = (/ bcnv(i-3,j,k,2)
c$$$     $            ,bcnv(i-2,j,k,2)
c$$$     $            ,bcnv(i-1,j,k,2)
c$$$     $            ,bcnv(i  ,j,k,2) /)
c$$$          call IntDriver1d (4,x(1:4),vec,1,xt,bbcnv(2:2),ord)
c$$$          vec = (/ bcnv(i-3,j,k,3)
c$$$     $            ,bcnv(i-2,j,k,3)
c$$$     $            ,bcnv(i-1,j,k,3)
c$$$     $            ,bcnv(i  ,j,k,3) /)
c$$$          call IntDriver1d (4,x(1:4),vec,1,xt,bbcnv(3:3),ord)
c ====================

          gsub = 0.5*(g_def%gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .               +g_def%gmetric%grid(igrid)%gsub(i+1,j,k,:,:))
          bbcov = matmul(gsub,bbcnv)

          jach = 0.5*(g_def%gmetric%grid(igrid)%jac(i  ,j,k)
     .               +g_def%gmetric%grid(igrid)%jac(i+1,j,k))

          b2 = dot_product(bbcov,bbcnv)/jach

          if (b2 > 0d0) then
cc            ietah = 2./(eta(i+1,j,k)+eta(i,j,k))
            ietah =0.5*(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

c =================
            grad_p
     .         = grad(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,hex=1)
c =================
cc            grad_p = 0.5*(eta (i  ,j,k)
cc     .                   *grad(i  ,j,k,nx,ny,nz,igrid,igrid,igrid,prs)
cc     $                   +eta (i+1,j,k)
cc     .                   *grad(i+1,j,k,nx,ny,nz,igrid,igrid,igrid,prs))
cc     .              *ietah
c =================
c$$$            vec = (/ prs(i-2,j,k)
c$$$     $              ,prs(i-1,j,k)
c$$$     $              ,prs(i  ,j,k)
c$$$     $              ,prs(i+1,j,k) /)
c$$$            call IntDriver1d (4,x(2:5),vec,1,xt,grad_p(1:1),ord,deriv=1)
c =================

            alpha = dot_product(E0,bbcnv)/jach*ietah/b2

            j0(i+1,j,k,:)=alpha*bbcnv + cross_product(bbcov,grad_p)/b2
          endif
        enddo
      enddo

c     End program

      end subroutine J_pinch_bc

      end module local_BCS_variables

