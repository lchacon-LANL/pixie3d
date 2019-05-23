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

        use xdraw_io

        use app_iosetup

        integer :: order_bc=1

        logical :: regularize_A = .false.  !Whether we regularize A @ SP

        logical :: no_eq_bc,save_tn,no_tn_bc,bc_vars_allcd

        !Resistive wall vars
        logical :: res_wall_bc   =.true.
     .            ,implicit_rw_bc=.false.
     .            ,rw_dbg        =.true.

        real(8),allocatable,dimension(:) :: Brn_m,Bin_m

        real(8) :: tau_w=1d0,r_v=1.5d0

      contains

c     allocBCVariables
c     ###################################################################
      subroutine allocBCVariables

c     -------------------------------------------------------------------
c     Allocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer :: alloc_stat

        integer :: nx,ny,nz
      
c     Begin program

        nx = gv%gparams%nxv(1) 
        ny = gv%gparams%nyv(1)
        nz = gv%gparams%nzv(1)

        allocate(E_bc(0:nx+1,0:ny+1,0:nz+1,3))

        E_bc = 0d0
      
        if(res_wall_bc.and.coords=='hel')allocate(Brn_m(nyd),Bin_m(nyd))

        bc_vars_allcd = .true.

c     End program

      end subroutine allocBCVariables

c     deallocBCVariables
c     ###################################################################
      subroutine deallocBCVariables

c     -------------------------------------------------------------------
c     Deallocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer :: alloc_stat

c     Begin program

        deallocate(E_bc,Brn_m,Bin_m,stat=alloc_stat)
      
        bc_vars_allcd = .false.

c     End program

      end subroutine deallocBCVariables

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

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX) = EQU

        if (gamma > 1d0)
     .       where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU

      case ('ppnsl','p3nsl','ohtor')           !No slip BCs, pinch flow

        pinch_flow = .true.

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = EQU

        !Dirichlet boundary conditions for tangential components
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) = EQU
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ) = EQU

        if (gamma > 1d0)
     .       where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = EQU !Do not change EQU

      case ('vmec','tor','stell')

cc        where (abs(bbcs(:,IBX:IBZ)) == NEU) bbcs(:,IBX:IBZ) = DEF !for testing

        !Dirichlet boundary conditions for tangential components
        where (abs(bbcs(:,IVY)) == NEU) bbcs(:,IVY) = DIR
        where (abs(bbcs(:,IVZ)) == NEU) bbcs(:,IVZ) = DIR

cc        where (abs(bbcs(:,IRHO)) == NEU) bbcs(:,IRHO) = EQU
        where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = DIR

      case ('efit')

cc        where (abs(bbcs(:,IBX:IBZ)) == NEU) bbcs(:,IBX:IBZ) = DEF !for testing

        !Dirichlet boundary conditions for tangential components
        where (abs(bbcs(:,IVY)) == NEU) bbcs(:,IVY) = DIR
        where (abs(bbcs(:,IVZ)) == NEU) bbcs(:,IVZ) = DIR

cc        where (abs(bbcs(:,IRHO)) == NEU) bbcs(:,IRHO) = EQU
cc        where (abs(bbcs(:,IPRS)) == NEU) bbcs(:,IPRS) = DIR

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
      varray%array_var(IPRS)%bc_dep_list(1,1) = ITMP_AUX*AUX
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

      integer :: i,j,k,ig,jg,kg,ord,nx,ny,nz

      real(8) :: b2,bbcnv(3),bbcov(3),dum(3),grad_p(3),etah,gsup(3,3)
     .          ,gsub(3,3),jac,E00(3)

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
            E00 = E_0(i,j,k,:) + E0
            
            v0(i,j,k,1)=(E00(2)*bbcov(3)-E00(3)*bbcov(2)
     .                 -etah*grad_p(1))/b2
            v0(i,j,k,2)=(E00(3)*bbcov(1)-E00(1)*bbcov(3)
     .                 -etah*grad_p(2))/b2
            v0(i,j,k,3)=(E00(1)*bbcov(2)-E00(2)*bbcov(1)
     .                 -etah*grad_p(3))/b2

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

      integer :: i,j,k,ig,jg,kg,ord,nx,ny,nz
      real(8) :: Jb2,b2,bbcnv(3),bbcov(3),grad_p(3),alpha
     .          ,gsub(3,3),ietah,jach,E00(3)

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
     $                             /(g_def%xx(ig) -g_def%xx(ig-1))
     $                               *(0.5*g_def%xx(ig+1)
     $                                +0.5*g_def%xx(ig)
     $                                -    g_def%xx(ig-1))

          gsub = 0.5*(g_def%gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .               +g_def%gmetric%grid(igrid)%gsub(i+1,j,k,:,:))
          bbcov = matmul(gsub,bbcnv)
c ====================
cc          !High-order extrapolation of B to boundary, arbitrary mesh (critical)
cc          vec = (/ bcnv(i-3,j,k,1)
cc     $            ,bcnv(i-2,j,k,1)
cc     $            ,bcnv(i-1,j,k,1)
cc     $            ,bcnv(i  ,j,k,1) /)
cc          call IntDriver1d (4,x(1:4),vec,1,xt,bbcnv(1:1),ord)
cc          vec = (/ bcnv(i-3,j,k,2)
cc     $            ,bcnv(i-2,j,k,2)
cc     $            ,bcnv(i-1,j,k,2)
cc     $            ,bcnv(i  ,j,k,2) /)
cc          call IntDriver1d (4,x(1:4),vec,1,xt,bbcnv(2:2),ord)
cc          vec = (/ bcnv(i-3,j,k,3)
cc     $            ,bcnv(i-2,j,k,3)
cc     $            ,bcnv(i-1,j,k,3)
cc     $            ,bcnv(i  ,j,k,3) /)
cc          call IntDriver1d (4,x(1:4),vec,1,xt,bbcnv(3:3),ord)
cc
cc          !High-order extrapolation of B to boundary, arbitrary mesh (critical)
cc          vec = (/ bcov(i-3,j,k,1)
cc     $            ,bcov(i-2,j,k,1)
cc     $            ,bcov(i-1,j,k,1)
cc     $            ,bcov(i  ,j,k,1) /)
cc          call IntDriver1d (4,x(1:4),vec,1,xt,bbcov(1:1),ord)
cc          vec = (/ bcov(i-3,j,k,2)
cc     $            ,bcov(i-2,j,k,2)
cc     $            ,bcov(i-1,j,k,2)
cc     $            ,bcov(i  ,j,k,2) /)
cc          call IntDriver1d (4,x(1:4),vec,1,xt,bbcov(2:2),ord)
cc          vec = (/ bcov(i-3,j,k,3)
cc     $            ,bcov(i-2,j,k,3)
cc     $            ,bcov(i-1,j,k,3)
cc     $            ,bcov(i  ,j,k,3) /)
cc          call IntDriver1d (4,x(1:4),vec,1,xt,bbcov(3:3),ord)
c ====================

          jach = 0.5*(g_def%gmetric%grid(igrid)%jac(i  ,j,k)
     .               +g_def%gmetric%grid(igrid)%jac(i+1,j,k))

          Jb2 = dot_product(bbcov,bbcnv)
          b2  = Jb2/jach

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

            E00 = E_0(i+1,j,k,:)+E0
            alpha = dot_product(E00,bbcnv)*ietah/Jb2

            if (PRESENT(dbcnv)) then
              bbcnv = dbcnv(i,j,k,:)
            endif

            j0(i+1,j,k,:)=alpha*bbcnv + cross_product(bbcov,grad_p)/b2
          endif
        enddo
      enddo

c     End program

      end subroutine J_pinch_bc

c     find_Ebc
c     #################################################################
      subroutine find_Ebc(g_def,bcov_n,bcov,E_0)

c     -----------------------------------------------------------------
c     Finds tangential cov electric field components according to
c     resistive wall BCs.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      real(8),dimension(0:,0:,0:,:) :: bcov_n,bcov,E_0

c     Local variables

c     Begin program

      if (.not.(implicit_rw_bc.or.save_tn)) return !Explicit RW BC
 
      if (implicit_rw_bc) then
cc        write(*,*) "bcov_n BC",sqrt(sum(bcov_n(nxd-1:nxd,1:nyd,1,:)**2))
cc        write(*,*) "bcov   BC",sqrt(sum(bcov  (nxd-1:nxd,1:nyd,1,:)**2))
        call rw_bc(g_def,bcov_n,bcov  ,E_0)
      else
cc        write(*,*) "bcov_n BC",sqrt(sum(bcov_n(nxd-1:nxd,1:nyd,1,:)**2))
        call rw_bc(g_def,bcov_n,bcov_n,E_0)
      endif
c$$$      write (*,*) "E_bc BC",sqrt(sum(E_0(nxd+1,1:nyd,1,:)**2))
c$$$     .             ,sqrt(sum(E0**2))

c     End program

      end subroutine find_Ebc
      
c     rw_bc
c     #################################################################
      subroutine rw_bc(g_def,bcov_n,bcov,E_0)

c     -----------------------------------------------------------------
c     Finds tangential cov electric field components according to
c     resistive wall BCs.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      real(8),dimension(0:,0:,0:,:) :: bcov_n,bcov,E_0

c     Local variables

      integer :: i,j,k,nx,ny,nz,nn,mm

      real(8) :: phinm  (0:size(bcov,2)-1
     .                  ,0:size(bcov,3)-1)
     .          ,bcovext(0:size(bcov,2)-1
     .                  ,0:size(bcov,3)-1,2:3)
     .          ,bcovint(0:size(bcov,2)-1
     .                  ,0:size(bcov,3)-1,2:3)

      real(8) :: dum,dh,bcov_avg(2:3),send_buf(2),rec_buf(2)
      
c     Begin program
     
      if (g_def%ihi(1) /= nxd) return
      
      nx = size(bcov,1)-2
      ny = size(bcov,2)-2
      nz = size(bcov,3)-2

      i = nx + 1

      select case(coords)
      case('hel')
        call rw_hel
cc      case('cyl')
cc      case('tor')
      case default
        call pstop('rw_bc'
     .       ,'Resistive wall BC not implemented for '//coords)
      end select

c     Find internal tangential magnetic field components, EXCLUDING m=n=0 mode

      !Extrapolate to face
      bcovint(:,:,2:3) = 1.5*bcov(i-1,:,:,2:3)
     .                  -0.5*bcov(i-2,:,:,2:3)

      !Perform angular average
      bcov_avg(2) = sum(bcovint(1:ny,1:nz,2))
      bcov_avg(3) = sum(bcovint(1:ny,1:nz,3))

#if defined(petsc)
      send_buf = bcov_avg(2:3)
      call MPI_Allreduce(send_buf,rec_buf,2,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,g_def%MPI_COMM_YZ,mpierr)
      bcov_avg = rec_buf
#endif

      bcov_avg = bcov_avg/(nyd*nzd)
      
      !Subtract average
      bcovint(:,:,2) = bcovint(:,:,2) - bcov_avg(2)
      bcovint(:,:,3) = bcovint(:,:,3) - bcov_avg(3)
      
c     Find electric field components @ r=1 (store in (nx+1))

      !Cnv from curl
      E_0(i,:,:,1) = 0d0
      E_0(i,:,:,2) =-(bcovext(:,:,3)-bcovint(:,:,3))/tau_w
      E_0(i,:,:,3) = (bcovext(:,:,2)-bcovint(:,:,2))/tau_w

      !Xform E_0 to Cov @ boundary
      do k=0,nz+1
        do j=0,ny+1
          E_0(i,j,k,:)
     .         =XformToCov_ijk_face(g_def,i-1,j,k,1,E_0(i,j,k,:),1)
        enddo
      enddo

c$$$      write (*,*) "Bext",sqrt(sum(bcovext(:,:,2:3)**2))
c$$$      write (*,*) "||Ebc||",time,sqrt(sum(E_0(i,1:ny,1:nz,2:3)**2))
c$$$      dum = sqrt(sum(E_0(i,:,:,2:3)**2))
c$$$      if (my_rank == 0) write (*,*) "||Ebc||",time,dum
c$$$      write (*,*) '==========================='
c$$$      write (*,*) 'Helical cov E'
c$$$      do j=1,nyd
c$$$        write (*,*) i,j,E_0(i,j,1,2:3)
c$$$      enddo
c$$$      write (*,*) '==========================='
c$$$      write (*,*) 'Inside helical cov B'
c$$$      do j=1,nyd
c$$$        write (*,*) i,j,bcovint(j,1,2:3)
c$$$      enddo
c$$$      write (*,*) '==========================='
c$$$      write (*,*) 'Outside helical cov B'
c$$$      do j=1,nyd
c$$$        write (*,*) i,j,bcovext(j,1,2:3)
c$$$      enddo
c$$$      stop

c     End program

      contains

c     rw_hel
c     ####################################################################
      subroutine rw_hel

      integer :: i_bc,m,n,jglobal,i,j,m1,mm2,mm1,my_local_rank
      real(8) :: theta(nyd)
     .          ,Bn(1,nyd,1),gBn(1,nyd,1)
     .          ,Bnl(1,ny,1),gBnl(1,ny,1)
     .          ,Br_m(nyd),Bi_m(nyd)
     .          ,gBr_m(nyd),gBi_m(nyd),mm(nyd)
     .          ,phir_m(nyd),phii_m(nyd)
     .          ,d2phir_m(nyd),d2phii_m(nyd)
     .          ,dBr(nyd),Bth(nyd),Bz(nyd)
      real(8) :: dImv,dIm,Im,d2Im,dKmv,dKm,Km,d2Km,pv,pw,RR,kk,rm,rbc
     .          ,efac,gamma,dum

      logical :: debug_io
c     Begin program
      
      i_bc = nx + 1

#if defined(petsc)
      call MPI_Comm_rank(g_def%MPI_COMM_Y,my_local_rank,mpierr)
#else
      my_local_rank = 0
#endif
      
      debug_io = old_time_eval.and.rw_dbg.and.my_local_rank==0
      
c     Global angular mesh
      
      theta = g_def%yg(1:nyd)

c     Process Bnormal (old time)

#if !defined(plotting)
      if (save_tn) then
#endif
        !Extrapolate to r=1 old-time Bn
c$$$        Bnl(1,:,1) = 1.5*bcov_n(i_bc-1,1:ny,1,1)
c$$$     .              -0.5*bcov_n(i_bc-2,1:ny,1,1)
        Bnl(1,:,1) = 0.5*bcov_n(i_bc  ,1:ny,1,1)
     .              +0.5*bcov_n(i_bc-1,1:ny,1,1)

        !Gather Bnormal in parallel along boundary
#if defined(petsc)
        call find_global_nobc(Bnl ,Bn,mpi_comm=g_def%MPI_COMM_Y)
#else
        call find_global_nobc(Bnl ,Bn)
#endif

        !Perform FFT
        call fft(nyd,Bn (1,:,1),theta,nyd,Brn_m ,Bin_m ,mm,0,0)

#if !defined(plotting)
      endif
#endif
      
c     Process gradient of Bnormal (new time)
      
      !Compute new-time gradient of Bn
      dh = g_def%xg(nxd)-g_def%xg(nxd-1)
      if (implicit_rw_bc) then
        gBnl(1,:,1) =(bcov(i_bc-1,1:ny,1,1)
     .               -bcov(i_bc-2,1:ny,1,1))/dh
      else
        gBnl(1,:,1) =(bcov_n(i_bc  ,1:ny,1,1)
     .               -bcov_n(i_bc-1,1:ny,1,1))/dh
      endif
      
      !Gather gradient in parallel along boundary
#if defined(petsc)
      call find_global_nobc(gBnl,gBn,mpi_comm=g_def%MPI_COMM_Y)
#else
      call find_global_nobc(gBnl,gBn)
#endif

      !Perform FFT
      call fft(nyd,gBn(1,:,1),theta,nyd,gBr_m,gBi_m,mm,0,0)

      !Find potential solution Phi outside of resistive wall
      kk = g_def%params(2)

      if (kk == 0d0) then
        do m=1,nyd
          m1   = abs(nint(mm(m)))
          if (m1 == 0) then
            phir_m(m) = 0d0 ; phii_m(m) = 0d0
          else
            rm = r_v**(2*mm(m))
            
            !Evolve normal B component
            if (.not.source_eval) then
              gamma = 1 + mm(m)*(rm+1d0)/(rm-1d0)
              efac  = exp(-gamma*dt/tau_w)

              Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
              Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)

              if(m1==nint(g_def%params(1)).and.debug_io) then
                if (Brn_m(m)*Br_m(m) > 0d0) then  
                  dum=log(Br_m(m)/Brn_m(m))*tau_w/dt
                else
                  dum=0d0
                endif
                write (*,*)"Br norm. decay k=0 =",dum
              endif
            else
              Br_m(m) = Brn_m(m)
              Bi_m(m) = Bin_m(m)
            endif

            !Find phi_nm
            dum = (rm+1d0)/((1d0-rm)*mm(m))
            phir_m(m) = Br_m(m)*dum
            phii_m(m) = Bi_m(m)*dum
          endif
        enddo
      else
        pw = abs(kk)
        pv = abs(kk)*r_v
        do m=1,nyd
          m1   = abs(nint(mm(m)))
          mm1  = abs(m1-1)
          mm2  = abs(m1-2)

          !Bessel function evals
          Km   = BESSK(m1 ,pw)                   
          dKm  =-BESSK(mm1,pw) - m1/pw*BESSK(m1,pw)
          d2Km =                   BESSK(mm2,pw)
     .         +  (2*m1-1)/ pw    *BESSK(mm1,pw)
     .         + m1*(m1+1)/(pw*pw)*BESSK(m1 ,pw)
          dKmv =-BESSK(mm1,pv) - m1/pv*BESSK(m1,pv)
          Im   = BESSI(m1 ,pw)
          dIm  = BESSI(mm1,pw) - m1/pw*BESSI(m1,pw)
          d2Im =                   BESSI(mm2,pw)
     .         -  (2*m1-1)/ pw    *BESSI(mm1,pw)
     .         + m1*(m1+1)/(pw*pw)*BESSI(m1 ,pw)
          dImv = BESSI(mm1,pv) - m1/pv*BESSI(m1,pv)

          !Evolve normal B component
          if (.not.source_eval) then
            gamma =-pw*(dImv*d2Km-dKmv*d2Im)/(dImv*dKm -dKmv*dIm)
            efac  = exp(-gamma*dt/tau_w)
            Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
            Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)

            if(m1==nint(g_def%params(1)).and.debug_io) then
              if (Brn_m(m)*Br_m(m) > 0d0) then
                dum=log(Br_m(m)/Brn_m(m))*tau_w/dt
              else
                dum=0d0
              endif
              write (*,*)"Br norm. decay k/=0 =",dum
            endif
          else
            Br_m(m) = Brn_m(m)
            Bi_m(m) = Bin_m(m)
          endif
          
          !Find phi_nm
          dum = (dImv*Km-dKmv*Im)/((dImv*dKm-dKmv*dIm)*pw)
          phir_m(m)=Br_m(m)*dum
          phii_m(m)=Bi_m(m)*dum
        enddo
      endif

      !Bth [multiply by -sqrt(-1)*m; minus sign from convention of iFFT]
      Br_m = mm*phii_m
      Bi_m =-mm*phir_m
      call ifft(nyd,Bth,theta,nyd,Br_m,Bi_m,mm,0)

      !Get LOCAL helical (cov) magnetic field components @ r=1 (store in (nx+1))
      rbc = 1d0

      do k = klom,khip
        do j = jlom,jhip
          jglobal = j + gv%gparams%jlo(1) - 1
          if (jglobal > nyd) jglobal = jglobal - nyd
          if (jglobal < 1  ) jglobal = nyd + jglobal
          bcovext(j,k,2) = rbc*Bth(jglobal)
          bcovext(j,k,3) = 0d0
        enddo
      enddo
    
      end subroutine rw_hel

c$$$c     rw_cyl
c$$$c     ####################################################################
c$$$      subroutine rw_cyl
c$$$
c$$$      integer :: i_bc,m,n,jglobal,i,j,m1,mm2,mm1
c$$$      real(8) :: theta(nyd)
c$$$     .          ,Bn(1,nyd,1),gBn(1,nyd,1)
c$$$     .          ,Bnl(1,ny,1),gBnl(1,ny,1)
c$$$     .          ,Br_m(nyd),Bi_m(nyd)
c$$$     .          ,gBr_m(nyd),gBi_m(nyd),mm(nyd)
c$$$     .          ,phir_m(nyd),phii_m(nyd)
c$$$     .          ,d2phir_m(nyd),d2phii_m(nyd)
c$$$     .          ,dBr(nyd),Bth(nyd),Bz(nyd)
c$$$      real(8) :: dImv,dIm,Im,d2Im,dKmv,dKm,Km,d2Km,pv,pw,RR,kk,rm,rbc
c$$$     .          ,efac,gamma,dum
c$$$
c$$$c     Begin program
c$$$      
c$$$      i_bc = nx + 1
c$$$
c$$$c     Global angular mesh
c$$$      
c$$$      theta = g_def%yg(1:nyd)
c$$$
c$$$c     Process Bnormal (old time)
c$$$
c$$$#if !defined(plotting)
c$$$      if ((old_time_eval.or.source_eval)) then
c$$$#endif
c$$$        if (.not.allocated(Brn_m)) allocate(Brn_m(nyd),Bin_m(nyd))
c$$$        
c$$$        !Extrapolate to r=1 old-time Bn
c$$$        Bnl(1,:,1) = 1.5*bcov_n(i_bc-1,1:ny,1,1)
c$$$     .              -0.5*bcov_n(i_bc-2,1:ny,1,1)
c$$$
c$$$        !Gather Bnormal in parallel along boundary
c$$$#if defined(petsc)
c$$$        call find_global_nobc(Bnl ,Bn ,mpi_comm=g_def%MPI_COMM_Y)
c$$$#else
c$$$        call find_global_nobc(Bnl ,Bn )
c$$$#endif
c$$$
c$$$        !Perform FFT
c$$$        call fft(nyd,Bn (1,:,1),theta,nyd,Brn_m ,Bin_m ,mm,0,0)
c$$$
c$$$#if !defined(plotting)
c$$$      endif
c$$$#endif
c$$$      
c$$$c     Process gradient of Bnormal (new time)
c$$$      
c$$$      !Compute new-time gradient of Bn
c$$$      dh = g_def%xg(nxd)-g_def%xg(nxd-1)
c$$$      gBnl(1,:,1) =(bcov(i_bc-1,1:ny,1,1)
c$$$     .             -bcov(i_bc-2,1:ny,1,1))/dh
c$$$
c$$$      !Gather gradient in parallel along boundary
c$$$#if defined(petsc)
c$$$      call find_global_nobc(gBnl,gBn,mpi_comm=g_def%MPI_COMM_Y)
c$$$#else
c$$$      call find_global_nobc(gBnl,gBn)
c$$$#endif
c$$$
c$$$      !Perform FFT
c$$$      call fft(nyd,gBn(1,:,1),theta,nyd,gBr_m,gBi_m,mm,0,0)
c$$$
c$$$      !Find potential solution Phi outside of resistive wall
c$$$      kk = g_def%params(2)
c$$$
c$$$      if (kk == 0d0) then
c$$$        do m=1,nyd
c$$$          m1   = abs(nint(mm(m)))
c$$$          if (mm(m) == 0) then
c$$$            phir_m(m) = 0d0 ; phii_m(m) = 0d0
c$$$          else
c$$$            rm = r_v**(2*mm(m))
c$$$            
c$$$            !Evolve normal B component
c$$$            if (.not.source_eval) then
c$$$              gamma = 1 + mm(m)*(rm+1d0)/(rm-1d0)
c$$$              efac  = exp(-gamma*dt/tau_w)
c$$$
c$$$              Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
c$$$              Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)
c$$$
c$$$              if(m1==nint(g_def%params(1)).and.old_time_eval.and.rw_dbg)
c$$$     .        then
c$$$                if (Brn_m(m) /= 0d0.and. Br_m(m) < 0d0) then  
c$$$                  dum=log(Br_m(m)/Brn_m(m))*tau_w/dt
c$$$                else
c$$$                  dum=0d0
c$$$                endif
c$$$                write (*,*)"Br norm. decay k=0 =",dum
c$$$              endif
c$$$            endif
c$$$
c$$$            !Find phi_nm
c$$$            phir_m(m) = Br_m(m)/mm(m)*(rm+1d0)/(1d0-rm)
c$$$            phii_m(m) = Bi_m(m)/mm(m)*(rm+1d0)/(1d0-rm)
c$$$          endif
c$$$        enddo
c$$$      else
c$$$        pw = abs(kk)
c$$$        pv = abs(kk)*r_v
c$$$        do m=1,nyd
c$$$          m1   = abs(nint(mm(m)))
c$$$          mm1  = abs(m1-1)
c$$$          mm2  = abs(m1-2)
c$$$
c$$$          !Bessel function evals
c$$$          Km   = BESSK(m1 ,pw)                   
c$$$          dKm  =-BESSK(mm1,pw) - m1/pw*BESSK(m1,pw)
c$$$          d2Km =                   BESSK(mm2,pw)
c$$$     .         +  (2*m1-1)/ pw    *BESSK(mm1,pw)
c$$$     .         + m1*(m1+1)/(pw*pw)*BESSK(m1 ,pw)
c$$$          dKmv =-BESSK(mm1,pv) - m1/pv*BESSK(m1,pv)
c$$$          Im   = BESSI(m1 ,pw)
c$$$          dIm  = BESSI(mm1,pw) - m1/pw*BESSI(m1,pw)
c$$$          d2Im =                   BESSI(mm2,pw)
c$$$     .         -  (2*m1-1)/ pw    *BESSI(mm1,pw)
c$$$     .         + m1*(m1+1)/(pw*pw)*BESSI(m1 ,pw)
c$$$          dImv = BESSI(mm1,pv) - m1/pv*BESSI(m1,pv)
c$$$
c$$$          !Evolve normal B component
c$$$          if (.not.source_eval) then
c$$$            gamma =-pw*(dImv*d2Km-dKmv*d2Im)/(dImv*dKm -dKmv*dIm)
c$$$            efac  = exp(-gamma*dt/tau_w)
c$$$            Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
c$$$            Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)
c$$$
c$$$            if(m1==nint(g_def%params(1)).and.old_time_eval.and.rw_dbg)
c$$$     .      then
c$$$              if (Brn_m(m) < 0d0.and.Br_m(m) < 0d0) then
c$$$                dum=log(Br_m(m)/Brn_m(m))*tau_w/dt
c$$$              else
c$$$                dum=0d0
c$$$              endif
c$$$              write (*,*)"Br norm. decay k/=0 =",dum
c$$$            endif
c$$$          endif
c$$$          
c$$$          !Find phi_nm
c$$$          phir_m(m)=Br_m(m)/pw*(dImv*Km-dKmv*Im)/(dImv*dKm-dKmv*dIm)
c$$$          phii_m(m)=Bi_m(m)/pw*(dImv*Km-dKmv*Im)/(dImv*dKm-dKmv*dIm)
c$$$c$$$          d2phir_m(m)=-Br_m(m)*gamma
c$$$c$$$          d2phii_m(m)=-Bi_m(m)*gamma
c$$$        enddo
c$$$      endif
c$$$
c$$$c$$$      !Reconstruct GLOBAL cyl magnetic field components from grad(Phi) outside plasma wall
c$$$c$$$      !Bth [multiply by -sqrt(-1)*m; minus sign from convention of iFFT]
c$$$c$$$      call ifft(nyd,dBr,theta,nyd,d2phir_m,d2phii_m,mm,0)
c$$$      
c$$$      !Bth [multiply by -sqrt(-1)*m; minus sign from convention of iFFT]
c$$$      Br_m = mm*phii_m
c$$$      Bi_m =-mm*phir_m
c$$$      call ifft(nyd,Bth,theta,nyd,Br_m,Bi_m,mm,0)
c$$$
c$$$      !Bz [multiply by -sqrt(-1)*k; minus sign from convention of iFFT]
c$$$      m    = g_def%params(1)
c$$$      Br_m = mm*kk/m*phii_m  !mm vector, m scalar; for helical only
c$$$      Bi_m =-mm*kk/m*phir_m  !mm vector, m scalar; for helical only
c$$$      call ifft(nyd,Bz,theta,nyd,Br_m,Bi_m,mm,0)
c$$$
c$$$      !Get LOCAL helical (cov) magnetic field components @ r=1 (store in (nx+1))
c$$$      rbc = 1d0
c$$$
c$$$      do k = klom,khip
c$$$        do j = jlom,jhip
c$$$          jglobal = j + gv%gparams%jlo(1) - 1
c$$$          if (jglobal > nyd) jglobal = jglobal - nyd
c$$$          if (jglobal < 1  ) jglobal = nyd + jglobal
c$$$c$$$          dbcnv1ext(j,k) = rbc*dBr(jglobal)
c$$$          bcovext(j,k,2) = rbc*Bth(jglobal)
c$$$          bcovext(j,k,3) = Bz (jglobal) - rbc*kk/m*Bth(jglobal)
c$$$        enddo
c$$$      enddo
c$$$    
c$$$      end subroutine rw_cyl

      end subroutine rw_bc

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
      bbcs(:,1) = bcsq

      bcsq = bbcs(:,2)
      where (bcsq(3:4) == DEF) bcsq(3:4) =  DIR
      where (bcsq(1:2) == DEF) bcsq(1:2) = -NEU  !On covariant components
      where (bcsq(5:6) == DEF) bcsq(5:6) = -NEU  !On covariant components
      bbcs(:,2) = bcsq

      bcsq = bbcs(:,3)
      where (bcsq(1:4) == DEF) bcsq(1:4) = -NEU !On covariant components
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
        where (bcsq(1:2) == DEF) bcsq(1:2) = DIR
        where (bcsq(3:6) == DEF) bcsq(3:6) =-NEU !On covariant components
        bbcs(:,1) = bcsq

        bcsq = bbcs(:,2)
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
        bbcs(:,2) = bcsq

        bcsq = bbcs(:,3)
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
        bbcs(:,3) = bcsq
c$$$      endif

      call set_vec_symm_bc(bbcs,-1)

c     End program

      end subroutine default_V_BCs

      end module local_BCS_variables

