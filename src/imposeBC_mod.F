c module local_BCS_variables
c####################################################################
      module local_BCS_variables

        use BCS_variables

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

c       User-specified (via input) BC configuration parameters

cc        logical :: lag_pinch_bc=.false.,save_Bn  !Whether to lag the finite-flow pinch BC in time

c       Hardwired BC configuration parameters

        !Order of BC interpolation
#if defined(vec_pot)
        integer :: bc_order=2
#else
        integer :: bc_order=1
cc        integer :: bc_order=2
#endif

        INTERFACE XferBCs
          module procedure XferBCs_scl, XferBCs_vec
        END INTERFACE 

      contains

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

c     Exceptions for specific equilibria

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

cc        !Dirichlet boundary conditions for tangential velocity components
cc        where (abs(bbcs(:,IVY )) == NEU) bbcs(:,IVY ) = DIR
cc        where (abs(bbcs(:,IVZ )) == NEU) bbcs(:,IVZ ) = DIR

        where (abs(bbcs(:,ITMP)) == DIR) bbcs(:,ITMP) = NEU
 
cc      case('sgem')  !Symmetric domain
cc
cc        !Dirichlet boundary conditions for tangential velocity components
cc        where (abs(bbcs(:,IVY )) == NEU) bbcs(:,IVY ) = DIR
cc        where (abs(bbcs(:,IVZ )) == NEU) bbcs(:,IVZ ) = DIR
cc
cc        where (abs(bbcs(:,ITMP)) == DIR) bbcs(:,ITMP) = NEU

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
c     Defines BC dependencies for dependent variables.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray

c     Local variables

      integer :: ieq

c     Begin program

      !Rho
      varray%array_var(IRHO)%bc_dep_list(1,1) = IRHO
      varray%array_var(IRHO)%bc_dep_list(1,2) = SCALAR
      varray%array_var(IRHO)%bc_dep_list(1,3) = 0 !trivial

      !Tmp
      varray%array_var(ITMP)%bc_dep_list(1,1) = ITMP
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
          varray%array_var(ieq)%bc_dep_list(1,1) = IRHO
          varray%array_var(ieq)%bc_dep_list(1,2) = SCALAR
          varray%array_var(ieq)%bc_dep_list(1,3) = 0 !trivial

          varray%array_var(ieq)%bc_dep_list(2,1) = IVCNV*AUX
          varray%array_var(ieq)%bc_dep_list(2,2) = VECTOR
          varray%array_var(ieq)%bc_dep_list(2,3) = 0 !trivial
        enddo
      endif

      !B
      do ieq=IBX,IBZ
        varray%array_var(ieq)%bc_dep_list(1,1) = IBCNV*AUX
        varray%array_var(ieq)%bc_dep_list(1,2) = VECTOR
        varray%array_var(ieq)%bc_dep_list(1,3) = 0 !trivial
      enddo

c     End subroutine

      end subroutine defineDepBC_Deps

c     defineDepDom_Deps
c     ####################################################################
      subroutine defineDepDom_Deps(varray)
c     --------------------------------------------------------------------
c     Defines domain dependencies for dependent variables.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray

c     Local variables

      integer :: offset,ieq

c     Begin program

      !Rho
      varray%array_var(IRHO)%dom_dep_list(1,1) = IRHO
      varray%array_var(IRHO)%dom_dep_list(1,2) = SCALAR
      varray%array_var(IRHO)%dom_dep_list(1,3) = 1 !non-trivial

      varray%array_var(IRHO)%dom_dep_list(2,1) = IVCNV*AUX
      varray%array_var(IRHO)%dom_dep_list(2,2) = VECTOR
      varray%array_var(IRHO)%dom_dep_list(2,3) = 1 !non-trivial

      !Tmp
      varray%array_var(ITMP)%dom_dep_list(1,1) = ITMP
      varray%array_var(ITMP)%dom_dep_list(1,2) = SCALAR
      varray%array_var(ITMP)%dom_dep_list(1,3) = 1 !non-trivial

      varray%array_var(ITMP)%dom_dep_list(2,1) = IVCNV*AUX
      varray%array_var(ITMP)%dom_dep_list(2,2) = VECTOR
      varray%array_var(ITMP)%dom_dep_list(2,3) = 1 !non-trivial

      !V
      do ieq=IVX,IVZ
        varray%array_var(ieq)%dom_dep_list(1,1) = IRHO
        varray%array_var(ieq)%dom_dep_list(1,2) = SCALAR
        varray%array_var(ieq)%dom_dep_list(1,3) = 1 !non-trivial

        varray%array_var(ieq)%dom_dep_list(2,1) = ITMP
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

      !B
      do ieq=IBX,IBZ
        varray%array_var(ieq)%dom_dep_list(1,1) = IENI*AUX
        varray%array_var(ieq)%dom_dep_list(1,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(1,3) = 1 !non-trivial

        varray%array_var(ieq)%dom_dep_list(2,1) = IENI_0*AUX
        varray%array_var(ieq)%dom_dep_list(2,2) = VECTOR
        varray%array_var(ieq)%dom_dep_list(2,3) = 1 !non-trivial
      enddo

c     End subroutine

      end subroutine defineDepDom_Deps

ccc     XferBCs
ccc     ###################################################################
cc      function XferBCs(arr1) result(arr2)
cc
ccc     -------------------------------------------------------------------
ccc     Routine to transfer ghost cell values form arr1 to arr2
ccc     -------------------------------------------------------------------
cc
ccc     Input variables
cc
cc      real(8) :: arr1(0:,0:,0:)
cc      real(8) :: arr2(0:size(arr1,1)-1
cc     .               ,0:size(arr1,2)-1
cc     .               ,0:size(arr1,3)-1)
cc
ccc     Local variables
cc
cc      integer :: nxp,nyp,nzp
cc
ccc     Begin program
cc
cc      nxp = size(arr1,1)-1
cc      nyp = size(arr1,2)-1
cc      nzp = size(arr1,3)-1
cc
cc      !X boundaries
cc      arr2(0  ,:,:) = arr1(0  ,:,:)
cc      arr2(nxp,:,:) = arr1(nxp,:,:)
cc
cc      !Y boundaries
cc      arr2(:,0  ,:) = arr1(:,0  ,:)
cc      arr2(:,nyp,:) = arr1(:,nyp,:)
cc
cc      !Z boundaries
cc      arr2(:,:,0  ) = arr1(:,:,0  )
cc      arr2(:,:,nzp) = arr1(:,:,nzp)
cc
ccc     End program
cc
cc      end function XferBCs

c     XferBCs_scl
c     ###################################################################
      subroutine XferBCs_scl(arr1,arr2)

      implicit none

c     -------------------------------------------------------------------
c     Routine to transfer ghost cell values form arr1 to arr2
c     -------------------------------------------------------------------

c     Input variables

      real(8) :: arr1(0:,0:,0:)
      real(8) :: arr2(0:,0:,0:)

c     Local variables

      integer :: nxp,nyp,nzp

c     Begin program

      nxp = size(arr1,1)-1
      nyp = size(arr1,2)-1
      nzp = size(arr1,3)-1

      !X boundaries
      arr2(0  ,:,:) = arr1(0  ,:,:)
      arr2(nxp,:,:) = arr1(nxp,:,:)

      !Y boundaries
      arr2(:,0  ,:) = arr1(:,0  ,:)
      arr2(:,nyp,:) = arr1(:,nyp,:)

      !Z boundaries
      arr2(:,:,0  ) = arr1(:,:,0  )
      arr2(:,:,nzp) = arr1(:,:,nzp)

c     End program

      end subroutine XferBCs_scl

c     XferBCs_vec
c     ###################################################################
      subroutine XferBCs_vec(vec1,vec2)

      implicit none
c     -------------------------------------------------------------------
c     Routine to transfer ghost cell values form vec1 to vec2
c     -------------------------------------------------------------------

c     Input variables

      real(8) :: vec1(0:,0:,0:,:)
      real(8) :: vec2(0:,0:,0:,:)

c     Local variables

      integer :: nxp,nyp,nzp

c     Begin program

      nxp = size(vec1,1)-1
      nyp = size(vec1,2)-1
      nzp = size(vec1,3)-1

      !X boundaries
      vec2(0  ,:,:,:) = vec1(0  ,:,:,:)
      vec2(nxp,:,:,:) = vec1(nxp,:,:,:)

      !Y boundaries
      vec2(:,0  ,:,:) = vec1(:,0  ,:,:)
      vec2(:,nyp,:,:) = vec1(:,nyp,:,:)

      !Z boundaries
      vec2(:,:,0  ,:) = vec1(:,:,0  ,:)
      vec2(:,:,nzp,:) = vec1(:,:,nzp,:)

c     End program

      end subroutine XferBCs_vec

c     setBCsched
c     ###################################################################
      subroutine setBCsched(itime)

c     -------------------------------------------------------------------
c     Sets BC schedule for dependent and auxiliary variables, which
c     are identified by a 2-dim integer array, bc_seq. The first
c     dimension is of size nbc_seq, and indicates order of treatment
c     lexicographically (i.e., first component of first dimension is
c     to be treated first, second component of second dimension is
c     to be treated second, and so on).
c
c     The second dimension is of size two, and has the following convention:
c       * First component of the second dimension identifies variable,
c         and its sign indicates whether the variable is dependent (>0) or
c         auxiliary (<0).
c       * Second component of the second dimension identifies scalars (=0)
c         or vectors (=1).
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: itime

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
        gv%bc_grp(4+offset)%nbc_seq = 3
        gv%bc_grp(5+offset)%nbc_seq = 2
        gv%bc_grp(6+offset)%nbc_seq = 2
        gv%bc_grp(7+offset)%nbc_seq = 5
        gv%bc_grp(8+offset)%nbc_seq = 4
        gv%bc_grp(9+offset)%nbc_seq = 2
#else
        gv%bc_grp(1+offset)%nbc_seq = 1
        gv%bc_grp(2+offset)%nbc_seq = 1
        gv%bc_grp(3+offset)%nbc_seq = 2
        gv%bc_grp(4+offset)%nbc_seq = 3
        gv%bc_grp(5+offset)%nbc_seq = 2
        gv%bc_grp(6+offset)%nbc_seq = 5
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
        gv%bc_grp(1+offset)%bc_seq(1,1) = IRHO        !Dep. Variable
        gv%bc_grp(1+offset)%bc_seq(1,2) = SCALAR      !Scalar

        !Group 2    
        gv%bc_grp(2+offset)%bc_seq(1,1) = ITMP        !Dep. variable
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
        gv%bc_grp(4+offset)%bc_seq(3,2) = VECTOR      !Vector

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
        gv%bc_grp(7+offset)%bc_seq(5,2) = VECTOR      !Vector

        !Group 8
        gv%bc_grp(8+offset)%bc_seq(1,1) = IVECNV*AUX  !Aux. variable
        gv%bc_grp(8+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(8+offset)%bc_seq(2,1) = IVECOV*AUX  !Aux. variable
        gv%bc_grp(8+offset)%bc_seq(2,2) = VECTOR      !Vector

        gv%bc_grp(8+offset)%bc_seq(3,1) = IVECOV_N*AUX!Aux. variable
        gv%bc_grp(8+offset)%bc_seq(3,2) = VECTOR      !Vector

cc        gv%bc_grp(8+offset)%bc_seq(4,1) = IVSCNV*AUX  !Aux. variable
cc        gv%bc_grp(8+offset)%bc_seq(4,2) = VECTOR      !Vector

        gv%bc_grp(8+offset)%bc_seq(4,1) = IVEFCNV*AUX !Aux. variable
        gv%bc_grp(8+offset)%bc_seq(4,2) = VECTOR      !Vector

        !Group 9
        gv%bc_grp(9+offset)%bc_seq(1,1) = IEH*AUX    !Aux. variable
        gv%bc_grp(9+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(9+offset)%bc_seq(2,1) = IENI*AUX    !Aux. variable
        gv%bc_grp(9+offset)%bc_seq(2,2) = VECTOR      !Vector
#else
        !Group 1
        gv%bc_grp(1+offset)%bc_seq(1,1) = IRHO        !Dep. Variable
        gv%bc_grp(1+offset)%bc_seq(1,2) = SCALAR      !Scalar

        !Group 2    
        gv%bc_grp(2+offset)%bc_seq(1,1) = ITMP        !Dep. variable
        gv%bc_grp(2+offset)%bc_seq(1,2) = SCALAR      !Scalar

        !Group 3
        gv%bc_grp(3+offset)%bc_seq(1,1) = IETA*AUX    !Aux. variable
        gv%bc_grp(3+offset)%bc_seq(1,2) = SCALAR      !Scalar

        gv%bc_grp(3+offset)%bc_seq(2,1) = INU*AUX    !Aux. variable
        gv%bc_grp(3+offset)%bc_seq(2,2) = SCALAR      !Scalar

        !Group 4
cc        gv%bc_grp(4+offset)%bc_seq(1,1) = IBCNV_N*AUX !Aux. variable
cc        gv%bc_grp(4+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(4+offset)%bc_seq(1,1) = IBCNV*AUX   !Aux. variable
        gv%bc_grp(4+offset)%bc_seq(1,2) = VECTOR      !Vector

        gv%bc_grp(4+offset)%bc_seq(2,1) = IBCOV*AUX   !Aux. variable
        gv%bc_grp(4+offset)%bc_seq(2,2) = VECTOR      !Vector

        gv%bc_grp(4+offset)%bc_seq(3,1) = IBX         !Dep. variable
        gv%bc_grp(4+offset)%bc_seq(3,2) = VECTOR      !Vector

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
        gv%bc_grp(6+offset)%bc_seq(5,2) = VECTOR      !Vector

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

cc          gv%bc_grp(8+offset)%bc_seq(4,1) = IVSCNV*AUX  !Aux. variable
cc          gv%bc_grp(8+offset)%bc_seq(4,2) = VECTOR      !Vector

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

cc          gv%bc_grp(7+offset)%bc_seq(4,1) = IVSCNV*AUX  !Aux. variable
cc          gv%bc_grp(7+offset)%bc_seq(4,2) = VECTOR      !Vector

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

            grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,he=1)

            gsuper = 0.5*(gmetric%grid(igrid)%gsup(i+1,j,k,:,:)
     .                   +gmetric%grid(igrid)%gsup(i  ,j,k,:,:))

            grad_p = matmul(gsuper,grad_p)

            !ExB velocity
            v0(i+1,j,k,1)=(E0(2)*bbcov(3)-E0(3)*bbcov(2)
     .                    -etah*grad_p(1))/b2
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

      integer :: i,j,k
      real(8) :: b2,bbcnv(3),bbcov(3),grad_p(3),alpha,jach
     .          ,gsub(3,3),ietah
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
cc        fo = .true.
      endif

      i = nx

      if (.not.fo) then
        do k=1,nz
          do j=1,ny
cc            jach  = 0.5*(gmetric%grid(igrid)%jac(i  ,j,k)
cc     .                  +gmetric%grid(igrid)%jac(i+1,j,k))
cc            bbcnv = 0.5*(bcnv(i+1,j,k,:)+bcnv(i,j,k,:))/jach
            bbcnv = 0.5*
     .             (bcnv(i+1,j,k,:)/gmetric%grid(igrid)%jac(i+1,j,k)
     .             +bcnv(i  ,j,k,:)/gmetric%grid(igrid)%jac(i  ,j,k))
            bbcov = 0.5*(matmul(gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .                         ,bcnv(i  ,j,k,:))
     .                  +matmul(gmetric%grid(igrid)%gsub(i+1,j,k,:,:)
     .                         ,bcnv(i+1,j,k,:)))

            b2 = dot_product(bbcov,bbcnv)

            if (b2 > 0d0) then
cc              ietah = 2./(eta(i+1,j,k)+eta(i,j,k))
              ietah =0.5*(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

              grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,he=1)

              alpha = dot_product(E0,bbcnv)*ietah

              j0(i+1,j,k,:)=(alpha*bbcnv+cross_product(bbcov,grad_p))/b2
            endif
          enddo
        enddo
      else
        do k=1,nz
          do j=1,ny
cc            jach = gmetric%grid(igrid)%jac(i,j,k)
            jach  = 0.5*(gmetric%grid(igrid)%jac(i  ,j,k)
     .                  +gmetric%grid(igrid)%jac(i+1,j,k))
            bbcnv = bcnv(i,j,k,:)/jach

cc            gsub = gmetric%grid(igrid)%gsub(i,j,k,:,:)
            gsub = 0.5*(gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .                 +gmetric%grid(igrid)%gsub(i+1,j,k,:,:))
            bbcov = matmul(gsub,bcnv(i,j,k,:))

            b2 = dot_product(bbcov,bbcnv)
            if (b2 > 0d0) then
cc              ietah = 2./(eta(i+1,j,k)+eta(i,j,k))
              ietah =0.5*(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

              grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,he=1)

              alpha = dot_product(E0,bbcnv)*ietah

              j0(i+1,j,k,:)=(alpha*bbcnv+cross_product(bbcov,grad_p))/b2
            endif
          enddo
        enddo
      endif

c     End program

      end subroutine J_pinch_bc

      end module local_BCS_variables

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

        call extrapolateBC(array(:,:,:,ivar),array0(:,:,:,ivar)
cc     .                    ,dim,loc,order)
     .                    ,dim,loc,1)  !Required for mass conservation

      case (IBX,IBY,IBZ,IJX,IJY,IJZ)

        if (ivar /= dim .or. bcSYM(ibc)) then

          call extrapolateBC(array (:,:,:,ivar),array0(:,:,:,ivar)
     .                      ,dim,loc,order)

        else  !Imposes divergence-free constraint on B, J

          call findBCLoopLimits(dim,loc
     .                         ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
cc     .                         ,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                         ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

          select case (ibc)
          case (1)
 
            i = 1

            array(i-1,:,:,dim) = array(i+1,:,:,dim)  !To cancel normal derivative

            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = 2.*grid_params%dxh(ig)
                dh(2) = 2.*grid_params%dyh(jg)
                dh(3) = 2.*grid_params%dzh(kg)

                hey = 0
                if (j < 1   ) hey =  1
                if (j > nybc) hey = -1

                hez = 0
                if (k < 1   ) hez =  1
                if (k > nzbc) hez = -1

                diver=jac0
     .               *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                   ,hey=hey,hez=hez)

                array(i-1,j,k,dim) = array(i+1,j,k,dim) + dh(dim)*diver

              enddo
            enddo
 
          case (2)

            i = nxbc

            array(i+1,:,:,dim) = array(i-1,:,:,dim)  !To cancel normal derivative

            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = 2.*grid_params%dxh(ig)
                dh(2) = 2.*grid_params%dyh(jg)
                dh(3) = 2.*grid_params%dzh(kg)

                hey = 0
                if (j < 1   ) hey =  1
                if (j > nybc) hey = -1

                hez = 0
                if (k < 1   ) hez =  1
                if (k > nzbc) hez = -1

                diver =jac0
     .                *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                    ,hey=hey,hez=hez)

                array(i+1,j,k,dim) = array(i-1,j,k,dim) - dh(dim)*diver

              enddo
            enddo

          case (3)

            j = 1

            array(:,j-1,:,dim) = array(:,j+1,:,dim)  !To cancel normal derivative

            do i=imin,imax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = 2.*grid_params%dxh(ig)
                dh(2) = 2.*grid_params%dyh(jg)
                dh(3) = 2.*grid_params%dzh(kg)

                hex = 0
                if (i < 1   ) hex =  1
                if (i > nxbc) hex = -1

                hez = 0
                if (k < 1   ) hez =  1
                if (k > nzbc) hez = -1

                diver =jac0
     .                *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                    ,hex=hex,hez=hez)

                array(i,j-1,k,dim) = array(i,j+1,k,dim) + dh(dim)*diver
                  
              enddo
            enddo

          case (4)

            j = nybc

            array(:,j+1,:,dim) = array(:,j-1,:,dim)  !To cancel normal derivative

            do i=imin,imax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = 2.*grid_params%dxh(ig)
                dh(2) = 2.*grid_params%dyh(jg)
                dh(3) = 2.*grid_params%dzh(kg)

                hex = 0
                if (i < 1   ) hex =  1
                if (i > nxbc) hex = -1

                hez = 0
                if (k < 1   ) hez =  1
                if (k > nzbc) hez = -1

                diver =jac0
     .                *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                    ,hex=hex,hez=hez)

                array(i,j+1,k,dim) = array(i,j-1,k,dim) - dh(dim)*diver
                  
              enddo
            enddo

          case (5)

            k = 1

            array(:,:,k-1,dim) = array(:,:,k+1,dim)   !To cancel normal derivative

            do i=imin,imax
              do j=jmin,jmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = 2.*grid_params%dxh(ig)
                dh(2) = 2.*grid_params%dyh(jg)
                dh(3) = 2.*grid_params%dzh(kg)

                hex = 0
                if (i < 1   ) hex =  1
                if (i > nxbc) hex = -1

                hey = 0
                if (j < 1   ) hey =  1
                if (j > nybc) hey = -1

                diver =jac0
     .                *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                    ,hex=hex,hey=hey)

                array(i,j,k-1,dim) = array(i,j,k+1,dim) + dh(dim)*diver

              enddo
            enddo

          case (6)

            k = nzbc

            array(:,:,k+1,dim) = array(:,:,k-1,dim)  !To cancel normal derivative

            do i=imin,imax
              do j=jmin,jmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = 2.*grid_params%dxh(ig)
                dh(2) = 2.*grid_params%dyh(jg)
                dh(3) = 2.*grid_params%dzh(kg)

                hex = 0
                if (i < 1   ) hex =  1
                if (i > nxbc) hex = -1

                hey = 0
                if (j < 1   ) hey =  1
                if (j > nybc) hey = -1

                diver =jac0
     .                *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                    ,hex=hex,hey=hey)

                array(i,j,k+1,dim) = array(i,j,k-1,dim) - dh(dim)*diver

              enddo
            enddo

          case default

            messg = 'Boundary'//int2char(ibc)//' non existent'
            call pstop('vectorDirichletBC',messg)

          end select

        endif

      case (IAX,IAY,IAZ)

        if (ivar /= dim) then

          call extrapolateBC(array (:,:,:,ivar),array0(:,:,:,ivar)
     .                      ,dim,loc,order)

        else  !Imposes divergence-free constraint on B, J

          call findBCLoopLimits(dim,loc
cc     .                         ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                         ,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                         ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

          select case (ibc)
          case (1)
 
            i = 1

            array(i-1,:,:,dim) = array(i+1,:,:,dim)  !To cancel normal derivative

            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = grid_params%dxh(ig)  !<=== Different from above
                dh(2) = grid_params%dyh(jg)
                dh(3) = grid_params%dzh(kg)

                hey = 0
                if (j < 1   ) hey =  1
                if (j > nybc) hey = -1

                hez = 0
                if (k < 1   ) hez =  1
                if (k > nzbc) hez = -1

                diver=jac0
     .               *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                   ,hey=hey,hez=hez)

                array(i-1,j,k,dim) = array(i,j,k,dim) + dh(dim)*diver  !<=== Different from above

              enddo
            enddo
 
          case (2)

            i = nxbc

            array(i+1,:,:,dim) = array(i-1,:,:,dim)  !To cancel normal derivative

            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = grid_params%dxh(ig)
                dh(2) = grid_params%dyh(jg)
                dh(3) = grid_params%dzh(kg)

                hey = 0
                if (j < 1   ) hey =  1
                if (j > nybc) hey = -1

                hez = 0
                if (k < 1   ) hez =  1
                if (k > nzbc) hez = -1

                diver =jac0
     .                *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                    ,hey=hey,hez=hez)

                array(i+1,j,k,dim) = array(i,j,k,dim) - dh(dim)*diver

              enddo
            enddo

          case (3)

            j = 1

            array(:,j-1,:,dim) = array(:,j+1,:,dim)  !To cancel normal derivative

            do i=imin,imax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = grid_params%dxh(ig)
                dh(2) = grid_params%dyh(jg)
                dh(3) = grid_params%dzh(kg)

                hex = 0
                if (i < 1   ) hex =  1
                if (i > nxbc) hex = -1

                hez = 0
                if (k < 1   ) hez =  1
                if (k > nzbc) hez = -1

                diver =jac0
     .                *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                    ,hex=hex,hez=hez)

                array(i,j-1,k,dim) = array(i,j,k,dim) + dh(dim)*diver
                  
              enddo
            enddo

          case (4)

            j = nybc

            array(:,j+1,:,dim) = array(:,j-1,:,dim)  !To cancel normal derivative

            do i=imin,imax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = grid_params%dxh(ig)
                dh(2) = grid_params%dyh(jg)
                dh(3) = grid_params%dzh(kg)

                hex = 0
                if (i < 1   ) hex =  1
                if (i > nxbc) hex = -1

                hez = 0
                if (k < 1   ) hez =  1
                if (k > nzbc) hez = -1

                diver =jac0
     .                *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                    ,hex=hex,hez=hez)

                array(i,j+1,k,dim) = array(i,j,k,dim) - dh(dim)*diver
                  
              enddo
            enddo

          case (5)

            k = 1

            array(:,:,k-1,dim) = array(:,:,k+1,dim)   !To cancel normal derivative

            do i=imin,imax
              do j=jmin,jmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = grid_params%dxh(ig)
                dh(2) = grid_params%dyh(jg)
                dh(3) = grid_params%dzh(kg)

                hex = 0
                if (i < 1   ) hex =  1
                if (i > nxbc) hex = -1

                hey = 0
                if (j < 1   ) hey =  1
                if (j > nybc) hey = -1

                diver =jac0
     .                *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                    ,hex=hex,hey=hey)

                array(i,j,k-1,dim) = array(i,j,k,dim) + dh(dim)*diver

              enddo
            enddo

          case (6)

            k = nzbc

            array(:,:,k+1,dim) = array(:,:,k-1,dim)  !To cancel normal derivative

            do i=imin,imax
              do j=jmin,jmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
                jac0 = gmetric%grid(igxbc)%jac(i,j,k)

                dh(1) = grid_params%dxh(ig)
                dh(2) = grid_params%dyh(jg)
                dh(3) = grid_params%dzh(kg)

                hex = 0
                if (i < 1   ) hex =  1
                if (i > nxbc) hex = -1

                hey = 0
                if (j < 1   ) hey =  1
                if (j > nybc) hey = -1

                diver =jac0
     .                *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array
     .                    ,hex=hex,hey=hey)

                array(i,j,k+1,dim) = array(i,j,k,dim) - dh(dim)*diver

              enddo
            enddo

          case default

            messg = 'Boundary'//int2char(ibc)//' non existent'
            call pstop('vectorDirichletBC',messg)

          end select

        endif

cc        if (ivar /= dim) then
cc
cc          call extrapolateBC(array (:,:,:,ivar),array0(:,:,:,ivar)
cc     .                      ,dim,loc,order)
cc
cc        else !Imposes divergence-free constraint on A (not exactly the same as above!!!!)
cc
cc          call findBCLoopLimits(dim,loc
cc     .                         ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
cc     .                         ,imin ,imax ,jmin ,jmax ,kmin ,kmax)
cc
cc
cc          select case(dim)
cc          case(1)
cc            allocate(rhs(0:nybc+1,0:nzbc+1))
cc          case(2)
cc            allocate(rhs(0:nxbc+1,0:nzbc+1))
cc          case(3)
cc            allocate(rhs(0:nxbc+1,0:nybc+1))
cc          end select
cc
cc          rhs = 0d0
cc
cc          do i=imin,imax
cc            do j=jmin,jmax
cc              do k=kmin,kmax
cc
cc                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
cc                jac0 = gmetric%grid(igxbc)%jac(i,j,k)
cc
cc                dh(1) = grid_params%dxh(ig)  !<=== Different from above
cc                dh(2) = grid_params%dyh(jg)
cc                dh(3) = grid_params%dzh(kg)
cc
cc                select case (ibc)
cc                case (1)
cc                  array(i-1,j,k,dim) = array(i+1,j,k,dim)
cc                  diver =jac0
cc     .               *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)
cc                  rhs(j,k) = array(i,j,k,dim) + dh(dim)*diver !<=== Different from above
cc                case (2)
cc                  array(i+1,j,k,dim) = array(i-1,j,k,dim)
cc                  diver =jac0
cc     .               *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)
cc                  rhs(j,k) = array(i,j,k,dim) - dh(dim)*diver
cc                case (3)
cc                  array(i,j-1,k,dim) = array(i,j+1,k,dim)
cc                  diver =jac0
cc     .               *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)
cc                  rhs(i,k) = array(i,j,k,dim) + dh(dim)*diver
cc                case (4)
cc                  array(i,j+1,k,dim) = array(i,j-1,k,dim)
cc                  diver =jac0
cc     .               *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)
cc                  rhs(i,k) = array(i,j,k,dim) - dh(dim)*diver
cc                case (5)
cc                  array(i,j,k-1,dim) = array(i,j,k+1,dim)
cc                  diver =jac0
cc     .               *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)
cc                  rhs(i,j) = array(i,j,k,dim) + dh(dim)*diver
cc                case (6)
cc                  array(i,j,k+1,dim) = array(i,j,k-1,dim)
cc                  diver =jac0
cc     .               *div(i,j,k,nxbc,nybc,nzbc,igxbc,igybc,igzbc,array)
cc                  rhs(i,j) = array(i,j,k,dim) - dh(dim)*diver
cc                end select
cc
cc              enddo
cc            enddo
cc          enddo
cc
ccc       Update BC ghost nodes
cc
cc          select case (ibc)
cc          case (1)              !x0
cc          array(0     ,jmin:jmax,kmin:kmax,dim)=rhs(jmin:jmax,kmin:kmax)
cc          case (2)              !x1
cc          array(nxbc+1,jmin:jmax,kmin:kmax,dim)=rhs(jmin:jmax,kmin:kmax)
cc          case (3)              !y0
cc          array(imin:imax,0     ,kmin:kmax,dim)=rhs(imin:imax,kmin:kmax)
cc          case (4)              !y1
cc          array(imin:imax,nybc+1,kmin:kmax,dim)=rhs(imin:imax,kmin:kmax)
cc          case (5)              !z0
cc          array(imin:imax,jmin:jmax,0     ,dim)=rhs(imin:imax,jmin:jmax)
cc          case (6)              !z1
cc          array(imin:imax,jmin:jmax,nzbc+1,dim)=rhs(imin:imax,jmin:jmax)
cc          case default
cc            messg = 'Boundary'//int2char(ibc)//' non existent'
cc            call pstop('vectorDirichletBC',messg)
cc          end select
cc
cc          deallocate(rhs)
cc
cc        endif

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
     .          ,igp,igm,jgp,jgm,kgp,kgm,i1,j1,k1
      real(8) :: x1,x2,x3,dh(3),jac0,ee,alp,vec(3),cnv(3)
     .          ,cov(3),J0(3),bnorm,jn,etal,B00(3)
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
      case (IVX,IVY,IVZ,IJX,IJY,IJZ,IAX,IAY,IAZ) !Velocity components, J and A components (for symmetry BCs)

        if (coords /= 'car') then

          call findBCLoopLimits(dim,loc
     .                         ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                         ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

                gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)

                hessian = gmetric%grid(igxbc)%Gamma(i,j,k,:,:,:)

cc                ip = min(i+1,nxbc)
cc                im = max(i-1,1   )
cc                jp = min(j+1,nybc)
cc                jm = max(j-1,1   )
cc                kp = min(k+1,nzbc)
cc                km = max(k-1,1   )
cc
cc                dh(1) = 2.*grid_params%dxh(ig)
cc                if (i == nxbc) dh(1) = grid_params%dx(ig-1)
cc                if (i == 1   ) dh(1) = grid_params%dx(ig)
cc
cc                dh(2) = 2.*grid_params%dyh(jg)
cc                if (j == nybc) dh(2) = grid_params%dy(jg-1)
cc                if (j == 1   ) dh(2) = grid_params%dy(jg)
cc
cc                dh(3) = 2.*grid_params%dzh(kg)
cc                if (k == nzbc) dh(3) = grid_params%dz(kg-1)
cc                if (k == 1   ) dh(3) = grid_params%dz(kg)

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

      case (IBX,IBY,IBZ) 

        if (ivar /= dim) then !Not symmetry BC

          call findBCLoopLimits(dim,loc
     .                         ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                         ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

                gsuper = gmetric%grid(igxbc)%gsup(i  ,j,k,:,:)
                gsupip = gmetric%grid(igxbc)%gsup(i+1,j,k,:,:)
                gsupim = gmetric%grid(igxbc)%gsup(i-1,j,k,:,:)
                gsupjp = gmetric%grid(igxbc)%gsup(i,j+1,k,:,:)
                gsupjm = gmetric%grid(igxbc)%gsup(i,j-1,k,:,:)
                gsupkp = gmetric%grid(igxbc)%gsup(i,j,k+1,:,:)
                gsupkm = gmetric%grid(igxbc)%gsup(i,j,k-1,:,:)

cc                ip = min(i+1,nxbc)
cc                im = max(i-1,1   )
cc                jp = min(j+1,nybc)
cc                jm = max(j-1,1   )
cc                kp = min(k+1,nzbc)
cc                km = max(k-1,1   )
cc
cc                dh(1) = 2.*grid_params%dxh(ig)
cc                if (i == nxbc) dh(1) = grid_params%dx(ig-1)
cc                if (i == 1   ) dh(1) = grid_params%dx(ig)
cc
cc                dh(2) = 2.*grid_params%dyh(jg)
cc                if (j == nybc) dh(2) = grid_params%dy(jg-1)
cc                if (j == 1   ) dh(2) = grid_params%dy(jg)
cc
cc                dh(3) = 2.*grid_params%dzh(kg)
cc                if (k == nzbc) dh(3) = grid_params%dz(kg-1)
cc                if (k == 1   ) dh(3) = grid_params%dz(kg)

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

                  !!!!!!!Finite flow at wall!!!!!!!
                  if (pinch_flow) then

                    ee = array0(i+1,j,k,icomp)

                  !!!!!!!Finite current at wall!!!!!!!
                  else
                    if (ibc == 1) then
                      ee = 0.5*(array0(i-1,j,k,icomp)
     .                         +array0(i  ,j,k,icomp))
                    elseif (ibc == 2) then
                      ee = 0.5*(array0(i+1,j,k,icomp)
     .                         +array0(i  ,j,k,icomp))
                    endif
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
                    ee = 0.5*(array0(i,j-1,k,icomp)
     .                       +array0(i,j  ,k,icomp))
                  elseif (ibc == 4) then
                    ee = 0.5*(array0(i,j+1,k,icomp)
     .                       +array0(i,j  ,k,icomp))
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
                    ee = 0.5*(array0(i,j,k-1,icomp)
     .                       +array0(i,j,k  ,icomp))
                  elseif (ibc == 6) then
                    ee = 0.5*(array0(i,j,k+1,icomp)
     .                       +array0(i,j,k  ,icomp))
                  endif

                  rhs(i,j) = dh(dim)*(rhs(i,j) + isgn*ee)

                endif

              enddo
            enddo
          enddo

        endif

      case default

        write (*,*) 'Variable=',ieq
        call pstop('vectorNeumanBC','Variable not implemented')

      end select

c     Assign value

      select case (ibc)
      case (1)
        array(0     ,jmin:jmax,kmin:kmax,ivar)
     $                    = array(1+offset,jmin:jmax,kmin:kmax,ivar)
     $                    - rhs(jmin:jmax,kmin:kmax)
      case (2)
        array(nxbc+1,jmin:jmax,kmin:kmax,ivar)
     $                    = array(nxbc-offset,jmin:jmax,kmin:kmax,ivar)
     $                    + rhs(jmin:jmax,kmin:kmax)
      case (3)
        array(imin:imax,0     ,kmin:kmax,ivar)
     $                    = array(imin:imax,1+offset,kmin:kmax,ivar)
     $                    - rhs(imin:imax,kmin:kmax)
      case (4)
        array(imin:imax,nybc+1,kmin:kmax,ivar)
     $                    = array(imin:imax,nybc-offset,kmin:kmax,ivar)
     $                    + rhs(imin:imax,kmin:kmax)
      case (5)
        array(imin:imax,jmin:jmax,0     ,ivar)
     $                    = array(imin:imax,jmin:jmax,1+offset,ivar)
     $                    - rhs(imin:imax,jmin:jmax)
      case (6)
        array(imin:imax,jmin:jmax,nzbc+1,ivar)
     $                    = array(imin:imax,jmin:jmax,nzbc-offset,ivar)
     $                    + rhs(imin:imax,jmin:jmax)
      end select

      deallocate(rhs)

c     End program

      end subroutine vectorNeumannBC
