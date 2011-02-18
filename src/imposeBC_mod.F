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

          varray%array_var(ieq)%bc_dep_list(2,1) = IVCOV*AUX
          varray%array_var(ieq)%bc_dep_list(2,2) = VECTOR
          varray%array_var(ieq)%bc_dep_list(2,3) = 0 !trivial
        enddo
      endif

      !B
      do ieq=IBX,IBZ
        varray%array_var(ieq)%bc_dep_list(1,1) = IBCOV*AUX
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

c     XferBCs_scl
c     ###################################################################
      subroutine XferBCs_scl(arr1,arr2,average)

      implicit none

c     -------------------------------------------------------------------
c     Routine to transfer ghost cell values form arr1 to arr2
c     -------------------------------------------------------------------

c     Input variables

      real(8) :: arr1(0:,0:,0:)
      real(8) :: arr2(0:,0:,0:)
      logical :: average

c     Local variables

      integer :: nxp,nyp,nzp

c     Begin program

      nxp = size(arr1,1)-1
      nyp = size(arr1,2)-1
      nzp = size(arr1,3)-1

      if (average) then
        arr2(nxp,:,:) = 0.5*(arr1(nxp  ,:,:)
     .                      +arr1(nxp-1,:,:))
        arr2(0  ,:,:) = 0.5*(arr1(1    ,:,:)
     .                      +arr1(0    ,:,:))

        arr2(:,nyp,:) = 0.5*(arr1(:,nyp  ,:)
     .                      +arr1(:,nyp-1,:))
        arr2(:,0  ,:) = 0.5*(arr1(:,1    ,:)
     .                      +arr1(:,0    ,:))

        arr2(:,:,nzp) = 0.5*(arr1(:,:,nzp  )
     .                      +arr1(:,:,nzp-1))
        arr2(:,:,0  ) = 0.5*(arr1(:,:,1    )
     .                      +arr1(:,:,0    ))
      else
        !X boundaries
        arr2(0  ,:,:) = arr1(0  ,:,:)
        arr2(nxp,:,:) = arr1(nxp,:,:)

        !Y boundaries
        arr2(:,0  ,:) = arr1(:,0  ,:)
        arr2(:,nyp,:) = arr1(:,nyp,:)

        !Z boundaries
        arr2(:,:,0  ) = arr1(:,:,0  )
        arr2(:,:,nzp) = arr1(:,:,nzp)
      endif

c     End program

      end subroutine XferBCs_scl

c     XferBCs_vec
c     ###################################################################
      subroutine XferBCs_vec(vec1,vec2,average)

      implicit none
c     -------------------------------------------------------------------
c     Routine to transfer ghost cell values form vec1 to vec2
c     -------------------------------------------------------------------

c     Input variables

      real(8) :: vec1(0:,0:,0:,:)
      real(8) :: vec2(0:,0:,0:,:)
      logical :: average

c     Local variables

      integer :: nxp,nyp,nzp

c     Begin program

      nxp = size(vec1,1)-1
      nyp = size(vec1,2)-1
      nzp = size(vec1,3)-1

      if (average) then
        vec2(nxp,:,:,:) = 0.5*(vec1(nxp  ,:,:,:)
     .                        +vec1(nxp-1,:,:,:))
        vec2(0  ,:,:,:) = 0.5*(vec1(1    ,:,:,:)
     .                        +vec1(0    ,:,:,:))

        vec2(:,nyp,:,:) = 0.5*(vec1(:,nyp  ,:,:)
     .                        +vec1(:,nyp-1,:,:))
        vec2(:,0  ,:,:) = 0.5*(vec1(:,1    ,:,:)
     .                        +vec1(:,0    ,:,:))

        vec2(:,:,nzp,:) = 0.5*(vec1(:,:,nzp  ,:)
     .                        +vec1(:,:,nzp-1,:))
        vec2(:,:,0  ,:) = 0.5*(vec1(:,:,1    ,:)
     .                        +vec1(:,:,0    ,:))
      else
        !X boundaries
        vec2(0  ,:,:,:) = vec1(0  ,:,:,:)
        vec2(nxp,:,:,:) = vec1(nxp,:,:,:)

        !Y boundaries
        vec2(:,0  ,:,:) = vec1(:,0  ,:,:)
        vec2(:,nyp,:,:) = vec1(:,nyp,:,:)

        !Z boundaries
        vec2(:,:,0  ,:) = vec1(:,:,0  ,:)
        vec2(:,:,nzp,:) = vec1(:,:,nzp,:)
      endif

c     End program

      end subroutine XferBCs_vec

c     setBCsched
c     ###################################################################
      subroutine setBCsched(varray,itime)

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

c     Local variables

      integer :: ierr,offset=0,igr

#if 1
      INTERFACE
        subroutine create_bc_schedule(u,aux,nbc_grp,bc_grp,
     &            inherit,interior_only,print_used,reverse)

        use variable_setup, only: var_array
        use variable_setup, only: aux_array
        use variable_setup, only: bcgrp
        type(var_array) :: u
        type(aux_array) :: aux
        logical :: inherit
        logical :: interior_only
        logical :: print_used
        logical :: reverse
        integer :: nbc_grp
        type(bcgrp), pointer :: bc_grp(:)
        end subroutine create_bc_schedule
      END INTERFACE

c     Begin program

#if defined(samrai)
      call create_bc_schedule(varray,gv%aux,gv%nbc_grp,gv%bc_grp
     .                       ,.false.,.false.,.false.,.false.)
#else
      call create_bc_schedule(varray,gv%aux,gv%nbc_grp,gv%bc_grp
     .                       ,.false.,.false.,.false.,.false.)
cc     .                       ,.false.,.false.,(itime==0),.false.)
#endif

#else

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
        gv%bc_grp(9+offset)%nbc_seq = 2
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

#endif

c     End program

      end subroutine setBCsched

ccc     setBCsched
ccc     ###################################################################
cc      subroutine setBCsched(varray,itime)
cc
ccc     -------------------------------------------------------------------
ccc     Sets BC schedule for dependent and auxiliary variables, which
ccc     are identified by a 2-dim integer array, bc_seq. The first
ccc     dimension is of size nbc_seq, and indicates order of treatment
ccc     lexicographically (i.e., first component of first dimension is
ccc     to be treated first, second component of first dimension is
ccc     to be treated second, and so on).
ccc
ccc     The second dimension is of size two, and has the following convention:
ccc       * First component of the second dimension identifies variable,
ccc         and its sign indicates whether the variable is dependent (>0) or
ccc         auxiliary (<0).
ccc       * Second component of the second dimension identifies whether
ccc         variable is a scalar (=0) or a vector (=1).
ccc     -------------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer :: itime
cc
cc      type(var_array),pointer :: varray
cc
ccc     Local variables
cc
cc      integer :: ierr,offset=0,igr
cc
cc#if 1
cc      INTERFACE
cc        subroutine create_bc_schedule(u,aux,nbc_grp,bc_grp,
cc     &            inherit,interior_only,print_used,reverse)
cc
cc        use variable_setup, only: var_array
cc        use variable_setup, only: aux_array
cc        use variable_setup, only: bcgrp
cc        type(var_array) :: u
cc        type(aux_array) :: aux
cc        logical :: inherit
cc        logical :: interior_only
cc        logical :: print_used
cc        logical :: reverse
cc        integer :: nbc_grp
cc        type(bcgrp), pointer :: bc_grp(:)
cc        end subroutine create_bc_schedule
cc      END INTERFACE
cc
ccc     Begin program
cc
cc      call create_bc_schedule(varray,gv%aux,gv%nbc_grp,gv%bc_grp
cc     .                       ,.false.,.false.,(itime==0),.false.)
cc#else
cc
ccc     Set number of BC groups
cc
cc      gv%nbc_grp = 3
cc
ccc     Check BC groups
cc
cc      if (     associated(gv%bc_grp)
cc     .    .and.size(gv%bc_grp) /= gv%nbc_grp) then
cc
cc        do igr=1,size(gv%bc_grp)
cc          deallocate(gv%bc_grp(igr)%bc_seq)
cc        enddo
cc        deallocate(gv%bc_grp)
cc
cc      endif
cc
ccc     Setup schedule
cc
cc      if (.not.associated(gv%bc_grp)) then
cc
cc        allocate(gv%bc_grp(gv%nbc_grp),STAT=ierr)
cc
cc        !Define group sizes
cc        gv%bc_grp(1)%nbc_seq = 4
cc        gv%bc_grp(2)%nbc_seq = 17
cc        gv%bc_grp(3)%nbc_seq = 6
cc
cc        !Allocate groups
cc        do igr=1,gv%nbc_grp
cc          allocate(gv%bc_grp(igr)
cc     .               %bc_seq(gv%bc_grp(igr)%nbc_seq,2),STAT=ierr)
cc
cc          gv%bc_grp(igr)%bc_seq = 0
cc        enddo
cc
cc        !Set groups
cc
cc        !Group 1
cc        gv%bc_grp(1)%bc_seq(1,1) = IRHO_0*AUX        !Dep. Variable
cc        gv%bc_grp(1)%bc_seq(1,2) = SCALAR      !Scalar
cc
cc        gv%bc_grp(1)%bc_seq(2,1) = ITMP_0*AUX        !Dep. variable
cc        gv%bc_grp(1)%bc_seq(2,2) = SCALAR      !Scalar
cc
cc        gv%bc_grp(1)%bc_seq(3,1) = IBCNV_0*AUX  !Aux. variable
cc        gv%bc_grp(1)%bc_seq(3,2) = VECTOR       !Vector
cc
cc        gv%bc_grp(1)%bc_seq(4,1) = IBCOV_0*AUX  !Aux. variable
cc        gv%bc_grp(1)%bc_seq(4,2) = VECTOR       !Vector
cc
cc        gv%bc_grp(1)%bc_seq(5,1) = IRHO  !Aux. Variable
cc        gv%bc_grp(1)%bc_seq(5,2) = SCALAR       !Scalar
cc
cc        !Group 2
cc        gv%bc_grp(2)%bc_seq(1,1) = IRHO  !Aux. Variable
cc        gv%bc_grp(2)%bc_seq(1,2) = SCALAR       !Scalar
cc
cc        gv%bc_grp(2)%bc_seq(2,1) = ITMP   !Aux. variable
cc        gv%bc_grp(2)%bc_seq(2,2) = SCALAR       !Scalar
cc
cc        gv%bc_grp(2)%bc_seq(3,1) = IETA_0*AUX   !Aux. variable
cc        gv%bc_grp(2)%bc_seq(3,2) = SCALAR       !Scalar
cc
cc        gv%bc_grp(2)%bc_seq(4,1) = INU*AUX    !Aux. variable
cc        gv%bc_grp(2)%bc_seq(4,2) = SCALAR      !Scalar
cc
cc        gv%bc_grp(2)%bc_seq(5,1) = IJCNV_0*AUX  !Aux. variable
cc        gv%bc_grp(2)%bc_seq(5,2) = VECTOR       !Vector
cc
cc        gv%bc_grp(2)%bc_seq(6,1) = IJCOV_0*AUX  !Aux. variable
cc        gv%bc_grp(2)%bc_seq(6,2) = VECTOR       !Vector
cc
cc        gv%bc_grp(2)%bc_seq(7,1) = IVCNV_0*AUX  !Aux. variable
cc        gv%bc_grp(2)%bc_seq(7,2) = VECTOR       !Vector
cc
cc        gv%bc_grp(2)%bc_seq(8,1) = IBCNV*AUX   !Aux. variable
cc        gv%bc_grp(2)%bc_seq(8,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(2)%bc_seq(9,1) = IBCOV*AUX   !Aux. variable
cc        gv%bc_grp(2)%bc_seq(9,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(2)%bc_seq(10,1) = IBX         !Dep. variable
cc        gv%bc_grp(2)%bc_seq(10,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(2)%bc_seq(11,1) = IBY         !Dep. variable
cc        gv%bc_grp(2)%bc_seq(11,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(2)%bc_seq(12,1) = IBZ         !Dep. variable
cc        gv%bc_grp(2)%bc_seq(12,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(2)%bc_seq(13,1) = IVCNV*AUX   !Aux. variable
cc        gv%bc_grp(2)%bc_seq(13,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(2)%bc_seq(14,1) = IVCOV*AUX   !Aux. variable
cc        gv%bc_grp(2)%bc_seq(14,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(2)%bc_seq(15,1) = IVX         !Dep. variable
cc        gv%bc_grp(2)%bc_seq(15,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(2)%bc_seq(16,1) = IVY         !Dep. variable
cc        gv%bc_grp(2)%bc_seq(16,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(2)%bc_seq(17,1) = IVZ         !Dep. variable
cc        gv%bc_grp(2)%bc_seq(17,2) = VECTOR      !Vector
cc
cc        !Group 3
cc        gv%bc_grp(3)%bc_seq(1,1) = IETA*AUX    !Aux. variable
cc        gv%bc_grp(3)%bc_seq(1,2) = SCALAR      !Scalar
cc
cc        gv%bc_grp(3)%bc_seq(2,1) = IJCNV*AUX   !Aux. variable
cc        gv%bc_grp(3)%bc_seq(2,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(3)%bc_seq(3,1) = IJCOV*AUX   !Aux. variable
cc        gv%bc_grp(3)%bc_seq(3,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(3)%bc_seq(4,1) = IENI_0*AUX   !Aux. variable
cc        gv%bc_grp(3)%bc_seq(4,2) = VECTOR       !Vector
cc
cc        gv%bc_grp(3)%bc_seq(5,1) = IENI*AUX    !Aux. variable
cc        gv%bc_grp(3)%bc_seq(5,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(3)%bc_seq(6,1) = IDIVPI*AUX  !Aux. variable
cc        gv%bc_grp(3)%bc_seq(6,2) = VECTOR      !Vector
cc
cc      endif
cc
cc#endif
cc
cccc      do igr=1,gv%nbc_grp
cccc        write (*,*) 'Group=',igr
cccc        do ierr=1,gv%bc_grp(igr)%nbc_seq
cccc          write (*,*) '  Var =',gv%bc_grp(igr)%bc_seq(ierr,1)
cccc          write (*,*) '  Type=',gv%bc_grp(igr)%bc_seq(ierr,2)
cccc        enddo
cccc      enddo
cccc      stop
cc
ccc     End program
cc
cc      end subroutine setBCsched

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
cc                  if (pinch_flow) then
cc
cc                    ee = array0(i+1,j,k,icomp)
cc
cc                  !!!!!!!Finite current at wall!!!!!!!
cc                  else
                    if (ibc == 1) then
                      ee = array0(i-1,j,k,icomp)
cc                      ee = 0.5*(array0(i-1,j,k,icomp)
cc     .                         +array0(i  ,j,k,icomp))
                    elseif (ibc == 2) then
                      ee = array0(i+1,j,k,icomp)
cc                      ee = 0.5*(array0(i+1,j,k,icomp)
cc     .                         +array0(i  ,j,k,icomp))
                    endif
cc                  endif

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
cc                    ee = 0.5*(array0(i,j-1,k,icomp)
cc     .                       +array0(i,j  ,k,icomp))
                  elseif (ibc == 4) then
                    ee = array0(i,j+1,k,icomp)
cc                    ee = 0.5*(array0(i,j+1,k,icomp)
cc     .                       +array0(i,j  ,k,icomp))
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
cc                    ee = 0.5*(array0(i,j,k-1,icomp)
cc     .                       +array0(i,j,k  ,icomp))
                  elseif (ibc == 6) then
                    ee = array0(i,j,k+1,icomp)
cc                    ee = 0.5*(array0(i,j,k+1,icomp)
cc     .                       +array0(i,j,k  ,icomp))
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

c create_bc_schedule
c######################################################################
      subroutine create_bc_schedule(u,aux,nbc_grp,bc_grp,
     &   inherit,interior_only,print_used,reverse)

      use variable_setup, only: var_array,aux_array,bcgrp

      INTERFACE
         function get_dep_index(dep_var,n_var)
            INTEGER :: get_dep_index
            INTEGER, INTENT(IN) :: dep_var(3)
            INTEGER, INTENT(IN) :: n_var(4)
         end function get_dep_index
         function check_dep_list(N,dep_list)
            INTEGER :: check_dep_list
            INTEGER :: dep_list(N,3)
            INTEGER i, j, N
         end function check_dep_list
         subroutine print_var_name(u,aux,k)
            use variable_setup, only: var_array, aux_array
            type(var_array) :: u
            type(aux_array) :: aux
            integer k
         end subroutine print_var_name
       END INTERFACE
 
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
      
      end subroutine create_bc_schedule

c check_dep_list
c######################################################################
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


c get_dep_index
c######################################################################
      integer function get_dep_index(dep_var,n_var)
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


c print_var_name
c######################################################################
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
cc#if defined(samrai)
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
cc#endif

