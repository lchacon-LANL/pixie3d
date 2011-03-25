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
     .                             ,res_ijk,vis_ijk,solve_prs

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

      case ('ppnsl','p3nsl','ohtor')           !No slip BCs, pinch flow

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX ) = EQU

        !Dirichlet boundary conditions for tangential components
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY ) = DIR
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ ) = DIR

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

#if !defined(petsc)
!#if 0
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
     .                       ,.true.,.false.,.false.,.false.)
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

            grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,hex=1)

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

              grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,hex=1)

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

              grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,hex=1)

              alpha = dot_product(E0,bbcnv)*ietah

              j0(i+1,j,k,:)=(alpha*bbcnv+cross_product(bbcov,grad_p))/b2
            endif
          enddo
        enddo
      endif

c     End program

      end subroutine J_pinch_bc

c     res_bc
c     #############################################################
      function res_bc(nx,ny,nz,igrid,field) result(rsvty)
c     -------------------------------------------------------------
c     This function computes the resistivity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,igrid
      real(8) :: field(0:nx+1,0:ny+1,0:nz+1)
     .          ,rsvty(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer :: i,j,k,imin,imax,jmin,jmax,kmin,kmax,dim,loc

c     Begin program

      do dim=1,3
        do loc=0,1

          call findBCLoopLimits(dim,loc
     .                     ,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                     ,imin ,imax ,jmin ,jmax ,kmin ,kmax )

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax
                rsvty(i,j,k)=res_ijk(i,j,k,igrid,field(i,j,k))
              enddo
            enddo
          enddo

        enddo
      enddo

c     End program

      end function res_bc

c     vis_bc
c     #############################################################
      function vis_bc(nx,ny,nz,igrid,field) result(vscty)
c     -------------------------------------------------------------
c     This function computes the viscosity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,igrid
      real(8) :: field(0:nx+1,0:ny+1,0:nz+1)
     .          ,vscty(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer :: i,j,k,imin,imax,jmin,jmax,kmin,kmax,dim,loc

c     Begin program

      do dim=1,3
        do loc=0,1

          call findBCLoopLimits(dim,loc
     .                     ,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                     ,imin ,imax ,jmin ,jmax ,kmin ,kmax )

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax
                vscty(i,j,k)=vis_ijk(i,j,k,igrid,field(i,j,k))
              enddo
            enddo
          enddo

        enddo
      enddo

c     End program

      end function vis_bc

      end module local_BCS_variables

