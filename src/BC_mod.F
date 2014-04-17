ccc     setBCsched_app
ccc     ###################################################################
cc      subroutine setBCsched_app
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
cc      use variables
cc
cc      use auxiliaryVariables
cc
cc      use nlfunction_setup
cc
cc      implicit none
cc
ccc     Call variables
cc
ccc     Local variables
cc
cc      integer :: ierr,offset=0,igr
cc
ccc     Begin program
cc
cccc      write (*,*) 'DIAG -- setBCsched_app',itime
cc
ccc     Select offset based on whether we are at the beginning of calculation
cc
cc      if (itime == 0) then
cc#if defined(vec_pot)
cc        offset = 8
cc#else
cc        offset = 7
cc#endif
cc      elseif (itime > 0) then
cc
cc        offset = 0
cc
cc      endif
cc
ccc     Set number of BC groups
cc
cc#if defined(vec_pot)
cc      gv%nbc_grp = 9+offset
cc#else
cc      if (fake_ve) then
cc        gv%nbc_grp = 9+offset
cc      else
cc        gv%nbc_grp = 8+offset
cc      endif
cc#endif
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
cc        if (itime == 0) then
cc          gv%bc_grp(1)%nbc_seq = 1
cc          gv%bc_grp(2)%nbc_seq = 1
cc          gv%bc_grp(3)%nbc_seq = 1
cc
cc          offset = 0
cc#if defined(vec_pot)
cc          gv%bc_grp(4)%nbc_seq = 1
cc          offset = 1
cc#endif
cc          gv%bc_grp(4+offset)%nbc_seq = 2
cc          gv%bc_grp(5+offset)%nbc_seq = 3
cc          gv%bc_grp(6+offset)%nbc_seq = 1
cc          gv%bc_grp(7+offset)%nbc_seq = 1
cc          offset = offset + 7
cc        endif
cc
cc#if defined(vec_pot)
cc        gv%bc_grp(1+offset)%nbc_seq = 1
cc        gv%bc_grp(2+offset)%nbc_seq = 0
cc        gv%bc_grp(3+offset)%nbc_seq = 2
cccc        gv%bc_grp(4+offset)%nbc_seq = 5
cc        gv%bc_grp(4+offset)%nbc_seq = 2
cc        gv%bc_grp(5+offset)%nbc_seq = 2
cc        gv%bc_grp(6+offset)%nbc_seq = 2
cccc        gv%bc_grp(7+offset)%nbc_seq = 7
cc        gv%bc_grp(7+offset)%nbc_seq = 4
cc        gv%bc_grp(8+offset)%nbc_seq = 4
cc        gv%bc_grp(9+offset)%nbc_seq = 4
cc#else
cc        gv%bc_grp(1+offset)%nbc_seq = 2
cc        gv%bc_grp(2+offset)%nbc_seq = 0
cc        gv%bc_grp(3+offset)%nbc_seq = 2
cccc        gv%bc_grp(4+offset)%nbc_seq = 5
cc        gv%bc_grp(4+offset)%nbc_seq = 2
cc        gv%bc_grp(5+offset)%nbc_seq = 2
cccc        gv%bc_grp(6+offset)%nbc_seq = 7
cc        gv%bc_grp(6+offset)%nbc_seq = 4
cc        if (fake_ve) then
cc          gv%bc_grp(7+offset)%nbc_seq = 3
cc          gv%bc_grp(8+offset)%nbc_seq = 4
cc          gv%bc_grp(9+offset)%nbc_seq = 1
cc        else
cc          gv%bc_grp(7+offset)%nbc_seq = 3
cc          gv%bc_grp(8+offset)%nbc_seq = 4
cc        endif
cc#endif
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
cc        if (itime == 0) then
cc          !Group 1
cc          gv%bc_grp(1)%bc_seq(1,1) = IRHO_0*AUX          !Aux. Variable
cc          gv%bc_grp(1)%bc_seq(1,2) = SCALAR              !Scalar
cc
cc          !Group 2   
cc          gv%bc_grp(2)%bc_seq(1,1) = ITMP_0*AUX          !Aux. variable
cc          gv%bc_grp(2)%bc_seq(1,2) = SCALAR              !Scalar
cc
cc          !Group 3
cc          gv%bc_grp(3)%bc_seq(1,1) = IETA_0*AUX          !Aux. variable
cc          gv%bc_grp(3)%bc_seq(1,2) = SCALAR              !Scalar
cc
cc          !Group 4
cc          offset = 0
cc#if defined(vec_pot)
cc          gv%bc_grp(4)%bc_seq(1,1) = IACOV_0*AUX         !Aux. variable
cc          gv%bc_grp(4)%bc_seq(1,2) = VECTOR              !Vector
cc          offset = 1
cc#endif
cc          !Group 4+
cc          gv%bc_grp(4+offset)%bc_seq(1,1) = IBCNV_0*AUX  !Aux. variable
cc          gv%bc_grp(4+offset)%bc_seq(1,2) = VECTOR       !Vector
cc
cc          gv%bc_grp(4+offset)%bc_seq(2,1) = IBCOV_0*AUX  !Aux. variable
cc          gv%bc_grp(4+offset)%bc_seq(2,2) = VECTOR       !Vector
cc
cc          !Group 5+
cc          gv%bc_grp(5+offset)%bc_seq(1,1) = IHETA*AUX    !Aux. variable
cc          gv%bc_grp(5+offset)%bc_seq(1,2) = SCALAR       !Scalar
cc
cc          gv%bc_grp(5+offset)%bc_seq(2,1) = IJCNV_0*AUX  !Aux. variable
cc          gv%bc_grp(5+offset)%bc_seq(2,2) = VECTOR       !Vector
cc
cc          gv%bc_grp(5+offset)%bc_seq(3,1) = IJCOV_0*AUX  !Aux. variable
cc          gv%bc_grp(5+offset)%bc_seq(3,2) = VECTOR       !Vector
cc
cc          !Group 6+
cc          gv%bc_grp(6+offset)%bc_seq(1,1) = IVCNV_0*AUX  !Aux. variable
cc          gv%bc_grp(6+offset)%bc_seq(1,2) = VECTOR       !Vector
cc
cc          !Group 7+
cccc          gv%bc_grp(7+offset)%bc_seq(1,1) = IVECNV_0*AUX !Aux. variable
cccc          gv%bc_grp(7+offset)%bc_seq(1,2) = VECTOR       !Vector
cc
cc          gv%bc_grp(7+offset)%bc_seq(1,1) = IENI_0*AUX   !Aux. variable
cc          gv%bc_grp(7+offset)%bc_seq(1,2) = VECTOR       !Vector
cc
cc          offset = offset + 7
cc        endif
cc
cc#if defined(vec_pot)
cc        !Group 1
cc        gv%bc_grp(1+offset)%bc_seq(1,1) = IRHO_AUX*AUX!Aux. Variable
cc        gv%bc_grp(1+offset)%bc_seq(1,2) = SCALAR      !Scalar
cc
cc        gv%bc_grp(1+offset)%bc_seq(2,1) = ITMP_AUX*AUX!Aux. variable
cc        gv%bc_grp(1+offset)%bc_seq(2,2) = SCALAR      !Scalar
cc
cc        !Group 3
cc        gv%bc_grp(3+offset)%bc_seq(1,1) = IETA*AUX    !Aux. variable
cc        gv%bc_grp(3+offset)%bc_seq(1,2) = SCALAR      !Scalar
cc
cc        gv%bc_grp(3+offset)%bc_seq(2,1) = INU*AUX    !Aux. variable
cc        gv%bc_grp(3+offset)%bc_seq(2,2) = SCALAR      !Scalar
cc
cc        !Group 4
cc        gv%bc_grp(4+offset)%bc_seq(1,1) = IACOV*AUX   !Aux. variable
cc        gv%bc_grp(4+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(4+offset)%bc_seq(2,1) = IACNV*AUX   !Aux. variable
cc        gv%bc_grp(4+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cccc        gv%bc_grp(4+offset)%bc_seq(3,1) = IAX         !Dep. variable
cccc        gv%bc_grp(4+offset)%bc_seq(3,2) = SCALAR      !Vector
cccc
cccc        gv%bc_grp(4+offset)%bc_seq(4,1) = IAY         !Dep. variable
cccc        gv%bc_grp(4+offset)%bc_seq(4,2) = SCALAR      !Vector
cccc
cccc        gv%bc_grp(4+offset)%bc_seq(5,1) = IAZ         !Dep. variable
cccc        gv%bc_grp(4+offset)%bc_seq(5,2) = SCALAR      !Vector
cc
cc        !Group 5
cccc        gv%bc_grp(5+offset)%bc_seq(1,1) = IBCNV_N*AUX !Aux. variable
cccc        gv%bc_grp(5+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(5+offset)%bc_seq(1,1) = IBCNV*AUX   !Aux. variable
cc        gv%bc_grp(5+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(5+offset)%bc_seq(2,1) = IBCOV*AUX   !Aux. variable
cc        gv%bc_grp(5+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc        !Group 6
cc        gv%bc_grp(6+offset)%bc_seq(1,1) = IJCNV*AUX   !Aux. variable
cc        gv%bc_grp(6+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cccc        gv%bc_grp(6+offset)%bc_seq(2,1) = IJCNV_N*AUX !Aux. variable
cccc        gv%bc_grp(6+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(6+offset)%bc_seq(2,1) = IJCOV*AUX   !Aux. variable
cc        gv%bc_grp(6+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc        !Group 7
cc        gv%bc_grp(7+offset)%bc_seq(1,1) = IVCNV*AUX   !Aux. variable
cc        gv%bc_grp(7+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(7+offset)%bc_seq(2,1) = IVCNV_N*AUX !Aux. variable
cc        gv%bc_grp(7+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(7+offset)%bc_seq(3,1) = IVCOV*AUX   !Aux. variable
cc        gv%bc_grp(7+offset)%bc_seq(3,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(7+offset)%bc_seq(4,1) = IVCNV_N*AUX !Aux. variable
cc        gv%bc_grp(7+offset)%bc_seq(4,2) = VECTOR      !Vector
cc
cccc        gv%bc_grp(7+offset)%bc_seq(5,1) = IVX         !Dep. variable
cccc        gv%bc_grp(7+offset)%bc_seq(5,2) = SCALAR      !Scalar
cccc
cccc        gv%bc_grp(7+offset)%bc_seq(6,1) = IVY         !Dep. variable
cccc        gv%bc_grp(7+offset)%bc_seq(6,2) = SCALAR      !Scalar
cccc
cccc        gv%bc_grp(7+offset)%bc_seq(7,1) = IVZ         !Dep. variable
cccc        gv%bc_grp(7+offset)%bc_seq(7,2) = SCALAR      !Scalar
cc
cc        !Group 8
cc        gv%bc_grp(8+offset)%bc_seq(1,1) = IVECNV*AUX  !Aux. variable
cc        gv%bc_grp(8+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(8+offset)%bc_seq(2,1) = IVECOV*AUX  !Aux. variable
cc        gv%bc_grp(8+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(8+offset)%bc_seq(3,1) = IVECOV_N*AUX!Aux. variable
cc        gv%bc_grp(8+offset)%bc_seq(3,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(8+offset)%bc_seq(4,1) = IVEFCNV*AUX !Aux. variable
cc        gv%bc_grp(8+offset)%bc_seq(4,2) = VECTOR      !Vector
cc
cc        !Group 9
cc        gv%bc_grp(9+offset)%bc_seq(1,1) = IDIVPI*AUX  !Aux. variable
cc        gv%bc_grp(9+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(9+offset)%bc_seq(2,1) = IDIVPE*AUX  !Aux. variable
cc        gv%bc_grp(9+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(9+offset)%bc_seq(3,1) = IEH*AUX     !Aux. variable
cc        gv%bc_grp(9+offset)%bc_seq(3,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(9+offset)%bc_seq(4,1) = IENI*AUX    !Aux. variable
cc        gv%bc_grp(9+offset)%bc_seq(4,2) = VECTOR      !Vector
cc 
cc#else
cc        !Group 1
cc        gv%bc_grp(1+offset)%bc_seq(1,1) = IRHO_AUX*AUX!Aux. Variable
cc        gv%bc_grp(1+offset)%bc_seq(1,2) = SCALAR      !Scalar
cc
cc        gv%bc_grp(1+offset)%bc_seq(2,1) = ITMP_AUX*AUX!Aux. variable
cc        gv%bc_grp(1+offset)%bc_seq(2,2) = SCALAR      !Scalar
cc
cc        !Group 3
cc        gv%bc_grp(3+offset)%bc_seq(1,1) = IETA*AUX    !Aux. variable
cc        gv%bc_grp(3+offset)%bc_seq(1,2) = SCALAR      !Scalar
cc
cc        gv%bc_grp(3+offset)%bc_seq(2,1) = INU*AUX     !Aux. variable
cc        gv%bc_grp(3+offset)%bc_seq(2,2) = SCALAR      !Scalar
cc
cc        !Group 4
cccc        gv%bc_grp(4+offset)%bc_seq(1,1) = IBCNV_N*AUX !Aux. variable
cccc        gv%bc_grp(4+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(4+offset)%bc_seq(1,1) = IBCNV*AUX   !Aux. variable
cc        gv%bc_grp(4+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(4+offset)%bc_seq(2,1) = IBCOV*AUX   !Aux. variable
cc        gv%bc_grp(4+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cccc        gv%bc_grp(4+offset)%bc_seq(3,1) = IBX         !Dep. variable
cccc        gv%bc_grp(4+offset)%bc_seq(3,2) = SCALAR      !Scalar
cccc
cccc        gv%bc_grp(4+offset)%bc_seq(4,1) = IBY         !Dep. variable
cccc        gv%bc_grp(4+offset)%bc_seq(4,2) = SCALAR      !Scalar
cccc
cccc        gv%bc_grp(4+offset)%bc_seq(5,1) = IBZ         !Dep. variable
cccc        gv%bc_grp(4+offset)%bc_seq(5,2) = SCALAR      !Scalar
cc
cc        !Group 5
cc        gv%bc_grp(5+offset)%bc_seq(1,1) = IJCNV*AUX   !Aux. variable
cc        gv%bc_grp(5+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cccc        gv%bc_grp(5+offset)%bc_seq(2,1) = IJCNV_N*AUX !Aux. variable
cccc        gv%bc_grp(5+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(5+offset)%bc_seq(2,1) = IJCOV*AUX   !Aux. variable
cc        gv%bc_grp(5+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc        !Group 6
cc        gv%bc_grp(6+offset)%bc_seq(1,1) = IVCNV*AUX   !Aux. variable
cc        gv%bc_grp(6+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(6+offset)%bc_seq(2,1) = IVCNV_N*AUX !Aux. variable
cc        gv%bc_grp(6+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(6+offset)%bc_seq(3,1) = IVCOV*AUX   !Aux. variable
cc        gv%bc_grp(6+offset)%bc_seq(3,2) = VECTOR      !Vector
cc
cc        gv%bc_grp(6+offset)%bc_seq(4,1) = IVCNV_N*AUX !Aux. variable
cc        gv%bc_grp(6+offset)%bc_seq(4,2) = VECTOR      !Vector
cc
cccc        gv%bc_grp(6+offset)%bc_seq(5,1) = IVX         !Dep. variable
cccc        gv%bc_grp(6+offset)%bc_seq(5,2) = SCALAR      !Scalar
cccc
cccc        gv%bc_grp(6+offset)%bc_seq(6,1) = IVY         !Dep. variable
cccc        gv%bc_grp(6+offset)%bc_seq(6,2) = SCALAR      !Scalar
cccc
cccc        gv%bc_grp(6+offset)%bc_seq(7,1) = IVZ         !Dep. variable
cccc        gv%bc_grp(6+offset)%bc_seq(7,2) = SCALAR      !Scalar
cc
cc        if (fake_ve) then
cc          !Group 7
cc          gv%bc_grp(7+offset)%bc_seq(1,1) = IDIVPI*AUX  !Aux. variable
cc          gv%bc_grp(7+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc          gv%bc_grp(7+offset)%bc_seq(2,1) = IEH*AUX    !Aux. variable
cc          gv%bc_grp(7+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc          gv%bc_grp(7+offset)%bc_seq(3,1) = IENI*AUX    !Aux. variable
cc          gv%bc_grp(7+offset)%bc_seq(3,2) = VECTOR      !Vector
cc
cc          !Group 8
cc          gv%bc_grp(8+offset)%bc_seq(1,1) = IVECNV*AUX  !Aux. variable
cc          gv%bc_grp(8+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc          gv%bc_grp(8+offset)%bc_seq(2,1) = IVECOV*AUX  !Aux. variable
cc          gv%bc_grp(8+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc          gv%bc_grp(8+offset)%bc_seq(3,1) = IVECOV_N*AUX!Aux. variable
cc          gv%bc_grp(8+offset)%bc_seq(3,2) = VECTOR      !Vector
cc
cc          gv%bc_grp(8+offset)%bc_seq(4,1) = IVEFCNV*AUX !Aux. variable
cc          gv%bc_grp(8+offset)%bc_seq(4,2) = VECTOR      !Vector
cc
cc          !Group 9
cc          gv%bc_grp(9+offset)%bc_seq(1,1) = IDIVPE*AUX  !Aux. variable
cc          gv%bc_grp(9+offset)%bc_seq(1,2) = VECTOR      !Vector
cc        else
cc          !Group 7
cc          gv%bc_grp(7+offset)%bc_seq(1,1) = IVECNV*AUX  !Aux. variable
cc          gv%bc_grp(7+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc          gv%bc_grp(7+offset)%bc_seq(2,1) = IVECOV*AUX  !Aux. variable
cc          gv%bc_grp(7+offset)%bc_seq(2,2) = VECTOR      !Vector
cc                     
cc          gv%bc_grp(7+offset)%bc_seq(3,1) = IVECOV_N*AUX!Aux. variable
cc          gv%bc_grp(7+offset)%bc_seq(3,2) = VECTOR      !Vector
cc
cc          !Group 8
cc          gv%bc_grp(8+offset)%bc_seq(1,1) = IDIVPI*AUX  !Aux. variable
cc          gv%bc_grp(8+offset)%bc_seq(1,2) = VECTOR      !Vector
cc
cc          gv%bc_grp(8+offset)%bc_seq(2,1) = IDIVPE*AUX  !Aux. variable
cc          gv%bc_grp(8+offset)%bc_seq(2,2) = VECTOR      !Vector
cc
cc          gv%bc_grp(8+offset)%bc_seq(3,1) = IEH*AUX    !Aux. variable
cc          gv%bc_grp(8+offset)%bc_seq(3,2) = VECTOR      !Vector
cc
cc          gv%bc_grp(8+offset)%bc_seq(4,1) = IENI*AUX    !Aux. variable
cc          gv%bc_grp(8+offset)%bc_seq(4,2) = VECTOR      !Vector
cc        endif
cc#endif
cc
cc      endif
cc
ccc     End program
cc
cc      end subroutine setBCsched_app

c module local_BCS_variables
c####################################################################
      module local_BCS_variables

        use variables

        use imposeBCinterface, ONLY: extrapolateBC

        use auxiliaryVariables

        use operators

        use timestepping

        use equilibrium

        use nlfunction_setup, ONLY: res,chi_perp,vis,E_fld_res
     .                             ,ion_hall,fake_ve,E_fld_hall
     .                             ,EOM_divPi,hres
     .                             ,nc_eom_v,nc_eom_jxb
     .                             ,save_tn,pinch_flow
     .                             ,res_ijk,vis_ijk,form_ve,form_jfake

#if defined(vec_pot)
        integer :: order_bc=1
#else
        integer :: order_bc=1
#endif
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

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX ) = EQU

        if (gamma > 1d0)
     .       where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = EQU
cc     .       where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = DIR

      case ('ppnsl','p3nsl','ohtor')           !No slip BCs, pinch flow

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX ) = EQU

        !Dirichlet boundary conditions for tangential components
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY ) = DIR
        where (bbcs(:,IVZ) == NEU) bbcs(:,IVZ ) = DIR

        if (gamma > 1d0)
     .       where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = EQU
cc     .       where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = DIR

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

      integer :: i,j,k,ig,jg,kg,ord

      real(8) :: b2,bbcnv(3),bbcov(3),dum(3),grad_p(3),etah,gsup(3,3)
     .          ,gsub(3,3),jach

c$$$      real(8) :: x(5),vec(4),xt(1)

c     Begin program

      if (.not.pinch_flow) return

      if (igrid /= 1) then
        write (*,*) 'ExB_flow_bc not ready for MG treatment'
        stop
      endif

      i = nx

c$$$      ord = 3
c$$$
c$$$      call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
c$$$
c$$$      x = (/grid_params%xx(ig-3)
c$$$     $     ,grid_params%xx(ig-2)
c$$$     $     ,grid_params%xx(ig-1)
c$$$     $     ,grid_params%xx(ig  )
c$$$     $     ,grid_params%xx(ig+1)/)
c$$$
c$$$      xt = 0.5*(grid_params%xx(ig+1)
c$$$     $         +grid_params%xx(ig  ))

      do k=1,nz
        do j=1,ny
          jach = 0.5*(gmetric%grid(igrid)%jac(i+1,j,k)
     .               +gmetric%grid(igrid)%jac(i  ,j,k))

          bbcnv=0.5*(bcnv(i+1,j,k,:)+bcnv(i  ,j,k,:))
          bbcov=0.5*(matmul(gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .                     ,bcnv(i  ,j,k,:))
     .              +matmul(gmetric%grid(igrid)%gsub(i+1,j,k,:,:)
     .                     ,bcnv(i+1,j,k,:)))

          b2 = dot_product(bbcnv,bbcov)/jach

          if (b2 > 0d0) then
            if (PRESENT(dbcov)) then  !Linearize ExB flow BC
              dum = 0.5*(dbcov(i,j,k,:)+dbcov(i+1,j,k,:))
              bbcov = dum - 2*dot_product(bbcnv,dum)/(jach*b2)*bbcov
            endif

            grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,hex=1)
c$$$            grad_p = 0.5*(grad(i  ,j,k,nx,ny,nz,igrid,igrid,igrid,prs)
c$$$     $                   +grad(i+1,j,k,nx,ny,nz,igrid,igrid,igrid,prs))
c$$$
c$$$            vec = (/ prs(i-2,j,k)
c$$$     $              ,prs(i-1,j,k)
c$$$     $              ,prs(i  ,j,k)
c$$$     $              ,prs(i+1,j,k) /)
c$$$            call IntDriver1d (4,x(2:5),vec,1,xt,grad_p(1:1),ord,deriv=1)

            gsup = 0.5*(gmetric%grid(igrid)%gsup(i+1,j,k,:,:)
     .                 +gmetric%grid(igrid)%gsup(i  ,j,k,:,:))

            grad_p = matmul(gsup,grad_p)

cc            etah = 2d0/(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))
            etah = 0.5*(eta(i+1,j,k)+eta(i,j,k))

            !ExB velocity
            v0(i+1,j,k,1)=(E0(2)*bbcov(3)-E0(3)*bbcov(2)
     .                    -etah*grad_p(1))/b2
            v0(i+1,j,k,2)=(E0(3)*bbcov(1)-E0(1)*bbcov(3)
     .                    -etah*grad_p(2))/b2
            v0(i+1,j,k,3)=(E0(1)*bbcov(2)-E0(2)*bbcov(1)
     .                    -etah*grad_p(3))/b2
          endif
        enddo
      enddo

c     End program

      end subroutine ExB_flow_bc

c     J_pinch_bc
c     #################################################################
      subroutine J_pinch_bc(igrid,nx,ny,nz,bcnv,prs,eta,j0)

c     -----------------------------------------------------------------
c     Finds current boundary condition when Eo is imposed at boundary
c     (returns covariant current).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,nx,ny,nz
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1)   :: prs,eta
      real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3) :: bcnv,j0

c     Local variables

      integer :: i,j,k,ig,jg,kg,ord
      real(8) :: b2,bbcnv(3),bbcov(3),grad_p(3),alpha
     .          ,gsub(3,3),ietah,jach

c$$$      real(8) :: x(5),vec(4),xt(1)

c     Begin program

      if (.not.pinch_flow) return

      if (igrid /= 1) then
        write (*,*) 'J_pinch_bc not ready for MG treatment'
        stop
      endif

      i = nx

c$$$      ord = 3

      call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

c$$$      x = (/grid_params%xx(ig-3)
c$$$     $     ,grid_params%xx(ig-2)
c$$$     $     ,grid_params%xx(ig-1)
c$$$     $     ,grid_params%xx(ig  )
c$$$     $     ,grid_params%xx(ig+1)/)
c$$$
c$$$      xt = 0.5*(grid_params%xx(ig+1)
c$$$     $         +grid_params%xx(ig  ))

      do k=1,nz
        do j=1,ny
          !First-order extrapolation of B to boundary (critical)
          bbcnv = bcnv(i-1,j,k,:) + (bcnv(i,j,k,:)-bcnv(i-1,j,k,:))
     $                               /(grid_params%xx(ig)
     $                                -grid_params%xx(ig-1))
     $                               *(0.5*grid_params%xx(ig+1)
     $                                +0.5*grid_params%xx(ig)
     $                                -    grid_params%xx(ig-1))

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

          gsub = 0.5*(gmetric%grid(igrid)%gsub(i  ,j,k,:,:)
     .               +gmetric%grid(igrid)%gsub(i+1,j,k,:,:))
          bbcov = matmul(gsub,bbcnv)

          jach = 0.5*(gmetric%grid(igrid)%jac(i  ,j,k)
     .               +gmetric%grid(igrid)%jac(i+1,j,k))

          b2 = dot_product(bbcov,bbcnv)/jach

          if (b2 > 0d0) then
            ietah = 2d0/(eta(i+1,j,k)+eta(i,j,k))
cc            ietah =0.5*(1d0/eta(i+1,j,k)+1d0/eta(i,j,k))

            grad_p = grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,prs,hex=1)
c$$$            grad_p = 0.5*(grad(i  ,j,k,nx,ny,nz,igrid,igrid,igrid,prs)
c$$$     $                   +grad(i+1,j,k,nx,ny,nz,igrid,igrid,igrid,prs))
c$$$
c$$$            vec = (/ prs(i-2,j,k)
c$$$     $              ,prs(i-1,j,k)
c$$$     $              ,prs(i  ,j,k)
c$$$     $              ,prs(i+1,j,k) /)
c$$$            call IntDriver1d (4,x(2:5),vec,1,xt,grad_p(1:1),ord,deriv=1)

            alpha = dot_product(E0,bbcnv)/jach*ietah/b2

            j0(i+1,j,k,:)=alpha*bbcnv + cross_product(bbcov,grad_p)/b2

          endif
        enddo
      enddo

c     End program

      end subroutine J_pinch_bc

      end module local_BCS_variables

