c createNonlinearFunction
c######################################################################
      subroutine createNonlinearFunction

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use parameters

      use grid_aliases

      use transport_params

      implicit none

c Call variables

c Local variables

c Begin program

      call allocateApplicationVariables

      call defineGridAliases

c End programs

      end subroutine createNonlinearFunction

c allocateApplicationVariables
c######################################################################
      subroutine allocateApplicationVariables

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use parameters

      use nlfunction_setup

      use precond_variables

      implicit none

c Call variables

c Local variables

c Begin program

        allocate (pcnv  (ilom:ihip,jlom:jhip,klom:khip,3)
     $           ,vecnv (ilom:ihip,jlom:jhip,klom:khip,3)
     $           ,vdummy(ilom:ihip,jlom:jhip,klom:khip,3))

        allocate (eeta  (ilom:ihip,jlom:jhip,klom:khip)
     .           ,nuu   (ilom:ihip,jlom:jhip,klom:khip))

        allocate (vcnv  (ilom:ihip,jlom:jhip,klom:khip,3)
     .           ,vcov  (ilom:ihip,jlom:jhip,klom:khip,3)
     .           ,vlap  (ilom:ihip,jlom:jhip,klom:khip,3)
     .           ,E_ni  (ilom:ihip,jlom:jhip,klom:khip,3)
     .           ,jcov_0(ilom:ihip,jlom:jhip,klom:khip,3))

#if defined(vec_pot)
        allocate (acnv  (ilom:ihip,jlom:jhip,klom:khip,3)
     $           ,acov  (ilom:ihip,jlom:jhip,klom:khip,3)
     $           ,vlap_a(ilom:ihip,jlom:jhip,klom:khip,3))
#else
        allocate (bcnv  (ilom:ihip,jlom:jhip,klom:khip,3)
     .           ,bcov  (ilom:ihip,jlom:jhip,klom:khip,3))
#endif

c SI operator
        allocate(b_n(ilom:ihip,jlom:jhip,klom:khip,3)
     .          ,p_n(ilom:ihip,jlom:jhip,klom:khip,1)
     .          ,vcov_n(ilom:ihip,jlom:jhip,klom:khip,3)
     .          ,dv_dt (ilom:ihip,jlom:jhip,klom:khip,3))
c SI operator

c Auxiliary variables

        call allocAuxVariables

c Preconditioner

        call allocPrecVariables

c End program

      end subroutine allocateApplicationVariables

c defineGridAliases
c####################################################################
      subroutine defineGridAliases
c--------------------------------------------------------------------
c     Sets aliases for grid quantities
c--------------------------------------------------------------------

      use grid

      use grid_aliases

      implicit none

c Call variables

c Local variables

c Begin program

      if (.not.associated(xx)) then
        xx => grid_params%xx
        yy => grid_params%yy
        zz => grid_params%zz

        dxh => grid_params%dxh
        dyh => grid_params%dyh
        dzh => grid_params%dzh

        dx => grid_params%dx
        dy => grid_params%dy
        dz => grid_params%dz
      endif

c End

      end subroutine defineGridAliases

c destroyNonlinearFunction
c######################################################################
      subroutine destroyNonlinearFunction

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use parameters

      use nlfunction_setup

      use precond_variables

      implicit none

c Call variables

c Local variables

c Begin program

        deallocate (pcnv,vecnv,vdummy)

        deallocate (eeta,nuu)

        deallocate (vcnv,vcov,vlap,E_ni,jcov_0)

#if defined(vec_pot)
        deallocate (acnv,acov,vlap_a)
#else
        deallocate (bcnv,bcov)
#endif

c SI operator
        deallocate(b_n,p_n,vcov_n,dv_dt)
c SI operator

        call deallocAuxVariables

        call deallocPrecVariables

c End programs

      end subroutine destroyNonlinearFunction
