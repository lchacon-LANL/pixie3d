c setupNonlinearFunction
c#################################################################
      subroutine setupNonlinearFunction(igx,igy,igz,varray)
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use parameters

      use var_setup

      use nlfunction_setup

      use equilibrium

      use timeStepping

      use error

      use imposeBCinterface

#if !defined(samrai)
      use nk
#endif
      use local_BCS_variables

      use xdraw_io

      implicit none

c Call variables

      integer :: igx,igy,igz
      type(var_array),pointer :: varray

c Local variables

      integer :: nx,ny,nz
      real(8) :: minrho,mintmp,lrmin,ltmin

c Interfaces

#if defined(flux_rhs)
      INTERFACE
        function set_nlf_flux_x(igr,nx,ny,nz) result(flx_x)
        use var_setup
        integer :: igr,nx,ny,nz
        real(8) :: flx_x(nx+1,ny,nz,neqd)
        end function set_nlf_flux_x

        function set_nlf_flux_y(igr,nx,ny,nz) result(flx_y)
        use var_setup
        integer :: igr,nx,ny,nz
        real(8) :: flx_y(nx,ny+1,nz,neqd)
        end function set_nlf_flux_y

        function set_nlf_flux_z(igr,nx,ny,nz) result(flx_z)
        use var_setup
        integer :: igr,nx,ny,nz
        real(8) :: flx_z(nx,ny,nz+1,neqd)
        end function set_nlf_flux_z

        function set_nlf_src(igr,nx,ny,nz) result(src)
        use var_setup
        integer :: igr,nx,ny,nz
        real(8) :: src(nx,ny,nz,neqd)
        end function set_nlf_src

      END INTERFACE
#endif

c Begin program

      nx = gv%gparams%nxv(igx)
      ny = gv%gparams%nyv(igy)
      nz = gv%gparams%nzv(igz)

c Impose boundary conditions and find auxiliary quantities

#if !defined (samrai)
      call applyBC(igx,varray,gv%aux,it=itime,bc_order=order_bc)
#else
      call defineAuxVariables(gv%aux)  !Set pointers
#endif

c Fill LOCAL auxiliary variables

      call fillLocalAuxVars(igx,varray,gv%aux)

c Safeguards

      if (old_time_eval) then
        call chk_neg(rho,'density')
        call chk_neg(prs,'pressure')
        call chk_neg(heat_src,'heat source')
      endif
      
c Setup fluxes

#if defined(flux_rhs)
      gv%flux_x = set_nlf_flux_x(igx,nx,ny,nz)
      gv%flux_y = set_nlf_flux_y(igx,nx,ny,nz)
      gv%flux_z = set_nlf_flux_z(igx,nx,ny,nz)
      gv%src    = set_nlf_src   (igx,nx,ny,nz)
#endif

c End program

      contains

c     chk_neg
c     ###############################################################
      subroutine chk_neg(array,desc)

      implicit none
c     ---------------------------------------------------------------
c     Checks for negative values in array, and if so terminates
c     with message.
c     ---------------------------------------------------------------

c     Call variables

      real(8) :: array(0:nx+1,0:ny+1,0:nz+1)
      character(*) :: desc

c     Local variables

      integer :: index(3)
      real(8) :: minv,ltmin,minvg

c     Begin program

      minv = minval(array(1:nx,1:ny,1:nz))

      minvg = pmin(minv)  !Parallel min
      
      if (minvg < -1d-14) then

        if (minv < -1d-14) then
          index = globalIndex(gv%gparams,igx
     .                       ,minloc(array(1:nx,1:ny,1:nz)))
          write(*,*) 'proc=',my_rank,' ',trim(desc),'=',minv
     .              ,' global pos=',index
        endif

        call pstop('setupNonlinearFunction',
     .   'Warning: negative values for '//trim(desc)//' are occurring')

      endif

      end subroutine chk_neg

      end subroutine setupNonlinearFunction

c killNonlinearFunction
c#################################################################
      subroutine killNonlinearFunction
c------------------------------------------------------------------
c     Frees memory
c------------------------------------------------------------------

      use nlfunction_setup

      implicit none

c Begin program

      call destroyAppBCs

c End program

      end subroutine killNonlinearFunction
