c setupNonlinearFunction
c#################################################################
      subroutine setupNonlinearFunction(igx,igy,igz,varray)
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use parameters

      use variable_setup

      use nlfunction_setup

      use equilibrium

      use timeStepping

      use error

      use imposeBCinterface

      use newtongm

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

#if !defined (samrai)
      INTERFACE
        subroutine applyBC (varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine applyBC
      END INTERFACE
#endif

c Begin program

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Set aliases

      rho => varray%array_var(IRHO)%array
      px  => varray%array_var(IVX )%array
      py  => varray%array_var(IVY )%array
      pz  => varray%array_var(IVZ )%array
#if !defined(vec_pot)
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
#else
      ax  => varray%array_var(IAX )%array
      ay  => varray%array_var(IAY )%array
      az  => varray%array_var(IAZ )%array
#endif
      tmp => varray%array_var(ITMP)%array

c Enslave density when appropriate

      if (.not.solve_rho) then
        if (adiabatic) then
          where (u_0%array_var(ITMP)%array /= 0d0)
            rho =     u_0%array_var(IRHO)%array
     .          *(tmp/u_0%array_var(ITMP)%array)**(1./gamma)
          end where
cc        else
cc          rho = 1d0
        endif
      endif

c Impose boundary conditions and find auxiliary quantities

#if !defined (samrai)
      save_tn = old_time_eval .or. source_eval
      call applyBC(varray,igx,igy,igz)
#else
      call allocAuxVariables(gv%aux)  !Set pointers
#endif

c Fill LOCAL auxiliary variables

      call fillLocalAuxVars(igx,igy,igz,varray,gv%aux)

c Safeguards

      call chk_neg(rho,'density')
      call chk_neg(tmp,'temperature')

cc      !Mass density positivity
cc      minrho = minval(rho(1:nx,1:ny,1:nz))
cc
cc#if defined(petsc)
cc      lrmin = minrho
cc      call MPI_Allreduce(lrmin,minrho,1
cc     .                  ,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
cc#endif
cc
cc      if (minrho < 0d0) then
cc        close (urecord)  !Close record file before exiting
cc        write(*,*)'proc=',my_rank,'rho=',minrho,'loc='
cc     $            ,minloc(rho(1:nx,1:ny,1:nz))
cc        call pstop('setupNonlinearFunction',
cc     .             'Warning: negative densities are occurring')
cc      endif
cc
cc      !Temperature positivity (enforce)
cc      mintmp = minval(tmp(1:nx,1:ny,1:nz))
cc
cc#if defined(petsc)
cc      ltmin = mintmp
cc      call MPI_Allreduce(ltmin,mintmp,1
cc     .                  ,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
cc#endif
cc
cc      if (mintmp < 0d0) then
cc        close (urecord) !Close record file before exiting
cc        write(*,*)'proc=',my_rank,'tmp=',mintmp,'loc='
cc     $            ,minloc(tmp(1:nx,1:ny,1:nz))
cc        call pstop('setupNonlinearFunction',
cc     .             'Warning: negative temperatures are occurring')
cc      endif

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

#if defined(petsc)
      call MPI_Allreduce(minv,minvg,1
     .                  ,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
#endif

      if (minvg < 0d0) then

        close (urecord) !Close record file before exiting

        if (minv < 0d0) then
          index = globalIndex(igx,minloc(array(1:nx,1:ny,1:nz)))
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

      call deallocLocalAuxVars

c End program

      end subroutine killNonlinearFunction
