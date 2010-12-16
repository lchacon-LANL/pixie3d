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

c Impose boundary conditions and find auxiliary quantities

#if !defined (samrai)
      save_tn = old_time_eval .or. source_eval
cc      save_Bn = (save_tn .and. lag_pinch_bc)  !Initialize Bn when processing source

      call applyBC(varray,igx,igy,igz)
#else
      call allocAuxVariables(gv%aux)  !Set pointers
#endif

c Fill LOCAL auxiliary variables

      call fillLocalAuxVars(igx,igy,igz,varray,gv%aux)

c Safeguards

      !Mass density positivity
      minrho = minval(rho(1:nx,1:ny,1:nz))

#if defined(petsc)
      lrmin = minrho
      call MPI_Allreduce(lrmin,minrho,1
     .                  ,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
#endif

      if (minrho < 0d0) then
        close (urecord)  !Close record file before exiting
        write(*,*)'proc=',my_rank,'rho=',minrho,'loc='
     $            ,minloc(rho(1:nx,1:ny,1:nz))
        call pstop('setupNonlinearFunction',
     .             'Warning: negative densities are occurring')
      endif

      !Temperature positivity (enforce)
      mintmp = minval(tmp(1:nx,1:ny,1:nz))

#if defined(petsc)
      ltmin = mintmp
      call MPI_Allreduce(ltmin,mintmp,1
     .                  ,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
#endif

      if (mintmp < 0d0) then
        close (urecord) !Close record file before exiting
        write(*,*)'proc=',my_rank,'tmp=',mintmp,'loc='
     $            ,minloc(tmp(1:nx,1:ny,1:nz))
        call pstop('setupNonlinearFunction',
     .             'Warning: negative temperatures are occurring')
      endif

c End program

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
