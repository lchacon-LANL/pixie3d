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

      integer    :: i,j,k,ig,jg,kg,nx,ny,nz

      integer    :: dim,loc,ibc,ieq,bctype,bcnd(6,3)

      real(8)    :: dh,vmax,cnv(3),cov(3),minrho,mintmp,lrmin,ltmin
     .             ,x1,y1,z1,bnorm

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
      call allocAuxVariables(gv%aux)
#endif

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

c Store B/|B| unit vector (for parallel heat conduction and/or boostrap current)

      if (chi_par /= chi .or. bootstrap) then
        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz
     .                                      ,ig,jg,kg,x1,y1,z1)

              bnorm = vectorNorm(i,j,k,igx,igy,igz
     .                     ,bcnv(i,j,k,1),bcnv(i,j,k,2),bcnv(i,j,k,3)
     .                     ,.false.)

              bhat(i,j,k,:) = bcnv(i,j,k,:)/sqrt(bnorm)
            enddo
          enddo
        enddo
      endif

c Other auxiliary variables

      !Momentum
      pcnv(:,:,:,1) = px
      pcnv(:,:,:,2) = py
      pcnv(:,:,:,3) = pz

      !Temperature transport velocity, v_star = v - di*j/rho/a_p
      if (slava) then
        vscnv = vcnv            !Also to recover old TM Hall performance results
      else
        vscnv(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv(:,:,:,1)/rho/a_p
        vscnv(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv(:,:,:,2)/rho/a_p
        vscnv(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv(:,:,:,3)/rho/a_p
      endif

      !Psi
      psi = 0d0

cc      !Previous time step B
cc      if (save_tn) bcnv_n = bcnv !Save previous time step
cc
cc      !Previous time step J
cc      if (save_tn .and. ion_hall) jcnv_n = jcnv
cc
cc      !Equilibrium electron velocity
cc      vecnv_0(:,:,:,1) = vcnv_0(:,:,:,1)-di*jcnv_0(:,:,:,1)/rho_0
cc      vecnv_0(:,:,:,2) = vcnv_0(:,:,:,2)-di*jcnv_0(:,:,:,2)/rho_0
cc      vecnv_0(:,:,:,3) = vcnv_0(:,:,:,3)-di*jcnv_0(:,:,:,3)/rho_0

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

c Call variables

c Local variables

c Begin program

cc      call deallocAuxVariables(gv%aux)

c End program

      end subroutine killNonlinearFunction
