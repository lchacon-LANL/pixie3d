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

      use constants

      use timeStepping

      use error

      use imposeBCinterface

      use newtongm

      use local_BCS_variables

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

      INTERFACE
        subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine imposeBoundaryConditions
      END INTERFACE

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

cc      save_Bn = (jit == -2 .and. rfx_pinch_bc)  !Initialize Bn when processing source
      save_Bn = (jit < 0 .and. rfx_pinch_bc)  !Initialize Bn when processing source

      call imposeBoundaryConditions(varray,igx,igy,igz)

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
cc        if (my_rank == 0) write (*,*) 'rho=',minrho
        write(*,*)'proc=',my_rank,'rho=',lrmin,'loc='
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

c Viscosity (resistivity found in BC routine)

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            nuu(i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

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

c Divergence of ion stress tensor

      do k=1,nz
        do j=1,ny
          do i=1,nx
            div_pi(i,j,k,:)=EOM_divPi(i,j,k,nx,ny,nz,igx,igy,igz,vcnv
     .                               ,diff=nuu)
          enddo
        enddo
      enddo

c Electron velocity

      if (di > 0d0) then
        vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv(:,:,:,1)/rho
        vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv(:,:,:,2)/rho
        vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv(:,:,:,3)/rho

        vecov(:,:,:,1) = vcov(:,:,:,1)-di*jcov(:,:,:,1)/rho
        vecov(:,:,:,2) = vcov(:,:,:,2)-di*jcov(:,:,:,2)/rho
        vecov(:,:,:,3) = vcov(:,:,:,3)-di*jcov(:,:,:,3)/rho
      else
        vecnv = vcnv
        vecov = vcov
      endif

      if (jit == -2) vecnv_0 = vecnv   !Store equilibrium info

c Divergence of electron stress tensor

      if (di*heta > 0d0) then
        vdummy = vecnv - vecnv_0
ctest        vdummy = vcnv

        do k=1,nz
          do j=1,ny
            do i=1,nx
              div_pe(i,j,k,:)=EOM_divPe(i,j,k,nx,ny,nz,igx,igy,igz
     .                                 ,vdummy)
            enddo
          enddo
        enddo
      else
        div_pe = 0d0
      endif

c Electric field (covariant vector)

      if (hall_new .and. jit < 0) vcov_n = vcov !Store time n info

      if (de > 0d0 .and. jit < 0) vecov_n = vecov !Store time n info

      call find_E(nx,ny,nz,igx,igy,igz,E_ni)

c Other auxiliary variables

      pcnv(:,:,:,1) = px
      pcnv(:,:,:,2) = py
      pcnv(:,:,:,3) = pz

c End program

      end subroutine setupNonlinearFunction

c killNonlinearFunction
c#################################################################
      subroutine killNonlinearFunction
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use nlfunction_setup

      implicit none

c Call variables

c Local variables

c Begin program

cc#if !defined(vec_pot)
cc      nullify(rho,bx,by,bz,px,py,pz,tmp)
cc#else
cc      nullify(rho,ax,ay,az,px,py,pz,tmp)
cc#endif

#if defined(samrai)
      call deallocAuxVariables
#endif

c End program

      end subroutine killNonlinearFunction
