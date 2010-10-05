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

      real(8),allocatable,dimension(:,:,:,:) :: vdum2

c Interfaces

#if !defined (samrai)
      INTERFACE
        subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine imposeBoundaryConditions
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
cc      save_tn = (jit < 0)
      save_tn = old_time_eval .or. source_eval
      save_Bn = (save_tn .and. lag_pinch_bc)  !Initialize Bn when processing source

      call imposeBoundaryConditions(varray,igx,igy,igz)
#else
      call allocAuxVariables
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

c Electron velocity

cc#if defined(vec_pot)
cc      if (di > 0d0) then
cc        vecnv(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv(:,:,:,1)/rho
cc        vecnv(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv(:,:,:,2)/rho
cc        vecnv(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv(:,:,:,3)/rho
cc
cc        vecov(:,:,:,1) = vcov(:,:,:,1)-di*jcov(:,:,:,1)/rho
cc        vecov(:,:,:,2) = vcov(:,:,:,2)-di*jcov(:,:,:,2)/rho
cc        vecov(:,:,:,3) = vcov(:,:,:,3)-di*jcov(:,:,:,3)/rho
cc      else
cc        vecnv = vcnv
cc        vecov = vcov
cc      endif
cc#endif

cc      if (ion_hall .and. jit < 0) vcov_n = vcov !Store time n info

cc      if (de > 0d0 .and. jit < 0) vecov_n = vecov !Store time n info

c Divergence of electron stress tensor

cc      if (di*heta > 0d0) then
cc        if (ion_hall) then
cc          if (jit == -2) vcnv_0 = vcnv !Store equilibrium info
cc
cc          !curlcurlv = v-di*curl(curl(-di*v))/rho
cc          allocate(vdum2(0:nx+1,0:ny+1,0:nz+1,3))
cc          vdum2 = -di*(vcnv-vcnv_0)
cc
cc          do k=1,nz
cc            do j=1,ny
cc              do i=1,nx
cc                curlcurlv(i,j,k,:) = curlcurl(i,j,k,nx,ny,nz
cc     .                                    ,igx,igy,igz,vdum2
cc     .                                    ,vol=.false.)
cc              enddo
cc            enddo
cc          enddo
cc
cc#if defined(vec_pot)
cc          bcnd(:,1) = u_n%array_var(IAX)%bconds
cc          bcnd(:,2) = u_n%array_var(IAY)%bconds
cc          bcnd(:,3) = u_n%array_var(IAZ)%bconds
cc
cccc          where (bcnd == -EQU) bcnd =-EXT
cc          where (bcnd == -EQU) bcnd = EXT
cc#else
cc          bcnd(:,1) = u_n%array_var(IBX)%bconds
cc          bcnd(:,2) = u_n%array_var(IBY)%bconds
cc          bcnd(:,3) = u_n%array_var(IBZ)%bconds
cccc          where (bcnd == -NEU) bcnd = EQU !Use covariant components for tangential BCs
cc          where (bcnd == -NEU) bcnd = EXT
cc#endif
cc
cc          call setBC(IAX,3,nx,ny,nz,curlcurlv,vdum2,vzeros,bcnd
cc     .              ,igx,igy,igz,iorder=2,is_cnv=.true.)
cc
cccc          curlcurlv(:,:,:,1) = vcnv(:,:,:,1)-di*curlcurlv(:,:,:,1)/rho
cccc          curlcurlv(:,:,:,2) = vcnv(:,:,:,2)-di*curlcurlv(:,:,:,2)/rho
cccc          curlcurlv(:,:,:,3) = vcnv(:,:,:,3)-di*curlcurlv(:,:,:,3)/rho
cc
cc          curlcurlv(:,:,:,1) = -di*curlcurlv(:,:,:,1)/rho
cc          curlcurlv(:,:,:,2) = -di*curlcurlv(:,:,:,2)/rho
cc          curlcurlv(:,:,:,3) = -di*curlcurlv(:,:,:,3)/rho
cc
cc          deallocate(vdum2)
cc        else
cc          if (jit == -2) vecnv_0 = vecnv !Store equilibrium info
cc
cc          curlcurlv = vecnv - vecnv_0
cc        endif
cc
cc        do k=1,nz
cc          do j=1,ny
cc            do i=1,nx
cc              div_pe(i,j,k,:)=EOM_divPe(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                                 ,curlcurlv)
cc            enddo
cc          enddo
cc        enddo
cc      else
cc        div_pe = 0d0
cc      endif
cc
ccc Electric field (covariant vector)
cc
cc#if defined(vec_pot)
cc      call find_E(nx,ny,nz,igx,igy,igz,E_ni)
cc#endif

cc      write (*,*) 'DIAG -- setupNLF'
cc        open(unit=110,file='debug.bin',form='unformatted'
cc     .     ,status='replace')
cc        call contour(E_ni(0:nx+1,0:ny+1,1,1),nx+2,ny+2,0d0,1d0
cc     .              ,0d0,1d0,0,110)
cc        call contour(E_ni(0:nx+1,0:ny+1,1,2),nx+2,ny+2,0d0,1d0
cc     .              ,0d0,1d0,1,110)
cc        call contour(E_ni(0:nx+1,0:ny+1,1,3),nx+2,ny+2,0d0,1d0
cc     .              ,0d0,1d0,1,110)
cc
cccc        stop
cc        close(110)

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
c     Frees memory
c------------------------------------------------------------------

      use nlfunction_setup

      implicit none

c Call variables

c Local variables

c Begin program

#if defined(samrai)
      call deallocAuxVariables
#endif

c End program

      end subroutine killNonlinearFunction
