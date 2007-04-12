c setupNonlinearFunction
c#################################################################
      subroutine setupNonlinearFunction(varray)
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use parameters

      use variables

      use nlfunction_setup

      use equilibrium

      use constants

      use timeStepping

      use error

      use imposeBCinterface

      implicit none

c Call variables

      type (var_array),target :: varray

c Local variables

      integer(4) :: i,j,k,ig,jg,kg

      integer(4) :: dim,loc,ibc,ieq,bctype,bcnd(6,3)

      real(8)    :: dh,vmax,cnv(3),minrho,mintmp
     .             ,pe(ilom:ihip,jlom:jhip,klom:khip)

c Begin program

      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Impose boundary conditions and find auxiliary quantities

      call imposeBoundaryConditions(varray,igx,igy,igz)

cc      write (*,*) 'DIAG -- setupNLF'

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

c Safeguards

#if defined(petsc)
      call MPI_Allreduce(minval(rho),minrho,1,MPI_REAL8,MPI_MIN
     .                    ,MPI_COMM_WORLD,mpierr)
      call MPI_Allreduce(minval(tmp),mintmp,1,MPI_REAL8,MPI_MIN
     .                    ,MPI_COMM_WORLD,mpierr)
#else
      minrho = minval(rho)
      mintmp = minval(tmp)
#endif

      if (minrho < 0d0) then
        close (urecord)  !Close record file before exiting
        if (my_rank == 0) write (*,*) 'rho=',minrho
        call pstop('setupNonlinearFunction',
     .             'Warning: negative densities are occurring')
      endif

      if (mintmp < 0d0) then
        close (urecord) !Close record file before exiting
        if (my_rank == 0) write (*,*) 'tmp=',mintmp
        call pstop('setupNonlinearFunction',
     .             'Warning: negative temperatures are occurring')
      endif

c Electron velocity

      if (di > 0d0) then
        vecnv(:,:,:,1) = vcnv(:,:,:,1) - di*jcnv(:,:,:,1)/rho
        vecnv(:,:,:,2) = vcnv(:,:,:,2) - di*jcnv(:,:,:,2)/rho
        vecnv(:,:,:,3) = vcnv(:,:,:,3) - di*jcnv(:,:,:,3)/rho
      else
        vecnv = vcnv
      endif

c Electron forces (covariant vector)

      if (di > 0d0) then
        pe = rho*tmp

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              !Electron pressure
              call grad(i,j,k,nx,ny,nz,igx,igy,igz,pe
     .                 ,grd(i,j,k,1),grd(i,j,k,2),grd(i,j,k,3))
              !Electron hall (normally included directly via ve in vexB tem)
cc              fx_cov(i,j,k) = fx_cov(i,j,k)
cc     .                    - (jy(i,j,k)*bz(i,j,k)
cc     .                      -jz(i,j,k)*by(i,j,k))
cc     .                      /gmetric%grid(igx)%jac(i,j,k)
cc              fy_cov(i,j,k) = fy_cov(i,j,k)
cc     .                    - (jz(i,j,k)*bx(i,j,k)
cc     .                      -jx(i,j,k)*bz(i,j,k))
cc     .                      /gmetric%grid(igx)%jac(i,j,k)
cc              fz_cov(i,j,k) = fz_cov(i,j,k)
cc     .                    - (jx(i,j,k)*by(i,j,k)
cc     .                      -jy(i,j,k)*bx(i,j,k))
cc     .                      /gmetric%grid(igx)%jac(i,j,k)
            enddo
          enddo
        enddo

        grd(:,:,:,1) = di*grd(:,:,:,1)/rho
        grd(:,:,:,2) = di*grd(:,:,:,2)/rho
        grd(:,:,:,3) = di*grd(:,:,:,3)/rho
      else
        grd = 0d0
      endif
        
c Viscosity (resistivity found in BC routine)

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            nuu(i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
cc            eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

c Auxiliary variables

#if !defined(vec_pot)
      if (solenoidal) then  !Add resistivity term to grad(pe) term
        grd(:,:,:,1) = eeta*jcov(:,:,:,1) - grd(:,:,:,1)
        grd(:,:,:,2) = eeta*jcov(:,:,:,2) - grd(:,:,:,2)
        grd(:,:,:,3) = eeta*jcov(:,:,:,3) - grd(:,:,:,3)
      endif
#else
      acov(:,:,:,1) = ax
      acov(:,:,:,2) = ay
      acov(:,:,:,3) = az

      !Covariant vec-laplacian/curlcurl on A (w/ BCs)
      if (heta > 0d0) then
cc        vlap_a = -jcov
        do k=1,nz
          do j=1,ny
            do i=1,nx
cc              vlap_a(i,j,k,:) =
cc       .           veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,acov,vol=.false.)
              cnv=-curlcurl(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
     .                     ,vol=.false.)
              call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,vlap_a(i,j,k,1),vlap_a(i,j,k,2),vlap_a(i,j,k,3)
     .             ,cnv(1),cnv(2),cnv(3),.false.)
            enddo
          enddo
        enddo

        bcnd(:,1) = varray%array_var(IAX)%bconds
        bcnd(:,2) = varray%array_var(IAY)%bconds
        bcnd(:,3) = varray%array_var(IAZ)%bconds

        call setBC(IJX,3,nx,ny,nz,vdummy,vlap_a,vzeros,bcnd,igx,igy,igz
     .            ,iorder=2,is_cnv=.false.)
      endif

#endif

      pcnv(:,:,:,1) = px
      pcnv(:,:,:,2) = py
      pcnv(:,:,:,3) = pz

c SI operator

      if (k_si > 0d0) then
        b_n(:,:,:,1) = u_n%array_var(IBX)%array
        b_n(:,:,:,2) = u_n%array_var(IBY)%array
        b_n(:,:,:,3) = u_n%array_var(IBZ)%array

        p_n(:,:,:,1) = 2*u_n%array_var(IRHO)%array
     .                  *u_n%array_var(ITMP)%array

        dv_dt(:,:,:,1) = (vcnv(:,:,:,1)-u_n%array_var(IVX)%array)/dt
        dv_dt(:,:,:,2) = (vcnv(:,:,:,2)-u_n%array_var(IVY)%array)/dt
        dv_dt(:,:,:,3) = (vcnv(:,:,:,3)-u_n%array_var(IVZ)%array)/dt

        max_dv_dt = maxval(dv_dt)
      endif

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

#if !defined(vec_pot)
      nullify(rho,bx,by,bz,px,py,pz,tmp)
#else
      nullify(rho,ax,ay,az,px,py,pz,tmp)
#endif

c End program

      end subroutine killNonlinearFunction
