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

cc      use debug

      use error

      implicit none

c Call variables

      type (var_array),target :: varray

c Local variables

      integer(4) :: i,j,k,ig,jg,kg

      integer(4) :: dim,loc,ibc,ieq,bctype

      real(8)    :: dh,vmax
     .             ,fx_cov(ilom:ihip,jlom:jhip,klom:khip)
     .             ,fy_cov(ilom:ihip,jlom:jhip,klom:khip)
     .             ,fz_cov(ilom:ihip,jlom:jhip,klom:khip)
     .             ,pe    (ilom:ihip,jlom:jhip,klom:khip)

c Begin program

      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Impose boundary conditions and find auxiliary quantities

      call imposeBoundaryConditions(varray,igx,igy,igz)

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

      if (minval(rho) < 0d0) then
        close (urecord)  !Close record file before exiting
        write (*,*) minval(rho)
        call pstop('setupNonlinearFunction',
     .             'Warning: negative densities are occurring')
      endif

      if (minval(tmp) < 0d0) then
        close (urecord) !Close record file before exiting
        write (*,*) minval(rho)
        call pstop('setupNonlinearFunction',
     .             'Warning: negative temperatures are occurring')
      endif

c Electron velocity

      if (di > 0d0) then
        vecnv(:,:,:,1) = vx -di*jx/rho
        vecnv(:,:,:,2) = vy -di*jy/rho
        vecnv(:,:,:,3) = vz -di*jz/rho
      else
        vecnv(:,:,:,1) = vx
        vecnv(:,:,:,2) = vy
        vecnv(:,:,:,3) = vz
      endif

c Electron forces (covariant vector)

      if (di > 0d0) then
        pe = rho*tmp

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              !Electron pressure
              call grad(i,j,k,nx,ny,nz,igx,igy,igz,pe
     .                 ,fx_cov(i,j,k),fy_cov(i,j,k),fz_cov(i,j,k))
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
      endif
        
c Viscosity (resistivity found in BC routine)

      if (nc_eom_v) then
        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              nuu  (i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
cc              eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz)
            enddo
          enddo
        enddo
      else
        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              nuu  (i,j,k) = rho(i,j,k)*vis(i,j,k,nx,ny,nz,igx,igy,igz)
cc              eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz)
            enddo
          enddo
        enddo
      endif

c Auxiliary variables

      if (di > 0d0) then    !Add electron pressure term
        ejx = - di*fx_cov/rho
        ejy = - di*fy_cov/rho
        ejz = - di*fz_cov/rho
      else
        ejx = 0d0
        ejy = 0d0
        ejz = 0d0
      endif

#if !defined(vec_pot)
      if (solenoidal) then  !Add resistivity term
        ejx = ejx + jx_cov*eeta
        ejy = ejy + jy_cov*eeta
        ejz = ejz + jz_cov*eeta
      endif

      bcnv(:,:,:,1) = bx
      bcnv(:,:,:,2) = by
      bcnv(:,:,:,3) = bz
#else
      acov(:,:,:,1) = ax
      acov(:,:,:,2) = ay
      acov(:,:,:,3) = az

      jcov(:,:,:,1) = jx_cov
      jcov(:,:,:,2) = jy_cov
      jcov(:,:,:,3) = jz_cov
#endif

      vcnv(:,:,:,1) = vx
      vcnv(:,:,:,2) = vy
      vcnv(:,:,:,3) = vz

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
