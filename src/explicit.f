c explicit
c####################################################################
      subroutine explicit(vn,vnp,iter,iout)

c--------------------------------------------------------------------
c     Performs first-order predictor-corrector explicit integration.
c     Elliptic coupling vorticity/stream-function is performed after
c     corrector iteration with MG-preconditioned CG.
c--------------------------------------------------------------------

      use parameters

      use grid

      use variable_setup

      use timeStepping

      use newtongm

      use constants

      use equilibrium
      implicit none

c Call variables

      integer          :: iter,iout
      type (var_array) :: vnp,vn

c Local variables

      integer          :: i,j,k,ii,ieq

      double precision :: ftemp(ntotd)

c Externals

      real(8) :: laplace
      external   laplace

c Begin program

      iter = 0

c Predictor

      ftemp = fold

cc      call evaluateNonlinearFunction(vnp,ftemp)

      do k = 1,nzd
        do j = 1,nyd
          do i = 1,nxd

            ii = neqd*(i-1 + nxd*(j-1) + nxd*nyd*(k-1))

            do ieq=1,neqd
              vnp%array_var(ieq)%array(i,j,k) = 
     .             vn%array_var(ieq)%array(i,j,k) 
     .             -dt*(ftemp(ii+ieq) - fsrc(ii+ieq))
     .                /volume(i,j,k,1,1,1)
cc              vnp%array_var(ieq)%array(i,j,k) = fsrc(ii+ieq)
            enddo

          enddo
        enddo
      enddo

c Set BCs

      call imposeBoundaryConditions(vnp)

c Corrector

      call evaluateNonlinearFunction(vnp,ftemp)

      do k = 1,nzd
        do j = 1,nyd
          do i = 1,nxd

            ii = neqd*(i-1 + nxd*(j-1) + nxd*nyd*(k-1))

            do ieq=1,neqd
              vnp%array_var(ieq)%array(i,j,k) = 
     .              vn%array_var(ieq)%array(i,j,k) 
     .             -dt*(ftemp(ii+ieq) - fsrc(ii+ieq))
     .                /volume(i,j,k,1,1,1)
            enddo

          enddo
        enddo
      enddo

c Set BCs

      call imposeBoundaryConditions(vnp)

c End program

 200  format (/,' New_it   Av_updt    Rel_res    Damping  GMRES')
 210  format (i4,3x,1p3e11.3,i4)
 220  format ('Max newton its. exceeded; relative residual: ',1p1e10.2)
 230  format ('    Relative residual =',f7.2,' > 1.0')
 320  format ('residual vor, sf, t, vpar, n =',1p5e10.2)

      end subroutine explicit

c findExplicitDt
c####################################################################
      subroutine findExplicitDt
c--------------------------------------------------------------------
c     This subroutine calculates the explicit CFL time step
c--------------------------------------------------------------------

      use variables

      use grid

      use transport_params

      use timeStepping

      use nlfunction_setup

      implicit none

c Call variables

c Local variables

      integer :: i,j,k,ig,jg,kg
      real(8) :: dxx,dyy,dzz,diffmax,eta2,nu2,beta,bnorm
      real(8) :: kk,kv_par,kb_par,dt_cfl,dt_cour,cs,tmp_max,v_alf
      real(8) :: x1,x2,x3

c Externals

      real(8) :: res,vis
      external   res,vis

c Begin program

      dxx = minval(grid_params%dx(1:nxd))
      dyy = minval(grid_params%dy(1:nyd))
      dzz = minval(grid_params%dz(1:nzd))

c Calculate maximum transport coefficient on grid

      diffmax = max(maxval(eeta+kdiv),maxval(rho*nuu),dd,chi)

c Calculate maximum sound speed on grid

      beta = 0d0
      bnorm = 0d0

      do k=1,nzd
        do j=1,nyd
          do i=1,nxd

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            x1 = grid_params%xx(ig)
            x2 = grid_params%yy(jg)
            x3 = grid_params%zz(kg)

            bnorm = max(bnorm,(bx_cov(i,j,k)*bx(i,j,k)
     .                        +by_cov(i,j,k)*by(i,j,k)
     .                        +bz_cov(i,j,k)*bz(i,j,k))
     .                        /jacobian(x1,x2,x3))

            beta = max(beta,2*rho(i,j,k)*tmp(i,j,k))

          enddo
        enddo
      enddo

      tmp_max= maxval(abs(tmp))

      cs = sqrt(2*gamma*tmp_max)

      v_alf = sqrt(bnorm)

c Calculate corresponding CFL

      kk = 2*sqrt(1./dxx**2 + 1./dyy**2 + 1./dzz**2)

      kb_par = bx_max/dxx + by_max/dyy + bz_max/dzz
      kv_par = vx_max/dxx + vy_max/dyy + vz_max/dzz

      dt_cfl  = kv_par + sqrt(cs**2 + v_alf**2)*kk
cc      dt_cfl  = kv_par + kb_par + sqrt(cs**2 + v_alf**2)*kk
cc      dt_cfl  = sqrt(kv_par**2 + (cs**2 + 1.)*kk**2)
      dt_cour = diffmax*kk**2

cc      write (*,*) 'Sound speed',cs,' Alfven speed',v_alf
cc     .           ,' Max. diff',diffmax
cc      stop

      if (dt_cfl <= dt_cour) then
        dt = 0.8/dt_cour
cc        write (*,*) 'Courant'
      else
        dt = 0.8*(2*dt_cour**2 + dt_cfl**2 - dt_cour*dt_cfl)/
     .             (dt_cour**3 + dt_cfl**3)
cc        dt = 0.8/dt_cfl
cc        write (*,*) 'CFL'
      endif

      dt = min(dt,dtbase)

      dtexp = dtexp + dt

c End program

      return
      end
