c explicit
c####################################################################
      subroutine explicit(vn,vnp,iter,iout)

c--------------------------------------------------------------------
c     Performs first-order predictor-corrector explicit integration.
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
            enddo

          enddo
        enddo
      enddo

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
            enddo

          enddo
        enddo
      enddo

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

      integer(4) :: i,j,k,ig,jg,kg
      real(8)    :: dxx,dyy,dzz,diffmax,eta2,nu2,beta,bnorm,norm,vnorm
      real(8)    :: kk,kv_par,kb_par,dt_cfl,dt_cour,cs,tmp_max,v_alf
      real(8)    :: x1,x2,x3,idx,idy,idz,vxx,vyy,vzz
      logical    :: cartsn

c Externals

      real(8) :: res,vis
      external   res,vis

c Begin program

c Calculate maximum transport coefficient on grid

      diffmax = max(maxval(eeta),maxval(rho*nuu),dd,chi)

c Calculate maximum sound speed on grid

      beta  = 0d0
      bnorm = 0d0
      vnorm = 0d0
      kk    = 0d0
      kv_par= 0d0

      do k=1,nzd
        do j=1,nyd
          do i=1,nxd

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            !Maximum magnetic field norm and maximum beta
            norm = vectorNorm(i,j,k,igx,igy,igz
     .                       ,bx(i,j,k),by(i,j,k),bz(i,j,k)
     .                        ,.false.)
            beta = max(beta,4*rho(i,j,k)*tmp(i,j,k)/norm)
            bnorm = max(bnorm,norm/rho(i,j,k))

            !Maximum kk
            idx  = 1./dx(ig)
            idy  = 1./dy(jg)
            idz  = 1./dz(kg)
            norm = vectorNorm(i,j,k,igx,igy,igz,idx,idy,idz,.true.)
            kk   = max(kk,norm)

            !Maximum k.v
            vxx = rvx(i,j,k)/rho(i,j,k)
            vyy = rvy(i,j,k)/rho(i,j,k)
            vzz = rvz(i,j,k)/rho(i,j,k)
            norm = scalarProduct(i,j,k,igx,igy,igz
     .                          ,idx,idy,idz,vxx,vyy,vzz)
            kv_par = max(kv_par,norm)

            !Maximum velocity field norm
            norm = vectorNorm(i,j,k,igx,igy,igz,vxx,vyy,vzz,.false.)
            vnorm = max(vnorm,norm)

          enddo
        enddo
      enddo

      tmp_max= maxval(abs(tmp))

      cs = sqrt(2*gamma*tmp_max)

      kk = 2*sqrt(kk)

c Calculate corresponding CFL

cc      dt_cfl  = kv_par + sqrt(cs**2 + bnorm)*kk
      dt_cfl  = sqrt(vnorm + cs**2 + bnorm)*kk
      dt_cour = diffmax*kk**2

      if (dt_cfl <= dt_cour) then
        dt = 0.8/dt_cour
cc        write (*,*) 'Courant'
      else
        dt = 0.8*(2*dt_cour**2 + dt_cfl**2 - dt_cour*dt_cfl)/
     .             (dt_cour**3 + dt_cfl**3)
cc        write (*,*) 'CFL'
      endif

cc      write (*,*) 'Sound speed',cs,' Alfven speed',v_alf
cc     .           ,' Max. diff',diffmax
cc      write (*,*) 'dx',dxx,'dy',dyy,'kk',kk
cc      write (*,*) 'dt',dt
cc      stop

      dt = min(dt,dtbase)

      dtexp = dtexp + dt

c End program

      end subroutine findExplicitDt
