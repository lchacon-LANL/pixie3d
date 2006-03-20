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
     .             -0.5*dt*(ftemp(ii+ieq) - fsrc(ii+ieq))
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
 210  format (i4,3x,1p,3e11.3,i4)
 220  format ('Max newton its. exceeded; relative residual: ',1p,1e10.2)
 230  format ('    Relative residual =',f7.2,' > 1.0')
 320  format ('residual vor, sf, t, vpar, n =',1p,5e10.2)

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
      real(8)    :: k2,k2max,kv_par,kb_par2,dt_cfl,dt_cour
     .             ,cs2,ca2,tmp_max
      real(8)    :: x1,x2,x3,idx,idy,idz,vxx,vyy,vzz
      logical    :: cartsn

c Begin program

c Calculate maximum transport coefficient on grid

      diffmax = max(maxval(eeta),maxval(nuu),dd,chi)

c Calculate CFL

cc      beta   = 0d0
cc      bnorm  = 0d0
cc      vnorm  = 0d0
cc      k2max  = 0d0
cc      kv_par = 0d0
cc      kb_par2= 0d0
      dt_cfl = 0d0

      do k=1,nzd
        do j=1,nyd
          do i=1,nxd

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            !Maximum magnetic field norm and maximum beta
            bnorm = vectorNorm(i,j,k,igx,igy,igz
     .                        ,bx(i,j,k),by(i,j,k),bz(i,j,k)
     .                        ,.false.)

            ca2 = bnorm/rho(i,j,k)

            !Maximum kk
            idx  = 1./dx(ig)
            if (nxd == 1) idx = 1d-2
            idy  = 1./dy(jg)
            if (nyd == 1) idy = 1d-2
            idz  = 1./dz(kg)
            if (nzd == 1) idz = 1d-2

            k2 = 4*vectorNorm(i,j,k,igx,igy,igz,idx,idy,idz,.true.)

            k2max = max(k2max,k2)

            !Maximum k.v
            vxx = rvx(i,j,k)/rho(i,j,k)
            vyy = rvy(i,j,k)/rho(i,j,k)
            vzz = rvz(i,j,k)/rho(i,j,k)
            kv_par = scalarProduct(i,j,k,igx,igy,igz
     .                            ,idx,idy,idz,vxx,vyy,vzz)

            !Maximum k.B
            kb_par2 = scalarProduct(i,j,k,igx,igy,igz,idx,idy,idz
     .                          ,bx(i,j,k),by(i,j,k),bz(i,j,k))**2
     .               /bnorm

            !Maximum velocity field norm
cc            norm = vectorNorm(i,j,k,igx,igy,igz,vxx,vyy,vzz,.false.)
cc            vnorm = max(vnorm,norm)

            !Sound speed
            cs2 = a_p*gamma*tmp(i,j,k)

            !Find CFL
            dt_cfl = max(0.5*(icfl(cs2,ca2,k2,kb_par2,di)+abs(kv_par))
     .                  ,dt_cfl)

          enddo
        enddo
      enddo

c Calculate courant number

      dt_cour = diffmax*k2max

c Calculate time step

cc      !Inverse CFL
cc      write (*,*) 'CFL',1/dt_cfl
cccc      dt_cfl = icfl(cs2,ca2,kk**2,kb_par**2,di)
cc
cc      stop
cccc      dt_cfl = kv_par + dt_cfl   !Add Doppler shift
cc
cccccc      dt_cfl  = kv_par + sqrt(cs2 + ca2)*kk
cccc      dt_cfl  = sqrt(vnorm + cs2 + ca2)*kk
cc
cc      !Inverse courant

      !Combine both
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

      contains

c     icfl
c     ##########################################################
      function icfl(cs2,ca2,k2,k2par,di)

c     ----------------------------------------------------------
c     Finds CFL frequency
c     ----------------------------------------------------------

      implicit none

c     Call variables

      real(8)   :: icfl,cs2,ca2,k2,k2par,di

c     Local variables

      real(8)   :: a,b,c,d,root(3),ckpar

c     Begin program

      ckpar = ca2*k2par

c     Solve cubic dispersion relation for omega^2

      a =  ckpar**2*cs2
      b = -ckpar*(ca2 + 2*cs2 + ca2*cs2*k2*di**2)
      c = (ca2+cs2) + ckpar*(1d0/k2 + ca2*di**2)
      d = -1d0/k2

      root = solve_cubic(a,b,c,d)

c     Find maximum real root

      where (root > 0d0) 
        root = sqrt(root)
      elsewhere
        root = 0d0
      end where

      icfl = maxval(root)

c     End program

      end function icfl

c     solve_cubic
c     ##########################################################
      function solve_cubic(a,b,c,d) result(root)

c     ----------------------------------------------------------
c     Solves for maximum of roots of cubic polynomial:
c        a + b x +c x^2 + d x^3 = 0
c     ----------------------------------------------------------

      implicit none

c     Call variables

      real(8)   :: root(3),a,b,c,d

c     Local variables

c     Begin program

      root(1) =
     -     -c/(3.*d) - (2**0.3333333333333333*(-c**2 + 3*b*d))/
     -   (3.*d*(-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -        Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -          (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -      0.3333333333333333) + 
     -  (-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -      Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -        (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -    0.3333333333333333/(3.*2**0.3333333333333333*d)

      root(2) =
     -     -c/(3.*d) + ((1 + (0,1)*Sqrt(3.))*(-c**2 + 3*b*d))/
     -   (3.*2**0.6666666666666666*d*
     -     (-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -        Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -          (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -      0.3333333333333333) - 
     -  ((1 - (0,1)*Sqrt(3.))*
     -     (-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -        Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -          (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -      0.3333333333333333)/(6.*2**0.3333333333333333*d)

      root(3) =
     .     -c/(3.*d) + ((1 - (0,1)*Sqrt(3.))*(-c**2 + 3*b*d))/
     -   (3.*2**0.6666666666666666*d*
     -     (-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -        Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -          (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -      0.3333333333333333) - 
     -  ((1 + (0,1)*Sqrt(3.))*
     -     (-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -        Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -          (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -      0.3333333333333333)/(6.*2**0.3333333333333333*d)

c     End program

      end function solve_cubic

      end subroutine findExplicitDt
