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

      integer    :: i,j,k,ig,jg,kg
      real(8)    :: dxx,dyy,dzz,diffmax,eta2,nu2,beta,bnorm,norm,vnorm
      real(8)    :: k2,k2max,kv_par,kb_par2,w_cfl,w_cour
     .             ,cs2,ca2,tmp_max
      real(8)    :: x1,x2,x3,idx,idy,idz,vxx,vyy,vzz
      logical    :: cartsn

c Begin program

c Calculate maximum transport coefficient on grid

      diffmax = max(maxval(eeta),maxval(nuu),dd,chi)

c Calculate CFL

      w_cfl = 0d0

      do k=1,nzd
        do j=1,nyd
          do i=1,nxd

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            !Maximum magnetic field norm and maximum beta
            bnorm = vectorNorm(i,j,k,igx,igy,igz
     .                        ,bcnv(i,j,k,1),bcnv(i,j,k,2),bcnv(i,j,k,3)
     .                        ,.false.)

            ca2 = bnorm/rho(i,j,k)

            !Maximum kk
            idx  = 1./dx(ig)
            if (nxd == 1) idx = 0d0
            idy  = 1./dy(jg)
            if (nyd == 1) idy = 0d0
            idz  = 1./dz(kg)
            if (nzd == 1) idz = 0d0

            k2 = 4*vectorNorm(i,j,k,igx,igy,igz,idx,idy,idz,.true.)

            k2max = max(k2max,k2)

            !Maximum k.v
            vxx = px(i,j,k)/rho(i,j,k)
            vyy = py(i,j,k)/rho(i,j,k)
            vzz = pz(i,j,k)/rho(i,j,k)
            kv_par = scalarProduct(i,j,k,igx,igy,igz
     .                            ,idx,idy,idz,vxx,vyy,vzz)

            !Maximum k.B
            kb_par2 = scalarProduct(i,j,k,igx,igy,igz,idx,idy,idz
     .                 ,bcnv(i,j,k,1),bcnv(i,j,k,2),bcnv(i,j,k,3))**2
     .               /bnorm

            !Sound speed
            cs2 = a_p*gamma*tmp(i,j,k)

            !Find CFL
            w_cfl = max((icfl(cs2,ca2,k2,kb_par2,di)+abs(kv_par))
     .                  ,w_cfl)

          enddo
        enddo
      enddo

c Calculate courant number

      w_cour = diffmax*k2max

c Calculate time step

      if (w_cfl <= w_cour) then
        dt = 1d0/w_cour
cc        write (*,*) 'Courant'
      else
        dt = 1d0*(2*w_cour**2 + w_cfl**2 - w_cour*w_cfl)/
     .             (w_cour**3 + w_cfl**3)
cc        write (*,*) 'CFL'
      endif

      dt = min(dt,dtbase)

c Accumulate explicit time step for average computation

      dtexp = dtexp + dt

c End program

      contains

c     icfl
c     ##########################################################
      function icfl(cs2,ca2,k2,k2par,di)

      use math

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
      b = -ckpar*(ca2 + 2*cs2 + cs2*k2*di**2)
      c = (ca2+cs2) + ckpar*(1d0/k2 + di**2)
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

      end subroutine findExplicitDt
