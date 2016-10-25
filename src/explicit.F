
#if !defined(samrai)
c set_explicit_dt
c####################################################################
      subroutine set_explicit_dt(dt_exp)
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

      real(8) :: dt_exp

c Local variables

      integer :: i,j,k,ig,jg,kg,igx,igy,igz,nx,ny,nz
      real(8) :: dxx,dyy,dzz,diffmax,eta2,nu2,bnorm,norm,vnorm
      real(8) :: k2,k2max,kv_par,kb_par2,kbp2max,w_cfl,w_cour
     .          ,cs2,ca2,tmp_max,roots(3)
      real(8) :: x1,x2,x3,idx,idy,idz,vxx,vyy,vzz,ldt

      character(7) :: exp_cfl_type

      INTERFACE
      function icfl(cs2,ca2,k2,k2par,di) result(root)
      real(8) :: cs2,ca2,k2,k2par,di,root(3)
      end function icfl
      END INTERFACE

c Begin program

c diag ****
cc      call scan_omega
c diag ****

      igx = 1
      igy = 1
      igz = 1

      nx = gv%gparams%nxv(igx)
      ny = gv%gparams%nyv(igy)
      nz = gv%gparams%nzv(igz)

c Calculate CFL

      w_cfl = 0d0
      k2max = 0d0
      kbp2max = 0d0

      do k=1,nz
        do j=1,ny
          do i=1,nx

            call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)

            !Maximum magnetic field norm and maximum beta
            bnorm=vectorNorm(gv%gparams,i,j,k,igx,bcnv(i,j,k,:),.false.)

            ca2 = bnorm/rho(i,j,k)

            !Maximum kk
            idx  = pi/gv%gparams%dx(ig)
            if (nx == 1) idx = 0d0
            idy  = pi/gv%gparams%dy(jg)
            if (ny == 1) idy = 0d0
            idz  = pi/gv%gparams%dz(kg)
            if (nz == 1) idz = 0d0

            k2 = vectorNorm(gv%gparams,i,j,k,igx,(/idx,idy,idz/),.true.)

            k2max = max(k2max,k2)

            !k.v
            kv_par = scalarProduct(gv%gparams,i,j,k,igx,(/idx,idy,idz/)
     .                            ,abs(vcnv(i,j,k,:)))

            !(k.B/B)^2
            if (bnorm > 0d0) then
              kb_par2=scalarProduct(gv%gparams,i,j,k,igx,(/idx,idy,idz/)
     .                             ,abs(bcnv(i,j,k,:)))**2
     .                 /bnorm
            else
              kb_par2 = 0d0
            endif

            kbp2max = max(kbp2max,kb_par2)

            !Sound speed gamma*T
            cs2 = gamma*tmp(i,j,k)

            !Find CFL frequency
            roots = icfl(cs2,ca2,k2,kb_par2,di)
            w_cfl=max(maxval(roots)+kv_par,w_cfl)

          enddo
        enddo
      enddo

c Calculate courant number

      if (.not.lagrangian) then
        w_cour = max(maxval(eeta)*k2max,maxval(nuu)*k2max,dd*k2max
     .          ,chi*k2max,di**2*k2max**2*maxval(h_eta),chi_par*kbp2max)
      else
        w_cour = max(maxval(eeta)*k2max,maxval(nuu)*k2max,dd*k2max
     .          ,chi*k2max,di**2*k2max**2*maxval(h_eta))
      endif

c Calculate time step

cc      if (w_cfl <= w_cour) then
cc        dt_exp = 2d0*pi/w_cour
cc        exp_cfl_type='Courant'
cc      else
cc        dt_exp = 0.5*pi*(2*w_cour**2 + w_cfl**2 - w_cour*w_cfl)
cc     .                 /(w_cour**3 + w_cfl**3)
cc        exp_cfl_type='CFL'
cc      endif
cc
cc      if (my_rank == 0.and.itime == 0)
cc     .     write (*,*) 'Fast time scale,'
cc     .                 ,trim(exp_cfl_type),'='
cc     .                 ,max(w_cour,w_cfl)

      dt_exp = 0.5*pi*(2*w_cour**2 + w_cfl**2 - w_cour*w_cfl)
     .           /(w_cour**3 + w_cfl**3)

c Synchronize time step among all processors

#if defined(petsc)
      ldt = dt_exp
      call MPI_Allreduce(ldt,dt_exp,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_MIN,MPI_COMM_WORLD,mpierr)
#endif

c Compare with input time step

      dt_exp = min(dt_exp,dtbase)

c End program

      contains

c     scan_omega
c     #############################################################
      subroutine scan_omega

c     -------------------------------------------------------------
c     Scans frequency for Hall MHD dispersion relation. On input:
c       * beta: plasma beta
c       * alpha = k/k||
c     -------------------------------------------------------------
      implicit none

c     Call variables

c     Local variables

      integer :: ndiv,i
      real(8) :: beta,theta,k_max,delta,rho_s
      real(8) :: cs2,ca2,k2,k2par,di,root(3)

      character(100) :: filename

c     Begin program

      di   = 1d0

cc      !Schnack's tests case (PoP 2006)
cc      beta = 0.2
cc      theta = 0.46  !in pi units

      !Large guide field case
      beta = 1d-4
      theta = 0.46  !in pi units

      ndiv = 500

      write(filename,'(a,1p,e8.2,a,0p,f4.2,a)')
     .     'omega_beta=',beta,'_theta=',theta,'pi.txt'

      open(unit=6666,file=trim(filename),status='unknown')

      cs2 = 0.5*beta

      rho_s = di*sqrt(cs2)

      k_max = 100.0/rho_s

      delta = k_max/ndiv

      do i=1,ndiv
        k2 = (i*delta)**2
        kb_par2 = k2*cos(theta*pi)**2
      
        root = icfl(cs2,1d0,k2,kb_par2,di)

        write (6666,*) sqrt(k2)*rho_s,root
      enddo

      close(6666)

      stop

      end subroutine scan_omega

      end subroutine set_explicit_dt

c set_user_dt
c####################################################################
      subroutine set_user_dt(dt)
c--------------------------------------------------------------------
c     This subroutine computes user-defined timestep, based on flow
c     CFL constraints.
c--------------------------------------------------------------------

      use variables

      use nlfunction_setup, ONLY: IVX,IVY,IVZ

      implicit none

c Call variables

      real(8) :: dt

c Local variables

      integer :: i,j,k,ig,jg,kg,igx,igy,igz,nx,ny,nz
      real(8) :: kv_par,w_cfl,dt_cfl
      real(8) :: ih(3),ldt,v(3)

c Begin program

      igx = 1
      igy = 1
      igz = 1

      nx = gv%gparams%nxv(igx)
      ny = gv%gparams%nyv(igy)
      nz = gv%gparams%nzv(igz)

c Calculate flow-based CFL

      w_cfl = 0d0

      do k=1,nz
        do j=1,ny
          do i=1,nx

            call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)

            !Maximum kk
            ih(1)  = pi/gv%gparams%dx(ig)
            if (nx == 1) ih(1) = 0d0
            ih(2)  = pi/gv%gparams%dy(jg)
            if (ny == 1) ih(2) = 0d0
            ih(3)  = pi/gv%gparams%dz(kg)
            if (nz == 1) ih(3) = 0d0

            !k.v
            v(1) = u_np%array_var(IVX)%array(i,j,k)
            v(2) = u_np%array_var(IVY)%array(i,j,k)
            v(3) = u_np%array_var(IVZ)%array(i,j,k)
            kv_par = scalarProduct(gv%gparams,i,j,k,igx,ih,v)

            !Find flow CFL
            w_cfl=max(kv_par,w_cfl)

          enddo
        enddo
      enddo

c Calculate flow-CFL time step

      if (w_cfl > 0d0) then
        dt_cfl = 0.15*pi/w_cfl
      else
        dt_cfl = 1d10
      endif

c Synchronize timestep among all processors

#if defined(petsc)
      ldt = dt_cfl
      call MPI_Allreduce(ldt,dt_cfl,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_MIN,MPI_COMM_WORLD,mpierr)
#endif

c Compare with input time step

      dt = min(dt,dt_cfl)

      if (my_rank==0.and.dt==dt_cfl)
     .     write (*,'(a,1p,1e11.4)')
     .     '   >>>> Flow-CFL-limited dt=',dt_cfl

c End program

      end subroutine set_user_dt

#else

c findExplicitDt
c####################################################################
       subroutine findExplicitDt(patch_var,dt_exp)
c--------------------------------------------------------------------
c     This subroutine calculates the explicit CFL time step
c--------------------------------------------------------------------

      use variables

      use grid

      use transport_params

      use timeStepping

      use auxiliaryVariables

      implicit none

c Call variables

      type(patch), TARGET :: patch_var

      real(8) :: dt_exp

c Local variables

      integer :: i,j,k,ig,jg,kg,igx,igy,igz,nx,ny,nz
      real(8) :: dxx,dyy,dzz,diffmax,eta2,nu2,bnorm,norm,vnorm
      real(8) :: k2,k2max,kv_par,kb_par2,kbp2max,w_cfl,w_cour
     .          ,cs2,ca2,tmp_max,roots(3)
      real(8) :: x1,x2,x3,idx,idy,idz,ldt

      INTERFACE
      function icfl(cs2,ca2,k2,k2par,di) result(root)
        real(8) :: cs2,ca2,k2,k2par,di,root(3)
      end function icfl
      END INTERFACE

c Begin program

      gv => patch_var
      tmp   => gv%aux%var_list(ITMP_AUX)%array
      rho   => gv%aux%var_list(IRHO_AUX)%array
      vcnv  => gv%aux%vec_list(IVCNV)%vec
      bcnv  => gv%aux%vec_list(IBCNV)%vec
      eeta  => gv%aux%var_list(IETA)%array
      h_eta => gv%aux%var_list(IHETA)%array
      nuu   => gv%aux%var_list(INU)%array

       igx = 1
       igy = 1
       igz = 1

       nx = gv%gparams%nxv(igx)
       ny = gv%gparams%nyv(igy)
       nz = gv%gparams%nzv(igz)

c Calculate CFL

       w_cfl = 0d0
       k2max = 0d0
       kbp2max = 0d0

       do k=1,nz
         do j=1,ny
           do i=1,nx

             call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)

             !Maximum magnetic field norm and maximum beta
             bnorm = vectorNorm(gv%gparams,i,j,k,igx,bcnv(i,j,k,:)
     .                         ,.false.)

             ca2 = bnorm/rho(i,j,k)

             !Maximum kk
             idx  = pi/gv%gparams%dx(ig)
             if (nx == 1) idx = 0d0
             idy  = pi/gv%gparams%dy(jg)
             if (ny == 1) idy = 0d0
             idz  = pi/gv%gparams%dz(kg)
             if (nz == 1) idz = 0d0

             k2 = vectorNorm(gv%gparams,i,j,k,igx,(/idx,idy,idz/)
     .                      ,.true.)

             k2max = max(k2max,k2)

             !k.v
             kv_par = scalarProduct(gv%gparams,i,j,k,igx,(/idx,idy,idz/)
     .                            ,abs(vcnv(i,j,k,:)))

             !(k.B/B)^2
             if (bnorm > 0d0) then
               kb_par2 = scalarProduct(gv%gparams,i,j,k,igx
     .                           ,(/idx,idy,idz/),abs(bcnv(i,j,k,:)))**2
     .                 /bnorm
             else
               kb_par2 = 0d0
             endif

             kbp2max = max(kbp2max,kb_par2)

             !Sound speed gamma*p/rho
             cs2 = gamma*tmp(i,j,k)

             !Find CFL frequency
             roots = icfl(cs2,ca2,k2,kb_par2,di)
             w_cfl=max(maxval(roots)+kv_par,w_cfl)

           enddo
         enddo
       enddo

c Calculate courant number

       if (di > 0d0) then
         w_cour = max(maxval(eeta)*k2max
     .               ,maxval(nuu+h_eta)*k2max
     .               ,dd*k2max
     .               ,chi*k2max
     .               ,maxval(di**2*h_eta*k2max**2)
     .               ,chi_par*kbp2max)
       else
         w_cour = max(maxval(eeta)*k2max,maxval(nuu)*k2max,dd*k2max
     .               ,chi*k2max,chi_par*kbp2max)
       endif

c Calculate time step

cc      if (w_cfl <= w_cour) then
cc        dt = 2d0*pi/w_cour
cc        exp_cfl_type='Courant'
cc      else
         dt_exp = 0.5*pi*(2*w_cour**2 + w_cfl**2 - w_cour*w_cfl)
     .              /(w_cour**3 + w_cfl**3)
cc        exp_cfl_type='CFL'
cc      endif

c End program

       end subroutine findExplicitDt

#endif
