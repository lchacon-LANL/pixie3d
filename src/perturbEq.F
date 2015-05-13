c     perturbEquilibrium
c     #################################################################
      subroutine perturbEquilibrium(varr)

c     -----------------------------------------------------------------
c     Perturbs equilibrium quantity in array0 with a sinusoidal
c     perturbation of magnitude perturb, and introduces it in array.
c     -----------------------------------------------------------------

      use grid

      use var_setup

      use timeStepping

      use app_iosetup

      use equilibrium

      use transport_params

      implicit none

c     Call variables

      type(var_array),pointer :: varr

c     Local variables

      integer :: i,j,k,ig,jg,kg,igx,igy,igz,ieq
      real(8) :: x1,y1,z1,jac,kx,ky,kk,eps
     .          ,vplus,vminus,vshear,cs2,kdotB,cth
     .          ,vxx,vyy,vzz,ivxx,ivyy,ivzz,xi,ixi,yy,iyy
     .          ,bxx,byy,bzz,ibxx,ibyy,ibzz,axx,ayy,azz
     .          ,w_fast,w_slow,w_shr,freq,roots(3),X,iX
      real(8) :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      integer :: ip,im,jp,jm,kp,km
      real(8) :: dxx,dyy,dzz,RR

      integer :: mm,nn

cc      real(8) :: kbesel, beseli, dbeseli, dbeselip, dbeselim

      logical :: check

      integer :: ierr
      character(50) :: command

      save  :: w_fast,w_slow,w_shr,freq

! Daniele, 04 May 2015
! to allow for multiple perturbed harmonics in the ohtor case
      integer :: npw,ipw
      integer,allocatable,dimension(:) :: mm_arr, nn_arr
      real(8),allocatable,dimension(:) :: br_pert_bc_arr
     .     ,br_pert_phase_arr

      INTERFACE
        function icfl(cs2,ca2,k2,k2par,di) result(root)
        real(8)   :: cs2,ca2,k2,k2par,di,root(3)
        end function icfl
      END INTERFACE

      INTERFACE
        subroutine perturbEquilibrium_ppnch_bc(varray,br_pert,br_phase)
        use var_setup
        type(var_array),pointer :: varray
        real(8) :: br_pert,br_phase
        end subroutine perturbEquilibrium_ppnch_bc
      END INTERFACE

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      select case(equil)
      case ('2dkai','3dkai')

        command = 'test -f '//trim(prt_file)
        ierr = shellcmd(trim(command))

        if (ierr == 0) then
          call perturbEquilibrium_kaitm
        else
          do ieq=1,varr%nvar
            call perturbEquilibrium_def(varr%array_var(ieq)%array
     .                                 ,varr%array_var(ieq)%bconds
     .                                 ,pert(ieq),.false.)
          enddo
        endif

      case ('ppnch','ppn3d','ppnsl','ppnst','p3nsl','2dhel')

        if (is_file(prt_file)) then
          call perturbEquilibrium_ppnch_ext
cc          do ieq=1,varr%nvar
cc            call perturbEquilibrium_ppnch(varr%array_var(ieq)%array
cc     .                                   ,varr%array_var(ieq)%bconds
cc     .                                   ,pert(ieq),ieq)
cc          enddo
        else
          do ieq=1,varr%nvar
            call perturbEquilibrium_def(varr%array_var(ieq)%array
     .                                 ,varr%array_var(ieq)%bconds
     .                                 ,pert(ieq),.false.)
          enddo

          if (perturb_rstrt) then
            !Substract old one
            call perturbEquilibrium_ppnch_bc(varr,-br_pert_bc_old
     .                                      ,br_pert_phase_old)

            !Add new one
            call perturbEquilibrium_ppnch_bc(varr,br_pert_bc
     .                                      ,br_pert_phase)
          else
cc            write(*,*) 'calling ppnch_br'
            call perturbEquilibrium_ppnch_bc(varr,br_pert_bc
     .                                      ,br_pert_phase)
          endif

          !Save current perturb state in input file for restart
          call write_real_field_to_inputf(br_pert_bc,'br_pert_bc_old')
          call write_real_field_to_inputf(br_pert_phase
     $                                   ,'br_pert_phase_old')

        endif

! Daniele, 03/02/2012. case 'ohtor' has been written from scratch
      case ('ohtor')

        do ieq=1,varr%nvar
          call perturbEquilibrium_def(varr%array_var(ieq)%array
     .                               ,varr%array_var(ieq)%bconds
     .                               ,pert(ieq),.false.)
        enddo

        !Perturb normal B-field component ensuring div(B)=0 numerically
        if (br_pert_bc /= 0d0) then

          RR = grid_params%params(1)
cc          mm = grid_params%params(4)
cc          nn = grid_params%params(5)
!          mm = nh2
!          nn =-nh3

! Daniele, 04 May 2015
! to allow for multiple perturbed harmonics in the ohtor case
           if (is_file(prt_file)) then

              if (my_rank == 0) then
                 write (*,*)
                 write (*,*) 'Reading external perturbation file '
     .                ,trim(prt_file)
              endif

              open(unit=u_prtf,file=trim(prt_file),form='formatted')

              read(u_prtf,*) npw

              if (my_rank == 0) write(*,*) 'Perturbed harmonics=',npw

              allocate(mm_arr(npw)
     .                ,nn_arr(npw)
     .                ,br_pert_bc_arr(npw)
     .                ,br_pert_phase_arr(npw))

              read(u_prtf,*) mm_arr
              read(u_prtf,*) nn_arr
              read(u_prtf,*) br_pert_bc_arr
              read(u_prtf,*) br_pert_phase_arr

              if (my_rank == 0) then
                 write(*,*) 'mm array=',mm_arr
                 write(*,*) 'nn array=',nn_arr
                 write(*,*) 'br_pert_bc array=',br_pert_bc_arr
                 write(*,*) 'br_pert_phase array=',br_pert_phase_arr
              endif

              close(u_prtf)

           else

              npw = 1

              allocate(mm_arr(npw)
     .                ,nn_arr(npw)
     .                ,br_pert_bc_arr(npw)
     .                ,br_pert_phase_arr(npw))

              mm_arr(1) = nh2
              nn_arr(1) =-nh3
              br_pert_bc_arr(1) = br_pert_bc
              br_pert_phase_arr(1) = br_pert_phase

           endif

           do ipw = 1,npw

              mm = mm_arr(ipw)
              nn = nn_arr(ipw)
              br_pert_bc = br_pert_bc_arr(ipw)
              br_pert_phase = br_pert_phase_arr(ipw)

          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                call getCurvilinearCoordinates(i,j,k,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)

                ip = ig+1
                im = ig-1
                jp = jg+1
                jm = jg-1
                kp = kg+1
                km = kg-1

                dxx = 2*grid_params%dxh(ig)
                dyy = 2*grid_params%dyh(jg)
                dzz = 2*grid_params%dzh(kg)

! Daniele, 19 Dec. 2014
! divergence-free perturbation with b^3=0 (and limited to the n=0 case)
! used for the n=0 case when the divergence-free perturbation with b^2=0 fails
! to be extended to generic br_pert_phase (now only particular case br_pert_phase=0.)
! to be possibly extended to generic n
! b^1=r^m instead of r^(m+2) similar to the helical/cylindrical case
                if (nn == 0) then

                   varr%array_var(IBX)%array(i,j,k) =
     .                  varr%array_var(IBX)%array(i,j,k)
     .                  +x1**mm*br_pert_bc*RR
     .                  *(cos(mm*grid_params%yy(jp))
     .                  -cos(mm*grid_params%yy(jm)))
     .                  /dyy/mm !Derivative done numerically to preserve div(B)=0

                   varr%array_var(IBY)%array(i,j,k) =
     .                  varr%array_var(IBY)%array(i,j,k)
     .                  -br_pert_bc*RR
     .                  *(grid_params%xx(ip)**mm
     .                  -grid_params%xx(im)**mm)
     .                  /dxx !Derivative done numerically to preserve div(B)=0
     .                  *cos(mm*y1)/mm

                   if (mm /= -1) then

                      varr%array_var(IBX)%array(i,j,k) =
     .                     varr%array_var(IBX)%array(i,j,k)
     .                     +0.5*x1**(mm+1)*br_pert_bc
     .                     *(cos((mm+1)*grid_params%yy(jp))
     .                     -cos((mm+1)*grid_params%yy(jm)))
     .                     /dyy/(mm+1) !Derivative done numerically to preserve div(B)=0

                      varr%array_var(IBY)%array(i,j,k) =
     .                     varr%array_var(IBY)%array(i,j,k)
     .                     -0.5*br_pert_bc
     .                     *(grid_params%xx(ip)**(mm+1)
     .                     -grid_params%xx(im)**(mm+1))
     .                     /dxx !Derivative done numerically to preserve div(B)=0
     .                     *cos((mm+1)*y1)/(mm+1)

                   endif

                   if (mm /= 1) then
                      varr%array_var(IBX)%array(i,j,k) =
     .                     varr%array_var(IBX)%array(i,j,k)
     .                     +0.5*x1**(mm+1)*br_pert_bc
     .                     *(cos((mm-1)*grid_params%yy(jp))
     .                     -cos((mm-1)*grid_params%yy(jm)))
     .                     /dyy/(mm-1) !Derivative done numerically to preserve div(B)=0

                      varr%array_var(IBY)%array(i,j,k) =
     .                     varr%array_var(IBY)%array(i,j,k)
     .                     -0.5*br_pert_bc
     .                     *(grid_params%xx(ip)**(mm+1)
     .                     -grid_params%xx(im)**(mm+1))
     .                     /dxx !Derivative done numerically to preserve div(B)=0
     .                     *cos((mm-1)*y1)/(mm-1)

                   endif

                else

! Daniele, 19 Dec. 2014
! divergence-free perturbation with b^3=0 (fails in the n=0 case)
                varr%array_var(IBX)%array(i,j,k) =
     .               varr%array_var(IBX)%array(i,j,k)
!     .               +x1**(mm+10)*br_pert_bc*(RR+x1*cos(y1)) ! only even powers work!
     .               +x1**(mm+2)*br_pert_bc*(RR+x1*cos(y1)) ! only even powers work!
!     .               +x1**mm*br_pert_bc*(RR+x1*cos(y1)) ! only even powers work!
     .                 *(cos(mm*y1+nn*grid_params%zz(kp)+br_pert_phase)
     .                  -cos(mm*y1+nn*grid_params%zz(km)+br_pert_phase))
     .                 /dzz/nn !Derivative done numerically to preserve div(B)=0
                varr%array_var(IBZ)%array(i,j,k) =
     .               varr%array_var(IBZ)%array(i,j,k)
     .               -br_pert_bc*
!     .               (grid_params%xx(ip)**(mm+10)
     .               (grid_params%xx(ip)**(mm+2)
!     .               (grid_params%xx(ip)**mm
     .                *(RR+grid_params%xx(ip)*cos(y1))
!     .               -grid_params%xx(im)**(mm+10)
     .               -grid_params%xx(im)**(mm+2)
!     .               -grid_params%xx(im)**mm
     .                *(RR+grid_params%xx(im)*cos(y1)))
     .                 /dxx !Derivative done numerically to preserve div(B)=0
     .               *cos(mm*y1+nn*z1+br_pert_phase)/nn

                endif

              enddo
            enddo
          enddo

          enddo

       endif

      case ('msw')

        kx = 2*pi*nh1/(xmax-xmin)
        ky = 2*pi*nh2/(ymax-ymin)

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              jac = 1d0

              varr%array_var(IRHO)%array(i,j,k) =
     .                     varr%array_var(IRHO)%array(i,j,k)
     .                     +jac*1d-3*cos(kx*(x1-xmin)+ky*(y1-ymin))

              varr%array_var(ITMP)%array(i,j,k) =
     .                     varr%array_var(ITMP)%array(i,j,k)
     .                     +a_p*jac*1d-3*cos(kx*(x1-xmin)+ky*(y1-ymin))

#if defined(vec_pot)
              varr%array_var(IAY)%array(i,j,k) =
     .                     varr%array_var(IAY)%array(i,j,k)
     .                     + jac*1d-3/kx
     .                          *sin(kx*(x1-xmin)+ky*(y1-ymin))
#else
              varr%array_var(IBZ)%array(i,j,k) =
     .                     varr%array_var(IBZ)%array(i,j,k)
     .                     + jac*1d-3*cos(kx*(x1-xmin)+ky*(y1-ymin))
#endif
            enddo
          enddo
        enddo

      case ('3wv')

        kx = 2*pi*nh1/(xmax-xmin)
        ky = 2*pi*nh2/(ymax-ymin)

        kdotB = kx*B0(1)+ky*B0(2)

        kk = sqrt(kx**2+ky**2)

        cs2 = 0.5*gamma*beta

        cth = kdotB/kk

        eps = 1d-5

        !Find MHD frequencies  (do this only once)
        if (ieq == 1) then
          roots = icfl(cs2,1d0,kk**2,kdotB**2,di)

          w_fast = maxval(roots)
          w_slow = minval(roots)
          w_shr  = fmed(roots(1),roots(2),roots(3))

          write (*,*) 'B angle =',grid_params%params(1)
          write (*,*)
          write (*,*) 'shear fq=',w_shr
          write (*,*) 'fast  fq=',w_fast
          write (*,*) 'slow  fq=',w_slow
          write (*,*)

          !Select frequency
          select case(nint(grid_params%params(3)))
          case(1)
            if (my_rank == 0)
     .           write (*,*) '===> Selecting fast freq...'
            freq = w_fast
          case(2)
            if (my_rank == 0)
     .           write (*,*) '===> Selecting slow freq...'
            freq = w_slow
          case(3)
            if (my_rank == 0)
     .           write (*,*) '===> Selecting shear freq...'
            freq = w_shr
          case default
            call pstop('perturbEquilibrium','Select frequency!')
          end select

          check = (dt <= pi/16./freq)

          if (.not.check) then
            if (my_rank == 0) write (*,*) 'Time step limit=',pi/16./freq
            call pstop('perturbEquilibrium'
     .                ,'Time step must be < than this!')
          endif

        endif

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              if ((cs2*kdotB*ky/freq**2-B0(2)) /= 0d0) then
                vxx = eps*cos(kx*(x1-xmin)+ky*(y1-ymin))
                vyy =-(cs2*kdotB*kx/freq**2-B0(1))*vxx
     .               /(cs2*kdotB*ky/freq**2-B0(2))
              else
                vyy = eps*cos(kx*(x1-xmin)+ky*(y1-ymin))
                vxx =-(cs2*kdotB*ky/freq**2-B0(2))*vyy
     .               /(cs2*kdotB*kx/freq**2-B0(1))
              endif

              xi = kx*vxx + ky*vyy
              X  = cs2*kdotB*xi/freq**2

              if ((cs2*kdotB*ky/freq**2-B0(2)) /= 0d0) then
                ivxx = eps*sin(kx*(x1-xmin)+ky*(y1-ymin))
                ivyy =-(cs2*kdotB*kx/freq**2-B0(1))*ivxx
     .                /(cs2*kdotB*ky/freq**2-B0(2))
              else
                ivyy = eps*sin(kx*(x1-xmin)+ky*(y1-ymin))
                ivxx =-(cs2*kdotB*ky/freq**2-B0(2))*ivyy
     .                /(cs2*kdotB*kx/freq**2-B0(1))
              endif

              ixi = kx*ivxx + ky*ivyy
              iX  = cs2*kdotB*ixi/freq**2

              yy  = iX*kdotB*di*(freq**2/cs2-kk**2)/(freq-kdotB**2/freq)
              iyy = -X*kdotB*di*(freq**2/cs2-kk**2)/(freq-kdotB**2/freq)

              vzz =( yy + B0(3)*(ky* vxx - kx* vyy))/(ky*B0(1)-kx*B0(2))
              ivzz=(iyy + B0(3)*(ky*ivxx - kx*ivyy))/(ky*B0(1)-kx*B0(2))

              bxx = (kdotB*vxx-xi*B0(1))/freq +di*ky*ivzz
              byy = (kdotB*vyy-xi*B0(2))/freq -di*kx*ivzz
              bzz = (kdotB*vzz-xi*B0(3))/freq +di*(kx*ivyy-ky*ivxx)

#if defined(vec_pot)
              ibxx = (kdotB*ivxx-ixi*B0(1))/freq -di*ky*vzz
              ibyy = (kdotB*ivyy-ixi*B0(2))/freq +di*kx*vzz
              ibzz = (kdotB*ivzz-ixi*B0(3))/freq -di*(kx*vyy-ky*vxx)

              if (kx /= 0d0) then
                axx = 0d0
                azz =-ibyy/kx
                ayy = ibzz/kx
              elseif (ky /= 0d0) then
                ayy = 0d0
                azz = ibxx/ky
                axx =-ibzz/ky
              endif
#endif

              varr%array_var(IRHO)%array(i,j,k) =
     .               varr%array_var(IRHO)%array(i,j,k) - xi/freq

              varr%array_var(ITMP)%array(i,j,k) =
     .               varr%array_var(ITMP)%array(i,j,k)
     .             - cs2*xi/freq !For total pressure

#if !defined(vec_pot)
              varr%array_var(IBX)%array(i,j,k)
     .               = varr%array_var(IBX)%array(i,j,k) + bxx

              varr%array_var(IBY)%array(i,j,k)
     .               = varr%array_var(IBY)%array(i,j,k) + byy

              varr%array_var(IBZ)%array(i,j,k)
     .               = varr%array_var(IBZ)%array(i,j,k) + bzz
#else
              varr%array_var(IAX)%array(i,j,k)
     .               = varr%array_var(IAX)%array(i,j,k) + axx

              varr%array_var(IAY)%array(i,j,k)
     .               = varr%array_var(IAY)%array(i,j,k) + ayy

              varr%array_var(IAZ)%array(i,j,k)
     .               = varr%array_var(IAZ)%array(i,j,k) + azz
#endif
              varr%array_var(IVX)%array(i,j,k)
     .               = varr%array_var(IVX)%array(i,j,k) + vxx

              varr%array_var(IVY)%array(i,j,k)
     .               = varr%array_var(IVY)%array(i,j,k) + vyy

              varr%array_var(IVZ)%array(i,j,k)
     .               = varr%array_var(IVZ)%array(i,j,k) + vzz
            enddo
          enddo
        enddo

      case ('whslr')

        kx = 2*pi*nh1/(xmax-xmin)
        ky = 2*pi*nh2/(ymax-ymin)

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              jac = gmetric%grid(igx)%jac(i,j,k)
#if !defined(vec_pot) 
              varr%array_var(IBY)%array(i,j,k) =
     .             varr%array_var(IBY)%array(i,j,k)
     .           + 1d-3*cos(kx*(x1-xmin)+ky*(y1-ymin))*jac
#else
              varr%array_var(IAZ)%array(i,j,k) =
     .             varr%array_var(IAZ)%array(i,j,k)
     .           - 1d-3/kx*sin(kx*(x1-xmin)+ky*(y1-ymin))
#endif
              varr%array_var(IVZ)%array(i,j,k) =
     .             varr%array_var(IVZ)%array(i,j,k)
     .           + 1d-3/kx/di*sin(kx*(x1-xmin)+ky*(y1-ymin))*jac
            enddo
          enddo
        enddo

      case ('kaw')

        kx = 2*pi*nh1/(xmax-xmin)
        ky = 2*pi*nh2/(ymax-ymin)
        kk = sqrt(kx**2+ky**2)

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              varr%array_var(IVY)%array(i,j,k)
     .               = varr%array_var(IVY)%array(i,j,k)
     .                       - 1d-7*kx/ky*(1.-(di*kk)**2)
     .                             *cos(kx*(x1-xmin)+ky*(y1-ymin))

              varr%array_var(IVX)%array(i,j,k)
     .               = varr%array_var(IVX)%array(i,j,k)
     .                       + 1d-7*cos(kx*(x1-xmin)+ky*(y1-ymin))
            enddo
          enddo
        enddo

      case ('ic','ic-hk','dtm','gem','tmcar','tmsin','kai1')

        do ieq=1,varr%nvar
          call perturbEquilibrium_psi(varr%array_var(ieq)%array
     .                               ,varr%array_var(ieq)%bconds
     .                               ,pert(ieq),ieq)
        enddo

cc      case ('3dpar')
cc
cc        call perturbEquilibrium_3dpar(array,bcs,pert(ieq),ieq)

      case default

        do ieq=1,varr%nvar
          call perturbEquilibrium_def(varr%array_var(ieq)%array
     .                               ,varr%array_var(ieq)%bconds
     .                               ,pert(ieq),.false.)
        enddo

      end select

c     End program

      contains

c     perturbEquilibrium_def
c     #################################################################
      subroutine perturbEquilibrium_def(array,bcs,perturb,cartesian)

c     -----------------------------------------------------------------
c     Perturbs equilibrium quantity in array0 with a sinusoidal
c     perturbation of magnitude perturb, and introduces it in array.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: bcs(6)
      real(8) :: perturb
      real(8) :: array (ilom:ihip,jlom:jhip,klom:khip)
      logical :: cartesian

c     Local variables

      integer :: i,j,k,ig,jg,kg,iglobal,ipw1,ipw2,ipw3,mode
      real(8) :: x1,y1,z1
      real(8) :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

c     Begin program

      if (cartesian) then

        do ipw1=0,npw1-1
          do ipw2=0,npw2-1
            do ipw3=0,npw3-1

              mode = nh1+sign(ipw1,nh1)
              do i = ilom,ihip
                call getCartesianCoordinates(i,jlo,klo,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)
                fx(i) = factor(xmin,xmax,x1,bcs(1:2),mode,odd(1))
              enddo

              mode = nh2+sign(ipw2,nh2)
              do j = jlom,jhip
                call getCartesianCoordinates(ilo,j,klo,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)
                fy(j) = factor(ymin,ymax,y1,bcs(3:4),mode,odd(2))
              enddo

              mode = nh3+sign(ipw3,nh3)
              do k = klom,khip
                call getCartesianCoordinates(ilo,jlo,k,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)
                fz(k) = factor(zmin,zmax,z1,bcs(5:6),mode,odd(3))
              enddo

              do k = klom,khip
                do j = jlom,jhip
                  do i = ilom,ihip
                    array(i,j,k) =array(i,j,k)+perturb*fx(i)*fy(j)*fz(k)
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo

      else

        do ipw1=0,npw1-1
          do ipw2=0,npw2-1
            do ipw3=0,npw3-1

              mode = nh1+sign(ipw1,nh1)
              do i = ilom,ihip
                call getCurvilinearCoordinates(i,jlo,klo,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)
                fx(i) = factor(xmin,xmax,x1,bcs(1:2),mode,odd(1))
              enddo

              mode = nh2+sign(ipw2,nh2)
              do j = jlom,jhip
                call getCurvilinearCoordinates(ilo,j,klo,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)
                fy(j) = factor(ymin,ymax,y1,bcs(3:4),mode,odd(2))
              enddo

              mode = nh3+sign(ipw3,nh3)
              do k = klom,khip
                call getCurvilinearCoordinates(ilo,jlo,k,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)
                fz(k) = factor(zmin,zmax,z1,bcs(5:6),mode,odd(3))
              enddo

              do k = klom,khip
                do j = jlom,jhip
                  do i = ilom,ihip
                    array(i,j,k) =array(i,j,k)+perturb*fx(i)*fy(j)*fz(k)
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo

      endif

c     End program

      end subroutine perturbEquilibrium_def

c     perturbEquilibrium_kaitm
c     #################################################################
      subroutine perturbEquilibrium_kaitm

c     -----------------------------------------------------------------
c     Perturbs equilibrium quantities according to eigenmodes computed
c     by Kai Germaschewski (note: eigenvectors are not correct, and
c     result in evolution by PIXIE3D to correct them).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      integer    :: i,j,k,ig,jg,kg,iglobal
      real(8)    :: x1,y1,z1,kk,mm
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      real(8) :: dum,rr(nxd,neqd),ii(nxd,neqd),pert(neqd)

c     Begin program

      mm = grid_params%params(1)
      kk = grid_params%params(2)

      if (my_rank == 0) write (*,*) 'Reading KAI TM perturbations'

      open(unit=u_prtf,file=trim(prt_file),status='old')

      do i=1,nxd
        read(u_prtf,*) dum,rr(i,1),rr(i,8),rr(i,2:7)
     .                 ,ii(i,1),ii(i,8),ii(i,2:7)
cc        write(*,*) dum,rr(1),rr(8),rr(2:7),ii(1),ii(8),ii(2:7)
      enddo

      do k = klom,khip
        do j = jlom,jhip
          do i = ilo,ihi
            call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

            iglobal = i + grid_params%ilo(igx) - 1

            pert = rr(iglobal,:)*cos(y1)+ii(iglobal,:)*sin(y1)

            pert(IRHO) = x1*pert(IRHO)

            pert(IVX) = x1*pert(IVX)
            if (equil == '2dkai') pert(IVY)=pert(IVY)+kk*x1/mm*pert(IVZ)
            if (alt__eom()) pert(IVY) = x1*pert(IVY)
            pert(IVZ) = x1*pert(IVZ)

            pert(IBX) = x1*pert(IBX)
            if (equil == '2dkai') pert(IBY)=pert(IBY)+kk*x1/mm*pert(IBZ)
            pert(IBZ) = x1*pert(IBZ)

            do ieq=1,varr%nvar-1
              varr%array_var(ieq)%array(i,j,k) =
     .               varr%array_var(ieq)%array(i,j,k) + pert(ieq)
            enddo

            pert(ITMP) = pert(ITMP)*a_p  !Transform to total temperature

            varr%array_var(ITMP)%array(i,j,k) =
     .               varr%array_var(ITMP)%array(i,j,k)*(1d0 
     .               + pert(IRHO)/varr%array_var(IRHO)%array(i,j,k))
     .               + pert(ITMP)*varr%array_var(IRHO)%array(i,j,k)

          enddo
        enddo
      enddo

      close (u_prtf)

c     End program

      end subroutine perturbEquilibrium_kaitm

c     perturbEquilibrium_ppnch_eig
c     #################################################################
      subroutine perturbEquilibrium_ppnch_eig(array,bcs,perturb,ieq)

c     -----------------------------------------------------------------
c     Perturbs equilibrium quantity in array0 with a sinusoidal
c     perturbation of magnitude perturb, and introduces it in array.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: bcs(6),ieq
      real(8)    :: perturb
      real(8)    :: array (ilom:ihip,jlom:jhip,klom:khip)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,iglobal,nxx
      real(8)    :: x1,y1,z1,kk,mm
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      real(8) :: dum,rr(nxd,neqd),ii(nxd,neqd),pert(neqd)

c     Begin program

      mm = grid_params%params(1)
      kk = grid_params%params(2)

c     Read perturbation file

      if (my_rank == 0) write (*,*) 'Reading PPNCH eigenmodes...'

      open(unit=u_prtf,file=trim(prt_file),status='old')

      read(u_prtf,*) nxx

      if (nxx /= nxd) then
        call pstop('perturbEquilibrium_ppnch','Grid sizes do not agree')
      endif

      do i=1,nxd
        read(u_prtf,*) dum,rr(i,1:8),ii(i,1:8)
      enddo

      close (u_prtf)

      ii = -ii

c     Find perturbations

      do i = ilo,ihi
        do k = klom,khip
          do j = jlom,jhip
            call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

            iglobal = i + grid_params%ilo(igx) - 1

            pert = rr(iglobal,:)*cos(y1)+ii(iglobal,:)*sin(y1)

            pert(2) = x1*pert(2)
            if (coords == 'hel') pert(3) = pert(3)+kk*x1/mm*pert(4)
            pert(4) = x1*pert(4)

            pert(5) = x1*pert(5)
            if (coords == 'hel') pert(6) = pert(6)+kk*x1/mm*pert(7)
            pert(7) = x1*pert(7)

cc            array(i,j,k) = array(i,j,k) + pert(ieq)
cc            array(i,j,k) = array(i,j,k) + perturb*pert(ieq)
            array(i,j,k) = array(i,j,k) + 1d-3*pert(ieq)
          enddo
        enddo
      enddo

c     End program

      end subroutine perturbEquilibrium_ppnch_eig

c     perturbEquilibrium_ppnch_ext
c     #################################################################
      subroutine perturbEquilibrium_ppnch_ext
c     -----------------------------------------------------------------
c     Marco, 14/11/2014
c     Reads fields from and externally supplied file.
c     Eigenfunctions of the fields are read from the external file
c     In the external file it is written
c     nxd
c     nyd
c     nzd
c     fft_br
c     fft_bt
c     fft_bz
c     fft_vr
c     fft_vt
c     fft_vz
c     fft_tmp
cc     -----------------------------------------------------------------

      use B_tools

cc      use app_iosetup, ONLY:it_divcl

      implicit none

c     Call variables

c     Local variables

      integer :: i,j,k,ig,jg,kg,iglobal,nr,ff,nx,ny,nz
      real(8) :: x1,y1,z1,kk,mm,y00

      complex(8), allocatable,dimension(:,:,:)  :: fft_br
     .                                            ,fft_bt
     .                                            ,fft_bz
     .                                            ,fft_vr
     .                                            ,fft_vz
     .                                            ,fft_vt          
     .                                            ,fft_tmp

      real(8), allocatable :: rin(:),rdum(:,:),idum(:,:)
 
      real(8),allocatable,dimension(:,:,:,:) :: db_cnv,dv_cnv
      real(8),allocatable,dimension(:,:,:)   :: dprs

c     Begin program

c     Read external file

      if (my_rank == 0) then
        write (*,*)
        write (*,*) 'Reading external perturbation file '
     .             ,trim(prt_file)
      endif

      open(unit=u_prtf,file=trim(prt_file),form='unformatted')

      read(u_prtf) nr
      read(u_prtf) nh2
      read(u_prtf) nh3

      if (my_rank == 0) write(*,*) 'Perturbation mesh=',nr,nh2,nh3

      allocate(fft_br(0:nr+1,0:nh2-1,0:nh3-1)
     .        ,fft_bt(0:nr+1,0:nh2-1,0:nh3-1)
     .        ,fft_bz(0:nr+1,0:nh2-1,0:nh3-1)
     .        ,fft_vr(0:nr+1,0:nh2-1,0:nh3-1)
     .        ,fft_vt(0:nr+1,0:nh2-1,0:nh3-1)
     .        ,fft_vz(0:nr+1,0:nh2-1,0:nh3-1)
     .        ,fft_tmp(0:nr+1,0:nh2-1,0:nh3-1)
     .        ,rin(0:nr+1),rdum(nxd,7),idum(nxd,7))

      read(u_prtf) fft_br
      read(u_prtf) fft_bt
      read(u_prtf) fft_bz
      read(u_prtf) fft_vr
      read(u_prtf) fft_vt
      read(u_prtf) fft_vz
      read(u_prtf) fft_tmp

      close (u_prtf)

c     Reconstruct input mesh

      call findGlobalGrid(xmin,xmax,bcond(1),bcond(2),rin,g_pack(1),1)

c     Reconstruct field on target mesh

      nx = ihi-ilo+1
      ny = jhi-jlo+1
      nz = khi-klo+1

      allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3))
      allocate(dv_cnv(0:nx+1,0:ny+1,0:nz+1,3))
      allocate(dprs  (0:nx+1,0:ny+1,0:nz+1))

      db_cnv = 0d0
      dv_cnv = 0d0
      dprs   = 0d0

      select case(coords)
      case('hel')

        mm = grid_params%params(1)
        kk = grid_params%params(2)

        !We are NOT including equilibrium component (ff=0)
        do ff = 1, nh2/2-1 
          
c#c magnetic field
          call IntDriver1d (nr+2,rin,dreal(fft_br(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),rdum(:,1),3,0)  !Cubic splines
          call IntDriver1d (nr+2,rin,dimag(fft_br(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),idum(:,1),3,0)  !Cubic splines

          call IntDriver1d (nr+2,rin,dreal(fft_bt(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),rdum(:,2),3,0)  !Cubic splines
          call IntDriver1d (nr+2,rin,dimag(fft_bt(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),idum(:,2),3,0)  !Cubic splines

          call IntDriver1d (nr+2,rin,dreal(fft_bz(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),rdum(:,3),3,0)  !Cubic splines
          call IntDriver1d (nr+2,rin,dimag(fft_bz(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),idum(:,3),3,0)  !Cubic splines
c#c velocity field
          call IntDriver1d (nr+2,rin,dreal(fft_vr(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),rdum(:,4),3,0)  !Cubic splines
          call IntDriver1d (nr+2,rin,dimag(fft_vr(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),idum(:,4),3,0)  !Cubic splines

          call IntDriver1d (nr+2,rin,dreal(fft_vt(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),rdum(:,5),3,0)  !Cubic splines
          call IntDriver1d (nr+2,rin,dimag(fft_vt(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),idum(:,5),3,0)  !Cubic splines

          call IntDriver1d (nr+2,rin,dreal(fft_vz(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),rdum(:,6),3,0)  !Cubic splines
          call IntDriver1d (nr+2,rin,dimag(fft_vz(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),idum(:,6),3,0)  !Cubic splines
c#c temperature field
          call IntDriver1d (nr+2,rin,dreal(fft_tmp(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),rdum(:,7),3,0)  !Cubic splines
          call IntDriver1d (nr+2,rin,dimag(fft_tmp(:,ff,0))
     .                     ,nxd,grid_params%xg(1:nxd),idum(:,7),3,0)  !Cubic splines

          do k = klo, khi
            do j = jlo, jhi
              do i = ilo, ihi
                call getCurvilinearCoordinates(i,j,k,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)

                iglobal = i + grid_params%ilo(igx) - 1

                pert(IBX) = + rdum(iglobal,1) * cos(ff*y1)
     .                      - idum(iglobal,1) * sin(ff*y1)
                pert(IBY) = + rdum(iglobal,2) * cos(ff*y1)
     .                      - idum(iglobal,2) * sin(ff*y1)
                pert(IBZ) = + rdum(iglobal,3) * cos(ff*y1)
     .                      - idum(iglobal,3) * sin(ff*y1)
                pert(IVX) = + rdum(iglobal,4) * cos(ff*y1)
     .                      - idum(iglobal,4) * sin(ff*y1)
                pert(IVY) = + rdum(iglobal,5) * cos(ff*y1)
     .                      - idum(iglobal,5) * sin(ff*y1)
                pert(IVZ) = + rdum(iglobal,6) * cos(ff*y1)
     .                      - idum(iglobal,6) * sin(ff*y1)
                pert(ITMP)= + rdum(iglobal,7) * cos(ff*y1)
     .                      - idum(iglobal,7) * sin(ff*y1)

                db_cnv(i,j,k,1) = db_cnv(i,j,k,1) + pert(IBX)*x1
                db_cnv(i,j,k,2) = db_cnv(i,j,k,2) + pert(IBY)
     .                                            + kk/mm * pert(IBZ)*x1
                db_cnv(i,j,k,3) = db_cnv(i,j,k,3) + pert(IBZ)*x1
                dv_cnv(i,j,k,1) = dv_cnv(i,j,k,1) + pert(IVX)*x1

                if (alt__eom()) then
                  dv_cnv(i,j,k,2) = dv_cnv(i,j,k,2)
     .                 + (pert(IVY) + kk/mm*pert(IVZ)*x1)*x1
                else
                  dv_cnv(i,j,k,2) = dv_cnv(i,j,k,2) + pert(IVY)
     .                                              + kk/mm*pert(IVZ)*x1
                endif

                dv_cnv(i,j,k,3) = dv_cnv(i,j,k,3) + pert(IVZ)*x1

                dprs(i,j,k)     = dprs(i,j,k)     + pert(ITMP)*x1
              enddo
            enddo
          enddo

        enddo

      case default

        call pstop('perturbEquilibrium_ppnch'
     .            ,'Coordinate not implemented')

      end select

      deallocate(fft_br
     .          ,fft_bt
     .          ,fft_bz
     .          ,fft_vr
     .          ,fft_vt
     .          ,fft_vz
     .          ,fft_tmp
     .          ,rin,rdum,idum)

c     Divergence clean perturbation

      if (my_rank == 0) then
        write (*,*)
        write (*,*) 'Divergence cleaning perturbation...'
        if (ilevel > 0) write (*,*)
      endif

      call vlap_div_clean(ilevel,nx,ny,nz,igx,igy,igz
     .                   ,db_cnv,w_pre_bcs=.true.
     .                   ,gm_it=it_divcl,rtol=0.1*rtol)

      if (my_rank == 0.and.ilevel == 0) then
        write (*,*) 'Divergence cleaning converged in',it_divcl,' its'
      endif

c     Add equilibrium

      varr%array_var(IBX)%array = db_cnv(:,:,:,1)
     .                 +varr%array_var(IBX)%array
      varr%array_var(IBY)%array = db_cnv(:,:,:,2)
     .                 +varr%array_var(IBY)%array
      varr%array_var(IBZ)%array = db_cnv(:,:,:,3)
     .                 +varr%array_var(IBZ)%array
      varr%array_var(IVX)%array = dv_cnv(:,:,:,1)
     .                 +varr%array_var(IVX)%array
      varr%array_var(IVY)%array = dv_cnv(:,:,:,2)
     .                 +varr%array_var(IVY)%array
      varr%array_var(IVZ)%array = dv_cnv(:,:,:,3)
     .                 +varr%array_var(IVZ)%array
      varr%array_var(ITMP)%array = dprs
     .                 +varr%array_var(ITMP)%array

      deallocate(db_cnv,dv_cnv,dprs)

c     End program

      end subroutine perturbEquilibrium_ppnch_ext

ccc     perturbEquilibrium_gem
ccc     #################################################################
cc      subroutine perturbEquilibrium_gem(array,bcs,perturb,ieq)
cc
ccc     -----------------------------------------------------------------
ccc     Perturbs equilibrium quantity in array0 with a sinusoidal
ccc     perturbation of magnitude perturb, and introduces it in array.
ccc     -----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer    :: bcs(6),ieq
cc      real(8)    :: perturb
cc      real(8)    :: array (ilom:ihip,jlom:jhip,klom:khip)
cc
ccc     Local variables
cc
cc      integer    :: i,j,k,ig,jg,kg,ip,im,jp,jm
cc      real(8)    :: x1,y1,z1,kk,mm,lx,ly
cc      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 
cc
cc      real(8)    :: psi(ilom:ihip,jlom:jhip,klom:khip)
cc
ccc     Begin program
cc
cc      if (ieq == IBX .or. ieq == IBY) then
cc
cc        lx = xmax-xmin
cc        if (bcond(1) == FSYM .or. bcond(2) == FSYM) lx = 2*(xmax-xmin)
cc        ly = ymax-ymin
cc        if (bcond(3) == FSYM .or. bcond(4) == FSYM) ly = 2*(ymax-ymin)
cc
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)
cc
cc              psi(i,j,k)=-sin(  pi*(x1-xmin)/lx)
cc     .                   *cos(2*pi*(y1-ymin)/ly)
cc            enddo
cc          enddo
cc        enddo
cc
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              ip = min(i+1,ihip)
cc              im = max(i-1,ilom)
cc              jp = min(j+1,jhip)
cc              jm = max(j-1,jlom)
cc
cc              if (ieq == IBX) then
cc                array(i,j,k) = array(i,j,k)
cc     .                   - perturb*(psi(i,jp,k)-psi(i,jm,k))
cc     .                              /(grid_params%yy(jg+jp-j)
cc     .                               -grid_params%yy(jg+jm-j))
cc              else
cc                array(i,j,k) = array(i,j,k)
cc     .                   + perturb*(psi(ip,j,k)-psi(im,j,k))
cc     .                              /(grid_params%xx(ig+ip-i)
cc     .                               -grid_params%xx(ig+im-i))
cc              endif
cc            enddo
cc          enddo
cc        enddo
cc
cccc      elseif (ieq == IVX .or. ieq == IVY) then
cccc
cccc        do k = klom,khip
cccc          do j = jlom,jhip
cccc            do i = ilom,ihip
cccc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cccc     .                                    ,x1,y1,z1)
cccc
cccc              psi(i,j,k)=sin(2*pi*(x1-xmin)/(xmax-xmin))
cccc     .                  *sin(2*pi*(y1-ymin)/(ymax-ymin))
cccc            enddo
cccc          enddo
cccc        enddo
cccc        
cccc        do k = klom,khip
cccc          do j = jlom,jhip
cccc            do i = ilom,ihip
cccc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cccc
cccc              ip = min(i+1,ihip)
cccc              im = max(i-1,ilom)
cccc              jp = min(j+1,jhip)
cccc              jm = max(j-1,jlom)
cccc
cccc              if (ieq == IVX) then
cccc                array(i,j,k) = array(i,j,k)
cccc     .                   - perturb*(psi(i,jp,k)-psi(i,jm,k))
cccc     .                              /(grid_params%yy(jg+jp-j)
cccc     .                               -grid_params%yy(jg+jm-j))
cccc              else
cccc                array(i,j,k) = array(i,j,k)
cccc     .                   + perturb*(psi(ip,j,k)-psi(im,j,k))
cccc     .                              /(grid_params%xx(ig+ip-i)
cccc     .                               -grid_params%xx(ig+im-i))
cccc              endif
cccc            enddo
cccc          enddo
cccc        enddo
cc
cc      else
cc
cc        call perturbEquilibrium_def(array,bcs,perturb)
cc
cc      endif
cc
ccc     End program
cc
cc      end subroutine perturbEquilibrium_gem
cc
ccc     perturbEquilibrium_dtm
ccc     #################################################################
cc      subroutine perturbEquilibrium_dtm(array,bcs,perturb,ieq)
cc
ccc     -----------------------------------------------------------------
ccc     Perturbs equilibrium quantity in array0 with a sinusoidal
ccc     perturbation of magnitude perturb, and introduces it in array.
ccc     -----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer    :: bcs(6),ieq
cc      real(8)    :: perturb
cc      real(8)    :: array (ilom:ihip,jlom:jhip,klom:khip)
cc
ccc     Local variables
cc
cc      integer    :: i,j,k,ig,jg,kg,ip,im,jp,jm
cc      real(8)    :: x1,y1,z1,kk,mm
cc      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 
cc
cc      real(8)    :: psi(ilom:ihip,jlom:jhip,klom:khip)
cc
ccc     Begin program
cc
cc      if (ieq == IBX .or. ieq == IBY) then
cc
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,x1,y1,z1)
cc
cc              psi(i,j,k)=-cos(2*pi*(x1-xmin)/(xmax-xmin))
cc     .                   *cos(2*pi*(y1-ymin)/(ymax-ymin))
cc            enddo
cc          enddo
cc        enddo
cc        
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              ip = min(i+1,ihip)
cc              im = max(i-1,ilom)
cc              jp = min(j+1,jhip)
cc              jm = max(j-1,jlom)
cc
cc              if (ieq == IBX) then
cc                array(i,j,k) = array(i,j,k)
cc     .                   - perturb*(psi(i,jp,k)-psi(i,jm,k))
cc     .                              /(grid_params%yy(jg+jp-j)
cc     .                               -grid_params%yy(jg+jm-j))
cc              else
cc                array(i,j,k) = array(i,j,k)
cc     .                   + perturb*(psi(ip,j,k)-psi(im,j,k))
cc     .                              /(grid_params%xx(ig+ip-i)
cc     .                               -grid_params%xx(ig+im-i))
cc              endif
cc            enddo
cc          enddo
cc        enddo
cc
cc      else
cc
cc        call perturbEquilibrium_def(array,bcs,perturb)
cc
cc      endif
cc
ccc     End program
cc
cc      end subroutine perturbEquilibrium_dtm
cc
ccc     perturbEquilibrium_3dpar
ccc     #################################################################
cc      subroutine perturbEquilibrium_3dpar(array,bcs,perturb,ieq)
cc
ccc     -----------------------------------------------------------------
ccc     Perturbs equilibrium quantity in array0 with a sinusoidal
ccc     perturbation of magnitude perturb, and introduces it in array.
ccc     -----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer    :: bcs(6),ieq
cc      real(8)    :: perturb
cc      real(8)    :: array (ilom:ihip,jlom:jhip,klom:khip)
cc
ccc     Local variables
cc
cc      integer    :: i,j,k,ig,jg,kg,igx,igy,igz,ip,im,jp,jm
cc     .             ,ipw1,ipw2,ipw3
cc      real(8)    :: x1,y1,z1,kk,mm,nn,rs,iq0,RR,sigma,r0
cc      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 
cc
cc      real(8)    :: psi(ilom:ihip,jlom:jhip,klom:khip)
ccc
ccc     Begin program
cc
cc      igx = 1
cc      igy = 1
cc      igz = 1
cc
cc      if (ieq == IBX .or. ieq == IBY) then
cc
cc        mm = grid_params%params(1)
cc        kk = grid_params%params(2)
cc        RR = grid_params%params(3)
cc        iq0= grid_params%params(4)
cc
cc        nn = abs(kk*RR)
cc
cccc        !mm = 2, nn=1
cccc        mm = 2
cccc        nn = 1
cccc
cccc        kk = -nn/RR
cccc
cccc        sigma = 0.1
cccc
cccc        rs = dlambda*sqrt(q0*mm/nn-1)
cccc
cccc        r0 = rs-mm*sigma**2/rs
cccc
cccc        write (*,*) 'here',rs,r0,sigma,mm,nn,q0
cccc
cccc        do k = klom,khip
cccc          do j = jlom,jhip
cccc            do i = ilom,ihip
cccc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cccc     .                                    ,x1,y1,z1)
cccc
cccc              psi(i,j,k)=(x1/dlambda)**mm
cccc     .                  *exp(-.5*((x1-r0)/sigma)**2)
cccc     .                  *cos(mm*y1+kk*z1)
cccc            enddo
cccc          enddo
cccc        enddo
cccc
cccc        !mm = 3, nn=2
cccc        mm = 3
cccc        nn = 2
cccc
cccc        kk = -nn/RR
cccc
cccc        sigma = 0.1
cccc
cccc        rs = dlambda*sqrt(q0*mm/nn-1)
cccc
cccc        r0 = rs-mm*sigma**2/rs
cccc
cccc        do k = klom,khip
cccc          do j = jlom,jhip
cccc            do i = ilom,ihip
cccc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cccc     .                                    ,x1,y1,z1)
cccc
cccc              psi(i,j,k)=psi(i,j,k) + (x1/dlambda)**mm
cccc     .                  *exp(-.5*((x1-r0)/sigma)**2)
cccc     .                  *cos(mm*y1+kk*z1)
cccc            enddo
cccc          enddo
cccc        enddo
cc
cc        sigma = 0.1
cc
cc        psi = 0d0
cc
cc        do ipw3=nh3,nh3+sign(npw3-1,nh3),sign(1,nh3)
cc          do ipw2=max(floor(ipw3/iq0),nh2),nh2+sign(npw2-1,nh2)
cc     .           ,sign(1,nh2)
cc
cc            mm = ipw2
cc            nn = ipw3
cc
cc            kk = -nn/RR
cc
cccc            if (mm == 3d0 .and. nn == 1d0) cycle
cc            if (iq0*mm/nn < 1d0) cycle
cc
cc            rs = dlambda*sqrt(iq0*mm/nn-1)
cc
cc            if (rs > 1d0 .or. rs == 0d0) cycle
cc
cc            write (*,*) 'n,m=',nn,mm
cc
cc            r0 = rs-mm*sigma**2/rs
cc
cc            do k = klom,khip
cc              do j = jlom,jhip
cc                do i = ilom,ihip
cc                  call getCurvilinearCoordinates(i,j,k,igx,igy,igz
cc     .                                          ,ig,jg,kg,x1,y1,z1)
cc
cc                  psi(i,j,k)=psi(i,j,k) + (x1/dlambda)**mm
cc     .                      *exp(-.5*((x1-r0)/sigma)**2)*(1-x1)**mm
cc     .                      *cos(mm*y1+kk*z1)
cc                enddo
cc              enddo
cc            enddo
cc
cc          enddo
cc        enddo
cc
cc        !Add perturbation
cc        
cc        do k = klo,khi
cc          do j = jlo,jhi
cc            do i = ilo,ihi
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              ip = min(i+1,ihip)
cc              im = max(i-1,ilom)
cc              jp = min(j+1,jhip)
cc              jm = max(j-1,jlom)
cc
cc              if (ieq == IBX) then
cc                array(i,j,k) = array(i,j,k)
cc     .                   - perturb*(psi(i,jp,k)-psi(i,jm,k))
cc     .                              /(grid_params%yy(jg+jp-j)
cc     .                               -grid_params%yy(jg+jm-j))
cc              else
cc                array(i,j,k) = array(i,j,k)
cc     .                   + perturb*(psi(ip,j,k)-psi(im,j,k))
cc     .                              /(grid_params%xx(ig+ip-i)
cc     .                               -grid_params%xx(ig+im-i))
cc              endif
cc            enddo
cc          enddo
cc        enddo
cc
cc      else
cc
cc        call perturbEquilibrium_def(array,bcs,perturb)
cc
cc      endif
cc
ccc     End program
cc
cc      end subroutine perturbEquilibrium_3dpar

c     perturbEquilibrium_psi
c     #################################################################
      subroutine perturbEquilibrium_psi(array,bcs,perturb,ieq)

c     -----------------------------------------------------------------
c     Perturbs equilibrium quantity in array0 with a sinusoidal
c     perturbation of magnitude perturb, and introduces it in array.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: bcs(6),ieq
      real(8) :: perturb
      real(8) :: array (ilom:ihip,jlom:jhip,klom:khip)

c     Local variables

      integer :: i,j,k,ig,jg,kg,ip,im,jp,jm,bcsl(6),nx,ny,nz
      real(8) :: x1,y1,z1,kk,mm,xle,yle,zle
      real(8) :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      real(8) :: psi(ilom:ihip,jlom:jhip,klom:khip),grd(3)

c     Begin program

      if (ieq == IBX .or. ieq == IBY) then

        !Find vector potential
        psi = 0d0

        !Set psi BCs (Neumann if symmetry, Dirichlet otherwise)
        bcsl = bcond
        where(bcsl == DEF) bcsl = DIR

        call perturbEquilibrium_def(psi,bcsl,perturb,.true.)

        !Find delta_B
        nx = ihip-ilom-1
        ny = jhip-jlom-1
        nz = khip-klom-1

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              grd = grad(i,j,k,nx,ny,nz,igx,igy,igz,psi)

              if (ieq == IBX) then
                array(i,j,k) = array(i,j,k) - grd(2)
              else
                array(i,j,k) = array(i,j,k) + grd(1)
              endif
            enddo
          enddo
        enddo

      else

        call perturbEquilibrium_def(array,bcs,perturb,.true.)

      endif

c     End program

      end subroutine perturbEquilibrium_psi

c     factor
c     ####################################################################
      function factor(xmin,xmax,x,bcs,nh,odd) result(ff)

        implicit none

        real(8)    :: xmin,xmax,x,period,ff
        integer    :: bcs(2),nh

        logical    :: neumann(2),dirichlet(2),spoint(2),odd

c     Begin program

        spoint    = (abs(bcs)==SP)
        neumann   = (abs(bcs)==NEU).or.(bcs== SYM).or.(bcs== FSYM)
        dirichlet = (abs(bcs)==DIR).or.(bcs==-SYM).or.(bcs==-FSYM)
     .          .or.(abs(bcs)==EQU)

        if (nh == 0 .and. spoint(1)) then
          call pstop('factor','m for r^m perturbation is zero')
        endif

        period = pi
        if (odd) period = 2*pi

        if (bcs(1) == PER) then
          ff = cos(nh*2*pi*(x-xmin)/(xmax-xmin))
        elseif (random) then
          call random_number(ff)
        elseif (neumann(1) .and. neumann(2)) then
          ff = cos(period*(x-xmin)/(xmax-xmin))
        elseif (neumann(1) .and. dirichlet(2)) then
          if (.not.odd) then
            period = period/2.
          else
            period = 3*period/4.
          endif
          ff = cos(period*(x-xmin)/(xmax-xmin))
        elseif (dirichlet(1) .and. neumann(2)) then
          if (.not.odd) then
            period = period/2.
          else
            period = 3*period/4.
          endif
          ff = sin(period*(x-xmin)/(xmax-xmin))
        elseif (spoint(1) .and. dirichlet(2)) then
          ff = (sin(period*(x-xmin)/(xmax-xmin)))**(nh+2) !To satisfy regularity at r=0 (r^m)
     .         *sign(1d0,sin(period*(x-xmin)/(xmax-xmin)))**(nh+1)
        elseif (spoint(1) .and. neumann(2)) then
          if (.not.odd) then
            period = period/2.
            ff = (sin(period*(x-xmin)/(xmax-xmin)))**(nh+2) !To satisfy regularity at r=0 (r^m)
          else
            period = 3*period/4.
            ff = (sin(period*(x-xmin)/(xmax-xmin)))**(nh+2) !To satisfy regularity at r=0 (r^m)
     .        *sign(1d0,sin(period*(x-xmin)/(xmax-xmin)))
          endif
        else
          ff = sin(period*(x-xmin)/(xmax-xmin))
        endif

      end function factor

      end subroutine perturbEquilibrium

c     icfl
c     ##########################################################
      function icfl(cs2,ca2,k2,k2par,di) result(root)

      use math

c     ----------------------------------------------------------
c     Finds CFL frequency
c     ----------------------------------------------------------

      implicit none

c     Call variables

      real(8)   :: cs2,ca2,k2,k2par,di,root(3)

c     Local variables

      real(8)   :: a,b,c,d,ckpar
      complex(8):: r(3)

c     Begin program

      ckpar = ca2*k2par

c     Solve cubic dispersion relation for omega^2 (Checked, Luis 6/3/08)

      a =  1d0
      b = -ckpar*(1+k2*di**2)-k2*(ca2 + cs2)
      c = ckpar*(k2**2*cs2*di**2+2*k2*cs2+k2*ca2)
      d = -ckpar**2*cs2*k2

      r = solve_cubic(a,b,c,d)

ccc     Check roots
cc
cc      write (*,*) 'Coefficients=',a,b,c,d
cc
cc      r = a*r**3 + b*r**2 + c*r + d
cc
cc      write (*,*) 'Dispersion relation=',r
cc      stop

c     Find CFL constraint

      root = sqrt(abs(r))

c     End program

      end function icfl

c     perturbEquilibrium_ppnch_bc
c     #################################################################
      subroutine perturbEquilibrium_ppnch_bc(varr,br_pert,br_phase)
c     -----------------------------------------------------------------
c     Perturb normal B-field component ensuring div(B)=0 numerically
c     -----------------------------------------------------------------

      use grid

      use var_setup

      use timeStepping

      use app_iosetup

      use equilibrium

      use transport_params

      implicit none

c     Call variables

      type(var_array),pointer :: varr

      real(8) :: br_pert,br_phase

c     Local variables

      integer :: i,j,k,ig,jg,kg,ip,im,jp,jm,kp,km,mm,igx,igy,igz
      real(8) :: x1,y1,z1,kk,dxx,dyy,dzz,RR

      real(8) :: kbesel,beseli,dbeseli,dbeselip, dbeselim,br_pert_l

c     Begin program
      
      igx = 1
      igy = 1
      igz = 1

      if (br_pert /= 0d0) then
 
        !Luis, 11/17/14. We have now reserved nh2, nh3 for perturbations
        select case(coords)
        case('hel')
          mm = nh2
          kk = 0d0
        case('cyl')
          RR = grid_params%params(3)
          mm = nh2
          kk = nh3*1d0/RR
! Daniele, 19 Dec. 2014
! if k=0 the vacuum solution is b^1=r^m
! see Fitzpatrick PoP 1999
! implemented in the irrot_br_pert=.false. case
          if (kk == 0d0) irrot_br_pert=.false.
        end select

        if (irrot_br_pert) then  !Irrotational magnetic perturbation (lap(psi)=0)

          !#! Marco: rescaling of br_pert_bc to give it back it's original br(a) meaning
          kbesel = -mm*grid_params%params(2)
          dbeseli = 0.5*(BESSI(mm+1,kbesel)
     .                 + BESSI(mm-1,kbesel))
!          br_pert_l = 2.d0 * br_pert / (dbeseli)
! Daniele, 19 Dec. 2014
          br_pert_l = br_pert/(mm*dbeseli)

#if defined(vec_pot)

          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                call getCurvilinearCoordinates(i,j,k,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)

cc                kbesel = -mm*grid_params%params(2)
                beseli = BESSI(mm,kbesel*x1)
                dbeseli = 0.5*(BESSI(mm+1,kbesel*x1)
     .                       + BESSI(mm-1,kbesel*x1))

                varr%array_var(IAX)%array(i,j,k) =
     .               varr%array_var(IAX)%array(i,j,k)
     .               -(mm/kbesel)**2.*br_pert_l/x1*beseli
     .             *sin(mm*y1+kk*z1+br_phase)

                varr%array_var(IAY)%array(i,j,k) =
     .               varr%array_var(IAY)%array(i,j,k)
     .               -(mm/kbesel)*br_pert_l*x1*dbeseli
     .             *cos(mm*y1+kk*z1+br_phase)

                if (coords.eq.'hel') then

                varr%array_var(IAZ)%array(i,j,k) =
     .               varr%array_var(IAZ)%array(i,j,k)
     .               -br_pert_l*x1*dbeseli
     .             *cos(mm*y1+kk*z1+br_phase)
                
                endif

              enddo
            enddo
          enddo

          !Redefine equilibrium to force correct BCs on perturbed A
!          acov_0(:,:,:,1)=varr%array_var(IAX)%array
!          acov_0(:,:,:,2)=varr%array_var(IAY)%array
!          acov_0(:,:,:,3)=varr%array_var(IAZ)%array
          u_0%array_var(IAX)%array=varr%array_var(IAX)%array          
          u_0%array_var(IAY)%array=varr%array_var(IAY)%array
          u_0%array_var(IAZ)%array=varr%array_var(IAZ)%array         
#else

          do k = klo,khi
            do j = jlo,jhi
              do i = ilo,ihi
                call getCurvilinearCoordinates(i,j,k,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)

cc                kbesel = -mm*grid_params%params(2)
                beseli = BESSI(mm,kbesel*x1)
                dbeseli = 0.5*(BESSI(mm+1,kbesel*x1)
     .                       + BESSI(mm-1,kbesel*x1))

                dxx = 2*grid_params%dxh(ig)
                dyy = 2*grid_params%dyh(jg)
                dzz = 2*grid_params%dzh(kg)

                ip = ig+1
                im = ig-1
                jp = jg+1
                jm = jg-1
                kp = kg+1
                km = kg-1

                dbeselip =0.5*(BESSI(mm+1,kbesel*grid_params%xx(ip))
     .                       + BESSI(mm-1,kbesel*grid_params%xx(ip)))
                dbeselim =0.5*(BESSI(mm+1,kbesel*grid_params%xx(im))
     .                       + BESSI(mm-1,kbesel*grid_params%xx(im)))

                varr%array_var(IBX)%array(i,j,k) =
     .               varr%array_var(IBX)%array(i,j,k)
     .               -br_pert_l*x1*dbeseli
     .             *(cos(mm*grid_params%yy(jp)+kk*z1+br_phase)
     .              -cos(mm*grid_params%yy(jm)+kk*z1+br_phase))
     .             /dyy !Derivative done numerically to preserve div(B)=0

                varr%array_var(IBY)%array(i,j,k) =
     .               varr%array_var(IBY)%array(i,j,k)
     .               +br_pert_l*
     .               (grid_params%xx(ip)*dbeselip
     .               -grid_params%xx(im)*dbeselim)
     .               /dxx !Derivative done numerically to preserve div(B)=0
     .               *cos(mm*y1+kk*z1+br_phase)

                varr%array_var(IBZ)%array(i,j,k) =
     .               varr%array_var(IBZ)%array(i,j,k)
     .               -mm*br_pert_l*x1*beseli
     .               *cos(mm*y1+kk*z1+br_phase)

              enddo
            enddo
          enddo

#endif

        else   !Standard magnetic perturbation from psi=r^(m+2)*cos(m*th + kk*z)

          do k = klo,khi
            do j = jlo,jhi
              do i = ilo,ihi
                call getCurvilinearCoordinates(i,j,k,igx,igy,igz
     .                                        ,ig,jg,kg,x1,y1,z1)

                ip = ig+1
                im = ig-1
                jp = jg+1
                jm = jg-1
                kp = kg+1
                km = kg-1

                dxx = 2*grid_params%dxh(ig)
                dyy = 2*grid_params%dyh(jg)
                dzz = 2*grid_params%dzh(kg)

! Daniele, 19 Dec. 2014
! b^1=r^m instead of r^(m+2)
! to be consistent with the vacuum solution in the kk=0 case, 
! see Fitzpatrick PoP 1999

                varr%array_var(IBX)%array(i,j,k) =
     .               varr%array_var(IBX)%array(i,j,k)
!     .               +x1**(mm+2.)*br_pert
     .               +x1**mm*br_pert
     .             *(cos(mm*grid_params%yy(jp)+kk*z1+br_phase)
     .              -cos(mm*grid_params%yy(jm)+kk*z1+br_phase))
     .             /dyy/mm !Derivative done numerically to preserve div(B)=0

                varr%array_var(IBY)%array(i,j,k) =
     .               varr%array_var(IBY)%array(i,j,k)
     .               -br_pert*
!     .               (grid_params%xx(ip)**(mm+2.)
!     .               -grid_params%xx(im)**(mm+2.))
     .               (grid_params%xx(ip)**mm
     .               -grid_params%xx(im)**mm)
     .               /dxx !Derivative done numerically to preserve div(B)=0
     .               *cos(mm*y1+kk*z1+br_phase)/mm

              enddo
            enddo
          enddo

        endif

      endif

c     End program

      end subroutine perturbEquilibrium_ppnch_bc
