c     perturbEquilibrium
c     #################################################################
      subroutine perturbEquilibrium(array,bcs,ieq)

c     -----------------------------------------------------------------
c     Perturbs equilibrium quantity in array0 with a sinusoidal
c     perturbation of magnitude perturb, and introduces it in array.
c     -----------------------------------------------------------------

      use grid

      use variable_setup

      use timeStepping

      use iosetup

      use equilibrium

      use transport_params

      use nlfunction_setup, ONLY:solve_prs

      implicit none

c     Call variables

      integer    :: bcs(6),ieq
      real(8)    :: array (ilom:ihip,jlom:jhip,klom:khip)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,igx,igy,igz
      real(8)    :: x1,y1,z1,jac,kx,ky,kk,eps
     .             ,vplus,vminus,vshear,cs2,kdotB,cth
     .             ,vxx,vyy,vzz,ivxx,ivyy,ivzz,xi,ixi,yy,iyy
     .             ,bxx,byy,bzz,ibxx,ibyy,ibzz,axx,ayy,azz
     .             ,w_fast,w_slow,w_shr,freq,roots(3),X,iX
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      logical    :: check

      integer    :: ierr
      character(50) :: command

      save  :: w_fast,w_slow,w_shr,freq

      INTERFACE
      function icfl(cs2,ca2,k2,k2par,di) result(root)
      real(8)   :: cs2,ca2,k2,k2par,di,root(3)
      end function icfl
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
          call perturbEquilibrium_kaitm(array,bcs,pert(ieq),ieq)
        else
          call perturbEquilibrium_def  (array,bcs,pert(ieq),ieq)
        endif

      case ('ppnch','ppn3d','ppnsl','ppnst','p3nsl')

        command = 'test -f '//trim(prt_file)
        ierr = shellcmd(trim(command))

        if (ierr == 0) then
          call perturbEquilibrium_ppnch(array,bcs,pert(ieq),ieq)
        else
          call perturbEquilibrium_def  (array,bcs,pert(ieq),ieq)
        endif

      case ('msw')

        kx = 2*pi*nh1/(xmax-xmin)
        ky = 2*pi*nh2/(ymax-ymin)

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)
cc              jac = gmetric%grid(igx)%jac(i,j,k)
cc              if (ieq == IRHO) jac = 1d0
              jac = 1d0
              if (ieq == IRHO) then
                array(i,j,k) = array(i,j,k)
     .                       + 1d-3*cos(kx*(x1-xmin)+ky*(y1-ymin))
#if defined(vec_pot) 
              elseif (ieq == IAY) then
                array(i,j,k) = array(i,j,k)
     .                     + jac*1d-3/kx
     .                          *sin(kx*(x1-xmin)+ky*(y1-ymin))
#else
              elseif (ieq == IBZ) then
                array(i,j,k) = array(i,j,k)
     .                     + jac*1d-3*cos(kx*(x1-xmin)+ky*(y1-ymin))
#endif
              endif
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

              if (ieq == IRHO) then
                array(i,j,k) = array(i,j,k) - xi/freq
              elseif (ieq == ITMP) then
                if (solve_prs) then
                  array(i,j,k) = array(i,j,k) - cs2*xi/(a_p*freq) !For e pressure
                else
                  array(i,j,k) = array(i,j,k)
     .                         -(cs2/a_p-array(i,j,k))*xi/freq    !For temperature
                endif
#if !defined(vec_pot)
              elseif (ieq == IBX) then
                array(i,j,k) = array(i,j,k) + bxx
              elseif (ieq == IBY) then
                array(i,j,k) = array(i,j,k) + byy
              elseif (ieq == IBZ) then
                array(i,j,k) = array(i,j,k) + bzz
#else
              elseif (ieq == IAX) then
                array(i,j,k) = array(i,j,k) + axx
              elseif (ieq == IAY) then
                array(i,j,k) = array(i,j,k) + ayy
              elseif (ieq == IAZ) then
                array(i,j,k) = array(i,j,k) + azz
#endif
              elseif (ieq == IVX) then
                array(i,j,k) = array(i,j,k) + vxx
              elseif (ieq == IVY) then
                array(i,j,k) = array(i,j,k) + vyy
              elseif (ieq == IVZ) then
                array(i,j,k) = array(i,j,k) + vzz
              endif
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
#if !defined(vec_pot) 
              if (ieq == IBY) then
                array(i,j,k) = array(i,j,k)
     .                       + 1d-3*cos(kx*(x1-xmin)+ky*(y1-ymin))
c$$$              elseif (ieq == IBZ) then   !This sets a traveling wave
c$$$                array(i,j,k) = array(i,j,k)
c$$$     .                       + 1d-3*sin(nh1*2*pi*(x1-xmin)/(xmax-xmin))
#else
              if (ieq == IAZ) then
                array(i,j,k) = array(i,j,k)
     .                       - 1d-3/kx*sin(kx*(x1-xmin)+ky*(y1-ymin))
cc     .                       - 1d-3*(xmax-xmin)/2./pi/nh1
cc     $                             *sin(nh1*2*pi*(x1-xmin)/(xmax-xmin)
cc     .                                 +nh2*2*pi*(y1-ymin)/(ymax-ymin))
#endif
              elseif (ieq == IVZ) then
                array(i,j,k) = array(i,j,k)
     .                   + 1d-3/kx/di*sin(kx*(x1-xmin)+ky*(y1-ymin))
              endif
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
              if (ieq == IVY) then
                array(i,j,k) = array(i,j,k)
     .                       - 1d-7*kx/ky*(1.-(di*kk)**2)
     .                             *cos(kx*(x1-xmin)+ky*(y1-ymin))
              elseif (ieq == IVX) then
                array(i,j,k) = array(i,j,k)
     .                       + 1d-7*cos(kx*(x1-xmin)+ky*(y1-ymin))
              endif
            enddo
          enddo
        enddo

      case ('gem','tmcar','tmsin','kai1')

        call perturbEquilibrium_gem(array,bcs,pert(ieq),ieq)

      case ('dtm')

        call perturbEquilibrium_dtm(array,bcs,pert(ieq),ieq)

      case ('ic','ic-hk')

        call perturbEquilibrium_ic(array,bcs,pert(ieq),ieq)

cc      case ('3dpar')
cc
cc        call perturbEquilibrium_3dpar(array,bcs,pert(ieq),ieq)

      case default

        call perturbEquilibrium_def(array,bcs,pert(ieq),ieq)

      end select

c     End program

      contains

c     perturbEquilibrium_def
c     #################################################################
      subroutine perturbEquilibrium_def(array,bcs,perturb,ieq)

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

      integer    :: i,j,k,ig,jg,kg,iglobal
     .             ,ipw1,ipw2,ipw3,mode
      real(8)    :: x1,y1,z1
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

c     Begin program

      do ipw1=0,npw1-1
        do ipw2=0,npw2-1
          do ipw3=0,npw3-1

            mode = nh1+sign(ipw1,nh1)
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,jlo,klo,igx,igy,igz
     .                                      ,ig,jg,kg,x1,y1,z1)
              fx(i) = factor(xmin,xmax,x1,bcs(1:2),mode,odd(1))
            enddo

            mode = nh2+sign(ipw2,nh2)
            do j = jlom,jhip
              call getCurvilinearCoordinates(ilo,j,klo,igx,igy,igz
     .                                      ,ig,jg,kg,x1,y1,z1)
              fy(j) = factor(ymin,ymax,y1,bcs(3:4),mode,odd(2))
            enddo

            mode = nh3+sign(ipw3,nh3)
            do k = klom,khip
              call getCurvilinearCoordinates(ilo,jlo,k,igx,igy,igz
     .                                      ,ig,jg,kg,x1,y1,z1)
              fz(k) = factor(zmin,zmax,z1,bcs(5:6),mode,odd(3))
            enddo

            do k = klom,khip
              do j = jlom,jhip
                do i = ilom,ihip
                  array(i,j,k) = array(i,j,k)+perturb*fx(i)*fy(j)*fz(k)
                enddo
              enddo
            enddo

          enddo
        enddo
      enddo

c     End program

      end subroutine perturbEquilibrium_def

c     perturbEquilibrium_kaitm
c     #################################################################
      subroutine perturbEquilibrium_kaitm(array,bcs,perturb,ieq)

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

            pert(2) = x1*pert(2)
            if (equil == '2dkai') pert(3) = pert(3)+kk*x1/mm*pert(4)
            pert(4) = x1*pert(4)

            pert(5) = x1*pert(5)
            if (equil == '2dkai') pert(6) = pert(6)+kk*x1/mm*pert(7)
            pert(7) = x1*pert(7)

            array(i,j,k) = array(i,j,k) + pert(ieq)
          enddo
        enddo
      enddo

      close (u_prtf)

c     End program

      end subroutine perturbEquilibrium_kaitm

c     perturbEquilibrium_ppnch
c     #################################################################
      subroutine perturbEquilibrium_ppnch(array,bcs,perturb,ieq)

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

      end subroutine perturbEquilibrium_ppnch

c     perturbEquilibrium_gem
c     #################################################################
      subroutine perturbEquilibrium_gem(array,bcs,perturb,ieq)

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

      integer    :: i,j,k,ig,jg,kg,ip,im,jp,jm
      real(8)    :: x1,y1,z1,kk,mm,lx,ly
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      real(8)    :: psi(ilom:ihip,jlom:jhip,klom:khip)

c     Begin program

      if (ieq == IBX .or. ieq == IBY) then

        lx = xmax-xmin
        if (bcond(1) == FSYM .or. bcond(2) == FSYM) lx = 2*(xmax-xmin)
        ly = ymax-ymin
        if (bcond(3) == FSYM .or. bcond(4) == FSYM) ly = 2*(ymax-ymin)

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              psi(i,j,k)=-sin(  pi*(x1-xmin)/lx)
     .                   *cos(2*pi*(y1-ymin)/ly)
            enddo
          enddo
        enddo

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              ip = min(i+1,ihip)
              im = max(i-1,ilom)
              jp = min(j+1,jhip)
              jm = max(j-1,jlom)

              if (ieq == IBX) then
                array(i,j,k) = array(i,j,k)
     .                   - perturb*(psi(i,jp,k)-psi(i,jm,k))
     .                              /(grid_params%yy(jg+jp-j)
     .                               -grid_params%yy(jg+jm-j))
              else
                array(i,j,k) = array(i,j,k)
     .                   + perturb*(psi(ip,j,k)-psi(im,j,k))
     .                              /(grid_params%xx(ig+ip-i)
     .                               -grid_params%xx(ig+im-i))
              endif
            enddo
          enddo
        enddo

cc      elseif (ieq == IVX .or. ieq == IVY) then
cc
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,x1,y1,z1)
cc
cc              psi(i,j,k)=sin(2*pi*(x1-xmin)/(xmax-xmin))
cc     .                  *sin(2*pi*(y1-ymin)/(ymax-ymin))
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
cc              if (ieq == IVX) then
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

      else

        call perturbEquilibrium_def(array,bcs,perturb,ieq)

      endif

c     End program

      end subroutine perturbEquilibrium_gem

c     perturbEquilibrium_dtm
c     #################################################################
      subroutine perturbEquilibrium_dtm(array,bcs,perturb,ieq)

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

      integer    :: i,j,k,ig,jg,kg,ip,im,jp,jm
      real(8)    :: x1,y1,z1,kk,mm
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      real(8)    :: psi(ilom:ihip,jlom:jhip,klom:khip)

c     Begin program

      if (ieq == IBX .or. ieq == IBY) then

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              psi(i,j,k)=-cos(2*pi*(x1-xmin)/(xmax-xmin))
     .                   *cos(2*pi*(y1-ymin)/(ymax-ymin))
            enddo
          enddo
        enddo
        
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              ip = min(i+1,ihip)
              im = max(i-1,ilom)
              jp = min(j+1,jhip)
              jm = max(j-1,jlom)

              if (ieq == IBX) then
                array(i,j,k) = array(i,j,k)
     .                   - perturb*(psi(i,jp,k)-psi(i,jm,k))
     .                              /(grid_params%yy(jg+jp-j)
     .                               -grid_params%yy(jg+jm-j))
              else
                array(i,j,k) = array(i,j,k)
     .                   + perturb*(psi(ip,j,k)-psi(im,j,k))
     .                              /(grid_params%xx(ig+ip-i)
     .                               -grid_params%xx(ig+im-i))
              endif
            enddo
          enddo
        enddo

      else

        call perturbEquilibrium_def(array,bcs,perturb,ieq)

      endif

c     End program

      end subroutine perturbEquilibrium_dtm

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
cc
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
cc        call perturbEquilibrium_def(array,bcs,perturb,ieq)
cc
cc      endif
cc
ccc     End program
cc
cc      end subroutine perturbEquilibrium_3dpar

c     perturbEquilibrium_ic
c     #################################################################
      subroutine perturbEquilibrium_ic(array,bcs,perturb,ieq)

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

      integer    :: i,j,k,ig,jg,kg,ip,im,jp,jm
      real(8)    :: x1,y1,z1,kk,mm
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      real(8)    :: psi(ilom:ihip,jlom:jhip,klom:khip)

c     Begin program

      if (ieq == IBX .or. ieq == IBY) then

cc        if (equil == 'ic') then
cc          do k = klom,khip
cc            do j = jlom,jhip
cc              do i = ilom,ihip
cc                call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)
cc
cc                psi(i,j,k)=sin(  pi*(x1-xmin)/(xmax-xmin))
cc     .                    *cos(2*pi*(y1-ymin)/(ymax-ymin))
cc     .                    *cos(2*pi*(z1-zmin)/(zmax-zmin))
cc              enddo
cc            enddo
cc          enddo
cc        elseif (equil == 'ic-hk') then
          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

                psi(i,j,k)= cos(0.5*pi*x1)*cos(pi*y1)*cos(pi*z1)
              enddo
            enddo
          enddo
cc        endif

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              ip = min(i+1,ihip)
              im = max(i-1,ilom)
              jp = min(j+1,jhip)
              jm = max(j-1,jlom)

              if (ieq == IBX) then
                array(i,j,k) = array(i,j,k)
     .                   - perturb*(psi(i,jp,k)-psi(i,jm,k))
     .                              /(grid_params%yy(jg+jp-j)
     .                               -grid_params%yy(jg+jm-j))
              else
                array(i,j,k) = array(i,j,k)
     .                   + perturb*(psi(ip,j,k)-psi(im,j,k))
     .                              /(grid_params%xx(ig+ip-i)
     .                               -grid_params%xx(ig+im-i))
              endif
            enddo
          enddo
        enddo

      else

        call perturbEquilibrium_def(array,bcs,perturb,ieq)

      endif

c     End program

      end subroutine perturbEquilibrium_ic

c     factor
c     ####################################################################
      function factor(xmin,xmax,x,bcs,nh,odd) result(ff)

        implicit none

        real(8)    :: xmin,xmax,x,period,ff
        integer    :: bcs(2),nh

        logical    :: neumann(2),dirichlet(2),spoint(2),odd

c     Begin program

        spoint    = (abs(bcs)==SP)
        neumann   = (abs(bcs)==NEU).or.(bcs== SYM)
        dirichlet = (abs(bcs)==DIR).or.(bcs==-SYM).or.(abs(bcs)==EQU)

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
