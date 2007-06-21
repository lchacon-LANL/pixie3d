c     perturbEquilibrium
c     #################################################################
      subroutine perturbEquilibrium(array,bcs,perturb,ieq)

c     -----------------------------------------------------------------
c     Perturbs equilibrium quantity in array0 with a sinusoidal
c     perturbation of magnitude perturb, and introduces it in array.
c     -----------------------------------------------------------------

      use icond

      use grid

      use variable_setup

      use timeStepping

      use constants

      use iosetup

      use equilibrium

      implicit none

c     Call variables

      integer    :: bcs(6),ieq
      real(8)    :: perturb
      real(8)    :: array (ilom:ihip,jlom:jhip,klom:khip)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,igx,igy,igz
      real(8)    :: x1,y1,z1,jac
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      integer    :: system,ierr
      character(50) :: command

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      select case(equil)
      case ('kaitm','kai3d')

        command = 'ls '//trim(prt_file)//' >& /dev/null'
        ierr = system(trim(command))

        if (ierr == 0) then
          call perturbEquilibrium_kaitm(array,bcs,perturb,ieq)
        else
          call perturbEquilibrium_def  (array,bcs,perturb,ieq)
        endif

      case ('ppnch','ppn3d','ppnsl','ppnst')

        command = 'ls '//trim(prt_file)//' >& /dev/null'
        ierr = system(trim(command))

        if (ierr == 0) then
          call perturbEquilibrium_ppnch(array,bcs,perturb,ieq)
        else
          call perturbEquilibrium_def  (array,bcs,perturb,ieq)
        endif

      case ('msw')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)
              jac = gmetric%grid(igx)%jac(i,j,k)
              if (ieq == 1) jac = 1d0
              array(i,j,k) = array(i,j,k)
     .                     + jac*perturb*cos(2*pi*(x1-xmin)/(xmax-xmin)
     .                                      +2*pi*(y1-ymin)/(ymax-ymin))
            enddo
          enddo
        enddo

      case ('gem')

        call perturbEquilibrium_gem(array,bcs,perturb,ieq)

      case default

        call perturbEquilibrium_def(array,bcs,perturb,ieq)

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

      integer    :: i,j,k,ig,jg,kg,igx,igy,igz,iglobal
      real(8)    :: x1,y1,z1,kk,mm
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      do i = ilom,ihip
        call getCurvilinearCoordinates(i,jlo,klo,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fx(i) = factor(xmin,xmax,x1,bcs(1:2),nh1)
      enddo

      do j = jlom,jhip
        call getCurvilinearCoordinates(ilo,j,klo,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fy(j) = factor(ymin,ymax,y1,bcs(3:4),nh2)
      enddo

      do k = klom,khip
        call getCurvilinearCoordinates(ilo,jlo,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fz(k) = factor(zmin,zmax,z1,bcs(5:6),nh3)
      enddo

      do k = klom,khip
        do j = jlom,jhip
          do i = ilom,ihip
            array(i,j,k) = array(i,j,k) + perturb*fx(i)*fy(j)*fz(k)
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

      integer    :: i,j,k,ig,jg,kg,igx,igy,igz,iglobal
      real(8)    :: x1,y1,z1,kk,mm
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      real(8) :: dum,rr(nxd,neqd),ii(nxd,neqd),pert(neqd)

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      mm = grid_params%params(1)
      kk = grid_params%params(2)

      write (*,*) 'Reading KAI TM perturbations'

      open(unit=111,file=trim(prt_file),status='old')

      do i=1,nxd
        read(111,*) dum,rr(i,1),rr(i,8),rr(i,2:7)
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
            if (equil == 'kaitm') pert(3) = pert(3)+kk*x1/mm*pert(4)
            pert(4) = x1*pert(4)

            pert(5) = x1*pert(5)
            if (equil == 'kaitm') pert(6) = pert(6)+kk*x1/mm*pert(7)
            pert(7) = x1*pert(7)

            array(i,j,k) = array(i,j,k) + pert(ieq)
          enddo
        enddo
      enddo

      close (111)

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

      integer    :: i,j,k,ig,jg,kg,igx,igy,igz,iglobal,nxx
      real(8)    :: x1,y1,z1,kk,mm
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      real(8) :: dum,rr(nxd,neqd),ii(nxd,neqd),pert(neqd)

c     Begin program

      if (ieq == IBX .or. ieq == IBY .or. ieq == IBZ) return

      igx = 1
      igy = 1
      igz = 1

      mm = grid_params%params(1)
      kk = grid_params%params(2)

c     Read perturbation file

      write (*,*) 'Reading PPNCH eigenmodes...'

      open(unit=111,file=trim(prt_file),status='old')

      read(111,*) nxx

      if (nxx /= nxd) then
        call pstop('perturbEquilibrium_ppnch','Grid sizes do not agree')
      endif

      do i=1,nxd
        read(111,*) dum,rr(i,1:8),ii(i,1:8)
      enddo

      close (111)

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

      integer    :: i,j,k,ig,jg,kg,igx,igy,igz,ip,im,jp,jm
      real(8)    :: x1,y1,z1,kk,mm
      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

      real(8)    :: psi(ilom:ihip,jlom:jhip,klom:khip)

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      if (ieq == IBX .or. ieq == IBY) then

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              psi(i,j,k)=sin(pi*x1/(xmax-xmin))*sin(2*pi*y1/(ymax-ymin))
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

      end subroutine perturbEquilibrium_gem

c     factor
c     ####################################################################
      function factor(xmin,xmax,x,bcs,nh) result(ff)

        implicit none

        real(8)    :: xmin,xmax,x,period,ff
        integer    :: bcs(2),nh

        logical    :: neumann(2),dirichlet(2),spoint(2)

        spoint    = (bcs == SP)
        neumann   = (abs(bcs) == NEU) .or. (bcs == SYM)
        dirichlet = (bcs == DIR) .or. (bcs ==-SYM) .or. (bcs == EQU)

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
