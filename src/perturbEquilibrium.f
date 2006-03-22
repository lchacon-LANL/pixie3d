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

      integer(4) :: bcs(6),ieq
      real(8)    :: perturb
      real(8)    :: array (0:nxdp,0:nydp,0:nzdp)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,igx,igy,igz
      real(8)    :: x1,y1,z1,jac
      real(8)    :: fx(0:nxdp),fy(0:nydp),fz(0:nzdp) 

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      select case(equil)
      case ('kaitm','kai3d')

        call perturbEquilibrium_kaitm(array,bcs,perturb,ieq)

      case ('msw')

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd
              array(i,j,k) = array(i,j,k) + perturb*fx(i)*fy(j)*fz(k)
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

        do i = 1,nxd
          call getCurvilinearCoordinates(i,1,1,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)
          fx(i) = factor(xmin,xmax,x1,bcs(1:2),nh1)
        enddo

        do j = 1,nyd
          call getCurvilinearCoordinates(1,j,1,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)
          fy(j) = factor(ymin,ymax,y1,bcs(3:4),nh2)
        enddo

        do k = 1,nzd
          call getCurvilinearCoordinates(1,1,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)
          fz(k) = factor(zmin,zmax,z1,bcs(5:6),nh3)
        enddo

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd
              array(i,j,k) = array(i,j,k) + perturb*fx(i)*fy(j)*fz(k)
            enddo
          enddo
        enddo

      end select

c     End program

      contains

c     perturbEquilibrium_kaitm
c     #################################################################
      subroutine perturbEquilibrium_kaitm(array,bcs,perturb,ieq)

c     -----------------------------------------------------------------
c     Perturbs equilibrium quantity in array0 with a sinusoidal
c     perturbation of magnitude perturb, and introduces it in array.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: bcs(6),ieq
      real(8)    :: perturb
      real(8)    :: array (0:nxdp,0:nydp,0:nzdp)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,igx,igy,igz
      real(8)    :: x1,y1,z1,kk,mm
      real(8)    :: fx(0:nxdp),fy(0:nydp),fz(0:nzdp) 

      real(8) :: dum,rr(neqd),ii(neqd),pert(neqd)

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      mm = grid_params%params(1)
      kk = grid_params%params(2)

      write (*,*) 'Reading KAI TM perturbations'

      open(unit=111,file='M64_000020.asc',status='old')

      k = 1
      do i=1,nxd
        read(111,*) dum,rr(1),rr(8),rr(2:7),ii(1),ii(8),ii(2:7)
cc        write(*,*) dum,rr(1),rr(8),rr(2:7),ii(1),ii(8),ii(2:7)
        do j=1,nyd
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)

          pert = rr*cos(y1)+ii*sin(y1)

          pert(2) = x1*pert(2)
          pert(3) = pert(3)+kk*x1/mm*pert(4)
          pert(4) = x1*pert(4)

          pert(5) = x1*pert(5)
          pert(6) = pert(6)+kk*x1/mm*pert(7)
          pert(7) = x1*pert(7)

          array(i,j,k) = array(i,j,k) + pert(ieq)
        enddo
cc        write (*,*) x1
      enddo

      close (111)
cc      stop

c     End program

      end subroutine perturbEquilibrium_kaitm

c     perturbEquilibrium_gem
c     #################################################################
      subroutine perturbEquilibrium_gem(array,bcs,perturb,ieq)

c     -----------------------------------------------------------------
c     Perturbs equilibrium quantity in array0 with a sinusoidal
c     perturbation of magnitude perturb, and introduces it in array.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: bcs(6),ieq
      real(8)    :: perturb
      real(8)    :: array (0:nxdp,0:nydp,0:nzdp)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,igx,igy,igz
      real(8)    :: x1,y1,z1,kk,mm
      real(8)    :: fx(0:nxdp),fy(0:nydp),fz(0:nzdp) 

      real(8)    :: dum,rr(neqd),ii(neqd),pert(neqd)

      real(8)    :: psi(0:nxdp,0:nydp,0:nzdp)

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      if (ieq == 5 .or. ieq == 6) then
        k = 1

cc        pi = acos(-1d0)
        do j=0,nydp
          do i=0,nxdp
            call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

            psi(i,j,k)=sin(pi*x1/(xmax-xmin))*sin(2*pi*y1/(ymax-ymin))
          enddo
        enddo

        do j=1,nyd
          do i=1,nxd
            call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

            if (ieq == 5) then
              array(i,j,k) = array(i,j,k)
     .                 - perturb*(psi(i,j+1,k)-psi(i,j-1,k))
     .                            /2./grid_params%dyh(jg)
            else
              array(i,j,k) = array(i,j,k)
     .                 + perturb*(psi(i+1,j,k)-psi(i-1,j,k))
     .                            /2./grid_params%dxh(ig)
            endif
          enddo
        enddo
      endif

c     End program

      end subroutine perturbEquilibrium_gem

c     factor
c     ####################################################################
      function factor(xmin,xmax,x,bcs,nh) result(ff)

        implicit none

        real(8)    :: xmin,xmax,x,period,ff
        integer(4) :: bcs(2),nh

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
