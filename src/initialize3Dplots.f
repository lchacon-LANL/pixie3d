c defineGraphics
c####################################################################
      subroutine defineGraphics

c--------------------------------------------------------------------
c     Set graphics files and dumping intervals
c--------------------------------------------------------------------

      use parameters

      use variables

      use graphics

      use nlfunction_setup

      implicit none

c Call variables

c Local variables

      integer  :: ieq

c Begin program

      call defineGraphicsIO

c Define graphics group #1: Contravariant variables

      graph(1)%cartesian=.false.
      graph(1)%descr='Cnv variables'

      do ieq = 1,neqd
        graph(1)%array_graph(ieq)%array => u_np%array_var(ieq)%array
        graph(1)%array_graph(ieq)%descr =  u_np%array_var(ieq)%descr
      enddo

      graph(1)%array_graph(neqd+1)%array => jx
      graph(1)%array_graph(neqd+1)%descr = 'J^1'

      graph(1)%array_graph(neqd+2)%array => jy
      graph(1)%array_graph(neqd+2)%descr = 'J^2'

      graph(1)%array_graph(neqd+3)%array => jz
      graph(1)%array_graph(neqd+3)%descr = 'J^3'

      graph(1)%array_graph(neqd+4)%array => vx
      graph(1)%array_graph(neqd+4)%descr = 'V^1'

      graph(1)%array_graph(neqd+5)%array => vy
      graph(1)%array_graph(neqd+5)%descr = 'V^2'

      graph(1)%array_graph(neqd+6)%array => vz
      graph(1)%array_graph(neqd+6)%descr = 'V^3'

      graph(1)%array_graph(neqd+7:ngraph)%descr = ''

c Define graphics group #2: Cov variables

      graph(2)%cartesian=.false.
      graph(2)%descr='Cov variables'

      graph(2)%array_graph(IRHO)%array => u_np%array_var(IRHO)%array
      graph(2)%array_graph(IRHO)%descr =  u_np%array_var(IRHO)%descr

      graph(2)%array_graph(ITMP)%array => u_np%array_var(ITMP)%array
      graph(2)%array_graph(ITMP)%descr =  u_np%array_var(ITMP)%descr

      graph(2)%array_graph(IVX)%array => vx_cov
      graph(2)%array_graph(IVX)%descr = 'V_1'

      graph(2)%array_graph(IVY)%array => vy_cov
      graph(2)%array_graph(IVY)%descr = 'V_2'

      graph(2)%array_graph(IVZ)%array => vz_cov
      graph(2)%array_graph(IVZ)%descr = 'V_3'

      graph(2)%array_graph(IBX)%array => bx_cov
      graph(2)%array_graph(IBX)%descr = 'B_1'

      graph(2)%array_graph(IBY)%array => by_cov
      graph(2)%array_graph(IBY)%descr = 'B_2'

      graph(2)%array_graph(IBZ)%array => bz_cov
      graph(2)%array_graph(IBZ)%descr = 'B_3'

      graph(2)%array_graph(neqd+1)%array => jx_cov
      graph(2)%array_graph(neqd+1)%descr = 'J_1'

      graph(2)%array_graph(neqd+2)%array => jy_cov
      graph(2)%array_graph(neqd+2)%descr = 'J_2'

      graph(2)%array_graph(neqd+3)%array => jz_cov
      graph(2)%array_graph(neqd+3)%descr = 'J_3'

      graph(2)%array_graph(neqd+4:ngraph)%descr = ''

c Define graphics group #3: Cartesian variables

      graph(3)%cartesian=.true.
      graph(3)%descr='Car variables'

      graph(3)%array_graph(IRHO)%array => u_np%array_var(IRHO)%array
      graph(3)%array_graph(IRHO)%descr =  u_np%array_var(IRHO)%descr

      graph(3)%array_graph(ITMP)%array => u_np%array_var(ITMP)%array
      graph(3)%array_graph(ITMP)%descr =  u_np%array_var(ITMP)%descr

      graph(3)%array_graph(IVX)%array => vx_car
      graph(3)%array_graph(IVX)%descr = 'Vx'

      graph(3)%array_graph(IVY)%array => vy_car
      graph(3)%array_graph(IVY)%descr = 'Vy'

      graph(3)%array_graph(IVZ)%array => vz_car
      graph(3)%array_graph(IVZ)%descr = 'Vz'

      graph(3)%array_graph(IBX)%array => bx_car
      graph(3)%array_graph(IBX)%descr = 'Bx'

      graph(3)%array_graph(IBY)%array => by_car
      graph(3)%array_graph(IBY)%descr = 'By'

      graph(3)%array_graph(IBZ)%array => bz_car
      graph(3)%array_graph(IBZ)%descr = 'Bz'

      graph(3)%array_graph(neqd+1)%array => jx_car
      graph(3)%array_graph(neqd+1)%descr = 'Jx'

      graph(3)%array_graph(neqd+2)%array => jy_car
      graph(3)%array_graph(neqd+2)%descr = 'Jy'

      graph(3)%array_graph(neqd+3)%array => jz_car
      graph(3)%array_graph(neqd+3)%descr = 'Jz'

      graph(3)%array_graph(neqd+4:ngraph)%descr = ''

c Define graphics group #4: Diagnostics

      graph(4)%cartesian=.true.
      graph(4)%descr='Diagnostics'

      graph(4)%array_graph(1)%array => nuu
      graph(4)%array_graph(1)%descr = 'nu'

      graph(4)%array_graph(2)%array => eeta
      graph(4)%array_graph(2)%descr = 'eta'

      graph(4)%array_graph(3)%array => divrgB
      graph(4)%array_graph(3)%descr = 'local div(B)'

      graph(4)%array_graph(4)%array => divrgV
      graph(4)%array_graph(4)%descr = 'local div(V)'

      graph(4)%array_graph(5)%array => divrgJ
      graph(4)%array_graph(5)%descr = 'local div(J)'

      graph(4)%array_graph(6)%array => Pflux
      graph(4)%array_graph(6)%descr = 'Poloidal flux'

      graph(4)%array_graph(7)%array => qfactor
      graph(4)%array_graph(7)%descr = 'q factor'

      graph(4)%array_graph(8)%array => p_tot
      graph(4)%array_graph(8)%descr = 'Total pressure'

      graph(4)%array_graph(9:ngraph)%descr = ''

c Define graphics group #5: Perturbations

      graph(5)%cartesian=.false.
      graph(5)%descr='Perturbations'

      do ieq = 1,neqd
        graph(5)%array_graph(ieq)%array => u_graph%array_var(ieq)%array
        graph(5)%array_graph(ieq)%descr =
     $       trim(u_graph%array_var(ieq)%descr) // ' pert'
      enddo

      graph(5)%array_graph(neqd+1:ngraph)%descr = ''

c End program

      contains

c     defineGraphicsIO
c     ##################################################################
      subroutine defineGraphicsIO

        implicit none

        graphfile(1) = 'cnv.bin'
        ugraph(1)    = 20
        graphfile(2) = 'cov.bin'
        ugraph(2)    = 21
        graphfile(3) = 'car.bin'
        ugraph(3)    = 22
        graphfile(4) = 'diag.bin'
        ugraph(4)    = 23
        graphfile(5) = 'pert.bin'
        ugraph(5)    = 24

        profilefile(1)= 'cnv-profiles.bin'
        uprofile(1)   = 30
        profilefile(2)= 'cov-profiles.bin'
        uprofile(2)   = 31
        profilefile(3)= 'car-profiles.bin'
        uprofile(3)   = 32
        profilefile(4)= 'diag-profiles.bin'
        uprofile(4)   = 33
        profilefile(5)= 'pert-profiles.bin'
        uprofile(5)   = 34

        debugfile    = 'debug.bin'
        udebug       = 40
        
        drawgraph(1) = 'drawcnv.in'
        drawgraph(2) = 'drawcov.in'
        drawgraph(3) = 'drawcar.in'
        drawgraph(4) = 'drawdiag.in'
        drawgraph(5) = 'drawpert.in'

        drawprof(1)  = 'drawcnvprof.in'
        drawprof(2)  = 'drawcovprof.in'
        drawprof(3)  = 'drawcarprof.in'
        drawprof(4)  = 'drawdiagprof.in'
        drawprof(5)  = 'drawpertprof.in'

      end subroutine defineGraphicsIO

      end subroutine defineGraphics

c prepareTimeStepPlots
c####################################################################
      subroutine prepareTimeStepPlots

c--------------------------------------------------------------------
c     Set graphics files and dumping intervals
c--------------------------------------------------------------------

      use variables

      use graphics

      use nlfunction_setup

      use timeStepping

      use constants

      implicit none

c Call variables

c Local variables

      integer(4) :: i,j,k,ig,jg,kg
      real(8)    :: mm,kk,RR,ll,x1,y1,z1
      logical    :: cartsn

c Begin program

c Allocate variables

cc        allocate (divrgB(0:nxdp,0:nydp,0:nzdp)
cc     .           ,divrgV(0:nxdp,0:nydp,0:nzdp)
cc     .           ,divrgJ(0:nxdp,0:nydp,0:nzdp)
cc     .           ,Pflux (0:nxdp,0:nydp,0:nzdp)
cc     .           ,p_tot (0:nxdp,0:nydp,0:nzdp)
cc     .           ,qfactor(0:nxdp,0:nydp,0:nzdp))

c Find perturbed quantities

      u_graph = u_np - u_0

cdiag ******
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
cc     .                         ,cartsn)
cc            jac = jacobian(x1,y1,z1,cartsn)
cc
cc            u_graph%array_var(IBZ)%array(i,j,k) =
cccc     .           u_graph%array_var(IBZ)%array(i,j,k)/jac
cccc     .          *u_0    %array_var(IBZ)%array(i,j,k)/jac
cc     .           u_np%array_var(IBZ)%array(i,j,k)
cc     .          *u_np%array_var(IBX)%array(i,j,k)/jac
cccc            u_graph%array_var(IRHO)%array(i,j,k) =
cccc     .           u_graph%array_var(IRHO)%array(i,j,k)/jac**2
cc          enddo
cc        enddo
cc      enddo
cc
cc      call imposeBoundaryConditions(u_graph)
cc      if (source) u_graph = fsrc
cdiag ******

c Poloidal flux diagnostics  (use graphics limits)

      do k = kmin,kmax
        do i = imin,imax
          j = jmin
          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
          if (i == imin) then
            Pflux(i,j,k) = 0d0
          else
            Pflux(i,j,k) = Pflux(i-1,j,k) - by(i,j,k)*dxh(ig)
          endif
          do j = jmin+1,jmax
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            Pflux(i,j,k) = Pflux(i,j-1,k) + bx(i,j,k)*dyh(jg)
          enddo
        enddo
      enddo

c Poloidally averaged q-factor (use graphics limits)

      if (coords == 'hel') then
        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        do k = kmin,kmax
          do i = imax,imin+1,-1
            qfactor(i,jmin-1,k) = 0d0
            ll = 0d0
            do j = jmin,jmax
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
              qfactor(i,j,k) = qfactor(i,j-1,k)
     .             + mm*bz(i,j,k)/RR/(by(i,j,k)-kk*bz(i,j,k))*dyh(jg)
              ll = ll + dyh(jg)
            enddo
            qfactor(i,:,k) = qfactor(i,jmax,k)/ll
          enddo
          qfactor(imin,:,k) = qfactor(imin+1,:,k)
        enddo
      else
        qfactor = 0d0
      endif

c Divergence diagnostics

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartsn)
            jac = jacobian(x1,y1,z1,cartsn)
            divrgJ(i,j,k) = divJ(i,j,k)/jac
            divrgB(i,j,k) = divB(i,j,k)/jac
            divrgV(i,j,k) = divV(i,j,k)/jac
          enddo
        enddo
      enddo

      if (bcond(1) == SP) then
      !Find div at singular point
        call FillGhostNodes(IRHO,1,0,SP,divrgJ,zeros)
        call FillGhostNodes(IRHO,1,0,SP,divrgB,zeros)
        call FillGhostNodes(IRHO,1,0,SP,divrgV,zeros)

      !Replace div at innermost radius with value at SP
        divrgJ(1,:,:) = divrgJ(0,:,:)
        divrgB(1,:,:) = divrgB(0,:,:)
cc        divrgV(1,:,:) = divrgV(0,:,:)
      endif

c Total pressure (use graphics limits)

      p_tot = 0d0
      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax
            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartsn)
            jac = jacobian(x1,y1,z1,cartsn)
            
            p_tot(i,j,k) = (2.*rho(i,j,k)*tmp(i,j,k))
     .                     +(bx(i,j,k)*bx_cov(i,j,k)
     .                      +by(i,j,k)*by_cov(i,j,k)
     .                      +bz(i,j,k)*bz_cov(i,j,k))/2./jac

          enddo
        enddo
      enddo

c End program

      end subroutine prepareTimeStepPlots
