c defineGraphics
c####################################################################
      subroutine defineGraphics(igx,igy,igz)

c--------------------------------------------------------------------
c     Set graphics files and dumping intervals
c--------------------------------------------------------------------

      use parameters

      use variables

      use graphics

      use auxiliaryVariables

      use equilibrium

      implicit none

c Call variables

      integer(4) :: igx,igy,igz

c Local variables

      integer(4) :: ieq,i,nx,ny,nz

c Begin program

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Allocate auxiliary graphics variables

      allocate (bx_car(0:nx+1,0:ny+1,0:nz+1)
     .         ,by_car(0:nx+1,0:ny+1,0:nz+1)
     .         ,bz_car(0:nx+1,0:ny+1,0:nz+1)
     .         ,jx_car(0:nx+1,0:ny+1,0:nz+1)
     .         ,jy_car(0:nx+1,0:ny+1,0:nz+1)
     .         ,jz_car(0:nx+1,0:ny+1,0:nz+1)
     .         ,vx_car(0:nx+1,0:ny+1,0:nz+1)
     .         ,vy_car(0:nx+1,0:ny+1,0:nz+1)
     .         ,vz_car(0:nx+1,0:ny+1,0:nz+1)
     .         ,divrgB(0:nx+1,0:ny+1,0:nz+1)
     .         ,divrgV(0:nx+1,0:ny+1,0:nz+1)
     .         ,divrgJ(0:nx+1,0:ny+1,0:nz+1)
     .         ,Pflux (0:nx+1,0:ny+1,0:nz+1)
     .         ,p_tot (0:nx+1,0:ny+1,0:nz+1)
     .         ,qfactor(0:nx+1,0:ny+1,0:nz+1))

c Define number of graphics groups

      ngroups = 5

c Allocate graphics arrays

      call allocateGraphicsVariables(ngroups)

c Define I/O

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

      sel_gr(1,:) = sel_graph

      prof_ivar(1,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (1,:) = 0        !No log scales
      prof_spline(1) = .false.  !No splines

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

      sel_gr(2,:) = sel_graph

      prof_ivar(2,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (2,:) = 0        !No log scales
      prof_spline(2) = .false.  !No splines

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

      sel_gr(3,:) = sel_graph

      prof_ivar(3,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (3,:) = 0        !No log scales
      prof_spline(3) = .false.  !No splines

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

      sel_gr(4,:) = (/ (i,i=1,8),0 /)

      prof_ivar(4,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (4,:) = 0        !No log scales
      prof_spline(4) = .false.  !No splines

c Define graphics group #5: Perturbations

      graph(5)%cartesian=.false.
      graph(5)%descr='Perturbations'

      do ieq = 1,neqd
        graph(5)%array_graph(ieq)%array => u_graph%array_var(ieq)%array
        graph(5)%array_graph(ieq)%descr =
     $       trim(u_graph%array_var(ieq)%descr) // ' pert'
      enddo

      graph(5)%array_graph(neqd+1:ngraph)%descr = ''

      sel_gr(5,:) = (/ (i,i=1,neqd),(0,i=neqd+1,9) /)

      prof_ivar(5,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (5,:) = 0        !No log scales
      prof_spline(5) = .false.  !No splines

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
      subroutine prepareTimeStepPlots(iigx,iigy,iigz)

c--------------------------------------------------------------------
c     Set graphics files and dumping intervals
c--------------------------------------------------------------------

      use variables

      use graphics

      use auxiliaryVariables

      use timeStepping

      use constants

      use operators

      use equilibrium

      implicit none

c Call variables

      integer(4) :: iigx,iigy,iigz

c Local variables

      integer(4) :: i,j,k,ig,jg,kg
      real(8)    :: mm,kk,RR,ll,x1,y1,z1,jac
      logical    :: cartsn

c Begin program

      igx = iigx
      igy = iigy
      igz = iigz

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Impose boundary conditions

      call imposeBoundaryConditions(u_np,igx,igy,igz)

      rho => u_np%array_var(IRHO)%array
      rvx => u_np%array_var(IVX )%array
      rvy => u_np%array_var(IVY )%array
      rvz => u_np%array_var(IVZ )%array
      bx  => u_np%array_var(IBX )%array
      by  => u_np%array_var(IBY )%array
      bz  => u_np%array_var(IBZ )%array
      tmp => u_np%array_var(ITMP)%array

c Find perturbed quantities

      u_graph = u_np - u_0

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
            divrgJ(i,j,k) = div(i,j,k,jx,jy,jz)/jac
            divrgB(i,j,k) = div(i,j,k,bx,by,bz)/jac
            divrgV(i,j,k) = div(i,j,k,vx,vy,vz)/jac
          enddo
        enddo
      enddo

      if (bcond(1) == SP) then
      !Find div at singular point
cc        call FillGhostNodes(IRHO,1,0,SP,divrgJ,zeros)
cc        call FillGhostNodes(IRHO,1,0,SP,divrgB,zeros)
cc        call FillGhostNodes(IRHO,1,0,SP,divrgV,zeros)

      !Replace div at innermost radius with value at SP
cc        divrgJ(1,:,:) = divrgJ(0,:,:)
cc        divrgB(1,:,:) = divrgB(0,:,:)
cc        divrgV(1,:,:) = divrgV(0,:,:)

        divrgJ(1,:,:) = 0d0
        divrgB(1,:,:) = 0d0
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
