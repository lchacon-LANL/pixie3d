c module diag_setup
c####################################################################
      module diag_setup

        use parameters

        use variables

        use timeStepping

        use equilibrium

        use graphics

        use nlfunction_setup

        implicit none

        integer(4) :: ulineplot
        character*(20) :: lineplotfile

        integer(4) :: ndiag

        real(8) :: Npar0,px0,py0,pz0,Ek0,Em0,Et0,Iz0,Tflux0
        real(8) :: Npar ,px ,py ,pz ,Em ,Ek ,Et ,Iz ,Tflux

      end module diag_setup

c initializeDiagnostics
c####################################################################
      subroutine initializeDiagnostics

c--------------------------------------------------------------------
c     Initializes diagnostics.
c--------------------------------------------------------------------

      use diag_setup

      implicit none

c Call variables

c Local variables

      logical :: initialize

c Begin program

c Define diagnostics

      call defineDiagnostics

c Initialize diagnostic references

      initialize = .true.
      call evaluateDiagnostics(u_0,initialize)

c Calculate initial diagnostics

      initialize = .false.
      if (.not.restart) then
        call evaluateDiagnostics(u_n,initialize)
      else
        call dumpDiagnostics
      endif

c Create draw*.in file

      ndiag = count (sel_diag /= 0)

      call createDrawInGfile(ndiag,lineplotfile,'Time histories'
     .           ,'time',sel_diag,diag_desc,diag_ivar,diag_log
     .           ,'drawgamma.in',.false.)

c End 

      end subroutine initializeDiagnostics

c finalizeDiagnostics
c####################################################################
      subroutine finalizeDiagnostics

c--------------------------------------------------------------------
c     Close diagnostic files
c--------------------------------------------------------------------

        use diag_setup

        implicit none

c Call variables

c Local variables

c Begin program

c Finalize time traces output (separator)

        write(ulineplot)

c Close time traces file

        close(unit=ulineplot)

c End program

      end subroutine finalizeDiagnostics

c defineDiagnostics
c####################################################################
      subroutine defineDiagnostics

c--------------------------------------------------------------------
c     Defines and sets up diagnostics
c--------------------------------------------------------------------

      use diag_setup

      implicit none

c Call variables

c Local variables

c Begin program

c Open files

      lineplotfile = 'timeplots.bin'
      ulineplot    = 10

      if (.not.restart) then
        open(unit=ulineplot,file=lineplotfile,form='unformatted'
     .      ,status='replace')
      else
        open(unit=ulineplot,file=lineplotfile,form='unformatted'
     .      ,status='old',position='append')
      endif

c Define diagnostics

      diag_desc(IRHO)   = '||drho||'       
      diag_desc(IVX)    = '||dpx||'        
      diag_desc(IVY)    = '||dpy||'        
      diag_desc(IVZ)    = '||dpz||'        
      diag_desc(IBX)    = '||dbx||'        
      diag_desc(IBY)    = '||dby||'        
      diag_desc(IBZ)    = '||dbz||'        
      diag_desc(ITMP)   = '||dtmp||'       
      diag_desc(neqd+1) = 'Magnetic energy'
      diag_desc(neqd+2) = 'Kinetic energy' 
      diag_desc(neqd+3) = 'Thermal energy' 
      diag_desc(neqd+4) = 'Total energy'   
      diag_desc(neqd+5) = 'Time step'      
      diag_desc(neqd+6) = 'Growth rate'    
      diag_desc(neqd+7) = '|div(B)|'
      diag_desc(neqd+8) = 'Conservation of flux'      
      diag_desc(neqd+9) = 'Total particles'
      diag_desc(neqd+10)= 'Total X momentum'
      diag_desc(neqd+11)= 'Total Y momentum'
      diag_desc(neqd+12)= 'Total Z momentum'
      diag_desc(neqd+13)= 'Vflux at boundaries'
      diag_desc(neqd+14)= 'Toroidal current Iz'
      diag_desc(neqd+15)= 'Toroidal flux'

      diag_desc(neqd+16:ngraph) = ''

c Define corresponding independent variables

      diag_ivar(:) = 0  !This means that the independent variable is time
                        !for all variables

c Specify log plots (log_plot = 1)

      diag_log(1:neqd) = 1  !Log scale for perturbations of dependent variables

      diag_log(neqd+1:ngraph) = 0

c End program

      end subroutine defineDiagnostics

c evaluateDiagnostics
c####################################################################
      subroutine evaluateDiagnostics(varray,init)

c--------------------------------------------------------------------
c     Calculates diagnostics.
c--------------------------------------------------------------------

      use diag_setup

      implicit none

c Call variables

      logical          :: init

      type (var_array) :: varray

c Local variables

      integer(4) :: i,j,k,ig,jg,kg,ieq
      real(8)    :: dpert(neqd),mag(neqd),dmag1,dmag2
      real(8)    :: array(0:nxdp,0:nydp,0:nzdp)

      real(8)    :: energy,Vflux,Bflux,diverB,x1,y1,z1
      logical    :: cartsn

c Externals

      real(8)    :: integral
      external   :: integral

c Begin program

      if (init) then
        Npar0 = 0d0
        px0   = 0d0
        py0   = 0d0
        pz0   = 0d0
        Em0   = 0d0
        Ek0   = 0d0
        Et0   = 0d0
        Iz0   = 0d0
        Tflux0= 0d0
      endif

c Transform auxiliary vector fields to cartesian coordinates

      call postProcessSolution(varray)

c Find maximum velocity and magnetic field components (cartesian)

      vx_max = maxval(abs(vx_car))
      vy_max = maxval(abs(vy_car))
      vz_max = maxval(abs(vz_car))

      bx_max = maxval(abs(bx_car))
      by_max = maxval(abs(by_car))
      bz_max = maxval(abs(bz_car))

c Particle diagnostics

      Npar = integral(nxd,nyd,nzd,rho,igx,igy,igz,.false.) - Npar0

c Momentum diagnostics (in Cartesian coords.)

      array = rho*vx_car
      px = integral(nxd,nyd,nzd,array,igx,igy,igz,.false.) - px0
                                                        
      array = rho*vy_car
      py = integral(nxd,nyd,nzd,array,igx,igy,igz,.false.) - py0
                                                        
      array = rho*vz_car
      pz = integral(nxd,nyd,nzd,array,igx,igy,igz,.false.) - pz0

c Energy diagnostics

      !Magnetic energy

      array = rho*(bx_car**2 + by_car**2 + bz_car**2)
      Em = 0.5*integral(nxd,nyd,nzd,array,igx,igy,igz,.false.) - Em0

      !Ion kinetic energy

      array = rho*(vx_car**2 + vy_car**2 + vz_car**2)
      Ek = 0.5*integral(nxd,nyd,nzd,array,igx,igy,igz,.false.) - Ek0

      !Thermal energy

      if (gamma /= 1d0) then
        array = 2*rho*tmp/(gamma-1.)
      else
        array = 0d0
      endif
      Et = integral(nxd,nyd,nzd,array,igx,igy,igz,.false.) - Et0

      !Total energy

      energy = Ek + Em + Et

c Current diagnostics

      Iz = integral(nxd,nyd,1,jz_car,igx,igy,igz,.false.) - Iz0

c Toroidal flux diagnostics

      Tflux = integral(nxd,nyd,1,bz_car,igx,igy,igz,.false.) - Tflux0

c Initial assignments

      if (init) then
        Npar0 = Npar
        px0   = px  
        py0   = py  
        pz0   = pz  
        Em0   = Em  
        Ek0   = Ek  
        Et0   = Et 
        Iz0   = Iz
        Tflux0= Tflux
        return
      endif

c Divergence diagnostics

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartsn)
            jac = jacobian(x1,y1,z1,cartsn)
            divrgB(i,j,k) = divB(i,j,k)/jac
            divrgV(i,j,k) = divV(i,j,k)/jac
          enddo
        enddo
      enddo

      !Total B divergence (conservation of flux)
      Bflux  = integral(nxd,nyd,nzd,divrgB,igx,igy,igz,.false.)

      !Local divergence (statement of numerical accuracy)
      divrgB = abs(divrgB)
      diverB = integral(nxd,nyd,nzd,divrgB,igx,igy,igz,.true.)

      !Total flow divergence (conservation of flow)
      Vflux  = integral(nxd,nyd,nzd,divrgV,igx,igy,igz,.false.)

c Growth rate diagnostics

      do ieq=1,neqd

        array = (varray%array_var(ieq)%array
     .          -u_0   %array_var(ieq)%array )**2

        dpert(ieq) = integral(nxd,nyd,nzd,array,igx,igy,igz,.true.)

cc        if (dpert(ieq).gt.0d0) dpert(ieq) = log(sqrt(dpert(ieq)))
        if (dpert(ieq).eq.0d0) then
          dpert(ieq) = 1d-20
        else
          dpert(ieq) = sqrt(dpert(ieq))
        endif

      enddo

c Calculation of local growth rate for CN

      do ieq=1,neqd

        array = (u_n%array_var(ieq)%array
     .          -u_0%array_var(ieq)%array )**2

        dmag1 = integral(nxd,nyd,nzd,array,igx,igy,igz,.true.)

        array = (varray%array_var(ieq)%array
     .          +u_n   %array_var(ieq)%array
     .       -2.*u_0   %array_var(ieq)%array )**2

        dmag2 = integral(nxd,nyd,nzd,array,igx,igy,igz,.true.)

        if (dpert(ieq) /= 0d0.and.dmag2 /= 0d0) then
cc          mag(ieq) = .5*dt*sqrt(dmag2)/(exp(dpert(ieq))-sqrt(dmag1))
          mag(ieq) = .5*dt*sqrt(dmag2)/(dpert(ieq)-sqrt(dmag1))
        else
          mag(ieq) = 1e30
        endif

      enddo

      gammat = 1./minval(abs(mag))

c Diagnostic assignment

      diagnostics(1:neqd) = dpert(:)
      diagnostics(neqd+1) = Em
      diagnostics(neqd+2) = Ek
      diagnostics(neqd+3) = Et
      diagnostics(neqd+4) = energy
      diagnostics(neqd+5) = dt
      diagnostics(neqd+6) = gammat
      diagnostics(neqd+7) = diverB
      diagnostics(neqd+8) = Bflux
      diagnostics(neqd+9) = Npar
      diagnostics(neqd+10)= px
      diagnostics(neqd+11)= py
      diagnostics(neqd+12)= pz
      diagnostics(neqd+13)= Vflux
      diagnostics(neqd+14)= Iz
      diagnostics(neqd+15)= Tflux

cc      do i=1,count(diag_desc /= "")
cc        write (*,*) diag_desc(i),diagnostics(i)
cc      enddo

c Dump diagnostics

      call dumpDiagnostics

c End 

      end subroutine evaluateDiagnostics

c dumpDiagnostics
c######################################################################
      subroutine dumpDiagnostics

c -------------------------------------------------------------------
c     Dumps time plots
c -------------------------------------------------------------------

      use diag_setup

      implicit none

c Call variables

c Local variables
 
c Begin program

      write(ulineplot) real(time),real(diagnostics)

c End program

      end subroutine dumpDiagnostics

c integral
c ####################################################################
      real(8) function integral(nx,ny,nz,array,igx,igy,igz,avg)

c -------------------------------------------------------------------
c     Integrates array(i,j) on domain (nx)x(ny) and returns result.
c -------------------------------------------------------------------

      use grid

      implicit none

c Call variables

      integer(4) :: igx,igy,igz,nx,ny,nz
      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1)
      logical    :: avg

c Local variables

      integer(4) :: i,j,k
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8)    :: tvolume,vol

c Begin program

c Find integral limits

      imin = 1
      imax = nx

      jmin = 1
      jmax = ny

      kmin = 1
      kmax = nz

c Integrate

      integral = 0d0
      tvolume  = 0d0

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax
            vol = volume(i,j,k,igx,igy,igz)
            integral = integral + array(i,j,k)*vol
            tvolume = tvolume + vol
          enddo
        enddo
      enddo

      if (avg) integral = integral/tvolume

c End 

      end function integral

c postProcessSolution
c ####################################################################
      subroutine postProcessSolution(varray)

c -------------------------------------------------------------------
c     Finds cartesian components of vector fields.
c -------------------------------------------------------------------

      use diag_setup

      implicit none

c Call variables

      type (var_array) :: varray

c Local variables

      integer(4) :: i,j,k
      logical :: covariant,to_cartsn,to_cnv

c Begin program

      call imposeBoundaryConditions(varray)

      to_cartsn = .true.

c Find cartesian velocity components

      covariant    = .false.

      vx_car = vx
      vy_car = vy
      vz_car = vz
      call transformVector(igx,igy,igz
     .                    ,0,nx+1,0,ny+1,0,nz+1
     .                    ,vx_car,vy_car,vz_car,covariant,to_cartsn)

c Find covariant velocity components

      to_cnv = .false.

      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,vx_cov(i,j,k),vy_cov(i,j,k),vz_cov(i,j,k)
     .             ,vx(i,j,k),vy(i,j,k),vz(i,j,k),to_cnv)
          enddo
        enddo
      enddo

c Find magnetic field components in cartesian coordinates

      covariant    = .false.

      bx_car = bx
      by_car = by
      bz_car = bz
      call transformVector(igx,igy,igz
     .                    ,0,nx+1,0,ny+1,0,nz+1
     .                    ,bx_car,by_car,bz_car,covariant,to_cartsn)

c Find current components in cartesian coordinates

      covariant    = .false.

      jx_car = jx
      jy_car = jy
      jz_car = jz
      call transformVector(igx,igy,igz
     .                    ,0,nx+1,0,ny+1,0,nz+1
     .                    ,jx_car,jy_car,jz_car,covariant,to_cartsn)

c End 

      end subroutine postProcessSolution

