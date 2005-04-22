c module diag_setup
c####################################################################
      module diag_setup

        use parameters

        use variables

        use timeStepping

        use equilibrium

        use graphics

        use icond

        use auxiliaryVariables

        use constants

        use transport_params

        use operators

        use generalPurposeFunctions

        use nlfunction_setup

        implicit none

        integer(4) :: ulineplot
        character*(20) :: lineplotfile

        character*(20) :: diag_desc(ngraph)

        integer(4)     :: diag_ivar(ngraph)
     .                   ,diag_log(ngraph)

        integer(4) :: ndiag

        real(8) :: diagnostics(ngraph)

        real(8) :: Npar0,px0,py0,pz0,Ek0,Em0,Et0,Iz0,Tflux0
        real(8) :: Npar ,px ,py ,pz ,Em ,Ek ,Et ,Iz ,Tflux

      end module diag_setup

c initializeDiagnostics
c####################################################################
      subroutine initializeDiagnostics(varray,iigx,iigy,iigz)

c--------------------------------------------------------------------
c     Initializes diagnostics.
c--------------------------------------------------------------------

      use diag_setup

      implicit none

c Call variables

      integer(4)       :: iigx,iigy,iigz
      type (var_array) :: varray

c Local variables

      logical :: initialize

c Begin program

c Define diagnostics and create draw*.in file

      if (isProc(0)) then
        call defineDiagnostics
        ndiag = count (sel_diag /= 0)

        call createDrawInGfile(ndiag,lineplotfile,'Time histories'
     .           ,'time',sel_diag,diag_desc,diag_ivar,diag_log
     .           ,'drawgamma.in',.false.)
      endif

c Initialize diagnostic references

      initialize = .true.
      call evaluateDiagnostics(varray,iigx,iigy,iigz,initialize)

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

        if (.not.isProc(0)) return

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

      lineplotfile = 'drawtplots.bin'
      ulineplot    = 10

      open(unit=ulineplot,file=lineplotfile,form='unformatted'
     .    ,status='replace')

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
      subroutine evaluateDiagnostics(varray,iigx,iigy,iigz,init)

c--------------------------------------------------------------------
c     Calculates diagnostics.
c--------------------------------------------------------------------

      use diag_setup

      implicit none

c Call variables

      integer(4)       :: iigx,iigy,iigz
      logical          :: init

      type (var_array) :: varray

c Local variables

      integer(4) :: i,j,k,ig,jg,kg,ieq
      real(8)    :: array(ilom:ihip,jlom:jhip,klom:khip)

      real(8)    :: energy,Vflux,Bflux,diverB,x1,y1,z1,dpert(neqd)
      logical    :: cartsn

c Begin program

      igx = iigx
      igy = iigy
      igz = iigz

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

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

      rho => varray%array_var(IRHO)%array
      rvx => varray%array_var(IVX )%array
      rvy => varray%array_var(IVY )%array
      rvz => varray%array_var(IVZ )%array
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
      tmp => varray%array_var(ITMP)%array

c Particle diagnostics

      Npar = integral(nx,ny,nz,rho,igx,igy,igz,.false.) - Npar0

c Momentum diagnostics (in Cartesian coords.)

      array = rho*vx_car
      px = integral(nx,ny,nz,array,igx,igy,igz,.false.) - px0
                                                        
      array = rho*vy_car
      py = integral(nx,ny,nz,array,igx,igy,igz,.false.) - py0
                                                        
      array = rho*vz_car
      pz = integral(nx,ny,nz,array,igx,igy,igz,.false.) - pz0

c Energy diagnostics

      !Magnetic energy

      array = rho*(bx_car**2 + by_car**2 + bz_car**2)
      Em = 0.5*integral(nx,ny,nz,array,igx,igy,igz,.false.) - Em0

      !Ion kinetic energy

      array = rho*(vx_car**2 + vy_car**2 + vz_car**2)
      Ek = 0.5*integral(nx,ny,nz,array,igx,igy,igz,.false.) - Ek0

      !Thermal energy

      if (gamma /= 1d0) then
        array = 2*rho*tmp/(gamma-1.)
      else
        array = 0d0
      endif
      Et = integral(nx,ny,nz,array,igx,igy,igz,.false.) - Et0

      !Total energy

      energy = Ek + Em + Et

c Current diagnostics

      Iz = integral(nx,ny,1,jz_car,igx,igy,igz,.false.) - Iz0

c Toroidal flux diagnostics

      Tflux = integral(nx,ny,1,bz_car,igx,igy,igz,.false.) - Tflux0

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

c Velocity divergence diagnostics

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            array(i,j,k) = div(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz)
          enddo
        enddo
      enddo

      !Total flow divergence (conservation of flow)
      Vflux  = integral(nx,ny,nz,array,igx,igy,igz,.false.)

c Magnetic divergence diagnostics

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            array(i,j,k) = div(i,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
          enddo
        enddo
      enddo

      !Total B divergence (conservation of flux)
      Bflux  = integral(nx,ny,nz,array,igx,igy,igz,.false.)

      !Local divergence (statement of numerical accuracy)
      array = abs(array)

      diverB = integral(nx,ny,nz,array,igx,igy,igz,.true.)

c Growth rate diagnostics

      do ieq=1,neqd

        array = (varray%array_var(ieq)%array
     .          -u_0   %array_var(ieq)%array )**2

        dpert(ieq) = integral(nx,ny,nz,array,igx,igy,igz,.true.)

cc        if (dpert(ieq).gt.0d0) dpert(ieq) = log(sqrt(dpert(ieq)))
        if (dpert(ieq).eq.0d0) then
          dpert(ieq) = 0.1*maxval(abs(pert))
        else
          dpert(ieq) = sqrt(dpert(ieq))
        endif

      enddo

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

      if (.not.isProc(0)) return

      write(ulineplot) real(time),real(diagnostics)
cc      write(*,*) real(time),real(diagnostics)

c End program

      end subroutine dumpDiagnostics
