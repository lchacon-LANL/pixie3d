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

        use auxPlotVariables

        use constants

        use transport_params

        use operators

        use generalPurposeFunctions

        use nlfunction_setup

        implicit none

        integer    :: ulineplot
        character*(20) :: lineplotfile

        character*(20) :: diag_desc(ngraph)

        integer        :: diag_ivar(ngraph)
     .                   ,diag_log(ngraph)

        integer    :: ndiag

        integer    :: drank,dnp,split_key,MPI_COMM_DIAG,iplt,jplt,kplt
        logical    :: isProcessor

        real(8) :: diagnostics(ngraph)

        real(8) :: Npar0,Ipx0,Ipy0,Ipz0,Ek0,Em0,Et0,Iz0,Tflux0
        real(8) :: Npar ,Ipx ,Ipy ,Ipz ,Em ,Ek ,Et ,Iz ,Tflux

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

      integer          :: iigx,iigy,iigz
      type(var_array),pointer :: varray

c Local variables

      logical :: initialize

c Interfaces

      INTERFACE
        subroutine evaluateDiagnostics(varray,iigx,iigy,iigz,init)
        use variable_setup
        integer    :: iigx,iigy,iigz
        logical    :: init
        type(var_array),pointer :: varray
        end subroutine evaluateDiagnostics
      END INTERFACE

c Begin program

      igx = iigx
      igy = iigy
      igz = iigz

c Define processor for local diagnostics

      !Identify local diagnostic proc
      isProcessor = inProc(iplot,jplot,kplot,igx,igy,igz)

#if defined(petsc)
      !Create local diagnostics communicator for line integral
      split_key = MPI_UNDEFINED

      if (inProc(iplot,jplot-nyd/2,kplot,igx,igy,igz)) then
        call fromGlobalToLocalLimits(iplot,jplot-nyd/2,kplot
     .                              ,iplt,jplt,kplt,igx,igy,igz)
        split_key = 1
      endif

      if (isProcessor) then
        call fromGlobalToLocalLimits(iplot,jplot,kplot,iplt,jplt,kplt
     .                              ,igx,igy,igz)
        split_key = 1
      endif

      call MPI_Comm_split(MPI_COMM_WORLD,split_key,my_rank
     .                   ,MPI_COMM_DIAG,mpierr)

      if (split_key /= MPI_UNDEFINED) then
        call MPI_Comm_rank(MPI_COMM_DIAG,drank,mpierr)
        call MPI_Comm_size(MPI_COMM_DIAG,dnp  ,mpierr)
      endif

cc      write (*,*) 'Proc:',drank,dnp
#else
      iplt = iplot
      jplt = jplot
      kplt = kplot
#endif

c Define diagnostics

      call defineDiagnostics

c Create draw*.in file

      if (isProcessor) then
        ndiag = count (sel_diag /= 0)

        call createDrawInGfile(ndiag,lineplotfile,'Time histories'
     .           ,'time',sel_diag,diag_desc,diag_ivar,diag_log
     .           ,'drawgamma.in',.false.)
      endif

c Initialize diagnostic value references

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

        character(30) :: ft_file
        integer :: uft=2345
        logical :: tor

c Begin program

        if (.not.isProcessor) return

c Finalize time traces output (separator)

        write(ulineplot)

c Close time traces file

        close(unit=ulineplot)

c Poincare output (field tracer)

#if !defined(petsc)
cc        write (*,*)
cc        write (*,*) 'Dumping Poincare file...'
cc        ft_file = 'poincare.bin'
cc        open(unit=uft,file=trim(ft_file)
cc     .      ,status='unknown',form='unformatted')
cc
cc        tor = (coords == 'vmc' .or. coords == 'tor')
cc
cc        write (uft) nxd,nyd,nzd,bcond,tor
cc        write (uft) grid_params%xx(1:nxd+2)
cc        write (uft) grid_params%yy(1:nyd+2)
cc        write (uft) grid_params%zz(1:nzd+2)
cc        write (uft) u_np%array_var(IBX)%array
cc        write (uft) u_np%array_var(IBY)%array
cc        write (uft) u_np%array_var(IBZ)%array
cc        write (uft) gmetric%grid(1)%jac
cc        write (uft) gmetric%grid(1)%car
cc
cc        close(uft)
#endif

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

      if (.not.isProcessor) return

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
#if !defined(vec_pot)
      diag_desc(IBX)    = '||dbx||'        
      diag_desc(IBY)    = '||dby||'        
      diag_desc(IBZ)    = '||dbz||'        
#else
      diag_desc(IAX)    = '||dax||'        
      diag_desc(IAY)    = '||day||'        
      diag_desc(IAZ)    = '||daz||'        
#endif
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
      diag_desc(neqd+16)= 'Ez'
      diag_desc(neqd+17)= 'Reconnected flux'

      diag_desc(neqd+18:ngraph) = ''

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

      integer          :: iigx,iigy,iigz
      logical          :: init

      type(var_array),pointer :: varray

c Local variables

      integer    :: i,j,k,ig,jg,kg,ieq
      real(8)    :: array (ilom:ihip,jlom:jhip,klom:khip)
     .             ,Pfluxn(ilom:ihip,jlom:jhip,klom:khip)

      real(8)    :: energy,Vflux,Bflux,diverB,x1,y1,z1,dpert(neqd)
     .             ,tn=0d0,ez,lrflux,rflux
      logical    :: cartsn

      INTERFACE
        subroutine postProcessSolution(varray,iigx,iigy,iigz)
        use variable_setup
        integer    :: iigx,iigy,iigz
        type(var_array),pointer :: varray
        end subroutine postProcessSolution
      END INTERFACE

c Begin program

      igx = iigx
      igy = iigy
      igz = iigz

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      if (init) then
        Npar0 = 0d0
        Ipx0  = 0d0
        Ipy0  = 0d0
        Ipz0  = 0d0
        Em0   = 0d0
        Ek0   = 0d0
        Et0   = 0d0
        Iz0   = 0d0
        Tflux0= 0d0
      endif

      rho => varray%array_var(IRHO)%array
      px  => varray%array_var(IVX )%array
      py  => varray%array_var(IVY )%array
      pz  => varray%array_var(IVZ )%array
#if !defined(vec_pot)
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
#else
      ax  => varray%array_var(IAX )%array
      ay  => varray%array_var(IAY )%array
      az  => varray%array_var(IAZ )%array
#endif
      tmp => varray%array_var(ITMP)%array

c Postprocess solution

      call postProcessSolution(varray,igx,igy,igz)

c Particle diagnostics

      Npar = integral(nx,ny,nz,rho,igx,igy,igz,.false.) - Npar0

c Momentum diagnostics (in Cartesian coords.)

      array = rho*vx_car
      Ipx = integral(nx,ny,nz,array,igx,igy,igz,.false.) - Ipx0
                                                        
      array = rho*vy_car
      Ipy = integral(nx,ny,nz,array,igx,igy,igz,.false.) - Ipy0
                                                        
      array = rho*vz_car
      Ipz = integral(nx,ny,nz,array,igx,igy,igz,.false.) - Ipz0

c Energy diagnostics

      !Magnetic energy

      array = (bx_car**2 + by_car**2 + bz_car**2)
      Em = 0.5*integral(nx,ny,nz,array,igx,igy,igz,.false.)

      !Ion kinetic energy

      array = rho*(vx_car**2 + vy_car**2 + vz_car**2)
      Ek = 0.5*integral(nx,ny,nz,array,igx,igy,igz,.false.)

      !Thermal energy

      array = rho*tmp

      if (gamma /= 1d0) then
        Et = a_p/(gamma-1)*integral(nx,ny,nz,array,igx,igy,igz,.false.)
      else
        Et = 0d0
      endif

      !Total energy

      energy = Ek + Em + Et - Ek0 - Em0 - Et0

c Current diagnostics

      Iz = integral(nx,ny,1,jz_car,igx,igy,igz,.false.) - Iz0

c Toroidal flux diagnostics

      Tflux = integral(nx,ny,1,bz_car,igx,igy,igz,.false.) - Tflux0

c Initial assignments

      if (init) then
        Npar0 = Npar
        Ipx0  = Ipx  
        Ipy0  = Ipy  
        Ipz0  = Ipz  
        Em0   = Em  
        Ek0   = Ek  
        Et0   = Et 
        Iz0   = Iz
        Tflux0= Tflux
        return
      endif

c Velocity divergence diagnostics

      !Total flow divergence (conservation of flow)
      Vflux  = integral(nx,ny,nz,divrgV,igx,igy,igz,.false.)

c Magnetic divergence diagnostics

      !Total B divergence (conservation of flux)
      Bflux  = integral(nx,ny,nz,divrgB,igx,igy,igz,.false.)

      !Local divergence (statement of numerical accuracy)
      array(1:nx,1:ny,1:nz) = abs(divrgB(1:nx,1:ny,1:nz))

      diverB = integral(nx,ny,nz,array,igx,igy,igz,.true.)

c Parallel electric field diagnostic

      if (isProcessor) then
        if (time == tn) then
          tn = time - 1d5
          Pfluxn(iplt,jplt,1) = Pflux(iplt,jplt,1)
        endif

        ez = (Pflux (iplt,jplt,1)
     .       -Pfluxn(iplt,jplt,1))/(time-tn)

        tn = time
        Pfluxn(iplt,jplt,1) = Pflux(iplt,jplt,1)
      endif

c Reconnected flux (GEM challenge)

#if defined(petsc)
      if (split_key /= MPI_UNDEFINED) then
        if (dnp == 1) then
          rflux = Pflux(iplt,jplt      ,1)
     .          - Pflux(iplt,jplt-nyd/2,1)
        else
          lrflux = -(-1)**drank*Pflux(iplt,jplt,1)

          call MPI_Allreduce(lrflux,rflux,1,MPI_DOUBLE_PRECISION
     .                      ,MPI_SUM,MPI_COMM_DIAG,mpierr)
        endif
      endif
#else
      rflux = minval(Pflux(iplt,1:nyd,1))
     .      - maxval(Pflux(iplt,1:nyd,1))
#endif

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
      diagnostics(neqd+10)= Ipx
      diagnostics(neqd+11)= Ipy
      diagnostics(neqd+12)= Ipz
      diagnostics(neqd+13)= Vflux
      diagnostics(neqd+14)= Iz
      diagnostics(neqd+15)= Tflux
      diagnostics(neqd+16)= Ez
      diagnostics(neqd+17)= rflux

cc      do i=1,count(diag_desc /= "")
cc        write (*,*) diag_desc(i),diagnostics(i)
cc      enddo

c Dump diagnostics

      if (isProcessor) then
        write(ulineplot) real(time),real(diagnostics)
cc      write(*,*) real(time),real(diagnostics)
      endif

c End 

      end subroutine evaluateDiagnostics
