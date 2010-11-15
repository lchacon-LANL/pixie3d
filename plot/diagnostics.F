c module diag_setup
c####################################################################
      module diag_setup

        use parameters

        use variables

        use timeStepping

        use equilibrium

        use graphics

        use auxiliaryVariables

        use auxPlotVariables

        use transport_params

        use operators

        use nlfunction_setup

        use app_iosetup

        implicit none

        integer    :: ulineplot
        character*(20) :: lineplotfile

        character*(20) :: diag_desc(ngraph)

        integer        :: diag_ivar(ngraph)
     .                   ,diag_log(0:ngraph)

        integer    :: ndiag

        integer    :: drank,ddnp,split_key,MPI_COMM_DIAG,iplt,jplt,kplt
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

      integer :: igx,igy,igz,nx,ny,nz,bcs(6,neqd),ierr,ifile
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

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Check for external map

      if (coords == 'ext') call read_ext_map

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
        call MPI_Comm_size(MPI_COMM_DIAG,ddnp ,mpierr)
      endif

cc      write (*,*) 'Proc:',drank,ddnp
#else
      iplt = iplot
      jplt = jplot
      kplt = kplot
#endif

c Define BCs

      call defineBoundaryConditions(neqd,bcs)

c Define diagnostics

      call defineDiagnostics

c Create draw*.in file

      if (isProcessor) then
        ndiag = count (sel_diag /= 0)

        call createDrawInGfile(ndiag,ngraph,lineplotfile
     .           ,'Time histories'
     .           ,'time',sel_diag,diag_desc,diag_ivar,diag_log
     .           ,'drawgamma.in',.true.,.false.)
      endif

c Initialize diagnostic value references

      initialize = .true.
      call evaluateDiagnostics(varray,iigx,iigy,iigz,initialize)

c End 

      contains

c     read_ext_map
c     ###############################################################
      subroutine read_ext_map

        implicit none

        real(8),pointer,dimension(:,:,:,:) :: ext_map

        integer,pointer,dimension(:) :: mumap

        integer :: nxg,nyg,nzg

c     Begin program

        nxg = grid_params%nxgl(iigx)
        nyg = grid_params%nygl(iigy)
        nzg = grid_params%nzgl(iigz)

        ierr = openFiles(my_rank,np,u_mapf,map_file,mumap)

        if (ierr /= 0) then
          call pstop('evaluateDiagnostics'
     .              ,'Problem reading map files')
        endif

c     Fill GLOBAL map

        allocate(ext_map(0:nxg+1,0:nyg+1,0:nzg+1,3))

        do ifile=1,size(mumap)
          ierr=readVectorFromFiles(mumap(ifile),ext_map,.false.)
        enddo

        call closeFiles(mumap)

c     Fill LOCAL map and related metrics

        call defineGridMetric(grid_params,igr=iigx
     .           ,xcar=ext_map(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:))

        deallocate(ext_map)

        gmetric => grid_params%gmetric

      end subroutine read_ext_map

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

        if (.not.isProcessor) return

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
      diag_desc(neqd+18)= 'RMS |J|'
      diag_desc(neqd+19)= 'Max Vx'
      diag_desc(neqd+20)= 'Max Vy'
      diag_desc(neqd+21)= 'Max Vz'

      diag_desc(neqd+22:ngraph) = ''

c Define corresponding independent variables

      diag_ivar(:) = 0  !This means that the independent variable is time
                        !for all variables

c Specify log plots (log_plot = 1)

      diag_log = 0

      diag_log(1:neqd) = 1  !Log scale for perturbations of dependent variables

c End program

      end subroutine defineDiagnostics

c evaluateDiagnostics
c####################################################################
      subroutine evaluateDiagnostics(varray,iigx,iigy,iigz,init)

c--------------------------------------------------------------------
c     Calculates diagnostics.
c--------------------------------------------------------------------

      use diag_setup

      use local_BCS_variables

      implicit none

c Call variables

      integer          :: iigx,iigy,iigz
      logical          :: init

      type(var_array),pointer :: varray

c Local variables

      integer    :: i,j,k,ig,jg,kg,ieq,igx,igy,igz,nx,ny,nz
      real(8)    :: array (ilom:ihip,jlom:jhip,klom:khip)
     .             ,vv    (ilom:ihip,jlom:jhip,klom:khip,3)

      real(8)    :: energy,Vflux,Bflux,diverB,x1,y1,z1,dpert(neqd)
     .             ,tn=0d0,ez,lrflux,rflux,rmsj,vmax(3),lvmax(3)

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

c Code for explicit pinch flow BC

      if (lag_pinch_bc.and.pinch_flow) then

#if defined(vec_pot)
        save_Bn = .true.

        acov(:,:,:,1) = u_n%array_var(IAX)%array
        acov(:,:,:,2) = u_n%array_var(IAY)%array
        acov(:,:,:,3) = u_n%array_var(IAZ)%array

        bcnv_n = curl(nx,ny,nz,igx,igy,igz,acov)
#else
        save_Bn = .false.  !No need to save Bn, since u_n contains ghost cells

        bcnv_n(:,:,:,1) = u_n%array_var(IBX)%array
        bcnv_n(:,:,:,2) = u_n%array_var(IBY)%array
        bcnv_n(:,:,:,3) = u_n%array_var(IBZ)%array
#endif
      endif

c Code for electric field BC

cc      save_tn = .false.
      save_tn = .true.

cc      vcov_n(:,:,:,1) = u_n%array_var(IVX)%array
cc      vcov_n(:,:,:,2) = u_n%array_var(IVY)%array
cc      vcov_n(:,:,:,3) = u_n%array_var(IVZ)%array

      alpha = 1.-cnfactor

c Postprocess solution

      call postProcessSolution(varray,igx,igy,igz)

c Particle diagnostics

      Npar = integral(igx,igy,igz,nx,ny,nz,rho,average=.false.) - Npar0

c Momentum diagnostics (in Cartesian coords.)

      array = rho*vx_car
      Ipx = integral(igx,igy,igz,nx,ny,nz,array,average=.false.) - Ipx0
                                                        
      array = rho*vy_car
      Ipy = integral(igx,igy,igz,nx,ny,nz,array,average=.false.) - Ipy0
                                                        
      array = rho*vz_car
      Ipz = integral(igx,igy,igz,nx,ny,nz,array,average=.false.) - Ipz0

c Energy diagnostics

      !Magnetic energy

      array = (bx_car**2 + by_car**2 + bz_car**2)
      Em = 0.5*integral(igx,igy,igz,nx,ny,nz,array,average=.false.)

      !Ion kinetic energy

      array = rho*(vx_car**2 + vy_car**2 + vz_car**2)
      Ek = 0.5*integral(igx,igy,igz,nx,ny,nz,array,average=.false.)

      !Thermal energy

      if (adiabatic) then
        array = tmp
      else
        array = rho*tmp
      endif

      if (gamma /= 1d0) then
        Et = a_p/(gamma-1)
     .      *integral(igx,igy,igz,nx,ny,nz,array,average=.false.)
      else
        Et = 0d0
      endif

      !Total energy error

      energy = Ek + Em + Et - Ek0 - Em0 - Et0

c Current diagnostics

      Iz = integral(igx,igy,igz,nx,ny,nz,jz_car,average=.false.) - Iz0

      array =(jcnv(:,:,:,1)*jcov(:,:,:,1)
     .       +jcnv(:,:,:,2)*jcov(:,:,:,2)
     .       +jcnv(:,:,:,3)*jcov(:,:,:,3))/gmetric%grid(igx)%jac

      rmsj = integral(igx,igy,igz,nx,ny,nz,array,average=.false.)

      rmsj = sqrt(rmsj)

c Toroidal flux diagnostics

      Tflux = integral(igx,igy,igz,nx,ny,nz,bz_car,average=.false.)
     .      - Tflux0

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

c Normalize conservation errors

       if (Npar0 /= 0d0) Npar  = Npar /Npar0 
cc       if (Ipx0  /= 0d0) Ipx   = Ipx  /Ipx0  
cc       if (Ipy0  /= 0d0) Ipy   = Ipy  /Ipy0  
cc       if (Ipz0  /= 0d0) Ipz   = Ipz  /Ipz0  
       if (Iz0   /= 0d0) Iz    = Iz   /Iz0   
       if (Tflux0/= 0d0) Tflux = Tflux/Tflux0

       energy = energy/(Ek0 + Em0 + Et0)

c Velocity divergence diagnostics

      !Total flow divergence (conservation of flow)
      Vflux  = integral(igx,igy,igz,nx,ny,nz,divrgV,average=.false.)

c Magnetic divergence diagnostics

      !Total B divergence (conservation of flux)
      Bflux  = integral(igx,igy,igz,nx,ny,nz,divrgB,average=.false.)

      !Local divergence (statement of numerical accuracy)
      array(1:nx,1:ny,1:nz) = abs(divrgB(1:nx,1:ny,1:nz))

      diverB = integral(igx,igy,igz,nx,ny,nz,array,average=.true.)

c Parallel electric field and reconnected flux diagnostic (assumes periodic direction in Y)

      if (isProcessor) then

        if (time == tn) then
cc          tn = time - 1d5
cc          Pfluxn = Pflux
cc        endif
          ez = 0d0 ; rflux = 0d0
          Pfluxn = Pflux
        else

          select case(trim(equil))
          case('dtm')
            ez = (maxval(Pflux (1,:,1))-minval(Pflux (1,:,1))
     .           -maxval(Pfluxn(1,:,1))+minval(Pfluxn(1,:,1)))/(time-tn)

            rflux = abs(maxval(Pflux(1,:,1))-minval(Pflux(1,:,1)))
          case default
            ez = (Pflux (iplt,jplt,1)
     .           -Pfluxn(iplt,jplt,1))/(time-tn)

            rflux = abs(maxval(Pflux(iplt,:,1))-minval(Pflux(iplt,:,1)))
          end select

          tn = time
          Pfluxn = Pflux

        endif

      endif

ccc Reconnected flux (GEM challenge)
cc
cc      if (trim(equil) == 'gem') then
cc
cc#if defined(petsc)
cc        if (split_key /= MPI_UNDEFINED) then
cc          if (ddnp == 1) then
cc            rflux = Pflux(iplt,jplt      ,1)
cc     .            - Pflux(iplt,jplt-nyd/2,1)
cc          else
cc            lrflux = -(-1)**drank*Pflux(iplt,jplt,1)
cc
cc            call MPI_Allreduce(lrflux,rflux,1,MPI_DOUBLE_PRECISION
cc     .                        ,MPI_SUM,MPI_COMM_DIAG,mpierr)
cc          endif
cc        endif
cc#else
cc        rflux = abs(maxval(Pflux(1,:,1))-minval(Pflux(1,:,1)))
cc#endif
cc      else
cc        rflux = 0d0
cc      endif

c Growth rate diagnostics

      do ieq=1,neqd

        array = (varray%array_var(ieq)%array
     .          -u_0   %array_var(ieq)%array )**2

        dpert(ieq) = integral(igx,igy,igz,nx,ny,nz,array,average=.true.)

        if (dpert(ieq).eq.0d0) then
          dpert(ieq) = 0.1*maxval(abs(pert))
        else
          dpert(ieq) = sqrt(dpert(ieq))
        endif

      enddo

c Maximum velocity diagnostics

      vmax = (/ maxval(vx_car),maxval(vy_car),maxval(vz_car) /)

#if defined(petsc)      
      lvmax = vmax
      call MPI_Reduce(lvmax,vmax,3,MPI_DOUBLE_PRECISION
     .               ,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#endif

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
      diagnostics(neqd+18)= rmsj
      diagnostics(neqd+19)= vmax(1)
      diagnostics(neqd+20)= vmax(2)
      diagnostics(neqd+21)= vmax(3)

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
