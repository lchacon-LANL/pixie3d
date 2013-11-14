c module diag_setup
c####################################################################
      module diag_setup

        use parameters

        use variables

        use ts_setup

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

        integer :: diag_ivar(ngraph)
     .            ,diag_log(0:ngraph)

        integer :: ndiag

#if defined(petsc)
        integer :: drank,ddnp,split_key,MPI_COMM_DIAG=MPI_COMM_NULL
#endif
        integer :: iplt,jplt,kplt
        logical :: isProcessor

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

      integer :: igx,igy,igz,nx,ny,nz,ierr,ifile
      logical :: initialize

c Interfaces

      INTERFACE
        subroutine evaluateDiagnostics(varray,iigx,iigy,iigz,init)
        use var_setup
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

c Create diagnostic communicator

      select case(trim(equil))
      case('dtm','gem')

#if defined(petsc)
        split_key = MPI_UNDEFINED

        if (inProc(iplot,grid_params%jlo(igy),kplot,igx,igy,igz)) then
          call fromGlobalToLocalLimits(iplot,grid_params%jlo(igy),kplot
     .                                ,iplt,jplt,kplt,igx,igy,igz)
          split_key = 1
        endif

        call MPI_Comm_split(MPI_COMM_WORLD,split_key,my_rank
     .                     ,MPI_COMM_DIAG,mpierr)

        if (split_key /= MPI_UNDEFINED) then
          call MPI_Comm_rank(MPI_COMM_DIAG,drank,mpierr)
          call MPI_Comm_size(MPI_COMM_DIAG,ddnp ,mpierr)
          isProcessor = (drank == 0)
        else
          isProcessor = .false.
        endif
#else
        isProcessor = (my_rank == 0)

        iplt = iplot
        jplt = jplot
        kplt = kplot
#endif

      case('ic-hk','ic')

      case default

        !Identify local diagnostic proc
        isProcessor = inProc(iplot,jplot,kplot,igx,igy,igz)

        if (isProcessor) then
          call fromGlobalToLocalLimits(iplot,jplot,kplot,iplt,jplt,kplt
     .                                ,igx,igy,igz)
        endif

      end select

c Define BCs

      pinch_flow = (    equil == 'ppnch' .or. equil == 'ppn3d'
     .              .or.equil == 'ppnsl' .or. equil == 'ppnst'
     .              .or.equil == '2dhel' .or. equil == 'p3nsl'
     .              .or.equil == 'ohtor')

cc      call defineDepBCs(neqd,bcs)

c Define diagnostics

      call defineDiagnostics

c Create draw*.in file

      if (isProcessor) then
        ndiag = count (sel_diag /= 0)

        call createDrawInGfile(ndiag,ngraph,(/lineplotfile/)
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

c Define GENERIC diagnostics

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
      diag_desc(neqd+16)= 'RMS |J|'
      diag_desc(neqd+17)= 'Max Vx'
      diag_desc(neqd+18)= 'Max Vy'
      diag_desc(neqd+19)= 'Max Vz'

      diag_desc(neqd+20)= 'Ez'
      diag_desc(neqd+21)= 'Reconnected flux'

      diag_desc(neqd+22)= 'delta/2'
      diag_desc(neqd+23)= 'w/2'

      diag_desc(neqd+24)= 'max(Bx) along y'
      diag_desc(neqd+25)= 'max(By) along x'

      diag_desc(neqd+26)= 'Max Vex'
      diag_desc(neqd+27)= 'Max Vey'
      diag_desc(neqd+28)= 'Max Vez'

      diag_desc(neqd+29:ngraph) = ''

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

      integer :: iigx,iigy,iigz
      logical :: init

      type(var_array),pointer :: varray

c Local variables

      integer :: i,j,k,ig,jg,kg,ieq,igx,igy,igz,nx,ny,nz,pos(3)
      real(8) :: array (ilom:ihip,jlom:jhip,klom:khip)
     .          ,vv    (ilom:ihip,jlom:jhip,klom:khip,3)

      real(8) :: energy,Vflux,Bflux,diverB,x1,y1,z1,dpert(neqd)
     .          ,told=0d0,ez,lrflux,rflux,rfluxn,rmsj,vmax(3),vemax(3)
     .          ,lvmax(3),ystar,lapjzm,lapjzp,dvm,dvp
     .          ,dl,ww=0d0,Bup=0d0,Bup0=1d0,Bdown=0d0
     .          ,x_xpt,y_xpt,z_xpt,x,y,z,xm,ym,zm,vi_out,ve_out

      logical :: cartsn

      save    :: rfluxn,told

      INTERFACE
        subroutine postProcessSolution(igrid,varray)
        use var_setup
        integer    :: igrid
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

c Code for electric field BC

cc      save_tn = .false.
      save_tn = .true.

      alpha = 1.-cnfactor

c Postprocess solution

      call postProcessSolution(igx,varray)

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

      array = rho*tmp

      if (gamma /= 1d0) then
        Et = 1d0/(gamma-1)
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

      Tflux = integral(igx,igy,igz,nx,ny,nz,bcnv(:,:,:,3)
     .                ,average=.false.,logical_dom=.true.)
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

c Maximum ION velocity diagnostics

      vmax = (/ maxval(vx_car),maxval(vy_car),maxval(vz_car) /)

#if defined(petsc)      
      lvmax = vmax
      call MPI_Reduce(lvmax,vmax,3,MPI_DOUBLE_PRECISION
     .               ,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#endif

c Maximum ELECTRON velocity diagnostics

      vemax = (/ maxval(vex),maxval(vey),maxval(vez) /)

#if defined(petsc)      
      lvmax = vemax
      call MPI_Reduce(lvmax,vemax,3,MPI_DOUBLE_PRECISION
     .               ,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#endif


c Problem specific diagnostics  <<<<<<<<<<<<<<<<<<<<<<<<<<<

c Parallel electric field and reconnected flux diagnostic (assumes periodic direction in Y)

      select case(trim(equil))
      case('dtm','gem')

#if defined(petsc)
        if (MPI_COMM_DIAG /= MPI_COMM_NULL) then
          lrflux = minval(Pflux(iplt,:,kplt))
          call MPI_Allreduce(lrflux,rflux,1,MPI_DOUBLE_PRECISION
     .                      ,MPI_MIN,MPI_COMM_DIAG,mpierr)

          lrflux = maxval(Pflux(iplt,:,kplt))-rflux
          call MPI_Allreduce(lrflux,rflux,1,MPI_DOUBLE_PRECISION
     .                      ,MPI_MAX,MPI_COMM_DIAG,mpierr)
        endif
#else
        rflux = maxval(Pflux(iplt,:,kplt))-minval(Pflux(iplt,:,kplt))
#endif

        if (time == told) then
          Ez = 0d0
        else
          Ez = (rflux-rfluxn)/(time-told)
        endif

        told = time
        rfluxn = rflux

      case('ic-hk','ic')

#if defined(petsc)
        if (isProcessor) then
          lrflux = maxval(Pflux(iplt,:,kplt))-Pflux(iplt,jplt,kplt)
        else
          lrflux = 0d0
        endif

        call MPI_Allreduce(lrflux,rflux,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_MAX,MPI_COMM_WORLD,mpierr)
#else
        rflux = maxval(Pflux(iplt,:,kplt))-Pflux(iplt,jplt,kplt)
#endif

        if (time == told) then
          Ez = 0d0
        else
          Ez = abs((rflux-rfluxn))/(time-told)
        endif

        told = time
        rfluxn = rflux

      case default

        Ez = 0d0
        rflux = 0d0

      end select

c Reconnection region diagnostics (assumes periodic direction in Y and X-point is at (iplt,jplt))

      select case(trim(equil))
      case('ic-hk','ic')!,'dtm','gem')

        i = iplt
        j = jplt
        k = 1

c       Current sheet thickness

        call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                              ,x_xpt,y_xpt,z_xpt)

        !Find mesh point where lapjz changes sign
        j = jplt+1
        lapjzm =lap(i,j,k,nx,ny,nz,igx,igy,igz,jz_car(:,:,:))
        do j = jplt+2,ny
          lapjzp = lap(i,j,k,nx,ny,nz,igx,igy,igz,jz_car(:,:,:))

          !Check whether sign changes
          if (lapjzp*lapjzm < 0d0) exit
          lapjzm = lapjzp
        enddo

        !find y*
        call getCartesianCoordinates(i,j  ,k,igx,igy,igz,ig,jg,kg
     .                              ,x,y ,z)
        call getCartesianCoordinates(i,j-1,k,igx,igy,igz,ig,jg,kg
     .                              ,x,ym,z)

        if (lapjzp /= lapjzm) then
          ystar = y - (y-ym)/(lapjzp-lapjzm)*lapjzm
        else
          ystar = 0d0
        endif

        !Compute delta/2 (half-width half-maximum)
        dl = 1.1768*ystar

ccc       Upstream magnetic field (at delta/2)
cc
cc        do j = jplt,ny
cc          call getCartesianCoordinates(i,j  ,k,igx,igy,igz,ig,jg,kg
cc     .                                ,x,y,z)
cc          call getCartesianCoordinates(i,j-1,k,igx,igy,igz,ig,jg,kg
cc     .                                ,x,ym,z)
cc          !Check whether sign of (y-dl) changes
cc          if((ym-dl)*(y-dl)<0d0) exit
cc        enddo
cc
cc        Bup = bx_car(iplt,j-1,kplt)
cc     .      +(bx_car(iplt,j  ,kplt)-bx_car(iplt,j-1,kplt))
cc     .      /(y-ym)*(dl-ym)

c       Maximum ION in-plane flows

        vmax(1) = maxval(abs(vx_car(:,jplt,kplt)))
        vmax(2) = maxval(abs(vy_car(iplt,:,kplt)))

c       Maximum ION in-plane flows

        vemax(1) = maxval(abs(vex(:,jplt,kplt)))
        vemax(2) = maxval(abs(vey(iplt,:,kplt)))

c       Maximum upstream magnetic field

cc        j = jplt+1
cc        dvm = bx_car(i,j,k)-bx_car(i,j-1,k)
cc        do j = jplt+2,ny
cc          dvp = bx_car(i,j,k)-bx_car(i,j-1,k)
cc          !Check whether sign of dbx/dy changes
cc          if(dvp*dvm<0d0) exit
cc          dvm = dvp
cc        enddo
cc
cc        Bup = bx_car(i,j,k)

        Bup = maxval(abs(bx_car(iplt,:,kplt)))

        if (itime == 0 .and. trim(equil) =='ic-hk') Bup0 = Bup

cc        call getCartesianCoordinates(i,j  ,k,igx,igy,igz,ig,jg,kg
cc     .                              ,x,y,z)
cc        call getCartesianCoordinates(i,j-1,k,igx,igy,igz,ig,jg,kg
cc     .                              ,x,ym,z)
cc
cc        Bup = bx(i,j-1,k)+(bx(i,j,k)-bx(i,j-1,k))/(y-ym)*(dl-ym)

c       Current sheet length  (assumes periodic direction in Y and X-point is at (iplt,jplt))

        j = jplt
        k = kplt

        !Find mesh point where X-component of Ve reaches extrema
        pos(1:1) = maxloc(abs(vex(:,jplt,kplt))) ; i = pos(1)

        call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,ww,y,z)

        ww = abs(ww)

        Bdown = by_car(i,jplt,kplt)

cc        i = iplt+1
cc        ve_out = 0d0
cc        dvm = vex(i  ,j,k)-vex(i-1,j,k)
cc        do i = iplt+2,nx
cc          dvp = vex(i+1,j,k)-vex(i  ,j,k)
cc
cc          !Check whether sign changes and assing current length
cc          if (dvp*dvm < 0d0.and. abs(vex(i,j,k)) > ve_out) then
cc            ve_out=max(abs(vex(i,j,k)),ve_out)
cc            call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,ww,y,z)
cc          endif
cc
cc          dvm = dvp
cc        enddo

cc        !Ion outflow
cc        i = iplt+1
cc        vi_out = 0d0
cc        dvm = vx_car(i  ,j,k)-vx_car(i-1,j,k)
cc        do i = iplt+2,nx
cc          dvp = vx_car(i+1,j,k)-vx_car(i  ,j,k)
cc
cc          !Check whether sign changes
cc          if (dvp*dvm < 0d0) vi_out=max(abs(vx_car(i,j,k)),vi_out)
cc          dvm = dvp
cc        enddo
cc
cc        !find w (half length)
cc        call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,ww,y,z)

c       Adjust for X-pt position

        dl = dl - y_xpt
        ww = ww - x_xpt

      case default

        dl = 0d0
        ww = 0d0
        Bup = 0d0
        Bdown = 0d0
        Bup0 = 1d0
        vi_out = 0d0
        ve_out = 0d0

      end select

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
      diagnostics(neqd+16)= rmsj
      diagnostics(neqd+17:neqd+19)= vmax

      diagnostics(neqd+20)= Ez/Bup0**2
      diagnostics(neqd+21)= rflux
      diagnostics(neqd+22) = dl
      diagnostics(neqd+23) = ww

      diagnostics(neqd+24) = Bup/Bup0
      diagnostics(neqd+25) = Bdown/Bup0
      diagnostics(neqd+26:neqd+28) = vemax

      !Compensate for symmetry BCs
      if (bcond(1) == FSYM .or. bcond(2) == FSYM) then
        diagnostics(neqd+1:neqd+3) = 2*diagnostics(neqd+1:neqd+3) !To compensate for symmetry
      endif

      if (bcond(3) == FSYM .or. bcond(4) == FSYM) then
        diagnostics(neqd+1:neqd+3) = 2*diagnostics(neqd+1:neqd+3) !To compensate for symmetry
      endif

      if (bcond(5) == FSYM .or. bcond(6) == FSYM) then
        diagnostics(neqd+1:neqd+3) = 2*diagnostics(neqd+1:neqd+3) !To compensate for symmetry
      endif

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
