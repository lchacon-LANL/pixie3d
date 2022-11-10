
c readGraphicsInput
c######################################################################
      subroutine readGraphicsInput

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use app_iosetup

      use graphics_io

      use equilibrium

      use acc_tests

      implicit none

c Call variables

c Local variables

      integer :: i,ioerr

      !For compatibility
      namelist /graphdef/ sel_diag,sel_graph,ndplot,dplot
     .                   ,hdf_plot,hdf_debug,car_diag_plots
     .                   ,prof_conf,cont_conf,clean
     .                   ,iplot,jplot,kplot,poincare,xdraw_plot
     .                   ,adios_plot,temp_acc,acc_ref_file
     .                   ,heat_flx_diag

      namelist /plots/    clean,poincare

c Begin program

c Graphics defaults

      ndplot = 0            ! Postprocessing interval (# of time steps; integer)
      dplot  = 0d0          !       "          "      (time interval; real)
      hdf_plot =.false.     ! Whether an HDF5 file is to be created
      hdf_debug=.false.     ! Debugging option in HDF5
      adios_plot =.false.   ! Whether an ADIOS-BP file is to be created
      xdraw_plot = .true.   ! Whether to dump XDRAW files
      car_diag_plots=.true. ! Whether diagnostic plots are in Cartesian geom or not
      heat_flx_diag =.false.! Whether to output heat-flux boundary plots
      
      !Local (point) diagnostic configuration (XDRAW)
      iplot = 1
      jplot = 1
      kplot = 1

      !Line profile configuration (XDRAW)
      prof_conf%line   = 1          ! Direction (1 -> x, 2 -> y, 3 -> z)
      prof_conf%label  ='x'         ! Label
      prof_conf%coords = (/1,1,1/)  ! Line coordinates 

      !Contour slice configuration (XDRAW)
      cont_conf%plane  = 3          ! Normal to cut plane (1 -> x, 2 -> y, 3 -> z)
      cont_conf%label  = (/'x','y'/)! Contour plot axes labels
      cont_conf%coords = (/1,1,1/)  ! Plane coordinates along normal

      sel_diag = 0   !Array of size "xdraw_cont_lim" (currently set to 16)
                     !  indicating time histories to be dumped (see
                     !  drawgamma.in --generated after postprocessing--
                     !  or diagnostics.F for information on available variables)

      sel_graph = 0  !Selects graphics to be shown with XDRAW contour plotter.

c Read computation initializations (external)

      call readInput

c Read graphics initialization parameters

      open(unit=uinput,file=in_file,status='old')

      read(uinput,graphdef,iostat=ioerr)

      close(unit=uinput)

      if (ioerr /= 0) then
        open(unit=uinput,file=in_file,status='old')
        read(uinput,plots,iostat=ioerr)
        close(unit=uinput)
        if (ioerr /= 0.and.my_rank == 0) then
          write (*,*) 'Problem reading namelist "plots"'
        endif

        ioerr=readAccTestConfig (in_file)
        if (ioerr /= 0.and.my_rank == 0) then
          write (*,*) 'Problem reading namelist "acc_test"'
        endif

        ioerr=readGraphicsConfig(in_file)
        if (ioerr /= 0.and.my_rank == 0) then
          write (*,*) 'Problem reading namelist "graphics"'
        endif
      endif

c Consistency check

      if (iplot > nxd .or. jplot > nyd .or. kplot > nzd) then
        call pstop('readGraphicsInput'
     .            ,'Plot indices out of range')
      endif

c End program

      end subroutine readGraphicsInput
