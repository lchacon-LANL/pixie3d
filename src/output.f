c output
c######################################################################
      subroutine output

c----------------------------------------------------------------------
c     Writes program output to standard output
c----------------------------------------------------------------------

      use equilibrium

      use timeStepping

      use precond_setup

      use counters

      use iosetup

      use grid

      use nlfunction_setup

      use icond

      implicit none

c Call variables

c Local variables

      integer(4)  :: ngrd(3)

      character(10):: ts

c Begin program

      ngrd = (/ grid_params%ngrdx,grid_params%ngrdy,grid_params%ngrdz /)

c Find maximum velocity and magnetic field components (cartesian)

      vx_max = maxval(abs(vx))
      vy_max = maxval(abs(vy))
      vz_max = maxval(abs(vz))

cc      bx_max = maxval(abs(bx_car))
cc      by_max = maxval(abs(by_car))
cc      bz_max = maxval(abs(bz_car))

c Output

      if (itime.lt.inewtime) then

        call warnings

        if     (sm_flag.eq.0) then
          ts = ' Theta'
        elseif (sm_flag.eq.1) then
          ts = ' Rannacher'
        elseif (sm_flag.eq.2) then
          ts = ' BDF2'
        endif

cc        write (*,10) nxd,nyd,nzd,ngrd,precon,timecorr,source
        write (*,10) nxd,nyd,nzd,ngrd,timecorr,source,numerical_grid,ts

        write (*,*) 
        write (*,100)

        if (restart) write (*,110) itime,time,dt,vx_max,vy_max,vz_max

      else

        if (ilevel.ge.1) then
          write (*,*) 
          write (*,100)
        endif

        write(*,110) itime,time,dt,vx_max,vy_max,vz_max,itgmres,itnewt
     .              ,cnfactor

      endif

c End program

 10   format (/,'  Grid mesh:................',i4,'x',i3,'x',i3
     .        /,'  Number of grid levels:....',i4,',',i2,',',i2
cc     .        /,'  Preconditioning method:...',a11,
     .        /,'  Time correction...........',l4,
     .        /,'  External source...........',l4,
     .        /,'  Numerical grid............',l4
     .        /,'  TS scheme.................',a)
 100  format ('   itime    time     dt    vx_max vy_max vz_max GMRES',
     .        '  Newton  CN factor')
 110  format (i7,f9.2,1x,1p,e8.1,1x,0p,3f7.3,2i6,3x,f7.3)

      contains

c     #################################################################
      subroutine warnings

c     -----------------------------------------------------------------
c     Writes pertinent warings about input variables at the beginning
c     of simulation
c     -----------------------------------------------------------------

      implicit none

      if (nc_eom_f .or. nc_eom_v .or. (.not.solenoidal)) then
        write (*,*)
        write (*,*) '****************** NOTICE ******************'

        if (nc_eom_f)
     .       write (*,*) 'Using non-conservative EM    part of EOM'
        if (nc_eom_v)
     .       write (*,*) 'Using non-conservative hydro part of EOM'
        if (.not.solenoidal)
     .       write (*,*) 'Using non-solenoidal version'
      
        write (*,*) '****************** NOTICE ******************'
      endif

      end subroutine warnings

      end subroutine
