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

c Begin program

      ngrd = (/ grid_params%ngrdx,grid_params%ngrdy,grid_params%ngrdz /)

ccc Find maximum velocity and magnetic field components (cartesian)
cc
cc      vx_max = maxval(abs(vx_car))
cc      vy_max = maxval(abs(vy_car))
cc      vz_max = maxval(abs(vz_car))
cc
cc      bx_max = maxval(abs(bx_car))
cc      by_max = maxval(abs(by_car))
cc      bz_max = maxval(abs(bz_car))

c Output

      if (itime.lt.inewtime) then

cc        write (*,10) nxd,nyd,nzd,ngrd,precon,timecorr,source
        write (*,10) nxd,nyd,nzd,ngrd,timecorr,source

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
     .        /,'  External source...........',l4)
 100  format ('   itime    time     dt    vx_max vy_max vz_max GMRES',
     .        '  Newton  CN factor')
 110  format (i7,f9.2,1x,1pe8.1,1x,0p3f7.3,2i6,3x,f7.3)

      end subroutine
