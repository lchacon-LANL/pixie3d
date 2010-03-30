c output
c######################################################################
      subroutine output(varr)

c----------------------------------------------------------------------
c     Writes program output to standard output
c----------------------------------------------------------------------

      use equilibrium

      use timeStepping

      use precond_setup

      use counters

      use app_iosetup

      use grid

      use nlfunction_setup

      implicit none

c Call variables

      type(var_array),pointer :: varr

c Local variables

      integer       :: ngrd(3)

      character(10) :: advective,ts

      real(8)       :: vmax(3),lvmax(3)

c Begin program

      ngrd = (/ grid_params%ngrdx,grid_params%ngrdy,grid_params%ngrdz /)

c Find maximum velocity and magnetic field components (cartesian)

      vx_max = maxval(abs(varr%array_var(IVX)%array))
      vy_max = maxval(abs(varr%array_var(IVY)%array))
      vz_max = maxval(abs(varr%array_var(IVZ)%array))

      vmax = (/ vx_max,vy_max,vz_max /)

#if defined(petsc)
      lvmax = vmax
      call MPI_Reduce(lvmax,vmax,3,MPI_DOUBLE_PRECISION
     .               ,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#endif

cc      bx_max = maxval(abs(bx_car))
cc      by_max = maxval(abs(by_car))
cc      bz_max = maxval(abs(bz_car))

c Output

      if (my_rank == 0) then

        if (itime.lt.inewtime) then

          if (.not.test) call warnings

          if     (advect.eq.1 .or. advect == 7) then
            advective = 'UPWIND'
          elseif (advect.eq.2) then
            advective = '   ZIP'
          elseif (advect.eq.3) then
            advective = ' QUICK'
          elseif (advect.eq.4) then
            advective = ' SMART'
          elseif (advect.eq.5) then
            advective = 'SSMART'
          elseif (advect.eq.6) then
            advective = 'HORDER'
cc          elseif (advect.eq.7) then
cc            advective = '   TVD'
cc          else
cc            advective = 'PBRKCT'
          endif

          if (cnfactor < 1d0) then
            if     (sm_flag.eq.0 .and. cnfactor > 0d0) then
              ts = ' Theta'
            elseif (sm_flag.eq.0 .and. cnfactor == 0d0) then
              ts = ' BDF1 (BE)'
            elseif (sm_flag.eq.1) then
              ts = ' Rannacher'
            elseif (sm_flag.eq.2) then
              ts = ' BDF2'
            endif
          elseif (cnfactor == 1d0) then
            ts = ' Explicit'
          else
            call pstop('output','cnfactor not within allowed limits')
          endif

          write (*,10)nxd,nyd,nzd,ngrd,pc_type,timecorr,source
     .               ,numerical_grid,advective,ts

          write (*,*) 
          write (*,100)

          if (restart) write (*,110) itime,time,dt,vmax

        else

          if (ilevel.ge.1) then
            write (*,*) 
            write (*,100)
          endif

          write(*,110) itime,time,dt,vmax,itgmres,itnewt,cnfactor

        endif

        call flush(6)

      endif

c End program

 10   format (/,'  Grid mesh:................',i5,'x',i5,'x',i5
     .        /,'  Number of grid levels:....',i5,',',i5,',',i5
     .        /,'  Preconditioning method:...',a4,
     .        /,'  Time correction...........',l4,
     .        /,'  External source...........',l4,
     .        /,'  Numerical grid............',l4
     .        /,'  Advection scheme..........',a
     .        /,'  TS scheme.................',a)
 100  format ('   itime    time     dt     vx_max   vy_max   vz_max',
     .        '  GMRES  Newton  CN factor')

cc#if defined(RFX)
 110  format (i7,f9.3,1x,1p,e8.1,1x,3e9.1,0p,2i6,3x,f7.3)
cc#else
cc 110  format (i7,f9.3,1x,1p,e8.1,1x,0p,3f7.3,2i6,3x,f7.3)
cc#endif
      contains

c     #################################################################
      subroutine warnings

c     -----------------------------------------------------------------
c     Writes initial message and pertinent warings about input variables
c     at the beginning of simulation
c     -----------------------------------------------------------------

      implicit none

      character (5) :: release

c     Initial message

      write (release,'(i1,a1,i1)') REL1,'.',REL2

      write (*,*)
      write (*,*) '***********************************************'
      write (*,*) '****** PIXIE3D Extended MHD code, V. '//
     $            trim(release) //' ******'
      write (*,*) '******       Developer: L. Chacon        ******'
      write (*,*) '******  Los Alamos National Laboratory   ******'
      write (*,*) '******       Los Alamos, NM 87545        ******'
      write (*,*) '******           LA-CC 07-005            ******'
      write (*,*) '***********************************************'

c     Warnings

      write (*,*)
      write (*,*) '****************** NOTICES ******************'

#if defined(vec_pot) 
      write (*,*) 'Using VECTOR POTENTIAL version'
#else
      write (*,*) 'Using MAGNETIC FIELD version'
#endif

      if (nc_eom_jxb)
     .     write (*,*) 'Using non-conservative JxB     part of EOM'
      if (nc_eom_gp)
     .     write (*,*) 'Using non-conservative grad(p) part of EOM'
      if (nc_eom_v)
     .     write (*,*) 'Using non-conservative div(vv) part of EOM'
      if (.not.solenoidal)
     .     write (*,*) 'Using non-solenoidal version'
      if (.not.solve_rho)
     .     write (*,*) 'Not solving for density: rho = 1'
      if (sym_st)
     .     write (*,*) 'Using symmetric viscous stress tensor in EOM'
      if (ion_hall)
     .     write (*,*) 'Solving ion-based Hall MHD formulation'
      if (adiabatic)    write (*,*) 'Using polytropic EOS'
      if (gamma == 1d0) write (*,*) 'Using isothermal EOS (gamma = 1)'

      write (*,*) '****************** NOTICES ******************'

      end subroutine warnings

      end subroutine
