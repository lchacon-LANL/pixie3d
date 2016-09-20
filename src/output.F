c output
c######################################################################
      subroutine output(varr)

c----------------------------------------------------------------------
c     Writes program output to standard output
c----------------------------------------------------------------------

      use equilibrium

      use timeStepping

      use precond_setup

      use app_iosetup

      use grid

      use nlfunction_setup

      implicit none

c Call variables

      type(var_array),pointer :: varr

c Local variables

      integer       :: ngrd(3),nx,ny,nz

      character(10) :: advective,vadvective,tsmeth

      real(8)       :: vmax(3),lvmax(3)

c Begin program

      ngrd = (/ gv%gparams%ngrdx,gv%gparams%ngrdy,gv%gparams%ngrdz /)

      nx = gv%gparams%nxv(1)
      ny = gv%gparams%nyv(1)
      nz = gv%gparams%nzv(1)

c Find maximum velocity and magnetic field components (cartesian)

      vmax(1) = maxval(abs(vcnv(1:nx,1:ny,1:nz,1)))
      vmax(2) = maxval(abs(vcnv(1:nx,1:ny,1:nz,2)))
      vmax(3) = maxval(abs(vcnv(1:nx,1:ny,1:nz,3)))

#if defined(petsc)
      lvmax = vmax
      call MPI_Reduce(lvmax,vmax,3,MPI_DOUBLE_PRECISION
     .               ,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
      where(vmax < 1d-20) vmax = 0d0
#endif

c Output

      if (my_rank == 0) then

        if (itime.lt.inewtime) then

          if (.not.test) call warnings

          select case(advect)
          case (1)
            advective = ' UPWIND'
          case (2)
            advective = '    ZIP'
          case (3)
            advective = '  QUICK'
          case (4)
            advective = '  SMART'
          case (5)
            advective = ' SSMART'
          case (6)
            advective = ' HORDER'
          case (7)
            advective = '  Gamma'
cc            advective = '   TVD'
          end select

          select case(v_advect)
          case (1)
            vadvective = ' UPWIND'
          case (2)
            vadvective = '    ZIP'
          case (3)
            vadvective = '  QUICK'
          case (4)
            vadvective = '  SMART'
          case (5)
            vadvective = ' SSMART'
          case (6)
            vadvective = ' HORDER'
          case (7)
            vadvective = '  Gamma'
cc            vadvective = '   TVD'
          end select

          if (cnfactor <= 5d-1) then
            if     (sm_flag.eq.0 .and. cnfactor > 0d0) then
              tsmeth = ' Theta'
            elseif (sm_flag.eq.0 .and. cnfactor == 0d0) then
              tsmeth = ' BDF1 (BE)'
            elseif (sm_flag.eq.1) then
              tsmeth = ' Rannacher'
            elseif (sm_flag.eq.2) then
              tsmeth = ' BDF2'
            endif
          elseif (cnfactor <= 1d0) then
            tsmeth = ' Explicit'
          else
            call pstop('output','cnfactor not within allowed limits')
          endif

          write (*,10)nxd,nyd,nzd,ngrd,pc_type,ts_adapt_dt,source
     .               ,numerical_grid,advective,vadvective,tsmeth

          write (*,*) 
          if (cnfactor > 5d-1) then
            write (*,105)
          else
            write (*,100)
          endif

          if (ts_restart) write (*,110) itime,time,dt,vmax

        else

          if (ilevel.ge.1) then
            write (*,*) 
            if (cnfactor == 1d0) then
              write (*,105)
            else
              write (*,100)
            endif
          endif

          if (cnfactor > 5d-1) then
            write(*,115) itime,time,dt,vmax
          else
            write(*,110) itime,time,dt,vmax,itgmres,itnewt,cnfactor
#if !defined(vec_pot)
     .                  ,it_divcl
            it_divcl = 0
#endif
          endif

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
     .        /,'  Scalar advection scheme...',a
     .        /,'  Velocity advection scheme.',a
     .        /,'  TS scheme.................',a)
 100  format ('   itime    time     dt     vx_max   vy_max   vz_max',
#if defined(vec_pot)
     .        '  GMRES  Newton  CN factor')
#else
     .        '  GMRES  Newton  CN factor  GMRES-DC')
#endif
 105  format ('   itime    time     dt     vx_max   vy_max   vz_max')

#if defined(vec_pot)
 110  format (i7,f9.3,1x,1p,e8.1,1x,3e9.1,0p,2i6,3x,f7.3)
#else
 110  format (i7,f9.3,1x,1p,e8.1,1x,3e9.1,0p,2i6,3x,f7.3,3x,i6)
#endif

 115  format (i7,f9.3,1x,1p,e8.1,1x,3e9.1)

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

      write (release,'(i1,a1,f3.1)') REL1,'.',REL2

      write (*,*)
      write (*,*) '****************************************************'
      write (*,*) '******  PIXIE3D Extended MHD code, V. '// 
     $                                       trim(release) //'   ******'
      write (*,*) '******                                        ******'
      write (*,*) '******         Developer: L. Chacon           ******'
      write (*,*) '****** Los Alamos National Laboratory (12-)   ******'
      write (*,*) '****** Oak Ridge  National Laboratory (08-12) ******'
      write (*,*) '****** Los Alamos National Laboratory (04-08) ******'
      write (*,*) '******             LA-CC 07-005               ******'
      write (*,*) '******                                        ******'
      write (*,*) '******        with contributions from:        ******'
      write (*,*) '******                                        ******'
      write (*,*) '******       Daniele Bonfiglio, CNR-RFX       ******'
      write (*,*) '******    Mark Berrill, Bobby Philip, ORNL    ******'
      write (*,*) '****************************************************'

c     Warnings

      write (*,*)
      write (*,*) '****************** NOTICES ******************'

#if defined(vec_pot) 
      write (*,*) 'Using VECTOR POTENTIAL formulation'
#else
      write (*,*) 'Using MAGNETIC FIELD formulation'
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
     .     write (*,*) 'Not solving for density'
      if (sym_st)
     .     write (*,*) 'Using symmetric viscous stress tensor in EOM'
      if (ion_hall .and. (di > 0d0))
     .     write (*,*) 'Solving ion-based Hall MHD formulation'
      if (slava) write (*,*) 'Using Slava Lukins Hall MHD formulation'
      if (fake_ve) write (*,*) 'Using fake Ve in Pi_e'
      if (adiabatic) write (*,*) 'Using polytropic EOS'
cc      if (use_p_eom) write (*,*) 'Solving for pressure instead of temp'
      if (gamma == 1d0) write (*,*) 'Using isothermal EOS (gamma = 1)'
      if (spitzer) write (*,*) 'Using Spitzer resistivity'
#if !defined(vec_pot)
      if (.not.pc_divclean) write (*,*) 'Div-cleaning postprocessed'
#endif
      write (*,*) '****************** NOTICES ******************'


      if (advect == 5) then
        write (*,*)
     .   'WARNING: Smooth smart does not work for non-uniform meshes'
      endif

      end subroutine warnings

      end subroutine
