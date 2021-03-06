c output
c######################################################################
      subroutine output(varr)

c----------------------------------------------------------------------
c     Writes program output to standard output
c----------------------------------------------------------------------

      use equilibrium

      use timeStepping

      use PC_setup

      use app_iosetup

      use grid

      use nlfunction_setup

      implicit none

c Call variables

      type(var_array),pointer :: varr

c Local variables

      integer       :: ngrd(3),nx,ny,nz

      character(10) :: advective,vadvective,tsmeth,pxo
      character(40) :: pc
      
      real(8)       :: vmax(3),lvmax(3),dt__exp,ncfl

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

c Compute explicit time step

      if (itime >= inewtime.or.ts_restart) then
        call set_explicit_dt(dt__exp)
      endif

c Output

      if (my_rank == 0) then

        if (itime < inewtime) then

          if (.not.test) call warnings

          select case(advect)
          case (1)
            advective = ' UPWIND'
          case (2)
            advective = ' ZIP   '
          case (3)
            advective = ' QUICK '
          case (4)
            advective = ' SMART '
          case (5)
            advective = ' SSMART'
          case (6)
            advective = ' HORDER'
          case (7)
            advective = ' Gamma '
cc            advective = '   TVD'
          end select

          select case(v_advect)
          case (1)
            vadvective = ' UPWIND'
          case (2)
            vadvective = ' ZIP   '
          case (3)
            vadvective = ' QUICK '
          case (4)
            vadvective = ' SMART '
          case (5)
            vadvective = ' SSMART'
          case (6)
            vadvective = ' HORDER'
          case (7)
            vadvective = ' Gamma '
cc            vadvective = '   TVD'
          end select

          if (cnfactor <= 5d-1) then
            if     (sm_flag.eq.0 .and. cnfactor > 0d0) then
              tsmeth = ' Theta'
            elseif (sm_flag.eq.0 .and. cnfactor == 0d0) then
              tsmeth = ' BDF1(BE)'
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

          select case(trim(pc_type))
          case('s1')
            pc = ' SI PC w/o flow'
          case('s2')
            pc = ' SI PC with flow, iterative'
          case('s3')
            pc = ' SI PC with flow, commutator'
          case('s4')
            pc = ' EMHD split PC'
          case ('id')
            pc = ' Identity'
          case default
            pc = ' Undefined!'
          end select

          select case(xport_order)
          case(2)
            pxo = ' 2nd      '
          case(4)
            pxo = ' 4th-7pt  '
          case(-4)
            pxo = ' 4th-5pt  '
          end select
          
#if defined(petsc) && defined(coarse_MG)
          write (*,10)nxd,nyd,nzd,npx,npy,npz,trim(pc),ts_adapt_dt
     .           ,source,numerical_grid,advective,vadvective,pxo,tsmeth
#else
          write (*,10)nxd,nyd,nzd,ngrd,trim(pc),ts_adapt_dt
     .           ,source,numerical_grid,advective,vadvective,pxo,tsmeth
#endif
          write (*,*) 
          if (cnfactor > 5d-1) then
            write (*,105)
          else
            write (*,100)
          endif

          if (ts_restart) write (*,110) itime,time,dt,0d0,vmax

        else

          ncfl = dt/dt__exp

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
            write(*,110) itime,time,dt,ncfl,vmax,itgmres,itnewt,cnfactor
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
#if defined(petsc) && defined(coarse_MG)
     .        /,'  Number of processors:.....',i5,',',i5,',',i5
#else
     .        /,'  Number of grid levels:....',i5,',',i5,',',i5
#endif
     .        /,'  Preconditioning method:...',a,
     .        /,'  Timestep adaption.........',l2,
     .        /,'  External source...........',l2,
     .        /,'  Numerical grid............',l2
     .        /,'  Scalar advection scheme...',a
     .        /,'  Velocity advection scheme.',a
     .        /,'  Anisotropic transport sch.',a
     .        /,'  TS scheme.................',a)
 100  format ('    itime      time     dt     #CFL   vx_max   vy_max',
#if defined(vec_pot)
     .        '   vz_max    GMRES Newton CN factor')
#else
     .        '   vz_max    GMRES Newton CN factor  GMRES-DC')
#endif
 105  format ('    itime      time     dt     vx_max   vy_max   vz_max')

#if defined(vec_pot)
 110  format (i10,f10.2,1x,1p,2e8.1,3e9.1,0p,2i6,3x,f7.3)
#else
 110  format (i10,f10.2,1x,1p,2e8.1,3e9.1,0p,2i6,3x,f7.3,i9)
#endif

 115  format (i10,f10.2,1x,1p,e8.1,1x,3e9.1)

      contains

c     #################################################################
      subroutine warnings

c     -----------------------------------------------------------------
c     Writes initial message and pertinent warings about input variables
c     at the beginning of simulation
c     -----------------------------------------------------------------

      implicit none

      character (50) :: release

c     Initial message

c$$$      write (release,'(a)') REL1
c$$$
c$$$      write (*,*)
c$$$      write (*,*) '****************************************************'
c$$$      write (*,*) '******       PIXIE3D Extended MHD code        ******'
c$$$      write (*,*) '******         V. '//release(1:28)//         '******'
c$$$      write (*,*) '******                                        ******'
c$$$      write (*,*) '******         Developer: L. Chacon           ******'
c$$$      write (*,*) '****** Los Alamos National Laboratory (12-)   ******'
c$$$      write (*,*) '****** Oak Ridge  National Laboratory (08-12) ******'
c$$$      write (*,*) '****** Los Alamos National Laboratory (04-08) ******'
c$$$      write (*,*) '******             LA-CC 07-005               ******'
c$$$      write (*,*) '******                                        ******'
c$$$      write (*,*) '******        with contributions from:        ******'
c$$$      write (*,*) '******                                        ******'
c$$$      write (*,*) '******       Daniele Bonfiglio, CNR-RFX       ******'
c$$$      write (*,*) '******    Mark Berrill, Bobby Philip, ORNL    ******'
c$$$      write (*,*) '****************************************************'

c     Warnings

      write (*,*)
      write (*,*) '****************** CONFIG SUMMARY ******************'

#if defined(vec_pot) 
      write (*,*) 'Using VECTOR POTENTIAL formulation'
#else
      write (*,*) 'Using MAGNETIC FIELD formulation'
#endif

      if (limit_dt_flow_cfl)
     .     write (*,*) 'Limiting timestep to flow CFL'
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
      if (solve_prs) then
        write (*,*) 'Solving for P=Pi+Pe'
      else
        write (*,*) 'Solving for Te'
      endif

      if (adiabatic_rho) write (*,*) 'Using polytropic EOS for density'
      if (adiabatic) write (*,*) 'Using adiabatic P,T update'
      if (sym_st)
     .     write (*,*) 'Using symmetric viscous stress tensor in EOM'
      if (ion_hall)
     .     write (*,*) 'Solving ion-based Hall MHD formulation'
      if (slava) write (*,*) 'Using Slava Lukins Hall MHD formulation'
      if (gamma == 1d0) write (*,*) 'Using isothermal EOS (gamma = 1)'
      if (spitzer) write (*,*) 'Using Spitzer resistivity'
      if (chi_bgnskii) write (*,*) 'Using Braginskii transport coeffs'

      if (use_p_eom)
     .     write (*,*) 'Using grad(p) instead of grad(nT) in EOM'
      if (use_p_diff)
     .     write (*,*) 'Using lap(p) instead of lap(T) in pressure eq.'
      if (use_p_classic_eq) then
        if (solve_prs) then
          write (*,*) 'Using p.div(v) form of pressure eq.'
        else
          write (*,*) 'Using (gamma-1)*Te.div(v) form of Te eq.'
        endif
      else
        if (solve_prs) then
          write (*,*) 'Using -v.grad(p) form of pressure eq.'
        else
          write (*,*) 'Using (2-gamma)*v.grad(Te) form of Te eq.'
        endif
      endif
#if !defined(vec_pot)
      if (.not.pc_divclean) write (*,*) 'Div-cleaning postprocessed'
#endif
      if (refine_eq) write (*,*) 'Running Equilibrium-checking mode'
      write (*,*) '****************************************************'

      if (advect == 5) then
        write (*,*)
     .   'WARNING: Smooth smart does not work for non-uniform meshes'
      endif

      end subroutine warnings

      end subroutine
