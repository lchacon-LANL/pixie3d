c setGraphicsGroups
c####################################################################
      subroutine setGraphicsGroups

c--------------------------------------------------------------------
c     Set number of graphics groups
c--------------------------------------------------------------------

      use graphics

      implicit none

c Call variables

c Local variables

c Begin program

      ngroups = 5   !Number of plotting groups

c End program

      end subroutine setGraphicsGroups

c defineGraphicsGroups
c####################################################################
      subroutine defineGraphicsGroups

c--------------------------------------------------------------------
c     Set graphics files and dumping intervals
c--------------------------------------------------------------------

      use parameters

      use variables

      use graphics

      use auxiliaryVariables

      use auxPlotVariables

      use equilibrium

      implicit none

c Call variables

c Local variables

      integer    :: ieq,i

c Begin program

c Define I/O

      call defineGraphicsIO

c Define graphics group #1: Contravariant variables

      graph(1)%cartesian=.false.
      graph(1)%descr='Cnv_variables'

#if !defined(vec_pot)
      do ieq = 1,neqd
        graph(1)%array_graph(ieq)%array => u_np%array_var(ieq)%array
        graph(1)%array_graph(ieq)%descr =  u_np%array_var(ieq)%descr
      enddo

      graph(1)%array_graph(neqd+1)%array => jx
      graph(1)%array_graph(neqd+1)%descr = 'J^1'

      graph(1)%array_graph(neqd+2)%array => jy
      graph(1)%array_graph(neqd+2)%descr = 'J^2'

      graph(1)%array_graph(neqd+3)%array => jz
      graph(1)%array_graph(neqd+3)%descr = 'J^3'

      graph(1)%array_graph(neqd+4)%array => vx
      graph(1)%array_graph(neqd+4)%descr = 'V^1'

      graph(1)%array_graph(neqd+5)%array => vy
      graph(1)%array_graph(neqd+5)%descr = 'V^2'

      graph(1)%array_graph(neqd+6)%array => vz
      graph(1)%array_graph(neqd+6)%descr = 'V^3'

      graph(1)%array_graph(neqd+7)%array => vex
      graph(1)%array_graph(neqd+7)%descr = 'Ve_1'

      graph(1)%array_graph(neqd+8)%array => vey
      graph(1)%array_graph(neqd+8)%descr = 'Ve_2'

      graph(1)%array_graph(neqd+9)%array => vez
      graph(1)%array_graph(neqd+9)%descr = 'Ve_3'

      graph(1)%array_graph(neqd+10:ngraph)%descr = ''
#else
      graph(1)%array_graph(IRHO)%array => u_np%array_var(IRHO)%array
      graph(1)%array_graph(IRHO)%descr =  u_np%array_var(IRHO)%descr

      graph(1)%array_graph(ITMP)%array => u_np%array_var(ITMP)%array
      graph(1)%array_graph(ITMP)%descr =  u_np%array_var(ITMP)%descr

      graph(1)%array_graph(IVX)%array => u_np%array_var(IVX)%array
      graph(1)%array_graph(IVX)%descr =  u_np%array_var(IVX)%descr
                                                           
      graph(1)%array_graph(IVY)%array => u_np%array_var(IVY)%array
      graph(1)%array_graph(IVY)%descr =  u_np%array_var(IVY)%descr
                                                           
      graph(1)%array_graph(IVZ)%array => u_np%array_var(IVZ)%array
      graph(1)%array_graph(IVZ)%descr =  u_np%array_var(IVZ)%descr

      graph(1)%array_graph(IAX)%array => ax_cnv
      graph(1)%array_graph(IAX)%descr = 'A^1'

      graph(1)%array_graph(IAY)%array => ay_cnv
      graph(1)%array_graph(IAY)%descr = 'A^2'

      graph(1)%array_graph(IAZ)%array => az_cnv
      graph(1)%array_graph(IAZ)%descr = 'A^3'

      graph(1)%array_graph(neqd+1)%array => bx
      graph(1)%array_graph(neqd+1)%descr = 'B^1'
                                 
      graph(1)%array_graph(neqd+2)%array => by
      graph(1)%array_graph(neqd+2)%descr = 'B^2'
                                 
      graph(1)%array_graph(neqd+3)%array => bz
      graph(1)%array_graph(neqd+3)%descr = 'B^3'

      graph(1)%array_graph(neqd+4)%array => jx
      graph(1)%array_graph(neqd+4)%descr = 'J^1'

      graph(1)%array_graph(neqd+5)%array => jy
      graph(1)%array_graph(neqd+5)%descr = 'J^2'

      graph(1)%array_graph(neqd+6)%array => jz
      graph(1)%array_graph(neqd+6)%descr = 'J^3'

      graph(1)%array_graph(neqd+7)%array => vx
      graph(1)%array_graph(neqd+7)%descr = 'V^1'

      graph(1)%array_graph(neqd+8)%array => vy
      graph(1)%array_graph(neqd+8)%descr = 'V^2'

      graph(1)%array_graph(neqd+9)%array => vz
      graph(1)%array_graph(neqd+9)%descr = 'V^3'

      graph(1)%array_graph(neqd+10:ngraph)%descr = ''
#endif

cc      sel_gr(1,:) = sel_graph
      sel_gr(1,:) = (/ (i,i=1,xdraw_cont_lim) /)

      prof_ivar(1,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (1,:) = 0        !No log scales
      prof_spline(1) = .false.  !No splines

#if !defined(RFX)
c Define graphics group #2: Cov variables

      graph(2)%cartesian=.false.
      graph(2)%descr='Cov_variables'

#if !defined(vec_pot)
      graph(2)%array_graph(IRHO)%array => u_np%array_var(IRHO)%array
      graph(2)%array_graph(IRHO)%descr =  u_np%array_var(IRHO)%descr

      graph(2)%array_graph(ITMP)%array => u_np%array_var(ITMP)%array
      graph(2)%array_graph(ITMP)%descr =  u_np%array_var(ITMP)%descr

      graph(2)%array_graph(IVX)%array => vx_cov
      graph(2)%array_graph(IVX)%descr = 'V_1'

      graph(2)%array_graph(IVY)%array => vy_cov
      graph(2)%array_graph(IVY)%descr = 'V_2'

      graph(2)%array_graph(IVZ)%array => vz_cov
      graph(2)%array_graph(IVZ)%descr = 'V_3'

      graph(2)%array_graph(IBX)%array => bx_cov
      graph(2)%array_graph(IBX)%descr = 'B_1'

      graph(2)%array_graph(IBY)%array => by_cov
      graph(2)%array_graph(IBY)%descr = 'B_2'

      graph(2)%array_graph(IBZ)%array => bz_cov
      graph(2)%array_graph(IBZ)%descr = 'B_3'

      graph(2)%array_graph(neqd+1)%array => jx_cov
      graph(2)%array_graph(neqd+1)%descr = 'J_1'

      graph(2)%array_graph(neqd+2)%array => jy_cov
      graph(2)%array_graph(neqd+2)%descr = 'J_2'

      graph(2)%array_graph(neqd+3)%array => jz_cov
      graph(2)%array_graph(neqd+3)%descr = 'J_3'

      graph(2)%array_graph(neqd+4:ngraph)%descr = ''
#else
      graph(2)%array_graph(IRHO)%array => u_np%array_var(IRHO)%array
      graph(2)%array_graph(IRHO)%descr =  u_np%array_var(IRHO)%descr

      graph(2)%array_graph(ITMP)%array => u_np%array_var(ITMP)%array
      graph(2)%array_graph(ITMP)%descr =  u_np%array_var(ITMP)%descr

      graph(2)%array_graph(IVX)%array => vx_cov
      graph(2)%array_graph(IVX)%descr = 'V_1'

      graph(2)%array_graph(IVY)%array => vy_cov
      graph(2)%array_graph(IVY)%descr = 'V_2'

      graph(2)%array_graph(IVZ)%array => vz_cov
      graph(2)%array_graph(IVZ)%descr = 'V_3'

      graph(2)%array_graph(IAX)%array => u_np%array_var(IAX)%array
      graph(2)%array_graph(IAX)%descr =  u_np%array_var(IAX)%descr
                                                           
      graph(2)%array_graph(IAY)%array => u_np%array_var(IAY)%array
      graph(2)%array_graph(IAY)%descr =  u_np%array_var(IAY)%descr
                                                           
      graph(2)%array_graph(IAZ)%array => u_np%array_var(IAZ)%array
      graph(2)%array_graph(IAZ)%descr =  u_np%array_var(IAZ)%descr

      graph(2)%array_graph(neqd+1)%array => bx_cov
      graph(2)%array_graph(neqd+1)%descr = 'B_1'
                                 
      graph(2)%array_graph(neqd+2)%array => by_cov
      graph(2)%array_graph(neqd+2)%descr = 'B_2'
                                 
      graph(2)%array_graph(neqd+3)%array => bz_cov
      graph(2)%array_graph(neqd+3)%descr = 'B_3'

      graph(2)%array_graph(neqd+4)%array => jx_cov
      graph(2)%array_graph(neqd+4)%descr = 'J_1'

      graph(2)%array_graph(neqd+5)%array => jy_cov
      graph(2)%array_graph(neqd+5)%descr = 'J_2'

      graph(2)%array_graph(neqd+6)%array => jz_cov
      graph(2)%array_graph(neqd+6)%descr = 'J_3'

      graph(2)%array_graph(neqd+7:ngraph)%descr = ''
#endif

cc      sel_gr(2,:) = sel_graph
      sel_gr(2,:) = (/ (i,i=1,xdraw_cont_lim) /)

      prof_ivar(2,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (2,:) = 0        !No log scales
      prof_spline(2) = .false.  !No splines

#else

c Define graphics group #2: Cyl variables

      graph(2)%cartesian=.false.
      graph(2)%descr='Cyl_variables'

#if !defined(vec_pot)
      graph(2)%array_graph(1)%array => vx_cov
      graph(2)%array_graph(1)%descr = 'Vr'

      graph(2)%array_graph(2)%array => vy_cov
      graph(2)%array_graph(2)%descr = 'Vth'

      graph(2)%array_graph(3)%array => vz_cov
      graph(2)%array_graph(3)%descr = 'Vz'

      graph(2)%array_graph(4)%array => bx_cov
      graph(2)%array_graph(4)%descr = 'Br'

      graph(2)%array_graph(5)%array => by_cov
      graph(2)%array_graph(5)%descr = 'Bth'

      graph(2)%array_graph(6)%array => bz_cov
      graph(2)%array_graph(6)%descr = 'Bz'

      graph(2)%array_graph(7)%array => jx_cov
      graph(2)%array_graph(7)%descr = 'Jr'

      graph(2)%array_graph(8)%array => jy_cov
      graph(2)%array_graph(8)%descr = 'Jth'

      graph(2)%array_graph(9)%array => jz_cov
      graph(2)%array_graph(9)%descr = 'Jz'

      graph(2)%array_graph(10:ngraph)%descr = ''
#else
      graph(2)%array_graph(1)%array => vx_cov
      graph(2)%array_graph(1)%descr = 'Vr'

      graph(2)%array_graph(2)%array => vy_cov
      graph(2)%array_graph(2)%descr = 'Vth'

      graph(2)%array_graph(3)%array => vz_cov
      graph(2)%array_graph(3)%descr = 'Vz'

      graph(2)%array_graph(4)%array => ax_cnv
      graph(2)%array_graph(4)%descr =  'Ar'
                                                        
      graph(2)%array_graph(5)%array => ay_cnv
      graph(2)%array_graph(5)%descr =  'Ath'
                                                        
      graph(2)%array_graph(6)%array => az_cnv
      graph(2)%array_graph(6)%descr =  'Az'

      graph(2)%array_graph(7)%array => bx_cov
      graph(2)%array_graph(7)%descr = 'Br'
                           
      graph(2)%array_graph(8)%array => by_cov
      graph(2)%array_graph(8)%descr = 'Bth'
                           
      graph(2)%array_graph(9)%array => bz_cov
      graph(2)%array_graph(9)%descr = 'Bz'

      graph(2)%array_graph(10)%array => jx_cov
      graph(2)%array_graph(10)%descr = 'Jr'

      graph(2)%array_graph(11)%array => jy_cov
      graph(2)%array_graph(11)%descr = 'Jth'

      graph(2)%array_graph(12)%array => jz_cov
      graph(2)%array_graph(12)%descr = 'Jz'

      graph(2)%array_graph(13:ngraph)%descr = ''
#endif

cc      sel_gr(2,:) = sel_graph
      sel_gr(2,:) = (/ (i,i=1,xdraw_cont_lim) /)

      prof_ivar(2,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (2,:) = 0        !No log scales
      prof_spline(2) = .false.  !No splines

#endif

c Define graphics group #3: Cartesian variables

      graph(3)%cartesian=.true.
      graph(3)%descr='Car_variables'

#if !defined(vec_pot)
      graph(3)%array_graph(IRHO)%array => u_np%array_var(IRHO)%array
      graph(3)%array_graph(IRHO)%descr =  u_np%array_var(IRHO)%descr

      graph(3)%array_graph(ITMP)%array => u_np%array_var(ITMP)%array
      graph(3)%array_graph(ITMP)%descr =  u_np%array_var(ITMP)%descr

      graph(3)%array_graph(IVX)%array => vx_car
      graph(3)%array_graph(IVX)%descr = 'Vx'

      graph(3)%array_graph(IVY)%array => vy_car
      graph(3)%array_graph(IVY)%descr = 'Vy'

      graph(3)%array_graph(IVZ)%array => vz_car
      graph(3)%array_graph(IVZ)%descr = 'Vz'

      graph(3)%array_graph(IVX:IVZ)%vector_name ='V' !Identifies vector comp.

      graph(3)%array_graph(IBX)%array => bx_car
      graph(3)%array_graph(IBX)%descr = 'Bx'

      graph(3)%array_graph(IBY)%array => by_car
      graph(3)%array_graph(IBY)%descr = 'By'

      graph(3)%array_graph(IBZ)%array => bz_car
      graph(3)%array_graph(IBZ)%descr = 'Bz'

      graph(3)%array_graph(IBX:IBZ)%vector_name ='B' !Identifies vector comp.

      graph(3)%array_graph(neqd+1)%array => jx_car
      graph(3)%array_graph(neqd+1)%descr = 'Jx'

      graph(3)%array_graph(neqd+2)%array => jy_car
      graph(3)%array_graph(neqd+2)%descr = 'Jy'

      graph(3)%array_graph(neqd+3)%array => jz_car
      graph(3)%array_graph(neqd+3)%descr = 'Jz'

      graph(3)%array_graph(neqd+1:neqd+3)%vector_name = 'J' !Identifies vector comp.

      graph(3)%array_graph(neqd+4:ngraph)%descr = ''
#else
      graph(3)%array_graph(IRHO)%array => u_np%array_var(IRHO)%array
      graph(3)%array_graph(IRHO)%descr =  u_np%array_var(IRHO)%descr

      graph(3)%array_graph(ITMP)%array => u_np%array_var(ITMP)%array
      graph(3)%array_graph(ITMP)%descr =  u_np%array_var(ITMP)%descr

      graph(3)%array_graph(IVX)%array => vx_car
      graph(3)%array_graph(IVX)%descr = 'Vx'

      graph(3)%array_graph(IVY)%array => vy_car
      graph(3)%array_graph(IVY)%descr = 'Vy'

      graph(3)%array_graph(IVZ)%array => vz_car
      graph(3)%array_graph(IVZ)%descr = 'Vz'

      graph(3)%array_graph(IVX:IVZ)%vector_name ='V' !Identifies vector comp.

      graph(3)%array_graph(IAX)%array => ax_car
      graph(3)%array_graph(IAX)%descr = 'Ax'

      graph(3)%array_graph(IAY)%array => ay_car
      graph(3)%array_graph(IAY)%descr = 'Ay'

      graph(3)%array_graph(IAZ)%array => az_car
      graph(3)%array_graph(IAZ)%descr = 'Az'

      graph(3)%array_graph(neqd+1)%array => bx_car
      graph(3)%array_graph(neqd+1)%descr = 'Bx'
                                 
      graph(3)%array_graph(neqd+2)%array => by_car
      graph(3)%array_graph(neqd+2)%descr = 'By'
                                 
      graph(3)%array_graph(neqd+3)%array => bz_car
      graph(3)%array_graph(neqd+3)%descr = 'Bz'

      graph(3)%array_graph(neqd+1:neqd+3)%vector_name ='B' !Identifies vector comp.

      graph(3)%array_graph(neqd+4)%array => jx_car
      graph(3)%array_graph(neqd+4)%descr = 'Jx'

      graph(3)%array_graph(neqd+5)%array => jy_car
      graph(3)%array_graph(neqd+5)%descr = 'Jy'

      graph(3)%array_graph(neqd+6)%array => jz_car
      graph(3)%array_graph(neqd+6)%descr = 'Jz'

      graph(3)%array_graph(neqd+4:neqd+6)%vector_name = 'J' !Identifies vector comp.

      graph(3)%array_graph(neqd+7:ngraph)%descr = ''
#endif

cc      sel_gr(3,:) = sel_graph
      sel_gr(3,:) = (/ (i,i=1,xdraw_cont_lim) /)

      prof_ivar(3,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (3,:) = 0        !No log scales
      prof_spline(3) = .false.  !No splines

c Define graphics group #4: Diagnostics

      graph(4)%cartesian=.true.
cc      graph(4)%cartesian=.false.
      graph(4)%descr='Diagnostics'

      graph(4)%array_graph(1)%array => nuu
      graph(4)%array_graph(1)%descr = 'nu'

      graph(4)%array_graph(2)%array => eeta
      graph(4)%array_graph(2)%descr = 'eta'

      graph(4)%array_graph(3)%array => divrgB
      graph(4)%array_graph(3)%descr = 'local div(B)'

      graph(4)%array_graph(4)%array => divrgV
      graph(4)%array_graph(4)%descr = 'local div(V)'

      graph(4)%array_graph(5)%array => divrgJ
      graph(4)%array_graph(5)%descr = 'local div(J)'

      graph(4)%array_graph(6)%array => Pflux
      graph(4)%array_graph(6)%descr = 'Poloidal flux'

      graph(4)%array_graph(7)%array => qfactor
      graph(4)%array_graph(7)%descr = 'q factor'

cc      graph(4)%array_graph(7)%array => qfactr2
cc      graph(4)%array_graph(7)%descr = 'avrgd q factor'

      graph(4)%array_graph(8)%array => lambda
      graph(4)%array_graph(8)%descr = 'lambda'

cc      graph(4)%array_graph(8)%array => lambda2
cc      graph(4)%array_graph(8)%descr = 'avrgd lambda'

      graph(4)%array_graph(9)%array => p_tot
      graph(4)%array_graph(9)%descr = 'Pressure'

      graph(4)%array_graph(10)%array => vpar
      graph(4)%array_graph(10)%descr = 'V||'

      graph(4)%array_graph(11)%array => vpsi
      graph(4)%array_graph(11)%descr = 'V.grad(psi)'

      graph(4)%array_graph(12)%array => vperp
      graph(4)%array_graph(12)%descr = 'V.(b x grad(psi))'

      graph(4)%array_graph(13)%array => jpar
      graph(4)%array_graph(13)%descr = 'E_1'

      graph(4)%array_graph(14)%array => jpsi
      graph(4)%array_graph(14)%descr = 'E_2'

      graph(4)%array_graph(15)%array => jperp
      graph(4)%array_graph(15)%descr = 'E_3'

      graph(4)%array_graph(16:ngraph)%descr = ''

cc      sel_gr(4,:) = (/ (i,i=6,14) /)
cc      sel_gr(4,:) = (/ (i,i=1,9) /)
      sel_gr(4,:) = (/ (i,i=1,xdraw_cont_lim) /)

      prof_ivar(4,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (4,:) = 0        !No log scales
      prof_spline(4) = .false.  !No splines

c Define graphics group #5: Perturbations

      graph(5)%cartesian=.false.
cc      graph(5)%cartesian=.true.
      graph(5)%descr='Perturbations'

      do ieq = 1,neqd
        graph(5)%array_graph(ieq)%array => u_graph%array_var(ieq)%array
        graph(5)%array_graph(ieq)%descr =
     $       trim(u_graph%array_var(ieq)%descr) // ' pert'
      enddo

      graph(5)%array_graph(neqd+1:ngraph)%descr = ''

      sel_gr(5,:) = (/ (i,i=1,neqd),(0,i=neqd+1,xdraw_cont_lim) /)

      prof_ivar(5,:) = 0        !All profiles have same (default) independent variable: x
      prof_log (5,:) = 0        !No log scales
      prof_spline(5) = .false.  !No splines

c End program

      end subroutine defineGraphicsGroups

c defineGraphicsIO
c ##################################################################
      subroutine defineGraphicsIO

c--------------------------------------------------------------------
c     Set graphics files
c--------------------------------------------------------------------

      use graphics

      implicit none

c Call variables

c Local variables

      hdf5_fname   = 'pixie3d.h5'
      adios_fname  = 'pixie3d.bp'

      if (np > 1) return

      graphfile(1) = 'drawcnv.bin'
      ugraph(1)    = 20
#if defined(RFX)
      graphfile(2) = 'drawcyl.bin'
#else
      graphfile(2) = 'drawcov.bin'
#endif
      ugraph(2)    = 21
      graphfile(3) = 'drawcar.bin'
      ugraph(3)    = 22
      graphfile(4) = 'drawdiag.bin'
      ugraph(4)    = 23
      graphfile(5) = 'drawpert.bin'
      ugraph(5)    = 24

      profilefile(1)= 'drawcnv-prof.bin'
      uprofile(1)   = 30
#if defined(RFX)
      profilefile(2)= 'drawcyl-prof.bin'
#else
      profilefile(2)= 'drawcov-prof.bin'
#endif
      uprofile(2)   = 31
      profilefile(3)= 'drawcar-prof.bin'
      uprofile(3)   = 32
      profilefile(4)= 'drawdiag-prof.bin'
      uprofile(4)   = 33
      profilefile(5)= 'drawpert-prof.bin'
      uprofile(5)   = 34

      debugfile    = 'debug.bin'
      udebug       = 40

      drawgraph(1) = 'drawcnv.in'
#if defined(RFX)
      drawgraph(2) = 'drawcyl.in'
#else
      drawgraph(2) = 'drawcov.in'
#endif
      drawgraph(3) = 'drawcar.in'
      drawgraph(4) = 'drawdiag.in'
      drawgraph(5) = 'drawpert.in'

      drawprof(1)  = 'drawcnv-prof.in'
#if defined(RFX)
      drawprof(2)  = 'drawcyl-prof.in'
#else
      drawprof(2)  = 'drawcov-prof.in'
#endif  
      drawprof(3)  = 'drawcar-prof.in'
      drawprof(4)  = 'drawdiag-prof.in'
      drawprof(5)  = 'drawpert-prof.in'

      end subroutine defineGraphicsIO
