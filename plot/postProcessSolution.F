#if defined(RFX)
#include "rfx_module"
#endif

c postProcessSolution
c####################################################################
      subroutine postProcessSolution(igrid,varray)

c--------------------------------------------------------------------
c     Postprocess solution in varray to find derived quantities 
c     for plotting.
c--------------------------------------------------------------------

      use auxiliaryVariables

      use auxPlotVariables

      use ts_io

      use nlfunction_setup

      use rw_bc

      use poincare_io

      use app_iosetup, ONLY: car_diag_plots,heat_flx_diag,fourier
     $                      ,poincare,hyperResistivity

      use re
      
#if defined(RFX)
      use rfx
#endif

      implicit none

c Call variables

      integer :: igrid,nn,my_rank_rad,np_rad

      type(var_array),pointer :: varray

c Local variables

      integer :: i,j,k,ig,jg,kg,ieq,bcs(6),igx,igy,igz,nx,ny,nz,ierr
      real(8) :: x1,y1,z1,ee(3),vr,vth,norm,lnorm,cov(3),cnv(3),gradp(3)

      logical :: cartsn,covariant,to_cartsn,to_cnv,tor
      
c Begin program

      igx = igrid
      igy = igrid
      igz = igrid

      nx = gv%gparams%nxv(igx)
      ny = gv%gparams%nyv(igy)
      nz = gv%gparams%nzv(igz)

      to_cartsn = .true.
      covariant = .false.
      to_cnv    = .false.

c Impose boundary conditions
c (finds all covariant and contravariant components of interest)

cc      write (*,*) 'Dumping perturbations of all quantities'
cc      call AXPYDerivedType(1d0,varray,-1d0,u_0,u_graph)
cc      call equateDerivedType(varray,u_graph)
cc      varray%array_var(IRHO)%array = u_0%array_var(IRHO)%array
cc
cc      if (chk_src) then
cc        varray%array_var(IRHO)%array = 1d0
cc      endif
cc
cc      call applyBCplot(igx,varray,gv%aux,it=itime,bc_order=order_bc)
      call applyBC(igx,varray,gv%aux,it=itime,bc_order=order_bc)

c Find perturbed quantities (u_pert = varray - u_ic)

      call AXPYDerivedType(1d0,varray,-1d0,u_ic,u_pert)

      !Remove Jacobian factors from perturbations
      u_pert%array_var(IRHO)%array = u_pert%array_var(IRHO)%array
     $                              /gv%gparams%gmetric%grid(igrid)%jac
      u_pert%array_var(IRHO)%descr = 'Rho'

      if (solve_prs) then
        u_pert%array_var(IPRS)%array = u_pert%array_var(IPRS)%array
     $                               /gv%gparams%gmetric%grid(igrid)%jac
        u_pert%array_var(IPRS)%descr = 'Prs'
      else
        u_pert%array_var(ITMPE)%array = u_pert%array_var(ITMPE)%array
     $                               /gv%gparams%gmetric%grid(igrid)%jac
        u_pert%array_var(ITMPE)%descr = 'Te'
      endif
      
      if (alt__eom()) then
        u_pert%array_var(IVY)%array = u_pert%array_var(IVY)%array
     $                               /gv%gparams%gmetric%grid(igrid)%jac
        u_pert%array_var(IVY)%descr = 'P^2'
      endif

c Assignments
      
      vx = vcnv(:,:,:,1)
      vy = vcnv(:,:,:,2)
      vz = vcnv(:,:,:,3)

      vx_cov = vcov(:,:,:,1)
      vy_cov = vcov(:,:,:,2)
      vz_cov = vcov(:,:,:,3)

#if defined(vec_pot)
      ax_cnv = acnv(:,:,:,1)
      ay_cnv = acnv(:,:,:,2)
      az_cnv = acnv(:,:,:,3)
#endif

c Deltas
c$$$      jx = jcnv(:,:,:,1)-jcnv_0(:,:,:,1)
c$$$      jy = jcnv(:,:,:,2)-jcnv_0(:,:,:,2)
c$$$      jz = jcnv(:,:,:,3)-jcnv_0(:,:,:,3)
c$$$
c$$$      bx_cov = bcov(:,:,:,1)-bcov_0(:,:,:,1)
c$$$      by_cov = bcov(:,:,:,2)-bcov_0(:,:,:,2)
c$$$      bz_cov = bcov(:,:,:,3)-bcov_0(:,:,:,3)
c$$$
c$$$      jx_cov = jcov(:,:,:,1)-jcov_0(:,:,:,1)
c$$$      jy_cov = jcov(:,:,:,2)-jcov_0(:,:,:,2)
c$$$      jz_cov = jcov(:,:,:,3)-jcov_0(:,:,:,3)
c$$$
c$$$#if defined(vec_pot)
c$$$      bx = bcnv(:,:,:,1)-bcnv_0(:,:,:,1)
c$$$      by = bcnv(:,:,:,2)-bcnv_0(:,:,:,2)
c$$$      bz = bcnv(:,:,:,3)-bcnv_0(:,:,:,3)
c$$$#endif

      jx = jcnv(:,:,:,1)
      jy = jcnv(:,:,:,2)
      jz = jcnv(:,:,:,3)

      bx_cov = bcov(:,:,:,1)
      by_cov = bcov(:,:,:,2)
      bz_cov = bcov(:,:,:,3)

      jx_cov = jcov(:,:,:,1)
      jy_cov = jcov(:,:,:,2)
      jz_cov = jcov(:,:,:,3)

#if defined(vec_pot)
      ax_cnv = acnv(:,:,:,1)
      ay_cnv = acnv(:,:,:,2)
      az_cnv = acnv(:,:,:,3)

      bx = bcnv(:,:,:,1)
      by = bcnv(:,:,:,2)
      bz = bcnv(:,:,:,3)
#endif

c Find Cartesian components of ALL vectors

      !Cartesian velocities
      vx_car = vcnv(:,:,:,1)
      vy_car = vcnv(:,:,:,2)
      vz_car = vcnv(:,:,:,3)
      call transformVector(gv%gparams,igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,vx_car,vy_car,vz_car,'cnv','car')

#if defined(vec_pot)
      !Cartesian A components
      covariant = .true.
      ax_car = acov(:,:,:,1)
      ay_car = acov(:,:,:,2)
      az_car = acov(:,:,:,3)
      call transformVector(gv%gparams,igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,ax_car,ay_car,az_car,'cov','car')
      covariant = .false.
#endif

      !Cartesian B components
      bx_car = bcnv(:,:,:,1)!-bcnv_0(:,:,:,1)
      by_car = bcnv(:,:,:,2)!-bcnv_0(:,:,:,2)
      bz_car = bcnv(:,:,:,3)!-bcnv_0(:,:,:,3)
      call transformVector(gv%gparams,igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,bx_car,by_car,bz_car,'cnv','car')

      !Find Cartesian J components
      jx_car = jcnv(:,:,:,1)
      jy_car = jcnv(:,:,:,2)
      jz_car = jcnv(:,:,:,3)
      call transformVector(gv%gparams,igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,jx_car,jy_car,jz_car,'cnv','car')

      jmag = sqrt(jx_car**2+jy_car**2+jz_car**2)
      
c Return if checking sources

cc      if (chk_src) return

c Electron velocity

      vex = vecnv(:,:,:,1)
      vey = vecnv(:,:,:,2)
      vez = vecnv(:,:,:,3)
      call transformVector(gv%gparams,igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,vex,vey,vez,'cnv','car')

c Poloidal flux diagnostics

#if !defined(vec_pot)
      call pol_flux_diag(Pflux)
#else
      Pflux(:,:,:) = acov(:,:,:,3)
#endif

c Poloidal velocity diagnostic

      if (pinch_flow) then
        norm = maxval(sqrt(vx_car**2+vy_car**2))

#if defined(petsc)
        lnorm = norm
        call MPI_Reduce(lnorm,norm,1,MPI_DOUBLE_PRECISION
     .                 ,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#endif 
        if (my_rank == 0) write (*,*) 'Maximum poloidal velocity=',norm
      endif

c Resistive wall diagnostics

      if (res_wall_bc) then
        if (my_rank == 0) write (*,*) '||Bn||',time,
     .   sqrt(integral(gv%gparams,igx,bcov(:,:,:,1)**2,average=.false.))
        if (my_rank == 0) write (*,*) '||Bt||',time,
     .   sqrt(integral(gv%gparams,igx,bcov(:,:,:,2)**2,average=.false.))
      endif

c Runaway electron diagnostics

      if (re_pkg) then
        if (.not.associated(nre)) then
          allocate(nre(0:nx+1,0:ny+1,0:nz+1)
     .            ,jre(0:nx+1,0:ny+1,0:nz+1,3))
        endif
        call re_qtys(eeta,bcnv,bcov,E_ni,nre,jre)

        jrex_car = jre(:,:,:,1)
        jrey_car = jre(:,:,:,2)
        jrez_car = jre(:,:,:,3)
        call transformVector(gv%gparams,igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                      ,jrex_car,jrey_car,jrez_car,'cov','car')
      endif
      
c Velocity and current intrinsic components (wrt B)

      do k = kming,kmaxg
        do j = jming,jmaxg
          do i = iming,imaxg
            call getCurvilinearCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .                                    ,ig,jg,kg,x1,y1,z1)

            !V_parallel
            norm = sqrt(vectorNorm(gv%gparams,i,j,k,igx,bcnv(i,j,k,:)
     .                            ,.false.))

            if (norm == 0d0) then
              vpar(i,j,k) = 0d0
            else
              vpar(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx
     .                                   ,bcov(i,j,k,:)
     .                                   ,vcnv(i,j,k,:))/norm
            endif

            !V.grad(psi)/|grad(psi)|
            cov = grad(gv%gparams,i,j,k,nx,ny,nz,igx,igy,igz,Pflux)
            norm = sqrt(vectorNorm(gv%gparams,i,j,k,igx,cov,.true.))

            if (norm == 0d0) then
              vpsi(i,j,k) = 0d0
            else
              vpsi(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx,cov
     .                                   ,vcnv(i,j,k,:))/norm
            endif

            !V.(b x grad(psi))/|(b x grad(psi))|
            cnv = crossProduct(gv%gparams,i,j,k,igx,bcov(i,j,k,:),cov
     .                        ,.false.)

            norm = sqrt(vectorNorm(gv%gparams,i,j,k,igx,cnv,.false.))

            if (norm == 0d0) then
              vperp(i,j,k) = 0d0
            else
              vperp(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx
     .                                    ,vcov(i,j,k,:),cnv)/norm
            endif
            
cc            !J_parallel
cc            norm = sqrt(vectorNorm(gv%gparams,i,j,k,igx,bcnv(i,j,k,:),.false.))
cc
cc            if (norm == 0d0) then
cc              jpar(i,j,k) = 0d0
cc            else
cc              jpar(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx
cc     .                                   ,bcov(i,j,k,:),jcnv(i,j,k,:))
cc     .                     /norm
cc            endif
cc
cc            !J.grad(psi)/|grad(psi)|
cc            cov =  grad(i,j,k,nx,ny,nz,igx,igy,igz,Pflux)
cc            norm = sqrt(vectorNorm(gv%gparams,i,j,k,igx,cov,.true.))
cc
cc            if (norm == 0d0) then
cc              jpsi(i,j,k) = 0d0
cc            else
cc              jpsi(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx,cov,jcnv(i,j,k,:))
cc     .                     /norm
cc            endif
cccc            vpsi(i,j,k) = cov(1)
cc
cc            !J.(b x grad(psi))/|(b x grad(psi))|
cc            cnv(1) = bcov(i,j,k,2)*cov(3)
cc     .              -bcov(i,j,k,3)*cov(2)
cc            cnv(2) = bcov(i,j,k,3)*cov(1)
cc     .              -bcov(i,j,k,1)*cov(3)
cc            cnv(3) = bcov(i,j,k,1)*cov(2)
cc     .              -bcov(i,j,k,2)*cov(1)
cc
cc            norm = sqrt(vectorNorm(gv%gparams,i,j,k,igx,cnv,.false.))
cc
cc            if (norm == 0) then
cc              jperp(i,j,k) = 0d0
cc            else
cc              jperp(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx,jcov(i,j,k,:),cnv)
cc     .                       /norm
cc            endif
          enddo
        enddo
      enddo

c Hall electric field components (covariant)

      if (di > 0d0) then
        E1  = di*E_h(:,:,:,1)
        E2  = di*E_h(:,:,:,2)
        E3  = di*E_h(:,:,:,3)
      else
        E1  = E_ni(:,:,:,1)
        E2  = E_ni(:,:,:,2)
        E3  = E_ni(:,:,:,3)
      endif

      if (car_diag_plots)
     .  call transformVector(gv%gparams,igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                      ,E1,E2,E3,'cov','car')
      
c Pressure diagnostics

      p_tot = rho*tmp

c Hyper-resistivity diagnostic (solve for hyperresistivity lambda)

      if (hyperResistivity.and.TS_chk_plot_dump(.false.))
     $     call hyperres_diag(hlambda)

c grad(P) diagnostic (relative to the magnetic field)

      do k = kming,kmaxg
        do j = jming,jmaxg
          do i = iming,imaxg
            gradp = grad(gv%gparams,i,j,k,nx,ny,nz,igx,igy,igz,p_tot)

            call getCurvilinearCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .                                    ,ig,jg,kg,x1,y1,z1)

            !grad(P)_parallel
            norm = sqrt(vectorNorm(gv%gparams,i,j,k,igx,bcnv(i,j,k,:)
     .                            ,.false.))

            if (norm == 0d0) then
              gppar(i,j,k) = 0d0
            else
              gppar(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx
     .                                   ,bcnv(i,j,k,:)
     .                                   ,gradp(:))/norm
            endif

            !grad(P).grad(psi)/|grad(psi)|
            cov = grad(gv%gparams,i,j,k,nx,ny,nz,igx,igy,igz,Pflux)
            norm = sqrt(vectorNorm(gv%gparams,i,j,k,igx,cov,.true.))
            cnv = XformToCnv_ijk(gv%gparams,i,j,k,igx,cov)

            if (norm == 0d0) then
              gppsi(i,j,k) = 0d0
            else
              gppsi(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx,cnv
     .                                   ,gradp(:))/norm
            endif

            !grad(P).(b x grad(psi))/|(b x grad(psi))|  
            cnv = crossProduct(gv%gparams,i,j,k,igx,bcov(i,j,k,:),cov
     .                        ,.false.)
            ! We now need grad(P) in cov:
            norm = sqrt(vectorNorm(gv%gparams,i,j,k,igx,cnv,.false.))

            if (norm == 0d0) then
              gpperp(i,j,k) = 0d0 
            else
              gpperp(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx
     .                                    ,gradp(:),cnv)/norm
            endif

          enddo
        enddo
      enddo

c JxB-g(P) diagnostic (cnv)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            if (nc_eom_jxb) then
              cnv = crossProduct(gv%gparams,i,j,k,igx
     $                          ,jcov(i,j,k,:),bcov(i,j,k,:),.false.)
            else
              cnv = EOM_jxb(gv%gparams,i,j,k,igx,bcnv,bcnv)
            endif

            if (use_p_eom) then
              cnv = cnv
     $             -EOM_gp (gv%gparams,i,j,k,igx,p_tot,nc_eom_gp)
            else
              cnv = cnv
     $             -EOM_gnT(gv%gparams,i,j,k,igx,rho,tmp,nc_eom_gp)
            endif

            F1(i,j,k) = cnv(1)
            F2(i,j,k) = cnv(2)
            F3(i,j,k) = cnv(3)
            if (alt__eom().and.(.not.car_diag_plots))
     $     F2(i,j,k) = F2(i,j,k)*gv%gparams%gmetric%grid(igx)%jac(i,j,k)
          enddo
        enddo
      enddo

      if (car_diag_plots)
     .  call transformVector(gv%gparams,igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,F1,F2,F3,'cov','car')

c Ohm's law diagnostic

cc      write (*,*) 'Current',J0
cc      i = nxd
cc      k = nzd
cc      do j = 1,nyd
cc        ee(1) = eta*(-J0(1)+0.5*(jx_cov(i,j,k)+jx_cov(i+1,j,k)))
cc     .             - (0.5*(vy(i,j,k)+vy(i+1,j,k))
cc     .               *0.5*(bz(i,j,k)+bz(i+1,j,k))
cc     .               -0.5*(vz(i,j,k)+vz(i+1,j,k))
cc     .               *0.5*(by(i,j,k)+by(i+1,j,k)))
cc     .              /0.5/(gv%gparams%gmetric%grid(1)%jac(i  ,j,k)
cc     .                   +gv%gparams%gmetric%grid(1)%jac(i+1,j,k))
cc        ee(2) = eta*(-J0(2)+0.5*(jy_cov(i,j,k)+jy_cov(i+1,j,k)))
cc     .             - (0.5*(vz(i,j,k)+vz(i+1,j,k))
cc     .               *0.5*(bx(i,j,k)+bx(i+1,j,k))
cc     .               -0.5*(vx(i,j,k)+vx(i+1,j,k))
cc     .               *0.5*(bz(i,j,k)+bz(i+1,j,k)))
cc     .              /0.5/(gv%gparams%gmetric%grid(1)%jac(i  ,j,k)
cc     .                   +gv%gparams%gmetric%grid(1)%jac(i+1,j,k))
cc        ee(3) = eta*(-J0(3)+0.5*(jz_cov(i,j,k)+jz_cov(i+1,j,k)))
cc     .             - (0.5*(vx(i,j,k)+vx(i+1,j,k))
cc     .               *0.5*(by(i,j,k)+by(i+1,j,k))
cc     .               -0.5*(vy(i,j,k)+vy(i+1,j,k))
cc     .               *0.5*(bx(i,j,k)+bx(i+1,j,k)))
cc     .              /0.5/(gv%gparams%gmetric%grid(1)%jac(i  ,j,k)
cc     .                   +gv%gparams%gmetric%grid(1)%jac(i+1,j,k))
cc        write (*,*) 'Ohms law at',i,j,k,' = ',ee
cc      enddo

c Toroidal diagnostics (quasi-linear q-profile, lambda profile)

      if (     coords == 'hel'
     .    .or. coords == 'cyl'
     .    .or. coords == 'tor'
     .    .or. trim(equil) == 'vmec'
     $    .or. trim(equil) == 'efit') then
        call tor_diag
      else
        qfactor = 0d0
        lambda  = 0d0
      endif

c Fourier output

      if (np == 1 .and. fourier.and.TS_chk_plot_dump(.false.)) then
        call fourier_dump
      elseif (np > 1 .and. fourier .and. my_rank == 0) then
        write (*,*) 'Fourier dump not available in parallel'
      endif

c Poincare output (field tracer)

      tor = (equil  == 'vmec'
     $  .or. coords == 'tor'
     $  .or. coords == 'sha'
     $  .or. coords == 'sh2'
     $  .or. coords == 'sh3'
     $  .or. coords == 'stl')

      if (poincare.and.TS_chk_plot_dump(.false.))
     . call poincare_plot_dump(gv%gparams,time,bcnv,bx_car,by_car,bz_car
     .                        ,tor)

c Divergence diagnostics

      divrgJ = div(gv%gparams,igx,jcnv)
      divrgB = div(gv%gparams,igx,bcnv)
      divrgV = div(gv%gparams,igx,vcnv)

c Heat source diagnostics

      !Internal heat source
      if ((.not.adiabatic).and.(gamma > 1d0)) then
        do k=1,nz
          do j=1,ny
            do i=1,nx
              !Internal heat source
              heat_src(i,j,k) =
     .             int_heat_src(gv%gparams,i,j,k,igx,eeta,nuu,h_eta
     .                         ,jcov,vcnv,vecnv)

              !External heat source
              heat_src(i,j,k) = heat_src(i,j,k)
     .        +ext_heat_src(gv%gparams,i,j,k,igx,rho(i,j,k),tmpe(i,j,k))
            enddo
          enddo
        enddo
      else
        heat_src = 0d0
      endif

#if defined(RFX)
c RFX dump
      if (np == 1.and.TS_chk_plot_dump(.false.)) then
        call rfx_dump_fields
      elseif (np > 1 .and. my_rank == 0) then
        write (*,*) 'RFX fields dump not available in parallel'
      endif
#endif

c Fields dump

cc      if (itime == 0) then
cc        write (*,*) 'Removing dump files...'
cc        ierr=rm_files('*.txt')
cc      else
cc        call vector_dump('b-field',bx_car,by_car,bz_car)
cc        call vector_dump('v-field',vx_car,vy_car,vz_car)
cc        call vector_dump('j-field',jx_car,jy_car,jz_car)
cc        call scalar_dump('rho-field',rho)
cc      endif

c End program

      contains

c     toroidal_average
c     ################################################################
      subroutine toroidal_average(TargetVarIn,AvgVarOut)

c     ----------------------------------------------------------------
c     Averages a provided scalar or vector component along the z-axis.
c     ----------------------------------------------------------------

      implicit none

c     Call variables

      real(8),dimension(0:nx+1,0:ny+1,0:nz+1),intent(in) :: TargetVarIn
      real(8),dimension(0:nx+1,0:ny+1),intent(out) :: AvgVarOut
      real(8),dimension(0:nx+1,0:ny+1) :: dummy
      
c     Local variables

c     Begin program

      AvgVarOut = sum(TargetVarIn(:,:,1:nz),3)
#if defined(petsc)
      dummy = AvgVarOut
      call MPI_Allreduce(dummy,AvgVarOut,size(AvgVarOut)
     .                  ,MPI_DOUBLE_PRECISION,MPI_SUM
     .                  ,gv%gparams%MPI_COMM_Z,mpierr)
#endif
      AvgVarOut = AvgVarOut/nzd
      
      end subroutine toroidal_average

c     pol_flux_diag
c     ################################################################
      subroutine pol_flux_diag(pflux)

c     ---------------------------------------------------------------
c     Calculates poloidal flux from B-field in parallel.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: pflux(0:,0:,0:)

c     Local variables

      real(8),allocatable,dimension(:,:,:) :: bcnvavg
#if defined(petsc)
      real(8),allocatable,dimension(:,:,:) ::dummy
      real(8),allocatable,dimension(:) :: mpivec
#endif

      integer :: np_mpi,my_rank_mpi
      
c     Begin program
      
      !Average fields toroidally
      allocate(bcnvavg(0:nx+1,0:ny+1,3))
      
      bcnvavg = sum(bcnv(:,:,1:nz,:),3)
#if defined(petsc)
      allocate(dummy(0:nx+1,0:ny+1,3))
      dummy = bcnvavg
      call MPI_Allreduce(dummy,bcnvavg,size(bcnvavg)
     .                  ,MPI_DOUBLE_PRECISION,MPI_SUM
     .                  ,gv%gparams%MPI_COMM_Z,mpierr)
      deallocate(dummy)
#endif
      bcnvavg = bcnvavg/nzd
      
      !Integrate poloidal flux
      Pflux(0,0,0) = 0d0   !Set value at origin

      do i = 1,nx+1
        call getMGmap(gv%gparams,i,1,1,igx,igy,igz,ig,jg,kg)
        Pflux(i,0,0) = Pflux(i-1,0,0)
     .                 - 0.5*(bcnvavg(i-1,0,2)
     .                       +bcnvavg(i  ,0,2))*gv%gparams%dx(ig-1)
      enddo

#if defined(petsc)
      !Communicate radial integration constants
      if (npx > 1) then
        nn = 1
        allocate(mpivec(npx*nn))
        call MPI_Allgather(Pflux(nx,0,0)
     .                          ,nn,MPI_DOUBLE_PRECISION
     .                   ,mpivec,nn,MPI_DOUBLE_PRECISION
     .                   ,gv%gparams%MPI_COMM_X,mpierr)

        !Sum integration constants
        do j=2,npx
          mpivec(1+nn*(j-1):nn*j) = mpivec(1+nn*(j-1):nn* j   )
     .                             +mpivec(1+nn*(j-2):nn*(j-1))
        enddo

        !Add integration constants
        call initMPI(gv%gparams%MPI_COMM_X,np_mpi,my_rank_mpi)

        j = my_rank_mpi
        if (j > 0) then
        do i=0,nx+1
            Pflux(i,0:0,0) = Pflux(i,0:0,0) + mpivec(1+nn*(j-1):nn*j)
          enddo
        endif

        deallocate(mpivec)
      endif
#endif

      do j = 1,ny+1
        call getMGmap(gv%gparams,1,j,1,igx,igy,igz,ig,jg,kg)
        Pflux(:,j,0) = Pflux(:,j-1,0)
     .               + 0.5*(bcnvavg(:,j-1,1)
     .                     +bcnvavg(:,j  ,1))*gv%gparams%dy(jg-1)
      enddo
 
c$$$#if defined(petsc)
c$$$      !Communicate poloidal integration constants (not needed because all start at SP)
c$$$      if (npy > 1.and.(.not.bcSP())) then
c$$$        nn = nx+2
c$$$        allocate(mpivec(npy*nn))
c$$$        call MPI_Allgather(Pflux(0:nx+1,ny,0)
c$$$     .                           ,nn,MPI_DOUBLE_PRECISION
c$$$     .                    ,mpivec,nn,MPI_DOUBLE_PRECISION
c$$$     .                    ,gv%gparams%MPI_COMM_Y,mpierr)
c$$$
c$$$        !Add integration constants
c$$$        call initMPI(gv%gparams%MPI_COMM_Y,np_mpi,my_rank_mpi)
c$$$
c$$$        !Sum integration constants
c$$$        do j=2,npy
c$$$          mpivec(1+nn*(j-1):nn*j) = mpivec(1+nn*(j-1):nn* j   )
c$$$     .                             +mpivec(1+nn*(j-2):nn*(j-1))
c$$$        enddo
c$$$
c$$$        i = my_rank_mpi 
c$$$        if (i > 0) then
c$$$          do j=0,ny+1
c$$$            Pflux(:,j,0) = Pflux(:,j,0) + mpivec(1+nn*(i-1):nn*i)
c$$$          enddo
c$$$        endif
c$$$
c$$$        deallocate(mpivec)
c$$$      endif
c$$$#endif

      do k=1,nz+1
        Pflux(:,:,k) = Pflux(:,:,0)
      enddo
        
      deallocate(bcnvavg)
      
      end subroutine pol_flux_diag

c     tor_diag
c     ################################################################
      subroutine tor_diag

c     ---------------------------------------------------------------
c     Calculates toroidal diagnostics q-profile and lambda.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

      if (     coords == 'hel'
     .    .or. coords == 'cyl'
     .    .or. coords == 'tor') call q_profile_simple

      if (heat_flx_diag) call heat_flx_tor_diag

      end subroutine tor_diag

c     q_profile_simple
c     ################################################################
      subroutine q_profile_simple

c     ---------------------------------------------------------------
c     Calculates toroidal diagnostics q-profile and lambda.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      real(8) ::  mm,kk,RR,ll,dv
      real(8),allocatable,dimension(:,:) :: b00,j00_cov,v00,q00,l00,ldum

c     Begin program

      allocate(b00(iming:imaxg,3),j00_cov(iming:imaxg,3)
     .        ,v00(iming:imaxg,3),q00    (iming:imaxg,1)
     .        ,l00(iming:imaxg,1))

      if (    coords == 'hel') then
        mm = gv%gparams%params(1)
        kk = gv%gparams%params(2)
        RR = gv%gparams%params(3)
      elseif (coords == 'cyl') then
        mm = 1d0
        kk = 0d0
        RR = gv%gparams%params(3)
      elseif (coords == 'tor') then
        mm = 1d0
        kk = 0d0
        RR = 1d0
      endif

      !Find mean fields
      do i = iming,imaxg
        b00    (i,:) = 0d0
        j00_cov(i,:) = 0d0
        v00    (i,:) = 0d0
        q00    (i,:) = 0d0
        l00    (i,:) = 0d0
        ll = 0d0
        do k = kming,kmaxg-1   !Careful with the periodic limits
          do j = jming,jmaxg-1 !Careful with the periodic limits
            call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)

            dv = gv%gparams%dyh(jg)*gv%gparams%dzh(kg)

            b00(i,1) = b00(i,1) + bcnv(i,j,k,1)*dv
            b00(i,2) = b00(i,2) + bcnv(i,j,k,2)*dv
            b00(i,3) = b00(i,3) + bcnv(i,j,k,3)*dv

            j00_cov(i,1) = j00_cov(i,1) + jcov(i,j,k,1)*dv
            j00_cov(i,2) = j00_cov(i,2) + jcov(i,j,k,2)*dv
            j00_cov(i,3) = j00_cov(i,3) + jcov(i,j,k,3)*dv

            v00(i,1) = v00(i,1) + vcov(i,j,k,1)*dv
            v00(i,2) = v00(i,2) + vcov(i,j,k,2)*dv
            v00(i,3) = v00(i,3) + vcov(i,j,k,3)*dv

cc               q00(i,1) = q00(i,1)
cc    .                + bcnv(i,j,k,3)/RR
cc    .                   /(bcnv(i,j,k,2)-kk/mm*bcnv(i,j,k,3))*dv
            l00(i,1) = l00(i,1)
     .               + scalarProduct(gv%gparams,i,j,k,igx,jcov(i,j,k,:)
     .                              ,bcnv(i,j,k,:))
     .                 /vectorNorm(gv%gparams,i,j,k,igx,bcnv(i,j,k,:)
     .                            ,.false.)
     .                *dv

            ll = ll + dv
          enddo
        enddo
      enddo

      !Communicate sums for each radial slice
#if defined(petsc)
      allocate(ldum(iming:imaxg,3))

      ldum = b00
      call MPI_Allreduce(ldum,b00,3*(imaxg-iming+1)
     .                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,gv%gparams%MPI_COMM_YZ,mpierr)

      ldum = j00_cov
      call MPI_Allreduce(ldum,j00_cov,3*(imaxg-iming+1)
     .                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,gv%gparams%MPI_COMM_YZ,mpierr)

      ldum = v00
      call MPI_Allreduce(ldum,v00,3*(imaxg-iming+1)
     .                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,gv%gparams%MPI_COMM_YZ,mpierr)

      ldum(:,1) = l00(:,1)
      call MPI_Allreduce(ldum(:,1),l00(:,1),(imaxg-iming+1)
     .                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,gv%gparams%MPI_COMM_YZ,mpierr)
      deallocate(ldum)

      ll = psum((/ll/),mpi_comm=gv%gparams%MPI_COMM_YZ)
c$$$      dv = ll
c$$$      call MPI_Allreduce(dv,ll,1
c$$$     .                  ,MPI_DOUBLE_PRECISION
c$$$     .                  ,MPI_SUM,gv%gparams%MPI_COMM_YZ,mpierr)
#endif

      !Average radial quantities
      b00     = b00    /ll
      j00_cov = j00_cov/ll
      v00     = v00    /ll
      q00     = q00    /ll
      l00     = l00    /ll

      !Find q-factor and lambda profile
      do i = iming,imaxg
        qfactor(i,:,:) = mm*b00(i,3)/RR/(b00(i,2)-kk*b00(i,3))
        lambda (i,:,:) = l00(i,1)
      enddo
      qfactor(iming,:,:)=qfactor(iming+1,:,:)
      lambda (iming,:,:)=lambda (iming+1,:,:)

c diag  : dump text data
cc          if (time == 0d0) then
cc            open(unit=1000,file='q-lambda.txt',status='unknown')
cc          else
cc            open(unit=1000,file='q-lambda.txt',status='unknown'
cc     .        ,position='append')
cc          endif
cc          write(1000,*) 'Time = ',time
cc          write(1000,*)'      Qlnr Q-prof    Avgd Q-prof    Qlnr Lambda'
cc     .                 ,'      Avgd Lambda         Vr'
cc          do i=iming,imaxg
cc            write (1000,*) qfactor(i,1,1),lambda (i,1,1),v00(i,1)
cc          enddo
cc          write (1000,*)
cc          close(1000)

#if defined(RFX)
      !RFX dump
      if (np == 1) then
        call rfx_dump_mean_profiles(b00,j00_cov,v00)
      elseif (np > 1 .and. my_rank == 0) then
        write (*,*) 'RFX mean-profile dump not available in parallel'
      endif
#endif

      deallocate(b00,j00_cov,v00,q00,l00)

      end subroutine q_profile_simple

c     heat_flx_tor_diag
c     ################################################################
      subroutine heat_flx_tor_diag

c     ---------------------------------------------------------------
c     Calculates heat flux diagnostic at radial boundary
c     ---------------------------------------------------------------

      use pixie3d_io
      
      implicit none

c     Call variables

c     Local variables

      integer :: nx,ny,nz,nxg,nyg,nzg,igr,funit,j,k,my_rank_l
      real(8),allocatable,dimension(:,:,:,:) :: heat_flx,heat_flx_g
     $                                         ,dum1,dum2,chi_tnsr

      character(100) :: fileout
      
c     Begin program

      igr = 1
      
      nx = gv%gparams%nxv(igr)
      ny = gv%gparams%nyv(igr)
      nz = gv%gparams%nzv(igr)

      nxg = gv%gparams%nxgl(igr)
      nyg = gv%gparams%nygl(igr)
      nzg = gv%gparams%nzgl(igr)

c     Allocate variables
      
      allocate(heat_flx  (0:nx+1,0:ny+1,0:nz+1,3)
     .        ,dum1      (0:nx+1,0:ny+1,0:nz+1,1)
     .        ,dum2      (0:nx+1,0:ny+1,0:nz+1,1)
     .        ,heat_flx_g(0:nx+1,0:nyg+1,0:nzg+1,3)
     .        ,chi_tnsr  (0:nx+1,0:ny+1,0:nz+1,6))

      heat_flx = 0d0 ; heat_flx_g = 0d0
      dum1 = 0d0 ; dum2 = 0d0; chi_tnsr = 0d0

c     Compute heat flux
      
      dum1(:,:,:,1) = cchi_par
      dum2(:,:,:,1) = cchi
      chi_tnsr = find_diff_tnsr(gv%gparams,igr,bcnv,dum1,dum2,6)
      heat_flx = tnsr_flx_mesh(gv%gparams,igr,tmpe
     $                        ,varray%array_var(IPRS)%bconds,chi_tnsr
     $                        ,xport_order,.true.,limit_mthd=SMART)

c     Find global array
      
      call find_global_nobc(heat_flx  (1:nx,1:ny ,1:nz ,:)
     $                     ,heat_flx_g(1:nx,1:nyg,1:nzg,:)
#if defined(petsc)
     $                     ,mpi_comm=gv%gparams%MPI_COMM_YZ)
#else
     $     )
#endif

c     Dump array (gnuplot format)
      
#if defined(petsc)
      call MPI_Comm_rank(gv%gparams%MPI_COMM_YZ,my_rank_l,mpierr)
#else
      my_rank_l = 0
#endif
      if (my_rank_l == 0.and.isBdry(gv%gparams,nx,igr,2)) then
        write (fileout,'(a,1pe9.3,a)') "heat_flx_time=",time,".dat"
        call info("Dumping "//trim(fileout))
        funit = find_unit(1234)
        open(unit=funit,file=trim(fileout),status="unknown")
        write (funit,*) "#Plot with gnuplot with:"
        write (funit,*) "#  set autoscale xfix"                                           
        write (funit,*) "#  set autoscale yfix"                                          
        write (funit,*) "#  set autoscale cbfix"                                          
        write (funit,*) "#  set pm3d map"                                         
        write (funit,*) '#  set xlabel "{/Symbol f}"'                                           
        write (funit,*) '#  set ylabel "{/Symbol q}"'                                             
        write (funit,*) '#  set title "'//trim(fileout)//'"'                                          
        write (funit,*) '#  splot "'//trim(fileout)
     $                   //'" matrix nonuniform notitle'
        write (funit,*) nzg,gv%gparams%zg(1:nzg)
        do j=1,nyg
          write(funit,*) gv%gparams%yg(j),heat_flx_g(nx,j,1:nzg,1)
        enddo
        close(funit)
      endif
        
      deallocate(heat_flx,heat_flx_g,dum1,dum2,chi_tnsr)

      end subroutine heat_flx_tor_diag
      
c     fourier_dump
c     ################################################################
      subroutine fourier_dump

c     ---------------------------------------------------------------
c     Integrates array(i,j,k) on periodic dimensions of domain
c     (nx)x(ny)x(nz) to find quasi-linear profiles.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      character(30) :: ft_file
      integer :: uft=2345

c     Begin program

c     Open file

      write (ft_file,'(1pe10.3)') 0.01*int(time*100)
      ft_file = 'fourier_t='//trim(ft_file)//'.txt'

      open(unit=uft,file=trim(ft_file),status='unknown')

c     Dump arrays

      write (*,*) 'Dumping Fourier file ',ft_file

      do i=1,nxd
        do j=1,nyd
          do k=1,nzd
            write (uft,*) bx_car(i,j,k),by_car(i,j,k),bz_car(i,j,k)
          enddo
        enddo
      enddo

      write (uft,*)

c     Close file

      close(uft)

c     End program

      end subroutine fourier_dump

c     vector_dump
c     ################################################################
      subroutine vector_dump(qty,vx,vy,vz)

c     ---------------------------------------------------------------
c     Dumps vector components to an ascii file
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      character(*) :: qty
      real(8) :: vx(0:nx+1,0:ny+1,0:nz+1)
     .          ,vy(0:nx+1,0:ny+1,0:nz+1)
     .          ,vz(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      character(30) :: file
      integer :: uft=2345

c     Begin program

c     Open file

cc      write (file,'(1pe9.3)') 0.01*int(time*100)
cc      file = trim(qty)//'_t='//trim(file)//'.txt'
      file = trim(qty)//'.txt'

      if (.not.is_file(file)) then
        open(unit=uft,file=trim(file),status='new')

        write (*,*) 'Creating vector file ',file

        write (uft,'(3i5)') nx,ny,nz
      else
        open(unit=uft,file=trim(file),status='old',position='append')
      endif

c     Dump arrays

      write (*,*) 'Dumping vector data on ',file

      write (uft,'(e12.3)') time
      do k=1,nz
        do j=1,ny
          do i=1,nx
            write(uft,'(3i5,3e12.3)')i,j,k,vx(i,j,k),vy(i,j,k),vz(i,j,k)
          enddo
        enddo
      enddo

      write (uft,*)

c     Close file

      close(uft)

c     End program

      end subroutine vector_dump

c     scalar_dump
c     ################################################################
      subroutine scalar_dump(qty,scl)

c     ---------------------------------------------------------------
c     Dumps scalar to an ascii file
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      character(*) :: qty
      real(8) :: scl(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      character(30) :: file
      integer :: uft=2345

c     Begin program

c     Open file

cc      write (file,'(1pe9.3)') 0.01*int(time*100)
cc      file = trim(qty)//'_t='//trim(file)//'.txt'
      file = trim(qty)//'.txt'

      if (.not.is_file(file)) then
        open(unit=uft,file=trim(file),status='new')

        write (*,*) 'Creating scalar file ',file

        write (uft,'(3i5)') nx,ny,nz
      else
        open(unit=uft,file=trim(file),status='old',position='append')
      endif

c     Dump arrays

      write (*,*) 'Dumping scalar data on ',file

      write (uft,'(e12.3)') time
      do k=1,nz
        do j=1,ny
          do i=1,nx
            write (uft,'(3i5,e12.3)') i,j,k,scl(i,j,k)
          enddo
        enddo
      enddo

      write (uft,*)

c     Close file

      close(uft)

c     End program

      end subroutine scalar_dump

c     hyperres_diag
c     ##################################################################
      subroutine hyperres_diag(lambda)

      use fixed_point_accelerator

      use B_tools, xi2 => xi
      
      integer :: igrid,ii,i,j,k,guess,nn,bc_lbda(6),bc_xi(6)
     $          ,order,it,bcsb(6,3),bcsa(6,3),bcs(6,3)
      
      real(8) :: vxB(3),BMagAvg(3),ResOhmsLaw(3),dvol,vol,lerror,l0,l1
     $          ,lambda(0:nx+1,0:ny+1,0:nz+1)
            
      real(8),allocatable,dimension(:,:,:)   :: JparAvg,RHSfluc,BAvg
     $                                         ,RHSsrc,JAvg,diffcf,xxi
     $                                         ,lambda_o,EAvg,dummy
      real(8),allocatable,dimension(:,:,:,:) :: gJpar,gJpar_cnv
     $                                         ,ecnv,ajcov,ajcov_o,ecov
      real(8),allocatable,dimension(:,:)     :: etaAvg,x,y
      real(8),allocatable,dimension(:)       :: lbda_err,xi,xa
      
      type (fpa_state) :: fpav
      
      igrid = igx

      allocate(BAvg     (0:nx+1,0:ny+1,3)
     .        ,JAvg     (0:nx+1,0:ny+1,3)
     .        ,EAvg     (0:nx+1,0:ny+1,3)
     .        ,etaAvg   (0:nx+1,0:ny+1)
     .        ,JparAvg  (0:nx+1,0:ny+1,0:nz+1)
     .        ,RHSfluc  (0:nx+1,0:ny+1,0:nz+1)
     .        ,RHSsrc   (0:nx+1,0:ny+1,0:nz+1)
     .        ,xxi      (0:nx+1,0:ny+1,0:nz+1)
     .        ,lambda_o (0:nx+1,0:ny+1,0:nz+1)
     .        ,diffcf   (0:nx+1,0:ny+1,0:nz+1)
     .        ,gJpar    (0:nx+1,0:ny+1,0:nz+1,3)
     .        ,gJpar_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,ecnv     (0:nx+1,0:ny+1,0:nz+1,3)
     .        ,ecov     (0:nx+1,0:ny+1,0:nz+1,3)
     .        ,ajcov    (0:nx+1,0:ny+1,0:nz+1,3)
     .        ,ajcov_o  (0:nx+1,0:ny+1,0:nz+1,3))

      nn = nx*ny*nz
      allocate(x(nn,1),y(nn,1),lbda_err(nn),xi(nn),xa(3*nn))

      ecnv = XformToCnv(gv%gparams,igx,E_ni)

      !Average each vector component toroidally
      do i = 1,3
        call toroidal_average(bcov(:,:,:,i),BAvg(:,:,i))
        call toroidal_average(jcnv(:,:,:,i),JAvg(:,:,i))
        call toroidal_average(ecnv(:,:,:,i),EAvg(:,:,i))
      end do
      call toroidal_average(eeta,etaAvg) ! Average scalar toroidally

      !Compute Ohm's law sources
      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            vxB = crossProduct(gv%gparams,i,j,k,igx
     .                        ,vcov(i,j,k,:),bcov(i,j,k,:),.false.)

            BMagAvg = BAvg(i,j,:)
     $              /vectorNorm(gv%gparams,i,j,k,igx,BAvg(i,j,:),.true.)

            JparAvg(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx
     .                                    ,BMagAvg(:)
     .                                    ,JAvg(i,j,:))

            ResOhmsLaw = ecnv(i,j,k,:) - EAvg(i,j,:) + vxB
     .                   - eeta(i,j,k)*(jcnv(i,j,k,:)-JAvg(i,j,:))
     .                   !+ etaAvg(i,j)*JAvg(i,j,:)

            RHSfluc(i,j,k) = scalarProduct(gv%gparams,i,j,k,igx
     .                                    ,BAvg(i,j,:)
     .                                    ,ResOhmsLaw(:))
          enddo
        enddo
      enddo

      lbda_err = reshape(RHSfluc(1:nx,1:ny,1:nz),shape(lbda_err))

      lerror = rms(lbda_err)

      if (lerror < 1d-12) return

      !Compute gradient of j_|| (in cnv representation)
      bcs(:,1) = bcond
      bcs(:,2) = bcond
      bcs(:,3) = bcond
      where(bcs == DEF) bcs =-EXT
      
      gJpar = grad(gv%gparams,igx,JparAvg)
      call setBC_ez(gv%gparams,igx,gJpar,.false.,bcs,iorder=order_bc)
      gJpar_cnv = XformToCnv(gv%gparams,igx,gJpar)

c     Find scalar potential xi

      bc_xi = bcond
      where(bc_xi == DEF) bc_xi = NEU !Homogeneous Neumann

      !Compute residual
      vol = 0d0
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            dvol = gv%gparams%gmetric%grid(igx)%dvol(i,j,k)
            y(ii,1) = RHSfluc(i,j,k)*dvol
            vol = vol + dvol
          enddo
        enddo
      enddo
      vol = psum((/vol/))
      vol = psum(y(:,1))/vol

      !Regularize residual by subtracting mean (needed by Neumann BCs)
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)
            dvol = gv%gparams%gmetric%grid(igx)%dvol(i,j,k)
            y(ii,1) = y(ii,1) - dvol*vol
          enddo
        enddo
      enddo
        
      !Solve for scalar potential
      if (my_rank == 0) write (*,*) " Hyper-res grad(xi) solve..."
      guess = 0
      call cSolver(gv%gparams,1,nn,y,x,bc_xi
     .            ,igrid,ilevel,guess,lapneu_mtvc,.true.
     .            ,tol           = rtol
     .            ,gm_driver     = .true. !pc_gm_driver
     .            ,mg_debug      = .false. !mg_debug
     .            ,mg_galerkin   = .false. !mg_glrkin
     .            ,mg_order_res  = 0 !mg_ores
     .            ,mg_order_prol = 1 !mg_oprol
     .            ,mg_gm_coarse_solve  = .false. !mg_gm_coarse
     .            ,mg_coarse_grid_level = 2
     .            ,mg_vcyc       = 100 !mg_vcyc
     .            ,mg_smooth     = "jb" !sm_type
     .            ,sm_it         = 5 !sm_iter
     .            ,sm_omega      = 0.7d0 !sm_omega
     .            ,sm_ncolors    = 4   !sm_ncolors
     .            ,sm_line_relax = .false. !sm_line_relax
     .            ,sm_zebra_relax= .false. !sm_zebra_relax
     .            ,sm_zebra_it   = sm_zebra_it
     .            ,sing_matrix   = .true.)

      !Find xi (with BCs)
      xxi(1:nx,1:ny,1:nz) = reshape(x,shape(xxi(1:nx,1:ny,1:nz)))
      call setBC_ez(gv%gparams,igrid,xxi,bc_xi,iorder=order_bc)

c     Solve for vector potential

      call default_A_BCs(bcsa)
      call default_B_BCs(bcsb)

      !Homogeneous Dirichlet (lambda=0 @ Bdry)
      bc_lbda = bcond
      where(bc_lbda == DEF) bc_lbda = DIR

      !Vector-potential/lambda iteration
      lambda = 0d0 ; x = 0d0 ; xi = 0d0 ; ajcov = 0d0 ; ajcov_o = 0d0

      bcs(:,1) = bcond
      bcs(:,2) = bcond
      bcs(:,3) = bcond
      where(bcs == DEF) bcs = EXT

      call fpa_create(fpav,xa,maxv=10)

      !Compute lambda (with BCs) and vector potential
      if (my_rank == 0) write (*,*) " Hyper-res curl(A_J) solve..."

      lerror = 1d0 ; l1 = 1d0
      
      do it=1,30

        !Reconstruct lambda*V_J vector
        ecov = grad(gv%gparams,igx,xxi)
        ecnv = XformToCnv(gv%gparams,igx,ecov)
        ecnv = ecnv + curl(gv%gparams,igx,ajcov)

        call setBC_ez(gv%gparams,igx,ecnv,.true.,bcs,iorder=order_bc)

        !Compute lambda (with BCs)
        lambda_o = lambda
        lambda = scalarProduct_mesh(gv%gparams,igx,gJpar,ecnv)
     $       /(scalarProduct_mesh(gv%gparams,igx,gJpar,gJpar_cnv))
c$$$        lambda = max(lambda,0d0)
        
        call setBC_ez(gv%gparams,igrid,lambda,bc_lbda,iorder=order_bc)

c$$$        !Compute error in lambda
c$$$        x(:,1) = reshape(lambda  (1:nx,1:ny,1:nz)
c$$$     $                  -lambda_o(1:nx,1:ny,1:nz),shape(x(:,1)))
c$$$        lerror = rms(x(:,1))
c$$$
c$$$        if (it == 1) l0 = lerror
c$$$        if (my_rank==0) write (*,*) "Lambda it error, it=",it
c$$$     $                             ,",error=",lerror/l0
c$$$        if (lerror < rtol*l0.and.it > 1) exit

        if (lerror < rtol*l1.and.it > 1) exit

        !Solve for vector potential A_J
c$$$        ecnv(:,:,:,1) = (lambda-lambda_o)*gJpar_cnv(:,:,:,1)
c$$$        ecnv(:,:,:,2) = (lambda-lambda_o)*gJpar_cnv(:,:,:,2)
c$$$        ecnv(:,:,:,3) = (lambda-lambda_o)*gJpar_cnv(:,:,:,3)
        ecnv(:,:,:,1) = lambda*gJpar_cnv(:,:,:,1)
        ecnv(:,:,:,2) = lambda*gJpar_cnv(:,:,:,2)
        ecnv(:,:,:,3) = lambda*gJpar_cnv(:,:,:,3)

        ajcov_o = ajcov
        call vlap_curl_inv(gv%gparams,igx,ecnv,ajcov,bcsb,bcsa,ilevel
     .                    ,rtol=5d-1) !try .2

        !Apply fixed-point accelerator on A_J update
c$$$        xa = -reshape(ajcov(1:nx,1:ny,1:nz,:),shape(xa))
        xa = -reshape(ajcov  (1:nx,1:ny,1:nz,:)
     $               -ajcov_o(1:nx,1:ny,1:nz,:),shape(xa))

        call fpa_correction(fpav,xa,dp=dot2)

        ajcov(1:nx,1:ny,1:nz,:) = ajcov_o(1:nx,1:ny,1:nz,:)
     $         - reshape(xa,shape(ajcov  (1:nx,1:ny,1:nz,:)))

        call setMGBC(gv%gparams,0,3,nx,ny,nz,igx,ajcov,bcsa
     $              ,icomp=(/IAX/),is_vec=.true.,is_cnv=.false.
     .              ,iorder=order_bc)

        !Compute error in A_J
        lerror = rms(reshape(xa,(/3*nn/)))
        if (it == 1) l1 = max(lerror,findRoundOff())
        if (my_rank==0) write (*,100) "  >>> A_J error, it=",it
     $                             ,",error=",lerror/l1
 100    format(a,i3,a,1pe12.3)

      enddo

c diag ****
      if (my_rank == 0) then
        call createDrawInCfile(15,'lambda.bin','j','t','x','y'
     $                    ,(/'xi  ','Jpar','gJp1','gJp2','gJp3'
     $                      ,'A_1 ','A_2 ','A_3 ','lbda','hlm1'
     $                      ,'hlm2','hlm3','hlm1','hlm2','hlm3'/)
     $                    ,'-c -X0 -L57','drawlambda.in')
      endif

      !WARNING: this debug code requires global communication
      call dumpGlobalField(gv%gparams,1,"lambda.bin",xxi,0)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",JparAvg,1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",gJpar(:,:,:,1),1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",gJpar(:,:,:,2),1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",gJpar(:,:,:,3),1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",aJcov(:,:,:,1),1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",aJcov(:,:,:,2),1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",aJcov(:,:,:,3),1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",lambda,1)

      ecov = grad(gv%gparams,igx,xxi)
      ecnv = XformToCnv(gv%gparams,igx,ecov)
      ecnv = ecnv + curl(gv%gparams,igx,ajcov)

      call dumpGlobalField(gv%gparams,1,"lambda.bin",ecnv(:,:,:,1),1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",ecnv(:,:,:,2),1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",ecnv(:,:,:,3),1)

      ecnv(:,:,:,1) = lambda*gJpar_cnv(:,:,:,1)
      ecnv(:,:,:,2) = lambda*gJpar_cnv(:,:,:,2)
      ecnv(:,:,:,3) = lambda*gJpar_cnv(:,:,:,3)

      call dumpGlobalField(gv%gparams,1,"lambda.bin",ecnv(:,:,:,1),1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",ecnv(:,:,:,2),1)
      call dumpGlobalField(gv%gparams,1,"lambda.bin",ecnv(:,:,:,3),1)

c$$$        pause
c diag ****
        
      !Destroy fixed-point accelerator
      call fpa_destroy(fpav)

      deallocate(BAvg,JAvg,EAvg,etaAvg,JparAvg,RHSfluc,RHSsrc,xxi,diffcf
     $          ,gJpar,gJpar_cnv,ecnv,ecov,ajcov,ajcov_o,lambda_o)

      deallocate(x,y,lbda_err,xi,xa)
      
      end subroutine hyperres_diag

c     lmbda_mtvc
c     ####################################################################
      subroutine lmbda_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
c     --------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A.x  matrix-free
c     for the continuity equation.
c     In call:
c      * gpos: vector index of position on the numerical grid
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D
c              , 27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c     --------------------------------------------------------------------

      use setMGBC_interface

      use local_BCS_variables, ONLY: order_bc,nc_eom_v

      implicit none

c     Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)
      type(grid_mg_def),pointer :: g_def

c     Local variables

      integer    :: isig,ijk,ijkg,nnx,nny,nnz,ip,im,jp,jm,kp,km
      integer    :: imin,imax,jmin,jmax,kmin,kmax,i,j,k,ig,jg,kg
 
      real(8),allocatable,dimension(:,:,:,:) :: drho
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8),pointer    ,dimension(:,:,:) :: jac

      real(8)    :: upwind,slap,dvol

      logical    :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c     Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nnx = g_def%nxv(igrid)
      nny = g_def%nyv(igrid)
      nnz = g_def%nzv(igrid)

      jac => g_def%gmetric%grid(igrid)%jac

c     Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c     Map vector x to array for processing

      allocate(drho(0:nnx+1,0:nny+1,0:nnz+1,neq))

      drho = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped
      !For finding the diagonal, gpos > 0
      call fillArray(g_def,igrid,neq,x,drho,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IRHO/))

c     Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk = i + nnx*(j-1) + nnx*nny*(k-1)

            if (MG_get_res_vol_wgt()) then
              dvol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)*jac(i,j,k)
            else
              dvol = jac(i,j,k)
            endif

            !Conservative advection
            upwind = flx_advec(g_def,i,j,k,nnx,nny,nnz
     $                        ,igrid,igrid,igrid
     $                        ,v0_cnv,drho(:,:,:,1),UPWD)

c$$$            y(ijk) = dvol*upwind
c$$$            y(ijk)=dvol*(0.01*drho(i,j,k,1) + upwind)
            y(ijk)=dvol*(cnp/dt*drho(i,j,k,1) + alpha*upwind)

          enddo
        enddo
      enddo

c     End program

      deallocate(drho)
      nullify(v0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine lmbda_mtvc

c     lapneu_mtvc
c     ##################################################################
      subroutine lapneu_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
c     ------------------------------------------------------------------
c     This subroutine is a matvec test, y = A.x.
c     In call:
c      * gpos: vector index of position on the numerical grid
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D
c              , 27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnd: boundary conditions on x vector.
c      -----------------------------------------------------------------

      use setMGBC_interface

      use local_BCS_variables, ONLY: order_bc,nc_eom_v

      implicit none

c     Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)
      type(grid_mg_def),pointer :: g_def

c     Local variables

      integer :: i,j,k,ig,jg,kg,isig,ijk,ijkg,nxx,nyy,nzz
      integer :: imin,imax,jmin,jmax,kmin,kmax,order

      real(8),allocatable,dimension(:,:,:,:) :: xarr

      logical :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c     Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

c     Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c     Map vector x to array for processing

      allocate(xarr(0:nxx+1,0:nyy+1,0:nzz+1,neq))

cc      xarr = 0d0

      order = 0

      call fillArray(g_def,igrid,neq,x,xarr,bcnd,order,gpos=gpos)

c     Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

            y(ijk) = lap(g_def,i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                  ,xarr(:,:,:,1),vol=MG_get_res_vol_wgt()
     $                  ,zeroflx=.true.)

          enddo
        enddo
      enddo

c     End program

      deallocate(xarr)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine lapneu_mtvc
      
      end subroutine postProcessSolution
