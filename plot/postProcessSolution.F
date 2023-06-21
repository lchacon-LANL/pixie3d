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
     $                      ,poincare

#if defined(RFX)
      use rfx
#endif

      implicit none

c Call variables

      integer :: igrid,nn,my_rank_rad,np_rad

      type(var_array),pointer :: varray

c Local variables

      integer :: i,j,k,ig,jg,kg,ieq,bcs(6),igx,igy,igz,nx,ny,nz,ierr
      real(8) :: x1,y1,z1,ee(3),vr,vth,norm,lnorm
     .          ,cov(3),cnv(3)
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

c JxB-g(P) diagnostice (cnv)

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
     .                         ,jcnv,jcov,vcnv,vecnv)

              !External heat source
              heat_src(i,j,k) = heat_src(i,j,k)
     .                        + ext_heat_src(gv%gparams,i,j,k,igx)
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
      real(8),allocatable,dimension(:) :: mpivec
#endif

      integer :: np_mpi,my_rank_mpi
      
c     Begin program
      
      !Average fields toroidally
      allocate(bcnvavg(0:nx+1,0:ny+1,3))
      
      bcnvavg = sum(bcnv(:,:,1:nz,:),3)
#if defined(petsc)
      call MPI_Allreduce(bcnvavg,bcnvavg,size(bcnvavg)
     .                  ,MPI_DOUBLE_PRECISION,MPI_SUM
     .                  ,gv%gparams%MPI_COMM_Z,mpierr)
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
      
c     Local variables

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

      end subroutine postProcessSolution
