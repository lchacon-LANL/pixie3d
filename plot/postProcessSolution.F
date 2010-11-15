#if defined(RFX)
#include "rfx_module"
#endif

c postProcessSolution
c####################################################################
      subroutine postProcessSolution(varray,iigx,iigy,iigz)

c--------------------------------------------------------------------
c     Postprocess solution in varray to find derived quantities 
c     for plotting.
c--------------------------------------------------------------------

      use variable_setup

      use graphics

      use auxiliaryVariables

      use auxPlotVariables

      use timeStepping

      use operators

      use equilibrium

      use nlfunction_setup

      use local_BCS_variables

      use diag_setup

#if defined(RFX)
      use rfx
#endif

      implicit none

c Call variables

      integer    :: iigx,iigy,iigz

      type(var_array),pointer :: varray

c Local variables

      integer    :: i,j,k,ig,jg,kg,ieq,bcs(6),igx,igy,igz,nx,ny,nz
      real(8)    :: x1,y1,z1,ee(3),vr,vth,norm,lnorm
     .             ,cov(3),cnv(3)
      logical    :: cartsn,covariant,to_cartsn,to_cnv

c Interfaces

      INTERFACE
        subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine imposeBoundaryConditions
      END INTERFACE

c Begin program

      igx = iigx
      igy = iigy
      igz = iigz

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      to_cartsn = .true.
      covariant = .false.
      to_cnv    = .false.

c Impose boundary conditions
c (finds all covariant and contravariant components of interest)

cc      write (*,*) 'Dumping perturbations of all quantities'
cc      call AXPYDerivedType(1d0,varray,-1d0,u_0,u_graph)
cc      varray = u_graph
cc      varray%array_var(1)%array = u_0%array_var(1)%array

cc      if (chk_src) then
cc        varray%array_var(IRHO)%array = 1d0
cc      endif

      call imposeBoundaryConditions(varray,igx,igy,igz)

      rho => varray%array_var(IRHO)%array
      px  => varray%array_var(IVX )%array
      py  => varray%array_var(IVY )%array
      pz  => varray%array_var(IVZ )%array
#if !defined(vec_pot)
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
#else
      ax  => varray%array_var(IAX )%array
      ay  => varray%array_var(IAY )%array
      az  => varray%array_var(IAZ )%array
#endif
      tmp => varray%array_var(ITMP)%array

#if defined(vec_pot)
      ax_cnv = acnv(:,:,:,1)
      ay_cnv = acnv(:,:,:,2)
      az_cnv = acnv(:,:,:,3)
#endif

c Find perturbed quantities (u_graph = varray - u_ic)

      call AXPYDerivedType(1d0,varray,-1d0,u_ic,u_graph)

c Assignments

      jx = jcnv(:,:,:,1)
      jy = jcnv(:,:,:,2)
      jz = jcnv(:,:,:,3)

      vx = vcnv(:,:,:,1)
      vy = vcnv(:,:,:,2)
      vz = vcnv(:,:,:,3)

      bx_cov = bcov(:,:,:,1)
      by_cov = bcov(:,:,:,2)
      bz_cov = bcov(:,:,:,3)

      jx_cov = jcov(:,:,:,1)
      jy_cov = jcov(:,:,:,2)
      jz_cov = jcov(:,:,:,3)

      vx_cov = vcov(:,:,:,1)
      vy_cov = vcov(:,:,:,2)
      vz_cov = vcov(:,:,:,3)

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
      call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,vx_car,vy_car,vz_car,'cnv','car')

#if defined(vec_pot)
      !Cartesian A components
      covariant = .true.
      ax_car = acov(:,:,:,1)
      ay_car = acov(:,:,:,2)
      az_car = acov(:,:,:,3)
      call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,ax_car,ay_car,az_car,'cov','car')
      covariant = .false.
#endif

      !Cartesian B components
      bx_car = bcnv(:,:,:,1)
      by_car = bcnv(:,:,:,2)
      bz_car = bcnv(:,:,:,3)
      call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,bx_car,by_car,bz_car,'cnv','car')

      !Find Cartesian J components
      jx_car = jcnv(:,:,:,1)
      jy_car = jcnv(:,:,:,2)
      jz_car = jcnv(:,:,:,3)
      call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,jx_car,jy_car,jz_car,'cnv','car')

c Return if checking sources

cc      if (chk_src) return

c Electron velocity

      vex = vecnv(:,:,:,1)
      vey = vecnv(:,:,:,2)
      vez = vecnv(:,:,:,3)

c Poloidal flux diagnostics

#if !defined(vec_pot)
      do k = 0,nz+1
        do i = 0,nx+1
          j = 0
          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
          if (i == 0) then
            Pflux(i,j,k) = 0d0
          else
            Pflux(i,j,k) = Pflux(i-1,j,k)
     .                   - 0.5*(bcnv(i-1,j,k,2)
     .                         +bcnv(i  ,j,k,2))*dx(ig-1)
          endif
          do j = 1,ny+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            Pflux(i,j,k) = Pflux(i,j-1,k)
     .                   + 0.5*(bcnv(i,j-1,k,1)
     .                         +bcnv(i,j  ,k,1))*dy(jg-1)
          enddo
        enddo
      enddo
#else
      Pflux = az
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

c Velocity and current intrinsic components (wrt B)

      do k = kming,kmaxg
        do j = jming,jmaxg
          do i = iming,imaxg
            call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

            !V_parallel
            norm = sqrt(vectorNorm(i,j,k,igx,igy,igz
     .                            ,bcnv(i,j,k,1)
     .                            ,bcnv(i,j,k,2)
     .                            ,bcnv(i,j,k,3),.false.))

            if (norm == 0d0) then
              vpar(i,j,k) = 0d0
            else
              vpar(i,j,k) = scalarProduct(i,j,k,igx,igy,igz
     .                     ,bcov(i,j,k,1),bcov(i,j,k,2),bcov(i,j,k,3)
     .                     ,vcnv(i,j,k,1),vcnv(i,j,k,2),vcnv(i,j,k,3))
     .                     /norm
            endif

            !V.grad(psi)/|grad(psi)|
            cov = grad(i,j,k,nx,ny,nz,igx,igy,igz,Pflux)
            norm = sqrt(vectorNorm(i,j,k,igx,igy,igz
     .                             ,cov(1),cov(2),cov(3),.true.))

            if (norm == 0d0) then
              vpsi(i,j,k) = 0d0
            else
              vpsi(i,j,k) = scalarProduct(i,j,k,igx,igy,igz
     .                     ,cov(1),cov(2),cov(3)
     .                     ,vcnv(i,j,k,1),vcnv(i,j,k,2),vcnv(i,j,k,3))
     .                     /norm
            endif
cc            vpsi(i,j,k) = cov(1)

            !V.(b x grad(psi))/|(b x grad(psi))|
            cnv(1) = bcov(i,j,k,2)*cov(3)
     .              -bcov(i,j,k,3)*cov(2)
            cnv(2) = bcov(i,j,k,3)*cov(1)
     .              -bcov(i,j,k,1)*cov(3)
            cnv(3) = bcov(i,j,k,1)*cov(2)
     .              -bcov(i,j,k,2)*cov(1)

            norm = sqrt(vectorNorm(i,j,k,igx,igy,igz
     .                            ,cnv(1),cnv(2),cnv(3),.false.))

            if (norm == 0d0) then
              vperp(i,j,k) = 0d0
            else
              vperp(i,j,k) = scalarProduct(i,j,k,igx,igy,igz
     .                       ,vcov(i,j,k,1),vcov(i,j,k,2),vcov(i,j,k,3)
     .                       ,cnv(1),cnv(2),cnv(3))
     .                       /norm
            endif
            
cc            !J_parallel
cc            norm = sqrt(vectorNorm(i,j,k,igx,igy,igz
cc     .                            ,bcnv(i,j,k,1)
cc     .                            ,bcnv(i,j,k,2)
cc     .                            ,bcnv(i,j,k,3),.false.))
cc
cc            if (norm == 0d0) then
cc              jpar(i,j,k) = 0d0
cc            else
cc              jpar(i,j,k) = scalarProduct(i,j,k,igx,igy,igz
cc     .                     ,bcov(i,j,k,1),bcov(i,j,k,2),bcov(i,j,k,3)
cc     .                     ,jcnv(i,j,k,1),jcnv(i,j,k,2),jcnv(i,j,k,3))
cc     .                     /norm
cc            endif
cc
cc            !J.grad(psi)/|grad(psi)|
cc            call grad(i,j,k,nx,ny,nz,igx,igy,igz,Pflux
cc     .               ,cov(1),cov(2),cov(3),0,0,0)
cc            norm = sqrt(vectorNorm(i,j,k,igx,igy,igz
cc     .                             ,cov(1),cov(2),cov(3),.true.))
cc
cc            if (norm == 0d0) then
cc              jpsi(i,j,k) = 0d0
cc            else
cc              jpsi(i,j,k) = scalarProduct(i,j,k,igx,igy,igz
cc     .                     ,cov(1),cov(2),cov(3)
cc     .                     ,jcnv(i,j,k,1),jcnv(i,j,k,2),jcnv(i,j,k,3))
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
cc            norm = sqrt(vectorNorm(i,j,k,igx,igy,igz
cc     .                            ,cnv(1),cnv(2),cnv(3),.false.))
cc
cc            if (norm == 0) then
cc              jperp(i,j,k) = 0d0
cc            else
cc              jperp(i,j,k) = scalarProduct(i,j,k,igx,igy,igz
cc     .                       ,jcov(i,j,k,1),jcov(i,j,k,2),jcov(i,j,k,3)
cc     .                       ,cnv(1),cnv(2),cnv(3))
cc     .                       /norm
cc            endif
          enddo
        enddo
      enddo

c Hall electric field components (covariant)

      if (di > 0d0) then
        jpar  = di*E_h(:,:,:,1)
        jpsi  = di*E_h(:,:,:,2)
        jperp = di*E_h(:,:,:,3)
      else
        jpar  = E_ni(:,:,:,1)
        jpsi  = E_ni(:,:,:,2)
        jperp = E_ni(:,:,:,3)
      endif

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
cc     .              /0.5/(gmetric%grid(1)%jac(i  ,j,k)
cc     .                   +gmetric%grid(1)%jac(i+1,j,k))
cc        ee(2) = eta*(-J0(2)+0.5*(jy_cov(i,j,k)+jy_cov(i+1,j,k)))
cc     .             - (0.5*(vz(i,j,k)+vz(i+1,j,k))
cc     .               *0.5*(bx(i,j,k)+bx(i+1,j,k))
cc     .               -0.5*(vx(i,j,k)+vx(i+1,j,k))
cc     .               *0.5*(bz(i,j,k)+bz(i+1,j,k)))
cc     .              /0.5/(gmetric%grid(1)%jac(i  ,j,k)
cc     .                   +gmetric%grid(1)%jac(i+1,j,k))
cc        ee(3) = eta*(-J0(3)+0.5*(jz_cov(i,j,k)+jz_cov(i+1,j,k)))
cc     .             - (0.5*(vx(i,j,k)+vx(i+1,j,k))
cc     .               *0.5*(by(i,j,k)+by(i+1,j,k))
cc     .               -0.5*(vy(i,j,k)+vy(i+1,j,k))
cc     .               *0.5*(bx(i,j,k)+bx(i+1,j,k)))
cc     .              /0.5/(gmetric%grid(1)%jac(i  ,j,k)
cc     .                   +gmetric%grid(1)%jac(i+1,j,k))
cc        write (*,*) 'Ohms law at',i,j,k,' = ',ee
cc      enddo

c Toroidal diagnostics (quasi-linear q-profile, lambda profile)

      if (     coords == 'hel'
     .    .or. coords == 'cyl'
     .    .or. coords == 'tor'
     .    .or. coords == 'vmc') then
        call tor_diag
      else
        qfactor = 0d0
        lambda  = 0d0
      endif

c Fourier output

      if (np == 1 .and. fourier) then
        call fourier_dump
      elseif (np > 1 .and. fourier .and. my_rank == 0) then
        write (*,*) 'Fourier dump not available in parallel'
      endif

c Poincare output (field tracer)

      if (poincare) call poincare_plot

c Divergence diagnostics

      divrgJ = div(igx,jcnv)
      divrgB = div(igx,bcnv)
      divrgV = div(igx,vcnv)

c Thermal pressure

      if (adiabatic) then
        p_tot = a_p*tmp     + 0.5*(bx_car**2+by_car**2+bz_car**2)
      else
        p_tot = a_p*rho*tmp + 0.5*(bx_car**2+by_car**2+bz_car**2)
      endif

c Transport coefficients

      do k = kming,kmaxg
        do j = jming,jmaxg
          do i = iming,imaxg
            nuu  (i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
            eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz,Pflux)
          enddo
        enddo
      enddo

#if defined(RFX)
c RFX dump
      if (np == 1) then
        call rfx_dump_fields
      elseif (np > 1 .and. my_rank == 0) then
        write (*,*) 'RFX fields dump not available in parallel'
      endif
#endif

c End program

      contains

c     tor_diag
c     ################################################################
      subroutine tor_diag

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
     .         ,v00(iming:imaxg,3),q00    (iming:imaxg,1)
     .         ,l00(iming:imaxg,1))

       if (    coords == 'hel') then
         mm = grid_params%params(1)
         kk = grid_params%params(2)
         RR = grid_params%params(3)
       elseif (coords == 'cyl') then
         mm = 1d0
         kk = 0d0
         RR = grid_params%params(3)
       else
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
             call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

             dv = dyh(jg)*dzh(kg)

             b00(i,1) = b00(i,1) + bcnv(i,j,k,1)*dv
             b00(i,2) = b00(i,2) + bcnv(i,j,k,2)*dv
             b00(i,3) = b00(i,3) + bcnv(i,j,k,3)*dv

             j00_cov(i,1) = j00_cov(i,1) + jcov(i,j,k,1)*dv
             j00_cov(i,2) = j00_cov(i,2) + jcov(i,j,k,2)*dv
             j00_cov(i,3) = j00_cov(i,3) + jcov(i,j,k,3)*dv

             v00(i,1) = v00(i,1) + vcov(i,j,k,1)*dv
             v00(i,2) = v00(i,2) + vcov(i,j,k,2)*dv
             v00(i,3) = v00(i,3) + vcov(i,j,k,3)*dv

cc                q00(i,1) = q00(i,1)
cc     .                + bcnv(i,j,k,3)/RR
cc     .                   /(bcnv(i,j,k,2)-kk/mm*bcnv(i,j,k,3))*dv
             l00(i,1) = l00(i,1)
     .                + scalarProduct(i,j,k,igx,igy,igz
     .                       ,jcov(i,j,k,1),jcov(i,j,k,2),jcov(i,j,k,3)
     .                       ,bcnv(i,j,k,1),bcnv(i,j,k,2),bcnv(i,j,k,3))
     .                  /vectorNorm(i,j,k,igx,igy,igz
     .                             ,bcnv(i,j,k,1)
     .                             ,bcnv(i,j,k,2)
     .                             ,bcnv(i,j,k,3),.false.)
     .                 *dv

             ll = ll + dv
           enddo
         enddo
       enddo

       !Communicate sums for each radial slice
#if defined(petsc)
       allocate(ldum(iming:imaxg,3))

       ldum = b00
       call MPI_Allreduce(ldum,b00,3*(imaxg-iming+1)
     .                   ,MPI_DOUBLE_PRECISION
     .                   ,MPI_SUM,MPI_COMM_YZ,mpierr)

       ldum = j00_cov
       call MPI_Allreduce(ldum,j00_cov,3*(imaxg-iming+1)
     .                   ,MPI_DOUBLE_PRECISION
     .                   ,MPI_SUM,MPI_COMM_YZ,mpierr)

       ldum = v00
       call MPI_Allreduce(ldum,v00,3*(imaxg-iming+1)
     .                   ,MPI_DOUBLE_PRECISION
     .                   ,MPI_SUM,MPI_COMM_YZ,mpierr)

       ldum(:,1) = l00(:,1)
       call MPI_Allreduce(ldum(:,1),l00(:,1),(imaxg-iming+1)
     .                   ,MPI_DOUBLE_PRECISION
     .                   ,MPI_SUM,MPI_COMM_YZ,mpierr)
       deallocate(ldum)

       dv = ll
       call MPI_Allreduce(dv,ll,1
     .                   ,MPI_DOUBLE_PRECISION
     .                   ,MPI_SUM,MPI_COMM_YZ,mpierr)
#endif

       !Average radial quantities
       b00     = b00    /ll
       j00_cov = j00_cov/ll
       v00     = v00    /ll
       q00     = q00    /ll
       l00     = l00    /ll

       !Find q-factor and lambda profile
       do i = iming,imaxg
         qfactor(i,:,:) = b00(i,3)/RR/(b00(i,2)-kk/mm*b00(i,3))
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

      end subroutine tor_diag

c     poincare_plot
c     ################################################################
      subroutine poincare_plot

c     ---------------------------------------------------------------
c     Integrates array(i,j,k) on periodic dimensions of domain
c     (nx)x(ny)x(nz) to find quasi-linear profiles.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      character(30) :: ft_file,ft_file_old=''
      integer :: uft=2345
      logical :: tor

      save :: ft_file_old

c     Begin program

c     Open file

      write (ft_file,'(1pe9.3)') time

      if (trim(ft_file) == trim(ft_file_old)) then
        return
      else
        ft_file_old = ft_file
      endif

      ft_file = 'poinc_t='//trim(ft_file)//'.bin'

#if defined(petsc)
      ft_file = trim(ft_file)//'_proc'//trim(int2char(my_rank))
#endif

      open(unit=uft,file=trim(ft_file)
     .    ,status='unknown',form='unformatted')

c     Dump arrays

      write (*,*) 'Dumping Poincare file ',ft_file

      tor = (equil == 'vmec' .or. coords == 'tor')

      write (uft) nxd,nyd,nzd,bcond,tor
cc      write (uft) grid_params%xx(1:nxl+2)
cc      write (uft) grid_params%yy(1:nyl+2)
cc      write (uft) grid_params%zz(1:nzl+2)
      write (uft) grid_params%xg  !Dump global position array
      write (uft) grid_params%yg  !"
      write (uft) grid_params%zg  !"

      if (np > 1) then
        write (uft) ilomg,ihipg,jlomg,jhipg,klomg,khipg
      endif

      write (uft) varray%array_var(IBX)%array
      write (uft) varray%array_var(IBY)%array
      write (uft) varray%array_var(IBZ)%array
      write (uft) bx_car
      write (uft) by_car
      write (uft) bz_car
      write (uft) gmetric%grid(1)%jac
      write (uft) gmetric%grid(1)%car

c     Close file

      close(uft)

c     End program

      end subroutine poincare_plot

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

      write (ft_file,'(1pe9.3)') 0.01*int(time*100)
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

      end subroutine postProcessSolution
