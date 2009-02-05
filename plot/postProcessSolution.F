#if defined(RFX)
      include "rfx_module"
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

      use constants

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

      integer    :: i,j,k,ig,jg,kg,ieq,bcs(6)
      real(8)    :: mm,kk,RR,ll,x1,y1,z1,ee(3),vr,vth,norm,lnorm
     .             ,cov(3),cnv(3)
      logical    :: cartsn,covariant,to_cartsn,to_cnv

      real(8),allocatable,dimension(:,:) :: b00,j00_cov,v00,q00,l00

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

      if(.not.chk_src) call imposeBoundaryConditions(varray,igx,igy,igz)

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

      vecnv(:,:,:,1) = vcnv(:,:,:,1) - di*jcnv(:,:,:,1)/rho
      vecnv(:,:,:,2) = vcnv(:,:,:,2) - di*jcnv(:,:,:,2)/rho
      vecnv(:,:,:,3) = vcnv(:,:,:,3) - di*jcnv(:,:,:,3)/rho

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

cc      bcs = bcond
cc      where (bcs == DEF) bcs = EXT

cc      call setBC(IRHO,nx,ny,nz,Pflux,zeros,bcs,igx,igy,igz)
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
            call grad(i,j,k,nx,ny,nz,igx,igy,igz,Pflux
     .               ,cov(1),cov(2),cov(3))
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
            
            !J_parallel
            norm = sqrt(vectorNorm(i,j,k,igx,igy,igz
     .                            ,bcnv(i,j,k,1)
     .                            ,bcnv(i,j,k,2)
     .                            ,bcnv(i,j,k,3),.false.))

            if (norm == 0d0) then
              jpar(i,j,k) = 0d0
            else
              jpar(i,j,k) = scalarProduct(i,j,k,igx,igy,igz
     .                     ,bcov(i,j,k,1),bcov(i,j,k,2),bcov(i,j,k,3)
     .                     ,jcnv(i,j,k,1),jcnv(i,j,k,2),jcnv(i,j,k,3))
     .                     /norm
            endif

            !J.grad(psi)/|grad(psi)|
            call grad(i,j,k,nx,ny,nz,igx,igy,igz,Pflux
     .               ,cov(1),cov(2),cov(3))
            norm = sqrt(vectorNorm(i,j,k,igx,igy,igz
     .                             ,cov(1),cov(2),cov(3),.true.))

            if (norm == 0d0) then
              jpsi(i,j,k) = 0d0
            else
              jpsi(i,j,k) = scalarProduct(i,j,k,igx,igy,igz
     .                     ,cov(1),cov(2),cov(3)
     .                     ,jcnv(i,j,k,1),jcnv(i,j,k,2),jcnv(i,j,k,3))
     .                     /norm
            endif
cc            vpsi(i,j,k) = cov(1)

            !J.(b x grad(psi))/|(b x grad(psi))|
            cnv(1) = bcov(i,j,k,2)*cov(3)
     .              -bcov(i,j,k,3)*cov(2)
            cnv(2) = bcov(i,j,k,3)*cov(1)
     .              -bcov(i,j,k,1)*cov(3)
            cnv(3) = bcov(i,j,k,1)*cov(2)
     .              -bcov(i,j,k,2)*cov(1)

            norm = sqrt(vectorNorm(i,j,k,igx,igy,igz
     .                            ,cnv(1),cnv(2),cnv(3),.false.))

            if (norm == 0) then
              jperp(i,j,k) = 0d0
            else
              jperp(i,j,k) = scalarProduct(i,j,k,igx,igy,igz
     .                       ,jcov(i,j,k,1),jcov(i,j,k,2),jcov(i,j,k,3)
     .                       ,cnv(1),cnv(2),cnv(3))
     .                       /norm
            endif
            
          enddo
        enddo
      enddo

c Velocity intrinsic components (wrt B)

      do k = kming,kmaxg
        do j = jming,jmaxg
          do i = iming,imaxg
            call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

          enddo
        enddo
      enddo

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

c Mean-field (poloidally averaged) q-factor (use graphics limits)

c FIX PARALLEL
      if (np == 1) then

        qfactor = 0d0
        if (    coords == 'hel'
     .     .or. coords == 'cyl'
     .     .or. coords == 'tor'
     .     .or. coords == 'vmc') then

          allocate(b00(iming:imaxg,3),j00_cov(iming:imaxg,3)
     .            ,v00(iming:imaxg,3),q00    (iming:imaxg,1)
     .            ,l00(iming:imaxg,1))

          if (    coords == 'hel'
     .       .or. coords == 'cyl'
     .       .or. coords == 'tor') then
            mm = grid_params%params(1)
            kk = grid_params%params(2)
            RR = grid_params%params(3)
          else !Define to avoid affecting q-profile
            mm = 1d0
            kk = 0d0
            RR = 1d0
          endif

          !Find mean fields
          do k = kming,kmaxg
            do i = iming,imaxg
              b00    (i,:) = 0d0
              j00_cov(i,:) = 0d0
              v00    (i,:) = 0d0
              q00    (i,:) = 0d0
              l00    (i,:) = 0d0
              ll = 0d0
              do j = jming,jmaxg-1  !Careful with the periodic limits
                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
                b00(i,1) = b00(i,1) + bcnv(i,j,k,1)*dyh(jg)
                b00(i,2) = b00(i,2) + bcnv(i,j,k,2)*dyh(jg)
                b00(i,3) = b00(i,3) + bcnv(i,j,k,3)*dyh(jg)
                j00_cov(i,1) = j00_cov(i,1) + jcov(i,j,k,1)*dyh(jg)
                j00_cov(i,2) = j00_cov(i,2) + jcov(i,j,k,2)*dyh(jg)
                j00_cov(i,3) = j00_cov(i,3) + jcov(i,j,k,3)*dyh(jg)
                v00(i,1) = v00(i,1) + vcov(i,j,k,1)*dyh(jg)
                v00(i,2) = v00(i,2) + vcov(i,j,k,2)*dyh(jg)
                v00(i,3) = v00(i,3) + vcov(i,j,k,3)*dyh(jg)
                if (coords == 'hel') then
                  q00(i,1) = q00(i,1)
     .                + bcnv(i,j,k,3)/RR
     .                   /(bcnv(i,j,k,2)-kk/mm*bcnv(i,j,k,3))*dyh(jg)
                else
                  q00(i,1) = q00(i,1)
     .                     + bcnv(i,j,k,3)/bcnv(i,j,k,2)*dyh(jg)
                endif
                l00(i,1) = l00(i,1)
     .                + scalarProduct(i,j,k,igx,igy,igz
     .                       ,jcov(i,j,k,1),jcov(i,j,k,2),jcov(i,j,k,3)
     .                       ,bcnv(i,j,k,1),bcnv(i,j,k,2),bcnv(i,j,k,3))
     .                  /vectorNorm(i,j,k,igx,igy,igz
     .                             ,bcnv(i,j,k,1)
     .                             ,bcnv(i,j,k,2)
     .                             ,bcnv(i,j,k,3),.false.)
     .                 *dyh(jg)

                ll = ll + dyh(jg)
              enddo
              b00    (i,:) = b00    (i,:)/ll
              j00_cov(i,:) = j00_cov(i,:)/ll
              v00    (i,:) = v00    (i,:)/ll
              q00    (i,:) = q00    (i,:)/ll
              l00    (i,:) = l00    (i,:)/ll
            enddo
          enddo

          !Find q-factor and lambda profile
          do k = kming,kmaxg
            do i = iming,imaxg
              do j = jming,jmaxg
                if (coords == 'hel') then
                  qfactor(i,j,k) = b00(i,3)/RR
     .                           /(b00(i,2)-kk/mm*b00(i,3))
                else
                  qfactor(i,j,k) = b00(i,3)/b00(i,2)
                endif
                lambda (i,j,k) = scalarProduct(i,j,k,igx,igy,igz
     .                         ,j00_cov(i,1),j00_cov(i,2),j00_cov(i,3)
     .                         ,b00    (i,1),b00    (i,2),b00    (i,3))
     .                          /vectorNorm(i,j,k,igx,igy,igz
     .                            ,b00(i,1),b00(i,2),b00(i,3),.false.)
                qfactr2(i,j,k) = q00(i,1)
                lambda2(i,j,k) = l00(i,1)
              enddo
            enddo
            qfactor(iming,jming:jmaxg,k)=qfactor(iming+1,jming:jmaxg,k)
            lambda (iming,jming:jmaxg,k)=lambda (iming+1,jming:jmaxg,k)
            qfactr2(iming,jming:jmaxg,k)=qfactr2(iming+1,jming:jmaxg,k)
            lambda2(iming,jming:jmaxg,k)=lambda2(iming+1,jming:jmaxg,k)
          enddo

c diag  : dump text data
          if (time == 0d0) then
            open(unit=1000,file='q-lambda.txt',status='unknown')
          else
            open(unit=1000,file='q-lambda.txt',status='unknown'
     .        ,position='append')
          endif
          write(1000,*) 'Time = ',time
          write(1000,*)'      Qlnr Q-prof    Avgd Q-prof    Qlnr Lambda'
     .                 ,'      Avgd Lambda         Vr'
          do i=iming,imaxg
            write (1000,*) qfactor(i,1,1),qfactr2(i,1,1)
     .                    ,lambda (i,1,1),lambda2(i,1,1),v00(i,1)
          enddo
          write (1000,*)
          close(1000)

#if defined(RFX)
          call rfx_dump_mean_profiles(b00,j00_cov,v00)
#endif
          deallocate(b00,j00_cov,v00,q00,l00)

        endif

      else
         if (my_rank == 0) then
           write (*,*) 'Q-profile and lambda-profile not'
     $                ,' available in parallel'
         endif
      endif
c FIX PARALLEL

c Fourier output

      if (np == 1 .and. fourier) then
        call fourier_dump
      elseif (np > 1 .and. fourier .and. my_rank == 0) then
        write (*,*) 'Fourier dump not available in parallel'
      endif

c Poincare output (field tracer)

      if (np == 1 .and. poincare) then
        call poincare_plot
      elseif (np > 1 .and. fourier .and. my_rank == 0) then
        write (*,*) 'Poincare dump not available in parallel'
      endif

c Divergence diagnostics

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            divrgJ(i,j,k) = div(i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,jcnv(:,:,:,1)
     .                         ,jcnv(:,:,:,2)
     .                         ,jcnv(:,:,:,3))
            divrgB(i,j,k) = div(i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,bcnv(:,:,:,1)
     .                         ,bcnv(:,:,:,2)
     .                         ,bcnv(:,:,:,3))
            divrgV(i,j,k) = div(i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,vcnv(:,:,:,1)
     .                         ,vcnv(:,:,:,2)
     .                         ,vcnv(:,:,:,3))
          enddo
        enddo
      enddo

c Thermal pressure (use graphics limits)

      p_tot = a_p*rho*tmp

c Transport coefficients

      do k = kming,kmaxg
        do j = jming,jmaxg
          do i = iming,imaxg
            nuu  (i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
            eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

#if defined(RFX)
c RFX dump
      call rfx_dump_fields
#endif

c End program

      contains

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

cc      write (ft_file,'(1pe9.3)') 0.01*int(time*100)
      write (ft_file,'(1pe9.3)') time

      if (trim(ft_file) == trim(ft_file_old)) then
        return
      else
        ft_file_old = ft_file
      endif

      ft_file = 'poinc_t='//trim(ft_file)//'.bin'

      open(unit=uft,file=trim(ft_file)
     .    ,status='unknown',form='unformatted')

c     Dump arrays

      write (*,*) 'Dumping Poincare file ',ft_file

      tor = (coords == 'vmc' .or. coords == 'tor')

      write (uft) nxd,nyd,nzd,bcond,tor
      write (uft) grid_params%xx(1:nxd+2)
      write (uft) grid_params%yy(1:nyd+2)
      write (uft) grid_params%zz(1:nzd+2)
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

ccc     ql_integral
ccc     ################################################################
cc      function ql_integral(nx,ny,nz,array,igx,igy,igz,avg) result(qlp)
cc
ccc     ---------------------------------------------------------------
ccc     Integrates array(i,j,k) on periodic dimensions of domain
ccc     (nx)x(ny)x(nz) to find quasi-linear profiles.
ccc     ---------------------------------------------------------------
cc
cc      implicit none
cc
ccc    Call variables
cc
cc      integer    :: igx,igy,igz,nx,ny,nz
cc      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1),integral
cc     .             
cc      logical    :: avg
cc
ccc     Local variables
cc
cc      integer    :: i,j,k
cc
cc      real(8)    :: tvolume,vol,lvolume,lintegral
cc
ccc     Begin program
cc
ccc     Integrate
cc
cc      lintegral = 0d0
cc      lvolume   = 0d0
cc
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            vol = volume(i,j,k,igx,igy,igz)
cc
cc            if (isSYM(i,igx,1,0)) vol = 0.5*vol
cc            if (isSYM(i,igx,1,1)) vol = 0.5*vol
cc            if (isSYM(j,igy,2,0)) vol = 0.5*vol
cc            if (isSYM(j,igy,2,1)) vol = 0.5*vol
cc            if (isSYM(k,igz,3,0)) vol = 0.5*vol
cc            if (isSYM(k,igz,3,1)) vol = 0.5*vol
cc
cc            lintegral = lintegral + array(i,j,k)*vol
cc            lvolume = lvolume + vol
cc          enddo
cc        enddo
cc      enddo
cc
cc#if defined(petsc)
cc      call MPI_Allreduce(lintegral,integral,1,MPI_DOUBLE_PRECISION
cc     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
cc      call MPI_Allreduce(lvolume  ,tvolume ,1,MPI_DOUBLE_PRECISION
cc     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
cc#else
cc      integral = lintegral
cc      tvolume  = lvolume
cc#endif
cc
cc      if (avg) integral = integral/tvolume
cc
ccc     End 
cc
cc      end function ql_integral

      end subroutine postProcessSolution
