c postProcessSolution
c####################################################################
      subroutine postProcessSolution(varray,iigx,iigy,iigz)

c--------------------------------------------------------------------
c     Postprocess solution in varray to find derived quantities 
c     for plotting.
c--------------------------------------------------------------------

      use variables

      use graphics

      use auxiliaryVariables

      use timeStepping

      use constants

      use operators

      use equilibrium

      use nlfunction_setup

      use vectorOps

      use imposeBCinterface

      implicit none

c Call variables

      integer(4) :: iigx,iigy,iigz

      type(var_array) :: varray

c Local variables

      integer(4) :: i,j,k,ig,jg,kg,ieq
      real(8)    :: mm,kk,RR,ll,x1,y1,z1,ee(3)
      logical    :: cartsn,covariant,to_cartsn,to_cnv

      real(8),allocatable,dimension(:,:) :: b00,j00_cov,v00,q00,l00

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
cc      call substractDerivedType(varray,u_0,u_graph)
cc      varray = u_graph
cc      varray%array_var(1)%array = u_0%array_var(1)%array

      call imposeBoundaryConditions(varray,igx,igy,igz)

      rho => varray%array_var(IRHO)%array
      px => varray%array_var(IVX )%array
      py => varray%array_var(IVY )%array
      pz => varray%array_var(IVZ )%array
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
      tmp => varray%array_var(ITMP)%array

c Find perturbed quantities (u_graph = varray - u_ic)

      call substractDerivedType(varray,u_ic,u_graph)

c Find Cartesian components of ALL vectors

      !Cartesian velocities
      vx_car = vx
      vy_car = vy
      vz_car = vz
      call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,vx_car,vy_car,vz_car,covariant,to_cartsn)

      !Cartesian B components
      bx_car = bx
      by_car = by
      bz_car = bz
      call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,bx_car,by_car,bz_car,covariant,to_cartsn)

      !Find Cartesian J components

      jx_car = jx
      jy_car = jy
      jz_car = jz
      call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,jx_car,jy_car,jz_car,covariant,to_cartsn)

c Electron velocity

      vex = vx - di*jx/rho
      vey = vy - di*jy/rho
      vez = vz - di*jz/rho

c Poloidal flux diagnostics  (use graphics limits)

      do k = kming,kmaxg
        do i = iming,imaxg
          j = jming
          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
          if (i == iming) then
            Pflux(i,j,k) = 0d0
          else
            Pflux(i,j,k) = Pflux(i-1,j,k) - by(i,j,k)*dxh(ig)
          endif
          do j = jming+1,jmaxg
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            Pflux(i,j,k) = Pflux(i,j-1,k) + bx(i,j,k)*dyh(jg)
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
      qfactor = 0d0
      if (coords == 'hel' .or. coords == 'cyl') then

        allocate(b00(iming:imaxg,3),j00_cov(iming:imaxg,3)
     .          ,v00(iming:imaxg,3),q00    (iming:imaxg,1)
     .          ,l00(iming:imaxg,1))

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)

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
              b00(i,1) = b00(i,1) + bx(i,j,k)*dyh(jg)
              b00(i,2) = b00(i,2) + by(i,j,k)*dyh(jg)
              b00(i,3) = b00(i,3) + bz(i,j,k)*dyh(jg)
              j00_cov(i,1) = j00_cov(i,1) + jx_cov(i,j,k)*dyh(jg)
              j00_cov(i,2) = j00_cov(i,2) + jy_cov(i,j,k)*dyh(jg)
              j00_cov(i,3) = j00_cov(i,3) + jz_cov(i,j,k)*dyh(jg)
              v00(i,1) = v00(i,1) + vx_cov(i,j,k)*dyh(jg)
              q00(i,1) = q00(i,1)
     .              + mm*bz(i,j,k)/RR/(by(i,j,k)-kk*bz(i,j,k))*dyh(jg)
              l00(i,1) = l00(i,1)
     .              + scalarProduct(i,j,k,igx,igy,igz
     .                     ,jx_cov(i,j,k),jy_cov(i,j,k),jz_cov(i,j,k)
     .                     ,bx    (i,j,k),by    (i,j,k),bz    (i,j,k))
     .                /vectorNorm(i,j,k,igx,igy,igz
     .                          ,bx(i,j,k),by(i,j,k),bz(i,j,k),.false.)
     .               *dyh(jg)

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
              qfactor(i,j,k) = mm*b00(i,3)/RR/(b00(i,2)-kk*b00(i,3))
              lambda (i,j,k) = scalarProduct(i,j,k,igx,igy,igz
     .                          ,j00_cov(i,1),j00_cov(i,2),j00_cov(i,3)
     .                          ,b00    (i,1),b00    (i,2),b00    (i,3))
     .                        /vectorNorm(i,j,k,igx,igy,igz
     .                          ,b00(i,1),b00(i,2),b00(i,3),.false.)
              qfactr2(i,j,k) = q00(i,1)
              lambda2(i,j,k) = l00(i,1)
            enddo
          enddo
          qfactor(iming,:,k) = qfactor(iming+1,:,k)
          lambda (iming,:,k) = lambda (iming+1,:,k)
          qfactr2(iming,:,k) = qfactr2(iming+1,:,k)
          lambda2(iming,:,k) = lambda2(iming+1,:,k)
        enddo

c diag: dump text data
        if (time == 0d0) then
          open(unit=1000,file='q-lambda.txt',status='unknown')
        else
          open(unit=1000,file='q-lambda.txt',status='unknown'
     .      ,access='append')
        endif
        write (1000,*) 'Time = ',time
        write (1000,*) '      Qlnr Q-prof    Avgd Q-prof    Qlnr Lambda'
     .                ,'      Avgd Lambda         Vr'
        do i=iming,imaxg
          write (1000,*) qfactor(i,1,1),qfactr2(i,1,1)
     .                  ,lambda (i,1,1),lambda2(i,1,1),v00(i,1)
        enddo
        write (1000,*)
        close(1000)
c diag: dump text data

        deallocate(b00,j00_cov,v00)

c diag: dump text data
cc        if (time == 0d0) then
          open(unit=1000,file='fourier.txt',status='unknown')
cc        else
cc          open(unit=1000,file='q-lambda.txt',status='unknown'
cc     .      ,access='append')
cc        endif
          write (*,*) 'Dumping fourier.txt'
        do i=1,nxd
          do j=1,nyd
            do k=1,nzd
              write (1000,*) bx_car(i,j,k),by_car(i,j,k),bz_car(i,j,k)
            enddo
          enddo
        enddo
        write (1000,*)
        close(1000)
c diag: dump text data
      endif
c FIX PARALLEL

c Divergence diagnostics

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
cc            jac = gmetric%grid(igx)%jac(i,j,k)
            jac = 1d0
            divrgJ(i,j,k) = jac*div(i,j,k,nx,ny,nz,igx,igy,igz,jx,jy,jz)
            divrgB(i,j,k) = jac*div(i,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
            divrgV(i,j,k) = jac*div(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz)
          enddo
        enddo
      enddo

      call setBC(IRHO,nx,ny,nz,divrgJ,zeros,bcond,igx,igy,igz)
      call setBC(IRHO,nx,ny,nz,divrgB,zeros,bcond,igx,igy,igz)
      call setBC(IRHO,nx,ny,nz,divrgV,zeros,bcond,igx,igy,igz)

c Total pressure (use graphics limits)

      p_tot = 0d0

      do k = kming,kmaxg
        do j = jming,jmaxg
          do i = iming,imaxg
            jac = gmetric%grid(igx)%jac(i,j,k)

            p_tot(i,j,k) = (a_p*rho(i,j,k)*tmp(i,j,k))
     .                     +(bx(i,j,k)*bx_cov(i,j,k)
     .                      +by(i,j,k)*by_cov(i,j,k)
     .                      +bz(i,j,k)*bz_cov(i,j,k))/2./jac

          enddo
        enddo
      enddo

c Transport coefficients

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            nuu  (i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
            eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

c End program

cc      contains
cc
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
cc      integer(4) :: igx,igy,igz,nx,ny,nz
cc      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1),integral
cc     .             
cc      logical    :: avg
cc
ccc     Local variables
cc
cc      integer(4) :: i,j,k
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
