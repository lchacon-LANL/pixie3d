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

      use imposeBCinterface

      implicit none

c Call variables

      integer(4) :: iigx,iigy,iigz

      type(var_array) :: varray

c Local variables

      integer(4) :: i,j,k,ig,jg,kg
      real(8)    :: mm,kk,RR,ll,x1,y1,z1
      logical    :: cartsn,covariant,to_cartsn,to_cnv

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

      call imposeBoundaryConditions(varray,igx,igy,igz)

      rho => varray%array_var(IRHO)%array
      rvx => varray%array_var(IVX )%array
      rvy => varray%array_var(IVY )%array
      rvz => varray%array_var(IVZ )%array
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
      tmp => varray%array_var(ITMP)%array

c Find perturbed quantities

      u_graph = varray - u_ic

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

c Poloidally averaged q-factor (use graphics limits)

      if (coords == 'hel') then
        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        do k = kming,kmaxg
          do i = imaxg,iming+1,-1
            qfactor(i,jming-1,k) = 0d0
            ll = 0d0
            do j = jming,jmaxg
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
              qfactor(i,j,k) = qfactor(i,j-1,k)
     .             + mm*bz(i,j,k)/RR/(by(i,j,k)-kk*bz(i,j,k))*dyh(jg)
              ll = ll + dyh(jg)
            enddo
            qfactor(i,:,k) = qfactor(i,jmaxg,k)/ll
          enddo
          qfactor(iming,:,k) = qfactor(iming+1,:,k)
        enddo
      else
        qfactor = 0d0
      endif

c Divergence diagnostics

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            divrgJ(i,j,k) = div(i,j,k,nx,ny,nz,igx,igy,igz,jx,jy,jz)
            divrgB(i,j,k) = div(i,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
            divrgV(i,j,k) = div(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz)
          enddo
        enddo
      enddo

      call setBC(IRHO,nx,ny,nz,divrgJ,zeros,bcond,igx,igy,igz)
      call setBC(IRHO,nx,ny,nz,divrgB,zeros,bcond,igx,igy,igz)
      call setBC(IRHO,nx,ny,nz,divrgV,zeros,bcond,igx,igy,igz)

      !Average divergence around singular point (we consider whole control volume)
      if (bcond(1) == SP) then
        do k = 1,nz
          divrgJ(1,:,k) = divrgJ(0,:,k)
          divrgB(1,:,k) = divrgB(0,:,k)
          divrgV(1,:,k) = divrgV(0,:,k)
        enddo
      endif

c Total pressure (use graphics limits)

      p_tot = 0d0

      do k = kming,kmaxg
        do j = jming,jmaxg
          do i = iming,imaxg
            jac = gmetric%grid(igx)%jac(i,j,k)

            p_tot(i,j,k) = (2.*rho(i,j,k)*tmp(i,j,k))
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

      end subroutine postProcessSolution

c transformVector
c######################################################################
      subroutine transformVector(igx,igy,igz
     .                          ,imin,imax,jmin,jmax,kmin,kmax
     .                          ,arr1,arr2,arr3,covariant,to_cartsn)

c----------------------------------------------------------------------
c     Transforms vectors components in arrays arr1,arr2,arr3
c     from Cartesian to curvilinear (to_cartesian=.false.)
c     and viceversa, either with covariant (covariant=.true.) or 
c     contravariant curvilinear vectors.
c----------------------------------------------------------------------

      use grid

      implicit none

c     Input variables

        integer(4) :: imin,imax,jmin,jmax,kmin,kmax
        integer(4) :: igx,igy,igz
        logical    :: covariant,to_cartsn
        real(8)    :: arr1(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr2(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr3(imin:imax,jmin:jmax,kmin:kmax)

c     Local variables

        integer(4) :: i,j,k
        real(8)    :: vec(3)

c     Begin program

        if (to_cartsn) then

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call transformVectorToCartesian
     .               (i,j,k,igx,igy,igz
     .               ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .               ,covariant
     .               ,vec(1),vec(2),vec(3))

                arr1(i,j,k) = vec(1)
                arr2(i,j,k) = vec(2)
                arr3(i,j,k) = vec(3)
                
              enddo
            enddo
          enddo

        else

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call transformVectorToCurvilinear
     .               (i,j,k,igx,igy,igz
     .               ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .               ,covariant
     .               ,vec(1),vec(2),vec(3))

                arr1(i,j,k) = vec(1)
                arr2(i,j,k) = vec(2)
                arr3(i,j,k) = vec(3)
                
              enddo
            enddo
          enddo

        endif

      end subroutine transformVector
