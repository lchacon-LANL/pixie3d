c postProcessSolution
c####################################################################
      subroutine postProcessSolution(varray,iigx,iigy,iigz)

c--------------------------------------------------------------------
c     Set graphics files and dumping intervals
c--------------------------------------------------------------------

      use variables

      use graphics

      use auxiliaryVariables

      use timeStepping

      use constants

      use operators

      use equilibrium

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

      u_graph = varray - u_0

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
            divrgJ(i,j,k) = div(i,j,k,nx,ny,nz,jx,jy,jz)
            divrgB(i,j,k) = div(i,j,k,nx,ny,nz,bx,by,bz)
            divrgV(i,j,k) = div(i,j,k,nx,ny,nz,vx,vy,vz)
          enddo
        enddo
      enddo

      call imposeBConScalar(IRHO,divrgJ,zeros,bcond)
      call imposeBConScalar(IRHO,divrgB,zeros,bcond)
      call imposeBConScalar(IRHO,divrgV,zeros,bcond)

      !Average divergence around singular point (we consider whole control volume)
      if (bcond(1) == SP) then
        do k = 1,nz
          divrgJ(1,:,k) = sum(divrgJ(1,1:ny,k))/ny
          divrgB(1,:,k) = sum(divrgB(1,1:ny,k))/ny
          divrgV(1,:,k) = sum(divrgV(1,1:ny,k))/ny
        enddo
      endif

c Total pressure (use graphics limits)

      p_tot = 0d0

      do k = kming,kmaxg
        do j = jming,jmaxg
          do i = iming,imaxg
            call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartsn)
            jac = jacobian(x1,y1,z1,cartsn)
            
            p_tot(i,j,k) = (2.*rho(i,j,k)*tmp(i,j,k))
     .                     +(bx(i,j,k)*bx_cov(i,j,k)
     .                      +by(i,j,k)*by_cov(i,j,k)
     .                      +bz(i,j,k)*bz_cov(i,j,k))/2./jac

          enddo
        enddo
      enddo

c End program

      end subroutine postProcessSolution
