c defineBoundaryConditions
c####################################################################
      subroutine defineBoundaryConditions (neq,bcs)
c--------------------------------------------------------------------
c     Defines boundary conditions of physical quantities.
c     On input:
c       * neq -> number of equations
c     On output:
c       * bcs -> real array of size (6,neq) containing BC setup:
c           + bcs(1) ---> at x0
c           + bcs(2) ---> at x1
c           + bcs(3) ---> at y0
c           + bcs(4) ---> at y1
c           + bcs(5) ---> at z0
c           + bcs(6) ---> at z1
c--------------------------------------------------------------------

      use equilibrium

      use grid

      use icond

      implicit none

c Call variables

      integer(4) :: neq,bcs(6,neq)

c Local variables

      integer(4) :: ieq,bcsq(6)

c Begin program

c Defaults

      bcsq = bcs(:,IRHO)
      where (bcsq == DEF) bcsq = NEU
      bcs(:,IRHO) = bcsq

      bcsq = bcs(:,IVX)
      where (bcsq == DEF) bcsq = DIR
      bcs(:,IVX) = bcsq

      bcsq = bcs(:,IVY)
      where (bcsq == DEF) bcsq = NEU
      bcs(:,IVY) = bcsq

      bcsq = bcs(:,IVZ)
      where (bcsq == DEF) bcsq = NEU
      bcs(:,IVZ) = bcsq

      bcsq = bcs(:,IBX)
      where (bcsq == DEF) bcsq = DIR
      bcs(:,IBX) = bcsq
cc      if (bcs(1,IBX) == NEU) bcs(1,IBX) = DIR
cc      if (bcs(2,IBX) == NEU) bcs(2,IBX) = DIR

      bcsq = bcs(:,IBY)
      where (bcsq == DEF) bcsq = NEU
      bcs(:,IBY) = bcsq
cc      if (bcs(3,IBY) == NEU) bcs(3,IBY) = DIR
cc      if (bcs(4,IBY) == NEU) bcs(4,IBY) = DIR

      bcsq = bcs(:,IBZ)
      where (bcsq == DEF) bcsq = NEU
      bcs(:,IBZ) = bcsq
cc      if (bcs(5,IBZ) == NEU) bcs(5,IBZ) = DIR
cc      if (bcs(6,IBZ) == NEU) bcs(6,IBZ) = DIR

      bcsq = bcs(:,ITMP)
      where (bcsq == DEF) bcsq = NEU !To allow isothermal case
      bcs(:,ITMP) = bcsq

c Exceptions

      select case (equil)

      case ('rfp')

cc        bcs(2,IBY) = EQU  !Imposed by equilibrium

cc        bcs(2,IBZ) = EQU  !Imposed by equilibrium

      end select

c End

      end subroutine

c imposeBoundaryConditions
c####################################################################
      subroutine imposeBoundaryConditions (varray)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use equilibrium

      use grid

      use variables

      use nlfunction_setup

      use icond

      use constants

      implicit none

c Call variables

      type (var_array) :: varray

c Local variables

      integer(4)       :: neq,ieq,i,j,k,ig,jg,kg
      integer(4)       :: dim,loc,bctype,ibc

c Begin program

      rho => varray%array_var(IRHO)%array
      rvx => varray%array_var(IVX )%array
      rvy => varray%array_var(IVY )%array
      rvz => varray%array_var(IVZ )%array
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
      tmp => varray%array_var(ITMP)%array

c Find covariant magnetic field

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,bx_cov(i,j,k),by_cov(i,j,k),bz_cov(i,j,k)
     .             ,bx    (i,j,k),by    (i,j,k),bz    (i,j,k),.false.)
          enddo
        enddo
      enddo

c Fill ghost nodes

      neq = varray%nvar

      do bctype=1,5            !Enforces a particular order in the BCs (see grid_mod.f)
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
            do ieq = 1,neq 
              if (varray%array_var(ieq)%bconds(ibc) == bctype) then
                if (bctype == EQU) then
                  call FillGhostNodes(ieq,dim,loc,bctype
     .                               ,varray%array_var(ieq)%array
     .                               ,u_0   %array_var(ieq)%array)
                else
                  call FillGhostNodes(ieq,dim,loc,bctype
     .                               ,varray%array_var(ieq)%array
     .                               ,zeros)
                endif
              endif
            enddo
          enddo
        enddo
      enddo

c Impose vector singular point BCs

      if (bcond(1) == SP) then
        call vectorSingularBC(rvx,rvy,rvz,.false.,2)
        call vectorSingularBC(bx ,by ,bz ,.false.,2)
      endif

c Synchronize periodic boundaries

      bctype=PER

      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)
          do ieq = 1,neq 
            if (varray%array_var(ieq)%bconds(ibc) == bctype) then
              call FillGhostNodes(ieq,dim,loc,bctype
     .                           ,varray%array_var(ieq)%array
     .                           ,zeros)
            endif
          enddo
        enddo
      enddo

c Synchronize velocity field

      call postProcessV

c Synchronize magnetic field

      call postProcessB

c Calculate current

      call postProcessJ

c End

      contains

c     postProcessV
c     #################################################################
      subroutine postProcessV
c     -----------------------------------------------------------------
c     Synchronizes velocity components (covariant, contravariant)
c     at boundaries.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      integer(4) :: i,j,k,dim,loc
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      real(8)    :: x1,x2,x3,gsub(3,3)
      logical    :: cartesian

      integer(4) :: ic,ig,jg,kg
      real(8)    :: x,y,z,xp,yp,zp
      real(8)    :: avg_q,avg_vol,cx,cy,cz
      real(8),allocatable,dimension(:) :: ax0,ay0,az0

c     Begin program

c     Find velocity field

      where (rho /= 0d0)
        vx = rvx/rho
        vy = rvy/rho
        vz = rvz/rho
      end where

c     End program

      end subroutine postProcessV

c     postProcessB
c     #################################################################
      subroutine postProcessB
c     -----------------------------------------------------------------
c     Synchronizes magnetic field components (covariant, contravariant)
c     at boundaries.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      integer(4) :: i,j,k,dim,loc
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      real(8)    :: x1,x2,x3,gsuper(3,3),gsub(3,3)
      logical    :: cartesian,cov_to_cnv

      integer(4) :: ic,ig,jg,kg
      real(8)    :: x,y,z,xp,yp,zp
      real(8)    :: avg_q,avg_vol,cx,cy,cz
      real(8),allocatable,dimension(:) :: ax0,ay0,az0

c     Begin program

c     Determine postprocessing operation

      cov_to_cnv = .false.
      do dim = 1,3
        do loc = 0,1
          ibc = (1+loc)+2*(dim-1)
          if (    varray%array_var(IBX)%bconds(ibc) == NEU
     .        .or.varray%array_var(IBY)%bconds(ibc) == NEU
     .        .or.varray%array_var(IBZ)%bconds(ibc) == NEU) then
            cov_to_cnv = .true.
          endif
        enddo
      enddo

c     FIND CONTRAVARIANT COMPONENTS AT BOUNDARY
      if (cov_to_cnv) then

c     Find covariant magnetic field NORMAL components at boundaries

      do dim = 1,3
        do loc = 0,1
          ibc = (1+loc)+2*(dim-1)

          if (dim == 1) then
            imin=1  +    loc *(nx-1)
            imax=nx + (1-loc)*(1-nx)
            jmin=1
            jmax=ny
            kmin=1
            kmax=nz
          elseif (dim == 2) then
            imin=1 
            imax=nx
            jmin=1  +    loc *(ny-1)
            jmax=ny + (1-loc)*(1-ny)
            kmin=1
            kmax=nz
          elseif (dim == 3) then
            imin=1 
            imax=nx
            jmin=1
            jmax=ny
            kmin=1  +    loc *(nz-1)
            kmax=nz + (1-loc)*(1-nz)
          endif

          select case (ibc)
          case (1)

            if (bcond(1) == SP) then

              call vectorSingularBC(bx_cov,by_cov,bz_cov,.true.,2)

            else
              do i=imin,imax
                do j=jmin,jmax
                  do k=kmin,kmax

                    call getCoordinates(i-1,j,k,igx,igy,igz,ig,jg,kg
     .                                 ,x1,x2,x3,cartesian)

                    gsuper = g_super(x1,x2,x3,cartesian)

                    bx_cov(i-1,j,k) = -(gsuper(1,2)*by_cov(i-1,j,k)
     .                                 +gsuper(1,3)*bz_cov(i-1,j,k)
     .                                 -bx(i-1,j,k) )/gsuper(1,1)

                  enddo
                enddo
              enddo
            endif

          case (2)

            do i=imin,imax
              do j=jmin,jmax
                do k=kmin,kmax

                  call getCoordinates(i+1,j,k,igx,igy,igz,ig,jg,kg
     .                               ,x1,x2,x3,cartesian)

                  gsuper = g_super(x1,x2,x3,cartesian)

                  bx_cov(i+1,j,k) = -(gsuper(1,2)*by_cov(i+1,j,k)
     .                               +gsuper(1,3)*bz_cov(i+1,j,k)
     .                               -bx(i+1,j,k) )/gsuper(1,1)
                enddo
              enddo
            enddo

          case (3)

            do i=imin,imax
              do j=jmin,jmax
                do k=kmin,kmax

                  call getCoordinates(i,j-1,k,igx,igy,igz,ig,jg,kg
     .                               ,x1,x2,x3,cartesian)

                  gsuper = g_super(x1,x2,x3,cartesian)

                  by_cov(i,j-1,k) = -(gsuper(2,1)*bx_cov(i,j-1,k)
     .                               +gsuper(2,3)*bz_cov(i,j-1,k)
     .                               -by(i,j-1,k) )/gsuper(2,2)

                enddo
              enddo
            enddo

          case (4)

            do i=imin,imax
              do j=jmin,jmax
                do k=kmin,kmax

                  call getCoordinates(i,j+1,k,igx,igy,igz,ig,jg,kg
     .                               ,x1,x2,x3,cartesian)

                  gsuper = g_super(x1,x2,x3,cartesian)

                  by_cov(i,j+1,k) = -(gsuper(2,1)*bx_cov(i,j+1,k)
     .                               +gsuper(2,3)*bz_cov(i,j+1,k)
     .                               -by(i,j+1,k) )/gsuper(2,2)

                enddo
              enddo
            enddo

          case (5)

            do i=imin,imax
              do j=jmin,jmax
                do k=kmin,kmax

                  call getCoordinates(i,j,k-1,igx,igy,igz,ig,jg,kg
     .                               ,x1,x2,x3,cartesian)

                  gsuper = g_super(x1,x2,x3,cartesian)

                  bz_cov(i,j,k-1) = -(gsuper(3,1)*bx_cov(i,j,k-1)
     .                               +gsuper(3,2)*by_cov(i,j,k-1)
     .                               -bz(i,j,k-1) )/gsuper(3,3)
                enddo
              enddo
            enddo

          case (6)

            do i=imin,imax
              do j=jmin,jmax
                do k=kmin,kmax

                  call getCoordinates(i,j,k+1,igx,igy,igz,ig,jg,kg
     .                               ,x1,x2,x3,cartesian)

                  gsuper = g_super(x1,x2,x3,cartesian)

                  bz_cov(i,j,k+1) = -(gsuper(3,1)*bx_cov(i,j,k+1)
     .                               +gsuper(3,2)*by_cov(i,j,k+1)
     .                               -bz(i,j,k+1) )/gsuper(3,3)
                enddo
              enddo
            enddo

          end select

        enddo
      enddo

c     Enforce PER BC on covariant magnetic field components

      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          bctype=PER

          ieq = IBX
          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
            call FillGhostNodes(ieq,dim,loc,bctype,bx_cov,zeros)
          endif

          ieq = IBY
          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
            call FillGhostNodes(ieq,dim,loc,bctype,by_cov,zeros)
          endif

          ieq = IBZ
          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
            call FillGhostNodes(ieq,dim,loc,bctype,bz_cov,zeros)
          endif

        enddo
      enddo

c     Find contravariant magnetic field at boundaries

      do dim = 1,3
        do loc = 0,1
          if (dim == 1) then
            imin=0 + loc*(nx+1)
            imax=0 + loc*(nx+1)
            jmin=0
            jmax=ny+1
            kmin=0
            kmax=nz+1
          elseif (dim == 2) then
            imin=0
            imax=nx+1
            jmin=0 + loc*(ny+1)
            jmax=0 + loc*(ny+1)
            kmin=0
            kmax=nz+1
          elseif (dim == 3) then
            imin=0
            imax=nx+1
            jmin=0
            jmax=ny+1
            kmin=0 + loc*(nz+1)
            kmax=0 + loc*(nz+1)
          endif

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax
                call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,bx_cov(i,j,k),by_cov(i,j,k),bz_cov(i,j,k)
     .             ,bx(i,j,k),by(i,j,k),bz(i,j,k),.true.)
              enddo
            enddo
          enddo

        enddo
      enddo

c     FIND COVARIANT COMPONENTS AT BOUNDARY
      else

      do dim = 1,3
        do loc = 0,1
          if (dim == 1) then
            imin=0 + loc*(nx+1)
            imax=0 + loc*(nx+1)
            jmin=0
            jmax=ny+1
            kmin=0
            kmax=nz+1
          elseif (dim == 2) then
            imin=0
            imax=nx+1
            jmin=0 + loc*(ny+1)
            jmax=0 + loc*(ny+1)
            kmin=0
            kmax=nz+1
          elseif (dim == 3) then
            imin=0
            imax=nx+1
            jmin=0
            jmax=ny+1
            kmin=0 + loc*(nz+1)
            kmax=0 + loc*(nz+1)
          endif

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax
                call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,bx_cov(i,j,k),by_cov(i,j,k),bz_cov(i,j,k)
     .             ,bx(i,j,k),by(i,j,k),bz(i,j,k),.false.)
              enddo
            enddo
          enddo

        enddo
      enddo

      endif

c     End program

      end subroutine postProcessB

c     postProcessJ
c     #################################################################
      subroutine postProcessJ
c     -----------------------------------------------------------------
c     Calculates current components (covariant, contravariant).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      integer(4) :: i,j,k,dim,loc,order
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      real(8)    :: x1,x2,x3,gsub(3,3)
      real(8),allocatable,dimension(:) :: xint,vint
      real(8)                          :: xo(2),vo(2)
      logical    :: cartesian

      integer(4) :: ic,ig,jg,kg
      real(8)    :: x,y,z,xp,yp,zp
      real(8)    :: avg_q,avg_vol,cx,cy,cz
      real(8),allocatable,dimension(:) :: ax0,ay0,az0

c     Externals

      real(8)    :: quad_int
      external   :: quad_int

c     Begin program

c     Contravariant current

c This seems to create an instability at the dirichlet boundary
ccheck      do k = 0,nz+1
ccheck        do j = 0,ny+1
ccheck         do i = 0,nx+1
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            jx(i,j,k) = curl(i,j,k,bx_cov,by_cov,bz_cov,1)
            jy(i,j,k) = curl(i,j,k,bx_cov,by_cov,bz_cov,2)
            jz(i,j,k) = curl(i,j,k,bx_cov,by_cov,bz_cov,3)
cc            jx(i,j,k) = curl2(i,j,k,bx_cov,by_cov,bz_cov,1)
cc            jy(i,j,k) = curl2(i,j,k,bx_cov,by_cov,bz_cov,2)
cc            jz(i,j,k) = curl2(i,j,k,bx_cov,by_cov,bz_cov,3)
          enddo
        enddo
      enddo

c     Enforce div(J)=0 on contravariant current components

      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          bctype=DIR

          ieq = IBX
          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
            call FillGhostNodes(-ieq,dim,loc,bctype,jx,zeros)
          endif

          ieq = IBY
          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
            call FillGhostNodes(-ieq,dim,loc,bctype,jy,zeros)
          endif

          ieq = IBZ
          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
            call FillGhostNodes(-ieq,dim,loc,bctype,jz,zeros)
          endif

        enddo
      enddo

c     Impose SP boundary conditions on contravariant components of J

      if (bcond(1) == SP) call vectorSingularBC(jx,jy,jz,.false.,2)

c     Enforce PER BC on contravariant current components

      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          bctype=PER

          ieq = IBX
          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
            call FillGhostNodes(ieq,dim,loc,bctype,jx,zeros)
          endif

          ieq = IBY
          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
            call FillGhostNodes(ieq,dim,loc,bctype,jy,zeros)
          endif

          ieq = IBZ
          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
            call FillGhostNodes(ieq,dim,loc,bctype,jz,zeros)
          endif

        enddo
      enddo

c     Find covariant current

      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,jx_cov(i,j,k),jy_cov(i,j,k),jz_cov(i,j,k)
     .             ,jx(i,j,k),jy(i,j,k),jz(i,j,k),.false.)
          enddo
        enddo
      enddo

c     End program

      end subroutine postProcessJ

c     vectorSingularBC
c     #################################################################
      subroutine vectorSingularBC(vec1,vec2,vec3,cov,order)
c     -----------------------------------------------------------------
c     Averages vector components around singular point and calculates
c     curvilinear components at singular point.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: order
      real(8)    :: vec1(0:nx+1,0:ny+1,0:nz+1)
     .             ,vec2(0:nx+1,0:ny+1,0:nz+1)
     .             ,vec3(0:nx+1,0:ny+1,0:nz+1)
      logical    :: cov

c     Local variables

      integer(4) :: i,j,k
      logical    :: cartesian

      integer(4) :: ic,ig,jg,kg,order1
      real(8)    :: x0,avg_q,avg_vol,vol
      real(8)    :: cx1,cy1,cz1
      real(8),allocatable,dimension(:) :: ax0,ay0,az0,cx,cy,cz

c     External

      real(8) :: quad_int
      external   quad_int

c     Begin program

      allocate(ax0(nz),ay0(nz),az0(nz))

c     Find average cartesian coordinates

      if (order == 3) then
        order1 = 3
      else
        order1 = order
      endif

      do k=1,nz
        ax0(k) = 0d0
        ay0(k) = 0d0
        az0(k) = 0d0
        avg_vol = 0d0
        do j=1,ny
          i = 1
          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
          vol = volume(i,j,k,igx,igy,igz)
          x0  = xx(ig-1)

          allocate(cx(order1+1),cy(order1+1),cz(order1+1))
          do i=1,order1+1
            call transformVectorToCartesian(i,j,k,igx,igy,igz
     .           ,vec1(i,j,k),vec2(i,j,k),vec3(i,j,k),cov
     .           ,cx(i),cy(i),cz(i))

cc          write (*,*) cov
cc          write (*,*) vec1(i,j,k),vec2(i,j,k),vec3(i,j,k)
cc          write (*,*) cx(i),cy(i),cz(i)
cc          write (*,*) 'Inverse'
cc          call transformVectorToCurvilinear(i,j,k,igx,igy,igz
cc     .               ,cx(i),cy(i),cz(i),cov
cc     .               ,vec1(i,j,k),vec2(i,j,k),vec3(i,j,k))
cc          write (*,*) vec1(i,j,k),vec2(i,j,k),vec3(i,j,k)
cc          write (*,*)

          enddo

          call IntDriver1d(order1+1,xx(ig),cx,1,x0,cx1,order)
          call IntDriver1d(order1+1,xx(ig),cy,1,x0,cy1,order)
          call IntDriver1d(order1+1,xx(ig),cz,1,x0,cz1,order)
          deallocate(cx,cy,cz)

          ax0(k)  = ax0(k)  + vol*cx1
          ay0(k)  = ay0(k)  + vol*cy1
          az0(k)  = az0(k)  + vol*cz1
          avg_vol = avg_vol + vol
        enddo
        ax0(k) = ax0(k)/avg_vol
        ay0(k) = ay0(k)/avg_vol
        az0(k) = az0(k)/avg_vol
cc        write (*,*) ax0(k),ay0(k),az0(k)
      enddo
      
c     Transform to curvilinear components at SP

      i = 0
      do k=1,nz
        do j=1,ny
          call transformVectorToCurvilinear(i,j,k,igx,igy,igz
     .               ,ax0(k),ay0(k),az0(k),cov
     .               ,vec1(i,j,k),vec2(i,j,k),vec3(i,j,k))
        enddo
      enddo

      deallocate(ax0,ay0,az0)

ccc     ALTERNATIVE WAY: EXTRAPOLATE CURVILINEAR COORDINATES DIRECTLY
cc
cc      if (order == 3) then
cc        order1 = 3
cc      else
cc        order1 = order
cc      endif
cc
cc      do k=1,nz
cc        do j=1,ny
cc          i = 1
cc          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc          x0  = xx(ig-1)
cc
cc          allocate(cx(order1+1),cy(order1+1),cz(order1+1))
cc          do i=1,order1+1
cc            cx(i) = vec1(i,j,k)
cc            cy(i) = vec2(i,j,k)
cc            cz(i) = vec3(i,j,k)
cc          enddo
cc          call IntDriver1d(order1+1,xx(ig),cx,1,x0,cx1,order)
cc          call IntDriver1d(order1+1,xx(ig),cy,1,x0,cy1,order)
cc          call IntDriver1d(order1+1,xx(ig),cz,1,x0,cz1,order)
cc          deallocate(cx,cy,cz)
cc
cc          vec1(0,j,k) = cx1
cc          vec2(0,j,k) = cy1
cc          vec3(0,j,k) = cz1
cc        enddo
cc      enddo

c     End program

      end subroutine vectorSingularBC

      end subroutine imposeBoundaryConditions

c FillGhostNodes
c####################################################################
      subroutine FillGhostNodes(ieq,dim,loc,bctype,array,array0)

c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array.
c
c     On input:
c       * ieq    -> equation identifier
c       * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c       * loc    -> location in dimension (0 -> right, 1 -> left)
c       * bctype -> type of BC (dirichlet, neumann, periodic, etc.)
c       * array  -> real array with ghost-nodes
c       * array0 -> auxiliary real array
c--------------------------------------------------------------------

      use grid

      use nlfunction_setup

      use icond

      implicit none       !For safe fortran

c Call variables

      integer(4) :: ieq,dim,loc,bctype
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1)

c Local variables

      integer(4) :: neq,ibc
      integer(4) :: i,j,k,ig,jg,kg

      real(8),allocatable,dimension(:,:) :: rhs

c Begin program

      ibc = (1+loc)+2*(dim-1)

      select case (ibc)

      case (1)

c     X0

      select case(bctype)
      case(PER)
        array(0,:,:) = array(nx,:,:)
      case(EQU)
        array(0,:,:) = array0(0,:,:)
      case(SP)
        call singularBC(ieq,dim,loc,2)
      case(DIR)
        call dirichletBC(ieq,dim,loc,1)
        array(0,:,:) = rhs(:,:)
        deallocate(rhs)
      case(NEU)
        call neumannBC(ieq,dim,loc)
        array(0,:,:) = array(1,:,:) - rhs(:,:)
        deallocate(rhs)
      case default
        write (*,*) 'BC not implemented'
        stop
      end select

      case(2)

c     X1

      select case(bctype)
      case(PER)
        array(nx+1,:,:) = array(1,:,:)
      case(EQU)
        array(nx+1,:,:) = array0(nx+1,:,:)
      case(DIR)
        call dirichletBC(ieq,dim,loc,1)
        array(nx+1,:,:) = rhs(:,:)
        deallocate(rhs)
      case(NEU)
        call neumannBC(ieq,dim,loc)
        array(nx+1,:,:) = array(nx,:,:) + rhs(:,:)
        deallocate(rhs)
      case default
        write (*,*) 'BC not implemented'
        stop
      end select

      case (3)

c     Y0

      select case(bctype)
      case(PER)
        array(:,0,:) = array(:,ny,:)
      case(EQU)
        array(:,0,:) = array0(:,0,:)
      case(DIR)
        call dirichletBC(ieq,dim,loc,1)
        array(:,0,:) = rhs(:,:)
        deallocate(rhs)
      case(NEU)
        call neumannBC(ieq,dim,loc)
        array(:,0,:) = array(:,1,:) - rhs(:,:)
        deallocate(rhs)
      case default
        write (*,*) 'BC not implemented'
        stop
      end select

      case (4)

c     Y1

      select case(bctype)
      case(PER)
        array(:,ny+1,:) = array(:,1,:)
      case(EQU)
        array(:,ny+1,:) = array0(:,ny+1,:)
      case(DIR)
        call dirichletBC(ieq,dim,loc,1)
        array(:,ny+1,:) = rhs(:,:)
        deallocate(rhs)
      case(NEU)
        call neumannBC(ieq,dim,loc)
        array(:,ny+1,:) = array(:,ny,:) + rhs(:,:)
        deallocate(rhs)
      case default
        write (*,*) 'BC not implemented'
        stop
      end select

      case (5)

c     Z0

      select case(bctype)
      case(PER)
        array(:,:,0) = array(:,:,nz)
      case(EQU)
        array(:,:,0) = array0(:,:,0)
      case(DIR)
        call dirichletBC(ieq,dim,loc,1)
        array(:,:,0) = rhs(:,:)
        deallocate(rhs)
      case(NEU)
        call neumannBC(ieq,dim,loc)
        array(:,:,0) = array(:,:,1) - rhs(:,:)
        deallocate(rhs)
      case default
        write (*,*) 'BC not implemented'
        stop
      end select

      case (6)

c     Z1

      select case(bctype)
      case(PER)
        array(:,:,nz+1) = array(:,:,1)
      case(EQU)
        array(:,:,nz+1) = array0(:,:,nz+1)
      case(DIR)
        call dirichletBC(ieq,dim,loc,1)
        array(:,:,nz+1) = rhs(:,:)
        deallocate(rhs)
      case(NEU)
        call neumannBC(ieq,dim,loc)
        array(:,:,nz+1) = array(:,:,nz) + rhs(:,:)
        deallocate(rhs)
      case default
        write (*,*) 'BC not implemented'
        stop
      end select

      end select

c End

      contains

c     singularBC
c     #################################################################
      subroutine singularBC(ieq,dim,loc,order)
c     -----------------------------------------------------------------
c     Imposes singular point BC. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     -----------------------------------------------------------------

      use grid

      implicit none

c     Call variables

      integer(4) :: ieq,dim,loc,order

c     Local variables

      integer(4) :: i,j,k,order1
      real(8)    :: avg_q,avg_vol,rho0,vol

c     Begin program

      select case (ieq)
      case (IRHO,ITMP)  !Find average at singular point and store it in ghost node

cc        do k=1,nz
cc          avg_q = 0d0
cc          avg_vol = 0d0
cc          do j=1,ny
cc            avg_q   = avg_q   + volume(1,j,k,igx,igy,igz)*array(1,j,k)
cc            avg_vol = avg_vol + volume(1,j,k,igx,igy,igz)
cc          enddo
cc          array(0,:,k) = avg_q/avg_vol
cc        enddo

      if (order == 3) then
        order1 = 3
      else
        order1 = order
      endif

      do k=1,nz
        avg_q   = 0d0
        avg_vol = 0d0
        do j=1,ny
          i = 1
          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
          vol = volume(i,j,k,igx,igy,igz)
          x0  = xx(ig-1)

          call IntDriver1d(order1+1,xx(ig),array(i:i+order1+1,j,k)
     .                    ,1,x0,rho0,order)

          avg_q   = avg_q   + vol*rho0
          avg_vol = avg_vol + vol
        enddo
        array(0,:,k) = avg_q/avg_vol
      enddo

      end select

c     End program

      end subroutine singularBC

c     neumannBC
c     #################################################################
      subroutine neumannBC(ieq,dim,loc)
c     -----------------------------------------------------------------
c     Imposes neumann BC for a scalar. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Neumann BC.
c     -----------------------------------------------------------------

      use grid

      implicit none

c     Call variables

      integer(4) :: ieq,dim,loc

c     Local variables

      integer(4) :: i,j,k,ip,im,jp,jm,kp,km,icomp
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      real(8)    :: x1,x2,x3,dh(3),nabla_v(3,3)
      logical    :: cartesian

c     Begin program

      if (dim == 1) then
        imin=1  +    loc *(nx-1)
        imax=nx + (1-loc)*(1-nx)
        jmin=1
        jmax=ny
        kmin=1
        kmax=nz
        allocate(rhs(0:ny+1,0:nz+1))
      elseif (dim == 2) then
        imin=1 
        imax=nx
        jmin=1  +    loc *(ny-1)
        jmax=ny + (1-loc)*(1-ny)
        kmin=1
        kmax=nz
        allocate(rhs(0:nx+1,0:nz+1))
      elseif (dim == 3) then
        imin=1 
        imax=nx
        jmin=1
        jmax=ny
        kmin=1  +    loc *(nz-1)
        kmax=nz + (1-loc)*(1-nz)
        allocate(rhs(0:nx+1,0:ny+1))
      endif

      rhs = 0d0

      select case (ieq)
      case (IRHO,ITMP)

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
     .                           ,cartesian)

              gsuper = g_super(x1,x2,x3,cartesian)

              ip = min(i+1,nx)
              im = max(i-1,1)
              jp = min(j+1,ny)
              jm = max(j-1,1)
              kp = min(k+1,nz)
              km = max(k-1,1)

              dh(1) = 2.*dxh(ig)
              if (i == nx) dh(1) = dx(ig-1)
              if (i == 1 ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == ny) dh(2) = dy(jg-1)
              if (j == 1 ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nz) dh(3) = dz(kg-1)
              if (k == 1 ) dh(3) = dz(kg)

              if (dim == 1) then
                rhs(j,k) = -dh(dim)
     .             *(gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3))
     .              /gsuper(dim,dim)
              elseif (dim == 2) then
                rhs(i,k) = -dh(dim)*
     .              (gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
     .              +gsuper(dim,3)*(array(i,j,kp)-array(i,j,km))/dh(3))
     .              /gsuper(dim,dim)
              elseif (dim == 3) then
                rhs(i,j) = -dh(dim)
     .             *(gsuper(dim,1)*(array(ip,j,k)-array(im,j,k))/dh(1)
     .              +gsuper(dim,2)*(array(i,jp,k)-array(i,jm,k))/dh(2))
     .              /gsuper(dim,dim)
              endif

            enddo
          enddo
        enddo

      case (IVX,IVY,IVZ) !Velocity components

        if (ieq == IVX) icomp = 1
        if (ieq == IVY) icomp = 2
        if (ieq == IVZ) icomp = 3

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
     .                           ,cartesian)

              gsuper = g_super(x1,x2,x3,cartesian)

cc              if (loc == 0) then
cc                if (dim==1) then
cc                  nabla_v = fnabla_v_bc(i-1,j,k,dim)
cc                elseif (dim==2) then
cc                  nabla_v = fnabla_v_bc(i,j-1,k,dim)
cc                else
cc                  nabla_v = fnabla_v_bc(i,j,k-1,dim)
cc                endif
cc              else
                nabla_v = fnabla_v_bc(i,j,k,x1,x2,x3,cartesian,dim)
cc              endif

              dh(1) = 2.*dxh(ig)
              if (i == nx) dh(1) = dx(ig-1)
              if (i == 1 ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == ny) dh(2) = dy(jg-1)
              if (j == 1 ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nz) dh(3) = dz(kg-1)
              if (k == 1 ) dh(3) = dz(kg)

              if (dim == 1) then
                rhs(j,k) = -dh(dim)
     .               *(gsuper(dim,1)*nabla_v(1,icomp)
     .                +gsuper(dim,2)*nabla_v(2,icomp)
     .                +gsuper(dim,3)*nabla_v(3,icomp))
     .                /gsuper(dim,dim)
              elseif (dim == 2) then
                rhs(i,k) = -dh(dim)
     .               *(gsuper(dim,1)*nabla_v(1,icomp)
     .                +gsuper(dim,2)*nabla_v(2,icomp)
     .                +gsuper(dim,3)*nabla_v(3,icomp))
     .                /gsuper(dim,dim)
              elseif (dim == 3) then
                rhs(i,j) = -dh(dim)
     .               *(gsuper(dim,1)*nabla_v(1,icomp)
     .                +gsuper(dim,2)*nabla_v(2,icomp)
     .                +gsuper(dim,3)*nabla_v(3,icomp))
     .                /gsuper(dim,dim)
              endif

            enddo
          enddo
        enddo

      case (IBX,IBY,IBZ) 

        if (ieq == IBX) icomp = 1
        if (ieq == IBY) icomp = 2
        if (ieq == IBZ) icomp = 3

        if (icomp == dim) then
          write (*,*) 'Error in B in neumannBC: icomp = dim=',dim
          stop
        endif

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
     .                           ,cartesian)

              gsuper = g_super(x1,x2,x3,cartesian)

              ip = min(i+1,nx)
              im = max(i-1,1)
              jp = min(j+1,ny)
              jm = max(j-1,1)
              kp = min(k+1,nz)
              km = max(k-1,1)

              dh(1) = 2.*dxh(ig)
              if (i == nx) dh(1) = dx(ig-1)
              if (i == 1 ) dh(1) = dx(ig)

              dh(2) = 2.*dyh(jg)
              if (j == ny) dh(2) = dy(jg-1)
              if (j == 1 ) dh(2) = dy(jg)

              dh(3) = 2.*dzh(kg)
              if (k == nz) dh(3) = dz(kg-1)
              if (k == 1 ) dh(3) = dz(kg)

              if (dim == 1) then

                jx(i,j,k) = (bz_cov(i,jp,k)-bz_cov(i,jm,k))/dh(2)
     .                     -(by_cov(i,j,kp)-by_cov(i,j,km))/dh(3)

                if (icomp == 2) then
                  rhs(j,k) = dh(dim)
     .               *( (bx_cov(i,jp,k)-bx_cov(i,jm,k))/dh(2)
     .                 +gsuper(dim,3)*jx(i,j,k)/gsuper(dim,dim) )
                  if (loc == 0) then
                    by_cov(i-1,j,k) = by_cov(i,j,k) - rhs(j,k)
                  else
                    by_cov(i+1,j,k) = by_cov(i,j,k) + rhs(j,k)
                  endif
                elseif (icomp == 3) then
                  rhs(j,k) = dh(dim)
     .               *( (bx_cov(i,j,kp)-bx_cov(i,j,km))/dh(3)
     .                 -gsuper(dim,2)*jx(i,j,k)/gsuper(dim,dim) )
                  if (loc == 0) then
                    bz_cov(i-1,j,k) = bz_cov(i,j,k) - rhs(j,k)
                  else
                    bz_cov(i+1,j,k) = bz_cov(i,j,k) + rhs(j,k)
                  endif
                endif

              elseif (dim == 2) then

                jy(i,j,k) = (bx_cov(i,j,kp)-bx_cov(i,j,km))/dh(3)
     .                     -(bz_cov(ip,j,k)-bz_cov(im,j,k))/dh(1)

                if (icomp == 3) then
                  rhs(i,k) = dh(dim)
     .               *( (by_cov(i,j,kp)-by_cov(i,j,km))/dh(3)
     .                 +gsuper(dim,1)*jy(i,j,k)/gsuper(dim,dim) )
                  if (loc == 0) then
                    bz_cov(i,j-1,k) = bz_cov(i,j,k) - rhs(i,k)
                  else
                    bz_cov(i,j+1,k) = bz_cov(i,j,k) + rhs(i,k)
                  endif
                elseif (icomp == 1) then
                  rhs(i,k) = dh(dim)
     .               *( (by_cov(ip,j,k)-by_cov(im,j,k))/dh(1)
     .                 -gsuper(dim,3)*jy(i,j,k)/gsuper(dim,dim) )
                  if (loc == 0) then
                    bx_cov(i,j-1,k) = bx_cov(i,j,k) - rhs(i,k)
                  else
                    bx_cov(i,j+1,k) = bx_cov(i,j,k) + rhs(i,k)
                  endif
                endif

              elseif (dim == 3) then

                jz(i,j,k) = (by_cov(ip,j,k)-by_cov(im,j,k))/dh(1)
     .                     -(bx_cov(i,jp,k)-bx_cov(i,jm,k))/dh(2)

                if (icomp == 1) then
                  rhs(i,j) = dh(dim)
     .               *( (bz_cov(ip,j,k)-bz_cov(im,j,k))/dh(1)
     .                 +gsuper(dim,2)*jz(i,j,k)/gsuper(dim,dim) )
                  if (loc == 0) then
                    bx_cov(i,j,k-1) = bx_cov(i,j,k) - rhs(i,j)
                  else
                    bx_cov(i,j,k+1) = bx_cov(i,j,k) + rhs(i,j)
                  endif
                elseif (icomp == 2) then
                  rhs(i,j) = dh(dim)
     .               *( (bz_cov(i,jp,k)-bz_cov(i,jm,k))/dh(2)
     .                 -gsuper(dim,1)*jz(i,j,k)/gsuper(dim,dim) )
                  if (loc == 0) then
                    by_cov(i,j,k-1) = by_cov(i,j,k) - rhs(i,j)
                  else
                    by_cov(i,j,k+1) = by_cov(i,j,k) + rhs(i,j)
                  endif
                endif

              endif

            enddo
          enddo
        enddo

      case default

        write (*,*) 'Error in neumannBC'
        stop

      end select

c     End program

      end subroutine neumannBC

c     dirichletBC
c     #################################################################
      subroutine dirichletBC(ieq,dim,loc,order)
c     -----------------------------------------------------------------
c     Imposes dirichlet BC. On input:
c        * ieq -> equation number (i.e., vector component)
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c        * order -> order of extrapolation (when used)
c     This routine fills up the bi-dimensional array rhs, which 
c     contains the right hand side of the Dirichlet BC.
c     -----------------------------------------------------------------

      use grid

      implicit none

c     Call variables

      integer(4) :: ieq,dim,loc,order

c     Local variables

      integer(4) :: i,j,k,ip,im,jp,jm,kp,km,icomp,ibc
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      real(8)    :: x1,x2,x3,dh(3),nabla_v(3,3),jac
      logical    :: cartesian

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      if (dim == 1) then
        imin=1  +    loc *(nx-1)
        imax=nx + (1-loc)*(1-nx)
        jmin=1
        jmax=ny
        kmin=1
        kmax=nz
        allocate(rhs(0:ny+1,0:nz+1))
      elseif (dim == 2) then
        imin=1 
        imax=nx
        jmin=1  +    loc *(ny-1)
        jmax=ny + (1-loc)*(1-ny)
        kmin=1
        kmax=nz
        allocate(rhs(0:nx+1,0:nz+1))
      elseif (dim == 3) then
        imin=1 
        imax=nx
        jmin=1
        jmax=ny
        kmin=1  +    loc *(nz-1)
        kmax=nz + (1-loc)*(1-nz)
        allocate(rhs(0:nx+1,0:ny+1))
      endif

      rhs = 0d0

      select case (ieq)
      case (IRHO,ITMP,IVX,IVY,IVZ)

        call interpolate(imin,imax,jmin,jmax,kmin,kmax,order)

      case (IBX,IBY,IBZ) !Imposes divergence-free constraint on B-field

        if (ieq == IBX) icomp = 1
        if (ieq == IBY) icomp = 2
        if (ieq == IBZ) icomp = 3

        if (icomp /= dim) then

          call interpolate(imin,imax,jmin,jmax,kmin,kmax,order)

        else

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              dh(1) = 2.*dxh(ig)
              dh(2) = 2.*dyh(jg)
              dh(3) = 2.*dzh(kg)

              select case (ibc)
              case (1)
                array(i-1,j,k) = array(i+1,j,k)
                rhs(j,k) = array(i+1,j,k) + dh(1)*divB(i,j,k)
              case (2)
                array(i+1,j,k) = array(i-1,j,k)
                rhs(j,k) = array(i-1,j,k) - dh(1)*divB(i,j,k)
              case (3)
                array(i,j-1,k) = array(i,j+1,k)
                rhs(i,k) = array(i,j+1,k) + dh(2)*divB(i,j,k)
              case (4)
                array(i,j+1,k) = array(i,j-1,k)
                rhs(i,k) = array(i,j-1,k) - dh(2)*divB(i,j,k)
              case (5)
                array(i,j,k-1) = array(i,j,k+1)
                rhs(i,j) = array(i,j,k+1) + dh(3)*divB(i,j,k)
              case (6)
                array(i,j,k+1) = array(i,j,k-1)
                rhs(i,j) = array(i,j,k-1) - dh(3)*divB(i,j,k)
              end select

              enddo
            enddo
          enddo

        endif

      case (-IBX,-IBY,-IBZ) !Imposes divergence-free constraint on current

        if (ieq == -IBX) icomp = 1
        if (ieq == -IBY) icomp = 2
        if (ieq == -IBZ) icomp = 3

        if (icomp /= dim) then

          call interpolate(imin,imax,jmin,jmax,kmin,kmax,order)

        else

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              dh(1) = 2.*dxh(ig)
              dh(2) = 2.*dyh(jg)
              dh(3) = 2.*dzh(kg)

              select case (ibc)
              case (1)
                array(i-1,j,k) = array(i+1,j,k)
                rhs(j,k) = array(i+1,j,k) + dh(1)*divJ(i,j,k)
              case (2)
                array(i+1,j,k) = array(i-1,j,k)
                rhs(j,k) = array(i-1,j,k) - dh(1)*divJ(i,j,k)
              case (3)
                array(i,j-1,k) = array(i,j+1,k)
                rhs(i,k) = array(i,j+1,k) + dh(2)*divJ(i,j,k)
              case (4)
                array(i,j+1,k) = array(i,j-1,k)
                rhs(i,k) = array(i,j-1,k) - dh(2)*divJ(i,j,k)
              case (5)
                array(i,j,k-1) = array(i,j,k+1)
                rhs(i,j) = array(i,j,k+1) + dh(3)*divJ(i,j,k)
              case (6)
                array(i,j,k+1) = array(i,j,k-1)
                rhs(i,j) = array(i,j,k-1) - dh(3)*divJ(i,j,k)
              end select

              enddo
            enddo
          enddo

        endif

      case default

        write (*,*) 'Error in dirichletBC'
        stop

      end select

c     End program

      end subroutine dirichletBC

c     interpolate
c     #######################################################################
      subroutine interpolate(imin,imax,jmin,jmax,kmin,kmax,order)

        implicit none

c     Call variables

        integer(4) :: imin,imax,jmin,jmax,kmin,kmax,order

c     Local variables

        integer(4) :: i,j,k,ig,jg,kg
        real(8),allocatable,dimension(:) :: xint,vint
        real(8)                          :: xo(1),vo(1)
        logical    :: quadratic

c     Externals

        real(8)    :: quad_int
        external   :: quad_int

c     Begin program

        quadratic = .false.

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              select case (ibc)
              case (1)

cc                allocate(xint(0:nx),vint(0:nx))
cc
cc              !Define input abcissae
cc                xint(0)    = xx(ig-1)+dxh(ig-1)
cc                xint(1:nx) = xx(ig:ig+nx-1)
cc
cc              !Define extrapolation point(s)
cc                xo(1)      = xx(ig-1)
cc
cc              !Define vector of values
cc                vint(0)    = array0(i-1,j,k)
cc                vint(1:nx) = array(1:nx,j,k)
cc                call IntDriver1d(nx+1,xint,vint,1,xo,vo,order)
cc                rhs(j,k) = vo(1)
cc
cc                deallocate(xint,vint)

                rhs(j,k) = 
     .             quad_int(xx(ig-1)+dxh(ig-1),xx(ig),xx(ig+1),xx(ig+2)
     .                     ,array0(i-1,j,k),array(i,j,k)
     .                     ,array (i+1,j,k),array(i+2,j,k)
     .                     ,xx(ig-1),order )
              case (2)

cc                allocate(xint(0:nx),vint(0:nx))
cc
cc              !Define input abcissae
cc                xint(nx)     = xx(ig+1)-dxh(ig+1)
cc                xint(0:nx-1) = xx(ig-nx+1:ig)
cc
cc              !Define extrapolation point(s)
cc                xo(1)        = xx(ig+1)
cc
cc              !Define vector of values
cc                vint(nx)     = array0(i+1,j,k)
cc                vint(0:nx-1) = array(1:nx,j,k)
cc                call IntDriver1d(nx+1,xint,vint,1,xo,vo,order)
cc                rhs(j,k) = vo(1)
cc
cc                deallocate(xint,vint)

                rhs(j,k) =
     .              quad_int(xx(ig+1)-dxh(ig+1),xx(ig),xx(ig-1),xx(ig-2)
     .                      ,array0(i+1,j,k),array(i,j,k)
     .                      ,array (i-1,j,k),array(i-2,j,k)
     .                      ,xx(ig+1),order )
              case (3)

cc                allocate(xint(0:ny),vint(0:ny))
cc
cc              !Define input abcissae
cc                xint(0)    = yy(jg-1)+dyh(jg-1)
cc                xint(1:ny) = yy(jg:jg+ny-1)
cc
cc              !Define extrapolation point(s)
cc                xo(1)      = yy(jg-1)
cc
cc              !Define vector of values
cc                vint(0)    = array0(i,j-1,k)
cc                vint(1:ny) = array(i,1:ny,k)
cc                call IntDriver1d(ny+1,xint,vint,1,xo,vo,order)
cc                rhs(i,k) = vo(1)
cc
cc                deallocate(xint,vint)

                rhs(i,k) =
     .              quad_int(yy(jg-1)+dyh(jg-1),yy(jg),yy(jg+1),yy(jg+2)
     .                      ,array0(i,j-1,k),array(i,j,k)
     .                      ,array (i,j+1,k),array(i,j+2,k)
     .                      ,yy(jg-1),order )
              case (4)

cc                allocate(xint(0:ny),vint(0:ny))
cc
cc              !Define input abcissae
cc                xint(ny)     = yy(jg+1)-dyh(jg+1)
cc                xint(0:ny-1) = yy(jg-ny+1:jg)
cc
cc              !Define extrapolation point(s)
cc                xo(1)        = yy(jg+1)
cc
cc              !Define vector of values
cc                vint(ny)     = array0(i,j+1,k)
cc                vint(0:ny-1) = array(i,1:ny,k)
cc                call IntDriver1d(ny+1,xint,vint,1,xo,vo,order)
cc                rhs(i,k) = vo(1)
cc
cc                deallocate(xint,vint)

                rhs(i,k) =
     .              quad_int(yy(jg+1)-dyh(jg+1),yy(jg),yy(jg-1),yy(jg-2)
     .                      ,array0(i,j+1,k),array(i,j,k)
     .                      ,array (i,j-1,k),array(i,j-2,k)
     .                      ,yy(jg+1),order )
              case (5)

cc                allocate(xint(0:nz),vint(0:nz))
cc
cc              !Define input abcissae
cc                xint(0)    = zz(kg-1)+dzh(kg-1)
cc                xint(1:nz) = zz(kg:kg+nz-1)
cc
cc              !Define extrapolation point(s)
cc                xo(1)      = zz(kg-1)
cc
cc              !Define vector of values
cc                vint(0)    = array0(i,j,k-1)
cc                vint(1:nz) = array(i,j,1:nz)
cc                call IntDriver1d(nz+1,xint,vint,1,xo,vo,order)
cc                rhs(i,j) = vo(1)
cc
cc                deallocate(xint,vint)

                rhs(i,j) =
     .              quad_int(zz(kg-1)+dzh(kg-1),zz(kg),zz(kg+1),zz(kg+2)
     .                      ,array0(i,j,k-1),array(i,j,k)
     .                      ,array (i,j,k+1),array(i,j,k+2)
     .                      ,zz(kg-1),order )
              case (6)

cc                allocate(xint(0:nz),vint(0:nz))
cc
cc              !Define input abcissae
cc                xint(nz)     = zz(kg+1)-dzh(kg+1)
cc                xint(0:nz-1) = zz(kg-nz+1:kg)
cc
cc              !Define extrapolation point(s)
cc                xo(1)        = zz(kg+1)
cc
cc              !Define vector of values
cc                vint(nz)     = array0(i,j,k+1)
cc                vint(0:nz-1) = array(i,j,1:nz)
cc                call IntDriver1d(nz+1,xint,vint,1,xo,vo,order)
cc                rhs(i,j) = vo(1)
cc
cc                deallocate(xint,vint)

                rhs(i,j) =
     .              quad_int(zz(kg+1)-dzh(kg+1),zz(kg),zz(kg-1),zz(kg-2)
     .                      ,array0(i,j,k+1),array(i,j,k)
     .                      ,array (i,j,k-1),array(i,j,k-2)
     .                      ,zz(kg+1),order )
              end select
            enddo
          enddo
        enddo

      end subroutine interpolate

      end subroutine FillGhostNodes

c quad_int
c #####################################################################
      real(8) function quad_int(x0,x1,x2,x3,y0,y1,y2,y3,x,order)
     .        result(y)
c ---------------------------------------------------------------------
c     Quadratic interpolation (extrapolation).
c ---------------------------------------------------------------------

      implicit none

c Call variables

      integer(4) :: order
      real(8)    :: x0,x1,x2,x3,y0,y1,y2,y3,x

c Local variables

c Begin program

      select case (order)
      case (3)
        y = y0*(x-x1)*(x-x2)*(x-x3)/(x0-x1)/(x0-x2)/(x0-x3)
     .     +y1*(x-x0)*(x-x2)*(x-x3)/(x1-x0)/(x1-x2)/(x1-x3)
     .     +y2*(x-x0)*(x-x1)*(x-x3)/(x2-x0)/(x2-x1)/(x2-x3)
     .     +y3*(x-x0)*(x-x1)*(x-x2)/(x3-x0)/(x3-x1)/(x3-x2)
      case (2)
        y = y0*(x-x1)*(x-x2)/(x0-x1)/(x0-x2)
     .     +y1*(x-x0)*(x-x2)/(x1-x0)/(x1-x2)
     .     +y2*(x-x0)*(x-x1)/(x2-x0)/(x2-x1)
      case (1)
        y = y0*(x-x1)/(x0-x1)
     .     +y1*(x-x0)/(x1-x0)
      case (0)
        y = y0
      end select

c End program

      end function quad_int

c imposeBConFluxes
c####################################################################
      subroutine imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                            ,flxkp,flxkm,bconds)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on fluxes
c--------------------------------------------------------------------

      use grid

      use nlfunction_setup

      implicit none

c Call variables

      integer(4) :: i,j,k,bconds(6)
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

c Local variables

c Begin program

      if (i == 1 .and. bconds(1) == SP) flxim = 0d0

cc      if (i == nx .and. bconds(2) == NEU) then
cc        flxip = 0d0
cc      elseif (i == 1 .and.(bconds(1) == NEU.or.bconds(1) == SP))then
cc        flxim = 0d0
cc      endif
cc
cc      if (j == ny .and. bconds(4) == NEU) then
cc        flxjp = 0d0
cc      elseif (j == 1 .and. bconds(3) == NEU) then
cc        flxjm = 0d0
cc      endif
cc
cc      if (k == nz .and. bconds(6) == NEU) then
cc        flxkp = 0d0
cc      elseif (k == 1 .and. bconds(5) == NEU) then
cc        flxkm = 0d0
cc      endif

c End

      end subroutine imposeBConfluxes
