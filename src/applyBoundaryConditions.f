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

      select case (coords)
      case ('car','scl')

        bcsq = bcs(:,IRHO)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IRHO) = bcsq

        bcsq = bcs(:,IVX)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IVX) = bcsq

        bcsq = bcs(:,IVY)
        where (bcsq == OTH) bcsq = DIR
        bcs(:,IVY) = bcsq

        bcsq = bcs(:,IVZ)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IVZ) = bcsq

        bcsq = bcs(:,IBX)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IBX) = bcsq

        bcsq = bcs(:,IBY)
        where (bcsq == OTH) bcsq = DIR
        bcs(:,IBY) = bcsq

        bcsq = bcs(:,IBZ)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IBZ) = bcsq

        bcsq = bcs(:,ITMP)
        where (bcsq == OTH) bcsq = NEU !To allow isothermal case
        bcs(:,ITMP) = bcsq

      case ('cyl')

        bcsq = bcs(:,IRHO)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IRHO) = bcsq

        bcsq = bcs(:,IVX)
        where (bcsq == OTH) bcsq = DIR
        bcs(:,IVX) = bcsq

        bcsq = bcs(:,IVY)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IVY) = bcsq

        bcsq = bcs(:,IVZ)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IVZ) = bcsq

        bcsq = bcs(:,IBX)
        where (bcsq == OTH) bcsq = DIR
        bcs(:,IBX) = bcsq

        bcsq = bcs(:,IBY)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IBY) = bcsq

        bcsq = bcs(:,IBZ)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IBZ) = bcsq

        bcsq = bcs(:,ITMP)
        where (bcsq == OTH) bcsq = NEU !To allow isothermal case
        bcs(:,ITMP) = bcsq

      case ('tor')

        bcsq = bcs(:,IRHO)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IRHO) = bcsq

        bcsq = bcs(:,IVX)
        where (bcsq == OTH) bcsq = DIR
        bcs(:,IVX) = bcsq

        bcsq = bcs(:,IVY)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IVY) = bcsq

        bcsq = bcs(:,IVZ)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IVZ) = bcsq

        bcsq = bcs(:,IBX)
        where (bcsq == OTH) bcsq = DIR
        bcs(:,IBX) = bcsq

        bcsq = bcs(:,IBY)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IBY) = bcsq

        bcsq = bcs(:,IBZ)
        where (bcsq == OTH) bcsq = NEU
        bcs(:,IBZ) = bcsq

        bcsq = bcs(:,ITMP)
        where (bcsq == OTH) bcsq = NEU !To allow isothermal case
        bcs(:,ITMP) = bcsq

      case default

        write (*,*) 'Could not define boundary conditions for grid type'
        write (*,*) 'Aborting...'
        stop

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

      call setupNonlinearFunction(varray)

c Fill ghost nodes

      neq = varray%nvar

      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)
          do bctype=1,4    !Enforces a particular order in the BCs: DIR,NEU,PER,SP
            do ieq = 1,neq 
              if (varray%array_var(ieq)%bconds(ibc) == bctype) then
                call FillGhostNodes(ieq,dim,loc,bctype
     .                             ,varray%array_var(ieq)%array
     .                             ,u_0   %array_var(ieq)%array)
              endif
            enddo
          enddo
        enddo
      enddo

c Postprocess magnetic field

      call postProcessB

c Take care of vector components at singular point

      if (varray%array_var(1)%bconds(1) == SP) then
        call singularBCvector(rvx,rvy,rvz)
        call singularBCvector(bx,by,bz)
      endif

c End

      contains

c     singularBCvector
c     #################################################################
      subroutine singularBCvector(ax,ay,az)
c     -----------------------------------------------------------------
c     Averages scalar quantity around singular point
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg
      real(8)    :: cov(3),cnv1(3),cnv2(3),cnv3(3)
      real(8)    :: x1,x2,x20,x3,jac
      real(8)    :: ix(0:ny+1),iy(0:ny+1)

c     Begin program


      i = 1
      ig = i + grid_params%istartx(igx)

      x1  = dxh(ig)/2
      x20 = yy(1 + grid_params%istartx(igy))

      do k = 1,nz
        do j = 2,ny
          jg = j + grid_params%istarty(igy)
          kg = k + grid_params%istartz(igz)

          x2 = yy(jg)
          x3 = zz(kg)

          cnv1 = contravariantVector(1,x1,x20,x3)
          cnv2 = contravariantVector(2,x1,x20,x3)
          cnv3 = contravariantVector(3,x1,x20,x3)

          jac = jacobian(x1,x2,x3)

          cov = covariantVector(1,x1,x2,x3)
          ix(j) = ax(i,1,k)*jac*dot_product(cnv1,cov)
     .           +ay(i,1,k)*jac*dot_product(cnv2,cov)
cc     .           +az(i,1,k)*jac*dot_product(cnv3,cov)
          
          cov = covariantVector(2,x1,x2,x3)
          iy(j) = ax(i,1,k)*jac*dot_product(cnv1,cov)
     .           +ay(i,1,k)*jac*dot_product(cnv2,cov)
cc     .           +az(i,1,k)*jac*dot_product(cnv3,cov)

        enddo
        ax(i,2:ny,k) = ix(2:ny)
        ay(i,2:ny,k) = iy(2:ny)
      enddo

c     End program

      end subroutine singularBCvector

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

c     Begin program

c     Find covariant magnetic field components at Neumann boundaries

      do dim = 1,3
        do loc = 0,1
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

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax

                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

                if (dim == 1) then

                  if (loc == 0) then
cc
cc                    by_cov(i-1,j,k) = by(i-1,j,k)
cc                    bz_cov(i-1,j,k) = bz(i-1,j,k)

                    gsuper = g_super(xx(ig-1),yy(jg),zz(kg))

                    bx_cov(i-1,j,k) = -(gsuper(1,2)*by_cov(i-1,j,k)
     .                                 +gsuper(1,3)*bz_cov(i-1,j,k)
     .                                 -bx(i-1,j,k) )
     .                                 /gsuper(1,1)
                  else
cc
cc                    by_cov(i+1,j,k) = by(i+1,j,k)
cc                    bz_cov(i+1,j,k) = bz(i+1,j,k)

                    gsuper = g_super(xx(ig+1),yy(jg),zz(kg))

                    bx_cov(i+1,j,k) = -(gsuper(1,2)*by_cov(i+1,j,k)
     .                                 +gsuper(1,3)*bz_cov(i+1,j,k)
     .                                 -bx(i+1,j,k) )
     .                                 /gsuper(1,1)
                  endif

                elseif (dim == 2) then

                  if (loc == 0) then
cc
cc                    bx_cov(i,j-1,k) = bx(i,j-1,k)
cc                    bz_cov(i,j-1,k) = bz(i,j-1,k)

                    gsuper = g_super(xx(ig),yy(jg-1),zz(kg))

                    by_cov(i,j-1,k) = -(gsuper(2,1)*bx_cov(i,j-1,k)
     .                                 +gsuper(2,3)*bz_cov(i,j-1,k)
     .                                 -by(i,j-1,k) )
     .                                 /gsuper(2,2)
                  else
cc
cc                    bx_cov(i,j+1,k) = bx(i,j+1,k)
cc                    bz_cov(i,j+1,k) = bz(i,j+1,k)

                    gsuper = g_super(xx(ig),yy(jg+1),zz(kg))

                    by_cov(i,j+1,k) = -(gsuper(2,1)*bx_cov(i,j+1,k)
     .                                 +gsuper(2,3)*bz_cov(i,j+1,k)
     .                                 -by(i,j+1,k) )
     .                                 /gsuper(2,2)
                  endif

                elseif (dim == 3) then

                  if (loc == 0) then
cc
cc                    bx_cov(i,j,k-1) = bx(i,j,k-1)
cc                    by_cov(i,j,k-1) = by(i,j,k-1)

                    gsuper = g_super(xx(ig),yy(jg),zz(kg-1))

                    bz_cov(i,j,k-1) = -(gsuper(3,1)*bx_cov(i,j,k-1)
     .                                 +gsuper(3,2)*by_cov(i,j,k-1)
     .                                 -bz(i,j,k-1) )
     .                                 /gsuper(3,3)
                  else
cc
cc                    bx_cov(i,j,k+1) = bx(i,j,k+1)
cc                    by_cov(i,j,k+1) = by(i,j,k+1)

                    gsuper = g_super(xx(ig),yy(jg),zz(kg+1))

                    bz_cov(i,j,k+1) = -(gsuper(3,1)*bx_cov(i,j,k+1)
     .                                 +gsuper(3,2)*by_cov(i,j,k+1)
     .                                 -bz(i,j,k+1) )
     .                                 /gsuper(3,3)
                  endif

                endif

              enddo
            enddo
          enddo

        enddo
      enddo

c     Enforce periodic BC on covariant magnetic field components

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

c     Find contravariant magnetic field

      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            gsuper = G_super(xx(ig),yy(jg),zz(kg))
            cov = (/ bx_cov(i,j,k),by_cov(i,j,k),bz_cov(i,j,k) /)
            cnv = matmul(gsuper,cov)
            bx(i,j,k) = cnv(1)
            by(i,j,k) = cnv(2)
            bz(i,j,k) = cnv(3)
          enddo
        enddo
      enddo

c     End program

      end subroutine postProcessB

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

      logical    :: quadratic

      real(8),allocatable,dimension(:,:) :: rhs

c Begin program

      quadratic = .false.

      ibc = (1+loc)+2*(dim-1)

      select case (ibc)

      case (1)

c     X0

      select case(bctype)
      case(PER)
        array(0,:,:) = array(nx,:,:)
      case(SP)
        call singularBC(array(1,:,:))
      case(DIR)
cc        array(0,:,:) = array0(0,:,:)
        i = 0
        do j=0,ny+1
          do k=0,nz+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            array(i,j,k) = quad_int(xx(ig)+dxh(ig),xx(ig+1),xx(ig+2)
     .                     ,array0(i,j,k),array(i+1,j,k),array(i+2,j,k)
     .                     ,xx(ig),quadratic )
          enddo
        enddo

      case(NEU)
cc        array(0,:,:) = array(1,:,:)        !Required for proper initialization
        call neumannBC(ieq,dim,loc)
c2nd        array(0,:,:) = array(2,:,:) - rhs(:,:)
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
      case(DIR)
cc        array(nx+1,:,:) = array0(nx+1,:,:)

        i = nx+1
        do j=0,ny+1
          do k=0,nz+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            array(i,j,k) = quad_int(xx(ig)-dxh(ig),xx(ig-1),xx(ig-2)
     .                     ,array0(i,j,k),array(i-1,j,k),array(i-2,j,k)
     .                     ,xx(ig),quadratic )
          enddo
        enddo

      case(NEU)
cc        array(nx+1,:,:) = array(nx,:,:)        !Required for proper initialization
        call neumannBC(ieq,dim,loc)
c2nd        array(nx+1,:,:) = array(nx-1,:,:) + rhs(:,:)
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
      case(DIR)
cc        array(:,0,:) = array0(:,0,:)
        j = 0
        do i=0,nx+1
          do k=0,nz+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            array(i,j,k) = quad_int(yy(jg)+dyh(jg),yy(jg+1),yy(jg+2)
     .                     ,array0(i,j,k),array(i,j+1,k),array(i,j+2,k)
     .                     ,yy(jg),quadratic )
          enddo
        enddo
      case(NEU)
cc        array(:,0,:) = array(:,1,:)        !Required for proper initialization
        call neumannBC(ieq,dim,loc)
c2nd        array(:,0,:) = array(:,2,:) - rhs(:,:)
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
      case(DIR)
cc        array(:,ny+1,:) = array0(:,ny+1,:)
        j = ny+1
        do i=0,nx+1
          do k=0,nz+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            array(i,j,k) = quad_int(yy(jg)-dyh(jg),yy(jg-1),yy(jg-2)
     .                     ,array0(i,j,k),array(i,j-1,k),array(i,j-2,k)
     .                     ,yy(jg),quadratic )
          enddo
        enddo
      case(NEU)
cc        array(:,ny+1,:) = array(:,ny,:)         !Required for proper initialization
        call neumannBC(ieq,dim,loc)
c2nd        array(:,ny+1,:) = array(:,ny-1,:) + rhs(:,:)
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
      case(DIR)
cc        array(:,:,0) = array0(:,:,0)
        k = 0
        do i=0,nx+1
          do j=0,ny+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            array(i,j,k) = quad_int(zz(kg)+dzh(kg),zz(kg+1),zz(kg+2)
     .                     ,array0(i,j,k),array(i,j,k+1),array(i,j,k+2)
     .                     ,zz(kg),quadratic )
          enddo
        enddo
      case(NEU)
cc        array(:,:,0) = array(:,:,1)         !Required for proper initialization
        call neumannBC(ieq,dim,loc)
c2nd        array(:,:,0) = array(:,:,2) - rhs(:,:)
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
      case(DIR)
cc        array(:,:,nz+1) = array0(:,:,nz+1)
        k = nz+1
        do i=0,nx+1
          do j=0,ny+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            array(i,j,k) = quad_int(zz(kg)-dzh(kg),zz(kg-1),zz(kg-2)
     .                     ,array0(i,j,k),array(i,j,k-1),array(i,j,k-2)
     .                     ,zz(kg),quadratic )
          enddo
        enddo
      case(NEU)
cc        array(:,:,nz+1) = array(:,:,nz)         !Required for proper initialization
        call neumannBC(ieq,dim,loc)
c2nd        array(:,:,nz+1) = array(:,:,nz-1) + rhs(:,:)
        array(:,:,nz+1) = array(:,:,nz) + rhs(:,:)
        deallocate(rhs)
      case default
        write (*,*) 'BC not implemented'
        stop
      end select

      end select

c End

      contains

c     quad_int
c     #################################################################
      real(8) function quad_int(x0,x1,x2,y0,y1,y2,x,quadratic) result(y)
c     -----------------------------------------------------------------
c     Quadratic interpolation (extrapolation).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: x0,x1,x2,y0,y1,y2,x
      logical    :: quadratic

c     Local variables

c     Begin program

      if (quadratic) then
        y = y0*(x-x1)*(x-x2)/(x0-x1)/(x0-x2)
     .     +y1*(x-x0)*(x-x2)/(x1-x0)/(x1-x2)
     .     +y2*(x-x0)*(x-x1)/(x2-x0)/(x2-x1)
      else
        y = y0*(x-x1)/(x0-x1)
     .     +y1*(x-x0)/(x1-x0)
      endif

c     End program

      end function quad_int

c     singularBC
c     #################################################################
      subroutine singularBC(arr)
c     -----------------------------------------------------------------
c     Averages scalar quantity around singular point
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: arr(0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: j,k
      real(8)    :: avg_q,avg_vol

c     Begin program

      do k=1,nz
        avg_q = 0d0
        avg_vol = 0d0
        do j=1,ny
          avg_q   = avg_q   + volume(1,j,k,igx,igy,igz)*arr(j,k)
          avg_vol = avg_vol + volume(1,j,k,igx,igy,igz)
        enddo
        arr(:,k) = avg_q/avg_vol
      enddo
          
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

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              x1 = xx(ig)
              x2 = yy(jg)
              x3 = zz(kg)

              gsuper = g_super(x1,x2,x3)

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

        if (icomp == dim) then
          write (*,*) 'Error in neumannBC: icomp = dim'
          stop
        endif

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              x1 = xx(ig)
              x2 = yy(jg)
              x3 = zz(kg)

              gsuper = g_super(x1,x2,x3)
cc
cc              if (loc == 0) then
cc                if (dim==1) then
cc                  nabla_v = fnabla_v_bc(i-1,j,k,dim,1d0)
cc                elseif (dim==2) then
cc                  nabla_v = fnabla_v_bc(i,j-1,k,dim,1d0)
cc                else
cc                  nabla_v = fnabla_v_bc(i,j,k-1,dim,1d0)
cc                endif
cc              else
                nabla_v = fnabla_v_bc(i,j,k,dim,1d0)
cc              endif

cc              nabla_v = fnabla_v_bc(i,j,k,0,0d0)

c2nd              nabla_v = fnabla_v(i,j,k,0,1d0)

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

                if (rhs(i,k) /= 0d0) then
cc                  write (*,*) i,j,k
cc                  write (*,*) dim,icomp
cc                  write (*,*) gsuper(dim,1),gsuper(dim,2),gsuper(dim,3)
cc     .                       ,gsuper(dim,dim)
cc                  write (*,*) nabla_v(1,icomp),nabla_v(2,icomp)
cc     .                       ,nabla_v(3,icomp)
cc                  write (*,*) rhs(i,k)
cc                  stop
                endif
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
          write (*,*) 'Error in neumannBC: icomp = dim'
          stop
        endif

        do i=imin,imax
          do j=jmin,jmax
            do k=kmin,kmax

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              x1 = xx(ig)
              x2 = yy(jg)
              x3 = zz(kg)

              gsuper = g_super(x1,x2,x3)
              gsub   = g_sub  (x1,x2,x3)

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

                jx(i,j,k) = (bz_cov(i,jp,k)-bz_cov(i,jp,k))/dh(2)
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
     .                 +gsuper(dim,2)*jy(i,j,k)/gsuper(dim,dim) )
                  if (loc == 0) then
                    bx_cov(i,j,k-1) = bx_cov(i,j,k) - rhs(i,j)
                  else
                    bx_cov(i,j,k+1) = bx_cov(i,j,k) + rhs(i,j)
                  endif
                elseif (icomp == 2) then
                  rhs(i,j) = dh(dim)
     .               *( (bz_cov(i,jp,k)-bz_cov(i,jm,k))/dh(2)
     .                 -gsuper(dim,1)*jy(i,j,k)/gsuper(dim,dim) )
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

      end subroutine FillGhostNodes

c imposeBConFluxes
c####################################################################
      subroutine imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                            ,flxkp,flxkm,bcond)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on fluxes
c--------------------------------------------------------------------

      use grid

      use nlfunction_setup

      implicit none

c Call variables

      integer(4) :: i,j,k,bcond(6)
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

c Local variables

c Begin program

      if (i == nx) then
        if (bcond(2) == NEU) then
          flxip = 0d0
        endif
      elseif (i == 1) then
        if (bcond(1) == NEU) then
          flxim = 0d0
        elseif (bcond(1) == SP) then
          flxim = 0d0
          flxjp = 0d0
          flxjm = 0d0
        endif
      endif

      if (j == ny) then
        if (bcond(4) == NEU) then
          flxjp = 0d0
        endif
      elseif (j == 1) then
        if (bcond(3) == NEU) then
          flxjm = 0d0
        endif
      endif

      if (k == nz) then
        if (bcond(6) == NEU) then
          flxkp = 0d0
        endif
      elseif (k == 1) then
        if (bcond(5) == NEU) then
          flxkm = 0d0
        endif
      endif

c End

      end subroutine imposeBConfluxes
