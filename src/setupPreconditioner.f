c setupPreconditioner
c####################################################################
      subroutine setupPreconditioner(ntot,x)

c--------------------------------------------------------------------
c     Finds velocity and magnetic field components in all grids
c     for matrix-free preconditioning.
c--------------------------------------------------------------------

      use grid

      use precond_variables

      use variables

      implicit none

c Call variables

      integer ::  ntot
      real*8      x(ntot)

c Local variables

cc      integer*4   i,j,ig,ieq,ii,iii,isig,ip,jp,jj,jjj,dd,neq

cc      integer,allocatable,dimension(:,:) :: bbcond

      type (var_array),target :: varray

c Externals

cc      external    bihar_mtvc,bihar_mtvc_cpl,psi_mtvc_cpl,bz_mtvc_cpl
cc     .           ,alf_mtvc

c Begin program

c Allocate variables

cc      allocate (bxx(ntotd2p),byy(ntotd2p))
cc      allocate (vxx(ntotd2p),vyy(ntotd2p))
cc      allocate (vxe(ntotd2p),vye(ntotd2p))
cc      allocate (pp(ndiagdp,ntotd2p,neqd))
cc      allocate (diag_mu(ntotd2p))
cc      allocate (idiagp(ndiagdp,10),bcs(4,neqd))
cc      allocate (ue(0:nxdp,0:nydp))
cc      allocate ( d_sc(1,  ntotd2p)
cc     .          ,d_bc(2,2*ntotd2p)
cc     .          ,d_pc(2,2*ntotd2p)
cc     .          ,d_bh(2,2*ntotd2p)
cc     .          ,d_alf(1, ntotd2p))
cc
ccc Unpack vector x, taking BCs from t=n solution
cc
cc      call mapVectorToStructure(x,varray)
cc
ccc Extract arrays and BC's
cc
cc      uu => varray%array_var(PHI)%array
cc      vv => varray%array_var(VOR)%array
cc      tt => varray%array_var(PSI)%array
cc      ww => varray%array_var(IVZ)%array
cc      nn => varray%array_var(IBZ)%array
cc
cc      do ieq=1,neqd
cc        bcs(:,ieq) = varray%array_var(PSI)%bconds
cc      enddo
cc
ccc Build preconditioning matrix pp
cc
ccc     Stencil mapping
cc
cc      do i = 1,ngrd
cc        idiagp(1,i) = -nxvp(i)
cc        idiagp(2,i) = -1
cc        idiagp(3,i) = 0
cc        idiagp(4,i) = 1
cc        idiagp(5,i) = nxvp(i)
cc      enddo
cc
ccc     Compute the diagonal blocks of the jacobian in all grids
cc
cc      call buildDiagBlocks(varray)
cc
ccc Calculate magnetic field components (for SI operators) 
ccc PROBLEM: Setting order=1 here destroys symmetry in whistler matvec.
cc
cc      call restrictVectorFieldFromSF(nxd,nyd,tt,bxx,byy,0,bcs(:,PSI))
cc
ccc Find ion velocity components in all grids
cc
cc      call restrictVectorFieldFromSF(nxd,nyd,uu,vxx,vyy,0,bcs(:,PHI))
cc
ccc Find electron velocity components in all grids
cc
cc      call restrictVectorFieldFromSF(nxd,nyd,ue,vxe,vye,0,bcs(:,PHI))
cc
ccc Store appropriate boundary conditions
cc
cc      bc_sc(:,1) = bcs(:,PSI)
cc
cc      bc_cpld(:,1) = bcs(:,PSI)
cc      bc_cpld(1,2) = 3          !Bottom: linear extrapolation
cc      bc_cpld(2,2) = 3          !Top
cc      bc_cpld(3,2) = 0          !Left: periodic
cc      bc_cpld(4,2) = 0          !Right
cc
ccc Find required diagonals in all grids
cc
cc      diag_mu(:) = pp(3,:,VOR)
cc
cc      if (di > 0d0) then
cc        neq = 1
cc        call find_mf_diag_neq(neq,neq*ntotdp,bihar_mtvc,ngrd,bc_sc,d_sc)
cc
cc
cc        if (precon == 'bc') then
cc          neq = 2
cc          call find_mf_diag_neq(neq,neq*ntotdp,bz_mtvc_cpl,ngrd,bc_cpld
cc     .                         ,d_bc)
cc          call find_mf_diag_neq(neq,neq*ntotdp,bihar_mtvc_cpl,ngrd
cc     .                         ,bc_cpld,d_bh)
cc
cc        elseif (precon == 'pc') then
cc          neq = 2
cc          call find_mf_diag_neq(neq,neq*ntotdp,psi_mtvc_cpl,ngrd
cc     .                         ,bc_cpld,d_pc)
cc          call find_mf_diag_neq(neq,neq*ntotdp,bihar_mtvc_cpl,ngrd
cc     .                         ,bc_cpld,d_bh)
cc        endif
cc      else
cc        neq = 1
cc        call find_mf_diag_neq(neq,neq*ntotdp,alf_mtvc,ngrd,bc_sc,d_alf)
cc      endif

c End program

      return
      end

ccc buildDiagBlocks
ccc###################################################################
cc      subroutine buildDiagBlocks(varray)
ccc-------------------------------------------------------------------
ccc    This subroutine computes the diagonal blocks of the Jacobian.
ccc
ccc     Boundary conditions are imposed via bcond(4):
ccc        * bcond(1)  -- > bottom
ccc        * bcond(2)  -- > top
ccc        * bcond(3)  -- > left
ccc        * bcond(4)  -- > right
ccc     If bcond = 0, BC's --> periodic
ccc     If bcond = 1, BC's --> dirichlet
ccc     If bcond = 2, BC's --> neumann (natural)
ccc-------------------------------------------------------------------
cc
cc      use grid
cc
cc      use parameters
cc
cc      use precond_variables
cc
cc      use variable_setup
cc
cc      use nlfunction_setup
cc
cc      implicit none
cc
ccc Call variables
cc
cc      type (var_array),target :: varray
cc
ccc Local variables
cc
cc      integer*4     i,j
cc
cc      real*8        diff2(0:nxd+1,0:nyd+1),sf(0:nxd+1,0:nyd+1)
cc
cc      integer, pointer,dimension(:) :: bcond
cc
ccc Begin program
cc
ccc Stream function block
cc
cc      diff2 = -1d0
cc      sf    = 0d0
cc
cc      call createMatrix(pp(:,:,PHI),1d30,1d0,sf,diff2,bcs(:,PHI),2d0)
cc
ccc Vorticity block
cc
cc      sf = uu
cc
cc      call createMatrix(pp(:,:,VOR),dt,alpha,sf,nuu,bcs(:,VOR),1d0)
cc
ccc Poloidal flux block
cc
cc      diff2 = eeta + de**2/dt/alpha
cc
cc      sf = ue
cc
cc      call createMatrix(pp(:,:,PSI),dt,alpha,sf,diff2,bcs(:,PSI),1d0)
cc
ccc vz block
cc
cc      sf = uu
cc
cc      call createMatrix(pp(:,:,IVZ),dt,alpha,sf,nuu,bcs(:,IVZ),1d0)
cc
ccc Bz block
cc
cc      call createMatrix(pp(:,:,IBZ),dt,alpha,sf,diff2,bcs(:,IBZ),1d0)
cc
ccc End program
cc
cc      return
cc      end
cc
ccc createMatrix
ccc####################################################################
cc      subroutine createMatrix(p1,dtt,coeff,phi,diff,bcond,rscale)
cc
ccc--------------------------------------------------------------------
ccc     Computes and stores discretization matrix for a convection
ccc     diffusion equation with stream function phi, diffusion
ccc     coefficient diff, and SI constant si, in all grids.
ccc
ccc     Boundary conditions are stored in bcond(4):
ccc        * bcond(1)  -- > bottom
ccc        * bcond(2)  -- > top
ccc        * bcond(3)  -- > left
ccc        * bcond(4)  -- > right
ccc     If bcond = 0, BC's --> periodic
ccc     If bcond = 1, BC's --> dirichlet
ccc     If bcond = 2, BC's --> neumann (natural)
ccc--------------------------------------------------------------------
cc
cc      use parameters
cc
cc      use grid
cc
cc      use timeStepping
cc
cc      implicit none    !For safe fortran
cc
ccc Call variables
cc
cc      real*8       dtt,coeff,rscale
cc      real*8       p1(ndiagdp,ntotd2p)
cc      real*8       phi(0:nxdp,0:nydp)
cc      real*8       diff(0:nxdp,0:nydp)
cc      integer*4    bcond(4)
cc
ccc Local variables
cc
cc      integer*4    i,j,ii,isig
cc      real*8       uvel,vvel,constx,consty
cc      real*8       nuvel,nvvel,theta_x,theta_y
cc
ccc Begin program
cc
cc      isig = istartp(ngrd) - 1
cc
cc      do j = 1,nyd
cc        do i = 1,nxd
cc
cc          ii  = i + nxd*(j-1) + isig
cc
ccc     Define uvel and vvel from the stream function
cc
cc          uvel = -(phi(i,j+1)-phi(i,j-1))/dy(ngrd)/2.0d0
cc          vvel =  (phi(i+1,j)-phi(i-1,j))/dx(ngrd)/2.0d0
cc
ccc     Center
cc
cc          constx = 2d0
cc          consty = 2d0
cc          if (bcond(1).eq.2.and.j.eq.1  ) constx = 1d0
cc          if (bcond(2).eq.2.and.j.eq.nyd) constx = 1d0
cc          if (bcond(3).eq.2.and.i.eq.1  ) consty = 1d0
cc          if (bcond(4).eq.2.and.i.eq.nxd) consty = 1d0
cc
cc          p1(3,ii) = dy(ngrd)*dx(ngrd)/coeff/dtt
cc     &             + diff(i,j)*constx*dx(ngrd)/dy(ngrd)
cc     &             + diff(i,j)*consty*dy(ngrd)/dx(ngrd)
cccc     &             + dx(ngrd)*abs(vvel) + dy(ngrd)*abs(uvel)
cc
cc          theta_x = 0d0
cc          if (uvel.ne.0d0)
cc     .       theta_x = min(dabs(p1(3,ii)*2.*dx(ngrd)/uvel),1d0)
cc
cc          theta_y = 0d0
cc          if (vvel.ne.0d0)
cc     .       theta_y = min(dabs(p1(3,ii)*2.*dy(ngrd)/vvel),1d0)
cc
cc          nuvel = (1.-theta_x)*uvel
cc          nvvel = (1.-theta_y)*vvel
cc
cc          p1(3,ii) = p1(3,ii)
cc     &             + dx(ngrd)*abs(nvvel) + dy(ngrd)*abs(nuvel)
cc
ccc     South
cc
cc          if( j .eq. 1.and.bcond(1).ne.0)then
cc            p1(1,ii) = 0.0d0
cc          else
cc            p1(1,ii) = -diff(i,j)*dx(ngrd)/dy(ngrd)
cc     &                 -dx(ngrd)*theta_y*vvel/2d0
cc     &                 -dx(ngrd)*max(0d0,nvvel)
cc          endif
cc
ccc     West
cc
cc          if( i .eq. 1.and.bcond(3).ne.0)then
cc            p1(2,ii) = 0.0d0
cc          else
cc            p1(2,ii) = -diff(i,j)*dy(ngrd)/dx(ngrd)
cc     &                 -dy(ngrd)*theta_x*uvel/2d0
cc     &                 -dy(ngrd)*max(0d0,nuvel)
cc          endif
cc
ccc     East
cc
cc          if( i .eq. nxd.and.bcond(4).ne.0)then
cc            p1(4,ii) = 0.0d0
cc          else
cc            p1(4,ii) = -diff(i,j)*dy(ngrd)/dx(ngrd)
cc     &                 +dy(ngrd)*theta_x*uvel/2d0
cc     &                 +dy(ngrd)*min(0d0,nuvel)
cc          endif
cc
ccc     North
cc
cc          if( j .eq. nyd.and.bcond(2).ne.0)then
cc            p1(5,ii) = 0.0d0
cc          else
cc            p1(5,ii) = -diff(i,j)*dx(ngrd)/dy(ngrd)
cc     &                 +dx(ngrd)*theta_y*vvel/2d0
cc     &                 +dx(ngrd)*min(0d0,nvvel)
cc          endif
cc
cc        enddo
cc      enddo
cc
ccc Normalization
cc
cc      p1 = coeff*p1
cc
ccc Piece-wise constant operator restriction to all grids.
cc
cc      call restrict5ptMatrix(p1,rscale,bcond)
cc
ccc End program
cc
cc      return
cc      end
cc
ccc restrict5ptMatrix
ccc##################################################################
cc      subroutine restrict5ptMatrix(p1,rscasf,bcond)
cc
ccc------------------------------------------------------------------
ccc     This is a routine to build coarse grid 5pt-stencil jacobians
ccc     from fine grid 5pt-stencil jacobians
ccc------------------------------------------------------------------
cc
cc      use parameters
cc
cc      use grid
cc
cc      implicit none    !For safe fortran
cc
ccc Call variables
cc
cc      real*8       p1(ndiagdp,ntotd2p),rscasf
cc      integer*4    bcond(4)
cc
ccc Local variables
cc
cc      integer*4    i,j,if,jf,ic,jc,ii,isig,isigf,isigc,ig
cc      integer*4    nxc,nxf,nyc,nyf
cc      integer*4    iuf4,iuf3,iuf2,iuf1,iuc
cc
ccc Begin program
cc
cc      do ig = ngrd-1,1,-1
cc
cc        nxc = nxvp(ig)
cc        nyc = nyvp(ig)
cc
cc        nxf = nxvp(ig+1)
cc        nyf = nyvp(ig+1)
cc
cc        isigc = istartp(ig  ) - 1
cc        isigf = istartp(ig+1) - 1
cc
cc        do ic = 1,nxc
cc          do jc = 1,nyc
cc
cc            if = 2*ic
cc            jf = 2*jc
cc
cc            iuc  = ic     + nxc*(jc-1) + isigc
cc            iuf4 = if     + nxf*(jf-1) + isigf
cc            iuf3 = if - 1 + nxf*(jf-1) + isigf
cc            iuf2 = if     + nxf*(jf-2) + isigf
cc            iuf1 = if - 1 + nxf*(jf-2) + isigf
cc
ccc        south
cc
cc            if( jc .eq. 1.and.bcond(1).ne.0)then
cc              p1(1,iuc) = 0.0d0
cc            else
cc              p1(1,iuc) = (p1(1,iuf2) + p1(1,iuf1))/rscasf
cc            endif
cc
ccc        west
cc
cc            if( ic .eq. 1.and.bcond(3).ne.0)then
cc              p1(2,iuc) = 0.0d0
cc            else
cc              p1(2,iuc) = (p1(2,iuf3) + p1(2,iuf1))/rscasf
cc            endif
cc
ccc        center
cc
cc            p1(3,iuc)=(p1(3,iuf1)+p1(3,iuf2)+p1(3,iuf3)+p1(3,iuf4)
cc     .               + p1(1,iuf4) + p1(2,iuf4)
cc     .               + p1(1,iuf3) + p1(4,iuf3)
cc     .               + p1(2,iuf2) + p1(5,iuf2)
cc     .               + p1(4,iuf1) + p1(5,iuf1))/rscasf
cc
ccc        east
cc
cc            if( ic .eq. nxc .and.bcond(4).ne.0)then
cc              p1(4,iuc) = 0.0d0
cc            else
cc              p1(4,iuc) = (p1(4,iuf2) + p1(4,iuf4))/rscasf
cc            endif
cc
ccc        north
cc
cc            if( jc .eq. nyc .and.bcond(2).ne.0)then
cc              p1(5,iuc) = 0.0d0
cc            else
cc              p1(5,iuc) = (p1(5,iuf3) + p1(5,iuf4))/rscasf
cc            endif
cc
cc          enddo
cc        enddo
cc      enddo
cc
ccc End 
cc
cc      return
cc      end
cc
ccc restrictVectorFieldFromSF
ccc####################################################################
cc      subroutine restrictVectorFieldFromSF(nx,ny,sf,vxx,vyy,order,bcond)
ccc-------------------------------------------------------------------
ccc     Finds vector components from stream function psi in all grids,
ccc     and maps them into MG vectors.
ccc-------------------------------------------------------------------
cc
cc      use grid
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer         nx,ny,order,bcond(4)
cc      real*8          sf(0:nx+1,0:ny+1)
cc     .               ,vx(0:nx+1,0:ny+1),vy(0:nx+1,0:ny+1)
cc      real*8          vxx(2*nx*ny),vyy(2*nx*ny)
cc
ccc Local variables
cc
cc      integer*4       nxf,nyf,nxc,nyc,igrid,igridf
cc
cc      real*8          dx1,dy1,coeff,coeff1,coeff2
cc
cc      double precision, allocatable,dimension(:,:)::sff,sfc,vxc,vyc
cc
ccc Begin program
cc
ccc Map vector components in finest grid onto MG vector
cc
cc      dx1 = dx(ngrd)
cc      dy1 = dy(ngrd)
cc
cc      nxf  = nxvp(ngrd)
cc      nyf  = nyvp(ngrd)
cc
cc      call vec_find(nxf,nyf,dx1,dy1,sf,vx,vy)
cc
cc      call mapArrayToMGVector(nxf,nyf,vx,vxx,ngrd)
cc      call mapArrayToMGVector(nxf,nyf,vy,vyy,ngrd)
cc
cc      allocate(sff(0:nxf+1,0:nyf+1))
cc
cc      sff = sf
cc
cc      igridf = ngrd
cc
ccc Restrict stream function and find vector components in coarser grids
cc
cc      do igrid = ngrd-1,2,-1
cc
ccc     Characterize coarse grid and define arrays
cc
cc        dx1 = dx(igrid)
cc        dy1 = dy(igrid)
cc
cc        nxc  = nxvp(igrid)
cc        nyc  = nyvp(igrid)
cc
cc        allocate(sfc(0:nxc+1,0:nyc+1))
cc        allocate(vxc (0:nxc+1,0:nyc+1))
cc        allocate(vyc (0:nxc+1,0:nyc+1))
cc
ccc     Restrict array sff -> sfc
cc
cc        call restrictArraytoArray(nxf,nyf,sff,igridf
cc     .                           ,nxc,nyc,sfc,igrid ,order,bcond)
cc
ccc     Form vector components in coarse grid
cc
cc        call vec_find(nxc,nyc,dx1,dy1,sfc,vxc,vyc)
cc
ccc     Map vector components onto MG vector
cc
cc        call mapArrayToMGVector(nxc,nyc,vxc,vxx,igrid)
cc        call mapArrayToMGVector(nxc,nyc,vyc,vyy,igrid)
cc
ccc     Transfer grid information
cc
cc        if (order.eq.0) then
cc          igridf = igrid
cc          nxf = nxc
cc          nyf = nyc
cc
cc          deallocate(sff)
cc          allocate(sff(0:nxf+1,0:nyf+1))
cc
cc          sff = sfc
cc        endif
cc
ccc     Deallocate variables
cc
cc        deallocate(sfc,vxc,vyc)
cc
cc      enddo
cc
ccc End program
cc
cc      deallocate(sff)
cc
cc      return
cc      end
