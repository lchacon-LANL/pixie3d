c  module matvec
c ###################################################################
      module matvec

        use grid

        use setMGBC_interface

        use mg_internal

        use constants

      contains

c     symm_test
c     ###############################################################
      subroutine symm_test(neq,igrid,matvec,bcnd)
c     ---------------------------------------------------------------
c     Performs symmetry test of matvec on grid "igrid".
c     ---------------------------------------------------------------

      implicit none     !For safe fortran

c     Call variables

      integer(4) :: neq,igrid,bcnd(6,neq)

      external      matvec

c     Local variables

      real(8),allocatable,dimension(:)::x1,dummy,dummy2

      real(8)    :: dd1,dd2,error
      integer(4) :: nx,ny,nz,nn,ii,jj,i1,j1,i2,j2,ix1,iy1,ix2,iy2,ieq

c     Begin program

c     Initialize variables

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      nn = neq*nx*ny*nz

      write (*,*) 'Performing symmetry test of system matrix ',
     .            'on grid:',nx,'x',ny,'x',nz,'...'

      allocate(x1(nn),dummy(nn),dummy2(nn))

      error = 0d0

      do ii = 1,nn
        x1    (ii) = 0d0
        dummy (ii) = 0d0
        dummy2(ii) = 0d0
      enddo

c     Check symmetry

      do ii = 1,nn/neq

        do ieq =1,neq

c       Find column vector ii

          call findBaseVector(ii,ieq,neq,nn,x1,1d0)

          call matvec(0,nn,x1,dummy,igrid,bcnd)

          call findBaseVector(ii,ieq,neq,nn,x1,0d0)

c       Compare column vector ii with corresponding row vector (intersect in
c       diagonal)

          do jj = ii,nn/neq

            call findBaseVector(jj,ieq,neq,nn,x1,1d0)

            call matvec(ii,nn,x1,dummy2,igrid,bcnd)

            call findBaseVector(jj,ieq,neq,nn,x1,0d0)

            dd1 = abs(dummy(jj) - dummy2(ii))
            if(abs(dummy(jj)).gt.1d-15.or.abs(dummy2(ii)).gt.1d-15) then
              write(*,15) jj,ii,dummy(jj),ii,jj,dummy2(ii),dd1
     .             ,100*dd1/max(abs(dummy(jj)),abs(dummy2(ii)))
              error = error + dd1
            endif

          enddo

        enddo
      enddo

      write (*,20) error

      stop

c     End program

      deallocate (x1,dummy,dummy2)

 15   format ('(',i3,',',i3,'):',1pe10.2,'; (',i3,',',i3,'):',e10.2,
     .        '  Error:',e10.2,'  %error:',0pf7.2)
 20   format (/,'Total relative error:',1pe10.3)

      end subroutine symm_test

      end module matvec

ccc module mgarraySetup
ccc ######################################################################
cc      module mgarraySetup
cc
cccc        use parameters
cccc
cccc        use grid
cc
cc        use mg_internal
cc
cc        use setMGBC_interface
cc
cccc        use equilibrium
cccc
cccc        use grid_aliases
cccc
cccc        use matvec
cc
cccc        integer(4) :: icmp   !Passed to setMGBC to define appropriate BC operation
cc
cc        type :: garray
cc          real(8),pointer,dimension(:,:,:,:) :: array
cc        end type garray
cc
cc        type :: mg_array
cc          type(garray),pointer,dimension(:) :: grid
cc        end type mg_array
cc
cc        logical :: is__cnv
cc
cc      contains
cc
ccc     allocateMGArray
ccc     #################################################################
cc      subroutine allocateMGarray(neq,mgarray)
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4)      :: neq
cc        type(mg_array)  :: mgarray
cc
ccc     Local variables
cc
cc        integer(4)      :: igrid,nxp,nyp,nzp
cc
ccc     Begin program
cc
cc        if (.not.associated(mgarray%grid)) then
cc          allocate(mgarray%grid(grid_params%ngrid))
cc          do igrid=1,grid_params%ngrid
cc            nxp = grid_params%nxv(igrid)+1
cc            nyp = grid_params%nyv(igrid)+1
cc            nzp = grid_params%nzv(igrid)+1
cc            allocate(mgarray%grid(igrid)%array(0:nxp,0:nyp,0:nzp,neq))
cc            mgarray%grid(igrid)%array = 0d0
cc          enddo
cc        endif
cc
cccc        do igrid=1,grid_params%ngrid
cccc          if (.not.associated(mgarray%grid(igrid)%array)) then
cccc            nxp = grid_params%nxv(igrid)+1
cccc            nyp = grid_params%nyv(igrid)+1
cccc            nzp = grid_params%nzv(igrid)+1
cccc            allocate(mgarray%grid(igrid)%array(0:nxp,0:nyp,0:nzp,neq))
cccc            mgarray%grid(igrid)%array = 0d0
cccc          endif
cccc        enddo
cc
ccc     End program
cc
cc      end subroutine allocateMGarray
cc
ccc     deallocateMGArray
ccc     #################################################################
cc      subroutine deallocateMGArray(mgarray)
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        type(mg_array)  :: mgarray
cc
ccc     Local variables
cc
cc        integer          :: igrid
cc
ccc     Begin program
cc
cc        if (associated(mgarray%grid)) then
cc          do igrid=1,grid_params%ngrid
cc            if (associated(mgarray%grid(igrid)%array)) then
cc              deallocate(mgarray%grid(igrid)%array)
cc            endif
cc          enddo
cc          deallocate(mgarray%grid)
cc        endif
cc
ccc     End program
cc
cc      end subroutine deallocateMGArray
cc
ccc     restrictMGArray
ccc     #################################################################
cc      subroutine restrictMGArray(icmp,neq,mgarray,bcnd,igrid,order
cc     .                          ,iscnv)
ccc     -----------------------------------------------------------------
ccc     Restricts MG array in all grids with ghost nodes.
ccc     -----------------------------------------------------------------
cc
cc      implicit none    !For safe fortran
cc
ccc     Call variables
cc
cc      integer(4)     :: neq,icmp,bcnd(6,neq),order,igrid
cc      type(mg_array) :: mgarray
cc      logical,optional,intent(IN) :: iscnv
cc
ccc     Local variables
cc
cc      integer(4)     :: igf,nxf,nyf,nzf,igc,nxc,nyc,nzc
cc
ccc     Begin program
cc
cc      if (PRESENT(iscnv)) then
cc        is__cnv = iscnv
cc      else
cc        is__cnv = .true.   !Contravariant representation by default
cc      endif
cc
ccc     Consistency check
cc
cc      if (size(mgarray%grid(igrid)%array,4) /= neq) then
cc        write (*,*) 'Cannot restrict MG array: ',
cc     .              'inconsistent number of vector components'
cc        write (*,*) neq,size(mgarray%grid)
cc        write (*,*) 'Aborting...'
cc        stop
cc      endif
cc
ccc     Restrict array
cc
cc      do igc=igrid+1,grid_params%ngrid
cc        igf = igc-1
cc
cc        nxf = grid_params%nxv(igf)
cc        nyf = grid_params%nyv(igf)
cc        nzf = grid_params%nzv(igf)
cc        nxc = grid_params%nxv(igc)
cc        nyc = grid_params%nyv(igc)
cc        nzc = grid_params%nzv(igc)
cc
cc        call restrictArrayToArray(icmp,neq
cc     .       ,igf,nxf,nyf,nzf,mgarray%grid(igf)%array
cc     .       ,igc,nxc,nyc,nzc,mgarray%grid(igc)%array
cc     .       ,order,.false.,bcnd)
cc      enddo
cc
ccc     Reset grid level quantities (modified in setMGBC within restrictArrayToArray)
cc
cccc      igx = igrid
cccc      igy = igrid
cccc      igz = igrid
cccc
cccc      nx = grid_params%nxv(igrid)
cccc      ny = grid_params%nyv(igrid)
cccc      nz = grid_params%nzv(igrid)
cc
ccc     End program
cc
cc      end subroutine restrictMGArray
cc
ccc     restrictArrayToArray
ccc     #################################################################
cc      subroutine restrictArrayToArray(icmp,neq,igf,nxf,nyf,nzf,arrayf
cc     .                                        ,igc,nxc,nyc,nzc,arrayc
cc     .                               ,order,volf,bcnd)
ccc     -----------------------------------------------------------------
ccc     Restricts array to array in all grids (with ghost nodes),
ccc     starting at grid igf.
ccc     -----------------------------------------------------------------
cc
cc      implicit none    !For safe fortran
cc
ccc     Call variables
cc
cc      integer(4) :: neq,igf,nxf,nyf,nzf,igc,nxc,nyc,nzc
cc     .             ,order,bcnd(6,neq),icmp
cc      real(8)    :: arrayf(0:nxf+1,0:nyf+1,0:nzf+1,neq)
cc     .             ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1,neq)
cc      logical    :: volf
cc
ccc     Local variables
cc
cc      integer(4) :: igridf,igridc,isigf,isigc,i,j,k,ii,ieq
cc     .             ,nxxf,nyyf,nzzf,nxxc,nyyc,nzzc,ntotc,ntotf
cc     .             ,bcmod(6,neq)
cc
cc      real(8),allocatable,dimension(:) :: vecf,vecc
cc      logical    :: fpointers
cc
ccc     Begin program
cc
cc      call allocPointers(neq,fpointers)
cc
ccc     Consistency check
cc
cc      nxxf = nxv(igf)
cc      nyyf = nyv(igf)
cc      nzzf = nzv(igf)
cc
cc      if (nxf /= nxxf .or. nyf /= nyyf .or. nzf /= nzzf) then
cc        write (*,*) 'Grid mismatch in restrictArrayToArray:'
cc        write (*,*) 'Aborting...'
cc        stop
cc      endif
cc
ccc     Allocate vectors
cc
cc      nxxc = nxxf
cc      nyyc = nyyf
cc      nzzc = nzzf
cc
cc      ntotf = neq*nxxf*nyyf*nzzf
cc      ntotc = ntotf
cc
cc      allocate(vecf(ntotf))
cc      allocate(vecc(ntotc))
cc
ccc     Map arrays onto MG vector
cc
cc      !Set igrid=1 since vecf is NOT a MG vector
cc      call mapArrayToMGVector(neq,nxxf,nyyf,nzzf,arrayf,vecc,1)
cc
ccc     Restrict MG vectors
cc
cc      do igridc = igf+1,igc
cc
cc        igridf = igridc-1
cc
ccc       Characterize coarse and fine grids
cc
cc        nxxf = grid_params%nxv(igridf)
cc        nyyf = grid_params%nyv(igridf)
cc        nzzf = grid_params%nzv(igridf)
cc
cc        nxxc = grid_params%nxv(igridc)
cc        nyyc = grid_params%nyv(igridc)
cc        nzzc = grid_params%nzv(igridc)
cc
cc        ntotf = neq*nxxf*nyyf*nzzf
cc        ntotc = neq*nxxc*nyyc*nzzc
cc
ccc       Allocate coarse mesh vector
cc
cc        deallocate(vecf)
cc        allocate(vecf(ntotf))
cc
cc        vecf = vecc
cc
cc        deallocate(vecc)
cc        allocate(vecc(ntotc))
cc
ccc       Restrict vector
cc
cc        call crestrict(neq,vecc,ntotc,nxxc,nyyc,nzzc
cc     .                    ,vecf,ntotf,nxxf,nyyf,nzzf
cc     .                ,order,igridf,volf)
cc
cc      enddo
cc
ccc     Map vector to array
cc
cc      call mapMGVectorToArray(0,neq,vecc,nxc,nyc,nzc,arrayc,igc,.false.)
cc
cc      if (icmp /= 0) then
cc        bcmod = bcnd
cc        where (bcnd == EQU)
cc          bcmod = EXT
cc        end where
cc        call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcmod,icomp=icmp
cc     .              ,is_cnv=is__cnv)
cc      endif
cc
ccc     Deallocate vectors
cc
cc      deallocate(vecf,vecc)
cc
cc      call deallocPointers(fpointers)
cc
ccc     End program
cc
cc      end subroutine restrictArrayToArray
cc
ccc     restrictArrayToMGVector
ccc     #################################################################
cc      subroutine restrictArrayToMGVector(neq,nx,ny,nz,array,mgvector
cc     .                                  ,igr0,order,volf)
ccc     -----------------------------------------------------------------
ccc     Restricts array to mgvector in all grids (without ghost nodes),
ccc     starting at grid igr0.
ccc     -----------------------------------------------------------------
cc
cc      implicit none    !For safe fortran
cc
ccc     Call variables
cc
cc      integer(4) :: neq,nx,ny,nz,order
cc      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1,neq)
cc      real(8)    :: mgvector(*)
cc      logical    :: volf
cc
ccc     Local variables
cc
cc      integer(4) :: ieq,nxf,nyf,nzf,nxc,nyc,nzc,igridc,igridf,igr0
cc     .             ,isigf,isigc,ntotc,ntotf
cc      logical    :: fpointers
cc
ccc     Begin program
cc
cc      call allocPointers(neq,fpointers)
cc
ccc     Consistency check
cc
cc      nxf = nxv(igr0)
cc      nyf = nyv(igr0)
cc      nzf = nzv(igr0)
cc
cc      if (nxf /= nx .or. nyf /= ny .or. nzf /= nz) then
cc        write (*,*) 'Grid mismatch in restrictArray:'
cc        write (*,*) 'Aborting...'
cc        stop
cc      endif
cc
ccc     Map array in initial grid onto MG vector
cc
cc      call mapArrayToMGVector(neq,nx,ny,nz,array,mgvector,igr0)
cc
ccc     Restrict array to coarser grids
cc
cc      do igridc = igr0+1,grid_params%ngrid
cc
cc        igridf = igridc-1
cc
ccc       Characterize coarse and fine grids
cc
cc        nxf = nxv(igridf)
cc        nyf = nyv(igridf)
cc        nzf = nzv(igridf)
cc
cc        nxc = nxv(igridc)
cc        nyc = nyv(igridc)
cc        nzc = nzv(igridc)
cc
cc        isigc = istart(igridc)
cc        isigf = istart(igridf)
cc
cc        ntotf = neq*nxf*nyf*nzf
cc        ntotc = neq*nxc*nyc*nzc
cc
ccc       Restrict MG vector
cc
cc        call crestrict(neq,mgvector(isigc),ntotc,nxc,nyc,nzc
cc     .                    ,mgvector(isigf),ntotf,nxf,nyf,nzf
cc     .                ,order,igridf,volf)
cc      enddo
cc
cc      call deallocPointers(fpointers)
cc
cc      end subroutine restrictArrayToMGVector
cc
cc      end module mgarraySetup

c module precond_variables
c ######################################################################
      module precond_variables

        use timeStepping

        use precond_setup

        use transport_params

        use equilibrium

        use parameters

        use constants

        use auxiliaryVariables

        use nlfunction_setup

        use operators

        use variables

        use grid

        use mgarraySetup

        integer(4) :: i,j,k,ig,jg,kg,ieq
        real(8)    :: x1,y1,z1
        logical    :: cartsn,covariant,to_cartsn,to_cnv

        integer(4) :: ntotd2p

        real(8), allocatable, dimension(:,:,:,:):: p0

        real(8), allocatable, dimension(:,:)    :: rho_diag,tmp_diag
     .                                            ,b_diag,v_diag

        real(8), allocatable, dimension(:,:)    :: mgj0cnv,mgadvdiffV0

        real(8), allocatable, dimension(:)      :: mgdivV0

        integer(4), allocatable, dimension(:,:) :: bcs

        type (var_array),target :: varray

        type (mg_array ),target :: gp0,gb0,gv0,grho0,gb0_cov

        logical :: form_diag=.true.,vol_wgt=.true.,gm_smooth=.false.

      contains

c     allocPrecVariables
c     ###################################################################
      subroutine allocPrecVariables

c     -------------------------------------------------------------------
c     Allocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer(4) :: alloc_stat

c     Begin program

        ntotd2p = 2*ntotd/neqd

        allocate (mgj0cnv    (ntotd2p,3))
        allocate (mgadvdiffV0(ntotd2p,3))
        allocate (mgdivV0    (ntotd2p))
        allocate (bcs(6,neqd+3))

        allocate (rho_diag(1,  ntotd2p)
     .           ,tmp_diag(1,  ntotd2p)
     .           ,  b_diag(3,3*ntotd2p)
     .           ,  v_diag(3,3*ntotd2p),STAT=alloc_stat)

        call allocateMGArray(1,gp0)
        call allocateMGArray(1,grho0)
        call allocateMGArray(3,gv0)
        call allocateMGArray(3,gb0)
        call allocateMGArray(3,gb0_cov)

c     End program

      end subroutine allocPrecVariables

c     deallocPrecVariables
c     ###################################################################
      subroutine deallocPrecVariables

c     -------------------------------------------------------------------
c     Deallocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

        deallocate (mgj0cnv)
        deallocate (mgdivV0)
        deallocate (mgadvdiffV0)
        deallocate (bcs)

        call deallocateMGArray(gp0)
        call deallocateMGArray(grho0)
        call deallocateMGArray(gb0)
        call deallocateMGArray(gb0_cov)
        call deallocateMGArray(gv0)

        call deallocateDerivedType(varray)

c     End program

      end subroutine deallocPrecVariables

c     findCoeffs
c     ###################################################################
      subroutine findCoeffs

c     -------------------------------------------------------------------
c     Finds coefficients for linearized systems in preconditioner
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer(4) :: order,nxx,nyy,nzz,igrid,ii,ivar,igr
        real(8)    :: dvol

        real(8), allocatable, dimension(:,:,:,:) :: vector,vel

c     Begin program

        nxx = grid_params%nxv(igx)
        nyy = grid_params%nyv(igy)
        nzz = grid_params%nzv(igz)

        nx = nxx
        ny = nyy
        nz = nzz

        igrid = igx

        order=0

c     Store density in all grids (w/o BCs)

cc        call restrictArrayToMGVector(1,nxx,nyy,nzz,rho,mgrho0,igrid,order
cc     .                            ,.false.)

        grho0%grid(igrid)%array(:,:,:,1) = rho

        call restrictMGArray(IRHO,1,grho0,bcs(:,IRHO),igrid,order)

c     Store magnetic field cnv components in all grids (w/ BCs)

        gb0%grid(igrid)%array(:,:,:,1) = bx
        gb0%grid(igrid)%array(:,:,:,2) = by
        gb0%grid(igrid)%array(:,:,:,3) = bz

        call restrictMGArray(IBX,3,gb0,bcs(:,IBX:IBZ),igrid,order)

c     Store magnetic field cov components in all grids (w/ BCs)

        gb0_cov%grid(igrid)%array(:,:,:,1) = bx_cov
        gb0_cov%grid(igrid)%array(:,:,:,2) = by_cov
        gb0_cov%grid(igrid)%array(:,:,:,3) = bz_cov

        !Cannot use restricMGArray, since this routine only works for
        !cnv vector components.
        do igr=igrid+1,grid_params%ngrid
          do k=0,grid_params%nzv(igr)+1
            do j=0,grid_params%nyv(igr)+1
              do i=0,grid_params%nxv(igr)+1
                call transformFromCurvToCurv(i,j,k,igr,igr,igr
     .                            ,gb0_cov%grid(igr)%array(i,j,k,1)
     .                            ,gb0_cov%grid(igr)%array(i,j,k,2)
     .                            ,gb0_cov%grid(igr)%array(i,j,k,3)
     .                            ,gb0    %grid(igr)%array(i,j,k,1)
     .                            ,gb0    %grid(igr)%array(i,j,k,2)
     .                            ,gb0    %grid(igr)%array(i,j,k,3)
     .                            ,.false.)
              enddo
            enddo
          enddo
        enddo

cc        call restrictMGArray(IBX,3,gb0_cov,bcs(:,IBX:IBZ),igrid,order)

c     Store ion velocity components in all grids (w/ BCs)

        gv0%grid(igrid)%array(:,:,:,1) = vx
        gv0%grid(igrid)%array(:,:,:,2) = vy
        gv0%grid(igrid)%array(:,:,:,3) = vz

        call restrictMGArray(IVX,3,gv0,bcs(:,IVX:IVZ),igrid,order)

c     Store current components in all grids (w/o BCs)

        call restrictArrayToMGVector(1,nxx,nyy,nzz,jx,mgj0cnv(:,1),igrid
     .                  ,order,.false.)
        call restrictArrayToMGVector(1,nxx,nyy,nzz,jy,mgj0cnv(:,2),igrid
     .                  ,order,.false.)
        call restrictArrayToMGVector(1,nxx,nyy,nzz,jz,mgj0cnv(:,3),igrid
     .                  ,order,.false.)

c     Find auxiliary quantities and store them in all grids

        !Velocity divergence (w/o BCs)
        allocate(divrgV(0:nxx+1,0:nyy+1,0:nzz+1))

        do k=1,nzz
          do j=1,nyy
            do i=1,nxx
              divrgV(i,j,k) = div(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                           ,vx,vy,vz)
            enddo
          enddo
        enddo

        call restrictArrayToMGVector(1,nxx,nyy,nzz,divrgV,mgdivV0,igrid
     .                            ,order,.false.)

        deallocate(divrgV)

        !pressure (w/ BCs)
        gp0%grid(igrid)%array(:,:,:,1) = 2.*rho*tmp

        call restrictMGArray(IRHO,1,gp0,bcs(:,IRHO),igrid,order)

        !v0/dt+theta(v0.grad(v0)-mu*veclap(v0))
        allocate(vector(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .          ,vel   (0:nxx+1,0:nyy+1,0:nzz+1,3))

        vel (:,:,:,1) = vx
        vel (:,:,:,2) = vy
        vel (:,:,:,3) = vz

        do k = 1,nzz
          do j = 1,nyy
            do i = 1,nxx
              ii  = i + nxx*(j-1) + nxx*nyy*(k-1)

              jac    = gmetric%grid(igrid)%jac(i,j,k)
              nabla_v= fnabla_v(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                         ,vx,vy,vz,0)

              !Eqn 1
              do ivar=1,3
                vector(i,j,k,ivar) = 
     .                            gv0%grid(igrid)%array(i,j,k,ivar)/dt
     .                          + alpha*vx(i,j,k)*nabla_v(1,ivar)/jac
     .                          + alpha*vy(i,j,k)*nabla_v(2,ivar)/jac
     .                          + alpha*vz(i,j,k)*nabla_v(3,ivar)/jac
     .                          - alpha*veclaplacian(i,j,k,nxx,nyy,nzz
     .                                              ,igrid,igrid,igrid
     .                                              ,vel,nuu,alt_eom
     .                                              ,ivar,vol=.false.)
cc     .                                              ,ivar)/dvol
              enddo

            enddo
          enddo
        enddo

        do ivar=1,3
          call restrictArrayToMGVector(1,nxx,nyy,nzz
     .                              ,vector(:,:,:,ivar)
     .                              ,mgadvdiffV0(:,ivar)
     .                              ,igrid,order,.false.)
        enddo

        deallocate(vector,vel)

c     End program

      end subroutine findCoeffs

c     find_curl_vxb
c     ###################################################################
      subroutine find_curl_vxb(i,j,k,nx,ny,nz,vv,bb,a1,a2,a3
     .                        ,half_elem,igrid)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of -curl(vv x bb) at the
c     grid node (i,j,k). One sided derivatives are employed when half_elem=1
c     (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nx,ny,nz,half_elem,igrid
        real(8)    :: a1,a2,a3,vv(0:nx+1,0:ny+1,0:nz+1,3)
     $                        ,bb(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq,igx,igy,igz
        integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

        real(8)    :: idhx,idhy,idhz
        real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .               ,jacp,jacm,jach,jac0

        real(8)    :: vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm

        real(8)    :: bxip,bxim,bxjp,bxjm,bxkp,bxkm
     .               ,byip,byim,byjp,byjm,bykp,bykm
     .               ,bzip,bzim,bzjp,bzjm,bzkp,bzkm

        logical    :: sing_point

c     Begin program

        igx = igrid
        igy = igrid
        igz = igrid

        sing_point = isSP(i,j,k,igx,igy,igz)

c     Defaults

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        idhx = 0.5/dxh(ig)
        idhy = 0.5/dyh(jg)
        idhz = 0.5/dzh(kg)

        jac  = gmetric%grid(igx)%jac(i,j,k)

c     Exceptions

        select case(half_elem)
        case (1)
          idhx = 1./dx(ig)
          im = i

          if (.not.isSP(ip,j,k,igx,igy,igz)) then  !i>0 and/or not a SP at i=0

            jacip  = gmetric%grid(igx)%jac(ip,j,k)
            jacim  = gmetric%grid(igx)%jac(i ,j,k)
            jacjp  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
     .                   +gmetric%grid(igx)%jac(i ,jp,k))
            jacjm  = 0.5*(gmetric%grid(igx)%jac(ip,jm,k)
     .                   +gmetric%grid(igx)%jac(i ,jm,k))
            jackp  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
     .                   +gmetric%grid(igx)%jac(i ,j,kp))
            jackm  = 0.5*(gmetric%grid(igx)%jac(ip,j,km)
     .                   +gmetric%grid(igx)%jac(i ,j,km))

            vxip = vv(ip,j,k,1)
            vxim = vv(i ,j,k,1)
            vyip = vv(ip,j,k,2)
            vyim = vv(i ,j,k,2)
            vzip = vv(ip,j,k,3)
            vzim = vv(i ,j,k,3)

            vxjp = (vv(ip,jp,k,1)+vv(i,jp,k,1))*0.5
            vxjm = (vv(ip,jm,k,1)+vv(i,jm,k,1))*0.5
            vyjp = (vv(ip,jp,k,2)+vv(i,jp,k,2))*0.5
            vyjm = (vv(ip,jm,k,2)+vv(i,jm,k,2))*0.5
            vzjp = (vv(ip,jp,k,3)+vv(i,jp,k,3))*0.5
            vzjm = (vv(ip,jm,k,3)+vv(i,jm,k,3))*0.5

            vxkp = (vv(ip,j,kp,1)+vv(i,j,kp,1))*0.5
            vxkm = (vv(ip,j,km,1)+vv(i,j,km,1))*0.5
            vykp = (vv(ip,j,kp,2)+vv(i,j,kp,2))*0.5
            vykm = (vv(ip,j,km,2)+vv(i,j,km,2))*0.5
            vzkp = (vv(ip,j,kp,3)+vv(i,j,kp,3))*0.5
            vzkm = (vv(ip,j,km,3)+vv(i,j,km,3))*0.5

            bxip = bb(ip,j,k,1)
            bxim = bb(i ,j,k,1)
            byip = bb(ip,j,k,2)
            byim = bb(i ,j,k,2)
            bzip = bb(ip,j,k,3)
            bzim = bb(i ,j,k,3)

            bxjp = (bb(ip,jp,k,1)+bb(i,jp,k,1))*0.5
            bxjm = (bb(ip,jm,k,1)+bb(i,jm,k,1))*0.5
            byjp = (bb(ip,jp,k,2)+bb(i,jp,k,2))*0.5
            byjm = (bb(ip,jm,k,2)+bb(i,jm,k,2))*0.5
            bzjp = (bb(ip,jp,k,3)+bb(i,jp,k,3))*0.5
            bzjm = (bb(ip,jm,k,3)+bb(i,jm,k,3))*0.5

            bxkp = (bb(ip,j,kp,1)+bb(i,j,kp,1))*0.5
            bxkm = (bb(ip,j,km,1)+bb(i,j,km,1))*0.5
            bykp = (bb(ip,j,kp,2)+bb(i,j,kp,2))*0.5
            bykm = (bb(ip,j,km,2)+bb(i,j,km,2))*0.5
            bzkp = (bb(ip,j,kp,3)+bb(i,j,kp,3))*0.5
            bzkm = (bb(ip,j,km,3)+bb(i,j,km,3))*0.5

          else  !i=0 and a SP at i=0: first order differences @ i=0

            jacip  = gmetric%grid(igx)%jac(ip,j,k)
            jacim  = gmetric%grid(igx)%jac(i ,j,k)
            jacjp  = gmetric%grid(igx)%jac(ip,jp,k)
            jacjm  = gmetric%grid(igx)%jac(ip,jm,k)
            jackp  = gmetric%grid(igx)%jac(ip,j,kp)
            jackm  = gmetric%grid(igx)%jac(ip,j,km)

            vxip = vv(ip,j,k,1)
            vxim = vv(i ,j,k,1)
            vyip = vv(ip,j,k,2)
            vyim = vv(i ,j,k,2)
            vzip = vv(ip,j,k,3)
            vzim = vv(i ,j,k,3)

            vxjp = vv(ip,jp,k,1)
            vxjm = vv(ip,jm,k,1)
            vyjp = vv(ip,jp,k,2)
            vyjm = vv(ip,jm,k,2)
            vzjp = vv(ip,jp,k,3)
            vzjm = vv(ip,jm,k,3)

            vxkp = vv(ip,j,kp,1)
            vxkm = vv(ip,j,km,1)
            vykp = vv(ip,j,kp,2)
            vykm = vv(ip,j,km,2)
            vzkp = vv(ip,j,kp,3)
            vzkm = vv(ip,j,km,3)

            bxip = bb(ip,j,k,1)
            bxim = bb(i ,j,k,1)
            byip = bb(ip,j,k,2)
            byim = bb(i ,j,k,2)
            bzip = bb(ip,j,k,3)
            bzim = bb(i ,j,k,3)

            bxjp = bb(ip,jp,k,1)
            bxjm = bb(ip,jm,k,1)
            byjp = bb(ip,jp,k,2)
            byjm = bb(ip,jm,k,2)
            bzjp = bb(ip,jp,k,3)
            bzjm = bb(ip,jm,k,3)

            bxkp = bb(ip,j,kp,1)
            bxkm = bb(ip,j,km,1)
            bykp = bb(ip,j,kp,2)
            bykm = bb(ip,j,km,2)
            bzkp = bb(ip,j,kp,3)
            bzkm = bb(ip,j,km,3)
          endif
        case (2)
          idhy = 1./dy(jg)
          jm = j

          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
     .                 +gmetric%grid(igx)%jac(ip,j ,k))
          jacim  = 0.5*(gmetric%grid(igx)%jac(im,jp,k)
     .                 +gmetric%grid(igx)%jac(im,j ,k))
          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
          jacjm  = gmetric%grid(igx)%jac(i,jm,k)
          jackp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
     .                 +gmetric%grid(igx)%jac(i,j ,kp))
          jackm  = 0.5*(gmetric%grid(igx)%jac(i,jp,km)
     .                 +gmetric%grid(igx)%jac(i,j ,km))

          if (sing_point) then !We avoid i=0 with first order differences
            idhx = 1./dx(ig)
            jacim= 0.5*(gmetric%grid(igx)%jac(i,jp,k)
     .                 +gmetric%grid(igx)%jac(i,j ,k))
            vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1))*0.5
            vxim = (vv(i ,j,k,1)+vv(i ,jp,k,1))*0.5
            vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2))*0.5
            vyim = (vv(i ,j,k,2)+vv(i ,jp,k,2))*0.5
            vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3))*0.5
            vzim = (vv(i ,j,k,3)+vv(i ,jp,k,3))*0.5

cc            vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1)
cc     .             +vv(ip,j,k,1)+vv(ip,jp,k,1))*0.25
cc            vxim = (vv(im,j,k,1)+vv(im,jp,k,1))
cc            vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2)
cc     .             +vv(ip,j,k,2)+vv(ip,jp,k,2))*0.25
cc            vyim = (vv(im,j,k,2)+vv(im,jp,k,2))
cc            vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3)
cc     .             +vv(ip,j,k,3)+vv(ip,jp,k,3))*0.25
cc            vzim = (vv(im,j,k,3)+vv(im,jp,k,3))
          else
            vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1))*0.5
            vxim = (vv(im,j,k,1)+vv(im,jp,k,1))*0.5
            vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2))*0.5
            vyim = (vv(im,j,k,2)+vv(im,jp,k,2))*0.5
            vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3))*0.5
            vzim = (vv(im,j,k,3)+vv(im,jp,k,3))*0.5
          endif

          vxjp = vv(i,jp,k,1)
          vxjm = vv(i,j ,k,1)
          vyjp = vv(i,jp,k,2)
          vyjm = vv(i,j ,k,2)
          vzjp = vv(i,jp,k,3)
          vzjm = vv(i,j ,k,3)

          vxkp = (vv(i,j,kp,1)+vv(i,jp,kp,1))*0.5
          vxkm = (vv(i,j,km,1)+vv(i,jp,km,1))*0.5
          vykp = (vv(i,j,kp,2)+vv(i,jp,kp,2))*0.5
          vykm = (vv(i,j,km,2)+vv(i,jp,km,2))*0.5
          vzkp = (vv(i,j,kp,3)+vv(i,jp,kp,3))*0.5
          vzkm = (vv(i,j,km,3)+vv(i,jp,km,3))*0.5

          if (sing_point) then !We avoid i=0 with first order differences
            bxip = (bb(ip,j,k,1)+bb(ip,jp,k,1))*0.5
            bxim = (bb(i ,j,k,1)+bb(i ,jp,k,1))*0.5
            byip = (bb(ip,j,k,2)+bb(ip,jp,k,2))*0.5
            byim = (bb(i ,j,k,2)+bb(i ,jp,k,2))*0.5
            bzip = (bb(ip,j,k,3)+bb(ip,jp,k,3))*0.5
            bzim = (bb(i ,j,k,3)+bb(i ,jp,k,3))*0.5

cc            bxip = (bb(ip,j,k,1)+bb(ip,jp,k,1)
cc     .             +bb(ip,j,k,1)+bb(ip,jp,k,1))*0.25
cc            bxim = (bb(im,j,k,1)+bb(im,jp,k,1))
cc            byip = (bb(ip,j,k,2)+bb(ip,jp,k,2)
cc     .             +bb(ip,j,k,2)+bb(ip,jp,k,2))*0.25
cc            byim = (bb(im,j,k,2)+bb(im,jp,k,2))
cc            bzip = (bb(ip,j,k,3)+bb(ip,jp,k,3)
cc     .             +bb(ip,j,k,3)+bb(ip,jp,k,3))*0.25
cc            bzim = (bb(im,j,k,3)+bb(im,jp,k,3))
          else
            bxip = (bb(ip,j,k,1)+bb(ip,jp,k,1))*0.5
            bxim = (bb(im,j,k,1)+bb(im,jp,k,1))*0.5
            byip = (bb(ip,j,k,2)+bb(ip,jp,k,2))*0.5
            byim = (bb(im,j,k,2)+bb(im,jp,k,2))*0.5
            bzip = (bb(ip,j,k,3)+bb(ip,jp,k,3))*0.5
            bzim = (bb(im,j,k,3)+bb(im,jp,k,3))*0.5
          endif

          bxjp = bb(i,jp,k,1)
          bxjm = bb(i,j ,k,1)
          byjp = bb(i,jp,k,2)
          byjm = bb(i,j ,k,2)
          bzjp = bb(i,jp,k,3)
          bzjm = bb(i,j ,k,3)

          bxkp = (bb(i,j,kp,1)+bb(i,jp,kp,1))*0.5
          bxkm = (bb(i,j,km,1)+bb(i,jp,km,1))*0.5
          bykp = (bb(i,j,kp,2)+bb(i,jp,kp,2))*0.5
          bykm = (bb(i,j,km,2)+bb(i,jp,km,2))*0.5
          bzkp = (bb(i,j,kp,3)+bb(i,jp,kp,3))*0.5
          bzkm = (bb(i,j,km,3)+bb(i,jp,km,3))*0.5

        case (3)
          idhz = 1./dz(kg)
          km = k

          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
     .                 +gmetric%grid(igx)%jac(ip,j,k ))
          jacim  = 0.5*(gmetric%grid(igx)%jac(im,j,kp)
     .                 +gmetric%grid(igx)%jac(im,j,k ))
          jacjp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
     .                 +gmetric%grid(igx)%jac(i,jp,k ))
          jacjm  = 0.5*(gmetric%grid(igx)%jac(i,jm,kp)
     .                 +gmetric%grid(igx)%jac(i,jm,k ))
          jackp  = gmetric%grid(igx)%jac(i,j,kp)
          jackm  = gmetric%grid(igx)%jac(i,j,km)

          if (sing_point) then !We avoid i=0 with first order differences
            idhx = 1./dx(ig)
            jacim= 0.5*(gmetric%grid(igx)%jac(i,j,kp)
     .                 +gmetric%grid(igx)%jac(i,j,k ))
            vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1))*0.5
            vxim = (vv(i ,j,k,1)+vv(i ,j,kp,1))*0.5
            vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2))*0.5
            vyim = (vv(i ,j,k,2)+vv(i ,j,kp,2))*0.5
            vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3))*0.5
            vzim = (vv(i ,j,k,3)+vv(i ,j,kp,3))*0.5

cc            vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1)
cc     .             +vv(i ,j,k,1)+vv(i ,j,kp,1))*0.25
cc            vxim = (vv(im,j,k,1)+vv(im,j,kp,1))
cc            vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2)
cc     .             +vv(i ,j,k,2)+vv(i ,j,kp,2))*0.25
cc            vyim = (vv(im,j,k,2)+vv(im,j,kp,2))
cc            vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3)
cc     .             +vv(i ,j,k,3)+vv(i ,j,kp,3))*0.25
cc            vzim = (vv(im,j,k,3)+vv(im,j,kp,3))

          else
            vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1))*0.5
            vxim = (vv(im,j,k,1)+vv(im,j,kp,1))*0.5
            vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2))*0.5
            vyim = (vv(im,j,k,2)+vv(im,j,kp,2))*0.5
            vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3))*0.5
            vzim = (vv(im,j,k,3)+vv(im,j,kp,3))*0.5
          endif

          vxjp = (vv(i,jp,k,1)+vv(i,jp,kp,1))*0.5
          vxjm = (vv(i,jm,k,1)+vv(i,jm,kp,1))*0.5
          vyjp = (vv(i,jp,k,2)+vv(i,jp,kp,2))*0.5
          vyjm = (vv(i,jm,k,2)+vv(i,jm,kp,2))*0.5
          vzjp = (vv(i,jp,k,3)+vv(i,jp,kp,3))*0.5
          vzjm = (vv(i,jm,k,3)+vv(i,jm,kp,3))*0.5

          vxkp = vv(i,j,kp,1)
          vxkm = vv(i,j,k ,1)
          vykp = vv(i,j,kp,2)
          vykm = vv(i,j,k ,2)
          vzkp = vv(i,j,kp,3)
          vzkm = vv(i,j,k ,3)

          if (sing_point) then !We avoid i=0 with first order differences
            bxip = (bb(ip,j,k,1)+bb(ip,j,kp,1))*0.5
            bxim = (bb(i ,j,k,1)+bb(i ,j,kp,1))*0.5
            byip = (bb(ip,j,k,2)+bb(ip,j,kp,2))*0.5
            byim = (bb(i ,j,k,2)+bb(i ,j,kp,2))*0.5
            bzip = (bb(ip,j,k,3)+bb(ip,j,kp,3))*0.5
            bzim = (bb(i ,j,k,3)+bb(i ,j,kp,3))*0.5

cc            bxip = (bb(ip,j,k,1)+bb(ip,j,kp,1)
cc     .             +bb(i ,j,k,1)+bb(i ,j,kp,1))*0.25
cc            bxim = (bb(im,j,k,1)+bb(im,j,kp,1))
cc            byip = (bb(ip,j,k,2)+bb(ip,j,kp,2)
cc     .             +bb(i ,j,k,2)+bb(i ,j,kp,2))*0.25
cc            byim = (bb(im,j,k,2)+bb(im,j,kp,2))
cc            bzip = (bb(ip,j,k,3)+bb(ip,j,kp,3)
cc     .             +bb(i ,j,k,3)+bb(i ,j,kp,3))*0.25
cc            bzim = (bb(im,j,k,3)+bb(im,j,kp,3))
          else
            bxip = (bb(ip,j,k,1)+bb(ip,j,kp,1))*0.5
            bxim = (bb(im,j,k,1)+bb(im,j,kp,1))*0.5
            byip = (bb(ip,j,k,2)+bb(ip,j,kp,2))*0.5
            byim = (bb(im,j,k,2)+bb(im,j,kp,2))*0.5
            bzip = (bb(ip,j,k,3)+bb(ip,j,kp,3))*0.5
            bzim = (bb(im,j,k,3)+bb(im,j,kp,3))*0.5
          endif

          bxjp = (bb(i,jp,k,1)+bb(i,jp,kp,1))*0.5
          bxjm = (bb(i,jm,k,1)+bb(i,jm,kp,1))*0.5
          byjp = (bb(i,jp,k,2)+bb(i,jp,kp,2))*0.5
          byjm = (bb(i,jm,k,2)+bb(i,jm,kp,2))*0.5
          bzjp = (bb(i,jp,k,3)+bb(i,jp,kp,3))*0.5
          bzjm = (bb(i,jm,k,3)+bb(i,jm,kp,3))*0.5

          bxkp = bb(i,j,kp,1)
          bxkm = bb(i,j,k ,1)
          bykp = bb(i,j,kp,2)
          bykm = bb(i,j,k ,2)
          bzkp = bb(i,j,kp,3)
          bzkm = bb(i,j,k ,3)

        case default

          jacip  = gmetric%grid(igx)%jac(ip,j,k)
          jacim  = gmetric%grid(igx)%jac(im,j,k)
          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
          jacjm  = gmetric%grid(igx)%jac(i,jm,k)
          jackp  = gmetric%grid(igx)%jac(i,j,kp)
          jackm  = gmetric%grid(igx)%jac(i,j,km)

          !Velocity
          if (sing_point) then
cc            vxip = vv(ip,j,k,1)+vv(i,j,k,1)
cc            vxim = 2.*vv(im,j,k,1)
cc            vyip = vv(ip,j,k,2)+vv(i,j,k,2)
cc            vyim = 2.*vv(im,j,k,2)
cc            vzip = vv(ip,j,k,3)+vv(i,j,k,3)
cc            vzim = 2.*vv(im,j,k,3)
            idhx = 1./dx(ig)
            jacim= gmetric%grid(igx)%jac(i,j,k)
            vxip = vv(ip,j,k,1)
            vxim = vv(i ,j,k,1)
            vyip = vv(ip,j,k,2)
            vyim = vv(i ,j,k,2)
            vzip = vv(ip,j,k,3)
            vzim = vv(i ,j,k,3)
          else
            vxip = vv(ip,j,k,1)
            vxim = vv(im,j,k,1)
            vyip = vv(ip,j,k,2)
            vyim = vv(im,j,k,2)
            vzip = vv(ip,j,k,3)
            vzim = vv(im,j,k,3)
          endif

          vxjp = vv(i,jp,k,1)
          vxjm = vv(i,jm,k,1)
          vyjp = vv(i,jp,k,2)
          vyjm = vv(i,jm,k,2)
          vzjp = vv(i,jp,k,3)
          vzjm = vv(i,jm,k,3)

          vxkp = vv(i,j,kp,1)
          vxkm = vv(i,j,km,1)
          vykp = vv(i,j,kp,2)
          vykm = vv(i,j,km,2)
          vzkp = vv(i,j,kp,3)
          vzkm = vv(i,j,km,3)

          !Magnetic field
          if (sing_point) then
cc            bxip = bb(ip,j,k,1)+bb(i,j,k,1)
cc            bxim = 2.*bb(im,j,k,1)
cc            byip = bb(ip,j,k,2)+bb(i,j,k,2)
cc            byim = 2.*bb(im,j,k,2)
cc            bzip = bb(ip,j,k,3)+bb(i,j,k,3)
cc            bzim = 2.*bb(im,j,k,3)
            bxip = bb(ip,j,k,1)
            bxim = bb(i ,j,k,1)
            byip = bb(ip,j,k,2)
            byim = bb(i ,j,k,2)
            bzip = bb(ip,j,k,3)
            bzim = bb(i ,j,k,3)
          else
            bxip = bb(ip,j,k,1)
            bxim = bb(im,j,k,1)
            byip = bb(ip,j,k,2)
            byim = bb(im,j,k,2)
            bzip = bb(ip,j,k,3)
            bzim = bb(im,j,k,3)
          endif

          bxjp = bb(i,jp,k,1)
          bxjm = bb(i,jm,k,1)
          byjp = bb(i,jp,k,2)
          byjm = bb(i,jm,k,2)
          bzjp = bb(i,jp,k,3)
          bzjm = bb(i,jm,k,3)

          bxkp = bb(i,j,kp,1)
          bxkm = bb(i,j,km,1)
          bykp = bb(i,j,kp,2)
          bykm = bb(i,j,km,2)
          bzkp = bb(i,j,kp,3)
          bzkm = bb(i,j,km,3)

        end select

c     Components

        !component 1

        flxjp = ( vyjp*bxjp-vxjp*byjp )/jacjp
        flxjm = ( vyjm*bxjm-vxjm*byjm )/jacjm

        flxkp = ( vzkp*bxkp-vxkp*bzkp )/jackp
        flxkm = ( vzkm*bxkm-vxkm*bzkm )/jackm

        a1 =  (flxjp-flxjm)*idhy
     .       +(flxkp-flxkm)*idhz

        !component 2

        flxip = ( vxip*byip-vyip*bxip )/jacip
        flxim = ( vxim*byim-vyim*bxim )/jacim

        flxkp = ( vzkp*bykp-vykp*bzkp )/jackp
        flxkm = ( vzkm*bykm-vykm*bzkm )/jackm

        a2 =  (flxip-flxim)*idhx
     .       +(flxkp-flxkm)*idhz

        !component 3

        flxip = ( vxip*bzip-vzip*bxip )/jacip
        flxim = ( vxim*bzim-vzim*bxim )/jacim

        flxjp = ( vyjp*bzjp-vzjp*byjp )/jacjp
        flxjm = ( vyjm*bzjm-vzjm*byjm )/jacjm

        a3 =  (flxip-flxim)*idhx
     .       +(flxjp-flxjm)*idhy

      end subroutine find_curl_vxb

c     cSolver
c #   #####################################################################
      subroutine cSolver(neq,ntotp,b,x,bcnd,igrid,out,guess
     $                  ,matvec,dg,ncolors,line_relax,gm_smth,cvrg_tst)
c     ---------------------------------------------------------------
c     This subroutine solves a coupled system of neq equations. 
c     In call sequence:
c       * neq: number of coupled equations
c       * ntotp: number of mesh points
c       * b: rhs
c       * x: solution
c       * bcnd: boundary condition defs.
c       * igrid: MG grid level (igrid=1 is finest level)
c       * out: level of output information
c       * guess: whether a non-trivial initial guess is provided
c               (iguess=1) or not (iguess=0)
c       * matvec (external): matrix-vector product definition.
c
c     Optional variables:
c       * dg: matrix neq*neq diagonal block (for stationary its).
c       * ncolors: number of colors in grid (for GS).
c       * line_relax: whether we want line relaxation
c       * gm_smth: whether we want GMRES as a smoother
c       * cvrg_tst: whether we want to perform a convergence test
c     ---------------------------------------------------------------

        use mlsolverSetup

        implicit none

c     Call variables

        integer(4) :: neq,ntotp,igrid,bcnd(6,neq),out,guess
        real(8)    :: x(ntotp,neq),b(ntotp,neq)

        real(8)   ,optional,intent(IN) :: dg(neq,2*neq*ntotp)
        integer(4),optional,intent(IN) :: ncolors
        logical   ,optional,intent(IN) :: line_relax,gm_smth,cvrg_tst

        external   matvec

c     Local variables

        integer(4) :: ntot
        real(8)    :: xi(ntotp*neq),bi(ntotp*neq)
        real(8), target :: ddg(neq,2*neq*ntotp)
        logical    :: gm_smooth,cvrg_test

c     Begin program

c     Process optional arguments

        if (PRESENT(dg)) ddg = dg

        if (PRESENT(gm_smth)) then
          gm_smooth = gm_smth
        else
          gm_smooth = .false.
        endif

        if (PRESENT(cvrg_tst)) then
          cvrg_test = cvrg_tst
        else
          cvrg_test = .false.
        endif

c     Interlace variables for coupled solve

        if (.not.cvrg_test) then
          do i=1,ntotp
            do ieq=1,neq
              xi(neq*(i-1)+ieq) = x(i,ieq)
              bi(neq*(i-1)+ieq) = b(i,ieq)
            enddo
          enddo
        else
          bi = 0d0
          call random_number(xi)
          guess = 1
        endif

c     Solve coupled MG

c       Initialize solver

        call solverInit

c       Upper_level solver options (MG)

        call solverOptionsInit

        solverOptions%tol      = mgtol
        solverOptions%vcyc     = maxvcyc
        solverOptions%igridmin = 2
        solverOptions%orderres = 2
        solverOptions%orderprol= 2
        solverOptions%mg_mu    = 1
        solverOptions%vol_res  = vol_wgt

        if (PRESENT(dg)) solverOptions%diag => ddg

        if (PRESENT(ncolors)) solverOptions%ncolors  = ncolors

        if (gm_smooth) then
        !With GMRES smoothing
          solverOptions%mg_coarse_solver_depth = 4 !GMRES, as defined below
          solverOptions%fdiag = .false.
        else
        !With JB smoothing
          solverOptions%mg_coarse_solver_depth = 3 !GMRES, as defined below
        endif

        !Vertex relaxation
cc        solverOptions%vertex_based_relax = .true.

        !Plane/line relaxation
        if (PRESENT(line_relax)) then
          solverOptions%mg_line_relax  = line_relax
          solverOptions%mg_line_nsweep = 1
          solverOptions%mg_line_vcyc   = 1
          solverOptions%mg_line_tol    = 1d-1
          solverOptions%mg_line_omega  = 1d0
          solverOptions%mg_line_coarse_solver_depth = 0
        endif

        call assembleSolverHierarchy('mg')

c       Next level solver (smoother)

        call solverOptionsInit

        solverOptions%iter    = nsweep
        solverOptions%tol     = mgtol

        if (.not.gm_smooth) then
cc          if (PRESENT(dg)) solverOptions%diag => ddg

cc          !GS
cc          solverOptions%omega   = 1d0
cc          if (PRESENT(ncolors)) solverOptions%ncolors = ncolors
cccc          solverOptions%vertex_based_relax = .true.
cc
cc          call assembleSolverHierarchy('gs')

          !JB
          if (solverOptions%mg_line_relax) then
            solverOptions%omega   = 0.6 !For line relax
          else
            solverOptions%omega   = 0.8 !For point relax
          endif

          call assembleSolverHierarchy('jb')

        else

          !GMRES
          solverOptions%krylov_subspace = nsweep
          solverOptions%stp_test        = 1 

          call assembleSolverHierarchy('gm')
          call assembleSolverHierarchy('id') !GMRES preconditioner

        endif

c       Coarsest grid solve for outer MG

        call solverOptionsInit

        solverOptions%tol             = 1d-5
        solverOptions%krylov_subspace = 1000
        solverOptions%iter            = 1000
        solverOptions%stp_test        = 1 

        call assembleSolverHierarchy('gm')
        call assembleSolverHierarchy('id') !GMRES preconditioner

c       Coarsest grid solve for inner line/plane MG

        call solverOptionsInit

        solverOptions%omega   = 1d0
        solverOptions%iter    = 100
        solverOptions%tol     = 1d-4

        call assembleSolverHierarchy('gs')

c       Invoke solver

        ntot=neq*ntotp
        call getSolver(neq,ntot,bi,xi,matvec,igrid,bcnd,guess,out,1)

c       Get output data

cc      call getSolverOptions(1)

c       Kill solver

        call solverKill

c     Unpack solution for output

        do i = 1,ntotp
          do ieq=1,neq
            x(i,ieq) = xi(neq*(i-1)+ieq)
          enddo
        enddo

c     End program

      end subroutine cSolver

      end module precond_variables
