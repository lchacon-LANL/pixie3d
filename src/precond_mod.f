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

c module mgarraySetup
c ######################################################################
      module mgarraySetup

        use parameters

        use grid

        use mg_internal

        use equilibrium

        use grid_aliases

        use matvec

cc        integer(4) :: icmp   !Passed to setMGBC to define appropriate BC operation

        type :: garray
          real(8),pointer,dimension(:,:,:,:) :: array
        end type garray

        type :: mg_array
          type(garray),pointer,dimension(:) :: grid
        end type mg_array

      contains

c     allocateMGArray
c     #################################################################
      subroutine allocateMGarray(neq,mgarray)

        implicit none

c     Call variables

        integer(4)      :: neq
        type(mg_array)  :: mgarray

c     Local variables

        integer(4)      :: igrid,nxp,nyp,nzp

c     Begin program

        if (.not.associated(mgarray%grid)) then
          allocate(mgarray%grid(grid_params%ngrid))
          do igrid=1,grid_params%ngrid
            nxp = grid_params%nxv(igrid)+1
            nyp = grid_params%nyv(igrid)+1
            nzp = grid_params%nzv(igrid)+1
            allocate(mgarray%grid(igrid)%array(0:nxp,0:nyp,0:nzp,neq))
            mgarray%grid(igrid)%array = 0d0
          enddo
        endif

cc        do igrid=1,grid_params%ngrid
cc          if (.not.associated(mgarray%grid(igrid)%array)) then
cc            nxp = grid_params%nxv(igrid)+1
cc            nyp = grid_params%nyv(igrid)+1
cc            nzp = grid_params%nzv(igrid)+1
cc            allocate(mgarray%grid(igrid)%array(0:nxp,0:nyp,0:nzp,neq))
cc            mgarray%grid(igrid)%array = 0d0
cc          endif
cc        enddo

c     End program

      end subroutine allocateMGarray

c     deallocateMGArray
c     #################################################################
      subroutine deallocateMGArray(mgarray)

        implicit none

c     Call variables

        type(mg_array)  :: mgarray

c     Local variables

        integer          :: igrid

c     Begin program

        if (associated(mgarray%grid)) then
          do igrid=1,grid_params%ngrid
            if (associated(mgarray%grid(igrid)%array)) then
              deallocate(mgarray%grid(igrid)%array)
            endif
          enddo
          deallocate(mgarray%grid)
        endif

c     End program

      end subroutine deallocateMGArray

c     restrictMGArray
c     #################################################################
      subroutine restrictMGArray(icmp,neq,mgarray,bcnd,igrid,order)
c     -----------------------------------------------------------------
c     Restricts MG array in all grids with ghost nodes.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer(4)     :: neq,icmp,bcnd(6,neq),order,igrid
      type(mg_array) :: mgarray

c     Local variables

      integer(4)     :: igf,nxf,nyf,nzf,igc,nxc,nyc,nzc

c     Begin program

c     Consistency check

      if (size(mgarray%grid(igrid)%array,4) /= neq) then
        write (*,*) 'Cannot restrict MG array: ',
     .              'inconsistent number of vector components'
        write (*,*) neq,size(mgarray%grid)
        write (*,*) 'Aborting...'
        stop
      endif

c     Restrict array

      do igc=igrid+1,grid_params%ngrid
        igf = igc-1

        nxf = grid_params%nxv(igf)
        nyf = grid_params%nyv(igf)
        nzf = grid_params%nzv(igf)
        nxc = grid_params%nxv(igc)
        nyc = grid_params%nyv(igc)
        nzc = grid_params%nzv(igc)

        call restrictArrayToArray(icmp,neq
     .       ,igf,nxf,nyf,nzf,mgarray%grid(igf)%array
     .       ,igc,nxc,nyc,nzc,mgarray%grid(igc)%array
     .       ,order,.false.,bcnd)
      enddo

c     Reset grid level quantities (modified in setMGBC within restrictArrayToArray)

      igx = igrid
      igy = igrid
      igz = igrid

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

c     End program

      end subroutine restrictMGArray

c     restrictArrayToArray
c     #################################################################
      subroutine restrictArrayToArray(icmp,neq,igf,nxf,nyf,nzf,arrayf
     .                                        ,igc,nxc,nyc,nzc,arrayc
     .                               ,order,volf,bcnd)
c     -----------------------------------------------------------------
c     Restricts array to array in all grids (with ghost nodes),
c     starting at grid igf.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer(4) :: neq,igf,nxf,nyf,nzf,igc,nxc,nyc,nzc
     .             ,order,bcnd(6,neq),icmp
      real(8)    :: arrayf(0:nxf+1,0:nyf+1,0:nzf+1,neq)
     .             ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1,neq)
      logical    :: volf

c     Local variables

      integer(4) :: igridf,igridc,isigf,isigc,i,j,k,ii,ieq
     .             ,nxxf,nyyf,nzzf,nxxc,nyyc,nzzc,ntotc,ntotf
     .             ,bcmod(6,neq)

      real(8),allocatable,dimension(:) :: vecf,vecc
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

c     Consistency check

      nxxf = nxv(igf)
      nyyf = nyv(igf)
      nzzf = nzv(igf)

      if (nxf /= nxxf .or. nyf /= nyyf .or. nzf /= nzzf) then
        write (*,*) 'Grid mismatch in restrictArrayToArray:'
        write (*,*) 'Aborting...'
        stop
      endif

c     Allocate vectors

      nxxc = nxxf
      nyyc = nyyf
      nzzc = nzzf

      ntotf = neq*nxxf*nyyf*nzzf
      ntotc = ntotf

      allocate(vecf(ntotf))
      allocate(vecc(ntotc))

c     Map arrays onto MG vector

      !Set igrid=1 since vecf is NOT a MG vector
      call mapArrayToMGVector(neq,nxxf,nyyf,nzzf,arrayf,vecc,1)

c     Restrict MG vectors

      do igridc = igf+1,igc

        igridf = igridc-1

c       Characterize coarse and fine grids

        nxxf = grid_params%nxv(igridf)
        nyyf = grid_params%nyv(igridf)
        nzzf = grid_params%nzv(igridf)

        nxxc = grid_params%nxv(igridc)
        nyyc = grid_params%nyv(igridc)
        nzzc = grid_params%nzv(igridc)

        ntotf = neq*nxxf*nyyf*nzzf
        ntotc = neq*nxxc*nyyc*nzzc

c       Allocate coarse mesh vector

        deallocate(vecf)
        allocate(vecf(ntotf))

        vecf = vecc

        deallocate(vecc)
        allocate(vecc(ntotc))

c       Restrict vector

        call crestrict(neq,vecc,ntotc,nxxc,nyyc,nzzc
     .                    ,vecf,ntotf,nxxf,nyyf,nzzf
     .                ,order,igridf,volf)

      enddo

c     Map vector to array

      call mapMGVectorToArray(0,neq,vecc,nxc,nyc,nzc,arrayc,igc,.false.)

      if (icmp /= 0) then
        bcmod = bcnd
        where (bcnd == EQU)
          bcmod = EXT
        end where
        call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcmod,icomp=icmp)
      endif

c     Deallocate vectors

      deallocate(vecf,vecc)

      call deallocPointers(fpointers)

c     End program

      end subroutine restrictArrayToArray

c     restrictArrayToMGVector
c     #################################################################
      subroutine restrictArrayToMGVector(neq,nx,ny,nz,array,mgvector
     .                                  ,igr0,order,volf)
c     -----------------------------------------------------------------
c     Restricts array to mgvector in all grids (without ghost nodes),
c     starting at grid igr0.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer(4) :: neq,nx,ny,nz,order
      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1,neq)
      real(8)    :: mgvector(*)
      logical    :: volf

c     Local variables

      integer(4) :: ieq,nxf,nyf,nzf,nxc,nyc,nzc,igridc,igridf,igr0
     .             ,isigf,isigc,ntotc,ntotf
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

c     Consistency check

      nxf = nxv(igr0)
      nyf = nyv(igr0)
      nzf = nzv(igr0)

      if (nxf /= nx .or. nyf /= ny .or. nzf /= nz) then
        write (*,*) 'Grid mismatch in restrictArray:'
        write (*,*) 'Aborting...'
        stop
      endif

c     Map array in initial grid onto MG vector

      call mapArrayToMGVector(neq,nx,ny,nz,array,mgvector,igr0)

c     Restrict array to coarser grids

      do igridc = igr0+1,grid_params%ngrid

        igridf = igridc-1

c       Characterize coarse and fine grids

        nxf = nxv(igridf)
        nyf = nyv(igridf)
        nzf = nzv(igridf)

        nxc = nxv(igridc)
        nyc = nyv(igridc)
        nzc = nzv(igridc)

        isigc = istart(igridc)
        isigf = istart(igridf)

        ntotf = neq*nxf*nyf*nzf
        ntotc = neq*nxc*nyc*nzc

c       Restrict MG vector

        call crestrict(neq,mgvector(isigc),ntotc,nxc,nyc,nzc
     .                    ,mgvector(isigf),ntotf,nxf,nyf,nzf
     .                ,order,igridf,volf)
      enddo

      call deallocPointers(fpointers)

      end subroutine restrictArrayToMGVector

      end module mgarraySetup

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

        type (mg_array ),target :: gp0,gb0,gv0,grho0

        logical :: form_diag=.true.,vol_wgt=.true.

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

        integer(4) :: order,nxx,nyy,nzz,igrid,ii,ivar
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

c     Store magnetic field components in all grids (w/ BCs)

        gb0%grid(igrid)%array(:,:,:,1) = bx
        gb0%grid(igrid)%array(:,:,:,2) = by
        gb0%grid(igrid)%array(:,:,:,3) = bz

        call restrictMGArray(IBX,3,gb0,bcs(:,IBX:IBZ),igrid,order)

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
            jacim= gmetric%grid(igx)%jac(i,jp,k)
            vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1))*0.5
            vxim = (vv(i ,j,k,1)+vv(i ,jp,k,1))*0.5
            vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2))*0.5
            vyim = (vv(i ,j,k,2)+vv(i ,jp,k,2))*0.5
            vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3))*0.5
            vzim = (vv(i ,j,k,3)+vv(i ,jp,k,3))*0.5
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
            jacim= gmetric%grid(igx)%jac(i,jp,k)
            vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1))*0.5
            vxim = (vv(i ,j,k,1)+vv(i ,j,kp,1))*0.5
            vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2))*0.5
            vyim = (vv(i ,j,k,2)+vv(i ,j,kp,2))*0.5
            vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3))*0.5
            vzim = (vv(i ,j,k,3)+vv(i ,j,kp,3))*0.5
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
            vxip = vv(ip,j,k,1)+vv(i,j,k,1)
            vxim = 2.*vv(im,j,k,1)
            vyip = vv(ip,j,k,2)+vv(i,j,k,2)
            vyim = 2.*vv(im,j,k,2)
            vzip = vv(ip,j,k,3)+vv(i,j,k,3)
            vzim = 2.*vv(im,j,k,3)
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
            bxip = bb(ip,j,k,1)+bb(i,j,k,1)
            bxim = 2.*bb(im,j,k,1)
            byip = bb(ip,j,k,2)+bb(i,j,k,2)
            byim = 2.*bb(im,j,k,2)
            bzip = bb(ip,j,k,3)+bb(i,j,k,3)
            bzim = 2.*bb(im,j,k,3)
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

ccc     find_curl_vxb
ccc     ###################################################################
cc      subroutine find_curl_vxb(i,j,k,nx,ny,nz,vv,bb,a1,a2,a3
cc     .                        ,half_elem,igrid)
cc
ccc     -------------------------------------------------------------------
ccc     Finds contravariant components (a1,a2,a3) of -curl(vv x bb) at the
ccc     grid node (i,j,k). One sided derivatives are employed when half_elem=1
ccc     (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
ccc     -------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,nx,ny,nz,half_elem,igrid
cc        real(8)    :: a1,a2,a3,vv(0:nx+1,0:ny+1,0:nz+1,3)
cc     $                        ,bb(0:nx+1,0:ny+1,0:nz+1,3)
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq,igx,igy,igz
cc        integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
cc
cc        real(8)    :: idhx,idhy,idhz
cc        real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc
cc        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
cc     .               ,jacp,jacm,jach,jac0
cc
cc        real(8)    :: vxip,vxim,vxjp,vxjm,vxkp,vxkm
cc     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
cc     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
cc
cc        real(8)    :: bxip,bxim,bxjp,bxjm,bxkp,bxkm
cc     .               ,byip,byim,byjp,byjm,bykp,bykm
cc     .               ,bzip,bzim,bzjp,bzjm,bzkp,bzkm
cc
cc        logical    :: sing_point
cc
ccc     Begin program
cc
cc        igx = igrid
cc        igy = igrid
cc        igz = igrid
cc
cc        sing_point = isSP(i,j,k,igx,igy,igz)
cc
ccc     Defaults
cc
cc        ip = i+1
cc        im = i-1
cc        jp = j+1
cc        jm = j-1
cc        kp = k+1
cc        km = k-1
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        idhx = 0.5/dxh(ig)
cc        idhy = 0.5/dyh(jg)
cc        idhz = 0.5/dzh(kg)
cc
cc        jac  = gmetric%grid(igx)%jac(i,j,k)
cc
ccc     Exceptions
cc
cc        select case(half_elem)
cc        case (1)
cc          idhx = 1./dx(ig)
cc          im = i
cc
cc          if (.not.isSP(ip,j,k,igx,igy,igz)) then  !i>0 and/or not a SP at i=0
cc
cc            jacip  = gmetric%grid(igx)%jac(ip,j,k)
cc            jacim  = gmetric%grid(igx)%jac(i ,j,k)
cc            jacjp  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
cc     .                   +gmetric%grid(igx)%jac(i ,jp,k))
cc            jacjm  = 0.5*(gmetric%grid(igx)%jac(ip,jm,k)
cc     .                   +gmetric%grid(igx)%jac(i ,jm,k))
cc            jackp  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
cc     .                   +gmetric%grid(igx)%jac(i ,j,kp))
cc            jackm  = 0.5*(gmetric%grid(igx)%jac(ip,j,km)
cc     .                   +gmetric%grid(igx)%jac(i ,j,km))
cc
cc            vxip = vv(ip,j,k,1)
cc            vxim = vv(i ,j,k,1)
cc            vyip = vv(ip,j,k,2)
cc            vyim = vv(i ,j,k,2)
cc            vzip = vv(ip,j,k,3)
cc            vzim = vv(i ,j,k,3)
cc
cc            vxjp = (vv(ip,jp,k,1)+vv(i,jp,k,1))*0.5
cc            vxjm = (vv(ip,jm,k,1)+vv(i,jm,k,1))*0.5
cc            vyjp = (vv(ip,jp,k,2)+vv(i,jp,k,2))*0.5
cc            vyjm = (vv(ip,jm,k,2)+vv(i,jm,k,2))*0.5
cc            vzjp = (vv(ip,jp,k,3)+vv(i,jp,k,3))*0.5
cc            vzjm = (vv(ip,jm,k,3)+vv(i,jm,k,3))*0.5
cc
cc            vxkp = (vv(ip,j,kp,1)+vv(i,j,kp,1))*0.5
cc            vxkm = (vv(ip,j,km,1)+vv(i,j,km,1))*0.5
cc            vykp = (vv(ip,j,kp,2)+vv(i,j,kp,2))*0.5
cc            vykm = (vv(ip,j,km,2)+vv(i,j,km,2))*0.5
cc            vzkp = (vv(ip,j,kp,3)+vv(i,j,kp,3))*0.5
cc            vzkm = (vv(ip,j,km,3)+vv(i,j,km,3))*0.5
cc
cc            bxip = bb(ip,j,k,1)
cc            bxim = bb(i ,j,k,1)
cc            byip = bb(ip,j,k,2)
cc            byim = bb(i ,j,k,2)
cc            bzip = bb(ip,j,k,3)
cc            bzim = bb(i ,j,k,3)
cc
cc            bxjp = (bb(ip,jp,k,1)+bb(i,jp,k,1))*0.5
cc            bxjm = (bb(ip,jm,k,1)+bb(i,jm,k,1))*0.5
cc            byjp = (bb(ip,jp,k,2)+bb(i,jp,k,2))*0.5
cc            byjm = (bb(ip,jm,k,2)+bb(i,jm,k,2))*0.5
cc            bzjp = (bb(ip,jp,k,3)+bb(i,jp,k,3))*0.5
cc            bzjm = (bb(ip,jm,k,3)+bb(i,jm,k,3))*0.5
cc
cc            bxkp = (bb(ip,j,kp,1)+bb(i,j,kp,1))*0.5
cc            bxkm = (bb(ip,j,km,1)+bb(i,j,km,1))*0.5
cc            bykp = (bb(ip,j,kp,2)+bb(i,j,kp,2))*0.5
cc            bykm = (bb(ip,j,km,2)+bb(i,j,km,2))*0.5
cc            bzkp = (bb(ip,j,kp,3)+bb(i,j,kp,3))*0.5
cc            bzkm = (bb(ip,j,km,3)+bb(i,j,km,3))*0.5
cc
cc          else  !i=0 and a SP at i=0: first order differences @ i=0
cc
cc            jacip  = gmetric%grid(igx)%jac(ip,j,k)
cc            jacim  = gmetric%grid(igx)%jac(i ,j,k)
cc            jacjp  = gmetric%grid(igx)%jac(ip,jp,k)
cc            jacjm  = gmetric%grid(igx)%jac(ip,jm,k)
cc            jackp  = gmetric%grid(igx)%jac(ip,j,kp)
cc            jackm  = gmetric%grid(igx)%jac(ip,j,km)
cc
cc            vxip = vv(ip,j,k,1)
cc            vxim = vv(i ,j,k,1)
cc            vyip = vv(ip,j,k,2)
cc            vyim = vv(i ,j,k,2)
cc            vzip = vv(ip,j,k,3)
cc            vzim = vv(i ,j,k,3)
cc
cc            vxjp = vv(ip,jp,k,1)
cc            vxjm = vv(ip,jm,k,1)
cc            vyjp = vv(ip,jp,k,2)
cc            vyjm = vv(ip,jm,k,2)
cc            vzjp = vv(ip,jp,k,3)
cc            vzjm = vv(ip,jm,k,3)
cc
cc            vxkp = vv(ip,j,kp,1)
cc            vxkm = vv(ip,j,km,1)
cc            vykp = vv(ip,j,kp,2)
cc            vykm = vv(ip,j,km,2)
cc            vzkp = vv(ip,j,kp,3)
cc            vzkm = vv(ip,j,km,3)
cc
cc            bxip = bb(ip,j,k,1)
cc            bxim = bb(i ,j,k,1)
cc            byip = bb(ip,j,k,2)
cc            byim = bb(i ,j,k,2)
cc            bzip = bb(ip,j,k,3)
cc            bzim = bb(i ,j,k,3)
cc
cc            bxjp = bb(ip,jp,k,1)
cc            bxjm = bb(ip,jm,k,1)
cc            byjp = bb(ip,jp,k,2)
cc            byjm = bb(ip,jm,k,2)
cc            bzjp = bb(ip,jp,k,3)
cc            bzjm = bb(ip,jm,k,3)
cc
cc            bxkp = bb(ip,j,kp,1)
cc            bxkm = bb(ip,j,km,1)
cc            bykp = bb(ip,j,kp,2)
cc            bykm = bb(ip,j,km,2)
cc            bzkp = bb(ip,j,kp,3)
cc            bzkm = bb(ip,j,km,3)
cc          endif
cc        case (2)
cc          idhy = 1./dy(jg)
cc          jm = j
cc
cc          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
cc     .                 +gmetric%grid(igx)%jac(ip,j ,k))
cc          jacim  = 0.5*(gmetric%grid(igx)%jac(im,jp,k)
cc     .                 +gmetric%grid(igx)%jac(im,j ,k))
cc          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
cc          jacjm  = gmetric%grid(igx)%jac(i,jm,k)
cc          jackp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
cc     .                 +gmetric%grid(igx)%jac(i,j ,kp))
cc          jackm  = 0.5*(gmetric%grid(igx)%jac(i,jp,km)
cc     .                 +gmetric%grid(igx)%jac(i,j ,km))
cc
cc          if (sing_point) then !We avoid i=0 with first order differences
cc            idhx = 1./dx(ig)
cc            jacim= gmetric%grid(igx)%jac(i,jp,k)
cc            vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1))*0.5
cc            vxim = (vv(i ,j,k,1)+vv(i ,jp,k,1))*0.5
cc            vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2))*0.5
cc            vyim = (vv(i ,j,k,2)+vv(i ,jp,k,2))*0.5
cc            vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3))*0.5
cc            vzim = (vv(i ,j,k,3)+vv(i ,jp,k,3))*0.5
cc          else
cc            vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1))*0.5
cc            vxim = (vv(im,j,k,1)+vv(im,jp,k,1))*0.5
cc            vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2))*0.5
cc            vyim = (vv(im,j,k,2)+vv(im,jp,k,2))*0.5
cc            vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3))*0.5
cc            vzim = (vv(im,j,k,3)+vv(im,jp,k,3))*0.5
cc          endif
cc
cc          vxjp = vv(i,jp,k,1)
cc          vxjm = vv(i,j ,k,1)
cc          vyjp = vv(i,jp,k,2)
cc          vyjm = vv(i,j ,k,2)
cc          vzjp = vv(i,jp,k,3)
cc          vzjm = vv(i,j ,k,3)
cc
cc          vxkp = (vv(i,j,kp,1)+vv(i,jp,kp,1))*0.5
cc          vxkm = (vv(i,j,km,1)+vv(i,jp,km,1))*0.5
cc          vykp = (vv(i,j,kp,2)+vv(i,jp,kp,2))*0.5
cc          vykm = (vv(i,j,km,2)+vv(i,jp,km,2))*0.5
cc          vzkp = (vv(i,j,kp,3)+vv(i,jp,kp,3))*0.5
cc          vzkm = (vv(i,j,km,3)+vv(i,jp,km,3))*0.5
cc
cc          if (sing_point) then !We avoid i=0 with first order differences
cc            bxip = (bb(ip,j,k,1)+bb(ip,jp,k,1))*0.5
cc            bxim = (bb(i ,j,k,1)+bb(i ,jp,k,1))*0.5
cc            byip = (bb(ip,j,k,2)+bb(ip,jp,k,2))*0.5
cc            byim = (bb(i ,j,k,2)+bb(i ,jp,k,2))*0.5
cc            bzip = (bb(ip,j,k,3)+bb(ip,jp,k,3))*0.5
cc            bzim = (bb(i ,j,k,3)+bb(i ,jp,k,3))*0.5
cc          else
cc            bxip = (bb(ip,j,k,1)+bb(ip,jp,k,1))*0.5
cc            bxim = (bb(im,j,k,1)+bb(im,jp,k,1))*0.5
cc            byip = (bb(ip,j,k,2)+bb(ip,jp,k,2))*0.5
cc            byim = (bb(im,j,k,2)+bb(im,jp,k,2))*0.5
cc            bzip = (bb(ip,j,k,3)+bb(ip,jp,k,3))*0.5
cc            bzim = (bb(im,j,k,3)+bb(im,jp,k,3))*0.5
cc          endif
cc
cc          bxjp = bb(i,jp,k,1)
cc          bxjm = bb(i,j ,k,1)
cc          byjp = bb(i,jp,k,2)
cc          byjm = bb(i,j ,k,2)
cc          bzjp = bb(i,jp,k,3)
cc          bzjm = bb(i,j ,k,3)
cc
cc          bxkp = (bb(i,j,kp,1)+bb(i,jp,kp,1))*0.5
cc          bxkm = (bb(i,j,km,1)+bb(i,jp,km,1))*0.5
cc          bykp = (bb(i,j,kp,2)+bb(i,jp,kp,2))*0.5
cc          bykm = (bb(i,j,km,2)+bb(i,jp,km,2))*0.5
cc          bzkp = (bb(i,j,kp,3)+bb(i,jp,kp,3))*0.5
cc          bzkm = (bb(i,j,km,3)+bb(i,jp,km,3))*0.5
cc
cc        case (3)
cc          idhz = 1./dz(kg)
cc          km = k
cc
cc          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
cc     .                 +gmetric%grid(igx)%jac(ip,j,k ))
cc          jacim  = 0.5*(gmetric%grid(igx)%jac(im,j,kp)
cc     .                 +gmetric%grid(igx)%jac(im,j,k ))
cc          jacjp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
cc     .                 +gmetric%grid(igx)%jac(i,jp,k ))
cc          jacjm  = 0.5*(gmetric%grid(igx)%jac(i,jm,kp)
cc     .                 +gmetric%grid(igx)%jac(i,jm,k ))
cc          jackp  = gmetric%grid(igx)%jac(i,j,kp)
cc          jackm  = gmetric%grid(igx)%jac(i,j,km)
cc
cc          if (sing_point) then !We avoid i=0 with first order differences
cc            idhx = 1./dx(ig)
cc            jacim= gmetric%grid(igx)%jac(i,jp,k)
cc            vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1))*0.5
cc            vxim = (vv(i ,j,k,1)+vv(i ,j,kp,1))*0.5
cc            vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2))*0.5
cc            vyim = (vv(i ,j,k,2)+vv(i ,j,kp,2))*0.5
cc            vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3))*0.5
cc            vzim = (vv(i ,j,k,3)+vv(i ,j,kp,3))*0.5
cc          else
cc            vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1))*0.5
cc            vxim = (vv(im,j,k,1)+vv(im,j,kp,1))*0.5
cc            vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2))*0.5
cc            vyim = (vv(im,j,k,2)+vv(im,j,kp,2))*0.5
cc            vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3))*0.5
cc            vzim = (vv(im,j,k,3)+vv(im,j,kp,3))*0.5
cc          endif
cc
cc          vxjp = (vv(i,jp,k,1)+vv(i,jp,kp,1))*0.5
cc          vxjm = (vv(i,jm,k,1)+vv(i,jm,kp,1))*0.5
cc          vyjp = (vv(i,jp,k,2)+vv(i,jp,kp,2))*0.5
cc          vyjm = (vv(i,jm,k,2)+vv(i,jm,kp,2))*0.5
cc          vzjp = (vv(i,jp,k,3)+vv(i,jp,kp,3))*0.5
cc          vzjm = (vv(i,jm,k,3)+vv(i,jm,kp,3))*0.5
cc
cc          vxkp = vv(i,j,kp,1)
cc          vxkm = vv(i,j,k ,1)
cc          vykp = vv(i,j,kp,2)
cc          vykm = vv(i,j,k ,2)
cc          vzkp = vv(i,j,kp,3)
cc          vzkm = vv(i,j,k ,3)
cc
cc          if (sing_point) then !We avoid i=0 with first order differences
cc            bxip = (bb(ip,j,k,1)+bb(ip,j,kp,1))*0.5
cc            bxim = (bb(i ,j,k,1)+bb(i ,j,kp,1))*0.5
cc            byip = (bb(ip,j,k,2)+bb(ip,j,kp,2))*0.5
cc            byim = (bb(i ,j,k,2)+bb(i ,j,kp,2))*0.5
cc            bzip = (bb(ip,j,k,3)+bb(ip,j,kp,3))*0.5
cc            bzim = (bb(i ,j,k,3)+bb(i ,j,kp,3))*0.5
cc          else
cc            bxip = (bb(ip,j,k,1)+bb(ip,j,kp,1))*0.5
cc            bxim = (bb(im,j,k,1)+bb(im,j,kp,1))*0.5
cc            byip = (bb(ip,j,k,2)+bb(ip,j,kp,2))*0.5
cc            byim = (bb(im,j,k,2)+bb(im,j,kp,2))*0.5
cc            bzip = (bb(ip,j,k,3)+bb(ip,j,kp,3))*0.5
cc            bzim = (bb(im,j,k,3)+bb(im,j,kp,3))*0.5
cc          endif
cc
cc          bxjp = (bb(i,jp,k,1)+bb(i,jp,kp,1))*0.5
cc          bxjm = (bb(i,jm,k,1)+bb(i,jm,kp,1))*0.5
cc          byjp = (bb(i,jp,k,2)+bb(i,jp,kp,2))*0.5
cc          byjm = (bb(i,jm,k,2)+bb(i,jm,kp,2))*0.5
cc          bzjp = (bb(i,jp,k,3)+bb(i,jp,kp,3))*0.5
cc          bzjm = (bb(i,jm,k,3)+bb(i,jm,kp,3))*0.5
cc
cc          bxkp = bb(i,j,kp,1)
cc          bxkm = bb(i,j,k ,1)
cc          bykp = bb(i,j,kp,2)
cc          bykm = bb(i,j,k ,2)
cc          bzkp = bb(i,j,kp,3)
cc          bzkm = bb(i,j,k ,3)
cc
cc        case default
cc
cc          jacip  = gmetric%grid(igx)%jac(ip,j,k)
cc          jacim  = gmetric%grid(igx)%jac(im,j,k)
cc          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
cc          jacjm  = gmetric%grid(igx)%jac(i,jm,k)
cc          jackp  = gmetric%grid(igx)%jac(i,j,kp)
cc          jackm  = gmetric%grid(igx)%jac(i,j,km)
cc
cc          !Velocity
cc          if (sing_point) then
cc            vxip = vv(ip,j,k,1)+vv(i,j,k,1)
cc            vxim = 2.*vv(im,j,k,1)
cc            vyip = vv(ip,j,k,2)+vv(i,j,k,2)
cc            vyim = 2.*vv(im,j,k,2)
cc            vzip = vv(ip,j,k,3)+vv(i,j,k,3)
cc            vzim = 2.*vv(im,j,k,3)
cc          else
cc            vxip = vv(ip,j,k,1)
cc            vxim = vv(im,j,k,1)
cc            vyip = vv(ip,j,k,2)
cc            vyim = vv(im,j,k,2)
cc            vzip = vv(ip,j,k,3)
cc            vzim = vv(im,j,k,3)
cc          endif
cc
cc          vxjp = vv(i,jp,k,1)
cc          vxjm = vv(i,jm,k,1)
cc          vyjp = vv(i,jp,k,2)
cc          vyjm = vv(i,jm,k,2)
cc          vzjp = vv(i,jp,k,3)
cc          vzjm = vv(i,jm,k,3)
cc
cc          vxkp = vv(i,j,kp,1)
cc          vxkm = vv(i,j,km,1)
cc          vykp = vv(i,j,kp,2)
cc          vykm = vv(i,j,km,2)
cc          vzkp = vv(i,j,kp,3)
cc          vzkm = vv(i,j,km,3)
cc
cc          !Magnetic field
cc          if (sing_point) then
cc            bxip = bb(ip,j,k,1)+bb(i,j,k,1)
cc            bxim = 2.*bb(im,j,k,1)
cc            byip = bb(ip,j,k,2)+bb(i,j,k,2)
cc            byim = 2.*bb(im,j,k,2)
cc            bzip = bb(ip,j,k,3)+bb(i,j,k,3)
cc            bzim = 2.*bb(im,j,k,3)
cc          else
cc            bxip = bb(ip,j,k,1)
cc            bxim = bb(im,j,k,1)
cc            byip = bb(ip,j,k,2)
cc            byim = bb(im,j,k,2)
cc            bzip = bb(ip,j,k,3)
cc            bzim = bb(im,j,k,3)
cc          endif
cc
cc          bxjp = bb(i,jp,k,1)
cc          bxjm = bb(i,jm,k,1)
cc          byjp = bb(i,jp,k,2)
cc          byjm = bb(i,jm,k,2)
cc          bzjp = bb(i,jp,k,3)
cc          bzjm = bb(i,jm,k,3)
cc
cc          bxkp = bb(i,j,kp,1)
cc          bxkm = bb(i,j,km,1)
cc          bykp = bb(i,j,kp,2)
cc          bykm = bb(i,j,km,2)
cc          bzkp = bb(i,j,kp,3)
cc          bzkm = bb(i,j,km,3)
cc
cc        end select
cc
ccc     Components
cc
cc        !component 1
cc
cc        flxjp = ( vyjp*bxjp-vxjp*byjp )/jacjp
cc        flxjm = ( vyjm*bxjm-vxjm*byjm )/jacjm
cc
cc        flxkp = ( vzkp*bxkp-vxkp*bzkp )/jackp
cc        flxkm = ( vzkm*bxkm-vxkm*bzkm )/jackm
cc
cc        a1 =  (flxjp-flxjm)*idhy
cc     .       +(flxkp-flxkm)*idhz
cc
cc        !component 2
cc
cc        flxip = ( vxip*byip-vyip*bxip )/jacip
cc        flxim = ( vxim*byim-vyim*bxim )/jacim
cc
cc        flxkp = ( vzkp*bykp-vykp*bzkp )/jackp
cc        flxkm = ( vzkm*bykm-vykm*bzkm )/jackm
cc
cc        a2 =  (flxip-flxim)*idhx
cc     .       +(flxkp-flxkm)*idhz
cc
cc        !component 3
cc
cc        flxip = ( vxip*bzip-vzip*bxip )/jacip
cc        flxim = ( vxim*bzim-vzim*bxim )/jacim
cc
cc        flxjp = ( vyjp*bzjp-vzjp*byjp )/jacjp
cc        flxjm = ( vyjm*bzjm-vzjm*byjm )/jacjm
cc
cc        a3 =  (flxip-flxim)*idhx
cc     .       +(flxjp-flxjm)*idhy
cc
cc      end subroutine find_curl_vxb

      end module precond_variables
