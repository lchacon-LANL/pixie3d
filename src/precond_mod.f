c  module matvec
c ###################################################################
      module matvec

        use grid

        use mg_internal

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

        integer(4) :: icomp   !Passed to setMGBC to define appropriate BC operation

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
        endif

        do igrid=1,grid_params%ngrid
          if (.not.associated(mgarray%grid(igrid)%array)) then
            nxp = grid_params%nxv(igrid)+1
            nyp = grid_params%nyv(igrid)+1
            nzp = grid_params%nzv(igrid)+1
            allocate(mgarray%grid(igrid)%array(0:nxp,0:nyp,0:nzp,neq))
            mgarray%grid(igrid)%array = 0d0
          endif
        enddo

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

      real(8),allocatable,dimension(:) :: vecf,vecc
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,grid_params,fpointers)

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

      call mapMGVectorToArray(0,neq,vecc,nxc,nyc,nzc,arrayc,igridc
     .                       ,.false.)

      if (icmp /= 0) then
        icomp=icmp              !Define icomp for BCs
        call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcnd)
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

      call allocPointers(neq,grid_params,fpointers)

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

        real(8), allocatable, dimension(:,:,:,:,:) :: nabla_v0

        real(8), allocatable, dimension(:,:,:,:)   :: p0

        real(8), allocatable, dimension(:,:,:)  :: mgnablaV0

        real(8), allocatable, dimension(:,:)    :: rho_diag,tmp_diag
     .                                            ,b_diag,v_diag

        real(8), allocatable, dimension(:,:)    :: mgj0cnv,mgadvdiffV0

        real(8), allocatable, dimension(:)      :: mgdivV0

        integer(4), allocatable, dimension(:,:) :: bcs

        type (var_array),target :: varray

        type (mg_array ),target :: gp0,gb0,gv0,grho0

      contains

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

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq
        integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

        real(8)    :: dhx,dhy,dhz
        real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

        real(8)    :: xim,yim,zim,xip,yip,zip
     .               ,xjm,yjm,zjm,xjp,yjp,zjp
     .               ,xkm,ykm,zkm,xkp,ykp,zkp
     .               ,x0,y0,z0

        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .               ,jacp,jacm,jach,jac0

        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm

        real(8)    :: bxx,byy,bzz
     .               ,bxip,bxim,bxjp,bxjm,bxkp,bxkm
     .               ,byip,byim,byjp,byjm,bykp,bykm
     .               ,bzip,bzim,bzjp,bzjm,bzkp,bzkm

        real(8),allocatable,dimension(:,:,:,:) :: b0_cnv

        logical    :: sing_point

c     Begin program

        sing_point = .false.
        if (i == 1 .and. bcond(1) == SP) sing_point = .true.

c     Defaults

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        jac = gmetric%grid(igx)%jac(i,j,k)

        dhx = 2.*dxh(ig)
        dhy = 2.*dyh(jg)
        dhz = 2.*dzh(kg)

c     Exceptions

cc        if (hex == 1) then
cc          im = i
cc          dhx = dx(ig)
cc        endif
cc
cc        if (hey == 1) then
cc          jm = j
cc          dhy = dy(jg)
cc        endif
cc
cc        if (hez == 1) then
cc          km = k
cc          dhz = dz(kg)
cc        endif

        select case(half_elem)
        case (1)
          dhx = dx(ig)
          im = i

          vxx = 0.5*(vv(i,j,k,1)+vv(ip,j,k,1))
          vyy = 0.5*(vv(i,j,k,2)+vv(ip,j,k,2))
          vzz = 0.5*(vv(i,j,k,3)+vv(ip,j,k,3))

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

          bxx = 0.5*(bb(i,j,k,1)+bb(ip,j,k,1))
          byy = 0.5*(bb(i,j,k,2)+bb(ip,j,k,2))
          bzz = 0.5*(bb(i,j,k,3)+bb(ip,j,k,3))

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
        case (2)
          dhy = dy(jg)
          jm = j

          vxx = 0.5*(vv(i,j,k,1)+vv(i,jp,k,1))
          vyy = 0.5*(vv(i,j,k,2)+vv(i,jp,k,2))
          vzz = 0.5*(vv(i,j,k,3)+vv(i,jp,k,3))

          if (sing_point) then
            vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1))*0.5
     .            +(vv(i ,j,k,1)+vv(i ,jp,k,1))*0.5
            vxim = (vv(im,j,k,1)+vv(im,jp,k,1))
            vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2))*0.5
     .            +(vv(i ,j,k,2)+vv(i ,jp,k,2))*0.5
            vyim = (vv(im,j,k,2)+vv(im,jp,k,2))
            vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3))*0.5
     .            +(vv(i ,j,k,3)+vv(i ,jp,k,3))*0.5
            vzim = (vv(im,j,k,3)+vv(im,jp,k,3))
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

          bxx = 0.5*(bb(i,j,k,1)+bb(i,jp,k,1))
          byy = 0.5*(bb(i,j,k,2)+bb(i,jp,k,2))
          bzz = 0.5*(bb(i,j,k,3)+bb(i,jp,k,3))

          if (sing_point) then
            bxip = (bb(ip,j,k,1)+bb(ip,jp,k,1))*0.5
     .            +(bb(i ,j,k,1)+bb(i ,jp,k,1))*0.5
            bxim = (bb(im,j,k,1)+bb(im,jp,k,1))
            byip = (bb(ip,j,k,2)+bb(ip,jp,k,2))*0.5
     .            +(bb(i ,j,k,2)+bb(i ,jp,k,2))*0.5
            byim = (bb(im,j,k,2)+bb(im,jp,k,2))
            bzip = (bb(ip,j,k,3)+bb(ip,jp,k,3))*0.5
     .            +(bb(i ,j,k,3)+bb(i ,jp,k,3))*0.5
            bzim = (bb(im,j,k,3)+bb(im,jp,k,3))
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
          dhz = dz(kg)
          km = k

          vxx = 0.5*(vv(i,j,k,1)+vv(i,j,kp,1))
          vyy = 0.5*(vv(i,j,k,2)+vv(i,j,kp,2))
          vzz = 0.5*(vv(i,j,k,3)+vv(i,j,kp,3))

          if (sing_point) then
            vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1))*0.5
     .            +(vv(i ,j,k,1)+vv(i ,j,kp,1))*0.5
            vxim = (vv(im,j,k,1)+vv(im,j,kp,1))
            vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2))*0.5
     .            +(vv(i ,j,k,2)+vv(i ,j,kp,2))*0.5
            vyim = (vv(im,j,k,2)+vv(im,j,kp,2))
            vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3))*0.5
     .            +(vv(i ,j,k,3)+vv(i ,j,kp,3))*0.5
            vzim = (vv(im,j,k,3)+vv(im,j,kp,3))
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

          bxx = 0.5*(bb(i,j,k,1)+bb(i,j,kp,1))
          byy = 0.5*(bb(i,j,k,2)+bb(i,j,kp,2))
          bzz = 0.5*(bb(i,j,k,3)+bb(i,j,kp,3))

          if (sing_point) then
            bxip = (bb(ip,j,k,1)+bb(ip,j,kp,1))*0.5
     .            +(bb(i ,j,k,1)+bb(i ,j,kp,1))*0.5
            bxim = (bb(im,j,k,1)+bb(im,j,kp,1))
            byip = (bb(ip,j,k,2)+bb(ip,j,kp,2))*0.5
     .            +(bb(i ,j,k,2)+bb(i ,j,kp,2))*0.5
            byim = (bb(im,j,k,2)+bb(im,j,kp,2))
            bzip = (bb(ip,j,k,3)+bb(ip,j,kp,3))*0.5
     .            +(bb(i ,j,k,3)+bb(i ,j,kp,3))*0.5
            bzim = (bb(im,j,k,3)+bb(im,j,kp,3))
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

          !Velocity
          vxx = vv(i,j,k,1)
          vyy = vv(i,j,k,2)
          vzz = vv(i,j,k,3)

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
          bxx = bb(i,j,k,1)
          byy = bb(i,j,k,2)
          bzz = bb(i,j,k,3)

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

c     Grid quantities

cc        jacip  = jacobian(xip,yip,zip,cartsn)
cc        jacim  = jacobian(xim,yim,zim,cartsn)
cc        jacjp  = jacobian(xjp,yjp,zjp,cartsn)
cc        jacjm  = jacobian(xjm,yjm,zjm,cartsn)
cc        jackp  = jacobian(xkp,ykp,zkp,cartsn)
cc        jackm  = jacobian(xkm,ykm,zkm,cartsn)

        jacip  = gmetric%grid(igx)%jac(ip,j,k)
        jacim  = gmetric%grid(igx)%jac(im,j,k)
        jacjp  = gmetric%grid(igx)%jac(i,jp,k)
        jacjm  = gmetric%grid(igx)%jac(i,jm,k)
        jackp  = gmetric%grid(igx)%jac(i,j,kp)
        jackm  = gmetric%grid(igx)%jac(i,j,km)

c     Components

        !component 1

        flxjp = ( vyjp*bxjp-vxjp*byjp )/(jac+jacjp)*2
        flxjm = ( vyjm*bxjm-vxjm*byjm )/(jac+jacjm)*2

        flxkp = ( vzkp*bxkp-vxkp*bzkp )/(jac+jackp)*2
        flxkm = ( vzkm*bxkm-vxkm*bzkm )/(jac+jackm)*2

        a1 =  (flxjp-flxjm)/dhy
     .       +(flxkp-flxkm)/dhz

        !component 2

        flxip = ( vxip*byip-vyip*bxip )/(jac+jacip)*2
        flxim = ( vxim*byim-vyim*bxim )/(jac+jacim)*2

        flxkp = ( vzkp*bykp-vykp*bzkp )/(jac+jackp)*2
        flxkm = ( vzkm*bykm-vykm*bzkm )/(jac+jackm)*2

        a2 =  (flxip-flxim)/dhx
     .       +(flxkp-flxkm)/dhz

        !component 3

        flxip = ( vxip*bzip-vzip*bxip )/(jac+jacip)*2
        flxim = ( vxim*bzim-vzim*bxim )/(jac+jacim)*2

        flxjp = ( vyjp*bzjp-vzjp*byjp )/(jac+jacjp)*2
        flxjm = ( vyjm*bzjm-vzjm*byjm )/(jac+jacjm)*2

        a3 =  (flxip-flxim)/dhx
     .       +(flxjp-flxjm)/dhy

      end subroutine find_curl_vxb

      end module precond_variables
