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

c module precond_variables
c ######################################################################
      module precond_variables

        use timeStepping

        use precond_setup

        use transport_params

        use equilibrium

cc        use parameters

        use constants

        use auxiliaryVariables

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

        real(8), allocatable, dimension(:,:,:)  :: mgnablaV0

        real(8), allocatable, dimension(:)      :: mgdivV0

        integer(4), allocatable, dimension(:,:) :: bcs

        type (var_array),target :: varray

        type (mg_array ),target :: gp0,gb0,gv0,grho0,gb0_cov

        real(8),pointer,dimension(:,:,:):: rho,rvx,rvy,rvz,bx,by,bz,tmp

        logical :: form_diag=.true.,vol_wgt=.true.,gm_smooth=.false.
     .            ,si_car=.true.,line_relax=.false.

        INTERFACE ASSIGNMENT (=)
          module procedure scatter,gather
        END INTERFACE

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
        allocate (mgnablaV0  (ntotd2p,3,3))
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
        deallocate (mgnablaV0)
        deallocate (bcs)

        call deallocateMGArray(gp0)
        call deallocateMGArray(grho0)
        call deallocateMGArray(gb0)
        call deallocateMGArray(gb0_cov)
        call deallocateMGArray(gv0)

        call deallocateDerivedType(varray)

c     End program

      end subroutine deallocPrecVariables

c     gather
c     ###################################################################
      subroutine gather(xout,xin)

c     -------------------------------------------------------------------
c     Deallocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        real(8),intent(IN)  :: xin(:,:)
        real(8),intent(OUT) :: xout(:)

c     Local variables

        integer(4) :: i,j,k,ii,iii,ieq
        real(8)    :: dvol

c     Begin program

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              do ieq=1,neqd
                iii = ieq + neqd*(ii-1)
                xout(iii) = xin(ii,ieq)
              enddo
            enddo
          enddo
        enddo

c     End program

      end subroutine gather

c     scatter
c     ###################################################################
      subroutine scatter(xout,xin)

c     -------------------------------------------------------------------
c     Deallocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        real(8),intent(IN)  :: xin(:)
        real(8),intent(OUT) :: xout(:,:)

c     Local variables

        integer(4) :: i,j,k,ii,iii,ieq
        real(8)    :: dvol

c     Begin program

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              if (.not.vol_wgt) then  !Residuals are volume-weighed by default
                dvol = 1d0/volume(i,j,k,igx,igy,igz)
              else
                dvol = 1d0
              endif

              do ieq=1,neqd
                iii = ieq + neqd*(ii-1)
                xout(ii,ieq) = xin(iii)*dvol
              enddo

            enddo
          enddo
        enddo

c     End program

      end subroutine scatter

c     findCoeffs
c     ###################################################################
      subroutine findCoeffs

c     -------------------------------------------------------------------
c     Finds coefficients for linearized systems in preconditioner
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer(4) :: order,nxx,nyy,nzz,igrid,ii,ivar,igr,i,j
        real(8)    :: dvol,veclap(3)

        real(8), allocatable, dimension(:,:,:,:) :: vector,vel
        real(8), allocatable, dimension(:,:,:,:,:) :: tensor

c     Begin program

        nxx = grid_params%nxv(igx)
        nyy = grid_params%nyv(igy)
        nzz = grid_params%nzv(igz)

        nx  = nxx
        ny  = nyy
        nz  = nzz

        igrid = igx

        order = 0

c     Store density in all grids (w/o BCs)

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

        call restrictMGArray(IBX,3,gb0_cov,bcs(:,IBX:IBZ),igrid,order
     .                      ,iscnv=.false.)

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

        !v0/dt+theta(v0.grad(v0)-mu*veclap(v0)) and nabla_v0
        allocate(vector(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .          ,vel   (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .          ,tensor(0:nxx+1,0:nyy+1,0:nzz+1,3,3))

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
              veclap = nuu(i,j,k)
     $                *veclaplacian(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                             ,vel,alt_eom,vol=.false.)

              do ivar=1,3
                vector(i,j,k,ivar) = 
     .                            gv0%grid(igrid)%array(i,j,k,ivar)/dt
     .                          + alpha*vx(i,j,k)*nabla_v(1,ivar)/jac
     .                          + alpha*vy(i,j,k)*nabla_v(2,ivar)/jac
     .                          + alpha*vz(i,j,k)*nabla_v(3,ivar)/jac
     .                          - alpha*veclap(ivar)
              enddo

              tensor(i,j,k,:,:) = nabla_v
            enddo
          enddo
        enddo

        do ivar=1,3
          call restrictArrayToMGVector(1,nxx,nyy,nzz
     .                              ,vector(:,:,:,ivar)
     .                              ,mgadvdiffV0(:,ivar)
     .                              ,igrid,order,.false.)
        enddo

        do i=1,3
          do j=1,3
            call restrictArrayToMGVector(1,nxx,nyy,nzz
     .                              ,tensor(:,:,:,i,j)
     .                              ,mgnablaV0(:,i,j)
     .                              ,igrid,order,.false.)
          enddo
        enddo

        deallocate(vector,vel,tensor)

c     End program

      end subroutine findCoeffs

c     find_curl_vxb
c     ###################################################################
      subroutine find_curl_vxb(i,j,k,nx,ny,nz,vv,bb,a1,a2,a3
     .                                  ,half_elem,igrid)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of -curl(vv x bb) at the
c     grid node (i,j,k). One sided derivatives are employed when half_elem=1
c     (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nx,ny,nz,half_elem,igrid
        real(8)    :: a1,a2,a3
        real(8),target :: vv(0:nx+1,0:ny+1,0:nz+1,3)
     $                   ,bb(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq,igx,igy,igz
        integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

        real(8)    :: idhx,idhy,idhz,a(3)
        real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .               ,jacp,jacm,jacph,jacmh,jach,jac0

        real(8)    :: vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm

        real(8)    :: bxip,bxim,bxjp,bxjm,bxkp,bxkm
     .               ,byip,byim,byjp,byjm,bykp,bykm
     .               ,bzip,bzim,bzjp,bzjm,bzkp,bzkm

        logical    :: sing_point,sing_coord

c     Begin program

        igx = igrid
        igy = igrid
        igz = igrid

        sing_point = isSP(i,j,k,igx,igy,igz)

        sing_coord =(i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx))
     .              .and. (bcond(1) == SP)
        sing_coord = .false.

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

          if (isSP(ip,j,k,igx,igy,igz)) then !i=0 and a SP at ip=1

            jacip  = gmetric%grid(igx)%jac(ip,j,k)
            jacim  = gmetric%grid(igx)%jac(i ,j,k)
            jacjp  = gmetric%grid(igx)%jac(i,jp,k)
            jacjm  = gmetric%grid(igx)%jac(i,jm,k)
            jackp  = gmetric%grid(igx)%jac(i,j,kp)
            jackm  = gmetric%grid(igx)%jac(i,j,km)

            vxip = vv(ip,j,k,1)
            vxim = vv(i ,j,k,1)
            vyip = vv(ip,j,k,2)
            vyim = vv(i ,j,k,2)
            vzip = vv(ip,j,k,3)
            vzim = vv(i ,j,k,3)

            vxjp = vv(i,jp,k,1)
            vxjm = vv(i,jm,k,1)
            vyjp = vv(i,jp,k,2)
            vyjm = vv(i,jm,k,2)
            vzkp = vv(i,j,kp,3)
            vzkm = vv(i,j,km,3)

            vxkp = vv(i,j,kp,1)
            vxkm = vv(i,j,km,1)
            vykp = vv(i,j,kp,2)
            vykm = vv(i,j,km,2)
            vzkp = vv(i,j,kp,3)
            vzkm = vv(i,j,km,3)

            bxip = bb(ip,j,k,1)
            bxim = bb(i ,j,k,1)
            byip = bb(ip,j,k,2)
            byim = bb(i ,j,k,2)
            bzip = bb(ip,j,k,3)
            bzim = bb(i ,j,k,3)

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

cc          elseif (sing_point) then !SP at i=1
          elseif (sing_coord) then

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

c$$$            vxjp = 0.5*(vv(ip,jp,k,1)+vv(i,jp,k,1))
c$$$            vxjm = 0.5*(vv(ip,jm,k,1)+vv(i,jm,k,1))
c$$$            vyjp = 0.5*(vv(ip,jp,k,2)+vv(i,jp,k,2))
c$$$            vyjm = 0.5*(vv(ip,jm,k,2)+vv(i,jm,k,2))
c$$$            vzjp = 0.5*(vv(ip,jp,k,3)+vv(i,jp,k,3))
c$$$            vzjm = 0.5*(vv(ip,jm,k,3)+vv(i,jm,k,3))

            vxjp = 0.5*(vv(ip,jp,k,1)/gmetric%grid(igx)%jac(ip,jp,k)
     .                 +vv(i ,jp,k,1)/gmetric%grid(igx)%jac(i ,jp,k))
     .                *jacjp
            vxjm = 0.5*(vv(ip,jm,k,1)/gmetric%grid(igx)%jac(ip,jm,k)
     .                 +vv(i ,jm,k,1)/gmetric%grid(igx)%jac(i ,jm,k))
     .                *jacjm
            vyjp = 0.5*(vv(ip,jp,k,2)+vv(i,jp,k,2))
            vyjm = 0.5*(vv(ip,jm,k,2)+vv(i,jm,k,2))
            vzjp = 0.5*(vv(ip,jp,k,3)/gmetric%grid(igx)%jac(ip,jp,k)
     .                 +vv(i ,jp,k,3)/gmetric%grid(igx)%jac(i ,jp,k))
     .                *jacjp
            vzjm = 0.5*(vv(ip,jm,k,3)/gmetric%grid(igx)%jac(ip,jm,k)
     .                 +vv(i ,jm,k,3)/gmetric%grid(igx)%jac(i ,jm,k))
     .                *jacjm

c$$$            vxkp = 0.5*(vv(ip,j,kp,1)+vv(i,j,kp,1))
c$$$            vxkm = 0.5*(vv(ip,j,km,1)+vv(i,j,km,1))
c$$$            vykp = 0.5*(vv(ip,j,kp,2)+vv(i,j,kp,2))
c$$$            vykm = 0.5*(vv(ip,j,km,2)+vv(i,j,km,2))
c$$$            vzkp = 0.5*(vv(ip,j,kp,3)+vv(i,j,kp,3))
c$$$            vzkm = 0.5*(vv(ip,j,km,3)+vv(i,j,km,3))

            vxkp = 0.5*(vv(ip,j,kp,1)/gmetric%grid(igx)%jac(ip,j,kp)
     .                 +vv(i ,j,kp,1)/gmetric%grid(igx)%jac(i ,j,kp))
     .                *jackp                                       
            vxkm = 0.5*(vv(ip,j,km,1)/gmetric%grid(igx)%jac(ip,j,km)
     .                 +vv(i ,j,km,1)/gmetric%grid(igx)%jac(i ,j,km))
     .                *jackm                                       
            vykp = 0.5*(vv(ip,j,kp,2)+vv(i,j,kp,2))                
            vykm = 0.5*(vv(ip,j,km,2)+vv(i,j,km,2))                
            vzkp = 0.5*(vv(ip,j,kp,3)/gmetric%grid(igx)%jac(ip,j,kp)
     .                 +vv(i ,j,kp,3)/gmetric%grid(igx)%jac(i ,j,kp))
     .                *jackp                                       
            vzkm = 0.5*(vv(ip,j,km,3)/gmetric%grid(igx)%jac(ip,j,km)
     .                 +vv(i ,j,km,3)/gmetric%grid(igx)%jac(i ,j,km))
     .                *jackm

            bxip = bb(ip,j,k,1)
            bxim = bb(i ,j,k,1)
            byip = bb(ip,j,k,2)
            byim = bb(i ,j,k,2)
            bzip = bb(ip,j,k,3)
            bzim = bb(i ,j,k,3)

c$$$            bxjp = 0.5*(bb(ip,jp,k,1)+bb(i,jp,k,1))
c$$$            bxjm = 0.5*(bb(ip,jm,k,1)+bb(i,jm,k,1))
c$$$            byjp = 0.5*(bb(ip,jp,k,2)+bb(i,jp,k,2))
c$$$            byjm = 0.5*(bb(ip,jm,k,2)+bb(i,jm,k,2))
c$$$            bzjp = 0.5*(bb(ip,jp,k,3)+bb(i,jp,k,3))
c$$$            bzjm = 0.5*(bb(ip,jm,k,3)+bb(i,jm,k,3))

            bxjp = 0.5*(bb(ip,jp,k,1)/gmetric%grid(igx)%jac(ip,jp,k)
     .                 +bb(i ,jp,k,1)/gmetric%grid(igx)%jac(i ,jp,k))
     .                *jacjp
            bxjm = 0.5*(bb(ip,jm,k,1)/gmetric%grid(igx)%jac(ip,jm,k)
     .                 +bb(i ,jm,k,1)/gmetric%grid(igx)%jac(i ,jm,k))
     .                *jacjm
            byjp = 0.5*(bb(ip,jp,k,2)+bb(i,jp,k,2))
            byjm = 0.5*(bb(ip,jm,k,2)+bb(i,jm,k,2))
            bzjp = 0.5*(bb(ip,jp,k,3)/gmetric%grid(igx)%jac(ip,jp,k)
     .                 +bb(i ,jp,k,3)/gmetric%grid(igx)%jac(i ,jp,k))
     .                *jacjp
            bzjm = 0.5*(bb(ip,jm,k,3)/gmetric%grid(igx)%jac(ip,jm,k)
     .                 +bb(i ,jm,k,3)/gmetric%grid(igx)%jac(i ,jm,k))
     .                *jacjm

c$$$            bxkp = 0.5*(bb(ip,j,kp,1)+bb(i,j,kp,1))
c$$$            bxkm = 0.5*(bb(ip,j,km,1)+bb(i,j,km,1))
c$$$            bykp = 0.5*(bb(ip,j,kp,2)+bb(i,j,kp,2))
c$$$            bykm = 0.5*(bb(ip,j,km,2)+bb(i,j,km,2))
c$$$            bzkp = 0.5*(bb(ip,j,kp,3)+bb(i,j,kp,3))
c$$$            bzkm = 0.5*(bb(ip,j,km,3)+bb(i,j,km,3))

            bxkp = 0.5*(bb(ip,j,kp,1)/gmetric%grid(igx)%jac(ip,j,kp)
     .                 +bb(i ,j,kp,1)/gmetric%grid(igx)%jac(i ,j,kp))
     .                *jackp                                       
            bxkm = 0.5*(bb(ip,j,km,1)/gmetric%grid(igx)%jac(ip,j,km)
     .                 +bb(i ,j,km,1)/gmetric%grid(igx)%jac(i ,j,km))
     .                *jackm                                       
            bykp = 0.5*(bb(ip,j,kp,2)+bb(i,j,kp,2))                
            bykm = 0.5*(bb(ip,j,km,2)+bb(i,j,km,2))                
            bzkp = 0.5*(bb(ip,j,kp,3)/gmetric%grid(igx)%jac(ip,j,kp)
     .                 +bb(i ,j,kp,3)/gmetric%grid(igx)%jac(i ,j,kp))
     .                *jackp                                       
            bzkm = 0.5*(bb(ip,j,km,3)/gmetric%grid(igx)%jac(ip,j,km)
     .                 +bb(i ,j,km,3)/gmetric%grid(igx)%jac(i ,j,km))
     .                *jackm

          else

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

            vxjp = 0.5*(vv(ip,jp,k,1)+vv(i,jp,k,1))
            vxjm = 0.5*(vv(ip,jm,k,1)+vv(i,jm,k,1))
            vyjp = 0.5*(vv(ip,jp,k,2)+vv(i,jp,k,2))
            vyjm = 0.5*(vv(ip,jm,k,2)+vv(i,jm,k,2))
            vzjp = 0.5*(vv(ip,jp,k,3)+vv(i,jp,k,3))
            vzjm = 0.5*(vv(ip,jm,k,3)+vv(i,jm,k,3))

            vxkp = 0.5*(vv(ip,j,kp,1)+vv(i,j,kp,1))
            vxkm = 0.5*(vv(ip,j,km,1)+vv(i,j,km,1))
            vykp = 0.5*(vv(ip,j,kp,2)+vv(i,j,kp,2))
            vykm = 0.5*(vv(ip,j,km,2)+vv(i,j,km,2))
            vzkp = 0.5*(vv(ip,j,kp,3)+vv(i,j,kp,3))
            vzkm = 0.5*(vv(ip,j,km,3)+vv(i,j,km,3))

            bxip = bb(ip,j,k,1)
            bxim = bb(i ,j,k,1)
            byip = bb(ip,j,k,2)
            byim = bb(i ,j,k,2)
            bzip = bb(ip,j,k,3)
            bzim = bb(i ,j,k,3)

            bxjp = 0.5*(bb(ip,jp,k,1)+bb(i,jp,k,1))
            bxjm = 0.5*(bb(ip,jm,k,1)+bb(i,jm,k,1))
            byjp = 0.5*(bb(ip,jp,k,2)+bb(i,jp,k,2))
            byjm = 0.5*(bb(ip,jm,k,2)+bb(i,jm,k,2))
            bzjp = 0.5*(bb(ip,jp,k,3)+bb(i,jp,k,3))
            bzjm = 0.5*(bb(ip,jm,k,3)+bb(i,jm,k,3))

            bxkp = 0.5*(bb(ip,j,kp,1)+bb(i,j,kp,1))
            bxkm = 0.5*(bb(ip,j,km,1)+bb(i,j,km,1))
            bykp = 0.5*(bb(ip,j,kp,2)+bb(i,j,kp,2))
            bykm = 0.5*(bb(ip,j,km,2)+bb(i,j,km,2))
            bzkp = 0.5*(bb(ip,j,kp,3)+bb(i,j,kp,3))
            bzkm = 0.5*(bb(ip,j,km,3)+bb(i,j,km,3))
          endif

        case (2)

          idhy = 1./dy(jg)
          jm = j

          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
     .                 +gmetric%grid(igx)%jac(ip,j ,k))
          jacim  = 0.5*(gmetric%grid(igx)%jac(im,jp,k)
     .                 +gmetric%grid(igx)%jac(im,j ,k))
          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
          jacjm  = gmetric%grid(igx)%jac(i,j ,k)
          jackp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
     .                 +gmetric%grid(igx)%jac(i,j ,kp))
          jackm  = 0.5*(gmetric%grid(igx)%jac(i,jp,km)
     .                 +gmetric%grid(igx)%jac(i,j ,km))

          if (sing_point) then
            idhx = 1./dxh(ig)

cc            jacip  =0.25*(gmetric%grid(igx)%jac(ip,jp,k)
cc     .                   +gmetric%grid(igx)%jac(i ,jp,k)
cc     .                   +gmetric%grid(igx)%jac(ip,j ,k)
cc     .                   +gmetric%grid(igx)%jac(i ,j ,k))
cc            jacim  = 0.5*(gmetric%grid(igx)%jac(im,jp,k)
cc     .                   +gmetric%grid(igx)%jac(im,j ,k))
cc
cc            vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1)
cc     .             +vv(i ,j,k,1)+vv(i ,jp,k,1))*0.25
cc            vxim = (vv(im,j,k,1)+vv(im,jp,k,1))*0.5
cc            vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2)
cc     .             +vv(i ,j,k,2)+vv(i ,jp,k,2))*0.25
cc            vyim = (vv(im,j,k,2)+vv(im,jp,k,2))*0.5
cc            vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3)
cc     .             +vv(i ,j,k,3)+vv(i ,jp,k,3))*0.25
cc            vzim = (vv(im,j,k,3)+vv(im,jp,k,3))*0.5

            jacip = 0.25*(gmetric%grid(igx)%jac(ip,jp,k)
     .                  +gmetric%grid(igx)%jac(i ,jp,k)
     .                  +gmetric%grid(igx)%jac(ip,j ,k)
     .                  +gmetric%grid(igx)%jac(i ,j ,k))
            jacim =  0.5*(gmetric%grid(igx)%jac(im,jp,k)
     .                   +gmetric%grid(igx)%jac(im,j ,k))

            vxip = 0.25*(vv(ip,j ,k,1)/gmetric%grid(igx)%jac(ip,j ,k)
     .                  +vv(ip,jp,k,1)/gmetric%grid(igx)%jac(ip,jp,k)
     .                  +vv(i ,j ,k,1)/gmetric%grid(igx)%jac(i ,j ,k)
     .                  +vv(i ,jp,k,1)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacip
            vxim = 0.5  *(vv(im,j,k,1)+vv(im,jp,k,1))
            vyip = 0.25 *(vv(ip,j,k,2)+vv(ip,jp,k,2)
     .                   +vv(i ,j,k,2)+vv(i ,jp,k,2))
            vyim = 0.5  *(vv(im,j,k,2)+vv(im,jp,k,2))
            vzip = 0.25*(vv(ip,j ,k,3)/gmetric%grid(igx)%jac(ip,j ,k) 
     .                  +vv(ip,jp,k,3)/gmetric%grid(igx)%jac(ip,jp,k) 
     .                  +vv(i ,j ,k,3)/gmetric%grid(igx)%jac(i ,j ,k) 
     .                  +vv(i ,jp,k,3)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacip
            vzim = 0.5  *(vv(im,j,k,3)+vv(im,jp,k,3))
          elseif (sing_coord) then
            idhx = 1./dxh(ig)

            jacip = 0.25*(gmetric%grid(igx)%jac(ip,jp,k)
     .                  +gmetric%grid(igx)%jac(i ,jp,k)
     .                  +gmetric%grid(igx)%jac(ip,j ,k)
     .                  +gmetric%grid(igx)%jac(i ,j ,k))
            jacim = 0.25*(gmetric%grid(igx)%jac(im,jp,k)
     .                  +gmetric%grid(igx)%jac(i ,jp,k)
     .                  +gmetric%grid(igx)%jac(im,j ,k)
     .                  +gmetric%grid(igx)%jac(i ,j ,k))

            vxip = 0.25*(vv(ip,j ,k,1)/gmetric%grid(igx)%jac(ip,j ,k)
     .                  +vv(ip,jp,k,1)/gmetric%grid(igx)%jac(ip,jp,k)
     .                  +vv(i ,j ,k,1)/gmetric%grid(igx)%jac(i ,j ,k)
     .                  +vv(i ,jp,k,1)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacip
            vxim = 0.25*(vv(im,j ,k,1)/gmetric%grid(igx)%jac(im,j ,k)
     .                  +vv(im,jp,k,1)/gmetric%grid(igx)%jac(im,jp,k)
     .                  +vv(i ,j ,k,1)/gmetric%grid(igx)%jac(i ,j ,k)
     .                  +vv(i ,jp,k,1)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacim

            vyip = 0.25 *(vv(ip,j,k,2)+vv(ip,jp,k,2)
     .                   +vv(i ,j,k,2)+vv(i ,jp,k,2))
            vyim = 0.25 *(vv(im,j,k,2)+vv(im,jp,k,2)
     .                   +vv(i ,j,k,2)+vv(i ,jp,k,2))

            vzip = 0.25*(vv(ip,j ,k,3)/gmetric%grid(igx)%jac(ip,j ,k) 
     .                  +vv(ip,jp,k,3)/gmetric%grid(igx)%jac(ip,jp,k) 
     .                  +vv(i ,j ,k,3)/gmetric%grid(igx)%jac(i ,j ,k) 
     .                  +vv(i ,jp,k,3)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacip
            vzim = 0.25*(vv(im,j ,k,3)/gmetric%grid(igx)%jac(im,j ,k) 
     .                  +vv(im,jp,k,3)/gmetric%grid(igx)%jac(im,jp,k) 
     .                  +vv(i ,j ,k,3)/gmetric%grid(igx)%jac(i ,j ,k) 
     .                  +vv(i ,jp,k,3)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacim
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

          if (sing_point) then
cc            bxip = (bb(ip,j,k,1)+bb(ip,jp,k,1)
cc     .             +bb(i ,j,k,1)+bb(i ,jp,k,1))*0.25
cc            bxim = (bb(im,j,k,1)+bb(im,jp,k,1))*0.5
cc            byip = (bb(ip,j,k,2)+bb(ip,jp,k,2)
cc     .             +bb(i ,j,k,2)+bb(i ,jp,k,2))*0.25
cc            byim = (bb(im,j,k,2)+bb(im,jp,k,2))*0.5
cc            bzip = (bb(ip,j,k,3)+bb(ip,jp,k,3)
cc     .             +bb(i ,j,k,3)+bb(i ,jp,k,3))*0.25
cc            bzim = (bb(im,j,k,3)+bb(im,jp,k,3))*0.5

            bxip = 0.25*(bb(ip,j ,k,1)/gmetric%grid(igx)%jac(ip,j ,k)
     .                  +bb(ip,jp,k,1)/gmetric%grid(igx)%jac(ip,jp,k)
     .                  +bb(i ,j ,k,1)/gmetric%grid(igx)%jac(i ,j ,k)
     .                  +bb(i ,jp,k,1)/gmetric%grid(igx)%jac(i ,jp,k))
     .                   *jacip
            bxim = 0.5  *(bb(im,j,k,1)+bb(im,jp,k,1))
            byip = 0.25 *(bb(ip,j,k,2)+bb(ip,jp,k,2)
     .                   +bb(i ,j,k,2)+bb(i ,jp,k,2))
            byim = 0.5  *(bb(im,j,k,2)+bb(im,jp,k,2))
            bzip = 0.25*(bb(ip,j ,k,3)/gmetric%grid(igx)%jac(ip,j ,k) 
     .                  +bb(ip,jp,k,3)/gmetric%grid(igx)%jac(ip,jp,k) 
     .                  +bb(i ,j ,k,3)/gmetric%grid(igx)%jac(i ,j ,k) 
     .                  +bb(i ,jp,k,3)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacip
            bzim = 0.5  *(bb(im,j,k,3)+bb(im,jp,k,3))

          elseif (sing_coord) then
            bxip = 0.25*(bb(ip,j ,k,1)/gmetric%grid(igx)%jac(ip,j ,k)
     .                  +bb(ip,jp,k,1)/gmetric%grid(igx)%jac(ip,jp,k)
     .                  +bb(i ,j ,k,1)/gmetric%grid(igx)%jac(i ,j ,k)
     .                  +bb(i ,jp,k,1)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacip
            bxim = 0.25*(bb(im,j ,k,1)/gmetric%grid(igx)%jac(im,j ,k)
     .                  +bb(im,jp,k,1)/gmetric%grid(igx)%jac(im,jp,k)
     .                  +bb(i ,j ,k,1)/gmetric%grid(igx)%jac(i ,j ,k)
     .                  +bb(i ,jp,k,1)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacim

            byip = 0.25 *(bb(ip,j,k,2)+bb(ip,jp,k,2)
     .                   +bb(i ,j,k,2)+bb(i ,jp,k,2))
            byim = 0.25 *(bb(im,j,k,2)+bb(im,jp,k,2)
     .                   +bb(i ,j,k,2)+bb(i ,jp,k,2))

            bzip = 0.25*(bb(ip,j ,k,3)/gmetric%grid(igx)%jac(ip,j ,k) 
     .                  +bb(ip,jp,k,3)/gmetric%grid(igx)%jac(ip,jp,k) 
     .                  +bb(i ,j ,k,3)/gmetric%grid(igx)%jac(i ,j ,k) 
     .                  +bb(i ,jp,k,3)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacip
            bzim = 0.25*(bb(im,j ,k,3)/gmetric%grid(igx)%jac(im,j ,k) 
     .                  +bb(im,jp,k,3)/gmetric%grid(igx)%jac(im,jp,k) 
     .                  +bb(i ,j ,k,3)/gmetric%grid(igx)%jac(i ,j ,k) 
     .                  +bb(i ,jp,k,3)/gmetric%grid(igx)%jac(i ,jp,k))
     .                  *jacim
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

          if (sing_point) then
            idhx = 1./dxh(ig)

cc            jacip  =0.25*(gmetric%grid(igx)%jac(ip,j,kp)
cc     .                   +gmetric%grid(igx)%jac(i ,j,kp)
cc     .                   +gmetric%grid(igx)%jac(ip,j,k )
cc     .                   +gmetric%grid(igx)%jac(i ,j,k ))
cc            jacim  = 0.5*(gmetric%grid(igx)%jac(im,j,kp)
cc     .                   +gmetric%grid(igx)%jac(im,j,k ))
cc
cc            vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1)
cc     .             +vv(i ,j,k,1)+vv(i ,j,kp,1))*0.25
cc            vxim = (vv(im,j,k,1)+vv(im,j,kp,1))*0.5
cc            vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2)
cc     .             +vv(i ,j,k,2)+vv(i ,j,kp,2))*0.25
cc            vyim = (vv(im,j,k,2)+vv(im,j,kp,2))*0.5
cc            vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3)
cc     .             +vv(i ,j,k,3)+vv(i ,j,kp,3))*0.25
cc            vzim = (vv(im,j,k,3)+vv(im,j,kp,3))*0.5

            jacip = 0.25*(gmetric%grid(igx)%jac(ip,j,kp)
     .                   +gmetric%grid(igx)%jac(i ,j,kp)
     .                   +gmetric%grid(igx)%jac(ip,j,k )
     .                   +gmetric%grid(igx)%jac(i ,j,k ))
            jacim = 0.5 *(gmetric%grid(igx)%jac(im,j,kp)
     .                   +gmetric%grid(igx)%jac(im,j,k ))

            vxip = 0.25*(vv(ip,j,k ,1)/gmetric%grid(igx)%jac(ip,j,k )
     .                  +vv(ip,j,kp,1)/gmetric%grid(igx)%jac(ip,j,kp)
     .                  +vv(i ,j,k ,1)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +vv(i ,j,kp,1)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacip
            vxim = 0.5*(vv(im,j,k,1)+vv(im,j,kp,1))
            vyip = 0.25*(vv(ip,j,k ,2)
     .                  +vv(ip,j,kp,2)
     .                  +vv(i ,j,k ,2)
     .                  +vv(i ,j,kp,2))
            vyim = 0.5*(vv(im,j,k,2)+vv(im,j,kp,2))
            vzip = 0.25*(vv(ip,j,k ,3)/gmetric%grid(igx)%jac(ip,j,k )
     .                  +vv(ip,j,kp,3)/gmetric%grid(igx)%jac(ip,j,kp)
     .                  +vv(i ,j,k ,3)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +vv(i ,j,kp,3)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacip
            vzim = 0.5*(vv(im,j,k,3)+vv(im,j,kp,3))

          elseif (sing_coord) then
            idhx = 1./dxh(ig)

            jacip = 0.25*(gmetric%grid(igx)%jac(ip,j,kp)
     .                   +gmetric%grid(igx)%jac(i ,j,kp)
     .                   +gmetric%grid(igx)%jac(ip,j,k )
     .                   +gmetric%grid(igx)%jac(i ,j,k ))

            vxip = 0.25*(vv(ip,j,k ,1)/gmetric%grid(igx)%jac(ip,j,k )
     .                  +vv(ip,j,kp,1)/gmetric%grid(igx)%jac(ip,j,kp)
     .                  +vv(i ,j,k ,1)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +vv(i ,j,kp,1)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacip

            vyip = 0.25*(vv(ip,j,k ,2)
     .                  +vv(ip,j,kp,2)
     .                  +vv(i ,j,k ,2)
     .                  +vv(i ,j,kp,2))
            vzip = 0.25*(vv(ip,j,k ,3)/gmetric%grid(igx)%jac(ip,j,k )
     .                  +vv(ip,j,kp,3)/gmetric%grid(igx)%jac(ip,j,kp)
     .                  +vv(i ,j,k ,3)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +vv(i ,j,kp,3)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacip

            jacim = 0.25*(gmetric%grid(igx)%jac(im,j,kp)
     .                   +gmetric%grid(igx)%jac(i ,j,kp)
     .                   +gmetric%grid(igx)%jac(im,j,k )
     .                   +gmetric%grid(igx)%jac(i ,j,k ))

            vxim = 0.25*(vv(im,j,k ,1)/gmetric%grid(igx)%jac(im,j,k )
     .                  +vv(im,j,kp,1)/gmetric%grid(igx)%jac(im,j,kp)
     .                  +vv(i ,j,k ,1)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +vv(i ,j,kp,1)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacim

            vyim = 0.25*(vv(im,j,k ,2)
     .                  +vv(im,j,kp,2)
     .                  +vv(i ,j,k ,2)
     .                  +vv(i ,j,kp,2))
            vzim = 0.25*(vv(im,j,k ,3)/gmetric%grid(igx)%jac(im,j,k )
     .                  +vv(im,j,kp,3)/gmetric%grid(igx)%jac(im,j,kp)
     .                  +vv(i ,j,k ,3)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +vv(i ,j,kp,3)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacim

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

          if (sing_point) then
cc            bxip = (bb(ip,j,k,1)+bb(ip,j,kp,1)
cc     .             +bb(i ,j,k,1)+bb(i ,j,kp,1))*0.25
cc            bxim = (bb(im,j,k,1)+bb(im,j,kp,1))*0.5
cc            byip = (bb(ip,j,k,2)+bb(ip,j,kp,2)
cc     .             +bb(i ,j,k,2)+bb(i ,j,kp,2))*0.25
cc            byim = (bb(im,j,k,2)+bb(im,j,kp,2))*0.5
cc            bzip = (bb(ip,j,k,3)+bb(ip,j,kp,3)
cc     .             +bb(i ,j,k,3)+bb(i ,j,kp,3))*0.25
cc            bzim = (bb(im,j,k,3)+bb(im,j,kp,3))*0.5

            bxip = 0.25*(bb(ip,j,k ,1)/gmetric%grid(igx)%jac(ip,j,k )
     .                  +bb(ip,j,kp,1)/gmetric%grid(igx)%jac(ip,j,kp)
     .                  +bb(i ,j,k ,1)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +bb(i ,j,kp,1)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacip
            bxim = 0.5*(bb(im,j,k,1)+bb(im,j,kp,1))
            byip = 0.25*(bb(ip,j,k ,2)
     .                  +bb(ip,j,kp,2)
     .                  +bb(i ,j,k ,2)
     .                  +bb(i ,j,kp,2))
            byim = 0.5*(bb(im,j,k,2)+bb(im,j,kp,2))
            bzip = 0.25*(bb(ip,j,k ,3)/gmetric%grid(igx)%jac(ip,j,k )
     .                  +bb(ip,j,kp,3)/gmetric%grid(igx)%jac(ip,j,kp)
     .                  +bb(i ,j,k ,3)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +bb(i ,j,kp,3)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacip
            bzim = 0.5*(bb(im,j,k,3)+bb(im,j,kp,3))

          elseif (sing_coord) then
            bxip = 0.25*(bb(ip,j,k ,1)/gmetric%grid(igx)%jac(ip,j,k )
     .                  +bb(ip,j,kp,1)/gmetric%grid(igx)%jac(ip,j,kp)
     .                  +bb(i ,j,k ,1)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +bb(i ,j,kp,1)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacip

            byip = 0.25*(bb(ip,j,k ,2)
     .                  +bb(ip,j,kp,2)
     .                  +bb(i ,j,k ,2)
     .                  +bb(i ,j,kp,2))

            bzip = 0.25*(bb(ip,j,k ,3)/gmetric%grid(igx)%jac(ip,j,k )
     .                  +bb(ip,j,kp,3)/gmetric%grid(igx)%jac(ip,j,kp)
     .                  +bb(i ,j,k ,3)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +bb(i ,j,kp,3)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacip

            bxim = 0.25*(bb(im,j,k ,1)/gmetric%grid(igx)%jac(im,j,k )
     .                  +bb(im,j,kp,1)/gmetric%grid(igx)%jac(im,j,kp)
     .                  +bb(i ,j,k ,1)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +bb(i ,j,kp,1)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacim

            byim = 0.25*(bb(im,j,k ,2)
     .                  +bb(im,j,kp,2)
     .                  +bb(i ,j,k ,2)
     .                  +bb(i ,j,kp,2))

            bzim = 0.25*(bb(im,j,k ,3)/gmetric%grid(igx)%jac(im,j,k )
     .                  +bb(im,j,kp,3)/gmetric%grid(igx)%jac(im,j,kp)
     .                  +bb(i ,j,k ,3)/gmetric%grid(igx)%jac(i ,j,k )
     .                  +bb(i ,j,kp,3)/gmetric%grid(igx)%jac(i ,j,kp))
     .                  *jacim

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
          
          vec1 => vv
          vec2 => bb

          a = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .                  ,btensor_x,btensor_y,btensor_z,vol=.false.)

          nullify(vec1,vec2)

          a1 = a(1)
          a2 = a(2)
          a3 = a(3)

          return

          jacip  = gmetric%grid(igx)%jac(ip,j,k)
          jacim  = gmetric%grid(igx)%jac(im,j,k)
          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
          jacjm  = gmetric%grid(igx)%jac(i,jm,k)
          jackp  = gmetric%grid(igx)%jac(i,j,kp)
          jackm  = gmetric%grid(igx)%jac(i,j,km)

          !Velocity
          if (sing_point) then
            idhx = 1./dxh(ig)
            jach = 0.5*(jacip+jac)

            vxip = 0.5*(vv(ip,j,k,1)/jacip+vv(i,j,k,1)/jac)*jach
            vxim = vv(im,j,k,1)
            vyip = 0.5*(vv(ip,j,k,2)      +vv(i,j,k,2))
            vyim = vv(im,j,k,2)
            vzip = 0.5*(vv(ip,j,k,3)/jacip+vv(i,j,k,3)/jac)*jach
            vzim = vv(im,j,k,3)
          elseif (sing_coord) then
            idhx = 1./dxh(ig)
            jacph = 0.5*(jacip+jac)
            jacmh = 0.5*(jacim+jac)

            vxip = 0.5*(vv(ip,j,k,1)/jacip+vv(i,j,k,1)/jac)*jacph
            vxim = 0.5*(vv(im,j,k,1)/jacim+vv(i,j,k,1)/jac)*jacmh
            vyip = 0.5*(vv(ip,j,k,2)      +vv(i,j,k,2)    )
            vyim = 0.5*(vv(im,j,k,2)      +vv(i,j,k,2)    )
            vzip = 0.5*(vv(ip,j,k,3)/jacip+vv(i,j,k,3)/jac)*jacph
            vzim = 0.5*(vv(im,j,k,3)/jacim+vv(i,j,k,3)/jac)*jacmh
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
            bxip = 0.5*(bb(ip,j,k,1)/jacip+bb(i,j,k,1)/jac)*jach
            bxim = bb(im,j,k,1)
            byip = 0.5*(bb(ip,j,k,2)      +bb(i,j,k,2))
            byim = bb(im,j,k,2)
            bzip = 0.5*(bb(ip,j,k,3)/jacip+bb(i,j,k,3)/jac)*jach
            bzim = bb(im,j,k,3)

            jacip = jach
          elseif (sing_coord) then
            bxip = 0.5*(bb(ip,j,k,1)/jacip+bb(i,j,k,1)/jac)*jacph
            bxim = 0.5*(bb(im,j,k,1)/jacim+bb(i,j,k,1)/jac)*jacmh
            byip = 0.5*(bb(ip,j,k,2)      +bb(i,j,k,2)    )
            byim = 0.5*(bb(im,j,k,2)      +bb(i,j,k,2)    )
            bzip = 0.5*(bb(ip,j,k,3)/jacip+bb(i,j,k,3)/jac)*jacph
            bzim = 0.5*(bb(im,j,k,3)/jacim+bb(i,j,k,3)/jac)*jacmh

            jacim = jacmh
            jacip = jacph
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
        logical    :: gm_smooth,cvrg_test,lrelax

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

        if (PRESENT(line_relax)) then
          lrelax = line_relax
        else
          lrelax = .false.
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

c diag ****
c$$$        call solverOptionsInit
c$$$
c$$$        solverOptions%iter    = 100
c$$$        solverOptions%tol     = mgtol
c$$$
c$$$        solverOptions%omega   = 0.8 !For point relax
c$$$
c$$$        if (PRESENT(dg)) solverOptions%diag => ddg
c$$$
c$$$        if (PRESENT(ncolors)) solverOptions%ncolors  = ncolors
c$$$
c$$$cc        call assembleSolverHierarchy('jb')
c$$$cc        call assembleSolverHierarchy('gs')
c$$$
c$$$        solverOptions%stp_test = 1 
c$$$        solverOptions%krylov_subspace = solverOptions%iter 
c$$$        call assembleSolverHierarchy('gm')
c$$$cc        call assembleSolverHierarchy('id')
c diag ****

c       Upper_level solver options (MG)

        call solverOptionsInit

        solverOptions%tol      = mgtol
        solverOptions%vcyc     = maxvcyc
        solverOptions%igridmin = 3
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
        if (lrelax) then
          solverOptions%igridmin        = 3
          solverOptions%orderres        = 2
          solverOptions%orderprol       = 2
          solverOptions%mg_mu           = 1    
          solverOptions%mg_line_relax   = .true.
          solverOptions%mg_line_nsweep  = 1
          solverOptions%mg_line_vcyc    = 100
          solverOptions%mg_line_tol     = 1d-1
          solverOptions%mg_line_omega   = 1d0
          solverOptions%mg_line_x       = .true.
          solverOptions%mg_line_y       = .true.
          solverOptions%mg_line_z       = .true.
          solverOptions%mg_line_solve   = "mg"
          solverOptions%mg_line_coarse_solver_depth = 0
        endif

        call assembleSolverHierarchy('mg')

c       Next level solver (smoother)

        call solverOptionsInit

        solverOptions%iter    = nsweep
        solverOptions%tol     = mgtol

        if (lrelax) then
          solverOptions%iter  = 2
        endif

        if (.not.gm_smooth) then
cc          if (PRESENT(dg)) solverOptions%diag => ddg

          !GS
cc          solverOptions%omega   = 1d0
cc          if (PRESENT(ncolors)) solverOptions%ncolors = ncolors
cccc          solverOptions%vertex_based_relax = .true.
cc
cc          call assembleSolverHierarchy('gs')

          !JB
          if (lrelax) then
            solverOptions%omega   = 0.67 !For line relax
          else
            solverOptions%omega   = 0.8 !For point relax
          endif

          call assembleSolverHierarchy('jb')

        else

          !GMRES
          solverOptions%krylov_subspace = nsweep
          solverOptions%stp_test        = 1 

          call assembleSolverHierarchy('gm')

          call solverOptionsInit
          call assembleSolverHierarchy('id') !GMRES preconditioner

        endif

c       Coarsest grid solve for outer MG

        call solverOptionsInit

        solverOptions%tol             = 1d-5
        solverOptions%krylov_subspace = 1000
        solverOptions%iter            = 1000
c diag ****
cc        solverOptions%iter            = 10
c diag ****
        solverOptions%stp_test        = 1 
        solverOptions%omega           = 1d0

c diag ****
cc        call assembleSolverHierarchy('jb')
c diag ****
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
