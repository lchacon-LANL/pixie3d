c module par_int
c #############################################################################
      module par_int

        use orbit

        use grid

        use timestepping,ONLY: dt,dtold

        !Configuration
        logical :: bwgt

        !Field line accumulators
        real(8) :: ct0,ct1  !chi*dt
        real(8) :: temp_s,kernel_s,length_s,b0_s

        !Private variables
        integer,private :: ss_order=1   !Interpolation order in surface of section
     .                    ,kmax=100     !Maximum number of Fourier modes in Fourier Green's f. form

        real(8),private :: orbit_atol=1e-12    !Orbit integrator absolute tolerance
     .                    ,orbit_rtol=1e-12    !Orbit integrator relative tolerance
     .                    ,orbit_dt  =1e-2     !Orbit integrator time step
cc        real(8),private :: orbit_atol=1e-10    !Orbit integrator absolute tolerance
cc     .                    ,orbit_rtol=1e-8     !Orbit integrator relative tolerance
cc     .                    ,orbit_dt  =1e-2     !Orbit integrator time step

        real(8),private :: alpha,dx0,dy0,dz0
     .                    ,t_tol,T_cutoff2,S_cutoff2!,s_max=0d0

        logical,private :: l2d,x_per,y_per,z_per,solenoidal

        real(8),pointer,dimension(:,:,:) :: ll,b0

        real(8),pointer,dimension(:,:,:,:) :: sp_coef

        INTERFACE lagrangian_2d
          module procedure lagrangian_2d_scl,lagrangian_2d_vec
        END INTERFACE

        INTERFACE lagrangian_3d
          module procedure lagrangian_3d_scl,lagrangian_3d_vec
        END INTERFACE

        private :: lagrangian_2d_scl,lagrangian_2d_vec
     .            ,lagrangian_3d_scl,lagrangian_3d_vec
     .            ,lagrangian_2d_init,lagrangian_3d_init

      contains

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     PRIVATE ROUTINES
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c     lagrangian_3d_scl
c     ##############################################################
      subroutine lagrangian_3d_scl(igrid,nx,ny,nz,tt
     .                            ,chi,iout,int_f,lsode_f)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian Green's function temperature integration
c     of tt and returns it in tt.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igrid,iout

        real(8) :: chi

        real(8) :: tt(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg

        real(8) :: x0,y0,z0

        real(8),pointer,dimension(:,:,:) :: tg,zeros

        logical :: have_lsode_f

c     Interfaces

        optional :: int_f

        INTERFACE
          function int_f(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: int_f,fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function int_f
        END INTERFACE

        optional :: lsode_f

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_f(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_f
        END INTERFACE

c     Begin program

        if (chi==0d0) return

        ct0   = chi*dt

        have_lsode_f = PRESENT(lsode_f)

c     Output

        if (my_rank == 0) then
          if (iout == 1) then
            write (*,FMT="(a)",ADVANCE="NO") 'Lagrangian 3D: '
          elseif (iout <= 0) then
            write (*,*) 'Lagrangian step in 3D...'
          endif
        endif

c     Assemble global temperature field

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

        allocate(tg(0:nxg+1,0:nyg+1,0:nzg+1))

        call find_global(tt,tg)

c     Spline temperature field

        call splineFld(tg)

        deallocate(tg)

c     Orbit integration

        do k=1,nz
          do j=1,ny
            do i=1,nx

              !Find location on global mesh

              call fromLocalToGlobalLimits(i,j,k,ig,jg,kg
     .                                    ,igrid,igrid,igrid)

              x0 = grid_params%xg(ig)
              y0 = grid_params%yg(jg)
              z0 = grid_params%zg(kg)

              temp_s   = 0d0    !Initialize temperature
              kernel_s = 0d0    !Initialize kernel integral (for normalization)

              call orbit_find(iout-2,x0,y0,z0,.true.
     .                       ,line_int=int_f,lsode_f=lsode_f) !Forward

              if (have_lsode_f) then
                temp_s   = temp_s   + lsode_out(1)
                kernel_s = kernel_s + lsode_out(2)
              endif

              call orbit_find(iout-2,x0,y0,z0,.false.
     .                       ,line_int=int_f,lsode_f=lsode_f) !Backward

              if (have_lsode_f) then
                temp_s   = temp_s   + lsode_out(1)
                kernel_s = kernel_s + lsode_out(2)
              endif

              if ((my_rank == 0).and.(iout == 1))
     .                    write (*,FMT="(a)",ADVANCE="NO") '+'

              !Fill local temperature
              tt(i,j,k) = temp_s/kernel_s
            enddo
          enddo
        enddo

        if ((my_rank == 0).and.(iout == 1)) write (*,*)

      end subroutine lagrangian_3d_scl

c     lagrangian_3d_vec
c     ##############################################################
      subroutine lagrangian_3d_vec(igrid,nx,ny,nz,vec_in,vec_out
     .                            ,chi,iout,lsode_f)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian Green's function temperature integration
c     of tt and returns it in tt.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igrid,iout

        real(8) :: chi

        real(8) :: vec_in(0:,0:,0:,:),vec_out(0:,0:,0:,:)

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_f(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_f
        END INTERFACE

c     Local variables

        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg,ieq,neq

        real(8) :: x0,y0,z0,temp_n,gkernel_dt,temp_nm,gkernel_2dt,src_dt
     .            ,pkernel_dt,src_2dt,pkernel_2dt

        real(8),pointer,dimension(:,:,:,:) :: vg

c     Begin program

        if (chi==0d0) return

c     Output

        if (my_rank == 0) then
          if (iout == 1) then
            write (*,FMT="(a)",ADVANCE="NO") 'Lagrangian 3D: '
          elseif (iout <= 0) then
            write (*,*) 'Lagrangian step in 3D...'
          endif
        endif

c     Assemble global temperature field

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

        neq = size(vec_in,4)

        if (neq /= 3) then
          call pstop('lagrangian_3d_vec'
     .              ,'Wrong number of input fields')
        endif

        if (size(vec_out,4) /= 4) then
          call pstop('lagrangian_3d_vec'
     .              ,'Wrong number of output fields')
        endif

        allocate(vg(0:nxg+1,0:nyg+1,0:nzg+1,neq)
     .          ,sp_coef (nxg+2,nyg+2,nzg+2,neq))

        call find_global(vec_in,vg)

c     Spline temperature field

        call splineFlds(vg,fldcoef=sp_coef)

        deallocate(vg)

c     Orbit integration

        do k=1,nz
          do j=1,ny
            do i=1,nx

              !Find location on global mesh

              call fromLocalToGlobalLimits(i,j,k,ig,jg,kg
     .                                    ,igrid,igrid,igrid)

              x0 = grid_params%xg(ig)
              y0 = grid_params%yg(jg)
              z0 = grid_params%zg(kg)

              ct0 = chi*dt/b0(i,j,k)
              ct1 = chi*(dt+dtold)/b0(i,j,k)

              call orbit_find(iout-2,x0,y0,z0,.true.,lsode_f=lsode_f) !Forward

              temp_n      = lsode_out(1)
              gkernel_dt  = lsode_out(2)
              temp_nm     = lsode_out(3)
              gkernel_2dt = lsode_out(4)
              src_dt      = lsode_out(5)
              pkernel_dt  = lsode_out(6)
              src_2dt     = lsode_out(7)
              pkernel_2dt = lsode_out(8)

              call orbit_find(iout-2,x0,y0,z0,.false.,lsode_f=lsode_f) !Backward

              temp_n      = temp_n      + lsode_out(1)
              gkernel_dt  = gkernel_dt  + lsode_out(2)
              temp_nm     = temp_nm     + lsode_out(3)
              gkernel_2dt = gkernel_2dt + lsode_out(4)
              src_dt      = src_dt      + lsode_out(5)
              pkernel_dt  = pkernel_dt  + lsode_out(6)
              src_2dt     = src_2dt     + lsode_out(7)
              pkernel_2dt = pkernel_2dt + lsode_out(8)

              if ((my_rank == 0).and.(iout == 1))
     .                    write (*,FMT="(a)",ADVANCE="NO") '+'

              vec_out(i,j,k,1) = temp_n /gkernel_dt
              vec_out(i,j,k,2) = temp_nm/gkernel_2dt
              vec_out(i,j,k,3) = src_dt /pkernel_dt 
              vec_out(i,j,k,4) = src_2dt/pkernel_2dt
            enddo
          enddo
        enddo

        if ((my_rank == 0).and.(iout == 1)) write (*,*)

c     Free up orbit integral work space

        deallocate(sp_coef)

      end subroutine lagrangian_3d_vec

c     lagrangian_2d_scl
c     ##############################################################
      subroutine lagrangian_2d_scl(igrid,nx,ny,nz,tt
     .                            ,chi,iout,int_f,lsode_f)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian temperature integration according to B
c     field (bx,by,bz).
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igrid,iout

        real(8) :: chi

        real(8) :: tt (0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg

        real(8) :: x0,y0,z0

        real(8),pointer,dimension(:,:,:) :: tg,zeros

        logical :: have_lsode_f

c     Interfaces

        optional :: int_f

        INTERFACE
          function int_f(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: int_f,fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function int_f
        END INTERFACE

        optional :: lsode_f

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_f(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_f
        END INTERFACE

c     Begin program

        if (maxval(tt) == minval(tt)) return

        if (chi==0d0) return

        ct0 = chi*dt

        have_lsode_f = PRESENT(lsode_f)

c     Output

        if (my_rank == 0) then
          if (iout == 1) then
            write (*,FMT="(a)",ADVANCE="NO") 'Lagrangian 2D: '
          elseif (iout <= 0) then
            if (solenoidal) then
              write (*,*) 'Solenoidal lagrangian step in 2D...'
            else
              write (*,*) 'LSODE lagrangian step in 2D...'
            endif
          endif
        endif

c     Assemble global temperature field

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

        allocate(tg(0:nxg+1,0:nyg+1,0:nzg+1))

        call find_global(tt,tg)

c     Spline temperature field

        call splineFld(tg)

        deallocate(tg)

c     Orbit integration

        if (have_lsode_f) then

          do k=1,nz
            do j=1,ny
              do i=1,nx

                call fromLocalToGlobalLimits(i,j,k,ig,jg,kg
     .                                      ,igrid,igrid,igrid)

                x0 = grid_params%xg(ig)
                y0 = grid_params%yg(jg)
                z0 = grid_params%zg(kg)

                temp_s   = 0d0        !Initialize temperature
                kernel_s = 0d0        !Initialize kernel integral (for normalization)
                length_s = ll(i,j,k)  !Initialize arc length

                call orbit_find(iout-2,x0,y0,z0,.true.,s_max=length_s
     .                         ,lsode_f=lsode_f)
                temp_s   = lsode_out(1)
                kernel_s = lsode_out(2)

                if ((my_rank == 0).and.(iout == 1))
     .                      write (*,FMT="(a)",ADVANCE="NO") '+'

                tt(i,j,k) = temp_s/kernel_s
              enddo
            enddo
          enddo

        else

          do k=1,nz
            do j=1,ny
              do i=1,nx

                call fromLocalToGlobalLimits(i,j,k,ig,jg,kg
     .                                      ,igrid,igrid,igrid)

                x0 = grid_params%xg(ig)
                y0 = grid_params%yg(jg)
                z0 = grid_params%zg(kg)

                temp_s   = 0d0        !Initialize temperature
                kernel_s = 0d0        !Initialize kernel integral (for normalization)
                length_s = ll(i,j,k)  !Initialize arc length

                call orbit_find(iout-2,x0,y0,z0,.true.,s_max=length_s
     .                         ,line_int=int_f)

                if ((my_rank == 0).and.(iout == 1))
     .                      write (*,FMT="(a)",ADVANCE="NO") '+'

                tt(i,j,k) = temp_s/kernel_s
              enddo
            enddo
          enddo

        endif

c     Free up orbit integral work space

        if ((my_rank == 0).and.(iout == 1)) write (*,*)

      end subroutine lagrangian_2d_scl

c     lagrangian_2d_vec
c     ##############################################################
      subroutine lagrangian_2d_vec(igrid,nx,ny,nz,vec_in,vec_out
     .                            ,chi,iout,lsode_f)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian integration of vector field vec_in according 
c     to B field (bx,by,bz). Result is stored on vector field vec_out:
c       * vec_in: number of fields to integrate (spline)
c       * vec_out: number of integrals to be performed on vec_in fields.
c     Thus, vec_in and vec_out can have different 4th dimension.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igrid,iout

        real(8) :: chi

        real(8) :: vec_in(0:,0:,0:,:),vec_out(0:,0:,0:,:)

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_f(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_f
        END INTERFACE

c     Local variables

        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg,ieq,neq

        real(8) :: x0,y0,z0,temp_n,gkernel_dt,temp_nm,gkernel_2dt,src_dt
     .            ,pkernel_dt,src_2dt,pkernel_2dt

        real(8),pointer,dimension(:,:,:,:) :: vg

c     Begin program

        if (chi==0d0) return

c     Output

        if (my_rank == 0) then
          if (iout == 1) then
            write (*,FMT="(a)",ADVANCE="NO") 'Lagrangian 2D: '
          elseif (iout <= 0) then
            write (*,*) 'LSODE lagrangian step in 2D...'
          endif
        endif

c     Assemble global temperature field

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

        neq = size(vec_in,4)

        if (neq /= 3) then
          call pstop('lagrangian_2d_vec'
     .              ,'Wrong number of input fields')
        endif

        if (size(vec_out,4) /= 4) then
          call pstop('lagrangian_2d_vec'
     .              ,'Wrong number of output fields')
        endif

        allocate(vg(0:nxg+1,0:nyg+1,0:nzg+1,neq)
     .          ,sp_coef (nxg+2,nyg+2,nzg+2,neq))

        call find_global(vec_in,vg)

c     Spline fields

        call splineFlds(vg,fldcoef=sp_coef)

        deallocate(vg)

c     LSODE orbit integration

        do k=1,nz
          do j=1,ny
            do i=1,nx

              call fromLocalToGlobalLimits(i,j,k,ig,jg,kg
     .                                    ,igrid,igrid,igrid)

              x0 = grid_params%xg(ig)
              y0 = grid_params%yg(jg)
              z0 = grid_params%zg(kg)

              length_s  = ll(i,j,k) !Initialize arc length

              ct0 = chi*dt/b0(i,j,k)
              ct1 = chi*(dt+dtold)/b0(i,j,k)

              call orbit_find(iout-2,x0,y0,z0,.true.
     .                       ,s_max=length_s
     .                       ,lsode_f=lsode_f)

              if ((my_rank == 0).and.(iout == 1))
     .                    write (*,FMT="(a)",ADVANCE="NO") '+'

              temp_n      = lsode_out(1)
              temp_nm     = lsode_out(3)
              src_dt      = lsode_out(5)
              src_2dt     = lsode_out(7)

              gkernel_dt  = lsode_out(2)
              gkernel_2dt = lsode_out(4)
              pkernel_dt  = lsode_out(6)
              pkernel_2dt = lsode_out(8)

              vec_out(i,j,k,1) = temp_n /gkernel_dt
              vec_out(i,j,k,2) = temp_nm/gkernel_2dt
              vec_out(i,j,k,3) = src_dt /pkernel_dt 
              vec_out(i,j,k,4) = src_2dt/pkernel_2dt
            enddo
          enddo
        enddo

        if ((my_rank == 0).and.(iout == 1)) write (*,*)

c     Free up orbit integral work space

        deallocate(sp_coef)

      end subroutine lagrangian_2d_vec

c     lagrangian_3d_init
c     ##############################################################
      subroutine lagrangian_3d_init(igx,igy,igz,nx,ny,nz,bb,alp
     .                             ,G_tol,spline_order,nints,B_wgt)

      implicit none

c     --------------------------------------------------------------
c     Initializes 3D Lagrangian step.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igx,igy,igz,spline_order,nints
 
        real(8) :: G_tol,alp
        real(8),INTENT(IN) :: bb(0:nx+1,0:ny+1,0:nz+1,3)

        logical :: B_wgt

c     Local variables

        integer :: nxg,nyg,nzg,i

        real(8),pointer,dimension(:,:,:  ) :: jacg
        real(8),pointer,dimension(:,:,:,:) :: bcar,bg,bcarg,carg

c     Begin program

        l2d = .false.   !Indicate 3D treatment

        alpha = alp

        nxg = grid_params%nxgl(igx)
        nyg = grid_params%nygl(igy)
        nzg = grid_params%nzgl(igz)

        bwgt = B_wgt
        if (bwgt) then
          call pstop('lagrangian_3d_init'
     .              ,'3D version cannot do B.ds yet: needs 1/<B^2>')
        endif

c     Compute Green's integral cut-offs

        t_tol = G_tol !Initialize Green's function tolerance

        T_cutoff2 = 1d0
        S_cutoff2 = 1d0
        do i=1,10
          T_cutoff2 = abs(log(  sqrt(pi)*t_tol*sqrt(T_cutoff2)))
          S_cutoff2 = abs(log(2*sqrt(pi)*t_tol*     S_cutoff2 ))
        enddo

c     Find local cartesian components

        allocate(bcar(0:nx+1,0:ny+1,0:nz+1,3))

        bcar = bb

        call transformVector(igx,igy,igz
     .                      ,0,nx+1,0,ny+1,0,nz+1
     .                      ,bcar(:,:,:,1)
     .                      ,bcar(:,:,:,2)
     .                      ,bcar(:,:,:,3),'cnv','car')

c     Fill global arrays

        allocate(bg   (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,bcarg(0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,carg (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,jacg (0:nxg+1,0:nyg+1,0:nzg+1))

        call find_global(bb  ,bg)
        call find_global(bcar,bcarg)

        call find_global(gmetric%grid(igx)%jac,jacg)
        call find_global(gmetric%grid(igx)%car,carg)

c     Orbit initialization on global array

        call orbit_setup(nxg+2,nyg+2,nzg+2
     .                  ,grid_params%xg
     .                  ,grid_params%yg
     .                  ,grid_params%zg
     .                  ,bg(:,:,:,1)
     .                  ,bg(:,:,:,2)
     .                  ,bg(:,:,:,3)
     .                  ,bcarg(:,:,:,1)
     .                  ,bcarg(:,:,:,2)
     .                  ,bcarg(:,:,:,3)
     .                  ,jacg
     .                  ,carg
     .                  ,bcond
     .                  ,dtime=orbit_dt
     .                  ,atol=orbit_atol  !Field line integration absolute tolerance
     .                  ,rtol=orbit_rtol  !Field line integration relative tolerance
     .                  ,spline_order=spline_order
     .                  ,solen=.false.
     .                  ,lsode_neq=5 + nints  !5 field eqs + number of field integrals
     .                  ,ds_is_Bds=bwgt   !Use B.ds as integration variable
     .                  )

c     Deallocate memory

        deallocate(bg,bcarg,jacg,carg,bcar)

c     End program

      end subroutine lagrangian_3d_init

c     lagrangian_2d_init
c     ##############################################################
      subroutine lagrangian_2d_init(igx,igy,igz,nx,ny,nz,bb,alp,iout
     .                             ,G_tol,spline_order,nints,B_wgt,psi)

      implicit none

c     --------------------------------------------------------------
c     Initializes 2D Lagrangian step.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igx,igy,igz,iout,spline_order,nints

        real(8),INTENT(IN) :: bb (0:nx+1,0:ny+1,0:nz+1,3)
        real(8) :: G_tol,alp

        logical :: B_wgt

        real(8),INTENT(IN),optional :: psi(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg

        real(8) :: x0,y0,z0

        real(8),pointer,dimension(:,:,:  ) :: jacg,psig,psil
        real(8),pointer,dimension(:,:,:,:) :: bcar,bg,bcarg,carg

c     Interfaces

        INTERFACE
          function l_int(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: l_int,fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function l_int
        END INTERFACE

c     Begin program

        l2d = .true.   !Indicate 2D treatment

        solenoidal = PRESENT(psi)

        alpha = alp

        nxg = grid_params%nxgl(igx)
        nyg = grid_params%nygl(igy)
        nzg = grid_params%nzgl(igz)

        bwgt = B_wgt

c     Compute Green's integral cut-offs

        t_tol = G_tol !Initialize Green's function tolerance

        T_cutoff2 = 1d0
        S_cutoff2 = 1d0
        do i=1,10
          T_cutoff2 = abs(log(  sqrt(pi)*t_tol*sqrt(T_cutoff2)))
          S_cutoff2 = abs(log(2*sqrt(pi)*t_tol*     S_cutoff2 ))
        enddo

c     Find local cartesian components

        allocate(bcar(0:nx+1,0:ny+1,0:nz+1,3))

        bcar = bb

        call transformVector(igx,igy,igz
     .                      ,0,nx+1,0,ny+1,0,nz+1
     .                      ,bcar(:,:,:,1)
     .                      ,bcar(:,:,:,2)
     .                      ,bcar(:,:,:,3),'cnv','car')

c     Fill global arrays

        allocate(bg   (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,bcarg(0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,carg (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,jacg (0:nxg+1,0:nyg+1,0:nzg+1))

        call find_global(bcar,bcarg)

        call find_global(gmetric%grid(igx)%jac,jacg)
        call find_global(gmetric%grid(igx)%car,carg)

c     Orbit initialization on global array

        if (solenoidal) then

          if (bwgt) then
            call pstop('lagrangian_2d_init'
     .                ,'Solen version cannot do B.ds yet')
          endif

          allocate(psil(0:nx +1,0:ny +1,0:nz +1)
     .            ,psig(0:nxg+1,0:nyg+1,0:nzg+1))

          psil = psi
          call find_global(psil,psig)

          bg(:,:,:,1:2) = 0d0
          bg(:,:,:,3  ) = -psig

          call orbit_setup(nxg+2,nyg+2,nzg+2
     .                    ,grid_params%xg
     .                    ,grid_params%yg
     .                    ,grid_params%zg
     .                    ,bg(:,:,:,1)
     .                    ,bg(:,:,:,2)
     .                    ,bg(:,:,:,3)
     .                    ,bcarg(:,:,:,1)
     .                    ,bcarg(:,:,:,2)
     .                    ,bcarg(:,:,:,3)
     .                    ,jacg
     .                    ,carg
     .                    ,bcond
     .                    ,dtime=orbit_dt
     .                    ,atol=orbit_atol  !Field line integration absolute tolerance
     .                    ,rtol=orbit_rtol  !Field line integration relative tolerance
     .                    ,spline_order=spline_order
     .                    ,solen=.true.
     .                    ,B_input_is_A=.true.
     .                    )

          deallocate(psig,psil)
        else

          call find_global(bb,bg)

          call orbit_setup(nxg+2,nyg+2,nzg+2
     .                    ,grid_params%xg
     .                    ,grid_params%yg
     .                    ,grid_params%zg
     .                    ,bg(:,:,:,1)
     .                    ,bg(:,:,:,2)
     .                    ,bg(:,:,:,3)
     .                    ,bcarg(:,:,:,1)
     .                    ,bcarg(:,:,:,2)
     .                    ,bcarg(:,:,:,3)
     .                    ,jacg
     .                    ,carg
     .                    ,bcond
     .                    ,dtime=orbit_dt
     .                    ,atol=orbit_atol       !Field line integration absolute tolerance
     .                    ,rtol=orbit_rtol       !Field line integration relative tolerance
     .                    ,spline_order=spline_order
     .                    ,solen=.false.
     .                    ,ds_is_Bds=bwgt        !Use B.ds as integration variable
     .                    ,lsode_dt_lim=.true.
     .                    ,lsode_neq = 3         !Number of equations (only integrate field line)
     .                    ,lsode_s_indep=.true.  !Use arclength as independent variable
     .                    )

        endif

c     Output

        if (my_rank == 0) then
          if (iout == 1) then
            write (*,FMT="(a)",ADVANCE="NO") '2D orbit integrals: '
          elseif (iout <= 0) then
            if (solenoidal) then
              write (*,*) 'Solenoidal 2D orbit integrals...'
            else
              write (*,*) 'LSODE 2D orbit integrals...'
            endif
          endif
        endif

c     Find periodic orbit lengths and <B^2>

        allocate(ll(0:nx+1,0:ny+1,0:nz+1)
     .          ,b0(0:nx+1,0:ny+1,0:nz+1))

        ll = 0d0
        b0 = 1d0  !Set beta to unity by default

        do k=1,nz
          do j=1,ny
            do i=1,nx

              call fromLocalToGlobalLimits(i,j,k,ig,jg,kg,igx,igy,igz)

              x0 = grid_params%xg(ig)
              y0 = grid_params%yg(jg)
              z0 = grid_params%zg(kg)

              length_s  = 0d0    !Initialize orbit length
              b0_s = 0d0

              call orbit_find(iout-2,x0,y0,z0,.true.,line_int=l_int)

              if ((my_rank == 0).and.(iout == 1))
     .                    write (*,FMT="(a)",ADVANCE="NO") '+'

              ll(i,j,k) = length_s                  !Find int(dlambda)=int(B.ds)
              if (b_wgt) b0(i,j,k) = b0_s/length_s  !Find beta
            enddo
          enddo
        enddo

c diag ****
cc        open(unit=110,file='per-orbits-qtys.bin',form='unformatted'
cc     .      ,status='replace')
cc        call contour(ll(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,0,110)
cc        call contour(b0(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc        close(110)
c diag ****

c     Setup periodic LSODE

        if (.not.solenoidal) then

          call orbit_destroy

          call orbit_setup(nxg+2,nyg+2,nzg+2
     .                    ,grid_params%xg
     .                    ,grid_params%yg
     .                    ,grid_params%zg
     .                    ,bg(:,:,:,1)
     .                    ,bg(:,:,:,2)
     .                    ,bg(:,:,:,3)
     .                    ,bcarg(:,:,:,1)
     .                    ,bcarg(:,:,:,2)
     .                    ,bcarg(:,:,:,3)
     .                    ,jacg
     .                    ,carg
     .                    ,bcond
     .                    ,dtime=orbit_dt
     .                    ,atol=orbit_atol      !Field line integration absolute tolerance
     .                    ,rtol=orbit_rtol      !Field line integration relative tolerance
     .                    ,spline_order=spline_order
     .                    ,solen=.false.
     .                    ,ds_is_Bds=bwgt       !Use B.ds as integration variable
     .                    ,lsode_neq = 3+nints  !3 fields + number of field integrals
     .                    ,lsode_to_smax=.true. !Use arclength as independent variable
cc     .                    ,lsode_neq = 3
cc     .                    ,lsode_s_indep=.true.
cc     .                    ,lsode_dt_lim=.true.
     .                    )

        endif

        deallocate(bg,bcarg,jacg,carg,bcar)

c     End program

      end subroutine lagrangian_2d_init

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     PUBLIC ROUTINES
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c     lagr_update
c     ###############################################################
      subroutine lagr_update(igrid,bb,tnm,tn,tnp,ch_par,ch_perp,init
     .                      ,E_diag,iout,G_tol,alph,sp_order,psi)

c     ---------------------------------------------------------------
c     Performs Lagrangian update of temperature at time n+1 (tnp)
c     to enforce parallel transport.
c     ---------------------------------------------------------------

      use timeStepping

      implicit none

c     Call variables

      integer :: igrid,nx,ny,nz,iout
      real(8) :: ch_par,ch_perp
      real(8),dimension(0:,0:,0:)   :: tnm,tn,tnp
      real(8),dimension(0:,0:,0:,:) :: bb
      real(8),dimension(:,:,:),optional :: psi

      real(8),optional :: G_tol,alph

      integer,optional :: sp_order

      logical :: E_diag,init

c     Local variables

      integer :: sorder

      real(8) :: chir,gamma,sgamma,Estar,Enp,dE,En,gtol,alp

      real(8),dimension(:,:,:,:),pointer :: flds,integrals

c     Begin program

      chir = (ch_par-ch_perp)  !Isotropic implementation

      if (chir == 0d0) return

      if (PRESENT(G_tol)) then
        gtol = G_tol
      else
        gtol = 1d-5
      endif

      if (PRESENT(G_tol)) then
        alp = alph
      else
        alp = 2d0
      endif

      if (PRESENT(sp_order)) then
        sorder = sp_order
      else
        sorder = 2
      endif

      nx = size(tnp,1)-2
      ny = size(tnp,2)-2
      nz = size(tnp,3)-2

c     Integrate T_star

      if (E_diag) then
        Estar = integral(igrid,igrid,igrid,nx,ny,nz,tnp)
      endif

c     Lagrangian init

      if (init) then
        call lagrangian_init(igrid,nx,ny,nz,bb,alp,iout
     .                      ,gtol,8,sorder,.true.,psi=psi)
      endif

c     Lagrangian step

      allocate(flds     (0:nx+1,0:ny+1,0:nz+1,3)
     .        ,integrals(0:nx+1,0:ny+1,0:nz+1,4))

      flds(:,:,:,1) = tn
      flds(:,:,:,2) = tnm
      flds(:,:,:,3) = cnp*tnp + cn*tn + cnm*tnm

      call lagrangian_xport(igrid,flds,integrals,chir,iout)

      if (bdf2) then
        sgamma = dt/(dt+dtold)
        gamma = sgamma*sgamma

        tnp =(        integrals(:,:,:,1)  !G(Tn,dt/eps)
     .        - gamma*integrals(:,:,:,2)  !G(Tnm,2dt/eps)
     .        +       integrals(:,:,:,3)  !P(dT/dt,dt/eps)
     .        -sgamma*integrals(:,:,:,4)  !P(dT/dt,2dt/eps)
     .       )/(1d0-gamma)
      else
        tnp = integrals(:,:,:,1) + integrals(:,:,:,3) !Theta scheme
      endif

c     Integrate T_np

      if (E_diag) then
        Enp = integral(igrid,igrid,igrid,nx,ny,nz,tnp)

        if (my_rank == 0) write (*,*) 'Lagrangian step, time=',time
     .                             ,'delta_E=',Enp-Estar
      endif

      deallocate(flds,integrals)

      call lagrangian_destroy

      end subroutine lagr_update

c     lagrangian_init
c     ##############################################################
      subroutine lagrangian_init(igr,nx,ny,nz,bb,alp,iout,G_tol,nints
     .                          ,s_order,B_wgt,psi)

      implicit none

c     --------------------------------------------------------------
c     Initializes Lagrangian step.
c     --------------------------------------------------------------

c     Call variables

      integer :: nx,ny,nz,igr,s_order,nints,iout
 
      real(8) :: G_tol,alp
      real(8),INTENT(IN) :: bb(0:nx+1,0:ny+1,0:nz+1,3)

      real(8),INTENT(IN),optional :: psi(0:nx+1,0:ny+1,0:nz+1)

      logical :: B_wgt

c     Local variables

c     Begin program

      if (nx > 1 .and. ny > 1 .and. nz > 1) then
        call lagrangian_3d_init(igr,igr,igr,nx,ny,nz,bb,alp
     .                       ,G_tol,s_order,nints,B_wgt)
      else
        call lagrangian_2d_init(igr,igr,igr,nx,ny,nz,bb,alp,iout
     .                       ,G_tol,s_order,nints,B_wgt,psi=psi)
      endif

c     End program

      end subroutine lagrangian_init

c     lagrangian_init_from_file
c     ##############################################################
      subroutine lagrangian_init_from_file(file,igx,igy,igz,nx,ny,nz,bb
     .                                    ,alp,G_tol,spline_order,B_wgt)

      use iosetup

      implicit none

c     --------------------------------------------------------------
c     Initializes Lagrangian step from PIXIE3D's poincare output file.
c     Provides bb as output
c     --------------------------------------------------------------

c     Call variables

        character(*) :: file

        integer,INTENT(IN)  :: nx,ny,nz,igx,igy,igz,spline_order

        real(8) :: G_tol,alp
cc        real(8),INTENT(OUT) :: bb(0:nx+1,0:ny+1,0:nz+1,3)
        real(8),pointer,dimension(:,:,:,:) :: bb

        logical :: B_wgt

c     Local variables

        integer :: upixie=12345,ierr,nfiles

        integer :: nxg,nyg,nzg,bcnd(6)

        real(8),pointer,dimension(:)       :: xx,yy,zz
        real(8),pointer,dimension(:,:,:)   :: jac,jacg
        real(8),pointer,dimension(:,:,:,:) :: bcar,bg,bcarg,carg,car

        logical :: tor

        integer,pointer,dimension(:) :: mupixie

c     Begin program

        t_tol = G_tol  !Initialize Green's function tolerance

        alpha = alp

        bwgt = B_wgt

c     Check file configuration

        ierr = openFiles(my_rank,np,upixie,trim(file),mupixie)

        if (ierr /= 0) call pstop('openGraphicsfiles'
     .                           ,'Cannot find files')

c     Read global data

        nfiles = size(mupixie)

        if (nfiles /= np) then
          call read_global_data_async
        else
          call read_global_data_sync
        endif

c     Close files

        call closeFiles(mupixie)

c     Orbit initialization on global array

        call orbit_setup(nxg+2,nyg+2,nzg+2
     .                  ,xx,yy,zz
     .                  ,bg(:,:,:,1)
     .                  ,bg(:,:,:,2)
     .                  ,bg(:,:,:,3)
     .                  ,bcarg(:,:,:,1)
     .                  ,bcarg(:,:,:,2)
     .                  ,bcarg(:,:,:,3)
     .                  ,jacg
     .                  ,carg
     .                  ,bcond
     .                  ,dtime=1d-1
     .                  ,atol=orbit_atol  !Field line integration absolute tolerance
     .                  ,rtol=orbit_rtol  !Field line integration relative tolerance
     .                  ,spline_order=spline_order
     .                  ,solen=.false.
     .                  ,ds_is_Bds=bwgt   !Use B.ds as integration variable
     .                  )

c     Deallocate memory

        deallocate(bg,bcarg,jacg,carg,xx,yy,zz)

c     End program

      contains

c     read_global_data_sync
c     ##################################################################
      subroutine read_global_data_sync

c     ------------------------------------------------------------------
c     Reads global data from PIXIE3D files in a synchronous manner (i.e.,
c     one file per proc).
c     ------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer :: nxl,nyl,nzl,ilmg,ihpg,jlmg,jhpg,klmg,khpg,funit

c     Begin program

        if (my_rank == 0) then
          write (*,*) 'Synchronous input in Lagrangian step'
        endif

        funit = mupixie(1)

        read (funit) nxg,nyg,nzg,bcnd,tor

        if (tor .and. (coords /= 'vmc' .and. coords /= 'tor')) then
          call closeFiles(mupixie)
          call pstop('lagrangian_init','Coordinates are not toroidal')
        endif

        if (sum(bcnd-bcond) /= 0) then
          call closeFiles(mupixie)
          call pstop('lagrangian_init','BC conf. does not agree')
        endif

        if (nxd /= nxg .or. nyd /= nyg .or. nzd /= nzg) then
          call closeFiles(mupixie)
          write (*,*) 'Grid size=',nxg,nyg,nzg
          call pstop('lagrangian_init'
     .              ,'Global grid size does not agree')
        endif

c     Read global position arrays

        allocate(xx(nxg+2),yy(nyg+2),zz(nzg+2))

        read (funit) xx
        read (funit) yy
        read (funit) zz

c     Local limits

        if (np > 1) then
          read (funit) ilmg,ihpg,jlmg,jhpg,klmg,khpg
        else
          ilmg = 0
          ihpg = nxg+1

          jlmg = 0
          jhpg = nyg+1

          klmg = 0
          khpg = nzg+1
        endif

        nxl = ihpg-ilmg-1
        nyl = jhpg-jlmg-1
        nzl = khpg-klmg-1

        if (nx /= nxl .or. nyl /= nyl .or. nz /= nzl) then
          write (*,*) 'Local grid size=',nxl,nyl,nzl
          call pstop('lagrangian_init'
     .              ,'Local grid size does not agree')
        endif

c     Allocate local arrays

        allocate(bcar(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,jac (0:nx+1,0:ny+1,0:nz+1)
     .          ,car (0:nx+1,0:ny+1,0:nz+1,3))

c     Fill local arrays

        read (funit) bb(:,:,:,1)
        read (funit) bb(:,:,:,2)
        read (funit) bb(:,:,:,3)
        read (funit) bcar(:,:,:,1)
        read (funit) bcar(:,:,:,2)
        read (funit) bcar(:,:,:,3)
        read (funit) jac
        read (funit) car

c     Allocate global arrays

        allocate(bg   (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,bcarg(0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,carg (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,jacg (0:nxg+1,0:nyg+1,0:nzg+1))

c     Fill global arrays

        call find_global(bb,bg)
        call find_global(bcar,bcarg)

        call find_global(jac,jacg)
        call find_global(car,carg)

c     Update grid metrics

        call defineGridMetric(grid_params,xcar=car)

        gmetric => grid_params%gmetric

c     Deallocate memory

        deallocate(bcar,jac,car)

      end subroutine read_global_data_sync

c     read_global_data_async
c     ##################################################################
      subroutine read_global_data_async

c     ------------------------------------------------------------------
c     Reads global data from PIXIE3D files in an asynchronous manner (i.e.,
c     global solution is found in all procs by reading multiple files).
c     ------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer :: ifl,ilmg,ihpg,jlmg,jhpg,klmg,khpg,nx,ny,nz

c     Begin program

        if (my_rank == 0) then
          write (*,*) 'Asynchronous input in Lagrangian step'
        endif

c     Read setup info

        do ifl=1,nfiles
          read (mupixie(ifl)) nxg,nyg,nzg,bcnd,tor
        enddo

        if (tor .and. (coords /= 'vmc' .and. coords /= 'tor')) then
          call closeFiles(mupixie)
          call pstop('lagrangian_init','Coordinates are not toroidal')
        endif

        if (sum(bcnd-bcond) /= 0) then
          call closeFiles(mupixie)
          call pstop('lagrangian_init','BC conf. does not agree')
        endif

        if (nxd /= nxg .or. nyd /= nyg .or. nzd /= nzg) then
          call closeFiles(mupixie)
          write (*,*) 'Grid size=',nxg,nyg,nzg
          call pstop('lagrangian_init'
     .              ,'Global grid size does not agree')
        endif

c     Allocate global arrays

        allocate(xx(nxg+2),yy(nyg+2),zz(nzg+2))

        allocate(bg   (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,bcarg(0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,carg (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,jacg (0:nxg+1,0:nyg+1,0:nzg+1  ))

c     Read global arrays

        do ifl=1,nfiles

          read (mupixie(ifl)) xx
          read (mupixie(ifl)) yy
          read (mupixie(ifl)) zz

          if (nfiles > 1) then
            read (mupixie(ifl)) ilmg,ihpg,jlmg,jhpg,klmg,khpg
          else
            ilmg = 0
            ihpg = nxg+1

            jlmg = 0
            jhpg = nyg+1

            klmg = 0
            khpg = nzg+1
          endif

          read (mupixie(ifl)) bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,1)
          read (mupixie(ifl)) bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,2)
          read (mupixie(ifl)) bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)
          read (mupixie(ifl)) bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,1)
          read (mupixie(ifl)) bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,2)
          read (mupixie(ifl)) bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)
          read (mupixie(ifl)) jacg (ilmg:ihpg,jlmg:jhpg,klmg:khpg)
          read (mupixie(ifl)) carg (ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)

        enddo

c     Assign local B arrays

        bb = bg(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)

c     Update grid metrics

        nxl = ihipg-ilomg-1
        nyl = jhipg-jlomg-1
        nzl = khipg-klomg-1

        call defineGridMetric(grid_params
     .           ,xcar=carg(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:))

        gmetric => grid_params%gmetric

c     Deallocate memory

      end subroutine read_global_data_async

      end subroutine lagrangian_init_from_file

c     lagrangian_destroy
c     ##############################################################
      subroutine lagrangian_destroy

      implicit none

c     --------------------------------------------------------------
c     Destroys storage for Lagrangian step
c     --------------------------------------------------------------

c     Call variables

c     Local variables

      integer :: alloc_stat

c     Begin program

        call orbit_destroy

        deallocate(ll,stat=alloc_stat)
        deallocate(b0,stat=alloc_stat)

c     End program

      end subroutine lagrangian_destroy

c     lagrangian_xport
c     ##############################################################
      subroutine lagrangian_xport(igrid,in_flds,out_flds,chi,iout)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian temperature integration according to B
c     field (bx,by,bz).
c     --------------------------------------------------------------

c     Call variables

        integer :: igrid,iout

        real(8) :: chi,in_flds(:,:,:,:),out_flds(:,:,:,:)

c     Local variables

        integer ::  nx,ny,nz

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_2d_bdf2_3eq(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_2d_bdf2_3eq
        END INTERFACE

        INTERFACE
          subroutine lsode_3d_bdf2(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_3d_bdf2
        END INTERFACE

c     Begin program

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        if (l2d) then
          call lagrangian_2d_vec(igrid,nx,ny,nz,in_flds,out_flds
     .                          ,chi,iout,lsode_f=lsode_2d_bdf2_3eq)
        else
          call lagrangian_3d_vec(igrid,nx,ny,nz,in_flds,out_flds
     .                          ,chi,iout,lsode_f=lsode_3d_bdf2)
        endif

      end subroutine lagrangian_xport

c     lagrangian_T
c     ##############################################################
      subroutine lagrangian_T(igrid,tt,chi,iout)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian temperature integration according to B
c     field (bx,by,bz).
c     --------------------------------------------------------------

c     Call variables

        integer :: igrid,iout

        real(8) :: chi,tt(:,:,:)

c     Local variables

        integer ::  nx,ny,nz

        INTERFACE
          function t_int_per(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: t_int_per,fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function t_int_per
        END INTERFACE

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_3d_T(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_3d_T
        END INTERFACE

        INTERFACE
          subroutine lsode_2d_T_3eq(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_2d_T_3eq
        END INTERFACE

c     Begin program

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        if (l2d) then
          if (solenoidal) then
            call lagrangian_2d(igrid,nx,ny,nz,tt
     .                        ,chi,iout,int_f=t_int_per)
          else
            call lagrangian_2d(igrid,nx,ny,nz,tt
     .                        ,chi,iout,lsode_f=lsode_2d_T_3eq)
          endif
        else
          call lagrangian_3d(igrid,nx,ny,nz,tt
     .                      ,chi,iout,lsode_f=lsode_3d_T)
        endif

      end subroutine lagrangian_T

c     lagrangian_S
c     ##############################################################
      subroutine lagrangian_S(igrid,ss,chi,iout)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian temperature integration according to B
c     field (bx,by,bz).
c     --------------------------------------------------------------

c     Call variables

        integer :: igrid,iout

        real(8) :: chi,ss(:,:,:)

c     Local variables

        integer :: nx,ny,nz

        INTERFACE
          function s_int_per(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: s_int_per,fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function s_int_per
        END INTERFACE

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_3d_S(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_3d_S
        END INTERFACE

        INTERFACE
          subroutine lsode_2d_S_3eq(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_2d_S_3eq
        END INTERFACE

c     Begin program

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        if (l2d) then
          if (solenoidal) then
            call lagrangian_2d(igrid,nx,ny,nz,ss
     .                        ,chi,iout,int_f=s_int_per)
          else
            call lagrangian_2d(igrid,nx,ny,nz,ss
     .                        ,chi,iout,lsode_f=lsode_2d_S_3eq)
          endif
        else
          call lagrangian_3d(igrid,nx,ny,nz,ss
     .                      ,chi,iout,lsode_f=lsode_3d_S)

        endif

      end subroutine lagrangian_S

c     tkernel_3d
c     ###############################################################
      function tkernel_3d(ss,ct0,ex_loop) result(tkernel)

      implicit none

c     Call variables

      real(8) :: tkernel,ss,ct0
      logical :: ex_loop

c     Local variables

c     Begin program

      if (alpha == 2d0) then
        tkernel = 0.5/sqrt(pi*ct0)*exp(-0.25*ss**2/ct0) !ct0=chi_par*t0

cc        ex_loop = (abs(ss) > sqrt(s_max**2 + 4*T_cutoff2*ct0))
        ex_loop = (abs(ss) > sqrt(4*T_cutoff2*ct0))
      elseif (alpha == 1d0) then
        tkernel = 1d0/(pi*ct0*(1+(abs(ss)/ct0)**2))

cc        ex_loop = (abs(ss) >  s_max + 2/pi/t_tol*ct0)
        ex_loop = (abs(ss) >  2d0/pi/t_tol*ct0)
      else
        write (*,*) 'Fractional order',alpha,' not implemented'
        stop
      endif

c     End program

      end function tkernel_3d

c     tkernel_per
c     ###############################################################
      function tkernel_per(ss,ct0) result(tkernel)

      implicit none

c     Call variables

      real(8) :: tkernel,ss,ct0

c     Local variables

      integer :: kk,k0
      real(8) :: twopi_L,dtk,eta,smax,ismax
      logical :: ex_loop,ex_loop2

c     Begin program

      smax = length_s
      ismax = 1d0/smax

      twopi_L = 2d0*pi*ismax

      eta = ct0*(twopi_L)**alpha

      if (eta < 1d0) then

        tkernel = tkernel_3d(ss,ct0,ex_loop)

        kk = 1
        do 
          tkernel = tkernel + tkernel_3d(ss-kk*smax,ct0,ex_loop)
     .                      + tkernel_3d(ss+kk*smax,ct0,ex_loop2)

          if (ex_loop .and. ex_loop2) exit
          kk = kk + 1
        enddo

      else

        !Compute iteration limit
        k0 = ceiling((abs(log(0.5*t_tol))/eta)**(1d0/alpha))

        !Accumulate periodic Green's function
        tkernel = 1d0

        do kk=1,k0
          dtk = exp(-eta*kk**alpha)
          tkernel = tkernel + 2d0*dtk*cos(twopi_L*kk*ss)
        enddo

        tkernel = tkernel*ismax

      endif

c     End program

      end function tkernel_per

c     skernel_3d
c     ###############################################################
      function skernel_3d(ss,ct0,ex_loop) result(skernel)

      implicit none

c     Call variables

      real(8) :: skernel,ss,ct0
      logical :: ex_loop

c     Local variables

      real(8) :: eta,z

#if defined(pgf90)
      real(8) :: erfc
#endif

c     Begin program

      z = abs(ss)

      if (alpha == 2d0) then
        eta = 0.5*z/sqrt(ct0) !ct0=chi_par/chi_perp*t0
        skernel = (exp(-eta**2)/sqrt(pi)-eta*erfc(eta))/sqrt(ct0)

cc        ex_loop = (z > sqrt(s_max**2 + 4*S_cutoff2*ct0))
        ex_loop = (z > sqrt(4*S_cutoff2*ct0))
cc      elseif (alpha == 1d0) then
cc        skernel = 1d0/(pi*ct0*(1+(abs(ss)/ct0)**2))
cc
cc        slim = s_max + 2/pi/t_tol*ct0
cc
cc        ex_loop = (abs(ss) > slim)
      else
        write (*,*) 'Fractional order not implemented'
        stop
      endif

c     End program

      end function skernel_3d

c     skernel_per
c     ###############################################################
      function skernel_per(ss,ct0) result(skernel)

      implicit none

c     Call variables

      real(8) :: skernel,ss,ct0

c     Local variables

      integer :: kk,k0
      real(8) :: twopi_L,dtk,eta,smax,ismax
      logical :: ex_loop,ex_loop2

c     Begin program

      smax = length_s
      ismax = 1d0/smax

      twopi_L = 2d0*pi*ismax

      eta = ct0*(twopi_L)**alpha

      if (eta < 1d0/sqrt(t_tol*abs(log(t_tol)))) then

        skernel = skernel_3d(ss,ct0,ex_loop)

        kk = 1
        do
          skernel = skernel + skernel_3d(ss-kk*smax,ct0,ex_loop)
     .                      + skernel_3d(ss+kk*smax,ct0,ex_loop2)

          if (ex_loop .and. ex_loop2) exit
          kk = kk + 1
        enddo

      else

        !Compute iteration limit
        k0 = ceiling((2d0/t_tol/eta)**(1d0/alpha))
        k0 = min(k0,kmax) !To avoid aliasing errors

        !Accumulate periodic Green's function
        skernel = 1d0

        do kk=1,k0
          eta = ct0*(twopi_L*kk)**alpha
          dtk = (1d0-exp(-eta))/eta

          skernel = skernel + 2d0*dtk*cos(twopi_L*kk*ss)
        enddo

        skernel = skernel*ismax

      endif

c     End program

      end function skernel_per

c     per_orbit_init
c     ###############################################################
      subroutine per_orbit_init(x0,x,y0,y,z0,z)

c     ---------------------------------------------------------------
c     Initializes periodic orbit check.
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: x0,x,y0,y,z0,z

c     Local variables

c     Begin program

        orb_init = .false.

        dx0 = x-x0 !Determines initial direction of orbit
        dy0 = y-y0 !Determines initial direction of orbit
        dz0 = z-z0 !Determines initial direction of orbit

cc        x_per =               (bcond(1)==PER.and.grid_params%nxgl(1)>1)
cc        if (.not.x_per) y_per=(bcond(3)==PER.and.grid_params%nygl(1)>1)
cc        if (.not.y_per) z_per=(bcond(5)==PER.and.grid_params%nzgl(1)>1)

cc        if ((.not.x_per).and.(.not.y_per).and.(.not.z_per)) then
        !Failsafe in case of non-periodic domains
        !Select periodic direction based on largest displacement from surface
          x_per = (abs(dx0) > abs(dy0)) .and. (abs(dx0) > abs(dz0))
          y_per = (abs(dy0) > abs(dx0)) .and. (abs(dy0) > abs(dz0))
          z_per = (abs(dz0) > abs(dy0)) .and. (abs(dz0) > abs(dx0))
cc        endif

cc        write (*,*) 'per_orbit=',x0,y0,z0,x_per,y_per,z_per

        call ss_init(ss_order,x_per,x0,y_per,y0,z_per,z0)

      end subroutine per_orbit_init

c     per_orbit_chk
c     ###############################################################
      subroutine per_orbit_chk(so,ss,xo,x,yo,y,zo,z,ts,xs,ys,zs)

c     ---------------------------------------------------------------
c     Check for end of periodic orbit. Returns ts <> 0d0 if end is
c     found.
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: so,ss,xo,x,yo,y,zo,z,xs,ys,zs,ts

c     Local variables

        real(8) :: xoo=0d0,yoo=0d0,zoo=0d0,soo=0d0

c     Begin program

        ts = 0d0

        if (    (x_per.and.(x-xo)*dx0 > 0d0)
     .      .or.(y_per.and.(y-yo)*dy0 > 0d0)
     .      .or.(z_per.and.(z-zo)*dz0 > 0d0)) then

          call sschk(ss_order
     .              ,soo,so,ss
     .              ,xoo,xo,x
     .              ,yoo,yo,y
     .              ,zoo,zo,z,xs,ys,zs,ts)  !--> Returns ts <> 0 if intersects

        endif

      end subroutine per_orbit_chk

      end module par_int

c t_int
c ###################################################################
      function t_int(fwd,xo,x,yo,y,zo,z,so,ss) result(ex_loop)

c -------------------------------------------------------------------
c     Performs line integral of temperature.
c -------------------------------------------------------------------

        use par_int

        implicit none

c Call variables

        logical :: ex_loop,fwd
        real(8) :: xo,x,yo,y,zo,z,so,ss

c Local variables

        integer :: ierror
        real(8) :: xh,yh,zh,sh,tt,tk,it_max=0d0
        save :: it_max

c Begin program

c Find middle points

        xh = 0.5*(x +xo)
        yh = 0.5*(y +yo)
        zh = 0.5*(z +zo)
        sh = 0.5*(ss+so)

c Find temperature

        tt = evalFld(xh,yh,zh,ierror)

c Find kernel

        tk = tkernel_3d(sh,ct0,ex_loop)

c Integrate temperature

        temp_s   = temp_s   + tt*tk*(ss-so)
        kernel_s = kernel_s +    tk*(ss-so)

      end function t_int

c t_int_per
c ###################################################################
      function t_int_per(fwd,xo,x,yo,y,zo,z,so,ss) result(exloop)

c -------------------------------------------------------------------
c     Performs line integral of temperature in periodic domain.
c -------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      logical :: exloop,fwd
      real(8) :: xo,x,yo,y,zo,z,so,ss

c Local variables

      integer :: ierror
      real(8) :: xh,yh,zh,sh,tt,tk,xs,ys,zs,ts

c Begin program

      exloop = .false.

ccc Check whether orbit should terminate (x0,y0,z0 is initial orbit position)
cc
cc      if (orb_init) call per_orbit_init(xo,x,yo,y,zo,z)
cc
ccc Check whether orbit should terminate
cc
cc      call per_orbit_chk(so,ss,xo,x,yo,y,zo,z,ts,xs,ys,zs)  !--> Returns ts <> 0 if intersects
cc
cc      exloop = (ts /= 0d0)
cc
ccc Update current position
cc
cc      if (exloop) then
cc        x = xs ; y = ys ; z = zs
cc      else
        ts = ss
cc      endif

c Find middle points

      xh = 0.5*(x +xo)
      yh = 0.5*(y +yo)
      zh = 0.5*(z +zo)
      sh = 0.5*(ts+so)

c Find temperature

      tt = evalFld(xh,yh,zh,ierror)

c Find kernel

      tk = tkernel_per(sh,ct0)

c Integrate temperature

      temp_s   = temp_s   + tt*tk*(ts-so)
      kernel_s = kernel_s +    tk*(ts-so)

      end function t_int_per

c s_int
c ###################################################################
      function s_int(fwd,xo,x,yo,y,zo,z,so,ss) result(ex_loop)

c -------------------------------------------------------------------
c     Performs line integral of temperature.
c -------------------------------------------------------------------

        use par_int

        implicit none

c Call variables

        logical :: ex_loop,fwd
        real(8) :: xo,x,yo,y,zo,z,so,ss

c Local variables

        integer :: ierror
        real(8) :: xh,yh,zh,sh,tt,tk,it_max=0d0
        save :: it_max

c Begin program

c Find middle points

        xh = 0.5*(x +xo)
        yh = 0.5*(y +yo)
        zh = 0.5*(z +zo)
        sh = 0.5*(ss+so)

c Find temperature

        tt = evalFld(xh,yh,zh,ierror)

c Find kernel

        tk = skernel_3d(sh,ct0,ex_loop)

c Integrate temperature

        temp_s   = temp_s   + tt*tk*(ss-so)
        kernel_s = kernel_s +    tk*(ss-so)

      end function s_int

c s_int_per
c ###################################################################
      function s_int_per(fwd,xo,x,yo,y,zo,z,so,ss) result(exloop)

c -------------------------------------------------------------------
c     Performs line integral of temperature in periodic domain.
c -------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      logical :: exloop,fwd
      real(8) :: xo,x,yo,y,zo,z,so,ss

c Local variables

      integer :: ierror
      real(8) :: xh,yh,zh,sh,tt,tk,xs,ys,zs,ts

c Begin program

      exloop = .false.

ccc Check whether orbit should terminate
cc
cc      if (orb_init) call per_orbit_init(xo,x,yo,y,zo,z)
cc
ccc Check whether orbit should terminate
cc
cc      call per_orbit_chk(so,ss,xo,x,yo,y,zo,z,ts,xs,ys,zs)  !--> Returns ts <> 0 if intersects
cc
cc      exloop = (ts /= 0d0)
cc
ccc Update current position
cc
cc      if (exloop) then
cc        x = xs ; y = ys ; z = zs
cc      else
        ts = ss
cc      endif

c Find middle points

      xh = 0.5*(x +xo)
      yh = 0.5*(y +yo)
      zh = 0.5*(z +zo)
      sh = 0.5*(ts+so)

c Find temperature

      tt = evalFld(xh,yh,zh,ierror)

c Find kernel

      tk = skernel_per(sh,ct0)

c Integrate temperature

      temp_s   = temp_s   + tt*tk*(ts-so)
      kernel_s = kernel_s +    tk*(ts-so)

      end function s_int_per

c l_int
c ###################################################################
      function l_int(fwd,xo,x,yo,y,zo,z,so,ss) result(exloop)

c -------------------------------------------------------------------
c     Performs integral of arc length. On input, old and new logical
c     positions and field line length
c -------------------------------------------------------------------

      use par_int

      use equilibrium

      implicit none

c Call variables

      logical :: exloop,fwd
      real(8) :: xo,x,yo,y,zo,z,so,ss

c Local variables

      integer :: ierr
      real(8) :: xh,yh,zh,sh,bb,bx,by,bz,xs,ys,zs,ts,b2,ds

c Begin program

      if (orb_init) call per_orbit_init(xo,x,yo,y,zo,z)

c Check whether orbit should terminate

      call per_orbit_chk(so,ss,xo,x,yo,y,zo,z,ts,xs,ys,zs) !--> Returns ts <> 0 if intersects

      exloop = (ts /= 0d0)

c Update current position

      if (exloop) then
        x = xs ; y = ys ; z = zs
      else
        ts = ss
      endif

c Add arc length

      ds = (ts-so)

      length_s = length_s + ds  !Integrate dlambda

c Add B-weighed arc-length

      xh = .5*(xo+x)
      yh = .5*(yo+y)
      zh = .5*(zo+z)

      call getB(xh,yh,zh,bx,by,bz,.false.,.true.,ierr)

      b2 = bx**2 + by**2 + bz**2

      b0_s = b0_s + ds/b2       !Integrate dlambda/B^2

      end function l_int

c lsode_3d_T
c #########################################################################
      subroutine lsode_3d_T(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      use equilibrium

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,tt,it_max=0d0
      save :: it_max

c Begin program

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f(5,t,y(1:5),ydot(1:5),xx1=x1,yy1=y1,zz1=z1)

c Call additional integrators

      s1 = abs(y(5))

      tt = evalFld(x1,y1,z1,ierror)

      ydot(7) = tkernel_3d(s1,ct0,exit_loop)*ydot(5)  !To normalize

      ydot(6) = tt*ydot(7)

c End program

      end subroutine lsode_3d_T

c lsode_2d_T_3eq
c #########################################################################
      subroutine lsode_2d_T_3eq(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      use equilibrium

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,tt,kernel,it_max=0d0
      save :: it_max

c Begin program

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),xx1=x1,yy1=y1,zz1=z1)

c Call additional integrators

      tt = evalFld(x1,y1,z1,ierror)
      kernel = tkernel_per(t,ct0)

      ydot(4) = tt*kernel

      ydot(5) = kernel  !To normalize

c End program

      end subroutine lsode_2d_T_3eq

c lsode_3d_S
c #########################################################################
      subroutine lsode_3d_S(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      use equilibrium

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,tt,it_max=0d0
      save :: it_max

c Begin program

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f(5,t,y(1:5),ydot(1:5),xx1=x1,yy1=y1,zz1=z1)

c Call additional integrators

      s1 = abs(y(5))

      tt = evalFld(x1,y1,z1,ierror)

      ydot(7) = skernel_3d(s1,ct0,exit_loop)*ydot(5)  !To normalize

      ydot(6) = tt*ydot(7)

c End program

      end subroutine lsode_3d_S

c lsode_2d_S_3eq
c #########################################################################
      subroutine lsode_2d_S_3eq(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      use equilibrium

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,tt,it_max=0d0
      save :: it_max

c Begin program

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),xx1=x1,yy1=y1,zz1=z1)

c Call additional integrators

      tt = evalFld(x1,y1,z1,ierror)
      ydot(5) = skernel_per(t,ct0)

      ydot(4) = tt*ydot(5)

c End program

      end subroutine lsode_2d_S_3eq

c lsode_3d_bdf2
c #########################################################################
      subroutine lsode_3d_bdf2(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      use equilibrium

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,flds(size(sp_coef,4)),B2
      logical :: exl1,exl2,exl3,exl4

c Begin program

      if (neq /= 13) then
        call pstop('lsode_3d_bdf2','Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f(5,t,y(1:5),ydot(1:5),xx1=x1,yy1=y1,zz1=z1
     .                  ,B2=B2)

      if (.not.bwgt) B2=1d0

c Call additional integrators

      s1 = abs(y(5))

      flds = evalFlds(x1,y1,z1,sp_coef)

      !Tn
      ydot(7) = tkernel_3d(s1,ct0,exl1)*ydot(5)

      ydot(6) = flds(1)*ydot(7)

      !Tnm
      ydot(9) = tkernel_3d(s1,ct1,exl2)*ydot(5)

      ydot(8) = flds(2)*ydot(9)

      !Src @ dt
      ydot(11) = skernel_3d(s1,ct0,exl3)*ydot(5)/B2

      ydot(10) = flds(3)*ydot(11)

      !Src @ dt+dtold
      ydot(13) = skernel_3d(s1,ct1,exl4)*ydot(5)/B2

      ydot(12) = flds(3)*ydot(13)

c Are we done?

      exit_loop = exl1.and.exl2.and.exl3.and.exl4

c End program

      end subroutine lsode_3d_bdf2

c lsode_2d_bdf2_3eq
c #########################################################################
      subroutine lsode_2d_bdf2_3eq(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      use equilibrium

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,it_max=0d0,flds(size(sp_coef,4)),B2
      save :: it_max

c Begin program

      if (neq /= 11) then
        call pstop('lsode_2d_bdf2_3eq','Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),xx1=x1,yy1=y1,zz1=z1
     .                      ,B2=B2)

      if (.not.bwgt) B2=1d0

c Call additional integrators

      flds = evalFlds(x1,y1,z1,sp_coef)

      !Tn
      ydot(5) = tkernel_per(t,ct0)

      ydot(4) = flds(1)*ydot(5)

      !Tnm
      ydot(7) = tkernel_per(t,ct1)

      ydot(6) = flds(2)*ydot(7)

      !Src @ dt
      ydot(9) = skernel_per(t,ct0)/B2

      ydot(8) = flds(3)*ydot(9)

      !Src @ dt+dtold
      ydot(11) = skernel_per(t,ct1)/B2

      ydot(10) = flds(3)*ydot(11)

c End program

      end subroutine lsode_2d_bdf2_3eq

c lsode_2d_cn_3eq
c #########################################################################
      subroutine lsode_2d_cn_3eq(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      use equilibrium

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,it_max=0d0,flds(size(sp_coef,4)),B2
      save :: it_max

c Begin program

      if (neq /= 7) then
        call pstop('lsode_2d_cn_3eq'
     .            ,'Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),xx1=x1,yy1=y1,zz1=z1
     .                      ,B2=B2)

      if (.not.bwgt) B2=1d0

c Call additional integrators

      flds = evalFlds(x1,y1,z1,sp_coef)

      !Tn
      ydot(5) = tkernel_per(t,ct0)

      ydot(4) = flds(1)*ydot(5)

      !Src @ dt
      ydot(7) = skernel_per(t,ct0)/B2

      ydot(6) = flds(3)*ydot(7)

c End program

      end subroutine lsode_2d_cn_3eq
