module rw_bim_mod
   use math, only: pi, ellipticK, ellipticE
   use elliptic
   use grid_def_st
#if defined(petsc)
   use grid_mpi, only: find_global_nobc,mpi_comm_rank,mpierr
#else
   use grid_mpi, only: find_global_nobc 
#endif
   implicit none
   private

   integer, parameter :: &
      REGULAR = 0, &
      SINGULAR = REGULAR + 1

   type :: bim_data_t
      integer :: n
      double precision :: dxi
      double precision, allocatable, dimension(:) :: R, Z, xi, Bn, phi
      double precision, allocatable, dimension(:,:) :: Jnorm, L, K, Kinv, KinvL
   end type bim_data_t

   ! for verification tests
   double precision :: r0,Rmaj,a,rho0,z0

   logical :: bim_dump
   
   type(bim_data_t) :: bd

   ! general 18-point gauss-legendre quadrature rule for singular integrands
   ! containing one or more singular and/or hypersingular terms
   ! A*P_n(x) + B*P_n(x)*ln(|x-x0|) + C*P_n(x)/(x0-x) + D*P_n(x)/(x0-x)^2
   ! nodes
   double precision, parameter :: &
      xglg_18(18) = (/&
      -9.915651684209309e-01,&
      -9.558239495713978e-01,&
      -8.926024664975557e-01,&
      -8.037049589725231e-01,&
      -6.916870430603532e-01,&
      -5.597708310739475e-01,&
      -4.117511614628426e-01,&
      -2.518862256915055e-01,&
      -8.477501304173529e-02,&
       8.477501304173529e-02,&
       2.518862256915055e-01,&
       4.117511614628426e-01,&
       5.597708310739475e-01,&
       6.916870430603532e-01,&
       8.037049589725231e-01,&
       8.926024664975557e-01,&
       9.558239495713978e-01,&
       9.915651684209309e-01 /)
   ! weights
   double precision, parameter :: &
      wglg_18(18) = (/&
      -5.069710800869371e+02,&
       1.385924948889776e+03,&
      -1.729512527512267e+03,&
       1.400159558133540e+03,&
      -7.916434221273678e+02,&
       3.122202458902182e+02,&
      -8.069146562278169e+01,&
       1.189080097391863e+01,&
      -3.770585380920011e-01,&
      -3.770585379486420e-01,&
       1.189080097108004e+01,&
      -8.069146560534483e+01,&
       3.122202458286709e+02,&
      -7.916434219820235e+02,&
       1.400159557890686e+03,&
      -1.729512527224754e+03,&
       1.385924948665997e+03,&
      -5.069710800063730e+02 /)

   ! general 36-point gauss-legendre quadrature rule for singular integrands
   ! containing one or more singular and/or hypersingular terms
   ! A*P_n(x) + B*P_n(x)*ln(|x-x0|) + C*P_n(x)/(x0-x) + D*P_n(x)/(x0-x)^2
   ! nodes
   double precision, parameter :: &
      xglg_36(36) = (/&
      -9.978304624840858e-01,&
      -9.885864789022123e-01,&
      -9.720276910496980e-01,&
      -9.482729843995076e-01,&
      -9.174977745156591e-01,&
      -8.799298008903971e-01,&
      -8.358471669924753e-01,&
      -7.855762301322066e-01,&
      -7.294891715935566e-01,&
      -6.680012365855210e-01,&
      -6.015676581359806e-01,&
      -5.306802859262452e-01,&
      -4.558639444334203e-01,&
      -3.776725471196892e-01,&
      -2.966849953440283e-01,&
      -2.135008923168656e-01,&
      -1.287361038093848e-01,&
      -4.301819847370861e-02,&
       4.301819847370861e-02,&
       1.287361038093848e-01,&
       2.135008923168656e-01,&
       2.966849953440283e-01,&
       3.776725471196892e-01,&
       4.558639444334203e-01,&
       5.306802859262452e-01,&
       6.015676581359806e-01,&
       6.680012365855210e-01,&
       7.294891715935566e-01,&
       7.855762301322066e-01,&
       8.358471669924753e-01,&
       8.799298008903971e-01,&
       9.174977745156591e-01,&
       9.482729843995076e-01,&
       9.720276910496980e-01,&
       9.885864789022123e-01,&
       9.978304624840858e-01 /)

    ! weights
    double precision, parameter :: &
       wglg_36(36) = (/&
       -2.370763065337607e+00,&
       -8.181667077901314e-02,&
        2.140169556187191e+00,&
        2.520675380752111e+00,&
        7.085661950016020e-01,&
       -1.819271801993064e+00,&
       -2.808185557676485e+00,&
       -1.118150194369488e+00,&
        1.974290056523884e+00,&
        3.424527802095672e+00,&
        1.141420230715159e+00,&
       -3.075465301896946e+00,&
       -3.763159219147007e+00,&
        2.186818189824050e+00,&
        6.420761841181085e+00,&
       -6.870579144433949e+00,&
        2.522531789704407e+00,&
       -1.323700863509908e-01,&
       -1.323700863457605e-01,&
        2.522531789647625e+00,&
       -6.870579144301985e+00,&
        6.420761841105238e+00,&
        2.186818189748193e+00,&
       -3.763159219113578e+00,&
       -3.075465301826242e+00,&
        1.141420230784354e+00,&
        3.424527801980273e+00,&
        1.974290056459363e+00,&
       -1.118150194342107e+00,&
       -2.808185557623396e+00,&
       -1.819271801932196e+00,&
        7.085661950133750e-01,&
        2.520675380642497e+00,&
        2.140169556134462e+00,&
       -8.181667067939005e-02,&
       -2.370763065351329e+00 /)


   ! standard 18-point gauss-legendre quadrature rule
   ! nodes
   double precision, parameter :: &
      xgls_18(18) = (/ &
      -9.915651684209309e-01,&
      -9.558239495713978e-01,&
      -8.926024664975557e-01,&
      -8.037049589725231e-01,&
      -6.916870430603532e-01,&
      -5.597708310739475e-01,&
      -4.117511614628426e-01,&
      -2.518862256915055e-01,&
      -8.477501304173529e-02,&
       8.477501304173529e-02,&
       2.518862256915055e-01,&
       4.117511614628426e-01,&
       5.597708310739475e-01,&
       6.916870430603532e-01,&
       8.037049589725231e-01,&
       8.926024664975557e-01,&
       9.558239495713978e-01,&
       9.915651684209309e-01 /)
   ! weights
   double precision, parameter :: &
      wgls_18(18) = (/ &
      2.161601352648413e-02,&
      4.971454889496922e-02,&
      7.642573025488925e-02,&
      1.009420441062870e-01,&
      1.225552067114784e-01,&
      1.406429146706506e-01,&
      1.546846751262652e-01,&
      1.642764837458327e-01,&
      1.691423829631436e-01,&
      1.691423829631436e-01,&
      1.642764837458327e-01,&
      1.546846751262652e-01,&
      1.406429146706506e-01,&
      1.225552067114784e-01,&
      1.009420441062870e-01,&
      7.642573025488925e-02,&
      4.971454889496922e-02,&
      2.161601352648413e-02/)

   ! standard 36-point gauss-legendre quadrature rule
   ! nodes
   double precision, parameter :: &
      xgls_36(36) = (/ &
      -9.97830462484085801e-01,&
      -9.88586478902212296e-01,&
      -9.72027691049697995e-01,&
      -9.48272984399507579e-01,&
      -9.17497774515659059e-01,&
      -8.79929800890397074e-01,&
      -8.35847166992475299e-01,&
      -7.85576230132206565e-01,&
      -7.29489171593556640e-01,&
      -6.68001236585521019e-01,&
      -6.01567658135980565e-01,&
      -5.30680285926245165e-01,&
      -4.55863944433420265e-01,&
      -3.77672547119689228e-01,&
      -2.96684995344028257e-01,&
      -2.13500892316865587e-01,&
      -1.28736103809384800e-01,&
      -4.30181984737086076e-02,&
       4.30181984737086076e-02,&
       1.28736103809384800e-01,&
       2.13500892316865587e-01,&
       2.96684995344028257e-01,&
       3.77672547119689228e-01,&
       4.55863944433420265e-01,&
       5.30680285926245165e-01,&
       6.01567658135980565e-01,&
       6.68001236585521019e-01,&
       7.29489171593556640e-01,&
       7.85576230132206565e-01,&
       8.35847166992475299e-01,&
       8.79929800890397074e-01,&
       9.17497774515659059e-01,&
       9.48272984399507579e-01,&
       9.72027691049697995e-01,&
       9.88586478902212296e-01,&
       9.97830462484085801e-01 /)

    ! weights
    double precision, parameter :: &
       wgls_36(36) = (/ &
       5.56571966424778374e-03,&
       1.29159472840641044e-02,&
       2.01815152977351739e-02,&
       2.72986214985683553e-02,&
       3.42138107703074817e-02,&
       4.08757509236452321e-02,&
       4.72350834902660471e-02,&
       5.32447139777596778e-02,&
       5.88601442453245485e-02,&
       6.40397973550154292e-02,&
       6.87453238357362967e-02,&
       7.29418850056530038e-02,&
       7.65984106458706265e-02,&
       7.96878289120715594e-02,&
       8.21872667043396510e-02,&
       8.40782189796617924e-02,&
       8.53466857393384987e-02,&
       8.59832756703946266e-02,&
       8.59832756703946266e-02,&
       8.53466857393384987e-02,&
       8.40782189796617924e-02,&
       8.21872667043396510e-02,&
       7.96878289120715594e-02,&
       7.65984106458706265e-02,&
       7.29418850056530038e-02,&
       6.87453238357362967e-02,&
       6.40397973550154292e-02,&
       5.88601442453245485e-02,&
       5.32447139777596778e-02,&
       4.72350834902660471e-02,&
       4.08757509236452321e-02,&
       3.42138107703074817e-02,&
       2.72986214985683553e-02,&
       2.01815152977351739e-02,&
       1.29159472840641044e-02,&
       5.56571966424778374e-03 /)

   public :: rw_bim_symm_init, rw_bim_symm_solve

contains

   subroutine rw_bim_symm_init(g_def,dump,test)
     implicit none
     
     type(grid_mg_def),pointer :: g_def
     logical :: dump,test

     integer :: i,j,k,nx,ny,nz,nyg,igr,my_rank_y
     real(8) :: xx,yy

     real(8),allocatable,dimension(:,:,:,:) :: Jn,Jng
     real(8),allocatable,dimension(:,:,:)   :: Bn,JJ,RR,ZZ,Bng,Rg,Zg
     real(8),allocatable,dimension(:)       :: phi,theta

     igr = 1

     !Select only mpi ranks at r=1 boundary
     if (g_def%ihi(igr) /= g_def%nxgl(igr)) return

     !Find local rank (for IO)
#if defined(petsc)
     call MPI_Comm_rank(g_def%MPI_COMM_Y,my_rank_y,mpierr)
#else
     my_rank_y = 0
#endif

     bim_dump = (my_rank_y==0.and.dump)

     if (test) then
        r0   = 1d0                       ! Minor radius
        Rmaj = g_def%params(1)           ! Major radius (from input deck)
        a    = sqrt((Rmaj+r0)*(Rmaj-r0)) ! Radius of focal ring of torus
        rho0 = acosh(Rmaj/r0)            ! Segura & Gil CPC 124(1) 2000, Pages 104-122, Eq. 19
        z0   = cosh(rho0)
     endif
     
     !Global angular mesh and associated qtys
     nx = g_def%nxv(igr) 
     ny = g_def%nyv(igr)
     nz = g_def%nzv(igr)

     nyg = g_def%nygl(igr)

     allocate(JJ(1,ny,1),Jn(1,ny,1,3),RR(1,ny,1),ZZ(1,ny,1))

     allocate(Rg(1,nyg,1),Zg(1,nyg,1),Jng(1,nyg,1,3),theta(nyg))

     if (test) allocate(Bn(1,ny,1),Bng(1,nyg,1))
     
     theta = g_def%yg(1:nyg)

     i = nx ; k = 1
     do j=1,ny
        !Cell position at (nx+1/2) boundary in (R,Z) coords
        ZZ(1,j,1) = 0.5d0*(g_def%gmetric%grid(igr)%car(i  ,j,k,3) &
                          +g_def%gmetric%grid(igr)%car(i+1,j,k,3))
        xx = 0.5d0*(g_def%gmetric%grid(igr)%car(i  ,j,k,1)        &
                   +g_def%gmetric%grid(igr)%car(i+1,j,k,1))
        yy = 0.5d0*(g_def%gmetric%grid(igr)%car(i  ,j,k,2)        &
                   +g_def%gmetric%grid(igr)%car(i+1,j,k,2))
        RR(1,j,1) = sqrt(xx*xx + yy*yy)

        !Cartesian J.n (n->normal vector pointing *inward*)
        ! at (nx+1/2) bdry at phi=0 poloidal plane
        Jn(1,j,1,:) =                                           &
            -0.25d0*(g_def%gmetric%grid(igr)%cov(i  ,j,k  ,1,:) &
                    *g_def%gmetric%grid(igr)%jac(i  ,j,k)       &
                    +g_def%gmetric%grid(igr)%cov(i+1,j,k  ,1,:) &
                    *g_def%gmetric%grid(igr)%jac(i+1,j,k)       &
                    +g_def%gmetric%grid(igr)%cov(i  ,j,k-1,1,:) &
                    *g_def%gmetric%grid(igr)%jac(i  ,j,k-1)     &
                    +g_def%gmetric%grid(igr)%cov(i+1,j,k-1,1,:) &
                    *g_def%gmetric%grid(igr)%jac(i+1,j,k-1))

        !Jacobian interpolated at (nx+1/2) bdry
        JJ(1,j,1) = 0.5d0*(g_def%gmetric%grid(igr)%jac(i  ,j,k) &
                          +g_def%gmetric%grid(igr)%jac(i+1,j,k))

        if (test) Bn(1,j,1) = -JJ(1,j,1)*get_Bn_analytic(RR(1,j,1),ZZ(1,j,1))
     enddo

      !Gather qtys in parallel along boundary
#if defined(petsc)
      if (test) call find_global_nobc(Bn,Bng,mpi_comm=g_def%MPI_COMM_Y)
      call find_global_nobc(RR,Rg ,mpi_comm=g_def%MPI_COMM_Y)
      call find_global_nobc(ZZ,Zg ,mpi_comm=g_def%MPI_COMM_Y)
      call find_global_nobc(Jn,Jng,mpi_comm=g_def%MPI_COMM_Y)
#else
      if (test) call find_global_nobc(Bn,Bng)
      call find_global_nobc(RR,Rg )
      call find_global_nobc(ZZ,Zg )
      call find_global_nobc(Jn,Jng)
#endif

     if (bim_dump) write(*,*) " BIM: setting BCs"

     if (bim_dump) write(*,*) " BIM: allocating vars"
     call setup_bim_data(bd,nyg,theta,Rg(1,:,1),Zg(1,:,1),Jng(1,:,1,:))

     if (bim_dump) write(*,*) " BIM: setting BCs"
     if (test) bd%Bn(1:bd%n) = Bn(1,:,1)
     call set_BCs(bd)

     if (bim_dump) write(*,*) " BIM: computing matrix"
     call compute_KinvL_matrix(bd,test)

     !Compute analytical test
     if (test) then
        if (bim_dump) then
           write (*,*) "BIM DIAG ","theta ","Bng ","Rg ","Zg ","|J.ng|"
           do i=1,nyg
              write (*,*) "BIM DIAG",theta(i),Bng(1,i,1),Rg(1,i,1),Zg(1,i,1) &
                                    ,sqrt(sum(Jng(1,i,1,:)**2))
           enddo
        endif

        bd%Bn(1:bd%n) = Bng(1,:,1)
        call rw_bim_symm_solve(Bng(1,:,1),bd%phi(1:bd%n))

        if (bim_dump) then
           write(*,*) "BIM: computing error"
           call compute_error(bd,test)
        endif
     endif
     
     deallocate(JJ,Jn,RR,ZZ,Jng,Rg,Zg)

     if (test) deallocate(Bn,Bng)
     
   end subroutine rw_bim_symm_init

   subroutine rw_bim_symm_solve(Bn,phi)
      implicit none
      ! pass
        
      double precision, intent(in), dimension(:) :: Bn
      double precision, intent(out), dimension(size(Bn)) :: phi
      ! local

      bd%Bn(1:bd%n) = Bn

      if (bim_dump) write(*,*) "BIM: solving system"
      bd%phi(1:bd%n) = matmul(bd%Kinv,matmul(bd%L,bd%Bn(1:bd%n)))

      phi = bd%phi(1:bd%n)

   end subroutine rw_bim_symm_solve

   subroutine compute_error(bd,test)
      use toroidal_harmonics, only: p0, p1, p2
      implicit none
      ! pass
      type(bim_data_t) :: bd
      logical, intent(in) :: test
      integer :: i, fh
      double precision :: ct, th, jacobian, hrho, exact, error

      error = 0.d0
      open(newunit=fh,file="rw_bim_solution.txt",action="write")
      do i = 1, bd%n
         th = toroidal_theta(bd%R(i),bd%Z(i))
         ct = cos(th)
         jacobian = abs(a**3*sinh(rho0)/(ct-z0)**3)
         hrho = a/(z0-cos(th))

         if (test) then
            ! m = 0,n=0 mode
            !exact = sqrt(z0-ct)*p0(z0)

            ! m = 1,n=0 mode
            exact = sqrt(z0-ct)*p1(z0)*exp(dcmplx(0,th))

            error = error + abs(exact - bd%phi(i))
            write(fh,"(10(es25.15e3))") bd%R(i), bd%Z(i), bd%Jnorm(i,1), bd%Jnorm(i,3), &
               jacobian, hrho, th, bd%Bn(i), bd%phi(i), exact
         else
            write(fh,"(9(es25.15e3))") bd%R(i), bd%Z(i), bd%Jnorm(i,1), bd%Jnorm(i,3), &
               jacobian, hrho, th, bd%Bn(i), bd%phi(i)
         end if
      end do
      close(fh)

      if (test) then
         write(*,*) "BIM error:",bd%n, error/bd%n
         open(newunit=fh,file="error.txt")
         write(fh,*) bd%n, error/bd%n
         close(fh)
      end if

   end subroutine compute_error

!!$   subroutine set_Bn_initial_condition(bd)
!!$      implicit none
!!$      ! pass
!!$      type(bim_data_t) :: bd
!!$      ! local
!!$      integer :: i
!!$
!!$      do i = 1,bd%n
!!$
!!$         bd%Bn(i) = get_Bn_analytic(bd%R(i),bd%Z(i))
!!$         
!!$      end do
!!$
!!$   end subroutine set_Bn_initial_condition

   function evaldRdrho(th) result(dRdrho)
      implicit none
      double precision, intent(in) :: th
      double precision :: dRdrho

      dRdrho = a*(1-cos(th)*z0)/(cos(th)-z0)**2

   end function evaldRdrho

   function evaldZdrho(th) result(dZdrho)
      implicit none
      double precision, intent(in) :: th
      double precision :: dZdrho

      dZdrho = -a*sin(th)*sinh(rho0)/(z0-cos(th))**2

   end function evaldZdrho

   function evaldRdtheta(th) result(dRdtheta)
      implicit none
      double precision, intent(in) :: th
      double precision :: dRdtheta

      dRdtheta = -a*sinh(rho0)*sin(th)/(z0-cos(th))**2

   end function evaldRdtheta

   function evaldZdtheta(th) result(dZdtheta)
      implicit none
      double precision, intent(in) :: th
      double precision :: dZdtheta

      dZdtheta = a*(z0*cos(th)-1)/(z0-cos(th))**2

   end function evaldZdtheta

   function get_Bn_analytic(Rin,Z) result(b)
      use toroidal_harmonics, only: p0, p1, p2
      implicit none
      double precision, intent(in) :: Rin, Z
      double precision :: R, b, ct, th, hrho, jacobian, n(2), norm

      R = Rin
      th = toroidal_theta(R, Z)
      ct = cos(th)
      jacobian = abs(a**3*sinh(rho0)/(ct-z0)**3)
      hrho = a/(z0-cos(th))

      ! m = 0,n=0 mode
      !b = -sqrt(z0-ct)/(2*a*sinh(rho0))*&
         !(  (z0*ct-1)*p0(z0) + (z0-ct)*p1(z0) )

      ! m = 1,n=0 mode
      b = -exp(dcmplx(0,th))*sqrt(z0-ct)/(2*a*sinh(rho0))*&
      (  (3*z0*ct-2*z0**2-1)*p1(z0) + 3*(z0-ct)*p2(z0) )

      ! unit e_n.B = B^n/(J*|n|) where J is jacobian and |n| is magnitude of normal vector
      !b = b/(jacobian*norm)
      !b = b/jacobian

   end function get_Bn_analytic

   !> compute toroidal polar coordinate theta (also called sigma)
   !> from cylindrical R and Z assuming phi = 0
   function toroidal_theta(R,Z) result(theta)
      implicit none
      double precision, intent(in) :: R,Z
      double precision :: theta
      
      theta = acos(z0 - a*sinh(rho0)/R)
      theta = merge(theta, 2*pi-theta, Z < 0.d0)
   end function toroidal_theta

!!$   subroutine correct_allR(R)
!!$      implicit none
!!$      double precision, intent(inout) :: R(:)
!!$      integer :: i
!!$
!!$      do i = 1, size(R)
!!$         R(i) = correct_one_R(R(i))
!!$      end do
!!$   end subroutine correct_allR
!!$
!!$   function correct_one_R(Rin) result(Rout)
!!$      implicit none
!!$      double precision, intent(in) :: Rin
!!$      double precision :: Rout
!!$      Rout = Rin-1.65d0+Rmaj
!!$   end function correct_one_R


   subroutine setup_bim_data(bd,n,xi,R,Z,Jnorm)
      implicit none
      type(bim_data_t) :: bd
      ! pass
      integer :: n
      double precision, dimension(n), intent(in) :: xi, R, Z
      double precision, dimension(n,3), intent(in) :: Jnorm
      ! local
      integer :: i
      double precision :: dxlast

      bd%n = n
      allocate(bd%xi(0:n+1))
      bd%xi(1:n) = xi
      allocate(bd%R(0:n+1))
      bd%R(1:n) = R
      allocate(bd%Z(0:n+1))
      bd%Z(1:n) = Z
      allocate(bd%Bn(0:n+1))
      bd%Bn = 0d0
      allocate(bd%phi(0:n+1))
      bd%phi = 0d0
      allocate(bd%Jnorm(0:n+1,3))
      bd%Jnorm(1:n,:) = Jnorm
      ! matrices
      allocate(bd%L(n,n))
      bd%L = 0d0
      allocate(bd%K(n,n))
      bd%K = 0d0
      allocate(bd%Kinv(n,n))
      bd%Kinv = 0d0
      allocate(bd%KinvL(n,n))
      bd%KinvL = 0d0

      ! compute mesh spacing
      bd%dxi = xi(2) - xi(1)

   end subroutine setup_bim_data

   subroutine set_BCs(bd)
      implicit none
      type(bim_data_t) :: bd
      integer :: n
      n = bd%n
      ! periodic boundaries
      ! lo side
      bd%xi(0) = bd%xi(1) - bd%dxi
      bd%R(0) = bd%R(n)
      bd%Z(0) = bd%Z(n)
      bd%Bn(0) = bd%Bn(n)
      bd%Jnorm(0,:) = bd%Jnorm(n,:)
      ! hi side
      bd%xi(bd%n+1) = bd%xi(n) + bd%dxi
      bd%R(n+1) = bd%R(1)
      bd%Z(n+1) = bd%Z(1)
      bd%Bn(n+1) = bd%Bn(1)
      bd%Jnorm(n+1,:) = bd%Jnorm(1,:)

   end subroutine set_BCs

   !--------------------------------------------------------
   ! MATRIX CONSTRUCTION
   !--------------------------------------------------------

   subroutine compute_KinvL_matrix(bd,test)
      implicit none
      type(bim_data_t) :: bd
      logical, intent(in) :: test

      call compute_L_matrix(bd,test)
      call compute_K_matrix(bd,test)
      call compute_Kinv_matrix(bd)

      bd%KinvL = matmul(bd%Kinv,bd%L)

   end subroutine compute_KinvL_matrix

   subroutine compute_Kinv_matrix(bd)
      implicit none
      type(bim_data_t) :: bd
      double precision :: work(bd%n)
      integer :: ipiv(bd%n), ierr

      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      bd%Kinv = bd%K

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(bd%n, bd%n, bd%Kinv, bd%n, ipiv, ierr)

      if (ierr /= 0) then
         stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(bd%n, bd%Kinv, bd%n, ipiv, work, bd%n, ierr)

   end subroutine compute_Kinv_matrix

   subroutine compute_L_matrix(bd,test)
      implicit none
      type(bim_data_t) :: bd
      logical, intent(in) :: test
      integer :: i, j

      do i = 1, bd%n
         do j = 1, bd%n
            bd%L(i,j) = L_matrix_element(bd,i,j,test)
         end do
      end do

   end subroutine compute_L_matrix

   subroutine compute_K_matrix(bd,test)
      implicit none
      type(bim_data_t) :: bd
      logical, intent(in) :: test
      integer :: i, j

      do i = 1, bd%n
         do j = 1, bd%n
            bd%K(i,j) = K_matrix_element(bd,i,j,test)
         end do
      end do

   end subroutine compute_K_matrix

   function L_matrix_element(bd,i,j,test) result(Lij)
      implicit none
      type(bim_data_t) :: bd
      integer, intent(in) :: i, j
      logical, intent(in) :: test
      integer :: flavor
      double precision :: Lij, delta

      flavor = merge(SINGULAR,REGULAR,i==j)
      Lij = -integrate_toroidal(bd, i, j, flavor, integrand_g, test)

   end function L_matrix_element

   function K_matrix_element(bd,i,j,test) result(Kij)
      implicit none
      type(bim_data_t) :: bd
      integer, intent(in) :: i, j
      logical, intent(in) :: test
      integer :: flavor
      double precision :: Kij, delta

      flavor = merge(SINGULAR,REGULAR,i==j)
      Kij = -integrate_toroidal(bd, i, j, flavor, integrand_dgdn, test)

      if (i==j) Kij = Kij + 0.5d0

   end function K_matrix_element

   ! TOROIDAL TORSOINAL GREEN'S FUNCTION AND DERIVATIVE

   function evalchi(Ri,Zi,Rj,Zj) result(result)
      implicit none
      double precision, intent(in) :: Ri, Zi, Rj, Zj
      double precision :: result

      result = (Ri*Ri + Rj*Rj + (Zi-Zj)**2)/(2*Ri*Rj)
   end function evalchi

   function evaldchidR(R,Z,RP,ZP) result(result)
      implicit none
      double precision, intent(in) :: R, Z, RP, ZP
      double precision :: result

      result = 1/RP - (R*R + RP*RP + (Z-ZP)**2)/(2*R**2*RP)
   end function evaldchidR

   function evaldchidZ(R,Z,RP,ZP) result(result)
      implicit none
      double precision, intent(in) :: R, Z, RP, ZP
      double precision :: result

      result = (Z-ZP)/(R*RP)
   end function evaldchidZ

   function evaldchidRj(Ri,Zi,Rj,Zj) result(result)
      implicit none
      double precision, intent(in) :: Ri, Zi, Rj, Zj
      double precision :: result

      result = 1/Ri - (Ri*Ri + Rj*Rj + (Zi-Zj)**2)/(2*Ri*Rj**2)
   end function evaldchidRj

   function evaldchidZj(Ri,Zi,Rj,Zj) result(result)
      implicit none
      double precision, intent(in) :: Ri, Zi, Rj, Zj
      double precision :: result

      result = (Zj-Zi)/(Ri*Rj)
   end function evaldchidZj

   function evalmu(chi) result(result)
      implicit none
      double precision :: chi, result

      result = sqrt(2.d0/(chi+1.d0))
   end function evalmu

   function Kelliptic(mu) result(result)
      implicit none
      double precision, intent(in) :: mu
      double precision :: mu2, mymu, result
      mu2 = mu*mu
      mymu = 1.d0 - mu2
      result = ceik(mymu)
      !result = ellipticK(mymu)
   end function Kelliptic

   function Eelliptic(mu) result(result)
      implicit none
      double precision, intent(in) :: mu
      double precision :: mu2, mymu, result
      mu2 = mu*mu
      mymu = 1.d0-mu2
      result = ceie(mymu)
      !result = ellipticE(mymu)
   end function Eelliptic

   function evaldKdmu(mu) result(result)
      implicit none
      double precision :: mu, result
      result = Eelliptic(mu)/(mu*(1.d0-mu**2)) - Kelliptic(mu)/mu
   end function evaldKdmu

   function greens(Ri,Zi,Rj,Zj) result(G)
      implicit none
      double precision, intent(in) :: Ri, Zi, Rj, Zj
      double precision :: chi, mu, K, G

      chi = evalchi(Ri,Zi,Rj,Zj)
      mu = evalmu(chi)
      K = Kelliptic(mu)

      G = -mu*K/(2.d0*pi*sqrt(Ri*Rj))
   end function greens

   ! derivative wrt Rj and evaluated at (Ri,Zi) and (Rj,Zj)
   function dgreensdRj(Ri,Zi,Rj,Zj) result(dGdRj)
      implicit none
      double precision, intent(in) :: Ri, Zi, Rj, Zj
      double precision :: chi, mu, K, dGdRj, dchidRj

      chi = evalchi(Ri,Zi,Rj,Zj)
      mu = evalmu(chi)
      K = Kelliptic(mu)
      dchidRj = evaldchidRj(Ri,Zi,Rj,Zj)

      dGdRj = -greens(Ri,Zi,Rj,Zj)/(2*Rj) + mu**3/(8*pi*sqrt(Ri*Rj))*&
         Eelliptic(mu)/(1-mu**2)*dchidRj

   end function dgreensdRj

   ! derivative wrt Zj and evaluated at (Ri,Zi) and (Rj,Zj)
   function dgreensdZj(Ri,Zi,Rj,Zj) result(dGdZj)
      implicit none
      double precision, intent(in) :: Ri, Zi, Rj, Zj
      double precision :: chi, mu, K, dGdZj, dchidZj

      chi = evalchi(Ri,Zi,Rj,Zj)
      mu = evalmu(chi)
      K = Kelliptic(mu)
      dchidZj = evaldchidZj(Ri,Zi,Rj,Zj)

      dGdZj = mu**3/(8*pi*sqrt(Ri*Rj))*Eelliptic(mu)/(1-mu**2)*dchidZj

   end function dgreensdZj

   ! PARABOLIC INTERPOLANT
   function interpolate(x,q,x_want) result(q_want)
      implicit none
      double precision, intent(in) :: x(3), q(3), x_want
      double precision :: q_want, dx, dqdx, d2qdx2, a, b, dp, dm, dc
      integer :: i = 2
      dx = x(i) - x(i-1)

      dp = q(i+1)-q(i)
      dm = q(i)-q(i-1)
      dc = 0.5d0*(q(i+1)-q(i-1))

      dqdx = dc/dx
      d2qdx2 = (q(i+1) - 2*q(i) + q(i-1))/dx**2

      ! taylor series expansion
      q_want = q(i) + dqdx*(x_want-x(i)) + 0.5d0*d2qdx2*(x_want-x(i))**2

   end function interpolate

   ! INTEGRANDS

   function integrand_g(bd,i,j,xi_j) result(integrand)
      implicit none
      ! pass
      type(bim_data_t), intent(in) :: bd
      integer, intent(in) :: i, j
      double precision, intent(in) :: xi_j
      ! internal
      double precision :: integrand
      double precision :: Ri, Zi, Rj, Zj

      ! interpolation
      Rj = interpolate(bd%xi(j-1:j+1),bd%R(j-1:j+1),xi_j)
      Zj = interpolate(bd%xi(j-1:j+1),bd%Z(j-1:j+1),xi_j)

      Ri = bd%R(i)
      Zi = bd%Z(i)

      integrand = greens(Ri,Zi,Rj,Zj)

   end function integrand_g

   function integrand_dgdn(bd,i,j,xi_j) result(integrand)
      implicit none
      ! pass
      type(bim_data_t), intent(in) :: bd
      integer, intent(in) :: i, j
      double precision, intent(in) :: xi_j
      ! internal
      double precision :: integrand

      double precision :: normlo(3), normhi(3), norm(3), Ri, Zi, Rj, Zj
      double precision :: Jnorm(3)
      integer :: k

      ! interpolation
      Rj = interpolate(bd%xi(j-1:j+1),bd%R(j-1:j+1),xi_j)
      Zj = interpolate(bd%xi(j-1:j+1),bd%Z(j-1:j+1),xi_j)

      ! interpolation
      do k = 1, 3
         Jnorm(k) = interpolate(bd%xi(j-1:j+1),bd%Jnorm(j-1:j+1,k),xi_j)
      end do

      Ri = bd%R(i)
      Zi = bd%Z(i)

      integrand = (Jnorm(1)*dgreensdRj(Ri,Zi,Rj,Zj) + Jnorm(3)*dgreensdZj(Ri,Zi,Rj,Zj))

   end function integrand_dgdn

   ! SPECIAL GAUSSIAN QUADRATURE

   function integrate_toroidal(bd, i, j, flavor, integrand, test) result(result)
      implicit none
      ! pass
      type(bim_data_t), intent(in) :: bd
      integer, intent(in) :: i, j, flavor
      double precision, external :: integrand
      logical, intent(in) :: test
      ! internal
      double precision :: xi_j, t0, g1, result, error
      double precision :: x18(18), w18(18)
      double precision :: x36(36), w36(36)
      double precision :: result18, result36
      integer :: k

      ! compute integral transform
      ! integral transform from (lima,limb) to (-1,1)
      t0 = bd%xi(j)
      g1 = 2.d0/(bd%dxi)

      ! select appropriate nodes and weights
      select case(flavor)
      case(REGULAR)
         x18 = xgls_18; w18 = wgls_18
         x36 = xgls_36; w36 = wgls_36
      case(SINGULAR)
         x18 = xglg_18; w18 = wglg_18
         x36 = xglg_36; w36 = wglg_36
      end select

      ! 18-point quadrature
      result18 = 0.d0
      do k = 1, size(x18)
         xi_j = t0 + x18(k)/g1
         result18 = result18 + w18(k) * integrand(bd, i, j, xi_j) / g1
      end do

      ! 36-point quadrature
      result36 = 0.d0
      do k = 1, size(x36)
         xi_j = t0 + x36(k)/g1
         result36 = result36 + w36(k) * integrand(bd, i, j, xi_j) / g1
      end do

      ! compute error and return result of 36-point quadrature
      result = result36

      if (test) then
         error = (result18-result36)/result36
         if (error > 1.d-3) then
            print *, "error large in quadrature rule:"
            print *, result18, result36, error
         end if
      end if

   end function integrate_toroidal

end module rw_bim_mod
