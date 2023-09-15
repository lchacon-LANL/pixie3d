module rw_bim_mod
   use math, only: pi, ellipticK, ellipticE
   use elliptic, only: ceik, ceie
   implicit none
   private

   integer, parameter :: &
      REGULAR = 0, &
      SINGULAR = REGULAR + 1

   double precision, target :: &
      xglg_18(18), wglg_18(18), &
      xglg_36(36), wglg_36(36), &
      xgls_18(18), wgls_18(18), &
      xgls_36(36), wgls_36(36)

   type :: bim_data_t
      integer :: n
      double precision :: dxi
      double precision, allocatable, dimension(:) :: R, Z, xi, Bn, phi
      double precision, allocatable, dimension(:,:) :: Jnorm, L, K, Kinv, KinvL
   end type bim_data_t

   ! for verification tests
   double precision, parameter :: r0 = 1.d0 ! actual minor radius we want
   !double precision, parameter :: Rmaj = 5.d2 ! actual major radius we want
   !double precision, parameter :: Rmaj = 1.d1 ! actual major radius we want
   double precision, parameter :: Rmaj = 1.65d0 ! actual major radius we want
   ! radius of focal ring of torus
   double precision, parameter :: a = sqrt((Rmaj+r0)*(Rmaj-r0))
   ! Segura & Gil CPC Volume 124, Issue 1, 15 January 2000, Pages 104-122, Eq. 19
   double precision, parameter :: rho0 = acosh(Rmaj/r0) ! toroidal "radial inverse" coordinate for rho0
   double precision, parameter :: z0 = cosh(rho0)

   public :: rw_bim_mod_init, solve_toroidally_symmetric_bim_RZ, get_Bn_analytic

contains

   subroutine rw_bim_mod_init()
      implicit none
      call gauss_init()
   end subroutine rw_bim_mod_init

   subroutine solve_toroidally_symmetric_bim_RZ(n,xi,R,Z,Jnorm,Bn,phi)
      implicit none
      ! pass
      integer, intent(in) :: n
      double precision, intent(in), dimension(n) :: xi, R, Z, Bn
      double precision, intent(in), dimension(n,3) :: Jnorm
      double precision, intent(out), dimension(n) :: phi
      ! local
      integer :: i
      type(bim_data_t) :: bd

      call rw_bim_mod_init()

      print *, "bim allocating"
      call allocate_and_store_bim_data(bd,n,xi,R,Z,Bn,Jnorm,phi)

      !call correct_allR(bd%R)

      print *, "setting BCs"
      call set_BCs(bd)

      print *, "bim computing matrix"
      call compute_KinvL_matrix(bd)

      !print *, "bim setting initial Bn condition"
      !call set_Bn_initial_condition(bd)

      print *, "bim solving system"
      bd%phi(1:bd%n) = matmul(bd%Kinv,matmul(bd%L,bd%Bn(1:bd%n)))

      print *, "bim computing error"
      call compute_error(bd)

      phi = bd%phi(1:bd%n)

   end subroutine solve_toroidally_symmetric_bim_RZ

   subroutine compute_error(bd)
      use toroidal_harmonics, only: p0, p1, p2
      implicit none
      ! pass
      type(bim_data_t) :: bd
      integer :: i, fh
      double precision :: ct, th, jacobian, hrho, exact, error

      error = 0.d0
      open(newunit=fh,file="solution.txt",action="write")
      do i = 1, bd%n
         th = toroidal_theta(bd%R(i),bd%Z(i))
         ct = cos(th)
         jacobian = abs(a**3*sinh(rho0)/(ct-z0)**3)
         hrho = a/(z0-cos(th))

         ! m = 0,n=0 mode
         !exact = sqrt(z0-ct)*p0(z0)

         ! m = 1,n=0 mode
         exact = sqrt(z0-ct)*p1(z0)*exp(dcmplx(0,th))

         error = error + abs(exact - bd%phi(i))
         write(fh,"(10(es25.15e3))") bd%R(i), bd%Z(i), bd%Jnorm(i,1), bd%Jnorm(i,3), jacobian, hrho, th, bd%Bn(i), bd%phi(i), exact
      end do
      close(fh)

      open(newunit=fh,file="error.txt")
      write(fh,*) bd%n, error/bd%n
      close(fh)

   end subroutine compute_error

   subroutine set_Bn_initial_condition(bd)
      implicit none
      ! pass
      type(bim_data_t) :: bd
      ! local
      integer :: i

      do i = 1,bd%n

         bd%Bn(i) = get_Bn_analytic(bd%R(i),bd%Z(i))
         
      end do

   end subroutine set_Bn_initial_condition

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
      !theta = acos(z0 - a*sinh(rho0)/(R-1.65d0+Rmaj))
      theta = merge(theta, 2*pi-theta, Z < 0.d0)
   end function toroidal_theta

   subroutine correct_allR(R)
      implicit none
      double precision, intent(inout) :: R(:)
      integer :: i

      do i = 1, size(R)
         R(i) = correct_one_R(R(i))
      end do
   end subroutine correct_allR

   function correct_one_R(Rin) result(Rout)
      implicit none
      double precision, intent(in) :: Rin
      double precision :: Rout
      Rout = Rin-1.65d0+Rmaj
   end function correct_one_R


   subroutine allocate_and_store_bim_data(bd,n,xi,R,Z,Bn,Jnorm,phi)
      implicit none
      type(bim_data_t) :: bd
      ! pass
      integer :: n
      double precision, dimension(n), intent(in) :: xi, R, Z, Bn
      double precision, dimension(n,3), intent(in) :: Jnorm
      double precision, dimension(n), intent(in) :: phi
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
      bd%Bn(1:n) = Bn
      allocate(bd%phi(0:n+1))
      allocate(bd%Jnorm(0:n+1,3))
      bd%Jnorm(1:n,:) = Jnorm
      ! matrices
      allocate(bd%L(n,n))
      allocate(bd%K(n,n))
      allocate(bd%Kinv(n,n))
      allocate(bd%KinvL(n,n))

      ! compute mesh spacing
      bd%dxi = xi(2) - xi(1)

   end subroutine allocate_and_store_bim_data

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

   subroutine compute_KinvL_matrix(bd)
      implicit none
      type(bim_data_t) :: bd

      call compute_L_matrix(bd)
      call compute_K_matrix(bd)
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

   subroutine compute_L_matrix(bd)
      implicit none
      type(bim_data_t) :: bd
      integer :: i, j

      do i = 1, bd%n
         do j = 1, bd%n
            bd%L(i,j) = L_matrix_element(bd,i,j)
         end do
      end do

   end subroutine compute_L_matrix

   subroutine compute_K_matrix(bd)
      implicit none
      type(bim_data_t) :: bd
      integer :: i, j

      do i = 1, bd%n
         do j = 1, bd%n
            bd%K(i,j) = K_matrix_element(bd,i,j)
         end do
      end do

   end subroutine compute_K_matrix

   function L_matrix_element(bd,i,j) result(Lij)
      implicit none
      type(bim_data_t) :: bd
      integer, intent(in) :: i, j
      integer :: flavor
      double precision :: Lij, delta

      flavor = merge(SINGULAR,REGULAR,i==j)
      Lij = -integrate_toroidal(bd, i, j, flavor, integrand_g)

   end function L_matrix_element

   function K_matrix_element(bd,i,j) result(Kij)
      implicit none
      type(bim_data_t) :: bd
      integer, intent(in) :: i, j
      integer :: flavor
      double precision :: Kij, delta

      flavor = merge(SINGULAR,REGULAR,i==j)
      Kij = -integrate_toroidal(bd, i, j, flavor, integrand_dgdn)

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

   function integrate_toroidal(bd, i, j, flavor, integrand) result(result)
      implicit none
      ! pass
      type(bim_data_t), intent(in) :: bd
      integer, intent(in) :: i, j, flavor
      double precision, external :: integrand
      ! internal
      double precision :: xi_j, t0, g1, result18, result36, result, error
      double precision, pointer :: x(:), w(:)
      integer :: k

      ! compute integral transform
      ! integral transform from (lima,limb) to (-1,1)
      t0 = bd%xi(j)
      g1 = 2.d0/(bd%dxi)

      select case(flavor)
      case(REGULAR)
         x=>xgls_18; w=>wgls_18
      case(SINGULAR)
         x=>xglg_18; w=>wglg_18
      end select

      result18 = 0.d0
      do k = 1, size(x)
         xi_j = t0 + x(k)/g1
         result18 = result18 + w(k) * integrand(bd, i, j, xi_j) / g1
      end do
      nullify(x); nullify(w)

      select case(flavor)
      case(REGULAR)
         x=>xgls_36; w=>wgls_36
      case(SINGULAR)
         x=>xglg_36; w=>wglg_36
      end select

      result36 = 0.d0
      do k = 1, size(x)
         xi_j = t0 + x(k)/g1
         result36 = result36 + w(k) * integrand(bd, i, j, xi_j) / g1
      end do
      nullify(x); nullify(w)

      error = abs((result18-result36)/result36)
      if (error > 1.d-3) then
         print *, "error large in quadrature rule:", flavor
         !print *, result18, result36, error
         write(*,"(a8,es12.3,a6,i4,a6,i4)") &
            "error = ", error, "  i = ", i, "  j = ", j
      end if

      result = result36

   end function integrate_toroidal

   subroutine gauss_init()
      implicit none

      call read_quad(xglg_18,wglg_18,"../../../src/gl_general_18.txt")
      call read_quad(xglg_36,wglg_36,"../../../src/gl_general_36.txt")
      call read_quad(xgls_18,wgls_18,"../../../src/gl_standard_18.txt")
      call read_quad(xgls_36,wgls_36,"../../../src/gl_standard_36.txt")

   end subroutine gauss_init

   subroutine read_quad(x,w,fn)
      implicit none
      double precision :: x(:), w(:)
      character(len=*) :: fn
      integer :: i, fh
      open(file=fn,newunit=fh,action="read")
         do i = 1, size(x)
            read(fh,*) x(i), w(i)
         end do
      close(fh)
   end subroutine read_quad

end module rw_bim_mod
