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
      double precision, allocatable, dimension(:) :: R, Z, xi, J, Bn, phi, nmag
      double precision, allocatable, dimension(:,:) :: norm, L, K, Kinv, KinvL
   end type bim_data_t

   ! for verification tests
   double precision, parameter :: r0 = 1.d0 ! actual minor radius we want
   !double precision, parameter :: Rmaj = 5.d2 ! actual major radius we want
   double precision, parameter :: Rmaj = 1.d1 ! actual major radius we want
   double precision, parameter :: a = sqrt((Rmaj+r0)*(Rmaj-r0))
   ! Segura & Gil CPC Volume 124, Issue 1, 15 January 2000, Pages 104-122, Eq. 19
   double precision, parameter :: rho0 = acosh(Rmaj/r0) ! toroidal "radial inverse" coordinate for rho0
   double precision, parameter :: z0 = cosh(rho0)

   public :: rw_bim_mod_init, solve_toroidally_symmetric_bim_RZ

contains

   subroutine rw_bim_mod_init()
      implicit none
      call gauss_init()
   end subroutine rw_bim_mod_init

   subroutine solve_toroidally_symmetric_bim_RZ(n,xi,R,Z,J,norm,Bn,phi)
      implicit none
      ! pass
      integer, intent(in) :: n
      double precision, intent(in), dimension(n) :: xi, R, Z, J, Bn
      double precision, intent(in), dimension(n,3) :: norm
      double precision, intent(out), dimension(n) :: phi
      ! local
      integer :: i
      type(bim_data_t) :: bd

      call rw_bim_mod_init()

      print *, "bim allocating"
      call allocate_and_store_bim_data(bd,n,xi,R,Z,J,Bn,norm,phi)

      print *, "setting BCs"
      call set_BCs(bd)

      print *, "bim computing matrix"
      call compute_KinvL_matrix(bd)

      print *, "bim setting initial Bn condition"
      call set_Bn_initial_condition(bd)

      print *, "bim solving system"
      bd%phi(1:bd%n) = matmul(bd%KinvL,bd%Bn(1:bd%n))

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
      double precision :: ct, th, exact, error

      error = 0.d0
      open(newunit=fh,file="solution.txt",action="write")
      do i = 1, bd%n
         th = toroidal_theta(bd%R(i),bd%Z(i))
         ct = cos(th)

         ! m = 0,n=0 mode
         exact = sqrt(z0-ct)*p0(z0)

         ! m = 1,n=0 mode
         !exact = sqrt(z0-ct)*p1(z0)*exp(dcmplx(0,th))

         error = error + abs(exact - bd%phi(i))
         write(fh,"(6(es25.15e3))") bd%R(i), bd%Z(i), th, bd%Bn(i), bd%phi(i), exact
      end do
      close(fh)

      open(newunit=fh,file="error.txt")
      write(fh,*) bd%n, error/bd%n
      close(fh)

   end subroutine compute_error

   subroutine set_Bn_initial_condition(bd)
      use toroidal_harmonics, only: p0, p1, p2
      implicit none
      ! pass
      type(bim_data_t) :: bd
      ! local
      integer :: i
      double precision :: ct, th

      do i = 1,bd%n

         th = toroidal_theta(bd%R(i),bd%Z(i))
         
         ! m = 0,n=0 mode
         ct = cos(th)
         bd%Bn(i) = -sqrt(z0-ct)/(2*a*sinh(rho0))*&
            (  (z0*ct-1)*p0(z0) + (z0-ct)*p1(z0) )

         ! m = 1,n=0 mode
         !ct = cos(thetac(i))
         !bd%Bn(i) = -exp(dcmplx(0,th)*sqrt(z0-ct)/(2*a*sinh(rho0))*&
            !(  (3*z0*ct-2*z0**2-1)*p1(z0) + 3*(z0-ct)*p2(z0) )


      end do

   end subroutine set_Bn_initial_condition

   !> compute toroidal polar coordinate theta (also called sigma)
   !> from cylindrical R and Z assuming phi = 0
   function toroidal_theta(R,Z) result(theta)
      implicit none
      double precision, intent(in) :: R,Z
      double precision :: theta
      
      !theta = acos(z0 - a*sinh(rho0)/R)
      theta = acos(z0 - a*sinh(rho0)/(R-1.65d0+Rmaj))
      theta = merge(theta, 2*pi-theta, Z < 0.d0)
   end function toroidal_theta

   subroutine allocate_and_store_bim_data(bd,n,xi,R,Z,J,Bn,norm,phi)
      implicit none
      type(bim_data_t) :: bd
      ! pass
      integer :: n
      double precision, dimension(n), intent(in) :: xi, R, Z, J, Bn
      double precision, dimension(n,3), intent(in) :: norm
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
      allocate(bd%J(0:n+1))
      bd%J(1:n) = J
      allocate(bd%Bn(0:n+1))
      bd%Bn(1:n) = Bn
      allocate(bd%phi(0:n+1))
      allocate(bd%norm(0:n+1,3))
      bd%norm(1:n,:) = norm
      ! matrices
      allocate(bd%L(n,n))
      allocate(bd%K(n,n))
      allocate(bd%Kinv(n,n))
      allocate(bd%KinvL(n,n))

      ! compute normal vector magnitudes
      allocate(bd%nmag(0:n+1))
      do i = 1, n
         bd%nmag(i) = sqrt(bd%norm(i,1)**2 + bd%norm(i,3)**2)
      end do

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
      bd%J(0) = bd%J(n)
      bd%Bn(0) = bd%Bn(n)
      bd%norm(0,:) = bd%norm(n,:)
      bd%nmag(0) = bd%nmag(n)
      ! hi side
      bd%xi(bd%n+1) = bd%xi(n) + bd%dxi
      bd%R(n+1) = bd%R(1)
      bd%Z(n+1) = bd%Z(1)
      bd%J(n+1) = bd%J(1)
      bd%Bn(n+1) = bd%Bn(1)
      bd%norm(n+1,:) = bd%norm(1,:)
      bd%nmag(n+1) = bd%nmag(1)

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
      double precision :: Lij, delta

      if (i==j) then
         Lij = integrate_toroidal(bd, i, j, SINGULAR, integrand_g)
      else
         Lij = integrate_toroidal(bd, i, j, REGULAR, integrand_g)
      end if

   end function L_matrix_element

   function K_matrix_element(bd,i,j) result(Kij)
      implicit none
      type(bim_data_t) :: bd
      integer, intent(in) :: i, j
      double precision :: Kij, delta

      if (i==j) then
         Kij = integrate_toroidal(bd, i, j, SINGULAR, integrand_dgdn)
         Kij = Kij - 0.5d0
      else
         Kij = integrate_toroidal(bd, i, j, REGULAR, integrand_dgdn)
      end if

   end function K_matrix_element

   ! TOROIDAL TORSOINAL GREEN'S FUNCTION AND DERIVATIVE

   function evalchi(R,Z,RP,ZP) result(result)
      implicit none
      double precision, intent(in) :: R, Z, RP, ZP
      double precision :: result

      result = (R*R + RP*RP + (Z-ZP)**2)/(2*R*RP)
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
      !write(*,"(a9,8(es23.15))") "greens:", mu, chi, K, pi, Ri, Rj, Zi, Zj
   end function greens

   function dgreensdR(Ri,Zi,Rj,Zj) result(dGdR)
      implicit none
      double precision, intent(in) :: Ri, Zi, Rj, Zj
      double precision :: chi, mu, K, dGdR, dchidR

      chi = evalchi(Ri,Zi,Rj,Zj)
      mu = evalmu(chi)
      K = Kelliptic(mu)
      dchidR = evaldchidR(Ri,Zi,Rj,Zj)

      !dGdR = -greens(Ri,Zi,Rj,Zj)/(2*Ri) + mu**2/(8*pi*sqrt(Ri*Rj))*&
         !(Eelliptic(mu)/(1-mu**2) - (1-mu)*Kelliptic(mu)) * &
         !dchidR

      dGdR = -greens(Ri,Zi,Rj,Zj)/(2*Ri) + mu**3/(8*pi*sqrt(Ri*Rj))*&
         Eelliptic(mu)/(1-mu**2)*dchidR

   end function dgreensdR

   function dgreensdZ(Ri,Zi,Rj,Zj) result(dGdZ)
      implicit none
      double precision, intent(in) :: Ri, Zi, Rj, Zj
      double precision :: chi, mu, K, dGdZ, dchidZ

      chi = evalchi(Ri,Zi,Rj,Zj)
      mu = evalmu(chi)
      K = Kelliptic(mu)
      dchidZ = evaldchidZ(Ri,Zi,Rj,Zj)

      !dGdZ = mu**2/(8*pi*sqrt(Ri*Rj))*&
         !(Eelliptic(mu)/(1-mu**2) - (1-mu)*Kelliptic(mu)) * &
         !dchidZ

      dGdZ = mu**3/(8*pi*sqrt(Ri*Rj))*Eelliptic(mu)/(1-mu**2)*dchidZ

   end function dgreensdZ

   ! INTEGRANDS

   function integrand_g(bd,i,j,xi_j) result(integrand)
      implicit none
      ! pass
      type(bim_data_t), intent(in) :: bd
      integer, intent(in) :: i, j
      double precision, intent(in) :: xi_j
      ! internal
      double precision :: integrand
      double precision :: dxlo, dxhi, Jacj, nmagi, nmagj, Ri, Zi, Rj, Zj
      integer :: lo, hi

      ! interpolant
      lo = merge(j-1,j,xi_j < bd%xi(j))
      hi = merge(j,j+1,xi_j < bd%xi(j))
      dxlo = (xi_j - bd%xi(lo))/bd%dxi
      dxhi = (bd%xi(hi) - xi_j)/bd%dxi
      dxlo = min(max(0.d0,dxlo),1.d0)
      dxhi = min(max(0.d0,dxhi),1.d0)

      ! interpolation
      Jacj = bd%J(lo)*dxhi + bd%J(hi)*dxlo
      nmagj = bd%nmag(lo)*dxhi + bd%nmag(hi)*dxlo
      Rj = bd%R(lo)*dxhi + bd%R(hi)*dxlo
      Zj = bd%Z(lo)*dxhi + bd%Z(hi)*dxlo
      ! no interpolation
      !Jacj = bd%J(j)
      !nmagj = bd%nmag(j)

      Ri = bd%R(i)
      Zi = bd%Z(i)

      integrand = Jacj*nmagj*greens(Ri,Zi,Rj,Zj) / (bd%J(i)*bd%nmag(i))

   end function integrand_g

   function integrand_dgdn(bd,i,j,xi_j) result(integrand)
      implicit none
      ! pass
      type(bim_data_t), intent(in) :: bd
      integer, intent(in) :: i, j
      double precision, intent(in) :: xi_j
      ! internal
      double precision :: integrand

      double precision :: Jacj, nmagi, nmagj, Ri, Zi, Rj, Zj
      double precision :: dxlo, dxhi, nidoter, nidotez
      integer :: lo, hi

      ! interpolant
      lo = merge(j-1,j,xi_j < bd%xi(j))
      hi = merge(j,j+1,xi_j < bd%xi(j))
      dxlo = (xi_j - bd%xi(lo))/bd%dxi
      dxhi = (bd%xi(hi) - xi_j)/bd%dxi
      dxlo = min(max(0.d0,dxlo),1.d0)
      dxhi = min(max(0.d0,dxhi),1.d0)

      ! interpolation
      Jacj = bd%J(lo)*dxhi + bd%J(hi)*dxlo
      nmagj = bd%nmag(lo)*dxhi + bd%nmag(hi)*dxlo
      Rj = bd%R(lo)*dxhi + bd%R(hi)*dxlo
      Zj = bd%Z(lo)*dxhi + bd%Z(hi)*dxlo
      ! no interpolation
      !Jacj = bd%J(j)
      !nmagj = bd%nmag(j)

      nidoter = bd%norm(i,1)/bd%nmag(i)
      nidotez = bd%norm(i,3)/bd%nmag(i)

      Ri = bd%R(i)
      Zi = bd%Z(i)

      integrand = Jacj*nmagj*&
         (nidoter*dgreensdR(Ri,Zi,Rj,Zj) + nidotez*dgreensdZ(Ri,Zi,Rj,Zj))

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
         print *, "error large in quadrature rule:"
         print *, result18, result36, error
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
