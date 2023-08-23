module toroidal_harmonics
   use dtorh
   implicit none

   ! vars for interfacing with dtorh
   integer, parameter :: mmax = 2
   integer :: nmax = 2
   integer, parameter :: mdim = 3
   integer, parameter :: ndim = 3
   double precision :: pl(0:mdim,0:ndim), ql(0:mdim,0:ndim)
   integer :: newm, newn(0:mdim)

   public
contains
   function q0(x) result(result)
      ! P^m_{n-1/2}
      ! m=0,n=0
      use dtorh
      implicit none
      double precision :: x,result
      call dtorh2(x,mdim,ndim,mmax,nmax,pl,ql,newm,newn)
      result = ql(0,0)
   end function q0

   function p0(x) result(result)
      ! P^m_{n-1/2}
      ! m=0,n=0
      use dtorh
      implicit none
      double precision :: x,result
      call dtorh2(x,mdim,ndim,mmax,nmax,pl,ql,newm,newn)
      result = pl(0,0)
   end function p0

   function q1(x) result(result)
      ! P^m_{n-1/2}
      ! m=0,n=1
      use dtorh
      implicit none
      double precision :: x,result
      call dtorh2(x,mdim,ndim,mmax,nmax,pl,ql,newm,newn)
      result = ql(0,1)
   end function q1

   function p1(x) result(result)
      ! P^m_{n-1/2}
      ! m=0,n=1
      use dtorh
      implicit none
      double precision :: x,result
      call dtorh2(x,mdim,ndim,mmax,nmax,pl,ql,newm,newn)
      result = pl(0,1)
   end function p1

   function q2(x) result(result)
      ! P^m_{n-1/2}
      ! m=0,n=2
      use dtorh
      implicit none
      double precision :: x,result
      call dtorh2(x,mdim,ndim,mmax,nmax,pl,ql,newm,newn)
      result = ql(0,2)
   end function q2

   function p2(x) result(result)
      ! m=0,n=2
      use dtorh
      implicit none
      double precision :: x,result
      call dtorh2(x,mdim,ndim,mmax,nmax,pl,ql,newm,newn)
      result = pl(0,2)
   end function p2
end module toroidal_harmonics
