c-----------------------------------------------------------------------
c     file sing.f.
c     computations relating to singular surfaces.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. sing_mod.
c     1. sing_scan.
c     2. sing_find.
c     3. sing_lim.
c     4. sing_vmat.
c     5. sing_mmat.
c     6. sing_solve.
c     7. sing_matmul.
c     8. sing_vmat_diagnose.
c     9. sing_get_ua.
c     10. sing_get_ca.
c     11. sing_der.
c     12. sing_ua_diagnose.
c-----------------------------------------------------------------------
c     subprogram 0. sing_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE sing_mod
      USE fourfit_mod
      IMPLICIT NONE

      INTEGER :: msol,sing_order=2
      INTEGER, DIMENSION(:), POINTER :: r1,r2,n1,n2
      COMPLEX(r8), DIMENSION(2,2) :: m0mat

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. sing_scan.
c     scans singular surfaces and prints information about them.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_scan

      INTEGER :: ising
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/2x,"i",5x,"psi",8x,"rho",9x,"q",10x,"q1",8x,"di0",9x,"di",
     $     8x,"err"/)
c-----------------------------------------------------------------------
c     scan over singular surfaces and print output.
c-----------------------------------------------------------------------
      WRITE(out_unit,'(/1x,a)')"Singular Surfaces:"
      WRITE(out_unit,10)
      msol=mpert
      DO ising=1,msing
         CALL sing_vmat(ising)
      ENDDO
      WRITE(out_unit,10)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_scan
c-----------------------------------------------------------------------
c     subprogram 2. sing_find.
c     finds positions of singular values of q.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_find

      INTEGER, PARAMETER :: itmax=200,nsing=1000
      INTEGER :: iex,ising,m,dm,it
      INTEGER, DIMENSION(nsing) :: m_sing
      REAL(r8) :: dq,psifac,psifac0,psifac1,singfac
      REAL(r8), DIMENSION(nsing) :: psising,qsing,q1sing
c-----------------------------------------------------------------------
c     start loop over extrema to find singular surfaces.
c-----------------------------------------------------------------------
      ising=0
      DO iex=1,mex
         dq=qex(iex)-qex(iex-1)
         m=nn*qex(iex-1)
         IF(dq > 0)m=m+1
         dm=SIGN(one,dq*nn)
c-----------------------------------------------------------------------
c     find singular surfaces by binary search.
c-----------------------------------------------------------------------
         DO
            IF((m-nn*qex(iex-1))*(m-nn*qex(iex)) > 0)EXIT
            it=0
            psifac0=psiex(iex-1)
            psifac1=psiex(iex)
            DO
               it=it+1
               psifac=(psifac0+psifac1)/2
               CALL spline_eval(sq,psifac,0)
               singfac=(m-nn*sq%f(4))*dm
               IF(singfac > 0)THEN
                  psifac0=psifac
                  psifac=(psifac+psifac1)/2
               ELSE
                  psifac1=psifac
                  psifac=(psifac+psifac0)/2
               ENDIF
               IF(ABS(singfac) <= 1e-12)EXIT
               IF(it > itmax)
     $              CALL program_stop("sing_find can't find root")
            ENDDO
c-----------------------------------------------------------------------
c     store singular surfaces.
c-----------------------------------------------------------------------
            ising=ising+1
            CALL spline_eval(sq,psifac,1)
            m_sing(ising)=m
            qsing(ising)=REAL(m,r8)/nn
            q1sing(ising)=sq%f1(4)
            psising(ising)=psifac
            m=m+dm
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     transfer to permanent storage.
c-----------------------------------------------------------------------
      msing=ising
      ALLOCATE(sing(msing))
      DO ising=1,msing
         sing(ising)%m=m_sing(ising)
         sing(ising)%psifac=psising(ising)
         sing(ising)%rho=SQRT(psising(ising))
         sing(ising)%q=qsing(ising)
         sing(ising)%q1=q1sing(ising)
      ENDDO
      DEALLOCATE(psiex,qex)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_find
c-----------------------------------------------------------------------
c     subprogram 3. sing_lim.
c     computes limiter values.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_lim

      INTEGER :: it,itmax=50
      INTEGER, DIMENSION(1) :: jpsi
      REAL(r8) :: dpsi,q,q1,eps=1e-10
c-----------------------------------------------------------------------
c     compute psilim and qlim.
c-----------------------------------------------------------------------
      qlim=qmax
      q1lim=sq%fs1(mpsi,4)
      psilim=psihigh
      IF(.NOT. sas_flag)RETURN
c-----------------------------------------------------------------------
c     normalze dmlim to interval [0,1).
c-----------------------------------------------------------------------
      DO
         IF(dmlim < 1)EXIT
         dmlim=dmlim-1
      ENDDO
      DO
         IF(dmlim >= 0)EXIT
         dmlim=dmlim+1
      ENDDO
c-----------------------------------------------------------------------
c     compute qlim.
c-----------------------------------------------------------------------
      qlim=(INT(nn*qlim)+dmlim)/nn
      DO
         IF(qlim <= qmax)EXIT
         qlim=qlim-1._r8/nn
      ENDDO
c-----------------------------------------------------------------------
c     use newton iteration to find psilim.
c-----------------------------------------------------------------------
      jpsi=MINLOC(ABS(sq%fs(:,4)-qlim))
      psilim=sq%xs(jpsi(1))
      it=0
      DO
         it=it+1
         CALL spline_eval(sq,psilim,1)
         q=sq%f(4)
         q1=sq%f1(4)
         dpsi=(qlim-q)/q1
         psilim=psilim+dpsi
         IF(ABS(dpsi) < eps*ABS(psilim) .OR. it > itmax)EXIT
      ENDDO
      q1lim=q1
c-----------------------------------------------------------------------
c     abort if not found.
c-----------------------------------------------------------------------
      IF(it > itmax)THEN
         CALL program_stop("Can't find psilim.")
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_lim
c-----------------------------------------------------------------------
c     subprogram 4. sing_vmat.
c     computes asymptotic behavior at the singular surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_vmat(ising)

      INTEGER, INTENT(IN) :: ising

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: ipert0,ipert,k
      REAL(r8) :: psifac,di,di0,q,q1,rho,dpsi
      REAL(r8), PARAMETER :: half=.5_r8
      COMPLEX(r8) :: det
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     allocate and initialize solution array.
c-----------------------------------------------------------------------
      IF(ising < 1 .OR. ising > msing)RETURN
      singp => sing(ising)
      ALLOCATE(singp%vmat(mpert,2*mpert,2,0:2*sing_order),
     $     singp%mmat(mpert,2*mpert,2,0:2*sing_order+2))
c-----------------------------------------------------------------------
c     zero di if out of range.
c-----------------------------------------------------------------------
      ipert0=NINT(nn*singp%q)-mlow+1
      IF(ipert0 <= 0 .OR. mlow > nn*q .OR. mhigh < nn*q)THEN
         singp%di=0
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     allocate and compute ranges.
c-----------------------------------------------------------------------
      ALLOCATE(sing(ising)%n1(mpert-1),sing(ising)%n2(2*mpert-2))
      q=singp%q
      singp%r1=(/ipert0/)
      singp%r2=(/ipert0,ipert0+mpert/)
      singp%n1=(/(ipert,ipert=1,ipert0-1),(ipert,ipert=ipert0+1,mpert)/)
      singp%n2=(/singp%n1,singp%n1+mpert/)
      r1 => singp%r1
      r2 => singp%r2
      n1 => singp%n1
      n2 => singp%n2
c-----------------------------------------------------------------------
c     interpolate local values.
c-----------------------------------------------------------------------
      psifac=singp%psifac
      CALL spline_eval(locstab,psifac,0)
      di0=locstab%f(1)/psifac
      q1=singp%q1
      rho=singp%rho
c-----------------------------------------------------------------------
c     compute Mercier criterion and singular power.
c-----------------------------------------------------------------------
      CALL sing_mmat(ising)
      m0mat=TRANSPOSE(singp%mmat(r1(1),r2,:,0))
      di=m0mat(1,1)*m0mat(2,2)-m0mat(2,1)*m0mat(1,2)
      singp%di=di
      singp%alpha=SQRT(-CMPLX(singp%di))
      ALLOCATE(singp%power(2*mpert))
      singp%power=0
      singp%power(ipert0)=-singp%alpha
      singp%power(ipert0+mpert)=singp%alpha
      WRITE(out_unit,'(i3,1p,7e11.3)')ising,psifac,rho,q,q1,di0,
     $     singp%di,singp%di/di0-1
c-----------------------------------------------------------------------
c     zeroth-order non-resonant solutions.
c-----------------------------------------------------------------------
      singp%vmat=0
      DO ipert=1,mpert
         singp%vmat(ipert,ipert,1,0)=1
         singp%vmat(ipert,ipert+mpert,2,0)=1
      ENDDO
c-----------------------------------------------------------------------
c     zeroth-order resonant solutions.
c-----------------------------------------------------------------------
      singp%vmat(ipert0,ipert0,1,0)=1
      singp%vmat(ipert0,ipert0+mpert,1,0)=1
      singp%vmat(ipert0,ipert0,2,0)
     $     =-(m0mat(1,1)+singp%alpha)/m0mat(1,2)
      singp%vmat(ipert0,ipert0+mpert,2,0)
     $     =-(m0mat(1,1)-singp%alpha)/m0mat(1,2)
      det=CONJG(singp%vmat(ipert0,ipert0,1,0))
     $     *singp%vmat(ipert0,ipert0+mpert,2,0)
     $     -CONJG(singp%vmat(ipert0,ipert0+mpert,1,0))
     $     *singp%vmat(ipert0,ipert0,2,0)
      singp%vmat(ipert0,:,:,0)=singp%vmat(ipert0,:,:,0)/SQRT(det)
c-----------------------------------------------------------------------
c     compute higher-order solutions.
c-----------------------------------------------------------------------
      DO k=1,2*sing_order
         CALL sing_solve(k,singp%power,singp%mmat,singp%vmat)
      ENDDO
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL sing_vmat_diagnose(ising)
         dpsi=-1e-2
         CALL sing_ua_diagnose(ising,dpsi)
         CALL program_stop("Termination by sing_vmat.")
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_vmat
c-----------------------------------------------------------------------
c     subprogram 5. sing_mmat.
c     computes series expansion of coefficient matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_mmat(ising)

      INTEGER, INTENT(IN) :: ising

      INTEGER :: ipert0,ipert,jpert,kpert,isol,iqty,info,msol,m,i,j,n,
     $     fac0,fac1
      INTEGER, DIMENSION(mpert) :: mvec
      REAL(r8) :: psifac
      REAL(r8), DIMENSION(0:3) :: q
      REAL(r8), DIMENSION(mpert,0:3) :: singfac
      COMPLEX(r8), PARAMETER :: one=1,half=.5_r8
      COMPLEX(r8), DIMENSION(mband+1,mpert,0:sing_order) :: f,ff,g
      COMPLEX(r8), DIMENSION(mband+1,mpert) :: f0
      COMPLEX(r8), DIMENSION(2*mband+1,mpert,0:sing_order) :: k
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: v
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2,0:sing_order) :: x
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: mmat

      LOGICAL, PARAMETER :: diagnose=.FALSE.
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"io",4x,"is",4x,"ip",2(4x,"re x",i1,6x,"im x",i1,2x)/)
 20   FORMAT(3i6,1p,4e11.3)
c-----------------------------------------------------------------------
c     evaluate cubic splines.
c-----------------------------------------------------------------------
      msol=2*mpert
      psifac=sing(ising)%psifac
      CALL spline_eval(sq,psifac,3)
      CALL cspline_eval(fmats,psifac,3)
      CALL cspline_eval(gmats,psifac,3)
      CALL cspline_eval(kmats,psifac,3)
c-----------------------------------------------------------------------
c     evaluate safety factor and its derivatives.
c-----------------------------------------------------------------------
      q(0)=sq%f(4)
      q(1)=sq%f1(4)
      q(2)=sq%f2(4)
      q(3)=sq%f3(4)
c-----------------------------------------------------------------------
c     evaluate singfac and its derivatives.
c-----------------------------------------------------------------------
      ipert0=sing(ising)%m-mlow+1
      mvec=(/(m,m=mlow,mhigh)/)
      singfac=0
      singfac(:,0)=mvec-nn*q(0)
      singfac(:,1)=-nn*q(1)
      singfac(:,2)=-nn*q(2)
      singfac(:,3)=-nn*q(3)
      singfac(ipert0,0)=-nn*q(1)
      singfac(ipert0,1)=-nn*q(2)/2
      singfac(ipert0,2)=-nn*q(3)/3
      singfac(ipert0,3)=0
c-----------------------------------------------------------------------
c     compute factored Hermitian banded matrix F and its derivatives.
c-----------------------------------------------------------------------
      f=0
      iqty=0
      DO jpert=1,mpert
         DO ipert=jpert,MIN(mpert,jpert+mband)
            iqty=iqty+1
            f(1+ipert-jpert,jpert,0)
     $           =singfac(ipert,0)*fmats%f(iqty)
            IF(sing_order < 1)CYCLE
            f(1+ipert-jpert,jpert,1)
     $           =singfac(ipert,0)*fmats%f1(iqty)
     $           +singfac(ipert,1)*fmats%f(iqty)
            IF(sing_order < 2)CYCLE
            f(1+ipert-jpert,jpert,2)
     $           =singfac(ipert,0)*fmats%f2(iqty)
     $           +2*singfac(ipert,1)*fmats%f1(iqty)
     $           +singfac(ipert,2)*fmats%f(iqty)
            IF(sing_order < 3)CYCLE
            f(1+ipert-jpert,jpert,3)
     $           =singfac(ipert,0)*fmats%f3(iqty)
     $           +3*singfac(ipert,1)*fmats%f2(iqty)
     $           +3*singfac(ipert,2)*fmats%f1(iqty)
     $           +singfac(ipert,3)*fmats%f(iqty)
            IF(sing_order < 4)CYCLE
            f(1+ipert-jpert,jpert,4)
     $           =4*singfac(ipert,1)*fmats%f3(iqty)
     $           +6*singfac(ipert,2)*fmats%f2(iqty)
     $           +4*singfac(ipert,3)*fmats%f1(iqty)
            IF(sing_order < 5)CYCLE
            f(1+ipert-jpert,jpert,5)
     $           =10*singfac(ipert,2)*fmats%f3(iqty)
     $           +10*singfac(ipert,3)*fmats%f2(iqty)
            IF(sing_order < 6)CYCLE
            f(1+ipert-jpert,jpert,5)
     $           =20*singfac(ipert,3)*fmats%f3(iqty)
         ENDDO
      ENDDO
      f0=f(:,:,0)
c-----------------------------------------------------------------------
c     compute product of factored Hermitian banded matrix F.
c-----------------------------------------------------------------------
      ff=0
      fac0=1
      DO n=0,sing_order
         fac1=1
         DO j=0,n
            DO jpert=1,mpert
               DO ipert=jpert,MIN(mpert,jpert+mband)
                  DO kpert=MAX(1,ipert-mband),jpert
                     ff(1+ipert-jpert,jpert,n)
     $                    =ff(1+ipert-jpert,jpert,n)
     $                    +fac1*f(1+ipert-kpert,kpert,j)
     $                    *CONJG(f(1+jpert-kpert,kpert,n-j))
                  ENDDO
               ENDDO
            ENDDO
            fac1=fac1*(n-j)/(j+1)
         ENDDO
         ff(:,:,n)=ff(:,:,n)/fac0
         fac0=fac0*(n+1)
      ENDDO
c-----------------------------------------------------------------------
c     compute non-Hermitian banded matrix K.
c-----------------------------------------------------------------------
      k=0
      iqty=0
      DO jpert=1,mpert
         DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
            iqty=iqty+1
            k(1+mband+ipert-jpert,jpert,0)
     $           =singfac(ipert,0)*kmats%f(iqty)
            IF(sing_order < 1)CYCLE
            k(1+mband+ipert-jpert,jpert,1)
     $           =singfac(ipert,0)*kmats%f1(iqty)
     $           +singfac(ipert,1)*kmats%f(iqty)
            IF(sing_order < 2)CYCLE
            k(1+mband+ipert-jpert,jpert,2)
     $           =singfac(ipert,0)*kmats%f2(iqty)/2
     $           +singfac(ipert,1)*kmats%f1(iqty)
     $           +singfac(ipert,2)*kmats%f(iqty)/2
            IF(sing_order < 3)CYCLE
            k(1+mband+ipert-jpert,jpert,3)
     $           =singfac(ipert,0)*kmats%f3(iqty)/6
     $           +singfac(ipert,1)*kmats%f2(iqty)/2
     $           +singfac(ipert,2)*kmats%f1(iqty)/2
     $           +singfac(ipert,3)*kmats%f(iqty)/6
            IF(sing_order < 4)CYCLE
            k(1+mband+ipert-jpert,jpert,4)
     $           =singfac(ipert,1)*kmats%f3(iqty)/6
     $           +singfac(ipert,2)*kmats%f2(iqty)/4
     $           +singfac(ipert,3)*kmats%f1(iqty)/6
            IF(sing_order < 5)CYCLE
            k(1+mband+ipert-jpert,jpert,5)
     $           =singfac(ipert,2)*kmats%f3(iqty)/12
     $           +singfac(ipert,3)*kmats%f2(iqty)/12
            IF(sing_order < 6)CYCLE
            k(1+mband+ipert-jpert,jpert,6)
     $           =singfac(ipert,3)*kmats%f3(iqty)/36
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute Hermitian banded matrix G.
c-----------------------------------------------------------------------
      g=0
      iqty=0
      DO jpert=1,mpert
         DO ipert=jpert,MIN(mpert,jpert+mband)
            iqty=iqty+1
            g(1+ipert-jpert,jpert,0)=gmats%f(iqty)
            IF(sing_order < 1)CYCLE
            g(1+ipert-jpert,jpert,1)=gmats%f1(iqty)
            IF(sing_order < 2)CYCLE
            g(1+ipert-jpert,jpert,2)=gmats%f2(iqty)/2
            IF(sing_order < 3)CYCLE
            g(1+ipert-jpert,jpert,3)=gmats%f2(iqty)/6
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute identity.
c-----------------------------------------------------------------------
      v=0
      DO ipert=1,mpert
         v(ipert,ipert,1)=1
         v(ipert,ipert+mpert,2)=1
      ENDDO
c-----------------------------------------------------------------------
c     compute zeroth-order x1.
c-----------------------------------------------------------------------
      x=0
      DO isol=1,msol
         x(:,isol,1,0)=v(:,isol,2)
         CALL zgbmv('N',mpert,mpert,mband,mband,-one,k(:,:,0),
     $        2*mband+1,v(:,isol,1),1,one,x(:,isol,1,0),1)
      ENDDO
      CALL zpbtrs('L',mpert,mband,msol,f0,mband+1,x(:,:,1,0),mpert,info)
c-----------------------------------------------------------------------
c     compute higher-order x1.
c-----------------------------------------------------------------------
      DO i=1,sing_order
         DO isol=1,msol
            DO j=1,i
               CALL zhbmv('L',mpert,mband,-one,ff(:,:,j),
     $              mband+1,x(:,isol,1,i-j),1,one,x(:,isol,1,i),1)
            ENDDO
            CALL zgbmv('N',mpert,mpert,mband,mband,-one,k(:,:,i),
     $           2*mband+1,v(:,isol,1),1,one,x(:,isol,1,i),1)
         ENDDO
         CALL zpbtrs('L',mpert,mband,msol,f0,mband+1,x(:,:,1,i),mpert,
     $        info)
      ENDDO
c-----------------------------------------------------------------------
c     compute x2.
c-----------------------------------------------------------------------
      DO i=0,sing_order
         DO isol=1,msol
            DO j=0,i
               CALL zgbmv('C',mpert,mpert,mband,mband,one,k(:,:,j),
     $              2*mband+1,x(:,isol,1,i-j),1,one,x(:,isol,2,i),1)
            ENDDO
            CALL zhbmv('L',mpert,mband,one,g(:,:,i),
     $           mband+1,v(:,isol,1),1,one,x(:,isol,2,i),1)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL ascii_open(debug_unit,"xmat.out","UNKNOWN")
         DO j=0,sing_order
            DO isol=1,2*mpert
               WRITE(debug_unit,10)(i,i,i=1,2)
               DO ipert=1,mpert
                  WRITE(debug_unit,20)j,isol,ipert,
     $                 x(ipert,isol,:,j)
               ENDDO
            ENDDO
         ENDDO
         WRITE(debug_unit,10)(i,i,i=1,2)
         CALL ascii_close(debug_unit)
      ENDIF
c-----------------------------------------------------------------------
c     principal terms of mmat.
c-----------------------------------------------------------------------
      mmat => sing(ising)%mmat
      mmat=0
      j=0
      DO i=0,sing_order
         mmat(r1,r2,:,j)=x(r1,r2,:,i)
         mmat(r1,n2,:,j+1)=x(r1,n2,:,i)
         mmat(n1,r2,:,j+1)=x(n1,r2,:,i)
         mmat(n1,n2,:,j+2)=x(n1,n2,:,i)
         j=j+2
      ENDDO
c-----------------------------------------------------------------------
c     shearing terms.
c-----------------------------------------------------------------------
      mmat(r1,r2(1),1,0)=mmat(r1,r2(1),1,0)+half
      mmat(r1,r2(2),2,0)=mmat(r1,r2(2),2,0)-half
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_mmat
c-----------------------------------------------------------------------
c     subprogram 6. sing_solve.
c     solves iteratively for the next order in the power series vmat.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_solve(k,power,mmat,vmat)

      INTEGER, INTENT(IN) :: k
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: power
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2,0:k), INTENT(IN) :: mmat
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2,0:k), INTENT(INOUT) :: vmat

      INTEGER :: l,isol
      REAL(r8), PARAMETER :: two=2
      COMPLEX(r8) :: det
      COMPLEX(r8), DIMENSION(2,2) :: a
      COMPLEX(r8), DIMENSION(2) :: x
c-----------------------------------------------------------------------
c     compute rhs.
c-----------------------------------------------------------------------
      DO l=1,k
         vmat(:,:,:,k)=vmat(:,:,:,k)
     $        +sing_matmul(mmat(:,:,:,l),vmat(:,:,:,k-l))
      ENDDO
c-----------------------------------------------------------------------
c     evaluate solutions.
c-----------------------------------------------------------------------
      DO isol=1,2*mpert
         a=m0mat
         a(1,1)=a(1,1)-k/two-power(isol)
         a(2,2)=a(2,2)-k/two-power(isol)
         det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
         x=-vmat(r1(1),isol,:,k)
         vmat(r1(1),isol,1,k)=(a(2,2)*x(1)-a(1,2)*x(2))/det
         vmat(r1(1),isol,2,k)=(a(1,1)*x(2)-a(2,1)*x(1))/det
         vmat(n1,isol,:,k)=vmat(n1,isol,:,k)/(power(isol)+k/two)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_solve
c-----------------------------------------------------------------------
c     subprogram 7. sing_matmul.
c     multiplies matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION sing_matmul(a,b) RESULT(c)

      COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: a,b
      COMPLEX(r8), DIMENSION(SIZE(a,1),SIZE(b,2),2) :: c

      CHARACTER(64) :: message
      INTEGER :: i,j,m,n
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      m=SIZE(b,1)
      n=SIZE(b,2)
      IF(SIZE(a,2) /= 2*m)THEN
         WRITE(message,'(2(a,i3))')
     $        "Sing_matmul: SIZE(a,2) = ",SIZE(a,2),
     $        " /= 2*SIZE(b,1) = ",2*SIZE(b,1)
         CALL program_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     main computations.
c-----------------------------------------------------------------------
      c=0
      DO i=1,n
         DO j=1,2
            c(:,i,j)=c(:,i,j)
     $           +MATMUL(a(:,1:m,j),b(:,i,1))
     $           +MATMUL(a(:,1+m:2*m,j),b(:,i,2))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION sing_matmul
c-----------------------------------------------------------------------
c     subprogram 8. sing_vmat_diagnose.
c     diagnoses asymptotic behavior at the singular surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_vmat_diagnose(ising)

      INTEGER, INTENT(IN) :: ising

      INTEGER :: i,isol,iorder,ipert
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"io",4x,"is",4x,"ip",2(4x,"re v",i1,6x,"im v",i1,2x)/)
 20   FORMAT(3i6,1p,4e11.3)
c-----------------------------------------------------------------------
c     set pointer and open file.
c-----------------------------------------------------------------------
      singp => sing(ising)
      CALL ascii_open(debug_unit,"vmat.out","UNKNOWN")
c-----------------------------------------------------------------------
c     write output.
c-----------------------------------------------------------------------
      DO iorder=0,2
         DO isol=1,2*mpert
            WRITE(debug_unit,10)(i,i,i=1,2)
            DO ipert=1,mpert
               WRITE(debug_unit,20)iorder,isol,ipert,
     $              singp%vmat(ipert,isol,:,iorder)
            ENDDO
         ENDDO
      ENDDO
      WRITE(debug_unit,10)(i,i,i=1,2)
c-----------------------------------------------------------------------
c     close files.
c-----------------------------------------------------------------------
      CALL ascii_close(debug_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_vmat_diagnose
c-----------------------------------------------------------------------
c     subprogram 9. sing_get_ua.
c     computes asymptotic series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_get_ua(ising,psifac,ua)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: ua

      INTEGER :: iorder
      INTEGER, DIMENSION(:), POINTER :: r1,r2
      COMPLEX(r8) :: dpsi,pfac,sqrtfac
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmat
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     set pointers.
c-----------------------------------------------------------------------
      singp => sing(ising)
      vmat => singp%vmat
      r1 => singp%r1
      r2 => singp%r2
c-----------------------------------------------------------------------
c     compute distance from singular surface and its powers.
c-----------------------------------------------------------------------
      dpsi=psifac-singp%psifac
      sqrtfac=SQRT(dpsi)
      pfac=ABS(dpsi)**singp%alpha
c-----------------------------------------------------------------------
c     compute power series by Horner's method.
c-----------------------------------------------------------------------
      ua=vmat(:,:,:,2*sing_order)
      DO iorder=2*sing_order-1,0,-1
         ua=ua*sqrtfac+vmat(:,:,:,iorder)
      ENDDO
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      ua(r1,:,1)=ua(r1,:,1)/sqrtfac
      ua(r1,:,2)=ua(r1,:,2)*sqrtfac
      ua(:,r2(1),:)=ua(:,r2(1),:)/pfac
      ua(:,r2(2),:)=ua(:,r2(2),:)*pfac
c-----------------------------------------------------------------------
c     renormalize.
c-----------------------------------------------------------------------
      IF(psifac < singp%psifac)THEN
         ua(:,r2(1),:)=ua(:,r2(1),:)
     $        *ABS(ua(r1(1),r2(1),1))/ua(r1(1),r2(1),1)
         ua(:,r2(2),:)=ua(:,r2(2),:)
     $        *ABS(ua(r1(1),r2(2),1))/ua(r1(1),r2(2),1)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_get_ua
c-----------------------------------------------------------------------
c     subprogram 10. sing_get_ca.
c     computes asymptotic coefficients.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_get_ca(ising,psifac,u,ca)

      INTEGER, INTENT(IN) :: ising
      REAL(r8) , INTENT(IN):: psifac
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: u
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: ca

      INTEGER :: info,msol
      INTEGER, DIMENSION(2*mpert) :: ipiv
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua
      COMPLEX(r8), DIMENSION(2*mpert,2*mpert) :: temp1
      COMPLEX(r8), DIMENSION(2*mpert,SIZE(u,2)) :: temp2
c-----------------------------------------------------------------------
c     compute asymptotic coefficients.
c-----------------------------------------------------------------------
      msol=SIZE(u,2)
      CALL sing_get_ua(ising,psifac,ua)
      temp1(1:mpert,:)=ua(:,:,1)
      temp1(mpert+1:2*mpert,:)=ua(:,:,2)
      CALL zgetrf(2*mpert,2*mpert,temp1,2*mpert,ipiv,info)
      temp2(1:mpert,:)=u(:,:,1)
      temp2(mpert+1:2*mpert,:)=u(:,:,2)
      CALL zgetrs('N',2*mpert,msol,temp1,2*mpert,ipiv,
     $     temp2,2*mpert,info)
      ca(:,1:msol,1)=temp2(1:mpert,:)
      ca(:,1:msol,2)=temp2(mpert+1:2*mpert,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_get_ca
c-----------------------------------------------------------------------
c     subprogram 11. sing_der.
c     evaluates the differential equations of DCON.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_der(neq,psifac,u,du)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(mpert,msol,2), INTENT(IN) :: u
      COMPLEX(r8), DIMENSION(mpert,msol,2), INTENT(OUT) :: du

      INTEGER :: ipert,jpert,isol,iqty,info
      REAL(r8) :: q
      REAL(r8), DIMENSION(mpert) :: singfac
      COMPLEX(r8), PARAMETER :: one=1
      COMPLEX(r8), DIMENSION(mband+1,mpert) :: fmatb,gmatb
      COMPLEX(r8), DIMENSION(2*mband+1,mpert) :: kmatb
c-----------------------------------------------------------------------
c     cubic spline evaluation.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psifac,0)
      CALL cspline_eval(fmats,psifac,0)
      CALL cspline_eval(gmats,psifac,0)
      CALL cspline_eval(kmats,psifac,0)
c-----------------------------------------------------------------------
c     define local scalars.
c-----------------------------------------------------------------------
      q=sq%f(4)
      singfac=mlow-nn*q+(/(ipert,ipert=0,mpert-1)/)
      singfac=1/singfac
c-----------------------------------------------------------------------
c     copy Hermitian banded matrices F and G.
c-----------------------------------------------------------------------
      fmatb=0
      gmatb=0
      iqty=1
      DO jpert=1,mpert
         DO ipert=jpert,MIN(mpert,jpert+mband)
            fmatb(1+ipert-jpert,jpert)=fmats%f(iqty)
            gmatb(1+ipert-jpert,jpert)=gmats%f(iqty)
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     copy non-Hermitian banded matrix K.
c-----------------------------------------------------------------------
      kmatb=0
      iqty=1
      DO jpert=1,mpert
         DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
            kmatb(1+mband+ipert-jpert,jpert)=kmats%f(iqty)
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute du1.
c-----------------------------------------------------------------------
      du=0
      DO isol=1,msol
         du(:,isol,1)=u(:,isol,2)*singfac
         CALL zgbmv('N',mpert,mpert,mband,mband,-one,kmatb,
     $        2*mband+1,u(:,isol,1),1,one,du(:,isol,1),1)
      ENDDO
      CALL zpbtrs('L',mpert,mband,msol,fmatb,mband+1,du,mpert,info)
c-----------------------------------------------------------------------
c     compute du2.
c-----------------------------------------------------------------------
      DO isol=1,msol
         CALL zhbmv('L',mpert,mband,one,gmatb,
     $        mband+1,u(:,isol,1),1,one,du(:,isol,2),1 )
         CALL zgbmv('C',mpert,mpert,mband,mband,one,kmatb,
     $        2*mband+1,du(:,isol,1),1,one,du(:,isol,2),1)
         du(:,isol,1)=du(:,isol,1)*singfac
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_der
c-----------------------------------------------------------------------
c     subprogram 12. sing_ua_diagnose.
c     diagnoses asymptotic solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_ua_diagnose(ising,dpsi)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: dpsi

      INTEGER :: isol,ipert
      REAL(r8) :: psifac
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"is",4x,"ip",6x,"im ua1",5x,"re ua1",
     $     7x,"im ua2",5x,"re ua2"/)
 20   FORMAT(2i6,2x,1p,2e11.3,2x,2e11.3)
c-----------------------------------------------------------------------
c     evaluate asymptotic solutions.
c-----------------------------------------------------------------------
      psifac=sing(ising)%psifac+dpsi
      CALL sing_get_ua(ising,psifac,ua)
c-----------------------------------------------------------------------
c     diagnose asymptotic solutions.
c-----------------------------------------------------------------------
      CALL ascii_open(debug_unit,"ua.out","UNKNOWN")
      WRITE(debug_unit,'(a,i2,a,i1,a,1p,e10.3,a,e9.3)')
     $     "ising = ",ising,", sing_order = ",sing_order,
     $     ", dpsi = ",dpsi,", alpha = ",REAL(sing(ising)%alpha)
      WRITE(debug_unit,10)
      DO isol=1,2*msol
         DO ipert=1,mpert
            WRITE(debug_unit,20)isol,ipert,ua(ipert,isol,:)
         ENDDO
         WRITE(debug_unit,10)
      ENDDO
      CALL ascii_close(debug_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_ua_diagnose
      END MODULE sing_mod
