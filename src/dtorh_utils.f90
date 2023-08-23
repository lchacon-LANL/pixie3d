module dtorh_utils
   implicit none
   public

contains

   SUBROUTINE FRAC(Z,M,NMAX,EPS,TINYSQ,FC) 
      implicit none
      INTEGER M,NMAX,N,MM          
      DOUBLE PRECISION Z,EPS,TINYSQ,FC,DZ2,DN0,DN1,DN2,DN3,DN4,&
         B,A,C0,D0,DELTA
      N=NMAX
      MM=0
      DZ2=Z+Z
      DN0=N+M
      DN1=N+1.D0
      DN2=DN0+0.5D0
      DN3=DN0-0.5D0
      DN4=1.D0-M-M
      B=2.D0*DN1*Z/DN2 
      A=1.D0
      FC=TINYSQ
      C0=FC
      D0=0.D0
   81       D0=B+A*D0
      IF (DABS(D0).LT.TINYSQ) D0=TINYSQ
      C0=B+A/C0
      IF (DABS(C0).LT.TINYSQ) C0=TINYSQ
      D0=1.D0/D0
      DELTA=C0*D0
      FC=FC*DELTA
      MM=MM+1
      A=-(1.D0+DN4/(DN3+MM))  
      B=DZ2*(DN1+MM)/(DN2+MM)
      IF (MM.LT.10000) THEN
         IF (DABS(DELTA-1.D0).GT.EPS) GOTO 81
      END IF
      IF (MM.EQ.10000) WRITE(*,*)'CF CONVERGENCE FAILS'
      RETURN
   END subroutine frac

   SUBROUTINE FRACPS(QZ,M,N,EPS,TINYSQ,FC)
      implicit none
      INTEGER M,N,MM         
      DOUBLE PRECISION QZ,EPS,TINYSQ,FC,DN2,DZ2,DM,B,A,C0,D0,DELTA
      MM=0
      DN2=N*N
      DZ2=QZ+QZ
      DM=M-0.5D0
      B=DZ2*M
      A=DN2-DM*DM
      FC=TINYSQ
      C0=FC
      D0=0.D0
   82       D0=B+A*D0
      IF (DABS(D0).LT.TINYSQ) D0=TINYSQ
      C0=B+A/C0
      IF (DABS(C0).LT.TINYSQ) C0=TINYSQ
      D0=1.D0/D0
      DELTA=C0*D0
      FC=FC*DELTA
      MM=MM+1
      A=DN2-(MM+DM)*(MM+DM)  
      B=DZ2*(MM+M)
      IF (MM.LT.10000) THEN
         IF (ABS(DELTA-1.D0).GT.EPS) GOTO 82
      END IF
      IF (MM.EQ.10000) WRITE(*,*)'CF CONVERGENCE FAILS'
      RETURN
      END subroutine fracps


      DOUBLE PRECISION FUNCTION ELLIP2(AK)
         implicit none
         DOUBLE PRECISION AK
         DOUBLE PRECISION Q!,DRD,DRF
         Q=(1.D0-AK)*(1.D0+AK)
         ELLIP2=DRF(Q)-(AK)**2*DRD(Q)/3.D0
      END function ellip2
      
      DOUBLE PRECISION FUNCTION ELLIP1(AK)
         implicit none
         DOUBLE PRECISION AK!,DRF
         ELLIP1=DRF((1.D0-AK)*(1.D0+AK))
      END function ellip1
      
     
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   !  WE USE SLATEC ROUTINES IN THE EVALUATION OF ELLIPTIC INTEGRALS         C
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  

   DOUBLE PRECISION FUNCTION DRF (Y)
      implicit none
      DOUBLE PRECISION EPSLON, ERRTOL
      DOUBLE PRECISION C1, C2, C3, E2, E3, LAMDA
      DOUBLE PRECISION MU, S, XN, XNDEV
      DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, ZN, ZNDEV, ZNROOT
      LOGICAL FIRST
      SAVE ERRTOL,C1,C2,C3,FIRST
      DATA FIRST /.TRUE./
      !
      !***FIRST EXECUTABLE STATEMENT  DRF
      !
      IF (FIRST) THEN
         ERRTOL = 1.D-8
         C1 = 1.0D0/24.0D0
         C2 = 3.0D0/44.0D0
         C3 = 1.0D0/14.0D0
      ENDIF
      FIRST = .FALSE.
      DRF = 0.0D0
      XN = 0.D0
      YN = Y
      ZN = 1.D0

   30 MU = (XN+YN+ZN)/3.0D0
      XNDEV = 2.0D0 - (MU+XN)/MU
      YNDEV = 2.0D0 - (MU+YN)/MU
      ZNDEV = 2.0D0 - (MU+ZN)/MU
      EPSLON = MAX(ABS(XNDEV),ABS(YNDEV),ABS(ZNDEV))
      IF (EPSLON.LT.ERRTOL) GO TO 40
      XNROOT = SQRT(XN)
      YNROOT = SQRT(YN)
      ZNROOT = SQRT(ZN)
      LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
      XN = (XN+LAMDA)*0.250D0
      YN = (YN+LAMDA)*0.250D0
      ZN = (ZN+LAMDA)*0.250D0
      GO TO 30

   40 E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
      E3 = XNDEV*YNDEV*ZNDEV
      S  = 1.0D0 + (C1*E2-0.10D0-C2*E3)*E2 + C3*E3
      DRF = S/SQRT(MU)

   END function drf


   DOUBLE PRECISION FUNCTION DRD (Y)
      implicit none
      DOUBLE PRECISION EPSLON, ERRTOL
      DOUBLE PRECISION C1, C2, C3, C4, EA, EB, EC, ED, EF, LAMDA
      DOUBLE PRECISION MU, POWER4, SIGMA, S1, S2, XN, XNDEV
      DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT,  ZN, ZNDEV, ZNROOT    
      LOGICAL FIRST
      SAVE ERRTOL, C1, C2, C3, C4, FIRST
      DATA FIRST /.TRUE./
      !
      !***FIRST EXECUTABLE STATEMENT  DRD
      IF (FIRST) THEN
         ERRTOL = 1.D-8
         C1 = 3.0D0/14.0D0
         C2 = 1.0D0/6.0D0
         C3 = 9.0D0/22.0D0
         C4 = 3.0D0/26.0D0
      ENDIF
      FIRST = .FALSE.
      !
      !         CALL ERROR HANDLER IF NECESSARY.
      !
      DRD = 0.0D0

      XN = 0.D0
      YN = Y
      ZN = 1.D0
      SIGMA = 0.0D0
      POWER4 = 1.0D0
      !      
   30 MU = (XN+YN+3.0D0*ZN)*0.20D0
      XNDEV = (MU-XN)/MU
      YNDEV = (MU-YN)/MU
      ZNDEV = (MU-ZN)/MU
      EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV))
      IF (EPSLON.LT.ERRTOL) GO TO 40
      XNROOT = SQRT(XN)
      YNROOT = SQRT(YN)
      ZNROOT = SQRT(ZN)
      LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
      SIGMA = SIGMA + POWER4/(ZNROOT*(ZN+LAMDA))
      POWER4 = POWER4*0.250D0
      XN = (XN+LAMDA)*0.250D0
      YN = (YN+LAMDA)*0.250D0
      ZN = (ZN+LAMDA)*0.250D0
      GO TO 30
      !
   40 EA = XNDEV*YNDEV
      EB = ZNDEV*ZNDEV
      EC = EA - EB
      ED = EA - 6.0D0*EB
      EF = ED + EC + EC
      S1 = ED*(-C1+0.250D0*C3*ED-1.50D0*C4*ZNDEV*EF)
      S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))
      DRD = 3.0D0*SIGMA + POWER4*(1.0D0+S1+S2)/(MU*SQRT(MU))
      !
   END function drd

   DOUBLE PRECISION FUNCTION GAMMAH(N,OVER)
      implicit none
      INTEGER N,I,J
      DOUBLE PRECISION OVER
      I=N
      J=2*I-1
      GAMMAH=1.D0
   85 IF ((J.GE.1).AND.(GAMMAH.LT.OVER)) THEN
         GAMMAH=GAMMAH*J/2.D0
         I=I-1
         J=2*I-1
         GOTO 85
      ENDIF  
      IF (J.GT.1) THEN 
         GAMMAH=0.D0
      END IF   
   END function gammah

   DOUBLE PRECISION FUNCTION FACTCO(N,PL,M)
      implicit none
      INTEGER N,M,J  
      DOUBLE PRECISION PL,X1,X2 
      FACTCO=1.D0/PL
      X1=M+0.5D0
      X2=-M+0.5D0
      J=N
   86 IF (J.GE.0) THEN
         FACTCO=FACTCO*(J+X1)/(J+X2)
         J=J-1
         GOTO 86
      ENDIF
   END function factco

   SUBROUTINE EXPAN(Z,MODE,IPRE,OVER,QARGU,M,PL)
      implicit none
      INTEGER MODE,IPRE,M,K,I 
      DOUBLE PRECISION Z,OVER,QARGU,PL,PI,PISQ,DB,FL,CC,DZ,GAMMA,AR, &
         DFAC,DF1,A0,Z2I,DA1,DA2,DELTA,SUM,DCC,DCCP,DKK!,GAMMAH,FACTOR,PSI
      DOUBLE PRECISION PRECI(2)
      PARAMETER(PI=3.14159265358979323D0)
      PRECI(1)=1.D-13
      PRECI(2)=1.D-9
      PISQ=DSQRT(PI)
      DB=2.D0*DLOG(2.D0)
      FL=M/2.
      CC=ABS(FLOAT(INT(FL))-FL)
      IF (CC.LT.0.4D0) THEN
         AR=1.d0
      ELSE
         AR=-1.d0
      END IF      
      DZ=1.D0
      DO I=1,M
         DZ=DZ/QARGU
      end do
      IF (MODE.EQ.0) THEN
         GAMMA=GAMMAH(M,OVER)*AR*PISQ 
      ELSE 
         GAMMA=AR
      END IF
      DFAC=2.D0/PI*DZ*GAMMA/PISQ
      DF1=DLOG(2.D0*Z)
      A0=1.D0/DSQRT(2.D0*Z)
      Z2I=1.D0/(Z*Z)
      DELTA=1.D0
      SUM=0.D0
      K=0
      DA2=FACTOR(M,0)
      DA1=DB+PSI(M,0)
   87 DELTA=(DF1+DA1)*DA2*A0
      SUM=SUM+DELTA
      DCC=0.5D0+M+2.D0*K
      DCCP=DCC+1.D0
      DKK=K+1.D0
      DA2=DA2*DCCP*DCC/(DKK*DKK)
      DA1=DA1+1.D0/DKK-1.D0/DCCP-1.D0/DCC
      K=K+1
      A0=A0*0.25D0*Z2I
      IF (ABS(DELTA/SUM).GT.PRECI(IPRE)) GOTO 87
      PL=SUM*DFAC
   END subroutine expan

   DOUBLE PRECISION FUNCTION FACTOR(M,K)
      implicit none
      INTEGER M,K,N,I,J
      DOUBLE PRECISION X1
      FACTOR=1.D0
      IF (K.GE.1) THEN
         X1=M+0.5D0
         N=2*K-1
         I=K
         J=N
 88      IF (J.GE.0) THEN
            FACTOR=FACTOR*(J+X1)/I/I
            J=J-1
            I=I-1
            IF (I.EQ.0) I=1
            GOTO 88
         ENDIF 
      END IF  
   END function factor

   DOUBLE PRECISION FUNCTION PSI(M,K)
      implicit none
      INTEGER M,K,N,J,I
      DOUBLE PRECISION FACTR1,FACTR2
      PSI=0.D0
      FACTR1=0.D0
      FACTR2=0.D0
      N=2*K+M
      J=1
      IF (K.GE.1) THEN
   89       IF (J.LE.K) THEN
            FACTR1=FACTR1+1.D0/J
            J=J+1
            GOTO 89
         ENDIF
      END IF
      I=1
   90 IF (I.LE.N) THEN
         FACTR2=FACTR2+1.D0/(2.D0*I-1.D0)
         I=I+1
         GOTO 90
      ENDIF
      IF (N.EQ.0) FACTR2=0.D0
      PSI=FACTR1-2.D0*FACTR2
   END function psi

end module dtorh_utils


