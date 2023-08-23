module dtorh
   use dtorh_utils
   implicit none
   public
   ! by segura & Gil 2000, computer physics communications 
   ! Volume 124, Issue 1, 15 January 2000, Pages 104-122

contains

   SUBROUTINE DTORH2(Z,MDIM,NDIM,MMAX,NMAX,PL,QL,NEWM,NEWN)
      implicit none
      !  INPUT :                                                            C      
      !    Z       ARGUMENT OF THE FUNCTIONS                                C                                                                         
      !    MDIM    M-DIMENSION OF THE ARRAYS: MDIM MUST BE GREATER THAN     C
      !            MMAX                                                     C                       
      !    NDIM    N-DIMENSION OF THE ARRAYS: NDIM MUST BE GREATER THAN     C
      !            NMAX                                                     C
      !    MMAX    MAXIMUM ORDER OF THE FUNCTIONS :                         C   
      !            WE CALCULATE FUNCTIONS OF ALL ORDERS BELOW               C   
      !            MIN(NEWM,MMAX). NEWM IS DEFINED BELOW.                   C   
      !    NMAX    MAXIMUM DEGREE OF THE FUNCTIONS :                        C   
      !            WE GET  FUNCTIONS OF ALL THE DEGREES BELOW               C   
      !            MIN(NEWN(M),NMAX). NEWN(M) IS DEFINED BELOW .            C                                                                           
      !  OUTPUT :                                                           C    
      !   *IF MODE IS EQUAL TO 0:                                           C                                                                            
      !    PL(M,N)                                                          C    
      !            THESE VALUES ARE KEPT IN AN ARRAY                        C   
      !    QL(M,N)                                                          C    
      !            THESE VALUES ARE KEPT IN AN ARRAY                        C   
      !                                                                     C    
      !    NEWM    MAXIMUM  ORDER  OF FUNCTIONS CALCULATED WHEN             C   
      !            QL (MMAX+1,0)   IS LARGER THAN 1/TINY                    C   
      !            (OVERFLOW LIMIT = 1/TINY, TINY IS DEFINED BELOW)         C   
      !    NEWN(M) MAXIMUM  DEGREE  OF FUNCTIONS CALCULATED FOR A           C   
      !            GIVEN ORDER M WHEN PL (M,NMAX+1) IS LARGER THAN          C
      !            1/TINY (OVERFLOW LIMIT = 1/TINY, TINY IS DEFINED BELOW)  C          
      !    NOTE1:  FOR A PRECISION OF 10**(-12), IF Z>5 AND (Z/M)>0.22 THE  C
      !            CODE USES A SERIES EXPANSION FOR PL(M,0).                C
      !            WHEN Z<20 AND (Z/M)<0.22 A CONTINUED FRACTION            C
      !            IS APPLIED.                                              C
      !    NOTE2:  FOR A PRECISION OF 10**(-8), IF Z>5 AND (Z/M)>0.12       C
      !            THE CODE USES A SERIES EXPANSION FOR PL(M,0).WHEN Z<20   C
      !            AND (Z/M)<0.12 A CONTINUED FRACTION IS APPLIED.          C                                                                                             
      !   *IF MODE IS EQUAL TO 1:                                           C    
      !      THE SET OF FUNCTIONS EVALUATED IS:                             C    
      !           PL(M,N)/GAMMA(M+1/2),QL(M,N)/GAMMA(M+1/2),                C    
      !      WHICH ARE RESPECTIVELY STORED IN THE ARRAYS PL(M,N),QL(M,N)    C    
      !      NEWM AND NEWN REFER TO THIS NEW SET OF FUNCTIONS               C    
      !      NOTE1 AND NOTE2 ALSO APPLY IN THIS CASE                        C    
      !   *IF MODE IS EQUAL TO 2:                                           C    
      !      THE CODE PERFORMS AS FOR MODE 1, BUT THE RESTRICTION Z<20      C
      !      FOR THE EVALUATION OF THE CONTINUED FRACTION IS NOT CONSIDERED C                                                                                    C
      !      WARNING: USE ONLY IF HIGH M'S FOR Z>20 ARE REQUIRED. THE       C   
      !      EVALUATION OF THE CF MAY FAIL TO CONVERGE FOR TOO HIGH Z'S     C   
      !  PARAMETERS:                                                        C    
      !   MODE: ENABLES THE ENLARGEMENT OF THE RANGE OF ORDERS AND DEGREES  C
      !         THAT CAN BE EVALUATED.                                      C        
      !   EPS:  CONTROLS THE ACCURACY OF THE CONTINUED FRACTIONS AND        C
      !         SERIES.                                                     C
      !   IPRE: REQUIRED PRECISION IN THE EVALUATION OF TOROIDAL HARMONICS. C    
      !           *IF IPRE=1, PRECISION=10**(-12) (TAKING EPS<10**(-12))    C                                                                             C
      !           *IF IPRE=2, PRECISION=10**(-8) (TAKING EPS<10**(-8))      C    
      !   TINY: SMALL PARAMETER NEAR THE UNDERFLOW LIMIT.                   C                
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !   DECLARATION OF VARIABLES  C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER MODE,IPRE,MDIM,NDIM,MMAX,NMAX,NEWM,MP,ICAL,M,L,NMAXOL,NP,N,I
      INTEGER NEWN(0:MDIM)      
      DOUBLE PRECISION Z,PI,EPS,TINY,OVER,TINYSQ,QZ,QDC1, &
         QARGU,ARGU1,PISQ,D1,FL,CC,AR,GAMMA,DFACQS,FCP,DD,PL0, &
         QM0,DFAC3,FC,DFACC,FC2,DFAC4!,GAMMAH,ELLIP1,ELLIP2,FACTCO
      DOUBLE PRECISION PL(0:MDIM,0:NDIM),QL(0:MDIM,0:NDIM),PR(2)
      PARAMETER(PI=3.14159265358979323D0,EPS=1.D-14,TINY=1.D-290, &
         MODE=0,IPRE=1)
      OVER=1.D0/TINY
      TINYSQ=DSQRT(TINY)
      IF ((IPRE.NE.1).AND.(IPRE.NE.2)) THEN
         WRITE(6,*)'IPRE MUST BE 1 OR 2'
         STOP
      END IF
      PR(1)=.22D0
      PR(2)=.12D0
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !  EPS: REQUIRED ACCURACY FOR THE CONTINUED FRACTION     C
      !        (MODIFIED LENTZ)                                C
      !  TINY: SMALL PARAMETER TO PREVENT OVERFLOWS IN THE CF  C
      !         (CLOSE TO THE UNDERFLOW LIMIT)                 C                              
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (Z.LE.1.D0) THEN
         WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
         STOP
      END IF
      QZ=Z
      QDC1=QZ*QZ-1.D0
      QARGU=QZ/DSQRT(QDC1)
      ARGU1=DSQRT(2.D0/(Z+1.D0))
      PISQ=DSQRT(PI)
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                                        C
      !   WE USE THE CODE IF NMAX IS GREATER THAN OR EQUAL TO 2 C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(NMAX.LT.2) NMAX=2
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !C   WE EVALUATE Q^{0}_{-1/2},Q^{1}_{-1/2}            C
      !C   USING SLATEC ROUTINES FOR ELLIPTIC FUNCTIONS     C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      QL(0,0)=ARGU1*ELLIP1(ARGU1)
      QL(1,0)=-1.D0/DSQRT(2.D0*(QZ-1.D0))*ELLIP2(ARGU1)
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !C    WE APPLY FORWARD RECURRENCE IN M FOR Q'S        C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      MP=1
      IF (MODE.EQ.0) THEN
         !1          IF ((MP.LE.MMAX).AND.(ABS(QL(MP,0)).LT.OVER)) THEN
         !QL(MP+1,0)=-2.D0*MP*QARGU*QL(MP,0) 
         !*       -(MP-0.5D0)*(MP-0.5D0)*QL(MP-1,0) 
         !MP=MP+1
         !GOTO 1
         !ENDIF    
         do while ((MP.LE.MMAX).AND.(ABS(QL(MP,0)).LT.OVER))
            QL(MP+1,0)=-2.D0*MP*QARGU*QL(MP,0)  &
                   -(MP-0.5D0)*(MP-0.5D0)*QL(MP-1,0) 
            MP=MP+1
         end do
         IF ((MP-1).LT.MMAX) MMAX=MP-1
         NEWM=MMAX
      ELSE
         QL(0,0)=QL(0,0)/PISQ
         QL(1,0)=QL(1,0)*2.D0/PISQ
         !2          IF ((MP.LE.MMAX).AND.(ABS(QL(MP,0)).LT.OVER)) THEN
         !D1=MP+0.5D0
         !QL(MP+1,0)=-2.D0*MP*QARGU*QL(MP,0)/D1
         !*      -(MP-0.5D0)*QL(MP-1,0)/D1
         !MP=MP+1
         !GOTO 2
         !ENDIF
         do while ((MP.LE.MMAX).AND.(ABS(QL(MP,0)).LT.OVER)) 
            D1=MP+0.5D0
            QL(MP+1,0)=-2.D0*MP*QARGU*QL(MP,0)/D1 &
                  -(MP-0.5D0)*QL(MP-1,0)/D1
            MP=MP+1
         END do
         IF ((MP-1).LT.MMAX) MMAX=MP-1
         NEWM=MMAX
      END IF
      FL=MMAX/2.D0
      CC=ABS(FLOAT(INT(FL))-FL)
      IF (CC.LT.0.4D0) THEN
         AR=1.
      ELSE
         AR=-1.
      END IF
      IF (MODE.EQ.0) THEN
         GAMMA=GAMMAH(MMAX,OVER)*AR*PISQ  
         IF (ABS(GAMMA).LT.TINY) THEN
            WRITE(6,*)'MMAX IS TOO LARGE FOR MODE=0'
            WRITE(6,*)'BETTER TRY MODE=1'
            STOP
         END IF              
      ELSE 
         GAMMA=AR
      END IF  
      DFACQS=-(GAMMA/QL(MMAX+1,0))*GAMMA/PI
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !C  EVALUATION OF PL(MMAX,0),PL(MMAX+1,0)           C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !C  WE CHOOSE EXPANSION OR CF FOR PL(MMAX,0)        C
      !C  DEPENDING ON THE VALUES OF Z,MMAX AND MODE      C                           
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ICAL=1
      IF ((Z/MMAX).GT.PR(IPRE)) ICAL=2
      IF (Z.LT.5.D0) THEN
         ICAL=1
      ELSE IF (Z.GT.20.D0) THEN
         IF ((MODE.NE.2).AND.(ICAL.EQ.1)) ICAL=0
      END IF
      IF (ICAL.EQ.0) THEN
         WRITE(6,*)'YOU MUST CHOOSE MODE=2'
         STOP
      END IF
      IF (ICAL.EQ.1) THEN
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
         !C   WE CALCULATE THE CF FOR P'S     C
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
         CALL FRACPS(QARGU,MMAX+1,0,EPS,TINYSQ,FCP) 
         DD=MMAX+0.5D0
         IF (MODE.NE.0) THEN
            FCP=FCP/DD
            DFACQS=DFACQS/DD
         END IF
         PL(MMAX,0)=DFACQS/DSQRT(QDC1)/(1.D0-FCP*QL(MMAX,0)/QL(MMAX+1,0))
         PL(MMAX+1,0)=FCP*PL(MMAX,0)
      ELSE
         CALL EXPAN(Z,MODE,IPRE,OVER,QARGU,MMAX,PL0)
         PL(MMAX,0)=PL0
         DD=MMAX+0.5D0
         IF (MODE.NE.0) THEN
            FCP=FCP/DD
            DFACQS=DFACQS/DD
         END IF
         PL(MMAX+1,0)=(QL(MMAX+1,0)/QL(MMAX,0))*(PL(MMAX,0)- DFACQS/DSQRT(QDC1))
      END IF
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !C   EVALUATION OF PL(MMAX,1)    C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      QM0=QL(MMAX,0)
      DFAC3=(GAMMA/QM0)*GAMMA/PI/(0.5D0-MMAX)
      CALL FRAC(Z,MMAX,0,EPS,TINYSQ,FC)
      PL(MMAX,1)=PL(MMAX,0)*FC+DFAC3
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !C   EVALUATION OF PL(MMAX+1,1)   C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      QM0=QL(MMAX+1,0)
      IF (MODE.EQ.0) THEN 
         DFACC=(0.5D0+MMAX)
      ELSE
         DFACC=1.D0/(0.5D0+MMAX)
      END IF      
      DFAC3=-(GAMMA/QM0)*GAMMA*DFACC/PI         
      CALL FRAC(Z,MMAX+1,0,EPS,TINYSQ,FC2)
      PL(MMAX+1,1)=PL(MMAX+1,0)*FC2+DFAC3
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !C  WE APPLY BACKWARD RECURRENCE OVER M TO GET THE  C
      !C  SET PL(M,0),PL(M,1)                             C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MODE.EQ.0) THEN
         DO I=1,MMAX
            MP=MMAX+1-I
            M=MP-1
            PL(M,0)=-(PL(MP+1,0)+2.D0*MP*QARGU*PL(MP,0)) &
               /((0.5D0-MP)*(0.5D0-MP))
         end do
         DO I=1,MMAX
            MP=MMAX+1-I
            M=MP-1
            PL(M,1)=(PL(MP+1,1)+2.D0*MP*QARGU*PL(MP,1)) &
               /((1.5D0-MP)*(0.5D0+MP))
         end do
      ELSE
         DO I=1,MMAX
            MP=MMAX+1-I
            M=MP-1
            PL(M,0)=((MP+0.5D0)*PL(MP+1,0)+2.D0*MP*QARGU*PL(MP,0)) &
               /(0.5D0-MP)
         end do
         DO I=1,MMAX
            MP=MMAX+1-I
            M=MP-1
            PL(M,1)=((MP+0.5D0)*PL(MP+1,1)+2.D0*MP*QARGU*PL(MP,1))* &
               (MP-0.5D0)/((1.5D0-MP)*(0.5D0+MP))
         end do
      END IF              
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !C  NOW, WE PERFORM THE EVALUATION OVER N FOR EACH M  C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO L=0,MMAX
         NMAXOL=NMAX
         M=L
         FL=M/2.D0
         CC=ABS(FLOAT(INT(FL))-FL)
         IF (CC.LT.0.4D0) THEN
            AR=1.
         ELSE
            AR=-1.
         END IF
         IF (MODE.EQ.0) THEN
            GAMMA=GAMMAH(M,OVER)*AR*PISQ
         ELSE 
            GAMMA=AR
         END IF
         NP=1
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                                         C
         !C   WE USE THE RECURRENCE RELATIONS FOR P'S          C 
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         !8        IF ((NP.LE.NMAX).AND.(ABS((NP-M+0.5D0)*PL(L,NP)).LT.OVER)) THEN
         !PL(L,NP+1)=(2.D0*NP*Z*PL(L,NP)-(NP+M-0.5D0)*PL(L,NP-1))/(NP-M+0.5D0)
         !NP=NP+1
         !GOTO 8
         ! ENDIF
         do while ((NP.LE.NMAX).AND.(ABS((NP-M+0.5D0)*PL(L,NP)).LT.OVER))
            PL(L,NP+1)=(2.D0*NP*Z*PL(L,NP)-(NP+M-0.5D0)*PL(L,NP-1))/(NP-M+0.5D0)
            NP=NP+1
         end do
         NMAX=NP-1
         NEWN(L)=NMAX         
         DFAC4=(FACTCO(NMAX,PL(L,NMAX+1),M)*GAMMA)*GAMMA/PI/(NMAX+M+0.5D0)
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         !C   WE EVALUATE THE C.F. FOR Q'S USING LENTZ-THOMPSON        C                                                                        C
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         CALL FRAC(Z,M,NMAX,EPS,TINYSQ,FC)         
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
         !C     EVALUATION OF QL(L,NMAX+1) AND QL(L,NMAX) USING        C
         !C     THE WRONSKIAN W{PL(L,NMAX),QL(L,NMAX)},                C
         !C     THE KNOWN VALUES OF PL(L,NMAX+1) AND PL(L,NMAX)        C
         !C     THE VALUE OF H = QL(L,NMAX+1)/QL(L,NMAX)               C
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         QL(L,NMAX)=DFAC4/(1.D0-FC*PL(L,NMAX)/PL(L,NMAX+1))
         QL(L,NMAX+1)=QL(L,NMAX)*FC      
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         !C   WE USE THE BACKWARD RECURRENCE RELATION FOR Q'S          C                                                      
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO I=1,NMAX
            NP=NMAX+1-I
            N=NP-1
            QL(L,N)=((NP+NP)*Z*QL(L,NP)-(NP-M+0.5D0) * QL(L,NP+1))/(NP+M-0.5D0)
         end do
         NMAX=NMAXOL 
      end do
   END subroutine dtorh2

end module dtorh
