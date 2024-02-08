      MODULE fourier
!     
!     WRITTEN 08-09-06 BY R. SANCHEZ AS PART OF THE VMEC++ PROJECT (c)
!     
!     PURPOSE: CONVERTS QUANTITIES FROM FOURIER SPACE TO REAL SPACE AND VICEVERSA
!
!       NOTE: CALL FIXARRAY must be used once before calling any of the Fourier subroutines
!       to calculate all necessary cosine/sine factors on fixed mesh of collocation angles
!       NS, NTOR, MPOL and NFP used as in the wout.file
!
      USE stel_kinds
      USE stel_constants
      USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,           &
     &             mpol=>mpol_i, ntor=>ntor_i, nfp=>nfp_i
!Luis      USE timer_mod
      IMPLICIT NONE
!
!     VARIABLE DECLARATIONS  
!
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: orthonorm
      CONTAINS
   
     
        SUBROUTINE TOIJSP(XMN, XUV, UFLAG, VFLAG, IPARITY, IHALF)
!
!       DESCRIPTION:  This subroutine moves a quantity X to real space by summing over its Fourier harmonics X_MN.
!       IF UFLAG = 1, the first poloidal derivative of that quantity is obtained
!       IF VFLAG = 1, the first toroidal derivative of that quantity is obtained
!       IF both flags are equal to one, it gives the joint poloidal, toroidal derivative
!       IPARITY = 0, means the quantity X (NOT any of its derivatives!) is COSINE (even); = 1, X is SINE (odd)
!       IHALF = 0, means quantity is on the radial full mesh; =1, on the half radial mesh
!
        USE stel_kinds
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER(iprec), INTENT(IN):: ihalf                         ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH     
        REAL(rprec), DIMENSION(1:ns,0:mpol,-ntor:ntor),            &
     &    INTENT(IN):: xmn
        REAL(rprec), DIMENSION(1:ns,1:ntheta,1:nzeta),             &
     &    intent(OUT):: xuv
        INTEGER(iprec), INTENT(IN):: uflag, vflag                  ! UFLAG/VFLAG = order of poloidal/toroidal derivative
        INTEGER(iprec), INTENT(IN):: iparity                       ! IPARITY = 0, cosine (EVEN); = 1, sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER(iprec):: m, n, jk, lk, js
        INTEGER(iprec):: isign, in 
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:):: work0,        &
     &    work1, work2
        REAL(rprec) :: ton, toff
!-----------------------------------------------
!Luis        CALL second0(ton)

        ALLOCATE(work0(1:ns,0:mpol,-ntor:ntor))                    ! RS (10/03/06): Added WORK0 for n=0/m=0 normalization

        work0 = xmn
        DO js = 1+ihalf,ns
           work0(js,:,:) = orthonorm(:,:)*work0(js,:,:)
        END DO
             
        ALLOCATE(work1(1:ns,1:ntheta,-ntor:ntor),                        &
     &    work2(1:ns,1:ntheta,-ntor:ntor))
        work1 = zero; work2 = zero 

        DO jk = 1, ntheta
          DO m = 0, mpol                                         ! First DO sum over poloidal modes
            IF (uflag == 0) then
              work1(1+ihalf:ns,jk,:) = work1(1+ihalf:ns,jk,:) +                            &
     &          work0(1+ihalf:ns,m,:)*cosmu(m, jk)
              work2(1+ihalf:ns,jk,:) = work2(1+ihalf:ns,jk,:) +                            &
     &          work0(1+ihalf:ns,m,:)*sinmu(m, jk)
            ELSEIF (uflag == 1) then                             ! Poloidal derivative requested: USE COSMUM, SINMUM
              work1(1+ihalf:ns,jk,:) = work1(1+ihalf:ns,jk,:) -                            &
     &          work0(1+ihalf:ns,m,:)*sinmum(m, jk)
              work2(1+ihalf:ns,jk,:) = work2(1+ihalf:ns,jk,:) +                            &
     &          work0(1+ihalf:ns,m,:)*cosmum(m, jk)
            ELSE
              STOP 'UFLAG > 1 in TOUVSP subroutine'
            ENDIF
          ENDDO
        ENDDO
        DEALLOCATE(work0)

        xuv = zero                                                 ! Make sure that SPATIAL variable is zeroed

        DO lk = 1, nzeta
          
          IF (iparity == 0) THEN                                        ! Do first N=0 mode
          
            IF (vflag == 0) THEN
               xuv(1+ihalf:ns,:,lk) = work1(1+ihalf:ns,:,0) 
            ELSEIF (vflag == 1) THEN
               xuv(1+ihalf:ns,:,lk) = 0.d0   
            ELSE        
               STOP 'VFLAG > 1 in TOUVSP subroutine'            
            ENDIF
          
          ELSE
           
            IF (vflag == 0) THEN
                xuv(1+ihalf:ns,:,lk) = work2(1+ihalf:ns,:,0) 
            ELSEIF (vflag == 1) THEN
                xuv(1+ihalf:ns,:,lk) = 0.d0           
            ELSE
               STOP 'VFLAG > 1 in TOUVSP subroutine' 
            ENDIF          
          
          ENDIF          
          
          DO n = 1, ntor                                          ! Then sum over N>0 and N<0 toroidal modes
            IF (iparity == 0) THEN                                ! COSINE series

              IF (vflag == 0) THEN
                xuv(1+ihalf:ns,:,lk) =  xuv(1+ihalf:ns,:,lk)                       &
                  + (work1(1+ihalf:ns,:,n) + work1(1+ihalf:ns,:,-n))*cosnv(n,lk)   &
                  - (work2(1+ihalf:ns,:,n) - work2(1+ihalf:ns,:,-n))*sinnv(n,lk)
              ELSEIF (vflag == 1) THEN                           ! First toroidal derivative requested
                xuv(1+ihalf:ns,:,lk) =  xuv(1+ihalf:ns,:,lk)                       &
     &            - (work1(1+ihalf:ns,:,n) + work1(1+ihalf:ns,:,-n))*sinnvn(n,lk)  &
     &            - (work2(1+ihalf:ns,:,n) - work2(1+ihalf:ns,:,-n))*cosnvn(n,lk)
              ELSE
                STOP 'VFLAG > 1 in TOUVSP subroutine'
              ENDIF

            ELSE                                                  ! SINE series

              IF (vflag == 0) THEN
                xuv(1+ihalf:ns,:,lk) =  xuv(1+ihalf:ns,:,lk)                      &
     &            + (work1(1+ihalf:ns,:,n) - work1(1+ihalf:ns,:,-n))*sinnv(n,lk)  &
     &            + (work2(1+ihalf:ns,:,n) + work2(1+ihalf:ns,:,-n))*cosnv(n,lk)
              ELSEIF (vflag == 1) THEN                            ! First toroidal derivative requested
                xuv(1+ihalf:ns,:,lk) =  xuv(1+ihalf:ns,:,lk)                      &
     &            + (work1(1+ihalf:ns,:,n) - work1(1+ihalf:ns,:,-n))*cosnvn(n,lk) &
     &            - (work2(1+ihalf:ns,:,n) + work2(1+ihalf:ns,:,-n))*sinnvn(n,lk)
              ELSE
                STOP 'VFLAG > 1 in TOUVSP subroutine'
              ENDIF
              
            ENDIF
          ENDDO
        ENDDO

        DEALLOCATE(work1, work2)

!Luis        CALL second0(toff)
!Luis        time_toijsp = time_toijsp + (toff-ton)

        END SUBROUTINE TOIJSP



        SUBROUTINE TOMNSP(XUV, XMN, IPARITY, IHALF)
!
!       Description: This subroutine moves a quantity X to Fourier space producing its armonics X_MN by
!         summing over all the presribed (toroidal and poloidal) collocation points:
!           theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!         where M = mpol - 1 and N = ntor.
!       IPARITY = 0, means the quantity X (NOT any of its derivatives) is even; = 1, X is odd
!       IHALF = 0, means quantity is on the radial full mesh; =1, on the half radial mesh
!
        USE stel_kinds
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER(iprec), INTENT(IN):: ihalf                         ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH
        REAL(rprec), DIMENSION(1:ns,1:ntheta,1:nzeta),             &
     &    INTENT(IN):: xuv
        REAL(rprec), DIMENSION(1:ns,0:mpol,-ntor:ntor),            &
     &    INTENT(OUT):: xmn
        INTEGER(iprec), INTENT(IN):: iparity                          ! IPARITY = 0, cosine (EVEN); = 1; sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER(iprec):: m, n, jk, lk, isign, in, js
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:):: work1, work2, &
     &    work3, work4
        REAL(rprec) :: ton, toff
!-----------------------------------------------
!Luis        CALL second0(ton)

        ALLOCATE(work1(1:ns,0:mpol,nzeta),                         & 
     &    work2(1:ns,0:mpol,nzeta))
        work1 =zero; work2 = zero

        DO m = 0, mpol
          DO jk = 1, ntheta                                            ! First, add over poloidal collocation angles
            work1(1+ihalf:ns,m,:) = work1(1+ihalf:ns,m,:) +        &
     &        xuv(1+ihalf:ns,jk,:)*cosmui(m, jk)                       ! Note USE of cosmuI/sinmuI; include norm. factors
            work2(1+ihalf:ns,m,:) = work2(1+ihalf:ns,m,:) +        &
     &        xuv(1+ihalf:ns,jk,:)*sinmui(m, jk)
          ENDDO
        ENDDO

        xmn = zero                                                  ! Make sure that FOURIER variable is zeroed
        ALLOCATE(work3(1:ns,0:mpol,nzeta),                         & 
     &    work4(1:ns,0:mpol,nzeta))
        work3 =zero; work4 = zero

        DO n = 0, ntor  
          DO lk = 1, nzeta                                      ! Then, add over toroidal collocation angles
            IF (iparity == 0) THEN                              ! COSINE series              
              work3(1+ihalf:ns,:,lk) = work1(1+ihalf:ns,:,lk)*cosnv(n,lk)
              work4(1+ihalf:ns,:,lk) = work2(1+ihalf:ns,:,lk)*sinnv(n,lk)
              xmn(1+ihalf:ns,:,n) =  xmn(1+ihalf:ns,:,n)                         &
     &          + work3(1+ihalf:ns,:,lk) - work4(1+ihalf:ns,:,lk)
              IF (n .NE. 0) THEN
                xmn(1+ihalf:ns,1:mpol,-n) = xmn(1+ihalf:ns,1:mpol,-n)            &
     &            + work3(1+ihalf:ns,1:mpol,lk) + work4(1+ihalf:ns,1:mpol,lk)
              ENDIF
            ELSE                                                ! SINE series
              work3(1+ihalf:ns,:,lk) = work1(1+ihalf:ns,:,lk)*sinnv(n,lk)
              work4(1+ihalf:ns,:,lk) = work2(1+ihalf:ns,:,lk)*cosnv(n,lk)
              xmn(1+ihalf:ns,:,n) =  xmn(1+ihalf:ns,:,n)                         &
     &          + work3(1+ihalf:ns,:,lk) + work4(1+ihalf:ns,:,lk)
              IF (n .NE. 0) THEN
                xmn(1+ihalf:ns,1:mpol,-n) = xmn(1+ihalf:ns,1:mpol,-n)            &
     &            - work3(1+ihalf:ns,1:mpol,lk) + work4(1+ihalf:ns,1:mpol,lk)
              ENDIF
            ENDIF
          ENDDO
        ENDDO 
 
        xmn(:,0,-ntor:-1) = zero              ! Redundant: To avoid counting (0,-n) for non-zero n.

        DO js = 1+ihalf,ns
           xmn(js,:,:) = orthonorm(:,:)*xmn(js,:,:)
        END DO
        
        DEALLOCATE(work1, work2, work3, work4)

!Luis        CALL second0(toff)
!Luis        time_tomnsp = time_tomnsp + (toff-ton)

        END SUBROUTINE TOMNSP
      

      END MODULE fourier
