c-----------------------------------------------------------------------
c     file ode.f.
c     sets up and integrates Euler-Lagrange differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. ode_mod.
c     1. ode_run.
c     2. ode_axis_init.
c     3. ode_sing_init.
c     4. ode_ideal_cross.
c     5. ode_resist_cross.
c     6. ode_step.
c     7. ode_unorm.
c     8. ode_fixup.
c     9. ode_test.
c     10. ode_test_fixup.
c-----------------------------------------------------------------------
c     subprogram 0. ode_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ode_mod
      USE ode_output_mod
      USE debug_mod
      IMPLICIT NONE

      LOGICAL :: new
      LOGICAL, PARAMETER :: diagnose_ca=.FALSE.
      CHARACTER(12) :: next
      INTEGER, PRIVATE :: istate,liw,lrw,iopt,itask,itol,mf,ix
      INTEGER, DIMENSION(20), PRIVATE :: iwork
      REAL(r8), PRIVATE :: psiout
      REAL(r8), PARAMETER, PRIVATE :: eps=1e-10
      REAL(r8), DIMENSION(:), POINTER, PRIVATE :: rwork
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: atol

      INTEGER, DIMENSION(:), POINTER :: index
      REAL(r8), DIMENSION(:), POINTER :: unorm,unorm0
      COMPLEX(r8), DIMENSION(:,:), POINTER :: fixfac

      CHARACTER(4) :: sort_type="absm"
      INTEGER :: flag_count

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ode_run.
c     integrates the differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_run

      CHARACTER(64) :: message
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      IF(sing_start <= 0)THEN
         CALL ode_axis_init
      ELSE IF(sing_start <= msing)THEN
         CALL ode_sing_init
      ELSE
         WRITE(message,'(2(a,i2))')
     $        "sing_start = ",sing_start," > msing = ",msing
         CALL program_stop(message)
      ENDIF
      flag_count=0
      CALL ode_output_open
c     CALL bin_open(bin_bal1_unit,"fixup.bin","UNKNOWN","REWIND","none")
      IF(diagnose_ca)CALL ascii_open(ca_unit,"ca.out","UNKNOWN")
c-----------------------------------------------------------------------
c     integrate.
c-----------------------------------------------------------------------
      DO
         DO
            IF(istep > 0)CALL ode_unorm(.FALSE.)
            CALL ode_output_step(unorm)
            IF(ode_test())EXIT
            CALL ode_step
         ENDDO
c-----------------------------------------------------------------------
c     re-initialize.
c-----------------------------------------------------------------------
         IF(ising == ksing)EXIT
         IF(next == "cross")THEN
            IF(res_flag)THEN
               CALL ode_resist_cross
            ELSE
               CALL ode_ideal_cross
            ENDIF
         ELSE
            EXIT
         ENDIF
         flag_count=0
      ENDDO
c-----------------------------------------------------------------------
c     finalize.
c-----------------------------------------------------------------------
      CALL ode_output_close
c     CALL bin_close(bin_bal1_unit)
      DEALLOCATE(rwork,atol,unorm0,unorm,index,fixfac)
      IF(diagnose_ca)CALL ascii_close(ca_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_run
c-----------------------------------------------------------------------
c     subprogram 2. ode_axis_init.
c     initializes variables near axis.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_axis_init
      
      INTEGER :: ipert
      REAL(r8), DIMENSION(mpert) :: key,m
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      new=.TRUE.
      ix=0
      psiout=1
      psifac=sq%xs(0)
c-----------------------------------------------------------------------
c     find next singular surface.
c-----------------------------------------------------------------------
      ising=0
      DO
         ising=ising+1
         IF(ising > msing .OR. psilim < sing(ising)%psifac)EXIT
         q=sing(ising)%q
         IF(mlow <= nn*q .AND. mhigh >= nn*q)EXIT
      ENDDO
      IF(ising > msing .OR. psilim < sing(ising)%psifac)THEN
         psimax=psilim*(1-eps)
         next="finish"
      ELSE
         psimax=sing(ising)%psifac-singfac_min/ABS(nn*sing(ising)%q1)
         next="cross"
      ENDIF
c-----------------------------------------------------------------------
c     allocate and sort solutions by increasing value of |m-ms1|.
c-----------------------------------------------------------------------
      ALLOCATE(u(mpert,mpert,2),du(mpert,mpert,2),u_save(mpert,mpert,2),
     $     unorm0(2*mpert),unorm(2*mpert),index(2*mpert))
      index(1:mpert)=(/(ipert,ipert=1,mpert)/)
      m=mlow-1+index(1:mpert)
      SELECT CASE(sort_type)
      CASE("absm")
         key=ABS(m)
      CASE("sing")
         key=m
         IF(msing > 0)key=key-sing(1)%m
         key=-ABS(key)
      CASE DEFAULT
         CALL program_stop("Cannot recognize sort_type = "
     $        //TRIM(sort_type))
      END SELECT
      CALL bubble(key,index,1,mpert)
c-----------------------------------------------------------------------
c     initialize solutions.
c-----------------------------------------------------------------------
      u=0
      DO ipert=1,mpert
         u(index(ipert),ipert,2)=1
      ENDDO
      msol=mpert
      neq=4*mpert*msol
      u_save=u
      psi_save=psifac
c-----------------------------------------------------------------------
c     initialize integrator parameters.
c-----------------------------------------------------------------------
      istep=0
      iopt=1
      itask=5
      itol=2
      mf=10
      istate=1
c-----------------------------------------------------------------------
c     compute conditions at next singular surface.
c-----------------------------------------------------------------------
      q=sq%fs(0,4)
      IF(msing > 0)THEN
         m1=NINT(nn*sing(ising)%q)
      ELSE
         m1=NINT(nn*qlim)+NINT(SIGN(one,nn*sq%fs1(mpsi,4)))
      ENDIF
      singfac=ABS(m1-nn*q)
c-----------------------------------------------------------------------
c     set up work arrays.
c-----------------------------------------------------------------------
      liw=SIZE(iwork)
      lrw=22+64*mpert*msol
      ALLOCATE(rwork(lrw),atol(mpert,msol,2),fixfac(msol,msol))
      iwork=0
      rwork=0
      rwork(1)=psimax
      rwork(5)=psifac*1e-3
      rwork(11)=rwork(5)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_axis_init
c-----------------------------------------------------------------------
c     subprogram 3. ode_sing_init.
c     re-initializes solution across singular surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_sing_init

      LOGICAL, PARAMETER :: diagnose=.FALSE.,old_init=.FALSE.
      CHARACTER(1), DIMENSION(mpert) :: star
      INTEGER :: ipert,isol,ipert0,ieq,m
      REAL(r8) :: dpsi
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(3x,"mlow",1x,"mhigh",1x,"mpert",5x,"q",8x,"psifac",6x,
     $     "dpsi",4x,"order"//3i6,1p,3e11.3,i6/)
 20   FORMAT(1x,a,i2,a,i3,": ",a)
 30   FORMAT(/5x,"i",
     $     2x,2(3x,"re u(",i1,")",4x,"im u(",i1,")",1x),
     $     2x,2(3x,"re du(",i1,")",3x,"im du(",i1,")")/)
 40   FORMAT((i6,1p,2(2x,4e11.3),a3))
c-----------------------------------------------------------------------
c     initialize position.
c-----------------------------------------------------------------------
      new=.TRUE.
      ising=sing_start
      dpsi=singfac_min/ABS(nn*sing(ising)%q1)*10
      psifac=sing(ising)%psifac+dpsi
      q=sing(ising)%q+dpsi*sing(ising)%q1
c-----------------------------------------------------------------------
c     allocate and initialize solutions.
c-----------------------------------------------------------------------
      ALLOCATE(u(mpert,mpert,2),du(mpert,mpert,2),u_save(mpert,mpert,2),
     $     unorm0(2*mpert),unorm(2*mpert),index(2*mpert))
      IF(old_init)THEN
         u=0
         DO ipert=1,mpert
            u(ipert,ipert,2)=1
         ENDDO
      ELSE
         CALL sing_get_ua(ising,psifac,ua)
         u=ua(:,mpert+1:2*mpert,:)
      ENDIF
      u_save=u
      psi_save=psifac
      msol=mpert
      neq=4*mpert*msol
c-----------------------------------------------------------------------
c     diagnose output.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL sing_der(neq,psifac,u,du)
         CALL ascii_open(init_out_unit,"init.out","UNKNOWN")
         WRITE(init_out_unit,'(1x,a/)')"Output from ode_sing_init"
         WRITE(init_out_unit,10)mlow,mhigh,mpert,q,psifac,dpsi,
     $        sing_order
         ipert0=sing(ising)%m-mlow+1
         star=" "
         star(ipert0)="*"
         m=mlow
         DO isol=1,msol
            WRITE(init_out_unit,20)"isol = ",isol,", m = ",m,star(isol)
            WRITE(init_out_unit,30)(ieq,ieq,ieq=1,2),(ieq,ieq,ieq=1,2),
     $           (ieq,ieq,ieq=1,2)
            WRITE(init_out_unit,40)(ipert,u(ipert,isol,:),
     $           du(ipert,isol,:),star(ipert),ipert=1,mpert)
            WRITE(init_out_unit,30)(ieq,ieq,ieq=1,2),(ieq,ieq,ieq=1,2),
     $           (ieq,ieq,ieq=1,2)
            m=m+1
         ENDDO
         CALL ascii_close(init_out_unit)
         CALL program_stop("Termination by ode_sing_init.")
      ENDIF
c-----------------------------------------------------------------------
c     compute conditions at next singular surface.
c-----------------------------------------------------------------------
      DO
         ising=ising+1
         IF(ising > msing .OR. psilim < sing(ising)%psifac)EXIT
         q=sing(ising)%q
         IF(mlow <= nn*q .AND. mhigh >= nn*q)EXIT
      ENDDO
      IF(ising > msing .OR. psilim < sing(ising)%psifac)THEN
         m1=NINT(nn*qlim)+NINT(SIGN(one,nn*sq%fs1(mpsi,4)))
         psimax=psilim*(1-eps)
         next="finish"
      ELSE
         m1=NINT(nn*sing(ising)%q)
         psimax=sing(ising)%psifac-singfac_min/ABS(nn*sing(ising)%q1)
         next="cross"
      ENDIF
c-----------------------------------------------------------------------
c     set up integrator parameters.
c-----------------------------------------------------------------------
      istep=0
      istate=1
      iopt=1
      itask=5
      itol=2
      mf=10
      istate=1
c-----------------------------------------------------------------------
c     set up work arrays.
c-----------------------------------------------------------------------
      liw=SIZE(iwork)
      lrw=22+64*mpert*msol
      ALLOCATE(rwork(lrw),atol(mpert,msol,2),fixfac(msol,msol))
      iwork=0
      rwork=0
      rwork(1)=psimax
      rwork(5)=dpsi
      rwork(11)=rwork(5)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_sing_init
c-----------------------------------------------------------------------
c     subprogram 4. ode_ideal_cross.
c     re-initializes ideal solutions across singular surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_ideal_cross

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: ipert0,jsol,ipert
      REAL(r8) :: dpsi
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1x,a,1p,e9.3,0p,a,f6.3)
 20   FORMAT(/3x,"ising",3x,"psi",9x,"q",10x,"di",6x,"re alpha",
     $     3x,"im alpha"//i6,1p,5e11.3/)
 30   FORMAT(/3x,"is",4x,"psifac",6x,"dpsi",8x,"q",7x,"singfac",5x,
     $     "eval1"/)
 40   FORMAT(a,1p,e11.3/)
 50   FORMAT(1x,2(a,i3))
 60   FORMAT(/2x,"i",4x,"re ca1",5x,"im ca1",4x,"abs ca1",5x,"re ca2",
     $     5x,"im ca2",4x,"abs ca2"/)
 70   FORMAT(i3,1p,6e11.3)
c-----------------------------------------------------------------------
c     fixup solution at singular surface.
c-----------------------------------------------------------------------
      WRITE(*,10)"psi = ",sing(ising)%psifac,", q = ",sing(ising)%q
      CALL ode_unorm(.TRUE.)
c-----------------------------------------------------------------------
c     diagnose solution before reinitialization.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL sing_get_ca(ising,psifac,u,ca)
         CALL ascii_open(init_out_unit,"reinit.out","UNKNOWN")
         WRITE(init_out_unit,40)
     $        "Output from ode_ideal_cross, singfac = ",singfac
         WRITE(init_out_unit,'(a/)')
     $        "Asymptotic coefficients matrix before reinit:"
         DO jsol=1,msol
            WRITE(init_out_unit,50)"jsol = ",jsol,", m = ",mlow+jsol-1
            WRITE(init_out_unit,60)
            WRITE(init_out_unit,70)(ipert,
     $           ca(ipert,jsol,1),ABS(ca(ipert,jsol,1)),
     $           ca(ipert,jsol,2),ABS(ca(ipert,jsol,2)),ipert=1,mpert)
            WRITE(init_out_unit,60)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     write asymptotic coefficients before reinit.
c-----------------------------------------------------------------------
      IF(bin_euler)THEN
         CALL sing_get_ca(ising,psifac,u,ca)
         WRITE(euler_bin_unit)4
         WRITE(euler_bin_unit)sing(ising)%psifac,sing(ising)%q,
     $        sing(ising)%q1
         WRITE(euler_bin_unit)msol
         WRITE(euler_bin_unit)ca
      ENDIF
c-----------------------------------------------------------------------
c     re-initialize.
c-----------------------------------------------------------------------
      ipert0=NINT(nn*sing(ising)%q)-mlow+1
      dpsi=sing(ising)%psifac-psifac
      psifac=sing(ising)%psifac+dpsi
      CALL sing_get_ua(ising,psifac,ua)
      u(ipert0,:,:)=0
      CALL sing_der(neq,psifac,u,du)
      u=u+du*2*dpsi
      u(ipert0,:,:)=0
      u(:,index(1),:)=ua(:,ipert0+mpert,:)
c-----------------------------------------------------------------------
c     diagnose solution after reinitialization.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL sing_get_ca(ising,psifac,u,ca)
         WRITE(init_out_unit,'(a/)')
     $        "Asymptotic coefficients matrix after reinit:"
         DO jsol=1,msol
            WRITE(init_out_unit,50)"jsol = ",jsol,", m = ",mlow+jsol-1
            WRITE(init_out_unit,60)
            WRITE(init_out_unit,70)(ipert,
     $           ca(ipert,jsol,1),ABS(ca(ipert,jsol,1)),
     $           ca(ipert,jsol,2),ABS(ca(ipert,jsol,2)),ipert=1,mpert)
            WRITE(init_out_unit,60)
         ENDDO
         CALL ascii_close(init_out_unit)
         CALL program_stop("Termination by ode_ideal_cross.")
      ENDIF
c-----------------------------------------------------------------------
c     write asymptotic coefficients after reinit.
c-----------------------------------------------------------------------
      IF(bin_euler)THEN
         CALL sing_get_ca(ising,psifac,u,ca)
         WRITE(euler_bin_unit)msol
         WRITE(euler_bin_unit)ca
         WRITE(euler_bin_unit)
     $        sing(ising)%restype%e,sing(ising)%restype%f,
     $        sing(ising)%restype%h,sing(ising)%restype%m,
     $        sing(ising)%restype%g,sing(ising)%restype%k,
     $        sing(ising)%restype%eta,sing(ising)%restype%rho,
     $        sing(ising)%restype%taua,sing(ising)%restype%taur
      ENDIF
c-----------------------------------------------------------------------
c     find next ising.
c-----------------------------------------------------------------------
      DO
         ising=ising+1
         IF(ising > msing .OR. psilim < sing(ising)%psifac)EXIT
         q=sing(ising)%q
         IF(mlow <= nn*q .AND. mhigh >= nn*q)EXIT
      ENDDO
c-----------------------------------------------------------------------
c     compute conditions at next singular surface.
c-----------------------------------------------------------------------
      IF(ising > msing .OR. psilim < sing(ising)%psifac)THEN
         psimax=psilim*(1-eps)
         m1=NINT(nn*qlim)+NINT(SIGN(one,nn*sq%fs1(mpsi,4)))
         next="finish"
      ELSE
         psimax=sing(ising)%psifac-singfac_min/ABS(nn*sing(ising)%q1)
         m1=NINT(nn*sing(ising)%q)
         WRITE(crit_out_unit,20)ising,sing(ising)%psifac,
     $        sing(ising)%q,sing(ising)%di,sing(ising)%alpha
      ENDIF
c-----------------------------------------------------------------------
c     restart ode solver.
c-----------------------------------------------------------------------
      istate=1
      istep=istep+1
      rwork(1)=psimax
      new=.TRUE.
      u_save=u
      psi_save=psifac
c-----------------------------------------------------------------------
c     write to files.
c-----------------------------------------------------------------------
      WRITE(crit_out_unit,30)
      IF(crit_break)WRITE(crit_bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_ideal_cross
c-----------------------------------------------------------------------
c     subprogram 5. ode_resist_cross.
c     re-initializes resistive solutions across singular surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_resist_cross

      LOGICAL, PARAMETER :: sing_flag=.TRUE.,diagnose=.FALSE.
c$$$      CHARACTER(16) :: message
      INTEGER :: ipert0,msol_old,jsol,ipert
      REAL(r8) :: dpsi
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1x,a,1p,e9.3,0p,a,f6.3)
 20   FORMAT(/3x,"ising",3x,"psi",9x,"q",10x,"di",6x,"re alpha",
     $     3x,"im alpha"//i6,1p,5e11.3/)
 30   FORMAT(/3x,"is",4x,"psifac",6x,"dpsi",8x,"q",7x,"singfac",5x,
     $     "eval1"/)
 40   FORMAT(a,1p,e11.3/)
 50   FORMAT(1x,2(a,i3))
 60   FORMAT(/2x,"i",4x,"re ca1",5x,"im ca1",4x,"abs ca1",5x,"re ca2",
     $     5x,"im ca2",4x,"abs ca2"/)
 70   FORMAT(i3,1p,6e11.3)
c-----------------------------------------------------------------------
c     compute asymptotic coefficients before reinit.
c-----------------------------------------------------------------------
      CALL sing_get_ca(ising,psifac,u,ca)
      ipert0=NINT(nn*sing(ising)%q)-mlow+1
c-----------------------------------------------------------------------
c     sort unorm.
c-----------------------------------------------------------------------
      index(1:msol)=(/(ipert,ipert=1,msol)/)
      CALL bubble(unorm,index,1,mpert)
      CALL bubble(unorm,index,mpert+1,msol)
c-----------------------------------------------------------------------
c     write to crit.out.
c-----------------------------------------------------------------------
      WRITE(*,10)"psi = ",sing(ising)%psifac,", q = ",sing(ising)%q
      WRITE(crit_out_unit,'(/1x,a,i6,a,1p,e9.3,0p,a,f6.3,a,i2)')
     $     "Gaussian Reduction at istep = ",istep,
     $     ", psi = ",sing(ising)%psifac,", q = ",sing(ising)%q,
     $     ", index1 = ",index(1)
c-----------------------------------------------------------------------
c     initialize fixfac.
c-----------------------------------------------------------------------
      fixfac=0
      DO ipert=1,mpert
         fixfac(ipert,ipert)=1
      ENDDO
c-----------------------------------------------------------------------
c     fixup solutions.
c-----------------------------------------------------------------------
      DO jsol=1,msol
         IF(jsol == index(1))CYCLE
         fixfac(index(1),jsol)=-ca(ipert0,jsol,1)/ca(ipert0,index(1),1)
         ca(:,jsol,:)=ca(:,jsol,:)
     $        +ca(:,index(1),:)*fixfac(index(1),jsol)
         ca(ipert0,jsol,1)=0
      ENDDO
c-----------------------------------------------------------------------
c     diagnose solution before reinitialization.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL ascii_open(init_out_unit,"reinit.out","UNKNOWN")
         WRITE(init_out_unit,40)
     $        "Output from ode_resist_cross, singfac = ",singfac
         WRITE(init_out_unit,'(a/)')
     $        "Asymptotic coefficients matrix before reinit:"
         DO jsol=1,msol
            WRITE(init_out_unit,50)"jsol = ",jsol,", m = ",mlow+jsol-1
            WRITE(init_out_unit,60)
            WRITE(init_out_unit,70)(ipert,
     $           ca(ipert,jsol,1),ABS(ca(ipert,jsol,1)),
     $           ca(ipert,jsol,2),ABS(ca(ipert,jsol,2)),ipert=1,mpert)
            WRITE(init_out_unit,60)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     write asymptotic coefficients before reinit.
c-----------------------------------------------------------------------
      IF(bin_euler)THEN
         WRITE(euler_bin_unit)2
         WRITE(euler_bin_unit)sing_flag,msol
         WRITE(euler_bin_unit)fixfac,index(1:msol)
         WRITE(euler_bin_unit)4
         WRITE(euler_bin_unit)sing(ising)%psifac,sing(ising)%q,
     $        sing(ising)%q1
         WRITE(euler_bin_unit)msol
         WRITE(euler_bin_unit)ca
c$$$         IF(diagnose_ca)THEN
c$$$            WRITE(message,'(a,i2.2)')"ca_l_",ising
c$$$            CALL debug1(ca,message,ca_unit)
c$$$         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     reallocate.
c-----------------------------------------------------------------------
      msol_old=msol
      ca_old=ca
      IF(msol < 2*mpert)THEN
         msol=msol+1
         IF(msol < 2*mpert)msol=msol+1
         neq=4*mpert*msol
         lrw=22+64*mpert*msol
         DEALLOCATE(u,du,ca,u_save,rwork,atol,fixfac)
         ALLOCATE(u(mpert,msol,2),du(mpert,msol,2),ca(mpert,msol,2),
     $        u_save(mpert,msol,2),rwork(lrw),atol(mpert,msol,2),
     $        fixfac(msol,msol))
         rwork=0
      ENDIF
c-----------------------------------------------------------------------
c     re-initialize.
c-----------------------------------------------------------------------
      dpsi=sing(ising)%psifac-psifac
      psifac=sing(ising)%psifac+dpsi
      CALL sing_get_ua(ising,psifac,ua)
      u(:,1:msol_old,:)=sing_matmul(ua,ca_old(:,1:msol_old,:))
      u(:,index(1),:)=ua(:,ipert0+mpert,:)
      u(:,msol,:)=ua(:,ipert0,:)
      IF(msol_old < 2*mpert-1)THEN
         u(:,msol-1:msol-1,:)
     $        =sing_matmul(ua,ca_old(:,index(1):index(1),:))
      ENDIF
c-----------------------------------------------------------------------
c     reallocate ca_old.
c-----------------------------------------------------------------------
      DEALLOCATE(ca_old)
      ALLOCATE(ca_old(mpert,msol,2))
c-----------------------------------------------------------------------
c     write asymptotic coefficients after reinit.
c-----------------------------------------------------------------------
      IF(bin_euler)THEN
         CALL sing_get_ca(ising,psifac,u,ca)
         WRITE(euler_bin_unit)msol
         WRITE(euler_bin_unit)ca
c$$$         IF(diagnose_ca)THEN
c$$$            WRITE(message,'(a,i2.2)')"ca_r_",ising
c$$$            CALL debug1(ca,message,ca_unit)
c$$$         ENDIF
         WRITE(euler_bin_unit)
     $        sing(ising)%restype%e,sing(ising)%restype%f,
     $        sing(ising)%restype%h,sing(ising)%restype%m,
     $        sing(ising)%restype%g,sing(ising)%restype%k,
     $        sing(ising)%restype%eta,sing(ising)%restype%rho,
     $        sing(ising)%restype%taua,sing(ising)%restype%taur
      ENDIF
c-----------------------------------------------------------------------
c     diagnose solution after reinitialization.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         WRITE(init_out_unit,'(a/)')
     $        "Asymptotic coefficients matrix after reinit:"
         DO jsol=1,msol
            WRITE(init_out_unit,50)"jsol = ",jsol,", m = ",mlow+jsol-1
            WRITE(init_out_unit,60)
            WRITE(init_out_unit,70)(ipert,
     $           ca(ipert,jsol,1),ABS(ca(ipert,jsol,1)),
     $           ca(ipert,jsol,2),ABS(ca(ipert,jsol,2)),ipert=1,mpert)
            WRITE(init_out_unit,60)
         ENDDO
         CALL ascii_close(init_out_unit)
         CALL program_stop("Termination by ode_resist_cross.")
      ENDIF
c-----------------------------------------------------------------------
c     find next ising.
c-----------------------------------------------------------------------
      DO
         ising=ising+1
         IF(ising > msing .OR. psilim < sing(ising)%psifac)EXIT
         q=sing(ising)%q
         IF(mlow <= nn*q .AND. mhigh >= nn*q)EXIT
      ENDDO
c-----------------------------------------------------------------------
c     compute conditions at next singular surface.
c-----------------------------------------------------------------------
      IF(ising > msing .OR. psilim < sing(ising)%psifac)THEN
         psimax=psilim*(1-eps)
         m1=NINT(nn*qlim)+NINT(SIGN(one,nn*sq%fs1(mpsi,4)))
         next="finish"
      ELSE
         psimax=sing(ising)%psifac-singfac_min/ABS(nn*sing(ising)%q1)
         m1=NINT(nn*sing(ising)%q)
         WRITE(crit_out_unit,20)ising,sing(ising)%psifac,
     $        sing(ising)%q,sing(ising)%di,sing(ising)%alpha
      ENDIF
c-----------------------------------------------------------------------
c     restart ode solver.
c-----------------------------------------------------------------------
      istate=1
      istep=istep+1
      rwork(1)=psimax
      new=.TRUE.
      u_save=u
      psi_save=psifac
c-----------------------------------------------------------------------
c     write to files.
c-----------------------------------------------------------------------
      WRITE(crit_out_unit,30)
      IF(crit_break)WRITE(crit_bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_resist_cross
c-----------------------------------------------------------------------
c     subprogram 6. ode_step.
c     takes a step of the integrator.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_step

      INTEGER :: ipert,isol,ieq,jac,errloc(3)
      REAL(r8) :: singfac,rtol,atol0,tol,dt,errmax,ewtmax
      REAL(r8), PARAMETER :: dpsimax=1e-3,dpsimin=1e-5,dpsifac=2e-2
      COMPLEX(r8) :: err(mpert,msol,2),ewt(mpert,msol,2)

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      CHARACTER(256) :: message
c-----------------------------------------------------------------------
c     compute relative tolerances.
c-----------------------------------------------------------------------
      singfac=HUGE(singfac)
      IF(ising <= msing)singfac=ABS(sing(ising)%m-nn*q)
      IF(ising > 1)singfac=MIN(singfac,ABS(sing(ising-1)%m-nn*q))
      IF(singfac < crossover)THEN
         tol=tol_r
      ELSE
         tol=tol_nr
      ENDIF
      rtol=tol
c-----------------------------------------------------------------------
c     compute absolute tolerances.
c-----------------------------------------------------------------------
      DO ieq=1,2
         DO isol=1,msol
c            atol0=MAXVAL(ABS(u(:,isol,ieq)))*tol*tol
            atol0=MAXVAL(ABS(u(:,isol,ieq)))*tol
            IF(atol0 == 0)atol0=HUGE(atol0)
            atol(:,isol,ieq)=DCMPLX(atol0,atol0)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     choose psiout.
c-----------------------------------------------------------------------
      IF(node_flag)THEN
         DO
            IF(psifac < sq%xs(ix))EXIT
            ix=ix+1
         ENDDO
         psiout=sq%xs(ix)
         psiout=MIN(psiout,psimax)
         rwork(1)=psiout
      ELSE
         psiout=psimax
      ENDIF
c-----------------------------------------------------------------------
c     advance differential equations.
c-----------------------------------------------------------------------
      istep=istep+1
      CALL lsode(sing_der,neq,u,psifac,psiout,itol,rtol,atol,
     $     itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
c-----------------------------------------------------------------------
c     diagnose error.
c-----------------------------------------------------------------------
      IF(rwork(11) < 1e-14 .AND. diagnose)THEN
         dt=rwork(11)
         CALL sing_der(neq,psifac,u,du)
         CALL dewset(neq,itol,rtol,atol,u,ewt)
         err=du*dt/ewt
         errmax=MAXVAL(ABS(err))
         errloc=MAXLOC(ABS(err))
         ewtmax=ewt(errloc(1),errloc(2),errloc(3))
         ipert=errloc(1)
         isol=errloc(2)
         ieq=errloc(3)
         WRITE(message,'(4(a,i3),1p,3(a,e9.3))')
     $        "Termination by ode_step"//CHAR(10)
     $        //" ipert = ",ipert,", ieq = ",ieq,", isol = ",isol,
     $        ", msol = ",msol,CHAR(10)
     $        //" errmax = ",errmax,", ewt = ",ewtmax,
     $        ", atol = ",ABS(atol(ipert,isol,ieq))
         CALL program_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_step
c-----------------------------------------------------------------------
c     subprogram 7. ode_unorm.
c     computes unorm, tests solution matrix for Gaussian reduction.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_unorm(sing_flag)

      LOGICAL, INTENT(IN) :: sing_flag

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      CHARACTER(64) :: message
      INTEGER, DIMENSION(1) :: jmax
      REAL(r8) :: uratio
c-----------------------------------------------------------------------
c     compute norms of first solution vectors, abort if any are zero.
c-----------------------------------------------------------------------
      unorm(1:mpert)=SQRT(SUM(ABS(u(:,1:mpert,1))**2,1))
      unorm(mpert+1:msol)=SQRT(SUM(ABS(u(:,mpert+1:msol,2))**2,1))
      IF(MINVAL(unorm(1:msol)) == 0)THEN
         jmax=MINLOC(unorm(1:msol))
         WRITE(message,'(a,i2,a)')"_unorm: unorm(1,",jmax(1),") = 0"
         CALL program_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     normalize unorm and perform Gaussian reduction if required.
c-----------------------------------------------------------------------
      IF(new)THEN
         new=.FALSE.
         unorm0(1:msol)=unorm(1:msol)
      ELSE
         unorm(1:msol)=unorm(1:msol)/unorm0(1:msol)
         uratio=MAXVAL(unorm(1:msol))/MINVAL(unorm(1:msol))
         IF(uratio > ucrit .OR. sing_flag)THEN
            CALL ode_fixup(sing_flag,.FALSE.)
            IF(diagnose)CALL ode_test_fixup
            new=.TRUE.
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_unorm
c-----------------------------------------------------------------------
c     subprogram 8. ode_fixup.
c     performs Gaussian reduction of solution matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_fixup(sing_flag,test)

      LOGICAL, INTENT(IN) :: sing_flag,test

      LOGICAL, PARAMETER :: diagnose=.FALSE.,secondary=.FALSE.
      LOGICAL, SAVE :: new=.TRUE.
      LOGICAL, DIMENSION(2,msol) :: mask
      INTEGER :: ipert,isol,jsol,kpert,ksol
      INTEGER, DIMENSION(1) :: jmax
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/3x,"is",4x,"psifac",6x,"dpsi",8x,"q",7x,"singfac",5x,
     $     "eval1"/)
 20   FORMAT(1x,a,i6,a,1p,e9.3,0p,a,f6.3)
 40   FORMAT(2(a,i3))
 30   FORMAT(3x,"mlow",1x,"mhigh",1x,"mpert",2x,"msol",3x,"psifac",
     $     6x,"q"//4i6,1p,2e11.3/)
 50   FORMAT(/1x,"i",2(3x,"re u(",i1,")",4x,"im u(",i1,")",1x)/)
 60   FORMAT(i2,1p,4e11.3)
c-----------------------------------------------------------------------
c     open output file and write initial values.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         IF(new)THEN
            CALL ascii_open(init_out_unit,"fixup.out","UNKNOWN")
            new=.FALSE.
         ENDIF
         WRITE(init_out_unit,30)mlow,mhigh,mpert,msol,psifac,q
         WRITE(init_out_unit,'(/a/)')"input values:"
         DO isol=1,msol
            WRITE(init_out_unit,40)"isol = ",isol,", m = ",mlow+isol-1
            WRITE(init_out_unit,50)(ipert,ipert,ipert=1,2)
            WRITE(init_out_unit,60)(ipert,
     $           u(ipert,isol,:)/MAXVAL(ABS(u(:,isol,1))),ipert=1,mpert)
            WRITE(init_out_unit,50)(ipert,ipert,ipert=1,2)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     initial output.
c-----------------------------------------------------------------------
      WRITE(crit_out_unit,10)
      WRITE(crit_out_unit,20)
     $     "Gaussian Reduction at istep = ",istep,
     $     ", psi = ",psifac,", q = ",q
      IF(.NOT. sing_flag)WRITE(crit_out_unit,10)
      istate=1
      flag_count=0
c-----------------------------------------------------------------------
c     initialize fixfac.
c-----------------------------------------------------------------------
      IF(.NOT. test)THEN
         fixfac=0
         DO isol=1,msol
            fixfac(isol,isol)=1
         ENDDO
c-----------------------------------------------------------------------
c     sort unorm.
c-----------------------------------------------------------------------
         index(1:msol)=(/(ipert,ipert=1,msol)/)
         CALL bubble(unorm,index,1,mpert)
         CALL bubble(unorm,index,mpert+1,msol)
      ENDIF
c-----------------------------------------------------------------------
c     triangularize primary solutions.
c-----------------------------------------------------------------------
      mask=.TRUE.
      DO isol=1,mpert
         ksol=index(isol)
         mask(2,ksol)=.FALSE.
         IF(.NOT. test)THEN
            jmax=MAXLOC(ABS(u(:,ksol,1)),mask(1,1:mpert))
            kpert=jmax(1)
            mask(1,kpert)=.FALSE.
         ENDIF
         DO jsol=1,msol
            IF(mask(2,jsol))THEN
               IF(.NOT. test)fixfac(ksol,jsol)
     $              =-u(kpert,jsol,1)/u(kpert,ksol,1)
               u(:,jsol,:)=u(:,jsol,:)+u(:,ksol,:)*fixfac(ksol,jsol)
               IF(.NOT. test)u(kpert,jsol,1)=0
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     triangularize secondary solutions.
c-----------------------------------------------------------------------
      IF(msol > mpert .AND. secondary)THEN
         mask=.TRUE.
         DO isol=mpert+1,msol
            ksol=index(isol)
            mask(2,ksol)=.FALSE.
            IF(.NOT. test)THEN
               jmax=MAXLOC(ABS(u(:,ksol,2)),mask(1,1:mpert))
               kpert=jmax(1)
               mask(1,kpert)=.FALSE.
            ENDIF
            DO jsol=mpert+1,msol
               IF(mask(2,jsol))THEN
                  IF(.NOT. test)fixfac(ksol,jsol)
     $                 =-u(kpert,jsol,2)/u(kpert,ksol,2)
                  u(:,jsol,:)=u(:,jsol,:)+u(:,ksol,:)*fixfac(ksol,jsol)
                  IF(.NOT. test)u(kpert,jsol,2)=0
               ENDIF
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     save fixfac to file.
c-----------------------------------------------------------------------
      IF(bin_euler. AND. .NOT. test)THEN
         WRITE(euler_bin_unit)2
         WRITE(euler_bin_unit)sing_flag,msol
         WRITE(euler_bin_unit)fixfac,index(1:msol)
      ENDIF
c-----------------------------------------------------------------------
c     write output values and close output file.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         WRITE(init_out_unit,'(/a/)')"output values:"
         DO isol=1,msol
            WRITE(init_out_unit,40)"isol = ",isol,", m = ",mlow+isol-1
            WRITE(init_out_unit,50)(ipert,ipert,ipert=1,2)
            WRITE(init_out_unit,60)(ipert,u(ipert,isol,:),ipert=1,mpert)
            WRITE(init_out_unit,50)(ipert,ipert,ipert=1,2)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_fixup
c-----------------------------------------------------------------------
c     subprogram 9. ode_test.
c     tests for optimum approach to singular surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     delcarations.
c-----------------------------------------------------------------------
      FUNCTION ode_test() RESULT(flag)

      LOGICAL :: flag

      INTEGER :: isol
      REAL(r8), SAVE :: singfac_old,powmax
      REAL(r8) :: dsingfac,norm,dnorm,powmax_old
      REAL(r8), DIMENSION(msol) :: power
      COMPLEX(r8), DIMENSION(mpert,msol,2) :: dca
c-----------------------------------------------------------------------
c     simple return.
c-----------------------------------------------------------------------
      flag = psifac == psimax .OR. istep == nstep .OR. istate < 0
      IF(.NOT. res_flag .OR. flag
     $     .OR. ising > msing .OR. singfac > singfac_max)RETURN
c-----------------------------------------------------------------------
c     update values.
c-----------------------------------------------------------------------
      ca_old=ca
      CALL sing_get_ca(ising,psifac,u,ca)
      dca=ca-ca_old
      dsingfac=ABS((singfac-singfac_old)/singfac)
      singfac_old=singfac
c-----------------------------------------------------------------------
c     update powers.
c-----------------------------------------------------------------------
      powmax_old=powmax
      DO isol=1,msol
         norm=ABS(SUM(CONJG(ca(:,isol,:))*ca(:,isol,:)))
         dnorm=ABS(SUM(CONJG(ca(:,isol,:))*dca(:,isol,:)))
         power(isol)=dnorm/(norm*dsingfac)
      ENDDO
      powmax=MAXVAL(power)
c-----------------------------------------------------------------------
c     compute derivatives and flag.
c-----------------------------------------------------------------------
      flag_count=flag_count+1
      IF(flag_count < 3)RETURN
      flag=flag .OR. singfac < singfac_max .AND. powmax > powmax_old
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(ising == ksing)THEN
         WRITE(err_unit)REAL(LOG10(singfac),4),REAL(LOG10(powmax),4)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION ode_test
c-----------------------------------------------------------------------
c     subprogram 10. ode_test_fixup.
c     tests fixup routine.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_test_fixup

      CHARACTER(16)filename
      INTEGER :: isol,ipert
      INTEGER, SAVE :: ifix=0
      REAL(r8), DIMENSION(:,:,:), POINTER :: re,im
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: u_old
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/2x,"isol",2x,"index"/)
 20   FORMAT(2i6)
 30   FORMAT(2(a,i3))
 40   FORMAT(/1x,"i",2(3x,"re u(",i1,")",4x,"im u(",i1,")",1x)/)
 50   FORMAT(i2,1p,4e11.3)
c-----------------------------------------------------------------------
c     create initial u.
c-----------------------------------------------------------------------
      ALLOCATE(u_old(mpert,msol,2),re(mpert,msol,2),im(mpert,msol,2))
      u_old=u
      CALL put_seed(0)
      CALL RANDOM_NUMBER(re)
      CALL RANDOM_NUMBER(im)
      u=CMPLX(re,im)
      DEALLOCATE(re,im)
c-----------------------------------------------------------------------
c     diagnose index.
c-----------------------------------------------------------------------
      ifix=ifix+1
      WRITE(filename,'(a,i2.2,a)')"fixup_",ifix,".out"
      CALL ascii_open(debug_unit,TRIM(filename),"UNKNOWN")
      WRITE(debug_unit,10)
      DO isol=1,msol
         WRITE(debug_unit,20)isol,index(isol)
      ENDDO
      WRITE(debug_unit,10)
c-----------------------------------------------------------------------
c     diagnose initial u.
c-----------------------------------------------------------------------
      WRITE(debug_unit,'(/a/)')"initial values:"
      DO isol=1,msol
         WRITE(debug_unit,30)"isol = ",isol,", m = ",mlow+isol-1
         WRITE(debug_unit,40)(ipert,ipert,ipert=1,2)
         WRITE(debug_unit,50)(ipert,u(ipert,isol,:),ipert=1,mpert)
         WRITE(debug_unit,40)(ipert,ipert,ipert=1,2)
      ENDDO
c-----------------------------------------------------------------------
c     diagnose final u.
c-----------------------------------------------------------------------
      CALL ode_fixup(.FALSE.,.TRUE.)
      WRITE(debug_unit,'(/a/)')"final values:"
      DO isol=1,msol
         WRITE(debug_unit,30)"isol = ",isol,", m = ",mlow+isol-1
         WRITE(debug_unit,40)(ipert,ipert,ipert=1,2)
         WRITE(debug_unit,50)(ipert,u(ipert,isol,:),ipert=1,mpert)
         WRITE(debug_unit,40)(ipert,ipert,ipert=1,2)
      ENDDO
      CALL ascii_close(debug_unit)
      u=u_old
      DEALLOCATE(u_old)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_test_fixup
      END MODULE ode_mod
