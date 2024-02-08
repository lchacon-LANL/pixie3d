      MODULE safe_open_mod
!
!     Module for performing a "safe" open of a file for
!     a Fortran read/write operation. Makes sure the requested file
!     unit number is not in use, and increments it until an unused
!     unit is found
!
      CONTAINS

      SUBROUTINE safe_open(iunit, istat, filename, filestat,                   &
     &           fileform, record_in, access_in)
!
!     Subroutine for performing a "safe" open of a file for
!     a Fortran read/write operation. Makes sure the requested file
!     unit number is not in use, and increments it until an unused
!     unit is found
!
!  Note that:
! 1)  the actual i/o unit number used is returned in the first argument.
! 2)  the status variable from the OPEN command is returned as the second
!     argument.

!  Here are some examples of usage:
!   To open an existing namelist input file:
!      CALL safe_open(iou,istat,nli_file_name,'old','formatted')
!   To create a file, in order to write to it:
!      CALL safe_open(iou,istat,my_output_file_name,'replace','formatted')

!  JDH 08-30-2004. 
!     From Steve Hirshman's LIBSTELL/Sources/Modules/safe_open_mod.f
!     Rearranged comments, continuation lines, some statement ordering.
!     Commented out duplicate coding. Added variable lscratch, for clarity.
!     Should be NO change in functionality.

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(inout) :: iunit
      INTEGER, INTENT(out) :: istat
      CHARACTER(LEN=*), INTENT(in) :: filename, filestat, fileform
      INTEGER, INTENT(in), OPTIONAL :: record_in
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: access_in

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: cdelim = "apostrophe",                    & 
     &     cform="formatted", cunform="unformatted",                           & 
     &     cscratch="scratch", cseq="sequential"
      CHARACTER(LEN=10) :: acc_type
      LOGICAL :: lopen, lexist, linvalid, lscratch

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!
!     Check that unit is not already opened
!     Increment iunit until find one that is not in use
!
      iunit = MAX(iunit,11)
      linvalid = .true.
      DO WHILE (linvalid)
         INQUIRE(iunit, exist=lexist, opened=lopen, iostat=istat)
         linvalid = (istat.ne.0 .or. .not.lexist) .or. lopen
         IF (.not.linvalid) EXIT
         iunit = iunit + 1
      END DO

      IF (PRESENT(access_in)) THEN
         acc_type = TRIM(access_in)
      ELSE
         acc_type = 'SEQUENTIAL'
      END IF

!  JDH 08-24-2004. From here to end, changed variable lexist to lscratch
      lscratch = (filestat(1:1).eq.'s') .or. (filestat(1:1).eq.'S')      !Scratch file

! Here are the actual OPEN commands. Eight different cases.
      SELECT CASE (fileform(1:1))
      CASE ('u', 'U')
         IF (PRESENT(record_in)) THEN
            IF (lscratch) THEN   ! unformatted, record length specified, scratch 
               OPEN(unit=iunit, form=cunform, status=cscratch,                 &
     &              recl=record_in, access=acc_type, iostat=istat)
            ELSE             ! unformatted, record length specified, non-scratch 
               OPEN(unit=iunit, file=TRIM(filename), form=cunform,             &
     &              status=TRIM(filestat), recl=record_in,                     &
     &              access=acc_type, iostat=istat)
            END IF
         ELSE
            IF (lscratch) THEN ! unformatted, record length unspecified, scratch 
               OPEN(unit=iunit, form=cunform, status=cscratch,                 &
     &              access=acc_type, iostat=istat)
            ELSE           ! unformatted, record length unspecified, non-scratch 
               OPEN(unit=iunit, file=TRIM(filename), form=cunform,             &
     &              status=TRIM(filestat), access=acc_type,iostat=istat)
            END IF
         END IF

      CASE DEFAULT
         IF (PRESENT(record_in)) THEN
            IF (lscratch) THEN     ! formatted, record length specified, scratch 
               OPEN(unit=iunit, form=cform, status=cscratch,                   &
     &              delim=TRIM(cdelim), recl=record_in, access=acc_type,       &
     &              iostat=istat)
            ELSE               ! formatted, record length specified, non-scratch 
               OPEN(unit=iunit, file=TRIM(filename), form=cform,               &
     &              status=TRIM(filestat), delim=TRIM(cdelim),                 &
     &              recl=record_in, access=acc_type, iostat=istat)
            END IF
         ELSE
            IF (lscratch) THEN   ! formatted, record length unspecified, scratch 
               OPEN(unit=iunit, form=cform, status=cscratch,                   &
     &              delim=TRIM(cdelim), access=acc_type, iostat=istat)
            ELSE             ! formatted, record length unspecified, non-scratch 
               OPEN(unit=iunit, file=TRIM(filename), form=cform,               &
     &             status=TRIM(filestat), delim=TRIM(cdelim),                  &
     &             access=acc_type, iostat=istat)
            END IF
         END IF

      END SELECT

      END SUBROUTINE safe_open

      END MODULE safe_open_mod
