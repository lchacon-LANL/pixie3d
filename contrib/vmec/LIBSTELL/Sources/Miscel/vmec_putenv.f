      SUBROUTINE vmec_putenv(ename, evalue, ierror)
      IMPLICIT NONE
      INTEGER :: ierror
      CHARACTER(LEN=*), INTENT(in) :: ename, evalue
      CHARACTER(LEN=len_trim(ename)+len_TRIM(evalue)+2) :: temp
      INTEGER, EXTERNAL :: putenv

      temp = TRIM(ename) // "=" // TRIM(evalue) // char(0)

      ierror = putenv(temp)

      END SUBROUTINE vmec_putenv
