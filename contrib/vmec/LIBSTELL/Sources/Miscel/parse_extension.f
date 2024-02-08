      SUBROUTINE parse_extension(file_to_parse, file_or_extension, lnc)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: file_or_extension
      CHARACTER(LEN=*), INTENT(inout) :: file_to_parse
      LOGICAL, INTENT(out) :: lnc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: index_path, index_nc, istat=0
      LOGICAL :: ltxt
      CHARACTER(len=LEN(file_to_parse)) :: path
      CHARACTER(len=LEN(file_to_parse)) :: temp      !!Assumes file_to_parse can store it all
      CHARACTER(LEN=1) :: ch_test
!-----------------------------------------------
!
!     FIRST CHECK IF FILE_OR_EXTENSION IS A FILENAME (FILE_TO_PARSE EMBEDDED) 
!     OR AN EXTENSION
!
      index_path = INDEX(file_or_extension, TRIM(file_to_parse))

      IF (index_path .gt. 0) THEN
!
!     MUST BE <FILENAME). OR <FILENAME>_
!
         index_path = index_path + LEN_TRIM(file_to_parse)
         ch_test = file_or_extension(index_path:index_path)
         IF ((ch_test.ne.'.') .and. (ch_test.ne.'_')) index_path = 0
      END IF

      IF (index_path .gt. 0) THEN
         file_to_parse = file_or_extension
!
!     CHECK FOR netcdf FILE EXTENSION (*.nc)
!
         index_nc = INDEX(file_to_parse,".nc",BACK=.TRUE.)
         lnc = (index_nc .eq. (LEN_TRIM(file_to_parse)-2))
!
!     MAY HAVE PASSED FILE NAME EXTENSION WITHOUT .nc; CHECK IF FILE_TO_PARSE EXISTS
         IF (.not.lnc) THEN
            INQUIRE(FILE=file_to_parse, EXIST=ltxt, iostat=istat)
            IF (istat.ne.0 .or. .not.ltxt) THEN
               file_to_parse = TRIM(file_to_parse) // ".nc"
               lnc = .true.
            END IF
         END IF

      ELSE
!
!     CHECK IF NETCDF FILE EXISTS
!
         path = file_to_parse
         IF (file_or_extension(1:1) == '.' .or. 
     1       file_or_extension(1:1) == '_') THEN
            temp = TRIM(path) // file_or_extension
         ELSE
            temp = TRIM(path) // '_' // file_or_extension
         END IF

!
!     FIRST LOOK FOR FILE WITH .nc SUFFIX IN file_or_extension
!
         file_to_parse = temp
         index_nc = INDEX(file_to_parse,".nc",BACK=.TRUE.)
         lnc = (index_nc .eq. (LEN_TRIM(file_to_parse)-2))
!
!     NEXT LOOK FOR .txt SUFFIX
!
         IF (.not.lnc) THEN
            index_nc = INDEX(file_to_parse,".txt",BACK=.TRUE.)
            ltxt = (index_nc .eq. (LEN_TRIM(file_to_parse)-3))
!
!     CHECK IF file_or_extension WAS GIVEN WITHOUT EXPLICIT .nc OR .txt SUFFIX
!
            IF (.not.ltxt) THEN
               file_to_parse = TRIM(temp) // '.nc'
               INQUIRE (FILE=file_to_parse, EXIST=lnc, iostat=istat)
               IF (istat.ne.0 .or. .not.lnc) THEN
                  file_to_parse = TRIM(path) // '.' 
     1                         // TRIM(file_or_extension) // '.nc'
                  INQUIRE (FILE=file_to_parse, EXIST=lnc, iostat=istat)
                  IF (istat.ne.0 .or. .not.lnc) THEN
                     file_to_parse = TRIM(temp) // '.txt'
                     INQUIRE (FILE=file_to_parse, EXIST=ltxt, 
     1                        iostat=istat)
                     IF (.not.ltxt) THEN
                        file_to_parse = TRIM(path) // '.' // 
     1                               TRIM(file_or_extension) // '.txt'
                        INQUIRE (FILE=file_to_parse, EXIST=ltxt, 
     1                          iostat=istat)
                     END IF
                  END IF
               END IF
            END IF
!
!     DEFAULT (OLD STYLE) FILE NAME WHEN NONE OF THE ABOVE EXIST
!      
            IF ((istat.ne.0 .or. .not.ltxt) .and. .not.lnc) THEN
               file_to_parse = TRIM(path) // '.' // file_or_extension
            END IF

         END IF
 
      END IF
         
      END SUBROUTINE parse_extension
