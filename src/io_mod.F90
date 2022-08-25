! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for time-whimy things
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module time_io

    ! variable declaration
    implicit none
    ! variable types
    integer, parameter :: long=8
    ! cpu start time
    integer(long), private :: time_start
    ! strcut for date and time information at start
    integer, dimension(8), private :: time_info_start

contains
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialize time module
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_runtime
        ! variable declarations
        implicit none

        ! get time at start up
        call system_clock(time_start)
        ! get date and time information
        call date_and_time(values=time_info_start)

    end subroutine init_runtime

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate runtime in s
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function calculate_runtime() result(runtime)
        ! variable declaration
        implicit none
        ! current cpu time
        integer(long) :: time_end
        ! ticks per second
        integer(long) :: count_rate
        ! time difference
        integer(long) :: runtime

        ! get current cpu time
        call system_clock(time_end, count_rate)
        ! calculate difference
        runtime = (time_end - time_start) / count_rate

    end function calculate_runtime

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! function to get start date/time
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function get_start_datetime() result(start_date)

        ! variable declaration
        implicit none

        ! string to containd date time
        character(len=19) :: start_date

        ! format start date/time
        write (start_date, '(i4,5(a,i2.2))') time_info_start(1), '-', time_info_start(2), '-', time_info_start(3), &
            ' ', time_info_start(5), ':', time_info_start(6), ':', time_info_start(7)

    end function get_start_datetime

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! function to get current time
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_datetime() result(datetime)

        ! variable declaration
        implicit none

        ! string to containd date time
        character(len=19) :: datetime

        ! format start date/time
        write (datetime, '(i4,5(a,i2.2))') time_info_start(1), '-', time_info_start(2), '-', time_info_start(3), &
            ' ', time_info_start(5), ':', time_info_start(6), ':', time_info_start(7)

    end function get_datetime

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! function to get runtime
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_runtime() result(runtime_str)

        ! variable declaration
        implicit none

        ! return string with runtime
        character(len=8) :: runtime_str
        ! runtime in seconds
        integer(long) :: runtime
        ! runtime info
        integer(long), dimension(3) :: runtime_info

        ! get runtime
        runtime = int(calculate_runtime(), long)
        ! calculate second
        runtime_info(3) = mod(runtime, 60_long)
        ! subtract seconds
        runtime = runtime - runtime_info(3)
        ! calculate mins
        runtime_info(2) = mod(runtime, 3600_long) / 60
        ! subtract mins
        runtime = runtime - 60 * runtime_info(2)
        ! calculate hours
        runtime_info(1) = runtime / 3600

        ! format runtime
        write (runtime_str, '(i2,a,i2.2,a,i2.2)') runtime_info(1), ':', runtime_info(2), ':', runtime_info(3)

    end function get_runtime

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write time information
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_time_header(unit, prefix, current_info)
        ! variable declaration
        implicit none

        ! parameter unit
        integer, intent(in), optional :: unit
        ! prefix character argument
        character(len=1), intent(in), optional :: prefix
        ! write current time information and runtime
        logical, intent(in), optional :: current_info

        ! io device unit
        integer :: p_unit
        ! prefix character
        character(len=2) :: p_prefix

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! get optional unit argument
        if (present(unit)) then
            ! assign file unit
            p_unit = unit
            ! nothing passed
        else ! p_unit
            ! default terminal buffer
            p_unit = 6
        endif ! p_unit

        ! get optional prefix
        if (present(prefix)) then
            ! assign prefix
            p_prefix = prefix
        else ! p_prefix
            ! default prefix
            p_prefix = '  '
        endif ! p_prefix

        ! write start time
        write (p_unit, '(a,a,a)') adjustl(p_prefix), 'Run started on ', get_start_datetime()
        ! check if flag is present
        if (present(current_info)) then
            ! check flag
            if (current_info) then

                ! write current info
                write (p_unit, '(a,a,a)') adjustl(p_prefix), 'Current time: ', get_datetime()
                ! write runtime
                write (p_unit, '(a,a,a)') adjustl(p_prefix), 'Runtime: ', get_runtime()

            endif ! flag
        endif ! present

    end subroutine write_time_header

end module time_io

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for for exception handling
!
! Allows the unit testing system to check for errors
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pixie3d_exception

    ! variable declaration
    implicit none

    ! throw function interface
    abstract interface
        subroutine exception_handler(function_name, message)
            ! name of calling function
            character(len=*), intent(in) :: function_name
            ! message for error
            character(len=*), intent(in) :: message
        end subroutine exception_handler
    end interface

    ! function pointer to handler function, use to raise and exception
    procedure(exception_handler), pointer, protected :: raise => default_exception_handler

contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! default exception handler
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine default_exception_handler(fct, msg_str)

        ! error handling
        use grid_mpi, only: pstop

        ! variable declaration
        implicit none

        ! function name
        character(len=*), intent(in) :: fct
        ! message string
        character(len=*), intent(in) :: msg_str

        ! call framework routine to abort
        call pstop(fct, msg_str)

    end subroutine default_exception_handler

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! change default exception halder
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine set_exception_handler(method)

        ! variable declaration
        implicit none

        ! pointer to handler method
        procedure(exception_handler), optional :: method

        if (present(method)) then
            ! pointer handler to new method
            raise => method
        else
            ! reset pointer to default
            raise => default_exception_handler
        end if ! present

    end subroutine set_exception_handler

end module pixie3d_exception

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! functions for compile time information
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module compile_info

    ! check for the definition of the macros and set default values if necessary
    ! git branch
#ifndef GIT_BRANCH
#define GIT_BRANCH "Unknown"
#endif

    ! git commit date
#ifndef GIT_COMMIT_DATE
#define GIT_COMMIT_DATE "Unknown"
#endif

    ! git commit hash
#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "Unknown"
#endif

#ifndef COMPILE_HOST
#define COMPILE_HOST "Unknown"
#endif

contains
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compile date
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function compile_date() result(date)
        ! variable declaration
        implicit none
        ! return string
        character(len=11) :: date

        ! assign compile date to return string
        date = __DATE__
    end function compile_date

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compile time
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function compile_time() result(time)
        ! variable declaration
        implicit none
        ! return string
        character(len=8) :: time

        ! assign compile time to return value
        time = __TIME__
    end function compile_time

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compiler name
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function compiler_name() result(compiler)
        ! variable declaration
        implicit none
        ! return string
        character(len=10) :: compiler

        ! assign compiler name to return value, depending on preprocessor values
#if __INTEL_COMPILER
        compiler = "intel"
#elif _CRAYFTN
        compiler = "cray"
#elif __GNUC__
        compiler = "gfortran"
#else
        compiler = "Unknown"
#endif

    end function compiler_name

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compiler version
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function compiler_version() result(version)
        ! variable declaration
        implicit none
        ! return string
        character(len=14) :: version

        ! assign compiler version to return value, depending on preprocessor values
#if  __INTEL_COMPILER
        write (version, '(i4)') __INTEL_COMPILER
#elif _CRAYFTN
        write (version, '(i0,a,i0,a,i0)') _RELEASE_MAJOR,".",_RELEASE_MINOR,".",_RELEASE_PATCHLEVEL
#elif __GNUC__
        version = __VERSION__
#else
        version = "Unknown"
#endif

    end function compiler_version

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! git branch
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function git_branch() result(branch)
        ! variable declaration
        implicit none
        ! return string
        character(len=30) :: branch

        ! assign the git branch name
        branch = GIT_BRANCH

    end function git_branch

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! git commit hash
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function git_commit_hash() result(hash)
        ! variable declaration
        implicit none
        ! return string
        character(len=40) :: hash

        ! assign commit hash to return value
        hash = GIT_COMMIT_HASH

    end function git_commit_hash

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! git commit date
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function git_commit_date() result(date)
        ! variable declaration
        implicit none
        ! return string
        character(len=25) :: date

        ! assign commit date to return value
        date = GIT_COMMIT_DATE

    end function git_commit_date

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! host name of machine compiled on
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function compile_host() result(host)
        ! variable declaration
        implicit none
        ! return string
        character(len=100) :: host

        ! assign host name to return
        host = COMPILE_HOST

    end function compile_host

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! print formated information
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_compile_header(unit, prefix)
        ! variable declaration
        implicit none

        ! parameter unit
        integer, intent(in), optional :: unit
        ! io device unit
        integer :: p_unit

        ! prefix character argument
        character(len=1), intent(in), optional :: prefix
        ! prefix character
        character(len=2) :: p_prefix

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! get optional unit argument
        if (present(unit)) then
            ! assign file unit
            p_unit = unit
            ! nothing passed
        else ! p_unit
            ! default terminal buffer
            p_unit = 6
        endif ! p_unit

        ! get optional prefix
        if (present(prefix)) then
            ! assign prefix
            p_prefix = prefix
        else ! p_prefix
            ! default prefix
            p_prefix = '  '
        endif ! p_prefix

        ! write compiler information
        write (p_unit, '(11a)') adjustl(p_prefix), "Compiled on ", compile_date(), " ", compile_time(), " with ", &
            trim(compiler_name()), " ", trim(compiler_version()), ' on ', trim(compile_host())
        ! write git information
        write (p_unit, '(7a)') adjustl(p_prefix), "Git branch: ", trim(git_branch()), &
            " on ", trim(git_commit_date()), "   commit: ", trim(git_commit_hash())

    end subroutine write_compile_header
end module compile_info

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for ouput operations
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module pixie3d_io

    use grid_mpi, ONLY: my_rank
  
    ! variable declaration
    implicit none

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! constants for color codes
    ! when used in front of output, will turn font in terminal to this color
    ! example write (*, *) COLOR_RED, 'error message'
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! excape character used to indicate begin of color code
    character(len=1), parameter :: escape_char = char(27)

    ! default color
    character(len=4), parameter :: COLOR_DEFAULT = escape_char//'[0m'
    ! black
    character(len=5), parameter :: COLOR_BLACK = escape_char//'[30m'
    ! red
    character(len=5), parameter :: COLOR_RED = escape_char//'[31m'
    ! greem
    character(len=5), parameter :: COLOR_GREEN = escape_char//'[32m'
    ! yellow
    character(len=5), parameter :: COLOR_YELLOW = escape_char//'[33m'
    ! blue
    character(len=5), parameter :: COLOR_BLUE = escape_char//'[34m'
    ! purple
    character(len=5), parameter :: COLOR_PURPLE = escape_char//'[35m'
    ! aqua
    character(len=5), parameter :: COLOR_AGUA = escape_char//'[36m'

    ! light gray
    character(len=5), parameter :: COLOR_LIGHT_GRAY = escape_char//'[90m'
    ! peach
    character(len=5), parameter :: COLOR_PEACH = escape_char//'[91m'
    ! light green
    character(len=5), parameter :: COLOR_LIGHT_GREEN = escape_char//'[92m'
    ! light yellow
    character(len=5), parameter :: COLOR_LIGHT_YELLOW = escape_char//'[93m'
    ! light blue
    character(len=5), parameter :: COLOR__LIGHT_BLUE = escape_char//'[94m'
    ! pink
    character(len=5), parameter :: COLOR__PINK = escape_char//'[95m'
    ! light agua
    character(len=5), parameter :: COLOR_LIGHT_AGUA = escape_char//'[96m'
    ! pearl white
    character(len=5), parameter :: COLOR_PEARL_WHITE = escape_char//'[97m'

    ! arrays used to convert string cases
    character(len=*), private, parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
    character(len=*), private, parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  generic interface to allow easy conversion of numbers to strings
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    interface to_str
!!$        ! integer conversion
!!$        module procedure :: short_to_str, &
!!$            int_to_str, &
!!$            long_to_str
!!$        ! real conversion
!!$        module procedure :: real_to_str, &
!!$            double_to_str, &
!!$            quad_to_str
!!$    end interface to_str

contains
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! convert C host name to fortran
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_hostname(hostname) result(error)

        ! C bindings
        use iso_c_binding, only: c_int, c_char

        ! variable declaration
        implicit none

        ! standard conform way to get the host name is to call a C function
        interface
            integer(kind=c_int) function get_c_hostname(hname, len) bind(C, name='gethostname')
                ! C bindings
                use iso_c_binding, only: c_int, c_char
                ! variable declaration
                implicit none
                ! host name
                character(kind=c_char) :: hname(*)
                ! string length
                integer(kind=c_int), VALUE :: len
            end function get_c_hostname
        end interface

        ! host name
        character(len=*), intent(out) :: hostname
        ! error code
        integer :: error
        ! string length
        integer(kind=c_int), parameter :: str_len = 100
        ! C string
        character(kind=c_char) :: c_hostname(str_len)
        ! conversion variable
        character :: c
        ! loop variable
        integer :: i

        ! initialize host name
        hostname = ''
        ! call C function
        error = get_c_hostname(c_hostname, str_len)
        ! check return value
        if (error == 0) then
            ! loop over C strin
            do i = 1, str_len
                ! get current letter
                c = c_hostname(i)
                ! check for end of string
                if (c == char(0)) exit
                ! store in hostname
                hostname(i:i) = c
            end do ! i

            ! remove whitespace
            hostname = trim(hostname)

        end if ! error

    end function get_hostname

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write a message to the terminal, optional settings for color and channel
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message(msg_str, color, error)
        ! standard file units
        use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
        ! MPI parameters
!!$        use pixie3d_mpi, only: mpi_info

        ! variable declaration
        implicit none

        ! message string
        character(len=*), intent(in) :: msg_str
        ! optional color code to be used
        character(len=5), intent(in), optional :: color
        ! flag to use error output
        logical, intent(in), optional :: error

        ! unit for output
        integer :: unit_num

        ! check MPI rank
        if (my_rank == 0) then
            ! check for error flag
            if (.not. present(error)) then
                ! set to default unit
                unit_num = output_unit
            elseif (error) then
                ! set to error unit
                unit_num = error_unit
            else
                ! set to default unit
                unit_num = output_unit
            end if ! error flag

            ! check for color
            if (present(color)) then
                ! write to standard out with color
                write (unit_num, '(a,a,a)') color, trim(msg_str), COLOR_DEFAULT
            else
                ! write to standard out, no color
                write (unit_num, '(a)') trim(msg_str)
            end if ! color

        end if ! MPI rank

    end subroutine message

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write a highlighted text
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine info(msg_str)

        ! variable declaration
        implicit none

        ! message string
        character(len=*), intent(in) :: msg_str

        ! write using color agua
        call message(msg_str, COLOR_AGUA)

    end subroutine info

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write a warning message to the terminal
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine warning(msg_str)

        ! variable declaration
        implicit none

        ! message string
        character(len=*), intent(in) :: msg_str

        ! write using yellow
        call message(msg_str, COLOR_YELLOW)

    end subroutine warning

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write an error message to the terminal
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error(msg_str)

        ! variable declaration
        implicit none

        ! message string
        character(len=*), intent(in) :: msg_str

        ! write using red to error unit
        call message(msg_str, COLOR_RED, .true.)

    end subroutine error

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! critical error, write error message and terminate calculation
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine critical(fct, msg_str)

        ! error handling
        use pixie3d_exception, only: raise

        ! variable declaration
        implicit none

        ! function name
        character(len=*), intent(in) :: fct
        ! message string
        character(len=*), intent(in) :: msg_str

        ! write using red to error unit
        call message(msg_str, COLOR_RED, .true.)
        ! raise an exception
        call raise(fct, msg_str)

    end subroutine critical

!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    ! convert to lower case
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    pure function to_lower(str) result(lower)
!!$
!!$        ! type definitions
!!$        use pixie3d_types, only: standard
!!$
!!$        ! variable declaration
!!$        implicit none
!!$
!!$        ! input string
!!$        character(len=*), intent(in) :: str
!!$        ! lower case string
!!$        character(len=len_trim(str)) :: lower
!!$
!!$        ! loop variable
!!$        integer(standard) :: i
!!$        ! character number
!!$        integer(standard) :: c
!!$
!!$        ! copy string
!!$        lower = str
!!$
!!$        ! loop over string
!!$        do i = 1, len(lower)
!!$            ! find location of letter in upper case array
!!$            c = index(UPPER_CASE, lower(i:i))
!!$            ! if character is upper case
!!$            if (c /= 0) then
!!$                ! und convert to lower case
!!$                lower(i:i) = LOWER_CASE(c:c)
!!$            end if ! c
!!$
!!$        end do ! i
!!$
!!$    end function to_lower
!!$
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    ! convert to upper case
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    pure function to_upper(str) result(upper)
!!$
!!$        ! type definitions
!!$        use pixie3d_types, only: standard
!!$
!!$        ! variable declaration
!!$        implicit none
!!$
!!$        ! input string
!!$        character(len=*), intent(in) :: str
!!$        ! lower case string
!!$        character(len=len_trim(str)) :: upper
!!$
!!$        ! loop variable
!!$        integer(standard) :: i
!!$        ! character number
!!$        integer(standard) :: c
!!$
!!$        ! copy string
!!$        upper = str
!!$
!!$        ! loop over string
!!$        do i = 1, len(upper)
!!$            ! find location of letter in lower case array
!!$            c = index(LOWER_CASE, upper(i:i))
!!$            ! if character is lower case
!!$            if (c /= 0) then
!!$                ! and convert to upper case
!!$                upper(i:i) = UPPER_CASE(c:c)
!!$            end if ! c
!!$
!!$        end do ! i
!!$
!!$    end function to_upper
!!$
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    ! formating an short to string
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    pure function short_to_str(i, fmt_str) result(str)
!!$
!!$        ! type definitions
!!$        use pixie3d_types, only: short
!!$
!!$        ! variable declaration
!!$        implicit none
!!$
!!$        ! integer value to be converted
!!$        integer(short), intent(in) :: i
!!$        ! format string
!!$        character(len=*), intent(in), optional :: fmt_str
!!$
!!$        ! output string
!!$        character(len=:), allocatable :: str
!!$
!!$        ! temporary string
!!$        character(len=50) :: tmp_str
!!$
!!$        ! check for format string
!!$        if (present(fmt_str)) then
!!$            ! format using argument
!!$            write (tmp_str, fmt_str) i
!!$        else
!!$            ! format using default
!!$            write (tmp_str, '(g0)') i
!!$        end if ! present
!!$
!!$        ! allocate return string to appropriate legnth and assign
!!$        ! relies on reallocation
!!$        str = trim(tmp_str)
!!$
!!$    end function short_to_str
!!$
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    ! formating an standard to string
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    pure function int_to_str(i, fmt_str) result(str)
!!$
!!$        ! type definitions
!!$        use pixie3d_types, only: standard
!!$
!!$        ! variable declaration
!!$        implicit none
!!$
!!$        ! integer value to be converted
!!$        integer(standard), intent(in) :: i
!!$        ! format string
!!$        character(len=*), intent(in), optional :: fmt_str
!!$
!!$        ! output string
!!$        character(len=:), allocatable :: str
!!$
!!$        ! temporary string
!!$        character(len=50) :: tmp_str
!!$
!!$        ! check for format string
!!$        if (present(fmt_str)) then
!!$            ! format using argument
!!$            write (tmp_str, fmt_str) i
!!$        else
!!$            ! format using default
!!$            write (tmp_str, '(g0)') i
!!$        end if ! present
!!$
!!$        ! allocate return string to appropriate legnth and assign
!!$        ! relies on reallocation
!!$        str = trim(tmp_str)
!!$
!!$    end function int_to_str
!!$
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    ! formating an long to string
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    pure function long_to_str(i, fmt_str) result(str)
!!$
!!$        ! type definitions
!!$        use pixie3d_types, only: long
!!$
!!$        ! variable declaration
!!$        implicit none
!!$
!!$        ! integer value to be converted
!!$        integer(long), intent(in) :: i
!!$        ! format string
!!$        character(len=*), intent(in), optional :: fmt_str
!!$
!!$        ! output string
!!$        character(len=:), allocatable :: str
!!$
!!$        ! temporary string
!!$        character(len=50) :: tmp_str
!!$
!!$        ! check for format string
!!$        if (present(fmt_str)) then
!!$            ! format using argument
!!$            write (tmp_str, fmt_str) i
!!$        else
!!$            ! format using default
!!$            write (tmp_str, '(g0)') i
!!$        end if ! present
!!$
!!$        ! allocate return string to appropriate legnth and assign
!!$        ! relies on reallocation
!!$        str = trim(tmp_str)
!!$
!!$    end function long_to_str
!!$
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    ! formating a real to string
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    pure function real_to_str(r, fmt_str) result(str)
!!$
!!$        ! type definitions
!!$        use pixie3d_types, only: single
!!$
!!$        ! variable declaration
!!$        implicit none
!!$
!!$        ! integer value to be converted
!!$        real(single), intent(in) :: r
!!$        ! format string
!!$        character(len=*), intent(in), optional :: fmt_str
!!$
!!$        ! output string
!!$        character(len=:), allocatable :: str
!!$
!!$        ! temporary string
!!$        character(len=50) :: tmp_str
!!$
!!$        ! check for format string
!!$        if (present(fmt_str)) then
!!$            ! format using argument
!!$            write (tmp_str, fmt_str) r
!!$        else
!!$            ! format using default
!!$            write (tmp_str, '(g0)') r
!!$        end if ! present
!!$
!!$        ! allocate return string to appropriate legnth and assign
!!$        ! relies on reallocation
!!$        str = trim(tmp_str)
!!$
!!$    end function real_to_str
!!$
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    ! formating a double to string
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    pure function double_to_str(r, fmt_str) result(str)
!!$
!!$        ! type definitions
!!$        use pixie3d_types, only: double
!!$
!!$        ! variable declaration
!!$        implicit none
!!$
!!$        ! integer value to be converted
!!$        real(double), intent(in) :: r
!!$        ! format string
!!$        character(len=*), intent(in), optional :: fmt_str
!!$
!!$        ! output string
!!$        character(len=:), allocatable :: str
!!$
!!$        ! temporary string
!!$        character(len=50) :: tmp_str
!!$
!!$        ! check for format string
!!$        if (present(fmt_str)) then
!!$            ! format using argument
!!$            write (tmp_str, fmt_str) r
!!$        else
!!$            ! format using default
!!$            write (tmp_str, '(g0)') r
!!$        end if ! present
!!$
!!$        ! allocate return string to appropriate legnth and assign
!!$        ! relies on reallocation
!!$        str = trim(tmp_str)
!!$
!!$    end function double_to_str
!!$
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    ! formating a quad to string
!!$    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    pure function quad_to_str(r, fmt_str) result(str)
!!$
!!$        ! type definitions
!!$        use pixie3d_types, only: quad
!!$
!!$        ! variable declaration
!!$        implicit none
!!$
!!$        ! integer value to be converted
!!$        real(quad), intent(in) :: r
!!$        ! format string
!!$        character(len=*), intent(in), optional :: fmt_str
!!$
!!$        ! output string
!!$        character(len=:), allocatable :: str
!!$
!!$        ! temporary string
!!$        character(len=50) :: tmp_str
!!$
!!$        ! check for format string
!!$        if (present(fmt_str)) then
!!$            ! format using argument
!!$            write (tmp_str, fmt_str) r
!!$        else
!!$            ! format using default
!!$            write (tmp_str, '(g0)') r
!!$        end if ! present
!!$
!!$        ! allocate return string to appropriate legnth and assign
!!$        ! relies on reallocation
!!$        str = trim(tmp_str)
!!$
!!$    end function quad_to_str

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write terminal header
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine output_header(unit, prefix)
        ! standard file units
        use, intrinsic :: iso_fortran_env, only: output_unit

        ! get compile info
        use compile_info, only: write_compile_header
        ! get time information
        use time_io, only: write_time_header,init_runtime

        ! variable declaration
        implicit none

        ! parameter unit
        integer, intent(in), optional :: unit
        ! io device unit
        integer :: p_unit

        ! prefix character argument
        character(len=1), intent(in), optional :: prefix
        ! prefix character
        character(len=2) :: p_prefix
        ! host name
        character(len=25) :: hostname
        ! color string
        character(len=5) :: color
        ! default color
        character(len=4) :: reset_color

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call init_runtime

        ! get optional unit argument
        if (present(unit)) then
            ! assign file unit
            p_unit = unit
            ! nothing passed
        else ! p_unit
            ! default terminal buffer
            p_unit = output_unit
        endif ! p_unit

        ! do nothing if terminal and not rank 0
        if (p_unit == output_unit .and. .not. my_rank == 0) then
            ! return focus
            return
        end if ! rank

        ! get optional prefix
        if (present(prefix)) then
            ! assign prefix
            p_prefix = prefix//' '
        else ! p_prefix
            ! default prefix
            p_prefix = '  '
        endif ! p_prefix

        ! only print color if we print to terminal
        if (p_unit == output_unit) then
            ! set color to agua
            color = COLOR_AGUA
            ! set reset string
            reset_color = COLOR_DEFAULT
        else
            ! no color string for log files
            color = ''
            ! no reset string for log files
            reset_color = ''
        end if

        write (p_unit, '(4a)') &
            adjustl(p_prefix), trim(color), '                                                               ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '       _/_/_/    _/  _/    _/   _/  _/_/_/     _/_/_/   _/_/_/ ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '      _/    _/  _/    _/ _/    _/  _/             _/   _/    _/', trim(reset_color), &
            adjustl(p_prefix), trim(color), '     _/_/_/    _/     _/      _/  _/_/_/       _/_/   _/    _/ ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '    _/        _/    _/ _/    _/  _/             _/   _/    _/  ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '   _/        _/  _/     _/  _/  _/_/_/     _/_/_/   _/_/_/     ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                                                               ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'Extended MHD plasma physics code in arbitrary curvilinear geom.', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                                                               ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                                                               ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                   Developer: L. Chacon                        ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '           Los Alamos National Laboratory (12-)                ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '           Oak Ridge  National Laboratory (08-12)              ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '           Los Alamos National Laboratory (04-08)              ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                       LA-CC 07-005                            ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                                                               ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                  with contributions from:                     ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                                                               ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                 Daniele Bonfiglio, CNR-RFX                    ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '              Mark Berrill, Bobby Philip, ORNL                 ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                                                               ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                                                               ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'This program was produced under U.S. Government     ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'contract 89233218CNA000001 for Los Alamos National  ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'Laboratory (LANL), which is operated by Triad       ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'National Security, LLC for the U.S. Department of   ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'Energy/National Nuclear Security Administration.    ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'All rights in the program are reserved by Triad     ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'National Security, LLC, and the U.S. Department of  ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'Energy/NationalNuclear Security Administration.     ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'The Government is granted for itself and others     ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'acting on its behalf a nonexclusive, paid-up,       ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'irrevocable worldwide license in this material      ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'to reproduce, prepare derivative  works,            ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'distribute copies to the public, perform publicly   ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'and display publicly, and to permit others to do so.', trim(reset_color), &
            adjustl(p_prefix), ''

        write (p_unit, '(4a)') &
            adjustl(p_prefix), trim(color), 'This program is open source under the BSD-3 License. ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'Redistribution and use in source and binary forms,   ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'with or without modification, are permitted provided ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'that the following conditions are met:               ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '1. Redistributions of source code must retain the    ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'above copyright notice, this list of conditions and  ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'the following disclaimer.                            ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                                                     ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '2. Redistributions in binary form must reproduce     ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'the above copyright notice, this list of conditions  ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'and the following disclaimer in the documentation    ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'and/or other materials provided with the distribution', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                                                     ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '3. Neither the name of the copyright holder nor the  ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'names of its contributors may be used to endorse or  ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'promote products derived from this software without  ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'specific prior written permission.                   ', trim(reset_color), &
            adjustl(p_prefix), trim(color), '                                                     ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS   ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED  ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE       ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), 'IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS    ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), 'FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), 'SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,       ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), 'EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT  ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), 'NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR   ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), 'SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), 'INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    ', trim(reset_color), &
            adjustl(p_prefix), trim(color), 'LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,    ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), 'OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING  ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), 'IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF  ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), 'ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.           ', trim(reset_color), &      
            adjustl(p_prefix), trim(color), '', trim(reset_color), &      
            adjustl(p_prefix), ''

        ! write compile information
        call write_compile_header(p_unit, p_prefix)

        ! write compile options
        ! test for petsc
#ifdef petsc
        write (p_unit, '(a,t6,a,a,a)') adjustl(p_prefix), trim(color), 'PETSc support', trim(reset_color)
#endif
        ! test for openMP
#ifdef _OPENMP
        write (p_unit, '(a,t6,a,a,a)') adjustl(p_prefix), trim(color), 'OpenMP support', trim(reset_color)
#endif

        ! blanck line
        write (p_unit, '(a)') adjustl(p_prefix)

        ! get host name of running machine
        if (get_hostname(hostname) > 0) then
            hostname = 'Unkown'
        endif
        ! write running host name
        write (p_unit, '(6a)') adjustl(p_prefix), trim(color), 'Running on ', hostname, trim(reset_color)

        ! write run time information
        call write_time_header(p_unit, p_prefix, p_unit /= 6)

        ! add blank lines
        write (p_unit, '(a)') adjustl(p_prefix), adjustl(p_prefix)

    end subroutine output_header

end module pixie3d_io
