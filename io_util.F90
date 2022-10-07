#include "config.h"
!-------------------------------------------------------------------------
! Module for utilities of input and output
!-------------------------------------------------------------------------
module io_util
  use string, only : CHARLEN
  implicit none
  private
  ! generic interface for a function readenv
  ! 
  ! bool = readenv(envkey, envvalue [,thead_safe])
  ! 
  ! Arguments:
  !   envkey   = The name of the environment variable. (type:character)
  !   envvalue = The value of the environment variable. (type: character, integer, or real(double))
  !   thread_safe = If 'thread_safe' is set and it is .TRUE., the
  !      envrinment varialble is read only at the primary rank, and
  !      the value is broadcasted to all the rank.  This routine
  !      becomes thread-safe.
  ! 
  ! Return value:
  !   If environment variable is obtaned, this function returns .TRUE.. 
  !   Otherwise, it returns .FALSE.
  ! 
  interface readenv
     module procedure readenv_char, readenv_int, readenv_double
  end interface readenv

  ! generic interface for a subroutine read_env
  ! 
  ! call read_env(envkey, envvalue [,thead_safe] [,stat])
  ! 
  ! Arguments:
  !   envkey   = The name of the environment variable. (type:character)
  !   envvalue = The value of the environment variable. (type: character, integer, or real(double))
  !   thread_safe = If 'thread_safe' is set and it is .TRUE., the
  !      envrinment varialble is read only at the primary rank, and
  !      the value is broadcasted to all the rank.  This routine
  !      becomes thread-safe.
  !  stat = If environment variable is obtaned, stat is set at .TRUE..
  !      Otherwise, it is set at .FALSE..
  ! 
  interface read_env
     module procedure read_env_char, read_env_int, read_env_double
  end interface read_env

  public :: readenv, read_env, wchar, print_msg
contains
  ! --------------------------------------------------------------------
  ! read environment variable (character). return ture if success and false if fail.
  ! --------------------------------------------------------------------
  function readenv_char(envkey, envval, thread_safe) result(bool)
    use mpilib
    character(len=*),intent(in) :: envkey
    character(len=*),intent(out) :: envval
    logical,optional,intent(in) :: thread_safe
    logical :: bool
    integer :: length, stat
    bool = .FALSE.

    call get_environment_variable(envkey, envval, length, stat, .true.)    
    if (stat == 0 .and. length /= 0) bool = .TRUE. ! read val w/o error

    if (present(thread_safe)) then
       if (thread_safe) then
          call mpi_bcast( bool, 1, MPI_LOGICAL, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
          call mpi_bcast( envval, len(envval), MPI_CHARACTER, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
       end if
    end if
  end function readenv_char
  ! --------------------------------------------------------------------
  ! read environment variable (integer). return ture if success and false if fail.
  ! --------------------------------------------------------------------
  function readenv_int(envkey, envval, thread_safe) result(bool)
    use string
    character(len=*),intent(in) :: envkey
    integer,intent(out) :: envval
    logical,optional,intent(in) :: thread_safe
    character(len=CHARLEN) :: envvalc
    logical :: bool
    bool =  readenv_char(envkey, envvalc, thread_safe)
    if (bool) read(envvalc,*) envval
  end function readenv_int
  ! --------------------------------------------------------------------
  ! read environment variable (double). return ture if success and false if fail.
  ! --------------------------------------------------------------------
  function readenv_double(envkey, envval, thread_safe) result(bool)
    use string
    character(len=*),intent(in) :: envkey
    real(kind=DBL_KIND),intent(out) :: envval
    logical,optional,intent(in) :: thread_safe
    character(len=CHARLEN) :: envvalc
    logical :: bool
    bool =  readenv_char(envkey, envvalc, thread_safe)
    if (bool) read(envvalc,*) envval
  end function readenv_double
  ! --------------------------------------------------------------------
  ! read_env for character
  ! --------------------------------------------------------------------
  subroutine read_env_char(envkey, envval, thread_safe, stat)
    character(len=*),intent(in) :: envkey
    character(len=*),intent(out) :: envval
    logical,optional,intent(in) :: thread_safe
    logical,optional,intent(out) :: stat
    logical :: bool
    bool = readenv_char(envkey, envval, thread_safe)
    if (present(stat)) stat = bool
    if (.not. bool) print *, '*** error in read_env ', trim(envkey), ' ', trim(envval)
  end subroutine read_env_char
  ! --------------------------------------------------------------------
  ! read_env for integer
  ! --------------------------------------------------------------------
  subroutine read_env_int(envkey, envval, thread_safe, stat)
    character(len=*),intent(in) :: envkey
    integer,intent(out) :: envval
    logical,optional,intent(in) :: thread_safe
    logical,optional,intent(out) :: stat
    logical :: bool
    bool = readenv_int(envkey, envval, thread_safe)
    if (present(stat)) stat = bool
    if (.not. bool) print *, '*** error in read_env ', trim(envkey), ' ', envval
  end subroutine read_env_int
  ! --------------------------------------------------------------------
  ! read_env for double
  ! --------------------------------------------------------------------
  subroutine read_env_double(envkey, envval, thread_safe, stat)
    character(len=*),intent(in) :: envkey
    real(kind=DBL_KIND),intent(out) :: envval
    logical,optional,intent(in) :: thread_safe
    logical,optional,intent(out) :: stat
    logical :: bool
    bool = readenv_double(envkey, envval, thread_safe)
    if (present(stat)) stat = bool
    if (.not. bool) print *, '*** error in read_env ', trim(envkey), ' ', envval
  end subroutine read_env_double
  ! --------------------------------------------------------------------
  ! write string without left and right white spaces
  ! --------------------------------------------------------------------
  subroutine wchar(ifh, cline)
    character(*) :: cline
    integer,intent(in) :: ifh !file handle
    write(ifh,'(A)') trim(adjustl(cline))
  end subroutine wchar
  ! --------------------------------------------------------------------
  ! print messave
  ! --------------------------------------------------------------------
  subroutine print_msg(cline)
    use mpilib
    character(*) :: cline
    if ( get_myrank() == PRIMARY_RANK )  &
         call wchar(6,cline)
  end subroutine print_msg
end module io_util
