module io_util
  use string, only : CHARLEN
  implicit none
  private
  interface readenv
     module procedure readenv_char, readenv_int, readenv_double
  end interface readenv
  interface read_env
     module procedure read_env_char, read_env_int, read_env_double
  end interface read_env
  public :: readenv, read_env, wchar, print_msg
contains
  function readenv_char(envkey, envval, thread_safe) result(bool)
    use mpilib
    character(len=*),intent(in) :: envkey
    character(len=*),intent(out) :: envval
    logical,optional,intent(in) :: thread_safe
    logical :: bool
    integer :: length, stat
    bool = .FALSE.
    call get_environment_variable(envkey, envval, length, stat, .true.)
    if (stat == 0 .and. length /= 0) bool = .TRUE. 
    if (present(thread_safe)) then
       if (thread_safe) then
          call mpi_bcast( bool, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
          call mpi_bcast( envval, len(envval), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
       end if
    end if
  end function readenv_char
  function readenv_int(envkey, envval, thread_safe) result(bool)
    use string
    character(len=*),intent(in) :: envkey
    integer,intent(out) :: envval
    logical,optional,intent(in) :: thread_safe
    character(len=CHARLEN) :: envvalc
    logical :: bool
    bool = readenv_char(envkey, envvalc, thread_safe)
    if (bool) read(envvalc,*) envval
  end function readenv_int
  function readenv_double(envkey, envval, thread_safe) result(bool)
    use string
    character(len=*),intent(in) :: envkey
    real(kind=8),intent(out) :: envval
    logical,optional,intent(in) :: thread_safe
    character(len=CHARLEN) :: envvalc
    logical :: bool
    bool = readenv_char(envkey, envvalc, thread_safe)
    if (bool) read(envvalc,*) envval
  end function readenv_double
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
  subroutine read_env_double(envkey, envval, thread_safe, stat)
    character(len=*),intent(in) :: envkey
    real(kind=8),intent(out) :: envval
    logical,optional,intent(in) :: thread_safe
    logical,optional,intent(out) :: stat
    logical :: bool
    bool = readenv_double(envkey, envval, thread_safe)
    if (present(stat)) stat = bool
    if (.not. bool) print *, '*** error in read_env ', trim(envkey), ' ', envval
  end subroutine read_env_double
  subroutine wchar(ifh, cline)
    character(*) :: cline
    integer,intent(in) :: ifh 
    write(ifh,'(A)') trim(adjustl(cline))
  end subroutine wchar
  subroutine print_msg(cline)
    use mpilib
    character(*) :: cline
    if ( get_myrank() == 0 ) &
         call wchar(6,cline)
  end subroutine print_msg
end module io_util
