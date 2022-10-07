! A module for system call.
! This module is a wrapper for non-standard subroutines and functions.
! 

! Select your compiler.
 #define INTEL_FORTRAN
!#define GNU_FORTRAN

! If a subroutine 'execute_command_line' implemented in your fortran,
! a macro FORTRAN2008_OR_LATER is defined.
! If 'execute_command_line' is not implemented, the function 'system' is used instead. 
#define FORTRAN2008_OR_LATER

module systemcall
#ifdef INTEL_FORTRAN
  use ifport
#endif !INTEL_FORTRAN
  private
  public :: systemcall_unlink, systemcall_rename, systemcall_command, systemcall_symlink
!!$  public :: systemcall_command_safe
contains
  !---------------------------------------------------------------------
  ! Remove/unlink/delete file
  !---------------------------------------------------------------------
  subroutine systemcall_unlink(file, status)
    character(len=*),intent(IN) :: file
    integer,intent(OUT),optional :: status
    integer :: stat
    stat = unlink(file)
    if (stat /= 0) print *, '*** error when unlink ', trim(file)
    if (present(status)) status = stat
  end subroutine systemcall_unlink
  !---------------------------------------------------------------------
  ! Move/rename file
  !---------------------------------------------------------------------
  subroutine systemcall_rename(from, to , status)
    character(len=*),intent(IN) :: from, to
    integer,intent(OUT),optional :: status
    integer :: stat
    stat = rename(from, to)
    if (stat /= 0) print *, '*** error when rename ', trim(from), ' ', trim(to)
    if (present(status)) status = stat
  end subroutine systemcall_rename
  !---------------------------------------------------------------------
  ! Creat symblic link of 'from' to 'to'
  !---------------------------------------------------------------------
  subroutine systemcall_symlink(from, to , status)
    character(len=*),intent(IN) :: from, to
    integer,intent(OUT),optional :: status
    integer :: stat
#ifdef GNU_FORTRAN
    call symlnk(from, to, stat)
#else !GNU_FORTRAN
    call systemcall_command("/bin/ln -s " // trim(from) // " " // trim(to), stat)
#endif !GNU_FORTRAN
    if (stat /= 0) print *, '*** error when rename ', trim(from), ' ', trim(to)
    if (present(status)) status = stat
  end subroutine systemcall_symlink
  !---------------------------------------------------------------------
  ! System call for a command
  !---------------------------------------------------------------------
  subroutine systemcall_command(command, status)
    character(len=*),intent(IN) :: command
    integer,intent(OUT),optional :: status
    integer :: cmdstat
    character(200) :: message
#ifdef FORTRAN2008_OR_LATER
    ! execute_command_line is standard for Fortran2008 or later
    call execute_command_line(command, cmdstat=cmdstat, cmdmsg=message)
    if (cmdstat /= 0) then
       print *, '*** error when system ', trim(command), ' ' ,cmdstat
       print *, 'cmdmsg', message
    endif
#else !FORTRAN2008_OR_LATER
    cmdstat = system(command)
    if (cmdstat /= 0) then
       print *, '*** error when system ', trim(command), ' ' ,cmdstat
    endif
#endif !FORTRAN2008_OR_LATER
    if (present(status)) status = cmdstat
  end subroutine systemcall_command
!!$  !---------------------------------------------------------------------
!!$  ! System call for a command (safe)
!!$  !---------------------------------------------------------------------
!!$  subroutine systemcall_command_safe(command, status)
!!$    character(len=*),intent(IN) :: command
!!$    integer,intent(OUT),optional :: status
!!$    integer :: cmdstat
!!$    integer,parameter :: unit=11
!!$    character(len=*),parameter :: shell_cmd = "syscmd.sh"
!!$    open(unit, file=shell_cmd)
!!$    write(unit, '(A)') command
!!$    write(*, '(A)') command
!!$    call flush(unit)
!!$    close(unit)
!!$    call systemcall_command("sh ./"//shell_cmd, cmdstat)
!!$    call systemcall_unlink(shell_cmd)
!!$    if (present(status)) status = cmdstat
!!$  end subroutine systemcall_command_safe
end module systemcall
! program main
!   use systemcall
!   integer :: status
!   call systemcall_command('echo hello')
!   call systemcall_command('touch hoge.txt')
!   call systemcall_command('mv hoge.txt hoo.txt')
!   call systemcall_unlink('ho.txt', status)
!   call systemcall_rename('hote.txt','hoge.txt', status)
!   print *, status
! end program main
