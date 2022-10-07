module systemcall
  use ifport
  private
  public :: systemcall_unlink, systemcall_rename, systemcall_command, systemcall_symlink
contains
  subroutine systemcall_unlink(file, status)
    character(len=*),intent(IN) :: file
    integer,intent(OUT),optional :: status
    integer :: stat
    stat = unlink(file)
    if (stat /= 0) print *, '*** error when unlink ', trim(file)
    if (present(status)) status = stat
  end subroutine systemcall_unlink
  subroutine systemcall_rename(from, to , status)
    character(len=*),intent(IN) :: from, to
    integer,intent(OUT),optional :: status
    integer :: stat
    stat = rename(from, to)
    if (stat /= 0) print *, '*** error when rename ', trim(from), ' ', trim(to)
    if (present(status)) status = stat
  end subroutine systemcall_rename
  subroutine systemcall_symlink(from, to , status)
    character(len=*),intent(IN) :: from, to
    integer,intent(OUT),optional :: status
    integer :: stat
    call systemcall_command("/bin/ln -s " // trim(from) // " " // trim(to), stat)
    if (stat /= 0) print *, '*** error when rename ', trim(from), ' ', trim(to)
    if (present(status)) status = stat
  end subroutine systemcall_symlink
  subroutine systemcall_command(command, status)
    character(len=*),intent(IN) :: command
    integer,intent(OUT),optional :: status
    integer :: cmdstat
    character(200) :: message
    call execute_command_line(command, cmdstat=cmdstat, cmdmsg=message)
    if (cmdstat /= 0) then
       print *, '*** error when system ', trim(command), ' ' ,cmdstat
       print *, 'cmdmsg', message
    endif
    if (present(status)) status = cmdstat
  end subroutine systemcall_command
end module systemcall
