program main
  use io
  use refine
  use multi_timestep
  use writeSnap
  use sinkParticle
  call init
  call restoredata 
  call sp_read
  call step_all_level
  call dumpdata
  call sp_write
  call finalize
end program main
subroutine init
  use mpilib
  use grid
  use parameter
  use modelParameter
  use unit
  use kinzoku
  call mpilib_init
  call hellomsg
  call grid_init
  call parameter_init
  call unit_init
  call init_kinzoku
  call modelParameter_init
  call eos_init
end subroutine init
subroutine hellomsg
  use mpilib
  if (get_myrank() == 0) then
     print *, 'SFUMATO (C) Tomoaki Matsumoto 2004-NOW'
     call flush(6)
  endif
end subroutine hellomsg
subroutine finalize
  use mpilib
  use grid
  call grid_finalize
  call mpi_finalize(ierr)
end subroutine finalize
