#include "config.h"
! SFUMATO-AMR
! (c) Tomoaki Matsumoto <matsu@i.hosei.ac.jp>

program main
  use io
  use refine
  use multi_timestep
  use writeSnap
#ifdef SINKPARTICLE
  use sinkParticle
#endif ! SINKPARTICLE
  call init
  call restoredata              ! array should be allocated!!
#ifdef SINKPARTICLE
  call sp_read
#endif !SINKPARTICLE
  call step_all_level
  call dumpdata
#ifdef SINKPARTICLE
  call sp_write
#endif !SINKPARTICLE
  call finalize
end program main
!-------------------------------------------------------------------------
! Initialize program
!-------------------------------------------------------------------------
subroutine init
  use mpilib
!!$  use boundary
  use grid
  use parameter
  use modelParameter
  use unit
#ifdef METAL ! HFADDED
  use kinzoku
#endif
  call mpilib_init
  call hellomsg
  call grid_init
  call parameter_init
  call unit_init
#ifdef METAL ! HFADDED
  call init_kinzoku
#endif
  call modelParameter_init
  call eos_init
!!$  call boundary_init
end subroutine init
!-------------------------------------------------------------------------
subroutine hellomsg
!!$  use systemcall
  use mpilib
  if (get_myrank() == PRIMARY_RANK) then
!!$     call systemcall_command('date')
     print *, 'SFUMATO (C) Tomoaki Matsumoto 2004-NOW'
     call flush(6)
  endif
!!$  call mpi_barrier( MPI_COMM_WORLD, ierr )
end subroutine hellomsg
!-------------------------------------------------------------------------
! Initialize program
!-------------------------------------------------------------------------
subroutine finalize
  use mpilib
  use grid
  call grid_finalize
  call mpi_finalize(ierr)
end subroutine finalize
