!-------------------------------------------------------------------------
! Macro for packarr.F90

!
! example:
!
!    call pkar_reset
!    myrank = get_myrank()
!    PACK_SEND3( 3dim_array, myrank, ranks, rankd) 
!            Note: this command should be executed by both ranks and rankd.
!    call pkar_sendrecv()
!    UNPACK_RECV3(3dim_array, rankds)
!
! When array size fixed only in ranks, use PACK_SEND*_SZ, e.g., 
!
!    PACK_SEND3_SZ( 3dim_array, myrank, ranks, rankd, sizeofarray) 
!
!-------------------------------------------------------------------------
! Scalar variable
!-------------------------------------------------------------------------
#define PACK_SEND0(A, MYRANK, SRC, DST) \
  if ((MYRANK) == (SRC)) call pkar_push(A, 1, kind(A), DST) ; \
  if ((MYRANK) == (DST)) call pkar_recvlen(1, kind(A), SRC)

#define UNPACK_RECV0(A, MYRANK, SRC, DST) \
  if ((MYRANK) == (DST)) call pkar_pop(A, 1, kind(A), SRC)

!-------------------------------------------------------------------------
! 1 dimensional array
!-------------------------------------------------------------------------
#define PACK_SEND1(A, MYRANK, SRC, DST) \
  if ((MYRANK) == (SRC)) call pkar_push(A(PTF1(A)), size(A), kind(A), DST) ; \
  if ((MYRANK) == (DST)) call pkar_recvlen(size(A), kind(A), SRC)

#define PACK_SEND1_SZ(A, MYRANK, SRC, DST, SZ) \
  if ((MYRANK) == (SRC)) call pkar_push(A(PTF1(A)), size(A), kind(A), DST) ; \
  if ((MYRANK) == (DST)) call pkar_recvlen(SZ, kind(A), SRC)

#define UNPACK_RECV1(A, MYRANK, SRC, DST) \
  if ((MYRANK) == (DST)) call pkar_pop(A(PTF1(A)), size(A), kind(A), SRC)

!-------------------------------------------------------------------------
! 2 dimensional array
!-------------------------------------------------------------------------
#define PACK_SEND2(A, MYRANK, SRC, DST) \
  if ((MYRANK) == (SRC)) call pkar_push(A(PTF2(A)), size(A), kind(A), DST) ; \
  if ((MYRANK) == (DST)) call pkar_recvlen(size(A), kind(A), SRC)

#define PACK_SEND2_SZ(A, MYRANK, SRC, DST, SZ) \
  if ((MYRANK) == (SRC)) call pkar_push(A(PTF2(A)), size(A), kind(A), DST) ; \
  if ((MYRANK) == (DST)) call pkar_recvlen(SZ, kind(A), SRC)

#define UNPACK_RECV2(A, MYRANK, SRC, DST) \
  if ((MYRANK) == (DST)) call pkar_pop(A(PTF2(A)), size(A), kind(A), SRC)

!-------------------------------------------------------------------------
! 3 dimensional array
!-------------------------------------------------------------------------
#define PACK_SEND3(A, MYRANK, SRC, DST) \
  if ((MYRANK) == (SRC)) call pkar_push(A(PTF3(A)), size(A), kind(A), DST) ; \
  if ((MYRANK) == (DST)) call pkar_recvlen(size(A), kind(A), SRC)

#define PACK_SEND3_SZ(A, MYRANK, SRC, DST, SZ) \
  if ((MYRANK) == (SRC)) call pkar_push(A(PTF3(A)), size(A), kind(A), DST) ; \
  if ((MYRANK) == (DST)) call pkar_recvlen(SZ, kind(A), SRC)

#define UNPACK_RECV3(A, MYRANK, SRC, DST) \
  if ((MYRANK) == (DST)) call pkar_pop(A(PTF3(A)), size(A), kind(A), SRC)

!-------------------------------------------------------------------------
! 4 dimensional array
!-------------------------------------------------------------------------
#define PACK_SEND4(A, MYRANK, SRC, DST) \
  if ((MYRANK) == (SRC)) call pkar_push(A(PTF4(A)), size(A), kind(A), DST) ; \
  if ((MYRANK) == (DST)) call pkar_recvlen(size(A), kind(A), SRC)

#define PACK_SEND4_SZ(A, MYRANK, SRC, DST, SZ) \
  if ((MYRANK) == (SRC)) call pkar_push(A(PTF4(A)), size(A), kind(A), DST) ; \
  if ((MYRANK) == (DST)) call pkar_recvlen(SZ, kind(A), SRC)

#define UNPACK_RECV4(A, MYRANK, SRC, DST) \
  if ((MYRANK) == (DST)) call pkar_pop(A(PTF4(A)), size(A), kind(A), SRC)


