!BOP
!    !ROUTINE: BELLI_SIZE.h
!    !INTERFACE:
! #include BELLI_SIZE.h

!    !DESCRIPTION: \bv
!     ==================================================================
!     BELLI_SIZE.h
!     ==================================================================
!     Contains BELLI source receiver array dimension
!     \ev
!EOP

!   nts    :: No. of time series points
!     ================================
!     Number of time series:
!     ================================
      INTEGER nts
#ifdef BELLI_MULTIPLE_TIMES
      PARAMETER ( nts=1080 )
#else
      PARAMETER ( nts=1 )
#endif


!   nsd    :: No. of sound sources at range of 0 m
!   nrd    :: No. of sound receivers at a single range
!   nrr    :: No. of sound receivers at a single depth
!     ================================
!     Number of Sources:
!     ================================
      INTEGER nsd
#ifdef BELLI_MULTIPLE_SOURCES
      PARAMETER ( nsd=10 )
#else
      PARAMETER ( nsd=1 )
#endif

!     Number of Receivers:
!     ================================
      INTEGER nrd
      INTEGER nrr
#ifdef BELLI_MULTIPLE_RECEIVER_DEPTHS
      PARAMETER ( nrd=30 )
#else
      PARAMETER ( nrd=1 )
#endif

#ifdef BELLI_MULTIPLE_RECEIVER_RANGES
      PARAMETER ( nrr=30 )
#else
      PARAMETER ( nrr=1 )
#endif

!     Number of interpolation points:
!     ================================
      INTEGER BELLI_MAX_NC_SIZE
      PARAMETER ( BELLI_MAX_NC_SIZE = 15 )
      INTEGER BELLI_MAX_RANGE
      PARAMETER( BELLI_MAX_RANGE = 22 )
      INTEGER BELLI_MAX_IDW
      PARAMETER( BELLI_MAX_IDW = 4 )



!     Cost function sizes
!     ================================
! NFILESMAX_BELLI      :: maximum number of input files
! NOBSMAX_BELLI        :: maximum number of observations per file per tile

#ifdef ALLOW_COST
      INTEGER NFILESMAX_BELLI
      PARAMETER ( NFILESMAX_BELLI=1 )

      INTEGER NOBSMAX_BELLI
      PARAMETER ( NOBSMAX_BELLI=10 )

#endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
