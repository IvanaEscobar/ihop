!BOP
!    !ROUTINE: IHOP_SIZE.h
!    !INTERFACE:
! #include IHOP_SIZE.h

!    !DESCRIPTION: \bv
!     ==================================================================
!     IHOP_SIZE.h
!     ==================================================================
!     Contains IHOP source receiver array dimension
!     \ev
!EOP

!   nts    :: No. of time series points
!     ================================
!     Number of time series:
!     ================================
      INTEGER nts
#ifdef IHOP_MULTIPLE_TIMES
      PARAMETER ( nts=10 )
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
#ifdef IHOP_MULTIPLE_SOURCES
      PARAMETER ( nsd=10 )
#else 
      PARAMETER ( nsd=1 )
#endif
    
!     Number of Receivers:
!     ================================
      INTEGER nrd
      INTEGER nrr
#ifdef IHOP_MULTIPLE_RECEIVER_DEPTHS
      PARAMETER ( nrd=30 )
#else 
      PARAMETER ( nrd=1 )
#endif

#ifdef IHOP_MULTIPLE_RECEIVER_RANGES
      PARAMETER ( nrr=30 )
#else 
      PARAMETER ( nrr=1 )
#endif

!     Number of interpolation points:
!     ================================
      INTEGER IHOP_MAX_NC_SIZE
      PARAMETER ( IHOP_MAX_NC_SIZE = 10 )
!      INTEGER IHOP_NPTS_RANGE 
!      PARAMETER( IHOP_NPTS_RANGE = 6 )
!      INTEGER IHOP_IDW_NPTS
!      PARAMETER( IHOP_IDW_NPTS = 4 )



!     Cost function sizes 
!     ================================
! NFILESMAX_ihop      :: maximum number of input files
! NOBSMAX_ihop        :: maximum number of observations per file per tile
! NSAMPLESMAX        :: maximum number of samples per file per tile
! NUM_INTERP_PTS_OBS :: number of points used in interpolation for 
!                       model sampling

#ifdef ALLOW_COST
      INTEGER NFILESMAX_IHOP
      PARAMETER ( NFILESMAX_IHOP=3 )

      INTEGER NOBSMAX_IHOP
      PARAMETER ( NOBSMAX_IHOP=1000 )

      INTEGER NSAMPLESMAX
      PARAMETER ( NSAMPLESMAX=100 )

      INTEGER NSAMP_PER_OBS_MAX
      PARAMETER ( NSAMP_PER_OBS_MAX=1 )
#endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
