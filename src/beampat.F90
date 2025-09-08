#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: beamPat
MODULE beamPat
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
! Loads a source beam pattern

! !USES:
  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC NSBPPts, SrcBmPat, SBPFlag, readPat, writePat
!=======================================================================

! == Module variables ==
  INTEGER       :: NSBPPts ! Number of source beam-pattern points
  REAL (KIND=_RL90), ALLOCATABLE :: SrcBmPat( :, : )
  CHARACTER*(1) :: SBPFlag ! '*' or 'O' to indicate a directional or omni pattern

! == Derived types == None
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R ReadPat
! S/R writePat
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! SUBROUTINE: ReadPat
! !INTERFACE:
SUBROUTINE ReadPat( myThid )
! !DESCRIPTION:
! Reads the source beam pattern from a file or sets to default omni pattern.

! !USES:
  USE ihop_mod, only: SBPFile

! !INPUT PARAMETERS:
! myThid  :: my thread ID
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! I      :: Loop index
! IAllocStat :: Allocation status
! IOStat :: I/O status
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: I, IAllocStat, IOStat
!EOP

  IF (ALLOCATED(SrcBmPat)) DEALLOCATE(SrcBmPat)

  IF ( SBPFlag=='*' ) THEN
    OPEN( UNIT=SBPFile, FILE=TRIM( IHOP_fileroot ) // '.sbp', &
          STATUS='OLD', IOSTAT=IOStat, ACTION='READ' )
    IF ( IOStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      ! In adjoint mode we do not write output besides on the first run
      IF (IHOP_dumpfreq.GE.0) &
        WRITE(msgBuf,'(2A)') 'BEAMPATTERN ReadPat: ', &
          'Unable to open source beampattern file'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadPat'
    ENDIF ! IF ( IOStat.NE.0 )

    READ( SBPFile, * ) NSBPPts

    ALLOCATE( SrcBmPat( NSBPPts, 2 ), Stat=IAllocStat )
    IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BEAMPATTERN ReadPat: ', &
        'Insufficient memory for beam pattern data: reduce # SBP points'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadPat'
    ENDIF ! IF ( IAllocStat.NE.0 )

    DO I = 1, NSBPPts
      READ( SBPFile, * ) SrcBmPat( I, : )
    END DO

  ELSE   ! no pattern given, use omni source pattern
    NSBPPts = 2
    ALLOCATE( SrcBmPat( NSBPPts, 2 ), Stat = IAllocStat )
    IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(A)') 'BEAMPATTERN ReadPat: Insufficient memory'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadPat'
    ENDIF ! IF ( IAllocStat.NE.0 )

    SrcBmPat( 1, : ) = [ -180.0, 0.0 ]
    SrcBmPat( 2, : ) = [  180.0, 0.0 ]

  ENDIF ! IF ( SBPFlag=='*' )

  ! convert dB to linear scale
  SrcBmPat( :, 2 ) = 10**( SrcBmPat( :, 2 ) / 20 )

  RETURN
  END !SUBROUTINE ReadPat

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! SUBROUTINE: writePat
! !INTERFACE:
  SUBROUTINE writePat( myThid )
      ! IESCO25: UNUSED for now
! !DESCRIPTION:
! Writes the source beam pattern to an output file.

!! !USES:
!  USE ihop_mod, only: PRTFile
!
! !INPUT PARAMETERS:
! myThid  :: my thread ID
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None
!
!! !LOCAL VARIABLES:
!! msgBuf :: Informational/error message buffer
!! I      :: Loop index
!  CHARACTER*(MAX_LEN_MBUF):: msgBuf
!  INTEGER :: I
!EOP
!
!  ! I/O on main thread only
!  _BEGIN_MASTER(myThid)
!
!#ifdef IHOP_WRITE_OUT
!  IF ( SBPFlag.EQ.'*' ) THEN
!    WRITE( PRTFile, * )
!    WRITE( PRTFile, * ) '______________________________'
!    WRITE( PRTFile, * ) 'Using source beam pattern file'
!    
!    WRITE( PRTFile, * ) 'Number of source beam pattern points', NSBPPts
!
!    WRITE( PRTFile, * )
!    WRITE( PRTFile, * ) ' Angle (degrees)  Power (dB)'
!
!    DO I = 1, NSBPPts
!      WRITE( PRTFile, FMT = "( 2G11.3 )" ) SrcBmPat( I, : )
!    END DO
!  END IF
!#endif /* IHOP_WRITE_OUT */
!
!  ! I/O on main thread only
!  _END_MASTER(myThid)
!  _BARRIER
!
  RETURN
  END !SUBROUTINE writePat

END !MODULE beamPattern
