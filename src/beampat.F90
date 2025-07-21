#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE beamPat
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

  ! Loads a source beam pattern

  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

  PRIVATE

! public interfaces
!=======================================================================
    public NSBPPts, SrcBmPat, SBPFlag, readPat, writePat
!=======================================================================

  INTEGER       :: NSBPPts ! Number of source beam-pattern points
  REAL (KIND=_RL90), ALLOCATABLE :: SrcBmPat( :, : )
  CHARACTER*(1) :: SBPFlag ! '*' or 'O' to indicate a directional or omni pattern

CONTAINS
! **************************************************************************** !
  SUBROUTINE ReadPat( myThid )
  USE ihop_mod, only: SBPFile

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    INTEGER                          :: I, IAllocStat, IOStat

    IF (ALLOCATED(SrcBmPat)) DEALLOCATE(SrcBmPat)

    IF ( SBPFlag=='*' ) THEN
      OPEN( UNIT=SBPFile, FILE=TRIM( IHOP_fileroot ) // '.sbp', &
            STATUS='OLD', IOSTAT=IOStat, ACTION='READ' )
      IF ( IOStat/=0 ) THEN
#ifdef IHOP_WRITE_OUT
        ! In adjoint mode we do not write output besides on the first run
        IF (IHOP_dumpfreq.GE.0) &
          WRITE(msgBuf,'(2A)') 'BEAMPATTERN ReadPat: ', &
                               'Unable to open source beampattern file'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadPat'
      END IF

      READ( SBPFile, * ) NSBPPts

      ALLOCATE( SrcBmPat( NSBPPts, 2 ), Stat=IAllocStat )
      IF ( IAllocStat/=0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'BEAMPATTERN ReadPat: ', &
        'Insufficient memory for beam pattern data: reduce # SBP points'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadPat'
      END IF

      DO I = 1, NSBPPts
        READ( SBPFile, * ) SrcBmPat( I, : )
      END DO

    ELSE   ! no pattern given, use omni source pattern
      NSBPPts = 2
      ALLOCATE( SrcBmPat( NSBPPts, 2 ), Stat = IAllocStat )
      IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(A)') 'BEAMPATTERN ReadPat: Insufficient memory'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadPat'
      END IF
      SrcBmPat( 1, : ) = [ -180.0, 0.0 ]
      SrcBmPat( 2, : ) = [  180.0, 0.0 ]

    ENDIF

    ! convert dB to linear scale
    SrcBmPat( :, 2 ) = 10**( SrcBmPat( :, 2 ) / 20 )

  RETURN
  END !SUBROUTINE ReadPat

! **************************************************************************** !
  SUBROUTINE writePat( myThid )
  USE ihop_mod, only: PRTFile

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
  INTEGER :: I

  ! I/O on main thread only
  _BEGIN_MASTER(myThid)

#ifdef IHOP_WRITE_OUT
  ! In adjoint mode we do not write output besides on the first run
  IF (IHOP_dumpfreq.GE.0) THEN

    IF ( SBPFlag == '*' ) THEN
      WRITE( PRTFile, * )
      WRITE( PRTFile, * ) '______________________________'
      WRITE( PRTFile, * ) 'Using source beam pattern file'
      
      WRITE( PRTFile, * ) 'Number of source beam pattern points', NSBPPts

      WRITE( PRTFile, * )
      WRITE( PRTFile, * ) ' Angle (degrees)  Power (dB)'

      DO I = 1, NSBPPts
        WRITE( PRTFile, FMT = "( 2G11.3 )" ) SrcBmPat( I, : )
      END DO
    END IF

  END IF ! don't write in adjoint mode
#endif /* IHOP_WRITE_OUT */
  ! I/O on main thread only
  _END_MASTER(myThid)
  _BARRIER



  RETURN
  END !SUBROUTINE writePat

END !MODULE beamPattern
