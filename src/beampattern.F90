#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE beamPattern
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

  ! Loads a source beam pattern

  USE ihop_mod, only: PRTFile, SBPFile

  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

  PRIVATE

! public interfaces
!=======================================================================
    public NSBPPts, SrcBmPat, SBPFlag, ReadPat
!=======================================================================

  INTEGER                    :: NSBPPts ! Number of source beam-pattern points
  REAL (KIND=_RL90), ALLOCATABLE :: SrcBmPat( :, : )
  CHARACTER (LEN=1)          :: SBPFlag ! '*' or 'O' to indicate a directional or omni pattern

CONTAINS
  SUBROUTINE ReadPat( myThid )

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    INTEGER                          :: I, IAllocStat, IOStat

    IF ( SBPFlag == '*' ) THEN
#ifdef IHOP_WRITE_OUT
       ! In adjoint mode we do not write output besides on the first run
       IF (IHOP_dumpfreq.GE.0) THEN
        WRITE( PRTFile, * )
        WRITE( PRTFile, * ) '______________________________'
        WRITE( PRTFile, * ) 'Using source beam pattern file'
       ENDIF
#endif /* IHOP_WRITE_OUT */

       OPEN( UNIT = SBPFile, FILE = TRIM( IHOP_fileroot ) // '.sbp', STATUS = 'OLD',&
            IOSTAT = IOStat, ACTION = 'READ' )
       IF ( IOstat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
            ! In adjoint mode we do not write output besides on the first run
            IF (IHOP_dumpfreq.GE.0) &
             WRITE( PRTFile, * ) 'SBPFile = ', TRIM( IHOP_fileroot ) // '.sbp'
            WRITE(msgBuf,'(2A)') 'BEAMPATTERN ReadPat: ', &
                                 'Unable to open source beampattern file'
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadPat'
       END IF

       READ(  SBPFile, * ) NSBPPts
#ifdef IHOP_WRITE_OUT
       ! In adjoint mode we do not write output besides on the first run
       IF (IHOP_dumpfreq.GE.0) &
        WRITE( PRTFile, * ) 'Number of source beam pattern points', NSBPPts
#endif /* IHOP_WRITE_OUT */

       IF (ALLOCATED(SrcBmPat)) DEALLOCATE(SrcBmPat)
       ALLOCATE( SrcBmPat( NSBPPts, 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'BEAMPATTERN ReadPat: ', &
            'Insufficient memory for beam pattern data: reduce # SBP points'
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadPat'
        END IF

#ifdef IHOP_WRITE_OUT
       ! In adjoint mode we do not write output besides on the first run
       IF (IHOP_dumpfreq.GE.0) &
       WRITE( PRTFile, * )
       ! In adjoint mode we do not write output besides on the first run
       IF (IHOP_dumpfreq.GE.0) &
       WRITE( PRTFile, * ) ' Angle (degrees)  Power (dB)'
#endif /* IHOP_WRITE_OUT */

       DO I = 1, NSBPPts
          READ(  SBPFile, * ) SrcBmPat( I, : )
#ifdef IHOP_WRITE_OUT
          ! In adjoint mode we do not write output besides on the first run
          IF (IHOP_dumpfreq.GE.0) &
            WRITE( PRTFile, FMT = "( 2G11.3 )" ) SrcBmPat( I, : )
#endif /* IHOP_WRITE_OUT */
       END DO

    ELSE   ! no pattern given, use omni source pattern
        NSBPPts = 2
        IF (ALLOCATED(SrcBmPat)) DEALLOCATE(SrcBmPat)
        ALLOCATE( SrcBmPat( 2, 2 ), Stat = IAllocStat )
        IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'BEAMPATTERN ReadPat: ', &
                                 'Insufficient memory'
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadPat'
        END IF
        SrcBmPat( 1, : ) = [ -180.0, 0.0 ]
        SrcBmPat( 2, : ) = [  180.0, 0.0 ]
    ENDIF

    SrcBmPat( :, 2 ) = 10**( SrcBmPat( :, 2 ) / 20 )  ! convert dB to linear scale

  RETURN
  END !SUBROUTINE ReadPat

END !MODULE beamPattern
