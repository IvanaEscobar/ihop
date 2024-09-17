#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE arr_mod
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

  USE ihop_mod,   only: rad2deg, ARRFile

! ! USES
  implicit none
!  == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

  PRIVATE

! public interfaces
!=======================================================================

    public WriteArrivalsASCII, WriteArrivalsBinary, &
           MaxNArr, NArr, Arr, AddArr, U
#ifdef IHOP_THREED
    public NArr3D, Arr3D
#endif /* IHOP_THREED */

!=======================================================================

  INTEGER               :: MaxNArr
  INTEGER, ALLOCATABLE  :: NArr( :, : )
  COMPLEX, ALLOCATABLE  :: U( :, : )
#ifdef IHOP_THREED
  INTEGER, ALLOCATABLE  :: NArr3D( :, :, : )
#endif /* IHOP_THREED */

  TYPE Arrival
     INTEGER :: NTopBnc, NBotBnc
     REAL    :: SrcDeclAngle, RcvrDeclAngle, A, Phase
#ifdef IHOP_THREED
     REAL    :: SrcAzimAngle, RcvrAzimAngle
#endif /* IHOP_THREED */
     COMPLEX :: delay
  END TYPE

  TYPE(Arrival), ALLOCATABLE :: Arr( :, :, : )
#ifdef IHOP_THREED
  TYPE(Arrival), ALLOCATABLE :: Arr3D( :, :, :, : )
#endif /* IHOP_THREED */

CONTAINS
  SUBROUTINE AddArr( afreq, iz, ir, Amp, Phase, delay, SrcDeclAngle, &
                     RcvrDeclAngle, NumTopBnc, NumBotBnc )

    ! ADDs the amplitude and delay for an ARRival into a matrix of same.
    ! Extra logic included to keep only the strongest arrivals.

    ! arrivals with essentially the same phase are grouped into one
    REAL,   PARAMETER                  :: PhaseTol = 0.05
    INTEGER,              INTENT( IN ) :: NumTopBnc, NumBotBnc, iz, ir
    REAL    (KIND=_RL90), INTENT( IN ) :: afreq, Amp, Phase, SrcDeclAngle, &
                                          RcvrDeclAngle
    COMPLEX (KIND=_RL90), INTENT( IN ) :: delay
    LOGICAL              :: NewRay
    INTEGER              :: iArr( 1 ), Nt
    REAL                 :: AmpTot, w1, w2

    Nt     = NArr( ir, iz )    ! # of arrivals
    NewRay = .TRUE.

    ! Is this the second bracketting ray of a pair?
    ! If so, we want to combine the arrivals to conserve space.
    ! (test this by seeing if the arrival time is close to the previous one)
    ! (also need that the phase is about the same to make sure surface and
    ! direct paths are not joined)

    IF ( Nt >= 1 ) THEN
       IF ( afreq * ABS( delay - Arr( Nt, ir, iz )%delay ) < PhaseTol .AND. &
           ABS( Arr( Nt, ir, iz )%phase - Phase ) < PhaseTol ) NewRay = .FALSE.
    END IF

    IF ( NewRay ) THEN
       IF ( Nt >= MaxNArr ) THEN       ! space available to add an arrival?
          iArr = MINLOC( Arr( :, ir, iz )%A )   ! no: replace weakest arrival
          IF ( Amp > Arr( iArr(1), ir, iz )%A ) THEN
             Arr( iArr(1), ir, iz )%A             = SNGL( Amp )       ! amplitude
             Arr( iArr(1), ir, iz )%Phase         = SNGL( Phase )     ! phase
             Arr( iArr(1), ir, iz )%delay         = CMPLX( delay )    ! delay time
             Arr( iArr(1), ir, iz )%SrcDeclAngle  = SNGL( SrcDeclAngle )  ! angle
             Arr( iArr(1), ir, iz )%RcvrDeclAngle = SNGL( RcvrDeclAngle ) ! angle
             Arr( iArr(1), ir, iz )%NTopBnc       = NumTopBnc         ! # top bounces
             Arr( iArr(1), ir, iz )%NBotBnc       = NumBotBnc         ! # bottom bounces
          ENDIF
       ELSE
          NArr( ir, iz       )               = Nt + 1              ! # arrivals
          Arr( Nt + 1, ir, iz)%A             = SNGL( Amp )         ! amplitude
          Arr( Nt + 1, ir, iz)%Phase         = SNGL( Phase )       ! phase
          Arr( Nt + 1, ir, iz)%delay         = CMPLX( delay )      ! delay time
          Arr( Nt + 1, ir, iz)%SrcDeclAngle  = SNGL( SrcDeclAngle )    ! angle
          Arr( Nt + 1, ir, iz)%RcvrDeclAngle = SNGL( RcvrDeclAngle )   ! angle
          Arr( Nt + 1, ir, iz)%NTopBnc       = NumTopBnc           ! # top bounces
          Arr( Nt + 1, ir, iz)%NBotBnc       = NumBotBnc           ! # bottom bounces
       ENDIF
    ELSE      ! not a new ray
       ! calculate weightings of old ray information vs. new, based on
       ! amplitude of the arrival
       AmpTot = Arr( Nt, ir, iz )%A + SNGL( Amp )
       w1     = Arr( Nt, ir, iz )%A / AmpTot
       w2     = REAL( Amp ) / AmpTot

       Arr( Nt, ir, iz)%delay         =  w1 * Arr( Nt, ir, iz )%delay &
                                       + w2 * CMPLX( delay ) ! weighted sum
       Arr( Nt, ir, iz)%A             =  AmpTot
       Arr( Nt, ir, iz)%SrcDeclAngle  =  w1 * Arr( Nt, ir, iz )%SrcDeclAngle &
                                       + w2 * SNGL( SrcDeclAngle  )
       Arr( Nt, ir, iz)%RcvrDeclAngle =  w1 * Arr( Nt, ir, iz )%RcvrDeclAngle &
                                        + w2 * SNGL( RcvrDeclAngle )
    ENDIF

  END !SUBROUTINE AddArr

  ! **********************************************************************!

  SUBROUTINE WriteArrivalsASCII( r, Nrz, Nrr, SourceType )

    ! Writes the arrival data (Amplitude, delay for each eigenray)
    ! ASCII output file

    INTEGER,           INTENT( IN ) :: Nrz, Nrr     ! NRz per range, NRr
    REAL (KIND=_RL90), INTENT( IN ) :: r( Nr )      ! Rr
    CHARACTER (LEN=1), INTENT( IN ) :: SourceType   ! Beam%RunType(4:4)
    INTEGER             :: ir, iz, iArr
    REAL (KIND=_RL90)   :: factor

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) RETURN

#ifdef IHOP_WRITE_OUT
    WRITE( ARRFile, * ) MAXVAL( NArr( 1:Nrr, 1:Nrz ) )
#endif /* IHOP_WRITE_OUT */

    DO iz = 1, Nrz
       DO ir = 1, Nrr
          IF ( SourceType == 'X' ) THEN   ! line source
             factor =  4.0 * SQRT( PI )
          ELSE                            ! point source: default
             IF ( r ( ir ) == 0 ) THEN
                factor = 1e5                   ! avoid /0 at origin
             ELSE
                factor = 1. / SQRT( r( ir ) )  ! cyl. spreading
             END IF
          END IF

#ifdef IHOP_WRITE_OUT
          WRITE( ARRFile, * ) NArr( ir, iz )
          DO iArr = 1, NArr( ir, iz )
             ! You can compress the output file a lot by putting in an explicit
             ! format statement here ...
             ! However, you'll need to make sure you keep adequate precision
             WRITE( ARRFile, * ) &
             SNGL( factor ) * Arr( iArr, ir, iz )%A,             &
             SNGL( rad2deg ) *Arr( iArr, ir, iz )%Phase,        &
                        REAL( Arr( iArr, ir, iz )%delay ),       &
                       AIMAG( Arr( iArr, ir, iz )%delay ),       &
                              Arr( iArr, ir, iz )%SrcDeclAngle,  &
                              Arr( iArr, ir, iz )%RcvrDeclAngle, &
                              Arr( iArr, ir, iz )%NTopBnc,       &
                              Arr( iArr, ir, iz )%NBotBnc
          END DO  ! next arrival
#endif /* IHOP_WRITE_OUT */
       END DO  ! next receiver depth
    END DO  ! next range

  RETURN
  END !SUBROUTINE WriteArrivalsASCII

  ! **********************************************************************!

  SUBROUTINE WriteArrivalsBinary( r, Nrz, Nrr, SourceType )

    ! Writes the arrival data (amplitude, delay for each eigenray)
    ! Binary output file

    INTEGER,           INTENT( IN ) :: Nrz, Nrr     ! NRz per range, NRrr
    REAL (KIND=_RL90), INTENT( IN ) :: r( Nr )      ! Rr
    CHARACTER (LEN=1), INTENT( IN ) :: SourceType   ! Beam%RunType(4:4)
    INTEGER                 :: ir, iz, iArr
    REAL     (KIND=_RL90)   :: factor

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) RETURN

#ifdef IHOP_WRITE_OUT
    WRITE( ARRFile ) MAXVAL( NArr( 1:Nrr, 1:Nrz ) )
#endif /* IHOP_WRITE_OUT */

    DO iz = 1, Nrz
       DO ir = 1, Nrr
          IF ( SourceType == 'X' ) THEN   ! line source
             factor = 4.0 * SQRT( PI )
          ELSE                            ! point source
             IF ( r ( ir ) == 0 ) THEN
                factor = 1e5                   ! avoid /0 at origin
             ELSE
                factor = 1. / SQRT( r( ir ) )  ! cyl. spreading
             END IF
          END IF

#ifdef IHOP_WRITE_OUT
          WRITE( ARRFile ) NArr( ir, iz )
          DO iArr = 1, NArr( ir, iz )
             ! integers written out as reals below for fast reading in Matlab
             WRITE( ARRFile ) &
            SNGL( factor * Arr( iArr, ir, iz )%A ),           &
            SNGL( rad2deg *Arr( iArr, ir, iz )%Phase ),       &
                           Arr( iArr, ir, iz )%delay,         &
                           Arr( iArr, ir, iz )%SrcDeclAngle,  &
                           Arr( iArr, ir, iz )%RcvrDeclAngle, &
                     REAL( Arr( iArr, ir, iz )%NTopBnc ),     &
                     REAL( Arr( iArr, ir, iz )%NBotBnc )
          END DO   ! next arrival
#endif /* IHOP_WRITE_OUT */
       END DO   ! next receiver depth
    END DO   ! next range

  RETURN
  END !SUBROUTINE WriteArrivalsBinary

  ! **********************************************************************!
END !MODULE arr_mod
