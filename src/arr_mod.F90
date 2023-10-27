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

    public WriteArrivalsASCII, WriteArrivalsBinary, MaxNArr, NArr, NArr3D, &
           Arr, Arr3D, AddArr

!=======================================================================

  INTEGER               :: MaxNArr
  INTEGER, ALLOCATABLE  :: NArr( :, : ), NArr3D( :, :, : )

  TYPE Arrival
     INTEGER :: NTopBnc, NBotBnc
     REAL    :: SrcDeclAngle, SrcAzimAngle, RcvrDeclAngle, RcvrAzimAngle, &
                A, Phase
     COMPLEX :: delay
  END TYPE

  TYPE(Arrival), ALLOCATABLE :: Arr( :, :, : ), Arr3D( :, :, :, : )

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
    
    Nt     = NArr( iz, ir )    ! # of arrivals
    NewRay = .TRUE.

    ! Is this the second bracketting ray of a pair?
    ! If so, we want to combine the arrivals to conserve space.
    ! (test this by seeing if the arrival time is close to the previous one)
    ! (also need that the phase is about the same to make sure surface and 
    ! direct paths are not joined)

    IF ( Nt >= 1 ) THEN
       IF ( afreq * ABS( delay - Arr( iz, ir, Nt )%delay ) < PhaseTol .AND. &
           ABS( Arr( iz, ir, Nt )%phase - Phase ) < PhaseTol ) NewRay = .FALSE.
    END IF

    IF ( NewRay ) THEN
       IF ( Nt >= MaxNArr ) THEN       ! space available to add an arrival?
          iArr = MINLOC( Arr( iz, ir, : )%A )   ! no: replace weakest arrival
          IF ( Amp > Arr( iz, ir, iArr( 1 ) )%A ) THEN
             Arr( iz, ir, iArr( 1 ) )%A             = SNGL( Amp )       ! amplitude
             Arr( iz, ir, iArr( 1 ) )%Phase         = SNGL( Phase )     ! phase
             Arr( iz, ir, iArr( 1 ) )%delay         = CMPLX( delay )    ! delay time
             Arr( iz, ir, iArr( 1 ) )%SrcDeclAngle  = SNGL( SrcDeclAngle )  ! angle
             Arr( iz, ir, iArr( 1 ) )%RcvrDeclAngle = SNGL( RcvrDeclAngle ) ! angle
             Arr( iz, ir, iArr( 1 ) )%NTopBnc       = NumTopBnc         ! # top bounces
             Arr( iz, ir, iArr( 1 ) )%NBotBnc       = NumBotBnc         ! # bottom bounces
          ENDIF
       ELSE
          NArr( iz, ir         )               = Nt + 1              ! # arrivals
          Arr(  iz, ir, Nt + 1 )%A             = SNGL( Amp )         ! amplitude
          Arr(  iz, ir, Nt + 1 )%Phase         = SNGL( Phase )       ! phase
          Arr(  iz, ir, Nt + 1 )%delay         = CMPLX( delay )      ! delay time
          Arr(  iz, ir, Nt + 1 )%SrcDeclAngle  = SNGL( SrcDeclAngle )    ! angle
          Arr(  iz, ir, Nt + 1 )%RcvrDeclAngle = SNGL( RcvrDeclAngle )   ! angle
          Arr(  iz, ir, Nt + 1 )%NTopBnc       = NumTopBnc           ! # top bounces
          Arr(  iz, ir, Nt + 1 )%NBotBnc       = NumBotBnc           ! # bottom bounces
       ENDIF
    ELSE      ! not a new ray
       ! calculate weightings of old ray information vs. new, based on 
       ! amplitude of the arrival
       AmpTot = Arr( iz, ir, Nt )%A + SNGL( Amp )
       w1     = Arr( iz, ir, Nt )%A / AmpTot
       w2     = REAL( Amp ) / AmpTot

       Arr( iz, ir, Nt )%delay         =  w1 * Arr( iz, ir, Nt )%delay &        
                                        + w2 * CMPLX( delay ) ! weighted sum
       Arr( iz, ir, Nt )%A             =  AmpTot
       Arr( iz, ir, Nt )%SrcDeclAngle  =  w1 * Arr( iz, ir, Nt )%SrcDeclAngle & 
                                        + w2 * SNGL( SrcDeclAngle  )
       Arr( iz, ir, Nt )%RcvrDeclAngle =  w1 * Arr( iz, ir, Nt )%RcvrDeclAngle & 
                                        + w2 * SNGL( RcvrDeclAngle )
    ENDIF

  END !SUBROUTINE AddArr

  ! **********************************************************************!

  SUBROUTINE WriteArrivalsASCII( r, Nrz, Nr, SourceType )

    ! Writes the arrival data (Amplitude, delay for each eigenray)
    ! ASCII output file

    INTEGER,           INTENT( IN ) :: Nrz, Nr      ! NRz per range, NRr
    REAL (KIND=_RL90), INTENT( IN ) :: r( Nr )      ! Rr
    CHARACTER (LEN=1), INTENT( IN ) :: SourceType   ! Beam%RunType(4:4)
    INTEGER             :: ir, iz, iArr
    REAL (KIND=_RL90)   :: factor

#ifdef IHOP_WRITE_OUT
    WRITE( ARRFile, * ) MAXVAL( NArr( 1 : Nrz, 1 : Nr ) )
#endif /* IHOP_WRITE_OUT */

    DO iz = 1, Nrz
       DO ir = 1, Nr
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
          WRITE( ARRFile, * ) NArr( iz, ir )
          DO iArr = 1, NArr( iz, ir )
             ! You can compress the output file a lot by putting in an explicit 
             ! format statement here ...
             ! However, you'll need to make sure you keep adequate precision
             WRITE( ARRFile, * ) &
             SNGL( factor ) * Arr( iz, ir, iArr )%A,             &
             SNGL( rad2deg ) * Arr( iz, ir, iArr )%Phase,         &
                        REAL( Arr( iz, ir, iArr )%delay ),       &
                       AIMAG( Arr( iz, ir, iArr )%delay ),       &
                              Arr( iz, ir, iArr )%SrcDeclAngle,  &
                              Arr( iz, ir, iArr )%RcvrDeclAngle, &
                              Arr( iz, ir, iArr )%NTopBnc,       &
                              Arr( iz, ir, iArr )%NBotBnc
          END DO  ! next arrival
#endif /* IHOP_WRITE_OUT */
       END DO  ! next receiver depth
    END DO  ! next range

  RETURN
  END !SUBROUTINE WriteArrivalsASCII

  ! **********************************************************************!

  SUBROUTINE WriteArrivalsBinary( r, Nrz, Nr, SourceType )

    ! Writes the arrival data (amplitude, delay for each eigenray)
    ! Binary output file

    INTEGER,           INTENT( IN ) :: Nrz, Nr      ! NRz per range, NRr
    REAL (KIND=_RL90), INTENT( IN ) :: r( Nr )      ! Rr
    CHARACTER (LEN=1), INTENT( IN ) :: SourceType   ! Beam%RunType(4:4)
    INTEGER                 :: ir, iz, iArr
    REAL     (KIND=_RL90)   :: factor

#ifdef IHOP_WRITE_OUT
    WRITE( ARRFile ) MAXVAL( NArr( 1 : Nrz, 1 : Nr ) )
#endif /* IHOP_WRITE_OUT */

    DO iz = 1, Nrz
       DO ir = 1, Nr
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
          WRITE( ARRFile ) NArr( iz, ir )
          DO iArr = 1, NArr( iz, ir )
             ! integers written out as reals below for fast reading in Matlab
             WRITE( ARRFile ) &
            SNGL( factor * Arr( iz, ir, iArr )%A ),           &
            SNGL( rad2deg * Arr( iz, ir, iArr )%Phase ),       &
                           Arr( iz, ir, iArr )%delay,         &
                           Arr( iz, ir, iArr )%SrcDeclAngle,  &
                           Arr( iz, ir, iArr )%RcvrDeclAngle, &
                     REAL( Arr( iz, ir, iArr )%NTopBnc ),     &
                     REAL( Arr( iz, ir, iArr )%NBotBnc )
          END DO   ! next arrival
#endif /* IHOP_WRITE_OUT */
       END DO   ! next receiver depth
    END DO   ! next range

  RETURN
  END !SUBROUTINE WriteArrivalsBinary

  ! **********************************************************************!
END !MODULE arr_mod
