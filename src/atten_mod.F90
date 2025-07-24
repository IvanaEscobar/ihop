#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: atten_mod
MODULE atten_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!  Routines to convert a sound speed and attenuation units
!  Includes a formula for volume attenuation

! !USES:
#ifdef IHOP_WRITE_OUT
  USE ihop_mod, only: PRTFile
#endif
  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

! !SCOPE: 
  PRIVATE
!=======================================================================
   PUBLIC CRCI, T, Salinity, pH, z_bar, iBio, NBioLayers, bio
!=======================================================================

! == Module variables ==
  INTEGER, PARAMETER :: MaxBioLayers = 200
  INTEGER            :: iBio, NBioLayers
  ! Francois-Garrison volume attenuation; temperature, salinity, ...
  REAL (KIND=_RL90)  :: T = 20, Salinity = 35, pH = 8, z_bar = 0, FG

! == Derived types ==
  TYPE bioStructure
    REAL (KIND=_RL90) :: Z1, Z2, f0, Q, a0
  END TYPE bioStructure

  TYPE( bioStructure ) :: bio( MaxBioLayers )
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! FXN CRCI
! FXN Franc_Garr
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: CRCI
! !INTERFACE:
  FUNCTION CRCI( z, c, alpha, AttenUnit, beta, fT, myThid )
! !DESCRIPTION:
! Converts real wave speed and attenuation to a single
! complex wave speed (with positive imaginary part)

! !INPUT PARAMETERS:
! z        :: depth in meters
! c        :: real part of sound speed in m/s
! alpha    :: imaginary part of sound speed in Nepers/m
! AttenUnit:: attenuation unit (2 characters)
! beta     :: power law exponent for frequency dependence
! fT       :: frequency for which the attenuation is specified
! myThid   :: my thread ID
  REAL (KIND=_RL90), INTENT( IN ) :: z, c, alpha, beta, fT
  CHARACTER*(2),     INTENT( IN ) :: AttenUnit
  INTEGER,           INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS:
! CRCI     :: complex sound speed
  COMPLEX (KIND=_RL90) :: CRCI

! !LOCAL VARIABLES:
! f2 :: frequency squared
! afreq :: angular frequency
! alphaT :: imaginary part of sound speed in Nepers/m
! Thorp :: Thorp attenuation in Nepers/m
! a :: biological attenuation in Nepers/m
! FG :: Francois-Garrison attenuation in Nepers/m
   REAL (KIND=_RL90) :: f2, afreq, alphaT, Thorp, a, FG
!EOP

! Options:
! AttenUnit
! 6 Cases:    N Nepers/meter
!             M dB/meter      (M for Meters)
!             m dB/meter with a power law
!             F dB/m-kHz      (F for frequency dependent)
!             W dB/wavelength (W for Wavelength)
!             Q Q
!             L Loss parameter
!
! second letter adds volume attenuation according to standard laws:
!             T for Thorp
!             F for Francois Garrison
!             B for biological
!
! freq is the current frequency
! freq0 is the reference frequency for which the dB/meter was specified
!  (used only for 'm')

  afreq = 2.0 * PI * IHOP_freq

  !  Convert to Nepers/m
  alphaT = 0.0
  SELECT CASE ( AttenUnit( 1:1 ) )
  CASE ( 'N' )
    alphaT = alpha
  CASE ( 'M' )   ! dB/m
    alphaT = alpha / 8.6858896D0
  CASE ( 'm' )   ! dB/m with power law
    alphaT = alpha / 8.6858896D0
    IF ( IHOP_freq.LT.fT ) THEN   ! frequency raised to the power beta
      alphaT = alphaT * ( IHOP_freq / IHOP_freq ) ** beta
    ELSE                    ! linear in frequency
      alphaT = alphaT * ( IHOP_freq / IHOP_freq ) * &
                ( fT / IHOP_freq ) ** ( beta - 1 )
   ENDIF

  CASE ( 'F' )   ! dB/(m kHz)
    alphaT = alpha * IHOP_freq / 8685.8896D0
  CASE ( 'W' )   ! dB/wavelength
    IF ( c.NE.0.0 ) alphaT = alpha * IHOP_freq / ( 8.6858896D0 * c )
   !        The following lines give f^1.25 frequency dependence
   !        FAC = SQRT( SQRT( IHOP_freq / 50.0 ) )
   !        IF ( c /= 0.0 ) alphaT=FAC*alpha*IHOP_freq / ( 8.6858896D0*c )
  CASE ( 'Q' )   ! Quality factor
    IF ( c * alpha.NE.0.0 ) alphaT = afreq / ( 2.0 * c * alpha )
  CASE ( 'L' )   ! loss parameter
    IF ( c.NE.0.0 ) alphaT = alpha * afreq / c
  CASE DEFAULT
    ! Do nothing
    alphaT = 0.0
    STOP '1 ABNORMAL END: S/R CRCI'
  END SELECT ! AttenUnit( 1:1 )

  ! volume attenuation
  SELECT CASE ( AttenUnit( 2:2 ) )
  CASE ( 'T' )
    f2 = ( IHOP_freq / 1000.0 ) ** 2

    ! Original formula from Thorp 1967
    ! Thorp = 40.0 * f2 / ( 4100.0+f2 ) + 0.1 * f2 / ( 1.0+f2 )! dB/kyard
    ! Thorp = Thorp / 914.4D0              ! dB / m
    ! Thorp = Thorp / 8.6858896D0          ! Nepers / m
    ! Updated formula from JKPS Eq. 1.34
    Thorp = 3.3d-3 + 0.11 * f2 / ( 1.0+f2 ) + 44.0 * f2 / ( 4100.0+f2 ) &
            + 3d-4 * f2       ! dB/km
    Thorp = Thorp / 8685.8896 ! Nepers / m
    alphaT = alphaT + Thorp

  CASE ( 'F' )   ! Francois-Garrison
    FG     = Franc_Garr( IHOP_freq / 1000 )   ! dB/km
    FG     = FG / 8685.8896              ! Nepers / m
    alphaT = alphaT + FG

  CASE ( 'B' )   ! biological attenuation per Orest Diachok
    DO iBio = 1, NBioLayers
      IF ( z.GE.bio( iBio )%Z1 .AND. z.LE.bio( iBio )%Z2 ) THEN
        a = bio( iBio )%a0 / ( ( 1.0-bio( iBio )%f0**2 &
               / IHOP_freq**2  )**2 &
               + 1.0 / bio( iBio )%Q**2 ) ! dB/km
        a = a / 8685.8896              ! Nepers / m
        alphaT = alphaT + a

      ENDIF
    ENDDO ! DO iBio

  CASE DEFAULT
    ! Do nothing
    alphaT = alphaT
  END SELECT ! AttenUnit( 2:2 )

  ! Convert Nepers/m to equivalent imaginary sound speed
  alphaT = alphaT * c * c / afreq
  CRCI   = CMPLX( c, alphaT, KIND=_RL90 )

  IF ( alphaT.GT.c ) THEN
#ifdef IHOP_WRITE_OUT
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) THEN
      WRITE( PRTFile, * ) 'Complex sound speed: ', CRCI
      WRITE( PRTFile, '(2A)' ) 'Usually this means you have an attenuation',&
        'that is way too high'
    ENDIF

    WRITE(errorMessageUnit,'(2A)') 'ATTENMOD CRCI: complex sound speed',&
      'has an imaginary part > real part'
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R CRCI'
  ENDIF ! IF ( alphaT.GT.c )

  RETURN
  END !FUNCTION CRCI

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: Franc_Garr
! !INTERFACE:
  FUNCTION Franc_Garr( f )
! !DESCRIPTION:
! Computes the Francois-Garrison attenuation
! Based on a Matlab version by D. Jackson APL-UW
! Verified using F-G Table IV

! !USES: None

! !INPUT PARAMETERS:
! f :: frequency in kHz
  REAL (KIND=_RL90), INTENT(IN) :: f
! !OUTPUT PARAMETERS:
! Franc_Garr :: volume attenuation in dB/km
  REAL (KIND=_RL90) :: Franc_Garr

! !LOCAL VARIABLES:
! c          :: sound speed in m/s
! A1, A2, A3 :: coefficients for boric acid, magnesium sulfate, viscosity
! P1, P2, P3 :: pressure correction factors
! f1, f2     :: frequency correction factors for boric acid and magnesium sulfate
  REAL (KIND=_RL90) :: c, A1, A2, A3, P1, P2, P3, f1, f2
!EOP

  c = 1412 + 3.21 * T + 1.19 * Salinity + 0.0167 * z_bar

  ! Boric acid contribution
  A1 = 8.86 / c * 10**( 0.78 * pH - 5 )
  P1 = 1
  f1 = 2.8 * sqrt( Salinity / 35 ) * 10**( 4 - 1245 / ( T + 273 ) )

  ! Magnesium sulfate contribution
  A2 = 21.44 * Salinity / c * ( 1 + 0.025 * T )
  P2 = 1 - 1.37D-4 * z_bar + 6.2D-9 * z_bar**2
  f2 = 8.17 * 10**( 8 - 1990 / ( T + 273 ) ) / ( 1 + 0.0018*( Salinity - 35 ) )

  ! Viscosity
  P3 = 1 - 3.83D-5 * z_bar + 4.9D-10 * z_bar**2
  IF ( T.LT.20 ) THEN
    A3 = 4.937D-4 - 2.59D-5 * T + 9.11D-7 * T**2 - 1.5D-8  * T**3
  ELSE
    A3 = 3.964D-4 -1.146D-5 * T + 1.45D-7 * T**2 - 6.5D-10 * T**3
  ENDIF

  Franc_Garr = A1 * P1 * ( f1 * f**2 ) / ( f1**2 + f**2 ) &
              + A2 * P2 * ( f2 * f**2 ) / ( f2**2 + f**2 ) &
              + A3 * P3 * f**2

  RETURN
  END !FUNCTION Franc_Garr

END !MODULE atten_mod
