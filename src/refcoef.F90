#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: refCoef
MODULE refCoef
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   Module for reading and interpolating reflection coefficients.
! **************************************************************************** !
! NOTE: To be able to read using direct access, we assume a fixed pad of
! characters per line. So each line in BRCFile/TRCFile will have a fixed length
! of 100 characters. Helps TAF create the adjoint model. IESCO24
! **************************************************************************** !

! !USES:
  USE ihop_mod,   only: PRTFile, BRCFile, TRCFile, IRCFile
  IMPLICIT
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC ReadReflectionCoefficient, InterpolateReflectionCoefficient, &
        ReflectionCoef, RBot, RTop, NBotPts, NTopPts, writeRefCoef
!=======================================================================

! == Module variables ==
  INTEGER               :: NBotPts, NTopPts
  INTEGER               :: NkTab
  INTEGER,              ALLOCATABLE :: iTab( : )
  REAL    (KIND=_RL90), ALLOCATABLE :: xTab( : )
  COMPLEX (KIND=_RL90), ALLOCATABLE :: fTab( : ), gTab( : )
  
! == Derived types ==
  TYPE ReflectionCoef
    REAL(KIND=_RL90) :: theta, R, phi
  END TYPE

  TYPE(ReflectionCoef), ALLOCATABLE :: RBot( : ), RTop( : )
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R ReadReflectionCoefficient
! S/R InterpolateReflectionCoefficient
! S/R InterpolateIRC
! S/R writeRefCoef
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ReadReflectionCoefficient
! !INTERFACE:
  SUBROUTINE ReadReflectionCoefficient( myThid )
! !DESCRIPTION:
!   Reads the reflection coefficient from a file or initializes it to zero.

! !USES:
  USE bdry_mod, only: Bdry

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! BotRC, TopRC :: Flags indicating if reflection coefficients are read from files
! line :: Buffer for reading lines from the reflection coefficient files
! itheta, ik :: Loop indices for reading angles and coefficients
! IOStat :: I/O status for file operations
! iAllocStat :: Memory allocation status
! rec :: Record number for reading from files
! freq :: Frequency for internal reflection coefficient file
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  CHARACTER (LEN=1)   :: BotRC, TopRC
  CHARACTER (LEN=100) :: line
  INTEGER             :: itheta, ik, IOStat, iAllocStat, rec
  REAL  (KIND=_RL90)  :: freq
!EOP

  ! initiate local parameters
  BotRC = Bdry%Bot%HS%Opt( 1:1 )
  TopRC = Bdry%Top%HS%Opt( 2:2 )

  IF ( ALLOCATED( RTop ) ) DEALLOCATE( RTop )
  IF ( ALLOCATED( RBot ) ) DEALLOCATE( RBot )
  IF ( ALLOCATED( xTab ) ) DEALLOCATE( xTab, fTab, gTab, iTab )

  ! == read any Reflection Coefficients ==
  IF ( BotRC.EQ.'F' ) THEN
    OPEN( FIlE=TRIM( IHOP_fileroot )//'.brc', UNIT=BRCFile, STATUS='OLD',&
          IOSTAT=IOStat, ACTION='READ', ACCESS='direct', RECL=100 )
    IF ( IOStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)')'REFCOEF ReadReflectionCoeffcient: ',&
        'Unable to open Bottom Reflection Coefficient file'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
    ENDIF

    READ( BRCFile, REC=1, IOSTAT=IOStat, FMT='(A)' ) line
    READ( line,'(I10)' ) NBotPts

    ALLOCATE( RBot( NBotPts ), Stat=IAllocStat )
    IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'REFCOEF ReadReflectionCoeffcient: ', &
        'Insufficient memory for bot. refl. coef.: reduce # points'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
    ENDIF

    DO itheta=1,NBotPts
      rec = itheta+1
      READ( BRCFile, REC=rec, IOSTAT=IOStat, FMT='(A)' ) line
      READ( line,'(3F5.3)' ) RBot(itheta)%theta, RBot(itheta)%R, &
                              RBot(itheta)%phi
    ENDDO

    CLOSE( BRCFile )
    RBot%phi = deg2rad * RBot%phi   ! convert to radians

  ELSE   ! should allocate something anyway, since variable is passed
    ALLOCATE(  RBot( 1 ), Stat=IAllocStat )
    NBotPts = 1         !RG
    RBot(1)%theta = 0.  !RG
    RBot(1)%R = 0.      !RG
    RBot(1)%phi = 0.    !RG

  ENDIF ! BotRC .EQ. ...

  ! Optionally read in top reflection coefficient
  IF ( TopRC.EQ.'F' ) THEN
    OPEN( FIlE=TRIM( IHOP_fileroot )//'.trc', UNIT=TRCFile, STATUS='OLD',&
          IOSTAT=IOStat, ACTION='READ', ACCESS='direct', RECL=100 )
    IF ( IOStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'REFCOEF ReadReflectionCoeffcient: ', &
        'Unable to open Top Reflection Coefficient file'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
    ENDIF

    READ( TRCFile, REC=1, IOSTAT=IOStat, FMT='(A)' ) line
    READ( line,'(I10)' ) NTopPts

    ALLOCATE( RTop( NTopPts ), Stat=IAllocStat )
    IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'REFCOEF ReadReflectionCoeffcient: ', &
        'Insufficient memory for top refl. coef.: reduce # points'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
    ENDIF

    DO itheta=1,NTopPts
      rec = itheta+1
      READ( TRCFile, REC=rec, IOSTAT=IOStat, FMT='(A)' ) line
      READ( line,'(3F5.3)' ) RTop(itheta)%theta, RTop(itheta)%R, &
                              RTop(itheta)%phi
    ENDDO

    CLOSE( TRCFile )
    RTop%phi = deg2rad *  RTop%phi   ! convert to radians

  ELSE   ! should allocate something anyway, since variable is passed
    ALLOCATE( RTop( 1 ), Stat=IAllocStat )
    NTopPts = 1         !RG
    RTop(1)%theta = 0.  !RG
    RTop(1)%R = 0.      !RG
    RTop(1)%phi = 0.    !RG

  ENDIF ! TopRC .EQ. ...


  ! Optionally read in internal reflection coefficient data
  IF ( BotRC.EQ.'P' ) THEN
    OPEN( FIlE=TRIM( IHOP_fileroot )//'.irc', UNIT=TRCFile, STATUS='OLD',&
          IOSTAT=IOStat, ACTION='READ', ACCESS='direct', RECL=100 )
    IF ( IOStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'REFCOEF ReadReflectionCoeffcient: ', &
        'Unable to open Internal Reflection Coefficient file'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
    ENDIF

    READ( IRCFile, REC=1, IOSTAT=IOStat, FMT='(A)') line
    READ( line,'(A80,F10.6)' ) line, freq
    READ( IRCFile, REC=2, IOSTAT=IOStat, FMT='(A)') line
    READ( line,'(I10)' ) NkTab

    ALLOCATE( xTab( NkTab ), fTab( NkTab ), gTab( NkTab ), iTab( NkTab ), &
              Stat=IAllocStat )
    IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'REFCOEF ReadReflectionCoeffcient: ', &
        'Too many points in reflection coefficient'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
    ENDIF

    DO ik=1,NkTab
      rec = ik+2
      READ( IRCFile, REC=rec, IOSTAT=IOStat, FMT='(A)') line
      READ( line,'( 5G15.7,I5 )' )  xTab( ik ), fTab( ik ), &
            gTab( ik ), iTab( ik )
    ENDDO

    CLOSE( IRCFile )

  ENDIF

  RETURN
  END !SUBROUTINE ReadReflectionCoefficient

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: InterpolateReflectionCoefficient
! !INTERFACE:
  SUBROUTINE InterpolateReflectionCoefficient( RInt, R, NPts )
! !DESCRIPTION:
!   Interpolates the reflection coefficient for a given angle of incidence.
! Given an angle RInt%ThetaInt, returns the magnitude and
! phase of the reflection coefficient (RInt%R, RInt%phi).

! Uses linear interpolation using the two nearest abscissas
! Assumes phi has been unwrapped so that it varies smoothly.
! I tried modifying it to allow a complex angle of incidence but
! stopped when I realized I needed to fuss with a complex atan2 routine

! !USES: None

! !INPUT PARAMETERS:
! RInt :: Reflection coefficient to be interpolated
! R    :: Reflection coefficient table
! NPts :: Number of points in the reflection coefficient table
  TYPE(ReflectionCoef), INTENT( INOUT ) :: RInt
  TYPE(ReflectionCoef), INTENT( IN    ) :: R( NPts )
  INTEGER,              INTENT( IN    ) :: NPts
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! iLeft, iRight, iMid, iLR :: Loop indices for searching the nearest abscissas
! alpha :: Interpolation factor
! thetaIntr :: Angle of incidence in radians
  INTEGER             :: iLeft, iRight, iMid, iLR
  REAL (KIND=_RL90)   :: alpha, thetaIntr
!EOP

  iLeft  = 1
  iRight = NPts
! thetaIntr should be unnecessary? probably used when I was doing complex angles
  thetaIntr = REAL( RInt%Theta )

  ! Three cases: ThetaInt left, in, or right of tabulated interval

  IF ( thetaIntr.LT.R( iLeft )%theta ) THEN
    !iRight = 2
    RInt%R   = 0.0     ! R( iLeft  )%R
    RInt%phi = 0.0     ! R( iLeft  )%phi
#ifdef IHOP_WRITE_OUT
    ! In adjoint mode we do not write output besides on the first run
    IF ( IHOP_dumpfreq.GE.0 ) THEN
    WRITE( PRTFile, '(2A)' ) 'Warning in InterpolateReflectionCoefficient : ',&
      'Ref. Coef. being set to 0 outside tabulated domain'
    WRITE( PRTFile, * ) 'angle = ', thetaIntr, 'lower limit = ', &
      R( iLeft )%theta
    ENDIF
#endif /* IHOP_WRITE_OUT */

  ELSEIF ( thetaIntr.GT.R( iRight )%theta ) THEN
    !iLeft = NPts - 1
    RInt%R   = 0.0     ! R( iRight )%R
    RInt%phi = 0.0     ! R( iRight )%phi

  ELSE
    ! Search for bracketting abscissas: Log2( NPts ) stabs required for a
    ! bracket
    DO iLR = 1, NPts
      IF ( iLeft.NE.iRight-1 ) THEN
        iMid = ( iLeft + iRight ) / 2
        IF ( R( iMid )%theta.GT.thetaIntr ) THEN
          iRight = iMid
        ELSE
          iLeft  = iMid
        ENDIF
      ENDIF
    ENDDO

    ! Linear interpolation for reflection coef
    alpha    = ( RInt%theta - R( iLeft )%theta ) &
            / ( R( iRight )%theta - R( iLeft )%theta )
    RInt%R   = ( 1 - alpha ) * R( iLeft )%R   + alpha * R( iRight )%R
    RInt%phi = ( 1 - alpha ) * R( iLeft )%phi + alpha * R( iRight )%phi

  ENDIF ! IF ( thetaIntr .LT. R( iLeft )%theta )

  RETURN
  END !SUBROUTINE InterpolateReflectionCoefficient

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: InterpolateIRC
! !INTERFACE:
  SUBROUTINE InterpolateIRC( x, f, g, iPower, xTab, fTab, gTab, iTab, NkTab )
! !DESCRIPTION:
! Interpolates the internal reflection coefficient for a given x.
! Uses polynomial interpolation to approximate the function between the
! tabulated values. The interpolation is done using a polynomial of order N.

! !USES:
    USE poly_mod, only: Poly

! !INPUT PARAMETERS:
! x :: Complex number for which the reflection coefficient is to be interpolated
! f :: Real part of the reflection coefficient
! g :: Imaginary part of the reflection coefficient
! iPower :: Power of 10 to scale the reflection coefficient
! xTab :: Abscissas of the tabulated reflection coefficients
! fTab :: Real parts of the tabulated reflection coefficients
! gTab :: Imaginary parts of the tabulated reflection coefficients
! iTab :: Integer powers associated with the tabulated reflection coefficients
! NkTab :: Number of points in the tabulated reflection coefficients
  INTEGER,                INTENT( IN  ) :: NkTab, iTab( NkTab )
  REAL    ( KIND=_RL90 ), INTENT( IN  ) :: xTab( NkTab )
  COMPLEX ( KIND=_RL90 ), INTENT( IN  ) :: fTab( NkTab ), gTab( NkTab ), x
  INTEGER,                INTENT( OUT ) :: iPower
  COMPLEX ( KIND=_RL90 ), INTENT( OUT ) :: f, g
! !OUTPUT PARAMETERS: iPower, f, g

! !LOCAL VARIABLES:
! N :: Order of the polynomial for interpolation
! i, j, iDel, iLeft, iMid, iRight, NAct :: Loop indices and counters
! xReal :: Real part of the complex number x
! xT, fT, gT :: Arrays for interpolation  
  INTEGER, PARAMETER :: N = 3
  INTEGER                :: i, j, iDel, iLeft, iMid, iRight, NAct
  REAL    ( KIND=_RL90 ) :: xReal
  COMPLEX ( KIND=_RL90 ) :: xT( N ), fT( N ), gT( N )
!EOP

  xReal  = DBLE( x )        ! taking the real part
  iLeft  = 1
  iRight = NkTab

  ! Three cases: x left, in, or right of tabulated interval

  IF (     xReal.LT.xTab( iLeft  ) ) THEN
    f      = fTab( iLeft  )
    g      = gTab( iLeft  )
    iPower = iTab( iLeft  )

  ELSEIF ( xReal.GT.xTab( iRight ) ) THEN
    f      = fTab( iRight )
    g      = gTab( iRight )
    iPower = iTab( iRight )

  ELSE
    ! Search for bracketting abscissas:
    ! Log base 2 (NPts) stabs required for a bracket
    DO WHILE ( iLeft.NE.iRight-1 )
      iMid = ( iLeft + iRight )/2
      IF ( xTab( iMid ).GT.xReal ) THEN
        iRight = iMid
      ELSE
        iLeft  = iMid
      ENDIF
    ENDDO

    ! Extract the subset for interpolation and scale
    iLeft  = MAX( iLeft  - ( N - 2 ) / 2, 1     )
    iRight = MiN( iRight + ( N - 1 ) / 2, NkTab )

    NAct = iRight - iLeft + 1
    DO i = 1, NAct
      j       = i + iLeft - 1
      iDel    = iTab( j ) - iTab( iLeft )
      xT( i ) = xTab( j )
      fT( i ) = fTab( j ) * 10.0D0 ** iDel
      gT( i ) = gTab( j ) * 10.0D0 ** iDel
    ENDDO

    ! Construct the interpolate
    f      = Poly( x, xT, fT, NAct )
    g      = Poly( x, xT, gT, NAct )
    iPower = iTab( iLeft )

  ENDIF ! IF ( xReal .LT. xTab( iLeft  ) )

  RETURN
  END !SUBROUTINE InterpolateIRC

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: writeRefCoef
! !INTERFACE:
  SUBROUTINE writeRefCoef( myThid )
! !DESCRIPTION:
!   Writes information about the reflection coefficients to the output file.

! !USES:
  USE bdry_mod, only: Bdry

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! BotRC, TopRC :: Flags indicating if reflection coefficients are read from files
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
  CHARACTER*(1) :: BotRC, TopRC
!EOP

  ! initiate local parameters
  BotRC = Bdry%Bot%HS%Opt( 1:1 )
  TopRC = Bdry%Top%HS%Opt( 2:2 )


  ! I/O on main thread only
  _BEGIN_MASTER(myThid)

#ifdef IHOP_WRITE_OUT
  IF ( BotRC.EQ.'F' ) THEN
    WRITE( PRTFile,'(A)' ) &
      '__________________________________________________________________________'
    WRITE( PRTFile,'(A)' )
    WRITE( PRTFile,'(A)' ) 'Using tabulated bottom reflection coef.'
    WRITE( PRTFile,'(2A)' ) 'BRCFile=', TRIM( IHOP_fileroot )//'.brc'

    WRITE( PRTFile,'(2A,I10)' ) 'Number of points in bottom reflection ', &
      'coefficient = ', NBotPts
  ENDIF

  IF ( TopRC.EQ.'F' ) THEN
    WRITE( PRTFile,'(A)' ) &
      '__________________________________________________________________________'
    WRITE( PRTFile,'(A)' )
    WRITE( PRTFile,'(A)' ) 'Using tabulated top reflection coef.'
    WRITE( PRTFile,'(2A)' ) 'TRCFile=', TRIM( IHOP_fileroot )//'.trc'

    WRITE( PRTFile,'(2A,I10)' ) 'Number of points in top reflection ', &
      'coefficient = ', NTopPts
  ENDIF

  IF ( BotRC.EQ.'P' ) THEN
    WRITE( PRTFile, '(A)' ) &
      '__________________________________________________________________________'
    WRITE( PRTFile, '(2A,I10)' ) 'Number of points in internal reflection ', &
      'coefficient = ', NkTab
  ENDIF
#endif /* IHOP_WRITE_OUT */

  ! I/O on main thread only
  _END_MASTER(myThid)
  _BARRIER

  RETURN
  END !SUBROUTINE writeRefCoef

END !MODULE refCoef