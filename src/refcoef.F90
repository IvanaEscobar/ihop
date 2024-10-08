#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE refCoef
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

  ! reflection coefficient data

  USE ihop_mod,   only: PRTFile, BRCFile, TRCFile, IRCFile

! ! USES
  implicit none
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
    public ReadReflectionCoefficient, InterpolateReflectionCoefficient, &
           ReflectionCoef, RBot, RTop, NBotPts, NTopPts, writeRefCoef
!=======================================================================

  INTEGER               :: NBotPts, NTopPts
  INTEGER               :: NkTab
  INTEGER,              ALLOCATABLE :: iTab( : )
  REAL    (KIND=_RL90), ALLOCATABLE :: xTab( : )
  COMPLEX (KIND=_RL90), ALLOCATABLE :: fTab( : ), gTab( : )

  TYPE ReflectionCoef
      REAL(KIND=_RL90) :: theta, R, phi
  END TYPE

  TYPE(ReflectionCoef), ALLOCATABLE :: RBot( : ), RTop( : )

CONTAINS
! **************************************************************************** !
! NOTE: To be able to read using direct access, we assume a fixed pad of
! characters per line. So each line in BRCFile/TRCFile will have a fixed length
! of 100 characters. Helps TAF create the adjoint model. IESCO24
! **************************************************************************** !
  SUBROUTINE ReadReflectionCoefficient( myThid )
    USE bdry_mod, only: Bdry

    ! Optionally read in reflection coefficient for Top or Bottom boundary

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
! flag set to 'F' if refl. coef. is to be read from a File
    CHARACTER (LEN=1)   :: BotRC, TopRC
    CHARACTER (LEN=100) :: line
    INTEGER             :: itheta, ik, IOStat, iAllocStat, rec
    REAL  (KIND=_RL90)  :: freq

    ! initiate local parameters
    BotRC = Bdry%Bot%HS%Opt( 1:1 )
    TopRC = Bdry%Top%HS%Opt( 2:2 )

    IF ( ALLOCATED( RTop ) ) DEALLOCATE( RTop )
    IF ( ALLOCATED( RBot ) ) DEALLOCATE( RBot )
    IF ( ALLOCATED( xTab ) ) DEALLOCATE( xTab, fTab, gTab, iTab )

    ! == read any Reflection Coefficients ==
    IF ( BotRC == 'F' ) THEN
      OPEN( FIlE = TRIM( IHOP_fileroot )//'.brc', UNIT=BRCFile, STATUS='OLD',&
            IOSTAT=IOStat, ACTION='READ',ACCESS='direct',RECL=100 )
      IF ( IOStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)')'REFCOEF ReadReflectionCoeffcient: ',&
                      'Unable to open Bottom Reflection Coefficient file'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
      END IF

      READ( BRCFile,REC=1,IOSTAT=IOStat,FMT='(A)') line
      READ( line,'(I10)' ) NBotPts

      ALLOCATE( RBot( NBotPts ), Stat = IAllocStat )
      IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'REFCOEF ReadReflectionCoeffcient: ', &
                    'Insufficient memory for bot. refl. coef.: reduce # points'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
      END IF

      DO itheta=1,NBotPts
        rec = itheta+1
        READ( BRCFile,REC=rec,IOSTAT=IOStat,FMT='(A)' ) line
        READ( line,'(3F5.3)' ) RBot(itheta)%theta, RBot(itheta)%R, &
                               RBot(itheta)%phi
      ENDDO
      CLOSE( BRCFile )
      RBot%phi = deg2rad * RBot%phi   ! convert to radians

    ELSE   ! should allocate something anyway, since variable is passed
      ALLOCATE(  RBot( 1 ), Stat = IAllocStat )
      NBotPts = 1         !RG
      RBot(1)%theta = 0.  !RG
      RBot(1)%R = 0.      !RG
      RBot(1)%phi = 0.    !RG

    ENDIF ! BotRC == ...


    ! Optionally read in top reflection coefficient
    IF ( TopRC == 'F' ) THEN
      OPEN( FIlE = TRIM( IHOP_fileroot )//'.trc', UNIT=TRCFile, STATUS='OLD',&
            IOSTAT=IOStat, ACTION='READ',ACCESS='direct',RECL=100 )
      IF ( IOStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'REFCOEF ReadReflectionCoeffcient: ', &
                             'Unable to open Top Reflection Coefficient file'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
      END IF

      READ( TRCFile,REC=1,IOSTAT=IOStat,FMT='(A)') line
      READ( line,'(I10)' ) NTopPts

      ALLOCATE( RTop( NTopPts ), Stat = IAllocStat )
      IF ( iAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'REFCOEF ReadReflectionCoeffcient: ', &
                    'Insufficient memory for top refl. coef.: reduce # points'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
      END IF

      DO itheta=1,NTopPts
        rec = itheta+1
        READ( TRCFile,REC=rec,IOSTAT=IOStat,FMT='(A)' ) line
        READ( line,'(3F5.3)' ) RTop(itheta)%theta, RTop(itheta)%R, &
                               RTop(itheta)%phi
      ENDDO
      CLOSE( TRCFile )
      RTop%phi = deg2rad *  RTop%phi   ! convert to radians

    ELSE   ! should allocate something anyway, since variable is passed
      ALLOCATE( RTop( 1 ), Stat = iAllocStat )
      NTopPts = 1         !RG
      RTop(1)%theta = 0.  !RG
      RTop(1)%R = 0.      !RG
      RTop(1)%phi = 0.    !RG

    ENDIF ! TopRC == ...


    ! Optionally read in internal reflection coefficient data
    IF ( BotRC == 'P' ) THEN
      OPEN( FIlE = TRIM( IHOP_fileroot )//'.irc', UNIT=TRCFile, STATUS='OLD',&
            IOSTAT=IOStat, ACTION='READ',ACCESS='direct',RECL=100 )
      IF ( IOStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'REFCOEF ReadReflectionCoeffcient: ', &
                        'Unable to open Internal Reflection Coefficient file'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
      END IF

      READ( IRCFile,REC=1,IOSTAT=IOStat,FMT='(A)') line
      READ( line,'(A80,F10.6)' ) line, freq
      READ( IRCFile,REC=2,IOSTAT=IOStat,FMT='(A)') line
      READ( line,'(I10)' ) NkTab

      ALLOCATE( xTab( NkTab ), fTab( NkTab ), gTab( NkTab ), iTab( NkTab ), &
                Stat = iAllocStat )
      IF ( iAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'REFCOEF ReadReflectionCoeffcient: ', &
                               'Too many points in reflection coefficient'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadReflectionCoefficient'
      END IF

      DO ik=1,NkTab
        rec = ik+2
        READ( IRCFile,REC=rec,IOSTAT=IOStat,FMT='(A)') line
        READ( line,'( 5G15.7,I5 )' )  xTab( ik ), fTab( ik ), &
             gTab( ik ), iTab( ik )
      ENDDO
      CLOSE( IRCFile )
    ENDIF

  RETURN
  END !SUBROUTINE ReadReflectionCoefficient

! **************************************************************************** !
  SUBROUTINE InterpolateReflectionCoefficient( RInt, R, NPts )

    ! Given an angle RInt%ThetaInt, returns the magnitude and
    ! phase of the reflection coefficient (RInt%R, RInt%phi).

    ! Uses linear interpolation using the two nearest abscissas
    ! Assumes phi has been unwrapped so that it varies smoothly.
    ! I tried modifying it to allow a complex angle of incidence but
    ! stopped when I realized I needed to fuss with a complex atan2 routine

    INTEGER,              INTENT( IN    ) :: NPts       ! # pts in refl. coef.
    TYPE(ReflectionCoef), INTENT( IN    ) :: R( NPts )  ! Reflection coefficient table
    TYPE(ReflectionCoef), INTENT( INOUT ) :: RInt       ! interpolated value of refl. coef.
    INTEGER             :: iLeft, iRight, iMid, iLR
    REAL (KIND=_RL90)   :: alpha, Thetaintr

    iLeft  = 1
    iRight = NPts
 ! thetaIntr should be unnecessary? probably used when I was doing complex angles
    thetaIntr = REAL( RInt%Theta )

    ! Three cases: ThetaInt left, in, or right of tabulated interval

    IF ( thetaIntr < R( iLeft )%theta ) THEN
       !iRight = 2
       RInt%R   = 0.0     ! R( iLeft  )%R
       RInt%phi = 0.0     ! R( iLeft  )%phi
#ifdef IHOP_WRITE_OUT
       ! In adjoint mode we do not write output besides on the first run
       IF (IHOP_dumpfreq.GE.0) THEN
        WRITE( PRTFile, * ) 'Warning in InterpolateReflectionCoefficient : ',&
                            'Ref. Coef. being set to 0 outside tabulated domain'
        WRITE( PRTFile, * ) 'angle = ', thetaintr, 'lower limit = ', &
                            R( iLeft)%theta
       ENDIF
#endif /* IHOP_WRITE_OUT */

    ELSE IF( thetaIntr > R( iRight )%theta ) THEN
       !iLeft = NPts - 1
       RInt%R   = 0.0     ! R( iRight )%R
       RInt%phi = 0.0     ! R( iRight )%phi

    ELSE
       ! Search for bracketting abscissas: Log2( NPts ) stabs required for a
       ! bracket

       DO iLR = 1, NPts
        IF ( iLeft /= iRight - 1 ) THEN
          iMid = ( iLeft + iRight ) / 2
          IF ( R( iMid )%theta > thetaIntr ) THEN
             iRight = iMid
          ELSE
             iLeft  = iMid
          ENDIF
        ENDIF
       END DO

       ! Linear interpolation for reflection coef

       alpha    = ( RInt%theta - R( iLeft )%theta ) &
                / ( R( iRight )%theta - R( iLeft )%theta )
       RInt%R   = ( 1 - alpha ) * R( iLeft )%R   + alpha * R( iRight )%R
       RInt%phi = ( 1 - alpha ) * R( iLeft )%phi + alpha * R( iRight )%phi

    ENDIF

  RETURN
  END !SUBROUTINE InterpolateReflectionCoefficient

! **************************************************************************** !
  SUBROUTINE InterpolateIRC( x, f, g, iPower, xTab, fTab, gTab, iTab, NkTab )

    ! Internal reflection coefficient interpolator.
    ! Returns f, g, iPower for given x using tabulated values.
    ! Uses polynomial interpolation to approximate the function between the
    ! tabulated values

    USE poly_mod, only: Poly

    INTEGER, PARAMETER :: N = 3     ! order of the polynomial for interpolation
    INTEGER,                INTENT( IN  ) :: NkTab, iTab( NkTab )
    REAL    ( KIND=_RL90 ), INTENT( IN  ) :: xTab( NkTab )
    COMPLEX ( KIND=_RL90 ), INTENT( IN  ) :: fTab( NkTab ), gTab( NkTab ), x
    INTEGER,                INTENT( OUT ) :: iPower
    COMPLEX ( KIND=_RL90 ), INTENT( OUT ) :: f, g
    INTEGER                :: i, j, iDel, iLeft, iMid, iRight, NAct
    REAL    ( KIND=_RL90 ) :: xReal
    COMPLEX ( KIND=_RL90 ) :: xT( N ), fT( N ), gT( N )

    xReal  = DBLE( x )        ! taking the real part
    iLeft  = 1
    iRight = NkTab

    ! Three cases: x left, in, or right of tabulated interval

    IF (     xReal < xTab( iLeft  ) ) THEN
       f      = fTab( iLeft  )
       g      = gTab( iLeft  )
       iPower = iTab( iLeft  )
    ELSE IF( xReal > xTab( iRight ) ) THEN
       f      = fTab( iRight )
       g      = gTab( iRight )
       iPower = iTab( iRight )
    ELSE

       ! Search for bracketting abscissas:
       ! Log base 2 (NPts) stabs required for a bracket

       DO WHILE ( iLeft /= iRight-1 )
          iMid = ( iLeft + iRight )/2
          IF ( xTab( iMid ) > xReal ) THEN
             iRight = iMid
          ELSE
             iLeft  = iMid
          ENDIF
       END DO

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
       END DO

       ! Construct the interpolate

       f      = Poly( x, xT, fT, NAct )
       g      = Poly( x, xT, gT, NAct )
       iPower = iTab( iLeft )
    ENDIF

  RETURN
  END !SUBROUTINE InterpolateIRC

! **************************************************************************** !
  SUBROUTINE writeRefCoef( myThid )
    USE bdry_mod, only: Bdry

    ! Optionally read in reflection coefficient for Top or Bottom boundary

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
! flag set to 'F' if refl. coef. is to be read from a File
    CHARACTER*(1) :: BotRC, TopRC

    ! initiate local parameters
    BotRC = Bdry%Bot%HS%Opt( 1:1 )
    TopRC = Bdry%Top%HS%Opt( 2:2 )


  ! I/O on main thread only
  _BEGIN_MASTER(myThid)

#ifdef IHOP_WRITE_OUT
  ! In adjoint mode we do not write output besides on the first run
  IF (IHOP_dumpfreq.GE.0) THEN

    IF ( BotRC == 'F' ) THEN
      WRITE( PRTFile,'(2A)' ) '_______________________________________', &
                              '___________________________________'
      WRITE( PRTFile,'(A)' )
      WRITE( PRTFile,'(A)' ) 'Using tabulated bottom reflection coef.'
      WRITE( PRTFile,'(2A)' ) 'BRCFile=', TRIM( IHOP_fileroot )//'.brc'

      WRITE( PRTFile,'(2A,I10)' ) 'Number of points in bottom reflection ', &
                           'coefficient = ', NBotPts
    END IF

    IF ( TopRC == 'F' ) THEN
      WRITE( PRTFile,'(2A)' ) '_______________________________________', &
                              '___________________________________'
      WRITE( PRTFile,'(A)' )
      WRITE( PRTFile,'(A)' ) 'Using tabulated top reflection coef.'
      WRITE( PRTFile,'(2A)' ) 'TRCFile=', TRIM( IHOP_fileroot )//'.trc'

      WRITE( PRTFile,'(2A,I10)' ) 'Number of points in top reflection ', &
                           'coefficient = ', NTopPts
    END IF

    IF ( BotRC == 'P' ) THEN
      WRITE( PRTFile, * )
      WRITE( PRTFile, * ) 'Number of points in internal reflection ', &
                          'coefficient = ', NkTab
    END IF

  END IF ! don't write in adjoint mode
#endif /* IHOP_WRITE_OUT */
  ! I/O on main thread only
  _END_MASTER(myThid)
  _BARRIER

  RETURN
  END !SUBROUTINE writeRefCoef

END MODULE refCoef
