#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: bdry_mod
MODULE bdry_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!  Loads altimetry (top bdry) and bathymetry (bottom bdry) data

! !USES:
  IMPLICIT NONE
! == Global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC  initATI, initBTY, GetTopSeg, GetBotSeg, Bot, Top, &
          iSegTop, iSegBot, rTopSeg, rBotSeg, &
          atiType, btyType, HSInfo, Bdry, writeBdry
!=======================================================================

! == Module variables ==
  INTEGER, PARAMETER :: Number_to_Echo = 21
  INTEGER, PROTECTED :: NBtyPts = 2, NAtiPts = 2
  INTEGER            :: iSegTop, iSegBot ! indices to current active segment

  ! range intervals defining the current active segment
  REAL (KIND=_RL90) :: rTopSeg( 2 ), rBotSeg( 2 )
  CHARACTER*(2)     :: atiType = 'LS', btyType = 'LS'

! == Derived types ==
  TYPE HSInfo
    ! This type contains the information about the halfspace
    ! compressional and shear wave speeds/attenuations in user units
    REAL     (KIND=_RL90) :: alphaR, alphaI, betaR, betaI
    REAL     (KIND=_RL90) :: rho, Depth  ! density, depth
    COMPLEX  (KIND=_RL90) :: cP, cS      ! P-wave, S-wave speeds
    CHARACTER*(1)         :: BC          ! Boundary condition type
    CHARACTER*(6)         :: Opt
  END TYPE

  TYPE BdryPt
    REAL (KIND=_RL90) :: x( 2 ), &     ! segment coordinate
                         t( 2 ), &     ! segment tangent
                         n( 2 )        ! segment outward
    REAL (KIND=_RL90) :: Len, Kappa    ! length and curvature of a segement
    ! For the curvilinear grid option
    REAL (KIND=_RL90) :: Nodet( 2 ), & ! tangent at the node
                         Noden( 2 )    ! normal at the node
    REAL (KIND=_RL90) :: Dx, Dxx, &    ! 1st, 2nd derivatives wrt depth
                         Dss           ! derivative along tangent
    TYPE( HSInfo )    :: HS
  END TYPE

  TYPE(BdryPt), ALLOCATABLE :: Top( : ), Bot( : )

  TYPE BdryPt2
    TYPE( HSInfo )   :: HS
  END TYPE

  TYPE BdryType
    TYPE( BdryPt2 )   :: Top, Bot
  END TYPE BdryType

  TYPE(BdryType) :: Bdry
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R initATI
! S/R initBTY
! S/R ComputeBdryTangentNormal
! S/R GetTopSeg
! S/R GetBotSeg
! S/R WriteBdry
! S/R WriteTopBot
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: initATI
! !INTERFACE:
  SUBROUTINE initATI( TopATI, DepthT, myThid )
! !DESCRIPTION:
!  Reads in the top altimetry

! !USES:
    USE ihop_mod,        only: ATIFile
    USE monotonic_mod,   only: monotonic
    USE atten_mod,       only: CRCI
  ! fT = 1000 ONLY for acousto-elastic halfspaces, I will have to pass this
  ! parameter in a different way after ssp_mod is split btwn fixed and varia
  !USE initenvihop, only: fT

! !INPUT PARAMETERS:
! TopATI :: Type of top altimetry
! DepthT :: Depth at which to read altimetry
! myThid :: my thread ID
  CHARACTER*(1),     INTENT( IN ) :: TopATI
  REAL (KIND=_RL90), INTENT( IN ) :: DepthT
  INTEGER,           INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf     :: Informational/error message buffer
! ii         :: Loop index for reading altimetry points
! IOStat     :: I/O status
! IAllocStat :: Allocation status
! x          :: Array of altimetry points
! bPower     :: Power for attenuation conversion
! fT         :: Frequency for attenuation conversion
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: ii
  INTEGER :: IOStat, IAllocStat
  REAL (KIND=_RL90), ALLOCATABLE :: x(:)
  REAL (KIND=_RL90) :: bPower, fT
!EOP

  ! IESCO24 fT init
  bPower = 1.0
  fT     = 1000.0
  IF (ALLOCATED(Top)) DEALLOCATE(Top)


  ! Either read from an ATIFile or set +/- infty
  SELECT CASE ( TopATI )
  CASE ( '~', '*' )
    OPEN( UNIT=ATIFile, FILE=TRIM(IHOP_fileroot) // '.ati', &
          STATUS='OLD', IOSTAT=IOStat, ACTION='READ' )
    IF ( IOSTAT.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
        'Unable to open altimetry file'
      ! In adjoint mode we do not write output besides on the first run
      IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R initATI'
    ENDIF ! IF ( IOSTAT.NE.0 )

    READ( ATIFile, * ) atiType

    SELECT CASE ( atiType( 1:1 ) )
    CASE ( 'C', 'L' )
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
        'Unknown option for selecting altimetry interpolation'
      ! In adjoint mode we do not write output besides on the first run
      IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R initATI'
    END SELECT ! atiType( 1:1 )

    READ( ATIFile, * ) NAtiPts
    ! Extending altimetry to infinity to the left and right
    NAtiPts = NAtiPts + 2

    ALLOCATE( Top( NAtiPts ), Stat=IAllocStat )
    IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
        'Insufficient memory for altimetry data: reduce # ati points'
      ! In adjoint mode we do not write output besides on the first run
      IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R initATI'
    ENDIF ! IF ( IAllocStat.NE.0 )

    ! Read in the altimetry
    atiPt: DO ii = 2, NAtiPts - 1

      SELECT CASE ( atiType( 2:2 ) )
      CASE ( 'S', '' )
        READ( ATIFile, * ) Top( ii )%x
        ! init to default variables
        Top( ii )%HS%alphaR = -1.
        Top( ii )%HS%alphaI = -1.
        Top( ii )%HS%betaR = -1.
        Top( ii )%HS%betaI = -1.
        Top( ii )%HS%rho = -1.

      CASE ( 'L' )
        READ(  ATIFile, * ) Top( ii )%x, Top( ii )%HS%alphaR, &
                            Top( ii )%HS%betaR, Top( ii )%HS%rho, &
                            Top( ii )%HS%alphaI, Top( ii )%HS%betaI

      CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
          'Unknown option for selecting altimetry option'
        ! In adjoint mode we do not write output besides on the first run
        IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R initATI'
      END SELECT ! atiType( 2:2 )

      IF ( Top( ii )%x( 2 ).LT.DepthT ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
          'Altimetry rises above highest point in the sound speed profile'
        ! In adjoint mode we do not write output besides on the first run
        IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R initATI'
      ENDIF ! IF ( Top( ii )%x( 2 ) < DepthT )

    ENDDO atiPt

    CLOSE( ATIFile )
    Top(:)%x(1) = 1000.0 * Top(:)%x(1)   ! Convert ranges in km to m

  CASE DEFAULT   ! no altimetry given, use Grid depth for flat top
    NAtiPts = 2
    ALLOCATE( Top( NAtiPts ), Stat = IAllocStat )
    IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
        'Insufficient memory for altimetry data'
      ! In adjoint mode we do not write output besides on the first run
      IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R initATI'
    ENDIF ! IF ( IAllocStat.NE.0 )

    DO ii=1, NAtiPts
      Top( ii )%x = [ -sqrt( huge( Top( 1 )%x( 1 ) ) ) / 1.0d5, DepthT ]
      ! init to default variables
      Top( ii )%HS%alphaR = -1.
      Top( ii )%HS%alphaI = -1.
      Top( ii )%HS%betaR = -1.
      Top( ii )%HS%betaI = -1.
      Top( ii )%HS%rho = -1.
    ENDDO ! DO ii=1, NAtiPts

    Top(NAtiPts)%x(1) = -1 * Top(1)%x(1)

  END SELECT ! TopATI

  ! Set Top(1:NAtiPts)%t,n,Len,Nodet,Noden,kappa,Dx,Dxx,Dss
  CALL ComputeBdryTangentNormal( Top, 'Top', NAtiPts )

  ! Check if monotonic
  IF ( ALLOCATED(x) ) DEALLOCATE(x)
  ALLOCATE( x(NAtiPts) )
  x = Top%x(1)
  IF ( .NOT.monotonic( x, NAtiPts ) ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
      'Altimetry ranges are not monotonically increasing'
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R initATI'
  ENDIF ! IF ( .NOT.monotonic( x, NAtiPts ) )

  ! Set Top%HS%cP,cS
  DO ii = 1, NAtiPts
    ! compressional wave speed
    Top( ii )%HS%cP = -999.
    ! shear wave speed
    Top( ii )%HS%cS = -999.
  ENDDO
  ! convert range-dependent geoacoustic parameters from user to program units
  ! W is dB/wavelength
  IF ( atiType( 2:2 ).EQ.'L' ) THEN
    DO ii = 1, NAtiPts
      ! compressional wave speed
      Top( ii )%HS%cP = CRCI( 1D20, Top( ii )%HS%alphaR, &
                              Top( ii )%HS%alphaI, 'W ', bPower, fT, myThid )
      ! shear wave speed
      Top( ii )%HS%cS = CRCI( 1D20, Top( ii )%HS%betaR,  &
                              Top( ii )%HS%betaI, 'W ', bPower, fT,  myThid )
    ENDDO
  ENDIF ! IF ( atiType( 2:2 ).EQ.'L' )

  ! Set defaults for rest of Top derived type
  Top%HS%Depth = -999.
  Top%HS%BC  = ''
  Top%HS%Opt = ''

  RETURN
  END !SUBROUTINE initATI

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: initBTY
! !INTERFACE:
  SUBROUTINE initBTY( BotBTY, DepthB, myThid )
! !DESCRIPTION:
! Initializes the bottom bathymetry data structure

! !USES:
    USE ihop_mod,      only: BTYFile
    USE monotonic_mod, only: monotonic
    USE atten_mod,     only: CRCI
  ! fT = 1000 ONLY for acousto-elastic halfspaces, I will have to pass this
  ! parameter in a different way after ssp_mod is split btwn fixed and varia
  !USE initenvihop, only: fT

! !INPUT PARAMETERS:
! BotBTY :: Type of bottom bathymetry
! DepthB :: Depth at which to read bathymetry
! myThid :: my thread ID
  CHARACTER*(1),     INTENT( IN ) :: BotBTY
  REAL (KIND=_RL90), INTENT( IN ) :: DepthB
  INTEGER,           INTENT( IN ) :: myThid
!!OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf     :: Informational/error message buffer
! i,j,bi,bj  :: Loop indices for reading bathymetry points
! ii         :: Loop index for reading bathymetry points
! IOStat     :: I/O status
! IAllocStat :: Allocation status
! gcmbathy   :: Array of bathymetry points from GCM
! gcmmin     :: Minimum bathymetry point from GCM
! gcmmax     :: Maximum bathymetry point from GCM
! x          :: Array of bathymetry points
! firstnonzero :: Flag to check if first non-zero bathymetry point
! bPower     :: Power for attenuation conversion
! fT         :: Frequency for attenuation conversion
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: i,j,bi,bj,ii
  INTEGER :: IOStat, IAllocStat
  REAL (KIND=_RL90) :: gcmbathy(sNx,sNy), gcmmin, gcmmax
  REAL (KIND=_RL90), ALLOCATABLE :: x(:)
  LOGICAL :: firstnonzero
  REAL (KIND=_RL90)   :: bPower, fT
!EOP

  ! IESCO24 fT init
  bPower = 1.0
  fT     = 1000.0
  IF (ALLOCATED(Bot)) DEALLOCATE(Bot)


  ! Either read from an BTYFile or set +/- infty
  SELECT CASE ( BotBTY )
  CASE ( '~', '*' )
    OPEN( UNIT=BTYFile, FILE=TRIM(IHOP_fileroot) // '.bty', &
          STATUS='OLD', IOSTAT=IOStat, ACTION='READ' )
    IF ( IOSTAT.NE.0 ) THEN
# ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
        'Unable to open bathymetry file'
      ! In adjoint mode we do not write output besides on the first run
      IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R initBTY'
    ENDIF ! IF ( IOSTAT.NE.0 )

    READ( BTYFile, * ) btyType

    SELECT CASE ( btyType( 1:1 ) )
    CASE ( 'C','L' )
    CASE DEFAULT
# ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
        'Unknown option for selecting bathymetry interpolation'
      ! In adjoint mode we do not write output besides on the first run
      IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R initBTY'
    END SELECT ! btyType( 1:1 )

    READ( BTYFile, * ) NBtyPts
    ! Extend bathymetry to infinity on both sides
    NBtyPts = NBtyPts + 2

    ALLOCATE( Bot( NBtyPts ), Stat=IAllocStat )
    IF ( IAllocStat.NE.0 ) THEN
# ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
        'Insufficient memory for bathymetry data: reduce # bty points'
      ! In adjoint mode we do not write output besides on the first run
      IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R initBTY'
    ENDIF ! IF ( IAllocStat.NE.0 )


    SELECT CASE ( btyType( 2:2 ) )
    CASE ( 'S', '','L' )
    CASE DEFAULT
# ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
        'Unknown option for selecting bathymetry interpolation'
      ! In adjoint mode we do not write output besides on the first run
      IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R initBTY'
    END SELECT ! btyType( 2:2 )

    ! R_low check of GCM depths, change to positive values
    firstnonzero=.true.
    DO bj=myByLo(myThid),myByHi(myThid)
      DO bi=myBxLo(myThid),myBxHi(myThid)
        gcmbathy(1:sNx,1:sNy) = rkSign*R_low(1:sNx,1:sNy,bi,bj)
        DO i=1,sNx
          DO j=1,sNy
            IF ( gcmbathy(i,j).NE.0.0 ) THEN
              IF (firstnonzero) THEN
                gcmmin = gcmbathy(i,j)
                gcmmax = gcmbathy(i,j)
                firstnonzero=.false.
              ELSE
                IF ( gcmbathy(i,j).LT.gcmmin ) gcmmin = gcmbathy(i,j)
                IF ( gcmbathy(i,j).GT.gcmmax ) gcmmax = gcmbathy(i,j)
              END IF
            END IF
          END DO
        END DO
      END DO
    END DO

    ! Read in the bathymetry
    btyPt: DO ii = 2, NBtyPts - 1

      SELECT CASE ( btyType( 2:2 ) )
      CASE ( 'S', '' )   ! short format
        READ( BTYFile, * ) Bot( ii )%x
        ! init to default variables
        Bot( ii )%HS%alphaR = -1.
        Bot( ii )%HS%alphaI = -1.
        Bot( ii )%HS%betaR = -1.
        Bot( ii )%HS%betaI = -1.
        Bot( ii )%HS%rho = -1.
# ifdef IHOP_WRITE_OUT
        IF (Bot(ii)%x(2).LT.gcmmin .OR. Bot(ii)%x(2).GT.gcmmax) THEN
          WRITE(msgBuf,'(2A,F10.2,A,F10.2)')   &
            '** Warning ** BDRYMOD initBTY: ', &
            'ihop and gcm bathymetry vary, ihop:', Bot(ii)%x(2), &
            'gcm:', gcmmax
          ! In adjoint mode we do not write output besides on first run
          IF (IHOP_dumpfreq.GE.0) &
            CALL PRINT_MESSAGE(msgBuf, errorMessageUnit, &
                                SQUEEZE_RIGHT, myThid)
        ENDIF
# endif /* IHOP_WRITE_OUT */

      CASE ( 'L' )       ! long format
        READ( BTYFile, * ) Bot( ii )%x, Bot( ii )%HS%alphaR,    &
                          Bot( ii )%HS%betaR, Bot( ii )%HS%rho, &
                          Bot( ii )%HS%alphaI, Bot( ii )%HS%betaI

      CASE DEFAULT
# ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
          'Unknown option for selecting bathymetry interpolation'
        ! In adjoint mode we do not write output besides on first run
        IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R initBTY'
      END SELECT ! btyType( 2:2 )

      IF ( Bot( ii )%x( 2 ).GT.DepthB ) THEN
# ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
          'Bathymetry drops below lowest point in the sound speed profile'
        ! In adjoint mode we do not write output besides on first run
        IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
          STOP 'ABNORMAL END: S/R initBTY'
      ENDIF ! IF ( Bot( ii )%x( 2 ).GT.DepthB )
    END DO btyPt

    CLOSE( BTYFile )
    Bot(:)%x(1) = 1000.0 * Bot(:)%x(1)   ! Convert ranges in km to m

  CASE DEFAULT   ! no bathymetry given, use Grid depth for flat bottom
    NBtyPts = 2
    ALLOCATE( Bot( NBtyPts ), Stat = IAllocStat )
    IF ( IAllocStat.NE.0 ) THEN
#  ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
        'Insufficient memory for bathymetry data'
      ! In adjoint mode we do not write output besides on first run
      IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#  endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R initBTY'
    ENDIF ! IF ( IAllocStat.NE.0 )

    DO ii=1, NBtyPts
      Bot( ii )%x = [ -sqrt( huge( Bot( 1 )%x( 1 ) ) ) / 1.0d5, DepthB ]
      ! init to default variables
      Bot( ii )%HS%alphaR = -1.
      Bot( ii )%HS%alphaI = -1.
      Bot( ii )%HS%betaR = -1.
      Bot( ii )%HS%betaI = -1.
      Bot( ii )%HS%rho = -1.
    END DO
    Bot(NBtyPts)%x(1) = -1 * Bot(1)%x(1)

  END SELECT ! BotBTY

  ! Set Bot(1:NBtyPts)%t,n,Len,Nodet,Noden,kappa,Dx,Dxx,Dss
  CALL ComputeBdryTangentNormal( Bot, 'Bot', NBtyPts )

  IF ( ALLOCATED(x) ) DEALLOCATE(x)
  ALLOCATE( x(NBtyPts) )
  x = Bot%x(1)
  IF ( .NOT.monotonic( x, NBtyPts ) ) THEN
# ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
      'Bathymetry ranges are not monotonically increasing'
    ! In adjoint mode we do not write output besides on first run
    IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R initBTY'
  ENDIF ! IF ( .NOT.monotonic( x, NBtyPts ) )

  ! Initiate Bot
  DO ii = 1, NBtyPts
    ! compressional wave speed
    Bot( ii )%HS%cP = -999.
    ! shear wave speed
    Bot( ii )%HS%cS = -999.
  ENDDO
  ! convert range-dependent geoacoustic parameters from user to program units
  ! W is dB/wavelength
  IF ( btyType( 2:2 ) == 'L' ) THEN
    DO ii = 1, NBtyPts
      ! compressional wave speed
      Bot( ii )%HS%cP = CRCI( 1D20, Bot( ii )%HS%alphaR, &
                              Bot( ii )%HS%alphaI, 'W ', bPower, fT, myThid )
      ! shear wave speed
      Bot( ii )%HS%cS = CRCI( 1D20, Bot( ii )%HS%betaR,  &
                              Bot( ii )%HS%betaI,  'W ', bPower, fT, myThid )
    ENDDO
  ENDIF

  ! Set defaults for rest of Bot derived type
  Bot%HS%Depth = -999.
  Bot%HS%BC  = ''
  Bot%HS%Opt = ''

  RETURN
  END !SUBROUTINE initBTY

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ComputeBdryTangentNormal
! !INTERFACE:
  SUBROUTINE ComputeBdryTangentNormal( Bdry, BotTop, NPts )
! !DESCRIPTION:
! Computes tangent and outward-pointing normal to each segment of the
! boundary defined by the user.
! The boundary is also extended with a constant depth to infinity to cover
! cases where a ray leaves the domain.

! !USES: None

! !INPUT PARAMETERS:
! Bdry    :: Array of boundary points
! BotTop  :: Flag indicating bottom or top reflection
! NPts    :: Number of points in the boundary
  TYPE(BdryPt),  INTENT( INOUT ) :: Bdry( : )
  CHARACTER*(3), INTENT( IN )    :: BotTop
  INTEGER,       INTENT( IN )    :: NPts
! !OUTPUT PARAMETERS: Bdry(?)

! !LOCAL VARIABLES:
! phi        :: Array of bathymetry points
! sss        :: Temporary variable for averaging normals
! CurvilinearFlag :: Flag indicating curvilinear grid option
! ii         :: Loop index for computing tangent and normal
! IAllocStat :: Allocation status
  REAL (KIND=_RL90), ALLOCATABLE :: phi( : )
  REAL (KIND=_RL90)              :: sss
  CHARACTER*(2)                  :: CurvilinearFlag
  INTEGER                        :: ii
  INTEGER                        :: IAllocStat
!EOP

  SELECT CASE ( BotTop )
  CASE ( 'Bot' )
    CurvilinearFlag = btyType
  CASE ( 'Top' )
    CurvilinearFlag = atiType
  CASE DEFAULT
    ! Do nothing
    CurvilinearFlag = '-'
    STOP 'ABNORMAL END: S/R computetangentnormal'
  END SELECT

  ! extend the bathymetry to +/- infinity in a piecewise constant fashion

  Bdry( 1    )%x( 1 ) = -sqrt( huge( Bdry( 1 )%x( 1 ) ) ) / 1.0d5
  Bdry( 1    )%x( 2 ) = Bdry( 2        )%x( 2 )
  Bdry( 1    )%HS     = Bdry( 2        )%HS
  Bdry( NPts )%x( 1 ) = +sqrt( huge( Bdry( 1 )%x( 1 ) ) ) / 1.0d5
  Bdry( NPts )%x( 2 ) = Bdry( NPts - 1 )%x( 2 )
  Bdry( NPts )%HS     = Bdry( NPts - 1 )%HS

  ! compute tangent and outward-pointing normal to each bottom segment
  ! tBdry( 1, : ) = xBdry( 1, 2:NPts ) - xBdry( 1, 1:NPts - 1 )
  ! tBdry( 2, : ) = xBdry( 2, 2:NPts ) - xBdry( 2, 1:NPts - 1 )
  ! above caused compiler problems

  BoundaryPt: DO ii = 1, NPts - 1
    Bdry( ii )%t   = Bdry( ii + 1 )%x      - Bdry( ii )%x
    Bdry( ii )%Dx  = Bdry( ii )%t( 2 ) / Bdry( ii )%t( 1 ) ! 1st derivative

    ! normalize the tangent vector
    IF ( ALL(Bdry(ii)%t.EQ.0.0) ) THEN
      Bdry( ii )%Len = 0.0
      Bdry( ii )%t   = 0.0 * Bdry( ii )%t
    ELSE
      Bdry( ii )%Len = NORM2( Bdry( ii )%t )
      Bdry( ii )%t   = Bdry( ii )%t / Bdry( ii )%Len
    ENDIF

    SELECT CASE ( BotTop )
    CASE ( 'Bot' )
      Bdry( ii )%n( 1 ) = -Bdry( ii )%t( 2 )
      Bdry( ii )%n( 2 ) = +Bdry( ii )%t( 1 )
    CASE ( 'Top' )
      Bdry( ii )%n( 1 ) = +Bdry( ii )%t( 2 )
      Bdry( ii )%n( 2 ) = -Bdry( ii )%t( 1 )
    CASE DEFAULT
      ! Do nothing
      Bdry( ii )%n = 0
      STOP 'ABNORMAL END: S/R computetangentnormal'
    END SELECT

  END DO BoundaryPt

  ! curvilinear option: compute tangent and normal at node by averaging
  ! normals on adjacent segments
  IF ( CurvilinearFlag( 1:1 ).EQ.'C' ) THEN
    ! averaging two centered differences is equivalent to forming a single
    ! centered difference of two steps ...
    DO ii = 2, NPts - 1
      sss = Bdry( ii-1 )%Len / ( Bdry( ii-1 )%Len + Bdry( ii )%Len )
      sss = 0.5
      Bdry( ii )%Nodet = ( 1.0 - sss ) * Bdry( ii-1 )%t + sss * Bdry( ii )%t
    END DO

    Bdry( 1    )%Nodet = [ 1.0, 0.0 ]   ! tangent left-end node
    Bdry( NPts )%Nodet = [ 1.0, 0.0 ]   ! tangent right-end node

    SELECT CASE ( BotTop )
    CASE ( 'Bot' )
      Bdry( : )%Noden( 1 ) = -Bdry( : )%Nodet( 2 )
      Bdry( : )%Noden( 2 ) = +Bdry( : )%Nodet( 1 )
    CASE ( 'Top' )
      Bdry( : )%Noden( 1 ) = +Bdry( : )%Nodet( 2 )
      Bdry( : )%Noden( 2 ) = -Bdry( : )%Nodet( 1 )
    CASE DEFAULT
      ! Do nothing
      Bdry( : )%Noden(1) = 0
      Bdry( : )%Noden(2) = 0
      STOP 'ABNORMAL END: S/R computetangentnormal'
    END SELECT

    ! compute curvature in each segment
    IF (ALLOCATED(phi)) DEALLOCATE(phi)
    ALLOCATE( phi( NPts ), Stat = IAllocStat )
    ! phi is the angle at each node
    phi = atan2( Bdry( : )%Nodet( 2 ), Bdry( : )%Nodet( 1 ) )

    DO ii = 1, NPts - 1
      ! this is curvature = dphi/ds
      Bdry( ii )%kappa = ( phi( ii+1 ) - phi( ii ) ) / Bdry( ii )%Len
      ! second derivative
      Bdry( ii )%Dxx   = ( Bdry( ii+1 )%Dx     - Bdry( ii )%Dx     ) / &
                         ( Bdry( ii+1 )%x( 1 ) - Bdry( ii )%x( 1 ) )
      ! derivative in direction of tangent
      Bdry( ii )%Dss   = Bdry( ii )%Dxx * Bdry( ii )%t( 1 )**3

      Bdry( ii )%kappa = Bdry( ii )%Dss   !over-ride kappa !!!!!
    END DO

  ELSE ! not 'C'
      ! set derived types to dummy variables
      DO ii = 1, NPts - 1
        Bdry(ii)%Nodet = -1.
        Bdry(ii)%Noden = -1.
        Bdry(ii)%kappa = 0
        Bdry(ii)%Dxx   = -999.
        Bdry(ii)%Dss   = -999.
      END DO

      Bdry( 1    )%Nodet = [ 1.0, 0.0 ]   ! tangent left-end node
      Bdry( NPts )%Nodet = [ 1.0, 0.0 ]   ! tangent right-end node

  ENDIF ! IF ( CurvilinearFlag( 1:1 ).EQ.'C' )

  RETURN
  END !SUBROUTINE ComputeBdryTangentNormal

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: GetTopSeg
! !INTERFACE:
  SUBROUTINE GetTopSeg( r, myThid )
! !DESCRIPTION:
! Get the Top segment info (index and range interval) for range, r

! !USES:
  USE ihop_mod, only: PRTFile

! !INPUT PARAMETERS:
! r      :: Range at which to get the segment info
! myThid :: my thread ID
  REAL (KIND=_RL90), INTENT( IN ) :: r
  INTEGER,           INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf     :: Informational/error message buffer
! iSegTopT   :: Temporary variable for finding the segment index
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: iSegTopT( 1 )
!EOP

  iSegTopT = MAXLOC( Top( : )%x( 1 ), Top( : )%x( 1 ) < r )

  ! iSegTop MUST LIE IN [ 1, NAtiPts-1 ]
  IF ( iSegTopT( 1 ).GT.0 .AND. iSegTopT( 1 ).LT.NAtiPts ) THEN
    iSegTop = iSegTopT( 1 )
    ! segment limits in range
    rTopSeg = [ Top( iSegTop )%x( 1 ), Top( iSegTop+1 )%x( 1 ) ]
  ELSE
#ifdef IHOP_WRITE_OUT
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
      WRITE(msgBuf,'(A,F10.4)') 'r = ', r
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F10.4)') 'rLeft  = ', Top( 1       )%x( 1 )
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F10.4)') 'rRight = ', Top( NAtiPts )%x( 1 )
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)') 'BDRYMOD GetTopSeg', &
        'Top altimetry undefined above the ray'
      CALL PRINT_ERROR( msgBuf,myThid )
    ENDIF
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R GetTopSeg'
  ENDIF ! IF ( iSegTopT( 1 ).GT.0 .AND. iSegTopT( 1 ).LT.NAtiPts )

  RETURN
  END !SUBROUTINE GetTopSeg

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: GetBotSeg
! !INTERFACE:
  SUBROUTINE GetBotSeg( r, myThid )
! !DESCRIPTION:
! Get the Bottom segment info (index and range interval) for range, r

! !USES:
  USE ihop_mod, only: PRTFile

! !INPUT PARAMETERS:
! r      :: Range at which to get the segment info
! myThid :: my thread ID
  REAL (KIND=_RL90), INTENT( IN ) :: r
  INTEGER,           INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf     :: Informational/error message buffer
! iSegBotT   :: Temporary variable for finding the segment index
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: iSegBotT( 1 )

  iSegBotT = MAXLOC( Bot( : )%x( 1 ), Bot( : )%x( 1 ) < r )
  ! iSegBot MUST LIE IN [ 1, NBtyPts-1 ]
  IF ( iSegBotT( 1 ).GT.0 .AND. iSegBotT( 1 ).LT.NBtyPts ) THEN
    iSegBot = iSegBotT( 1 )
    ! segment limits in range
    rBotSeg = [ Bot( iSegBot )%x( 1 ), Bot( iSegBot + 1 )%x( 1 ) ]
  ELSE
#ifdef IHOP_WRITE_OUT
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
      WRITE(msgBuf,'(A,F10.4)') 'r = ', r
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F10.4)') 'rLeft  = ', Bot( 1       )%x( 1 )
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F10.4)') 'rRight = ', Bot( NBtyPts )%x( 1 )
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)') 'BDRYMOD GetBotSeg', &
        'Bottom bathymetry undefined below the source'
      CALL PRINT_ERROR( msgBuf,myThid )
    ENDIF
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R GetBotSeg'
  ENDIF ! IF ( iSegBotT( 1 ).GT.0 .AND. iSegBotT( 1 ).LT.NBtyPts )

  RETURN
  END !SUBROUTINE GetBotSeg

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteBdry
! !INTERFACE:
  SUBROUTINE WriteBdry( myThid )
! !DESCRIPTION:
! Write Top and Bottom boundary info to PRTFile

! !USES: None

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES: None

  ! Write Top
  CALL WriteTopBot( Top, 'Top', NAtiPts, myThid )
  ! Write Bot
  CALL WriteTopBot( Bot, 'Bot', NBtyPts, myThid )

  RETURN
  END ! SUBROUTINE WriteBdry

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteTopBot
! !INTERFACE:
  SUBROUTINE WriteTopBot( locBdry, BotTop, NPts, myThid )
! !DESCRIPTION:
! Write Top or Bottom boundary info to PRTFile

! !USES:
  USE ihop_mod, only: PRTFile

! !INPUT PARAMETERS:
! locBdry :: Array of boundary points to write
! BotTop  :: Flag indicating bottom or top reflection
! NPts    :: Number of points in the boundary
! myThid  :: my thread ID
  TYPE(BdryPt)                :: locBdry( : )
  CHARACTER*(3), INTENT( IN ) :: BotTop
  INTEGER,       INTENT( IN ) :: NPts
  INTEGER,       INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: locBdry(?) 

! !LOCAL VARIABLES:
! msgBuf     :: Informational/error message buffer
! BdryType   :: Type of boundary interpolation (curvilinear or linear)
! ReadFile   :: Flag indicating whether to read from a file or not
! ii         :: Loop index for writing boundary points
! ranges     :: Array of ranges in km
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
  CHARACTER*(2)        :: BdryType
  CHARACTER*(1)        :: ReadFile
  INTEGER              :: ii
  REAL(KIND=_RL90)     :: ranges
!EOP

  SELECT CASE ( BotTop )
  CASE ( 'Bot' )
    BdryType = btyType
    ReadFile = Bdry%Bot%HS%Opt( 2:2 )
  CASE ( 'Top' )
    BdryType = atiType
    ReadFile = Bdry%Top%HS%Opt( 5:5 )
  CASE DEFAULT
    ! Do nothing
    BdryType = '-'
    STOP 'ABNORMAL END: S/R WriteTopBot'
  END SELECT

  !   Only do I/O in the main thread
  _BARRIER
  _BEGIN_MASTER(myThid)

#ifdef IHOP_WRITE_OUT
  WRITE(msgBuf,'(A)') &
    '___________________________________________________________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  SELECT CASE ( ReadFile )
  CASE ( '~', '*' ) ! read from a file
    IF (BotTop.EQ.'Top') THEN
      WRITE(msgBuf,'(2A)') BotTop, ': Using top-altimetry file'
    ELSE ! 'Bot'
      WRITE(msgBuf,'(2A)') BotTop, ': Using bottom-bathymetry file'
    ENDIF
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    IF (BotTop.EQ.'Top') THEN
      WRITE(msgBuf,'(3A)') '  ATIFile = ', TRIM( IHOP_fileroot ), '.ati'
    ELSE ! 'Bot'
      WRITE(msgBuf,'(3A)') '  BTYFile = ', TRIM( IHOP_fileroot ), '.bty'
    ENDIF
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    SELECT CASE ( BdryType( 1:1 ) )
    CASE ( 'C' )
      WRITE(msgBuf,'(A)') '  Curvilinear Interpolation'
    CASE ( 'L' )
      WRITE(msgBuf,'(A)') '  Piecewise linear interpolation'
    CASE DEFAULT
      STOP 'ABNORMAL END: S/R writeTopBot'
    END SELECT
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    IF (BotTop.EQ.'Top') THEN
      WRITE(msgBuf,'(A,I10)') '  Number of altimetry points = ', NPts-2
    ELSE ! 'Bot'
      WRITE(msgBuf,'(A,I10)') '  Number of bathymetry points = ', NPts-2
    ENDIF
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    IF (BotTop.EQ.'Top') THEN
      WRITE(msgBuf,'(A)') ' Range (km)  Depth (m)'
    ELSE ! 'Bot'
      SELECT CASE ( BdryType( 2:2 ) )
      CASE ( 'S', '' )
        WRITE(msgBuf,'(A)') 'Short format (bathymetry only)'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') ' Range (km)  Depth (m)'

      CASE ( 'L' )
        WRITE(msgBuf,'(A)') 'Long format (bathymetry and geoacoustics)'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(2A)') &
            ' Range (km)  Depth (m)  alphaR (m/s)  betaR  rho (g/cm^3)',&
            'alphaI     betaI'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)')

      CASE DEFAULT
        STOP 'ABNORMAL END: S/R WriteTopBot'
      END SELECT
    ENDIF ! IF (BotTop.EQ.'Top')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    ! Write through all NPts
    boundaryPt: DO ii = 2, NPts - 1
      ! Print ranges in km
      ranges = locBdry( ii )%x(1) / 1000.

      SELECT CASE ( BdryType( 2:2 ) )
      CASE ( 'S', '' )
        IF ( ii.LT.Number_to_Echo .OR. ii.EQ.NPts ) THEN
          WRITE( msgBuf,"(2G11.3)" ) ranges, locBdry( ii )%x(2)
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        ENDIF

      CASE ( 'L' )
        IF ( ii.LT.Number_to_Echo .OR. ii.EQ.NPts ) THEN
          WRITE( msgBuf,'(2F10.2,3X,2F10.2,3X,F6.2,3X,2F10.4)' ) &
            ranges, locBdry(ii)%x(2), &
            locBdry(ii)%HS%alphaR, locBdry(ii)%HS%betaR, &
            locBdry(ii)%HS%rho, &
            locBdry(ii)%HS%alphaI, locBdry(ii)%HS%betaI
            CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        ENDIF

      CASE DEFAULT
        STOP 'ABNORMAL END: S/R writeTopBot'
      END SELECT ! BdryType( 2:2 )

    ENDDO boundaryPt

  CASE DEFAULT ! not reading from file
    IF (BotTop.EQ.'Top') THEN
      WRITE(msgBuf,'(3A)') BotTop, ': No ATIFile. Assume flat top.',&
        ' Set +/- infty.'
    ELSE ! 'Bot'
      WRITE(msgBuf,'(3A)') BotTop, ': No BTYFile. Assume flat bottom.',&
        ' Set +/- infty.'
    END IF
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  END SELECT ! ReadFile
# endif /* IHOP_WRITE_OUT */

  !   Only do I/O in the main thread
  _END_MASTER(myThid)

  RETURN
  END !SUBROUTINE WriteTopBot

END !MODULE bdry_mod
