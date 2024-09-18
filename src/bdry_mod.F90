#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE bdry_mod
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

  ! Loads altimetry (top bdry) and bathymetry (bottom bdry) data
  ! IEsco22: want to rname this for clarity to readTopBot since it takes user
  ! input and populates various structures used in the code

  IMPLICIT NONE
! == Global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

   PRIVATE

! public interfaces
!=======================================================================

   public   initATI, initBTY, GetTopSeg, GetBotSeg, Bot, Top, &
            IsegTop, IsegBot, rTopSeg, rBotSeg, &
            atiType, btyType, HSInfo, Bdry

!=======================================================================

   INTEGER, PARAMETER :: Number_to_Echo = 21
   INTEGER, PROTECTED :: NBtyPts = 2, NAtiPts = 2
   INTEGER            :: IsegTop, IsegBot ! indices to current active segment

   ! range intervals defining the current active segment
   REAL (KIND=_RL90) :: rTopSeg( 2 ), rBotSeg( 2 )
   CHARACTER*(2)     :: atiType = 'LS', btyType = 'LS'

   ! ***Halfspace properties***
   TYPE HSInfo
      ! compressional and shear wave speeds/attenuations in user units
      REAL     (KIND=_RL90)   :: alphaR, alphaI, betaR, betaI
      REAL     (KIND=_RL90)   :: rho, Depth  ! density, depth
      COMPLEX  (KIND=_RL90)   :: cP, cS      ! P-wave, S-wave speeds
      CHARACTER(LEN=1)        :: BC          ! Boundary condition type
      CHARACTER(LEN=6)        :: Opt
   END TYPE

   ! ***Boundary properties***
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

CONTAINS
  SUBROUTINE initATI( TopATI, DepthT, myThid )
    USE ihop_mod,        only: ATIFile
    USE monotonic_mod,   only: monotonic
    USE atten_mod,       only: CRCI
    ! Reads in the top altimetry
  ! IESCO24
  ! fT = 1000 ONLY for acousto-elastic halfspaces, I will have to pass this
  ! parameter in a different way after ssp_mod is split btwn fixed and varia
  !USE initenvihop, only: fT

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    CHARACTER (LEN= 1), INTENT( IN ) :: TopATI
    INTEGER :: ii
    INTEGER :: IOStat, IAllocStat
    REAL (KIND=_RL90),  INTENT( IN ) :: DepthT
    REAL (KIND=_RL90),  ALLOCATABLE  :: x(:)
    REAL (KIND=_RL90)   :: bPower, fT

    ! IESCO24 fT init
    bPower = 1.0
    fT     = 1000.0
    IF (ALLOCATED(Top)) DEALLOCATE(Top)


    ! Either read from an ATIFile or set +/- infty
    SELECT CASE ( TopATI )
      CASE ( '~', '*' )
        OPEN( UNIT = ATIFile,   FILE = TRIM( IHOP_fileroot ) // '.ati', &
              STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
        IF ( IOSTAT /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
          WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
                               'Unable to open altimetry file'
          ! In adjoint mode we do not write output besides on the first run
          IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
          STOP 'ABNORMAL END: S/R initATI'
        END IF

        READ(  ATIFile, * ) atiType

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
        END SELECT


        READ(  ATIFile, * ) NAtiPts
        ! Extending altimetry to infinity to the left and right
        NAtiPts = NAtiPts + 2

        ALLOCATE( Top( NAtiPts ), Stat = IAllocStat )
        IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
          WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
              'Insufficient memory for altimetry data: reduce # ati points'
          ! In adjoint mode we do not write output besides on the first run
          IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
          STOP 'ABNORMAL END: S/R initATI'
        END IF

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
          END SELECT

          IF ( Top( ii )%x( 2 ) < DepthT ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
                'Altimetry rises above highest point in the sound speed profile'
            ! In adjoint mode we do not write output besides on the first run
            IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R initATI'
          END IF
        END DO atiPt

        CLOSE( ATIFile )

        Top(:)%x(1) = 1000.0 * Top(:)%x(1)   ! Convert ranges in km to m

      CASE DEFAULT   ! no altimetry given, use SSP depth for flat top
        NAtiPts = 2
        ALLOCATE( Top( NAtiPts ), Stat = IAllocStat )
        IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
          WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
                               'Insufficient memory for altimetry data'
          ! In adjoint mode we do not write output besides on the first run
          IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
          STOP 'ABNORMAL END: S/R initATI'
        END IF

        DO ii=1, NAtiPts
          Top( ii )%x = [ -sqrt( huge( Top( 1 )%x( 1 ) ) ) / 1.0d5, DepthT ]
          ! init to default variables
          Top( ii )%HS%alphaR = -1.
          Top( ii )%HS%alphaI = -1.
          Top( ii )%HS%betaR = -1.
          Top( ii )%HS%betaI = -1.
          Top( ii )%HS%rho = -1.
        END DO
        Top(NAtiPts)%x(1) = -1 * Top(1)%x(1)

    END SELECT

    ! Set Top(1:NAtiPts)%t,n,Len,Nodet,Noden,kappa,Dx,Dxx,Dss
    CALL ComputeBdryTangentNormal( Top, 'Top', NAtiPts )

    ! Check if monotonic
    IF ( ALLOCATED(x) ) DEALLOCATE(x)
    ALLOCATE( x(NAtiPts) )
    x = Top%x(1)
    IF ( .NOT. monotonic( x, NAtiPts ) ) THEN
#ifdef IHOP_WRITE_OUT
     WRITE(msgBuf,'(2A)') 'BDRYMOD initATI', &
        'Altimetry ranges are not monotonically increasing'
     ! In adjoint mode we do not write output besides on the first run
     IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
     STOP 'ABNORMAL END: S/R initATI'
    END IF

    ! Set Top%HS%cP,cS
    DO ii = 1, NAtiPts
       ! compressional wave speed
       Top( ii )%HS%cP = -999.
       ! shear wave speed
       Top( ii )%HS%cS = -999.
    END DO
    ! convert range-dependent geoacoustic parameters from user to program units
    ! W is dB/wavelength
    IF ( atiType( 2:2 ) == 'L' ) THEN
      DO ii = 1, NAtiPts
        ! compressional wave speed
        Top( ii )%HS%cP = CRCI( 1D20, Top( ii )%HS%alphaR, &
                                Top( ii )%HS%alphaI, 'W ', bPower, fT, myThid )
        ! shear wave speed
        Top( ii )%HS%cS = CRCI( 1D20, Top( ii )%HS%betaR,  &
                                Top( ii )%HS%betaI, 'W ', bPower, fT,  myThid )
      END DO
    END IF

    ! Set defaults for rest of Top derived type
    Top%HS%Depth = -999.
    Top%HS%BC  = ''
    Top%HS%Opt = ''

    ! Write to PRTFile
    CALL WriteTopBot( Top, 'Top', NAtiPts, myThid ) 

  RETURN
  END !SUBROUTINE initATI

! **********************************************************************!

  SUBROUTINE initBTY( BotBTY, DepthB, myThid )
    USE ihop_mod, only: BTYFile
    USE monotonic_mod,   only: monotonic
    USE atten_mod,       only: CRCI

    ! Reads in the bottom bathymetry
  ! IESCO24
  ! fT = 1000 ONLY for acousto-elastic halfspaces, I will have to pass this
  ! parameter in a different way after ssp_mod is split btwn fixed and varia
  !USE initenvihop, only: fT

   !     == Routine Arguments ==
   !     myThid :: Thread number. Unused by IESCO
   !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

   !     == Local Variables ==
    LOGICAL :: firstnonzero
    CHARACTER (LEN= 1), INTENT( IN ) :: BotBTY
    INTEGER :: i,j,bi,bj,ii
    INTEGER :: IOStat, IAllocStat
    REAL (KIND=_RL90),  INTENT( IN ) :: DepthB
    REAL (KIND=_RL90) :: gcmbathy(sNx,sNy), gcmmin, gcmmax
    REAL (KIND=_RL90), ALLOCATABLE :: x(:)
    REAL (KIND=_RL90)   :: bPower, fT

    ! IESCO24 fT init
    bPower = 1.0
    fT     = 1000.0
    IF (ALLOCATED(Bot)) DEALLOCATE(Bot)


    ! Either read from an BTYFile or set +/- infty
    SELECT CASE ( BotBTY )
      CASE ( '~', '*' )
        OPEN( UNIT = BTYFile, FILE = TRIM( IHOP_fileroot ) // '.bty', &
              STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
        IF ( IOSTAT /= 0 ) THEN
# ifdef IHOP_WRITE_OUT
          WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
                               'Unable to open bathymetry file'
          ! In adjoint mode we do not write output besides on the first run
          IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
          STOP 'ABNORMAL END: S/R initBTY'
        END IF

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
        END SELECT


        READ( BTYFile, * ) NBtyPts
        ! Extend bathymetry to infinity on both sides
        NBtyPts = NBtyPts + 2

        ALLOCATE( Bot( NBtyPts ), Stat = IAllocStat )
        IF ( IAllocStat /= 0 ) THEN
# ifdef IHOP_WRITE_OUT
          WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
             'Insufficient memory for bathymetry data: reduce # bty points'
          ! In adjoint mode we do not write output besides on the first run
          IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
          STOP 'ABNORMAL END: S/R initBTY'
        END IF


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
        END SELECT

        ! R_low check of gcm depths, change to positive values
        firstnonzero=.true.
        DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
              gcmbathy(1:sNx,1:sNy) = rkSign*R_low(1:sNx,1:sNy,bi,bj)
              DO i=1,sNx
                 DO j=1,sNy
                    IF ( gcmbathy(i,j) .ne. 0.0 ) THEN
                       IF (firstnonzero) THEN
                          gcmmin = gcmbathy(i,j)
                          gcmmax = gcmbathy(i,j)
                          firstnonzero=.false.
                       ELSE
                          IF ( gcmbathy(i,j) < gcmmin ) gcmmin = gcmbathy(i,j)
                          IF ( gcmbathy(i,j) > gcmmax ) gcmmax = gcmbathy(i,j)
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
              IF (Bot(ii)%x(2) .lt. gcmmin .or. Bot(ii)%x(2) .gt. gcmmax) THEN
                WRITE(msgBuf,'(2A,F10.2,A,F10.2)') &
                 '** Warning ** BDRYMOD initBTY: ', &
                 'ihop and gcm bathymetry vary, ihop:', Bot(ii)%x(2), &
                 'gcm:', gcmmax
                ! In adjoint mode we do not write output besides on first run
                IF (IHOP_dumpfreq.GE.0) &
                  CALL PRINT_MESSAGE(msgBuf, errorMessageUnit, &
                                     SQUEEZE_RIGHT, myThid)
              END IF
# endif /* IHOP_WRITE_OUT */

            CASE ( 'L' )       ! long format
              READ(  BTYFile, * ) Bot( ii )%x, Bot( ii )%HS%alphaR, &
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
          END SELECT

          IF ( Bot( ii )%x( 2 ) > DepthB ) THEN
# ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
               'Bathymetry drops below lowest point in the sound speed profile'
            ! In adjoint mode we do not write output besides on first run
            IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R initBTY'
          END IF
        END DO btyPt

        CLOSE( BTYFile )

        Bot(:)%x(1) = 1000.0 * Bot(:)%x(1)   ! Convert ranges in km to m

    CASE DEFAULT   ! no bathymetry given, use SSP depth for flat bottom
      NBtyPts = 2
      ALLOCATE( Bot( NBtyPts ), Stat = IAllocStat )
      IF ( IAllocStat /= 0 ) THEN
#  ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
                             'Insufficient memory for bathymetry data'
        ! In adjoint mode we do not write output besides on first run
        IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
#  endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R initBTY'
      END IF

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

    END SELECT

    ! Set Bot(1:NBtyPts)%t,n,Len,Nodet,Noden,kappa,Dx,Dxx,Dss
    CALL ComputeBdryTangentNormal( Bot, 'Bot', NBtyPts )

    IF ( ALLOCATED(x) ) DEALLOCATE(x)
    ALLOCATE( x(NBtyPts) )
    x = Bot%x(1)
    IF ( .NOT. monotonic( x, NBtyPts ) ) THEN
# ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'BDRYMOD initBTY: ', &
          'Bathymetry ranges are not monotonically increasing'
      ! In adjoint mode we do not write output besides on first run
      IF (IHOP_dumpfreq.GE.0) CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R initBTY'
    END IF


    ! Initiate Bot
    DO ii = 1, NBtyPts
      ! compressional wave speed
      Bot( ii )%HS%cP = -999.
      ! shear wave speed
      Bot( ii )%HS%cS = -999.
    END DO
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
     END DO
   END IF

    ! Set defaults for rest of Bot derived type
    Bot%HS%Depth = -999.
    Bot%HS%BC  = ''
    Bot%HS%Opt = ''

    ! Write to PRTFile
    CALL WriteTopBot( Bot, 'Bot', NBtyPts, myThid )

   RETURN
  END !SUBROUTINE initBTY

! **********************************************************************!

  SUBROUTINE ComputeBdryTangentNormal( Bdry, BotTop, NPts )

    ! Does some pre-processing on the boundary points to pre-compute segment
    ! lengths  (%Len),
    ! tangents (%t, %nodet),
    ! normals  (%n, %noden), and
    ! curvatures (%kappa)
    !
    ! The boundary is also extended with a constant depth to infinity to cover
    ! cases where the ray leaves the domain defined by the user

    INTEGER, INTENT(IN)              :: NPts
    REAL (KIND=_RL90), ALLOCATABLE   :: phi( : )
    REAL (KIND=_RL90)                :: sss
    TYPE(BdryPt)                     :: Bdry( : )
    CHARACTER (LEN=3),  INTENT( IN ) :: BotTop ! Flag indicating bottom or top reflection
    CHARACTER (LEN=2)                :: CurvilinearFlag
    INTEGER :: ii
    INTEGER :: IAllocStat

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
       Bdry( ii )%Dx  = Bdry( ii )%t( 2 ) / Bdry( ii )%t( 1 )   ! 1st derivative

       ! normalize the tangent vector
       Bdry( ii )%Len = NORM2( Bdry( ii )%t )
       Bdry( ii )%t   = Bdry( ii )%t / Bdry( ii )%Len

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
    IF ( CurvilinearFlag( 1:1 ) == 'C' ) THEN
       ! averaging two centered differences is equivalent to forming a single
       ! centered difference of two steps ...
       DO ii = 2, NPts - 1
          sss = Bdry( ii - 1 )%Len / ( Bdry( ii - 1 )%Len + Bdry( ii )%Len )
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
         Bdry(ii)%Dxx = -999.
         Bdry(ii)%Dss = -999.
       END DO
       Bdry( 1    )%Nodet = [ 1.0, 0.0 ]   ! tangent left-end node
       Bdry( NPts )%Nodet = [ 1.0, 0.0 ]   ! tangent right-end node

    END IF

  RETURN
  END !SUBROUTINE ComputeBdryTangentNormal

  ! **********************************************************************!

  SUBROUTINE GetTopSeg( r, myThid )
    USE ihop_mod, only: PRTFile

    ! Get the Top segment info (index and range interval) for range, r

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    INTEGER IsegTopT( 1 )
    REAL (KIND=_RL90), INTENT( IN ) :: r

    IsegTopT = MAXLOC( Top( : )%x( 1 ), Top( : )%x( 1 ) < r )

    ! IsegTop MUST LIE IN [ 1, NAtiPts-1 ]
    IF ( IsegTopT( 1 ) > 0 .AND. IsegTopT( 1 ) < NAtiPts ) THEN
       IsegTop = IsegTopT( 1 )
       ! segment limits in range
       rTopSeg = [ Top( IsegTop )%x( 1 ), Top( IsegTop+1 )%x( 1 ) ]
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
    ENDIF

  RETURN
  END !SUBROUTINE GetTopSeg

  ! **********************************************************************!

  SUBROUTINE GetBotSeg( r, myThid )
    USE ihop_mod, only: PRTFile

    ! Get the Bottom segment info (index and range interval) for range, r
  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    INTEGER IsegBotT( 1 )
    REAL (KIND=_RL90), INTENT( IN ) :: r

    IsegBotT = MAXLOC( Bot( : )%x( 1 ), Bot( : )%x( 1 ) < r )
    ! IsegBot MUST LIE IN [ 1, NBtyPts-1 ]
    IF ( IsegBotT( 1 ) > 0 .AND. IsegBotT( 1 ) < NBtyPts ) THEN
       IsegBot = IsegBotT( 1 )
       ! segment limits in range
       rBotSeg = [ Bot( IsegBot )%x( 1 ), Bot( IsegBot + 1 )%x( 1 ) ]
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
    ENDIF

  RETURN
  END !SUBROUTINE GetBotSeg

  ! **********************************************************************!

  SUBROUTINE WriteTopBot( locBdry, BotTop, NPts, myThid )
    USE ihop_mod, only: PRTFile

  ! Write initATI, initBTY to PRTFile

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    INTEGER, INTENT(IN)              :: NPts
    REAL(KIND=_RL90)                 :: ranges
    TYPE(BdryPt)                     :: locBdry( : )
    CHARACTER (LEN=3),  INTENT( IN ) :: BotTop
    CHARACTER (LEN=2)                :: BdryType
    CHARACTER (LEN=1)                :: ReadFile 
    INTEGER                          :: ii

    SELECT CASE ( BotTop )
      CASE ( 'Bot' )
         BdryType = btyType
         ReadFile = Bdry%Bot%HS%Opt(2:2)
      CASE ( 'Top' )
         BdryType = atiType
         ReadFile = Bdry%Top%HS%Opt(5:5)
      CASE DEFAULT
         ! Do nothing
         BdryType = '-'
         STOP 'ABNORMAL END: S/R WrtieTopBot'
    END SELECT

    !   Only do I/O in the main thread
    _BARRIER
    _BEGIN_MASTER(myThid)

#ifdef IHOP_WRITE_OUT
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
      WRITE(msgBuf,'(2A)') '____________________________________________', &
                           '_______________'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )


      SELECT CASE ( ReadFile )
        CASE ( '~', '*' ) ! read from a file
          IF (BotTop.eq.'Top') THEN
            WRITE(msgBuf,'(2A)') BotTop, ': Using top-altimetry file'
          ELSE ! 'Bot'
            WRITE(msgBuf,'(2A)') BotTop, ': Using bottom-bathymetry file'
          ENDIF
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  
          IF (BotTop.eq.'Top') THEN
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

          IF (BotTop.eq.'Top') THEN
            WRITE(msgBuf,'(A,I10)') '  Number of altimetry points = ', NPts-2
          ELSE ! 'Bot'
            WRITE(msgBuf,'(A,I10)') '  Number of bathymetry points = ', NPts-2
          ENDIF
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(A)')
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

          IF (BotTop.eq.'Top') THEN
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
          ENDIF
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )


          ! Write through all NPts
          boundaryPt: DO ii = 2, NPts - 1
            ! Print ranges in km
            ranges = locBdry( ii )%x(1) / 1000.

            SELECT CASE ( BdryType( 2:2 ) )
              CASE ( 'S', '' )
                IF ( ii < Number_to_Echo .OR. ii == NPts ) THEN

                  WRITE( msgBuf,"(2G11.3)" ) ranges, locBdry( ii )%x(2)
                  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
                END IF

              CASE ( 'L' )
                IF ( ii < Number_to_Echo .OR. ii == NPts ) THEN
                  WRITE( msgBuf,'(2F10.2,3X,2F10.2,3X,F6.2,3X,2F10.4)' ) &
                      ranges, locBdry(ii)%x(2), &
                      locBdry(ii)%HS%alphaR, locBdry(ii)%HS%betaR, &
                      locBdry(ii)%HS%rho, &
                      locBdry(ii)%HS%alphaI, locBdry(ii)%HS%betaI
                   CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
                END IF
              CASE DEFAULT
                STOP 'ABNORMAL END: S/R writeTopBot'
            END SELECT

          END DO boundaryPt

        CASE DEFAULT ! not reading from file
          IF (BotTop.eq.'Top') THEN
            WRITE(msgBuf,'(3A)') BotTop, ': No ATIFile. Assume flat top.',&
                                         ' Set +/- infty.'
          ELSE ! 'Bot'
            WRITE(msgBuf,'(3A)') BotTop, ': No BTYFile. Assume flat bottom.',&
                                         ' Set +/- infty.'
          END IF
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      END SELECT

      WRITE(msgBuf,'(2A)')'_____________________________________________', &
                          '______________'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    END IF ! IHOP_dumpfreq.GE.0
# endif /* IHOP_WRITE_OUT */

    !   Only do I/O in the main thread
    _END_MASTER(myThid)

  RETURN
  END !SUBROUTINE WriteTopBot

END !MODULE bdry_mod
