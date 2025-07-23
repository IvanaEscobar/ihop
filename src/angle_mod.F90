#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: angle_mod
MODULE angle_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   Contains the subroutines to read the ray elevation angles

! !USES:
  USE subTab_mod,   only: SubTab
  USE srPos_mod,    only: Pos
  USE sort_mod,     only: Sort
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
  PUBLIC ReadRayElevationAngles, Angles, ialpha
#ifdef IHOP_THREED
  PUBLIC ReadRayBearingAngles
#endif /* IHOP_THREED */
!=======================================================================

! == Module variables ==
  INTEGER, PARAMETER :: Number_to_Echo = 10
  INTEGER            :: ialpha
#ifdef IHOP_THREED
  INTEGER            :: ibeta
#endif /* IHOP_THREED */
  INTEGER, PRIVATE   :: iAllocStat
  REAL (KIND=_RL90), PRIVATE, PARAMETER :: c0 = 1500.0

! == Derived types ==
  TYPE AnglesStructure
    INTEGER                        :: Nalpha = 0, iSingle_alpha = 0
    REAL (KIND=_RL90)              :: Dalpha
    REAL (KIND=_RL90), ALLOCATABLE :: aDeg( : )
    REAL (KIND=_RL90), ALLOCATABLE :: aRad( : )
#ifdef IHOP_THREED
    INTEGER                        :: Nbeta = 1, iSingle_beta = 0
    REAL (KIND=_RL90)              :: Dbeta
    REAL (KIND=_RL90), ALLOCATABLE :: beta( : )
#endif /* IHOP_THREED */
  END TYPE AnglesStructure

  TYPE( AnglesStructure ) :: Angles
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R ReadRayElevationAngles
! S/R ReadRayBearingAngles
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ReadRayElevationAngles
! !INTERFACE:
  SUBROUTINE ReadRayElevationAngles( Depth, TopOpt, RunType, myThid )
! !DESCRIPTION:
!   Read ray elevation angles from the environment file.

! !USES: None

! !INPUT PARAMETERS:
! Depth   :: Water depth at the receiver
! TopOpt  :: Options for ray trace
! RunType :: Type of iHOP run
! myThid  :: my thread ID
  REAL (KIND=_RL90), INTENT( IN ) :: Depth
  CHARACTER*(6),     INTENT( IN ) :: TopOpt, RunType
  INTEGER,           INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! d_theta_recommended :: Recommended angular spacing
! msgBuf              :: Informational/error message buffer
  REAL (KIND=_RL90)   :: d_theta_recommended
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
!EOP

  IF ( TopOpt( 6:6 ).EQ.'I' ) THEN ! option to trace a single beam
    Angles%Nalpha = 0
  ELSE
    Angles%Nalpha = IHOP_nalpha
  ENDIF

  IF ( Angles%Nalpha.EQ.0 ) THEN   ! automatically estimate Nalpha to use
    IF ( RunType( 1:1 ).EQ.'R' ) THEN
      ! For a ray trace plot, we don't want too many rays ...
      Angles%Nalpha = 50
    ELSE
      ! Letting ME choose? OK: ideas based on an isospeed ocean
      ! limit based on phase of adjacent beams at maximum range
      Angles%Nalpha = MAX( INT( 0.3*Pos%RR( Pos%nRR )*IHOP_freq/c0 ), 300 )

      ! limit based on having beams that are thin with respect to the water
      ! depth assumes also a full 360 degree angular spread of rays should
      ! check which Depth is used here, in case where there is a variable
      ! bathymetry
      d_theta_recommended = ATAN( Depth / ( 10.0*Pos%RR( Pos%nRR ) ) )
      Angles%Nalpha = MAX( INT( PI / d_theta_recommended ), Angles%Nalpha )
    ENDIF
  ENDIF ! IF ( Angles%Nalpha.EQ.0 )

  ALLOCATE( Angles%arad( MAX( 3, Angles%Nalpha ) ), &
            Angles%adeg( MAX( 3, Angles%Nalpha ) ), STAT=iAllocStat )
  IF ( iAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadRayElevationAngles:', &
      'Insufficient memory to store beam angles'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadRayElevationAngles'
  ENDIF ! IF ( iAllocStat.NE.0 )

  ! init Angles%adeg,arad
  Angles%adeg = 0.0
  Angles%arad = 0.0
  Angles%adeg(1:2) = IHOP_alpha
  IF ( Angles%Nalpha.GT.2 ) Angles%adeg(3) = -999.9

  ! set intermediate values of alpha
  CALL SubTab( Angles%adeg, Angles%Nalpha )
  CALL Sort(   Angles%adeg, Angles%Nalpha )

  ! full 360-degree sweep? remove duplicate beam
  IF ( Angles%Nalpha.GT.1 &
     .AND. ABS( MOD( Angles%adeg(Angles%Nalpha)-Angles%adeg(1), 360.0 ) ) &
     .LT.10.0*TINY(1.0D0) ) &
    Angles%Nalpha = Angles%Nalpha - 1

  IF ( Angles%Nalpha.GT.1 &
     .AND. Angles%adeg(Angles%Nalpha).EQ.Angles%adeg(1) ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadRayElevationAngles:', &
      'First and last beam take-off angle are identical'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadRayElevationAngles'
  ENDIF ! IF ( Angles%Nalpha>1 &

  IF ( TopOpt( 6:6 ).EQ.'I' ) THEN
    IF ( Angles%iSingle_alpha.LT.1 &
       .OR. Angles%iSingle_alpha.GT.Angles%Nalpha ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadRayElevationAngles:', &
      'Selected beam, iSingl not in [ 1, Angles%Nalpha ]'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadRayElevationAngles'
    ENDIF
  ENDIF ! IF ( TopOpt( 6:6 ).EQ.'I' )

  ! Angles in radians
  Angles%arad = Angles%adeg * deg2rad

  ! set uniform angle spacing
  IF ( Angles%Nalpha.GT.1 ) THEN
    Angles%Dalpha = ( Angles%arad( Angles%Nalpha )-Angles%arad( 1 ) ) &
                    / ( Angles%Nalpha - 1 )  ! angular spacing between beams
  ELSE
    Angles%Dalpha = 0.0
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadRayElevationAngles: ', &
      'Required: Nalpha>1, else add iSingle_alpha'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadRayElevationAngles'
  ENDIF ! IF ( Angles%Nalpha.GT.1 )

  RETURN
  END !SUBROUTINE ReadRayElevationAngles

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
#ifdef IHOP_THREED
!BOP
! !ROUTINE: ReadRayBearingAngles
! !INTERFACE:
  SUBROUTINE ReadRayBearingAngles( TopOpt, RunType, myThid )
! !DESCRIPTION:
! Read ray bearing angles from the environment file.

! !USES:
  USE ihop_mod,     only: PRTFile

! !INPUT PARAMETERS:
! TopOpt  :: Options for ray trace
! RunType :: Type of iHOP run
! myThid  :: my thread ID
  CHARACTER*(6), INTENT( IN ) :: TopOpt, RunType
  INTEGER,       INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf              :: Informational/error message buffer
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
!EOP

#ifdef IHOP_WRITE_OUT
  WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadBearingElevationAngles:', &
    '3D rays not supported in ihop'
  CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
  STOP 'ABNORMAL END: S/R ReadBearingElevationAngles'

  IF ( TopOpt( 6:6 ).EQ.'I' ) THEN
    ! option to trace a single beam
    !READ( ENVFile, * ) Angles%Nbeta, Angles%iSingle_beta
  ELSE
    !READ( ENVFile, * ) Angles%Nbeta
  ENDIF

  IF ( Angles%Nbeta.EQ.0 ) THEN   ! automatically estimate Nbeta to use
    IF ( RunType( 1:1 ).EQ.'R' ) THEN
      ! For a ray trace plot, we don't want too many rays ...
      Angles%Nbeta = 50
    ELSE
      Angles%Nbeta = MAX( INT( 0.1*Pos%RR( Pos%nRR )*IHOP_freq / c0 ), 300 )
    ENDIF
  ENDIF ! IF ( Angles%Nbeta.EQ.0 )

  ALLOCATE( Angles%beta( MAX( 3, Angles%Nbeta ) ), STAT=iAllocStat )
  IF ( iAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadBearingElevationAngles:', &
      'Insufficient memory to store beam angles'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadBearingElevationAngles'
  ENDIF ! IF ( iAllocStat.NE.0 )

  IF ( Angles%Nbeta.GT.2 ) Angles%beta( 3 ) = -999.9
  !READ( ENVFile, * ) Angles%beta

  CALL SubTab( Angles%beta, Angles%Nbeta )
  CALL Sort(   Angles%beta, Angles%Nbeta )

  ! full 360-degree sweep? remove duplicate beam
  IF ( Angles%Nbeta.GT.1 &
     .AND. ABS( MOD( Angles%beta( Angles%Nbeta )-Angles%beta( 1 ), 360.0D0 ) ).LT.10.0 * TINY( 1.0D0 ) ) &
    Angles%Nbeta = Angles%Nbeta - 1

  ! Nx2D CASE: beams must lie on rcvr radials--- replace beta with theta
  IF ( RunType( 6:6 ).EQ.'2' .AND. RunType( 1:1 ).NE.'R' ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)')
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) &
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)') 'Replacing beam take-off angles, beta, with ', &
      'receiver bearing lines, theta'
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) &
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    DEALLOCATE( Angles%beta )

    Angles%Nbeta = Pos%nTheta
    ALLOCATE( Angles%beta( MAX( 3, Angles%Nbeta ) ), STAT=iAllocStat )
    IF ( iAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadBearingElevationAngles:', &
        'Insufficient memory to store beam angles'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadBearingElevationAngles'
    ENDIF ! IF ( iAllocStat.NE.0 )

    ! Nbeta should = Ntheta
    Angles%beta( 1 : Angles%Nbeta ) = Pos%theta( 1 : Pos%nTheta )
  ENDIF ! IF ( RunType( 6:6 ).EQ.'2' .AND. RunType( 1:1 ).NE.'R' )

#ifdef IHOP_WRITE_OUT
  ! In adjoint mode we do not write output besides on the first run
  IF (IHOP_dumpfreq.GE.0) THEN
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A,I10)') 'Number of beams in bearing   = ', Angles%Nbeta
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    IF ( Angles%iSingle_beta.GT.0 ) THEN
      WRITE(msgBuf,'(A,I10)') 'Trace only beam number ', &
        Angles%iSingle_beta
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    ENDIF 

    WRITE(msgBuf,'(A)') 'Beam take-off angles (degrees)'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    IF ( Angles%Nbeta.GE.1 ) THEN
      WRITE(msgBuf,'(5G14.6)') &
        Angles%beta( 1:MIN(Angles%Nbeta,Number_to_Echo) )
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    ENDIF

    IF ( Angles%Nbeta.GT.Number_to_Echo ) THEN
      WRITE(msgBuf,'(A,G14.6)') ' ... ', Angles%beta( Angles%Nbeta )
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    ENDIF

  ENDIF ! IF (IHOP_dumpfreq.GE.0)
#endif /* IHOP_WRITE_OUT */

  IF ( Angles%Nbeta.GT.1 &
     .AND. Angles%beta( Angles%Nbeta ).EQ.Angles%beta(1) ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadBearingElevationAngles:', &
      'First and last beam take-off angle are identical'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadBearingElevationAngles'
  ENDIF ! IF ( Angles%Nbeta.GT.1 &

  IF ( TopOpt( 6:6 ).EQ.'I' ) THEN
    IF ( Angles%iSingle_beta.LT.1 &
       .OR. Angles%iSingle_beta.GT.Angles%Nbeta ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadBearingElevationAngles:', &
        'Selected beam, iSingle not in [ 1, Angles%Nbeta ]'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadBearingElevationAngles'
    ENDIF
  ENDIF ! IF ( TopOpt( 6:6 ).EQ.'I' )

  ! Angles in radians
  Angles%beta = Angles%beta * deg2rad

  Angles%Dbeta = 0.0
  IF ( Angles%Nbeta.NE.1 ) Angles%Dbeta = ( Angles%beta( Angles%NBeta ) - &
                              Angles%beta( 1 ) ) / ( Angles%Nbeta - 1 )

  RETURN
  END !SUBROUTINE ReadRayBearingAngles
#endif /* IHOP_THREED */

END !MODULE angle_mod