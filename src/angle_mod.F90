#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE angle_mod
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

  USE ihop_mod,     only: PRTFile
  USE subTab_mod,   only: SubTab
  USE srPos_mod,    only: Pos, Number_to_Echo
  USE sort_mod,     only: Sort

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

  public ReadRayElevationAngles, ReadRayBearingAngles, Angles, ialpha

!=======================================================================

  INTEGER          :: ialpha, ibeta
  INTEGER, PRIVATE :: iAllocStat
  REAL (KIND=_RL90), PRIVATE, PARAMETER :: c0 = 1500.0

  TYPE AnglesStructure
     INTEGER                        :: Nalpha = 0, Nbeta = 1, &
                                       iSingle_alpha = 0, iSingle_beta = 0
     REAL (KIND=_RL90)              :: Dalpha, Dbeta
     REAL (KIND=_RL90), ALLOCATABLE :: alpha( : ), beta( : )
  END TYPE AnglesStructure

  Type( AnglesStructure ) :: Angles

CONTAINS
  SUBROUTINE ReadRayElevationAngles( freq, Depth, TopOpt, RunType, myThid )

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    REAL (KIND=_RL90),  INTENT( IN  ) :: freq, Depth
    CHARACTER (LEN= 6), INTENT( IN  ) :: TopOpt, RunType
    REAL (KIND=_RL90)   :: d_theta_recommended

    IF ( TopOpt( 6 : 6 ) == 'I' ) THEN ! option to trace a single beam
       !READ( ENVFile, * ) Angles%Nalpha, Angles%iSingle_alpha 
    ELSE
       Angles%Nalpha = IHOP_nalpha
    END IF

    IF ( Angles%Nalpha == 0 ) THEN   ! automatically estimate Nalpha to use
       IF ( RunType( 1 : 1 ) == 'R' ) THEN
          ! For a ray trace plot, we don't want too many rays ...
          Angles%Nalpha = 50         
       ELSE
          ! you're letting ME choose? OK: ideas based on an isospeed ocean
          ! limit based on phase of adjacent beams at maximum range
          Angles%Nalpha = MAX( INT( 0.3 * Pos%Rr( Pos%NRr ) * freq / c0 ), 300 )

          ! limit based on having beams that are thin with respect to the water 
          ! depth assumes also a full 360 degree angular spread of rays should 
          ! check which Depth is used here, in case where there is a variable 
          ! bathymetry
          d_theta_recommended = ATAN( Depth / ( 10.0 * Pos%Rr( Pos%NRr ) ) )
          Angles%Nalpha = MAX( INT( PI / d_theta_recommended ), Angles%Nalpha )
       END IF
    END IF

    ALLOCATE( Angles%alpha( MAX( 3, Angles%Nalpha ) ), STAT = iAllocStat )
    IF ( iAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadRayElevationAngles:', &
        'Insufficient memory to store beam angles'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadRayElevationAngles'
    END IF

    Angles%alpha(1:2) = IHOP_alpha
    IF ( Angles%Nalpha > 2 ) Angles%alpha( 3 ) = -999.9

    CALL SubTab( Angles%alpha, Angles%Nalpha )
    CALL Sort(   Angles%alpha, Angles%Nalpha )

    ! full 360-degree sweep? remove duplicate beam
    IF ( Angles%Nalpha > 1 .AND. ABS( MOD( Angles%alpha( Angles%Nalpha ) - & 
                    Angles%alpha( 1 ), 360.0D0 ) ) < 10.0 * TINY( 1.0D0 ) ) &
                    Angles%Nalpha = Angles%Nalpha - 1

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)')'_________________________________________________', &
                        '__________'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A,I10)') 'Number of beams in elevation   = ', Angles%Nalpha
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    IF ( Angles%iSingle_alpha > 0 ) THEN
        WRITE(msgBuf,'(A,I10)') 'Trace only beam number ', Angles%iSingle_alpha
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF
    WRITE(msgBuf,'(A)') 'Beam take-off angles (degrees)'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    IF ( Angles%Nalpha >= 1 ) THEN
        WRITE(msgBuf,'(10F12.6)') Angles%alpha( 1:MIN(Angles%Nalpha,Number_to_Echo) )
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF
    IF ( Angles%Nalpha > Number_to_Echo ) THEN
        WRITE(msgBuf,'(A,F12.6)') ' ... ', Angles%alpha( Angles%Nalpha )
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF
#endif /* IHOP_WRITE_OUT */

    IF ( Angles%Nalpha>1 .AND. &
    Angles%alpha(Angles%Nalpha) == Angles%alpha(1) ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadRayElevationAngles:', &
        'First and last beam take-off angle are identical'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadRayElevationAngles'
    END IF

    IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
        IF ( Angles%iSingle_alpha<1 .OR. &
           Angles%iSingle_alpha>Angles%Nalpha ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadRayElevationAngles:', &
            'Selected beam, iSingl not in [ 1, Angles%Nalpha ]'
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadRayElevationAngles'
        END IF
    END IF
  
  RETURN
  END !SUBROUTINE ReadRayElevationAngles

  !**********************************************************************!

  SUBROUTINE ReadRayBearingAngles( freq, TopOpt, RunType, myThid )

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN ) :: freq
    CHARACTER (LEN= 6), INTENT( IN ) :: TopOpt, RunType

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadBearingElevationAngles:', &
                 '3D rays not supported in ihop'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadBearingElevationAngles'

    IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
       ! option to trace a single beam
       !READ( ENVFile, * ) Angles%Nbeta, Angles%iSingle_beta 
    ELSE
       !READ( ENVFile, * ) Angles%Nbeta
    END IF

    IF ( Angles%Nbeta == 0 ) THEN   ! automatically estimate Nbeta to use
       IF ( RunType( 1 : 1 ) == 'R' ) THEN
          ! For a ray trace plot, we don't want too many rays ...
          Angles%Nbeta = 50         
       ELSE
          Angles%Nbeta = MAX( INT( 0.1 * Pos%rr( Pos%NRr ) * freq / c0 ), 300 )
       END IF
    END IF

    ALLOCATE( Angles%beta( MAX( 3, Angles%Nbeta ) ), STAT = iAllocStat )
    IF ( iAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadBearingElevationAngles:', &
                        'Insufficient memory to store beam angles'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadBearingElevationAngles'
    END IF

    IF ( Angles%Nbeta > 2 ) Angles%beta( 3 ) = -999.9
    !READ( ENVFile, * ) Angles%beta

    CALL SubTab( Angles%beta, Angles%Nbeta )
    CALL Sort(   Angles%beta, Angles%Nbeta )

    ! full 360-degree sweep? remove duplicate beam
    IF ( Angles%Nbeta > 1 .AND. ABS( MOD( Angles%beta( Angles%Nbeta ) &
                    - Angles%beta( 1 ), 360.0D0 ) ) < 10.0 * TINY( 1.0D0 ) ) &
         Angles%Nbeta = Angles%Nbeta - 1

    ! Nx2D CASE: beams must lie on rcvr radials--- replace beta with theta
    IF ( RunType( 6 : 6 ) == '2' .AND. RunType( 1 : 1 ) /= 'R' ) THEN
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)')
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(A)') 'Replacing beam take-off angles, beta, with ', &
                           'receiver bearing lines, theta'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
       DEALLOCATE( Angles%beta )

       Angles%Nbeta = Pos%Ntheta
       ALLOCATE( Angles%beta( MAX( 3, Angles%Nbeta ) ), STAT = iAllocStat )
        IF ( iAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadBearingElevationAngles:', &
                            'Insufficient memory to store beam angles'
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadBearingElevationAngles'
        END IF

       ! Nbeta should = Ntheta
       Angles%beta( 1 : Angles%Nbeta ) = Pos%theta( 1 : Pos%Ntheta )  
    END IF

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A,I10)') 'Number of beams in bearing   = ', Angles%Nbeta
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    IF ( Angles%iSingle_beta > 0 ) THEN
        WRITE(msgBuf,'(A,I10)') 'Trace only beam number ', Angles%iSingle_beta
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF
    WRITE(msgBuf,'(A)') 'Beam take-off angles (degrees)'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    IF ( Angles%Nbeta >= 1 ) THEN
        WRITE(msgBuf,'(5G14.6)') Angles%beta( 1:MIN(Angles%Nbeta,Number_to_Echo) )
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF
    IF ( Angles%Nbeta > Number_to_Echo ) THEN
        WRITE(msgBuf,'(G14.6)') ' ... ', Angles%beta( Angles%Nbeta )
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF
#endif /* IHOP_WRITE_OUT */

    IF ( Angles%Nbeta>1 .AND. &
        Angles%beta( Angles%Nbeta )==Angles%beta(1) ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadBearingElevationAngles:', &
         'First and last beam take-off angle are identical'
        CALL PRINT_ERROR( msgBuf,myThid )

#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadBearingElevationAngles'
    END IF

    IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
        IF ( Angles%iSingle_beta < 1 .OR. &
            Angles%iSingle_beta > Angles%Nbeta ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'ANGLEMOD ReadBearingElevationAngles:', &
            'Selected beam, iSingl not in [ 1, Angles%Nbeta ]'
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadBearingElevationAngles'
        END IF
    END IF
    Angles%beta  = deg2rad * Angles%beta   ! convert to radians

    Angles%Dbeta = 0.0
    IF ( Angles%Nbeta /= 1 ) Angles%Dbeta = ( Angles%beta( Angles%NBeta ) - &
        Angles%beta( 1 ) ) / ( Angles%Nbeta - 1 )

  RETURN
  END !SUBROUTINE ReadRayBearingAngles

END !MODULE angle_mod
