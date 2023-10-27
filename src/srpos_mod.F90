#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE srpos_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>

  ! Reads in source depths, receiver depths, receiver ranges, and receiver bearings

  USE subTab_mod,       only: SubTab
  USE monotonic_mod,    only: monotonic
  USE sort_mod,         only: Sort
  USE ihop_mod,         only: PRTFile

! ! USES
  IMPLICIT NONE
!  == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

  PRIVATE

! public interfaces
!=======================================================================

    public Pos, Number_to_Echo, Nfreq, freqVec, ReadSxSy, ReadSzRz, &
           ReadRcvrRanges, &
#ifdef IHOP_THREED
           ReadRcvrBearings, &
#endif /* IHOP_THREED */
           ReadFreqVec

!=======================================================================

  INTEGER, PARAMETER    :: Number_to_Echo = 10
  INTEGER, PRIVATE      :: IAllocStat     ! Status after memory allocation
  INTEGER               :: Nfreq = 1      ! number of frequencies
  REAL (KIND=_RL90), ALLOCATABLE  :: freqVec( : )   ! frequency vector for braodband runs

  TYPE Position
     INTEGER              :: NSx = 1, NSy = 1, NSz = 1, & ! # of x,y,z coords
                             NRz = 1, NRr = 1, Ntheta = 1 ! # of z,r,theta coord`s
     REAL                 :: Delta_r, Delta_theta
     INTEGER,           ALLOCATABLE :: iSz( : ), iRz( : )
     REAL (KIND=_RL90), ALLOCATABLE :: Sx( : ), Sy( : ), Sz( : )          ! Source x, y, z coordinates
     REAL (KIND=_RL90), ALLOCATABLE :: Rr( : ), Rz( : ), ws( : ), wr( : ) ! Receiver r, z coordinates and weights for interpolation
     REAL (KIND=_RL90), ALLOCATABLE :: theta( : )                         ! Receiver bearings
  END TYPE Position

  TYPE ( Position ) :: Pos ! structure containing source and receiver positions

CONTAINS
  SUBROUTINE ReadfreqVec( freq0, BroadbandOption, myThid )

    ! Optionally reads a vector of source frequencies for a broadband run
    ! If the broadband option is not selected, then the input freq (a scalar) 
    ! is stored in the frequency vector
    
  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    REAL (KIND=_RL90),  INTENT( IN ) :: freq0   ! Source frequency
    CHARACTER,          INTENT( IN ) :: BroadbandOption*( 1 )
    INTEGER :: ifreq


    ! Broadband run?
    IF ( BroadbandOption == 'B' ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)')'________________________________________________',&
                            '___________'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') 
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') 
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,I10)') 'Number of frequencies =', Nfreq
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
        IF ( Nfreq <= 0 ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadfreqVec: ', &
                                 'Number of frequencies must be positive'
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadfreqVec'
        END IF
    END IF

    IF ( ALLOCATED( freqVec ) ) DEALLOCATE( freqVec )
    ALLOCATE( freqVec( MAX( 3, Nfreq ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadfreqVec: ', &
                             'Number of frequencies too large'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadfreqVec'
    END IF

    IF ( BroadbandOption == 'B' ) THEN
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Frequencies (Hz)'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
       freqVec( 3 ) = -999.9
       !READ(  ENVFile, * ) freqVec( 1 : Nfreq )
       CALL SubTab( freqVec, Nfreq )

#ifdef IHOP_WRITE_OUT
       WRITE( msgBuf, '(5G14.6)' ) ( freqVec( ifreq ), ifreq = 1, &
           MIN( Nfreq, Number_to_Echo ) )
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
       IF ( Nfreq > Number_to_Echo ) &
           WRITE( msgBuf,'(G14.6)' ) ' ... ', freqVec( Nfreq )
           CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    ELSE
       freqVec( 1 ) = freq0
    END IF

  RETURN
  END !SUBROUTINE ReadfreqVec

  !********************************************************************!

  SUBROUTINE ReadSxSy( myThid )

    ! Reads source x-y coordinates
    
  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==

#ifdef IHOP_THREED
    CALL ReadVector( Pos%NSx, Pos%Sx, 'source   x-coordinates, Sx', 'km', &
                    myThid )
    CALL ReadVector( Pos%NSy, Pos%Sy, 'source   y-coordinates, Sy', 'km', &
                    myThid )
#else /* IHOP_THREED */
    ALLOCATE( Pos%Sx( 1 ), Pos%Sy( 1 ), Stat=IAllocStat)
    IF (IAllocStat/=0) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf, *) 'allocation failed', IAllocStat
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSxSy'
    END IF

    Pos%Sx( 1 ) = 0.
    Pos%Sy( 1 ) = 0.
#endif /* IHOP_THREED */
  RETURN
  END !SUBROUTINE ReadSxSy

  !********************************************************************!

  SUBROUTINE ReadSzRz( zMin, zMax, myThid )

    ! Reads source and receiver z-coordinates (depths)
    ! zMin and zMax are limits for those depths; sources and receivers are 
    ! shifted to be within those limits

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    REAL,    INTENT( IN ) :: zMin, zMax

    CALL ReadVector( Pos%NSz, Pos%Sz, 'Source   depths, Sz', 'm', &
                    myThid )
    CALL ReadVector( Pos%NRz, Pos%Rz, 'Receiver depths, Rz', 'm', &
                    myThid )

    IF ( ALLOCATED( Pos%ws ) ) DEALLOCATE( Pos%ws, Pos%iSz )
    ALLOCATE( Pos%ws( Pos%NSz ), Pos%iSz( Pos%NSz ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadSzRz: ', &
                             'Too many sources'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSzRz'
    END IF

    IF ( ALLOCATED( Pos%wr ) ) DEALLOCATE( Pos%wr, Pos%iRz )
    ALLOCATE( Pos%wr( Pos%NRz ), Pos%iRz( Pos%NRz ), Stat = IAllocStat  )
    IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadSzRz: ', &
                             'Too many receivers'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSzRz'
    END IF

    ! *** Check for Sz/Rz in water column ***

#ifdef IHOP_WRITE_OUT
    IF ( ANY( Pos%Sz( 1 : Pos%NSz ) < zMin ) ) THEN
       WHERE ( Pos%Sz < zMin ) Pos%Sz = zMin
       WRITE(msgBuf,'(2A)') 'Warning in ReadSzRz : Source above or too ',&
                           'near the top bdry has been moved down'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF

    IF ( ANY( Pos%Sz( 1 : Pos%NSz ) > zMax ) ) THEN
       WHERE( Pos%Sz > zMax ) Pos%Sz = zMax
       WRITE(msgBuf,'(2A)') 'Warning in ReadSzRz : Source below or too ',&
                           'near the bottom bdry has been moved up'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF

    IF ( ANY( Pos%Rz( 1 : Pos%NRz ) < zMin ) ) THEN
       WHERE( Pos%Rz < zMin ) Pos%Rz = zMin
       WRITE(msgBuf,'(2A)') 'Warning in ReadSzRz : Receiver above or too ',&
           'near the top bdry has been moved down'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF

    IF ( ANY( Pos%Rz( 1 : Pos%NRz ) > zMax ) ) THEN
       WHERE( Pos%Rz > zMax ) Pos%Rz = zMax
       WRITE(msgBuf,'(2A)') 'Warning in ReadSzRz : Receiver below or too ',&
                           'near the bottom bdry has been moved up'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE ReadSzRz

  !********************************************************************!

  SUBROUTINE ReadRcvrRanges( myThid )

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==

    ! IESCO22: assuming receiver positions are equally spaced
    CALL ReadVector( Pos%NRr, Pos%Rr, 'Receiver ranges, Rr', 'km', myThid )

    ! calculate range spacing
    Pos%delta_r = 0.0
    IF ( Pos%NRr /= 1 ) Pos%delta_r = Pos%Rr( Pos%NRr ) - Pos%Rr( Pos%NRr-1 )

    IF ( .NOT. monotonic( Pos%Rr, Pos%NRr ) ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadRcvrRanges: ', &
                             'Receiver ranges are not monotonically increasing'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadRcvrRanges'
    END IF 
 
    RETURN
  END !SUBROUTINE ReadRcvrRanges

  !********************************************************************!
  
#ifdef IHOP_THREED
  SUBROUTINE ReadRcvrBearings( myThid )   ! for 3D bellhop

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==

! IEsco23: NOT SUPPORTED IN ihop
    CALL ReadVector( Pos%Ntheta, Pos%theta, 'receiver bearings, theta', &
        'degrees', myThid )

    ! full 360-degree sweep? remove duplicate angle
    IF ( Pos%Ntheta > 1 ) THEN
       IF ( ABS( MOD( Pos%theta( Pos%Ntheta ) - Pos%theta( 1 ), 360.0 ) ) &
           < 10.0 * TINY( 1.0D0 ) ) &
          Pos%Ntheta = Pos%Ntheta - 1
    END IF

    ! calculate angular spacing
    Pos%Delta_theta = 0.0
    IF ( Pos%Ntheta /= 1 ) Pos%Delta_theta = Pos%theta( Pos%Ntheta ) &
                                           - Pos%theta( Pos%Ntheta - 1 )

    IF ( .NOT. monotonic( Pos%theta, Pos%Ntheta ) ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadRcvrBearings: ', &
                            'Receiver bearings are not monotonically increasing'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadRcvrBearings'
    END IF 
 
    RETURN
  END !SUBROUTINE ReadRcvrBearings
#endif /* IHOP_THREED */
  !********************************************************************!

  SUBROUTINE ReadVector( Nx, x, Description, Units, myThid )

    ! Read a vector x
    ! Description is something like 'receiver ranges'
    ! Units       is something like 'km'

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    INTEGER,                        INTENT( IN ) :: Nx
    REAL (KIND=_RL90), ALLOCATABLE, INTENT( INOUT ) :: x( : )
    CHARACTER,                      INTENT( IN ) :: Description*( * ), &
                                                    Units*( * )
    INTEGER :: ix
   
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)')'__________________________________________________', &
                        '_________'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A,I10)') 'Number of ' // Description // ' = ', Nx
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    IF ( Nx <= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadVector: ', &
                             'Number of ' // Description // 'must be positive'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadVector'
    END IF

    IF ( .NOT. ALLOCATED( x ) ) THEN 
        ALLOCATE( x( MAX( 3, Nx ) ), Stat = IAllocStat )
        IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadVector: ', &
                                'Too many ' // Description
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadVector'
        END IF
    END IF

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') Description // ' (' // Units // ')'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    CALL SubTab( x, Nx )
    CALL Sort(   x, Nx )

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(5G14.6)') ( x( ix ), ix = 1, MIN( Nx, Number_to_Echo ) )
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    IF ( Nx > Number_to_Echo ) THEN
        WRITE(msgBuf,'(G14.6)') ' ... ', x( Nx )
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF

    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    ! Vectors in km should be converted to m for internal use
    IF ( LEN_TRIM( Units ) >= 2 ) THEN
       IF ( Units( 1 : 2 ) == 'km' ) x = 1000.0 * x
    END IF

  RETURN
  END !SUBROUTINE ReadVector

END MODULE srpos_mod
