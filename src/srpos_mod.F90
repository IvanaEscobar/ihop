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
  USE ihop_mod,         only: PRTFile

! ! USES
  IMPLICIT NONE
!  == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

  PRIVATE

! public interfaces
!=======================================================================

    public Pos, Nfreq, freqVec, &
           ReadSxSy,  ReadSzRz,  ReadRcvrRanges,  ReadFreqVec, &
           WriteSxSy, WriteSzRz, WriteRcvrRanges, WriteFreqVec
#ifdef IHOP_THREED
    public ReadRcvrBearings, WriteRcvrBearings
#endif /* IHOP_THREED */

!=======================================================================

  INTEGER, PARAMETER    :: Number_to_Echo = 10
  INTEGER, PRIVATE      :: IAllocStat     ! Status after memory allocation
  INTEGER               :: Nfreq = 1      ! number of frequencies
  REAL (KIND=_RL90), ALLOCATABLE  :: freqVec( : )   ! frequency vector for braodband runs

  TYPE Position
     INTEGER              :: NSx = 1, NSy = 1, NSz = 1, & ! # of x,y,z coords
                             NRz = 1, NRr = 1, nTheta = 1 ! # of z,r,theta coord`s
     REAL                 :: Delta_r, Delta_theta
     INTEGER,           ALLOCATABLE :: iSz( : ), iRz( : )
     REAL (KIND=_RL90), ALLOCATABLE :: Sx( : ), Sy( : ), Sz( : )          ! Source x, y, z coordinates
     REAL (KIND=_RL90), ALLOCATABLE :: Rr( : ), Rz( : ), ws( : ), wr( : ) ! Receiver r, z coordinates and weights for interpolation
     REAL (KIND=_RL90), ALLOCATABLE :: theta( : )                         ! Receiver bearings
  END TYPE Position

  TYPE ( Position ) :: Pos ! structure containing source and receiver positions

CONTAINS
  SUBROUTINE ReadfreqVec( BroadbandOption, myThid )

    ! Optionally reads a vector of source frequencies for a broadband run
    ! If the broadband option is not selected, then the input freq (a scalar)
    ! is stored in the frequency vector
    ! IHOP_freq is source frequency

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    CHARACTER,          INTENT( IN ) :: BroadbandOption*( 1 )
    INTEGER :: ifreq


    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
        ! Broadband run?
        IF ( BroadbandOption == 'B' ) THEN
            IF ( Nfreq <= 0 ) THEN
#ifdef IHOP_WRITE_OUT
                WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadfreqVec: ', &
                                 'Number of frequencies must be positive'
                CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
                STOP 'ABNORMAL END: S/R ReadfreqVec'
            END IF
        END IF
    ENDIF

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

    ! set default values
    freqVec = 0.0

    IF ( BroadbandOption == 'B' ) THEN
       freqVec(3) = -999.9
       !READ(  ENVFile, * ) freqVec( 1 : Nfreq )
       CALL SubTab( freqVec, Nfreq )

    ELSE
       freqVec(1) = IHOP_freq
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
    CALL ReadVector( Pos%nSX, Pos%SX, 'source   x-coordinates, Sx', 'km', &
                    myThid )
    CALL ReadVector( Pos%nSY, Pos%SY, 'source   y-coordinates, Sy', 'km', &
                    myThid )
#else /* IHOP_THREED */
    ALLOCATE( Pos%SX( 1 ), Pos%SY( 1 ), Stat=IAllocStat)
    IF (IAllocStat/=0) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf, *) 'allocation failed', IAllocStat
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSxSy'
    END IF

    Pos%SX( 1 ) = 0.
    Pos%SY( 1 ) = 0.
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
    REAL(KIND=_RL90),    INTENT( IN ) :: zMin, zMax
    INTEGER :: posShift

    CALL ReadVector( Pos%nSZ, Pos%SZ, 'Source   depths, Sz', 'm', &
                    myThid )
    CALL ReadVector( Pos%nRZ, Pos%RZ, 'Receiver depths, Rz', 'm', &
                    myThid )

    ! *** Check for Sz/Rz in water column ***
    posShift = 0
    IF ( ANY( Pos%SZ( 1:Pos%nSZ ) < zMin ) ) THEN
      WHERE ( Pos%SZ < zMin ) Pos%SZ = zMin
      posShift = 1
    END IF

    IF ( ANY( Pos%SZ( 1:Pos%nSZ ) > zMax ) ) THEN
      WHERE( Pos%SZ > zMax ) Pos%SZ = zMax
      posShift = 2
    END IF

    IF ( ANY( Pos%RZ( 1:Pos%nRZ ) < zMin ) ) THEN
      WHERE( Pos%RZ < zMin ) Pos%RZ = zMin
      posShift = 3
    END IF

    IF ( ANY( Pos%RZ( 1:Pos%nRZ ) > zMax ) ) THEN
      WHERE( Pos%RZ > zMax ) Pos%RZ = zMax
      posShift = 4
    END IF

#ifdef IHOP_WRITE_OUT
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
      IF ( posShift.eq.1 ) THEN
        WRITE(msgBuf,'(2A)') 'Warning in WriteSzRz : Source above or too ',&
                            'near the top bdry has been moved down'
        CALL PRINT_ERROR( msgBuf,myThid )
      END IF

      IF ( posShift.eq.2 ) THEN
        WRITE(msgBuf,'(2A)') 'Warning in WriteSzRz : Source below or too ',&
                            'near the bottom bdry has been moved up'
        CALL PRINT_ERROR( msgBuf,myThid )
      END IF

      IF ( posShift.eq.3 ) THEN
        WRITE(msgBuf,'(2A)') 'Warning in WriteSzRz : Receiver above or too ',&
            'near the top bdry has been moved down'
        CALL PRINT_ERROR( msgBuf,myThid )
      END IF

      IF ( posShift.eq.4 ) THEN
        WRITE(msgBuf,'(2A)') 'Warning in WriteSzRz : Receiver below or too ',&
                            'near the bottom bdry has been moved up'
        CALL PRINT_ERROR( msgBuf,myThid )
      END IF
    ENDIF
#endif /* IHOP_WRITE_OUT */

    IF ( ALLOCATED( Pos%ws ) ) DEALLOCATE( Pos%ws, Pos%iSz )
    ALLOCATE( Pos%ws( Pos%nSZ ), Pos%iSz( Pos%nSZ ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(A)') 'SRPOSITIONS ReadSzRz: Too many sources'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSzRz'
    END IF

    IF ( ALLOCATED( Pos%wr ) ) DEALLOCATE( Pos%wr, Pos%iRz )
    ALLOCATE( Pos%wr( Pos%nRZ ), Pos%iRz( Pos%nRZ ), Stat = IAllocStat  )
    IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(A)') 'SRPOSITIONS ReadSzRz: Too many receivers'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSzRz'
    END IF

    ! Set default values
    Pos%ws = 0
    Pos%isz = 0
    Pos%wr = 0
    Pos%irz = 0

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
    CALL ReadVector( Pos%nRR, Pos%RR, 'Receiver ranges, Rr', 'km', myThid )

    ! calculate range spacing
    Pos%delta_r = 0.0
    ! IESCO24: Assumes uniform spacing in ranges btwn receivers
    IF ( Pos%nRR /= 1 ) Pos%delta_r = Pos%RR( Pos%nRR ) - Pos%RR( Pos%nRR-1 )

    IF ( .NOT. monotonic( Pos%RR, Pos%nRR ) ) THEN
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
  SUBROUTINE ReadRcvrBearings( myThid )   ! for 3D models

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==

! IEsco23: NOT SUPPORTED IN ihop
    CALL ReadVector( Pos%nTheta, Pos%theta, 'receiver bearings, theta', &
        'degrees', myThid )

    ! full 360-degree sweep? remove duplicate angle
    IF ( Pos%nTheta > 1 ) THEN
       IF ( ABS( MOD( Pos%theta( Pos%nTheta ) - Pos%theta( 1 ), 360.0 ) ) &
           < 10.0 * TINY( 1.0D0 ) ) &
          Pos%nTheta = Pos%nTheta - 1
    END IF

    ! calculate angular spacing
    Pos%Delta_theta = 0.0
    IF ( Pos%nTheta /= 1 ) Pos%Delta_theta = Pos%theta( Pos%nTheta ) &
                                           - Pos%theta( Pos%nTheta - 1 )

    IF ( .NOT. monotonic( Pos%theta, Pos%nTheta ) ) THEN
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

    USE sort_mod,         only: Sort

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    INTEGER,                        INTENT( IN )    :: Nx
    REAL (KIND=_RL90), ALLOCATABLE, INTENT( INOUT ) :: x( : )
    CHARACTER,                      INTENT( IN )    :: Description*( * ), &
                                                       Units*( * )
    INTEGER :: ix

    IF ( Nx <= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SRPOS_MOD ReadVector: ', &
                             'Number of ' // Description // 'must be positive'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadVector'
    END IF

    IF ( .NOT. ALLOCATED( x ) ) THEN
        ALLOCATE( x( MAX( 3, Nx ) ), Stat = IAllocStat )
        IF ( IAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'SRPOS_MOD ReadVector: ', &
                                'Too many ' // Description
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadVector'
        END IF
    END IF

    CALL SubTab( x, Nx )
    CALL Sort(   x, Nx )

    ! Vectors in km should be converted to m for internal use
    IF ( LEN_TRIM( Units ) >= 2 ) THEN
       IF ( Units( 1:2 ) == 'km' ) x = 1000.0 * x
    END IF

  RETURN
  END !SUBROUTINE ReadVector

! ==============================================================================
! ==============================================================================
  SUBROUTINE writeFreqVec( BroadbandOption, myThid )

    ! Writes a vector of source frequencies for a broadband run

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    CHARACTER*(1), INTENT( IN ) :: BroadbandOption
    INTEGER :: ifreq


#ifdef IHOP_WRITE_OUT
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
      ! Broadband run?
      IF ( BroadbandOption == 'B' ) THEN
        WRITE(msgBuf,'(2A)')'___________________________________________',&
                            '________________'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)')
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)')
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,I10)') 'Number of frequencies =', Nfreq
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

        WRITE(msgBuf,'(A)') 'Frequencies (Hz)'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

        WRITE( msgBuf, '(5G14.6)' ) ( freqVec( ifreq ), ifreq = 1, &
            MIN( Nfreq, Number_to_Echo ) )
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        IF ( Nfreq > Number_to_Echo ) THEN
          WRITE( msgBuf,'(G14.6)' ) ' ... ', freqVec( Nfreq )
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        ENDIF
      ENDIF ! Broadband run
    ENDIF !adjoint run?
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE writeFreqVec

  !********************************************************************!

  SUBROUTINE WriteSxSy( myThid )

    ! Writes source x-y coordinates

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==

#ifdef IHOP_THREED
#ifdef IHOP_WRITE_OUT
    CALL WriteVector( Pos%nSX, Pos%SX, 'source   x-coordinates, Sx', 'km', &
                    myThid )
    CALL WriteVector( Pos%nSY, Pos%SY, 'source   y-coordinates, Sy', 'km', &
                    myThid )
#endif /* IHOP_WRITE_OUT */
#endif /* IHOP_THREED */

  RETURN
  END !SUBROUTINE WriteSxSy

  !********************************************************************!

  SUBROUTINE WriteSzRz( myThid )

    ! Writes source and receiver z-coordinates (depths)

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==

#ifdef IHOP_WRITE_OUT
    CALL WriteVector( Pos%nSZ, Pos%SZ, 'Source depths,   Sz', 'm', &
                    myThid )
    CALL WriteVector( Pos%nRZ, Pos%RZ, 'Receiver depths, Rz', 'm', &
                    myThid )

#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE WriteSzRz

  !********************************************************************!

  SUBROUTINE WriteRcvrRanges( myThid )

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

    REAL(KIND=_RL90) :: x(SIZE(Pos%RR))

    x = Pos%RR / 1000.0

#ifdef IHOP_WRITE_OUT
    ! IESCO22: assuming receiver positions are equally spaced
    CALL WriteVector( Pos%nRR, x, 'Receiver ranges, Rr', 'km', myThid )
#endif

    RETURN
  END !SUBROUTINE WriteRcvrRanges

  !********************************************************************!

#ifdef IHOP_THREED
  SUBROUTINE WriteRcvrBearings( myThid )   ! for 3D models

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==

! IEsco23: NOT SUPPORTED IN ihop
    CALL WriteVector( Pos%nTheta, Pos%theta, 'receiver bearings, theta', &
        'degrees', myThid )

    RETURN
  END !SUBROUTINE WriteRcvrBearings
#endif /* IHOP_THREED */
  !********************************************************************!

  SUBROUTINE WriteVector( Nx, x, Description, Units, myThid )

    ! Read a vector x
    ! Description is something like 'receiver ranges'
    ! Units       is something like 'km'

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    INTEGER,           INTENT( IN )    :: Nx
    REAL (KIND=_RL90), INTENT( INOUT ) :: x( : )
    CHARACTER,         INTENT( IN )    :: Description*( * ), Units*( * )
    INTEGER :: ix

#ifdef IHOP_WRITE_OUT
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)')'______________________________________________', &
                          '_____________'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,I10)') 'Number of ' // Description // ' = ', Nx
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

      WRITE(msgBuf,'(A)') Description // ' [' // Units // ']'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )


      WRITE(msgBuf,'(5G14.6)') ( x( ix ), ix = 1, MIN( Nx, Number_to_Echo ) )
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      IF ( Nx > Number_to_Echo ) THEN
        WRITE(msgBuf,'(G14.6)') ' ... ', x( Nx )
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      END IF

      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    ENDIF
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE WriteVector

END MODULE srpos_mod
