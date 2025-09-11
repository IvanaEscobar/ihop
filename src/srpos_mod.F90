#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: srpos_mod
MODULE srpos_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   This module contains subroutines to read and write source and receiver
!   positions, including source depths, receiver depths, receiver ranges,
!   and receiver bearings.

! !USES:
  USE sort_mod,         only: Sort
  USE subTab_mod,       only: SubTab
  USE monotonic_mod,    only: monotonic
  USE ihop_mod,         only: PRTFile
  IMPLICIT NONE
!  == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC Pos, Nfreq, freqVec, &
        ReadSxSy,  ReadSzRz,  ReadRcvrRanges,  ReadFreqVec, &
        WriteSxSy, WriteSzRz, WriteRcvrRanges, WriteFreqVec
#ifdef IHOP_THREED
  PUBLIC ReadRcvrBearings, WriteRcvrBearings
#endif /* IHOP_THREED */
!=======================================================================

! == Module variables ==
  INTEGER, PARAMETER    :: Number_to_Echo = 10
  INTEGER, PRIVATE      :: IAllocStat     ! Status after memory allocation
  INTEGER               :: Nfreq = 1      ! number of frequencies
  REAL (KIND=_RL90), ALLOCATABLE  :: freqVec( : )   ! frequency vector for braodband runs

! == Derived types ==
  TYPE Position
    INTEGER              :: nSX = 1, nSY = 1, nSZ = 1, & ! # of x,y,z coords
                            nRZ = 1, nRR = 1, nTheta = 1 ! # of z,r,theta coord`s
    REAL                 :: Delta_r, Delta_theta
    INTEGER,           ALLOCATABLE :: iSz( : ), iRz( : )
    REAL (KIND=_RL90), ALLOCATABLE :: Sx( : ), Sy( : ), Sz( : )          ! Source x, y, z coordinates
    REAL (KIND=_RL90), ALLOCATABLE :: Rr( : ), Rz( : ), ws( : ), wr( : ) ! Receiver r, z coordinates and weights for interpolation
    REAL (KIND=_RL90), ALLOCATABLE :: theta( : )                         ! Receiver bearings
  END TYPE Position

  TYPE ( Position ) :: Pos ! structure containing source and receiver positions
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R ReadfreqVec
! S/R ReadSxSy
! S/R ReadSzRz
! S/R ReadRcvrRanges
! S/R ReadRcvrBearings
! S/R ReadVector
! S/R writeFreqVec
! S/R WriteSxSy
! S/R WriteSzRz
! S/R WriteRcvrRanges
! S/R WriteRcvrBearings
! S/R WriteVector
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ReadfreqVec
! !INTERFACE:
  SUBROUTINE ReadfreqVec( BroadbandOption, myThid )
! !DESCRIPTION:
!   Reads a vector of source frequencies for a broadband run
!   If the broadband option is not selected, then the input freq (a scalar)
!   is stored in the frequency vector
!   IHOP_freq is source frequency

! !USES: None

! !INPUT PARAMETERS:
! BroadbandOption :: Character variable indicating whether the run is broadband
! myThid :: my thread ID
  CHARACTER, INTENT( IN ) :: BroadbandOption*( 1 )
  INTEGER,   INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! ifreq :: Frequency index
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: ifreq
!EOP

  ! In adjoint mode we do not write output besides on the first run
  IF ( IHOP_dumpfreq.GE.0 ) THEN
    ! Broadband run?
    IF ( BroadbandOption.EQ.'B' ) THEN
      IF ( Nfreq.LE.0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadfreqVec: ', &
          'Number of frequencies must be positive'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadfreqVec'
      ENDIF
    ENDIF
  ENDIF

  IF ( ALLOCATED( freqVec ) ) DEALLOCATE( freqVec )
  ALLOCATE( freqVec( MAX( 3, Nfreq ) ), Stat=IAllocStat )
  IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadfreqVec: ', &
      'Number of frequencies too large'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadfreqVec'
  ENDIF

  ! set default values
  freqVec = 0.0

  IF ( BroadbandOption.EQ.'B' ) THEN
    freqVec(3) = -999.9
    !READ(  ENVFile, * ) freqVec( 1 : Nfreq )
    CALL SubTab( freqVec, Nfreq )

  ELSE
    freqVec(1) = IHOP_freq

  ENDIF

  RETURN
  END !SUBROUTINE ReadfreqVec

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ReadSxSy
! !INTERFACE:
  SUBROUTINE ReadSxSy( myThid )
! !DESCRIPTION:
!   Reads source x-y coordinates.

! !USES: None

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
!EOP

#ifdef IHOP_THREED
  CALL ReadVector( Pos%nSX, Pos%SX, 'source   x-coordinates, Sx', 'km', &
                  myThid )
  CALL ReadVector( Pos%nSY, Pos%SY, 'source   y-coordinates, Sy', 'km', &
                  myThid )

#else /* IHOP_THREED */
  ALLOCATE( Pos%SX( 1 ), Pos%SY( 1 ), Stat=IAllocStat )
  IF ( IAllocStat.NE.0 ) THEN
# ifdef IHOP_WRITE_OUT
    WRITE(msgBuf, *) 'allocation failed', IAllocStat
    CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadSxSy'
  ENDIF

  Pos%SX( 1 ) = 0.
  Pos%SY( 1 ) = 0.
#endif /* IHOP_THREED */

  RETURN
  END !SUBROUTINE ReadSxSy

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ReadSzRz
! !INTERFACE:
  SUBROUTINE ReadSzRz( zMin, zMax, myThid )
! !DESCRIPTION:
!   Reads source and receiver z-coordinates (depths).
!   zMin and zMax are limits for those depths; sources and receivers are
!   shifted to be within those limits.

! !USES: None

! !INPUT PARAMETERS:
! zMin :: Minimum depth (m)
! zMax :: Maximum depth (m)
! myThid :: my thread ID
  REAL (KIND=_RL90), INTENT( IN ) :: zMin, zMax
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! posShift :: Variable to indicate if any positions were shifted
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: posShift
!EOP

  CALL ReadVector( Pos%nSZ, Pos%SZ, 'Source   depths, Sz', 'm', &
                  myThid )
  CALL ReadVector( Pos%nRZ, Pos%RZ, 'Receiver depths, Rz', 'm', &
                  myThid )

  ! *** Check for Sz/Rz in water column ***
  posShift = 0
  IF ( ANY( Pos%SZ( 1:Pos%nSZ ).LT.zMin ) ) THEN
    WHERE ( Pos%SZ.LT.zMin ) Pos%SZ = zMin
    posShift = 1
  ENDIF

  IF ( ANY( Pos%SZ( 1:Pos%nSZ ).GT.zMax ) ) THEN
    WHERE( Pos%SZ .GT. zMax ) Pos%SZ = zMax
    posShift = 2
  ENDIF

  IF ( ANY( Pos%RZ( 1:Pos%nRZ ).LT.zMin ) ) THEN
    WHERE( Pos%RZ.LT.zMin ) Pos%RZ = zMin
    posShift = 3
  ENDIF

  IF ( ANY( Pos%RZ( 1:Pos%nRZ ).GT.zMax ) ) THEN
    WHERE( Pos%RZ.GT.zMax ) Pos%RZ = zMax
    posShift = 4
  ENDIF

#ifdef IHOP_WRITE_OUT
  ! In adjoint mode we do not write output besides on the first run
  IF ( IHOP_dumpfreq.GE.0 ) THEN
    IF ( posShift.EQ.1 ) THEN
      WRITE(msgBuf,'(2A)') 'Warning in WriteSzRz : Source above or too ',&
        'near the top bdry has been moved down'
      CALL PRINT_ERROR( msgBuf,myThid )
    ENDIF

    IF ( posShift.EQ.2 ) THEN
      WRITE(msgBuf,'(2A)') 'Warning in WriteSzRz : Source below or too ',&
        'near the bottom bdry has been moved up'
      CALL PRINT_ERROR( msgBuf,myThid )
    ENDIF

    IF ( posShift.EQ.3 ) THEN
      WRITE(msgBuf,'(2A)') 'Warning in WriteSzRz : Receiver above or too ',&
        'near the top bdry has been moved down'
      CALL PRINT_ERROR( msgBuf,myThid )
    ENDIF

    IF ( posShift.EQ.4 ) THEN
      WRITE(msgBuf,'(2A)') 'Warning in WriteSzRz : Receiver below or too ',&
        'near the bottom bdry has been moved up'
      CALL PRINT_ERROR( msgBuf,myThid )
    ENDIF

  ENDIF
#endif /* IHOP_WRITE_OUT */

  IF ( ALLOCATED( Pos%ws ) ) DEALLOCATE( Pos%ws, Pos%iSz )
  ALLOCATE( Pos%ws( Pos%nSZ ), Pos%iSz( Pos%nSZ ), Stat=IAllocStat )
  IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 'SRPOSITIONS ReadSzRz: Too many sources'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadSzRz'
  ENDIF

  IF ( ALLOCATED( Pos%wr ) ) DEALLOCATE( Pos%wr, Pos%iRz )
  ALLOCATE( Pos%wr( Pos%nRZ ), Pos%iRz( Pos%nRZ ), Stat=IAllocStat )
  IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(A)') 'SRPOSITIONS ReadSzRz: Too many receivers'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadSzRz'
  ENDIF

  ! Set default values
  Pos%ws = 0
  Pos%isz = 0
  Pos%wr = 0
  Pos%irz = 0

  RETURN
  END !SUBROUTINE ReadSzRz

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ReadRcvrRanges
! !INTERFACE:
  SUBROUTINE ReadRcvrRanges( myThid )
! !DESCRIPTION:
!   Reads receiver ranges (distances from the source).

! !USES: None

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
!EOP

  ! IESCO22: assuming receiver positions are equally spaced
  CALL ReadVector( Pos%nRR, Pos%RR, 'Receiver ranges, Rr', 'km', myThid )

  ! calculate range spacing
  Pos%delta_r = 0.0
  ! IESCO24: Assumes uniform spacing in ranges btwn receivers
  IF ( Pos%nRR.NE.1 ) Pos%delta_r = Pos%RR( Pos%nRR ) - Pos%RR( Pos%nRR-1 )

  IF ( .NOT.monotonic( Pos%RR, Pos%nRR ) ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadRcvrRanges: ', &
      'Receiver ranges are not monotonically increasing'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadRcvrRanges'
  ENDIF

  RETURN
  END !SUBROUTINE ReadRcvrRanges

#ifdef IHOP_THREED
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ReadRcvrBearings
! !INTERFACE:
  SUBROUTINE ReadRcvrBearings( myThid )
! !DESCRIPTION:
!   Reads receiver bearings (angles in degrees) for 3D models.

! !USES: None

! !INPUT PARAMETERS:
! myThid :: my thread ID
INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
!EOP

! IEsco23: 3D NOT SUPPORTED IN ihop
  CALL ReadVector( Pos%nTheta, Pos%theta, 'receiver bearings, theta', &
    'degrees', myThid )

  ! full 360-degree sweep? remove duplicate angle
  IF ( Pos%nTheta.GT.1 ) THEN
      IF ( ABS( MOD( Pos%theta( Pos%nTheta ) - Pos%theta( 1 ), 360.0 ) ) &
          .LT.10.0*TINY( 1.0D0 ) ) &
        Pos%nTheta = Pos%nTheta - 1
  ENDIF

  ! calculate angular spacing
  Pos%Delta_theta = 0.0
  IF ( Pos%nTheta.NE.1 ) Pos%Delta_theta = Pos%theta( Pos%nTheta ) &
                                          - Pos%theta( Pos%nTheta-1 )

  IF ( .NOT.monotonic( Pos%theta, Pos%nTheta ) ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SRPOSITIONS ReadRcvrBearings: ', &
      'Receiver bearings are not monotonically increasing'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadRcvrBearings'
  ENDIF

  RETURN
  END !SUBROUTINE ReadRcvrBearings

#endif /* IHOP_THREED */
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ReadVector
! !INTERFACE:
  SUBROUTINE ReadVector( Nx, x, Description, Units, myThid )
! !DESCRIPTION:
!   Reads a vector x for source/receiver positions.

! !USES:
  USE sort_mod, only: Sort

! !INPUT PARAMETERS:
! Nx :: Number of elements in the vector
! x  :: Vector to be read
! Description :: Description of the vector (e.g., 'receiver ranges')
! Units :: Units of the vector (e.g., 'km')
! myThid :: my thread ID
  INTEGER,                        INTENT( IN )    :: Nx
  REAL (KIND=_RL90), ALLOCATABLE, INTENT( INOUT ) :: x( : )
  CHARACTER,                      INTENT( IN )    :: Description*( * ), &
                                                     Units*( * )
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: x

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! ix :: Index variable
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: ix
!EOP

  IF ( Nx.LE.0 ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SRPOS_MOD ReadVector: ', &
      'Number of ' // Description // 'must be positive'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadVector'
  ENDIF

  IF ( .NOT.ALLOCATED( x ) ) THEN
    ALLOCATE( x( MAX( 3, Nx ) ), Stat=IAllocStat )
    IF ( IAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'SRPOS_MOD ReadVector: ', &
        'Too many ' // Description
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ReadVector'
    ENDIF
  ENDIF

  CALL SubTab( x, Nx )
  CALL Sort(   x, Nx )

  ! Vectors in km should be converted to m for internal use
  IF ( LEN_TRIM( Units ).GE.2 ) THEN
      IF ( Units( 1:2 ).EQ.'km' ) x = 1000.0 * x
  ENDIF

  RETURN
  END !SUBROUTINE ReadVector

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: writeFreqVec
! !INTERFACE:
  SUBROUTINE writeFreqVec( BroadbandOption, myThid )
! !DESCRIPTION:
!   Writes a vector of source frequencies for a broadband run.

! !USES: None

! !INPUT PARAMETERS:
! BroadbandOption :: Character variable indicating whether the run is broadband
! myThid :: my thread ID
  CHARACTER*(1), INTENT( IN ) :: BroadbandOption
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! ifreq :: Frequency index
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: ifreq
!EOP

#ifdef IHOP_WRITE_OUT
  ! In adjoint mode we do not write output besides on the first run
  IF ( IHOP_dumpfreq.GE.0 ) THEN
    ! Broadband run?
    IF ( BroadbandOption.EQ.'B' ) THEN
      WRITE(msgBuf,'(A)') &
        '___________________________________________________________'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,I10)') 'Number of frequencies =', Nfreq
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

      WRITE(msgBuf,'(A)') 'Frequencies (Hz)'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

      WRITE( msgBuf, '(5G14.6)' ) ( freqVec( ifreq ), ifreq=1, &
          MIN( Nfreq, Number_to_Echo ) )
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      IF ( Nfreq .GT. Number_to_Echo ) THEN
        WRITE( msgBuf,'(G14.6)' ) ' ... ', freqVec( Nfreq )
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      ENDIF

    ENDIF ! Broadband run
    
  ENDIF ! adjoint run?
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE writeFreqVec

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteSxSy
! !INTERFACE:
  SUBROUTINE WriteSxSy( myThid )
! !DESCRIPTION:
!  Writes source x-y coordinates.

! !USES: None

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
!EOP

#ifdef IHOP_THREED
# ifdef IHOP_WRITE_OUT
  CALL WriteVector( Pos%nSX, Pos%SX, 'source   x-coordinates, Sx', 'km', &
                  myThid )
  CALL WriteVector( Pos%nSY, Pos%SY, 'source   y-coordinates, Sy', 'km', &
                  myThid )
# endif /* IHOP_WRITE_OUT */
#endif /* IHOP_THREED */

  RETURN
  END !SUBROUTINE WriteSxSy

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteSzRz
! !INTERFACE:
  SUBROUTINE WriteSzRz( myThid )
! !DESCRIPTION:
!  Writes source z-coordinates (depths).

! !USES: None

! !INPUT PARAMETERS:
! myThid :: my thread ID
    INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
!EOP

#ifdef IHOP_WRITE_OUT
  CALL WriteVector( Pos%nSZ, Pos%SZ, 'Source depths,   Sz', 'm', &
                  myThid )
  CALL WriteVector( Pos%nRZ, Pos%RZ, 'Receiver depths, Rz', 'm', &
                  myThid )
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE WriteSzRz

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteRcvrRanges
! !INTERFACE:
  SUBROUTINE WriteRcvrRanges( myThid )
! !DESCRIPTION:
!  Writes receiver ranges.

! !USES: None

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! x :: Vector of receiver ranges in km
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
  REAL(KIND=_RL90) :: x( SIZE(Pos%RR) )
!EOP

  x = Pos%RR / 1000.0
#ifdef IHOP_WRITE_OUT
  ! IESCO22: assuming receiver positions are equally spaced
  CALL WriteVector( Pos%nRR, x, 'Receiver ranges, Rr', 'km', myThid )
#endif

  RETURN
  END !SUBROUTINE WriteRcvrRanges

#ifdef IHOP_THREED
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteRcvrBearings
! !INTERFACE:
  SUBROUTINE WriteRcvrBearings( myThid )   ! for 3D models
! !DESCRIPTION:
!  Writes receiver bearings.

! !USES: None

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
!EOP

! IEsco23: NOT SUPPORTED IN ihop
  CALL WriteVector( Pos%nTheta, Pos%theta, 'receiver bearings, theta', &
    'degrees', myThid )

  RETURN
  END !SUBROUTINE WriteRcvrBearings

#endif /* IHOP_THREED */
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteVector
! !INTERFACE:
  SUBROUTINE WriteVector( Nx, x, Description, Units, myThid )
! !DESCRIPTION:
!   Writes a vector x for source/receiver positions.

! !USES: None

! !INPUT PARAMETERS:
! Nx :: Number of elements in the vector
! x  :: Vector to be written
! Description :: Description of the vector (e.g., 'receiver ranges')
! Units :: Units of the vector (e.g., 'km')
! myThid :: my thread ID
  INTEGER,           INTENT( IN )    :: Nx
  REAL (KIND=_RL90), INTENT( INOUT ) :: x( : )
  CHARACTER,         INTENT( IN )    :: Description*( * ), &
                                        Units*( * )
  INTEGER,           INTENT( IN )    :: myThid
! !OUTPUT PARAMETERS: x

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! ix :: Index variable
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
    INTEGER :: ix
!EOP

#ifdef IHOP_WRITE_OUT
  ! In adjoint mode we do not write output besides on the first run
  IF (IHOP_dumpfreq.GE.0) THEN
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') &
      '___________________________________________________________'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A,I10)') 'Number of ' // Description // ' = ', Nx
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    WRITE(msgBuf,'(A)') Description // ' [' // Units // ']'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )


    WRITE(msgBuf,'(5G14.6)') ( x( ix ), ix=1, MIN( Nx, Number_to_Echo ) )
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    IF ( Nx .GT. Number_to_Echo ) THEN
      WRITE(msgBuf,'(G14.6)') ' ... ', x( Nx )
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF

    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  ENDIF
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE WriteVector

END !MODULE srpos_mod
