#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE initenvihop
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>

  ! mbp 12/2018, based on much older subroutine

  USE ihop_mod,  only: PRTFile, RAYFile, DELFile, ARRFile, SHDFile, Title, Beam
  USE ssp_mod,   only: initSSP, SSP
  USE bdry_mod,  only: Bdry, HSInfo
  USE atten_mod, only: CRCI

! ! USES
  implicit none
!  == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef ALLOW_CAL
#include "cal.h"
#endif /* ALLOW_CAL */

  PRIVATE

!   == Public Interfaces ==
!=======================================================================
  public    initEnv, OpenOutputFiles, resetMemory
!=======================================================================

! INPUT/OUTPUT PARAMETERS:

! == External Functions ==
  INTEGER  ILNBLNK
  EXTERNAL ILNBLNK

CONTAINS
  SUBROUTINE initEnv( myTime, myIter, myThid )

    ! I/O routine for acoustic fixed inputS

    USE angle_mod,  only: ReadRayElevationAngles
    USE srpos_mod,  only: Pos, ReadSxSy, ReadSzRz, ReadRcvrRanges, ReadFreqVec
#ifdef IHOP_THREED
    USE angle_mod,  only: ReadRayBearingAngles
    USE srpos_mod,  only: ReadRcvrBearings
#endif /* IHOP_THREED */

  ! == Routine Arguments ==
  ! myTime  :: Current time in simulation
  ! myIter  :: Current time-step number
  ! myThid  :: my Thread Id number
  ! msgBuf  :: Used to build messages for printing.
    _RL, INTENT(IN)     ::  myTime
    INTEGER, INTENT(IN) ::  myIter, myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  ! == Local Variables ==
    REAL (KIND=_RL90),  PARAMETER :: c0 = 1500.0
    REAL (KIND=_RL90)  :: x(2), c, cimag, gradc(2), crz, czz, rho, Depth
    CHARACTER (LEN= 2) :: AttenUnit
    CHARACTER (LEN=10) :: PlotType
    
!$TAF init initEnv1 = 'initenvihop_initenv'

    ! init local variables
    AttenUnit = ''
    PlotType  = ''

    !RG
    Bdry%Bot%HS = HSInfo(0.,0.,0.,0., 0.,0. , (0.,0.),(0.,0.), '', '' )
    Bdry%Top%HS = HSInfo(0.,0.,0.,0., 0.,0. , (0.,0.),(0.,0.), '', '' )
 
    ! *** ihop info to PRTFile ***
    CALL openPRTFile( myTime, myIter, myThid )

    ! *** Top Boundary ***
    Bdry%Top%HS%Opt = IHOP_topopt
    Bdry%Top%HS%Depth = 0 !initiate to dummy value

    CALL ReadTopOpt( Bdry%Top%HS%BC, AttenUnit, myThid )
    CALL TopBot( AttenUnit, Bdry%Top%HS, myThid )

    ! *** Bottom Boundary ***
    Bdry%Bot%HS%Opt = IHOP_botopt 
    IF ( IHOP_depth.NE.0 ) THEN
        Bdry%Bot%HS%Depth = IHOP_depth
    ELSE
        ! Extend by 5 wavelengths
        Bdry%Bot%HS%Depth = rkSign*rF( Nr+1 ) + 5*c0/IHOP_freq 
    END IF

#ifdef IHOP_WRITE_OUT
    !   Only do I/O in the main thread
    _BEGIN_MASTER(myThid)

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
     WRITE(msgBuf,'(A)') 
     CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
     WRITE(msgBuf,'(A,F10.2,A)' ) 'Depth = ',Bdry%Bot%HS%Depth,' m'
     CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
     WRITE(msgBuf,'(2A)') 'Top options: ', Bdry%Top%HS%Opt
     CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
     WRITE(msgBuf,'(A)') 
     CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
     WRITE(msgBuf,'(2A)') 'Bottom options: ', Bdry%Bot%HS%Opt
     CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

     SELECT CASE ( Bdry%Bot%HS%Opt( 2 : 2 ) )
       CASE ( '~', '*' )
#ifdef IHOP_WRITE_OUT
         WRITE(msgBuf,'(A)') '    Bathymetry file selected'
         CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
       CASE( ' ' )
       CASE DEFAULT
#ifdef IHOP_WRITE_OUT
         WRITE(msgBuf,'(2A)') 'INITENVIHOP initEnv: ', & 
                          'Unknown bottom option letter in second position'
         CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
         STOP 'ABNORMAL END: S/R initEnv'
       END SELECT
    ENDIF
    !   Only do I/O in the main thread
    _END_MASTER(myThid)
#endif /* IHOP_WRITE_OUT */

    Bdry%Bot%HS%BC = Bdry%Bot%HS%Opt( 1:1 )
    CALL TopBot( AttenUnit, Bdry%Bot%HS, myThid )

! IESCO24: Write derived type with allocatable memory by type: Bdry from bdry_mod
! Scalar components
!$TAF store bdry%top%hs%depth,bdry%top%hs%bc = initEnv1

! IESCO24: Write derived type with allocatable memory by type: SSP from ssp_mod
! Scalar components
! Fixed arrays
! Allocatable arrays
!$TAF store ssp%czmat,ssp%seg%r,ssp%seg%x,ssp%seg%y,ssp%seg%z = initEnv1

! IESCO24: Write derived type with allocatable memory by type: Pos from srpos_mod
! Scalar components
! Allocatable arrays
!$TAF store pos%rr = initEnv1
!$TAF store pos%theta,pos%wr,pos%ws = initEnv1

    ! *** Ocean SSP ***
    x = [ 0.0 _d 0, Bdry%Bot%HS%Depth ]   ! tells SSP Depth to read to
    CALL initSSP( x, myThid )

    Bdry%Top%HS%Depth = SSP%z( 1 )   ! first SSP point is top depth
    Depth = Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth ! water column depth

    ! *** Source locations ***
    CALL ReadSxSy( myThid ) ! Read source/receiver x-y coordinates

    Pos%NSz = IHOP_nsd
    Pos%NRz = IHOP_nrd

    CALL AllocatePos( Pos%NSz, Pos%Sz, IHOP_sd )
    CALL AllocatePos( Pos%NRz, Pos%Rz, IHOP_rd )
    CALL ReadSzRz( Bdry%Top%HS%Depth, Bdry%Bot%HS%Depth, myThid )

    ! *** Receiver locations ***
    Pos%NRr = IHOP_nrr
    CALL AllocatePos( Pos%NRr, Pos%Rr, IHOP_rr )
    CALL ReadRcvrRanges( myThid )
#ifdef IHOP_THREED
    CALL ReadRcvrBearings( myThid )
#endif /* IHOP_THREED */

    ! *** Broadband frequencies ***
    CALL ReadfreqVec( Bdry%Top%HS%Opt( 6:6 ), myThid )

    ! *** Run type ***
    Beam%RunType = IHOP_runopt
    CALL ReadRunType( Beam%RunType, PlotType, myThid )

    CALL ReadRayElevationAngles( Depth, Bdry%Top%HS%Opt, Beam%RunType, myThid )
#ifdef IHOP_THREED
    CALL ReadRayBearingAngles( Bdry%Top%HS%Opt, Beam%RunType, myThid )
#endif /* IHOP_THREED */


    ! *** Acoustic grid ***

    ! Limits for tracing beams
#ifdef IHOP_THREED
    WRITE(msgBuf,'(2A)') 'INITENVIHOP initEnv: ', & 
                         '3D not supported in ihop'
    CALL PRINT_ERROR( msgBuf,myThid )
    STOP 'ABNORMAL END: S/R initEnv'
    !READ(  ENVFile, * ) Beam%deltas, Beam%Box%x, Beam%Box%y, Beam%Box%z
    Beam%Box%x = 1000.0 * Beam%Box%x   ! convert km to m
    Beam%Box%y = 1000.0 * Beam%Box%y   ! convert km to m

    ! Automatic step size selection
    IF ( Beam%deltas == 0.0 ) Beam%deltas = &
        ( Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth ) / 10.0   
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
        WRITE(msgBuf,'(A,A,G11.4,A)') &
            ' Step length,', ' deltas = ', Beam%deltas, ' m'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') 
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,G11.4,A)') &
            ' Maximum ray x-range, Box%X = ', Beam%Box%X,' m'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,G11.4,A)') &
            ' Maximum ray y-range, Box%Y = ', Beam%Box%Y,' m'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,G11.4,A)') &
            ' Maximum ray z-range, Box%Z = ', Beam%Box%Z,' m'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    ENDIF
#endif /* IHOP_WRITE_OUT */
#else /* not IHOP_THREED */
    ! Step size in meters [m]
    Beam%deltas = IHOP_step
    
    IF ( Beam%deltas == 0.0 ) THEN ! Automatic step size option
        Beam%deltas = ( Depth ) / 10.   
    END IF
#ifdef IHOP_WRITE_OUT
    !   Only do I/O if in the main thread
    _BEGIN_MASTER(myThid)
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
        WRITE(msgBuf,'(A)') 
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(2A)')'__________________________________________', &
                            '_________________'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') 
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,G11.4,A)') &
            ' Step length, deltas = ', Beam%deltas, ' m'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    ENDIF
    !   Only do I/O in the main thread
    _END_MASTER(myThid)
#endif /* IHOP_WRITE_OUT */

    ! Domain size
    Beam%Box%Z = Bdry%Bot%HS%Depth ! in [m]
    ! Extend beam box by a single step size forward
    Beam%Box%R = IHOP_rr(nrd) + Beam%deltas/1000. ! in [km]

#ifdef IHOP_WRITE_OUT
    !   Only do I/O if in the main thread
    _BEGIN_MASTER(myThid)
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN

        WRITE(msgBuf,'(A)') 
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,G11.4,A)') &
            ' Maximum ray range, Box%R = ', Beam%Box%R,' km'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,G11.4,A)') &
            ' Maximum ray depth, Box%Z = ', Beam%Box%Z,' m'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    ENDIF
    !   Only do I/O in the main thread
    _END_MASTER(myThid)
#endif /* IHOP_WRITE_OUT */

    Beam%Box%R = Beam%Box%R*1000.   ! convert km to m
#endif /* IHOP_THREED */


    ! *** Beam characteristics ***
    Beam%Type( 4:4 ) = Beam%RunType( 7:7 )   ! selects beam shift option
          
#ifdef IHOP_WRITE_OUT
    !   Only do I/O if in the main thread
    _BEGIN_MASTER(myThid)
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
        SELECT CASE ( Beam%Type( 4:4 ) )
        CASE ( 'S' )
           WRITE(msgBuf,'(A)') ' Beam shift in effect'
           CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        CASE DEFAULT
           WRITE(msgBuf,'(A)') ' No beam shift in effect'
           CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        END SELECT
        WRITE(msgBuf,'(A)')
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    ENDIF
    !   Only do I/O in the main thread
    _END_MASTER(myThid)
#endif /* IHOP_WRITE_OUT */

    ! don't worry about the beam type if this is a ray trace run
    ! using 'e' requires Beam%Type to be set
    IF ( Beam%RunType( 1:1 ) /= 'R' .OR. Beam%RunType( 1:1 ) /= 'E' ) THEN 

       ! Beam%Type( 1 : 1 ) is
       !   'G' or '^' Geometric hat beams in Cartesian coordinates
       !   'g' Geometric hat beams in ray-centered coordinates
       !   'B' Geometric Gaussian beams in Cartesian coordinates
       !   'b' Geometric Gaussian beams in ray-centered coordinates
       !   'S' Simple Gaussian beams
       ! Beam%Type( 2 : 2 ) controls the setting of the beam width
       !   'F' space Filling
       !   'M' minimum width
       !   'W' WKB beams
       ! Beam%Type( 3 : 3 ) controls curvature changes on boundary reflections
       !   'D' Double
       !   'S' Single
       !   'Z' Zero
       ! Beam%Type( 4 : 4 ) selects whether beam shifts are implemented on 
       ! boundary reflection
       !   'S' yes
       !   'N' no

       ! Curvature change can cause overflow in grazing case
       ! Suppress by setting BeamType( 3 : 3 ) = 'Z'

       Beam%Type( 1:1 ) = Beam%RunType( 2:2 )

        SELECT CASE ( Beam%Type( 1:1 ) )
        CASE ( 'G', 'g' , '^', 'B', 'b', 'S' )   
        CASE DEFAULT
#ifdef IHOP_WRITE_OUT
            !   Only do I/O if in the main thread
            _BEGIN_MASTER(myThid)
            ! In adjoint mode we do not write output besides on the first run
            IF (IHOP_dumpfreq.GE.0) THEN

                WRITE(msgBuf,'(2A)') 'INITENVIHOP initEnv: ', & 
                                'Unknown beam type (second letter of run type)'
                CALL PRINT_ERROR( msgBuf,myThid )
            ENDIF ! No output in adjoint mode
            !   Only do I/O in the main thread
            _END_MASTER(myThid)
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R initEnv'
        END SELECT

    END IF

  RETURN
  END !SUBROUTINE initEnv

  !**********************************************************************!

  SUBROUTINE ReadTopOpt( BC, AttenUnit, myThid )
    USE atten_mod, only: T, Salinity, pH, z_bar, iBio, NBioLayers, bio

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    CHARACTER (LEN= 1), INTENT( OUT )   :: BC ! Boundary condition type
    CHARACTER (LEN= 2), INTENT( INOUT ) :: AttenUnit

    SSP%Type  = IHOP_TopOpt( 1:1 )
    BC        = IHOP_TopOpt( 2:2 )
    AttenUnit = IHOP_TopOpt( 3:4 )
    SSP%AttenUnit = AttenUnit

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) RETURN

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    ! SSP approximation options
    SELECT CASE ( SSP%Type )
    CASE ( 'N' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    N2-linear approximation to SSP'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'C' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    C-linear approximation to SSP'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'P' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    PCHIP approximation to SSP'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'S' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Spline approximation to SSP'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'Q' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Quad approximation to SSP'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'A' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Analytic SSP option'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INITENVIHOP ReadTopOpt: ', & 
                             'Unknown option for SSP approximation'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    ! Attenuation options

    SELECT CASE ( AttenUnit( 1:1 ) )
    CASE ( 'N' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Attenuation units: nepers/m'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'F' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Attenuation units: dB/mkHz'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'M' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Attenuation units: dB/m'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'W' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Attenuation units: dB/wavelength'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'Q' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Attenuation units: Q'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'L' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Attenuation units: Loss parameter'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INITENVIHOP ReadTopOpt: ', & 
                             'Unknown attenuation units'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    ! optional addition of volume attenuation using standard formulas

    SELECT CASE ( AttenUnit( 2:2 ) )
    CASE ( 'T' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    THORP volume attenuation added'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'F' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Francois-Garrison volume attenuation added'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
       WRITE( PRTFile, &
              "( ' T = ', G11.4, 'degrees   S = ', G11.4, ' psu   pH = ', G11.4, ' z_bar = ', G11.4, ' m' )" ) &
            T, Salinity, pH, z_bar
#endif /* IHOP_WRITE_OUT */
    CASE ( 'B' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Biological attenaution'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(A,I10)') '      Number of Bio Layers = ', NBioLayers
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

       DO iBio = 1, NBioLayers
          !READ( ENVFile, *  ) bio( iBio )%Z1, bio( iBio )%Z2, bio( iBio )%f0, &
          !                    bio( iBio )%Q, bio( iBio )%a0
          WRITE(msgBuf, '(A,F10.4,A)') '      Top    of layer = ', bio( iBio )%Z1, ' m'
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf, '(A,F10.4,A)') '      Bottom of layer = ', bio( iBio )%Z2, ' m'
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf, '(A,F10.4,A)') '      Resonance frequency = ', bio( iBio )%f0, &
                              ' Hz'
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf, '(A,F10.4)') '      Q  = ', bio( iBio )%Q
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf, '(A,F10.4)') '      a0 = ', bio( iBio )%a0
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
       END DO
#endif /* IHOP_WRITE_OUT */
    CASE ( ' ' )
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INITENVIHOP ReadTopOpt: ', & 
                             'Unknown top option letter in fourth position'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    SELECT CASE ( IHOP_TopOpt( 5:5 ) )
    CASE ( '~', '*' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Altimetry file selected'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( '-', '_', ' ' )
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INITENVIHOP ReadTopOpt: ', & 
                             'Unknown top option letter in fifth position'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    SELECT CASE ( IHOP_TopOpt( 6:6 ) )
    CASE ( 'I' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Development options enabled'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( ' ' )
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INITENVIHOP ReadTopOpt: ', & 
                             'Unknown top option letter in sixth position'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

  RETURN
  END !SUBROUTINE ReadTopOpt

  !**********************************************************************!

  SUBROUTINE ReadRunType( RunType, PlotType, myThid )

    ! Read the RunType variable and print to .prt file

    USE srPos_mod, only: Pos

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    CHARACTER (LEN= 7), INTENT( INOUT ) :: RunType
    CHARACTER (LEN=10), INTENT( INOUT ) :: PlotType

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) RETURN

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    SELECT CASE ( RunType( 1:1 ) )
    CASE ( 'R' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Ray trace run'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'E' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Eigenray trace run'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'I' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Incoherent TL calculation'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'S' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Semi-coherent TL calculation'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'C' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Coherent TL calculation'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'A' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Arrivals calculation, ASCII  file output'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'a' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Arrivals calculation, binary file output'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'e' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Eigenrays + Arrivals run, ASCII file output'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INITENVIHOP ReadRunType: ', & 
            'Unknown RunType selected'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadRunType'
    END SELECT

    SELECT CASE ( RunType( 2:2 ) )
    CASE ( 'C' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Cartesian beams'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'R' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Ray centered beams'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'S' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Simple gaussian beams'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'b' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Geometric gaussian beams in ray-centered coordinates'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'B' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Geometric gaussian beams in Cartesian coordinates'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'g' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Geometric hat beams in ray-centered coordinates'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE DEFAULT
       RunType( 2:2 ) = 'G'
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Geometric hat beams in Cartesian coordinates'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    END SELECT

    SELECT CASE ( RunType( 4:4 ) )
    CASE ( 'R' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Point source (cylindrical coordinates)'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'X' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Line source (Cartesian coordinates)'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE DEFAULT
       RunType( 4:4 ) = 'R'
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Point source (cylindrical coordinates)'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    END SELECT

    SELECT CASE ( RunType( 5:5 ) )
    CASE ( 'R' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(2A)') 'Rectilinear receiver grid: Receivers at', &
                           ' ( Rr( ir ), Rz( ir ) ) )'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
       PlotType = 'rectilin  '
    CASE ( 'I' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Irregular grid: Receivers at Rr( : ) x Rz( : )'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
       IF ( Pos%NRz /= Pos%NRr ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'INITENVIHOP ReadRunType: ', & 
                    'Irregular grid option selected with NRz not equal to Nr'
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadRunType'
       END IF
       PlotType = 'irregular '
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(2A)') 'Rectilinear receiver grid: Receivers at', &
                           ' Rr( : ) x Rz( : )'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
       RunType( 5:5 ) = 'R'
       PlotType = 'rectilin  '
    END SELECT

    SELECT CASE ( RunType( 6:6 ) )
    CASE ( '2' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'N x 2D calculation (neglects horizontal refraction)'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( '3' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '3D calculation'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE DEFAULT
       RunType( 6:6 ) = '2'
    END SELECT

  RETURN
  END !SUBROUTINE ReadRunType

  !**********************************************************************!

  SUBROUTINE TopBot( AttenUnit, HS, myThid )
    USE ssp_mod, only: rhoR, alphaR, betaR, alphaI, betaI

    ! Handles top and bottom boundary conditions

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    CHARACTER (LEN=2), INTENT( IN    ) :: AttenUnit
    TYPE ( HSInfo ),   INTENT( INOUT ) :: HS
    REAL (KIND=_RL90) :: Mz, vr, alpha2_f     ! values related to grain size
    REAL (KIND=_RL90) :: ztemp, bPower, fT

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
    ! Echo to PRTFile user's choice of boundary condition

    SELECT CASE ( HS%BC )
    CASE ( 'V' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Surface modeled as a VACUUM'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'R' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Perfectly RIGID'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'A' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    ACOUSTO-ELASTIC half-space'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'G' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Grain size to define half-space'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'F' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    FILE used for reflection loss'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'W' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Writing an IRC file'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'P' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    reading PRECALCULATED IRC'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INITENVIHOP TopBot: ', & 
                             'Unknown boundary condition type'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R TopBot'
    END SELECT
    ENDIF ! no output on adjoint runs

    ! ****** Read in BC parameters depending on particular choice ******
    HS%cp  = 0.0
    HS%cs  = 0.0
    HS%rho = 0.0

    ! RG recommends resetting to the default values from ssp_mod.F90
    bPower = 1.0
    fT     = 1D20
    rhoR   = 1.0
    
    SELECT CASE ( HS%BC )
    CASE ( 'A' )                  ! *** Half-space properties ***
       ! IEsco23: MISSING IF BOTTOM BC CHECK
       zTemp    = HS%Depth
       alphaR   = IHOP_bcsound
       betaR    = IHOP_bcsoundshear
       rhoR     = IHOP_brho
       alphaI   = IHOP_bcsoundI
       betaI    = IHOP_bcsoundshearI
#ifdef IHOP_WRITE_OUT
       ! In adjoint mode we do not write output besides on the first run
       IF (IHOP_dumpfreq.GE.0) THEN
        WRITE(msgBuf,'(A)') 
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') &
         '   z (m)     alphaR (m/s)   betaR  rho (g/cm^3)  alphaI     betaI'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )' ) &
             zTemp, alphaR, betaR, rhoR, alphaI, betaI
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
       ENDIF
#endif /* IHOP_WRITE_OUT */
       ! dummy parameters for a layer with a general power law for attenuation
       ! these are not in play because the AttenUnit for this is not allowed yet
       fT            = 1000.0

       HS%cp  = CRCI( zTemp, alphaR, alphaI, AttenUnit, bPower, fT, myThid )
       HS%cs  = CRCI( zTemp, betaR,  betaI,  AttenUnit, bPower, fT, myThid )

       HS%rho = rhoR
    CASE ( 'G' )            ! *** Grain size (formulas from UW-APL HF Handbook)
       ! These formulas are from the UW-APL Handbook
       ! The code is taken from older Matlab and is unnecesarily verbose
       ! vr   is the sound speed ratio
       ! rhoR is the density ratio
       !READ(  ENVFile, *    ) zTemp, Mz
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'( F10.2, 3X, F10.2 )' ) zTemp, Mz
       ! In adjoint mode we do not write output besides on the first run
       IF (IHOP_dumpfreq.GE.0) &
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

       IF ( Mz >= -1 .AND. Mz < 1 ) THEN
          vr   = 0.002709 * Mz**2 - 0.056452 * Mz + 1.2778
          rhoR = 0.007797 * Mz**2 - 0.17057  * Mz + 2.3139
       ELSE IF ( Mz >= 1 .AND. Mz < 5.3 ) THEN
          vr   = -0.0014881 * Mz**3 + 0.0213937 * Mz**2 - 0.1382798 * Mz &
               + 1.3425
          rhoR = -0.0165406 * Mz**3 + 0.2290201 * Mz**2 - 1.1069031 * Mz &
               + 3.0455
       ELSE
          vr   = -0.0024324 * Mz + 1.0019
          rhoR = -0.0012973 * Mz + 1.1565
       END IF

       IF ( Mz >= -1 .AND. Mz < 0 ) THEN
          alpha2_f = 0.4556
       ELSE IF ( Mz >= 0 .AND. Mz < 2.6 ) THEN
          alpha2_f = 0.4556 + 0.0245 * Mz
       ELSE IF( Mz >= 2.6 .AND. Mz < 4.5 ) THEN
          alpha2_f = 0.1978 + 0.1245 * Mz
       ELSE IF( Mz >= 4.5 .AND. Mz < 6.0 ) THEN
          alpha2_f = 8.0399 - 2.5228 * Mz + 0.20098 * Mz ** 2
       ELSE IF( Mz >= 6.0 .AND. Mz < 9.5 ) THEN
          alpha2_f = 0.9431 - 0.2041 * Mz + 0.0117 * Mz ** 2
       ELSE
          alpha2_f =  0.0601
       END IF

       ! AttenUnit = 'L'   ! loss parameter
!!! following uses a reference sound speed of 1500 ???
!!! should be sound speed in the water, just above the sediment
       ! the term vr / 1000 converts vr to units of m per ms 
       alphaR = vr * 1500.0
       ! loss parameter Sect. IV., Eq. (4) of handbook
       alphaI = alpha2_f * ( vr / 1000 ) * 1500.0 * log( 10.0 ) / ( 40.0*PI )

       HS%cp  = CRCI( zTemp, alphaR, alphaI, 'L ', bPower, fT, myThid )
       HS%cs  = 0.0
       HS%rho = rhoR
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A,2F10.2,3X,A,F10.2,3X,A,F10.2)') &
           'Converted sound speed =', HS%cp, 'density = ', rhoR, &
           'loss parm = ', alphaI
       ! In adjoint mode we do not write output besides on the first run
       IF (IHOP_dumpfreq.GE.0) &
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    END SELECT

  RETURN
  END !SUBROUTINE TopBot

  ! **********************************************************************!

  SUBROUTINE OpenOutputFiles( fName, myTime, myIter, myThid )
    ! Write appropriate header information

    USE angle_mod,  only: Angles
    USE srPos_mod,  only: Pos

    ! == Routine Arguments ==
    !  myTime   :: Current time in simulation
    !  myIter   :: Current time-step number
    !  myThid   :: my Thread Id number
    _RL, INTENT(IN)     ::  myTime
    INTEGER, INTENT(IN) ::  myIter, myThid

    ! == Local variables ==
    CHARACTER*(MAX_LEN_FNAM), INTENT(IN) :: fName
    CHARACTER*(MAX_LEN_FNAM) :: fullName
    INTEGER             :: IL
    REAL               :: atten
    CHARACTER (LEN=10) :: PlotType

    ! add time step to filename
    IF (myIter.GE.0) THEN
        IL=ILNBLNK( fName )
        WRITE(fullName, '(2A,I10.10)') fName(1:IL),'.',myIter
    ELSE
        fullName = fName
    ENDIF

    SELECT CASE ( Beam%RunType( 1:1 ) )
    CASE ( 'R', 'E' )   ! Ray trace or Eigenrays
#ifdef IHOP_WRITE_OUT
       OPEN ( FILE = TRIM( fullName ) // '.ray', UNIT = RAYFile, &
              FORM = 'FORMATTED' )
       WRITE( RAYFile, * ) '''', Title( 1 : 50 ), ''''
       WRITE( RAYFile, * ) IHOP_freq
       WRITE( RAYFile, * ) Pos%NSx, Pos%NSy, Pos%NSz
       WRITE( RAYFile, * ) Angles%Nalpha
       WRITE( RAYFile, * ) Bdry%Top%HS%Depth
       WRITE( RAYFile, * ) Bdry%Bot%HS%Depth

#ifdef IHOP_THREED
       WRITE( RAYFile, * ) Angles%Nalpha, Angles%Nbeta
       WRITE( RAYFile, * ) '''xyz'''
#else /* IHOP_THREED */
       WRITE( RAYFile, * ) '''rz'''
#endif /* IHOP_THREED */
       FLUSH( RAYFile )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'e' )        ! eigenrays + arrival file in ascii format
#ifdef IHOP_WRITE_OUT
       OPEN ( FILE = TRIM( fullName ) // '.arr', UNIT = ARRFile, &
              FORM = 'FORMATTED' )

# ifdef IHOP_THREED
       WRITE( ARRFile, * ) '''3D'''
# else /* IHOP_THREED */
       WRITE( ARRFile, * ) '''2D'''
# endif /* IHOP_THREED */

       WRITE( ARRFile, * ) IHOP_freq

       ! write source locations
# ifdef IHOP_THREED
       WRITE( ARRFile, * ) Pos%NSx,    Pos%Sx(    1 : Pos%NSx )
       WRITE( ARRFile, * ) Pos%NSy,    Pos%Sy(    1 : Pos%NSy )
       WRITE( ARRFile, * ) Pos%NSz,    Pos%Sz(    1 : Pos%NSz )
# else /* IHOP_THREED */
       WRITE( ARRFile, * ) Pos%NSz,    Pos%Sz(    1 : Pos%NSz )
# endif /* IHOP_THREED */

       ! write receiver locations
       WRITE( ARRFile, *    ) Pos%NRz,    Pos%Rz(    1 : Pos%NRz )
       WRITE( ARRFile, *    ) Pos%NRr,    Pos%Rr(    1 : Pos%NRr )
# ifdef IHOP_THREED
       WRITE( ARRFile, * ) Pos%Ntheta, Pos%theta( 1 : Pos%Ntheta )
# endif /* IHOP_THREED */

       ! IEsco22: add erays to arrivals output
       OPEN ( FILE = TRIM( fullName ) // '.ray', UNIT = RAYFile, &
              FORM = 'FORMATTED' )
       WRITE( RAYFile, * ) '''', Title( 1 : 50 ), ''''
       WRITE( RAYFile, * ) IHOP_freq
       WRITE( RAYFile, * ) Pos%NSx, Pos%NSy, Pos%NSz
       WRITE( RAYFile, * ) Angles%Nalpha
       WRITE( RAYFile, * ) Bdry%Top%HS%Depth
       WRITE( RAYFile, * ) Bdry%Bot%HS%Depth

# ifdef IHOP_THREED
       WRITE( RAYFile, * ) Angles%Nalpha, Angles%Nbeta
       WRITE( RAYFile, * ) '''xyz'''
# else /* IHOP_THREED */
       WRITE( RAYFile, * ) '''rz'''
# endif /* IHOP_THREED */

       IF (writeDelay) THEN
        OPEN ( FILE = TRIM( fullName ) // '.delay', UNIT = DELFile, &
               FORM = 'FORMATTED' )
        WRITE( DELFile, * ) '''', Title( 1 : 50 ), ''''
        WRITE( DELFile, * ) IHOP_freq
        WRITE( DELFile, * ) Pos%NSx, Pos%NSy, Pos%NSz
        WRITE( DELFile, * ) Angles%Nalpha
        WRITE( DELFile, * ) Bdry%Top%HS%Depth
        WRITE( DELFile, * ) Bdry%Bot%HS%Depth

#ifdef IHOP_THREED
        WRITE( DELFile, * ) Angles%Nalpha, Angles%Nbeta
        WRITE( DELFile, * ) '''xyz'''
# else /* IHOP_THREED */
        WRITE( DELFile, * ) '''rz'''
# endif /* IHOP_THREED */
       ENDIF

       FLUSH( RAYFile )
       IF (writeDelay) FLUSH( DELFile )
       FLUSH( ARRFile )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'A' )        ! arrival file in ascii format
#ifdef IHOP_WRITE_OUT
       OPEN ( FILE = TRIM( fullName ) // '.arr', UNIT = ARRFile, &
              FORM = 'FORMATTED' )

# ifdef IHOP_THREED
       WRITE( ARRFile, * ) '''3D'''
# else /* IHOP_THREED */
       WRITE( ARRFile, * ) '''2D'''
# endif /* IHOP_THREED */

       WRITE( ARRFile, * ) IHOP_freq

       ! write source locations
# ifdef IHOP_THREED
       WRITE( ARRFile, * ) Pos%NSx,    Pos%Sx(    1 : Pos%NSx )
       WRITE( ARRFile, * ) Pos%NSy,    Pos%Sy(    1 : Pos%NSy )
# endif /* IHOP_THREED */
       WRITE( ARRFile, * ) Pos%NSz,    Pos%Sz(    1 : Pos%NSz )

       ! write receiver locations
       WRITE( ARRFile, *    ) Pos%NRz,    Pos%Rz(    1 : Pos%NRz )
       WRITE( ARRFile, *    ) Pos%NRr,    Pos%Rr(    1 : Pos%NRr )
# ifdef IHOP_THREED
       WRITE( ARRFile, * ) Pos%Ntheta, Pos%theta( 1 : Pos%Ntheta )
# endif /* IHOP_THREED */
       FLUSH( ARRFile )
#endif /* IHOP_WRITE_OUT */
    CASE ( 'a' )        ! arrival file in binary format
#ifdef IHOP_WRITE_OUT
       OPEN ( FILE = TRIM( fullName ) // '.arr', UNIT = ARRFile, &
              FORM = 'UNFORMATTED' )

# ifdef IHOP_THREED
       WRITE( ARRFile ) '''3D'''
# else /* IHOP_THREED */
       WRITE( ARRFile ) '''2D'''
# endif /* IHOP_THREED */

       WRITE( ARRFile    ) SNGL( IHOP_freq )

       ! write source locations
# ifdef IHOP_THREED
       WRITE( ARRFile    ) Pos%NSx,    Pos%Sx(    1 : Pos%NSx )
       WRITE( ARRFile    ) Pos%NSy,    Pos%Sy(    1 : Pos%NSy )
       WRITE( ARRFile    ) Pos%NSz,    Pos%Sz(    1 : Pos%NSz )
# else /* IHOP_THREED */
       WRITE( ARRFile    ) Pos%NSz,    Pos%Sz(    1 : Pos%NSz )
# endif /* IHOP_THREED */

       ! write receiver locations
       WRITE( ARRFile       ) Pos%NRz,    Pos%Rz(    1 : Pos%NRz )
       WRITE( ARRFile       ) Pos%NRr,    Pos%Rr(    1 : Pos%NRr )
# ifdef IHOP_THREED
       WRITE( ARRFile    ) Pos%Ntheta, Pos%theta( 1 : Pos%Ntheta )
# endif /* IHOP_THREED */
       FLUSH( ARRFile )
#endif /* IHOP_WRITE_OUT */
    CASE DEFAULT
       atten = 0.0

       ! following to set PlotType has alread been done in READIN if that was 
       ! used for input
       SELECT CASE ( Beam%RunType( 5 : 5 ) )
       CASE ( 'R' )
          PlotType = 'rectilin  '
       CASE ( 'I' )
          PlotType = 'irregular '
       CASE DEFAULT
          PlotType = 'rectilin  '
       END SELECT

       CALL WriteSHDHeader( TRIM( fullName ) // '.shd', Title, REAL( IHOP_freq ), &
                         atten, PlotType )
    END SELECT

  RETURN
  END !SUBROUTINE OpenOutputFiles

  !**********************************************************************!

  SUBROUTINE WriteSHDHeader( FileName, Title, freq0, Atten, PlotType )

    USE srPos_mod,  only: Pos, Nfreq, freqVec

    ! Write header to disk file

    REAL,      INTENT( IN ) :: freq0, Atten      ! Nominal frequency, stabilizing attenuation (for wavenumber integration only)
    CHARACTER, INTENT( IN ) :: FileName*( * )    ! Name of the file (could be a shade file or a Green's function file)
    CHARACTER, INTENT( IN ) :: Title*( * )       ! Arbitrary title
    CHARACTER, INTENT( IN ) :: PlotType*( 10 )   ! 
    INTEGER :: LRecl

    ! receiver bearing angles
    IF ( .NOT. ALLOCATED( Pos%theta ) ) THEN
       ALLOCATE( Pos%theta( 1 ) )
       Pos%theta( 1 ) = 0   ! dummy bearing angle
       Pos%Ntheta     = 1
    END IF

    ! source x-coordinates
    IF ( .NOT. ALLOCATED( Pos%Sx ) ) THEN
       ALLOCATE( Pos%Sx( 1 ) )
       Pos%sx( 1 ) = 0      ! dummy x-coordinate
       Pos%NSx     = 1
    END IF

    ! source y-coordinates
    IF ( .NOT. ALLOCATED( Pos%Sy ) ) THEN
       ALLOCATE( Pos%Sy( 1 ) )
       Pos%sy( 1 ) = 0      ! dummy y-coordinate
       Pos%NSy     = 1
    END IF

    IF ( PlotType( 1 : 2 ) /= 'TL' ) THEN
       ! MAX( 41, ... ) below because Title is already 40 words (or 80 bytes)
 ! words/record (NRr doubled for complex pressure storage)
       LRecl = MAX( 41, 2 * Nfreq, Pos%Ntheta, Pos%NSx, Pos%NSy, Pos%NSz, &
                    Pos%NRz, 2 * Pos%NRr )  

       OPEN ( FILE = FileName, UNIT = SHDFile, STATUS = 'REPLACE', &
              ACCESS = 'DIRECT', RECL = 4 * LRecl, FORM = 'UNFORMATTED')
       WRITE( SHDFile, REC = 1  ) LRecl, Title( 1 : 80 )
       WRITE( SHDFile, REC = 2  ) PlotType
       WRITE( SHDFile, REC = 3  ) Nfreq, Pos%Ntheta, Pos%NSx, Pos%NSy, Pos%NSz,& 
                                  Pos%NRz, Pos%NRr, freq0, atten
       WRITE( SHDFile, REC = 4  ) freqVec(   1 : Nfreq )
       WRITE( SHDFile, REC = 5  ) Pos%theta( 1 : Pos%Ntheta )

       WRITE( SHDFile, REC = 6  ) Pos%Sx( 1 : Pos%NSx )
       WRITE( SHDFile, REC = 7  ) Pos%Sy( 1 : Pos%NSy )
       WRITE( SHDFile, REC = 8  ) Pos%Sz( 1 : Pos%NSz )

       WRITE( SHDFile, REC = 9  ) Pos%Rz( 1 : Pos%NRz )
       WRITE( SHDFile, REC = 10 ) Pos%Rr( 1 : Pos%NRr )

    ELSE   ! compressed format for TL from FIELD3D
  ! words/record (NR doubled for complex pressure storage)
       LRecl = MAX( 41, 2 * Nfreq, Pos%Ntheta, Pos%NSz, Pos%NRz, 2 * Pos%NRr ) 

       OPEN ( FILE = FileName, UNIT = SHDFile, STATUS = 'REPLACE', &
              ACCESS = 'DIRECT', RECL = 4 * LRecl, FORM = 'UNFORMATTED')
       WRITE( SHDFile, REC = 1  ) LRecl, Title( 1 : 80 )
       WRITE( SHDFile, REC = 2  ) PlotType
       WRITE( SHDFile, REC = 3  ) Nfreq, Pos%Ntheta, Pos%NSx, Pos%NSy, Pos%NSz,& 
                                  Pos%NRz, Pos%NRr, freq0, atten
       WRITE( SHDFile, REC = 4  ) freqVec(   1 : Nfreq )
       WRITE( SHDFile, REC = 5  ) Pos%theta( 1 : Pos%Ntheta )

       WRITE( SHDFile, REC = 6  ) Pos%Sx( 1 ), Pos%Sx( Pos%NSx )
       WRITE( SHDFile, REC = 7  ) Pos%Sy( 1 ), Pos%Sy( Pos%NSy )
       WRITE( SHDFile, REC = 8  ) Pos%Sz( 1 : Pos%NSz )

       WRITE( SHDFile, REC = 9  ) Pos%Rz( 1 : Pos%NRz )
       WRITE( SHDFile, REC = 10 ) Pos%Rr( 1 : Pos%NRr )
    END IF

  RETURN
  END !SUBROUTINE WriteSHDHeader

  !**********************************************************************!

  SUBROUTINE WriteSHDField( P, NRz, NRr, IRec )

    ! Write the field to disk

    INTEGER, INTENT( IN )    :: NRz, NRr      ! # of receiver depths, ranges
    COMPLEX, INTENT( IN )    :: P( NRz, NRr ) ! Pressure field
    INTEGER, INTENT( INOUT ) :: iRec          ! last record read
    INTEGER                  :: iRz

    DO iRz = 1, NRz
       iRec = iRec + 1
       WRITE( SHDFile, REC = iRec ) P( iRz, : )
    END DO

  RETURN
  END !SUBROUTINE WriteSHDField

  !**********************************************************************!

  SUBROUTINE AllocatePos( Nx, x_out, x_in )

    ! Allocate and populate Pos structure from data.ihop

    INTEGER,          INTENT( IN  ) :: Nx    
    REAL(KIND=_RL90), INTENT( IN  ) :: x_in(:)
    REAL(KIND=_RL90), ALLOCATABLE, INTENT( OUT ) :: x_out(:)
    INTEGER                         :: i

    IF ( ALLOCATED(x_out) ) DEALLOCATE(x_out)
    ALLOCATE( x_out(MAX(3, Nx)) )

    ! set default values
    x_out    = 0.0
    x_out(3) = -999.9

    DO i = 1, Nx
        x_out(i) = x_in(i)
    END DO

  RETURN
  END !SUBROUTINE AllocatePos

  !**********************************************************************!
  SUBROUTINE openPRTFile ( myTime, myIter, myThid )

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     myTime :: Current time in simulation
  !     myIter :: Current time-step number
  !     msgBuf :: Used to build messages for printing.
    _RL, INTENT(IN)     ::  myTime
    INTEGER, INTENT(IN) ::  myIter, myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Arguments ==
#ifdef ALLOW_CAL
    INTEGER :: mydate(4)
#endif

    !   Only do I/O in the main thread
    _BARRIER
    _BEGIN_MASTER(myThid)

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) RETURN

#ifdef IHOP_WRITE_OUT
    WRITE(msgbuf,'(A)') 'iHOP Print File'
    CALL PRINT_MESSAGE( msgBuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgbuf,'(A)') 
    CALL PRINT_MESSAGE( msgBuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    ! *** TITLE ***
#ifdef IHOP_THREED
    WRITE(msgBuf,'(2A)') 'INITENVIHOP openPRTFile: ', & 
                         '3D not supported in ihop'
    CALL PRINT_ERROR( msgBuf,myThid )
    STOP 'ABNORMAL END: S/R initEnv'
    Title( 1 :11 ) = 'iHOP3D - '
    Title( 12:80 ) = IHOP_title
#else /* not IHOP_THREED */
    Title( 1 : 9 ) = 'iHOP - '
    Title( 10:80 ) = IHOP_title
#endif /* IHOP_THREED */

#ifdef IHOP_WRITE_OUT
    WRITE(msgbuf,'(A)') Title ! , ACHAR(10)
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)')'___________________________________________________', &
                        '________'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgbuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE( msgBuf, '(A,I10,A,F20.2,A)') 'GCM iter ', myIter,' at time = ', &
        myTime,' [sec]'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
# ifdef ALLOW_CAL
    CALL CAL_GETDATE( myIter,myTime,mydate,myThid )
    WRITE (msgBuf,'(A,I8,I6,I3,I4)') 'GCM cal date ', mydate(1), mydate(2), &
                                                       mydate(3), mydate(4)
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
# endif /* ALLOW_CAL */
    WRITE(msgbuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE( msgBuf, '(A,F11.4,A)' )'Frequency ', IHOP_freq, ' [Hz]'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)')'___________________________________________________', &
                        '________'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    !   Only do I/O in the main thread
    _END_MASTER(myThid)

  END !SUBROUTINE openPRTFile


  !**********************************************************************!

  SUBROUTINE resetMemory()
    USE srpos_mod,  only: Pos
    USE bdry_mod,   only: Top,Bot
    USE angle_mod,  only: Angles
    USE arr_mod,    only: Narr, Arr
    USE ihop_mod,   only: ray2D, MaxN, iStep
    
    ! From angle_mod
    IF (ALLOCATED(Angles%alpha))DEALLOCATE(Angles%alpha)
#ifdef IHOP_THREED
    IF (ALLOCATED(Angles%beta)) DEALLOCATE(Angles%beta)
#endif /* IHOP_THREED */
    ! From bdry_mod
    IF (ALLOCATED(Top))         DEALLOCATE(Top)
    IF (ALLOCATED(Bot))         DEALLOCATE(Bot)
    ! From bellhop
    IF (ALLOCATED(Pos%theta))   DEALLOCATE(Pos%theta)
    IF (ALLOCATED(Arr))         DEALLOCATE(Arr)
    IF (ALLOCATED(NArr))        DEALLOCATE(NArr)
    ! from initenvihop
    IF (ALLOCATED(Pos%Sx))      DEALLOCATE(Pos%Sx)
    IF (ALLOCATED(Pos%Sy))      DEALLOCATE(Pos%Sy)
    IF (ALLOCATED(Pos%Sz))      DEALLOCATE(Pos%Sz)
    ! From srpos_mod
    IF (ALLOCATED(Pos%ws))      DEALLOCATE(Pos%ws)
    IF (ALLOCATED(Pos%isz))     DEALLOCATE(Pos%isz)
    IF (ALLOCATED(Pos%wr))      DEALLOCATE(Pos%wr)
    IF (ALLOCATED(Pos%irz))     DEALLOCATE(Pos%irz)
    IF (ALLOCATED(Pos%rr))      DEALLOCATE(Pos%rr)
    IF (ALLOCATED(Pos%rz))      DEALLOCATE(Pos%rz)
    ! From ssp_mod
    IF (ALLOCATED(SSP%cMat))    DEALLOCATE(SSP%cMat)
    IF (ALLOCATED(SSP%czMat))   DEALLOCATE(SSP%czMat)
#ifdef IHOP_THREED
    IF (ALLOCATED(SSP%cMat3))   DEALLOCATE(SSP%cMat3)
    IF (ALLOCATED(SSP%czMat3))  DEALLOCATE(SSP%czMat3)
#endif /* IHOP_THREED */
    IF (ALLOCATED(SSP%Seg%r))   DEALLOCATE(SSP%Seg%r)
    ! From ihop_mod
    DO iStep = 1,MaxN
        ray2D(iStep)%x = [zeroRL, zeroRL]
        ray2D(iStep)%t = [zeroRL, zeroRL]
        ray2D(iStep)%p = [zeroRL, zeroRL]
        ray2D(iStep)%q = [zeroRL, zeroRL]
        ray2D(iStep)%c = zeroRL
        ray2D(iStep)%Amp = zeroRL
        ray2D(iStep)%Phase = zeroRL
        ray2D(iStep)%tau = (zeroRL, zeroRL)
    END DO

  END !SUBROUTINE resetMemory

END MODULE initenvihop
