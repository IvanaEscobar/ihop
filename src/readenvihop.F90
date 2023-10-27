#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE readEnviHop
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>

  ! mbp 12/2018, based on much older subroutine

  USE ihop_mod,     only: PRTFile, RAYFile, DELFile, ARRFile, SHDFile, &
                          Title, Beam
  USE ssp_mod,      only: EvaluateSSP, HSInfo, Bdry, SSP, zTemp, alphaR, betaR,&
                          alphaI, betaI, rhoR, betaPowerLaw, fT
  USE atten_mod,    only: CRCI, T, Salinity, pH, z_bar, iBio, NBioLayers, bio

! ! USES
  implicit none
!  == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

  PRIVATE

! public interfaces
!=======================================================================

  public    ReadEnvironment, OpenOutputFiles
  
!=======================================================================

CONTAINS
  SUBROUTINE ReadEnvironment( FileRoot, myThid )

    ! I/O routine for acoustic fixed inputS

    USE angle_mod,  only: ReadRayElevationAngles, ReadRayBearingAngles
    USE srpos_mod,  only: Pos, ReadSxSy, ReadSzRz, ReadRcvrRanges,         &
#ifdef IHOP_THREED
                              ReadRcvrBearings, &
#endif /* IHOP_THREED */
                              ReadFreqVec

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    REAL (KIND=_RL90),  PARAMETER   :: c0 = 1500.0
    CHARACTER (LEN=80), INTENT(IN ) :: FileRoot
    REAL               :: ZMin, ZMax
    REAL (KIND=_RL90)  :: x( 2 ), c, cimag, gradc( 2 ), crr, crz, czz, rho, &
                          Depth
    CHARACTER (LEN= 2) :: AttenUnit
    CHARACTER (LEN=10) :: PlotType

    !   Only do I/O if in the main thread
    _BARRIER
    _BEGIN_MASTER(myThid)

#ifdef IHOP_WRITE_OUT
    WRITE(msgbuf,'(A)') 'iHOP Print File'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgbuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    ! *** TITLE ***
#ifdef IHOP_THREED
    WRITE(msgBuf,'(2A)') 'READENVIHOP ReadEnvironment: ', & 
                         '3D not supported in ihop'
    CALL PRINT_ERROR( msgBuf,myThid )
    STOP 'ABNORMAL END: S/R ReadEnvironment'
    Title( 1 :11 ) = 'iHOP3D - '
    Title( 12:80 ) = IHOP_title
#else /* not IHOP_THREED */
    Title( 1 : 9 ) = 'iHOP - '
    Title( 10:80 ) = IHOP_title
#endif /* IHOP_THREED */

#ifdef IHOP_WRITE_OUT
    WRITE(msgbuf,'(A)') Title
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgbuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE( msgBuf, '(A,G11.4,A)' )'Frequency ', IHOP_freq, 'Hz'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    ! *** Top Boundary ***
    Bdry%Top%HS%Opt = IHOP_topopt
    Bdry%Top%HS%Depth = 0 !initiate to dummy value
    CALL ReadTopOpt( Bdry%Top%HS%Opt, Bdry%Top%HS%BC, AttenUnit, FileRoot, &
                     myThid )

    CALL TopBot( IHOP_freq, AttenUnit, Bdry%Top%HS, myThid )

    ! *** Ocean SSP ***
    IF ( IHOP_depth.NE.0 ) THEN
        Bdry%Bot%HS%Depth = IHOP_depth
    ELSE
        ! Extend by 5 wavelengths
        Bdry%Bot%HS%Depth = rkSign*rF( Nr+1 ) + 5*1500.0/IHOP_freq 
    END IF
    x = [ 0.0 _d 0, Bdry%Bot%HS%Depth ]   ! tells SSP Depth to read to

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A,F10.2,A)' ) 'Depth = ',Bdry%Bot%HS%Depth,' m'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)') 'Top options: ', Bdry%Top%HS%Opt
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    !   Only do I/O if in the main thread
    _END_MASTER(myThid)

    CALL EvaluateSSP( x, c, cimag, gradc, crr, crz, czz, rho, IHOP_freq, 'INI', myThid )

    !   Only do I/O if in the main thread
    _BEGIN_MASTER(myThid)

    Bdry%Top%HS%Depth = SSP%z( 1 )   ! first SSP point is top depth

    ! *** Bottom Boundary ***
    ! bottom depth should perhaps be set the same way?
    Bdry%Bot%HS%Opt = IHOP_botopt 
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)') 'Bottom options: ', Bdry%Bot%HS%Opt
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    SELECT CASE ( Bdry%Bot%HS%Opt( 2 : 2 ) )
    CASE ( '~', '*' )
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(A)') '    Bathymetry file selected'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE( ' ' )
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'READENVIHOP ReadEnvironment: ', & 
                             'Unknown bottom option letter in second position'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadEnvironment'
    END SELECT

    Bdry%Bot%HS%BC = Bdry%Bot%HS%Opt( 1 : 1 )
    CALL TopBot( IHOP_freq, AttenUnit, Bdry%Bot%HS, myThid )

    ! *** source and receiver locations ***
    CALL ReadSxSy( myThid ) ! Read source/receiver x-y coordinates
    ZMin = SNGL( Bdry%Top%HS%Depth )
    ZMax = SNGL( Bdry%Bot%HS%Depth )

    Pos%NSz = IHOP_nsd
    Pos%NRz = IHOP_nrd
    CALL AllocateSR( Pos%NSz, Pos%Sz, IHOP_sd )
    CALL AllocateSR( Pos%NRz, Pos%Rz, IHOP_rd )
    CALL ReadSzRz( ZMin, ZMax, myThid )

    Pos%NRr = IHOP_nrr
    CALL AllocateSR( Pos%NRr, Pos%Rr, IHOP_rr )
    CALL ReadRcvrRanges( myThid )

#ifdef IHOP_THREED
    CALL ReadRcvrBearings( myThid )
#endif /* IHOP_THREED */
    CALL ReadfreqVec( IHOP_freq,  Bdry%Top%HS%Opt( 6:6 ), myThid )

    ! *** run type ***
    Beam%RunType = IHOP_runopt
    CALL ReadRunType( Beam%RunType, PlotType, myThid )

    Depth = Zmax - Zmin   ! water depth
    CALL ReadRayElevationAngles( IHOP_freq, Depth, Bdry%Top%HS%Opt, &
        Beam%RunType, myThid )
#ifdef IHOP_THREED
    CALL ReadRayBearingAngles( IHOP_freq, Bdry%Top%HS%Opt, Beam%RunType, myThid )
#endif /* IHOP_THREED */

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)')'___________________________________________________', &
                        '________'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    ! Limits for tracing beams
#ifdef IHOP_THREED
    WRITE(msgBuf,'(2A)') 'READENVIHOP ReadEnvironment: ', & 
                         '3D not supported in ihop'
    CALL PRINT_ERROR( msgBuf,myThid )
    STOP 'ABNORMAL END: S/R ReadEnvironment'
    !READ(  ENVFile, * ) Beam%deltas, Beam%Box%x, Beam%Box%y, Beam%Box%z
    Beam%Box%x = 1000.0 * Beam%Box%x   ! convert km to m
    Beam%Box%y = 1000.0 * Beam%Box%y   ! convert km to m

    ! Automatic step size selection
    IF ( Beam%deltas == 0.0 ) Beam%deltas = &
        ( Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth ) / 10.0   
#ifdef IHOP_WRITE_OUT
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
#endif /* IHOP_WRITE_OUT */
#else /* not IHOP_THREED */
    ! Step size in meters [m]
    Beam%deltas = IHOP_step
    IF ( Beam%deltas == 0.0 ) THEN ! Automatic step size option
        Beam%deltas = ( Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth ) / 10.0   
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(A)') 
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,A,G11.4,A)') &
            ' Step length,', ' deltas = ', Beam%deltas, ' m (automatic step)'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    ELSE
        WRITE(msgBuf,'(A,A,G11.4,A)') &
            ' Step length,', ' deltas = ', Beam%deltas, ' m'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    END IF

    ! Domain size
    IF ( IHOP_zbox.NE.0 ) THEN
        Beam%Box%z = IHOP_zbox
    ELSE
        Beam%Box%z = Bdry%Bot%HS%Depth ! in m
    END IF
    IF ( IHOP_rbox.NE.0 ) THEN
        Beam%Box%r = IHOP_rbox
    ELSE
        ! Extend beam box by a single step size forward
        Beam%Box%r = ihop_rr(nrd) + Beam%deltas/1000. ! in [km]
    END IF
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A,G11.4,A)') &
        ' Maximum ray range, Box%R = ', Beam%Box%R,' km'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A,G11.4,A)') &
        ' Maximum ray depth, Box%Z = ', Beam%Box%Z,' m'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    Beam%Box%r = Beam%Box%r*1000.   ! convert km to m
#endif /* IHOP_THREED */

    ! *** Beam characteristics ***
       Beam%Type( 4 : 4 ) = Beam%RunType( 7 : 7 )   ! selects beam shift option
          
#ifdef IHOP_WRITE_OUT
       SELECT CASE ( Beam%Type( 4 : 4 ) )
       CASE ( 'S' )
          WRITE(msgBuf,'(A)') ' Beam shift in effect'
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
       CASE DEFAULT
          WRITE(msgBuf,'(A)') ' No beam shift in effect'
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
       END SELECT
#endif /* IHOP_WRITE_OUT */

       ! no worry about the beam type if this is a ray trace run
       IF ( Beam%RunType( 1:1 ) /= 'R' .OR. Beam%RunType( 1:1 ) /= 'E' ) THEN 

       ! Beam%Type( 1 : 1 ) is
       !   'G' or '^' Geometric hat beams in Cartesian coordinates
       !   'g' Geometric hat beams in ray-centered coordinates
       !   'B' Geometric Gaussian beams in Cartesian coordinates
       !   'b' Geometric Gaussian beams in ray-centered coordinates
       !   'S' Simple Gaussian beams
       !   'C' Cerveny Gaussian beams in Cartesian coordinates
       !   'R' Cerveny Gaussian beams in Ray-centered coordinates
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

       Beam%Type( 1 : 1 ) = Beam%RunType( 2 : 2 )
       SELECT CASE ( Beam%Type( 1 : 1 ) )
! geometric hat beams, geometric Gaussian beams, or simple Gaussian beams
       CASE ( 'G', 'g' , '^', 'B', 'b', 'S' )   
! Cerveny Gaussian Beams; read extra lines to specify the beam options
       CASE ( 'R', 'C' )   
          !READ(  ENVFile, * ) Beam%Type( 2 : 3 ), Beam%epsMultiplier, Beam%rLoop
#ifdef IHOP_WRITE_OUT
          WRITE(msgBuf,'(A)') 
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(A)') 
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(2A)') 'Type of beam = ', Beam%Type( 1:1 )
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

          SELECT CASE ( Beam%Type( 3 : 3 ) )
          CASE ( 'D' )
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(A)') 'Curvature doubling invoked'
            CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
          CASE ( 'Z' )
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(A)') 'Curvature zeroing invoked'
            CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
          CASE ( 'S' )
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(A)') 'Standard curvature condition'
            CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
          CASE DEFAULT
#ifdef IHOP_WRITE_OUT
                WRITE(msgBuf,'(2A)') 'READENVIHOP ReadEnvironment: ', & 
                                     'Unknown curvature condition'
                CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
                STOP 'ABNORMAL END: S/R ReadEnvironment'
          END SELECT

#ifdef IHOP_WRITE_OUT
          WRITE(msgBuf,'(A,F1.2)') 'UNUSED epsMultiplier', Beam%epsMultiplier
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(A,F10.2)') 'Range for choosing beam width', Beam%rLoop
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

          ! Images, windows
          WRITE(msgBuf,'(A)') 
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(A,I10)') 'Number of images, Nimage  = ', Beam%Nimage
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(A,I10)') 'Beam windowing parameter  = ', Beam%iBeamWindow
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(A)') 'Component                 = ', Beam%Component
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
       CASE DEFAULT
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A)') 'READENVIHOP ReadEnvironment: ', & 
                                'Unknown beam type (second letter of run type)'
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadEnvironment'
       END SELECT
    END IF

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    !-- Only do I/O if in the main thread
    _END_MASTER( myThid)
  RETURN
  END !SUBROUTINE ReadEnvironment

  !**********************************************************************!

  SUBROUTINE ReadTopOpt( TopOpt, BC, AttenUnit, FileRoot, myThid )

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    CHARACTER (LEN= 6), INTENT( OUT ) :: TopOpt
    CHARACTER (LEN= 1), INTENT( OUT ) :: BC         ! Boundary condition type
    CHARACTER (LEN= 2), INTENT( OUT ) :: AttenUnit
    CHARACTER (LEN=80), INTENT( IN  ) :: FileRoot

    TopOpt = IHOP_topopt
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    SSP%Type  = TopOpt( 1 : 1 )
    BC        = TopOpt( 2 : 2 )
    AttenUnit = TopOpt( 3 : 4 )
    SSP%AttenUnit = AttenUnit

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
        WRITE(msgBuf,'(2A)') 'READENVIHOP ReadTopOpt: ', & 
                             'Unknown option for SSP approximation'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    ! Attenuation options

    SELECT CASE ( AttenUnit( 1 : 1 ) )
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
        WRITE(msgBuf,'(2A)') 'READENVIHOP ReadTopOpt: ', & 
                             'Unknown attenuation units'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    ! optional addition of volume attenuation using standard formulas

    SELECT CASE ( AttenUnit( 2 : 2 ) )
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
        WRITE(msgBuf,'(2A)') 'READENVIHOP ReadTopOpt: ', & 
                             'Unknown top option letter in fourth position'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    SELECT CASE ( TopOpt( 5 : 5 ) )
    CASE ( '~', '*' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Altimetry file selected'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( '-', '_', ' ' )
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'READENVIHOP ReadTopOpt: ', & 
                             'Unknown top option letter in fifth position'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    SELECT CASE ( TopOpt( 6 : 6 ) )
    CASE ( 'I' )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') '    Development options enabled'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    CASE ( ' ' )
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'READENVIHOP ReadTopOpt: ', & 
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
    CHARACTER (LEN=10), INTENT( OUT ) :: PlotType

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    SELECT CASE ( RunType( 1 : 1 ) )
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
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'READENVIHOP ReadRunType: ', & 
            'Unknown RunType selected'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadRunType'
    END SELECT

    SELECT CASE ( RunType( 2 : 2 ) )
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
       RunType( 2 : 2 ) = 'G'
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Geometric hat beams in Cartesian coordinates'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    END SELECT

    SELECT CASE ( RunType( 4 : 4 ) )
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
       RunType( 4 : 4 ) = 'R'
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Point source (cylindrical coordinates)'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    END SELECT

    SELECT CASE ( RunType( 5 : 5 ) )
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
            WRITE(msgBuf,'(2A)') 'READENVIHOP ReadRunType: ', & 
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
       RunType( 5 : 5 ) = 'R'
       PlotType = 'rectilin  '
    END SELECT

    SELECT CASE ( RunType( 6 : 6 ) )
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
       RunType( 6 : 6 ) = '2'
    END SELECT

  RETURN
  END !SUBROUTINE ReadRunType

  !**********************************************************************!

  SUBROUTINE TopBot( freq, AttenUnit, HS, myThid )

    ! Handles top and bottom boundary conditions

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN    ) :: freq  ! frequency
    CHARACTER (LEN=2), INTENT( IN    ) :: AttenUnit
    TYPE ( HSInfo ),   INTENT( INOUT ) :: HS
    REAL (KIND=_RL90) :: Mz, vr, alpha2_f     ! values related to grain size

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
        WRITE(msgBuf,'(2A)') 'READENVIHOP TopBot: ', & 
                             'Unknown boundary condition type'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R TopBot'
    END SELECT

    ! ****** Read in BC parameters depending on particular choice ******

    HS%cp  = 0.0
    HS%cs  = 0.0
    HS%rho = 0.0

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
       WRITE(msgBuf,'(A)') 
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(A)') &
        '   z (m)     alphaR (m/s)   betaR  rho (g/cm^3)  alphaI     betaI'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )' ) &
            zTemp, alphaR, betaR, rhoR, alphaI, betaI
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
       ! dummy parameters for a layer with a general power law for attenuation
       ! these are not in play because the AttenUnit for this is not allowed yet
       !freq0         = freq
       betaPowerLaw  = 1.0
       ft            = 1000.0

       HS%cp  = CRCI( zTemp, alphaR, alphaI, freq, freq, AttenUnit, &
                      betaPowerLaw, ft, myThid )
       HS%cs  = CRCI( zTemp, betaR,  betaI,  freq, freq, AttenUnit, &
                      betaPowerLaw, ft, myThid )

       HS%rho = rhoR
    CASE ( 'G' )            ! *** Grain size (formulas from UW-APL HF Handbook)
       ! These formulas are from the UW-APL Handbook
       ! The code is taken from older Matlab and is unnecesarily verbose
       ! vr   is the sound speed ratio
       ! rhor is the density ratio
       !READ(  ENVFile, *    ) zTemp, Mz
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'( F10.2, 3X, F10.2 )' ) zTemp, Mz
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

       IF ( Mz >= -1 .AND. Mz < 1 ) THEN
          vr   = 0.002709 * Mz**2 - 0.056452 * Mz + 1.2778
          rhor = 0.007797 * Mz**2 - 0.17057  * Mz + 2.3139
       ELSE IF ( Mz >= 1 .AND. Mz < 5.3 ) THEN
          vr   = -0.0014881 * Mz**3 + 0.0213937 * Mz**2 - 0.1382798 * Mz &
               + 1.3425
          rhor = -0.0165406 * Mz**3 + 0.2290201 * Mz**2 - 1.1069031 * Mz &
               + 3.0455
       ELSE
          vr   = -0.0024324 * Mz + 1.0019
          rhor = -0.0012973 * Mz + 1.1565
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
       alphaI = alpha2_f * ( vr / 1000 ) * 1500.0 * log( 10.0 ) / ( 40.0 * PI )

       HS%cp  = CRCI( zTemp, alphaR, alphaI, freq, freq, 'L ', &
                      betaPowerLaw, ft, myThid )
       HS%cs  = 0.0
       HS%rho = rhoR
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A,2F10.2,3X,A,F10.2,3X,A,F10.2)') &
           'Converted sound speed =', HS%cp, 'density = ', rhor, &
           'loss parm = ', alphaI
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    END SELECT

  RETURN
  END !SUBROUTINE TopBot

  ! **********************************************************************!

  SUBROUTINE OpenOutputFiles( FileRoot )
    ! Write appropriate header information

    USE angle_mod,  only: Angles
    USE srPos_mod,  only: Pos

    ! == Routine Arguments ==

    ! == Local variables ==
    CHARACTER (LEN=80), INTENT( IN ) :: FileRoot
    REAL               :: atten
    CHARACTER (LEN=10) :: PlotType

    SELECT CASE ( Beam%RunType( 1 : 1 ) )
    CASE ( 'R', 'E' )   ! Ray trace or Eigenrays
#ifdef IHOP_WRITE_OUT
       OPEN ( FILE = TRIM( FileRoot ) // '.ray', UNIT = RAYFile, &
              FORM = 'FORMATTED' )
       WRITE( RAYFile, * ) '''', Title( 1 : 50 ), ''''
       WRITE( RAYFile, * ) IHOP_freq
       WRITE( RAYFile, * ) Pos%NSx, Pos%NSy, Pos%NSz
       WRITE( RAYFile, * ) Angles%Nalpha, Angles%Nbeta
       WRITE( RAYFile, * ) Bdry%Top%HS%Depth
       WRITE( RAYFile, * ) Bdry%Bot%HS%Depth

#ifdef IHOP_THREED
       WRITE( RAYFile, * ) '''xyz'''
#else /* IHOP_THREED */
       WRITE( RAYFile, * ) '''rz'''
#endif /* IHOP_THREED */
#endif /* IHOP_WRITE_OUT */

    CASE ( 'A' )        ! arrival file in ascii format
#ifdef IHOP_WRITE_OUT
       OPEN ( FILE = TRIM( FileRoot ) // '.arr', UNIT = ARRFile, &
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

       ! IEsco22: add to arrivals output
       OPEN ( FILE = TRIM( FileRoot ) // '.ray', UNIT = RAYFile, &
              FORM = 'FORMATTED' )
       WRITE( RAYFile, * ) '''', Title( 1 : 50 ), ''''
       WRITE( RAYFile, * ) IHOP_freq
       WRITE( RAYFile, * ) Pos%NSx, Pos%NSy, Pos%NSz
       WRITE( RAYFile, * ) Angles%Nalpha, Angles%Nbeta
       WRITE( RAYFile, * ) Bdry%Top%HS%Depth
       WRITE( RAYFile, * ) Bdry%Bot%HS%Depth

# ifdef IHOP_THREED
       WRITE( RAYFile, * ) '''xyz'''
# else /* IHOP_THREED */
       WRITE( RAYFile, * ) '''rz'''
# endif /* IHOP_THREED */

       OPEN ( FILE = TRIM( FileRoot ) // '.delay', UNIT = DELFile, &
              FORM = 'FORMATTED' )
       WRITE( DELFile, * ) '''', Title( 1 : 50 ), ''''
       WRITE( DELFile, * ) IHOP_freq
       WRITE( DELFile, * ) Pos%NSx, Pos%NSy, Pos%NSz
       WRITE( DELFile, * ) Angles%Nalpha, Angles%Nbeta
       WRITE( DELFile, * ) Bdry%Top%HS%Depth
       WRITE( DELFile, * ) Bdry%Bot%HS%Depth

#ifdef IHOP_THREED
       WRITE( DELFile, * ) '''xyz'''
# else /* IHOP_THREED */
       WRITE( DELFile, * ) '''rz'''
# endif /* IHOP_THREED */
#endif /* IHOP_WRITE_OUT */
    CASE ( 'a' )        ! arrival file in binary format
#ifdef IHOP_WRITE_OUT
       OPEN ( FILE = TRIM( FileRoot ) // '.arr', UNIT = ARRFile, &
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

       CALL WriteSHDHeader( TRIM( FileRoot ) // '.shd', Title, REAL( IHOP_freq ), &
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

  SUBROUTINE AllocateSR( Nx, x_out, x_in )

    ! Allocate and populate Pos structure from data.ihop

    INTEGER,          INTENT( IN  ) :: Nx    
    REAL(KIND=_RL90), INTENT( IN  ) :: x_in(:)
    REAL(KIND=_RL90), ALLOCATABLE, INTENT( OUT ) :: x_out(:)
    INTEGER                         :: i

    IF ( ALLOCATED(x_out) ) DEALLOCATE(x_out)
    ALLOCATE( x_out(MAX(3, Nx)) )
    x_out(3) = -999.9

    DO i = 1, Nx
        x_out(i) = x_in(i)
    END DO

  RETURN
  END !SUBROUTINE AllocateSR

  !**********************************************************************!

END MODULE readEnviHop
