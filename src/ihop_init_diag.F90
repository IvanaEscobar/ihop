#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE IHOP_INIT_DIAG
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>

  ! mbp 12/2018, based on much older subroutine

  USE ssp_mod,   only: SSP
  USE bdry_mod,  only: Bdry, HSInfo
  USE atten_mod, only: CRCI

! ! USES
  implicit none
!  == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef ALLOW_CAL
#include "cal.h"
#endif /* ALLOW_CAL */

  PRIVATE

!   == Public Interfaces ==
!=======================================================================
  public    initPRTFile, OpenOutputFiles, resetMemory
!=======================================================================

! INPUT/OUTPUT PARAMETERS:

! == External Functions ==
  INTEGER  ILNBLNK
  EXTERNAL ILNBLNK

CONTAINS
  SUBROUTINE initPRTFile( myTime, myIter, myThid )
    USE ihop_mod,  only: PRTFile, Beam
    USE angle_mod, only: Angles
    USE srpos_mod, only: WriteSxSy, WriteSzRz, WriteRcvrRanges, WriteFreqVec

    ! I/O routine for acoustic fixed inputS

  ! == Routine Arguments ==
  ! myTime  :: Current time in simulation
  ! myIter  :: Current time-step number
  ! myThid  :: my Thread Id number
  ! msgBuf  :: Used to build messages for printing.
    _RL, INTENT(IN)     ::  myTime
    INTEGER, INTENT(IN) ::  myIter, myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  ! == Local Variables ==
    INTEGER, PARAMETER :: Number_to_Echo = 10
!    REAL (KIND=_RL90),  PARAMETER :: c0 = 1500.0
!    REAL (KIND=_RL90)  :: x(2), c, cimag, gradc(2), crz, czz, rho, Depth
    REAL (KIND=_RL90)  :: ranges
!    CHARACTER (LEN=10) :: PlotType

    ! *** ihop info to PRTFile ***
    CALL openPRTFile( myTime, myIter, myThid )

#ifdef IHOP_WRITE_OUT
    !   Only do I/O in the main thread
    _BEGIN_MASTER(myThid)

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
      CALL WriteRunType( Beam%RunType, myThid )

      CALL WriteTopOpt( myThid )

      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F10.2,A)' ) 'Depth = ',Bdry%Bot%HS%Depth,' m'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)') 'Top options: ', Bdry%Top%HS%Opt
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

      CALL WriteTopBot( Bdry%Top%HS, myThid )

      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)') 'Bottom options: ', Bdry%Bot%HS%Opt
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

      SELECT CASE ( Bdry%Bot%HS%Opt( 2 : 2 ) )
        CASE ( '~', '*' )
          WRITE(msgBuf,'(A)') '    Bathymetry file selected'
          CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        CASE( ' ' )
        CASE DEFAULT
      END SELECT

      CALL WriteTopBot( Bdry%Bot%HS, myThid )


      CALL WriteSxSy( myThid )
      CALL WriteSzRz( Bdry%Top%HS%Depth, Bdry%Bot%HS%Depth, myThid )
      CALL WriteRcvrRanges( myThid )
#ifdef IHOP_THREED
      CALL WriteRcvrBearings( myThid )
#endif
      CALL WriteFreqVec( Bdry%Top%HS%Opt( 6:6 ), myThid )


      WRITE(msgBuf,'(2A)')'_____________________________________________', &
                          '______________'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,I10)') 'Number of beams in elevation   = ', &
                              Angles%Nalpha
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      IF ( Angles%iSingle_alpha > 0 ) THEN
        WRITE(msgBuf,'(A,I10)') 'Trace only beam number ', &
                                Angles%iSingle_alpha
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      END IF

      WRITE(msgBuf,'(A)') 'Beam take-off angles (degrees)'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

      IF ( Angles%Nalpha >= 1 ) THEN
        WRITE(msgBuf,'(10F12.3)') &
            Angles%alpha( 1:MIN(Angles%Nalpha,Number_to_Echo) )
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      END IF
      IF ( Angles%Nalpha > Number_to_Echo ) THEN
        WRITE(msgBuf,'(A,F12.6)') ' ... ', Angles%alpha( Angles%Nalpha )
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      END IF

      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)')'______________________________________________', &
                          '_____________'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,G11.4,A)') &
          ' Step length, deltas = ', Beam%deltas, ' m'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

#ifdef IHOP_THREED
      WRITE(msgBuf,'(A,G11.4,A)') &
          ' Maximum ray x-range, Box%X = ', Beam%Box%X,' m'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,G11.4,A)') &
          ' Maximum ray y-range, Box%Y = ', Beam%Box%Y,' m'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,G11.4,A)') &
          ' Maximum ray z-range, Box%Z = ', Beam%Box%Z,' m'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#else /* not IHOP_THREED */
      ranges = Beam%Box%R / 1000.0
      WRITE(msgBuf,'(A,G11.4,A)') &
          ' Maximum ray range, Box%R = ', ranges,' km'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,G11.4,A)') &
          ' Maximum ray depth, Box%Z = ', Beam%Box%Z,' m'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* not IHOP_THREED */

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

    ENDIF ! adjoint run check

    !   Only do I/O in the main thread
    _END_MASTER(myThid)
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE initPRTFile

  !**********************************************************************!

  SUBROUTINE WriteTopOpt( myThid )
    USE ihop_mod,  only: PRTFile
    USE ssp_mod,   only: SSP
    USE atten_mod, only: T, Salinity, pH, z_bar, iBio, NBioLayers, bio

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    CHARACTER (LEN= 1) :: BC ! Boundary condition type

    BC = IHOP_TopOpt( 2:2 )

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) RETURN

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 'Interior options: '
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    ! SSP approximation options
    SELECT CASE ( SSP%Type )
      CASE ( 'N' )
        WRITE(msgBuf,'(A)') '    N2-linear approximation to SSP'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'C' )
        WRITE(msgBuf,'(A)') '    C-linear approximation to SSP'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'P' )
        WRITE(msgBuf,'(A)') '    PCHIP approximation to SSP'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'S' )
        WRITE(msgBuf,'(A)') '    Spline approximation to SSP'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'Q' )
        WRITE(msgBuf,'(A)') '    Quad approximation to SSP'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'A' )
        WRITE(msgBuf,'(A)') '    Analytic SSP option'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE DEFAULT
    END SELECT

    ! Attenuation options
    SELECT CASE ( SSP%AttenUnit( 1:1 ) )
      CASE ( 'N' )
        WRITE(msgBuf,'(A)') '    Attenuation units: nepers/m'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'F' )
        WRITE(msgBuf,'(A)') '    Attenuation units: dB/mkHz'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'M' )
        WRITE(msgBuf,'(A)') '    Attenuation units: dB/m'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'W' )
        WRITE(msgBuf,'(A)') '    Attenuation units: dB/wavelength'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'Q' )
        WRITE(msgBuf,'(A)') '    Attenuation units: Q'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'L' )
        WRITE(msgBuf,'(A)') '    Attenuation units: Loss parameter'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE DEFAULT
    END SELECT

    ! optional addition of volume attenuation using standard formulas
    SELECT CASE ( SSP%AttenUnit( 2:2 ) )
      CASE ( 'T' )
        WRITE(msgBuf,'(A)') '    THORP volume attenuation added'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'F' )
        WRITE(msgBuf,'(A)') '    Francois-Garrison volume attenuation added'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE( PRTFile, &
          "( ' T = ', G11.4, 'degrees   S = ', G11.4, ' psu   pH = ', G11.4, ' z_bar = ', G11.4, ' m' )" ) &
          T, Salinity, pH, z_bar
      CASE ( 'B' )
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
      CASE ( ' ' )
      CASE DEFAULT
    END SELECT

    SELECT CASE ( IHOP_TopOpt( 5:5 ) )
      CASE ( '~', '*' )
        WRITE(msgBuf,'(A)') '    Altimetry file selected'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( '-', '_', ' ' )
      CASE DEFAULT
    END SELECT

    SELECT CASE ( IHOP_TopOpt( 6:6 ) )
      CASE ( 'I' )
        WRITE(msgBuf,'(A)') '    Development options enabled'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( ' ' )
      CASE DEFAULT
    END SELECT
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE WriteTopOpt

  !**********************************************************************!

  SUBROUTINE WriteRunType( RunType, myThid )
    USE ihop_mod,  only: PRTFile

  ! Write the RunType variable to .prt file

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    CHARACTER*(7), INTENT( IN ) :: RunType

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) RETURN

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    SELECT CASE ( RunType( 1:1 ) )
      CASE ( 'R' )
        WRITE(msgBuf,'(A)') 'Ray trace run'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'E' )
        WRITE(msgBuf,'(A)') 'Eigenray trace run'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'I' )
        WRITE(msgBuf,'(A)') 'Incoherent TL calculation'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'S' )
        WRITE(msgBuf,'(A)') 'Semi-coherent TL calculation'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'C' )
        WRITE(msgBuf,'(A)') 'Coherent TL calculation'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'A' )
        WRITE(msgBuf,'(A)') 'Arrivals calculation, ASCII  file output'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'a' )
        WRITE(msgBuf,'(A)') 'Arrivals calculation, binary file output'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'e' )
        WRITE(msgBuf,'(A)') 'Eigenrays + Arrivals run, ASCII file output'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE DEFAULT
    END SELECT

    SELECT CASE ( RunType( 2:2 ) )
      CASE ( 'C' )
        WRITE(msgBuf,'(A)') 'Cartesian beams'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'R' )
        WRITE(msgBuf,'(A)') 'Ray centered beams'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'S' )
        WRITE(msgBuf,'(A)') 'Simple gaussian beams'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'b' )
        WRITE(msgBuf,'(2A)') 'Geometric gaussian beams in ray-centered ', &
                             'coordinates'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'B' )
        WRITE(msgBuf,'(A)') 'Geometric gaussian beams in Cartesian coordinates'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'g' )
        WRITE(msgBuf,'(A)') 'Geometric hat beams in ray-centered coordinates'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE DEFAULT
    END SELECT

    SELECT CASE ( RunType( 4:4 ) )
      CASE ( 'R' )
        WRITE(msgBuf,'(A)') 'Point source (cylindrical coordinates)'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'X' )
        WRITE(msgBuf,'(A)') 'Line source (Cartesian coordinates)'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE DEFAULT
    END SELECT

    SELECT CASE ( RunType( 5:5 ) )
      CASE ( 'R' )
        WRITE(msgBuf,'(2A)') 'Rectilinear receiver grid: Receivers at', &
                             ' ( Rr( ir ), Rz( ir ) ) )'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'I' )
        WRITE(msgBuf,'(A)') 'Irregular grid: Receivers at Rr( : ) x Rz( : )'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE DEFAULT
    END SELECT

    SELECT CASE ( RunType( 6:6 ) )
      CASE ( '2' )
        WRITE(msgBuf,'(2A)') 'N x 2D calculation (neglects ', &
                             'horizontal refraction)'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( '3' )
        WRITE(msgBuf,'(A)') '3D calculation'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE DEFAULT
    END SELECT

    WRITE(msgBuf,'(2A)')'__________________________________________', &
                        '_________________'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE WriteRunType

  !**********************************************************************!

  SUBROUTINE WriteTopBot( HS, myThid )
    USE ihop_mod,  only: PRTFile
    USE ssp_mod, only: rhoR, alphaR, betaR, alphaI, betaI

    ! Handles top and bottom boundary conditions

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    TYPE ( HSInfo ),   INTENT( IN ) :: HS
    REAL (KIND=_RL90) :: Mz ! values related to grain size

#ifdef IHOP_WRITE_OUT
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
    ! Echo to PRTFile user's choice of boundary condition

    SELECT CASE ( HS%BC )
      CASE ( 'V' )
        WRITE(msgBuf,'(A)') '    Surface modeled as a VACUUM'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'R' )
        WRITE(msgBuf,'(A)') '    Perfectly RIGID'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'A' )
        WRITE(msgBuf,'(A)') '    ACOUSTO-ELASTIC half-space'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)')
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') &
         '   z [m]     alphaR [m/s]   betaR  rho [g/cm^3]  alphaI     betaI'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )' ) &
             HS%Depth, alphaR, betaR, rhoR, alphaI, betaI
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'G' )  ! *** Grain size (formulas from UW-APL HF Handbook)
        WRITE(msgBuf,'(A)') '    Grain size to define half-space'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)')
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'( F10.2, 3X, F10.2 )' ) HS%Depth, Mz
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,2F10.2,3X,A,F10.2,3X,A,F10.2)') &
            'Converted sound speed =', HS%cp, 'density = ', rhoR, &
            'loss parm = ', alphaI
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'F' )
        WRITE(msgBuf,'(A)') '    FILE used for reflection loss'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'W' )
        WRITE(msgBuf,'(A)') '    Writing an IRC file'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE ( 'P' )
        WRITE(msgBuf,'(A)') '    reading PRECALCULATED IRC'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      CASE DEFAULT
    END SELECT

    END IF ! if in adjoint mode
#endif /* IHOP_WRITE_OUT */
  RETURN
  END !SUBROUTINE WriteTopBot

  ! **********************************************************************!

  SUBROUTINE OpenOutputFiles( fName, myTime, myIter, myThid )
    USE ihop_mod,  only: RAYFile, DELFile, ARRFile, SHDFile, Title, Beam
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

       CALL WriteSHDHeader( TRIM( fullName ) // '.shd', Title, &
                            REAL( IHOP_freq ), atten, PlotType )
    END SELECT

  RETURN
  END !SUBROUTINE OpenOutputFiles

  !**********************************************************************!

  SUBROUTINE WriteSHDHeader( FileName, Title, freq0, Atten, PlotType )

    USE srPos_mod,  only: Pos, Nfreq, freqVec
    USE ihop_mod,   only: SHDFile

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
    USE ihop_mod,   only: SHDFile

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
    USE ihop_mod, only: PRTFile, Title

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     myTime :: Current time in simulation
  !     myIter :: Current time-step number
  !     msgBuf :: Used to build messages for printing.
    _RL, INTENT(IN)     ::  myTime
    INTEGER, INTENT(IN) ::  myIter, myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf


  !     == Local Arguments ==
    INTEGER                     :: iostat, ierr
  ! For MPI writing: inspo from eeboot_minimal.F
    CHARACTER*(MAX_LEN_FNAM)    :: fNam
    CHARACTER*(6)               :: fmtStr
    INTEGER                     :: mpiRC, IL
#ifdef ALLOW_CAL
    INTEGER :: mydate(4)
#endif

    ! Open the print file: template from eeboot_minimal.F
#ifdef IHOP_WRITE_OUT
    IF ( .NOT.usingMPI ) THEN
        WRITE(myProcessStr, '(I10.10)') myIter
        IL=ILNBLNK( myProcessStr )
        WRITE(fNam,'(4A)') TRIM(IHOP_fileroot),'.',myProcessStr(1:IL),'.prt'
        IF ( IHOP_dumpfreq.GE.0 ) &
         OPEN( PRTFile, FILE = fNam, STATUS = 'UNKNOWN', IOSTAT = iostat )
#ifdef ALLOW_USE_MPI
    ELSE ! using MPI
        CALL MPI_COMM_RANK( MPI_COMM_MODEL, mpiMyId, mpiRC )
        myProcId = mpiMyId
        IL = MAX(4,1+INT(LOG10(DFLOAT(nPx*nPy))))
        WRITE(fmtStr,'(2(A,I1),A)') '(I',IL,'.',IL,')'
        WRITE(myProcessStr,fmtStr) myProcId
        IL = ILNBLNK( myProcessStr )
        mpiPidIo = myProcId
        pidIO    = mpiPidIo

        IF( mpiPidIo.EQ.myProcId ) THEN
#  ifdef SINGLE_DISK_IO
         IF( myProcId.eq.0 ) THEN
#  endif
            IF (myIter.GE.0) THEN
                WRITE(fNam,'(4A,I10.10,A)') &
                    TRIM(IHOP_fileroot),'.',myProcessStr(1:IL),'.',myIter,'.prt'
            ELSE
                WRITE(fNam,'(4A)') &
                    TRIM(IHOP_fileroot),'.',myProcessStr(1:IL),'.prt'
            ENDIF

            IF ( IHOP_dumpfreq .GE. 0) &
             OPEN(PRTFile, FILE=fNam, STATUS='UNKNOWN', IOSTAT=iostat )
            IF ( iostat /= 0 ) THEN
                WRITE(*,*) 'ihop: IHOP_fileroot not recognized, ', &
                    TRIM(IHOP_fileroot)
                WRITE(msgBuf,'(A)') 'IHOP_INIT: Unable to recognize file'
                CALL PRINT_ERROR( msgBuf, myThid )
                STOP 'ABNORMAL END: S/R IHOP_INIT'
            END IF
#  ifdef SINGLE_DISK_IO
         END IF
#  endif
        END IF
# endif /* ALLOW_USE_MPI */
    END IF
#endif /* IHOP_WRITE_OUT */

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
    WRITE(msgBuf,'(2A)') 'IHOP_INIT_DIAG openPRTFile: ', &
                         '3D not supported in ihop'
    CALL PRINT_ERROR( msgBuf,myThid )
    STOP 'ABNORMAL END: S/R openPRTFile'
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
    USE angle_mod,  only: Angles
    USE arr_mod,    only: Narr, Arr, U
    USE ihop_mod,   only: ray2D, MaxN, iStep

    ! From ihop
    IF (ALLOCATED(Pos%theta))   DEALLOCATE(Pos%theta)
    IF (ALLOCATED(U))           DEALLOCATE(U)
    IF (ALLOCATED(Arr))         DEALLOCATE(Arr)
    IF (ALLOCATED(NArr))        DEALLOCATE(NArr)
    ! From ssp_mod
    IF (useSSPFile) THEN
        ! don't reset values, they've been read in from a file -_-
    ELSE
        SSP%cmat   = 1.0
        SSP%czmat  = 1.0
#ifdef IHOP_THREED
        SSP%cmat3  = 1.0
        SSP%czmat3 = 1.0
#endif /* IHOP_THREED */
    ENDIF
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

END MODULE IHOP_INIT_DIAG
