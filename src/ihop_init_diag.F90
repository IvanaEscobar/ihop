#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: ihop_init_diag
MODULE IHOP_INIT_DIAG
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION: 
!   Initialization of iHOP diagnostic files.

! !USES:
  USE bdry_mod,  only: Bdry, HSInfo
  USE atten_mod, only: CRCI
  IMPLICIT NONE
! == Global variables ==
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

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC initPRTFile, openOutputFiles, resetMemory
!=======================================================================

! == Module variables == None
! == External Functions ==
  INTEGER  ILNBLNK
  EXTERNAL ILNBLNK
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R initPRTFile
! S/R openPRTFile
! S/R WriteTopOpt
! S/R WriteRunType
! S/R WriteBdry
! S/R openOutputFiles
! S/R WriteSHDHeader
! S/R WriteSHDField
! S/R AllocatePos
! S/R resetMemory
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: initPRTFile
! !INTERFACE:
  SUBROUTINE initPRTFile( myTime, myIter, myThid )
! !DESCRIPTION:
!   Initializes the iHOP print file.

! !USES:
  USE ihop_mod,  only: PRTFile, Beam
  USE angle_mod, only: Angles
  USE srpos_mod, only: WriteSxSy, WriteSzRz, WriteRcvrRanges, WriteFreqVec

! !INPUT PARAMETERS:
! myTime  :: Current time in simulation
! myIter  :: Current time-step number
! myThid  :: my thread ID
  _RL,     INTENT( IN ) ::  myTime
  INTEGER, INTENT( IN ) ::  myIter, myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf          :: Informational/error message buffer
! Number_to_Echo  :: Number of angles to echo in the print file
! ranges          :: Range values for printing
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
  INTEGER, PARAMETER :: Number_to_Echo = 10
  REAL (KIND=_RL90)  :: ranges

  ! In adjoint mode we do not write output besides on the first run
  IF ( IHOP_dumpfreq.LT.0 ) RETURN

  ! *** ihop info to PRTFile ***
  CALL openPRTFile( myTime, myIter, myThid )

#ifdef IHOP_WRITE_OUT
  !   Only do I/O in the main thread
  _BEGIN_MASTER(myThid)

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

  CALL WriteBdry( Bdry%Top%HS, myThid )

  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(2A)') 'Bottom options: ', Bdry%Bot%HS%Opt
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  SELECT CASE ( Bdry%Bot%HS%Opt( 2:2 ) )
  CASE ( '~', '*' )
    WRITE(msgBuf,'(A)') '    Bathymetry file selected'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  CASE( ' ' )
  CASE DEFAULT
  END SELECT

  CALL WriteBdry( Bdry%Bot%HS, myThid )

  CALL WriteSxSy( myThid )
  CALL WriteSzRz( myThid )
  CALL WriteRcvrRanges( myThid )
# ifdef IHOP_THREED
  CALL WriteRcvrBearings( myThid )
# endif
  CALL WriteFreqVec( Bdry%Top%HS%Opt( 6:6 ), myThid )


  WRITE(msgBuf,'(A)') &
    '___________________________________________________________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A,I5)') 'Number of beams in elevation   = ', &
   Angles%nAlpha
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  IF ( Angles%iSingle_alpha.GT.0 ) THEN
    WRITE(msgBuf,'(A,I10)') 'Trace only beam number ', &
      Angles%iSingle_alpha
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  ENDIF

  WRITE(msgBuf,'(A)') 'Beam take-off angles (degrees)'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  IF ( Angles%nAlpha.GE.1 ) THEN
    WRITE(msgBuf,'(10F12.3)') &
      Angles%adeg( 1:MIN(Angles%nAlpha,Number_to_Echo) )
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  ENDIF

  IF ( Angles%nAlpha.GT.Number_to_Echo ) THEN
    WRITE(msgBuf,'(A,F12.6)') ' ... ', Angles%adeg( Angles%nAlpha )
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  ENDIF

  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') &
    '___________________________________________________________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A,G11.4,A)') &
    ' Step length, deltas = ', Beam%deltas, ' m'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

# ifdef IHOP_THREED
  WRITE(msgBuf,'(A,G11.6,A)') &
    ' Maximum ray x-range, Box%X = ', Beam%Box%X,' m'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A,G11.6,A)') &
    ' Maximum ray y-range, Box%Y = ', Beam%Box%Y,' m'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A,G11.6,A)') &
    ' Maximum ray z-range, Box%Z = ', Beam%Box%Z,' m'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
# else /* not IHOP_THREED */
  ranges = Beam%Box%R / 1000.0
  WRITE(msgBuf,'(A,G11.6,A)') &
    ' Maximum ray range, Box%R = ', ranges,' km'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A,G11.6,A)') &
    ' Maximum ray depth, Box%Z = ', Beam%Box%Z,' m'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
# endif /* not IHOP_THREED */

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

  !   Only do I/O in the main thread
  _END_MASTER(myThid)
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE initPRTFile

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: openPRTFile
! !INTERFACE:
  SUBROUTINE openPRTFile ( myTime, myIter, myThid )
! !DESCRIPTION:
!   Opens the iHOP print file for writing.

! !USES:
  USE ihop_mod, only: PRTFile, Title

! !INPUT PARAMETERS:
! myTime :: Current time in simulation
! myIter :: Current time-step number
! myThid :: my thread ID
  _RL,     INTENT( IN ) ::  myTime
  INTEGER, INTENT( IN ) ::  myIter, myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! iostat :: I/O status variable
! fNam   :: Name of the print file
! fmtStr :: Format string for writing
! mpiRC  :: MPI return code
! IL     :: Length of the process string
! mydate :: Date in pkg/cal format
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
  INTEGER                  :: iostat
! For MPI writing: inspo from eeboot_minimal.F
  CHARACTER*(MAX_LEN_FNAM)    :: fNam
  CHARACTER*(6)               :: fmtStr
  INTEGER                     :: mpiRC, IL
#ifdef ALLOW_CAL
  INTEGER :: mydate(4)
#endif

  ! In adjoint mode we do not write output besides on the first run
  IF ( IHOP_dumpfreq.LT.0 ) RETURN

  ! Open the print file: template from eeboot_minimal.F
#ifdef IHOP_WRITE_OUT
  IF ( .NOT.usingMPI ) THEN
    WRITE(myProcessStr, '(I10.10)') myIter
    IL=ILNBLNK( myProcessStr )
    WRITE(fNam,'(4A)') TRIM(IHOP_fileroot),'.',myProcessStr(1:IL),'.prt'
    OPEN( PRTFile, FILE=fNam, STATUS='UNKNOWN', IOSTAT=iostat )
    IF ( iostat.NE.0 ) THEN
      WRITE(*,*) 'ihop: IHOP_fileroot not recognized, ', &
        TRIM(IHOP_fileroot)
      WRITE(msgBuf,'(A)') 'IHOP_INIT: Unable to recognize file'
      CALL PRINT_ERROR( msgBuf, myThid )
      STOP 'ABNORMAL END: S/R IHOP_INIT'
    ENDIF

# ifdef ALLOW_USE_MPI
  ELSE ! using MPI
    IF ( myIter.GE.0 ) THEN
      WRITE(fNam,'(2A,I10.10,A)') &
        TRIM(IHOP_fileroot),'.',myIter,'.prt'
    ELSE
      WRITE(msgBuf,'(A,I)') 'IHOP_INIT: myIter is ', myIter
      CALL PRINT_ERROR( msgBuf, myThid )
      STOP 'ABNORMAL END: S/R IHOP_INIT'
    ENDIF

    OPEN(PRTFile, FILE=fNam, STATUS='UNKNOWN', IOSTAT=iostat )
    IF ( iostat.NE.0 ) THEN
      WRITE(*,*) 'ihop: IHOP_fileroot not recognized, ', &
        TRIM(IHOP_fileroot)
      WRITE(msgBuf,'(A)') 'IHOP_INIT: Unable to recognize file'
      CALL PRINT_ERROR( msgBuf, myThid )
      STOP 'ABNORMAL END: S/R IHOP_INIT'
    ENDIF
# endif /* ALLOW_USE_MPI */

  ENDIF ! IF ( .NOT.usingMPI )
#endif /* IHOP_WRITE_OUT */

  !   Only do I/O in the main thread
  _BARRIER
  _BEGIN_MASTER(myThid)

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
  Title( 1:9   ) = 'iHOP - '
  Title( 10:80 ) = IHOP_title
#endif /* IHOP_THREED */

#ifdef IHOP_WRITE_OUT
  WRITE(msgbuf,'(A)') TRIM(Title) ! , ACHAR(10)
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') &
    '___________________________________________________________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgbuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE( msgBuf, '(A,I10,A,F20.2,A)') 'GCM iter ', myIter, &
    ' at time = ', myTime,' [sec]'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
# ifdef ALLOW_CAL
  CALL CAL_GETDATE( myIter,myTime,mydate,myThid )
  WRITE (msgBuf,'(A,I8,I6,I3,I4)') 'GCM cal date ', mydate(1), &
    mydate(2), mydate(3), mydate(4)
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
# endif /* ALLOW_CAL */
  WRITE(msgbuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE( msgBuf, '(A,F11.4,A)' ) 'Frequency ', IHOP_freq, ' [Hz]'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') &
  '___________________________________________________________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

  !   Only do I/O in the main thread
  _END_MASTER(myThid)

  RETURN
  END !SUBROUTINE openPRTFile

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteTopOpt
! !INTERFACE:
  SUBROUTINE WriteTopOpt( myThid )
! !DESCRIPTION:
!  Write the top options for iHOP model.

! !USES:
  USE ihop_mod,  only: PRTFile
  USE ssp_mod,   only: Grid
  USE atten_mod, only: T, Salinity, pH, z_bar, iBio, NBioLayers, bio

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! BC     :: Boundary condition type
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  CHARACTER*(1) :: BC
!EOP

  BC = IHOP_TopOpt( 2:2 )

  ! In adjoint mode we do not write output besides on the first run
  IF (IHOP_dumpfreq.LT.0) RETURN

#ifdef IHOP_WRITE_OUT
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') 'Interior options: '
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  ! SSP_mod approximation options
  SELECT CASE ( Grid%Type )
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
  SELECT CASE ( Grid%AttenUnit( 1:1 ) )
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
  SELECT CASE ( Grid%AttenUnit( 2:2 ) )
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
    ENDDO

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

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteRunType
! !INTERFACE:
  SUBROUTINE WriteRunType( RunType, myThid )
! !DESCRIPTION:
! Write the RunType variable to .prt file

! !USES:
  USE ihop_mod, only: PRTFile

! !INPUT PARAMETERS:
! RunType :: String describing the run type
! myThid  :: my thread ID
  CHARACTER*(7), INTENT( IN ) :: RunType
  INTEGER,       INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
!EOP

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
    WRITE(msgBuf,'(A)') 'Irregular grid: Receivers at Rr(:) x Rz(:)'
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

  WRITE(msgBuf,'(A)') &
    '___________________________________________________________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE WriteRunType

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteBdry
! !INTERFACE:
  SUBROUTINE WriteBdry( HS, myThid )
! !DESCRIPTION:
!   Writes the boundary condition information to the print file.

! !USES:
  USE ihop_mod,only: PRTFile
  USE ssp_mod, only: rhoR, alphaR, betaR, alphaI, betaI

! !INPUT PARAMETERS:
! HS     :: Half-space information
! myThid :: my thread ID
  TYPE ( HSInfo ), INTENT( IN ) :: HS
  INTEGER,         INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! Mz     :: Grain size value
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  REAL (KIND=_RL90) :: Mz ! values related to grain size
!EOP

  ! In adjoint mode we do not write output besides on the first run
  IF (IHOP_dumpfreq.LT.0) RETURN

#ifdef IHOP_WRITE_OUT
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
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE writeBdry

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: openOutputFiles
! !INTERFACE:
  SUBROUTINE openOutputFiles( fName, myTime, myIter, myThid )
! !DESCRIPTION:
!   Opens the unique output files for iHOP model.

! !USES:
  USE ihop_mod,  only: RAYFile, DELFile, ARRFile, SHDFile, &
                       Title, Beam
  USE angle_mod, only: Angles
  USE srPos_mod, only: Pos

! !INPUT PARAMETERS:
! fName   :: Name of the output file
! myTime  :: Current time in simulation
! myIter  :: Current time-step number
! myThid  :: my thread ID
  CHARACTER*(MAX_LEN_FNAM), INTENT( IN ) :: fName
  _RL,                      INTENT( IN ) :: myTime
  INTEGER,                  INTENT( IN ) :: myIter, myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! fullName  :: Full name of the output file
! PlotType  :: Type of plot to generate
! fmtstr    :: Format string for output
! IL        :: Length of input filename
! atten     :: Attenuation value
! isOpen    :: Logical variable to check if file is open
  CHARACTER*(MAX_LEN_FNAM) :: fullName
  CHARACTER*(10) :: PlotType
  CHARACTER*(20) :: fmtstr
  INTEGER        :: IL
  REAL           :: atten
  LOGICAL        :: isOpen
!EOP

  ! In adjoint mode we do not write output besides on the first run
  IF (IHOP_dumpfreq.LT.0) RETURN

  IF (myIter.GE.0) THEN
    ! add time step to filename
    IL=ILNBLNK( fName )
    WRITE(fullName, '(2A,I10.10)') fName( 1:IL ), '.', myIter
  ELSE
    fullName = fName
  ENDIF

  SELECT CASE ( Beam%RunType( 1:1 ) )
  CASE ( 'R', 'E' )   ! Ray trace or Eigenrays
#ifdef IHOP_WRITE_OUT
    INQUIRE(UNIT=RAYFile, OPENED=isOpen)
    IF (isOpen) CLOSE(RAYFile)
    OPEN ( FILE=TRIM( fullName ) // '.ray', UNIT=RAYFile, &
          FORM='FORMATTED' )
    WRITE( RAYFile, '(3A)'    ) '''', Title( 1:50 ), ''''
    WRITE( RAYFile, '(F10.3)' ) IHOP_freq
    WRITE( RAYFile, '(3I6)'   ) Pos%nSX, Pos%nSY, Pos%nSZ
    WRITE( RAYFile, '(I5)'    ) Angles%nAlpha
    WRITE( RAYFile, '(F10.4)' ) Bdry%Top%HS%Depth
    WRITE( RAYFile, '(F10.4)' ) Bdry%Bot%HS%Depth
#ifdef IHOP_THREED
    WRITE( RAYFile, '(2I5)' ) Angles%nAlpha, Angles%nBeta
    WRITE( RAYFile, '(A)'   ) '''xyz'''
#else /* IHOP_THREED */
    WRITE( RAYFile, '(A)'  ) '''rz'''
#endif /* IHOP_THREED */

    FLUSH( RAYFile )
#endif /* IHOP_WRITE_OUT */

    IF (writeDelay) THEN
      INQUIRE(UNIT=DELFile, OPENED=isOpen)
      IF (isOpen) CLOSE(DELFile)
      OPEN ( FILE=TRIM( fullName ) // '.delay', UNIT=DELFile, &
            FORM='FORMATTED' )
      WRITE( DELFile, '(3A)'    ) '''', Title( 1:50 ), ''''
      WRITE( DELFile, '(F10.3)' ) IHOP_freq
      WRITE( DELFile, '(3I6)'   ) Pos%nSX, Pos%nSY, Pos%nSZ
      WRITE( DELFile, '(I5)'    ) Angles%nAlpha
      WRITE( DELFile, '(F10.4)' ) Bdry%Top%HS%Depth
      WRITE( DELFile, '(F10.4)' ) Bdry%Bot%HS%Depth
#ifdef IHOP_THREED
      WRITE( DELFile, '(2I5)' ) Angles%nAlpha, Angles%nBeta
      WRITE( DELFile, '(A)'   ) '''xyz'''
# else /* IHOP_THREED */
      WRITE( DELFile, '(A)'  ) '''rz'''
# endif /* IHOP_THREED */

      FLUSH( DELFile )

    ENDIF ! IF (writeDelay)

  CASE ( 'e' ) ! eigenrays + arrival file in ascii format
#ifdef IHOP_WRITE_OUT
    INQUIRE(UNIT=ARRFile, OPENED=isOpen)
    IF (isOpen) CLOSE(ARRFile)
    OPEN ( FILE=TRIM( fullName ) // '.arr', UNIT=ARRFile, &
          FORM='FORMATTED' )
# ifdef IHOP_THREED
    WRITE( ARRFile, '(A)' ) '''3D'''
# else /* IHOP_THREED */
    WRITE( ARRFile, '(A)' ) '''2D'''
# endif /* IHOP_THREED */
    WRITE( ARRFile, '(F10.4)' ) IHOP_freq

    ! write source locations
# ifdef IHOP_THREED
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSX, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSX, Pos%SX( 1:Pos%nSX )
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSY, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSY, Pos%SY( 1:Pos%nSY )
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSZ, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSZ, Pos%SZ( 1:Pos%nSZ )
# else /* IHOP_THREED */
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSZ, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSZ, Pos%SZ( 1:Pos%nSZ )
# endif /* IHOP_THREED */

    ! write receiver locations
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nRZ, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nRZ, Pos%RZ( 1:Pos%nRZ )
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nRR, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nRR, Pos%RR( 1:Pos%nRR )
# ifdef IHOP_THREED
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nTheta, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nTheta, Pos%theta( 1:Pos%nTheta )
# endif /* IHOP_THREED */

    FLUSH( ARRFile )

    ! IEsco22: add erays to arrivals output
    INQUIRE(UNIT=RAYFile, OPENED=isOpen)
    IF (isOpen) CLOSE(RAYFile)
    OPEN ( FILE=TRIM( fullName ) // '.ray', UNIT=RAYFile, &
          FORM='FORMATTED' )
    WRITE( RAYFile, '(3A)'    ) '''', Title( 1:50 ), ''''
    WRITE( RAYFile, '(F10.3)' ) IHOP_freq
    WRITE( RAYFile, '(3I6)'   ) Pos%nSX, Pos%nSY, Pos%nSZ
    WRITE( RAYFile, '(I5)'    ) Angles%nAlpha
    WRITE( RAYFile, '(F10.4)' ) Bdry%Top%HS%Depth
    WRITE( RAYFile, '(F10.4)' ) Bdry%Bot%HS%Depth
# ifdef IHOP_THREED
    WRITE( RAYFile, '(2I5)' ) Angles%nAlpha, Angles%nBeta
    WRITE( RAYFile, '(A)'   ) '''xyz'''
# else /* IHOP_THREED */
    WRITE( RAYFile, '(A)'  ) '''rz'''
# endif /* IHOP_THREED */

    FLUSH( RAYFile )

    IF (writeDelay) THEN
      INQUIRE(UNIT=DELFile, OPENED=isOpen)
      IF (isOpen) CLOSE(DELFile)
      OPEN ( FILE=TRIM( fullName ) // '.delay', UNIT=DELFile, &
            FORM='FORMATTED' )
      WRITE( DELFile, '(3A)'    ) '''', Title( 1:50 ), ''''
      WRITE( DELFile, '(F10.3)' ) IHOP_freq
      WRITE( DELFile, '(3I6)'   ) Pos%nSX, Pos%nSY, Pos%nSZ
      WRITE( DELFile, '(I5)'    ) Angles%nAlpha
      WRITE( DELFile, '(F10.4)' ) Bdry%Top%HS%Depth
      WRITE( DELFile, '(F10.4)' ) Bdry%Bot%HS%Depth
#ifdef IHOP_THREED
      WRITE( DELFile, '(2I5)' ) Angles%nAlpha, Angles%nBeta
      WRITE( DELFile, '(A)'   ) '''xyz'''
# else /* IHOP_THREED */
      WRITE( DELFile, '(A)'  ) '''rz'''
# endif /* IHOP_THREED */

      FLUSH( DELFile )

    ENDIF ! IF (writeDelay)
#endif /* IHOP_WRITE_OUT */

  CASE ( 'A' )        ! arrival file in ascii format
#ifdef IHOP_WRITE_OUT
    INQUIRE(UNIT=ARRFile, OPENED=isOpen)
    IF (isOpen) CLOSE(ARRFile)
    OPEN ( FILE=TRIM( fullName ) // '.arr', UNIT=ARRFile, &
          FORM='FORMATTED' )

# ifdef IHOP_THREED
    WRITE( ARRFile, '(A)' ) '''3D'''
# else /* IHOP_THREED */
    WRITE( ARRFile, '(A)' ) '''2D'''
# endif /* IHOP_THREED */
    WRITE( ARRFile, '(F10.4)' ) IHOP_freq

      ! write source locations
# ifdef IHOP_THREED
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSX, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSX, Pos%SX( 1:Pos%nSX )
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSY, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSY, Pos%SY( 1:Pos%nSY )
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSZ, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSZ, Pos%SZ( 1:Pos%nSZ )
# else /* IHOP_THREED */
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSZ, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSZ, Pos%SZ( 1:Pos%nSZ )
# endif /* IHOP_THREED */

    ! write receiver locations
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nRZ, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nRZ, Pos%RZ( 1:Pos%nRZ )
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nRR, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nRR, Pos%RR( 1:Pos%nRR )
# ifdef IHOP_THREED
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nTheta, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nTheta, Pos%theta( 1:Pos%nTheta )
# endif /* IHOP_THREED */

    FLUSH( ARRFile )
#endif /* IHOP_WRITE_OUT */

  CASE ( 'a' )        ! arrival file in binary format
#ifdef IHOP_WRITE_OUT
    INQUIRE(UNIT=ARRFile, OPENED=isOpen)
    IF (isOpen) CLOSE(ARRFile)
    OPEN ( FILE=TRIM( fullName ) // '.arr', UNIT=ARRFile, &
          FORM='UNFORMATTED' )

# ifdef IHOP_THREED
    WRITE( ARRFile, '(A)' ) '''3D'''
# else /* IHOP_THREED */
    WRITE( ARRFile, '(A)' ) '''2D'''
# endif /* IHOP_THREED */
    WRITE( ARRFile, '(F10.4)' ) IHOP_freq

    ! write source locations
# ifdef IHOP_THREED
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSX, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSX, Pos%SX( 1:Pos%nSX )
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSY, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSY, Pos%SY( 1:Pos%nSY )
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSZ, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSZ, Pos%SZ( 1:Pos%nSZ )
# else /* IHOP_THREED */
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nSZ, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nSZ, Pos%SZ( 1:Pos%nSZ )
# endif /* IHOP_THREED */

    ! write receiver locations
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nRZ, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nRZ, Pos%RZ( 1:Pos%nRZ )
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nRR, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nRR, Pos%RR( 1:Pos%nRR )
# ifdef IHOP_THREED
    WRITE( fmtstr, '(A,I0,A)' ) '(I4,', Pos%nTheta, 'G16.10)'
    WRITE( ARRFile, fmtstr    ) Pos%nTheta, Pos%theta( 1:Pos%nTheta )
# endif /* IHOP_THREED */

    FLUSH( ARRFile )
#endif /* IHOP_WRITE_OUT */

  CASE DEFAULT
    atten = 0.0

    ! following to set PlotType has alread been done in READIN if that was
    ! used for input
    SELECT CASE ( Beam%RunType( 5:5 ) )
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
  END !SUBROUTINE openOutputFiles

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteSHDHeader
! !INTERFACE:
  SUBROUTINE WriteSHDHeader( FileName, Title, freq0, Atten, PlotType )
! !DESCRIPTION:
!   Writes the header of a SHD file to disk. 

! !USES:
  USE srPos_mod,  only: Pos, Nfreq, freqVec
  USE ihop_mod,   only: SHDFile

! !INPUT PARAMETERS:
! FileName :: Name of the SHD file to write
! Title    :: Title of the SHD file
! freq0    :: Nominal frequency for the SHD file
! Atten    :: Stabilizing attenuation for wavenumber integration
! PlotType :: Type of plot to generate
  REAL,      INTENT( IN ) :: freq0, Atten
  CHARACTER, INTENT( IN ) :: FileName*( * )
  CHARACTER, INTENT( IN ) :: Title*( * )
  CHARACTER, INTENT( IN ) :: PlotType*( 10 )
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! LRecl :: Length of the record in the SHD file
  INTEGER :: LRecl
!EOP

  ! receiver bearing angles
  IF ( .NOT. ALLOCATED( Pos%theta ) ) THEN
    ALLOCATE( Pos%theta( 1 ) )
    Pos%theta( 1 ) = 0   ! dummy bearing angle
    Pos%nTheta     = 1
  END IF

  ! source x-coordinates
  IF ( .NOT. ALLOCATED( Pos%SX ) ) THEN
    ALLOCATE( Pos%SX( 1 ) )
    Pos%SX( 1 ) = 0      ! dummy x-coordinate
    Pos%nSX     = 1
  END IF

  ! source y-coordinates
  IF ( .NOT. ALLOCATED( Pos%SY ) ) THEN
    ALLOCATE( Pos%SY( 1 ) )
    Pos%SY( 1 ) = 0      ! dummy y-coordinate
    Pos%nSY     = 1
  END IF

  IF ( PlotType( 1 : 2 ) /= 'TL' ) THEN
    ! MAX( 41, ... ) below because Title is already 40 words (or 80 bytes)
    ! words/record (nRR doubled for complex pressure storage)
    LRecl = MAX( 41, 2 * Nfreq, Pos%nTheta, Pos%nSX, Pos%nSY, Pos%nSZ, &
                Pos%nRZ, 2 * Pos%nRR )

    OPEN ( FILE=FileName, UNIT=SHDFile, STATUS='REPLACE', &
          ACCESS='DIRECT', RECL=4 * LRecl, FORM='UNFORMATTED' )
    WRITE( SHDFile, REC=1  ) LRecl, Title( 1 : 80 )
    WRITE( SHDFile, REC=2  ) PlotType
    WRITE( SHDFile, REC=3  ) Nfreq, Pos%nTheta, Pos%nSX, Pos%nSY, &
      Pos%nSZ, Pos%nRZ, Pos%nRR, freq0, atten
    WRITE( SHDFile, REC=4  ) freqVec(   1 : Nfreq )
    WRITE( SHDFile, REC=5  ) Pos%theta( 1 : Pos%nTheta )

    WRITE( SHDFile, REC=6  ) Pos%SX( 1 : Pos%nSX )
    WRITE( SHDFile, REC=7  ) Pos%SY( 1 : Pos%nSY )
    WRITE( SHDFile, REC=8  ) Pos%SZ( 1 : Pos%nSZ )

    WRITE( SHDFile, REC=9  ) Pos%RZ( 1 : Pos%nRZ )
    WRITE( SHDFile, REC=10 ) Pos%RR( 1 : Pos%nRR )

  ELSE   ! compressed format for TL from FIELD3D
    ! words/record (NR doubled for complex pressure storage)
    LRecl = MAX( 41, 2 * Nfreq, Pos%nTheta, Pos%nSZ, Pos%nRZ, 2 * Pos%nRR )

    OPEN ( FILE=FileName, UNIT=SHDFile, STATUS='REPLACE', &
          ACCESS='DIRECT', RECL=4 * LRecl, FORM='UNFORMATTED')
    WRITE( SHDFile, REC=1  ) LRecl, Title( 1 : 80 )
    WRITE( SHDFile, REC=2  ) PlotType
    WRITE( SHDFile, REC=3  ) Nfreq, Pos%nTheta, Pos%nSX, Pos%nSY, &
      Pos%nSZ, Pos%nRZ, Pos%nRR, freq0, atten
    WRITE( SHDFile, REC=4  ) freqVec(   1 : Nfreq )
    WRITE( SHDFile, REC=5  ) Pos%theta( 1 : Pos%nTheta )

    WRITE( SHDFile, REC=6  ) Pos%SX( 1 ), Pos%SX( Pos%nSX )
    WRITE( SHDFile, REC=7  ) Pos%SY( 1 ), Pos%SY( Pos%nSY )
    WRITE( SHDFile, REC=8  ) Pos%SZ( 1 : Pos%nSZ )

    WRITE( SHDFile, REC=9  ) Pos%RZ( 1 : Pos%nRZ )
    WRITE( SHDFile, REC=10 ) Pos%RR( 1 : Pos%nRR )

  END IF

  RETURN
  END !SUBROUTINE WriteSHDHeader

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteSHDField
! !INTERFACE:
  SUBROUTINE WriteSHDField( P, nRZ, nRR, iRec )
! !DESCRIPTION:
!   Writes the pressure field to the SHD file.

! !USES:
  USE ihop_mod, only: SHDFile
    
! !INPUT PARAMETERS:
! P     :: Pressure field to write
! nRZ   :: Number of receiver depths
! nRR   :: Number of ranges
! iRec  :: Last record read in the SHD file (to be incremented)
  INTEGER, INTENT( IN )    :: nRZ, nRR
  COMPLEX, INTENT( IN )    :: P( nRZ, nRR )
  INTEGER, INTENT( INOUT ) :: iRec
! !OUTPUT PARAMETERS: iRec

! !LOCAL VARIABLES:
! iRz   :: Loop index for receiver depths
  INTEGER :: iRz
!EOP

  DO iRz = 1, nRZ
    iRec = iRec + 1
    WRITE( SHDFile, REC=iRec ) P( iRz, : )
  END DO

  RETURN
  END !SUBROUTINE WriteSHDField

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: AllocatePos
! !INTERFACE:
  SUBROUTINE AllocatePos( Nx, x_out, x_in )
! !DESCRIPTION:
!   Allocates and populates the Pos structure from the input data.

! !USES: None

! !INPUT PARAMETERS:
! Nx     :: Number of elements in the input array
! x_out  :: Output array of positions (allocated)
! x_in   :: Input array of positions
  INTEGER,          INTENT( IN  ) :: Nx
  REAL(KIND=_RL90), ALLOCATABLE, INTENT( OUT ) :: x_out(:)
  REAL(KIND=_RL90), INTENT( IN  ) :: x_in(:)
! !OUTPUT PARAMETERS: x_out

! !LOCAL VARIABLES:
! i :: Loop index
  INTEGER :: i
!EOP

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

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: resetMemory
! !INTERFACE:
  SUBROUTINE resetMemory()
! !DESCRIPTION:
!   Resets memory of iHOP model to default values for each new time step.

! !USES:
  USE ssp_mod,    only: SSP
  USE arr_mod,    only: nArrival, Arr, U
  USE ihop_mod,   only: ray2D, nMax, iStep

  ! From arr_mod.f90
  U                         = 0.
  nArrival(:,:)             = 0
  Arr(:,:,:)%NTopBnc        = -1
  Arr(:,:,:)%NBotBnc        = -1
  Arr(:,:,:)%SrcDeclAngle   = -999.
  Arr(:,:,:)%RcvrDeclAngle  = -999.
  Arr(:,:,:)%A              = -999.
  Arr(:,:,:)%Phase          = -999.
  Arr(:,:,:)%delay          = -999.
  Arr(:,:,:)%delayR         = -999.

  ! From ssp_mod.f90
  IF (useSSPFile) THEN
    ! don't reset values, they've been read in from a file -_-
  ELSE
    SSP%cmat   = -99.0
    SSP%czmat  = -99.0
#ifdef IHOP_THREED
    SSP%cmat3  = -99.0
    SSP%czmat3 = -99.0
#endif /* IHOP_THREED */
  ENDIF ! IF (useSSPFile)

  ! From ihop_mod.f90
  DO iStep = 1,nMax
    ray2D(iStep)%x = [zeroRL, zeroRL]
    ray2D(iStep)%t = [zeroRL, zeroRL]
    ray2D(iStep)%p = [zeroRL, zeroRL]
    ray2D(iStep)%q = [zeroRL, zeroRL]
    ray2D(iStep)%c = zeroRL
    ray2D(iStep)%Amp = zeroRL
    ray2D(iStep)%Phase = zeroRL
    ray2D(iStep)%tau = (zeroRL, zeroRL)
  END DO

  RETURN
  END !SUBROUTINE resetMemory

END !MODULE IHOP_INIT_DIAG
