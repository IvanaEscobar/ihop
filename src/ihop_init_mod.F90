#include "IHOP_OPTIONS.h"
MODULE ihop_init_mod

IMPLICIT NONE
PRIVATE

public init_fixed_env

CONTAINS
  SUBROUTINE INIT_FIXED_ENV ( myThid )
  ! Initiate fixed variable for ihop time series. Note: NO IHOP_THREED here
  ! ===========================================================================

  ! USES
    USE bdry_mod,  only: Bdry, HSInfo, initATI, initBTY
    USE srpos_mod, only: Pos, ReadSXSY, ReadSZRz, ReadRcvrRanges, ReadFreqVec
    USE ssp_mod,   only: SSP, initSSP, alphar
    USE ihop_mod,  only: Beam, rxyz, nRz_per_range
    USE angle_mod, only: Angles, ReadRayElevationAngles
#ifdef IHOP_THREED
    USE angle_mod, only: ReadRayBearingAngles
#endif /* IHOP_THREED */
    USE refcoef,   only: ReadReflectionCoefficient
    USE beampat,   only: SBPFlag, readPat
    USE arr_mod,   only: initArr

  ! ===========================================================================
  !     == Global Variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "IHOP.h"
#endif

  ! ===========================================================================
  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

!!  ! ===========================================================================
!!  !     == Local Variables ==
    REAL (KIND=_RL90), PARAMETER :: c0 = 1500.0
    CHARACTER (LEN=2)  :: AttenUnit
    CHARACTER (LEN=10) :: PlotType
    REAL (KIND=_RL90)  :: Depth
!!
!!  ! ===========================================================================

  !IESCO24: some notes while I noodle
    ! Use data.ihop, set time series invariant parameters. These are fixed
    ! parameters that do not depend on which time step you run ihop in.
    ! Primarily, the parameters are related to the acoustic grid:
    ! - From initenvihop.F90:initEnv
    !   - Bdry%Top, Bdry%Bot,
    !     SSP%AttenUnit,Type,nR,nZ,z,SSP%Seg%R,
    !     Pos%SX,SY,nSZ,nRZ,SZ,Rz,Ws,ISZ,Wr,Irz,nRR,Rr,Delta_r,
    !     Beam%RunType,Deltas,Nimage,iBeamWindow,Component,Multiplier,rloop,
    !     Beam%Box%R,Box%Z,Type,
    !     Angles%nAlpha,alpha,
    ! - From bdry_mod.F90:initATI
    !   - Top%Natipts,x,
    ! - From bdry_mod.F90:initBTY
    !   - Bot%Natipts,x,
    ! This subroutine will set parameters that shouldn't need to be modified
    ! throughout the MITgcm model run

    ! === Set local parameters ===
    AttenUnit = ''
    PlotType = ''
    Depth = -1.

    ! === Set nonallocatable derived type components from other modules ===
    Bdry%Bot%HS = HSInfo(0.,0.,0.,0., 0.,0. , (0.,0.),(0.,0.), '', '' )
    Bdry%Top%HS = HSInfo(0.,0.,0.,0., 0.,0. , (0.,0.),(0.,0.), '', '' )

    Pos%nSX = -1
    Pos%nSY = -1
    Pos%nSZ = -1
    Pos%nRZ = -1
    Pos%nRR = -1
    Pos%nTheta = 1
    Pos%Delta_r = -999.
    Pos%Delta_theta = -999.

    SSP%nPts = -1.
    SSP%nR = -1.
    SSP%nX = -1.
    SSP%nY = -1.
    SSP%nZ = -1.
    SSP%Z = -1.
    SSP%rho = -1.
    SSP%c = -1.
    SSP%Type = ''
    SSP%AttenUnit = ''

    Beam%nBeams = -1
    Beam%nImage = -1
    Beam%nSteps = -1
    Beam%iBeamWindow = -1
    Beam%deltas = -1.
    Beam%epsMultiplier = 1.
    Beam%rLoop = -1.
    Beam%Component = ''
    Beam%Type = 'G S '
    Beam%RunType = ''
    Beam%Box = rxyz(0.,0.,0.,0.)

    Angles%nAlpha = 0
    Angles%iSingle_alpha = 0
    Angles%Dalpha = -1.

    ! === From initenvihop.f90:initEnv ===
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

    Bdry%Bot%HS%BC = Bdry%Bot%HS%Opt( 1:1 )
    CALL TopBot( AttenUnit, Bdry%Bot%HS, myThid )

    SELECT CASE ( Bdry%Bot%HS%Opt( 2:2 ) )
      CASE( '~', '*', ' ' )
      CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'IHOP_INIT_MOD init_fixed_env: ',&
            'Unknown Bdry%Bot%HS%Opt(2)'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R init_fixed_env'
    END SELECT


    ! *** SSP parameters ***
    CALL initSSP( myThid )

    ! set Bdry%Top%HS%Depth from first SSP%Z
    Bdry%Top%HS%Depth = SSP%Z(1)
    ! set water column depth
    Depth = Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth


    ! *** Source locations ***
    CALL ReadSXSY( myThid ) ! Read source/receiver x-y coordinates

    Pos%nSZ = IHOP_nsd
    Pos%nRZ = IHOP_nrd

    CALL AllocatePos( Pos%nSZ, Pos%SZ, IHOP_sd, myThid )
    CALL AllocatePos( Pos%nRZ, Pos%RZ, IHOP_rd, myThid )
    CALL ReadSZRz( Bdry%Top%HS%Depth, Bdry%Bot%HS%Depth, myThid )


    ! *** Receiver locations ***
    Pos%nRR = IHOP_nRR
    CALL AllocatePos( Pos%nRR, Pos%RR, IHOP_rr, myThid )
    CALL ReadRcvrRanges( myThid )
    ! set dummy for receiver bearing
    CALL AllocatePos( Pos%nTheta, Pos%theta, IHOP_rr( 1:1 ), myThid )
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
    ! Step size in meters [m]
    Beam%deltas = IHOP_step

    ! Automatic step size option
    IF ( Beam%deltas == 0.0 ) THEN
        Beam%deltas = ( Depth ) / 10.
    END IF

    ! Domain size
    Beam%Box%Z = Bdry%Bot%HS%Depth ! in [m]
    ! Extend beam box by a single step size forward
    Beam%Box%R = IHOP_rr(nrd)*1000. + Beam%deltas ! in [m]


    ! *** Beam characteristics ***
    Beam%Type( 4:4 ) = Beam%RunType( 7:7 )   ! selects beam shift option

    ! don't worry about the beam type if this is a ray trace run
    ! IESCO23: using 'e' requires Beam%Type to be set
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
          WRITE(msgBuf,'(2A)') 'IHOP_INIT_MOD init_fixed_env: ', &
              'Unknown beam type (second letter of Beam%Type)'
          CALL PRINT_ERROR( msgBuf,myThid )
          !   Only do I/O in the main thread
          _END_MASTER(myThid)
#endif /* IHOP_WRITE_OUT */
          STOP 'ABNORMAL END: S/R init_fixed_env'
      END SELECT

    END IF ! Beam%RunType( 1:1 ) /= 'R' ...


! =================
! ================= from IHOP.F90: S/R IHOP_MAIN
! =================
    ! AlTImetry: OPTIONAL, default is no ATIFile
    CALL initATI( Bdry%Top%HS%Opt( 5:5 ), Bdry%Top%HS%Depth, myThid )
    ! BaThYmetry: OPTIONAL, default is BTYFile
    CALL initBTY( Bdry%Bot%HS%Opt( 2:2 ), Bdry%Bot%HS%Depth, myThid )
    ! (top and bottom): OPTIONAL
    CALL readReflectionCoefficient( myThid )


    ! Source Beam Pattern: OPTIONAL, default is omni source pattern
    SBPFlag = Beam%RunType( 3:3 )
    CALL readPat( myThid )


! Allocate arrival and U variables on all MPI processes
    SELECT CASE ( Beam%RunType( 5:5 ) )
      CASE ( 'I' )
        nRz_per_range = 1         ! irregular grid
      CASE DEFAULT
        nRz_per_range = Pos%nRZ   ! rectilinear grid
    END SELECT
!
    CALL initArr( myThid )

! =================
! ================= from IHOP.F90: S/R IHOP_MAIN
! =================

  RETURN
  END !SUBROUTINE IHOP_INIT_FIXED_ENV

  ! **********************************************************************!
  SUBROUTINE ReadTopOpt( BC, AttenUnit, myThid )
    USE atten_mod, only: T, Salinity, pH, z_bar, iBio, NBioLayers, bio
    USE ssp_mod, only: SSP

  ! ===========================================================================
  !     == Global Variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "IHOP.h"
#endif

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

    ! SSP approximation options
    SELECT CASE ( SSP%Type )
      CASE ( 'N','C','P','S','Q','A' )
      CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INIT_FIXED_ENV ReadTopOpt: ', &
                             'Unknown option for SSP approximation'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    ! Attenuation options
    SELECT CASE ( AttenUnit( 1:1 ) )
      CASE ( 'N','F','M','W','Q','L' )
      CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INIT_FIXED_ENV ReadTopOpt: ', &
                             'Unknown attenuation units'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    ! optional addition of volume attenuation using standard formulas
    SELECT CASE ( AttenUnit( 2:2 ) )
      CASE ( 'T','F','B',' ' )
      CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INIT_FIXED_ENV ReadTopOpt: ', &
                             'Unknown top option letter in fourth position'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    SELECT CASE ( IHOP_TopOpt( 5:5 ) )
      CASE ( '~', '*' )
      CASE ( '-', '_', ' ' )
      CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INIT_FIXED_ENV ReadTopOpt: ', &
                             'Unknown top option letter in fifth position'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

    SELECT CASE ( IHOP_TopOpt( 6:6 ) )
      CASE ( 'I' )
      CASE ( ' ' )
      CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INIT_FIXED_ENV ReadTopOpt: ', &
                             'Unknown top option letter in sixth position'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadTopOpt'
    END SELECT

  RETURN
  END !SUBROUTINE ReadTopOpt

  !**********************************************************************!
  SUBROUTINE TopBot( AttenUnit, HS, myThid )
    ! Handles top and bottom boundary conditions
    USE atten_mod, only: CRCI
    USE bdry_mod,  only: HSInfo
    USE ssp_mod,   only: alphaR, betaR, alphaI, betaI, rhoR

  ! ===========================================================================
  !     == Global Variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "IHOP.h"
#endif

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

!    ! ****** Read in BC parameters depending on particular choice ******
!    HS%cp  = 0.0
!    HS%cs  = 0.0
!    HS%rho = 0.0

    ! RG recommends resetting to the default values from ssp_mod.F90
    bPower = 1.0
    fT     = 1D20
    rhoR   = 1.0

    SELECT CASE ( HS%BC )
      CASE ( 'V','R','A','G','F','W','P' )
      CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INIT_FIXED_ENV TopBot: ', &
                             'Unknown boundary condition type'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R TopBot'
    END SELECT

    SELECT CASE ( HS%BC )
      CASE ( 'A' )                  ! *** Half-space properties ***
       ! IEsco23: MISSING IF BOTTOM BC CHECK
       zTemp    = HS%Depth
       alphaR   = IHOP_bcsound
       betaR    = IHOP_bcsoundshear
       rhoR     = IHOP_brho
       alphaI   = IHOP_bcsoundI
       betaI    = IHOP_bcsoundshearI

       ! dummy parameters for a layer with a general power law for attenuation
       ! these are not in play because the AttenUnit for this is not allowed yet
       fT            = 1000.0

       HS%cp  = CRCI( zTemp, alphaR, alphaI, AttenUnit, bPower, fT, myThid )
       HS%cs  = CRCI( zTemp, betaR,  betaI,  AttenUnit, bPower, fT, myThid )

       HS%rho = rhoR
      CASE DEFAULT
    END SELECT

  RETURN
  END !SUBROUTINE TopBot

! **************************************************************************** !
  SUBROUTINE AllocatePos( Nx, x_out, x_in, myThid )
  ! Allocate and populate Pos structure from data.ihop

  !     == Global Variables ==
!#include "SIZE.h"
!#include "GRID.h"
#include "EEPARAMS.h"

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  ! == Local Variables ==
    INTEGER,          INTENT( IN ) :: Nx
    REAL(KIND=_RL90), INTENT( IN ) :: x_in(:)
    REAL(KIND=_RL90), ALLOCATABLE, INTENT( OUT ) :: x_out(:)
    INTEGER           :: i, iAllocStat

    IF ( ALLOCATED(x_out) ) DEALLOCATE(x_out)
    ALLOCATE( x_out(MAX(3, Nx)), STAT=iAllocStat )
    IF ( iAllocStat/=0 ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'IHOP ALLOCATEPOS: failed allocation Pos'
      CALL PRINT_ERROR( msgBuf, myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ALLOCATEPOS'
    ENDIF

    ! set default values
    x_out    = 0.0
    x_out(3) = -999.9

    DO i = 1, Nx
        x_out(i) = x_in(i)
    END DO

  RETURN
  END !SUBROUTINE AllocatePos

  !**********************************************************************!
  SUBROUTINE ReadRunType( RunType, PlotType, myThid )

    ! Read the RunType variable and print to .prt file
    USE srPos_mod, only: Pos

  ! ===========================================================================
  !     == Global Variables ==
#include "EEPARAMS.h"

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    CHARACTER (LEN= 7), INTENT( INOUT ) :: RunType
    CHARACTER (LEN=10), INTENT( INOUT ) :: PlotType

    SELECT CASE ( RunType( 1:1 ) )
      CASE ( 'R','E','I','S','C','A','a','e' )
      CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'INIT_FIXED_ENV ReadRunType: ', &
            'Unknown RunType selected'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadRunType'
    END SELECT

    SELECT CASE ( RunType( 2:2 ) )
      CASE ( 'C','R','S','b','B','g' )
      CASE DEFAULT
        RunType( 2:2 ) = 'G'
    END SELECT

    SELECT CASE ( RunType( 4:4 ) )
      CASE ( 'R','X' )
      CASE DEFAULT
        RunType( 4:4 ) = 'R'
    END SELECT

    SELECT CASE ( RunType( 5:5 ) )
      CASE ( 'R' )
        PlotType = 'rectilin  '
      CASE ( 'I' )
        IF ( Pos%nRZ /= Pos%nRR ) THEN
#ifdef IHOP_WRITE_OUT
          WRITE(msgBuf,'(2A)') 'INIT_FIXED_ENV ReadRunType: ', &
                  'Irregular grid option selected with nRZ != nRR'
          CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
          STOP 'ABNORMAL END: S/R ReadRunType'
        END IF
        PlotType = 'irregular '
      CASE DEFAULT
        RunType( 5:5 ) = 'R'
        PlotType = 'rectilin  '
    END SELECT

    SELECT CASE ( RunType( 6:6 ) )
      CASE ( '2','3' )
      CASE DEFAULT
        RunType( 6:6 ) = '2'
    END SELECT

  RETURN
  END !SUBROUTINE ReadRunType
!**********************************************************************!
END MODULE ihop_init_mod
