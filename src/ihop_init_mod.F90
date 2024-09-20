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
    USE bdry_mod,  only: Bdry, HSInfo
    USE srpos_mod, only: Pos, ReadSxSy, ReadSzRz, ReadRcvrRanges, ReadFreqVec
    USE ssp_mod,   only: SSP, initSSP, alphar
    USE ihop_mod,  only: Beam, rxyz
    USE angle_mod, only: Angles, ReadRayElevationAngles

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
    !     SSP%AttenUnit,Type,Nr,Nz,z,SSP%Seg%r,
    !     Pos%Sx,Sy,Nsz,Nrz,Sz,Rz,Ws,Isz,Wr,Irz,Nrr,Rr,Delta_r,
    !     Beam%RunType,Deltas,Nimage,iBeamWindow,Component,Multiplier,rloop,
    !     Beam%Box%r,Box%z,Type,
    !     Angles%Nalpha,alpha,
    ! - From bdry_mod.F90:initATI
    !   - Top%Natipts,x,
    ! - From bdry_mod.F90:initBTY
    !   - Bot%Natipts,x,
    ! This subroutine will set parameters that shouldn't need to be modified
    ! throughout the MITgcm model run

    ! === Set local parameters ===
    AttenUnit = ''

    ! === Set nonallocatable derived type components from other modules ===
    Bdry%Bot%HS = HSInfo(0.,0.,0.,0., 0.,0. , (0.,0.),(0.,0.), '', '' )
    Bdry%Top%HS = HSInfo(0.,0.,0.,0., 0.,0. , (0.,0.),(0.,0.), '', '' )

    Pos%NSx = -1
    Pos%NSy = -1
    Pos%NSz = -1
    Pos%NRz = -1
    Pos%NRr = -1
    Pos%Ntheta = -1
    Pos%Delta_r = -999.
    Pos%Delta_theta = -999.

    SSP%NPts = -1.
    SSP%Nr = -1.
    SSP%Nx = -1.
    SSP%Ny = -1.
    SSP%Nz = -1.
    SSP%z = -1.
    SSP%rho = -1.
    SSP%c = -1.
    SSP%Type = ''
    SSP%AttenUnit = ''

    Beam%NBeams = -1
    Beam%NImage = -1
    Beam%NSteps = -1
    Beam%iBeamWindow = -1
    Beam%deltas = -1.
    Beam%epsMultiplier = 1.
    Beam%rLoop = -1.
    Beam%Component = ''
    Beam%Type = 'G S '
    Beam%RunType = ''
    Beam%Box = rxyz(0.,0.,0.,0.)

    Angles%Nalpha = 0
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

    ! set Bdry%Top%HS%Depth from first SSP%z
    Bdry%Top%HS%Depth = SSP%z(1)
    ! set water column depth
    Depth = Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth


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
! ================= Below from IHOP.F90: S/R IHOP_MAIN
! =================
!
!
!
!
!    ! AlTImetry: OPTIONAL, default is no ATIFile
!    CALL initATI( Bdry%Top%HS%Opt( 5:5 ), Bdry%Top%HS%Depth, myThid )
!    ! BaThYmetry: OPTIONAL, default is BTYFile
!    CALL initBTY( Bdry%Bot%HS%Opt( 2:2 ), Bdry%Bot%HS%Depth, myThid )
!    ! (top and bottom): OPTIONAL
!    CALL readReflectionCoefficient( Bdry%Bot%HS%Opt( 1:1 ), &
!                                    Bdry%Top%HS%Opt( 2:2 ), myThid )
!    ! Source Beam Pattern: OPTIONAL, default is omni source pattern
!    SBPFlag = Beam%RunType( 3:3 )
!    CALL readPat( myThid )
!    Pos%Ntheta = 1
!    ALLOCATE( Pos%theta( Pos%Ntheta ), Stat = IAllocStat )
!    IF ( IAllocStat/=0 ) THEN
!#ifdef IHOP_WRITE_OUT
!        WRITE(msgBuf,'(2A)') 'IHOP IHOP_INIT: failed allocation Pos%theta'
!        CALL PRINT_ERROR( msgBuf, myThid )
!#endif /* IHOP_WRITE_OUT */
!        STOP 'ABNORMAL END: S/R  IHOP_INIT'
!    ENDIF
!    Pos%theta( 1 ) = 0.
!
!
!! Allocate arrival and U variables on all MPI processes
!    SELECT CASE ( Beam%RunType( 5:5 ) )
!    CASE ( 'I' )
!       NRz_per_range = 1         ! irregular grid
!    CASE DEFAULT
!       NRz_per_range = Pos%NRz   ! rectilinear grid
!    END SELECT
!
!    IF ( ALLOCATED( U ) ) DEALLOCATE( U )
!     SELECT CASE ( Beam%RunType( 1:1 ) )
!     ! for a TL calculation, allocate space for the pressure matrix
!     CASE ( 'C', 'S', 'I' )        ! TL calculation
!          ALLOCATE ( U( NRz_per_range, Pos%NRr ), Stat = iAllocStat )
!          IF ( iAllocStat/=0 ) THEN
!#ifdef IHOP_WRITE_OUT
!              WRITE(msgBuf,'(2A)') 'IHOP IHOP_INIT: ', &
!                             'Insufficient memory for TL matrix: reduce Nr*NRz'
!              CALL PRINT_ERROR( msgBuf,myThid )
!#endif /* IHOP_WRITE_OUT */
!              STOP 'ABNORMAL END: S/R IHOP_INIT'
!          END IF
!          U = 0.0                                    ! init default value
!     CASE ( 'A', 'a', 'R', 'E', 'e' )   ! Arrivals calculation
!          ALLOCATE ( U( 1,1 ), Stat = iAllocStat )   ! open a dummy variable
!          U( 1,1 ) = 0.                              ! init default value
!     CASE DEFAULT
!          ALLOCATE ( U( 1,1 ), Stat = iAllocStat )   ! open a dummy variable
!          U( 1,1 ) = 0.                              ! init default value
!     END SELECT
!
!     ! for an arrivals run, allocate space for arrivals matrices
!     SELECT CASE ( Beam%RunType( 1:1 ) )
!     CASE ( 'A', 'a', 'e' )
!          ! allow space for at least MinNArr arrivals
!          MaxNArr = MAX( ArrivalsStorage / ( NRz_per_range * Pos%NRr ), &
!                         MinNArr )
!          ALLOCATE ( Arr( MaxNArr, Pos%NRr, NRz_per_range ), &
!                     NArr( Pos%NRr, NRz_per_range ), Stat = iAllocStat )
!          IF ( iAllocStat /= 0 ) THEN
!#ifdef IHOP_WRITE_OUT
!              WRITE(msgBuf,'(2A)') 'IHOP IHOP_INIT: ', &
!               'Not enough allocation for Arr; reduce ArrivalsStorage'
!              CALL PRINT_ERROR( msgBuf,myThid )
!#endif /* IHOP_WRITE_OUT */
!              STOP 'ABNORMAL END: S/R IHOP_INIT'
!          END IF
!     CASE DEFAULT
!          MaxNArr = 1
!          ALLOCATE ( Arr( 1, NRz_per_range, Pos%NRr ), &
!                     NArr( Pos%NRr, NRz_per_range ), Stat = iAllocStat )
!     END SELECT
!
!     ! init Arr, Narr
!     ! Arr = something
!     NArr( 1:Pos%NRr, 1:NRz_per_range ) = 0 ! IEsco22 unnecessary? NArr = 0 below
!
!#ifdef IHOP_WRITE_OUT
!     WRITE(msgBuf,'(A)')
!     ! In adjoint mode we do not write output besides on the first run
!     IF (IHOP_dumpfreq.GE.0) &
!       CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
!#endif /* IHOP_WRITE_OUT */
!
!
!
!! open all output files
!    IF ( IHOP_dumpfreq .GE. 0 ) &
!     CALL OpenOutputFiles( IHOP_fileroot, myTime, myIter, myThid )
!
!    ! Run ihop solver on a single processor
!    if (numberOfProcs.gt.1) then
!! Use same single processID as IHOP COST package
!!        if(myProcId.eq.(numberOfProcs-1)) then
!        if(myProcId.eq.0) then
!            CALL CPU_TIME( Tstart )
!            CALL ihopCore(myThid)
!            CALL CPU_TIME( Tstop )
!! Alternitavely, we can broadcast relevant info to all mpi processes Ask P.
!!#ifdef ALLOW_COST
!!            ! Broadcast info to all MPI procs for COST function accumulation
!!            CALL MPI_BCAST(i, 1, MPI_COMPLEX, myProcId, MPI_COMM_MODEL, ierr)
!!
!!#endif /* ALLOW_COST */
!        endif
!    else
!        CALL CPU_TIME( Tstart )
!        CALL ihopCore(myThid)
!        CALL CPU_TIME( Tstop )
!    endif
!
!#ifdef IHOP_WRITE_OUT
!    IF ( IHOP_dumpfreq.GE.0 ) THEN
!        ! print run time
!        if (numberOfProcs.gt.1) then
!            if(myProcId.ne.(numberOfProcs-1)) then
!                WRITE(msgBuf,'(A,I4,A)') 'NOTE: Proc ',myProcId, &
!                    " didn't run ihop"
!                CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
!            endif
!        endif
!        WRITE(msgBuf, '(A)' )
!        CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
!        WRITE(msgBuf, '(A,G15.3,A)' ) 'CPU Time = ', Tstop-Tstart, 's'
!        CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
!
!        ! close all files
!        IF ( IHOP_dumpfreq .GE. 0) THEN
!            SELECT CASE ( Beam%RunType( 1:1 ) )
!            CASE ( 'C', 'S', 'I' )  ! TL calculation
!               CLOSE( SHDFile )
!            CASE ( 'A', 'a' )       ! arrivals calculation
!               CLOSE( ARRFile )
!            CASE ( 'R', 'E' )       ! ray and eigen ray trace
!               CLOSE( RAYFile )
!            CASE ( 'e' )
!               CLOSE( RAYFile )
!               CLOSE( ARRFile )
!               IF ( writeDelay ) CLOSE( DELFile )
!            END SELECT
!
!            if (numberOfProcs.gt.1) then
!                ! Erase prtfiles that aren't on procid = 0
!                if(myProcId.ne.0) then
!                    CLOSE(PRTFile, STATUS='DELETE')
!                else
!                    CLOSE(PRTFile)
!                endif
!            else
!                CLOSE(PRTFile)
!            endif
!        ENDIF
!    ENDIF
!#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE

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

  ! **********************************************************************!
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
        IF ( Pos%NRz /= Pos%NRr ) THEN
#ifdef IHOP_WRITE_OUT
          WRITE(msgBuf,'(2A)') 'INIT_FIXED_ENV ReadRunType: ', &
                  'Irregular grid option selected with NRz not equal to Nr'
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
