#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: ihop_mod
MODULE ihop_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
! Defines modules used in iHOP

! !USES:
  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_USE_MPI
# include "EESUPPORT.h"
#endif

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC  rad2deg, oneCMPLX, &
          PRTFile, RAYFile, DELFile, SHDFile, ARRFile, SSPFile, &
          ATIFile, BTYFile, BRCFile, TRCFile, IRCFile, SBPFile, &
          nMax, nRz_per_range, iStep, afreq, SrcDeclAngle,      &
          Title, Beam, ray2D, ray2DPt, iSmallStepCtr, rxyz
#ifdef ALLOW_USE_MPI
    PUBLIC free_ihop_ray, BcastRay
#endif /* ALLOW_USE_MPI */
!=======================================================================

! == Module variables ==
  ! *** fixed parameters useful for ONLY ihop ***
  REAL(KIND=_RL90),     PARAMETER :: rad2deg = 180.D0 / PI
  COMPLEX (KIND=_RL90), PARAMETER :: oneCMPLX = ( 0.0D0, 1.0D0 )

  INTEGER, PARAMETER :: PRTFile = 61, &    ! standard output file
                        RAYFile = 21, &    ! ray paths file
                        DELFile = 22, &    ! ray paths file
                        SHDFile = 25, &    ! TL calc output file
                        ARRFile = 37, &    ! Arrivals calc output file
                        SSPFile = 40, &    ! optional 2D/3D SSP file
                        ATIFile = 41, &    ! optional 2D/3D altimetry
                        BTYFile = 42, &    ! optional 2D/3D bathymetry
                        BRCFile = 38, TRCFile = 39, IRCFile = 16, &
                        SBPFile = 50, &
                        nMax = 40000

  ! *** varying parameters for ihop ***
  INTEGER            :: iSmallStepCtr = 0
  INTEGER            :: nRz_per_range, iStep
  REAL (KIND=_RL90)  :: afreq, SrcDeclAngle, SrcAzimAngle
  CHARACTER*(80)     :: Title

#ifdef ALLOW_USE_MPI
  INTEGER :: MPI_IHOP_RAY = MPI_DATATYPE_NULL
  LOGICAL :: RAY_TYPE_COMMITTED = .false.
#endif /* ALLOW_USE_MPI */

! == Derived types ==
  ! *** Beam structure ***
  TYPE rxyz
    REAL (KIND=_RL90) :: r, x, y, z
  END TYPE rxyz

  TYPE BeamStructure
    INTEGER           :: nSteps, iBeamWindow
    REAL (KIND=_RL90) :: deltas, epsMultiplier = 1, rLoop
    CHARACTER*(1)     :: Component ! Pressure or displacement
    CHARACTER*(4)     :: Type = 'G S '
    CHARACTER*(7)     :: RunType
    TYPE( rxyz )      :: Box
  END TYPE BeamStructure

  TYPE( BeamStructure ) :: Beam

  ! *** ray structure ***
  TYPE ray2DPt
    INTEGER              :: nTopBnc, nBotBnc, nTurnPt
    REAL    (KIND=_RL90) :: x( 2 ), t( 2 ), p( 2 ), q( 2 )
    REAL    (KIND=_RL90) :: c, Amp, Phase
    COMPLEX (KIND=_RL90) :: tau
  END TYPE ray2DPt

  TYPE( ray2DPt )        :: ray2D( nMax )
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R BcastRay
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: BcastRay
! !INTERFACE:
  SUBROUTINE BcastRay( root, comm )
! !DESCRIPTION:
! Broadcasts the ray2D and Beam data, eg. tau, nSteps

! !USES: None

! !INPUT PARAMETERS:
! root :: MPI root
! comm :: MPI_COMM_WORLD; pre mpi_f08 is an INTEGER
  INTEGER, INTENT( IN ) :: root
  INTEGER, INTENT( IN ) :: comm
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! arrSize    :: indices for range and depth
! ierr       :: MPI return code
  INTEGER :: arrSize
  INTEGER :: ierr
!EOP  

#ifdef ALLOW_USE_MPI
  IF ( MPI_IHOP_RAY.EQ.MPI_DATATYPE_NULL ) THEN
    CALL RayTypeInit( ray2D(1) )
  ENDIF

  ! Broadcast MPI ray2D to all ranks, and free storage
  CALL MPI_Bcast( ray2D, nMax, MPI_IHOP_RAY, root, comm, ierr )

  ! We are on MPI rank 0
  CALL MPI_Bcast( Beam%nSteps, 1, MPI_INTEGER, root, comm, ierr )
#endif /* ALLOW_USE_MPI */

  RETURN
  END !SUBROUTINE BcastRay

#ifdef ALLOW_USE_MPI
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: free_ihop_ray
! !INTERFACE:
  SUBROUTINE free_ihop_ray(myThid)
! !DESCRIPTION:
! Free arrival datatype
    INTEGER, INTENT( IN ) :: myThid
    INTEGER :: ierr

    IF (RAY_TYPE_COMMITTED) THEN
      CALL MPI_Type_free(MPI_IHOP_RAY, ierr)
      MPI_IHOP_RAY = MPI_DATATYPE_NULL
      RAY_TYPE_COMMITTED = .false.
    ENDIF

  RETURN
  END !SUBROUTINE free_ihop_ray

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: RayTypeInti
! !INTERFACE:
  SUBROUTINE RayTypeInit( singleRayPt )
! !DESCRIPTION:
! Initializes the ray MPI Datatype

! !USES: None

! !INPUT PARAMETERS:
  TYPE( ray2DPt ), INTENT( IN ) :: singleRayPt
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! disp       :: address for new MPI datatype ray2dpt parameter
! base, addr :: base and address for current ray2dpt datatype and parameters
! ty         :: datatype of each ray2dpt parameter
! MPI_RL, MPI_CL :: MPI type depending on _RL90
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(11), base, addr(11)
  INTEGER :: ierr, n, bl(11), ty(11)
  INTEGER :: MPI_RL, MPI_CL
!EOP

  IF (RAY_TYPE_COMMITTED) RETURN

  ! Build Ray2Dpt MPI type
  IF (STORAGE_SIZE( real(singleRayPt%tau) ).EQ.64) THEN
    MPI_RL = MPI_DOUBLE_PRECISION
    MPI_CL = MPI_DOUBLE_COMPLEX
  ELSEIF (STORAGE_SIZE( real(singleRayPt%tau) ).EQ.32) THEN
    MPI_RL = MPI_REAL
    MPI_CL = MPI_COMPLEX
  ELSE
    STOP "ABNORMAL END MPI_DATATYPE: Unsupported _RL90 size for MPI"
  ENDIF

  ! Initiate all bl to 1 since all ray2D parameters are scalars
  bl=1
  CALL MPI_Get_address(singleRayPt, base, ierr)

  n=1
  CALL MPI_Get_address(singleRayPt%nTopBnc, addr(n), ierr)
  ty(n)=MPI_INTEGER

!  n=n+1
!  CALL MPI_Get_address(singleRayPt%nBotBnc, addr(n), ierr)
!  ty(n)=MPI_INTEGER
!
!  n=n+1
!  CALL MPI_Get_address(singleRayPt%nTurnPt, addr(n), ierr)
!  ty(n)=MPI_INTEGER
!
!  n=n+1
!  CALL MPI_Get_address(singleRayPt%x(1), addr(n), ierr)
!  ty(n)=MPI_RL
!  bl(n)=2
!
!  n=n+1
!  CALL MPI_Get_address(singleRayPt%t(1), addr(n), ierr)
!  ty(n)=MPI_RL
!  bl(n)=2
!
!  n=n+1
!  CALL MPI_Get_address(singleRayPt%p(1), addr(n), ierr)
!  ty(n)=MPI_RL
!  bl(n)=2
!
!  n=n+1
!  CALL MPI_Get_address(singleRayPt%q(1), addr(n), ierr)
!  ty(n)=MPI_RL
!  bl(n)=2
!
!  n=n+1
!  CALL MPI_Get_address(singleRayPt%c, addr(n), ierr)
!  ty(n)=MPI_RL
!
!  n=n+1
!  CALL MPI_Get_address(singleRayPt%Amp, addr(n), ierr)
!  ty(n)=MPI_RL
!
!  n=n+1
!  CALL MPI_Get_address(singleRayPt%Phase, addr(n), ierr)
!  ty(n)=MPI_RL

  n=n+1
  CALL MPI_Get_address(singleRayPt%tau, addr(n), ierr)
  ty(n)=MPI_CL

  disp(1:n)=addr(1:n)-base

  CALL MPI_Type_create_struct(n, bl, disp, ty, MPI_IHOP_RAY, ierr)
  CALL MPI_Type_commit(MPI_IHOP_RAY, ierr)
  RAY_TYPE_COMMITTED = .true.

RETURN
END
#endif /* ALLOW_USE_MPI */

END !MODULE ihop_mod
