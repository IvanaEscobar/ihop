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
    PUBLIC BcastRay
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
    REAL    (KIND=_RL90) :: x(2), t(2), p(2), q(2), h(2)
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
#ifdef ALLOW_USE_MPI
  INTEGER :: ierr
#endif /* ALLOW_USE_MPI */
#ifdef ALLOW_USE_MPI
  INTEGER :: MPI_IHOP_RAY2D = MPI_DATATYPE_NULL
  LOGICAL :: RAY2D_TYPE_COMMITTED = .false.
#endif /* ALLOW_USE_MPI */

!EOP  

#ifdef ALLOW_USE_MPI
! build custom MPI datatype for ray2d
  IF ( MPI_IHOP_RAY2D.EQ.MPI_DATATYPE_NULL ) THEN
    CALL ray2dPtTypeInit( ray2D(1), ray2D(2), &
      MPI_IHOP_RAY2D, RAY2D_TYPE_COMMITTED ) 
  ENDIF

  ! We are on MPI rank 0
  CALL MPI_Bcast( Beam%nSteps, 1, MPI_INTEGER, root, comm, ierr )

  !CALL MPI_Bcast( ray2d(beam%nsteps)%q(1), 1, MPI_RL, root, comm, ierr)
  CALL MPI_Bcast( ray2D, nMax, MPI_IHOP_RAY2D, root, comm, ierr )
#endif /* ALLOW_USE_MPI */

  RETURN
  END !SUBROUTINE BcastRay

#ifdef ALLOW_USE_MPI
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: free_ihop_ray2d
! !INTERFACE:
  SUBROUTINE free_ihop_ray2d(myThid)
! !DESCRIPTION:
! Free ray2d datatype
    INTEGER, INTENT( IN ) :: myThid
!    INTEGER :: ierr
!
!    IF (RAY2D_TYPE_COMMITTED) THEN
!      CALL MPI_Type_free(MPI_IHOP_RAY2D, ierr)
!      MPI_IHOP_RAY2D = MPI_DATATYPE_NULL
!      RAY2D_TYPE_COMMITTED = .false.
!    ENDIF

  RETURN
  END !SUBROUTINE free_ihop_arrival

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: Ray2DPtTypeInti
! !INTERFACE:
  SUBROUTINE Ray2DPtTypeInit( singleRay2D, nextRay2D, myType, TYPE_COMMITTED )

! !INPUT PARAMETERS:
  TYPE( ray2dpt ), INTENT( IN ) :: singleRay2D
  TYPE( ray2dpt ), INTENT( IN ) :: nextRay2D
  INTEGER, INTENT( INOUT ) :: MYTYPE
  LOGICAL, INTENT( INOUT ) :: TYPE_COMMITTED

  ! LOCAL VARIABLES:
  INTEGER, PARAMETER :: typeSize = 3
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(typeSize), addr(typeSize)
  INTEGER(KIND=MPI_ADDRESS_KIND) :: base, extent, a1, a2
  INTEGER :: n, bl(typeSize), ty(typeSize)
  INTEGER :: ierr, MPI_CL, MPI_RL, tmpType

  IF (TYPE_COMMITTED) RETURN

  ! Build partial ray2D MPI type
  IF (STORAGE_SIZE( REAL(singleRay2D%tau) ).EQ.64) THEN
    MPI_CL = MPI_DOUBLE_COMPLEX
    MPI_RL = MPI_DOUBLE_PRECISION
  ELSE
    MPI_CL = MPI_COMPLEX
    MPI_RL = MPI_REAL
  ENDIF

  bl=1
  CALL MPI_Get_address(singleRay2D, base, ierr)

  n=1
  CALL MPI_Get_address(singleRay2D%x(1), addr(n), ierr)
  bl(n)=2
  disp(n) = addr(n) - base
  ty(n)=MPI_RL

  n=n+1
  CALL MPI_Get_address(singleRay2D%q(1), addr(n), ierr)
  bl(n)=2
  disp(n) = addr(n) - base
  ty(n)=MPI_RL

  n=n+1
  CALL MPI_Get_address(singleRay2D%tau, addr(n), ierr)
  disp(n) = addr(n) - base
  ty(n)=MPI_CL

  CALL MPI_Type_create_struct(n, bl, disp, ty, tmpType, ierr)

  CALL MPI_get_address(singleRay2D, a1, ierr)
  CALL MPI_get_address(nextRay2D, a2, ierr)
  extent = a2-a1

  CALL MPI_Type_create_resized( tmpType, 0_MPI_ADDRESS_KIND, extent, &
                               myType, ierr )
  CALL MPI_Type_commit( myType, ierr )
  CALL MPI_Type_free( tmpType, ierr )
  TYPE_COMMITTED = .true.

  RETURN
  END
#endif /* ALLOW_USE_MPI */

END !MODULE ihop_mod
