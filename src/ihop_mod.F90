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

#ifdef ALLOW_USE_MPI
  INTEGER :: MPI_IHOP_RAY2D = MPI_DATATYPE_NULL
  LOGICAL :: RAY2D_TYPE_COMMITTED = .false.
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
  COMPLEX (KIND=_RL90) :: tauBuf(nMax)
  REAL (KIND=_RL90) :: qBuf(nMax), xBuf(2, nMax)
  INTEGER :: i
  TYPE( ray2DPt ) :: singleRay2D
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(7), base, addr(7)
  INTEGER :: ierr, n, bl(7), ty(7)
  INTEGER :: ierr, MPI_CL, MPI_RL
!EOP  

#ifdef ALLOW_USE_MPI
! build custom MPI datatype for ray2d
  IF ( MPI_IHOP_RAY2D.EQ.MPI_DATATYPE_NULL .AND. &
       .NOT.RAY2D_TYPE_COMMITTED ) THEN
    singleRay2D = ray2D(1)
    IF (STORAGE_SIZE( REAL(ray2d(1)%tau) ).EQ.64) THEN
      MPI_CL = MPI_DOUBLE_COMPLEX
      MPI_RL = MPI_DOUBLE_PRECISION
    ELSE
      MPI_CL = MPI_COMPLEX
      MPI_RL = MPI_REAL
    ELSE
      STOP "ABNORMAL END RAY MPI_DATATYPE: Unsupported _RL90 size for MPI"
    ENDIF

    bl=1
    CALL MPI_Get_address(singleRay2D, base, ierr)

    n=1
    CALL MPI_Get_address(singleRay2D%x, addr(n), ierr)
    ty(n)=MPI_RL

    n=n+1
    CALL MPI_Get_address(singleRay2D%q, addr(n), ierr)
    ty(n)=MPI_RL

    n=n+1
    CALL MPI_Get_address(singleRay2D%tau, addr(n), ierr)
    ty(n)=MPI_CL

    disp(1:n)=addr(1:n)-base

    CALL MPI_Type_create_struct(n, bl, disp, ty, MPI_IHOP_RAY2D, ierr)
    CALL MPI_Type_commit(MPI_IHOP_RAY2D, ierr)
    RAY2D_TYPE_COMMITTED = .true.

  ENDIF

  ! Build Ray2Dpt MPI type
  IF (STORAGE_SIZE( REAL(ray2d(1)%tau) ).EQ.64) THEN
    MPI_CL = MPI_DOUBLE_COMPLEX
    MPI_RL = MPI_DOUBLE_PRECISION
  ELSE
    MPI_CL = MPI_COMPLEX
    MPI_RL = MPI_REAL
  ENDIF

!  IF (myProcID.eq.root) THEN
!    DO i=1,nMax
!      tauBuf(i) = ray2D(i)%tau
!      qBuf(i) = ray2D(i)%q(1)
!      xBuf(1,i) = ray2D(i)%x(1)
!      xBuf(2,i) = ray2D(i)%x(2)
!    ENDDO
!  ENDIF
!  CALL MPI_Bcast( tauBuf, nMax, MPI_CL, root, comm, ierr )
!  CALL MPI_Bcast( qBuf, nMax, MPI_RL, root, comm, ierr )
!  CALL MPI_Bcast( xBuf, 2*nMax, MPI_RL, root, comm, ierr )
!
!  IF (myProcID.ne.root) THEN
!    DO i=1,nMax
!      ray2D(i)%tau  = tauBuf(i)
!      ray2D(i)%q(1) = qBuf(i)
!      ray2D(i)%x(1) = xBuf(1,i)
!      ray2D(i)%x(2) = xBuf(2,i)
!    ENDDO
!  ENDIF

  ! We are on MPI rank 0
  CALL MPI_Bcast( Beam%nSteps, 1, MPI_INTEGER, root, comm, ierr )
#endif /* ALLOW_USE_MPI */

  RETURN
  END !SUBROUTINE BcastRay

END !MODULE ihop_mod
