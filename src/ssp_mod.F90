#include "IHOP_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: ssp_mod
MODULE ssp_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
! This module contains the subroutines and functions to read, set, and evaluate
! SSP (Sound Speed Profile) data.

! !USES:
  USE ihop_mod, only: PRTFile
  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
!#ifdef ALLOW_USE_MPI
# include "EESUPPORT.h"
!#endif
#include "PARAMS.h"
#include "GRID.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC init_fixed_SSP, init_varia_SSP, setSSP, evalSSP, &
         Grid, SSP, &
         alphaR, betaR, alphaI, betaI, rhoR, iSegz, iSegr
!=======================================================================

! !FUNCTIONS:
  _RL CHEN_MILLERO
  EXTERNAL CHEN_MILLERO

! == Module Variables ==
  INTEGER, PARAMETER     :: MaxSSP = 201
  INTEGER                :: iSegr = 1, iSegz = 1
  LOGICAL                :: foundr, foundz
#ifdef IHOP_THREED
  INTEGER                :: iSegx = 1, iSegy = 1
  LOGICAL                :: foundx, foundy
#endif /* IHOP_THREED */
  INTEGER                :: iostat, iallocstat
  INTEGER,           PRIVATE :: iz
  REAL (KIND=_RL90), PRIVATE :: Depth, W
  ! DEFAULT values, IHOP only modifies alphaR
  REAL (KIND=_RL90)      :: alphaR = 1500, betaR = 0, alphaI = 0, &
                            betaI = 0, rhoR = 1
  ! SSP interpolation parameters, only used in ssp_mod
  COMPLEX (KIND=_RL90) :: n2(MaxSSP), n2z(MaxSSP)
  COMPLEX (KIND=_RL90) :: cSpln( 4, MaxSSP ), cCoef( 4, MaxSSP )

! == Derived types ==
  TYPE rxyz_vector
    REAL (KIND=_RL90), ALLOCATABLE :: r(:)
#ifdef IHOP_THREED
    REAL (KIND=_RL90), ALLOCATABLE :: x(:), y(:), z(:)
#endif
  END TYPE rxyz_vector

  TYPE SSPGrid
    INTEGER              :: nPts, nR, nX, nY, nZ
    REAL ( KIND=_RL90 )  :: z( MaxSSP ), rho( MaxSSP )
    TYPE ( rxyz_vector ) :: Seg
    CHARACTER*(1)        :: Type
    CHARACTER*(2)        :: AttenUnit
  END TYPE SSPGrid

  TYPE SSPVariable
    COMPLEX (KIND=_RL90)              :: c( MaxSSP ), cz( MaxSSP )
    REAL    (KIND=_RL90), ALLOCATABLE :: cMat( :,: ), czMat( :,: )
#ifdef IHOP_THREED
    REAL    (KIND=_RL90), ALLOCATABLE :: cMat3( :,:,: ), czMat3( :,:,: )
#endif /* IHOP_THREED */
  END TYPE SSPVariable

  TYPE( SSPGrid )     :: Grid
  TYPE( SSPVariable ) :: SSP
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R init_fixed_SSP
! S/R ReadSSP
! S/R init_fixed_Grid
! S/R init_varia_SSP

! S/R setSSP
! S/R evalSSP
! S/R n2Linear
! S/R cLinear
! S/R cPCHIP
! S/R cCubic
! S/R Quad
! S/R gcmSSP
! S/R writeSSP
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: init_fixed_SSP
! !INTERFACE:
  SUBROUTINE init_fixed_SSP( myThid )
! !DESCRIPTION:
!   Initialize the SSP derived type and set the SSP structures.

! !USES: None

! !INPUT PARAMETERS:
! myThid :: my thread ID 
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES: None
!EOP

  ! Check if SSPFile exists
  IF (useSSPFile) THEN
    CALL ReadSSP( myThid )
  ELSE
    CALL init_fixed_Grid( myThid )
  END IF

  RETURN
  END !SUBROUTINE init_fixed_SSP

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ReadSSP (LEGACY)
! !INTERFACE:
  SUBROUTINE ReadSSP( myThid )
! !DESCRIPTION:
!   Read SSP [m/s] from file and populate SSPStructure: SSP.

! !USES:
  USE atten_mod, only: CRCI
  USE ihop_mod,  only: SSPFile
  USE bdry_mod,  only: Bdry
! IESCO24
! fT = 1000 ONLY for acousto-elastic halfspaces, I will have to pass this
! parameter in a different way after ssp_mod is split btwn fixed and varia
! USE initenvihop, only: fT

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! Depth :: Depth of the bottom boundary
! bPower :: Power for attenuation calculation
! fT :: Frequency for attenuation calculation
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  REAL (KIND=_RL90) :: Depth
  REAL (KIND=_RL90) :: bPower, fT
!EOP

  ! IESCO24 fT init
  bPower = 1.0
  fT = 1000.0
  Depth = Bdry%Bot%HS%Depth

  ! I/O on main thread only
  _BEGIN_MASTER(myThid)

  ! OPEN SSPFile to read
  OPEN ( FILE=TRIM(IHOP_fileroot) // '.ssp', UNIT=SSPFile, &
    FORM='FORMATTED', STATUS='OLD', IOSTAT=iostat )
  IF ( IOSTAT.NE.0 ) THEN   ! successful open?
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 'SSPMOD ReadSSP: Unable to open the SSP file'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadSSP'
  ENDIF

  READ( SSPFile,  * ) Grid%nR, Grid%nZ

  ALLOCATE( SSP%cMat( Grid%nZ, Grid%nR ), &
            SSP%czMat( Grid%nZ-1, Grid%nR ), &
            Grid%Seg%R( Grid%nR ), &
            STAT=iallocstat )
  IF ( iallocstat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SSPMOD ReadSSP: ', &
      'Insufficient memory to store SSP'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadSSP'
  ENDIF

  ! Initiate to nonsense
  SSP%cMat  = -99.0 _d 0
  SSP%czMat = -99.0 _d 0

  ! Set SSP contents
  READ( SSPFile,  * ) Grid%Seg%R( 1:Grid%nR )
  Grid%Seg%R = 1000.0 * Grid%Seg%R   ! convert km to m

  READ( SSPFile,  * ) Grid%Z( 1:Grid%nZ )

  DO iz = 1, Grid%nZ
    READ(  SSPFile, * ) SSP%cMat( iz, : )
  ENDDO

  CLOSE( SSPFile )

  Grid%nPts = 1
  DO iz = 1, MaxSSP
    alphaR = SSP%cMat( iz, 1 )

    SSP%c(iz) = CRCI( Grid%Z(iz), alphaR, alphaI, Grid%AttenUnit, bPower, fT, &
                      myThid )
    Grid%rho(iz) = rhoR !IEsco22: set to a default value of 1

    ! verify depths are monotone increasing
    IF ( iz.GT.1 ) THEN
      IF ( Grid%Z( iz ).LE.Grid%Z( iz-1 ) ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A,F10.2)') 'SSPMOD ReadSSP: ', &
          'The depths in the SSP must be monotone increasing', Grid%Z(iz)
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSSP'
      ENDIF
    ENDIF

    ! compute gradient, cz
    IF ( iz.GT.1 ) SSP%cz( iz-1 ) = ( SSP%c( iz ) - SSP%c( iz-1 ) ) / &
                                    ( Grid%Z( iz ) - Grid%Z( iz-1 ) )

    ! Did we read the last point?
    IF ( ABS( Grid%Z( iz ) - Depth ).LT.100.*EPSILON( 1.0e0 ) ) THEN
      IF ( Grid%nPts.EQ.1 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SSPMOD ReadSSP: ', &
          'The SSP must have at least 2 points'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadSSP'
      ENDIF

      RETURN

    ENDIF ! IF ( ABS( Grid%Z( iz ) - Depth ).LT.100.*EPSILON( 1.0e0 ) )

    Grid%nPts = Grid%nPts + 1

  ENDDO ! DO iz = 1, Grid%nZ

  ! Fall through means too many points in the profile
#ifdef IHOP_WRITE_OUT
  WRITE(msgBuf,'(2A)') 'SSPMOD ReadSSP: ', &
    'Number of SSP points exceeds limit'
  CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
  STOP 'ABNORMAL END: S/R ReadSSP'

  ! I/O on main thread only
  _END_MASTER(myThid)

  RETURN
  END !SUBROUTINE ReadSSP

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: init_fixed_Grid
! !INTERFACE:
  SUBROUTINE init_fixed_Grid( myThid )
! !DESCRIPTION:
!   Initialize SSP Grid parameters that do not change within a time series.
! Sets Grid%nR,nZ,Seg%R, and ihop_sumweights

! !USES:
  USE bdry_mod, only: Bdry

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! skip_range :: Flag to skip range interpolation
! iallocstat :: Allocation status
! bi, bj, i, j, k, ii, jj :: GCM domain loop indices
! nii, njj :: Number of IDW points in range
! tolerance :: Tolerance for IDW interpolation
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  LOGICAL :: skip_range
  INTEGER :: iallocstat
  INTEGER :: bi,bj, i,j,k, ii,jj
  INTEGER :: nii(IHOP_NPTS_RANGE), njj(IHOP_NPTS_RANGE)
  REAL (KIND=_RL90) :: tolerance

  ! init local vars
  skip_range =.false.
  nii(:)    = 0
  njj(:)    = 0
  tolerance  = 5 _d -5

  ! init default Grid values (only fixed memory vars)
  Grid%nPts = -1
  Grid%Z    = -999.0
  Grid%rho  = -999.0

  ! set ihop SSP Grid size
  Grid%nZ = Nr+2 ! add z=0 z=Depth layers to GCM Nr
  Grid%nR = IHOP_NPTS_RANGE
  Grid%nPts = Grid%nZ

  ! set Grid%Z from rC, rkSign=-1 used bc ihop uses +ive depths
  Grid%Z( 1 )            = 0.0 _d 0
  Grid%Z( 2:(Grid%nZ-1) ) = rkSign*rC( 1:Nr )
  Grid%Z( Grid%nZ )       = Bdry%Bot%HS%Depth ! rkSign*rF(Nr+1)*1.05

  ! set Grid%Seg%R from data.ihop -> ihop_ranges
  !IF (ALLOCATED(Grid%Seg%R)) DEALLOCATE(Grid%Seg%R)
  ALLOCATE( Grid%Seg%R( Grid%nR ), STAT=iallocstat )
  IF ( iallocstat.NE.0 ) THEN
# ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SSPMOD init_fixed_Grid: ', &
      'Insufficient memory to store Grid%Seg%R'
    CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R init_fixed_Grid'
  ENDIF

  Grid%Seg%R( 1:Grid%nR ) = ihop_ranges( 1:Grid%nR )
  ! Modify from [m] to [km]
  Grid%Seg%R = 1000.0 * Grid%Seg%R

  !==================================================
  ! IDW Interpolate: COMPARING with LAT LONs (xC, yC)
  !==================================================
  ! Sum IDW weights
  DO ii = 1, Grid%nR
    ihop_sumweights(ii,:) = sum(ihop_idw_weights(ii,:))
  ENDDO

  ! Adapt IDW interpolation by bathymetry
  DO bj=myByLo(myThid),myByHi(myThid)
    DO bi=myBxLo(myThid),myBxHi(myThid)
      DO j=1,sNy
        DO i=1,sNx
          DO ii=1, Grid%nR
            skip_range = .FALSE.

            DO jj=1,IHOP_npts_idw
              IF ( ABS(xC(i, j, bi, bj) - ihop_xc(ii, jj)).LE.tolerance .AND. &
                   ABS(yC(i, j, bi, bj) - ihop_yc(ii, jj)).LE.tolerance ) THEN
                DO k=1,Nr
    ! No IDW interpolation on xc, yc centered ranges
    IF ( nii(ii).EQ.1 .AND. k.GT.njj(ii) ) THEN
      skip_range = .TRUE.
    ENDIF

    IF ( .NOT.skip_range ) THEN
      IF ( hFacC(i, j, k, bi, bj).EQ.0. ) THEN
        ihop_sumweights(ii, k) = &
          ihop_sumweights(ii, k) - ihop_idw_weights(ii, jj)

        ! No interpolation on xc, yc centered ranges
        IF ( ihop_idw_weights(ii, jj).EQ.0. ) THEN
          ihop_sumweights(ii, k:Nr) = 0.0
          nii(ii) = 1
          njj(ii) = k
        ENDIF

      ENDIF ! IF ( hFacC(i, j, k, bi, bj).EQ.0. )

      ! Set TINY and negative values to 0.0
      IF (ihop_sumweights(ii, k).LT.1D-13) ihop_sumweights(ii, k) = 0.0

    ENDIF ! IF ( .NOT.skip_range )
                ENDDO !k
              ENDIF
            ENDDO !jj
          ENDDO !ii
        ENDDO !i
      ENDDO !j
    ENDDO !bi
  ENDDO !bj

  RETURN
  END !SUBROUTINE init_fixed_Grid

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: init_varia_SSP
! !INTERFACE:
  SUBROUTINE init_varia_SSP( myThid )
! !DESCRIPTION:
!   Initialize the fixed SSP parameters that do not change within a time series.
! Sets SSP%c,cz,cMat,czMat

! !USES: None

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! iallocstat :: Allocation status
#ifdef IHOP_WRITE_OUT
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
#endif /* IHOP_WRITE_OUT */
  INTEGER :: iallocstat

  ! ONLY ALLOCATE cmat and czmat, to be filled per ihop run
  ALLOCATE( SSP%cMat( Grid%nZ, Grid%nR ), &
            SSP%czMat( Grid%nZ-1, Grid%nR ), &
            STAT=iallocstat )
  IF ( iallocstat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SSPMOD init_fixed_SSP: ', &
      'Insufficient memory to store SSP%cMat, SSP%czMat'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R init_fixed_SSP'
  ENDIF

  ! Initiate to nonsense
  SSP%cMat  = -99.0 _d 0
  SSP%czMat = -99.0 _d 0

  ! init default SSP values (only fixed memory vars)
  SSP%c    = -99.0 _d 0
  SSP%cz   = -99.0 _d 0

  RETURN
  END !SUBROUTINE init_varia_SSP

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: setSSP
! !INTERFACE:
  SUBROUTINE setSSP( myThid )
! !DESCRIPTION:
!   Set the SSP derived type based on the SSP interpolation scheme.

! !USES:
  USE pchip_mod,  only: PCHIP
  USE splinec_mod,only: cspline

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! ir, iz :: Loop indices
! mpiRC  :: MPI return code
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: ir, iz
#ifdef ALLOW_USE_MPI
  INTEGER :: mpiRC
  mpiRC  = 0
#endif
!EOP

!$TAF init setssp1 = 'ssp_mod_setssp'

! IESCO24: Write derived type with allocatable memory by type: SSP from ssp_mod
! Scalar components
! Fixed arrays
!$TAF store ssp%c = setssp1
! Allocatable arrays
!$TAF store ssp%cmat,ssp%czmat,grid%seg%r = setssp1

  ! init defaults for ssp_mod scoped arrays
  n2    = (-1.,-1.)
  n2z   = (-1.,-1.)
  cSpln = (-1.,-1.)
  cCoef = (-1.,-1.)

  ! Extract gcm SSP field
  IF ( .NOT.useSSPFile ) CALL gcmSSP( myThid )

  ! Write to PRTFile
  IF ( .NOT.usingMPI ) THEN
    CALL writeSSP( myThid )
#ifdef ALLOW_USE_MPI
  ELSE ! using MPI
    CALL MPI_COMM_RANK( MPI_COMM_MODEL, mpiMyId, mpiRC )
    myProcId = mpiMyId

    ! Hard coded write on single processor
    IF ( myProcId.EQ.0 ) CALL writeSSP( myThid )
#endif /* ALLOW_USE_MPI */
  ENDIF ! IF ( .NOT.usingMPI )

  ! Populate rest of SSP derived type based on SSP interpolation scheme
  SELECT CASE ( Grid%Type )
  CASE ( 'N' )  !  N2-linear profile option
    n2( 1:Grid%nPts ) = 1.0 / SSP%c( 1:Grid%nPts )**2
    !IEsco23 Test this: n2(  1:Grid%nZ ) = 1.0 / SSP%c( 1:Grid%nZ )**2

    ! compute gradient, n2z
    DO iz = 2, Grid%nPts
      n2z( iz-1 ) = (  n2(   iz ) - n2(    iz-1 ) ) / &
                    ( Grid%Z( iz ) - Grid%Z( iz-1 ) )
    END DO

  CASE ( 'C' )  !  C-linear profile option
  CASE ( 'P' )  !  monotone PCHIP ACS profile option
    !                                                               2      3
    ! compute coefficients of std cubic polynomial: c0 + c1*x + c2*x + c3*x
    !
    CALL PCHIP( Grid%Z, SSP%c, Grid%nPts, cCoef, cSpln )
!IEsco23 Test this:
!     CALL PCHIP( Grid%Z, SSP%c, Grid%nZ, cCoef, cSpln )

  CASE ( 'S' )  !  Cubic spline profile option
    cSpln( 1, 1:Grid%nPts ) = SSP%c( 1:Grid%nPts )
!IEsco23 Test this:
!     cSpln( 1, 1 : Grid%nZ ) = SSP%c( 1 : Grid%nZ )

    ! Compute spline coefs
    CALL cSpline( Grid%Z, cSpln( 1, 1 ), Grid%nPts, 0, 0, Grid%nPts )
!IEsco23 Test this:
!     CALL CSpline( Grid%Z, cSpln( 1,1 ), Grid%nZ,iBCBeg, iBCEnd, Grid%nZ )

  CASE ( 'Q' )
    ! calculate cz
    DO ir = 1, Grid%nR
      DO iz = 2, Grid%nZ
        ! delta_z = ( Grid%Z( iz2 ) - Grid%Z( iz2-1 ) )
        SSP%czMat( iz-1, ir ) = ( SSP%cMat( iz,   ir ) - &
                                  SSP%cMat( iz-1, ir ) ) / &
                                ( Grid%Z( iz ) - Grid%Z( iz-1 ) )
      END DO
    END DO

  CASE DEFAULT
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 'SSPMOD setSSP: Invalid SSP profile option'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R setSSP'
  END SELECT

  RETURN
  END !SUBROUTINE setSSP

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: evalSSP
! !INTERFACE:
  SUBROUTINE evalSSP( x, c, cimag, gradc, crr, crz, czz, rho, myThid )
! !DESCRIPTION:
!   Evaluate SSP at a given point x = (r, z) and return the sound speed

! !USES: None

! !INPUT PARAMETERS:
! x      :: r-z SSP evaluation point
! c      :: sound speed at x
! cimag  :: imaginary part of sound speed at x
! gradc  :: gradient of sound speed at x
! crr    :: radial derivative of sound speed at x
! crz    :: radial-z derivative of sound speed at x
! czz    :: z-z derivative of sound speed at x
! rho    :: density at x
! myThid :: my thread ID
  REAL (KIND=_RL90), INTENT( IN  ) :: x(2)
  REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, rho
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: c, cimag, gradc, crr, crz, czz, rho

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
  CHARACTER*(MAX_LEN_MBUF):: msgBuf

  SELECT CASE ( Grid%Type )
  CASE ( 'N' )  !  N2-linear profile option
    CALL n2Linear( x, c, cimag, gradc, crr, crz, czz, rho, myThid )
  CASE ( 'C' )  !  C-linear profile option
    CALL cLinear(  x, c, cimag, gradc, crr, crz, czz, rho, myThid )
  CASE ( 'P' )  !  monotone PCHIP ACS profile option
    CALL cPCHIP(   x, c, cimag, gradc, crr, crz, czz, rho, myThid )
  CASE ( 'S' )  !  Cubic spline profile option
    CALL cCubic(   x, c, cimag, gradc, crr, crz, czz, rho, myThid )
  CASE ( 'Q' )
    CALL Quad(     x, c, cimag, gradc, crr, crz, czz, rho, myThid )
  CASE DEFAULT
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'Profile option: ', Grid%Type
    ! In adjoint mode we do not write output besides on the first run
    IF ( IHOP_dumpfreq.GE.0 ) &
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 'SSPMOD evalSSP: Invalid SSP profile option'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R evalSSP'
    c = 0.         !RG
    cimag = 0.     !RG
    gradc = 0.     !RG
    crr   = 0.     !RG
    crz   = 0.     !RG
    czz   = 0.     !RG
    rho   = 0.     !RG

  END SELECT

  RETURN
  END !SUBROUTINE evalSSP

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: n2Linear
! !INTERFACE:
  SUBROUTINE n2Linear( x, c, cimag, gradc, crr, crz, czz, rho, myThid )
! !DESCRIPTION:
!   N2-linear interpolation of SSP data.

! !USES: None

! !INPUT PARAMETERS:
! x      :: r-z SSP evaluation point
! c      :: sound speed at x
! cimag  :: imaginary part of sound speed at x
! gradc  :: gradient of sound speed at x
! crr    :: radial derivative of sound speed at x
! crz    :: radial-z derivative of sound speed at x
! czz    :: z-z derivative of sound speed at x
! rho    :: density at x
! myThid :: my thread ID
  REAL (KIND=_RL90), INTENT( IN  ) :: x(2)
  REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, rho
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: c, cimag, gradc, crr, crz, czz, rho

! !LOCAL VARIABLES: None
!EOP

  iSegz = 1 !RG
  IF ( x( 2 ).LT.Grid%Z( iSegz ) .OR. x( 2 ).GT.Grid%Z( iSegz+1 ) ) THEN
    foundz=.false.
!IEsco23 Test this:
!     DO iz = 2, Grid%nZ   ! Search for bracketting Depths
    DO iz = 2, Grid%nPts   ! Search for bracketting Depths
      IF ( x( 2 ).LT.Grid%Z( iz ) .AND. .NOT.foundz ) THEN
        iSegz  = iz - 1
        foundz = .true.
      ENDIF

    ENDDO

  ENDIF ! IF ( x( 2 ).LT.Grid%Z( iSegz ) .OR. x( 2 ).GT.Grid%Z( iSegz+1 ) )

  W = ( x( 2 ) - Grid%Z( iSegz ) ) / ( Grid%Z( iSegz+1 ) - Grid%Z( iSegz ) )

  c     = REAL(  1.0D0 / SQRT( ( 1.0D0-W ) * n2( iSegz ) &
          + W * n2( iSegz+1 ) ) )
  cimag = AIMAG( 1.0D0 / SQRT( ( 1.0D0-W ) * n2( iSegz ) &
          + W * n2( iSegz+1 ) ) )

  gradc = [ 0.0D0, -0.5D0 * c * c * c * REAL( n2z( iSegz ) ) ]
  crr   = 0.0d0
  crz   = 0.0d0
  czz   = 3.0d0 * gradc( 2 ) * gradc( 2 ) / c

  rho   = ( 1.0D0-W ) * Grid%rho( iSegz ) + W * Grid%rho( iSegz+1 )

  RETURN
  END !SUBROUTINE n2Linear

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: cLinear
! !INTERFACE:
  SUBROUTINE cLinear( x, c, cimag, gradc, crr, crz, czz, rho, myThid  )
! !DESCRIPTION:
!   C-linear interpolation of SSP data.

! !USES: None

! !INPUT PARAMETERS:
! x      :: r-z SSP evaluation point
! c      :: sound speed at x
! cimag  :: imaginary part of sound speed at x
! gradc  :: gradient of sound speed at x
! crr    :: radial derivative of sound speed at x
! crz    :: radial-z derivative of sound speed at x
! czz    :: z-z derivative of sound speed at x
! rho    :: density at x
! myThid :: my thread ID
  REAL (KIND=_RL90), INTENT( IN  ) :: x(2)
  REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, rho
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: c, cimag, gradc, crr, crz, czz, rho

! !LOCAL VARIABLES: None
!EOP

  iSegz = 1                   !RG
  IF ( x(2).LT.Grid%Z( iSegz ) .OR. x(2).GT.Grid%Z( iSegz+1 ) ) THEN
      foundz=.false.
!IEsco23 Test this:
!      DO iz = 2, Grid%nZ   ! Search for bracketting Depths
      DO iz = 2, Grid%nPts   ! Search for bracketting Depths
        IF ( x(2).LT.Grid%Z( iz ) .AND. .NOT.foundz ) THEN
            iSegz  = iz - 1
            foundz = .true.
        ENDIF

      ENDDO

  ENDIF ! IF ( x(2).LT.Grid%Z( iSegz ) .OR. x(2).GT.Grid%Z( iSegz+1 ) )

  c     = REAL(  SSP%c( iSegz ) + ( x(2) - Grid%Z( iSegz ) ) * SSP%cz( iSegz ) )
  cimag = AIMAG( SSP%c( iSegz ) + ( x(2) - Grid%Z( iSegz ) ) * SSP%cz( iSegz ) )
  gradc = [ 0.0D0, REAL( SSP%cz( iSegz ) ) ]
  crr   = 0.0d0
  crz   = 0.0d0
  czz   = 0.0d0

  W     = ( x(2) - Grid%Z( iSegz ) ) / ( Grid%Z( iSegz+1 ) - Grid%Z( iSegz ) )
  rho   = ( 1.0D0-W ) * Grid%rho( iSegz ) + W * Grid%rho( iSegz+1 )

  RETURN
  END !SUBROUTINE cLinear

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: cPCHIP
! !INTERFACE:
  SUBROUTINE cPCHIP( x, c, cimag, gradc, crr, crz, czz, rho, myThid )
! !DESCRIPTION:
!   Piecewise cubic Hermite interpolation (PCHIP) for the sound speed
!   and its derivatives.

! !USES: None

! !INPUT PARAMETERS:
! x      :: r-z SSP evaluation point
! c      :: sound speed at x
! cimag  :: imaginary part of sound speed at x
! gradc  :: gradient of sound speed at x
! crr    :: radial derivative of sound speed at x
! crz    :: radial-z derivative of sound speed at x
! czz    :: z-z derivative of sound speed at x
! rho    :: density at x
! myThid :: my thread ID
  REAL (KIND=_RL90), INTENT( IN  ) :: x(2)
  REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, rho
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: c, cimag, gradc, crr, crz, czz, rho

! !LOCAL VARIABLES:
! xt     :: x(2) - Grid%Z(iSegz)
! c_cmplx :: complex sound speed at x
  REAL    (KIND=_RL90) :: xt
  COMPLEX (KIND=_RL90) :: c_cmplx
!EOP
    
  iSegz = 1                   !RG
  IF ( x(2).LT.Grid%Z( iSegz ) .OR. x(2).GT.Grid%Z( iSegz+1 ) ) THEN
      foundz=.false.
!IEsco23 Test this:
!      DO iz = 2, Grid%nZ   ! Search for bracketting Depths
      DO iz = 2, Grid%nPts   ! Search for bracketting Depths
        IF ( x(2).LT.Grid%Z( iz ) .AND. .NOT.foundz ) THEN
            iSegz  = iz - 1
            foundz = .true.
        ENDIF

      ENDDO

  ENDIF ! IF ( x(2).LT.Grid%Z( iSegz ) .OR. x(2).GT.Grid%Z( iSegz+1 ) )

  xt = x( 2 ) - Grid%Z( iSegz )
  c_cmplx = cCoef( 1, iSegz ) &
        + ( cCoef( 2, iSegz ) &
        + ( cCoef( 3, iSegz ) &
        +   cCoef( 4, iSegz ) * xt ) * xt ) * xt

  c     = REAL(  c_cmplx )
  cimag = AIMAG( c_cmplx )

  gradc = [ 0.0D0, &
            REAL( cCoef( 2, iSegz ) + ( 2.0D0 * cCoef( 3, iSegz ) &
                  + 3.0D0 * cCoef( 4, iSegz ) * xt ) * xt ) ]

  crr   = 0.0D0
  crz   = 0.0D0
  czz   = REAL( 2.0D0 * cCoef( 3, iSegz ) + &
                6.0D0 * cCoef( 4, iSegz ) * xt )   ! dgradc(2)/dxt

  W     = ( x(2) - Grid%Z( iSegz ) ) / ( Grid%Z( iSegz+1 ) - Grid%Z( iSegz ) )
  ! linear interp of density
  rho   = ( 1.0D0-W ) * Grid%rho( iSegz ) + W * Grid%rho( iSegz+1 )

  RETURN
  END !SUBROUTINE cPCHIP

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: cCubic
! !INTERFACE:
  SUBROUTINE cCubic( x, c, cimag, gradc, crr, crz, czz, rho, myThid  )
! !DESCRIPTION:
!   Cubic spline interpolation for the sound speed and its derivatives.

! !USES:
USE splinec_mod,  only: splineall

! !INPUT PARAMETERS:
! x      :: r-z SSP evaluation point
! c      :: sound speed at x
! cimag  :: imaginary part of sound speed at x
! gradc  :: gradient of sound speed at x
! crr    :: radial derivative of sound speed at x
! crz    :: radial-z derivative of sound speed at x
! czz    :: z-z derivative of sound speed at x
! rho    :: density at x
! myThid :: my thread ID
  REAL (KIND=_RL90), INTENT( IN  ) :: x(2)
  REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, rho
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: c, cimag, gradc, crr, crz, czz, rho

! !LOCAL VARIABLES:
! hSpline :: Depth offset from Grid%Z(iSegz)
! c_cmplx, cz_cmplx, czz_cmplx :: Complex sound speed and its derivatives
  REAL     (KIND=_RL90)   :: hSpline
  COMPLEX  (KIND=_RL90)   :: c_cmplx, cz_cmplx, czz_cmplx
!EOP

  ! *** Section to return SSP info ***

  iSegz = 1                   !RG
    IF ( x(2).LT.Grid%Z( iSegz ) .OR. x(2).GT.Grid%Z( iSegz+1 ) ) THEN
      foundz=.false.
!IEsco23 Test this:
!       DO iz = 2, Grid%nZ   ! Search for bracketting Depths
      DO iz = 2, Grid%nPts   ! Search for bracketting Depths
          IF ( x(2).LT.Grid%Z( iz ) .AND. .NOT.foundz ) THEN
            iSegz  = iz - 1
            foundz = .true.
          ENDIF

      ENDDO

    ENDIF ! IF ( x(2).LT.Grid%Z( iSegz ) .OR. x(2).GT.Grid%Z( iSegz+1 ) )

  hSpline = x(2) - Grid%Z( iSegz )

  CALL SplineALL( cSpln( 1, iSegz ), hSpline, c_cmplx, cz_cmplx, czz_cmplx )

  c     = DBLE(  c_cmplx )
  cimag = AIMAG( c_cmplx )
  gradc = [ 0.0D0, DBLE( cz_cmplx ) ]
  czz   = DBLE( czz_cmplx )
  crr   = 0.0d0
  crz   = 0.0d0

  ! linear interpolation for density
  W   = ( x(2) - Grid%Z( iSegz ) ) / ( Grid%Z( iSegz+1 ) - Grid%Z( iSegz ) )
  rho = ( 1.0D0-W ) * Grid%rho( iSegz ) + W * Grid%rho( iSegz+1 )

  RETURN
  END !SUBROUTINE cCubic

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: Quad
! !INTERFACE:
  SUBROUTINE Quad( x, c, cimag, gradc, crr, crz, czz, rho, myThid )
! !DESCRIPTION:
!   Bilinear quadrilateral interpolation of SSP data in 2D, Grid%Type = 'Q'

! !USES: None

! !INPUT PARAMETERS:
! x      :: r-z SSP evaluation point
! c      :: sound speed at x
! cimag  :: imaginary part of sound speed at x
! gradc  :: gradient of sound speed at x
! crr    :: radial derivative of sound speed at x
! crz    :: radial-z derivative of sound speed at x
! czz    :: z-z derivative of sound speed at x
! rho    :: density at x
! myThid :: my thread ID
  REAL (KIND=_RL90), INTENT( IN  ) :: x(2)
  REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, rho
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: c, cimag, gradc, crr, crz, czz, rho

! !LOCAL VARIABLES: None
! msgBuf :: Informational/error message buffer
! irT, iz2, isegzold :: Loop indices
! c1, c2, cz1, cz2 :: Sound speed and its derivatives at the segment ends
! cr, cz, s1, s2 :: Interpolation parameters
! delta_r, delta_z :: Range and depth segment lengths
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
  INTEGER             :: irT, iz2, isegzold
  REAL (KIND=_RL90)   :: c1, c2, cz1, cz2, cr, cz, s1, s2, delta_r, delta_z
!EOP

  ! INIT variables
  c1 = 0.
  c2 = 0.

  ! *** Section to return SSP info ***

  ! IESCO22: iSegz is the depth segment index containing x depth
  ! find depth-layer where x(2) in ( Grid%Z( iSegz ), Grid%Z( iSegz+1 ) )
  iSegz = 1                   !RG
  IF ( x(2).LT.Grid%Z( iSegz ) .OR. x(2).GT.Grid%Z( iSegz+1 ) ) THEN
      foundz=.false.
      DO iz = 2, Grid%nZ   ! Search for bracketting Depths
        IF ( x(2).LT.Grid%Z( iz ) .AND. .NOT.foundz ) THEN
            iSegz  = iz - 1
            foundz = .true.
        ENDIF

      ENDDO

  ENDIF ! IF ( x(2).LT.Grid%Z( iSegz ) .OR. x(2).GT.Grid%Z( iSegz+1 ) )

  ! Check that x is in SSP box range
  IF ( x(1).LT.Grid%Seg%R( 1 ) .OR. x(1).GT.Grid%Seg%R( Grid%nR ) ) THEN
#ifdef IHOP_WRITE_OUT
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) THEN
      WRITE(msgBuf,'(2A)') 'ray is outside the box where ocean ',&
        'soundspeed is defined'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,2F13.4)') ' x = ( r, z ) = ', x
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)') 'SSPMOD Quad: ', &
        'ray is outside the box where the soundspeed is defined'
      CALL PRINT_ERROR( msgBuf,myThid )
    ENDIF
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R Quad'
  ENDIF ! IF ( x(1).LT.Grid%Seg%R( 1 ) .OR. x(1).GT.Grid%Seg%R( Grid%nR ) )

  ! Range segment where x(1) in [ Grid%Seg%R( iSegr ), Grid%Seg%R( iSegr+1 ) )
  iSegr = 1           !RG
  IF ( x(1).LT.Grid%Seg%R( iSegr ) .OR. x(1).GE.Grid%Seg%R( iSegr+1 ) ) THEN
    foundr=.false.
    DO irT = 2, Grid%nR   ! Search for bracketting segment ranges
      IF ( x(1).LT.Grid%Seg%R( irT ) .AND. .NOT.foundr ) THEN
        iSegr = irT - 1
        foundr=.true.
      ENDIF

    ENDDO

  ENDIF ! IF ( x(1).LT.Grid%Seg%R( iSegr ) .OR. x(1).GE.Grid%Seg%R( iSegr+1 ) )

  !IESCO22: s2 is distance btwn field point, x(2), and ssp depth @ iSegz
  s2      = x(2)             - Grid%Z( iSegz )
  delta_z = Grid%Z( iSegz+1 ) - Grid%Z( iSegz )
  IF ( delta_z.LE.0 .OR. s2.GT.delta_z ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf, *) delta_z, s2, iSegz, Grid%Z(iSegz)
    CALL PRINT_ERROR( msgBuf,myThid )
    WRITE(msgBuf,'(2A)') 'SSPMOD Quad: ', &
      'depth is not monotonically increasing in Grid%Z'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R Quad'
  ENDIF ! IF ( delta_z.LE.0 .OR. s2.GT.delta_z )

  ! for depth, x(2), get the sound speed at both ends of range segment
  cz1 = SSP%czMat( iSegz, iSegr   )
  cz2 = SSP%czMat( iSegz, iSegr+1 )

  c1 = SSP%cMat( iSegz, iSegr   ) + s2*cz1
  c2 = SSP%cMat( iSegz, iSegr+1 ) + s2*cz2

  IF ((c1.EQ.0) .OR. (c2.EQ.0)) STOP "ABNORMAL END: S/R QUAD"

  ! s1 = proportional distance of x(1) in range
  delta_r = Grid%Seg%R( iSegr+1 ) - Grid%Seg%R( iSegr )
  s1 = ( x(1) - Grid%Seg%R( iSegr ) ) / delta_r
  ! piecewise constant extrapolation for ranges outside SSP box
  s1 = MIN( s1, 1.0D0 )
  s1 = MAX( s1, 0.0D0 )

  c = ( 1.0D0-s1 )*c1 + s1*c2 ! c @ x

  ! interpolate attenuation !!!! SSP in ENVFile needs to match first column of SSPFile
  s2    = s2 / delta_z   ! normalize depth layer
  ! volume attenuation is taken from the single c(z) profile
  cimag = AIMAG( ( 1.0D0-s2 )*SSP%c( iSegz ) + s2*SSP%c( iSegz+1 ) )

  cz  = ( 1.0D0-s1 )*cz1 + s1*cz2 ! cz @ x

  cr  = ( c2  - c1  ) / delta_r ! SSP grid cr
  crz = ( cz2 - cz1 ) / delta_r ! SSP grid crz

  gradc = [ cr, cz ]
  crr   = 0.0
  czz   = 0.0

  ! linear interpolation for density
  W   = ( x(2) - Grid%Z( iSegz ) ) / ( Grid%Z( iSegz+1 ) - Grid%Z( iSegz ) )
  rho = ( 1.0D0-W ) * Grid%rho( iSegz ) + W * Grid%rho( iSegz+1 )

  !IESCO22: for thesis, czz=crr=0, and rho=1 at all times
  RETURN
  END !SUBROUTINE Quad

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: gcmSSP
! !INTERFACE:
  SUBROUTINE gcmSSP( myThid )
! !DESCRIPTION:
!   Interpolate SSP from GCM grid to iHOP grid using adaptive IDW interpolation.

! !USES:
  USE atten_mod, only: CRCI
  USE bdry_mod, only: Bdry
! == Global variables ==
#ifdef ALLOW_AUTODIFF
#include "tamc.h"
#endif

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT(IN)     :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! fmtstr :: Format string for output messages
! interp_finished :: is adaptive IDW interp finished at this range point?
! iallocstat :: Allocation status
! bi, bj, i, j, k, ii, jj :: GCM domain loop indices
! njj :: number interpolation points per depth level
! dcdz, tolerance :: Tolerance for IDW interpolation
! tileSSP, tmpSSP, globSSP :: Arrays for storing interpolated SSP
! bPower, fT :: Parameters for attenuation calculation
! tkey, ijkey, hkey, lockey :: Keys for storing TAF tapes
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  CHARACTER*(80)          :: fmtstr
  LOGICAL :: interp_finished
  INTEGER :: iallocstat
  INTEGER :: bi,bj, i,j,k, ii,jj
  INTEGER :: njj(IHOP_NPTS_RANGE)
  INTEGER :: nZnR_size
  REAL (KIND=_RL90)             :: dcdz, tolerance
  REAL (KIND=_RL90), ALLOCATABLE:: tileSSP(:,:,:,:), tmpSSP(:,:,:), globSSP(:)
! IESCO24
! fT = 1000 ONLY for acousto-elastic halfspaces, I will have to pass this
! parameter in a different way after ssp_mod is split btwn fixed and varia
  REAL (KIND=_RL90)             :: bPower, fT
#ifdef ALLOW_AUTODIFF_TAMC
  INTEGER tkey, ijkey, hkey, lockey
!$TAF init loctape_ihop_gcmssp_bibj_ij_iijj_k = STATIC, nSx*nSy*sNx*sNy*IHOP_MAX_RANGE*IHOP_MAX_NC_SIZE*(Nr + 2)
#endif

  ! IESCO24 fT init
  bPower = 1.0
  fT = 1000.0

  ! init local vars
  interp_finished =.false.
  njj       = 0
  dcdz      = 0.0 _d 0
  tolerance = 5 _d -5
  nZnR_size = Grid%nZ*Grid%nR

  IF(ALLOCATED(tileSSP)) DEALLOCATE(tileSSP)
  IF(ALLOCATED(globSSP)) DEALLOCATE(globSSP)
  ALLOCATE( tileSSP(Grid%nZ,Grid%nR,nSx,nSy), &
            tmpSSP(nSx,nSy,nZnR_size), &
            globSSP(nZnR_size), STAT=iallocstat )
  IF ( iallocstat.NE.0 ) THEN
# ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SSPMOD gcmSSP: ', &
      'Insufficient memory to store tileSSP and/or globSSP'
    CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R gcmSSP'
  ENDIF

  ! Initiate to ceros
  tileSSP    = 0.0 _d 0
  tmpSSP     = 0.0 _d 0
  globSSP    = 0.0 _d 0

  ! interpolate SSP with adaptive IDW from gcm grid to ihop grid
  DO bj=myByLo(myThid),myByHi(myThid)
    DO bi=myBxLo(myThid),myBxHi(myThid)
#ifdef ALLOW_AUTODIFF_TAMC
      tkey = bi + (bj-1)*nSx + (ikey_dynamics-1)*nSx*nSy
#endif

! IESCO24: don't worry about gcm overlaps right now
      DO j=1,sNy
        DO i=1,sNx
#ifdef ALLOW_AUTODIFF_TAMC
          ijkey = i + (j-1)*sNx + (tkey-1)*sNx*sNy
!$TAF store njj(ii) = comlev1_bibj_ij_ihop, key=ijkey, kind=isbyte
#endif
          DO ii=1,IHOP_npts_range
            interp_finished = .FALSE.
            DO jj=1,IHOP_npts_idw
#ifdef ALLOW_AUTODIFF_TAMC
!$TAF STORE interp_finished = comlev1_bibj_ij_ihop, key=ijkey, kind=isbyte
#endif
              ! Interpolate from GCM grid cell centers
              IF ( ABS(xC(i, j, bi, bj) - ihop_xc(ii, jj)).LE.tolerance .AND. &
                   ABS(yC(i, j, bi, bj) - ihop_yc(ii, jj)).LE.tolerance .AND. &
                  .NOT.interp_finished ) THEN
                njj(ii) = njj(ii) + 1

                DO iz = 1, Grid%nZ - 1
#ifdef ALLOW_AUTODIFF_TAMC
                  hkey = jj + (ii-1)*IHOP_npts_idw + (ijkey-1)*sNy*sNx*nSy*nSx
                  lockey = iz + (hkey-1)*(Grid%nZ-1)*IHOP_npts_idw*IHOP_npts_range*sNx*sNy*nSx*nSy
!                  + ((jj-1) + ((ii-1) + (ijkey-1)*IHOP_npts_range)*IHOP_npts_idw)*(Grid%nZ-1)
!$TAF store njj(ii) = loctape_ihop_gcmssp_bibj_ij_iijj_k, key=lockey, kind=isbyte
#endif

    IF ( iz.EQ.1 ) THEN
      ! Top vlevel zero depth
      tileSSP(1, ii, bi, bj) = tileSSP(1, ii, bi, bj) + &
        CHEN_MILLERO(i, j, 0, bi, bj, myThid) * &
        ihop_idw_weights(ii, jj) / ihop_sumweights(ii, iz)

    ELSE ! 2:(Grid%nZ-1)
      ! Middle depth layers, only when not already underground
      IF ( ihop_sumweights(ii, iz-1).GT.0. ) THEN
        ! isolate njj increments for TAF 
        IF ( ihop_idw_weights(ii, jj).EQ.0. ) njj(ii) = IHOP_npts_idw + 1

        ! Exactly on a cell center, ignore interpolation
        IF ( ihop_idw_weights(ii, jj).EQ.0. ) THEN
          tileSSP(iz, ii, bi, bj) = ihop_ssp(i, j, iz-1, bi, bj)

        ! Apply IDW interpolation
        ELSEIF ( njj(ii).LE.IHOP_npts_idw ) THEN
          tileSSP(iz, ii, bi, bj) = tileSSP(iz, ii, bi, bj) + &
            ihop_ssp(i, j, iz-1, bi, bj) * &
            ihop_idw_weights(ii, jj) / ihop_sumweights(ii, iz-1)

        ELSE
          ! do nothing TAF NEEDS THIS LINE
          tileSSP(iz, ii, bi, bj) = tileSSP(iz, ii, bi, bj)

        ENDIF ! IF ( ihop_idw_weights(ii, jj).EQ.0. )

      ENDIF ! IF ( ihop_sumweights(ii, iz-1).GT.0. )

      ! Extrapolate through bathymetry; don't interpolate
      IF ( iz.EQ.Grid%nZ-1 .OR. ihop_sumweights(ii, iz-1).EQ.0.0 ) THEN
        k = iz

        IF ( njj(ii).GE.IHOP_npts_idw ) THEN
          ! Determine if you are at the last vlevel
          IF ( iz.EQ.Grid%nZ-1 .AND. ihop_sumweights(ii, iz-1).NE.0. ) &
            k = k + 1

          ! Calc depth gradient
          dcdz = (tileSSP(k-1, ii, bi, bj) - tileSSP(k-2, ii, bi, bj)) / &
                 (Grid%Z(k-1) - Grid%Z(k-2))
          ! Extrapolate
          tileSSP(k:Grid%nZ, ii, bi, bj) = &
            tileSSP(k-1, ii, bi, bj) + dcdz * Grid%Z(k:Grid%nZ)
          ! Move to next range point, ii
          interp_finished = .TRUE.

        ELSE
          interp_finished = .FALSE.

        ENDIF ! IF ( njj(ii).GE.IHOP_npts_idw )

      ENDIF ! IF ( iz.EQ.Grid%nZ-1 .OR. ihop_sumweights(ii, iz-1).EQ.0.0 )

    ENDIF ! IF ( iz.EQ.1 )
                ENDDO !iz
              ENDIF ! IF ( ABS(xC(i, j, bi, bj) - ihop_xc(ii, jj)).LE.tolerance .AND. &
                   ! ABS(yC(i, j, bi, bj) - ihop_yc(ii, jj)).LE.tolerance .AND. &
                   ! .NOT.interp_finished )
            ENDDO !jj
          ENDDO !ii
        ENDDO !i
      ENDDO !j
    ENDDO !bi
  ENDDO !bj

! IESCO24: MITgcm checkpoint69a uses a new global sum subroutine...
  DO bj=myByLo(myThid),myByHi(myThid)
    DO bi=myBxLo(myThid),myBxHi(myThid)
      k = 1
      DO jj = 1,Grid%nR
        DO ii = 1,Grid%nZ
          tmpSSP(bi,bj,k) = tileSSP(ii,jj,bi,bj)
          k = k + 1
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  CALL GLOBAL_SUM_VECTOR_RL(nZnR_size,tmpSSP,globSSP,myThid)

  ! reshape ssp to matrix
  k = 1
  DO jj = 1,Grid%nR
    DO ii = 1,Grid%nZ
      SSP%cMat(ii,jj) = globSSP(k)
      k = k + 1
    ENDDO
  ENDDO

  IF(ALLOCATED(tileSSP)) DEALLOCATE(tileSSP)
  IF(ALLOCATED(tmpSSP))  DEALLOCATE(tmpSSP)
  IF(ALLOCATED(globSSP)) DEALLOCATE(globSSP)

  !==================================================
  ! END IDW Interpolate
  !==================================================

  ! set vector structured c, rho, and cz for first range point
  IF ( .NOT.useSSPFile ) THEN ! if usessp, these have already been set
    DO iz = 1,Grid%nZ
      alphaR = SSP%cMat( iz, 1 )

      SSP%c(iz) = CRCI( Grid%Z(iz), alphaR, alphaI, Grid%AttenUnit, bPower, fT, &
                          myThid )
      Grid%rho(iz) = rhoR

      IF ( iz.GT.1 ) THEN
        IF ( Grid%Z( iz ).LE.Grid%Z( iz-1 ) ) THEN
#ifdef IHOP_WRITE_OUT
          WRITE( msgBuf,'(2A)' ) 'SSPMOD gcmSSP: ', &
            'The depths in the SSP must be monotone increasing'
          CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
          STOP 'ABNORMAL END: S/R gcmSSP'
        ENDIF
      ENDIF

      ! Compute gradient, cz
      IF ( iz.GT.1 ) SSP%cz( iz-1 ) = ( SSP%c( iz ) - SSP%c( iz-1 ) ) / &
                                   ( Grid%Z( iz ) - Grid%Z( iz-1 ) )

    ENDDO ! DO iz = 1,Grid%nZ

  ENDIF ! IF ( .NOT.useSSPFile )

RETURN
END !SUBROUTINE gcmSSP

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: writeSSP
! !INTERFACE:
  SUBROUTINE writeSSP( myThid )
! !DESCRIPTION:
!   Write the sound speed profile (SSP) to the output file.

! !USES:
  USE ihop_mod, only: PRTFile

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! fmtstr :: Format string for output messages
! ssptmp :: Temporary array for sound speed values
! iz :: Loop index for depth levels
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  CHARACTER*(80)    :: fmtstr
  REAL (KIND=_RL90) :: ssptmp(Grid%nR)
  INTEGER           :: iz

  ! In adjoint mode we do not write output besides on the first run
  IF ( IHOP_dumpfreq.LT.0 ) RETURN

  ! init local vars
  ssptmp = 0.0

  ! I/O on main thread only
  _BEGIN_MASTER(myThid)

#ifdef IHOP_WRITE_OUT
  ! Write relevant diagnostics
  WRITE(msgBuf,'(A)') &
    '___________________________________________________________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') "Sound Speed Field"
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  WRITE(msgBuf,'(2A)') 'Profile option: ', Grid%Type
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  IF ( Grid%nR.GT.1 ) THEN
    WRITE(msgBuf,'(A)') 'Using range-dependent sound speed'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  END IF
  IF ( Grid%nR.EQ.1 ) THEN
    WRITE(msgBuf,'(A)') 'Using range-independent sound speed'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  END IF

  WRITE(msgBuf,'(A,I10)') 'Number of SSP ranges = ', Grid%nR
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A,I10)') 'Number of SSP depths = ', Grid%nZ
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') 'Profile ranges [km]:'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(fmtStr,'(A,I10,A)') '(T11,',Grid%nR, 'F10.2)'
  ssptmp = Grid%Seg%R( 1:Grid%nR ) / 1000.0
  WRITE(msgBuf,fmtStr) ssptmp
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') 'Sound speed matrix:'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') ' Depth [m ]     Soundspeed [m/s]'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  ssptmp = 0.0
  WRITE(fmtStr,'(A,I10,A)') '(',Grid%nR+1, 'F10.2)'
  DO iz = 1, Grid%nZ
    ssptmp = ssp%cMat( iz,: )
    WRITE(msgBuf,fmtStr) Grid%Z( iz ), ssptmp
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  ENDDO

  IF ( useSSPFile ) THEN
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 'Sound speed profile:'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)')'      z         alphaR      betaR     rho      ', &
                        '  alphaI     betaI'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)')'     [m]         [m/s]      [m/s]   [g/cm^3]   ', &
                        '   [m/s]     [m/s]'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    WRITE(msgBuf,'(A)') &
      '___________________________________________________________'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    DO iz = 1, Grid%nPts
       WRITE(msgBuf,'( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4)') &
          Grid%Z( iz ), SSP%cMat(iz,1), betaR, rhoR, alphaI, betaI
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    ENDDO

  ENDIF ! IF ( useSSPFile )

  WRITE(msgBuf,'(A)') &
    '___________________________________________________________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

  ! I/O on main thread only
  _END_MASTER(myThid)
  _BARRIER

  RETURN
  END !SUBROUTINE writeSSP

END MODULE ssp_mod
