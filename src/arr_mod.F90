#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: arr_mod
MODULE arr_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   Contains the subroutines to read the ray arrival angles

! !USES:
  USE ihop_mod,   only: rad2deg, ARRFile
  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_USE_MPI
# include "EESUPPORT.h"
#endif
#include "IHOP_SIZE.h"
#include "IHOP.h"

! !SCOPE: 
  PRIVATE
!=======================================================================
    PUBLIC WriteArrivalsASCII, WriteArrivalsBinary, &
           initArr, nMaxArr, nArrival, Arr, AddArr, U
#ifdef IHOP_THREED
    PUBLIC nArr3D, Arr3D
#endif /* IHOP_THREED */
#ifdef ALLOW_USE_MPI
    PUBLIC free_ihop_arrival, BcastArr
#endif /* ALLOW_USE_MPI */
!=======================================================================

! == Module variables ==
  INTEGER               :: nMaxArr
  INTEGER, ALLOCATABLE  :: nArrival( :, : )
  COMPLEX, ALLOCATABLE  :: U( :, : )
#ifdef IHOP_THREED
  INTEGER, ALLOCATABLE  :: nArr3D( :, :, : )
#endif /* IHOP_THREED */

#ifdef ALLOW_USE_MPI
  INTEGER :: MPI_IHOP_ARRIVAL = MPI_DATATYPE_NULL
  LOGICAL :: ARRIVAL_TYPE_COMMITTED = .false.
#endif /* ALLOW_USE_MPI */

! == Derived types ==
  TYPE Arrival
    INTEGER              :: nTopBnc, nBotBnc
    REAL (KIND=_RL90)    :: SrcDeclAngle, RcvrDeclAngle, A, Phase
#ifdef IHOP_THREED
    REAL (KIND=_RL90)    :: SrcAzimAngle, RcvrAzimAngle
#endif /* IHOP_THREED */
    COMPLEX (KIND=_RL90) :: delay
  END TYPE

  TYPE(Arrival), ALLOCATABLE :: Arr( :, :, : )
#ifdef IHOP_THREED
  TYPE(Arrival), ALLOCATABLE :: Arr3D( :, :, :, : )
#endif /* IHOP_THREED */
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R initArr
! S/R AddArr
! S/R WriteArrivalsASCII
! S/R WriteArrivalsBinary
! S/R BcastArr
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: initArr
! !INTERFACE:
! **************************************************************************** !
  SUBROUTINE initArr( myThid )
! !DESCRIPTION:
! Initializes the arrival and U variables

! !USES:
  USE srPos_mod, only: Pos
  USE ihop_mod,  only: Beam, nRz_per_range

! !INPUT PARAMETERS:
! myThid  :: my thread ID
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf              :: Informational/error message buffer
! iAllocStat          :: Allocation status
! x, y                :: Dimensions of the U matrix
! arrStorage          :: Storage for arrivals
! minnArr             :: Minimum number of arrivals to allocate
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER             :: iAllocStat
  INTEGER             :: x, y
  INTEGER, PARAMETER  :: arrStorage = 100, minnArr = 10
!EOP

  ! reset memory
  IF (ALLOCATED(U))           DEALLOCATE(U)
  IF (ALLOCATED(Arr))         DEALLOCATE(Arr)
  IF (ALLOCATED(nArrival))    DEALLOCATE(nArrival)


  SELECT CASE ( Beam%RunType( 1:1 ) )
    CASE ( 'C', 'S', 'I' ) ! TL calculation
      ! Allocate space for the pressure matrix
      x = nRz_per_range
      y = Pos%nRR
    CASE ( 'A', 'a', 'R', 'E', 'e' )  ! Arrivals calculation
      x = 1
      y = 1
    CASE DEFAULT
      x = 1
      y = 1
  END SELECT ! Beam%RunType( 1:1 )

  ALLOCATE( U( x,y ), Stat=iAllocStat )
  IF ( iAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'ARR_MOD INITARR: ', &
      'Insufficient memory for TL matrix: reduce Nr*NRz'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R INITARR'
  ENDIF

  ! for an arrivals run, allocate space for arrivals matrices
  SELECT CASE ( Beam%RunType( 1:1 ) )
  CASE ( 'A', 'a', 'e' )
    ! allow space for at least minnArr arrivals
    nMaxArr = MAX( arrStorage / ( nRz_per_range * Pos%nRR ), &
                  minnArr )
  CASE DEFAULT
    nMaxArr = 1
  END SELECT ! Beam%RunType( 1:1 )

  ! init Arr, nArrival
  ALLOCATE( Arr( nMaxArr, Pos%nRR, nRz_per_range ), &
            nArrival(Pos%nRR, nRz_per_range), STAT=iAllocStat )
  IF ( iAllocStat.NE.0 ) THEN
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'ARR_MOD INITARR: ', &
      'Not enough allocation for Arr; reduce arrStorage'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R IHOP_MAIN'
  ENDIF

#ifdef ALLOW_USE_MPI
! create MPI DerivedType
  CALL ArrivalTypeInit( Arr(1,1,1) )

#endif /* ALLOW_USE_MPI */
  ! init default values
  U                         = 0.0
  nArrival(:,:)             = 0
  Arr(:,:,:)%NTopBnc        = -1
  Arr(:,:,:)%NBotBnc        = -1
  Arr(:,:,:)%SrcDeclAngle   = -999.
  Arr(:,:,:)%RcvrDeclAngle  = -999.
  Arr(:,:,:)%A              = -999.
  Arr(:,:,:)%Phase          = -999.
  Arr(:,:,:)%delay          = -999.

  END !SUBROUTINE initArr( myThid )

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: AddArr
! !INTERFACE:
  SUBROUTINE AddArr( afreq, iz, ir, Amp, Phase, delay, &
                     RcvrDeclAngle, nTopBnc, nBotBnc )
! !DESCRIPTION:
! Adds the amplitude and delay for an arrival into Arr.

! !USES:
  USE ihop_mod, only: SrcDeclAngle


! !INPUT PARAMETERS:
! afreq         :: frequency of the arrival
! iz            :: receiver depth index
! ir            :: range index
! Amp           :: amplitude of the arrival
! Phase         :: phase of the arrival
! delay         :: delay time of the arrival
! RcvrDeclAngle :: receiver declination angle
! nTopBnc       :: number of top bounces
! nBotBnc       :: number of bottom bounces
  INTEGER,              INTENT( IN ) :: iz, ir, nTopBnc, nBotBnc
  REAL    (KIND=_RL90), INTENT( IN ) :: afreq, Amp, Phase, &
                                        RcvrDeclAngle
  COMPLEX (KIND=_RL90), INTENT( IN ) :: delay
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
  REAL, PARAMETER :: PhaseTol = 0.05
  LOGICAL         :: NewRay
  INTEGER         :: iArr( 1 ), Nt
  REAL            :: AmpTot, w1, w2
!EOP

  Nt     = nArrival( ir, iz )    ! # of arrivals
  NewRay = .TRUE.

  ! Is this the second bracketting ray of a pair?
  ! If so, we want to combine the arrivals to conserve space.
  ! (test this by seeing if the arrival time is close to the previous one)
  ! (also need that the phase is about the same to make sure surface and
  ! direct paths are not joined)

  IF ( Nt.GE.1 ) THEN
    IF ( afreq * ABS( delay-Arr( Nt, ir, iz )%delay ).LT.PhaseTol &
       .AND. ABS( Arr( Nt, ir, iz )%phase-Phase ).LT.PhaseTol ) &
      NewRay = .FALSE.
  ENDIF

  IF ( NewRay ) THEN
    IF ( Nt.GE.nMaxArr ) THEN       ! space available to add an arrival?
      iArr = MINLOC( Arr( :, ir, iz )%A )   ! no: replace weakest arrival
      IF ( Amp.GT.Arr( iArr(1), ir, iz )%A ) THEN
        Arr( iArr(1), ir, iz)%A             = SNGL( Amp )       ! amplitude
        Arr( iArr(1), ir, iz)%Phase         = SNGL( Phase )     ! phase
        Arr( iArr(1), ir, iz)%delay         = CMPLX( delay )    ! delay time
        Arr( iArr(1), ir, iz)%SrcDeclAngle  = SNGL( SrcDeclAngle) ! angle
        Arr( iArr(1), ir, iz)%RcvrDeclAngle = SNGL(RcvrDeclAngle) ! angle
        Arr( iArr(1), ir, iz)%nTopBnc       = nTopBnc         ! top bounces
        Arr( iArr(1), ir, iz)%nBotBnc       = nBotBnc         ! bottom bounces
      ENDIF
    ELSE
      Nt                             = Nt+1                ! # arrivals
      Arr( Nt, ir, iz)%A             = SNGL( Amp )         ! amplitude
      Arr( Nt, ir, iz)%Phase         = SNGL( Phase )       ! phase
      Arr( Nt, ir, iz)%delay         = CMPLX( delay )      ! delay time
      Arr( Nt, ir, iz)%SrcDeclAngle  = SNGL( SrcDeclAngle )  ! angle
      Arr( Nt, ir, iz)%RcvrDeclAngle = SNGL( RcvrDeclAngle ) ! angle
      Arr( Nt, ir, iz)%nTopBnc       = nTopBnc           ! top bounces
      Arr( Nt, ir, iz)%nBotBnc       = nBotBnc           ! bottom bounces
    ENDIF !IF ( Nt.GE.nMaxArr )

  ELSE ! not a new ray
    ! calculate weights of old ray information vs. new
    AmpTot = Arr( Nt, ir, iz )%A + SNGL( Amp )
    w1     = Arr( Nt, ir, iz )%A / AmpTot
    w2     = REAL( Amp ) / AmpTot

    Arr( Nt, ir, iz)%delay         =  w1 * Arr( Nt, ir, iz )%delay &
                                    + w2 * CMPLX( delay ) ! weighted sum
    Arr( Nt, ir, iz)%A             =  AmpTot
    Arr( Nt, ir, iz)%SrcDeclAngle  =  w1 * Arr( Nt, ir, iz )%SrcDeclAngle &
                                    + w2 * SNGL( SrcDeclAngle  )
    Arr( Nt, ir, iz)%RcvrDeclAngle =  w1 * Arr( Nt, ir, iz )%RcvrDeclAngle &
                                    + w2 * SNGL( RcvrDeclAngle )
  ENDIF ! IF ( NewRay )

  ! Pass Nt to global Narrival
  nArrival(ir,iz) = Nt

  RETURN
  END !SUBROUTINE AddArr

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteArrivalsASCII
! !INTERFACE:
  SUBROUTINE WriteArrivalsASCII( r, Nrz, nRR, SourceType )
! !DESCRIPTION:
! Writes arrival data, eg. Amplitude, delay, for each eigenray in ASCII file

! !USES: None

! !INPUT PARAMETERS:
! r                  :: range vector
! Nrz                :: number of receiver depths
! nRR                :: number of ranges
! SourceType         :: type of source (point or line)
  REAL (KIND=_RL90), INTENT( IN ) :: r( Nr )
  INTEGER,           INTENT( IN ) :: Nrz, nRR
  CHARACTER (LEN=1), INTENT( IN ) :: SourceType
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! ir, iz, iArr :: indices for range and depth
! factor       :: factor for amplitude scaling
! arrFMT      :: format for writing arrivals
  INTEGER           :: ir, iz, iArr
  REAL (KIND=_RL90) :: factor
  CHARACTER*(36)    :: arrFMT
!EOP

  ! In adjoint mode we do not write output besides on the first run
  IF (IHOP_dumpfreq.LT.0) RETURN

  arrFMT='(G14.6,F10.2,F12.4,G10.2,2F14.6,2I6)'

#ifdef IHOP_WRITE_OUT
  WRITE( ARRFile,'(I0)' ) MAXVAL( nArrival( 1:nRR, 1:Nrz ) )
#endif /* IHOP_WRITE_OUT */

  DO iz = 1, Nrz
    DO ir = 1, nRR
      IF ( SourceType.EQ.'X' ) THEN   ! line source
        factor =  4.0 * SQRT( PI )
      ELSE                            ! point source: default
        IF ( r ( ir ).EQ.0 ) THEN
          factor = 1e5                   ! avoid /0 at origin
        ELSE
          factor = 1. / SQRT( r( ir ) )  ! cyl. spreading
        ENDIF

      ENDIF ! IF ( SourceType.EQ.'X' )

#ifdef IHOP_WRITE_OUT
      WRITE( ARRFile, '(I0)' ) nArrival( ir, iz )
      DO iArr = 1, nArrival( ir, iz )
        WRITE( ARRFile, arrFMT ) &
          SNGL( factor )  *Arr( iArr, ir, iz )%A,     &
          SNGL( rad2deg ) *Arr( iArr, ir, iz )%Phase, &
          REAL( Arr( iArr, ir, iz )%delay ),          &
          AIMAG( Arr( iArr, ir, iz )%delay ),         &
          Arr( iArr, ir, iz )%SrcDeclAngle,           &
          Arr( iArr, ir, iz )%RcvrDeclAngle,          &
          Arr( iArr, ir, iz )%NTopBnc,                &
          Arr( iArr, ir, iz )%NBotBnc

      ENDDO  ! next arrival

#endif /* IHOP_WRITE_OUT */
    ENDDO  ! DO ir: next range
  ENDDO  ! DO iz: next depth

  RETURN
  END !SUBROUTINE WriteArrivalsASCII

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteArrivalsBinary
! !INTERFACE:
  SUBROUTINE WriteArrivalsBinary( r, Nrz, nRR, SourceType )
! !DESCRIPTION:
! Writes the arrival data, eg. Amplitude, delay, for each eigenray in binary

! !USES: None

! !INPUT PARAMETERS:
! r                  :: range vector
! Nrz                :: number of receiver depths
! nRR                :: number of ranges
! SourceType         :: type of source (point or line)
  REAL (KIND=_RL90), INTENT( IN ) :: r( Nr )
  INTEGER,           INTENT( IN ) :: Nrz, nRR
  CHARACTER*(1),     INTENT( IN ) :: SourceType ! Beam%RunType(4:4)
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! ir, iz, iArr :: indices for range and depth
! factor       :: factor for amplitude scaling
  INTEGER           :: ir, iz, iArr
  REAL (KIND=_RL90) :: factor
!EOP  

  ! In adjoint mode we do not write output besides on the first run
  IF (IHOP_dumpfreq.LT.0) RETURN

#ifdef IHOP_WRITE_OUT
  WRITE( ARRFile, '(I0)' ) MAXVAL( nArrival( 1:nRR, 1:Nrz ) )
#endif /* IHOP_WRITE_OUT */

  DO iz = 1, Nrz
    DO ir = 1, nRR
      IF ( SourceType.EQ.'X' ) THEN   ! line source
        factor = 4.0 * SQRT( PI )
      ELSE                            ! point source
        IF ( r ( ir ).EQ.0 ) THEN
          factor = 1e5                   ! avoid /0 at origin
        ELSE
          factor = 1. / SQRT( r( ir ) )  ! cyl. spreading
        ENDIF
      ENDIF ! IF ( SourceType.EQ.'X' )

#ifdef IHOP_WRITE_OUT
      WRITE( ARRFile, '(I0)' ) nArrival( ir, iz )
      DO iArr = 1, nArrival( ir, iz )
        ! integers written out as reals below for fast reading in Matlab
        WRITE( ARRFile ) &
          SNGL( factor * Arr( iArr, ir, iz )%A ),      &
          SNGL( rad2deg * Arr( iArr, ir, iz )%Phase ), &
          Arr( iArr, ir, iz )%delay,                   &
          Arr( iArr, ir, iz )%SrcDeclAngle,            &
          Arr( iArr, ir, iz )%RcvrDeclAngle,           &
          REAL( Arr( iArr, ir, iz )%NTopBnc ),         &
          REAL( Arr( iArr, ir, iz )%NBotBnc )

      ENDDO   ! next arrival
#endif /* IHOP_WRITE_OUT */
    ENDDO ! DO ir: next range
  ENDDO ! DO iz: next depth

  RETURN
  END !SUBROUTINE WriteArrivalsBinary

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: BcastArr
! !INTERFACE:
  SUBROUTINE BcastArr( root, comm )
! !DESCRIPTION:
! Broadcasts the arrival data, eg. Amplitude, delay, to all MPI ranks

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
  IF ( MPI_IHOP_ARRIVAL.EQ.MPI_DATATYPE_NULL ) THEN
    CALL ArrivalTypeInit( Arr(1,1,1) )
  ENDIF

  ! We are on MPI rank 0
  arrSize = SIZE(nArrival)
  CALL MPI_Bcast( nArrival, arrSize, MPI_INTEGER, root, comm, ierr )

  ! Broadcast MPI Arrival to all ranks, and free storage
  arrSize = SIZE(Arr)
  CALL MPI_Bcast( Arr, arrSize, MPI_IHOP_ARRIVAL, root, comm, ierr )
#endif /* ALLOW_USE_MPI */

  RETURN
  END !SUBROUTINE BcastArr

#ifdef ALLOW_USE_MPI
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: free_ihop_arrival
! !INTERFACE:
  SUBROUTINE free_ihop_arrival(myThid)
! !DESCRIPTION:
! Free arrival datatype
    INTEGER, INTENT( IN ) :: myThid
    INTEGER :: ierr

    IF (ARRIVAL_TYPE_COMMITTED) THEN
      CALL MPI_Type_free(MPI_IHOP_ARRIVAL, ierr)
      MPI_IHOP_ARRIVAL = MPI_DATATYPE_NULL
      ARRIVAL_TYPE_COMMITTED = .false.
    ENDIF

  RETURN
  END !SUBROUTINE free_ihop_arrival

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ArrivalTypeInti
! !INTERFACE:
  SUBROUTINE ArrivalTypeInit( singleArrival )
! !DESCRIPTION:
! Initializes the arrival MPI Datatype

! !USES: None

! !INPUT PARAMETERS:
  TYPE( Arrival ), INTENT( IN ) :: singleArrival
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! disp       :: address for new MPI datatype Arrival parameter
! base, addr :: base and address for current Arrival datatype and parameters
! ty         :: datatype of each Arrival parameter
! MPI_RL, MPI_CL :: MPI type depending on _RL90
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(7), base, addr(7)
  INTEGER :: ierr, n, bl(7), ty(7)
  INTEGER :: MPI_RL, MPI_CL
!EOP

  IF (ARRIVAL_TYPE_COMMITTED) RETURN

  ! Build Arrival MPI type
  IF (STORAGE_SIZE( singleArrival%Phase ).EQ.64) THEN
    MPI_RL = MPI_DOUBLE_PRECISION
    MPI_CL = MPI_DOUBLE_COMPLEX
  ELSEIF (STORAGE_SIZE( singleArrival%Phase ).EQ.32) THEN
    MPI_RL = MPI_REAL
    MPI_CL = MPI_COMPLEX
  ELSE
    STOP "ABNORMAL END MPI_DATATYPE: Unsupported _RL90 size for MPI"
  ENDIF

  ! Initiate all bl to 1 since all Arr parameters are scalars
  bl=1
  CALL MPI_Get_address(singleArrival, base, ierr)

  n=1
  CALL MPI_Get_address(singleArrival%nTopBnc, addr(n), ierr)
  ty(n)=MPI_INTEGER

  n=n+1
  CALL MPI_Get_address(singleArrival%nBotBnc, addr(n), ierr)
  ty(n)=MPI_INTEGER

  n=n+1
  CALL MPI_Get_address(singleArrival%SrcDeclAngle, addr(n), ierr)
  ty(n)=MPI_RL

  n=n+1
  CALL MPI_Get_address(singleArrival%RcvrDeclAngle, addr(n), ierr)
  ty(n)=MPI_RL

  n=n+1
  CALL MPI_Get_address(singleArrival%A, addr(n), ierr)
  ty(n)=MPI_RL

  n=n+1
  CALL MPI_Get_address(singleArrival%Phase, addr(n), ierr)
  ty(n)=MPI_RL

  n=n+1
  CALL MPI_Get_address(singleArrival%delay, addr(n), ierr)
  ty(n)=MPI_CL

  disp(1:n)=addr(1:n)-base

  CALL MPI_Type_create_struct(n, bl, disp, ty, MPI_IHOP_ARRIVAL, ierr)
  CALL MPI_Type_commit(MPI_IHOP_ARRIVAL, ierr)
  ARRIVAL_TYPE_COMMITTED = .true.

RETURN
END
#endif /* ALLOW_USE_MPI */

END !MODULE arr_mod
