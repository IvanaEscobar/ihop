#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: influence
MODULE influence
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   This module contains routines for computing the influence of a beam on
!   the pressure field.

! !USES:
  USE ihop_mod,  only: rad2deg, oneCMPLX, PRTFile, maxN, Beam, ray2D, &
                       SrcDeclAngle, nRz_per_range
  USE srPos_mod, only: Pos
  IMPLICIT NONE
!  == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC InfluenceGeoHatRayCen, InfluenceGeoGaussianCart, &
           InfluenceGeoHatCart, ScalePressure
!=======================================================================

! == Module variables ==
  INTEGER,              PRIVATE :: iz, ir, iS
  REAL    (KIND=_RL90), PRIVATE :: Ratio1 = 1.0D0 ! scale factor for a line source
  REAL    (KIND=_RL90), PRIVATE :: W, s, n, Amp, phase, phaseInt, &
                                   q0, q, qold, RcvrDeclAngle, rA, rB
  COMPLEX (KIND=_RL90), PRIVATE :: delay
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R InfluenceGeoHatRayCen
! S/R InfluenceGeoHatCart
! S/R InfluenceGeoGaussianCart
! S/R ApplyContribution
! S/R ScalePressure
! FXN Hermite
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: InfluenceGeoHatRayCen
! !INTERFACE:
  SUBROUTINE InfluenceGeoHatRayCen( U, myThid )
! !DESCRIPTION:
!   Geometrically-spreading beams with a hat-shaped beam in ray-centered
!   coordinates

! !USES:
!  USE arr_mod, only: nArrival, Arr         !RG
  USE angle_mod, only: Angles

! !INPUT PARAMETERS:
! U :: complex pressure field
! myThid :: my thread ID
  COMPLEX, INTENT( INOUT ) :: U( nRz_per_range, Pos%nRR )
  INTEGER, INTENT( IN )    :: myThid
! !OUTPUT PARAMETERS: U

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! irA, irB :: indices of receivers
! II :: index for stepping through receivers
! nA, nB :: normalized distances to ray
! zr, L :: ray-centered coordinates
! dq :: step length along ray
! RcvrDeclAngleV :: vector of receiver declination angles
! dtau :: delay along ray
! skip_step :: logical to skip steps with no influence
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
  INTEGER              :: irA, irB, II
  REAL (KIND=_RL90)    :: nA, nB, zr, L, dq( Beam%Nsteps - 1 )
  REAL (KIND=_RL90)    :: znV( Beam%Nsteps ), rnV( Beam%Nsteps ), &
                           RcvrDeclAngleV ( Beam%Nsteps )
  COMPLEX (KIND=_RL90) :: dtau( Beam%Nsteps-1 )
  LOGICAL :: skip_step
!EOP

!$TAF init iRayCen0  = 'influence_iraycen'
!$TAF init iRayCen1  = static, (Beam%Nsteps-1)*nRz_per_range
!$TAF init iRayCen2  = static, (Beam%Nsteps-1)*ihop_nRR*nRz_per_range
!$TAF store ray2d = iRayCen0

  q0   = ray2D( 1 )%c / Angles%Dalpha   ! Reference for J = q0 / q

  dq   = ray2D( 2:Beam%Nsteps )%q( 1 ) - ray2D( 1:Beam%Nsteps-1 )%q( 1 )
  dtau = ray2D( 2:Beam%Nsteps )%tau    - ray2D( 1:Beam%Nsteps-1 )%tau

  ! Set the ray-centered coordinates (znV, rnV)
  ! pre-calculate ray normal based on tangent with c(s) scaling
  znV = -ray2D( 1:Beam%Nsteps )%t( 1 ) * ray2D( 1:Beam%Nsteps )%c
  rnV =  ray2D( 1:Beam%Nsteps )%t( 2 ) * ray2D( 1:Beam%Nsteps )%c

  DO ii=1,Beam%Nsteps
    IF ( ALL(ray2D(ii)%t==0.0) ) THEN
      RcvrDeclAngleV(ii) = 0.0
    ELSE
      RcvrDeclAngleV(ii) = rad2deg * ATAN2( ray2D(ii)%t(2), ray2D(ii)%t(1) )
    ENDIF
  ENDDO

  ! During reflection imag(q) is constant and adjacent normals cannot bracket
  ! a segment of the TL line, so no special treatment is necessary

  ! point source (cylindrical coordinates): default behavior
  Ratio1 = 1.0d0          !RG
  IF ( Beam%RunType( 4:4 ).EQ.'R' ) &
    Ratio1 = SQRT( ABS( COS( SrcDeclAngle / rad2deg ) ) )

  ray2D( 1:Beam%Nsteps )%Amp = Ratio1 * SQRT( ray2D( 1:Beam%Nsteps )%c ) &
                  * ray2D( 1:Beam%Nsteps )%Amp   ! pre-apply some scaling

  RcvrDepths: DO iz = 1, nRz_per_range
    zR = Pos%RZ( iz )

    phase = 0.0
    qOld  = ray2D( 1 )%q( 1 ) ! used to track KMAH index

    ! If normal is parallel to horizontal receiver line
    IF ( ABS( znV( 1 ) ).LT.1D-6 ) THEN
      nA  = 1D10
      rA  = 1D10
      irA = 1
    ELSE
      nA  = ( zR - ray2D( 1 )%x( 2 ) ) / znV( 1 )
      rA  = ray2D( 1 )%x( 1 ) + nA*rnV( 1 )
      !!! Find index of receiver: assumes uniform spacing in Pos%RR
      irA = MAX( MIN( INT( ( rA - Pos%RR( 1 ) ) / Pos%Delta_r )+1,  &
                     Pos%nRR ),                                    &
               1 )
    ENDIF

    Stepping: DO iS = 2, Beam%Nsteps
!$TAF store ira,irb,na,nb,phase,qold,ra = iRayCen1
      skip_step = .FALSE.

      ! Compute ray-centered coordinates, (znV, rnV)

      ! If normal is parallel to TL-line, skip to the next step on ray
      IF ( ABS(znV(iS)).LT.1D-10 ) THEN
        skip_step = .TRUE.
        rB = 1D10

      ELSE
        nB = (zR - ray2D(iS)%x(2)) / znV(iS)
        rB = ray2D(iS)%x(1) + nB * rnV(iS)

        ! Find index of receiver: assumes uniform spacing in Pos%RR
        irB = MAX(MIN(INT((rB - Pos%RR(1)) / Pos%Delta_r) + 1, &
                     Pos%nRR), &
                  1)

        ! Detect and skip duplicate points (happens at boundary reflection)
        IF ( ABS(ray2D(iS)%x(1) - ray2D(iS - 1)%x(1)).LT. &
             1.0D3 * SPACING(ray2D(iS)%x(1)) &
             .OR. irA.EQ.irB ) THEN
          rA = rB
          nA = nB
          irA = irB
          skip_step = .TRUE.
        ENDIF

      ENDIF ! IF ( ABS(znV(iS)).LT.1D-10 )

      IF ( .NOT.skip_step ) THEN
        !!! this should be pre-computed
        q = ray2D(iS - 1)%q(1)
        ! if phase shifts at caustics
        IF ((q.LE.0.0D0 .AND. qOld.GT.0.0D0) .OR. &
            (q.GE.0.0D0 .AND. qOld.LT.0.0D0)) &
          phase = phase + PI / 2.0D0

        qOld = q

        RcvrDeclAngle = RcvrDeclAngleV(iS)

        ! *** Compute contributions to bracketted receivers ***
        II = 0
        IF ( irB.LE.irA ) II = 1   ! going backwards in range

        ! Compute influence for each receiver
        DO ir = irA + 1 - II, irB + II, SIGN(1, irB - irA)
!!$TAF store Arr(:,ir,iz),nArrival(ir,iz) = iRayCen2
          W = (Pos%RR(ir) - rA) / (rB - rA)  ! relative range between rR
          n = ABS(nA + W * (nB - nA))
          q = ray2D(iS - 1)%q(1) + W * dq(iS - 1)  ! interpolated amplitude
          L = ABS(q) / q0   ! beam radius

          IF ( n.LT.L ) THEN  ! in beam window: update delay, Amp, phase
!$TAF store w = iRayCen2
            delay = ray2D(iS - 1)%tau + W * dtau(iS - 1)
            Amp = ray2D(iS)%Amp / SQRT(ABS(q))
            W = (L - n) / L  ! hat function: 1 on center, 0 on edge
            Amp = Amp * W
            phaseInt = ray2D(iS - 1)%Phase + phase
            !!! this should be precomputed
            IF ((q.LE.0.0D0 .AND. qOld.GT.0.0D0) .OR. &
                (q.GE.0.0D0 .AND. qOld.LT.0.0D0)) &
              phaseInt = phase + PI / 2.0D0  ! phase shifts at caustics

            CALL ApplyContribution(U(iz, ir))

          ENDIF ! IF ( n.LT.L ) THEN

        ENDDO ! DO ir

      ENDIF ! IF ( .NOT.skip_step )

      rA = rB
      nA = nB
      irA = irB

    ENDDO Stepping

  ENDDO RcvrDepths

  RETURN
  END !SUBROUTINE InfluenceGeoHatRayCen

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: InfluenceGeoHatCart
! !INTERFACE:
  SUBROUTINE InfluenceGeoHatCart( U, myThid )
! !DESCRIPTION:
!   Geometrically-spreading beams with a hat-shaped beam in Cartesian
!   coordinates. When Beam%Type(1:1)='G' or '^'

! !USES:
!   USE arr_mod, only: nArrival, Arr !RG
  USE angle_mod, only: Angles

! !INPUT PARAMETERS:
! U :: complex pressure field
! myThid :: my thread ID
  COMPLEX, INTENT( INOUT ) :: U( nRz_per_range, Pos%nRR )
  INTEGER, INTENT( IN )    :: myThid
! !OUTPUT PARAMETERS: U

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! irT, irTT :: indices of receivers
! x_ray, rayt, rayn :: vectors for ray coordinates
! x_rcvr :: coordinates of receivers
! rLen, RadiusMax, zMin, zMax, dqds :: beam variables
! dtauds :: delay along ray
! inRcvrRanges :: logical to skip steps with no influence
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER              :: irT(1), irTT ! irT needs size of 1, see MINLOC
  REAL (KIND=_RL90)    :: x_ray( 2 ), rayt( 2 ), rayn( 2 ), &
                        x_rcvr( 2, nRz_per_range ), rLen, RadiusMax, &
                        zMin, zMax, dqds
  COMPLEX (KIND=_RL90) :: dtauds
  LOGICAL              :: inRcvrRanges

!$TAF init iiitape1 = static, (Beam%Nsteps-1)
!$TAF init iiitape2 = static, (Beam%Nsteps-1)*ihop_nRR
!$TAF init iiitape3 = static, (Beam%Nsteps-1)*ihop_nRR*nRz_per_range

  q0           = ray2D( 1 )%c / Angles%Dalpha   ! Reference for J = q0 / q
  phase        = 0.0
  qOld         = ray2D( 1 )%q( 1 )       ! old KMAH index
  rA           = ray2D( 1 )%x( 1 )       ! range at start of ray, typically 0

  ! find index of first receiver to the right of rA
  irT = MINLOC( Pos%RR( 1 : Pos%nRR ), MASK = Pos%RR( 1 : Pos%nRR ).GT.rA )
  ir  = irT( 1 )
  ! if ray is left-traveling, get the first receiver to the left of rA
  IF ( ray2D( 1 )%t( 1 ).LT.0.0d0 .AND. ir.GT.1 ) ir = ir - 1

  ! point source: the default option
  Ratio1 = 1.0d0          !RG
  IF ( Beam%RunType( 4:4 ).EQ.'R' ) &
    Ratio1 = SQRT( ABS( COS( SrcDeclAngle / rad2deg ) ) )

  Stepping: DO iS = 2, Beam%Nsteps
!$TAF store phase,qold,ra = iiitape1
    rB     = ray2D( iS   )%x( 1 )
    x_ray  = ray2D( iS-1 )%x

    ! compute normalized tangent (we need to measure the step length)
    rayt = ray2D( iS )%x - x_ray
    IF ( ALL(rayt.EQ.0.0) ) THEN
      rlen = 0.0
    ELSE
      rlen = NORM2( rayt )
    ENDIF

    ! if duplicate point in ray, skip to next step along the ray
    IF ( rlen.GE.1.0D3*SPACING( ray2D( iS )%x( 1 ) ) ) THEN
!$TAF store rlen,rayt= iiitape1
      rayt = rayt / rlen                    ! unit tangent of ray @ A
      rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal  of ray @ A
      IF ( ALL(rayt.EQ.0.0) ) THEN
        RcvrDeclAngle = 0.0
      ELSE
        RcvrDeclAngle = rad2deg * ATAN2( rayt( 2 ), rayt( 1 ) )
      ENDIF

      q      = ray2D( iS-1 )%q( 1 )
      dqds   = ray2D( iS   )%q( 1 ) - q
      dtauds = ray2D( iS   )%tau    - ray2D( iS-1 )%tau

      !IESCO22: q only changes signs on direct paths, no top/bot bounces
      IF( q.LE.0. .AND. qOld.GT.0. .OR. q.GE.0. .AND. qOld.LT.0. ) &
        phase = phase + PI / 2.   ! phase shifts at caustics
      qOld = q

      ! Radius calc from beam radius projected onto vertical line
      RadiusMax = MAX( ABS( q ), ABS( ray2D( iS )%q( 1 ) ) ) &
               / q0 / ABS( rayt( 1 ) ) ! IESCO24: AKA rayn( 2 )

      ! depth limits of beam; IESCO22: a large range of about 1/2 box depth
      IF ( ABS( rayt( 1 ) ).GT.0.5 ) THEN   ! shallow angle ray
        zmin = min( x_ray( 2 ), ray2D( iS )%x( 2 ) ) - RadiusMax
        zmax = max( x_ray( 2 ), ray2D( iS )%x( 2 ) ) + RadiusMax
      ELSE                                  ! steep angle ray
        zmin = -HUGE( zmin )
        zmax = +HUGE( zmax )
      ENDIF

      ! compute beam influence for this segment of the ray
      inRcvrRanges=.TRUE.
      RcvrRanges: DO ir = 1, ihop_nRR
!$TAF store inrcvrranges = iiitape2
        ! is Rr( ir ) contained in [ rA, rB )? Then compute beam influence
        IF ( Pos%RR( ir ).GE.MIN( rA, rB ) &
            .AND. Pos%RR( ir ).LT.MAX( rA, rB ) &
            .AND. inRcvrRanges ) THEN
          x_rcvr( 1, 1:nRz_per_range ) = Pos%RR( ir )

          IF ( Beam%RunType( 5 : 5 ).EQ.'I' ) THEN ! irregular grid
            x_rcvr( 2, 1 ) = Pos%RZ( ir )
          ELSE ! default: rectilinear grid
            x_rcvr( 2, 1:nRz_per_range ) = Pos%RZ( 1:nRz_per_range )
          ENDIF

          RcvrDepths: DO iz = 1, nRz_per_range
!!$TAF store Arr(:,ir,iz),nArrival(ir,iz),q = iiitape3
            ! is x_rcvr( 2, iz ) contained in ( zmin, zmax )?
            IF (      x_rcvr( 2, iz ).GE.zmin &
                .AND. x_rcvr( 2, iz ).LE.zmax ) THEN
              ! normalized proportional distance along ray
              s = DOT_PRODUCT( x_rcvr( :, iz ) - x_ray, rayt ) / rlen
              ! normal distance to ray
              n = ABS( DOT_PRODUCT( x_rcvr( :, iz ) - x_ray, rayn ) )
              ! interpolated amplitude in [meters]
              q = q + s*dqds
              ! beam radius; IESCO22 smaller then previous RadiusMax
              RadiusMax = ABS( q / q0 )

              IF ( n.LT.RadiusMax ) THEN
#ifdef IHOP_WRITE_OUT
                WRITE(msgBuf,'(A,F10.2)') &
                  "Influence: Eigenray w RadiusMax = ", RadiusMax
                IF ( IHOP_dumpfreq.GE.0 ) &
                  CALL PRINT_MESSAGE( msgbuf, PRTFile, &
                                    SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
                ! interpolated delay
                delay    = ray2D( iS-1 )%tau + s*dtauds
                Amp      = Ratio1 * SQRT( ray2D( iS )%c / ABS( q ) ) &
                            * ray2D( iS )%Amp
                ! hat function: 1 on center, 0 on edge
                W        = ( RadiusMax - n ) / RadiusMax
                Amp      = Amp*W

                IF (      q.LE.0.0d0 .AND. qOld.GT.0.0d0 &
                     .OR. q.GE.0.0d0 .AND. qOld.LT.0.0d0 ) THEN
                  phaseInt = phase + PI / 2.   ! phase shifts at caustics
                  ! IESCO22: shouldn't this be = phaseInt + PI/2
                ELSE
                  phaseInt = ray2D( iS-1 )%Phase + phase
                ENDIF

                CALL ApplyContribution( U( iz, ir ) )

              ENDIF ! IF ( n.LT.RadiusMax )

            ENDIF ! is x_rcvr( 2, iz ) contained in ( zmin, zmax )?

          ENDDO RcvrDepths

        ENDIF ! IF ( Pos%RR( ir ).GE.MIN( rA, rB ) &

        ! bump receiver index, ir, towards rB
        IF ( Pos%RR( ir ).LT.rB ) THEN
          IF ( ir.GE.Pos%nRR        ) inRcvrRanges=.FALSE. ! to next step on ray
          irTT = ir + 1                     ! bump right
          IF ( Pos%RR( irTT ).GE.rB ) inRcvrRanges=.FALSE.
        !ELSE
        !   IF ( ir.LE.1              ) inRcvrRanges=.FALSE. ! to next step on ray
        !   irTT = ir - 1                     ! bump left
        !   IF ( Pos%RR( irTT ).LE.rB ) inRcvrRanges=.FALSE.
        !
        ENDIF

        !ir = irTT
      ENDDO RcvrRanges

      rA = rB

    ENDIF ! if duplicate point in ray, skip to next step along the ray

  ENDDO Stepping

  RETURN
  END !SUBROUTINE InfluenceGeoHatCart

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: InfluenceGeoGaussianCart
! !INTERFACE:
  SUBROUTINE InfluenceGeoGaussianCart( U, myThid )
! !DESCRIPTION:
!   Geometrically-spreading beams with a Gaussian beam in Cartesian
!   coordinates. When Beam%Type(1:1)='B'
!   beam window: kills beams outside e**(-0.5 * ibwin**2 )

! !USES:
!    USE arr_mod, only: nArrival, Arr         !RG
  USE angle_mod, only: Angles

! !INPUT PARAMETERS:
! U :: complex pressure field
! myThid :: my thread ID
  COMPLEX, INTENT( INOUT ) :: U( nRz_per_range, Pos%nRR )
  INTEGER, INTENT( IN )    :: myThid
! !OUTPUT PARAMETERS: U

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! BeamWindow :: beam window size in range
! irT, irTT :: indices of receivers
! x_ray, rayt, rayn :: vectors for ray coordinates
! x_rcvr :: coordinates of receivers
! rLen, RadiusMax, zMin, zMax, sigma, lambda, A, dqds :: beam variables
! dtauds :: delay along ray
! inRcvrRanges :: logical to skip steps with no influence
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
  INTEGER, PARAMETER   :: BeamWindow = 4
  INTEGER              :: irT( 1 ), irTT
  REAL (KIND=_RL90)    :: x_ray( 2 ), rayt( 2 ), rayn( 2 ), x_rcvr( 2 ), &
                        rLen, RadiusMax, zMin, zMax, sigma, lambda, A, dqds
  COMPLEX (KIND=_RL90) :: dtauds
  LOGICAL              :: inRcvrRanges
!EOP

!$TAF init iGauCart1 = static, (Beam%Nsteps-1)
!$TAF init iGauCart2 = static, (Beam%Nsteps-1)*ihop_nRR
!$TAF init iGauCart3 = static, (Beam%Nsteps-1)*ihop_nRR*nRz_per_range
  q0     = ray2D( 1 )%c / Angles%Dalpha   ! Reference for J = q0 / q
  phase  = 0
  qOld   = ray2D( 1 )%q( 1 )       ! used to track KMAH index
  rA     = ray2D( 1 )%x( 1 )       ! range at start of ray

  ! what if never satistified?
  ! what if there is a single receiver (ir = 0 possible)

  ! irT: find index of first receiver to the right of rA
  irT = MINLOC( Pos%RR( 1:Pos%nRR ), MASK=Pos%RR( 1:Pos%nRR ).GT.rA )
  ir  = irT( 1 )

  ! if ray is left-traveling, get the first receiver to the left of rA
  IF ( ray2D( 1 )%t( 1 ).LT.0.0d0 .AND. ir.GT.1 ) ir = ir - 1

  ! sqrt( 2 * PI ) represents a sum of Gaussians in free space
  IF ( Beam%RunType( 4:4 ).EQ.'R' ) THEN   ! point source
    Ratio1 = SQRT( ABS( COS( SrcDeclAngle / rad2deg ) ) ) / SQRT( 2.*PI )
  ELSE    ! line  source
    Ratio1 = 1 / SQRT( 2.*PI )
  ENDIF

  Stepping: DO iS = 2, Beam%Nsteps
!$TAF store phase,qold,ra = iGauCart1
    rB    = ray2D( iS   )%x( 1 )
    x_ray = ray2D( iS-1 )%x

    ! compute normalized tangent (compute it because we need to measure the
    ! step length)
    rayt = ray2D( iS )%x - ray2D( iS-1 )%x
    IF ( ALL(rayt.EQ.0.) ) THEN
      rlen = 0.0
    ELSE
      rlen = NORM2( rayt )
    ENDIF

    ! if duplicate point in ray, skip to next step along the ray
    IF ( rlen.GE.1.0D3 * SPACING( ray2D( iS )%x( 1 ) ) ) THEN

!$TAF store rlen,rayt = iGauCart1
    rayt = rayt / rlen
    rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal to ray
    IF ( ALL(rayt.EQ.0.0) ) THEN
      RcvrDeclAngle = 0.0
    ELSE
      RcvrDeclAngle = rad2deg * ATAN2( rayt( 2 ), rayt( 1 ) )
    ENDIF

    q      = ray2D( iS-1 )%q( 1 )
    dqds   = ray2D( iS )%q( 1 ) - q
    dtauds = ray2D( iS )%tau    - ray2D( iS-1 )%tau

    !IESCO22: q only changes signs on direct paths, no top/bot bounces
    IF ( q.LE.0.0 .AND. qOld.GT.0.0 .OR. q.GE.0.0 .AND. qOld.LT.0.0 ) &
      phase = phase + PI / 2.   ! phase shifts at caustics

    qOld = q

    ! calculate beam width beam radius projected onto vertical line
    lambda    = ray2D( iS-1 )%c / IHOP_freq
    sigma     = MAX( ABS( q ), ABS( ray2D( iS )%q( 1 ) ) ) &
               / q0 / ABS( rayt( 1 ) ) ! IESCO24: AKA rayn( 2 )
    sigma     = MAX( sigma, &
                     MIN( 0.2*IHOP_freq*REAL( ray2D( iS )%tau ), &
                        PI*lambda ) )
    ! Note on min: "Weinberg and Keenan suggest limiting a beam to a
    !               point, by imposing a minimum beam width of pilambda."
    !               - Jensen, Comp OA 2011
    ! default is 2 standard deviations of coverage of the Gaussian curve
    RadiusMax = BeamWindow*sigma

    ! depth limits of beam; IESCO22: a large range of about 1/2 box depth
    IF ( ABS( rayt( 1 ) ).GT.0.5 ) THEN   ! shallow angle ray
      zmin   = min( ray2D( iS-1 )%x( 2 ), ray2D( iS )%x( 2 ) ) - RadiusMax
      zmax   = max( ray2D( iS-1 )%x( 2 ), ray2D( iS )%x( 2 ) ) + RadiusMax
    ELSE                                 ! steep angle ray
      zmin = -HUGE( zmin )
      zmax = +HUGE( zmax )
    ENDIF

    ! compute beam influence for this segment of the ray
    inRcvrRanges=.TRUE.
    RcvrRanges: DO ir = 1, ihop_nRR
!$TAF store inrcvrranges = iGauCart2
      ! is Rr( ir ) contained in [ rA, rB )? Then compute beam influence
      IF ( Pos%RR( ir ).GE.MIN( rA, rB ) .AND. &
           Pos%RR( ir ).LT.MAX( rA, rB ) .AND. &
           inRcvrRanges ) THEN

        RcvrDepths: DO iz = 1, nRz_per_range
!!$TAF store arr,nArrival,q = iGauCart3
          IF ( Beam%RunType( 5 : 5 ).EQ.'I' ) THEN
            x_rcvr = [ Pos%RR( ir ), Pos%RZ( ir ) ]   ! irregular   grid
          ELSE
            x_rcvr = [ Pos%RR( ir ), Pos%RZ( iz ) ]   ! rectilinear grid
          ENDIF

          ! is x_rcvr( 2, iz ) contained in ( zmin, zmax )?
          IF ( x_rcvr( 2 ).GE.zmin .AND. x_rcvr( 2 ).LE.zmax ) THEN
            ! proportional distance along ray
            s = DOT_PRODUCT( x_rcvr - x_ray, rayt ) / rlen
            ! normal distance to ray
            n = ABS( DOT_PRODUCT( x_rcvr - x_ray, rayn ) )
            ! interpolated amplitude in [meters]
            q = q + s*dqds
            ! beam radius; IESCO22 smaller then previous RadiusMax
            sigma = ABS( q / q0 )
            sigma = MAX( sigma, &
                        MIN( 0.2*IHOP_freq*REAL( ray2D( iS )%tau ), &
                              PI*lambda ) )

            IF ( n.LT.BeamWindow*sigma ) THEN   ! Within beam window?
#ifdef IHOP_WRITE_OUT
              WRITE(msgBuf,'(A,F10.2)') &
                "Influence: Eigenray w RadiusMax = ", RadiusMax
               IF ( IHOP_dumpfreq .GE. 0) &
                 CALL PRINT_MESSAGE( msgbuf, PRTFile, &
                   SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
              ! interpolated delay
              A        = ABS( q0 / q )
              delay    = ray2D( iS-1 )%tau + s*dtauds
              Amp      = Ratio1 * SQRT( ray2D( iS )%c / ABS( q ) ) &
                          * ray2D( iS )%Amp
              ! W : Gaussian decay
              W        = EXP( -0.5*( n / sigma )**2 ) / ( sigma*A )
              Amp      = Amp*W
              phaseInt = ray2D( iS )%Phase + phase
              IF ( q.LE.0.0d0 .AND. qOld.GT.0.0d0 .OR. &
                   q.GE.0.0d0 .AND. qOld.LT.0.0d0 ) &
                phaseInt = phase + PI / 2.  ! phase shifts at caustics

              CALL ApplyContribution( U( iz, ir ) )

            ENDIF ! IF ( n.LT.BeamWindow*sigma )

          ENDIF ! IF ( x_rcvr( 2, iz ).GE.zmin .AND. x_rcvr( 2, iz ).LE.zmax )
        
        ENDDO RcvrDepths

      ENDIF ! IF ( Pos%RR( ir ).GE.MIN( rA, rB ) .AND. Pos%RR( ir ).LT.MAX( rA, rB ) .AND. inRcvrRanges )

      ! bump receiver index, ir, towards rB
      IF ( Pos%RR( ir ).LT.rB ) THEN
        IF ( ir.GE.Pos%nRR        ) inRcvrRanges=.FALSE. ! to next step on ray
        irTT = ir + 1                     ! bump right
        IF ( Pos%RR( irTT ).GE.rB ) inRcvrRanges=.FALSE.
      ELSE
        IF ( ir.LE.1              ) inRcvrRanges=.FALSE. ! to next step on ray
        irTT = ir - 1                     ! bump left
        IF ( Pos%RR( irTT ).LE.rB ) inRcvrRanges=.FALSE.
      ENDIF
      !ir = irTT
    ENDDO RcvrRanges

    rA = rB
    ENDIF ! if duplicate point in ray, skip to next step along the ray
  ENDDO Stepping

  RETURN
  END !SUBROUTINE InfluenceGeoGaussianCart

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ApplyContribution
! !INTERFACE:
  SUBROUTINE ApplyContribution( U )
! !DESCRIPTION:
!   Apply the contribution of the current ray to the pressure field U.

! !USES:
  USE writeray, only: WriteRayOutput
  USE arr_mod,  only: AddArr
  USE ihop_mod, only: afreq, RAYFile, DELFile, maxN

! !INPUT PARAMETERS:
! U :: complex pressure field
  COMPLEX, INTENT( INOUT ) :: U
! !OUTPUT PARAMETERS: U

! !LOCAL VARIABLES:
! tmpDelay :: temporary array for delay values
  INTEGER :: i
  REAL(KIND=_RL90) :: tmpX(maxN), tmpY(maxN), tmpDelay(maxN)
!EOP

! initiate temporary arrays
  tmpX = 0.
  tmpY = 0.
  tmpDelay = 0.
  DO i = 1,iS
    tmpX(i) = ray2D(i)%x(1)
    tmpY(i) = ray2D(i)%x(2)
    tmpDelay(i) = REAL(ray2D(i)%tau)
  ENDDO

! write rays and delays
  SELECT CASE( Beam%RunType( 1:1 ) )
  CASE ( 'E' )                ! eigenrays
    U=U
!    tmpDelay = 0.
    CALL WriteRayOutput( RAYFile, iS,         &
      tmpX(1:iS), tmpY(1:iS),                 &
      ray2D(iS)%NumTopBnc, ray2D(iS)%NumBotBnc )

  CASE ( 'e' )                ! eigenrays AND arrivals
    U=U
!    tmpDelay = 0.
    CALL WriteRayOutput( RAYFile, iS,         &
      tmpX(1:iS), tmpY(1:iS),                 &
      ray2D(iS)%NumTopBnc, ray2D(iS)%NumBotBnc )
    IF (writeDelay) THEN
      CALL WriteRayOutput( DELFile, iS,       &
        tmpDelay(1:iS), tmpY(1:iS),           &
        ray2D(iS)%NumTopBnc, ray2D(iS)%NumBotBnc )
    ENDIF

    CALL AddArr( afreq, iz, ir, Amp, phaseInt, delay,  &
      RcvrDeclAngle, ray2D( iS )%NumTopBnc, ray2D( iS )%NumBotBnc )

  CASE ( 'A', 'a' )           ! arrivals
    U=U
!    tmpDelay = 0.
    CALL AddArr( afreq, iz, ir, Amp, phaseInt, delay,  &
      RcvrDeclAngle, ray2D( iS )%NumTopBnc, ray2D( iS )%NumBotBnc )

  CASE ( 'C' )                ! coherent TL
!    tmpDelay = 0.
    U = U + CMPLX( Amp * EXP( -oneCMPLX * ( afreq * delay - phaseInt ) ) )

  CASE ( 'S', 'I' )                ! incoherent/semicoherent TL
!    tmpDelay = 0.
    IF ( Beam%Type( 1:1 ).EQ.'B' ) THEN   ! Gaussian beam
      U = U + SNGL( SQRT( 2. * PI ) &
            * ( Amp * EXP( AIMAG( afreq * delay ) ) )**2 )
    ELSE
      U = U + SNGL( &
               ( Amp * EXP( AIMAG( afreq * delay ) ) )**2 )
    ENDIF

  CASE DEFAULT                ! incoherent/semicoherent TL
!    tmpDelay = 0.
    IF ( Beam%Type( 1:1 ).EQ.'B' ) THEN   ! Gaussian beam
      U = U + SNGL( SQRT( 2. * PI ) &
            * ( Amp * EXP( AIMAG( afreq * delay ) ) )**2 )
    ELSE
      U = U + SNGL( &
               ( Amp * EXP( AIMAG( afreq * delay ) ) )**2 )
    ENDIF

  END SELECT ! SELECT CASE( Beam%RunType( 1:1 ) )

  RETURN
  END !SUBROUTINE ApplyContribution

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ScalePressure
! !INTERFACE:
  SUBROUTINE ScalePressure( c, r, U, nRz, nR, RunType, freq )
! !DESCRIPTION:
!   Scale the pressure field U according to the run type and range.

! !USES:
  USE angle_mod, only: Angles

! !INPUT PARAMETERS:
! c :: nominal sound speed, m/s
! r :: Rr ranges, m
! U :: complex pressure field
! nRz :: number of depths in the pressure field
! nR :: number of ranges in the pressure field
! RunType :: run type, e.g. 'C' for Cerveny Gaussian beams in Cartesian coordinates
! freq :: source frequency, Hz
  REAL (KIND=_RL90), INTENT( IN    ) :: c, r( nR )
  COMPLEX,           INTENT( INOUT ) :: U( nRz, nR )
  INTEGER,           INTENT( IN    ) :: nRz, nR
  CHARACTER*(5),     INTENT( IN    ) :: RunType
  REAL (KIND=_RL90), INTENT( IN    ) :: freq
! !OUTPUT PARAMETERS: U

! !LOCAL VARIABLES:
! const :: scale factor for field
! factor :: scale factor for cylindrical spreading
  REAL (KIND=_RL90) :: const, factor
!EOP

  ! Compute scale factor for field
  SELECT CASE ( RunType( 2:2 ) )
  CASE ( 'C' )   ! Cerveny Gaussian beams in Cartesian coordinates
    const = -Angles%Dalpha * SQRT( freq ) / c
  CASE ( 'R' )   ! Cerveny Gaussian beams in Ray-centered coordinates
    const = -Angles%Dalpha * SQRT( freq ) / c
  CASE DEFAULT
    const = -1.0
  END SELECT

  ! If incoherent run, convert intensity to pressure
  IF ( RunType( 1:1 ).NE.'C' ) U = SQRT( REAL( U ) )

  ! scale and/or incorporate cylindrical spreading
  Ranges: DO ir = 1, nR
    IF ( RunType( 4:4 ).EQ.'X' ) THEN   ! line source
      factor = -4.0 * SQRT( PI ) * const
    ELSE                                  ! point source
      IF ( r( ir ).EQ.0 ) THEN
        factor = 0.0D0 ! avoid /0 at origin, return pressure = 0
      ELSE
        factor = const / SQRT( ABS( r( ir ) ) )
      ENDIF
    ENDIF
    U( :, ir ) = SNGL( factor ) * U( :, ir )

  ENDDO Ranges

  RETURN
  END !SUBROUTINE ScalePressure

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: Hermite
! !INTERFACE:
  REAL (KIND=_RL90) FUNCTION Hermite( x, x1, x2 )
! !DESCRIPTION:
!   Hermite cubic taper function for smoothing influence functions.
! Calculates a smoothing function based on the h0 hermite cubic
! x is the point where the function is to be evaluated
! returns:
! [  0, x1  ] = 1
! [ x1, x2  ] = cubic taper from 1 to 0
! [ x2, inf ] = 0

! !USES: None

! !INPUT PARAMETERS:
! x :: point where the function is to be evaluated
! x1 :: first point of the taper
! x2 :: second point of the taper
  REAL (KIND=_RL90), INTENT( IN ) :: x, x1, x2
! !OUTPUT PARAMETERS: Hermite

! !LOCAL VARIABLES:
! Ax :: absolute value of x
! u :: normalized distance from x1 to x2
  REAL (KIND=_RL90) :: Ax, u
!EOP


  Ax  = ABS( x  )

  IF ( Ax.LE.x1 ) THEN
    Hermite = 1.0d0
  ELSEIF ( Ax.GE.x2 ) THEN
    Hermite = 0.0d0
  ELSE
    u       = ( Ax - x1 ) / ( x2 - x1 )
    Hermite = ( 1.0d0 + 2.0d0 * u ) * ( 1.0d0 - u ) ** 2
  ENDIF

  RETURN
  END !FUNCTION Hermite

END !MODULE influence
