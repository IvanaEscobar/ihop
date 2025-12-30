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
  USE ihop_mod,  only: rad2deg, oneCMPLX, PRTFile, nMax, Beam, ray2D, &
                       SrcDeclAngle, nRz_per_range
  USE srPos_mod, only: Pos
  IMPLICIT NONE
!  == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef ALLOW_COST
# include "IHOP_COST.h"
#endif

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC calculateInfluence, ScalePressure
!=======================================================================

!! == Module variables ==
!  INTEGER,              PRIVATE :: iz, ir
!  REAL    (KIND=_RL90), PRIVATE :: Ratio1 = 1.0D0 ! scale factor for a line source
!  REAL    (KIND=_RL90), PRIVATE :: W, s, n, Amp, phase, phaseInt, &
!                                   q0, q, qOld, RcvrDeclAngle, rA, rB
!  COMPLEX (KIND=_RL90), PRIVATE :: delay
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R calculateInfluence
! S/R InfluenceGeoHatRayCen
! S/R InfluenceGeoHatCart
! S/R InfluenceGeoGaussianCart
! S/R ApplyContribution
! S/R ScalePressure
! FXN Hermite
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: calculateInfluence
! !INTERFACE:
  SUBROUTINE calculateInfluence( myThid )
! !DESCRIPTION:
!   determine eigenrays and apply pressure field, U
  INTEGER, INTENT ( IN ) :: myThid

  SELECT CASE ( Beam%Type( 1:1 ) )
  CASE ( 'g' )
    CALL InfluenceGeoHatRayCen( myThid )
  CASE ( 'B' )
    CALL InfluenceGeoGaussianCart( myThid )
  CASE ( 'G','^' )
    CALL InfluenceGeoHatCart( myThid )
  CASE DEFAULT !IEsco22: thesis is in default behavior
    CALL InfluenceGeoHatCart( myThid )
  END SELECT

  RETURN
  END !SUBROUTINE calculateInfluence

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: InfluenceGeoHatRayCen
! !INTERFACE:
  SUBROUTINE InfluenceGeoHatRayCen( myThid )
! !DESCRIPTION:
!   Geometrically-spreading beams with a hat-shaped beam in ray-centered
!   coordinates

! !USES:
!  USE arr_mod, only: nArrival, Arr         !RG
  USE angle_mod, only: Angles

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN )    :: myThid
! !OUTPUT PARAMETERS: 

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
  INTEGER :: iH, iZ, iR
  REAL    (KIND=_RL90) :: Ratio1 ! scale factor for a line source
  REAL    (KIND=_RL90) :: W, s, n, phase, rA, rB, &
                          q0, q, qOld
  REAL    (KIND=_RL90) :: Amp, phaseInt, RcvrDeclAngle
  COMPLEX (KIND=_RL90) :: delay
  INTEGER              :: irA, irB, II
  REAL (KIND=_RL90)    :: nA, nB, zr, L, dq( Beam%nSteps-1 )
  REAL (KIND=_RL90)    :: znV( Beam%nSteps ), rnV( Beam%nSteps ), &
                           RcvrDeclAngleV ( Beam%nSteps )
  COMPLEX (KIND=_RL90) :: dtau( Beam%nSteps-1 )
  LOGICAL :: skip_step
!EOP
  REAL (KIND=_RL90) :: qAtmp, qBtmp, rayntmp, qtmp

!!$TAF init iRayCen0  = 'influence_iraycen'
!!$TAF init iRayCen1  = static, (Beam%nSteps-1)*nRz_per_range
!!$TAF init iRayCen2  = static, (Beam%nSteps-1)*ihop_nRR*nRz_per_range
!!$TAF store ray2d = iRayCen0
! init local variables
  W = 0
  s = 0
  n = 0
  rB = rA
  q = 0
  Amp = 0
  phaseint = 0
  rcvrdeclangle = 180
  delay = 0

  q0   = ray2D( 1 )%c / Angles%Dalpha   ! Reference for J = q0 / q

  dq   = ray2D( 2:Beam%nSteps )%q( 1 ) - ray2D( 1:Beam%nSteps-1 )%q( 1 )
  dtau = ray2D( 2:Beam%nSteps )%tau    - ray2D( 1:Beam%nSteps-1 )%tau

  ! Set the ray-centered coordinates (znV, rnV)
  ! pre-calculate ray normal based on tangent with c(s) scaling
  znV = -ray2D( 1:Beam%nSteps )%t( 1 ) * ray2D( 1:Beam%nSteps )%c
  rnV =  ray2D( 1:Beam%nSteps )%t( 2 ) * ray2D( 1:Beam%nSteps )%c

  DO ii=1,Beam%nSteps
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

  ! Scale amplitude
  ray2D( 1:Beam%nSteps )%Amp = Ratio1 * SQRT( ray2D( 1:Beam%nSteps )%c ) &
                  * ray2D( 1:Beam%nSteps )%Amp

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

    Stepping: DO iH = 2, Beam%nSteps
!$TAF store ira,irb,na,nb,phase,qold,ra = iRayCen1
      skip_step = .FALSE.

      ! Compute ray-centered coordinates, (znV, rnV)

      ! If normal is parallel to TL-line, skip to the next step on ray
      IF ( ABS(znV(iH)).LT.1D-10 ) THEN
        skip_step = .TRUE.
        rB = 1D10

      ELSE
        nB = (zR - ray2D(iH)%x(2)) / znV(iH)
        rB = ray2D(iH)%x(1) + nB * rnV(iH)

        ! Find index of receiver: assumes uniform spacing in Pos%RR
        irB = MAX(MIN(INT((rB - Pos%RR(1)) / Pos%Delta_r) + 1, &
                     Pos%nRR), &
                  1)

        ! Detect and skip duplicate points (happens at boundary reflection)
        IF ( ABS(ray2D(iH)%x(1) - ray2D( iH-1 )%x(1)).LT. &
             1.0D3*SPACING(ray2D(iH)%x(1) ) &
             .OR. irA.EQ.irB ) THEN
          rA = rB
          nA = nB
          irA = irB
          skip_step = .TRUE.
        ENDIF

      ENDIF ! IF ( ABS(znV(iH)).LT.1D-10 )

      IF ( .NOT.skip_step ) THEN
        !!! this should be pre-computed
        q = ray2D( iH-1 )%q(1)
        ! if phase shifts at caustics
        IF ((q.LE.0.0D0 .AND. qOld.GT.0.0D0) .OR. &
            (q.GE.0.0D0 .AND. qOld.LT.0.0D0)) &
          phase = phase + PI / 2.0D0

        qOld = q

        RcvrDeclAngle = RcvrDeclAngleV(iH)

        ! *** Compute contributions to bracketted receivers ***
        II = 0
        IF ( irB.LE.irA ) II = 1   ! going backwards in range

        ! Compute influence for each receiver
        DO ir = irA+1-II, irB+II, SIGN(1, irB-irA)
!!$TAF store Arr(:,ir,iz),nArrival(ir,iz) = iRayCen2
          W = (Pos%RR(ir) - rA) / (rB - rA)  ! relative range between rR
          n = ABS(nA + W*(nB - nA))
          q = ray2D( iH-1 )%q(1) + W*dq( iH-1 )  ! interpolated amplitude
          L = ABS(q) / q0   ! beam radius

          IF ( n.LT.L ) THEN  ! in beam window: update delay, Amp, phase
!!$TAF store w = iRayCen2
            delay = ray2D(iH - 1)%tau + W * dtau(iH - 1)
            Amp = ray2D(iH)%Amp / SQRT(ABS(q))
            W = (L - n) / L  ! hat function: 1 on center, 0 on edge
            Amp = Amp * W
            phaseInt = ray2D(iH - 1)%Phase + phase
            !!! this should be precomputed
            IF ((q.LE.0.0D0 .AND. qOld.GT.0.0D0) .OR. &
                (q.GE.0.0D0 .AND. qOld.LT.0.0D0)) &
              phaseInt = phase + PI / 2.0D0  ! phase shifts at caustics

            CALL ApplyContribution( iH, iR, iZ, delay, &
                                    amp, phaseint, rcvrdeclangle )

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
  SUBROUTINE InfluenceGeoHatCart( myThid )
! !DESCRIPTION:
!   Geometrically-spreading beams with a hat-shaped beam in Cartesian
!   coordinates. When Beam%Type(1:1)='G' or '^'

! !USES:
!   USE arr_mod, only: nArrival, Arr !RG
  USE angle_mod, only: Angles

! !INPUT PARAMETERS:
! myThid :: my thread ID
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
  INTEGER :: iH, iZ, iR
  REAL    (KIND=_RL90) :: Ratio1 ! scale factor for a line source
  REAL    (KIND=_RL90) :: W, s, n, phase, rA, xB(2), &
                          q0, q, qOld
  REAL    (KIND=_RL90) :: Amp, phaseInt, RcvrDeclAngle
  COMPLEX (KIND=_RL90) :: delay

!  INTEGER              :: irT(1), ! irT needs size of 1, see MINLOC
  INTEGER              :: irTT
  REAL (KIND=_RL90)    :: xA( 2 ), rayt( 2 ), rayn( 2 ), &
                          x_rcvr( 2, nRz_per_range ), rLen, RadiusMax, &
                          zMin, zMax, dqds
  REAL (KIND=_RL90) :: qAtmp, qBtmp, rayntmp, qtmp
  COMPLEX (KIND=_RL90) :: dtauds
  LOGICAL              :: inRcvrRanges

!!$TAF init iiitape1 = static, (Beam%nSteps-1)
!!$TAF init iiitape2 = static, (Beam%nSteps-1)*ihop_nRR
!!$TAF init iiitape3 = static, (Beam%nSteps-1)*ihop_nRR*nRz_per_range

  ! init local variables
  q0           = ray2D( 1 )%c / Angles%Dalpha   ! Reference for J = q0 / q
  phase        = 0.0
  qOld         = ray2D( 1 )%q( 1 )       ! old KMAH index
  rA           = ray2D( 1 )%x( 1 )       ! range at start of ray, typically 0
  !IESCO25: TAF TLM inits
  W = 0
  s = 0
  n = 0
  xB = [rA, zeroRL]
  q = 0
  Amp = 0
  phaseint = 0
  rcvrdeclangle = 180
  delay = 0
  xA = [zeroRL, zeroRL]
  rayt = [zeroRL, zeroRL]
  rayn = [zeroRL, zeroRL]
  rLen = 0
  RadiusMax = 0
  zMin = 0
  zMax = 0
  dqds = 0
  qatmp = 0
  qbtmp = 0
  rayntmp = 0
  qtmp = 0
  dtauds = (0,0)


!  ! find index of first receiver to the right of rA
!  irT = MINLOC( Pos%RR( 1:Pos%nRR ), MASK = Pos%RR( 1:Pos%nRR ).GT.rA )
!  ir  = irT( 1 )
!  ! if ray is left-traveling, get the first receiver to the left of rA
!  IF ( ray2D( 1 )%t( 1 ).LT.0.0d0 .AND. ir.GT.1 ) ir = ir - 1

  ! point source: the default option
  ! IESCO25: in application, +-90deg is never used. May need to hard code 
  ! defenses for TAF diff
  Ratio1 = 1.0 _d 0          !RG
  IF ( Beam%RunType( 4:4 ).EQ.'R' ) &
    Ratio1 = SQRT( ABS( COS( SrcDeclAngle / rad2deg ) ) )

  !IESCO25 for TAF TLM
  DO iH=1,35604
    geninfluence(iH) = 0.0
  ENDDO

  Stepping: DO iH = 2, Beam%nSteps
!!$TAF store phase,qold,ra = iiitape1
    xA  = ray2D( iH-1 )%x
    xB  = ray2D( iH   )%x
!    rB  = ray2D( iH   )%x( 1 )

    ! compute normalized tangent (we need to measure the step length)
    rayt = xB - xA !IESCO25: dxds
    IF ( ALL(rayt.EQ.0.0) ) THEN
      rlen = 0.0
    ELSE
      rlen = NORM2( rayt )
    ENDIF

    ! if duplicate point in ray, skip to next step along the ray
    IF ( rlen.GE.1.0D3*SPACING( xB(1) ) ) THEN
!!$TAF store rlen,rayt= iiitape1
      rayt = rayt / rlen                    ! unit tangent of ray @ A
      rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal  of ray @ A
      
      IF ( ALL(rayt.EQ.0.0) ) THEN
        RcvrDeclAngle = 0.0
      ELSE
        RcvrDeclAngle = rad2deg * ATAN2( rayt( 2 ), rayt( 1 ) )
      ENDIF

      q      = ray2D( iH-1 )%q( 1 )
      dqds   = ray2D( iH   )%q( 1 ) - q
      dtauds = ray2D( iH   )%tau    - ray2D( iH-1 )%tau

      !IESCO22: q only changes signs on direct paths, no top/bot bounces
      IF( q.LE.0. .AND. qOld.GT.0. .OR. q.GE.0. .AND. qOld.LT.0. ) &
        phase = phase + PI / 2.   ! phase shifts at caustics
      qOld   = ray2D( iH-1 )%q( 1 )

      ! Radius calc from beam radius projected onto vertical line
      ! IESCO25: make this TAF TLM friendly, remove ABS() calls
      qAtmp = ABS(qOld)
      qBtmp = ABS(ray2D(iH)%q(1))
      rayntmp = ABS(rayn(2))

      !IESCO25: RadusMax = MAX( qAtmp, qBtmp )
      IF (qAtmp.GE.qBtmp) THEN
        RadiusMax = qAtmp
      ELSE
        RadiusMax = qBtmp
      ENDIF
      RadiusMax = RadiusMax * rayntmp / q0

      ! depth limits of beam; IESCO22: a large range of about 1/2 box depth
      IF ( rayt(1).GT.0.5 .OR. rayt(1).LT.-0.5 ) THEN   ! shallow angle ray
        zmin = MIN( xA( 2 ), xB( 2 ) ) - RadiusMax
        zmax = MAX( xA( 2 ), xB( 2 ) ) + RadiusMax
      ELSE                                  ! steep angle ray
        zmin = -HUGE( zmin )
        zmax = +HUGE( zmax )
      ENDIF

      ! compute beam influence for this segment of the ray
      inRcvrRanges=.TRUE.
      RcvrRanges: DO ir = 1, ihop_nRR
!!$TAF store inrcvrranges = iiitape2
        ! is Rr( ir ) contained in [ rA, rB )? Then compute beam influence
        IF ( Pos%RR( ir ).GE.MIN( rA, xB(1) ) .AND. &
             Pos%RR( ir ).LT.MAX( rA, xB(1) ) .AND. &
             inRcvrRanges ) THEN
          x_rcvr( 1, 1:nRz_per_range ) = Pos%RR( ir )

          IF ( Beam%RunType( 5:5 ).EQ.'I' ) THEN ! irregular grid
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
              s = DOT_PRODUCT( x_rcvr( :,iz ) - xA, rayt ) 
              s = s / rlen
              ! normal distance to ray
              n = DOT_PRODUCT( x_rcvr( :,iz ) - xA, rayn )
              n = ABS( n )
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
                delay    = ray2D( iH-1 )%tau + s*dtauds
!      !IESCO25: test some parts of influence, send geninfluence to ihop_cost_modval
!              geninfluence(iH) = delay

                qtmp = ABS(q)
                IF ( ray2D( iH )%c.GT.zeroRL ) THEN
                  Amp = SQRT( ray2D( iH )%c / qtmp )
                ELSE
                  Amp = 0.0
                ENDIF
                Amp = Ratio1 * Amp * ray2D( iH )%Amp
                ! hat function: 1 on center, 0 on edge
                IF ( RadiusMax.NE.0 ) THEN
                  W = ( RadiusMax - n ) / RadiusMax
                ELSE
                  W = 0
                ENDIF
                Amp = Amp*W

                IF (      q.LE.0.0d0 .AND. qOld.GT.0.0d0 &
                     .OR. q.GE.0.0d0 .AND. qOld.LT.0.0d0 ) THEN
                  phaseInt = phase + PI / 2.   ! phase shifts at caustics
                  ! IESCO22: shouldn't this be = phaseInt + PI/2
                ELSE
                  phaseInt = ray2D( iH-1 )%Phase + phase
                ENDIF

                CALL ApplyContribution( iH, iR, iZ, delay, &
                                        amp, phaseint, rcvrdeclangle )

              ENDIF ! IF ( n.LT.RadiusMax )

            ENDIF ! is x_rcvr( 2, iz ) contained in ( zmin, zmax )?

          ENDDO RcvrDepths

        ENDIF ! IF ( Pos%RR( ir ).GE.MIN( rA, rB ) &

        ! bump receiver index, ir, towards rB
        IF ( Pos%RR( ir ).LT.xB(1) ) THEN
          IF ( ir.GE.Pos%nRR ) inRcvrRanges=.FALSE. ! to next step on ray
          irTT = ir + 1                     ! bump right
          IF ( Pos%RR( irTT ).GE.xB(1) ) inRcvrRanges=.FALSE.
        !ELSE
        !   IF ( ir.LE.1              ) inRcvrRanges=.FALSE. ! to next step on ray
        !   irTT = ir - 1                     ! bump left
        !   IF ( Pos%RR( irTT ).LE.rB ) inRcvrRanges=.FALSE.
        !
        ENDIF

        !ir = irTT
      ENDDO RcvrRanges

      rA = xB(1)

    ENDIF ! if duplicate point in ray, skip to next step along the ray

  ENDDO Stepping

  RETURN
  END !SUBROUTINE InfluenceGeoHatCart

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: InfluenceGeoGaussianCart
! !INTERFACE:
  SUBROUTINE InfluenceGeoGaussianCart( myThid )
! !DESCRIPTION:
!   Geometrically-spreading beams with a Gaussian beam in Cartesian
!   coordinates. When Beam%Type(1:1)='B'
!   beam window: kills beams outside e**(-0.5 * ibwin**2 )

! !USES:
!    USE arr_mod, only: nArrival, Arr         !RG
  USE angle_mod, only: Angles

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN )    :: myThid
! !OUTPUT PARAMETERS:

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
  INTEGER :: iH, iZ, iR
  REAL    (KIND=_RL90) :: Ratio1 ! scale factor for a line source
  REAL    (KIND=_RL90) :: W, s, n, phase, rA, xB(2),&
                          q0, q, qOld
  REAL    (KIND=_RL90) :: Amp, phaseInt, RcvrDeclAngle
  COMPLEX (KIND=_RL90) :: delay
  INTEGER, PARAMETER   :: BeamWindow = 4
  INTEGER              :: irT( 1 ), irTT
  REAL (KIND=_RL90)    :: xA( 2 ), rayt( 2 ), rayn( 2 ), x_rcvr( 2 ), &
                        rLen, RadiusMax, zMin, zMax, sigma, lambda, A, dqds
  COMPLEX (KIND=_RL90) :: dtauds
  LOGICAL              :: inRcvrRanges
  REAL (KIND=_RL90) :: qAtmp, qBtmp, rayntmp, qtmp, sigmatmp
!EOP

!$TAF init iGauCart1 = static, (Beam%nSteps-1)
!$TAF init iGauCart2 = static, (Beam%nSteps-1)*ihop_nRR
!$TAF init iGauCart3 = static, (Beam%nSteps-1)*ihop_nRR*nRz_per_range
  q0     = ray2D( 1 )%c / Angles%Dalpha   ! Reference for J = q0 / q
  phase  = 0
  qOld   = ray2D( 1 )%q( 1 )       ! used to track KMAH index
  rA     = ray2D( 1 )%x( 1 )       ! range at start of ray
! init local variables
  W = 0
  s = 0
  n = 0
  xB = [rA, zeroRL]
  q = 0
  Amp = 0
  phaseint = 0
  rcvrdeclangle = 180
  delay = 0

  ! what if never satistified?
  ! what if there is a single receiver (ir = 0 possible)

!  ! irT: find index of first receiver to the right of rA
!  irT = MINLOC( Pos%RR( 1:Pos%nRR ), MASK=Pos%RR( 1:Pos%nRR ).GT.rA )
!  ir  = irT( 1 )
!  ! if ray is left-traveling, get the first receiver to the left of rA
!  IF ( ray2D( 1 )%t( 1 ).LT.0.0d0 .AND. ir.GT.1 ) ir = ir - 1

  ! sqrt( 2 * PI ) represents a sum of Gaussians in free space
  Ratio1 = 1.0 _d 0
  IF ( Beam%RunType( 4:4 ).EQ.'R' ) THEN   ! point source
    Ratio1 = SQRT( ABS( COS( SrcDeclAngle / rad2deg ) ) ) / SQRT( 2.*PI )
  ELSE    ! line  source
    Ratio1 = 1 / SQRT( 2.*PI )
  ENDIF

  Stepping: DO iH = 2, Beam%nSteps
!!$TAF store phase,qold,ra = iGauCart1
    xA = ray2D( iH-1 )%x
    xB = ray2D( iH   )%x

    ! compute normalized tangent (compute it because we need to measure the
    ! step length)
    rayt = xB - xA ! dxds
    IF ( ALL(rayt.EQ.0.) ) THEN
      rlen = 0.0
    ELSE
      rlen = NORM2( rayt )
    ENDIF

    ! if duplicate point in ray, skip to next step along the ray
    IF ( rlen.GE.1.0D3*SPACING( xB(1) ) ) THEN
!!$TAF store rlen,rayt = iGauCart1
      rayt = rayt / rlen
      rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal to ray
      IF ( ALL(rayt.EQ.0.0) ) THEN
        RcvrDeclAngle = 0.0
      ELSE
        RcvrDeclAngle = rad2deg * ATAN2( rayt( 2 ), rayt( 1 ) )
      ENDIF

      q      = ray2D( iH-1 )%q( 1 )
      dqds   = ray2D( iH )%q( 1 ) - q
      dtauds = ray2D( iH )%tau    - ray2D( iH-1 )%tau

      !IESCO22: q only changes signs on direct paths, no top/bot bounces
      IF ( q.LE.0. .AND. qOld.GT.0. .OR. q.GE.0. .AND. qOld.LT.0. ) &
        phase = phase + PI / 2.   ! phase shifts at caustics
      qOld = q

      ! calculate beam width beam radius projected onto vertical line
      qAtmp = ABS(q)
      qBtmp = ABS(ray2D(iH)%q(1))
      rayntmp = ABS(rayn(2))

      IF (qAtmp.GE.qBtmp) THEN
        RadiusMax = qAtmp
      ELSE
        RadiusMax = qBtmp
      ENDIF
      RadiusMax = RadiusMax / q0 / rayntmp
      
      
      lambda   = ray2D( iH-1 )%c / IHOP_freq
      sigmatmp = MIN( 0.2*IHOP_freq*REAL( ray2D( iH )%tau ), &
                      PI*lambda )
      RadiusMax = MAX( RadiusMax, sigmatmp ) 
      ! Note on min: "Weinberg and Keenan suggest limiting a beam to a
      !               point, by imposing a minimum beam width of pilambda."
      !               - Jensen, Comp OA 2011
      ! default is 2 standard deviations of coverage of the Gaussian curve
      RadiusMax = BeamWindow*RadiusMax

      ! depth limits of beam; IESCO22: a large range of about 1/2 box depth
      IF ( rayt(1).GT.0.5 .OR. rayt(1).LT.-0.5 ) THEN   ! shallow angle ray
        zmin   = MIN( xA(2), xB(2) ) - RadiusMax
        zmax   = MAX( xA(2), xB(2) ) + RadiusMax
      ELSE                                 ! steep angle ray
        zmin = -HUGE( zmin )
        zmax = +HUGE( zmax )
      ENDIF

      ! compute beam influence for this segment of the ray
      inRcvrRanges=.TRUE.
      RcvrRanges: DO ir = 1, ihop_nRR
!!$TAF store inrcvrranges = iGauCart2
        ! is Rr( ir ) contained in [ rA, rB )? Then compute beam influence
        IF ( Pos%RR( ir ).GE.MIN( rA, xB(1) ) .AND. &
             Pos%RR( ir ).LT.MAX( rA, xB(1) ) .AND. &
             inRcvrRanges ) THEN
          RcvrDepths: DO iz = 1, nRz_per_range
!!$TAF store arr,nArrival,q = iGauCart3
            IF ( Beam%RunType( 5:5 ).EQ.'I' ) THEN
              x_rcvr = [ Pos%RR( ir ), Pos%RZ( ir ) ]   ! irregular   grid
            ELSE
              x_rcvr = [ Pos%RR( ir ), Pos%RZ( iz ) ]   ! rectilinear grid
            ENDIF

            ! is x_rcvr( 2, iz ) contained in ( zmin, zmax )?
            IF ( x_rcvr(2).GE.zmin .AND. x_rcvr(2).LE.zmax ) THEN
              ! proportional distance along ray
              s = DOT_PRODUCT( x_rcvr - xA, rayt )
              s = s / rlen
              ! normal distance to ray
              n = DOT_PRODUCT( x_rcvr - xA, rayn )
              n = ABS( n )
              ! interpolated amplitude in [meters]
              q = q + s*dqds

              ! beam radius; IESCO22 smaller then previous RadiusMax
              sigma = ABS( q / q0 )
              RadiusMax = MAX( sigma, sigmatmp )
              RadiusMax = BeamWindow*RadiusMax

              IF ( n.LT.RadiusMax ) THEN   ! Within beam window?
#ifdef IHOP_WRITE_OUT
                WRITE(msgBuf,'(A,F10.2)') &
                  "Influence: Eigenray w RadiusMax = ", RadiusMax
                 IF ( IHOP_dumpfreq .GE. 0) &
                   CALL PRINT_MESSAGE( msgbuf, PRTFile, &
                     SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
                ! interpolated delay
                delay    = ray2D( iH-1 )%tau + s*dtauds

                qtmp = ABS(q)
                IF ( ray2D( iH )%c.GT.zeroRL ) THEN
                  Amp = SQRT( ray2D( iH )%c / qtmp )
                ELSE
                  Amp = 0.0
                ENDIF
                Amp      = Ratio1 * Amp * ray2D( iH )%Amp
                ! W : Gaussian decay
                IF ( sigma.NE.0 ) THEN
                  W = EXP( -0.5*( n / sigma )**2 ) / ( sigma**2 )
                ELSE
                  W = 0
                ENDIF
                Amp = Amp*W

                IF ( q.LE.0.0d0 .AND. qOld.GT.0.0d0 .OR. &
                     q.GE.0.0d0 .AND. qOld.LT.0.0d0 ) THEN
                  phaseInt = phase + PI / 2.  ! phase shifts at caustics
                ELSE
                  phaseInt = ray2D( iH-1 )%Phase + phase
                ENDIF

                CALL ApplyContribution( iH, iR, iZ, delay, &
                                        amp, phaseint, rcvrdeclangle )

              ENDIF ! IF ( n.LT.RadiusMax )

            ENDIF ! IF ( x_rcvr( 2, iz ).GE.zmin .AND. x_rcvr( 2, iz ).LE.zmax )
        
          ENDDO RcvrDepths

        ENDIF ! IF ( Pos%RR( ir ).GE.MIN( rA, rB ) ...

        ! bump receiver index, ir, towards rB
        IF ( Pos%RR( ir ).LT.xB(1) ) THEN
          IF ( ir.GE.Pos%nRR ) inRcvrRanges=.FALSE. ! to next step on ray
          irTT = ir + 1                     ! bump right
          IF ( Pos%RR( irTT ).GE.xB(1) ) inRcvrRanges=.FALSE.
        !ELSE
        !  IF ( ir.LE.1              ) inRcvrRanges=.FALSE. ! to next step on ray
        !  irTT = ir - 1                     ! bump left
        !  IF ( Pos%RR( irTT ).LE.rB ) inRcvrRanges=.FALSE.
        ENDIF

        !ir = irTT
      ENDDO RcvrRanges

      rA = xB(1)

    ENDIF ! if duplicate point in ray, skip to next step along the ray

  ENDDO Stepping

  RETURN
  END !SUBROUTINE InfluenceGeoGaussianCart

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ApplyContribution
! !INTERFACE:
  SUBROUTINE ApplyContribution( iH, iR, iZ, tau, Amp, phaseInt, arrAngle )
! !DESCRIPTION:
!   Apply the contribution of the current ray to the pressure field U.

! !USES:
  USE writeray, only: WriteRayOutput
  USE arr_mod,  only: AddArr, U
  USE ihop_mod, only: afreq, RAYFile, DELFile, nMax

! !INPUT PARAMETERS:
  COMPLEX (KIND=_RL90), INTENT( IN ) :: tau
  REAL (KIND=_RL90), INTENT( IN ) :: Amp, phaseInt, arrAngle
  INTEGER, INTENT( IN ) :: iH, iR, iZ
! !OUTPUT PARAMETERS:

! !LOCAL VARIABLES:
! tmpDelay :: temporary array for delay values
! uPoint :: complex pressure field
  COMPLEX  :: uPoint
  INTEGER :: i
  REAL(KIND=_RL90) :: tmpX(nMax), tmpY(nMax), tmpDelay(nMax)
!EOP

! initiate temporary arrays
  tmpX = 0.
  tmpY = 0.
  tmpDelay = 0.
  DO i = 1,iH
    tmpX(i) = ray2D(i)%x(1)
    tmpY(i) = ray2D(i)%x(2)
    tmpDelay(i) = REAL(ray2D(i)%tau)
  ENDDO

! write rays and delays
  SELECT CASE( Beam%RunType( 1:1 ) )
  CASE ( 'E' )                ! eigenrays
    uPoint=uPoint
    CALL WriteRayOutput( RAYFile, iH,         &
      tmpX(1:iH), tmpY(1:iH),                 &
      ray2D(iH)%nTopBnc, ray2D(iH)%nBotBnc )

  CASE ( 'e' )                ! eigenrays AND arrivals
    uPoint=uPoint
    CALL WriteRayOutput( RAYFile, iH,         &
      tmpX(1:iH), tmpY(1:iH),                 &
      ray2D(iH)%nTopBnc, ray2D(iH)%nBotBnc )
    IF (writeDelay) THEN
      CALL WriteRayOutput( DELFile, iH,       &
        tmpDelay(1:iH), tmpY(1:iH),           &
        ray2D(iH)%nTopBnc, ray2D(iH)%nBotBnc )
    ENDIF

    CALL AddArr( afreq, iZ, iR, Amp, phaseInt, tau,  &
      arrAngle, ray2D(iH)%nTopBnc, ray2D(iH)%nBotBnc )

  CASE ( 'A', 'a' )           ! arrivals
    uPoint=uPoint
    CALL AddArr( afreq, iZ, iR, Amp, phaseInt, tau,  &
      arrAngle, ray2D(iH)%nTopBnc, ray2D(iH)%nBotBnc )

  CASE ( 'C' )                ! coherent TL
    uPoint = uPoint + CMPLX( Amp * EXP( -oneCMPLX * ( afreq * tau - phaseInt ) ) )

  CASE ( 'S', 'I' )                ! incoherent/semicoherent TL
    IF ( Beam%Type( 1:1 ).EQ.'B' ) THEN   ! Gaussian beam
      uPoint = uPoint + SNGL( SQRT( 2. * PI ) &
            * ( Amp * EXP( AIMAG( afreq * tau ) ) )**2 )
    ELSE
      uPoint = uPoint + SNGL( &
               ( Amp * EXP( AIMAG( afreq * tau ) ) )**2 )
    ENDIF

  CASE DEFAULT                ! incoherent/semicoherent TL
    IF ( Beam%Type( 1:1 ).EQ.'B' ) THEN   ! Gaussian beam
      uPoint = uPoint + SNGL( SQRT( 2. * PI ) &
            * ( Amp * EXP( AIMAG( afreq * tau ) ) )**2 )
    ELSE
      uPoint = uPoint + SNGL( &
               ( Amp * EXP( AIMAG( afreq * tau ) ) )**2 )
    ENDIF

  END SELECT ! SELECT CASE( Beam%RunType( 1:1 ) )

  ! upload pressure field
  U(iZ, iR) = uPoint

  RETURN
  END !SUBROUTINE ApplyContribution

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ScalePressure
! !INTERFACE:
  SUBROUTINE ScalePressure( c, r, U, nRz, nRcvr )
! !DESCRIPTION:
!   Scale the pressure field U according to the run type and range.

! !USES:
  USE ihop_mod,  only: Beam
  USE angle_mod, only: Angles

! !INPUT PARAMETERS:
! c :: nominal sound speed, m/s
! r :: Rr ranges, m
! U :: complex pressure field
! nRz :: number of depths in the pressure field
! nRcvr :: number of ranges in the pressure field
! RunType :: run type, e.g. 'C' for Cerveny Gaussian beams in Cartesian coordinates
! freq :: source frequency, Hz
  INTEGER,           INTENT( IN    ) :: nRz, nRcvr
  REAL (KIND=_RL90), INTENT( IN    ) :: c, r( nRcvr )
  COMPLEX,           INTENT( INOUT ) :: U( nRz, nRcvr )
! !OUTPUT PARAMETERS: U

! !LOCAL VARIABLES:
! const :: scale factor for field
! factor :: scale factor for cylindrical spreading
  REAL (KIND=_RL90) :: const, factor
  INTEGER :: iR
!EOP

  ! Compute scale factor for field
  SELECT CASE ( Beam%RunType( 2:2 ) )
  CASE ( 'C' )   ! Cerveny Gaussian beams in Cartesian coordinates
    const = -Angles%Dalpha * SQRT( IHOP_freq ) / c
  CASE ( 'R' )   ! Cerveny Gaussian beams in Ray-centered coordinates
    const = -Angles%Dalpha * SQRT( IHOP_freq ) / c
  CASE DEFAULT
    const = -1.0
  END SELECT

  ! If incoherent run, convert intensity to pressure
  IF ( Beam%RunType( 1:1 ).NE.'C' ) U = SQRT( REAL( U ) )

  ! scale and/or incorporate cylindrical spreading
  Ranges: DO iR = 1, nRcvr
    IF ( Beam%RunType( 4:4 ).EQ.'X' ) THEN   ! line source
      factor = -4.0 * SQRT( PI ) * const
    ELSE                                  ! point source
      IF ( r( iR ).EQ.0 ) THEN
        factor = 0.0D0 ! avoid /0 at origin, return pressure = 0
      ELSE
        IF ( r(iR).GT.0 ) THEN 
          factor = const / SQRT( r( iR ) )
        ELSE
          factor = zeroRL
        ENDIF
      ENDIF
    ENDIF
    U( :, iR ) = SNGL( factor ) * U( :, iR )

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


  Ax  = ABS( x )

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
