#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE influence
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

  ! Compute the beam influence, i.e. the contribution of a single beam to the
  ! complex pressure
  ! mbp 12/2018, based on much older subroutines

  USE ihop_mod,     only: rad2deg, oneCMPLX, PRTFile, maxN, Beam, ray2D, &
                          SrcDeclAngle, nRz_per_range
  USE srPos_mod,    only: Pos
! sspMod used to construct image beams in the Cerveny style beam routines

! ! USES
  implicit none
!  == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

  PRIVATE

! public interfaces
!=======================================================================

    public InfluenceGeoHatRayCen, InfluenceGeoGaussianCart, &
           InfluenceGeoHatCart, ScalePressure!, InfluenceSGB

!=======================================================================

  INTEGER,              PRIVATE :: iz, ir, iS
  REAL    (KIND=_RL90), PRIVATE :: Ratio1 = 1.0D0 ! scale factor for a line source
  REAL    (KIND=_RL90), PRIVATE :: W, s, n, Amp, phase, phaseInt, &
                                   q0, q, qold, RcvrDeclAngle, rA, rB
  COMPLEX (KIND=_RL90), PRIVATE :: delay

CONTAINS
  SUBROUTINE InfluenceGeoHatRayCen( U, dalpha, myThid )
    use arr_mod, only: nArr, Arr         !RG

    ! Geometrically-spreading beams with a hat-shaped beam in ray-centered
    ! coordinates

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN    ) :: dalpha ! take-off angle radians
    COMPLEX,           INTENT( INOUT ) :: U( nRz_per_range, Pos%nRR )   ! complex pressure field
    INTEGER              :: irA, irB, II
    REAL (KIND=_RL90)    :: nA, nB, zr, L, dq( Beam%Nsteps - 1 )
    REAL (KIND=_RL90)    :: znV( Beam%Nsteps ), rnV( Beam%Nsteps ), &
                            RcvrDeclAngleV ( Beam%Nsteps )
    COMPLEX (KIND=_RL90) :: dtau( Beam%Nsteps-1 )
    LOGICAL :: skip_step

    !!! need to add logic related to nRz_per_range

!$TAF init iRayCen0  = 'influence_iraycen'
!$TAF init iRayCen1  = static, (Beam%Nsteps-1)*nRz_per_range
!$TAF init iRayCen2  = static, (Beam%Nsteps-1)*ihop_nRR*nRz_per_range

!$TAF store ray2d = iRayCen0

    q0   = ray2D( 1 )%c / Dalpha   ! Reference for J = q0 / q

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
    IF ( Beam%RunType( 4:4 ) == 'R' ) &
        Ratio1 = SQRT( ABS( COS( SrcDeclAngle / rad2deg ) ) )

    ray2D( 1:Beam%Nsteps )%Amp = Ratio1 * SQRT( ray2D( 1:Beam%Nsteps )%c ) &
                        * ray2D( 1:Beam%Nsteps )%Amp   ! pre-apply some scaling

    RcvrDepths: DO iz = 1, nRz_per_range
       zR = Pos%RZ( iz )

       phase = 0.0
       qOld  = ray2D( 1 )%q( 1 ) ! used to track KMAH index

       ! If normal is parallel to horizontal receiver line
       IF ( ABS( znV( 1 ) ) < 1D-6 ) THEN
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
       END IF

       Stepping: DO iS = 2, Beam%Nsteps
!$TAF store ira,irb,na,nb,phase,qold,ra = iRayCen1
         skip_step = .FALSE.

         ! Compute ray-centered coordinates, (znV, rnV)

         ! If normal is parallel to TL-line, skip to the next step on ray
         IF (ABS(znV(iS)) < 1D-10) THEN
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
           IF (ABS(ray2D(iS)%x(1) - ray2D(iS - 1)%x(1)) < &
               1.0D3 * SPACING(ray2D(iS)%x(1)) &
               .OR. irA == irB) THEN
             rA = rB
             nA = nB
             irA = irB
             skip_step = .TRUE.
           END IF
         END IF

         IF (.NOT. skip_step) THEN
           !!! this should be pre-computed
           q = ray2D(iS - 1)%q(1)
           ! if phase shifts at caustics
           IF ((q <= 0.0D0 .AND. qOld > 0.0D0) .OR. &
               (q >= 0.0D0 .AND. qOld < 0.0D0)) &
               phase = phase + PI / 2.0D0
           qOld = q

           RcvrDeclAngle = RcvrDeclAngleV(iS)

           ! *** Compute contributions to bracketted receivers ***
           II = 0
           IF (irB <= irA) II = 1   ! going backwards in range

           ! Compute influence for each receiver
           DO ir = irA + 1 - II, irB + II, SIGN(1, irB - irA)
!$TAF store Arr(:,ir,iz),nArr(ir,iz) = iRayCen2
             W = (Pos%RR(ir) - rA) / (rB - rA)  ! relative range between rR
             n = ABS(nA + W * (nB - nA))
             q = ray2D(iS - 1)%q(1) + W * dq(iS - 1)  ! interpolated amplitude
             L = ABS(q) / q0   ! beam radius

             IF (n < L) THEN  ! in beam window: update delay, Amp, phase
!$TAF store w = iRayCen2
               delay = ray2D(iS - 1)%tau + W * dtau(iS - 1)
               Amp = ray2D(iS)%Amp / SQRT(ABS(q))
               W = (L - n) / L  ! hat function: 1 on center, 0 on edge
               Amp = Amp * W
               phaseInt = ray2D(iS - 1)%Phase + phase
               !!! this should be precomputed
               IF ((q <= 0.0D0 .AND. qOld > 0.0D0) .OR. &
                   (q >= 0.0D0 .AND. qOld < 0.0D0)) &
                   phaseInt = phase + PI / 2.0D0  ! phase shifts at caustics

               CALL ApplyContribution(U(iz, ir))
             END IF
           END DO
         END IF

         rA = rB
         nA = nB
         irA = irB
       END DO Stepping
    END DO RcvrDepths

  RETURN
  END !SUBROUTINE InfluenceGeoHatRayCen

  ! **********************************************************************!

  SUBROUTINE InfluenceGeoHatCart( U, Dalpha, myThid )
    use arr_mod, only: nArr, Arr         !RG
    ! Geometric, hat-shaped beams in Cartesisan coordinates

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN    ) :: & ! take-off angle, radians
                                          Dalpha   ! angular spacing
    COMPLEX,           INTENT( INOUT ) :: U( nRz_per_range, Pos%nRR ) ! complex pressure field
    INTEGER              :: irT(1), irTT ! irT needs size of 1, see MINLOC
    REAL (KIND=_RL90)    :: x_ray( 2 ), rayt( 2 ), rayn( 2 ), &
                            x_rcvr( 2, nRz_per_range ), rLen, RadiusMax, &
                            zMin, zMax, dqds
    COMPLEX (KIND=_RL90) :: dtauds
    LOGICAL              :: inRcvrRanges

!$TAF init iiitape1 = static, (Beam%Nsteps-1)
!$TAF init iiitape2 = static, (Beam%Nsteps-1)*ihop_nRR
!$TAF init iiitape3 = static, (Beam%Nsteps-1)*ihop_nRR*nRz_per_range

    q0           = ray2D( 1 )%c / Dalpha   ! Reference for J = q0 / q
    phase        = 0.0
    qOld         = ray2D( 1 )%q( 1 )       ! old KMAH index
    rA           = ray2D( 1 )%x( 1 )       ! range at start of ray, typically 0

    ! find index of first receiver to the right of rA
    irT = MINLOC( Pos%RR( 1 : Pos%nRR ), MASK = Pos%RR( 1 : Pos%nRR ) > rA )
    ir  = irT( 1 )
    ! if ray is left-traveling, get the first receiver to the left of rA
    IF ( ray2D( 1 )%t( 1 ) < 0.0d0 .AND. ir > 1 ) ir = ir - 1

    ! point source: the default option
    Ratio1 = 1.0d0          !RG
    IF ( Beam%RunType( 4 : 4 ) == 'R' ) &
        Ratio1 = SQRT( ABS( COS( SrcDeclAngle / rad2deg ) ) )

    Stepping: DO iS = 2, Beam%Nsteps

!$TAF store phase,qold,ra = iiitape1

       rB     = ray2D( iS   )%x( 1 )
       x_ray  = ray2D( iS-1 )%x

       ! compute normalized tangent (we need to measure the step length)
       rayt = ray2D( iS )%x - x_ray
       IF ( ALL(rayt==0.0) ) THEN
          rlen = 0.0
       ELSE
          rlen = NORM2( rayt )
       ENDIF

       ! if duplicate point in ray, skip to next step along the ray
       IF ( rlen .GE. 1.0D3 * SPACING( ray2D( iS )%x( 1 ) ) ) THEN

!$TAF store rlen,rayt= iiitape1

        rayt = rayt / rlen                    ! unit tangent of ray @ A
        rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal  of ray @ A
        IF ( ALL(rayt==0.0) ) THEN
           RcvrDeclAngle = 0.0
        ELSE
           RcvrDeclAngle = rad2deg * ATAN2( rayt( 2 ), rayt( 1 ) )
        ENDIF

        q      = ray2D( iS-1 )%q( 1 )
        dqds   = ray2D( iS   )%q( 1 ) - q
        dtauds = ray2D( iS   )%tau    - ray2D( iS-1 )%tau

        !IESCO22: q only changes signs on direct paths, no top/bot bounces
        IF( q <= 0.0 .AND. qOld > 0.0 .OR. q >= 0.0 .AND. qOld < 0.0 )&
            phase = phase + PI / 2.   ! phase shifts at caustics
        qOld = q

        ! Radius calc from beam radius projected onto vertical line
        RadiusMax = MAX( ABS( q ), ABS( ray2D( iS )%q( 1 ) ) ) &
                    / q0 / ABS( rayt( 1 ) ) ! IESCO24: AKA rayn( 2 )

        ! depth limits of beam; IESCO22: a large range of about 1/2 box depth
        IF ( ABS( rayt( 1 ) ) > 0.5 ) THEN   ! shallow angle ray
           zmin   = min( x_ray( 2 ), ray2D( iS )%x( 2 ) ) - RadiusMax
           zmax   = max( x_ray( 2 ), ray2D( iS )%x( 2 ) ) + RadiusMax
        ELSE                                 ! steep angle ray
           zmin = -HUGE( zmin )
           zmax = +HUGE( zmax )
        END IF

        ! compute beam influence for this segment of the ray
        inRcvrRanges=.TRUE.
        RcvrRanges: DO ir = 1, ihop_nRR
!$TAF store inrcvrranges = iiitape2
           ! is Rr( ir ) contained in [ rA, rB )? Then compute beam influence
           IF ( Pos%RR( ir ) >= MIN( rA, rB ) &
                .AND. Pos%RR( ir ) < MAX( rA, rB ) &
                .AND. inRcvrRanges ) THEN

              x_rcvr( 1, 1:nRz_per_range ) = Pos%RR( ir )
              IF ( Beam%RunType( 5 : 5 ) == 'I' ) THEN ! irregular grid
                 x_rcvr( 2, 1 ) = Pos%RZ( ir )
              ELSE ! default: rectilinear grid
                 x_rcvr( 2, 1:nRz_per_range ) = Pos%RZ( 1:nRz_per_range )
              END IF

              RcvrDepths: DO iz = 1, nRz_per_range
!$TAF store Arr(:,ir,iz),nArr(ir,iz),q = iiitape3
                 ! is x_rcvr( 2, iz ) contained in ( zmin, zmax )?
                 IF (       x_rcvr( 2, iz ) .GE. zmin &
                      .AND. x_rcvr( 2, iz ) .LE. zmax ) THEN
                    ! normalized proportional distance along ray
                    s = DOT_PRODUCT( x_rcvr( :, iz ) - x_ray, rayt ) / rlen
                    ! normal distance to ray
                    n = ABS( DOT_PRODUCT( x_rcvr( :, iz ) - x_ray, rayn ) )
                    ! interpolated amplitude in [meters]
                    q = q + s*dqds
                    ! beam radius; IESCO22 smaller then previous RadiusMax
                    RadiusMax = ABS( q / q0 )

                    IF ( n < RadiusMax ) THEN
#ifdef IHOP_WRITE_OUT
                        WRITE(msgBuf,'(A,F10.2)') &
                            "Influence: Eigenray w RadiusMax = ", RadiusMax
                        IF ( IHOP_dumpfreq .GE. 0) &
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

                       IF (      q <= 0.0d0 .AND. qOld > 0.0d0 &
                            .OR. q >= 0.0d0 .AND. qOld < 0.0d0 ) THEN
                        phaseInt = phase + PI / 2.   ! phase shifts at caustics
                        ! IESCO22: shouldn't this be = phaseInt + PI/2
                       ELSE
                        phaseInt = ray2D( iS-1 )%Phase + phase
                       END IF
                       CALL ApplyContribution( U( iz, ir ) )
                    END IF
                 END IF ! is x_rcvr( 2, iz ) contained in ( zmin, zmax )?
              END DO RcvrDepths
           END IF

           ! bump receiver index, ir, towards rB
           IF ( Pos%RR( ir ) < rB ) THEN
              IF ( ir >= Pos%nRR        ) inRcvrRanges=.FALSE. ! to next step on ray
              irTT = ir + 1                     ! bump right
              IF ( Pos%RR( irTT ) >= rB ) inRcvrRanges=.FALSE.
           !ELSE
           !   IF ( ir <= 1              ) inRcvrRanges=.FALSE. ! to next step on ray
           !   irTT = ir - 1                     ! bump left
           !   IF ( Pos%RR( irTT ) <= rB ) inRcvrRanges=.FALSE.
           END IF
           !ir = irTT
        END DO RcvrRanges

        rA = rB
       END IF ! if duplicate point in ray, skip to next step along the ray
    END DO Stepping

  RETURN
  END !SUBROUTINE InfluenceGeoHatCart

  ! **********************************************************************!

  SUBROUTINE InfluenceGeoGaussianCart( U, Dalpha, myThid )
    use arr_mod, only: nArr, Arr         !RG

    ! Geometric, Gaussian beams in Cartesian coordintes

    ! beam window: kills beams outside e**(-0.5 * ibwin**2 )
  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    INTEGER,           PARAMETER       :: BeamWindow = 4
    REAL (KIND=_RL90), INTENT( IN    ) :: dalpha ! take-off angle, angular spacing
    COMPLEX,           INTENT( INOUT ) :: U( nRz_per_range, Pos%nRR )  ! complex pressure field
    INTEGER              :: irT( 1 ), irTT
    REAL (KIND=_RL90)    :: x_ray( 2 ), rayt( 2 ), rayn( 2 ), x_rcvr( 2 ), &
                            rLen, RadiusMax, zMin, zMax, sigma, lambda, A, dqds
    COMPLEX (KIND=_RL90) :: dtauds
    LOGICAL              :: inRcvrRanges

!$TAF init iGauCart1 = static, (Beam%Nsteps-1)
!$TAF init iGauCart2 = static, (Beam%Nsteps-1)*ihop_nRR
!$TAF init iGauCart3 = static, (Beam%Nsteps-1)*ihop_nRR*nRz_per_range

    q0           = ray2D( 1 )%c / Dalpha   ! Reference for J = q0 / q
    phase        = 0
    qOld         = ray2D( 1 )%q( 1 )       ! used to track KMAH index
    rA           = ray2D( 1 )%x( 1 )       ! range at start of ray

    ! what if never satistified?
    ! what if there is a single receiver (ir = 0 possible)

    ! irT: find index of first receiver to the right of rA
    irT = MINLOC( Pos%RR( 1 : Pos%nRR ), MASK = Pos%RR( 1 : Pos%nRR ) > rA )
    ir  = irT( 1 )

    ! if ray is left-traveling, get the first receiver to the left of rA
    IF ( ray2D( 1 )%t( 1 ) < 0.0d0 .AND. ir > 1 ) ir = ir - 1

    ! sqrt( 2 * PI ) represents a sum of Gaussians in free space
    IF ( Beam%RunType( 4 : 4 ) == 'R' ) THEN   ! point source
       Ratio1 = SQRT( ABS( COS( SrcDeclAngle / rad2deg ) ) ) / SQRT( 2. * PI )
    ELSE    ! line  source
       Ratio1 = 1 / SQRT( 2. * PI )
    END IF

    Stepping: DO iS = 2, Beam%Nsteps
!$TAF store phase,qold,ra = iGauCart1
       rB    = ray2D( iS     )%x( 1 )
       x_ray = ray2D( iS - 1 )%x

       ! compute normalized tangent (compute it because we need to measure the
       ! step length)
       rayt = ray2D( iS )%x - ray2D( iS - 1 )%x
       IF ( ALL(rayt==0.0) ) THEN
          rlen = 0.0
       ELSE
          rlen = NORM2( rayt )
       ENDIF
       ! if duplicate point in ray, skip to next step along the ray
       IF ( rlen .GE. 1.0D3 * SPACING( ray2D( iS )%x( 1 ) ) ) THEN

!$TAF store rlen,rayt = iGauCart1

        rayt = rayt / rlen
        rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal to ray
        IF ( ALL(rayt==0.0) ) THEN
           RcvrDeclAngle = 0.0
        ELSE
           RcvrDeclAngle = rad2deg * ATAN2( rayt( 2 ), rayt( 1 ) )
        ENDIF

        q      = ray2D( iS-1 )%q( 1 )
        dqds   = ray2D( iS )%q( 1 ) - q
        dtauds = ray2D( iS )%tau    - ray2D( iS-1 )%tau

        !IESCO22: q only changes signs on direct paths, no top/bot bounces
        IF ( q <= 0.0 .AND. qOld > 0.0 .OR. q >= 0.0 .AND. qOld < 0.0 ) &
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
        IF ( ABS( rayt( 1 ) ) > 0.5 ) THEN   ! shallow angle ray
           zmin   = min( ray2D( iS-1 )%x( 2 ), ray2D( iS )%x( 2 ) ) - RadiusMax
           zmax   = max( ray2D( iS-1 )%x( 2 ), ray2D( iS )%x( 2 ) ) + RadiusMax
        ELSE                                 ! steep angle ray
           zmin = -HUGE( zmin )
           zmax = +HUGE( zmax )
        END IF

        ! compute beam influence for this segment of the ray
        inRcvrRanges=.TRUE.
        RcvrRanges: DO ir = 1, ihop_nRR
!$TAF store inrcvrranges = iGauCart2
           ! is Rr( ir ) contained in [ rA, rB )? Then compute beam influence
           IF ( Pos%RR( ir ) >= MIN( rA, rB ) &
                .AND. Pos%RR( ir ) < MAX( rA, rB ) &
                .AND. inRcvrRanges ) THEN

              RcvrDepths: DO iz = 1, nRz_per_range
!$TAF store arr,nArr,q = iGauCart3
                 IF ( Beam%RunType( 5 : 5 ) == 'I' ) THEN
                    x_rcvr = [ Pos%RR( ir ), Pos%RZ( ir ) ]   ! irregular   grid
                 ELSE
                    x_rcvr = [ Pos%RR( ir ), Pos%RZ( iz ) ]   ! rectilinear grid
                 END IF
                 ! is x_rcvr( 2, iz ) contained in ( zmin, zmax )?
                 IF (       x_rcvr( 2 ) .GE. zmin &
                      .AND. x_rcvr( 2 ) .LE. zmax ) THEN

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

                    IF ( n < BeamWindow*sigma ) THEN   ! Within beam window?
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
                       IF ( q <= 0.0d0 .AND. qOld > 0.0d0 &
                            .OR. q >= 0.0d0 .AND. qOld < 0.0d0 ) &
                        phaseInt = phase + PI / 2.  ! phase shifts at caustics

                       CALL ApplyContribution( U( iz, ir ) )
                    END IF
                 END IF ! is x_rcvr( 2, iz ) contained in ( zmin, zmax )?
              END DO RcvrDepths
           END IF

            ! bump receiver index, ir, towards rB
            IF ( Pos%RR( ir ) < rB ) THEN
               IF ( ir >= Pos%nRR        ) inRcvrRanges=.FALSE. ! to next step on ray
               irTT = ir + 1                     ! bump right
               IF ( Pos%RR( irTT ) >= rB ) inRcvrRanges=.FALSE.
            ELSE
               IF ( ir <= 1              ) inRcvrRanges=.FALSE. ! to next step on ray
               irTT = ir - 1                     ! bump left
               IF ( Pos%RR( irTT ) <= rB ) inRcvrRanges=.FALSE.
            END IF
            !ir = irTT
        END DO RcvrRanges

        rA = rB
       END IF ! if duplicate point in ray, skip to next step along the ray
    END DO Stepping

  RETURN
  END !SUBROUTINE InfluenceGeoGaussianCart

! **************************************************************************** !
  SUBROUTINE ApplyContribution( U )
    USE writeray, only: WriteRayOutput
    USE arr_mod,  only: AddArr
    USE ihop_mod, only: afreq, RAYFile, DELFile, maxN

    COMPLEX, INTENT( INOUT ) :: U
    REAL(KIND=_RL90) :: tmpDelay(maxN)

    tmpDelay = 0.
    SELECT CASE( Beam%RunType( 1:1 ) )
      CASE ( 'E' )                ! eigenrays
        U=U
        tmpDelay = 0.
        CALL WriteRayOutput( RAYFile, iS, ray2D%x(1), ray2D%x(2), &
            ray2D(iS)%NumTopBnc, ray2D(iS)%NumBotBnc )
      CASE ( 'e' )                ! eigenrays AND arrivals
        U=U
        tmpDelay = 0.
        CALL WriteRayOutput( RAYFile, iS, ray2D%x(1), ray2D%x(2), &
            ray2D(iS)%NumTopBnc, ray2D(iS)%NumBotBnc )
        IF (writeDelay) THEN
          tmpDelay = REAL(ray2D%tau)
          CALL WriteRayOutput( DELFile, iS, tmpDelay, ray2D%x(2), &
              ray2D(iS)%NumTopBnc, ray2D(iS)%NumBotBnc )
        END IF

        CALL AddArr( afreq, iz, ir, Amp, phaseInt, delay, SrcDeclAngle, &
                     RcvrDeclAngle, ray2D( iS )%NumTopBnc, &
                     ray2D( iS )%NumBotBnc )
      CASE ( 'A', 'a' )           ! arrivals
        U=U
        tmpDelay = 0.
        CALL AddArr( afreq, iz, ir, Amp, phaseInt, delay, SrcDeclAngle, &
                     RcvrDeclAngle, ray2D( iS )%NumTopBnc, &
                     ray2D( iS )%NumBotBnc )
      CASE ( 'C' )                ! coherent TL
        tmpDelay = 0.
        U = U + CMPLX( Amp * EXP( -oneCMPLX * ( afreq * delay - phaseInt ) ) )
      CASE ( 'S', 'I' )                ! incoherent/semicoherent TL
        tmpDelay = 0.
        IF ( Beam%Type( 1:1 ) == 'B' ) THEN   ! Gaussian beam
          U = U + SNGL( SQRT( 2. * PI ) &
                  * ( Amp * EXP( AIMAG( afreq * delay ) ) )**2 )
        ELSE
          U = U + SNGL( &
                    ( Amp * EXP( AIMAG( afreq * delay ) ) )**2 )
        END IF
      CASE DEFAULT                ! incoherent/semicoherent TL
        tmpDelay = 0.
        IF ( Beam%Type( 1:1 ) == 'B' ) THEN   ! Gaussian beam
          U = U + SNGL( SQRT( 2. * PI ) &
                  * ( Amp * EXP( AIMAG( afreq * delay ) ) )**2 )
        ELSE
          U = U + SNGL( &
                    ( Amp * EXP( AIMAG( afreq * delay ) ) )**2 )
        END IF
    END SELECT

  RETURN
  END !SUBROUTINE ApplyContribution

! **********************************************************************!

  SUBROUTINE ScalePressure( Dalpha, c, r, U, NRz, Nr, RunType, freq )

    ! Scale the pressure field

    INTEGER,           INTENT( IN    ) :: NRz, Nr
    REAL (KIND=_RL90), INTENT( IN    ) :: r( Nr )   ! Rr ranges
    REAL (KIND=_RL90), INTENT( IN    ) :: Dalpha, freq, c ! angular spacing between rays, source frequency, nominal sound speed
    COMPLEX,           INTENT( INOUT ) :: U( NRz, Nr )    ! Pressure field
    CHARACTER (LEN=5), INTENT( IN    ) :: RunType
    REAL (KIND=_RL90)                  :: const, factor

    ! Compute scale factor for field
    SELECT CASE ( RunType( 2 : 2 ) )
    CASE ( 'C' )   ! Cerveny Gaussian beams in Cartesian coordinates
       const = -Dalpha * SQRT( freq ) / c
    CASE ( 'R' )   ! Cerveny Gaussian beams in Ray-centered coordinates
       const = -Dalpha * SQRT( freq ) / c
    CASE DEFAULT
       const = -1.0
    END SELECT

    ! If incoherent run, convert intensity to pressure
    IF ( RunType( 1 : 1 ) /= 'C' ) U = SQRT( REAL( U ) )

    ! scale and/or incorporate cylindrical spreading
    Ranges: DO ir = 1, Nr
       IF ( RunType( 4 : 4 ) == 'X' ) THEN   ! line source
          factor = -4.0 * SQRT( PI ) * const
       ELSE                                  ! point source
          IF ( r ( ir ) == 0 ) THEN
             factor = 0.0D0         ! avoid /0 at origin, return pressure = 0
          ELSE
             factor = const / SQRT( ABS( r( ir ) ) )
          END IF
       END IF
       U( :, ir ) = SNGL( factor ) * U( :, ir )
    END DO Ranges

  RETURN
  END !SUBROUTINE ScalePressure

  ! **********************************************************************!

  REAL (KIND=_RL90) FUNCTION Hermite( x, x1, x2 )

    ! Calculates a smoothing function based on the h0 hermite cubic
    ! x is the point where the function is to be evaluated
    ! returns:
    ! [  0, x1  ] = 1
    ! [ x1, x2  ] = cubic taper from 1 to 0
    ! [ x2, inf ] = 0

    REAL (KIND=_RL90 ), INTENT( IN  ) :: x, x1, x2
    REAL (KIND=_RL90 )                :: Ax, u

    Ax  = ABS( x  )

    IF ( Ax <= x1 ) THEN
       Hermite = 1.0d0
    ELSE IF ( Ax >= x2 ) THEN
       Hermite = 0.0d0
    ELSE
       u       = ( Ax - x1 ) / ( x2 - x1 )
       Hermite = ( 1.0d0 + 2.0d0 * u ) * ( 1.0d0 - u ) ** 2
    END IF

  RETURN
  END !FUNCTION Hermite

END !MODULE influence
