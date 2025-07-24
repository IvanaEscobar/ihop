#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: step
MODULE step
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   This module implements 2D step along a single ray.

! !USES:
  USE ihop_mod, only: Beam, ray2DPt
  USE ssp_mod,  only: evalSSP, SSP, iSegz, iSegr
  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef ALLOW_CTRL
# include "CTRL_FIELDS.h"
#endif

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC Step2D
!=======================================================================

! == Module variables == None
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R Step2D
! S/R ReduceStep2D
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: Step2D
! !INTERFACE:
  SUBROUTINE Step2D( ray0, ray2, Topx, Topn, Botx, Botn, myThid )
! !DESCRIPTION:
!   March a single step along the ray in 2D.

! !USES: None

! !INPUT PARAMETERS:
! ray0 :: the ray at the beginning of the step
! Topx, Topn :: the top boundary coordinate and normal
! Botx, Botn :: the bottom boundary coordinate and normal
! myThid :: my thread ID
  TYPE( ray2DPt ),     INTENT( INOUT ) :: ray0, ray2
  REAL ( KIND=_RL90 ), INTENT( IN )    :: Topx( 2 ), Topn( 2 ), &
                                          Botx( 2 ), Botn( 2 )
  INTEGER, INTENT(IN) :: myThid
! !OUTPUT PARAMETERS: ray2

! !LOCAL VARIABLES:
! ray1 :: the ray at the end of the first phase of the step
! iSegz0, iSegr0 :: the segment indices for the ray
! gradc0, gradc1, gradc2 :: the gradient of the sound speed at the ray points
! c0, cimag0, crr0, crz0, czz0, csq0, cnn0_csq0 :: sound speed and its derivatives at ray0
! c1, cimag1, crr1, crz1, czz1, csq1, cnn1_csq1 :: sound speed and its derivatives at ray1
! c2, cimag2, crr2, crz2, czz2 :: sound speed and its derivatives at ray2
! urayt0, urayt1 :: unit tangent vectors at ray0 and ray1
! h, halfh, hw0, hw1 :: step size and its components
! ray2n :: outward unit normal at ray2
! RM, RN :: reflection coefficients
! gradcjump :: the jump in the gradient of the sound speed across an interface
! cnjump, csjump :: the jumps in the normal and tangential components of the sound speed gradient
! w0, w1 :: blending weights for the two phases of the step
! rho :: density at the ray points
  TYPE( ray2DPt ) :: ray1
  INTEGER         :: iSegz0, iSegr0
  REAL ( KIND=_RL90 ) :: gradc0( 2 ), gradc1( 2 ), gradc2( 2 ), &
                        c0, cimag0, crr0, crz0, czz0, csq0, cnn0_csq0, &
                        c1, cimag1, crr1, crz1, czz1, csq1, cnn1_csq1, &
                        c2, cimag2, crr2, crz2, czz2, &
                        urayt0( 2 ), urayt1( 2 ), &
                        h, halfh, hw0, hw1, ray2n( 2 ), RM, RN, &
                        gradcjump( 2 ), cnjump, csjump, w0, w1, rho
!EOP

  ! IESCO25: notes
  ! Does a single step along the ray
  ! x denotes the ray coordinate, (r,z)
  ! t denotes the scaled tangent to the ray (previously (rho, zeta))
  ! c * t would be the unit tangent

  ! The numerical integrator used here is a version of the polygon (a.k.a.
  ! midpoint, leapfrog, or Box) method, and similar
  ! to the Heun (second order Runge-Kutta method).
  ! However, it's modified to allow for a dynamic step change, while
  ! preserving second-order accuracy.

  ! *** Phase 1 (an Euler step)
  CALL evalSSP( ray0%x, c0, cimag0, gradc0, crr0, crz0, czz0, rho, myThid )

  csq0      = c0 * c0
  cnn0_csq0 = crr0*ray0%t( 2 )**2 - 2.0*crz0*ray0%t( 1 )*ray0%t( 2 ) &
              + czz0*ray0%t( 1 )**2 !IESCO22: chain rule map n deriv to r,z
  iSegz0    = iSegz     ! make note of current layer
  iSegr0    = iSegr

  h = Beam%deltas       ! initially set the step h, to the basic one, deltas
  urayt0 = c0 * ray0%t  ! unit tangent

  ! reduce h to land on boundary
  CALL ReduceStep2D( ray0%x, urayt0, iSegz0, iSegr0, Topx, Topn, Botx, &
                    Botn, h )
  halfh = 0.5 * h   ! first step of the modified polygon method is a half step

  ! Euler march forward
  ray1%x = ray0%x + halfh * urayt0
  ray1%t = ray0%t - halfh * gradc0 / csq0
  ray1%p = ray0%p - halfh * cnn0_csq0 * ray0%q
  ray1%q = ray0%q + halfh * c0        * ray0%p !IESCO22: q /= 0 for 'G' beam

  ! *** Phase 2 (update step size, and Polygon march forward)
  CALL evalSSP( ray1%x, c1, cimag1, gradc1, crr1, crz1, czz1, rho, myThid )
  csq1      = c1 * c1
  cnn1_csq1 = crr1*ray1%t( 2 )**2 - 2.0*crz1*ray1%t( 1 )*ray1%t( 2 ) &
              + czz1*ray1%t( 1 )**2

  ! BUG: The Munk test case with a horizontally launched ray caused problems.
  ! The ray vertexes on an interface and can ping-pong around that interface.
  ! Have to be careful in that case about big changes to the stepsize (that
  ! invalidate the leap-frog scheme) in phase II.
  ! A modified Heun or Box method could also work.

  urayt1 = c1 * ray1%t   ! unit tangent

  CALL ReduceStep2D( ray0%x, urayt1, iSegz0, iSegr0, Topx, Topn, Botx, &
                    Botn, h ) ! reduce h to stay in domain

  ! use blend of f' based on proportion of a full step used.
  w1  = h / ( 2.0d0 * halfh ) ! h/h_old
  w0  = 1.0d0 - w1
  hw0 = h * w0
  hw1 = h * w1                ! hw0 + hw1 = h

  ! blended march forward
  ray2%x   = ray0%x   + hw0 * urayt0              + hw1 * urayt1
  ray2%t   = ray0%t   - hw0 * gradc0 / csq0       - hw1 * gradc1 / csq1
  ray2%p   = ray0%p   - hw0 * cnn0_csq0 * ray0%q  - hw1 * cnn1_csq1 * ray1%q
  ray2%q   = ray0%q   + hw0 * c0        * ray0%p  + hw1 * c1        * ray1%p
  ray2%tau = ray0%tau + hw0 / CMPLX( c0, cimag0, KIND=_RL90 ) &
             + hw1 / CMPLX( c1, cimag1, KIND=_RL90 )

  ray2%Amp       = ray0%Amp
  ray2%Phase     = ray0%Phase
  ray2%NumTopBnc = ray0%NumTopBnc
  ray2%NumBotBnc = ray0%NumBotBnc

  ! If we crossed an interface, apply linear jump condition
  CALL evalSSP( ray2%x, c2, cimag2, gradc2, crr2, crz2, czz2, rho, myThid )
  ray2%c = c2

  IF ( iSegz.NE.iSegz0 .OR. iSegr.NE.iSegr0 ) THEN
    gradcjump =  gradc2 - gradc0
    ray2n     = [ -ray2%t( 2 ), ray2%t( 1 ) ]   ! outward unit normal / c

    cnjump    = DOT_PRODUCT( gradcjump, ray2n  )
    csjump    = DOT_PRODUCT( gradcjump, ray2%t )

    IF ( iSegz.NE.iSegz0 ) THEN         ! crossing in depth
      ! RM is tan( alpha ) where alpha is the angle of incidence
      RM = +ray2%t( 1 ) / ray2%t( 2 )
    ELSE ! iSegr /= iSegr0              ! crossing in range
      ! RM is tan( alpha ) where alpha is the angle of incidence
      RM = -ray2%t( 2 ) / ray2%t( 1 )
    ENDIF

    RN     = RM * ( 2*cnjump - RM*csjump ) / c2
    ray2%p = ray2%p - ray2%q*RN

  ENDIF ! IF ( iSegz.NE.iSegz0 .OR. iSegr.NE.iSegr0 )

  RETURN
  END !SUBROUTINE Step2D

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: ReduceStep2D
! !INTERFACE:
  SUBROUTINE ReduceStep2D( x0, urayt, iSegz0, iSegr0, Topx, Topn, &
                           Botx, Botn, h )
! !DESCRIPTION:
!   Calculate a reduced step size, h, that lands on any points where the
!   environment leaves water, or on the top or bottom boundary.

! !USES:
  USE ihop_mod, only: iSmallStepCtr
  USE bdry_mod, only: rTopSeg, rBotSeg

! !INPUT PARAMETERS:
! x0 :: the ray coordinate at the beginning of the step
! urayt :: the unit tangent vector at the ray coordinate
! iSegz0, iSegr0 :: SSP segment indices for the ray
! Topx, Topn :: the top boundary coordinate and normal
! Botx, Botn :: the bottom boundary coordinate and normal
  REAL (KIND=_RL90), INTENT( IN    ) :: x0( 2 ), urayt( 2 )
  INTEGER,           INTENT( IN    ) :: iSegz0, iSegr0
  REAL (KIND=_RL90), INTENT( IN    ) :: Topx( 2 ), Topn( 2 )
  REAL (KIND=_RL90), INTENT( IN    ) :: Botx( 2 ), Botn( 2 )
  REAL (KIND=_RL90), INTENT( INOUT ) :: h ! reduced step size
! !OUTPUT PARAMETERS: h

! !LOCAL VARIABLES:
! hInt, hTop, hBot, hSeg, hBoxr, hBoxz :: trial step sizes
! x, d, d0, rSeg :: vectors used to calculate step sizes
  REAL (KIND=_RL90) :: hInt, hTop, hBot, hSeg, hBoxr, hBoxz
  REAL (KIND=_RL90) :: x( 2 ), d( 2 ), d0( 2 ), rSeg( 2 )
!EOP

!$TAF init reducestep2d = static, 50
!$TAF store h = reducestep2d

  x = x0 + h * urayt ! take a trial Euler step

  ! interface crossing in depth
  hInt = huge( hInt ) ! Largest _RL90 number available
  IF ( ABS( urayt( 2 ) ).GT.EPSILON( hInt ) ) THEN
    IF       ( SSP%Z( iSegz0   ).GT.x(  2 ) ) THEN
      hInt = ( SSP%Z( iSegz0   ) - x0( 2 ) ) / urayt( 2 )
    ELSEIF   ( SSP%Z( iSegz0+1 ).LT.x(  2 ) ) THEN
      hInt = ( SSP%Z( iSegz0+1 ) - x0( 2 ) ) / urayt( 2 )
    ELSE
      ! Do nothing
      hInt = hInt
    ENDIF
  ENDIF

  ! top crossing
  hTop = huge( hTop )
  d = x - Topx             ! vector from top to ray
  IF ( DOT_PRODUCT( Topn, d ).GT.EPSILON( hTop ) ) THEN
    d0 = x0 - Topx         ! vector from top    node to ray origin
    hTop = -DOT_PRODUCT( d0, Topn ) / DOT_PRODUCT( urayt, Topn )
  ENDIF

  ! bottom crossing
  hBot = huge( hBot )
  d = x - Botx             ! vector from bottom to ray
  IF ( DOT_PRODUCT( Botn, d ).GT.EPSILON( hBot ) ) THEN
    d0 = x0 - Botx         ! vector from bottom node to ray origin
    hBot = -DOT_PRODUCT( d0, Botn ) / DOT_PRODUCT( urayt, Botn )
  ENDIF

  ! top or bottom segment crossing in range
  rSeg( 1 ) = MAX( rTopSeg( 1 ), rBotSeg( 1 ) )
  rSeg( 2 ) = MIN( rTopSeg( 2 ), rBotSeg( 2 ) )

!$TAF store rseg = reducestep2d
  IF ( SSP%Type.EQ.'Q' ) THEN ! Quad: 2D range-dependent SSP
    rSeg( 1 ) = MAX( rSeg( 1 ), SSP%Seg%R( iSegr0     ) )
    rSeg( 2 ) = MIN( rSeg( 2 ), SSP%Seg%R( iSegr0 + 1 ) )
  ENDIF

  ! interface crossing in range
  hSeg = huge( hSeg )
  IF ( ABS( urayt( 1 ) ).GT.EPSILON( hSeg ) ) THEN
    IF        ( x(  1 ).LT.rSeg( 1 ) ) THEN
      hSeg = -( x0( 1 ) - rSeg( 1 ) ) / urayt( 1 )
    ELSEIF    ( x(  1 ).GT.rSeg( 2 ) ) THEN
      hSeg = -( x0( 1 ) - rSeg( 2 ) ) / urayt( 1 )
    ELSE
      ! Do nothing
      hSeg = hSeg
    ENDIF
  ENDIF

  ! ray mask using a box centered at ( 0, 0 )
  hBoxr    = huge( hBoxr )
  hBoxz    = huge( hBoxz )

  IF ( ABS( x( 1 ) ).GT.Beam%Box%R ) hBoxr = ( Beam%Box%R - ABS( x0( 1 ) ) ) / ABS( urayt( 1 ) )
  IF ( ABS( x( 2 ) ).GT.Beam%Box%Z ) hBoxz = ( Beam%Box%Z - ABS( x0( 2 ) ) ) / ABS( urayt( 2 ) )

  h = MIN( h, hInt, hTop, hBot, hSeg, hBoxr, hBoxz )  ! take limit set by shortest distance to a crossing
  IF ( h.LT.1.0d-4 * Beam%deltas ) THEN   ! is it taking an infinitesimal step?
    h = 1.0d-4 * Beam%deltas            ! make sure we make some motion
    iSmallStepCtr = iSmallStepCtr + 1   ! keep a count of the number of sequential small steps
  ELSE
    iSmallStepCtr = 0   ! didn't do a small step so reset the counter
  ENDIF

  RETURN
  END !SUBROUTINE ReduceStep2D

END !MODULE step