#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: subTab_mod
MODULE subTab_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   Module for organizing ihop output

! If x(3) = -999.9 then subtabulation is performed
! i.e., a vector is generated with Nx points in [ x(1), x(2) ]
! If x(2) = -999.9 then x(1) is repeated into x(2)
! tabulates info in outputfile
! mbp 1/2015

! !USES:
  IMPLICIT NONE

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC SubTab
!=======================================================================

! == Module variables ==
  INTEGER, PRIVATE :: ix

! == Interfaces ==
  INTERFACE SubTab
     MODULE PROCEDURE SubTab_sngl, SubTab_dble
  END INTERFACE SubTab

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R SubTab_sngl
! S/R SubTab_dble
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: SubTab_sngl
! !INTERFACE:
  SUBROUTINE SubTab_sngl( x, Nx )
! !DESCRIPTION:
!   Subroutine to tabulate a single precision vector x of length Nx.

! !USES: None

! !INPUT PARAMETERS:
!   x  :: single precision vector of length Nx.
!   Nx :: The length of the vector x.
  REAL,    INTENT( INOUT ) :: x( Nx )
  INTEGER, INTENT( IN )    :: Nx
! !OUTPUT PARAMETERS: x

! !LOCAL VARIABLES:
! deltax :: increment between points in the vector x.
  REAL :: deltax
!EOP

  IF ( Nx.GE.3 ) THEN
    IF ( x( 3 ).EQ.-999.9 ) THEN   ! testing for equality here is dangerous
      IF ( x( 2 ).EQ.-999.9 ) x( 2 ) = x( 1 )
      deltax      = ( x( 2 ) - x( 1 ) ) / ( Nx - 1 )
      x( 1 : Nx ) = x( 1 ) + [ ( ix, ix = 0, Nx - 1 ) ] * deltax
    ENDIF
  ENDIF

  RETURN
  END SUBROUTINE !SubTab_sngl

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: SubTab_dble
! !INTERFACE:
  SUBROUTINE SubTab_dble( x, Nx )
! !DESCRIPTION:
!   Subroutine to tabulate a double precision vector x of length Nx.

! !USES: None

! !INPUT PARAMETERS:
!   x  :: double precision vector of length Nx.
!   Nx :: The length of the vector x.
  REAL (KIND=_RL90),  INTENT( INOUT ) :: x( Nx )
  INTEGER,            INTENT( IN    ) :: Nx
! !OUTPUT PARAMETERS: x
! !LOCAL VARIABLES:
! deltax :: increment between points in the vector x.
  REAL (KIND=_RL90)   :: deltax
!EOP

  IF ( Nx.GE.3 ) THEN
    IF ( x( 3 ).EQ.-999.9 ) THEN   ! testing for equality here is dangerous
      IF ( x( 2 ).EQ.-999.9 ) x( 2 ) = x( 1 )
      deltax      = ( x( 2 ) - x( 1 ) ) / ( Nx - 1 )
      x( 1 : Nx ) = x( 1 ) + [ ( ix, ix = 0, Nx - 1 ) ] * deltax
    ENDIF
  ENDIF

  RETURN
  END SUBROUTINE ! SubTab_dble

END MODULE ! subTab_mod