#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_INIT_VARIA

C !INTERFACE:
      SUBROUTINE IHOP_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize IHOP variables

C     !USES:
      IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

C     !INPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
C     !OUTPUT PARAMETERS: None

#ifdef ALLOW_IHOP
C     !LOCAL VARIABLES:
      INTEGER i,j,k,bi,bj
CEOP

      DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)

C     Initialize arrays in common blocks :
C     3-dim.
#ifdef IHOP_3D_STATE
          DO k=1,Nr
            DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
                ihop_ssp(i,j,k,bi,bj) = 0. _d 0
              ENDDO
            ENDDO
          ENDDO
#endif /* IHOP_3D_STATE */

C     2-dim.
#ifdef IHOP_2D_STATE
          DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
              ihop_sld(i,j,bi,bj) = 0. _d 0
            ENDDO
          ENDDO
#endif /* IHOP_2D_STATE */

C- end bi,bj loops
        ENDDO
      ENDDO

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN
#ifdef IHOP_3D_STATE
#endif /* IHOP_3D_STATE */

#ifdef IHOP_2D_STATE
#endif /* IHOP_2D_STATE */

      ELSE
C-    restart from a pickup:
        CALL IHOP_READ_PICKUP( nIter0, myThid )

C-    end start-from-iter-zero if/else block
      ENDIF

#endif /* ALLOW_IHOP */

      RETURN
      END