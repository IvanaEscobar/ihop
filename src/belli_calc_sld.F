#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_CALC_SLD
C !INTERFACE:
      SUBROUTINE IHOP_CALC_SLD( bi, bj, myThid )

C     !DESCRIPTION:
C     S/R IHOP_CALC_SLD
C     o Diagnose the Oceanic Sonic Layer depth
C     Note: output "hMixLayer" is in "r" unit, i.e., in Pa
C           when using P-coordinate.

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_GMREDI
# include "GMREDI.h"
#endif
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "IHOP.h"
#endif

C     !INPUT PARAMETERS:
C     == Routine Arguments ==
C     bi,bj     :: tile indices
C     myThid    :: my Thread Id number
      INTEGER bi, bj
      INTEGER myThid
C     !OUTPUT PARAMETERS: None

C     !FUNCTIONS:
#ifdef ALLOW_DIAGNOSTICS
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
#endif /* ALLOW_DIAGNOSTICS */

C     !LOCAL VARIABLES:
C     i, j, k   :: Loop counters
C     kDir      :: downward k increment
C     kU, kL    :: loop ranges, Up and Low
      INTEGER i, j, k
      INTEGER kDir, deltaK
      INTEGER kU, kL
      INTEGER kMin, kMax
      INTEGER ihop_kMin(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER ihop_kMax(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      LOGICAL calcSonicLayerDepth
      _RL     min_val, max_val
CEOP

      calcSonicLayerDepth = .FALSE.
#ifdef ALLOW_GMREDI
      IF ( useGMRedi .AND. .NOT.useKPP ) THEN
       calcSonicLayerDepth = GM_useSubMeso 
     &       .OR. GM_taper_scheme.EQ.'fm07'
     &       .OR. GM_useBatesK3d
      ENDIF
#endif
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics.AND. .NOT.calcSonicLayerDepth ) THEN
        calcSonicLayerDepth = DIAGNOSTICS_IS_ON('SLD     ',myThid)
      ENDIF
#endif
      IF ( calcSonicLayerDepth ) THEN

        IF ( usingPCoords ) THEN
         kDir    = -1
         deltaK  =  1
        ELSE
         kDir    =  1
         deltaK  =  0
        ENDIF

C--   First (only) method :
C     where the global maximum above the global minimum of sound
C     sound speed is, Z[sup(inf(ssp))]
C     (see Helber, Barron, et al.  JGR 2008 for default criterion)

C-    Z-coord: (kU,kL,kDir) = (2,Nr,1) ; P-coord: (kU,kL,kDir)= (Nr-1,1,-1)
        kU = 2 + deltaK*(Nr-3)
        kL = Nr - deltaK*(Nr-1)

C step 1: determine depth of global minimum
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          min_val = 9999.9
          kMin = -1*kDir

          DO k=kU,kL,kDir
           IF ( ihop_ssp(i,j,k,bi,bj).LT.min_val ) THEN
               min_val = ihop_ssp(i,j,k,bi,bj)
               kMin = k
           ENDIF
          ENDDO
          ihop_kMin (i,j) = kMin

         ENDDO
        ENDDO

C step 2: find vertical level of the truncated global maximum
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          max_val = -9999.9
          kMax = -1*kDir
          kMin = ihop_kMin(i,j)

          DO k=kU,kMin,kDir
           IF (ihop_ssp(i,j,k,bi,bj).GT.max_val ) THEN
               max_val = ihop_ssp(i,j,k,bi,bj)
               kMax = k
           ENDIF
          ENDDO
          ihop_kMax (i,j) = kMax

         ENDDO
        ENDDO

C step 3: save depth where maximum ssp exists, making SLD
C         a grid dependent diagnostic 
C IESCO25: update to determine depth of maximum SSP between
C         k levels... see calc_oce_mxlayer.F
C a) calc ssp(rC_UP/DOWN)
C b) line search max ssp UP/DOWN
C c) store depth of max ssp
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          kMax = ihop_kMax(i,j)
          ihop_sld(i,j,bi,bj) = gravitySign * rC(kMax) 
     &                          * maskC(i,j,kMax,bi,bj)

         ENDDO
        ENDDO

#ifdef ALLOW_DIAGNOSTICS
       IF ( useDiagnostics ) THEN
        CALL DIAGNOSTICS_FILL( ihop_sld, 'SLD     ',
     &                         0, 1, 1, bi, bj, myThid )
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */

C--   end if calcSonicLayerDepth
      ENDIF

      RETURN
      END