#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_SOUND_SPEED

C !INTERFACE:
      SUBROUTINE IHOP_SOUND_SPEED(
     I                       myThid )

C     !DESCRIPTION: 
C     o Diagnose speed of sound in seawater
C       from the algorithm by Chen and Millero (1977).
C       This is NOT the sound-speed that can be derived from
C       the equation of state (EOS). It is independent of
C       the model setup specific EOS.
C     
C     *==========================================================*
C     | o Reference:
C     | C. Chen and F. J. Millero, "Speed of sound in seawater at
C     | high pressures,"
C     | J. Acoust. Soc. Am. 672.5, 1129-1135 (1977).
C     *==========================================================*

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EOS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

C     !INPUT PARAMETERS:
C     myThid :: Thread number for this instance of the routine.
      INTEGER myThid
C     !OUTPUT PARAMETERS: None

C     !FUNCTIONS:
#ifdef ALLOW_IHOP

      _RL CHEN_MILLERO
      EXTERNAL CHEN_MILLERO

#ifdef ALLOW_DIAGNOSTICS
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
#endif /* ALLOW_DIAGNOSTICS */

C     !LOCAL VARIABLES:
      INTEGER bi,bj
      INTEGER i,j,k
c     LOGICAL calc_soundSpeed
CEOP

C      calc_soundSpeed = .FALSE.

C--   switch on this flag if Sound-Speed needed in any cost fxn
c#ifdef ALLOW_CSOUND_COST
c     calc_soundSpeed = .TRUE.
c#endif

#ifdef ALLOW_DIAGNOSTICS
c     IF ( useDiagnostics ) THEN
c       calc_soundSpeed = calc_soundSpeed
c    &               .OR. DIAGNOSTICS_IS_ON( 'SVEL    ', myThid )
c    &               .OR. DIAGNOSTICS_IS_ON( 'SLD     ', myThid )
c     ENDIF
#endif /* ALLOW_DIAGNOSTICS */

c     IF ( calc_soundSpeed ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            ihop_ssp(i,j,k,bi,bj) = 0.0 _d 0
           ENDDO
          ENDDO
         ENDDO
         DO k=1,Nr
          DO j=1,sNy
           DO i=1,sNx
            IF ( maskC(i,j,k,bi,bj).EQ.oneRS ) THEN
             ihop_ssp(i,j,k,bi,bj) = CHEN_MILLERO(i,j,k,bi,bj,myThid)
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO

C-    end-if calc_soundSpeed
c     ENDIF

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics
     &     .AND.DIAGNOSTICS_IS_ON( 'SLD     ', myThid) ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         CALL IHOP_CALC_SLD(bi,bj,myThid)

        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_IHOP */

      RETURN
      END
