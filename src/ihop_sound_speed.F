#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#include "IHOP_OPTIONS.h"

CBOP
C     !ROUTINE: IHOP_SOUND_SPEED
C     !INTERFACE:
      SUBROUTINE IHOP_SOUND_SPEED(
     I                       myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R IHOP_SOUND_SPEED
C     | o Diagnose speed of sound in seawater
C     |   from the algorithm by Chen and Millero (1977).
C     |   This is NOT the sound-speed that can be derived from
C     |   the equation of state (EOS). It is independent of
C     |   the model setup specific EOS.
C     |
C     | o Reference:
C     | C. Chen and F. J. Millero, "Speed of sound in seawater at 
C     | high pressures,"
C     | J. Acoust. Soc. Am. 672.5, 1129-1135 (1977).
C     *==========================================================*
C     \ev
C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
C#include "SURFACE.h"
#include "DYNVARS.h"
#include "EOS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine Arguments ==
C     myThid :: Thread number for this instance of the routine.
      INTEGER myThid

#ifdef ALLOW_IHOP 
C     !FUNCTIONS:
C      _RL CHEN_MILLERO
C      EXTERNAL CHEN_MILLERO
C
      _RL CHEN_MILLERO
      EXTERNAL CHEN_MILLERO

#ifdef ALLOW_DIAGNOSTICS
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
#endif /* ALLOW_DIAGNOSTICS */

C     !LOCAL VARIABLES:
C     == Local variables ==
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER bi,bj
      INTEGER i,j,k
      INTEGER ii, jj, kk, bii,bjj
      _RL c0,cw,a,b,d
      _RL pres,sal,temp
      LOGICAL calc_soundSpeed
CEOP

C      calc_soundSpeed = .FALSE.

C--   switch on this flag if Sound-Speed needed in any cost fxn
c#ifdef ALLOW_CSOUND_COST
      calc_soundSpeed = .TRUE.
c#endif

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        calc_soundSpeed = calc_soundSpeed
     &               .OR. DIAGNOSTICS_IS_ON( 'ihop_ssp', myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

      IF ( calc_soundSpeed ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k=1,Nr
          DO j=1,sNy
           DO i=1,sNx
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

C      IF ( calc_soundSpeed ) THEN
C       c0   = 1402.388 _d 0
C       cw   = 0. _d 0
C       a    = 0. _d 0
C       b    = 0. _d 0
C       d    = 0. _d 0
C       pres = 0. _d 0
C       sal  = 0. _d 0
C       temp = 0. _d 0
C       DO bj=myByLo(myThid),myByHi(myThid)
C        DO bi=myBxLo(myThid),myBxHi(myThid)
C         DO k=1,Nr
C          DO j=1,sNy
C           DO i=1,sNx
C            ihop_ssp(i,j,k,bi,bj) = 0.0 _d 0
C           ENDDO
C          ENDDO
C         ENDDO
C         DO k=1,Nr
C          DO j=1,sNy
C           DO i=1,sNx
C            IF ( maskC(i,j,k,bi,bj).EQ.oneRS ) THEN
CC pressure in dbar (for SW_TEMP)
C             pres = rhoConst*( totPhiHyd(i,j,k,bi,bj)
C     &               - rC(k)*gravity )*SItodBar
C             temp = SW_TEMP( SALT(i,j,k,bi,bj),
C     &                 THETA(i,j,k,bi,bj), pres, zeroRL )
C             sal  = SALT(i,j,k,bi,bj)
CC convert pressure to bar for Chen and Millero algorithm
C             pres = pres/(1. _d 1)
C             cw   = c0                      + 5.0383 _d 0*temp      +
C     &              -5.81090 _d -2*temp**2  + 3.3432 _d -4*temp**3  + 
C     &              -1.47797 _d -6*temp**4  + 3.1419 _d -9*temp**5  +
C     &              (1.53563 _d -1          + 6.8999 _d -4*temp     +
C     &               -8.1829 _d -6*temp**2  + 1.3632 _d -7*temp**3  +
C     &               -6.1260 _d -10*temp**4)*pres                   + 
C     &              (3.12600 _d -5          - 1.7111 _d -6*temp     +
C     &               2.59860 _d -8*temp**2  - 2.5353 _d -10*temp**3 +
C     &               1.04150 _d -12*temp**4)*pres**2                +
C     &              (-9.7729 _d -9          + 3.8513 _d -10*temp    +  
C     &               -2.3654 _d -12*temp**2)*pres**3
C             a    = 1.389 _d 0              - 1.2620 _d -2*temp     + 
C     &              7.166 _d -5*temp**2     + 2.0080 _d -6*temp**3  +
C     &              -3.21 _d -8*temp**4                             +
C     &              (9.47420 _d -5          - 1.2583 _d -5*temp     +
C     &               -6.4928 _d -8*temp**2  + 1.0515 _d -8*temp**3  +
C     &               -2.0142 _d -10*temp**4)*pres                   +
C     &              (-3.9064 _d -7          + 9.1061 _d -9*temp     +
C     &               -1.6009 _d -10*temp**2                         + 
C     &               7.99400 _d -12*temp**3)*pres**2                +
C     &              (1.10000 _d -10         + 6.6510 _d -12*temp    +
C     &               -3.3910 _d -13*temp**2)*pres**3
C             b    = -1.9220 _d -2           - 4.4200 _d -5*temp     + 
C     &              (7.3637 _d -5           + 1.7950 _d -7*temp)*pres
C             d    = 1.727 _d -3             - 7.9836 _d -6*pres 
C             ihop_ssp(i,j,k,bi,bj) = cw + a*sal + b*sal**(3./2.)    +
C     &                               d*sal**2
C            ENDIF
C           ENDDO
C          ENDDO
C         ENDDO
C        ENDDO
C       ENDDO

C-    end-if calc_soundSpeed
      ENDIF

#endif /* ALLOW_IHOP */ 

      RETURN
      END
