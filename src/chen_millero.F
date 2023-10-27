#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: CHEN_MILLERO
C     !INTERFACE:
      _RL FUNCTION CHEN_MILLERO(i,j,k,bi,bj,
     I                       myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | FXN: CHEN_MILLERO
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
#include "DYNVARS.h"
#include "EOS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine Arguments ==
C     myThid :: Thread number for this instance of the routine.
      INTEGER myThid

#ifdef ALLOW_IHOP 
C     !FUNCTIONS:
      _RL SW_TEMP
      EXTERNAL SW_TEMP

C     !LOCAL VARIABLES:
C     == Local variables ==
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER bi,bj
      INTEGER i,j,k
      _RL c0,cw,a,b,d
      _RL pres,sal,temp
CEOP

C initiate parameters
      c0   = 1402.388 _d 0
      cw   = 0. _d 0
      a    = 0. _d 0
      b    = 0. _d 0
      d    = 0. _d 0
      pres = 0. _d 0
      sal  = 0. _d 0
      temp = 0. _d 0

C pressure in dbar (for SW_TEMP)
C     let k=0 define the surface, a required level for ihop
      IF (k.EQ.0) THEN
       temp = SW_TEMP( SALT(i,j,1,bi,bj),
     &          THETA(i,j,1,bi,bj), pres, zeroRL )
       sal  = SALT(i,j,1,bi,bj)
      ELSE IF (k.GT.Nr) THEN
       pres = rhoConst*( totPhiHyd(i,j,Nr,bi,bj)
     &        - ( rF(Nr+1)*1.05 )*gravity )*SItodBar
       temp = SW_TEMP( SALT(i,j,Nr,bi,bj),
     &          THETA(i,j,Nr,bi,bj), pres, zeroRL )
       sal  = SALT(i,j,Nr,bi,bj)
      ELSE
       pres = rhoConst*( totPhiHyd(i,j,k,bi,bj)
     &        - rC(k)*gravity )*SItodBar
       temp = SW_TEMP( SALT(i,j,k,bi,bj),
     &          THETA(i,j,k,bi,bj), pres, zeroRL )
       sal  = SALT(i,j,k,bi,bj)
      ENDIF
C convert pressure to bar for Chen and Millero algorithm
      pres = pres/(1. _d 1)
      cw   = c0                      + 5.0383 _d 0*temp      +
     &       -5.81090 _d -2*temp**2  + 3.3432 _d -4*temp**3  + 
     &       -1.47797 _d -6*temp**4  + 3.1419 _d -9*temp**5  +
     &       (1.53563 _d -1          + 6.8999 _d -4*temp     +
     &        -8.1829 _d -6*temp**2  + 1.3632 _d -7*temp**3  +
     &        -6.1260 _d -10*temp**4)*pres                   + 
     &       (3.12600 _d -5          - 1.7111 _d -6*temp     +
     &        2.59860 _d -8*temp**2  - 2.5353 _d -10*temp**3 +
     &        1.04150 _d -12*temp**4)*pres**2                +
     &       (-9.7729 _d -9          + 3.8513 _d -10*temp    +  
     &        -2.3654 _d -12*temp**2)*pres**3
      a    = 1.389 _d 0              - 1.2620 _d -2*temp     + 
     &       7.166 _d -5*temp**2     + 2.0080 _d -6*temp**3  +
     &       -3.21 _d -8*temp**4                             +
     &       (9.47420 _d -5          - 1.2583 _d -5*temp     +
     &        -6.4928 _d -8*temp**2  + 1.0515 _d -8*temp**3  +
     &        -2.0142 _d -10*temp**4)*pres                   +
     &       (-3.9064 _d -7          + 9.1061 _d -9*temp     +
     &        -1.6009 _d -10*temp**2                         + 
     &        7.99400 _d -12*temp**3)*pres**2                +
     &       (1.10000 _d -10         + 6.6510 _d -12*temp    +
     &        -3.3910 _d -13*temp**2)*pres**3
      b    = -1.9220 _d -2           - 4.4200 _d -5*temp     + 
     &       (7.3637 _d -5           + 1.7950 _d -7*temp)*pres
      d    = 1.727 _d -3             - 7.9836 _d -6*pres 
      CHEN_MILLERO  = cw + a*sal + b*sal**(3./2.) + d*sal**2

#endif /* ALLOW_IHOP */ 

      RETURN
      END
