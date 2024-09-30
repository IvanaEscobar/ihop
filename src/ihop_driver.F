#include "IHOP_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: IHOP_DRIVER

C !INTERFACE: ==========================================================
      SUBROUTINE IHOP_DRIVER( myTime, myIter, myThid )

C !DESCRIPTION:
C Calculate ray travel times, tau

C !FORTRAN90 USE MODULES: ===============================================
      USE ihop, only: ihop_main
#ifdef ALLOW_AUTODIFF_TAMC
      USE arr_mod, only: arr
      USE bdry_mod, only: bdry
      USE ssp_mod, only: ssp
      USE srpos_mod, only: pos
      use angle_mod, only: angles
      use ihop_mod, only: beam, ray2d
#endif

      IMPLICIT NONE
C !GLOBAL VARIABLES: ===============================================
#ifdef ALLOW_IHOP
!   == Global variables ==
# include "SIZE.h"
# include "GRID.h"
# include "EEPARAMS.h"
# include "EESUPPORT.h"
# include "PARAMS.h"
# include "IHOP_SIZE.h"
# include "IHOP.h"
#endif /* ALLOW_IHOP */
#ifdef ALLOW_AUTODIFF_TAMC
# include "tamc.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myTime   :: Current time in simulation
C  myIter   :: Current time-step number
C  myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid

C !OUTPUT PARAMETERS: ==================================================

C !LOCAL VARIABLES: ====================================================
      INTEGER t
#ifdef ALLOW_AUTODIFF_TAMC
      integer tkey
#endif

CEOP

#ifdef ALLOW_IHOP
      DO t=1,nts
CML# ifdef ALLOW_AUTODIFF_TAMC
CML       tkey = nts*(ikey_dynamics-1) + t
CML!$TAF STORE arr, bdry, ray2d   = comlev1_ihop_nts, key = tkey
CML!$TAF STORE angles%dalpha      = comlev1_ihop_nts, key = tkey
CML!$TAF STORE bdry%bot, bdry%top = comlev1_ihop_nts, key = tkey
CML!$TAF STORE beam%nsteps        = comlev1_ihop_nts, key = tkey
CML!$TAF STORE pos%rz, pos%sz     = comlev1_ihop_nts, key = tkey
CML!$TAF STORE ssp%c              = comlev1_ihop_nts, key = tkey
CML# endif
        IF ( IHOP_iter(t).GE.0 ) THEN
            IF ( IHOP_iter(t).EQ.myIter ) THEN
                CALL IHOP_MAIN( myTime, myIter, myThid )
            ENDIF
        ENDIF
      ENDDO
#endif /* ALLOW_IHOP */

      RETURN
      END
