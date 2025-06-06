#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef ALLOW_OPENAD
# include "OPENAD_OPTIONS.h"
#endif
c#ifdef ALLOW_AUTODIFF
c# include "AUTODIFF_OPTIONS.h"
c#endif
#ifdef ALLOW_GENERIC_ADVDIFF
# include "GAD_OPTIONS.h"
#endif
#ifdef ALLOW_GGL90
# include "GGL90_OPTIONS.h"
#endif
#ifdef ALLOW_GMREDI
# include "GMREDI_OPTIONS.h"
#endif
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif
#ifdef ALLOW_GCHEM
# include "GCHEM_OPTIONS.h"
#endif
#ifdef ALLOW_DIC
# include "DIC_OPTIONS.h"
#endif
#ifdef ALLOW_IHOP
# include "IHOP_OPTIONS.h"
#endif
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif
#ifdef ALLOW_STREAMICE
# include "STREAMICE_OPTIONS.h"
#endif
#ifdef ALLOW_COST
# include "COST_OPTIONS.h"
#endif
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
c#ifdef ALLOW_ECCO
c# include "ECCO_OPTIONS.h"
c#endif
#ifdef ALLOW_TAPENADE
# ifdef ALLOW_SALT_PLUME
# include "SALT_PLUME_OPTIONS.h"
# endif
#ifdef ALLOW_ECCO
# include "ECCO_OPTIONS.h"
#endif
c# ifdef ALLOW_THSICE
c# include "THSICE_OPTIONS.h"
c# endif
#endif
#ifdef ALLOW_BLING
# include "BLING_OPTIONS.h"
#endif
CBOP
C     !ROUTINE: MAIN_DO_LOOP
C     !INTERFACE:
      SUBROUTINE MAIN_DO_LOOP( iloop, myTime, myIter, myThid )

C     !DESCRIPTION:
C     *================================================================*
C     | SUBROUTINE MAIN_DO_LOOP
C     | o Interface S/R that calls FORWARD_STEP
C     *================================================================*
C     | This S/R contains the model main "do-loop" for OpenAD built
C     | (i.e., ALLOW_OPENAD defined), with a call to FORWARD_STEP
C     | inside the loop.
C     | But without OpenAD, this S/R just calls directly FORWARD_STEP
C     | (in this case the main do-loop is done in S/R THE_MAIN_LOOP).
C     *================================================================*

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

c**************************************
#if (defined ALLOW_OPENAD || defined ALLOW_TAPENADE)

# ifdef ALLOW_CTRL
#  include "CTRL_SIZE.h"
#  include "CTRL.h"
#  include "CTRL_DUMMY.h"
#  include "CTRL_GENARR.h"
#  include "CTRL_OBCS.h"
# endif
# ifdef ALLOW_COST
#  include "cost.h"
# endif

# ifdef ALLOW_PTRACERS
#  include "PTRACERS_SIZE.h"
#  include "PTRACERS_FIELDS.h"
#  include "PTRACERS_START.h"
# endif
# ifdef ALLOW_GCHEM
#  include "GCHEM_SIZE.h"
#  include "GCHEM_FIELDS.h"
# endif
# ifdef ALLOW_CFC
#  include "CFC.h"
# endif
# ifdef ALLOW_DIC
#  include "DIC_VARS.h"
#  include "DIC_LOAD.h"
#  include "DIC_ATMOS.h"
#  include "DIC_CTRL.h"
#  include "DIC_COST.h"
# endif
# ifdef ALLOW_OBCS
#  include "OBCS_PARAMS.h"
#  include "OBCS_FIELDS.h"
#  include "OBCS_SEAICE.h"
#  ifdef ALLOW_PTRACERS
#   include "OBCS_PTRACERS.h"
#  endif
# endif
# ifdef ALLOW_SHELFICE
#  include "SHELFICE.h"
#  include "SHELFICE_COST.h"
# endif
# ifdef ALLOW_STREAMICE
#  include "STREAMICE.h"
#  include "STREAMICE_ADV.h"
#  include "STREAMICE_BDRY.h"
#  include "STREAMICE_CG.h"
# endif

#endif /* ALLOW_OPENAD or ALLOW_TAPENADE */

#ifdef ALLOW_TAPENADE
# ifdef ALLOW_EXF
#  include "EXF_PARAM.h"
#  include "EXF_INTERP_SIZE.h"
#  include "EXF_INTERP_PARAM.h"
#  include "EXF_FIELDS.h"
# endif
# ifdef ALLOW_CD_CODE
#  include "CD_CODE_VARS.h"
# endif
# ifdef ALLOW_GMREDI
#  include "GMREDI.h"
#  include "GMREDI_TAVE.h"
# endif
# ifdef ALLOW_TIMEAVE
#  include "TIMEAVE_STATV.h"
# endif
# include "FFIELDS.h"
# include "DYNVARS.h"
# include "SURFACE.h"
# include "GRID.h"
# ifdef ALLOW_SEAICE
#  include "SEAICE_SIZE.h"
#  include "SEAICE.h"
# endif
# ifdef ALLOW_KPP
#  include "KPP.h"
c#  include "KPP_PARAMS.h"
# endif
# ifdef ALLOW_DOWN_SLOPE
#  include "DWNSLP_SIZE.h"
#  include "DWNSLP_VARS.h"
# endif
# ifdef ALLOW_THSICE
#  include "THSICE_VARS.h"
# endif
# ifdef ALLOW_SALT_PLUME
#  include "SALT_PLUME.h"
# endif
# ifdef ALLOW_ECCO
#  include "ECCO_SIZE.h"
#  include "ECCO.h"
# endif
# ifdef ALLOW_BLING
#  include "BLING_VARS.h"
#  include "BLING_LOAD.h"
# endif
#endif /* ALLOW_TAPENADE */

c**************************************

C     !INPUT/OUTPUT PARAMETERS:
C     myTime :: time counter for this thread
C     myIter :: iteration counter for this thread
C     myThid :: thread number for this instance of the routine.
      INTEGER iloop
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
CEOP

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('MAIN_DO_LOOP',myThid)
#endif

#ifdef ALLOW_OPENAD
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>   LOOP   <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>  STARTS  <<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ifdef ALLOW_OPENAD_DIVA
      DO iloop = 1, nTimeSteps_l2
# else
      DO iloop = 1, nTimeSteps
# endif
#endif /* ALLOW_OPENAD */

#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('FORWARD_STEP',myThid)
#endif

#ifndef ALLOW_OPENAD
# ifdef ALLOW_ATM2D
        CALL TIMER_START('FORWARD_STEP_ATM2D  [MAIN_DO_LOOP]',myThid)
        CALL FORWARD_STEP_ATM2D( iloop, myTime, myIter, myThid )
        CALL TIMER_STOP ('FORWARD_STEP_ATM2D  [MAIN_DO_LOOP]',myThid)
# else
        CALL TIMER_START('FORWARD_STEP        [MAIN_DO_LOOP]',myThid)
        CALL FORWARD_STEP( iloop, myTime, myIter, myThid )
        CALL TIMER_STOP ('FORWARD_STEP        [MAIN_DO_LOOP]',myThid)
# endif
#else /* ALLOW_OPENAD */
# ifdef ALLOW_OPENAD_DIVA
        CALL TIMER_START('INNER_DO_LOOP       [MAIN_DO_LOOP]',myThid)
        nTimeSteps_l2 = 2
        CALL INNER_DO_LOOP( iloop, myTime, myIter, myThid )
        CALL TIMER_STOP ('INNER_DO_LOOP       [MAIN_DO_LOOP]',myThid)
# else
        CALL TIMER_START('FORWARD_STEP        [MAIN_DO_LOOP]',myThid)
        nTimeSteps_l2 = 2
        CALL FORWARD_STEP( iloop, myTime, myIter, myThid )
        CALL TIMER_STOP ('FORWARD_STEP        [MAIN_DO_LOOP]',myThid)
# endif
#endif /* ALLOW_OPENAD */

#ifdef ALLOW_OPENAD
      ENDDO
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>   LOOP   <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>  STOPS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#endif /* ALLOW_OPENAD */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('MAIN_DO_LOOP',myThid)
#endif

      RETURN
      END
