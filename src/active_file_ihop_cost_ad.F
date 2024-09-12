#include "IHOP_OPTIONS.h"

C--  File active_file_ihop_cost_ad.F:
C--   Contents
C--   o adactive_read_ihop_tile
C--   o adactive_write_ihop_tile
C--   o adactive_read_ihop_glob
C--   o adactive_write_ihop_glob

CBOP
C     !ROUTINE: ADACTIVE_READ_IHOP_TILE
C     !INTERFACE:
      SUBROUTINE ADACTIVE_READ_IHOP_TILE(
     I                              active_num_file,
     I                              adactive_var,
     I                              irec,
     I                              lAdInit,
     I                              myIter,
     I                              bi,
     I                              bj,
     I                              myThid )

C     !DESCRIPTION:
C     ==================================================================
C     | Adjoint of active_read_ihop_tile
C     ==================================================================

C !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C !INPUT PARAMETERS:
C     active_num_file: file number
C     adactive_var:    array
C     irec:            record number
C     myIter:          number of optimization iteration (default: 0)
C     myThid:          thread number for this instance
C     lAdInit:         initialisation of corresponding adjoint
C                      variable and write to active file
      INTEGER active_num_file
      _RL     adactive_var
      INTEGER irec
      INTEGER myIter
      INTEGER bi, bj, myThid
      LOGICAL lAdInit
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL active_read_ihop_tile_rl( fidadj_obs(active_num_file,bi,bj),
     &     active_num_file,
     &     adactive_var, lAdInit,
     &     irec, ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     REVERSE_SIMULATION, myIter, bi, bj, myThid)

#endif

      END

C     ==================================================================

CBOP
C     !ROUTINE: ADACTIVE_WRITE_IHOP_TILE
C     !INTERFACE:
      SUBROUTINE ADACTIVE_WRITE_IHOP_TILE(
     I                                  active_num_file,
     I                                  adactive_var,
     I                                  irec,
     I                                  myIter,
     I                                  bi,
     I                                  bj,
     I                                  myThid,
     I                                  dummy )

C     !DESCRIPTION:
C     ==================================================================
C     | Adjoint of active_write_ihop_tile
C     ==================================================================

C !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C !INPUT PARAMETERS:
C     active_num_file: file number
C     adactive_var:    array
C     irec:            record number
C     myIter:          number of optimization iteration (default: 0)
C     myThid:          thread number for this instance
C     lAdInit:         initialisation of corresponding adjoint
C                      variable and write to active file
      INTEGER active_num_file
      _RL     adactive_var
      INTEGER irec
      INTEGER myIter
      INTEGER bi, bj, myThid
      _RL     dummy
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL active_write_ihop_tile_rl( fidadj_obs(active_num_file,bi,bj),
     &     active_num_file,
     &     adactive_var, 
     &     irec, ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     REVERSE_SIMULATION, myIter, bi, bj, myThid)

#endif

      END

C     ==================================================================


CBOP
C     !ROUTINE: ADACTIVE_READ_IHOP_GLOB
C     !INTERFACE:
      SUBROUTINE ADACTIVE_READ_IHOP_GLOB(
     I                              active_num_file,
     I                              adactive_var,
     I                              irec,
     I                              lAdInit,
     I                              myIter,
     I                              myThid )

C     !DESCRIPTION:
C     ==================================================================
C     | SUBROUTINE ADACTIVE_READ_IHOP_GLOB
C     | Adjoint of ACTIVE_READ_IHOP_GLOB.F
C     ==================================================================

C !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C !INPUT PARAMETERS:
C     active_num_file: file number
C     adactive_var:    array
C     irec:            record number
C     myIter:          number of optimization iteration (default: 0)
C     myThid:          thread number for this instance
C     lAdInit:         initialisation of corresponding adjoint
C                      variable and write to active file
      INTEGER active_num_file
      _RL     adactive_var
      INTEGER irec
      INTEGER myIter
      INTEGER myThid
      LOGICAL lAdInit
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL active_read_ihop_glob_rl( fidadglobal(active_num_file),
     &     active_num_file,
     &     adactive_var, lAdInit,
     &     irec, ihopObs_ind_glob(active_num_file,irec),
     &     REVERSE_SIMULATION, myIter, myThid)

#endif

      END

C     ==================================================================

CBOP
C     !ROUTINE: ADACTIVE_WRITE_IHOP_GLOB
C     !INTERFACE:
      SUBROUTINE ADACTIVE_WRITE_IHOP_GLOB(
     I                                  active_num_file,
     I                                  adactive_var,
     I                                  irec,
     I                                  myIter,
     I                                  myThid,
     I                                  dummy )

C     !DESCRIPTION:
C     ==================================================================
C     | SUBROUTINE ADACTIVE_WRITE_IHOP_GLOB
C     | Adjoint of ACTIVE_WRITE_IHOP_GLOB.F
C     ==================================================================

C !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C !INPUT PARAMETERS:
C     active_num_file: file number
C     adactive_var:    array
C     irec:            record number
C     myIter:          number of optimization iteration (default: 0)
C     myThid:          thread number for this instance
C     lAdInit:         initialisation of corresponding adjoint
C                      variable and write to active file
      INTEGER active_num_file
      _RL     adactive_var
      INTEGER irec
      INTEGER myIter
      INTEGER myThid
      _RL     dummy
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL active_write_ihop_glob_rl( fidadglobal(active_num_file),
     &     active_num_file,
     &     adactive_var,
     &     irec, ihopObs_ind_glob(active_num_file,irec),
     &     REVERSE_SIMULATION, myIter,myThid)

#endif

      END

C     ==================================================================

