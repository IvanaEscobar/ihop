#include "IHOP_OPTIONS.h"

C--  File active_file_ihop_cost_ad.F:
C--   Contents
C--   o ftlactive_read_ihop_tile
C--   o ftlactive_write_ihop_tile
C--   o ftlactive_read_ihop_glob
C--   o ftlactive_write_ihop_glob

CBOP
C     !ROUTINE: FTLACTIVE_READ_IHOP_TILE
C     !INTERFACE:
      SUBROUTINE FTLACTIVE_READ_IHOP_TILE(
     I                              active_num_file,
     O                              active_var,
     I                              irec,
     I                              lAdInit,
     I                              myOptimIter,
     I                              bi,
     I                              bj,
     I                              myThid,
     I                              dummy,
     I                              ftlactive_var )

C     !DESCRIPTION:
C     ==================================================================
C     | SUBROUTINE FTLACTIVE_READ_IHOP_TILE
C     | Read an active record from an ihop .equi. tiled file
C     | for tangent linear calculations
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
C     active_var:      array
C     irec:            record number
C     myOptimIter:     number of optimization iteration (default: 0)
C     myThid:          thread number for this instance
C     lAdInit:         initialisation of corresponding adjoint
C                      variable and write to active file
      _RL     active_var
      _RL     ftlactive_var
      INTEGER irec, active_num_file
      INTEGER myOptimIter
      INTEGER bi, bj, myThid
      LOGICAL lAdInit
      _RL     dummy
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL active_read_ihop_tile_rl( fidfwd_obs(active_num_file,bi,bj),
     &     active_num_file,
     &     active_var, lAdInit,
     &     irec, ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     TANGENT_SIMULATION, myOptimIter, bi, bj, myThid)

      CALL active_read_ihop_tile_rl( fidftl_obs(active_num_file,bi,bj),
     &     active_num_file,
     &     ftlactive_var, lAdInit,
     &     irec, ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     TANGENT_SIMULATION, myOptimIter, bi, bj, myThid)
#endif

      RETURN
      END

C     ==================================================================

CBOP
C     !ROUTINE: FTLACTIVE_WRITE_IHOP_TILE
C     !INTERFACE:
      SUBROUTINE FTLACTIVE_WRITE_IHOP_TILE(
     I                              active_num_file,
     I                              active_var,
     I                              irec,
     I                              myOptimIter,
     I                              bi,
     I                              bj,
     I                              myThid,
     I                              dummy,
     I                              ftlactive_var,
     I                              ftldummy )

C     !DESCRIPTION:
C     ==================================================================
C     | SUBROUTINE FTLACTIVE_READ_IHOP_TILE
C     | Read an active record from an ihop .equi. tiled file
C     | for tangent linear calculations
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
C     active_var:      array
C     irec:            record number
C     myOptimIter:     number of optimization iteration (default: 0)
C     myThid:          thread number for this instance
      INTEGER active_num_file
      _RL     active_var
      _RL     ftlactive_var
      INTEGER irec
      INTEGER myOptimIter
      INTEGER bi, bj, myThid
      _RL     dummy
      _RL     ftldummy
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL active_write_ihop_tile_rl( fidfwd_obs(active_num_file,bi,bj),
     &     active_num_file,
     &     active_var,
     &     irec, ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     TANGENT_SIMULATION, myOptimIter, bi, bj, myThid)

      CALL active_write_ihop_tile_rl( fidftl_obs(active_num_file,bi,bj),
     &     active_num_file,
     &     ftlactive_var,
     &     irec, ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     TANGENT_SIMULATION, myOptimIter, bi, bj, myThid)

#endif

      RETURN
      END

C     ==================================================================


CBOP
C     !ROUTINE: FTLACTIVE_READ_IHOP_GLOB
C     !INTERFACE:
      SUBROUTINE FTLACTIVE_READ_IHOP_GLOB(
     I                              active_num_file,
     O                              active_var,
     I                              irec,
     I                              lAdInit,
     I                              myOptimIter,
     I                              myThid,
     I                              dummy,
     I                              ftlactive_var )

C     !DESCRIPTION:
C     ==================================================================
C     | SUBROUTINE FTLACTIVE_READ_IHOP_GLOB
C     | Read an active record from an ihop .equi. tiled file
C     | for tangent linear calculations
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
C     active_var:      array
C     irec:            record number
C     myOptimIter:     number of optimization iteration (default: 0)
C     myThid:          thread number for this instance
C     lAdInit:         initialisation of corresponding adjoint
C                      variable and write to active file
      _RL     active_var
      _RL     ftlactive_var
      INTEGER irec, active_num_file
      INTEGER myOptimIter
      INTEGER myThid
      LOGICAL lAdInit
      _RL     dummy
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL active_read_ihop_glob_rl( fidglobal(active_num_file),
     &     active_num_file,
     &     active_var, lAdInit,
     &     irec, ihopObs_ind_glob(active_num_file,irec),
     &     TANGENT_SIMULATION, myOptimIter, myThid)

      CALL active_read_ihop_glob_rl( fidftlglobal(active_num_file),
     &     active_num_file,
     &     ftlactive_var, lAdInit,
     &     irec, ihopObs_ind_glob(active_num_file,irec),
     &     TANGENT_SIMULATION, myOptimIter, myThid)

#endif

      RETURN
      END

C     ==================================================================

CBOP
C     !ROUTINE: FTLACTIVE_WRITE_IHOP_GLOB
C     !INTERFACE:
      SUBROUTINE FTLACTIVE_WRITE_IHOP_GLOB(
     I                           active_num_file,
     I                           active_var,
     I                           irec,
     I                           myOptimIter,
     I                           myThid,
     I                           dummy,
     I                           ftlactive_var,
     I                           ftldummy )

C     !DESCRIPTION:
C     ==================================================================
C     | SUBROUTINE FTLACTIVE_READ_IHOP_GLOB
C     | Read an active record from an ihop .equi. tiled file
C     | for tangent linear calculations
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
C     active_var:      array
C     irec:            record number
C     myOptimIter:     number of optimization iteration (default: 0)
C     myThid:          thread number for this instance
      integer active_num_file
      _RL     active_var
      _RL     ftlactive_var
      integer irec
      integer myOptimIter
      integer myThid
      _RL     dummy
      _RL     ftldummy
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL active_write_ihop_glob_rl( fidglobal(active_num_file),
     &     active_num_file,
     &     active_var,
     &     irec, ihopObs_ind_glob(active_num_file,irec),
     &     TANGENT_SIMULATION, myOptimIter,myThid)

      CALL active_write_ihop_glob_rl( fidftlglobal(active_num_file),
     &     active_num_file,
     &     ftlactive_var,
     &     irec, ihopObs_ind_glob(active_num_file,irec),
     &     TANGENT_SIMULATION, myOptimIter,myThid)

#endif

      RETURN
      END

C     ==================================================================
