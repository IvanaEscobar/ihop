#include "IHOP_OPTIONS.h"

C   ==================================================================
C   active_file_ihop.F: Routines that handle I/O of
C                          tiled and global active 
C                          variables for adjoint/tangent
C                          linear calculations, related to 
C                          NetCDF iHOP cost data files
C
C   Routines:
C   o active_read_ihop_tile  - Read an active 1D record from file
C   o active_write_ihop_tile - Write an active 1D record to file
C   o active_read_ihop_glob  - Read an active 1D record from file
C   o active_write_ihop_glob - Write an active 1D record to file
C   ==================================================================

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ACTIVE_READ_IHOP_TILE

C !INTERFACE:
      SUBROUTINE ACTIVE_READ_IHOP_TILE(
     I                               active_num_file,
     O                               active_var,
     I                               irec,
     I                               lAdInit,
     I                               myOptimIter,
     I                               bi,
     I                               bj,
     I                               myThid,
     I                               dummy )

C     !DESCRIPTION:
C     Read an active record from an iHOP .equi tiled NetCDF or 
C     Binary file

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C     !INPUT PARAMETERS:
C     active_num_file :: file ID 
C     irec            :: record number
C     lAdInit         :: initialisation of corresponding adjoint
C                        variable and write to active file
C     myOptimIter     :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER  active_num_file
      INTEGER  irec
      LOGICAL  lAdInit
      INTEGER  myOptimIter
      INTEGER  bi, bj, myThid
      _RL      dummy
C     !OUTPUT PARAMETERS:
C     active_var      :: active variable array
      _RL      active_var
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_READ_IHOP_TILE_RL(
     &     ncidFWD(active_num_file,bi,bj),
     &     active_num_file,
     &     active_var, lAdInit, irec,
     &     ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     FORWARD_SIMULATION, myOptimIter, bi, bj, myThid )

#endif

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ACTIVE_WRITE_IHOP_TILE

C !INTERFACE:
      SUBROUTINE ACTIVE_WRITE_IHOP_TILE(
     I                           active_num_file,
     I                           active_var,
     I                           irec,
     I                           myOptimIter,
     I                           bi,
     I                           bj,
     I                           myThid,
     O                           dummy )

C     !DESCRIPTION:
C     Write an active record from an iHOP .equi tiled NetCDF or 
C     Binary file

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C     !INPUT PARAMETERS:
C     active_num_file :: file ID 
C     active_var      :: active variable array
C     irec            :: record number
C     myOptimIter     :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER  active_num_file
      _RL      active_var
      INTEGER  irec
      INTEGER  myOptimIter
      INTEGER  bi, bj, myThid
C     !OUTPUT PARAMETERS:
      _RL      dummy
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_WRITE_IHOP_TILE_RL(
     &     ncidFWD(active_num_file,bi,bj),
     &     active_num_file,
     &     active_var, irec,
     &     ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     FORWARD_SIMULATION, myOptimIter, bi, bj, myThid )

#endif

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ACTIVE_READ_IHOP_GLOB

C !INTERFACE:
      SUBROUTINE ACTIVE_READ_IHOP_GLOB(
     I                               active_num_file,
     O                               active_var,
     I                               irec,
     I                               lAdInit,
     I                               myOptimIter,
     I                               myThid,
     I                               dummy )

C     !DESCRIPTION:
C     Read an active record from an iHOP .equi global NetCDF

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C     !INPUT PARAMETERS:
C     active_num_file :: file ID 
C     irec            :: record number
C     lAdInit         :: initialisation of corresponding adjoint
C                        variable and write to active file
C     myOptimIter     :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER  active_num_file
      INTEGER  irec
      LOGICAL  lAdInit
      INTEGER  myOptimIter
      INTEGER  bi, bj, myThid
      _RL      dummy
C     !OUTPUT PARAMETERS:
C     active_var      :: active variable array
      _RL      active_var
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_READ_IHOP_GLOB_RL(
     &     ncidGLOB(active_num_file),
     &     active_num_file,
     &     active_var, lAdInit, irec,
     &     ihopObs_ind_glob(active_num_file,irec),
     &     FORWARD_SIMULATION, myOptimIter, myThid )

#endif

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ACTIVE_WRITE_IHOP_GLOB

C !INTERFACE:
      SUBROUTINE ACTIVE_WRITE_IHOP_GLOB(
     I                           active_num_file,
     I                           active_var,
     I                           irec,
     I                           myOptimIter,
     I                           myThid,
     O                           dummy )

C     !DESCRIPTION:
C     Write an active record from an iHOP .equi global NetCDF

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C     !INPUT PARAMETERS:
C     active_num_file :: file ID 
C     active_var      :: active variable array
C     irec            :: record number
C     myOptimIter     :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER  active_num_file
      _RL      active_var
      INTEGER  irec
      INTEGER  myOptimIter
      INTEGER  bi, bj, myThid
C     !OUTPUT PARAMETERS:
      _RL      dummy
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_WRITE_IHOP_GLOB_RL(
     &     ncidGLOB(active_num_file),
     &     active_num_file,
     &     active_var, irec,
     &     ihopObs_ind_glob(active_num_file,irec),
     &     FORWARD_SIMULATION, myOptimIter, myThid )

#endif

      RETURN
      END
