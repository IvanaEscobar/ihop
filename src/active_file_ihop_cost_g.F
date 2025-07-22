#include "IHOP_OPTIONS.h"

C   ==================================================================
C   active_file_cost_ihop_g.F: Routines that handle I/O of
C                              tiled and global active 
C                              variables for adjoint
C                              calculations, related to NetCDF
C                              iHOP data files
C
C   Routines:
C   o g_active_read_ihop_tile  - Read a TL active 1D record from file
C   o g_active_write_ihop_tile - Write a TL active 1D record to file
C   o g_active_read_ihop_glob  - Read a TL active 1D record from file
C   o g_active_write_ihop_glob - Write a TL active 1D record to file
C   ==================================================================

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: G_ACTIVE_READ_IHOP_TILE

C !INTERFACE:
      SUBROUTINE G_ACTIVE_READ_IHOP_TILE(
     I                              active_num_file,
     O                              active_var,
     I                              g_active_var,
     I                              irec,
     I                              lAdInit,
     I                              myOptimIter,
     I                              bi,
     I                              bj,
     I                              myThid,
     I                              dummy )

C     !DESCRIPTION:
C     Tangent linear of active_read_ihop_tile
C     Read an active record and its TL from iHOP .equi tiled NetCDF or
C     Binary files

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
C     g_active_var    :: TL active variable array
C     irec            :: record number
C     lAdInit         :: initialisation of corresponding adjoint
C                        variable and write to active file
C     myOptimIter     :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER active_num_file
      _RL     g_active_var
      INTEGER irec
      LOGICAL lAdInit
      INTEGER myOptimIter
      INTEGER bi, bj, myThid
      _RL     dummy
C     !OUTPUT PARAMETERS:
C     active_var      :: active variable array
      _RL     active_var
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_READ_IHOP_TILE_RL(
     &     ncidFWD(active_num_file,bi,bj),
     &     active_num_file,
     &     active_var, lAdInit, irec,
     &     ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     FORWARD_SIMULATION, myOptimIter, bi, bj, myThid )

      CALL ACTIVE_READ_IHOP_TILE_RL(
     &     ncidTL(active_num_file,bi,bj),
     &     active_num_file,
     &     g_active_var, lAdInit, irec,
     &     ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     TANGENT_SIMULATION, myOptimIter, bi, bj, myThid )
#endif

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: G_ACTIVE_WRITE_IHOP_TILE

C !INTERFACE:
      SUBROUTINE G_ACTIVE_WRITE_IHOP_TILE(
     I                              active_num_file,
     I                              active_var,
     I                              g_active_var,
     I                              irec,
     I                              myOptimIter,
     I                              bi,
     I                              bj,
     I                              myThid,
     I                              dummy,
     O                              g_dummy )

C     !DESCRIPTION:
C     Tangent linear of active_write_ihop_tile
C     Write an active record and its TL from iHOP .equi tiled NetCDF
C     or Binary files

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
C     g_active_var    :: TL active variable array
C     irec            :: record number
C     myOptimIter     :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER active_num_file
      _RL     active_var
      _RL     g_active_var
      INTEGER irec
      INTEGER myOptimIter
      INTEGER bi, bj, myThid
      _RL     dummy
C     !OUTPUT PARAMETERS:
      _RL     g_dummy
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_WRITE_IHOP_TILE_RL(
     &     ncidFWD(active_num_file,bi,bj),
     &     active_num_file,
     &     active_var, irec,
     &     ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     FORWARD_SIMULATION, myOptimIter, bi, bj, myThid )

      CALL ACTIVE_WRITE_IHOP_TILE_RL(
     &     ncidTL(active_num_file,bi,bj),
     &     active_num_file,
     &     g_active_var, irec, 
     &     ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     TANGENT_SIMULATION, myOptimIter, bi, bj, myThid )

#endif

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: G_ACTIVE_READ_IHOP_GLOB

C !INTERFACE:
      SUBROUTINE G_ACTIVE_READ_IHOP_GLOB(
     I                              active_num_file,
     O                              active_var,
     I                              g_active_var,
     I                              irec,
     I                              lAdInit,
     I                              myOptimIter,
     I                              myThid,
     I                              dummy )

C     !DESCRIPTION:
C     Tangent linear of active_read_ihop_glob
C     Read an active record and its TL from iHOP .equi global NetCDF

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
C     g_active_var    :: TL active variable array
C     irec            :: record number
C     lAdInit         :: initialisation of corresponding adjoint
C                        variable and write to active file
C     myOptimIter     :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER active_num_file
      _RL     g_active_var
      INTEGER irec
      INTEGER myThid
      LOGICAL lAdInit
      INTEGER myOptimIter
      _RL     dummy
C     !OUTPUT PARAMETERS:
      _RL     active_var
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_READ_IHOP_GLOB_RL(
     &     ncidGLOB(active_num_file),
     &     active_num_file,
     &     active_var, lAdInit, irec, 
     &     ihopObs_ind_glob(active_num_file,irec),
     &     FORWARD_SIMULATION, myOptimIter, myThid )

      CALL ACTIVE_READ_IHOP_GLOB_RL(
     &     ncidTLGLOB(active_num_file),
     &     active_num_file,
     &     g_active_var, lAdInit, irec, 
     &     ihopObs_ind_glob(active_num_file,irec),
     &     TANGENT_SIMULATION, myOptimIter, myThid )

#endif

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: G_ACTIVE_WRITE_IHOP_GLOB

C !INTERFACE:
      SUBROUTINE G_ACTIVE_WRITE_IHOP_GLOB(
     I                           active_num_file,
     I                           active_var,
     I                           g_active_var,
     I                           irec,
     I                           myOptimIter,
     I                           myThid,
     I                           dummy,
     O                           g_dummy )

C     !DESCRIPTION:
C     Tangent linear of active_write_ihop_glob
C     Write an active record and its TL from iHOP .equi global NetCDF

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
C     g_active_var    :: TL active variable array
C     irec            :: record number
C     myOptimIter     :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER active_num_file
      _RL     active_var
      _RL     g_active_var
      INTEGER irec
      INTEGER myOptimIter
      INTEGER myThid
      _RL     dummy
C     !OUTPUT PARAMETERS:
      _RL     g_dummy
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_WRITE_IHOP_GLOB_RL(
     &     ncidGLOB(active_num_file),
     &     active_num_file,
     &     active_var, irec, 
     &     ihopObs_ind_glob(active_num_file,irec),
     &     FORWARD_SIMULATION, myOptimIter, myThid )

      CALL ACTIVE_WRITE_IHOP_GLOB_RL(
     &     ncidTLGLOB(active_num_file),
     &     active_num_file,
     &     g_active_var, irec, 
     &     ihopObs_ind_glob(active_num_file,irec),
     &     TANGENT_SIMULATION, myOptimIter, myThid )

#endif

      RETURN
      END