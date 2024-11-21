#include "IHOP_OPTIONS.h"

C   ==================================================================
C   active_file_cost_ihop_ad.F: Routines that handle I/O of
C                               tiled and global active 
C                               variables for adjoint
C                               calculations, related to NetCDF
C                               iHOP data files
C
C   Routines:
C   o adactive_read_ihop_tile  - Read an AD active 1D record from file
C   o adactive_write_ihop_tile - Write an AD active 1D record to file
C   o adactive_read_ihop_glob  - Read an AD active 1D record from file
C   o adactive_write_ihop_glob - Write an AD active 1D record to file
C   ==================================================================

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ADACTIVE_READ_IHOP_TILE

C !INTERFACE:
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
C     Adjoint of active_read_ihop_tile
C     Read an adjoint active record from an iHOP .equi tiled NetCDF or
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
C     adactive_var    :: active variable array
C     irec            :: record number
C     lAdInit         :: initialisation of corresponding adjoint
C                        variable and write to active file
C     myIter          :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER active_num_file
      _RL     adactive_var
      INTEGER irec
      LOGICAL lAdInit
      INTEGER myIter
      INTEGER bi, bj, myThid
C     !OUTPUT PARAMETERS:
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_READ_IHOP_TILE_RL(
     &     ncidAD(active_num_file,bi,bj),
     &     active_num_file,
     &     adactive_var, lAdInit, irec,
     &     ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     REVERSE_SIMULATION, myIter, bi, bj, myThid )

#endif

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ADACTIVE_WRITE_IHOP_TILE

C !INTERFACE:
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
C     Adjoint of active_write_ihop_tile
C     Write an adjoint active record from an iHOP .equi tiled NetCDF
C     or Binary file

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
C     adactive_var    :: active variable array
C     irec            :: record number
C     myIter          :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER active_num_file
      _RL     adactive_var
      INTEGER irec
      INTEGER myIter
      INTEGER bi, bj, myThid
      _RL     dummy
C     !OUTPUT PARAMETERS:
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_WRITE_IHOP_TILE_RL(
     &     ncidAD(active_num_file,bi,bj),
     &     active_num_file,
     &     adactive_var, irec,
     &     ihopObs_ind_glob_tiled(active_num_file,irec,bi,bj),
     &     REVERSE_SIMULATION, myIter, bi, bj, myThid )

#endif

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ADACTIVE_READ_IHOP_GLOB

C !INTERFACE:
      SUBROUTINE ADACTIVE_READ_IHOP_GLOB(
     I                              active_num_file,
     I                              adactive_var,
     I                              irec,
     I                              lAdInit,
     I                              myIter,
     I                              myThid )

C     !DESCRIPTION:
C     Adjoint of active_read_ihop_glob
C     Read an adjoint active record from an iHOP .equi global NetCDF

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
C     adactive_var    :: active variable array
C     irec            :: record number
C     lAdInit         :: initialisation of corresponding adjoint
C                        variable and write to active file
C     myIter          :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER active_num_file
      _RL     adactive_var
      INTEGER irec
      LOGICAL lAdInit
      INTEGER myIter
      INTEGER myThid
C     !OUTPUT PARAMETERS:
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_READ_IHOP_GLOB_RL(
     &     ncidADGLOB(active_num_file),
     &     active_num_file,
     &     adactive_var, lAdInit, irec,
     &     ihopObs_ind_glob(active_num_file,irec),
     &     REVERSE_SIMULATION, myIter, myThid )

#endif

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ADACTIVE_WRITE_IHOP_GLOB

C !INTERFACE:
      SUBROUTINE ADACTIVE_WRITE_IHOP_GLOB(
     I                                  active_num_file,
     I                                  adactive_var,
     I                                  irec,
     I                                  myIter,
     I                                  myThid,
     I                                  dummy )

C     !DESCRIPTION:
C     Adjoint of active_write_ihop_glob
C     Write an adjoint active record from an iHOP .equi global NetCDF

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
C     adactive_var    :: active variable array
C     irec            :: record number
C     myIter          :: optimization iteration (default: 0)
C     myThid          :: my thread ID
      INTEGER active_num_file
      _RL     adactive_var
      INTEGER irec
      INTEGER myIter
      INTEGER myThid
      _RL     dummy
C     !OUTPUT PARAMETERS:
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)

      CALL ACTIVE_WRITE_IHOP_GLOB_RL(
     &     ncidADGLOB(active_num_file),
     &     active_num_file,
     &     adactive_var, irec,
     &     ihopObs_ind_glob(active_num_file,irec),
     &     REVERSE_SIMULATION, myIter, myThid )

#endif

      RETURN
      END
