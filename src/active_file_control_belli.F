#include "IHOP_OPTIONS.h"

C     ================================================================
C     active_file_control_ihop.F: Routines that handle I/O of
C                                 tiled and global active 
C                                 variables for adjoint/tangent
C                                 linear calculations, related to 
C                                 NetCDF iHOP cost data files
C
C     Routines:
C     o active_read_ihop_tile_rl  - Read an active 1D record from file
C     o active_write_ihop_tile_rl - Write an active 1D record to file
C     o active_read_ihop_glob_rl  - Read an active 1D record from file
C     o active_write_ihop_glob_rl - Write an active 1D record to file
C     ================================================================

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ACTIVE_READ_IHOP_TILE_RL

C !INTERFACE:
      SUBROUTINE ACTIVE_READ_IHOP_TILE_RL(
     I                                  fid,
     I                                  active_num_file,
     O                                  active_var,
     I                                  lAdInit,
     I                                  irec,
     I                                  irecglob,
     I                                  theSimulationMode,
     I                                  myOptimIter,
     I                                  bi,
     I                                  bj,
     I                                  myThid )

C     !DESCRIPTION:
C     Read an active record from an iHOP .equi tiled NetCDF or 
C     Binary file

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
# include "netcdf.inc"
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER  fid
      INTEGER  active_num_file
      LOGICAL  lAdInit
      INTEGER  irec, irecglob
      INTEGER  theSimulationMode
      INTEGER  myOptimIter
      INTEGER  bi, bj, myThid
C     !OUTPUT PARAMETERS:
      _RL      active_var
CEOP

C     !LOCAL VARIABLES:
#if (defined ALLOW_IHOP) && (defined ALLOW_COST)
      INTEGER  err, varId
      INTEGER  vec_start, vec_count
      _RL      active_data_t
      REAL*8   vec_tmp(2)

      IF (ihopDoNcOutput) THEN
        vec_start = irec
        vec_count = 1

        err = NF_INQ_VARID( fid, ihop_nameequi, varId )
        CALL IHOP_COST_NF_ERROR('act_file_c 1',err,bi,bj,myThid)

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> FORWARD RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. FORWARD_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        IF (ihopDoNcOutput) THEN
          err = NF_GET_VARA_DOUBLE( fid, varId, 
     &          vec_start, vec_count, active_var )
          CALL IHOP_COST_NF_ERROR('act_file_c 2',err,bi,bj,myThid)

        ELSE
          READ( fid, REC = (irec-1)*2+1 ) vec_tmp
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( 2, vec_tmp )
#endif
          active_var = vec_tmp(1)

        ENDIF

        _END_MASTER( myThid )

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> ADJOINT RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. REVERSE_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        IF (ihopDoNcOutput) THEN
          err = NF_GET_VARA_DOUBLE( fid, varId, 
     &          vec_start, vec_count, active_data_t )
          CALL IHOP_COST_NF_ERROR('act_file_c 3',err,bi,bj,myThid)

C Add active_var from appropriate location to data.
          active_data_t = active_data_t + active_var

C Store the result on disk.
          err = NF_PUT_VARA_DOUBLE(fid, varId, 
     &          vec_start, vec_count, active_data_t )
          CALL IHOP_COST_NF_ERROR('act_file_c 4',err,bi,bj,myThid)

C Set active_var to zero.
          active_var = 0. _d 0

        ELSE
          READ( fid, REC = (irec-1)*2+1 ) vec_tmp
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( 2, vec_tmp )
#endif

C Add active_var from appropriate location to data.
          active_data_t = vec_tmp(1)
          active_data_t = active_data_t + active_var

C Store the result on disk.
          vec_tmp(1) = active_data_t
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( 2, vec_tmp )
#endif
          WRITE( fid, REC = (irec-1)*2+1 ) vec_tmp

C Set active_var to zero.
          active_var = 0. _d 0

        ENDIF

        _END_MASTER( myThid )

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> TANGENT RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. TANGENT_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        IF (ihopDoNcOutput) THEN
          err = NF_GET_VARA_DOUBLE( fid, varId,
     &          vec_start, vec_count, active_var )
          CALL IHOP_COST_NF_ERROR('act_file_c 5',err,bi,bj,myThid)

        ELSE
          READ( fid, REC = (irec-1)*2+1 ) vec_tmp
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( 2, vec_tmp )
#endif
          active_var = vec_tmp(1)

        ENDIF

        _END_MASTER( myThid )

      ENDIF

#endif /* ALLOW_IHOP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ACTIVE_WRITE_IHOP_TILE_RL

C !INTERFACE:
      SUBROUTINE ACTIVE_WRITE_IHOP_TILE_RL(
     I                                 fid,
     I                                 active_num_file,
     I                                 active_var,
     I                                 irec,
     I                                 irecglob,
     I                                 theSimulationMode,
     I                                 myOptimIter,
     I                                 bi,
     I                                 bj,
     I                                 myThid )


C     !DESCRIPTION:
C     Write an active record from an iHOP .equi tiled NetCDF or 
C     Binary file

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
# include "netcdf.inc"
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER  fid
      INTEGER  active_num_file
      _RL      active_var
      INTEGER  irec, irecglob
      INTEGER  theSimulationMode
      INTEGER  myOptimIter
      INTEGER  bi, bj, myThid
C     !OUTPUT PARAMETERS:
CEOP

C     !LOCAL VARIABLES:
#if (defined ALLOW_IHOP) && (defined ALLOW_COST)
      INTEGER  err, varId1, varId2
      INTEGER  vec_start, vec_count
      _RL      active_data_t
      REAL*8   vec_tmp(2)

      IF (ihopDoNcOutput) THEN
        vec_start = irec
        vec_count = 1

        err = NF_INQ_VARID( fid, ihop_nameequi, varId1 )
        CALL IHOP_COST_NF_ERROR('act_file_c 6',err,bi,bj,myThid)

      ENDIF


C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> FORWARD RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. FORWARD_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        IF (ihopDoNcOutput) THEN
          err = NF_PUT_VARA_DOUBLE( fid, varId1, 
     &          vec_start, vec_count, active_var )
          CALL IHOP_COST_NF_ERROR('act_file_c 7',err,bi,bj,myThid)

          err = NF_INQ_VARID( fid, 'ihopObs_ind_glob_tiled', varId2 )
          CALL IHOP_COST_NF_ERROR('act_file_c 8',err,bi,bj,myThid)
          err = NF_PUT_VAR1_INT( fid, varId2, vec_start, irecglob )
          CALL IHOP_COST_NF_ERROR('act_file_c 9',err,bi,bj,myThid)

        ELSE
          vec_tmp(1) = active_var
          vec_tmp(2) = irecglob

#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( 2, vec_tmp )
#endif
          WRITE( fid, REC = (irec-1)*2+1 ) vec_tmp

        ENDIF

        _END_MASTER( myThid )

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> ADJOINT RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. REVERSE_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        IF (ihopDoNcOutput) THEN
          err = NF_GET_VARA_DOUBLE( fid, varId1,
     &          vec_start, vec_count, active_data_t )
          CALL IHOP_COST_NF_ERROR('act_file_c10',err,bi,bj,myThid)

C Add active_var to data.
          active_var = active_var + active_data_t
          active_data_t = 0. _d 0

          err = NF_PUT_VARA_DOUBLE( fid, varId1,
     &          vec_start, vec_count, active_data_t )
          CALL IHOP_COST_NF_ERROR('act_file_c11',err,bi,bj,myThid)

        ELSE
          READ( fid, REC = (irec-1)*2+1 ) vec_tmp
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( 2, vec_tmp )
#endif

C Add active_var from appropriate location to data.
          active_var = active_var + active_data_t
          active_data_t = 0. _d 0

C Store the result on disk.
          vec_tmp(1) = active_data_t
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( 2, vec_tmp )
#endif
          WRITE( fid, REC = (irec-1)*2+1 ) vec_tmp

        ENDIF

        _END_MASTER( myThid )

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> TANGENT RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. TANGENT_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        IF (ihopDoNcOutput) THEN
          err = NF_PUT_VARA_DOUBLE( fid, varId1,
     &          vec_start, vec_count, active_var )
          CALL IHOP_COST_NF_ERROR('act_file_c12',err,bi,bj,myThid)

        ELSE
          vec_tmp(1)=active_var
          vec_tmp(2)=irecglob
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( 2, vec_tmp )
#endif  
          WRITE( fid, REC = (irec-1)*2+1 ) vec_tmp

        ENDIF

        _END_MASTER( myThid )

      ENDIF

#endif /* ALLOW_IHOP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ACTIVE_READ_IHOP_GLOB_RL

C !INTERFACE:
      SUBROUTINE ACTIVE_READ_IHOP_GLOB_RL(
     I                                  fid,
     I                                  active_num_file,
     O                                  active_var,
     I                                  lAdInit,
     I                                  irec,
     I                                  irecglob,
     I                                  theSimulationMode,
     I                                  myOptimIter,
     I                                  myThid )

C     !DESCRIPTION:
C     Read an active record from an iHOP .equi global NetCDF file

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
# include "netcdf.inc"
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER  fid
      INTEGER  active_num_file
      LOGICAL  lAdInit
      INTEGER  irec, irecglob
      INTEGER  theSimulationMode
      INTEGER  myOptimIter
      INTEGER  myThid
C     !OUTPUT PARAMETERS:
      _RL      active_var
CEOP

C     !LOCAL VARIABLES:
#if (defined ALLOW_IHOP) && (defined ALLOW_COST)
      INTEGER  err, varId
      INTEGER  vec_start, vec_count
      _RL      active_data_t
      REAL*8   vec_tmp(2)

      IF (ihopDoNcOutput) THEN
        vec_start = irecglob
        vec_count = 1

        err = NF_INQ_VARID( fid, ihop_nameequi, varId )
        CALL IHOP_COST_NF_ERROR('act_file_c13',err,1,1,myThid)

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> FORWARD RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. FORWARD_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        err = NF_GET_VARA_DOUBLE( fid, varId, 
     &        vec_start, vec_count, active_var )
        CALL IHOP_COST_NF_ERROR('act_file_c14',err,1,1,myThid)

        _END_MASTER( myThid )

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> ADJOINT RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. REVERSE_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        err = NF_GET_VARA_DOUBLE( fid, varId,
     &        vec_start, vec_count, active_data_t )
        CALL IHOP_COST_NF_ERROR('act_file_c15',err,1,1,myThid)

C Add active_var from appropriate location to data.
        active_data_t = active_data_t + active_var
C Store the result on disk.
        vec_start=irec

        err = NF_PUT_VARA_DOUBLE( fid, varId,
     &        vec_start, vec_count, active_data_t )
        CALL IHOP_COST_NF_ERROR('act_file_c16',err,1,1,myThid)

C Set active_var to zero.
        active_var = 0. _d 0

        _END_MASTER( myThid )

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> TANGENT RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. TANGENT_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        err = NF_GET_VARA_DOUBLE( fid, varId, 
     &        vec_start, vec_count, active_var )
        CALL IHOP_COST_NF_ERROR('act_file_c17',err,1,1,myThid)

        _END_MASTER( myThid )

      ENDIF

#endif /* ALLOW_IHOP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ACTIVE_WRITE_IHOP_GLOB_RL

C !INTERFACE:
      SUBROUTINE ACTIVE_WRITE_IHOP_GLOB_RL(
     I                                 fid,
     I                                 active_num_file,
     I                                 active_var,
     I                                 irec,
     I                                 irecglob,
     I                                 theSimulationMode,
     I                                 myOptimIter,
     I                                 myThid )

C     !DESCRIPTION:
C     Write an active record from an iHOP .equi global NetCDF file

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
# include "netcdf.inc"
# include "IHOP_SIZE.h"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif

C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER  fid
      INTEGER  active_num_file
      _RL      active_var
      INTEGER  irec, irecglob
      INTEGER  theSimulationMode
      INTEGER  myOptimIter
      INTEGER  myThid
C     !OUTPUT PARAMETERS:
CEOP

C     !LOCAL VARIABLES:
#if (defined ALLOW_IHOP) && (defined ALLOW_COST)
      INTEGER  err, varId
      INTEGER  vec_start, vec_count
      _RL      active_data_t
      REAL*8   vec_tmp(2)

      IF (ihopDoNcOutput) THEN
        vec_start = irecglob
        vec_count = 1

        err = NF_INQ_VARID( fid, ihop_nameequi, varId )
        CALL IHOP_COST_NF_ERROR('act_file_c18',err,1,1,myThid)

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> FORWARD RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. FORWARD_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        err = NF_PUT_VARA_DOUBLE( fid, varId, 
     &        vec_start, vec_count, active_var )
        CALL IHOP_COST_NF_ERROR('act_file_c19',err,1,1,myThid)

        _END_MASTER( myThid )

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> ADJOINT RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. REVERSE_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        vec_start = irec
        err = NF_GET_VARA_DOUBLE( fid, varId,
     &        vec_start, vec_count, active_data_t )
        CALL IHOP_COST_NF_ERROR('act_file_c20',err,1,1,myThid)

C       Add active_var to data.
        active_var = active_var + active_data_t
        active_data_t = 0. _d 0

        vec_start = irecglob
        err = NF_PUT_VARA_DOUBLE( fid, varId,
     &        vec_start, vec_count, active_data_t )
        CALL IHOP_COST_NF_ERROR('act_file_c21',err,1,1,myThid)

        _END_MASTER( myThid )

      ENDIF

C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>> TANGENT RUN <<<<<<<<<<<<<<<<<<<
C     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      IF (theSimulationMode .EQ. TANGENT_SIMULATION) THEN

        _BEGIN_MASTER( myThid )

        err = NF_PUT_VARA_DOUBLE( fid, varId, 
     &        vec_start, vec_count, active_var )
        CALL IHOP_COST_NF_ERROR('act_file_c22',err,1,1,myThid)

        _END_MASTER( myThid )

      ENDIF

#endif /* ALLOW_IHOP */

      RETURN
      END
