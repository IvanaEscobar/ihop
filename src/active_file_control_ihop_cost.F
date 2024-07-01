#include "IHOP_OPTIONS.h"

C--  File active_file_control_ihop_tile_rl.F:
C--   Contents
C--   o active_read_ihop_tile_rl
C--   o active_write_ihop_tile_rl

CBOP 0
C     !ROUTINE: ACTIVE_READ_ihop_TILE_RL
C     !INTERFACE:
      SUBROUTINE ACTIVE_READ_ihop_TILE_RL(
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
C     ==================================================================
C     | SUBROUTINE ACTIVE_READ_ihop_TILE_RL
C     | Read an active record from an IHOP .equi. tiled file
C     | (can be netcdf or binary)  
C     ==================================================================

C !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
# include "netcdf.inc"
# include "IHOP_SIZE.h"
# include "IHOP_COST.h"
#endif

C !INPUT PARAMETERS:
      integer  fid
      integer  active_num_file
      integer  irec, irecglob
      integer  theSimulationMode
      integer  myOptimIter
      integer  bi, bj, myThid
      _RL      active_var
      logical  lAdInit

CEOP

#ifdef ALLOW_IHOP
C !LOCAL VARIABLES:
      integer  err, varid1, varid2
      integer  vec_start, vec_count
      integer  ivar, i
      _RL      active_data_t
      real*8   vec_tmp(2)

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> FORWARD RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. FORWARD_SIMULATION) then

        _BEGIN_MASTER( myThid )

      if (ihopDoNcOutput) then
      
      vec_start=irec
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_GET_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_var)


      else

      read(fid,rec=(irec-1)*2+1 )
     & vec_tmp
#ifdef _BYTESWAPIO
            call MDS_BYTESWAPR8(2,vec_tmp)
#endif
      active_var=vec_tmp(1)
      read(fid,rec=(irec-1)*2+2 )
     & vec_tmp
#ifdef _BYTESWAPIO
            call MDS_BYTESWAPR8(2,vec_tmp)
#endif

      endif

        _END_MASTER( myThid )

      endif

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> ADJOINT RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. REVERSE_SIMULATION) then

        _BEGIN_MASTER( myThid )

      if (ihopDoNcOutput) then

      vec_start=irec
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_GET_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_data_t)

c     Add active_var from appropriate location to data.
      active_data_t = active_data_t + active_var
c     Store the result on disk.
      vec_start=irec
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_PUT_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_data_t)

c     Set active_var to zero.
      active_var = 0. _d 0

      else

      read(fid,rec=(irec-1)*2+1 )
     & vec_tmp
#ifdef _BYTESWAPIO
            call MDS_BYTESWAPR8(2,vec_tmp)
#endif
      active_data_t=vec_tmp(1)

c     Add active_var from appropriate location to data.
      active_data_t = active_data_t + active_var

c     Store the result on disk.
      vec_tmp(1)=active_data_t
#ifdef _BYTESWAPIO
            call MDS_BYTESWAPR8(2,vec_tmp)
#endif
      write(fid,rec=(irec-1)*2+1 )
     & vec_tmp

c     Set active_var to zero.
      active_var = 0. _d 0

      endif

        _END_MASTER( myThid )

      endif

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> TANGENT RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. TANGENT_SIMULATION) then

        _BEGIN_MASTER( myThid )

      if (ihopDoNcOutput) then

      vec_start=irec
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_GET_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_var)

      else

      read(fid,rec=(irec-1)*2+1 )
     & vec_tmp
#ifdef _BYTESWAPIO
            call MDS_BYTESWAPR8(2,vec_tmp)
#endif
      active_var=vec_tmp(1)

      endif

        _END_MASTER( myThid )

      endif

#endif /* ALLOW_IHOP */

      RETURN
      END

C     =================================================================

CBOP 0
C     !ROUTINE: ACTIVE_WRITE_ihop_TILE_RL
C     !INTERFACE:
      SUBROUTINE ACTIVE_WRITE_ihop_TILE_RL(
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
C     ==========================================================
C     | SUBROUTINE ACTIVE_WRITE_ihop_TILE_RL
C     | Write an active record to an ihop .equi. tiled file
C     | (can be netcdf or binary)  
C     ==========================================================

C !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
# include "netcdf.inc"
# include "IHOP_SIZE.h"
# include "IHOP_COST.h"
#endif

C !INPUT PARAMETERS:
      integer  fid
      integer  active_num_file
      integer  irec, irecglob
      integer  theSimulationMode
      integer  myOptimIter
      integer  bi, bj, myThid
      _RL      active_var
CEOP

#ifdef ALLOW_IHOP
C !LOCAL VARIABLES:
      integer  err, varid1, varid2, varid3
      integer  vec_start, vec_count
      integer  ivar, i
      _RL      active_data_t
      real*8   vec_tmp(2)

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> FORWARD RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. FORWARD_SIMULATION) then

        _BEGIN_MASTER( myThid )

      if (ihopDoNcOutput) then

      vec_start=irec
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_PUT_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_var)


      err = NF_INQ_VARID(fid,'ihopObs_ind_glob_tiled', varid3 )
      err = NF_PUT_VAR1_INT(fid, varid3, vec_start,
     & irecglob)

      else

      vec_tmp(1)=active_var
      vec_tmp(2)=irecglob
#ifdef _BYTESWAPIO
            call MDS_BYTESWAPR8(2,vec_tmp)
#endif
      write(fid,rec=(irec-1)*2+1 )
     & vec_tmp
      vec_tmp(2)=irecglob
#ifdef _BYTESWAPIO
            call MDS_BYTESWAPR8(2,vec_tmp)
#endif
      write(fid,rec=(irec-1)*2+2 )
     & vec_tmp

      endif

        _END_MASTER( myThid )

      endif

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> ADJOINT RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. REVERSE_SIMULATION) then

        _BEGIN_MASTER( myThid )

      if (ihopDoNcOutput) then

      vec_start=irec
      vec_count=1
      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_GET_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_data_t)

c     Add active_var to data.
      active_var = active_var + active_data_t
      active_data_t = 0. _d 0

      vec_start=irec
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_PUT_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_data_t)

      else

      read(fid,rec=(irec-1)*2+1 )
     & vec_tmp
#ifdef _BYTESWAPIO
            call MDS_BYTESWAPR8(2,vec_tmp)
#endif
      active_data_t=vec_tmp(1)

c     Add active_var from appropriate location to data.
      active_var = active_var + active_data_t
      active_data_t = 0. _d 0

c     Store the result on disk.
      vec_tmp(1)=active_data_t
#ifdef _BYTESWAPIO
            call MDS_BYTESWAPR8(2,vec_tmp)
#endif
      write(fid,rec=(irec-1)*2+1 )
     & vec_tmp

      endif

        _END_MASTER( myThid )

      endif

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> TANGENT RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. TANGENT_SIMULATION) then

        _BEGIN_MASTER( myThid )

      if (ihopDoNcOutput) then

      vec_start=irec
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_PUT_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_var)

      else

      vec_tmp(1)=active_var
      vec_tmp(2)=irecglob
#ifdef _BYTESWAPIO
            call MDS_BYTESWAPR8(2,vec_tmp)
#endif
      write(fid,rec=(irec-1)*2+1 )
     & vec_tmp

      endif

        _END_MASTER( myThid )

      endif

#endif /* ALLOW_IHOP */

      RETURN
      END

C     ==================================================================


CBOP 0
C     !ROUTINE: ACTIVE_READ_ihop_GLOB_RL
C     !INTERFACE:
      SUBROUTINE ACTIVE_READ_ihop_GLOB_RL(
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
C     ==================================================================
C     | SUBROUTINE ACTIVE_READ_ihop_GLOB_RL
C     | Read an active record from an IHOP .equi. global file
C     ==================================================================

C !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
# include "netcdf.inc"
# include "IHOP_SIZE.h"
# include "IHOP_COST.h"
#endif

C !INPUT PARAMETERS:
      integer  fid
      integer  active_num_file
      integer  irec, irecglob
      integer  theSimulationMode
      integer  myOptimIter
      integer  myThid
      _RL      active_var
      logical  lAdInit

CEOP

#ifdef ALLOW_IHOP
C !LOCAL VARIABLES:
      integer  err, varid1, varid2
      integer  vec_start, vec_count
      integer  ivar, i
      _RL      active_data_t
      real*8   vec_tmp(2)

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> FORWARD RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. FORWARD_SIMULATION) then

        _BEGIN_MASTER( myThid )
      
      vec_start=irecglob
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_GET_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_var)

        _END_MASTER( myThid )

      endif

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> ADJOINT RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. REVERSE_SIMULATION) then

        _BEGIN_MASTER( myThid )

      vec_start=irecglob
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_GET_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_data_t)

c     Add active_var from appropriate location to data.
      active_data_t = active_data_t + active_var
c     Store the result on disk.
      vec_start=irec
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_PUT_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_data_t)

c     Set active_var to zero.
      active_var = 0. _d 0

        _END_MASTER( myThid )

      endif

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> TANGENT RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. TANGENT_SIMULATION) then

        _BEGIN_MASTER( myThid )

      vec_start=irecglob
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_GET_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_var)

        _END_MASTER( myThid )

      endif

#endif /* ALLOW_IHOP */

      RETURN
      END

C     =================================================================

CBOP 0
C     !ROUTINE: ACTIVE_WRITE_ihop_GLOB_RL
C     !INTERFACE:
      SUBROUTINE ACTIVE_WRITE_ihop_GLOB_RL(
     I                                 fid,
     I                                 active_num_file,
     I                                 active_var,
     I                                 irec,
     I                                 irecglob,
     I                                 theSimulationMode,
     I                                 myOptimIter,
     I                                 myThid )

C     !DESCRIPTION: 
C     ==========================================================
C     | SUBROUTINE ACTIVE_WRITE_ihop_GLOB_RL
C     | Write an active record to an ihop .equi. global file
C     ==========================================================

C !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
# include "netcdf.inc"
# include "IHOP_SIZE.h"
# include "IHOP_COST.h"
#endif

C !INPUT PARAMETERS:
      integer  fid
      integer  active_num_file
      integer  irec, irecglob
      integer  theSimulationMode
      integer  myOptimIter
      integer  myThid
      _RL      active_var
CEOP

#ifdef ALLOW_IHOP
C !LOCAL VARIABLES:
      integer  err, varid1, varid2
      integer  vec_start, vec_count
      integer  ivar, i
      _RL      active_data_t
      real*8   vec_tmp(2)

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> FORWARD RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. FORWARD_SIMULATION) then

        _BEGIN_MASTER( myThid )

      vec_start=irecglob
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_PUT_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_var)

        _END_MASTER( myThid )

      endif

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> ADJOINT RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. REVERSE_SIMULATION) then

        _BEGIN_MASTER( myThid )

      vec_start=irec
      vec_count=1
      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_GET_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_data_t)

c     Add active_var to data.
      active_var = active_var + active_data_t
      active_data_t = 0. _d 0

      vec_start=irecglob
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_PUT_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_data_t)

        _END_MASTER( myThid )

      endif

c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>> TANGENT RUN <<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>             <<<<<<<<<<<<<<<<<<<

      if (theSimulationMode .eq. TANGENT_SIMULATION) then

        _BEGIN_MASTER( myThid )

      vec_start=irecglob
      vec_count=1

      err = NF_INQ_VARID(fid, ihop_nameequi, varid1)
      err = NF_PUT_VARA_DOUBLE(fid, varid1, vec_start, vec_count,
     & active_var)

        _END_MASTER( myThid )

      endif

#endif /* ALLOW_IHOP */

      RETURN
      END

C     ==================================================================

