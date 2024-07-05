#include "IHOP_OPTIONS.h"

CBOP
C !ROUTINE: IHOP_COST_INIT_EQUIFILES

C !INTERFACE: 
      SUBROUTINE IHOP_COST_INIT_EQUIFILES( 
     I                               num_file,
     I                               fid1,
     I                               file2,
     O                               fid2,
     I                               length,
     I                               bi,
     I                               bj,
     I                               myThid )

C !DESCRIPTION:
C Initialization of model counterparts files             
C for ihop obs data 

C !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "IHOP_COST.h"
# include "netcdf.inc"
#endif
CEOP

C !INPUT PARAMETERS:
      INTEGER num_file, fid1, fid2
      INTEGER bi, bj, myThid, length
      CHARACTER*(MAX_LEN_FNAM) file2

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)
C !LOCAL VARIABLES:
      INTEGER dimid,varid(3),varid0,err
      INTEGER irec
      Real*8  tmp_vec(2)

C Create a netcdf file
      err = NF_CREATE(file2,NF_CLOBBER,fid2)
      err = NF_DEF_DIM(fid2,'iOBS',length,dimid)
      err = NF_DEF_VAR(fid2,'ihopObs_ind_glob_tiled',NF_INT,1,
     &                  dimid,varid(1))
      err = NF_PUT_ATT_INT(fid2,varid(1),'_FillValue',NF_INT,1,-9999)

! Check if 'ihop_val' exists in the observations data file
      err = NF_INQ_VARID(fid1,ihop_nameval,varid0)

      IF (err.EQ.NF_NOERR) THEN

       err = NF_DEF_VAR(fid2,ihop_nameequi,
     &  NF_DOUBLE,1,dimid,varid(2))
       err = NF_PUT_ATT_DOUBLE(fid2,varid(2),
     &  '_FillValue',NF_DOUBLE,1,-9999. _d 0)
       err = NF_DEF_VAR(fid2,ihop_namemask,
     &  NF_DOUBLE,1,dimid,varid(3))
       err = NF_PUT_ATT_DOUBLE(fid2,varid(3),
     &  '_FillValue',NF_DOUBLE,1,-9999. _d 0)

      ENDIF

      PRINT *, 'escobar: EQUIFILES err 1', err

      err=NF_ENDDEF(fid2)
      err=NF_CLOSE(fid2)

      PRINT *, 'escobar: EQUIFILES err 2', err

      err = NF_OPEN(file2,NF_WRITE,fid2)

      PRINT *, 'escobar: EQUIFILES err 3', err
#endif

      END
