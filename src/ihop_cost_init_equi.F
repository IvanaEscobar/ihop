#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
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

C    !DESCRIPTION:
C    Initialization of model equivalent NetCDF files for iHOP obs data

C    !USES:
      IMPLICIT NONE
C == Global Variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "netcdf.inc"
# ifdef ALLOW_COST
# include "IHOP_COST.h"
# endif
#endif

C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER num_file, fid1
      CHARACTER*(*) file2
      INTEGER length
      INTEGER bi, bj, myThid
C     !OUTPUT PARAMETERS:
      INTEGER fid2
CEOP

C     !LOCAL VARIABLES:
#if (defined ALLOW_IHOP) && (defined ALLOW_COST)
      INTEGER dimId, varId(3), varId0
      INTEGER irec, err
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      _RL tmp_vec(2)

      IF (ihopDoNcOutput) THEN
C Create a netcdf file
        err = NF_CREATE( file2, NF_CLOBBER, fid2 )
        err = NF_DEF_DIM( fid2, 'iOBS', length, dimId )

        err = NF_DEF_VAR( fid2, 'ihopObs_ind_glob_tiled',
     &        NF_INT, 1, dimId,varId(1) )
        err = NF_PUT_ATT_INT( fid2, varId(1),
     &        '_FillValue', NF_INT, 1, -9999 )

! If 'ihop_val' exists in data file, then add info to equi file
        err = NF_INQ_VARID( fid1, ihop_nameval, varId0 )
        IF (err.EQ.NF_NOERR) THEN
          err = NF_DEF_VAR( fid2, ihop_nameequi,
     &          NF_DOUBLE, 1, dimId, varId(2) )
          err = NF_PUT_ATT_DOUBLE( fid2, varId(2),
     &          '_FillValue', NF_DOUBLE, 1, -9999. _d 0 )

          err = NF_DEF_VAR( fid2, ihop_namemask,
     &          NF_DOUBLE, 1, dimId, varId(3) )
          err = NF_PUT_ATT_DOUBLE(fid2,varId(3),
     &          '_FillValue', NF_DOUBLE, 1, -9999. _d 0 )

        ELSE
         WRITE(msgBuf,'(2A)')
     &        'S/R IHOP_COST_INIT_EQUIFILES: no variable ',
     &        ihop_nameval
         CALL PRINT_MESSAGE( msgBuf,
     &        standardMessageUnit, SQUEEZE_RIGHT, myThid )

        ENDIF

! Transition NetCDF file to data mode
        err = NF_ENDDEF( fid2 )

! Clear buffer and write data to file
        err = NF_CLOSE( fid2 )
        err = NF_OPEN( file2, NF_WRITE, fid2 )

      ELSE !IF (ihopDoNcOutput)
        PRINT *, "IHOP doesn't support binary model equivalents 2024"
        OPEN( fid2, FILE = file2, FORM = 'unformatted', 
     &   STATUS = 'unknown', ACCESS = 'direct', 
     &   RECL = 2*WORDLENGTH*2 )

        DO irec = 1, 2
          tmp_vec(irec) = 0
        ENDDO
# ifdef _BYTESWAPIO
        CALL MDS_BYTESWAPR8( 2, temp_vec )
# endif

        DO irec = length, 1, -1
          WRITE(fid2, REC = (irec*2-1)) tmp_vec
          WRITE(fid2, REC = (irec*2)) tmp_vec
        ENDDO

      ENDIF !IF (ihopDoNcOutput)
#endif /* ALLOW_IHOP && ALLOW_COST */

      RETURN
      END
