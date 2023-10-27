#include "IHOP_OPTIONS.h"

C     *==========================================================*
C     | subroutine ihop_init_fixed
C     | o initialization for netcdf ihop data
C     | started: Ivana Escobar July-2023
C     *==========================================================*
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: IHOP_INIT_FIXED

C     !INTERFACE:
      SUBROUTINE IHOP_INIT_FIXED( myThid )

C     !USES:
      IMPLICIT NONE
C ==================== Global Variables ===========================
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef HAVE_NETCDF
#include "netcdf.inc"
#endif

C ==================== Routine Arguments ==========================
C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C ==================== Routine Variables ==========================

#ifdef ALLOW_IHOP
C     !EXTERNAL VARIABLES:
      CHARACTER*(MAX_LEN_MBUF) :: msgBuf
C     !LOCAL VARIABLES:
      INTEGER :: iUnit
      INTEGER :: IL, err, ncid, dimid1, dimid2, varid1, varid2, varid3
      INTEGER :: k, i, j
      _RS :: data1vec( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      _RS :: data2vec( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      _RL :: data3vec( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      CHARACTER*(80) :: ihop_file

CC ==================== External Functions ==========================
C      INTEGER ILNBLNK
C      EXTERNAL ILNBLNK
CC---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C Only do I/O if in the main thread
       _BEGIN_MASTER( myThid )
 
CC Read NetCDF input: setting the acoustic domain
C      IL  = ILNBLNK( IHOP_interpfile )
C      IF (IL.NE.0) THEN
C          WRITE(ihop_file(1:80), '(1A)') IHOP_interpfile(1:IL)
C      ELSE
C          ihop_file = IHOP_interpfile
C      END IF
C      err = NF_OPEN(ihop_file, NF_NOWRITE, ncid)
C      IF (err .NE. NF_NOERR) THEN
C        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: ',
C     &  'NetCDF IHOP input file NOT found'
C        CALL PRINT_ERROR( msgBuf, myThid )
C        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
C      END IF
C
CC Get dimension IDs and lengths
C      err = NF_INQ_DIMID( ncid, 'nrange', dimid1 )
C      err = NF_INQ_DIMLEN( ncid, dimid1, IHOP_npts_range )
C
C      err = NF_INQ_DIMID( ncid, 'nidw', dimid2 )
C      err = NF_INQ_DIMLEN( ncid, dimid2, IHOP_npts_idw )
C
C      IF ( IHOP_npts_range>IHOP_MAX_NC_SIZE .OR. 
C     &     IHOP_npts_idw>IHOP_MAX_NC_SIZE ) THEN
C        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: reading dims',
C     &  'NetCDF Dimensions exceed maximum allowable size'
C        CALL PRINT_ERROR( msgBuf, myThid )
C        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
C      END IF
C
CC Get data variable IDs
C      err = NF_INQ_VARID( ncid, 'ihop_xc', varid1 )
C      err = NF_INQ_VARID( ncid, 'ihop_yc', varid2 )
C      err = NF_INQ_VARID( ncid, 'ihop_idw_weights', varid3 )
C
C      IF (err .NE. NF_NOERR) THEN
C        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: ',
C     &  'NetCDF IHOP data variable names are incorrect'
C        CALL PRINT_ERROR( msgBuf, myThid )
C        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
C      END IF
C
CC get data variables
C      err = NF_GET_VAR( ncid, varid1, data1vec )
C      err = NF_GET_VAR( ncid, varid2, data2vec )
C      err = NF_GET_VAR( ncid, varid3, data3vec )
C
C      IF (err .NE. NF_NOERR) THEN
C        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: ',
C     &  'NetCDF IHOP data variables are incorrect'
C        CALL PRINT_ERROR( msgBuf, myThid )
C        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
C      END IF
C
CC 'reshape' data variables to 2D array
C      k=1
C      DO i=1,IHOP_npts_range
C       DO j=1,IHOP_npts_idw
C        ihop_xc(i,j) = data1vec(k)
C        ihop_yc(i,j) = data2vec(k)
C        ihop_idw_weights(i,j) = data3vec(k)
C        k = k+1
C       END DO
C      END DO
C
CCC Close NetCDF file
C      err = NF_CLOSE( ncid )
C      IF (err .NE. NF_NOERR) THEN
C        WRITE(msgBuf, '(3A)') 'S/R IHOP_INIT_FIXED: ',
C     &  'NetCDF ERROR: Cannot close file ', IHOP_interpfile(1:IL)
C        CALL PRINT_ERROR( msgBuf, myThid )
C        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
C      END IF


C 'reshape' ihop_v_* arrays to 2D arrays
      k=1
      DO i=1,IHOP_npts_range
       DO j=1,IHOP_npts_idw
        ihop_xc(i,j) = ihop_v_xc(k)
        ihop_yc(i,j) = ihop_v_yc(k)
        ihop_idw_weights(i,j) = ihop_v_weight(k)
        k = k+1
       END DO
      END DO
C Only do I/O if in the main thread
      _END_MASTER( myThid )

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_MNC
      IF (useMNC) THEN
C       Initialize the MNC variable types for IHOP
        CALL IHOP_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL IHOP_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_IHOP */
      RETURN
      END
