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
      INTEGER :: IL, err, ncid, dimid1, dimid2
      INTEGER :: varid1, varid2, varid3, varid4
      INTEGER :: k, i, j
      _RS :: data1vec( IHOP_MAX_NC_SIZE )
      _RS :: data2vec( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      _RS :: data3vec( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      _RL :: data4vec( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      CHARACTER*(80) :: ihop_file

C ==================== External Functions ==========================
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C Only do I/O if in the main thread
       _BEGIN_MASTER( myThid )
 
C Read NetCDF input: setting the acoustic domain
      IL  = ILNBLNK( IHOP_interpfile )
      IF (IL.NE.0) THEN
          WRITE(ihop_file,'(2a)') ihop_interpfile(1:IL),'.nc'
      ELSE
          ihop_file = IHOP_interpfile
      END IF
      err = NF_OPEN(ihop_file, NF_NOWRITE, ncid)
      IF (err .NE. NF_NOERR) THEN
        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: ',
     &  'NetCDF IHOP input file NOT found'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
      END IF

C Get dimension IDs and lengths
      err = NF_INQ_DIMID( ncid, 'nrange', dimid1 )
      err = NF_INQ_DIMLEN( ncid, dimid1, IHOP_npts_range )

      err = NF_INQ_DIMID( ncid, 'nidw', dimid2 )
      err = NF_INQ_DIMLEN( ncid, dimid2, IHOP_npts_idw )

      IF ( IHOP_npts_range>IHOP_MAX_NC_SIZE .OR. 
     &     IHOP_npts_idw>IHOP_MAX_NC_SIZE ) THEN
        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: reading dims',
     &  'NetCDF Dimensions exceed maximum allowable size'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
      END IF

C Get data variable IDs
      err = NF_INQ_VARID( ncid, 'ihop_ranges', varid1 )
      err = NF_INQ_VARID( ncid, 'ihop_xc', varid2 )
      err = NF_INQ_VARID( ncid, 'ihop_yc', varid3 )
      err = NF_INQ_VARID( ncid, 'ihop_idw_weights', varid4 )

      IF (err .NE. NF_NOERR) THEN
        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: ',
     &  'NetCDF IHOP data variable names are incorrect'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
      END IF

C get data variables
      err = NF_GET_VAR( ncid, varid1, data1vec )
      err = NF_GET_VAR( ncid, varid2, data2vec )
      err = NF_GET_VAR( ncid, varid3, data3vec )
      err = NF_GET_VAR( ncid, varid4, data4vec )

      IF (err .NE. NF_NOERR) THEN
        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: ',
     &  'NetCDF IHOP data variables are incorrect'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
      END IF

C 'reshape' data variables to 2D array
      k=1
      DO i=1,IHOP_npts_range
       ihop_ranges(i) = data1vec(i)
       DO j=1,IHOP_npts_idw
        ihop_xc(i,j) = data2vec(k)
        ihop_yc(i,j) = data3vec(k)
        ihop_idw_weights(i,j) = data4vec(k)
        k = k+1
       END DO
      END DO

CC Close NetCDF file
      err = NF_CLOSE( ncid )
      IF (err .NE. NF_NOERR) THEN
        WRITE(msgBuf, '(3A)') 'S/R IHOP_INIT_FIXED: ',
     &  'NetCDF ERROR: Cannot close file ', IHOP_interpfile(1:IL)
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
      END IF

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
