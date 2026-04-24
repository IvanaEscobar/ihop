#include "BELLI_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: BELLI_INIT_FIXED
C !INTERFACE:
      SUBROUTINE BELLI_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     o initialization for netcdf ihop acoustic domain; pass to f90 code
C     started: Ivana Escobar July-2023

C     !USES:
      USE belli_init_mod, only: INIT_FIXED_ENV
      IMPLICIT NONE
C == Global Variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "BELLI_SIZE.h"
#include "BELLI.h"
#ifdef HAVE_NETCDF
#include "netcdf.inc"
#endif

C     !INPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
C     !OUTPUT PARAMETERS: None


#ifdef ALLOW_BELLI
C     !FUNCTIONS:
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_MBUF) :: msgBuf
      INTEGER :: iUnit
      INTEGER :: IL, err, ncid, dimId1, dimId2
      INTEGER :: varId1, varId2, varId3, varId4
      INTEGER :: k, i, j
      _RS :: data1vec( BELLI_MAX_NC_SIZE )
      _RS :: data2vec( BELLI_MAX_NC_SIZE*BELLI_MAX_NC_SIZE )
      _RS :: data3vec( BELLI_MAX_NC_SIZE*BELLI_MAX_NC_SIZE )
      _RL :: data4vec( BELLI_MAX_NC_SIZE*BELLI_MAX_NC_SIZE )
      CHARACTER*(MAX_LEN_FNAM) :: ihop_file
CEOP

C Only do I/O if in the main thread
       _BEGIN_MASTER( myThid )

C Read NetCDF input: setting the acoustic domain
      IL  = ILNBLNK( BELLI_interpfile )
      IF (IL.NE.0) THEN
          WRITE(ihop_file,'(2A)') ihop_interpfile(1:IL),'.nc'
      ELSE
          ihop_file = BELLI_interpfile
      END IF
      err = NF_OPEN(ihop_file, NF_NOWRITE, ncid)
      CALL BELLI_COST_NF_ERROR('INIT_FIXED NF_OPEN',err,1,1,myThid)
      IF (err .NE. NF_NOERR) THEN
        WRITE(msgBuf, '(2A)') 'S/R BELLI_INIT_FIXED: ',
     &    'NetCDF BELLI input file NOT found'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BELLI_INIT_FIXED'
      ENDIF

C Get dimension IDs and lengths
      err = NF_INQ_DIMID( ncid, 'nrange', dimId1 )
      CALL BELLI_COST_NF_ERROR('init_fixed 2',err,1,1,myThid)
      err = NF_INQ_DIMLEN( ncid, dimId1, BELLI_npts_range )

      CALL BELLI_COST_NF_ERROR('init_fixed 3',err,1,1,myThid)
      err = NF_INQ_DIMID( ncid, 'nidw', dimId2 )
      CALL BELLI_COST_NF_ERROR('init_fixed 4',err,1,1,myThid)
      err = NF_INQ_DIMLEN( ncid, dimId2, BELLI_npts_idw )
      CALL BELLI_COST_NF_ERROR('init_fixed 5',err,1,1,myThid)

      IF ( BELLI_npts_range>BELLI_MAX_RANGE .OR.
     &     BELLI_npts_idw>BELLI_MAX_IDW ) THEN
        WRITE(msgBuf, '(2A)') 'S/R BELLI_INIT_FIXED: reading dims ',
     &    'NetCDF Dimensions exceed maximum allowable size'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BELLI_INIT_FIXED'
      ENDIF

C Get data variable IDs
      err = NF_INQ_VARID( ncid, 'belli_ranges', varId1 )
      CALL BELLI_COST_NF_ERROR('init_fixed 6',err,1,1,myThid)
      err = NF_INQ_VARID( ncid, 'belli_xc', varId2 )
      CALL BELLI_COST_NF_ERROR('init_fixed 7',err,1,1,myThid)
      err = NF_INQ_VARID( ncid, 'belli_yc', varId3 )
      CALL BELLI_COST_NF_ERROR('init_fixed 8',err,1,1,myThid)
      err = NF_INQ_VARID( ncid, 'belli_idw_weights', varId4 )
      CALL BELLI_COST_NF_ERROR('init_fixed 9',err,1,1,myThid)

      IF (err .NE. NF_NOERR) THEN
        WRITE(msgBuf, '(2A)') 'S/R BELLI_INIT_FIXED: ',
     &    'NetCDF BELLI data variable names are incorrect'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BELLI_INIT_FIXED'
      ENDIF

C get data variables
      err = NF_GET_VAR( ncid, varId1, data1vec )
      CALL BELLI_COST_NF_ERROR('init_fixed10',err,1,1,myThid)
      err = NF_GET_VAR( ncid, varId2, data2vec )
      CALL BELLI_COST_NF_ERROR('init_fixed11',err,1,1,myThid)
      err = NF_GET_VAR( ncid, varId3, data3vec )
      CALL BELLI_COST_NF_ERROR('init_fixed12',err,1,1,myThid)
      err = NF_GET_VAR( ncid, varId4, data4vec )
      CALL BELLI_COST_NF_ERROR('init_fixed13',err,1,1,myThid)

      IF (err .NE. NF_NOERR) THEN
        WRITE(msgBuf, '(2A)') 'S/R BELLI_INIT_FIXED: ',
     &    'NetCDF BELLI data variables are incorrect'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BELLI_INIT_FIXED'
      ENDIF

C 'reshape' data variables to 2D arrays
      k=1
      DO i=1,BELLI_npts_range
        belli_ranges(i) = data1vec(i)
        DO j=1,BELLI_npts_idw
          belli_xc(i,j) = data2vec(k)
          belli_yc(i,j) = data3vec(k)
          belli_idw_weights(i,j) = data4vec(k)
          k = k+1
        END DO
      END DO

C Close NetCDF file
      err = NF_CLOSE( ncid )
      CALL BELLI_COST_NF_ERROR('init_fixed14',err,1,1,myThid)
      IF (err.NE.NF_NOERR) THEN
        WRITE(msgBuf, '(3A)') 'S/R BELLI_INIT_FIXED: ',
     &    'NetCDF ERROR: Cannot close file ', BELLI_interpfile(1:IL)
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BELLI_INIT_FIXED'
      ENDIF

C Only do I/O if in the main thread
      _END_MASTER( myThid )

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      DO k = 1,BELLI_MAX_NC_SIZE !Nr
        DO i=1,BELLI_MAX_RANGE !BELLI_npts_range
          belli_sumweights(i,k) = -1.0
        END DO
      END DO

C Pass fixed info to f90 ihop domain grid
      CALL INIT_FIXED_ENV( myThid )
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_MNC
      IF (useMNC) THEN
C       Initialize the MNC variable types for BELLI
        CALL BELLI_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL BELLI_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_BELLI */
      RETURN
      END
