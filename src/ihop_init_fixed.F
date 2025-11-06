#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_INIT_FIXED
C !INTERFACE:
      SUBROUTINE IHOP_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     o initialization for netcdf ihop acoustic domain; pass to f90 code
C     started: Ivana Escobar July-2023

C     !USES:
      USE ihop_init_mod, only: INIT_FIXED_ENV
      IMPLICIT NONE
C == Global Variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef HAVE_NETCDF
#include "netcdf.inc"
#endif

C     !INPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
C     !OUTPUT PARAMETERS: None


#ifdef ALLOW_IHOP
C     !FUNCTIONS:
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_MBUF) :: msgBuf
      INTEGER :: iUnit
      INTEGER :: IL, err, ncid, dimId1, dimId2
      INTEGER :: varId1, varId2, varId3, varId4
      INTEGER :: k, i, j
      _RS :: data1vec( IHOP_MAX_NC_SIZE )
      _RS :: data2vec( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      _RS :: data3vec( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      _RL :: data4vec( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      CHARACTER*(MAX_LEN_FNAM) :: ihop_file
CEOP

C Only do I/O if in the main thread
       _BEGIN_MASTER( myThid )

C Read NetCDF input: setting the acoustic domain
      IL  = ILNBLNK( IHOP_interpfile )
      IF (IL.NE.0) THEN
          WRITE(ihop_file,'(2A)') ihop_interpfile(1:IL),'.nc'
      ELSE
          ihop_file = IHOP_interpfile
      END IF
      err = NF_OPEN(ihop_file, NF_NOWRITE, ncid)
      CALL IHOP_COST_NF_ERROR('INIT_FIXED NF_OPEN',err,1,1,myThid)
      IF (err .NE. NF_NOERR) THEN
        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: ',
     &    'NetCDF IHOP input file NOT found'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
      ENDIF

C Get dimension IDs and lengths
      err = NF_INQ_DIMID( ncid, 'nrange', dimId1 )
      CALL IHOP_COST_NF_ERROR('init_fixed 2',err,1,1,myThid)
      err = NF_INQ_DIMLEN( ncid, dimId1, IHOP_npts_range )

      CALL IHOP_COST_NF_ERROR('init_fixed 3',err,1,1,myThid)
      err = NF_INQ_DIMID( ncid, 'nidw', dimId2 )
      CALL IHOP_COST_NF_ERROR('init_fixed 4',err,1,1,myThid)
      err = NF_INQ_DIMLEN( ncid, dimId2, IHOP_npts_idw )
      CALL IHOP_COST_NF_ERROR('init_fixed 5',err,1,1,myThid)

      IF ( IHOP_npts_range>IHOP_MAX_RANGE .OR.
     &     IHOP_npts_idw>IHOP_MAX_IDW ) THEN
        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: reading dims ',
     &    'NetCDF Dimensions exceed maximum allowable size'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
      ENDIF

C Get data variable IDs
      err = NF_INQ_VARID( ncid, 'ihop_ranges', varId1 )
      CALL IHOP_COST_NF_ERROR('init_fixed 6',err,1,1,myThid)
      err = NF_INQ_VARID( ncid, 'ihop_xc', varId2 )
      CALL IHOP_COST_NF_ERROR('init_fixed 7',err,1,1,myThid)
      err = NF_INQ_VARID( ncid, 'ihop_yc', varId3 )
      CALL IHOP_COST_NF_ERROR('init_fixed 8',err,1,1,myThid)
      err = NF_INQ_VARID( ncid, 'ihop_idw_weights', varId4 )
      CALL IHOP_COST_NF_ERROR('init_fixed 9',err,1,1,myThid)

      IF (err .NE. NF_NOERR) THEN
        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: ',
     &    'NetCDF IHOP data variable names are incorrect'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
      ENDIF

C get data variables
      err = NF_GET_VAR( ncid, varId1, data1vec )
      CALL IHOP_COST_NF_ERROR('init_fixed10',err,1,1,myThid)
      err = NF_GET_VAR( ncid, varId2, data2vec )
      CALL IHOP_COST_NF_ERROR('init_fixed11',err,1,1,myThid)
      err = NF_GET_VAR( ncid, varId3, data3vec )
      CALL IHOP_COST_NF_ERROR('init_fixed12',err,1,1,myThid)
      err = NF_GET_VAR( ncid, varId4, data4vec )
      CALL IHOP_COST_NF_ERROR('init_fixed13',err,1,1,myThid)

      IF (err .NE. NF_NOERR) THEN
        WRITE(msgBuf, '(2A)') 'S/R IHOP_INIT_FIXED: ',
     &    'NetCDF IHOP data variables are incorrect'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
      ENDIF

C 'reshape' data variables to 2D arrays
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

C Close NetCDF file
      err = NF_CLOSE( ncid )
      CALL IHOP_COST_NF_ERROR('init_fixed14',err,1,1,myThid)
      IF (err.NE.NF_NOERR) THEN
        WRITE(msgBuf, '(3A)') 'S/R IHOP_INIT_FIXED: ',
     &    'NetCDF ERROR: Cannot close file ', IHOP_interpfile(1:IL)
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R IHOP_INIT_FIXED'
      ENDIF

C Only do I/O if in the main thread
      _END_MASTER( myThid )

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      DO k = 1,IHOP_MAX_NC_SIZE !Nr
        DO i=1,IHOP_MAX_RANGE !IHOP_npts_range
          ihop_sumweights(i,k) = -1.0
        END DO
      END DO

C Pass info to f90 ihop domain grid
      CALL init_fixed_env( myThid )
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
