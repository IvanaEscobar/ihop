#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_COST_INIT_EQUI

C !INTERFACE:
      SUBROUTINE IHOP_COST_INIT_EQUI(
     I                               num_file,
     I                               fID1,
     I                               file2,
     O                               fID2,
     I                               length,
     I                               bi,
     I                               bj,
     I                               myThid )

C    !DESCRIPTION:
C    Initialization of tiled model equivalent NetCDF files for iHOP obs data

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
      INTEGER num_file, fID1
      CHARACTER*(*) file2
      INTEGER length
      INTEGER bi, bj, myThid
C     !OUTPUT PARAMETERS:
      INTEGER fID2

C     !LOCAL VARIABLES:
#if (defined ALLOW_IHOP) && (defined ALLOW_COST)
      INTEGER dimID, varID(3), varID0
      INTEGER irec, err
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      _RL tmp_vec(2)
CEOP

      IF (ihopDoNcOutput) THEN
C Create a netcdf file
        err = NF_CREATE( file2, NF_CLOBBER, fID2 )
        CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_CREATE fID2',
     &       err,bi,bj,myThid )
        err = NF_DEF_DIM( fID2, 'iOBS', length, dimID )
        CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_DEF_DIM iOBS',
     &       err,bi,bj,myThid )

        err = NF_DEF_VAR( fID2, 'ihopObs_ind_glob_tiled',
     &        NF_INT, 1, dimID, varID(1) )
        CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_DEF_VAR ind_glob',
     &       err,bi,bj,myThid )
        err = NF_PUT_ATT_INT( fID2, varID(1),
     &        '_FillValue', NF_INT, 1, 0 )
        CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_PUT_ATT_INT iOBS',
     &       err,bi,bj,myThid )

! If 'ihop_val' exists in data file, then add info to equi file
        err = NF_INQ_VARID( fID1, ihop_nameval, varID0 )
        CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_INQ_VARID ihop_val',
     &       err,bi,bj,myThid )

        IF (err.EQ.NF_NOERR) THEN
          err = NF_DEF_VAR( fID2, ihop_nameequi,
     &          NF_DOUBLE, 1, dimID, varID(2) )
          CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_DEF_VAR val',
     &         err,bi,bj,myThid )
          err = NF_PUT_ATT_DOUBLE( fID2, varID(2),
     &          '_FillValue', NF_DOUBLE, 1, 0. _d 0 )
          CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_PUT_ATT_D val',
     &         err,bi,bj,myThid )

          err = NF_DEF_VAR( fID2, ihop_namemask,
     &          NF_DOUBLE, 1, dimID, varID(3) )
          CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_DEF_VAR mask',
     &         err,bi,bj,myThid )
          err = NF_PUT_ATT_DOUBLE(fID2,varID(3),
     &          '_FillValue', NF_DOUBLE, 1, 0. _d 0 )
          CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_PUT_ATT_D mask',
     &         err,bi,bj,myThid )

        ELSE
          WRITE(msgBuf,'(2A)')
     &         'S/R IHOP_COST_INIT_EQUI: no variable ihop_val',
     &         ihop_nameval
          CALL PRINT_MESSAGE( msgBuf,
     &         standardMessageUnit, SQUEEZE_RIGHT, myThid )

        ENDIF

! Transition NetCDF file to data mode
        err = NF_ENDDEF( fID2 )
        CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_ENDDEF fID2',
     &       err,bi,bj,myThid )

! Clear buffer and write data to file
        err = NF_CLOSE( fID2 )
        CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_CLOSE fID2',
     &       err,bi,bj,myThid )
        err = NF_OPEN( file2, NF_WRITE, fID2 )
        CALL IHOP_COST_NF_ERROR( 'INIT_EQUI: NF_OPEN'//file2//'',
     &       err,bi,bj,myThid )

      ELSE !IF (ihopDoNcOutput)
        WRITE(msgBuf,'(2A)') "**ERROR** IHOP COST doesn't support ",
     &    "binary model equivalents est.2024"
        CALL PRINT_MESSAGE( msgBuf,
     &       errorMessageUnit, SQUEEZE_RIGHT, myThid )
        STOP "ABNORMAL END: S/R IHOP_COST_INIT_EQUI"

        OPEN( fID2, FILE = file2, FORM = 'unformatted', 
     &   STATUS = 'unknown', ACCESS = 'direct', 
     &   RECL = 2*WORDLENGTH*2 )

        DO irec = 1, 2
          tmp_vec(irec) = 0
        ENDDO
# ifdef _BYTESWAPIO
        CALL MDS_BYTESWAPR8( 2, temp_vec )
# endif

        DO irec = length, 1, -1
          WRITE(fID2, REC = (irec*2-1)) tmp_vec
          WRITE(fID2, REC = (irec*2)) tmp_vec
        ENDDO

      ENDIF !IF (ihopDoNcOutput)
#endif /* ALLOW_IHOP && ALLOW_COST */

      RETURN
      END
