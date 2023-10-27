#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: IHOP_READ_PICKUP

C     !INTERFACE:
      SUBROUTINE IHOP_READ_PICKUP( myIter, myThid )

C     !DESCRIPTION:
C     Reads current state of IHOP from a pickup file

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

C     !INPUT PARAMETERS:
C     myIter            :: time-step number
C     myThid            :: thread number
      INTEGER myIter
      INTEGER myThid

#ifdef ALLOW_IHOP
#if (defined IHOP_3D_STATE) || (defined IHOP_2D_STATE)

C     !LOCAL VARIABLES:
C     fn          :: character buffer for creating filename
C     fp          :: precision of pickup files
C     filePrec    :: pickup-file precision (read from meta file)
C     nbFields    :: number of fields in pickup file (read from meta file)
C     missFldList :: List of missing fields   (attempted to read but not found)
C     missFldDim  :: Dimension of missing fields list array: missFldList
C     nMissing    :: Number of missing fields (attempted to read but not found)
C     j           :: loop index
C     nj          :: record number
C     ioUnit      :: temp for writing msg unit
C     msgBuf      :: Informational/error message buffer
      INTEGER fp
      INTEGER filePrec, nbFields
      INTEGER missFldDim, nMissing
      INTEGER j, nj, ioUnit
      PARAMETER( missFldDim = 12 )
      CHARACTER*(10) suff
      CHARACTER*(MAX_LEN_FNAM) fn
      CHARACTER*(8) missFldList(missFldDim)
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C       IF ( pickupSuff.EQ.' ' ) THEN
C        IF ( rwSuffixType.EQ.0 ) THEN
C          WRITE(fn,'(A,I10.10)') 'pickup_ihop.', myIter
C        ELSE
C          CALL RW_GET_SUFFIX( suff, startTime, myIter, myThid )
C          WRITE(fn,'(A,A)') 'pickup_ihop.', suff
C        ENDIF
C       ELSE
C        WRITE(fn,'(A,A10)') 'pickup_ihop.', pickupSuff
C       ENDIF
C       fp = precFloat64
C
C       CALL READ_MFLDS_SET(
C     I                      fn,
C     O                      nbFields, filePrec,
C     I                      Nr, myIter, myThid )
C       _BEGIN_MASTER( myThid )
C       IF ( nbFields.GE.0 .AND. filePrec.NE.fp ) THEN
C         WRITE(msgBuf,'(2A,I4)') 'IHOP_READ_PICKUP: ',
C     &    'pickup-file binary precision do not match !'
C         CALL PRINT_ERROR( msgBuf, myThid )
C         WRITE(msgBuf,'(A,2(A,I4))') 'IHOP_READ_PICKUP: ',
C     &    'file prec.=', filePrec, ' but expecting prec.=', fp
C         CALL PRINT_ERROR( msgBuf, myThid )
C         CALL ALL_PROC_DIE( 0 )
C         STOP 'ABNORMAL END: S/R IHOP_READ_PICKUP (data-prec Pb)'
C       ENDIF
C       _END_MASTER( myThid )
C
C       IF ( nbFields.LE.0 ) THEN
CC-      No meta-file or old meta-file without List of Fields
C        ioUnit = errorMessageUnit
C        IF ( pickupStrictlyMatch ) THEN
C          WRITE(msgBuf,'(4A)') 'IHOP_READ_PICKUP: ',
C     &      'no field-list found in meta-file',
C     &      ' => cannot check for strick-matching'
C          CALL PRINT_ERROR( msgBuf, myThid )
C          WRITE(msgBuf,'(4A)') 'IHOP_READ_PICKUP: ',
C     &      'try with " pickupStrictlyMatch=.FALSE.,"',
C     &      ' in file: "data", NameList: "PARM03"'
C          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
C          CALL ALL_PROC_DIE( myThid )
C          STOP 'ABNORMAL END: S/R IHOP_READ_PICKUP'
C        ELSE
C          WRITE(msgBuf,'(4A)') 'WARNING >> IHOP_READ_PICKUP: ',
C     &      ' no field-list found'
C          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
C         IF ( nbFields.EQ.-1 ) THEN
CC-      No meta-file
C          WRITE(msgBuf,'(4A)') 'WARNING >> ',
C     &      ' try to read pickup as currently written'
C          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
C         ELSE
C          WRITE(msgBuf,'(4A)') 'IHOP_READ_PICKUP: ',
C     &      'no field-list found in meta-file'
C          CALL PRINT_ERROR( msgBuf, myThid )
C          CALL ALL_PROC_DIE( myThid )
C          STOP 'ABNORMAL END: S/R IHOP_READ_PICKUP'
C         ENDIF
C        ENDIF
C       ENDIF
C
CC---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C
C       IF ( nbFields.EQ.0 ) THEN
CC---   Old way to read pickup:
C
C       ELSE
CC---   New way to read IHOP pickup:
C        nj = 0
CC---    read IHOP 3-D fields for restart
C#ifdef IHOP_3D_STATE
C        CALL READ_MFLDS_3D_RL( 'ihop_ssp', ihop_ssp,
C     &                                 nj, fp, Nr, myIter, myThid )
C#endif /* IHOP_3D_STATE */
C        nj = nj*Nr
CC---    read IHOP 2-D fields for restart
C#ifdef IHOP_2D_STATE
C        CALL READ_MFLDS_3D_RL( 'myPaSur1', myPa_Surf1,
C     &                                 nj, fp, 1 , myIter, myThid )
C        CALL READ_MFLDS_3D_RL( 'myPaSur2', myPa_Surf2,
C     &                                 nj, fp, 1 , myIter, myThid )
C#endif /* IHOP_2D_STATE */
C
CC--    end: new way to read pickup file
C       ENDIF
C
CC--    Check for missing fields:
C       nMissing = missFldDim
C       CALL READ_MFLDS_CHECK(
C     O                     missFldList,
C     U                     nMissing,
C     I                     myIter, myThid )
C       IF ( nMissing.GT.missFldDim ) THEN
C         WRITE(msgBuf,'(2A,I4)') 'IHOP_READ_PICKUP: ',
C     &     'missing fields list has been truncated to', missFldDim
C         CALL PRINT_ERROR( msgBuf, myThid )
C         CALL ALL_PROC_DIE( myThid )
C         STOP 'ABNORMAL END: S/R IHOP_READ_PICKUP (list-size Pb)'
C       ENDIF
C       IF ( nMissing.GE.1 ) THEN
C        ioUnit = errorMessageUnit
C        DO j=1,nMissing
C         WRITE(msgBuf,'(4A)') 'IHOP_READ_PICKUP: ',
C     &       'cannot restart without field "',missFldList(nj),'"'
C         CALL PRINT_ERROR( msgBuf, myThid )
C        ENDDO
C        CALL ALL_PROC_DIE( myThid )
C        STOP 'ABNORMAL END: S/R IHOP_READ_PICKUP'
C       ENDIF
C
CC--    Update overlap regions:
C#ifdef IHOP_3D_STATE
C        CALL EXCH_3D_RL( ihop_ssp, Nr, myThid )
C#endif /* IHOP_3D_STATE */
C#ifdef IHOP_2D_STATE
C        CALL EXCH_XY_RL( myPa_Surf1, myThid )
C        CALL EXCH_XY_RL( myPa_Surf2, myThid )
C#endif /* IHOP_2D_STATE */

#endif /* IHOP_3D_STATE or IHOP_2D_STATE */
#endif /* ALLOW_IHOP */

      RETURN
      END
