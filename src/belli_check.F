#include "BELLI_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: BELLI_CHECK

C !INTERFACE:
      SUBROUTINE BELLI_CHECK( myThid )

C     !DESCRIPTION:
C     Check consistency with model configuration

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "BELLI_SIZE.h"
#include "BELLI.h"

C     !INPUT PARAMETERS:
C     myThid   :: my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_BELLI
C     !LOCAL VARIABLES:
C     msgBuf   :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER errCount
CEOP

      _BEGIN_MASTER(myThid)

       WRITE(msgBuf,'(A)') 'BELLI_CHECK: #define BELLI'
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , 1)

C- print out some key parameters :
       CALL WRITE_0D_C( TRIM(IHOP_fileroot), 0, INDEX_NONE,
     &  'IHOP_fileroot=',
     &  ' /* Name of file root */')

      _END_MASTER(myThid)

#endif /* ALLOW_BELLI */

      RETURN
      END
