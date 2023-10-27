#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: IHOP_CHECK

C !INTERFACE:
      SUBROUTINE IHOP_CHECK( myThid )

C     !DESCRIPTION:
C     Check consistency with model configuration

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

C     !INPUT PARAMETERS:
C     myThid   :: my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_IHOP
C     !LOCAL VARIABLES:
C     msgBuf   :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER errCount
CEOP

      _BEGIN_MASTER(myThid)

       WRITE(msgBuf,'(A)') 'IHOP_CHECK: #define IHOP'
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , 1)

C- print out some key parameters :
       CALL WRITE_0D_C( IHOP_fileroot, 0, INDEX_NONE,
     &  'IHOP_fileroot=',
     &  ' /* Name of file root */')

      _END_MASTER(myThid)

#endif /* ALLOW_IHOP */

      RETURN
      END
