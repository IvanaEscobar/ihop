#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C     !ROUTINE: IHOP_MNC_INIT

C     !INTERFACE:
      SUBROUTINE IHOP_MNC_INIT( myThid )

C     !DESCRIPTION:
C     Initialize MNC lookup tables for Mypackage.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_IHOP
#ifdef ALLOW_MNC

      IF (useMNC) THEN

        CALL MNC_CW_ADD_VNAME('ihop_soundSpeed',
     &      'Cen_xy_Hn__L__t', 4, 5, myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('ihop_soundSpeed',
     &      'units','[m/s]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('ihop_soundSpeed', 'description',
     &      'speed of sound in sea water', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('ihop_soundSpeed', 'long_name',
     &      'sound speed', myThid)

      ENDIF

#endif /* ALLOW_MNC */
#endif /* ALLOW_IHOP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
