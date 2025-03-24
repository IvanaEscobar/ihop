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

C#ifdef ALLOW_IHOP
C#ifdef ALLOW_MNC
C
C      IF (useMNC) THEN
C
C        CALL MNC_CW_ADD_VNAME('ihop_soundSpeed',
C     &      'Cen_xy_Hn__L__t', 4, 5, myThid)
C        CALL MNC_CW_ADD_VATTR_TEXT('ihop_soundSpeed',
C     &      'units','[m/s]', myThid)
C        CALL MNC_CW_ADD_VATTR_TEXT('ihop_soundSpeed', 'description',
C     &      'speed of sound in sea water', myThid)
C        CALL MNC_CW_ADD_VATTR_TEXT('ihop_soundSpeed', 'long_name',
C     &      'sound speed', myThid)
C
C      ENDIF
C
C#endif /* ALLOW_MNC */
C#endif /* ALLOW_IHOP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
