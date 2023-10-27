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

        CALL MNC_CW_ADD_VNAME('ihop_ssp',
     &      'Cen_xy_Hn__L__t', 4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('ihop_ssp',
     &      'units','[m/s]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('ihop_ssp','description',
     &      'scalar state variable ihop_ssp',myThid)

        CALL MNC_CW_ADD_VNAME('myPa_Surf1',
     &      'Cen_xy_Hn__-__t', 3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('myPa_Surf1',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('myPa_Surf1','description',
     &       'Surface scalar field myPa_Surf1',myThid)

        CALL MNC_CW_ADD_VNAME('myPa_Surf2',
     &      'Cen_xy_Hn__-__t', 3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('myPa_Surf2',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('myPa_Surf2','description',
     &       'Surface scalar field myPa_Surf2',myThid)

      ENDIF

#endif /* ALLOW_MNC */
#endif /* ALLOW_IHOP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
