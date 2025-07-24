#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_DIAGNOSTICS_INIT

C !INTERFACE:
      SUBROUTINE IHOP_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Initialize list of all available diagnostics

C     !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     === Local variables ===
C     diagNum   :: diagnostics number in list of available diag
C     diagMate  :: diag. mate number in list of available diag
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of diags
C     diagUnits :: local string (16c): physical units of a diag
C     diagTitle :: local string (80c): description of field in diag
      INTEGER       diagNum
      INTEGER       diagMate
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle
CEOP

#ifdef IHOP_3D_STATE
        diagName  = 'SVEL    '
        diagTitle = 'Sound speed from pkg/ihop'
        diagUnits = 'm/s             '
        diagCode  = 'SMR     MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* IHOP_3D_STATE */

        diagName  = 'SLD     '
        diagTitle = 'Sonic layer depth (>0) from pkg/ihop'
        diagUnits = 'm               '
        diagCode  = 'SM      M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
c     ENDIF

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END