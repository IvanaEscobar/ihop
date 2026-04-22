#include "BELLI_OPTIONS.h"
#ifdef ALLOW_COST
#include "COST_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: BELLI_COST_FINAL
C     !INTERFACE:
      SUBROUTINE BELLI_COST_FINAL( ifc, optimcycle, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE BELLI_COST_FINAL
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_COST
# include "cost.h"
#endif
#include "BELLI_SIZE.h"
#include "BELLI.h"
#include "BELLI_COST.h"

C     !INPUT/OUTPUT PARAMETERS:
C     ifc        :: file unit for costfunction_belli.XXXX
C     optimcycle :: cycle number of the off-line optimization.
C     myThid     :: my Thread Id number
      INTEGER ifc
      INTEGER optimcycle
      INTEGER myThid

#ifdef ALLOW_COST
C     ! FUNCTIONS:
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
      INTEGER num_file_belli
      _RL f_belli(NFILESMAX_BELLI)
      _RL no_belli(NFILESMAX_BELLI)

      CHARACTER*24 cfname
      INTEGER IL
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

c--   Each process has calculated the global part for itself.
      DO num_file_belli=1,NFILESMAX_BELLI
        glofc = glofc + mult_belli(num_file_belli)
     &                 *objf_belli(num_file_belli)
      ENDDO

C     This is just a copy since belli does not have the tile logic of
C     other part of the code
      DO num_file_belli=1,NFILESMAX_BELLI
       f_belli(num_file_belli) = objf_belli(num_file_belli)
       no_belli(num_file_belli) = num_belli(num_file_belli)
      ENDDO

CML I do not understand, why we do not need this or something like this here
CMLc--   Do global summation for each part of the cost function
CML      DO num_var=1,NVARMAX
CML       DO num_file=1,NFILESPROFMAX
CML        _GLOBAL_SUM_RL(f_belli(num_file_belli), myThid )
CML        _GLOBAL_SUM_RL(no_belli(num_file_belli), myThid )
CML       ENDDO
CML      ENDDO

C     start printing to STDOUT
      DO num_file_belli=1,NFILESMAX_BELLI
       IF ( no_belli(num_file_belli).GT.zeroRL ) THEN
        WRITE(msgBuf,'(A,1PE22.14,I2,1X,1PE9.2)')
     &    ' --> f_belli =', f_belli(num_file_belli), num_file_belli,
     &    mult_belli(num_file_belli)
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
       ENDIF
      ENDDO

c--   Each process has calculated the global part for itself.

C     only master thread of master CPU open and write to file
      IF ( ifc .NE. -1 ) THEN

       WRITE(cfname,'(A,I4.4)') 'costfunction_belli.',optimcycle
       OPEN(unit=ifc,file=cfname)
       WRITE(msgBuf,'(A,A)')
     &      'Writing belli cost function info to ', cfname
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )

       DO num_file_belli=1,NFILESMAX_BELLI
        IF ( no_belli(num_file_belli).GT.zeroRL ) THEN
         IL = ILNBLNK( B_Obs_Files(num_file_belli) )
         IL = MAX(IL,30)
         WRITE(ifc,'(4A,1PE22.14,1PE22.14,1X,1PE9.2)')
     &     B_Obs_Files(num_file_belli)(1:IL), ' ', belli_nameval,
     &     ' = ', f_belli(num_file_belli), no_belli(num_file_belli),
     &     mult_belli(num_file_belli)
        ENDIF
       ENDDO
       CLOSE(ifc)

      ENDIF

#endif /* ALLOW_COST */

      RETURN
      END
