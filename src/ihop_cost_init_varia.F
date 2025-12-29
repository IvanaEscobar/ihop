#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_COST_INIT_VARIA

C !INTERFACE:
      SUBROUTINE IHOP_COST_INIT_VARIA( myThid )

C    !DESCRIPTION:
C    Initialise variable cost function part

C    !USES:
      IMPLICIT NONE
C == Global Variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#ifdef ALLOW_COST
# include "IHOP_SIZE.h"
# include "IHOP_COST.h"
#endif
#ifdef ALLOW_GRDCHK
# include "netcdf.inc"
#endif

C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER myThid
C     !OUTPUT PARAMETERS:

C     !FUNCTIONS:
#ifdef ALLOW_GRDCHK
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK
#endif /* ALLOW_GRDCHK */

C     !LOCAL VARIABLES:
#ifdef ALLOW_COST
      INTEGER bi, bj
      INTEGER num_file
# ifdef ALLOW_GRDCHK
      INTEGER err
# endif /* ALLOW_GRDCHK */
CEOP

      DO bj=1,35604
        geninfluence(bj) = 0. _d 0
      ENDDO

C Initialize tiled cost fxn contributions.
      DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
          DO num_file = 1, NFILESMAX_IHOP
            ihop_dummy(num_file,bi,bj) = 0. _d 0
# ifdef ALLOW_GRDCHK
            IF ( obsNo_tiled(num_file,bi,bj).GT.0 .AND.
     &           ihopDoNcOutput ) THEN
C Sync tiled .equi file so data in buffer is not lost
              err = NF_SYNC( ncidFWD(num_file,bi,bj) )
              CALL IHOP_COST_NF_ERROR(
     &             'INIT_VARIA: NF_SYNC ncidFWD',
     &             err,1,1,myThid )
#  ifdef ALLOW_TANGENTLINEAR_RUN
          IF (TANGENT_SIMULATION) THEN
            err = NF_SYNC( ncidTL(num_filei,bi,bj) ) 
            CALL IHOP_COST_NF_ERROR(
     &           'INIT_VARIA: NF_SYNC ncidTL',
     &           err,1,1,myThid )
          ENDIF
#  endif /* ALLOW_TANGENTLINEAR_RUN */

            ENDIF !IF obsNo_tiled
# endif /* ALLOW_GRDCHK */

          ENDDO !DO num_file
        ENDDO !DO bi
      ENDDO !DO bj

C Initialize global cost fxn contributions.
      _BEGIN_MASTER(myThid)

      DO num_file = 1, NFILESMAX_IHOP
        objf_ihop(num_file) = 0. _d 0
        num_ihop(num_file)  = 0. _d 0
        ihop_globaldummy(num_file) = 0. _d 0
      ENDDO !DO num_file

      _END_MASTER(myThid)

C# ifdef ALLOW_GRDCHK
C      DO num_file = 1, NFILESMAX_IHOP
C        IF (obsNo(num_file).GT.0 .AND. ihopDoNcOutput) THEN
CC Sync global .equi file so data in buffer is not lost
C          err = NF_SYNC( ncidGLOB(num_file) )
C          CALL IHOP_COST_NF_ERROR(
C     &         'INIT_VARIA: NF_SYNC ncidGLOB',
C     &         err,1,1,myThid )
C#  ifdef ALLOW_TANGENTLINEAR_RUN
C          IF (TANGENT_SIMULATION) THEN
C            err = NF_SYNC( ncidTLGLOB(num_file) ) 
C            CALL IHOP_COST_NF_ERROR(
C     &           'INIT_VARIA: NF_SYNC ncidTLGLOB',
C     &           err,1,1,myThid )
C          ENDIF
C#  endif /* ALLOW_TANGENTLINEAR_RUN */
C        ENDIF
C      ENDDO !DO num_file
C
C# endif /* ALLOW_GRDCHK */

CC Initialize global variables for ihop_read_obs.F
C      ihop_curfile_buff = 0
C      DO l = 1, 1000
C        ihop_data_buff(l) = 0. _d 0
C        ihop_uncert_buff(l) = 0. _d 0
C      ENDDO

      _BARRIER
#endif /* ALLOW_COST */

      RETURN
      END
