#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_INIT_VARIA

C !INTERFACE:
      SUBROUTINE IHOP_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize IHOP variables

C     !USES:
      USE ssp_mod, only: SSP, Grid
      IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

C     !INPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
C     !OUTPUT PARAMETERS: None

#ifdef ALLOW_IHOP
C     !LOCAL VARIABLES:
      INTEGER i,j,k,bi,bj
CEOP

      DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)

C     Initialize arrays in common blocks :
C     3-dim.
#ifdef IHOP_3D_STATE
          DO k=1,Nr
            DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
                ihop_ssp(i,j,k,bi,bj) = 0. _d 0
              ENDDO
            ENDDO
          ENDDO

#endif /* IHOP_3D_STATE */

C     2-dim.
#ifdef IHOP_2D_STATE
          DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
              ihop_sld(i,j,bi,bj) = 0. _d 0
            ENDDO
          ENDDO
#endif /* IHOP_2D_STATE */

C- end bi,bj loops
        ENDDO
      ENDDO

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN
#ifdef IHOP_3D_STATE
#endif /* IHOP_3D_STATE */

#ifdef IHOP_2D_STATE
#endif /* IHOP_2D_STATE */

      ELSE
C-    restart from a pickup:
        CALL IHOP_READ_PICKUP( nIter0, myThid )

C-    end start-from-iter-zero if/else block
      ENDIF


C- GLOBAL reset SSP allocatable arrays (F90 feature)
      IF (ALLOCATED( SSP%cMat )) DEALLOCATE( SSP%cMat )
      IF (ALLOCATED( SSP%czMat )) DEALLOCATE( SSP%czMat )
      ALLOCATE( SSP%cMat( Grid%nZ, Grid%nR ), 
     &          SSP%czMat( Grid%nZ-1, Grid%nR ), 
     &          STAT=i )
      IF ( i.NE.0 ) THEN
C#ifdef IHOP_WRITE_OUT
C        WRITE(msgBuf,'(2A)') 'SSP_MOD::init_varia_SSP: ', &
C          'Insufficient memory to store SSP%cMat, SSP%czMat'
C        CALL PRINT_ERROR( msgBuf,myThid )
C#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R init_varia_SSP'
      ENDIF

      ! Init REAL arrays
      ssp%cmat = -99.d0
      ssp%czmat = -99.d0
      ! Init COMPLEX arrays
      ssp%c = -99.d0
      ssp%cz = -99.d0
C- END GLOBAL reset SSP allocatable arrays (F90 feature)
#endif /* ALLOW_IHOP */

      RETURN
      END
