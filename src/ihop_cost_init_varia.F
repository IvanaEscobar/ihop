#include "IHOP_OPTIONS.h"

      SUBROUTINE IHOP_COST_INIT_VARIA( myThid )

C     ==================================================================
C     SUBROUTINE IHOP_COST_INIT_VARIA
C     ==================================================================
C     o Initialise the variable cost function part.
C     added ihop term: ivana@oden.utexas.edu May-2024
C
C     ==================================================================

      IMPLICIT NONE

C     == global variables ==
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

C     == routine arguments ==
      INTEGER myThid

#ifdef ALLOW_COST
C     == local variables ==
      INTEGER bi,bj
      INTEGER l,num_file
# ifdef ALLOW_GRDCHK
      integer err, iG, jG, IL, JL
      character*(max_len_fnam) iFile
      character*(max_len_fnam) fnameequinc, fnameequincglo
# endif /* ALLOW_GRDCHK */

C     == external functions ==
# ifdef ALLOW_GRDCHK
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK
# endif /* ALLOW_GRDCHK */
C     == end of interface ==

C--   Initialize tiled cost fxn contributions.
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO num_file=1,NFILESMAX_IHOP
         ihop_dummy(num_file,bi,bj)=0. _d 0
# ifdef ALLOW_GRDCHK
         IF (obsNo_tiled(num_file,bi,bj).GT.0) THEN
C Close .equi file so data in buffer is not lost when run finishes
          err = NF_CLOSE(fidfwd_obs(num_file,bi,bj))

C Determine model equivalent file name
          iG=bi+(myXGlobalLo-1)/sNx ! Kludge until unstructered tiles
          jG=bj+(myYGlobalLo-1)/sNy ! Kludge until unstructered tiles
          IL  = ILNBLNK( ihopObsFiles(num_file) )
          WRITE(iFile,'(A)')
     &    ihopObsFiles(num_file)(1:IL)
          IL  = ILNBLNK( iFile )
          JL  = ILNBLNK( ihopObsDir )
          WRITE(fnameequinc,'(3A,2(I3.3,A))')
     &    ihopObsDir(1:JL),iFile(1:IL),'.',
     &    iG,'.',jG,'.equi.nc'

c re-open to reset write permissions
          JL  = ILNBLNK( fnameequinc )
          err = NF_OPEN( fnameequinc(1:JL),NF_WRITE,
     &    fidfwd_obs(num_file,bi,bj) )
         ENDIF
# endif /* ALLOW_GRDCHK */
        ENDDO
       ENDDO
      ENDDO

C--   Initialize global cost fxn contributions.
      DO num_file=1,NFILESMAX_IHOP
       num_ihop(num_file)  =0. _d 0
       objf_ihop(num_file) =0. _d 0
       ihop_globaldummy(num_file) =0. _d 0
# ifdef ALLOW_GRDCHK
         IF (obsNo(num_file).GT.0) THEN
C Close .equi file so data in buffer is not lost when run finishes
          err = NF_CLOSE(fidglobal(num_file))

C Determine model equivalent file name
          IL  = ILNBLNK( ihopObsFiles(num_file) )
          WRITE(iFile,'(A)')
     &    ihopObsFiles(num_file)(1:IL)
          IL  = ILNBLNK( iFile )
          JL  = ILNBLNK( ihopObsDir )
          WRITE(fnameequincglo,'(3A,2(I3.3,A))')
     &    ihopObsDir(1:JL),iFile(1:IL),'.equi.nc'

c re-open to reset write permissions
          JL  = ILNBLNK( fnameequincglo )
          err = NF_OPEN( fnameequincglo(1:JL),NF_WRITE,
     &    fidglobal(num_file) )
         ENDIF
# endif /* ALLOW_GRDCHK */
      ENDDO

C Initialize global variables for ihop_read_obs.F
      ihop_curfile_buff=0
      DO l=1,1000
       ihop_data_buff(l)=0. _d 0
       ihop_uncert_buff(l)=0. _d 0
      ENDDO

      _BARRIER
#endif /* ALLOW_COST */

      RETURN
      END
