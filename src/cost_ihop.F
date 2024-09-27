#include "IHOP_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: COST_IHOP

C !INTERFACE:
      SUBROUTINE COST_IHOP(
     I                        myIter, myTime, myThid )

C !DESCRIPTION:
C ==================================================================
C | Computes cost function for ihop data
C ==================================================================

C !USES:
      IMPLICIT NONE
C == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_CAL
# include "cal.h"
#endif
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "IHOP.h"
# include "netcdf.inc"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif
#ifdef ALLOW_CTRL
! >=c68u
# include "OPTIMCYCLE.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "tamc.h"
#endif

C !INPUT PARAMETERS:
C myTime    :: Current time in simulation
C myIter    :: Current time-step number
C myThid    :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_IHOP
C !LOCAL VARIABLES:
C bi,bj     :: Tile indices
      integer K,num_file,num_var,m
      integer bi,bj,iG,jG
      _RL ihop_modval, ihop_modvaltmp
      _RL ihop_data, ihop_uncert, ihop_weight
#ifndef ALLOW_CTRL
      integer optimcycle
#endif
      character*(max_len_mbuf) msgbuf
      character*(max_len_fnam) ihopObsFile, fnameequinc
      integer IL, JL, err
      integer irec, ii, varid1
      integer samples_cnt, obs_num, num_valid_samples
      logical exst
      _RL  objf_ihop_glo
      _RL  num_ihop_glo
      _RL  tmpgs
      INTEGER tmpgsi
      INTEGER ObsNo2,fid,dimid
      INTEGER nmodmean, nobsmean
      _RL offset, mod_mean, ihop_mean
#ifdef ALLOW_COST
      _RL  ihopObs_buff(NOBSMAX_IHOP)
      _RL  samples_mask_buff(NOBSMAX_IHOP)
      _RL  ihopObs_modval_glob(NOBSMAX_IHOP)
      _RL  samples_mask_glob(NOBSMAX_IHOP)
#endif

C !FUNCTIONS
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C ==================================================================

#ifndef ALLOW_CTRL
      optimcycle = 0
#endif

      write(msgbuf,'(A)') ' '
      CALL print_message( msgbuf,
     &  standardmessageunit,SQUEEZE_RIGHT, myThid)
      write(msgbuf,'(A)') '== cost_ihop: begin =='
      CALL print_message( msgbuf,
     &  standardmessageunit,SQUEEZE_RIGHT, myThid)

        _BEGIN_MASTER( myThid )

#ifdef ALLOW_COST
      JL  = ILNBLNK( ihopObsDir )

      DO num_file=1,NFILESMAX_IHOP
! Init default buffer values
       DO m=1,NOBSMAX_IHOP
        ihopObs_buff(m) = 0.0
       ENDDO


       DO bj=1,nSy
        DO bi=1,nSx

C File maintenance
         IF ( (ObsNo_tiled(num_file,bi,bj).GT.0).AND.
     &         (ihopDoNcOutput) ) THEN
c need to close the file so that the data is not lost when run finishes
           err = NF_CLOSE(fidfwd_obs(num_file,bi,bj))
c then re-open it to compute cost function
           IL  = ILNBLNK( ihopObsFiles(num_file) )
           write(ihopObsFile,'(A)') ihopObsFiles(num_file)(1:IL)
           IL  = ILNBLNK( ihopObsFile )

! Open tiled files
           iG=bi+(myXGlobalLo-1)/sNx ! Kludge until unstructered tiles
           jG=bj+(myYGlobalLo-1)/sNy ! Kludge until unstructered tiles
           write(fnameequinc,'(3A,I3.3,A,I3.3,A)')
     &      ihopObsDir(1:JL),ihopObsFile(1:IL),'.',iG,'.',jG,'.equi.nc'
c
           err = NF_OPEN(fnameequinc,NF_NOWRITE,
     &     fidfwd_obs(num_file,bi,bj))
         ENDIF

        ENDDO !bi
       ENDDO !bj

       DO bj=1,nSy
        DO bi=1,nSx

C Open tiled files and read to buffer
         DO m=1,NOBSMAX_IHOP
          IF (m.LE.ObsNo_tiled(num_file,bi,bj)) THEN

           ihop_modvaltmp=0.

           CALL active_read_ihop_tile(num_file,
     &         ihop_modvaltmp,
     &         m,.false.,optimcycle,bi,bj,myThid,
     &         ihop_dummy(num_file,bi,bj))

C Save model equi (of samples) and masks in buffer
C Combine all threads here
           irec = ihopObs_ind_glob_tiled(num_file,m,bi,bj)
           ihopObs_buff(irec) = ihopObs_buff(irec)
     &                               +ihop_modvaltmp

          ENDIF !if (m.LE.ObsNo_tiled(num_file,bi,bj)) then
         ENDDO !do m=..

        ENDDO !bj
       ENDDO !bi

C Combine MPI processes
       DO ii=1,NOBSMAX_IHOP
        tmpgs = ihopObs_buff(ii)
        _GLOBAL_SUM_RL(tmpgs, myThid)
        ihopObs_modval_glob(ii) = tmpgs
       ENDDO

       IF ( myProcId .eq. 0 ) THEN
C Loop over global obs
        DO m=1,NOBSMAX_IHOP
         ihop_globaldummy(num_file) = 0           !RG
         IF (m.LE.ObsNo(num_file)) THEN

          ihop_modval=0.

C Extract model equi: Keep in case you want NP avg in future dev
          ihop_modval = ihop_modval + ihopObs_modval_glob(m)

c Write to global netcdf file
          CALL active_write_ihop_glob(num_file,
     &     ihop_modval,
     &     m,optimcycle,myThid,
     &     ihop_globaldummy(num_file))

         ENDIF !if (m.LE.ObsNo(num_file)) then
        ENDDO !do m=..
C END Loop over global obs

cc       IF ( myProcId .eq. 0 ) THEN
        err = NF_CLOSE(fidglobal(num_file))
        IL  = ILNBLNK( ihopObsFiles(num_file) )
        write(ihopObsFile,'(A)') ihopObsFiles(num_file)(1:IL)
        IL  = ILNBLNK( ihopObsFile )
        write(fnameequinc,'(3A)')
     &   ihopObsDir(1:JL),ihopObsFile(1:IL),'.equi.nc'
        err = NF_OPEN(fnameequinc,NF_NOWRITE,fidglobal(num_file))
cc       ENDIF

C Loop over global obs
        DO m=1,NOBSMAX_IHOP
         IF (m.LE.ObsNo(num_file)) THEN

#ifdef ALLOW_AUTODIFF_TAMC
!$TAF STORE ihop_globaldummy = comlev1, key=ikey_dynamics, byte=isbyte
#endif

          ihop_modval=0.
          ihop_data=0.
          ihop_uncert=0.

C Read observation and uncertainty
          CALL ihop_cost_read_obs(num_file,1,
     &          ihopObs_ind_glob(num_file,m),
     &          ihop_data,myThid)

          CALL ihop_cost_read_obs(num_file,-1,
     &          ihopObs_ind_glob(num_file,m),
     &          ihop_uncert,myThid)

C Read model equivalent from global file
cc        IF ( myProcId .eq. 0 ) THEN
          CALL active_read_ihop_glob(num_file,
     &     ihop_modval,
     &     m,.false.,optimcycle,myThid,
     &     ihop_globaldummy(num_file))
cc        ENDIF

          IF (ihop_uncert.GT.0.) THEN
           ihop_weight=1/(ihop_uncert*ihop_uncert)
           objf_ihop(num_file)=
     &       objf_ihop(num_file)
     &       +ihop_weight
     &       *(ihop_modval-ihop_data)
     &       *(ihop_modval-ihop_data)

           num_ihop(num_file) = num_ihop(num_file) + 1
          ENDIF

         ENDIF !if (m.LE.ObsNo(num_file)) then
        ENDDO !do m...

        ENDIF !if myprocid = 0
      ENDDO !do num_file=1,NFILESMAX_IHOP

#endif /* ALLOW_COST */
        _END_MASTER( myThid )

#ifdef ALLOW_COST
c Print cost function values
      DO num_file=1,NFILESMAX_IHOP
c
          objf_ihop_glo =
     &             objf_ihop(num_file)
          num_ihop_glo =
     &             num_ihop(num_file)
c
      write(msgbuf,'(A,I2,A,2D12.5)')
     &  ' cost_ihop(', num_file, ' )= ',
     &  objf_ihop_glo, num_ihop_glo

      IF ( num_ihop_glo .GT. 0. ) CALL print_message( msgbuf,
     &  standardmessageunit,SQUEEZE_RIGHT, myThid)
c
      ENDDO
#endif /* ALLOW_COST */

      write(msgbuf,'(A)') '== cost_ihop: end   =='
      CALL print_message( msgbuf,
     &  standardmessageunit,SQUEEZE_RIGHT, myThid)
      write(msgbuf,'(A)') ' '
      CALL print_message( msgbuf,
     &  standardmessageunit,SQUEEZE_RIGHT, myThid)

#endif

C     ==================================================================

      RETURN
      END
