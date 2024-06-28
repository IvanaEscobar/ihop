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
! c68s
# include "optim.h"
! >=c68u
!# include "optimcycle.h"
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
      integer K,num_file,num_var,sample_num
      integer bi,bj,iG,jG
      _RL sample_modval, sample_modvaltmp
      _RL obs_modval
      _RL obs_data, obs_uncert, obs_weight
#ifndef ALLOW_CTRL
      integer optimcycle
#endif
      character*(max_len_mbuf) msgbuf
      character*(max_len_fnam) ihopfile, fnameequinc
      integer IL, JL, err
      integer irec, ii, varid1
      integer samples_cnt, obs_num, num_valid_samples
      _RL sample_mask_sum
      logical exst
      _RL  objf_ihop_glo
      _RL  num_ihop_glo
      _RL  samples_buff(NSAMPLESMAX)
      _RL  samples_mask_buff(NSAMPLESMAX)
      _RL  samples_modval_glob(NSAMPLESMAX)
      _RL  samples_mask_glob(NSAMPLESMAX)
      _RL  tmpgs
      INTEGER tmpgsi
      INTEGER ObsNo2,fid,dimid
      INTEGER nmodmean, nobsmean
      _RL offset, mod_mean, obs_mean
 
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
      DO num_file=1,NFILESMAX_IHOP

C File maintenance
       DO bj=1,nSy
        DO bi=1,nSx

         IF ( (ObsNo_tiled(num_file,bi,bj).GT.0).AND.
     &         (ihopDoNcOutput) ) THEN
c need to close the file so that the data is not lost when run finishes
           err = NF_CLOSE(fidfwd_obs(num_file,bi,bj))
c then re-open it to compute cost function
           iG=bi+(myXGlobalLo-1)/sNx ! Kludge until unstructered tiles
           jG=bj+(myYGlobalLo-1)/sNy ! Kludge until unstructered tiles
           IL  = ILNBLNK( ihopobsfiles(num_file) )
           write(ihopfile(1:80),'(1A)') ihopobsfiles(num_file)(1:IL)
           IL  = ILNBLNK( ihopfile )
           JL  = ILNBLNK( ihopObsDir )
           write(fnameequinc(1:80),'(3A,I3.3,A,I3.3,A)')
     &     ihopObsDir(1:JL),ihopfile(1:IL),'.',iG,'.',jG,'.equi.nc'
c
           err = NF_OPEN(fnameequinc,NF_NOWRITE,
     &     fidfwd_obs(num_file,bi,bj))
         ENDIF

        ENDDO
       ENDDO

c Loop over variables
       DO num_var=1,NOBSMAX_IHOP

c Loop over samples
         DO ii=1,NSAMPLESMAX
           samples_buff(ii) = 0.
           samples_mask_buff(ii) = 0.
         ENDDO

         DO bj=1,nSy
          DO bi=1,nSx

C Open tiled files and read to buffer
           DO sample_num=1,NSAMPLESMAX
            IF (sample_num.LE.ObsNo_tiled(num_file,bi,bj)) THEN

             sample_modvaltmp=0.

             CALL active_read_ihop_tile(num_file,
     &           sample_modvaltmp,
     &           sample_num,.false.,optimcycle,bi,bj,myThid,
     &           ihop_dummy(num_file,bi,bj))

C Save model equi (of samples) and masks in buffer
C Combine all threads here
             irec = ihopObs_ind_glob_tiled(num_file,sample_num,bi,bj)
             samples_buff(irec) = samples_buff(irec)
     &                                 +sample_modvaltmp
             samples_mask_buff(irec) = samples_mask_buff(irec)
     &                                 +sample_modmask(bi,bj)

            ENDIF !if (sample_num.LE.ObsNo_tiled(num_file,bi,bj)) then
           ENDDO !do sample_num=..
                
          ENDDO !bj
         ENDDO !bi

C Combine all processes
         DO ii=1,NSAMPLESMAX
          tmpgs = samples_buff(ii)
          _GLOBAL_SUM_RL (tmpgs, myThid) 
          samples_modval_glob(ii) = tmpgs
          tmpgs = samples_mask_buff(ii)
          _GLOBAL_SUM_RL (tmpgs, myThid)
          samples_mask_glob(ii) = tmpgs
         ENDDO     

         IF ( myProcId .eq. 0 ) THEN

C Loop over obs
         DO obs_num=1,NOBSMAX_IHOP
          IF (obs_num.LE.ObsNo(num_file)) THEN

           obs_modval=0.
           sample_mask_sum=0
           num_valid_samples=0

C Calculate model equi of each obs by averaging NP samples               
           DO sample_num=1,NSAMP_PER_OBS_MAX
            IF (sample_num.LE.ihopObs_np(num_file,obs_num)) THEN

              irec = ihopObs_sample1_ind(num_file,obs_num)
     &                + sample_num-1
               
              obs_modval = obs_modval + samples_modval_glob(irec)
     &               *samples_mask_glob(irec)
              sample_mask_sum = sample_mask_sum
     &               +samples_mask_glob(irec)
              IF (samples_mask_glob(irec) .GT. 0. _d 0) THEN
               num_valid_samples = num_valid_samples+1
              ENDIF

            ENDIF
           ENDDO
           
!IE: ignoring any time averaging for now
!C Time averaging
!           IF (ihopOperation(num_file).eq.1) THEN
!            obs_modval = obs_modval/ihopObs_delT(num_file, obs_num)
!           ENDIF
 
C Spatial averaging
           IF (num_valid_samples.gt.0) THEN
             obs_modval = obs_modval/sample_mask_sum
             ihopObs_modmask = 1
           ELSE
             obs_modval = -9999.
             ihopObs_modmask = 0
           ENDIF

c Write to global netcdf file
              CALL active_write_ihop_glob(num_file,
     &         obs_modval,
     &         obs_num,optimcycle,myThid,
     &         ihop_globaldummy(num_file))

          ENDIF !if (obs_num.LE.ObsNo(num_file)) then
         ENDDO !do obs_num=..

cc         IF ( myProcId .eq. 0 ) THEN
           err = NF_CLOSE(fidglobal(num_file))
           IL  = ILNBLNK( ihopobsfiles(num_file) )
           write(ihopfile(1:80),'(1a)')
     &      ihopobsfiles(num_file)(1:IL)
           IL  = ILNBLNK( ihopfile )
           JL  = ILNBLNK( ihopObsDir )
           write(fnameequinc(1:80),'(3a)')
     &      ihopObsDir(1:JL),ihopfile(1:IL),'.equi.nc'
           err = NF_OPEN(fnameequinc,NF_NOWRITE,fidglobal(num_file))
cc         ENDIF
         
C Loop over obs
         DO obs_num=1,NOBSMAX_IHOP
          IF (obs_num.LE.ObsNo(num_file)) THEN

           obs_modval=0.
           obs_data=0.
           obs_uncert=0.

C Read observation and uncertainty
           CALL ihop_cost_read_obs(num_file,1,
     &           ihopObs_ind_glob(num_file,obs_num),
     &           obs_data,myThid)

           CALL ihop_cost_read_obs(num_file,-1,
     &           ihopObs_ind_glob(num_file,obs_num),
     &           obs_uncert,myThid)

C Read model equivalent from global file
cc           IF ( myProcId .eq. 0 ) THEN
              CALL active_read_ihop_glob(num_file,
     &         obs_modval,
     &         obs_num,.false.,optimcycle,myThid,
     &         ihop_globaldummy(num_file))
cc           ENDIF

           IF (obs_uncert.GT.0.
     &            ) THEN
                 obs_weight=1/(obs_uncert*obs_uncert)
                 objf_ihop(num_file,bi,bj)=
     &             objf_ihop(num_file,bi,bj)
     &             +obs_weight*ihopObs_modmask
     &             *(obs_modval-obs_data)
     &             *(obs_modval-obs_data)
                 num_ihop(num_file,num_var)=
     &               num_ihop(num_file,num_var)
     &               +ihopObs_modmask
           ENDIF

          ENDIF !if (obs_num.LE.ObsNo(num_file)) then
         ENDDO !do obs_num...

         ENDIF !if myprocid = 0

       ENDDO !do num_var...
      ENDDO !do num_file=1,NFILESMAX_IHOP
  
#endif /* ALLOW_COST */
        _END_MASTER( myThid )

#ifdef ALLOW_COST
c Print cost function values
      DO num_file=1,NFILESMAX_IHOP
       DO num_var=1,NOBSMAX_IHOP
c
          objf_ihop_glo =
     &             objf_ihop(num_file,bi,bj)
          num_ihop_glo =
     &             num_ihop(num_file,num_var)
c
      write(msgbuf,'(A,I2,A,I2,A,2D12.5)')
     &  ' cost_ihop(',num_file,',',num_var,')= ',
     &  objf_ihop_glo,num_ihop_glo

      IF ( num_ihop_glo .GT. 0. ) CALL print_message( msgbuf,
     &  standardmessageunit,SQUEEZE_RIGHT, myThid)
c
       ENDDO
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
