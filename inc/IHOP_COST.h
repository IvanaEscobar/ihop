!     *==========================================================*
!     | IHOP_COST.h                                              |
!     | o ihop cost terms.                                       |
!     *==========================================================*

!    ihop logical parameters
      LOGICAL ihopDoNcOutput

      COMMON /ihop_cost_l/ ihopDoNcOutput

!    IHOP cost integer parameters 
!     cost_ihop_flag  :: cost ihop flag (see ihop_cost_test.F)

      INTEGER cost_ihop_flag
      INTEGER obs_ind_glob(NFILESMAX_IHOP,NOBSMAX_IHOP)
      INTEGER ihopOperation(NFILESMAX_IHOP)
      INTEGER sample_ind_glob(NFILESMAX_IHOP,NSAMPLESMAX,nsx,nsy)
      INTEGER obs_np(NFILESMAX_IHOP,NOBSMAX_IHOP)
      INTEGER ObsNo(NFILESMAX_IHOP)
      INTEGER sampleNo(NFILESMAX_IHOP,nsx,nsy)
      INTEGER fidfwd_obs(NFILESMAX_IHOP,nsx,nsy)
      INTEGER fidadj_obs(NFILESMAX_IHOP,nsx,nsy)
      INTEGER fiddata_obs(NFILESMAX_IHOP)
      INTEGER fidglobal(NFILESMAX_IHOP)
      INTEGER fidadglobal(NFILESMAX_IHOP)
      INTEGER obs_sample1_ind(NFILESMAX_IHOP,NOBSMAX_IHOP)
      COMMON /ihop_cost_i/                                                                                                          &
     &                  cost_ihop_flag, obs_ind_glob, ihopOperation, obs_np,                                                        &
     &                  fidfwd_obs, fidadj_obs, fidglobal, fidadglobal, 
     &                  fiddata_obs, sampleNo,                                                                                                   &
     &                  sample_ind_glob, ObsNo, obs_sample1_ind
      
!    IHOP buffers
      _RL ihop_data_buff(1000)
      _RL ihop_uncert_buff(1000)
      INTEGER ihop_minind_buff
      INTEGER ihop_maxind_buff
      INTEGER ihop_curfile_buff
      
      COMMON /IHOP_BUFF_R/ ihop_data_buff, ihop_uncert_buff
      COMMON /IHOP_BUFF_I/                                                                                                            &
     & ihop_minind_buff, ihop_maxind_buff, ihop_curfile_buff
     
!    IHOP cost real parameters
!     objf_ihop     :: ihop travel times
!     num_ihop      :: number of observations 
!     mult_ihop     :: multiplier applied to all cost terms

      _RL  sample_modmask(nsx,nsy)
      _RL  objf_ihop (nSx,nSy)
      _RL  num_ihop  (nSx,nSy)
      _RL  mult_ihop (NFILESMAX_IHOP)
      _RL obs_modmask
      _RL obs_delT(NFILESMAX_IHOP,NOBSMAX_IHOP)
      COMMON /IHOP_COST_R/                                                                                                          &
     &                sample_modmask,                                                                                               &
     &                objf_ihop,                                                                                                    &
     &                num_ihop,                                                                                                     &
     &                mult_ihop, obs_modmask, obs_delT

!    IHOP cost filenames
!     ihopObsDir    :: directory where ihop observations are found
!     ihopObsFile   :: file name for ihop observations 

      CHARACTER*(MAX_LEN_FNAM) ihopObsDir
      CHARACTER*(MAX_LEN_FNAM) ihopObsFiles(NFILESMAX_IHOP)
      CHARACTER*(8)  ihop_nameval
      CHARACTER*(12) ihop_namemask
      CHARACTER*(14) ihop_nameuncert
      CHARACTER*(8)  ihop_nameequi
      COMMON /IHOP_COST_C/                                                                                                          &
     &                ihopObsDir, ihopObsFiles,                                                                                     &
     &                ihop_nameval, ihop_namemask,                                                                                  &
     &                ihop_nameuncert, ihop_nameequi

      _RL ihop_dummy(NFILESMAX_IHOP,nsx,nsy)
      _RL ihop_globaldummy(NFILESMAX_IHOP)
      COMMON /IHOP_CTRL_DUMMY/
     &                ihop_dummy,
     &                ihop_globaldummy

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
