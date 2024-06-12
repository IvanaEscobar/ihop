!     *==========================================================*
!     | IHOP_COST.h                                              |
!     | o ihop cost terms.                                       |
!     *==========================================================*

!    OBSFIT logical parameters
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
      INTEGER fidglobal(NFILESMAX_IHOP)
      INTEGER obs_sample1_ind(NFILESMAX_IHOP,NOBSMAX_IHOP)
      COMMON /ihop_cost_i/                                                                                                          &
     &                  cost_ihop_flag, obs_ind_glob, ihopOperation, obs_np,                                                        &
     &                  fidfwd_obs, fidglobal, sampleNo,
     &                  sample_ind_glob, ObsNo, obs_sample1_ind
      
!    IHOP cost real parameters
!     objf_ihop     :: ihop travel times
!     num_ihop      :: number of observations 
!     mult_ihop     :: multiplier applied to all cost terms

      _RL  sample_modmask(nsx,nsy)
      _RL  objf_ihop (nSx,nSy)
      _RL  num_ihop  (nSx,nSy)
      _RL  mult_ihop (nSx,nSy)
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
      COMMON /IHOP_COST_C/                                                                                                          &
     &                ihopObsDir, ihopObsFiles

      _RL ihop_dummy(NFILESMAX_IHOP,nsx,nsy)
      _RL ihop_globaldummy(NFILESMAX_IHOP)
      COMMON /IHOP_CTRL_DUMMY/
     &                ihop_dummy,
     &                ihop_globaldummy

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
