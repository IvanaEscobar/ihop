!     *==========================================================*
!     | IHOP_COST.h                                              |
!     | o ihop cost terms.                                       |
!     *==========================================================*

!    ihop logical parameters
      LOGICAL ihopDoNcOutput

      COMMON /ihop_cost_l/ ihopDoNcOutput

!    IHOP cost integer parameters
!     ObsNo* :: No. of observations in a single ihop obs file, iOBS
!     ihopObs_ind_glob* :: MITgcm global index of each obs datum

      INTEGER ObsNo(NFILESMAX_IHOP)
      INTEGER ObsNo_tiled(NFILESMAX_IHOP,nsx,nsy)
      INTEGER ihopObs_ind_glob(NFILESMAX_IHOP,NOBSMAX_IHOP)
      INTEGER ihopObs_ind_glob_tiled(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      INTEGER fidfwd_obs(NFILESMAX_IHOP,nsx,nsy)
      INTEGER fidftl_obs(NFILESMAX_IHOP,nsx,nsy)
      INTEGER fidadj_obs(NFILESMAX_IHOP,nsx,nsy)
      INTEGER fiddata_obs(NFILESMAX_IHOP)
      INTEGER fidglobal(NFILESMAX_IHOP)
      INTEGER fidftlglobal(NFILESMAX_IHOP)
      INTEGER fidadglobal(NFILESMAX_IHOP)
      INTEGER ihopObs_i_tiled(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      INTEGER ihopObs_j_tiled(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      INTEGER ihopObs_k_tiled(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      INTEGER ihopObs_sample1_ind(NFILESMAX_IHOP,NOBSMAX_IHOP)
      COMMON /ihop_cost_i/                                                                                                          &
     &                  ObsNo, ObsNo_tiled,                                                                                         &
     &                  ihopObs_ind_glob, ihopObs_ind_glob_tiled,                                                                   &
     &                  fidfwd_obs, fidadj_obs, fidftl_obs,                                                                         &
     &                  fidglobal, fidadglobal, fidftlglobal,                                                                       &
     &                  fiddata_obs,                                                                                                &
     &                  ihopObs_i_tiled, ihopObs_j_tiled, ihopObs_k_tiled,                                                          &
     &                  ihopObs_sample1_ind

!    IHOP buffers
      _RL ihop_data_buff(1000)
      _RL ihop_uncert_buff(1000)
      INTEGER ihop_minind_buff
      INTEGER ihop_maxind_buff
      INTEGER ihop_curfile_buff

      COMMON /IHOP_BUFF_R/ ihop_data_buff, ihop_uncert_buff
      COMMON /IHOP_BUFF_I/                                                                                                          &
     & ihop_minind_buff, ihop_maxind_buff, ihop_curfile_buff

!    IHOP cost real parameters
!     objf_ihop     :: ihop travel times
!     num_ihop      :: number of observations
!     mult_ihop     :: multiplier applied to all cost terms
!     ihopObs_time  :: obs. start time

      _RL  objf_ihop (NFILESMAX_IHOP)
      _RL  num_ihop  (NFILESMAX_IHOP)
      _RL  mult_ihop (NFILESMAX_IHOP)
      _RL  ihopObs_time(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      _RL  ihopObs_lon(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      _RL  ihopObs_lat(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      _RL  ihopObs_depth(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      _RL  ihopObs_uncert(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      _RL  ihopObs_modmask
      _RL  ihopObs_modmask_tiled(nsx,nsy)
      COMMON /IHOP_COST_R/                                                                                                          &
     &                objf_ihop,                                                                                                    &
     &                num_ihop,                                                                                                     &
     &                mult_ihop,                                                                                                    &
     &                ihopObs_time, ihopObs_lat, ihopObs_lon, ihopObs_depth,                                                        &
     &                ihopObs_uncert,                                                                                               &
     &                ihopObs_modmask, ihopObs_modmask_tiled

!    IHOP cost filenames
!     ihopObsDir    :: directory where ihop observations are found
!     ihopObsFiles  :: file name for ihop observations

      CHARACTER*(MAX_LEN_FNAM) ihopObsDir
      CHARACTER*(MAX_LEN_FNAM) ihopObsFiles(NFILESMAX_IHOP)
      CHARACTER*(8)  ihop_nameval
      CHARACTER*(8)  ihop_namemask
      CHARACTER*(11) ihop_nameuncert
      CHARACTER*(7)  ihop_nameequi
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
