!     *==========================================================*
!     | IHOP_COST.h                                              |
!     | o ihop cost terms.                                       |
!     *==========================================================*

!    ihop logical parameters
      LOGICAL ihopDoNcOutput

      COMMON /IHOP_COST_L/ ihopDoNcOutput

!    IHOP cost integer parameters
!     ObsNo* :: No. of observations in a single ihop obs file, iOBS
!     ihopObs_ind_glob* :: MITgcm global index of each obs datum

      INTEGER ObsNo(NFILESMAX_IHOP)
      INTEGER ObsNo_tiled(NFILESMAX_IHOP,nsx,nsy)
      INTEGER ihopObs_ind_glob(NFILESMAX_IHOP,NOBSMAX_IHOP)
      INTEGER ihopObs_ind_glob_tiled(                                                                                               &
     &                        NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      INTEGER ncidFWD(NFILESMAX_IHOP,nsx,nsy)
      INTEGER ncidTL(NFILESMAX_IHOP,nsx,nsy)
      INTEGER ncidAD(NFILESMAX_IHOP,nsx,nsy)
      INTEGER ncidData(NFILESMAX_IHOP)
      INTEGER ncidGLOB(NFILESMAX_IHOP)
      INTEGER ncidTLGLOB(NFILESMAX_IHOP)
      INTEGER ncidADGLOB(NFILESMAX_IHOP)
      INTEGER ihopObs_i_tiled(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      INTEGER ihopObs_j_tiled(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      INTEGER ihopObs_k_tiled(NFILESMAX_IHOP,NOBSMAX_IHOP,nsx,nsy)
      INTEGER ihopObs_sample1_ind(NFILESMAX_IHOP,NOBSMAX_IHOP)
      COMMON /IHOP_COST_I/                                                                                                          &
     &                  ObsNo, ObsNo_tiled,                                                                                         &
     &                  ihopObs_ind_glob, ihopObs_ind_glob_tiled,                                                                   &
     &                  ncidFWD, ncidAD, ncidTL,                                                                                    &
     &                  ncidGLOB, ncidADGLOB, ncidTLGLOB,                                                                           &
     &                  ncidData,                                                                                                   &
     &                  ihopObs_i_tiled, ihopObs_j_tiled,                                                                           &
     &                  ihopObs_k_tiled,                                                                                            &
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
     &                ihopObs_time, ihopObs_lat, ihopObs_lon,                                                                       &
     &                ihopObs_depth,                                                                                                &
     &                ihopObs_uncert,                                                                                               &
     &                ihopObs_modmask, ihopObs_modmask_tiled

!    IHOP cost filenames
!     IHOPObs_Dir   :: directory where ihop observations are found
!     IHOPObs_Files :: file name for ihop observations

      CHARACTER*(MAX_LEN_FNAM) IHOPObs_Dir
      CHARACTER*(MAX_LEN_FNAM) IHOPObs_Files(NFILESMAX_IHOP)
      CHARACTER*(8)  ihop_nameval
      CHARACTER*(8)  ihop_namemask
      CHARACTER*(11) ihop_nameuncert
      CHARACTER*(7)  ihop_nameequi
      COMMON /IHOP_COST_C/                                                                                                          &
     &                IHOPObs_Dir, IHOPObs_Files,                                                                                   &
     &                ihop_nameval, ihop_namemask,                                                                                  &
     &                ihop_nameuncert, ihop_nameequi

      _RL ihop_dummy(NFILESMAX_IHOP,nsx,nsy)
      _RL ihop_globaldummy(NFILESMAX_IHOP)
      COMMON /IHOP_CTRL_DUMMY/
     &                ihop_dummy,
     &                ihop_globaldummy

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
