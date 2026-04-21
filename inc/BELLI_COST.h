!     *==========================================================*
!     | BELLI_COST.h                                              |
!     | o belli cost terms.                                       |
!     *==========================================================*

!    belli logical parameters
      LOGICAL belliDoNcOutput

      COMMON /BELLI_COST_L/ belliDoNcOutput

!    BELLI cost integer parameters
!     ObsNo* :: No. of observations in a single belli obs file, iOBS
!     B_Obs_ind_glob* :: MITgcm global index of each obs datum

      INTEGER ObsNo(NFILESMAX_BELLI)
      INTEGER ObsNo_tiled(NFILESMAX_BELLI,nsx,nsy)
      INTEGER B_Obs_ind_glob(NFILESMAX_BELLI,NOBSMAX_BELLI)
      INTEGER B_Obs_ind_glob_tiled(                                                                                                 &
     &                        NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER ncidFWD(NFILESMAX_BELLI,nsx,nsy)
      INTEGER ncidTL(NFILESMAX_BELLI,nsx,nsy)
      INTEGER ncidAD(NFILESMAX_BELLI,nsx,nsy)
      INTEGER ncidData(NFILESMAX_BELLI)
      INTEGER ncidGLOB(NFILESMAX_BELLI)
      INTEGER ncidTLGLOB(NFILESMAX_BELLI)
      INTEGER ncidADGLOB(NFILESMAX_BELLI)
      INTEGER B_Obs_i_tiled(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER B_Obs_j_tiled(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER B_Obs_k_tiled(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER B_Obs_sample1_ind(NFILESMAX_BELLI,NOBSMAX_BELLI)
      COMMON /BELLI_COST_I/                                                                                                         &
     &                  ObsNo, ObsNo_tiled,                                                                                         &
     &                  B_Obs_ind_glob, B_Obs_ind_glob_tiled,                                                                       &
     &                  ncidFWD, ncidAD, ncidTL,                                                                                    &
     &                  ncidGLOB, ncidADGLOB, ncidTLGLOB,                                                                           &
     &                  ncidData,                                                                                                   &
     &                  B_Obs_i_tiled, B_Obs_j_tiled,                                                                               &
     &                  B_Obs_k_tiled,                                                                                              &
     &                  B_Obs_sample1_ind

!    BELLI buffers
      _RL belli_data_buff(1000)
      _RL belli_uncert_buff(1000)
      INTEGER belli_minind_buff
      INTEGER belli_maxind_buff
      INTEGER belli_curfile_buff

      COMMON /BELLI_BUFF_R/ belli_data_buff, belli_uncert_buff
      COMMON /BELLI_BUFF_I/                                                                                                          &
     & belli_minind_buff, belli_maxind_buff, belli_curfile_buff

!    BELLI cost real parameters
!     objf_belli     :: belli travel times
!     num_belli      :: number of observations
!     mult_belli     :: multiplier applied to all cost terms
!     B_Obs_time  :: obs. start time

      _RL  objf_belli (NFILESMAX_BELLI)
      _RL  num_belli  (NFILESMAX_BELLI)
      _RL  mult_belli (NFILESMAX_BELLI)
      _RL  B_Obs_time(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  B_Obs_lon(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  B_Obs_lat(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  B_Obs_depth(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  B_Obs_uncert(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  B_Obs_modmask
      _RL  B_Obs_modmask_tiled(nsx,nsy)
      _RL geninfluence(35604)
      COMMON /BELLI_COST_R/                                                                                                         &
     &                objf_belli,                                                                                                   &
     &                num_belli,                                                                                                    &
     &                mult_belli,                                                                                                   &
     &                B_Obs_time, B_Obs_lat, B_Obs_lon,                                                                             &
     &                B_Obs_depth,                                                                                                  &
     &                B_Obs_uncert,                                                                                                 &
     &                B_Obs_modmask, B_Obs_modmask_tiled,                                                                           &
     &                geninfluence

!    BELLI cost filenames
!     B_Obs_Dir   :: directory where belli observations are found
!     B_Obs_Files :: file name for belli observations

      CHARACTER*(MAX_LEN_FNAM) B_Obs_Dir
      CHARACTER*(MAX_LEN_FNAM) B_Obs_Files(NFILESMAX_BELLI)
      CHARACTER*(8)  belli_nameval
      CHARACTER*(8)  belli_namemask
      CHARACTER*(11) belli_nameuncert
      CHARACTER*(7)  belli_nameequi
      COMMON /BELLI_COST_C/                                                                                                         &
     &                B_Obs_Dir, B_Obs_Files,                                                                                       &
     &                belli_nameval, belli_namemask,                                                                                &
     &                belli_nameuncert, belli_nameequi

      _RL belli_dummy(NFILESMAX_BELLI,nsx,nsy)
      _RL belli_globaldummy(NFILESMAX_BELLI)
      COMMON /BELLI_CTRL_DUMMY/                                                                                                     &
     &                belli_dummy,                                                                                                  &
     &                belli_globaldummy

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
