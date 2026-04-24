!     *==========================================================*
!     | BELLI_COST.h                                              |
!     | o belli cost terms.                                       |
!     *==========================================================*

!    belli logical parameters
      LOGICAL belliDoNcOutput

      COMMON /BELLI_COST_L/ belliDoNcOutput

!    BELLI cost integer parameters
!     bObsNo* :: No. of observations in a single belli obs file, iOBS
!     bObs_ind_glob* :: MITgcm global index of each obs datum

      INTEGER bObsNo(NFILESMAX_BELLI)
      INTEGER bObsNo_tiled(NFILESMAX_BELLI,nsx,nsy)
      INTEGER bObs_ind_glob(NFILESMAX_BELLI,NOBSMAX_BELLI)
      INTEGER bObs_ind_glob_tiled(                                                                                                  &
     &                        NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER ncidFWD(NFILESMAX_BELLI,nsx,nsy)
      INTEGER ncidTL(NFILESMAX_BELLI,nsx,nsy)
      INTEGER ncidAD(NFILESMAX_BELLI,nsx,nsy)
      INTEGER ncidData(NFILESMAX_BELLI)
      INTEGER ncidGLOB(NFILESMAX_BELLI)
      INTEGER ncidTLGLOB(NFILESMAX_BELLI)
      INTEGER ncidADGLOB(NFILESMAX_BELLI)
      INTEGER bObs_i_tiled(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER bObs_j_tiled(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER bObs_k_tiled(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER bObs_sample1_ind(NFILESMAX_BELLI,NOBSMAX_BELLI)
      COMMON /BELLI_COST_I/                                                                                                         &
     &                  ObsNo, ObsNo_tiled,                                                                                         &
     &                  bObs_ind_glob, bObs_ind_glob_tiled,                                                                         &
     &                  ncidFWD, ncidAD, ncidTL,                                                                                    &
     &                  ncidGLOB, ncidADGLOB, ncidTLGLOB,                                                                           &
     &                  ncidData,                                                                                                   &
     &                  bObs_i_tiled, bObs_j_tiled,                                                                                 &
     &                  bObs_k_tiled,                                                                                               &
     &                  bObs_sample1_ind

!    BELLI buffers
      _RL belli_data_buff(1000)
      _RL belli_uncert_buff(1000)
      INTEGER belli_minind_buff
      INTEGER belli_maxind_buff
      INTEGER belli_curfile_buff

      COMMON /BELLI_BUFF_R/ belli_data_buff, belli_uncert_buff
      COMMON /BELLI_BUFF_I/                                                                                                         &
     & belli_minind_buff, belli_maxind_buff, belli_curfile_buff

!    BELLI cost real parameters
!     objf_belli     :: belli travel times
!     num_belli      :: number of observations
!     mult_belli     :: multiplier applied to all cost terms
!     bObs_time  :: obs. start time

      _RL  objf_belli (NFILESMAX_BELLI)
      _RL  num_belli  (NFILESMAX_BELLI)
      _RL  mult_belli (NFILESMAX_BELLI)
      _RL  bObs_time(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  bObs_lon(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  bObs_lat(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  bObs_depth(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  bObs_uncert(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  bObs_modmask
      _RL  bObs_modmask_tiled(nsx,nsy)
      _RL geninfluence(35604)
      COMMON /BELLI_COST_R/                                                                                                         &
     &                objf_belli,                                                                                                   &
     &                num_belli,                                                                                                    &
     &                mult_belli,                                                                                                   &
     &                bObs_time, bObs_lat, bObs_lon,                                                                                &
     &                bObs_depth,                                                                                                   &
     &                bObs_uncert,                                                                                                  &
     &                bObs_modmask, bObs_modmask_tiled,                                                                             &
     &                geninfluence

!    BELLI cost filenames
!     bObs_Dir   :: directory where belli observations are found
!     bObs_Files :: file name for belli observations

      CHARACTER*(MAX_LEN_FNAM) bObs_Dir
      CHARACTER*(MAX_LEN_FNAM) bObs_Files(NFILESMAX_BELLI)
      CHARACTER*(8)  belli_nameval
      CHARACTER*(8)  belli_namemask
      CHARACTER*(11) belli_nameuncert
      CHARACTER*(7)  belli_nameequi
      COMMON /BELLI_COST_C/                                                                                                         &
     &                bObs_Dir, bObs_Files,                                                                                         &
     &                belli_nameval, belli_namemask,                                                                                &
     &                belli_nameuncert, belli_nameequi

      _RL belli_dummy(NFILESMAX_BELLI,nsx,nsy)
      _RL belli_globaldummy(NFILESMAX_BELLI)
      COMMON /BELLI_CTRL_DUMMY/                                                                                                     &
     &                belli_dummy,                                                                                                  &
     &                belli_globaldummy

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
