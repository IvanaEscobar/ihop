!     *==========================================================*
!     | BELLI_COST.h                                              |
!     | o belli cost terms.                                       |
!     *==========================================================*

!    belli logical parameters
      LOGICAL belliDoNcOutput

      COMMON /BELLI_COST_L/ belliDoNcOutput

!    BELLI cost integer parameters
!     belliObsNo* :: No. of observations in a single belli obs file, iOBS
!     belliObs_ind_glob* :: MITgcm global index of each obs datum

      INTEGER belliObsNo(NFILESMAX_BELLI)
      INTEGER belliObsNo_tiled(NFILESMAX_BELLI,nsx,nsy)
      INTEGER belliObs_ind_glob(NFILESMAX_BELLI,NOBSMAX_BELLI)
      INTEGER belliObs_ind_glob_tiled(                                                                                                  &
     &                        NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER ncidFWD(NFILESMAX_BELLI,nsx,nsy)
      INTEGER ncidTL(NFILESMAX_BELLI,nsx,nsy)
      INTEGER ncidAD(NFILESMAX_BELLI,nsx,nsy)
      INTEGER ncidData(NFILESMAX_BELLI)
      INTEGER ncidGLOB(NFILESMAX_BELLI)
      INTEGER ncidTLGLOB(NFILESMAX_BELLI)
      INTEGER ncidADGLOB(NFILESMAX_BELLI)
      INTEGER belliObs_i_tiled(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER belliObs_j_tiled(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER belliObs_k_tiled(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      INTEGER belliObs_sample1_ind(NFILESMAX_BELLI,NOBSMAX_BELLI)
      COMMON /BELLI_COST_I/                                                                                                         &
     &                  belliObsNo, belliObsNo_tiled,                                                                                       &
     &                  belliObs_ind_glob, belliObs_ind_glob_tiled,                                                                         &
     &                  ncidFWD, ncidAD, ncidTL,                                                                                    &
     &                  ncidGLOB, ncidADGLOB, ncidTLGLOB,                                                                           &
     &                  ncidData,                                                                                                   &
     &                  belliObs_i_tiled, belliObs_j_tiled,                                                                                 &
     &                  belliObs_k_tiled,                                                                                               &
     &                  belliObs_sample1_ind

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
!     belliObs_time  :: obs. start time

      _RL  objf_belli (NFILESMAX_BELLI)
      _RL  num_belli  (NFILESMAX_BELLI)
      _RL  mult_belli (NFILESMAX_BELLI)
      _RL  belliObs_time(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  belliObs_lon(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  belliObs_lat(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  belliObs_depth(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  belliObs_uncert(NFILESMAX_BELLI,NOBSMAX_BELLI,nsx,nsy)
      _RL  belliObs_modmask
      _RL  belliObs_modmask_tiled(nsx,nsy)
      _RL geninfluence(35604)
      COMMON /BELLI_COST_R/                                                                                                         &
     &                objf_belli,                                                                                                   &
     &                num_belli,                                                                                                    &
     &                mult_belli,                                                                                                   &
     &                belliObs_time, belliObs_lat, belliObs_lon,                                                                                &
     &                belliObs_depth,                                                                                                   &
     &                belliObs_uncert,                                                                                                  &
     &                belliObs_modmask, belliObs_modmask_tiled,                                                                             &
     &                geninfluence

!    BELLI cost filenames
!     belliObs_Dir   :: directory where belli observations are found
!     belliObs_Files :: file name for belli observations

      CHARACTER*(MAX_LEN_FNAM) belliObs_Dir
      CHARACTER*(MAX_LEN_FNAM) belliObs_Files(NFILESMAX_BELLI)
      CHARACTER*(8)  belli_nameval
      CHARACTER*(8)  belli_namemask
      CHARACTER*(11) belli_nameuncert
      CHARACTER*(7)  belli_nameequi
      COMMON /BELLI_COST_C/                                                                                                         &
     &                belliObs_Dir, belliObs_Files,                                                                                         &
     &                belli_nameval, belli_namemask,                                                                                &
     &                belli_nameuncert, belli_nameequi

      _RL belli_dummy(NFILESMAX_BELLI,nsx,nsy)
      _RL belli_globaldummy(NFILESMAX_BELLI)
      COMMON /BELLI_CTRL_DUMMY/                                                                                                     &
     &                belli_dummy,                                                                                                  &
     &                belli_globaldummy

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
