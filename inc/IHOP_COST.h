!     *==========================================================*
!     | IHOP_COST.h                                              |
!     | o ihop cost terms.                                       |
!     *==========================================================*

!    IHOP cost flag
!     cost_ihop_flag  :: cost ihop flag (see ihop_cost_test.F)

      INTEGER cost_ihop_flag
      COMMON /ihop_cost_i/ cost_ihop_flag
      
!    IHOP cost
!     objf_ihop     :: ihop travel times
!     num_ihop      :: number of observations 
!     mult_ihop     :: multiplier applied to all cost terms

      _RL  objf_ihop (nSx,nSy)
      _RL  num_ihop  (nSx,nSy)
      _RL  mult_ihop (nSx,nSy)
      COMMON /IHOP_COST__R/                                                                                                         &
     &                objf_ihop,                                                                                                    &
     &                num_ihop,                                                                                                     &
     &                mult_ihop

!    IHOP cost filenames
!     ihopObsDir    :: directory where ihop observations are found
!     ihopObsFile   :: file name for ihop observations 

      CHARACTER*(MAX_LEN_FNAM) ihopObsDir
      CHARACTER*(MAX_LEN_FNAM) ihopObsFile
      COMMON /IHOP_COST_C/                                                                                                          &
     &                ihopObsDir, ihopObsFile

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
