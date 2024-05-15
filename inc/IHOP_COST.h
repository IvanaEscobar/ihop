C     *==========================================================*
C     | IHOP_COST.h                                            |
C     | o ihop cost terms.                                    |
C     *==========================================================*

C     objf_ihop   :: ihop travel times
      COMMON /ihop_cost_objf/
     &                objf_ihop,
      _RL  objf_ihop        (nSx,nSy)

C-   IHOP_cutoff_area & _heff :: only used in pkg/ecco GET_EXCONC_DECONC S/R
      COMMON /ihop_cost_aux_r/
     &                num_ihop,
     &                mult_ihop,
     &                IHOP_cutoff_area,
      _RL  num_ihop  (nSx,nSy)
      _RL  mult_ihop
      _RL  IHOP_cutoff_area

      COMMON /ihop_cost_data_aux_i/
     &                           costIhopStart1,
     &                           costIhopStart2,
     &                           costIhopEnd1,
     &                           costIhopEnd2
      INTEGER costIhopStart1
      INTEGER costIhopStart2
      INTEGER costIhopEnd1
      INTEGER costIhopEnd2

      COMMON /ihop_cost_data_times_r/ costIhopStart, costIhopEnd
      _RL costIhopStart
      _RL costIhopEnd

C     cost_ihop_flag  :: cost_ihop flag (see ihop_cost_test.F)
      COMMON /ihop_cost_i/ cost_ihop_flag
      INTEGER cost_ihop_flag

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
