#ifdef ALLOW_BELLI
!BOP
!     !ROUTINE: BELLI.h
!     !INTERFACE:
!     #include BELLI.h
!
!     !DESCRIPTION:
!     *================================================================*
!     | BELLI.h
!     | o Header file defining "belli" parameters and variables
!     *================================================================*
!EOP

!     Package flag
      LOGICAL BELLI_MNC
      LOGICAL BELLI_MDSIO

      COMMON /BELLI_PACKAGE/                                                                                                        &
     &                      BELLI_MNC, BELLI_MDSIO

!     BELLI parameters
!     ===============
!--   COMMON /BELLI_PARAMS_L/ BELLI logical-type parameters:
!     writeDelay    :: true if delay is a desired output
!     useSSPFile    :: true if *.ssp is used instead MITgcm SSP

      LOGICAL writeDelay
      LOGICAL useSSPFile

      COMMON /BELLI_PARAMS_L/                                                                                                       &
     &      writeDelay, useSSPFile

!-- COMMON /BELLI_PARAMS_C/ BELLI Character-type parameters:
!   BELLI_fileroot   :: File name for reading in an environment
!   BELLI_title      :: Title name for writing into output files
!   BELLI_interpfile :: File name for reading NetCDF inputs
!   BELLI_topopt     :: SSP interpolation, top boundary type
!   BELLI_botopt     :: bottom boundary type
!   BELLI_runopt     :: run type [R/E/A]

      CHARACTER*(MAX_LEN_FNAM) BELLI_fileroot
      CHARACTER*(MAX_LEN_FNAM) BELLI_title
      CHARACTER*(MAX_LEN_FNAM) BELLI_interpfile
      CHARACTER*(6) BELLI_topopt
      CHARACTER*(2) BELLI_botopt
      CHARACTER*(7) BELLI_runopt

      COMMON /BELLI_PARAMS_C/                                                                                                       &
     &      BELLI_fileroot, BELLI_title,                                                                                            &
     &      BELLI_topopt, BELLI_botopt, BELLI_runopt,                                                                               &
     &      BELLI_interpfile

!-- COMMON /BELLI_PARAMS_I/ BELLI Integer-type parameters:
!   BELLI_nalpha :: No. of rays to propagate
!   BELLI_nts    :: No. of sample times
!   BELLI_nsd    :: No. of source depths [m]
!   BELLI_nrd    :: No. of receiver depths [m]
!   BELLI_nrr    :: No. of receiver ranges [km]
!   BELLI_iter   :: GCM iteration to run belli

      INTEGER BELLI_nalpha
      INTEGER BELLI_nts
      INTEGER BELLI_nsd
      INTEGER BELLI_nrd
      INTEGER BELLI_nrr
      INTEGER belli_iter(nts)
      INTEGER BELLI_npts_range
      INTEGER BELLI_npts_idw

      COMMON /BELLI_PARAMS_I/                                                                                                       &
     &      BELLI_nts, BELLI_nsd,                                                                                                   &
     &      BELLI_nrd, BELLI_nrr,                                                                                                   &
     &      BELLI_npts_range, BELLI_npts_idw,                                                                                       &
     &      BELLI_nalpha, belli_iter

!-- COMMON /BELLI_PARAMS_R/ BELLI Real-type parameters:
!   BELLI_dumpfreq       :: frequency of output dump to run directory
!   BELLI_freq           :: frequency [Hz]
!   BELLI_depth          :: depth of bottom [m]
!   BELLI_bcsound        :: bottom sound speed [m/s]
!   BELLI_bcsoundshear   :: shear bottom sound speed [m/s]
!   BELLI_bcsoundI       :: IMAG bottom sound speed [m/s]
!   BELLI_bcsoundshearI  :: IMAG shear bottom sound speed [m/s]
!   BELLI_brho           :: bottom density [kg/m^3]
!   BELLI_sd             :: source depth [m]
!   BELLI_rd             :: receiver depth [m]
!   BELLI_rr             :: receiver ranges [km]
!   BELLI_alpha          :: bearing launch angles [degrees]
!   BELLI_step           :: step length [m]

      _RL BELLI_dumpfreq
      _RL BELLI_freq
      _RL BELLI_depth
      _RL BELLI_bcsound
      _RL BELLI_bcsoundshear
      _RL BELLI_bcsoundI
      _RL BELLI_bcsoundshearI
      _RL BELLI_brho
      _RL BELLI_sd (nsd)
      _RL BELLI_rd (nrd)
      _RL BELLI_rr (nrr)
      _RL BELLI_alpha (2)
      _RL BELLI_step
      _RL belli_idw_weights ( BELLI_MAX_RANGE, BELLI_MAX_NC_SIZE )
      _RS belli_xc ( BELLI_MAX_RANGE, BELLI_MAX_NC_SIZE )
      _RS belli_yc ( BELLI_MAX_RANGE, BELLI_MAX_NC_SIZE )
      _RL belli_ranges ( BELLI_MAX_RANGE )
      _RL belli_sumweights ( BELLI_MAX_RANGE, BELLI_MAX_NC_SIZE )

      COMMON /BELLI_PARAMS_R/                                                                                                       &
     &      BELLI_dumpfreq,                                                                                                         &
     &      BELLI_freq, BELLI_depth, BELLI_bcsound, BELLI_bcsoundshear,                                                             &
     &      belli_brho, BELLI_bcsoundI, BELLI_bcsoundshearI,                                                                        &
     &      BELLI_sd, BELLI_rd, BELLI_rr, BELLI_alpha, BELLI_step,                                                                  &
     &      belli_yc, belli_xc, belli_idw_weights, belli_ranges,                                                                    & 
     &      belli_sumweights


#ifdef BELLI_3D_STATE
!C     BELLI 3-dim. fields
!   BELLI_ssp            :: speed of sound [m/s]
!                          (for diagnostic + belli AD model)
      _RL belli_ssp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      COMMON /BELLI_STATE_3D/                                                                                                       &
     &    belli_ssp
#endif /* BELLI_3D_STATE */

#ifdef BELLI_2D_STATE
!C     BELLI 2-dim. fields
!   BELLI_sld            :: sonic layer depth [m]
!                          (for diagnostic)
      _RL belli_sld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /BELLI_STATE_2D/                                                                                                       &
     &    belli_sld
#endif /* BELLI_2D_STATE */

#ifdef BELLI_TENDENCY
#endif /* BELLI_TENDENCY */

#endif /* ALLOW_BELLI */

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
