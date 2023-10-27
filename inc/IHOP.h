#ifdef ALLOW_IHOP
!BOP
!     !ROUTINE: IHOP.h
!     !INTERFACE:
!     #include IHOP.h
!
!     !DESCRIPTION:
!     *================================================================*
!     | IHOP.h
!     | o Header file defining "ihop" parameters and variables
!     *================================================================*
!EOP

!     Package flag
      LOGICAL IHOP_MNC
      LOGICAL IHOP_MDSIO

      COMMON /IHOP_PACKAGE/                                                                                                             &
     &                      IHOP_MNC, IHOP_MDSIO

!     IHOP parameters
!     ===============
!--   COMMON /IHOP_PARAMS_L/ IHOP logical-type parameters:
!     IHOP_bellOn   :: true if bellhop driver needs to run (UNUSED)
!     useSSPFile    :: true if *.ssp is used instead MITgcm SSP

      LOGICAL IHOP_bellOn
      LOGICAL useSSPFile

      COMMON /IHOP_PARAMS_L/                                                                                                             &
     &      IHOP_bellOn, useSSPFile

!-- COMMON /IHOP_PARAMS_C/ IHOP Character-type parameters:
!   IHOP_fileroot   :: File name for reading in an environment
!   IHOP_title      :: Title name for writing into output files 
!   IHOP_interpfile :: File name for reading NetCDF inputs
!   IHOP_topopt     :: SSP interpolation, top boundary type
!   IHOP_botopt     :: bottom boundary type
!   IHOP_runopt     :: run type (R/E/A)

      CHARACTER*(MAX_LEN_FNAM) IHOP_fileroot
      CHARACTER*(MAX_LEN_FNAM) IHOP_title
      CHARACTER*(MAX_LEN_FNAM) IHOP_interpfile
      CHARACTER*(6) IHOP_topopt
      CHARACTER*(2) IHOP_botopt
      CHARACTER*(7) IHOP_runopt

      COMMON /IHOP_PARAMS_C/                                                                                                             &
     &      IHOP_fileroot, IHOP_title,                                                                                                   &
     &      IHOP_topopt, IHOP_botopt, IHOP_runopt,                                                                                       &
     &      IHOP_interpfile

!-- COMMON /IHOP_PARAMS_I/ IHOP Integer-type parameters:
!   IHOP_nalpha :: No. of rays to propagate
!   IHOP_nsd    :: No. of source depths (m)
!   IHOP_nrd    :: No. of receiver depths (m)
!   IHOP_nrr    :: No. of receiver ranges (km)

      INTEGER IHOP_nalpha
      INTEGER IHOP_nsd
      INTEGER IHOP_nrd
      INTEGER IHOP_nrr
      INTEGER ihop_iter
      INTEGER IHOP_npts_range
      INTEGER IHOP_npts_idw

      COMMON /IHOP_PARAMS_I/                                                                                                            &
     &      IHOP_nsd,                                                                                                                   &
     &      IHOP_nrd, IHOP_nrr,                                                                                                         &
     &      IHOP_npts_range, IHOP_npts_idw,                                                                                             &
     &      IHOP_nalpha, ihop_iter

!-- COMMON /IHOP_PARAMS_R/ IHOP Real-type parameters:
!   IHOP_freq           :: frequency (Hz)
!   IHOP_depth          :: depth of bottom (m)
!   IHOP_bcsound        :: bottom sound speed (m/s) 
!   IHOP_bcsoundshear   :: shear bottom sound speed (m/s) 
!   IHOP_bcsoundI       :: IMAG bottom sound speed (m/s) 
!   IHOP_bcsoundshearI  :: IMAG shear bottom sound speed (m/s) 
!   IHOP_brho           :: bottom density (kg/m^3)
!   IHOP_sd             :: source depth (m)
!   IHOP_rd             :: receiver depth (m)
!   IHOP_rr             :: receiver ranges (km)
!   IHOP_alpha          :: bearing launch angles (degrees)
!   IHOP_step           :: step length (m)
!   IHOP_zbox           :: acoustic domain depth (m)
!   IHOP_rbox           :: acoustic domain range (km)

      _RL IHOP_freq
      _RL IHOP_depth
      _RL IHOP_bcsound
      _RL IHOP_bcsoundshear 
      _RL IHOP_bcsoundI
      _RL IHOP_bcsoundshearI 
      _RL IHOP_brho
      _RL IHOP_sd (nsd)
      _RL IHOP_rd (nrd)
      _RL IHOP_rr (nrr)
      _RL IHOP_alpha (2)
      _RL IHOP_step
      _RL IHOP_zbox
      _RL IHOP_rbox
      _RL ihop_v_weight ( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      _RS ihop_v_xc ( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      _RS ihop_v_yc ( IHOP_MAX_NC_SIZE*IHOP_MAX_NC_SIZE )
      _RL ihop_idw_weights ( IHOP_MAX_NC_SIZE, IHOP_MAX_NC_SIZE )
      _RS ihop_xc ( IHOP_MAX_NC_SIZE, IHOP_MAX_NC_SIZE )
      _RS ihop_yc ( IHOP_MAX_NC_SIZE, IHOP_MAX_NC_SIZE )
      _RL ihop_ranges ( IHOP_MAX_NC_SIZE )

      COMMON /IHOP_PARAMS_R/                                                                                                            &
     &      IHOP_freq, IHOP_depth, IHOP_bcsound, IHOP_bcsoundshear,                                                                     &
     &      ihop_brho, IHOP_bcsoundI, IHOP_bcsoundshearI,                                                                               &
     &      IHOP_sd, IHOP_rd, IHOP_rr, IHOP_alpha, IHOP_step,                                                                           &
     &      ihop_yc, ihop_xc, ihop_idw_weights, ihop_ranges,                                                                            &
     &      ihop_v_yc, ihop_v_xc, ihop_v_weight, IHOP_zbox, IHOP_rbox


#ifdef IHOP_3D_STATE
!C     IHOP 3-dim. fields
      _RL ihop_ssp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      COMMON /IHOP_STATE_3D/                                                                                                             &
     &    ihop_ssp 
#endif /* IHOP_3D_STATE */

#ifdef IHOP_2D_STATE
#endif /* IHOP_2D_STATE */

#ifdef IHOP_TENDENCY
#endif /* IHOP_TENDENCY */

#endif /* ALLOW_IHOP */

!EH3 ;;; Local Variables: ***
!EH3 ;;; mode:fortran ***
!EH3 ;;; End: ***
