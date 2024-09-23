!BOP
!     !ROUTINE: tamc.h
!     !INTERFACE:
!     #include "tamc.h"

!     !DESCRIPTION:
!     *================================================================*
!     | tamc.h
!     | o Header file defining parameters and variables for the use of
!     |   the Tangent Linear and Adjoint Model Compiler (TAMC)
!     |   or the Transformations in Fortran tool (TAF).
!     |
!     | started: Christian Eckert eckert@mit.edu  04-Feb-1999
!     | changed: Patrick Heimbach heimbach@mit.edu 06-Jun-2000
!     | cleanup: Martin Losch Martin.Losch@awi.de Nov-2022
!     *================================================================*
!EOP
#ifdef ALLOW_AUTODIFF_TAMC

!     TAMC checkpointing parameters:
!     ==============================
!
!     The checkpointing parameters have to be consistent with other model
!     parameters and variables. This has to be checked before the model is
!     run.
!

#ifdef ALLOW_TAMC_CHECKPOINTING

!     nchklev_1 :: length of inner loop (=size of storage in memory)
!     nchklev_2 :: length of second loop (stored on disk)
!     nchklev_3 :: length of outer loop of 3-level checkpointing
      INTEGER    nchklev_1
      PARAMETER( nchklev_1 =   2 )
      INTEGER    nchklev_2
      PARAMETER( nchklev_2 =   5 )
      INTEGER    nchklev_3
      PARAMETER( nchklev_3 =  50 )
#ifdef AUTODIFF_4_LEVEL_CHECKPOINT
!     nchklev_4 :: length of outer loop of 4-level checkpointing
      INTEGER    nchklev_4
      PARAMETER( nchklev_4 =   1 )
#endif

!--   Note always check for the correct sizes of the common blocks!
!     The product of the nchklev_X needs to be at least equal to
!     nTimeSteps.

#else /* ALLOW_TAMC_CHECKPOINTING undefined */

!     Without ALLOW_TAMC_CHECKPOINTING, nchklev_1 needs to be at least
!     equal to nTimeSteps. This (arbitrary) setting would accommodate a
!     short run (e.g., 10.d with deltaT=10.mn)
      INTEGER    nchklev_1
      PARAMETER( nchklev_1 = 1500 )

#endif /* ALLOW_TAMC_CHECKPOINTING */

!     TAMC keys:
!     ==========
!
!     The keys are used for storing and reading data of the reference
!     trajectory. Currently there is only one global key.
!     ikey_dynamics :: key for main time stepping loop

      COMMON /TAMC_KEYS_I/ ikey_dynamics
      INTEGER ikey_dynamics

!     isbyte :: precision of tapes (both memory and disk).
!               For smaller tapes replace 8 by 4.
      INTEGER    isbyte
      PARAMETER( isbyte    = 8 )

!     maxpass :: maximum number of (active + passive) tracers
!                Note: defined in PTRACERS_SIZE.h if compiling pkg/ptracers
#ifndef ALLOW_PTRACERS
      INTEGER    maxpass
      PARAMETER( maxpass   = 3 )
#endif
!     maxcube :: for Multi-Dim advection, max number of horizontal directions
      INTEGER    maxcube
      PARAMETER( maxcube   = 2 )

#ifdef ALLOW_CG2D_NSA
!     Parameter that is needed for the tape complev_cg2d_iter
!     cannot be smaller than the allowed number of iterations in cg2d
!     (numItersMax >= cg2dMaxIters in data-file)
      INTEGER numItersMax
      PARAMETER ( numItersMax = 100 )
#endif

#endif /* ALLOW_AUTODIFF_TAMC */
!     ================================================================
!     END OF HEADER TAMC
!     ================================================================
