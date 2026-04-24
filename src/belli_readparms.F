#include "BELLI_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_READPARMS

C !INTERFACE:
      SUBROUTINE IHOP_READPARMS( myThid )

C     !DESCRIPTION:
C     Initialize IHOP variables and constants, read in data.ihop

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_MNC
# include "MNC_PARAMS.h"
#endif
#include "BELLI_SIZE.h"
#include "BELLI.h"
#ifdef ALLOW_COST
# include "BELLI_COST.h"
#endif

C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER myThid
C     !OUTPUT PARAMETERS:

C    !FUNCTIONS:
      INTEGER  ILBLNK
      EXTERNAL ILBLNK
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
#ifdef ALLOW_BELLI
C     msgBuf :: Informational/error message buffer
C     iUnit  :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(128) fname
      LOGICAL exst
      INTEGER iUnit, IL, num_file
      INTEGER i, j

      NAMELIST /IHOP_PARM01/
     &      BELLI_fileroot,
     &      BELLI_title,
     &      IHOP_dumpFreq,
     &      BELLI_MNC

      NAMELIST /IHOP_PARM02/
     &      BELLI_topopt,
     &      BELLI_botopt,
     &      BELLI_runopt,
#ifdef BELLI_MULTIPLE_SOURCES
     &      BELLI_nsd,
#endif
#ifdef BELLI_MULTIPLE_RECEIVER_DEPTHS
     &      BELLI_nrd,
#endif
#ifdef BELLI_MULTIPLE_RECEIVER_RANGES
     &      BELLI_nrr,
#endif
     &      BELLI_nalpha,
     &      BELLI_freq,
     &      BELLI_depth,
     &      BELLI_bcsound,
     &      BELLI_bcsoundshear,
     &      IHOP_brho,
     &      BELLI_bcsoundI,
     &      BELLI_bcsoundshearI,
     &      BELLI_sd,
     &      BELLI_rd,
     &      BELLI_rr,
     &      BELLI_alpha,
     &      BELLI_step

      NAMELIST /IHOP_PARM03/
     &      writeDelay,
     &      useSSPFile,
     &      BELLI_interpfile,
     &      belli_iter,
     &      BELLI_npts_range,
     &      BELLI_npts_idw,
     &      belli_xc,
     &      belli_yc,
     &      belli_idw_weights,
     &      belli_ranges

#ifdef ALLOW_COST
      NAMELIST /IHOP_COST_NML/
     &      ihopdoncoutput,
     &      bObs_Dir,
     &      bObs_Files,
     &      mult_belli
#endif
CEOP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( .NOT.useIHOP ) THEN
C-    pkg IHOP is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.ihop is found
         CALL PACKAGES_UNUSED_MSG( 'useIHOP', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C     Open and read data.ihop file
      WRITE(msgBuf,'(A)') 'IHOP_READPARMS: opening data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                     'data.ihop', 'IHOP_READPARMS',
     O                     iUnit,
     I                     myThid )

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Default values for IHOP_PARM01
      BELLI_fileroot = ''
      BELLI_title    = 'ihop'
      BELLI_MNC      = useMNC
      IHOP_dumpFreq = dumpFreq

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Read parameters from open data file: IHOP_PARM01
      READ(UNIT=iUnit,NML=IHOP_PARM01)
      WRITE(msgBuf,'(A)')
     &    'IHOP_READPARMS: finished reading IHOP_PARM01 in data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C     Check paramters of IHOP_PARM01
      IF ( BELLI_fileroot.EQ.'' ) THEN
       WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "BELLI_fileroot" ',
     &  'Provide fileroot for ihop package'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF
      IF ( BELLI_title.EQ.'ihop' ) THEN
       WRITE(msgBuf, '(2A)') 
     &  '**WARNING** IHOP_READPARMS: "BELLI_title"',
     &  ' using default title name'
       CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Default values for IHOP_PARM02
      BELLI_topopt   = '      '
      BELLI_botopt   = '  '
      BELLI_runopt   = ' '
      BELLI_nsd      = 1
      BELLI_nrd      = 1
      BELLI_nrr      = 1
      BELLI_nalpha   = 1
      BELLI_freq     = 0. _d 0
      BELLI_depth    = 0. _d 0
      BELLI_bcsound  = 0. _d 0
      IHOP_brho     = 0. _d 0
      DO i=1,nsd
        BELLI_sd(i)  = 0. _d 0
      ENDDO
      DO i=1,nrd
        BELLI_rd(i)  = 0. _d 0
      ENDDO
      DO i=1,nrr
        BELLI_rr(i)  = 0. _d 0
      ENDDO
      BELLI_alpha    = 0. _d 0
      BELLI_step     = 0. _d 0
      BELLI_bcsoundI         = 0. _d 0
      BELLI_bcsoundshearI    = 0. _d 0
      BELLI_bcsoundshear     = 0. _d 0

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Read parameters from open data file: IHOP_PARM02
      READ(UNIT=iUnit,NML=IHOP_PARM02)
      WRITE(msgBuf,'(A)')
     &    'IHOP_READPARMS: finished reading IHOP_PARM02 in data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C     Check paramters of IHOP_PARM02
      IF ( BELLI_topopt(1:1) .NE. 'N' .AND.
     &     BELLI_topopt(1:1) .NE. 'C' .AND.
     &     BELLI_topopt(1:1) .NE. 'P' .AND.
     &     BELLI_topopt(1:1) .NE. 'S' .AND.
     &     BELLI_topopt(1:1) .NE. 'Q' .AND.
     &     BELLI_topopt(1:1) .NE. 'A'      ) THEN

       WRITE(msgBuf, '(2A)') 
     &  'S/R IHOP_READPARMS: "BELLI_topopt(1:1)" ',
     &  'unknown option for SSP approximation'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

      IF ( BELLI_topopt(2:2) .NE. 'V' .AND.
     &     BELLI_topopt(2:2) .NE. 'R' .AND.
     &     BELLI_topopt(2:2) .NE. 'A' .AND.
     &     BELLI_topopt(2:2) .NE. 'G' .AND.
     &     BELLI_topopt(2:2) .NE. 'F' .AND.
     &     BELLI_topopt(2:2) .NE. 'W' .AND.
     &     BELLI_topopt(2:2) .NE. 'P'      ) THEN

       WRITE(msgBuf, '(2A)') 
     &  'S/R IHOP_READPARMS: "BELLI_topopt(2:2)" ',
     &  'unknown top boundary condition type'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

      IF ( BELLI_topopt(3:3) .NE. 'N' .AND.
     &     BELLI_topopt(3:3) .NE. 'F' .AND.
     &     BELLI_topopt(3:3) .NE. 'M' .AND.
     &     BELLI_topopt(3:3) .NE. 'W' .AND.
     &     BELLI_topopt(3:3) .NE. 'Q' .AND.
     &     BELLI_topopt(3:3) .NE. 'L'      ) THEN

       WRITE(msgBuf, '(2A)')
     &  'S/R IHOP_READPARMS: "BELLI_topopt(3:3)" ',
     &  'unknown attenuation units'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

      IF ( BELLI_botopt(1:1) .NE. 'V' .AND.
     &     BELLI_botopt(1:1) .NE. 'R' .AND.
     &     BELLI_botopt(1:1) .NE. 'A' .AND.
     &     BELLI_botopt(1:1) .NE. 'G' .AND.
     &     BELLI_botopt(1:1) .NE. 'F' .AND.
     &     BELLI_botopt(1:1) .NE. 'W' .AND.
     &     BELLI_botopt(1:1) .NE. 'P'      ) THEN

       WRITE(msgBuf, '(2A)') 
     &  'S/R IHOP_READPARMS: "BELLI_botopt(1:1)" ',
     &  'unknown bottom boundary condition type'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

      IF ( BELLI_botopt(2:2) .NE. ' ' .AND.
     &     BELLI_botopt(2:2) .NE. '*' .AND.
     &     BELLI_botopt(2:2) .NE. '~'      ) THEN

       WRITE(msgBuf, '(2A)')
     &  'S/R IHOP_READPARMS: "BELLI_botopt(2:2)" ',
     &  'either select or deselect bathymetry file'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Default values for IHOP_PARM03 some are populated in
C     ihop_init_fixed.f
      writeDelay        = .TRUE.
      useSSPFile        = .FALSE.
      BELLI_interpfile   = ''
      DO i=1,nts
       belli_iter(i)     = -1
      ENDDO
      BELLI_npts_range   = 0
      BELLI_npts_idw     = 0
      DO i=1,BELLI_MAX_RANGE
       DO j=1,BELLI_MAX_NC_SIZE
        belli_xc(i,j)            = 0. _d 0
        belli_yc(i,j)            = 0. _d 0
        belli_idw_weights(i,j)   = 0. _d 0
       ENDDO
       belli_ranges(i)  = 0. _d 0
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Read parameters from open data file: IHOP_PARM03
      READ(UNIT=iUnit,NML=IHOP_PARM03)
      WRITE(msgBuf,'(A)')
     &    'IHOP_READPARMS: finished reading IHOP_PARM03 in data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C     Check paramters of IHOP_PARM03
      IF ( BELLI_interpfile.EQ.'' ) THEN
       WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "BELLI_interpfile" ',
     &  'Provide interpolation file name for ihop package'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

#ifdef ALLOW_COST
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Default values for IHOP_COST_NML
      ihopdoncoutput = .true.
      bObs_Dir    = ' '
      DO num_file = 1, NFILESMAX_BELLI
       bObs_Files(num_file)     = ' '
       mult_belli(num_file)         = 1. _d 0
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Read parameters from open data file: IHOP_COST_NML
      READ(UNIT=iUnit,NML=IHOP_COST_NML)
      WRITE(msgBuf,'(A)')
     & 'IHOP_READPARMS: finished reading IHOP_COST_NML in data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C     Check paramters of IHOP_COST_NML
      DO num_file=1,NFILESMAX_BELLI
       IF ( bObs_Files(num_file).NE.' ' ) THEN
        IL = ILNBLNK( bObs_Files(num_file) )
        fname = bObs_Files(num_file)(1:IL)//'.nc'
        INQUIRE( FILE=fname, EXIST=exst )
        IF ( .NOT.exst ) THEN
         WRITE(msgBuf,'(2A)') 'S/R IHOP_READPARMS: "bObs_Files " ',
     &    'Provide observations file name for ihop package'
         CALL PRINT_ERROR( msgBuf, myThid )
         STOP 'ABNORMAL END: S/R IHOP_READPARMS'
        ENDIF
       ENDIF
      ENDDO
#endif /* ALLOW_COST */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C     Make sure that we locally honor the global MNC on/off flag
      BELLI_MNC = BELLI_MNC .AND. useMNC
#ifndef ALLOW_MNC
C     Fix to avoid running without getting any output:
      BELLI_MNC = .FALSE.
#endif
      BELLI_MDSIO = (.NOT. BELLI_MNC) .OR. outputTypesInclusive

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_BELLI */
      RETURN
      END