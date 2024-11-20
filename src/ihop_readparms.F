#include "IHOP_OPTIONS.h"

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
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef ALLOW_COST
# include "IHOP_COST.h"
#endif

C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER myThid
C     !OUTPUT PARAMETERS:
CEOP

C    !FUNCTIONS:
      INTEGER  ILBLNK
      EXTERNAL ILBLNK
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
#ifdef ALLOW_IHOP
C     msgBuf :: Informational/error message buffer
C     iUnit  :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(128) fname
      LOGICAL exst
      INTEGER iUnit, IL, num_file
      INTEGER i, j

      NAMELIST /IHOP_PARM01/
     &      IHOP_fileroot,
     &      IHOP_title,
     &      IHOP_dumpFreq,
     &      IHOP_MNC

      NAMELIST /IHOP_PARM02/
     &      IHOP_topopt,
     &      IHOP_botopt,
     &      IHOP_runopt,
#ifdef IHOP_MULTIPLE_SOURCES
     &      IHOP_nsd,
#endif
#ifdef IHOP_MULTIPLE_RECEIVER_DEPTHS
     &      IHOP_nrd,
#endif
#ifdef IHOP_MULTIPLE_RECEIVER_RANGES
     &      IHOP_nrr,
#endif
     &      IHOP_nalpha,
     &      IHOP_freq,
     &      IHOP_depth,
     &      IHOP_bcsound,
     &      IHOP_bcsoundshear,
     &      IHOP_brho,
     &      IHOP_bcsoundI,
     &      IHOP_bcsoundshearI,
     &      IHOP_sd,
     &      IHOP_rd,
     &      IHOP_rr,
     &      IHOP_alpha,
     &      IHOP_step

      NAMELIST /IHOP_PARM03/
     &      writeDelay,
     &      useSSPFile,
     &      IHOP_interpfile,
     &      ihop_iter,
     &      IHOP_npts_range,
     &      IHOP_npts_idw,
     &      ihop_xc,
     &      ihop_yc,
     &      ihop_idw_weights,
     &      ihop_ranges

#ifdef ALLOW_COST
      NAMELIST /IHOP_COST_NML/
     &      ihopdoncoutput,
     &      IHOPObs_Dir,
     &      IHOPObs_Files,
     &      mult_ihop
#endif

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
      IHOP_fileroot = ''
      IHOP_title    = 'ihop'
      IHOP_MNC      = .false.
      IHOP_dumpFreq = dumpFreq

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Read parameters from open data file: IHOP_PARM01
      READ(UNIT=iUnit,NML=IHOP_PARM01)
      WRITE(msgBuf,'(A)')
     &    'IHOP_READPARMS: finished reading IHOP_PARM01 in data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C     Check paramters of IHOP_PARM01
      IF ( IHOP_fileroot.EQ.'' ) THEN
       WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "IHOP_fileroot" ',
     &  'Provide fileroot for ihop package'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF
      IF ( IHOP_title.EQ.'ihop' ) THEN
       WRITE(msgBuf, '(2A)') '**WARNING** IHOP_READPARMS: "IHOP_title"',
     &  ' using default title name'
       CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Default values for IHOP_PARM02
      IHOP_topopt   = '      '
      IHOP_botopt   = '  '
      IHOP_runopt   = ' '
      IHOP_nsd      = 1
      IHOP_nrd      = 1
      IHOP_nrr      = 1
      IHOP_nalpha   = 1
      IHOP_freq     = 0. _d 0
      IHOP_depth    = 0. _d 0
      IHOP_bcsound  = 0. _d 0
      IHOP_brho     = 0. _d 0
      DO i=1,nsd
        IHOP_sd(i)  = 0. _d 0
      ENDDO
      DO i=1,nrd
        IHOP_rd(i)  = 0. _d 0
      ENDDO
      DO i=1,nrr
        IHOP_rr(i)  = 0. _d 0
      ENDDO
      IHOP_alpha    = 0. _d 0
      IHOP_step     = 0. _d 0
      IHOP_bcsoundI         = 0. _d 0
      IHOP_bcsoundshearI    = 0. _d 0
      IHOP_bcsoundshear     = 0. _d 0

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Read parameters from open data file: IHOP_PARM02
      READ(UNIT=iUnit,NML=IHOP_PARM02)
      WRITE(msgBuf,'(A)')
     &    'IHOP_READPARMS: finished reading IHOP_PARM02 in data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C     Check paramters of IHOP_PARM02
      IF ( IHOP_topopt(1:1) .NE. 'N' .AND.
     &     IHOP_topopt(1:1) .NE. 'C' .AND.
     &     IHOP_topopt(1:1) .NE. 'P' .AND.
     &     IHOP_topopt(1:1) .NE. 'S' .AND.
     &     IHOP_topopt(1:1) .NE. 'Q' .AND.
     &     IHOP_topopt(1:1) .NE. 'A'      ) THEN

       WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "IHOP_topopt(1:1)" ',
     &  'unknown option for SSP approximation'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

      IF ( IHOP_topopt(2:2) .NE. 'V' .AND.
     &     IHOP_topopt(2:2) .NE. 'R' .AND.
     &     IHOP_topopt(2:2) .NE. 'A' .AND.
     &     IHOP_topopt(2:2) .NE. 'G' .AND.
     &     IHOP_topopt(2:2) .NE. 'F' .AND.
     &     IHOP_topopt(2:2) .NE. 'W' .AND.
     &     IHOP_topopt(2:2) .NE. 'P'      ) THEN

       WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "IHOP_topopt(2:2)" ',
     &  'unknown top boundary condition type'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

      IF ( IHOP_topopt(3:3) .NE. 'N' .AND.
     &     IHOP_topopt(3:3) .NE. 'F' .AND.
     &     IHOP_topopt(3:3) .NE. 'M' .AND.
     &     IHOP_topopt(3:3) .NE. 'W' .AND.
     &     IHOP_topopt(3:3) .NE. 'Q' .AND.
     &     IHOP_topopt(3:3) .NE. 'L'      ) THEN

       WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "IHOP_topopt(3:3)" ',
     &  'unknown attenuation units'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

      IF ( IHOP_botopt(1:1) .NE. 'V' .AND.
     &     IHOP_botopt(1:1) .NE. 'R' .AND.
     &     IHOP_botopt(1:1) .NE. 'A' .AND.
     &     IHOP_botopt(1:1) .NE. 'G' .AND.
     &     IHOP_botopt(1:1) .NE. 'F' .AND.
     &     IHOP_botopt(1:1) .NE. 'W' .AND.
     &     IHOP_botopt(1:1) .NE. 'P'      ) THEN

       WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "IHOP_botopt(1:1)" ',
     &  'unknown bottom boundary condition type'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

      IF ( IHOP_botopt(2:2) .NE. ' ' .AND.
     &     IHOP_botopt(2:2) .NE. '*' .AND.
     &     IHOP_botopt(2:2) .NE. '~'      ) THEN

       WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "IHOP_botopt(2:2)" ',
     &  'either select or deselect bathymetry file'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Default values for IHOP_PARM03 some are populated in
C     ihop_init_fixed.f
      writeDelay        = .TRUE.
      useSSPFile        = .FALSE.
      IHOP_interpfile   = ''
      DO i=1,nts
       ihop_iter(i)     = -1
      ENDDO
      IHOP_npts_range   = 0
      IHOP_npts_idw     = 0
      DO i=1,IHOP_MAX_RANGE
       DO j=1,IHOP_MAX_NC_SIZE
        ihop_xc(i,j)            = 0. _d 0
        ihop_yc(i,j)            = 0. _d 0
        ihop_idw_weights(i,j)   = 0. _d 0
       ENDDO
       ihop_ranges(i)  = 0. _d 0
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Read parameters from open data file: IHOP_PARM03
      READ(UNIT=iUnit,NML=IHOP_PARM03)
      WRITE(msgBuf,'(A)')
     &    'IHOP_READPARMS: finished reading IHOP_PARM03 in data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C     Check paramters of IHOP_PARM03
      IF ( IHOP_interpfile.EQ.'' ) THEN
       WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "IHOP_interpfile" ',
     &  'Provide interpolation file name for ihop package'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
      ENDIF

#ifdef ALLOW_COST
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Default values for IHOP_COST_NML
      ihopdoncoutput = .true.
      IHOPObs_Dir    = ' '
      DO num_file = 1, NFILESMAX_IHOP
       IHOPObs_Files(num_file)     = ' '
       mult_ihop(num_file)         = 1. _d 0
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Read parameters from open data file: IHOP_COST_NML
      READ(UNIT=iUnit,NML=IHOP_COST_NML)
      WRITE(msgBuf,'(A)')
     &    'IHOP_READPARMS: finished reading IHOP_COST_NML in data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C     Check paramters of IHOP_COST_NML
      DO num_file=1,NFILESMAX_IHOP
       IF ( IHOPObs_Files(num_file).NE.' ' ) THEN
        IL = ILNBLNK( IHOPObs_Files(num_file) )
        fname = IHOPObs_Files(num_file)(1:IL)//'.nc'
        INQUIRE( FILE=fname, EXIST=exst )
        IF ( .NOT.exst ) THEN
         WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "IHOPObs_Files " ',
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

CC     Make sure that we locally honor the global MNC on/off flag
C      IHOP_MNC = IHOP_MNC .AND. useMNC
C#ifndef ALLOW_MNC
CC     Fix to avoid running without getting any output:
C      IHOP_MNC = .FALSE.
C#endif
C      IHOP_MDSIO = (.NOT. IHOP_MNC) .OR. outputTypesInclusive

      _END_MASTER(myThid)
C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_IHOP */

      RETURN
      END
