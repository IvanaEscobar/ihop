#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 
C     !ROUTINE: IHOP_READPARMS

C     !INTERFACE:
      SUBROUTINE IHOP_READPARMS( myThid )

C     !DESCRIPTION:
C     Initialize IHOP variables and constants, read in data.ihop

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef ALLOW_MNC
#include "MNC_PARAMS.h"
#endif

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_IHOP
C     !LOCAL VARIABLES:
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit
      INTEGER i, j, v

      NAMELIST /IHOP_PARM01/
     &      IHOP_fileroot,
     &      IHOP_title,
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
     &      IHOP_step,
     &      IHOP_zbox,
     &      IHOP_rbox

      NAMELIST /IHOP_PARM03/
     &      useSSPFile,
     &      IHOP_interpfile,
     &      ihop_iter,
     &      IHOP_npts_range,
     &      IHOP_npts_idw,
     &      ihop_v_xc,
     &      ihop_v_yc,
     &      ihop_v_weight,
     &      ihop_xc,
     &      ihop_yc,
     &      ihop_idw_weights,
     &      ihop_ranges

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

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Default values for IHOP_PARM01
      IHOP_fileroot = ''
      IHOP_title    = 'ihop'
      IHOP_MNC      = .false.

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Open and read data.ihop file
      WRITE(msgBuf,'(A)') 'IHOP_READPARMS: opening data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                     'data.ihop', 'IHOP_READPARMS',
     O                     iUnit,
     I                     myThid )

C     Read parameters from open data file
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

      _END_MASTER(myThid)
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      _BEGIN_MASTER(myThid)
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
      IHOP_zbox     = 0. _d 0
      IHOP_rbox     = 0. _d 0
      IHOP_bcsoundI         = 0. _d 0
      IHOP_bcsoundshearI    = 0. _d 0
      IHOP_bcsoundshear     = 0. _d 0

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Read parameters from open data file
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
      useSSPFile        = .TRUE.
      IHOP_interpfile   = ''
      ihop_iter         = 0
      IHOP_npts_range   = 0
      IHOP_npts_idw     = 0
      v=1
      DO i=1,IHOP_MAX_NC_SIZE
       DO j=1,IHOP_MAX_NC_SIZE
        ihop_v_xc(v)            = 0. _d 0
        ihop_v_yc(v)            = 0. _d 0
        ihop_v_weight(v)        = 0. _d 0

        ihop_xc(i,j)            = 0. _d 0
        ihop_yc(i,j)            = 0. _d 0
        ihop_idw_weights(i,j)   = 0. _d 0
        v=v+1
       ENDDO
       ihop_ranges(i)  = 0. _d 0
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Read parameters from open data file
      READ(UNIT=iUnit,NML=IHOP_PARM03)
      WRITE(msgBuf,'(A)')
     &    'IHOP_READPARMS: finished reading IHOP_PARM03 in data.ihop'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

CC     Check paramters of IHOP_PARM03
C      IF ( IHOP_interpfile.EQ.'' ) THEN
C       WRITE(msgBuf, '(2A)') 'S/R IHOP_READPARMS: "IHOP_interpfile" ',
C     &  'Provide interpolation file name for ihop package'  
C       CALL PRINT_ERROR( msgBuf, myThid )
C       STOP 'ABNORMAL END: S/R IHOP_READPARMS'
C      ENDIF

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
