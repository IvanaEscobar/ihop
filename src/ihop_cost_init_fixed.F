#include "IHOP_OPTIONS.h"
#include "AD_CONFIG.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_COST_INIT_FIXED

C !INTERFACE:
      SUBROUTINE IHOP_COST_INIT_FIXED( myThid )

C    !DESCRIPTION:
C    Set contributions to cost function and the cost function
C    itself to zero. Both are defined in header file "IHOP_COST.h".

C    !USES:
      IMPLICIT NONE
C == Global Variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"
#include "IHOP_SIZE.h"
#ifdef ALLOW_COST
# include "IHOP_COST.h"
#endif
#ifdef ALLOW_CAL
# include "cal.h"
#endif
#ifdef HAVE_NETCDF
#include "netcdf.inc"
#endif

C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER myThid
C     !OUTPUT PARAMETERS:
CBOP

C     !FUNCTIONS:
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
#if ( defined (ALLOW_COST) && defined (ALLOW_IHOP) )
      CHARACTER*(MAX_LEN_MBUF) :: msgBuf
      CHARACTER*(MAX_LEN_FNAM) IHOPObs_File, fnamedatanc
      CHARACTER*(MAX_LEN_FNAM) fnameequinc
      CHARACTER*(MAX_LEN_FNAM) adfnameequinc, tlfnameequinc
      CHARACTER*(MAX_LEN_FNAM) fnameequincglo
      CHARACTER*(MAX_LEN_FNAM) adfnameequincglo, tlfnameequincglo
      LOGICAL exst
      INTEGER :: iUnit
      INTEGER :: JL, IL, err, ncid, dimid1, dimid2
      INTEGER :: varid1a, varid1b, varid2, varid3
      INTEGER :: varid4a, varid4b, varid4c
      INTEGER :: varid5a, varid5b, varidEquiv(2)
      INTEGER :: weighIhopObs, obsIsInRunTime
      INTEGER :: ObsNo_valid, ihopObsInTile
      INTEGER :: tmpdate(4), tmpdiff(4)
      INTEGER :: num_file, hh, obsno_hh
      INTEGER :: bi, bj, k, i, j, l, m
      INTEGER :: iG, jG
      _RL :: yymmddMin,yymmddMax
      _RL :: hhmmssMin,hhmmssMax
      _RL :: tmpyymmdd(NOBSMAX_IHOP), tmphhmmss(NOBSMAX_IHOP), diffsecs
      _RL :: tmp_ihopObs_time(NOBSMAX_IHOP)
      _RL :: tmp_lon, tmp_lon2(NOBSMAX_IHOP), tmp_lat2(NOBSMAX_IHOP)
      _RL :: lon_cur, lat_cur
      _RL :: tmp_depth2(NOBSMAX_IHOP), depth_cur
      _RL :: tmp_uncert2(NOBSMAX_IHOP)
      _RL :: ihopObs_i
      _RL :: ihopObs_j
      _RL :: ihopObs_k1
      _RL :: ihopObs_k2
      _RL :: lon_1
      _RL :: lon_2
      _RL :: lat_1
      _RL :: lat_2
      _RL :: depth_1
      _RL :: depth_2
      _RL :: lon_tmp1, lon_tmp2, lat_tmp1, lat_tmp2
      _RL :: lat_fac, lon_fac, depth_fac

      WRITE(msgBuf,'(A)') ''// CHAR(10) //
     & '                    iHOP COST general parameters :'
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )

      JL = ILNBLNK( IHOPObs_Dir )
      IF (JL.NE.0) THEN
        WRITE(msgBuf,'(2A)') '  IHOPObs_Dir ',IHOPObs_Dir(1:JL)
      ELSE
        WRITE(msgBuf,'(2A)') '  IHOPObs_Dir ','./'
      ENDIF
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )

      WRITE(msgBuf,'(A,L5)') '  ihopDoNcOutput ', ihopDoNcOutput
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )


!====================================
!       INIT IHOP OBS VARS
!====================================
      _BEGIN_MASTER( myThid )

      yymmddMin=modelstartdate(1)
      yymmddMax=modelenddate(1)
      hhmmssMin=modelstartdate(2)
      hhmmssMax=modelenddate(2)

      DO num_file = 1, NFILESMAX_IHOP
        ObsNo_hh = 0

        IHOPObs_File=' '
        IL = ILNBLNK( IHOPObs_Files(num_file) )
        IF (IL.NE.0) THEN
          WRITE(IHOPObs_File,'(A)') IHOPObs_Files(num_file)(1:IL)
          WRITE(msgBuf,'(A,I3,2A)') ''// CHAR(10) //
     &     '                    ihop file #', num_file,
     &     ' is ', IHOPObs_File(1:IL)
          CALL PRINT_MESSAGE( msgBuf,
     &         standardMessageUnit, SQUEEZE_RIGHT, myThid )
        ENDIF

        IL = ILNBLNK( IHOPObs_File )
        IF (IL.NE.0) THEN
!===========================================================
!       READ IHOP OBS NETCDF IDs
!===========================================================
          WRITE(fnamedatanc,'(2A)') IHOPObs_File(1:IL),'.nc'
          err = NF_OPEN( fnamedatanc, NF_NOWRITE, 
     &          fiddata_obs(num_file) )

!1) Read number of observations and available dimensions:
          ncid = fiddata_obs(num_file)
          err = NF_INQ_DIMID( ncid,'iOBS', dimid1 )
          err = NF_INQ_DIMLEN( ncid, dimid1, obsNo(num_file) )

          WRITE(msgBuf,'(A,I9)')
     &       '  # of observations in file                   =',
     &       obsNo(num_file)
          CALL PRINT_MESSAGE( msgBuf,
     &         standardMessageUnit, SQUEEZE_RIGHT, myThid )

          err =       NF_INQ_VARID( ncid, 'ihop_YYYYMMDD', varid1a )
          err = err + NF_INQ_VARID( ncid, 'ihop_HHMMSS', varid1b )
          err = err + NF_INQ_VARID( ncid, 'ihop_iter', varid2 )

          IF (err.NE.NF_NOERR) THEN
            WRITE(msgBuf,'(4A)')
     &       'IHOP_COST_INIT_FIXED: file ', IHOPObs_File(1:IL),
     &       '.nc is not in the pkg/ihop cost format: ',
     &       'check iter,uncert'
            CALL PRINT_ERROR( msgBuf, myThid )
          ENDIF
! Assume all ihop observations are instantaneous

!2 Weight for each ihop obs (if missing, assume equal weights):
          weighIhopObs = 1

          err = NF_INQ_VARID( ncid, 'ihop_uncert', varid3 )
          IF (err.NE.NF_NOERR) THEN
            WRITE(msgBuf,'(2A)')
     &        '      ihop obs file does not have uncertainties;',
     &        '      assume all obs datum are weighed evenly '
            CALL PRINT_MESSAGE( msgBuf,
     &           standardMessageUnit, SQUEEZE_RIGHT, myThid )

            weighIhopObs = 0

          ENDIF

!3 Read ihop observation positions and values:
          err =       NF_INQ_VARID( ncid, 'ihop_x', varid4a )
          err = err + NF_INQ_VARID( ncid, 'ihop_y', varid4b )
          err = err + NF_INQ_VARID( ncid, 'ihop_z', varid4c )

          IF (err.NE.NF_NOERR) THEN
            WRITE(msgBuf,'(3A)')
     &       'IHOP_COST_INIT_FIXED: file ', IHOPObs_File(1:IL),
     &       '.nc is not in the pkg/ihop cost format, check x,y,z'
            CALL PRINT_ERROR( msgBuf, myThid )
          ENDIF

          err =       NF_INQ_VARID( ncid, 'ihop_val', varid5a )

          IF (err.NE.NF_NOERR) THEN
           WRITE(msgBuf,'(3A)')
     &      'IHOP_COST_INIT_FIXED: file ', IHOPObs_File(1:IL),
     &      '.nc is not in the pkg/ihop cost format, check ihop_val'
           CALL PRINT_ERROR( msgBuf, myThid)
          ENDIF

!===========================================================
!       SET IHOP GLOBAL EQUIV NETCDF VARIABLES
!===========================================================
! Variable names for model equivalent files
          ihop_nameequi = 'mod_val'
          ihop_namemask = 'mod_mask'
          ihop_nameval  = 'ihop_val'
          ihop_nameuncert = 'ihop_uncert'

!4 Create global file for ihop model equivalent
          IF (myProcId.EQ.0) THEN
            JL  = ILNBLNK( IHOPObs_Dir )

            WRITE(fnameequincglo,'(3A)')
     &       IHOPObs_Dir(1:JL),
     &       IHOPObs_File(1:IL),'.equi.nc'

# ifdef ALLOW_ADJOINT_RUN
            WRITE(adfnameequincglo,'(4A)')
     &       IHOPObs_Dir(1:JL),'ad',
     &       IHOPObs_File(1:IL),'.equi.nc'
# endif
# ifdef ALLOW_TANGENTLINEAR_RUN
            WRITE(tlfnameequincglo,'(4A)')
     &       IHOPObs_Dir(1:JL),'tl',
     &       IHOPObs_File(1:IL),'.equi.nc'
# endif

            JL = ILNBLNK( fnameequincglo )
            INQUIRE( FILE = fnameequincglo(1:JL), EXIST = exst )
            IF (.NOT.exst) THEN
! Create NetCDF file
              err = NF_CREATE( fnameequincglo(1:JL), NF_CLOBBER,
     &              fidglobal(num_file) )
              err = NF_DEF_DIM( fidglobal(num_file),'iOBS',
     &              ObsNo(num_file), dimid1 )

              err = NF_DEF_VAR( fidglobal(num_file), 'mod_val',
     &              NF_DOUBLE, 1, dimid1, varidEquiv(1) )
              err = NF_PUT_ATT_DOUBLE( fidglobal(num_file),
     &              varidEquiv(1), '_FillValue',
     &              NF_DOUBLE, 1, -9999. _d 0 )

              err = NF_DEF_VAR( fidglobal(num_file), 'mod_mask',
     &              NF_DOUBLE, 1, dimid1, varidEquiv(2) )
              err = NF_PUT_ATT_DOUBLE( fidglobal(num_file),
     &              varidEquiv(2), '_FillValue',
     &              NF_DOUBLE, 1, -9999. _d 0 )

! Write buffer to file, reopen file with WRITE permissions
              err = NF_ENDDEF( fidglobal(num_file) )
              err = NF_CLOSE( fidglobal(num_file) )
              err = NF_OPEN( fnameequincglo(1:JL), NF_WRITE, 
     &              fidglobal(num_file) )

            ELSE
! Obtain existing NetCDF file id
              err = NF_OPEN( fnameequincglo(1:JL), NF_WRITE,
     &              fidglobal(num_file) )
            ENDIF

# ifdef ALLOW_ADJOINT_RUN
            JL = ILNBLNK( adfnameequincglo )
            INQUIRE( FILE = adfnameequincglo(1:JL), EXIST = exst )
            IF (.NOT.exst) THEN
! Create NetCDF AD file
              err = NF_CREATE( adfnameequincglo(1:JL), NF_CLOBBER,
     &              fidadglobal(num_file) )
              err = NF_DEF_DIM( fidadglobal(num_file), 'iOBS',
     &              ObsNo(num_file), dimid1 )

              err = NF_DEF_VAR( fidadglobal(num_file), 'mod_val',
     &              NF_DOUBLE, 1, dimid1, varidEquiv(1) )
              err = NF_PUT_ATT_DOUBLE( fidadglobal(num_file),
     &              varidEquiv(1), '_FillValue',
     &              NF_DOUBLE, 1, 0. _d 0 )

! Write buffer to file, reopen file with WRITE permissions
              err = NF_ENDDEF( fidadglobal(num_file) )
              err = NF_CLOSE( fidadglobal(num_file) )
              err = NF_OPEN( adfnameequincglo(1:JL), NF_WRITE,
     &              fidadglobal(num_file) )

            ELSE
! Obtain existing AD NetCDF file id
              err = NF_OPEN( adfnameequincglo(1:JL), NF_WRITE,
     &              fidadglobal(num_file) )
            ENDIF

# endif
# ifdef ALLOW_TANGENTLINEAR_RUN
            JL = ILNBLNK( tlfnameequincglo )
            INQUIRE( FILE = tlfnameequincglo(1:JL), EXIST = exst )
            IF (.NOT.exst) THEN
! Create NetCDF TL file
              err = NF_CREATE( tlfnameequincglo(1:JL), NF_CLOBBER,
     &              fidftlglobal(num_file) )
              err = NF_DEF_DIM( fidftlglobal(num_file), 'iOBS',
     &              ObsNo(num_file), dimid1 )

              err = NF_DEF_VAR( fidftlglobal(num_file), 'mod_val',
     &              NF_DOUBLE, 1, dimid1, varidEquiv(1) )
              err = NF_PUT_ATT_DOUBLE( fidftlglobal(num_file),
     &              varidEquiv(1), '_FillValue',
     &              NF_DOUBLE, 1, 0. _d 0 )

! Write buffer to file, reopen file with WRITE permissions
              err = NF_ENDDEF( fidftlglobal(num_file) )
              err = NF_CLOSE( fidftlglobal(num_file) )
              err = NF_OPEN( tlfnameequincglo(1:JL), NF_WRITE,
     &              fidftlglobal(num_file) )

            ELSE
! Obtain existing TL NetCDF file id
              err = NF_OPEN( tlfnameequincglo(1:JL), NF_WRITE,
     &              fidftlglobal(num_file) )
            ENDIF
# endif

          ENDIF !IF (myProcId.EQ.0)

!===========================================================
!       READ IHOP OBS NETCDF VARIABLES
!===========================================================
! NOTE: not applying reading in chunks
!       see pkg/profiles or pkg/obsfit for how to chunk.

! Default values
          DO k = 1, NOBSMAX_IHOP
            ihopObs_ind_glob(num_file,k) = 0
          ENDDO
          ObsNo_valid = 0

! Read time
          err = NF_GET_VARA_DOUBLE( ncid, varid1a, 
     &          1, ObsNo(num_file), tmpyymmdd )
          err = NF_GET_VARA_DOUBLE(ncid,varid1b,
     &          1, ObsNo(num_file), tmphhmmss )

          IF (err.NE.NF_NOERR) THEN
            WRITE(msgBuf,'(3A)')
     &       'IHOP_COST_INIT_FIXED: file ', IHOPObs_File(1:IL),
     &       '.nc was not read properly: check yyyymmdd, hhmmss.'
            CALL PRINT_ERROR( msgBuf, myThid)
          ENDIF

! Loop through number of observations, in obsfit this is a 'chunk'
          DO m = 1, ObsNo(num_file)
!AV if obs starts before model run, will be ignored even if it ends after
!AV model starts - ok for now
            IF (( (tmpyymmdd(m).GT.yymmddMin).OR.
     &            ((tmpyymmdd(m).EQ.yymmddMin) .AND.
     &             (tmphhmmss(m).GE.hhmmssMin)) ) .AND.
     &          ( (tmpyymmdd(m).LT.yymmddMax).OR.
     &            ((tmpyymmdd(m).EQ.yymmddMax) .AND.
     &             (tmphhmmss(m).LE.hhmmssMax)) )) THEN
              hh = INT(tmphhmmss(m)) / 10000

              IF (hh.LT.hoursPerDay) THEN
                obsIsInRunTime = 1
#ifdef ALLOW_CAL
                CALL CAL_FULLDATE( INT(tmpyymmdd(m)),
     &               INT(tmphhmmss(m)), tmpdate, myThid )
                CALL CAL_TIMEPASSED( modelstartdate, tmpdate,
     &               tmpdiff, myThid )
                CALL CAL_TOSECONDS( tmpdiff, diffsecs, myThid )
#endif /* ALLOW_CAL */

                diffsecs = diffsecs+nIter0*deltaTclock

              ELSE
! If tmphhmmss is out of range then disregard ihop obs
                obsIsInRunTime = 0
                diffsecs = -deltaTclock
                ObsNo_hh = ObsNo_hh+1

              ENDIF

            ELSE
             obsIsInRunTime = 0
             diffsecs = -deltaTclock

            ENDIF !IF (( (tmpyymmdd(m).GT.yymmddMin)

            IF (obsIsInRunTime.EQ.1) THEN
! If yes then store obs index position:
              ObsNo_valid = ObsNo_valid+1
              ihopObs_ind_glob(num_file,ObsNo_valid) = m

            ENDIF !IF (obsIsInRunTime.EQ.1)

! Regardless of whether ihop obs is in run time or not,
! Keep track of ihop obs index position and assign time to ihop obs
            tmp_ihopObs_time(m) = diffsecs

! Check that maximum size was not reached:
            IF (ObsNo(num_file).GT.NOBSMAX_IHOP) THEN
              WRITE(msgBuf,'(3A)')
     &         'IHOP_COST_INIT_FIXED: file ', IHOPObs_File(1:IL),
     &         '.nc was not read properly: increase NOBSMAX_IHOP.'
              CALL PRINT_ERROR( msgBuf, myThid )
            ENDIF

          ENDDO !DO m

! Calibrate number of VALID obs datum
          ObsNo(num_file) = ObsNo_valid

! Output ihop cost init info
          WRITE(msgBuf,'(A,I9)')
     &     '  # of ihop obs with erroneous HHMMSS values =', ObsNo_hh
          CALL PRINT_MESSAGE( msgBuf,
     &         standardMessageUnit, SQUEEZE_RIGHT, myThid )

          WRITE(msgBuf,'(A,I9)')
     &     '  # of ihop obs within time period  =', ObsNo(num_file)
          CALL PRINT_MESSAGE( msgBuf,
     &         standardMessageUnit, SQUEEZE_RIGHT, myThid )

! ==============================================================================
! Initiate tiled defaults for ihopObs before reading in from netcdf
          DO bj = 1, nSy
            DO bi = 1, nSx
              ObsNo_valid = 0
              ObsNo_tiled(num_file,bi,bj) = ObsNo(num_file)

              DO k = 1, NOBSMAX_IHOP
                ihopObs_time(num_file,k,bi,bj) =     -999. _d 0
                ihopObs_lon(num_file,k,bi,bj) =      -999. _d 0
                ihopObs_lat(num_file,k,bi,bj) =      -999. _d 0
                ihopObs_depth(num_file,k,bi,bj) =    -999. _d 0
                ihopObs_uncert(num_file,k,bi,bj) =      1. _d 0
                ihopObs_ind_glob_tiled(num_file,k,bi,bj) =    0
              ENDDO !DO k

! Read ALL info using respective varid and save to tmp, NO CHUNKS!
              IF (weighIhopObs.EQ.1) THEN
                err = NF_GET_VARA_DOUBLE( ncid, varid3,
     &                1, ObsNo_tiled(num_file,bi,bj), tmp_uncert2 )
              ENDIF

              err = NF_GET_VARA_DOUBLE( ncid, varid4a,
     &              1, ObsNo_tiled(num_file,bi,bj), tmp_lon2 )
              err = NF_GET_VARA_DOUBLE( ncid, varid4b,
     &              1, ObsNo_tiled(num_file,bi,bj), tmp_lat2 )
              err = NF_GET_VARA_DOUBLE( ncid, varid4c,
     &              1, ObsNo_tiled(num_file,bi,bj), tmp_depth2 )

              IF (err.NE.NF_NOERR) THEN
                WRITE(msgBuf,'(3A)')
     &           'IHOP_COST_INIT_FIXED: file ', IHOPObs_File(1:IL),
     &           '.nc was not read properly: check x,y,z,uncert.'
                CALL PRINT_ERROR( msgBuf, myThid )
              ENDIF

! Find obs datum in MITgcm spatial grid
! Loop through observations, NOTE: in obsfit this is a 'chunk'
              DO m = 1, ObsNo_tiled(num_file,bi,bj)
                ihopObsInTile = 1

! In a valid time?
                IF ( (tmp_ihopObs_time(m).LT.0. _d 0) .OR.
     &               (tmp_ihopObs_time(m).GT.modelend-modelstart) )
     &          THEN
                  ihopObsInTile = 0
                ENDIF

! Determine whether obs datum is in current tile domain
                lon_cur = tmp_lon2(m)
                lat_cur = tmp_lat2(m)
                depth_cur = tmp_depth2(m)

                IF (ihopObsInTile.EQ.1) THEN
                  IF (xC(sNx+1,1,bi,bj).LT.xC(1,1,bi,bj)) THEN
                    tmp_lon = xC(sNx+1,1,bi,bj)+360. _d 0
                  ELSE
                    tmp_lon = xC(sNx+1,1,bi,bj)
                  ENDIF

                  IF ( (xC(1,1,bi,bj).LE.lon_cur) .AND.
     &                 (tmp_lon.GT.lon_cur) .AND.
     &                 (yC(1,1,bi,bj).LE.lat_cur) .AND.
     &                 (yC(1,sNy+1,bi,bj).GT.lat_cur) ) THEN
! Do nothing

                  ELSEIF ( (xC(sNx+1,1,bi,bj).LT.xC(1,1,bi,bj)) .AND.
     &                     (xC(1,1,bi,bj).LE.lon_cur+360. _d 0) .AND.
     &                     (tmp_lon.GT.lon_cur+360. _d 0) .AND.
     &                     (yC(1,1,bi,bj).LE.lat_cur) .AND.
     &                     (yC(1,sNy+1,bi,bj).GT.lat_cur) ) THEN
                    lon_cur = lon_cur+360. _d 0

                  ELSE
! Not in this tile
                    ihopObsInTile = 0

                  ENDIF !IF ( (xC(1,1,bi,bj).LE.lon_cur)

! Determine i,j,k to the south-west-down of obs datum
! Set default values
                  ihopObs_i = -10
                  ihopObs_j = -10
                  ihopObs_k1 = -10
                  ihopObs_k2 = -10
                  lon_1 = -10
                  lon_2 = -10
                  lat_1 = -10
                  lat_2 = -10
                  depth_1 = -10
                  depth_2 = -10

                  IF (ihopObsInTile.EQ.1) THEN
! Value of k, below the data point:
! Case 1: above first depth level
                    IF (-rC(1).GT.depth_cur) THEN
                      ihopObs_k1 = 1
                      ihopObs_k2 = 1
                      depth_fac = 1. _d 0
                    ENDIF
! Case 2: below last depth level
                    IF (-rC(Nr).LE.depth_cur) THEN
                      ihopObs_k1 = Nr
                      ihopObs_k2 = Nr
                      depth_fac = 1. _d 0
                    ENDIF
! Case 3: between two depth levels
                    DO k = 1, Nr-1
                      IF ( (-rC(k).LE.depth_cur).AND.
     &                     (-rC(k+1).GT.depth_cur)) THEN
                        ihopObs_k1 = k
                        ihopObs_k2 = k+1
                        depth_1 = -rC(k)
                        depth_2 = -rC(k+1)
                        depth_fac =
     &                   (depth_cur-depth_1)/(depth_2-depth_1)
                      ENDIF
                    ENDDO

                    DO j = 1, sNy+1
                      DO i = 1, sNx+1
                        lon_tmp1 = xC(i,j,bi,bj)
                        lon_tmp2 = xC(i+1,j,bi,bj)
                        lat_tmp1 = yC(i,j,bi,bj)
                        lat_tmp2 = yC(i,j+1,bi,bj)

! Value of j, south of obs datum:
                        IF ( (lat_tmp1.LE.lat_cur) .AND.
     &                       (lat_tmp2.GT.lat_cur) ) THEN
                          ihopObs_j = j
                          lat_1 = lat_tmp1
                          lat_2 = lat_tmp2
                        ENDIF

! Value of i, west of obs datum:
                        IF (xC(i+1,j,bi,bj).LT.xC(1,j,bi,bj)) THEN
                          lon_tmp2=lon_tmp2+360
                        ENDIF
                        IF (xC(i,j,bi,bj).LT.xC(1,j,bi,bj)) THEN
                          lon_tmp1=lon_tmp1+360
                        ENDIF

                        IF ( (lon_tmp1.LE.lon_cur) .AND.
     &                       (lon_tmp2.GT.lon_cur) ) THEN
                          ihopObs_i = i
                          lon_1 = lon_tmp1
                          lon_2 = lon_tmp2
                        ENDIF

                      ENDDO !DO i
                    ENDDO !DO j

                  ENDIF !IF (ihopObsInTile.EQ.1)

                  IF ( ((ihopObs_i.EQ.-10).OR.(ihopObs_j.EQ.-10)) .OR.
     &                 ((ihopObs_k1.EQ.-10).OR.(ihopObs_k2.EQ.-10)) )
     &            THEN
                    ihopObsInTile = 0
                  ENDIF

! 'Reshape' observed datum to tiled structures
                  IF (ihopObsInTile.EQ.1) THEN
                    ObsNo_valid = ObsNo_valid+1
                    ihopObs_ind_glob_tiled(num_file,
     &               ObsNo_valid,bi,bj) = m

                    ihopObs_time(num_file,ObsNo_valid,bi,bj) =
     &               tmp_ihopObs_time(m)

                    IF (weighIhopObs.EQ.1) THEN
                      ihopObs_uncert(num_file,ObsNo_valid,bi,bj) =
     &                 tmp_uncert2(m)
                    ENDIF

                    ihopObs_lon(num_file,ObsNo_valid,bi,bj) = lon_cur
                    ihopObs_lat(num_file,ObsNo_valid,bi,bj) = lat_cur
                    ihopObs_depth(num_file,ObsNo_valid,bi,bj) =
     &               depth_cur

                    ihopObs_i_tiled(num_file,ObsNo_valid,bi,bj) = 
     &               ihopObs_i
                    ihopObs_j_tiled(num_file,ObsNo_valid,bi,bj) =
     &               ihopObs_j
                    ihopObs_k_tiled(num_file,ObsNo_valid,bi,bj) =
     &               ihopObs_k1

                  ENDIF !IF (ihopObsInTile.EQ.1) 3rd time

                ENDIF !IF (ihopObsInTile.EQ.1) 1st time

! ==============================================================================
! Check that maximum size was not reached:
                IF (m.GE.NOBSMAX_IHOP) THEN
                  WRITE(msgBuf,'(3A)')
     &             'IHOP_COST_INIT_FIXED: file ', IHOPObs_File(1:IL),
     &             '.nc was not read properly: increase NOBSMAX_IHOP.'
                  CALL PRINT_ERROR( msgBuf, myThid)
                ENDIF

              ENDDO !DO m

! Number of valid ihop obs in the file
              ObsNo_tiled(num_file,bi,bj) = ObsNo_valid

              WRITE(msgBuf,'(2(A,I4))')
     &         '  current tile is bi,bj                      =',
     &         bi,',',bj
              CALL PRINT_MESSAGE( msgBuf,
     &             standardMessageUnit, SQUEEZE_RIGHT, myThid )

              WRITE(msgBuf,'(A,I9)')
     &         '  # of ihop obs within tile and time period  =',
     &         ObsNo_tiled(num_file,bi,bj)
              CALL PRINT_MESSAGE( msgBuf,
     &             standardMessageUnit, SQUEEZE_RIGHT, myThid )

!===========================================================
! Create files for model counterparts to ihop obs
!===========================================================
              IF (obsNo_tiled(num_file,bi,bj).GT.0) THEN
                iG=bi+(myXGlobalLo-1)/sNx ! Kludge until unstructered tiles
                jG=bj+(myYGlobalLo-1)/sNy ! Kludge until unstructered tiles

                JL = ILNBLNK( IHOPObs_Dir )
! ALWAYS Write to nc file
                IF (ihopDoNcOutput) THEN
                  WRITE(fnameequinc,'(3A,2(I3.3,A))')
     &             IHOPObs_Dir(1:JL),IHOPObs_File(1:IL),
     &             '.',iG,'.',jG,'.equi.nc'
# ifdef ALLOW_ADJOINT_RUN
                  WRITE(adfnameequinc,'(4A,2(I3.3,A))')
     &             IHOPObs_Dir(1:JL),'ad',IHOPObs_File(1:IL),
     &             '.',iG,'.',jG,'.equi.nc'
# endif
# ifdef ALLOW_TANGENTLINEAR_RUN
                  WRITE(tlfnameequinc,'(4A,2(I3.3,A))')
     &             IHOPObs_Dir(1:JL),'ad',IHOPObs_File(1:IL),
     &             '.',iG,'.',jG,'.equi.nc'
# endif

                  JL = ILNBLNK( fnameequinc )
                  INQUIRE( FILE = fnameequinc(1:JL), EXIST = exst )
                  IF (.NOT.exst) THEN
                    CALL IHOP_COST_INIT_EQUI( num_file,
     &                   fiddata_obs(num_file), fnameequinc(1:JL),
     &                   fidfwd_obs(num_file,bi,bj),
     &                   ObsNo_tiled(num_file,bi,bj),
     &                   bi, bj, myThid )
                  ELSE
! Obtain existing NetCDF file id
                    err = NF_OPEN( fnameequinc(1:JL), NF_WRITE,
     &                    fidfwd_obs(num_file,bi,bj) )
                  ENDIF

# ifdef ALLOW_ADJOINT_RUN
                  JL = ILNBLNK( adfnameequinc )
                  INQUIRE( FILE = adfnameequinc(1:JL), EXIST = exst )
                  IF (.NOT.exst) THEN
                    CALL IHOP_COST_INIT_EQUI( num_file,
     &                   fiddata_obs(num_file), adfnameequinc(1:JL), 
     &                   fidadj_obs(num_file,bi,bj),
     &                   ObsNo_tiled(num_file,bi,bj),
     &                   bi, bj, myThid )
                  ELSE
! Obtain existing AD NetCDF file id
                    err = NF_OPEN( adfnameequinc(1:JL), NF_WRITE,
     &                    fidadj_obs(num_file,bi,bj) )
                  ENDIF

# endif /* ALLOW_ADJOINT_RUN */
# ifdef ALLOW_TANGENTLINEAR_RUN
                  JL = ILNBLNK( tlfnameequinc )
                  INQUIRE( FILE = tlfnameequinc(1:JL), EXIST = exst )
                  IF (.NOT.exst) THEN
                    CALL IHOP_COST_INIT_EQUI( num_file,
     &                   fiddata_obs(num_file),tlfnameequinc(1:JL),
     &                   fidftl_obs(num_file,bi,bj),
     &                   ObsNo_tiled(num_file,bi,bj),
     &                   bi, bj, myThid )
                  ELSE
! Obtain existing TL NetCDF file id
                    err = NF_OPEN( tlfnameequinc(1:JL), NF_WRITE,
     &                    fidftl_obs(num_file,bi,bj) )
                  ENDIF
# endif /* ALLOW_TANGENTLINEAR_RUN */
                ELSE
! Write to binary files...
                ENDIF !IF (ihopDoNcOutput)

              ENDIF !IF (obsNo_tiled(num_file,bi,bj).GT.0)

            ENDDO !DO bi
          ENDDO !DO bj
!===========================================================
        ELSE !IF (IL.NE.0)
          DO bj = 1, nSy
            DO bi = 1, nSx
              ObsNo(num_file) = 0
              ObsNo_tiled(num_file,bi,bj) = 0

              DO k = 1, NOBSMAX_IHOP
                ihopObs_lon(num_file,k,bi,bj) =      -999. _d 0
                ihopObs_lat(num_file,k,bi,bj) =      -999. _d 0
                ihopObs_depth(num_file,k,bi,bj) =    -999. _d 0
                ihopObs_uncert(num_file,k,bi,bj) =      1. _d 0
                ihopObs_ind_glob_tiled(num_file,k,bi,bj) =    0
              ENDDO

            ENDDO !DO bi
          ENDDO !DO bj

        ENDIF !IF (IL.NE.0)

      ENDDO !DO num_file

      _END_MASTER( myThid )
      _BARRIER

!===========================================================
      WRITE(msgBuf,'(A)') ' '// CHAR(10) //
     & '                    '//
     & '// ======================================================='
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     & '// ihop cost read obs data file >>> END <<<'
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     & '// ======================================================='
     & // CHAR(10) //'                    '
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )

#endif /* ALLOW_COST && ALLOW_IHOP*/

      RETURN
      END
