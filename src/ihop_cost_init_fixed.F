#include "IHOP_OPTIONS.h"

      subroutine ihop_cost_init_fixed( myThid )

C     ==================================================================
C     SUBROUTINE ihop_cost_init_fixed
C     ==================================================================
C
C     o Set contributions to the cost function and the cost function
C       itself to zero. The cost function and the individual contribu-
C       tions are defined in the header file "IHOP_COST.h".
C
C     ==================================================================
C     SUBROUTINE ihop_cost_init_fixed
C     ==================================================================

      implicit none

C     == global variables ==
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

!     == routine arguments ==
      integer myThid

#if ( defined (ALLOW_COST) && defined (ALLOW_IHOP) )
!     !EXTERNAL VARIABLES:
      CHARACTER*(MAX_LEN_MBUF) :: msgBuf
!     !LOCAL VARIABLES:
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
      CHARACTER*(MAX_LEN_FNAM) ihopObsfile, fnamedatanc
      CHARACTER*(MAX_LEN_FNAM) fnameequinc, adfnameequinc
      CHARACTER*(MAX_LEN_FNAM) fnameequincglo, adfnameequincglo
      LOGICAL exst

C ==================== External Functions ==========================
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     == end of interface ==

      WRITE(msgbuf,'(a)') ' '
      CALL print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT, myThid)
      WRITE(msgbuf,'(a)') 'ihop general parameters :'
      JL  = ILNBLNK( ihopObsDir )
      IF (JL.NE.0) THEN
        WRITE(msgbuf,'(a,a)') '  ihopObsDir ',ihopObsDir(1:JL)
      ELSE
        WRITE(msgbuf,'(a,a)') '  ihopObsDir ','./'
      ENDIF
      CALL print_message(
     &   msgbuf, standardmessageunit, SQUEEZE_RIGHT, myThid)
      WRITE(msgbuf,'(a,l5)') '  ihopDoNcOutput ',ihopDoNcOutput
      CALL print_message(
     &   msgbuf, standardmessageunit, SQUEEZE_RIGHT, myThid)
      WRITE(msgbuf,'(a)') ' '
      CALL print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT, myThid)


!====================================
!       INIT IHOP OBS VARS
!====================================
      _BEGIN_MASTER( myThid )

      yymmddMin=modelstartdate(1)
      yymmddMax=modelenddate(1)
      hhmmssMin=modelstartdate(2)
      hhmmssMax=modelenddate(2)

      DO num_file=1,NFILESMAX_IHOP

       ObsNo_hh=0

       ihopObsFile=' '
       IL  = ILNBLNK( ihopObsFiles(num_file) )
       IF (IL.NE.0) THEN
        WRITE(ihopObsFile,'(A)')
     &     ihopObsFiles(num_file)(1:IL)
        WRITE(msgbuf,'(A)') ' '
        CALL print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT, myThid)
        WRITE(msgbuf,'(A,I3,2A)')
     &     'ihop file #',num_file,' is ', ihopObsFile(1:IL)
        CALL print_message(
     &     msgbuf, standardmessageunit, SQUEEZE_RIGHT, myThid)
       ENDIF

       IL  = ILNBLNK( ihopObsFile )
       IF (IL.NE.0) THEN

!===========================================================
!       READ IHOP OBS NETCDF IDs 
!===========================================================

        WRITE(fnamedatanc,'(2A)') ihopObsFile(1:IL),'.nc'
        err = NF_OPEN(fnamedatanc, 0, fiddata_obs(num_file))

! Read number of observations:
        ncid = fiddata_obs(num_file)
        err = NF_INQ_DIMID(ncid,'iOBS', dimid1 )
        err = NF_INQ_DIMLEN(ncid, dimid1, obsNo(num_file) )
      
        WRITE(msgbuf,'(A,I9)')
     &     '  # of observations in file                   =',
     &     obsNo(num_file)
        CALL print_message(
     &     msgbuf, standardmessageunit, SQUEEZE_RIGHT, myThid)

! Read meta info for observations:
        err = NF_INQ_VARID(ncid,'ihop_YYYYMMDD', varid1a )
        err = err + NF_INQ_VARID(ncid,'ihop_HHMMSS', varid1b )
        err = err + NF_INQ_VARID(ncid,'ihop_iter', varid2 )

        IF (err.NE.NF_NOERR) THEN
         IL  = ILNBLNK( ihopObsFile )
         WRITE(msgBuf,'(3A)')
     &    'IHOP_COST_INIT_FIXED: file ', ihopObsFile(1:IL),
     &    '.nc is not in the pkg/ihop cost format: check iter,uncert'
         CALL PRINT_ERROR( msgBuf, myThid )
        ENDIF
! Assume all ihop observations are instantaneous

! Weight for each ihop obs (if missing, assume equal weights):
        weighIhopObs=1
        err = NF_INQ_VARID(ncid,'ihop_uncert', varid3 )
        IF (err.NE.NF_NOERR) THEN
         weighIhopObs=0

         WRITE(msgbuf,'(2A)')
     &     '      ihop obs file does not have uncertainties;',
     &     '      assume all obs datum are weighed evenly '
         CALL print_message(
     &      msgbuf, standardmessageunit, SQUEEZE_RIGHT, myThid)

        ENDIF

! Read the ihop observations' positions:
        err =       NF_INQ_VARID(ncid,'ihop_x', varid4a )
        err = err + NF_INQ_VARID(ncid,'ihop_y', varid4b )
        err = err + NF_INQ_VARID(ncid,'ihop_z', varid4c )
        IF (err.NE.NF_NOERR) THEN
         IL  = ILNBLNK( ihopObsFile )
         WRITE(msgBuf,'(3A)')
     &    'IHOP_COST_INIT_FIXED: file ', ihopObsFile(1:IL),
     &    '.nc is not in the pkg/ihop cost format, check x,y,z'
         CALL PRINT_ERROR( msgBuf, myThid)
        ENDIF

! Read the ihop observations' values:
        err =       NF_INQ_VARID(ncid,'ihop_tau', varid5a )
!        err = err + NF_INQ_VARID(ncid,'ihop_ang', varid5b )
        IF (err.NE.NF_NOERR) THEN
         IL  = ILNBLNK( ihopObsFile )
         WRITE(msgBuf,'(3A)')
     &    'IHOP_COST_INIT_FIXED: file ', ihopObsFile(1:IL),
     &    '.nc is not in the pkg/ihop cost format, check tau,ang'
         CALL PRINT_ERROR( msgBuf, myThid)
        ENDIF

!===========================================================
!       SET IHOP EQUIV NETCDF VARIABLES
!===========================================================

! Variables names for model equivalent files
        ihop_nameequi = 'mod_val'
        ihop_namemask = 'mod_mask'
        ihop_nameval  = 'obs_val'
        ihop_nameuncert = 'obs_uncert'

! Create new netcdf global file for ihop model equivalent 
        IF ( myProcId .eq. 0 ) THEN
         JL  = ILNBLNK( ihopObsDir )

         WRITE(fnameequincglo,'(3A)')
     &    ihopObsDir(1:JL),ihopObsFile(1:IL),'.equi.nc'
         WRITE(adfnameequincglo,'(4A)')
     &    ihopObsDir(1:JL),'ad',
     &    ihopObsFile(1:IL),'.equi.nc'

         inquire( file=fnameequincglo, exist=exst )
         IF (.NOT.exst) THEN

          err = NF_CREATE(fnameequincglo,NF_CLOBBER,fidglobal(num_file))
          err = NF_DEF_DIM(fidglobal(num_file),'iOBS',ObsNo(num_file),
     &           dimid1)

          err = NF_DEF_VAR(fidglobal(num_file),
     &         'mod_val',NF_DOUBLE,
     &         1,dimid1,varidEquiv(1))

          err = NF_PUT_ATT_DOUBLE(fidglobal(num_file),
     &         varidEquiv(1),'_FillValue',
     &         NF_DOUBLE,1, -9999. _d 0 )

          err = NF_DEF_VAR(fidglobal(num_file),
     &         'mod_mask',NF_DOUBLE,
     &         1,dimid1,varidEquiv(2))

          err = NF_PUT_ATT_DOUBLE(fidglobal(num_file),
     &         varidEquiv(2),'_FillValue',
     &         NF_DOUBLE,1, -9999. _d 0 )

          err = NF_ENDDEF(fidglobal(num_file))
          err = NF_CLOSE(fidglobal(num_file))
          err = NF_OPEN(fnameequincglo, NF_WRITE, fidglobal(num_file))

         ELSE
          err = NF_OPEN(fnameequincglo, NF_WRITE, fidglobal(num_file))
         ENDIF

#ifdef ALLOW_ADJOINT_RUN
         inquire( file=adfnameequincglo, exist=exst )
         IF (.NOT.exst) THEN

          err = NF_CREATE(adfnameequincglo,NF_CLOBBER,
     &           fidadglobal(num_file))
          err = NF_DEF_DIM(fidadglobal(num_file),'iOBS',
     &           ObsNo(num_file),dimid1)

          err = NF_DEF_VAR(fidadglobal(num_file),
     &         'mod_val',NF_DOUBLE,
     &         1,dimid1,varidEquiv(1))

          err = NF_PUT_ATT_DOUBLE(fidadglobal(num_file),
     &         varidEquiv(1),'_FillValue',
     &         NF_DOUBLE,1, 0. _d 0 )

          err = NF_DEF_VAR(fidglobal(num_file),
     &         'mod_mask',NF_DOUBLE,
     &         1,dimid1,varidEquiv(2))

          err = NF_PUT_ATT_DOUBLE(fidglobal(num_file),
     &         varidEquiv(2),'_FillValue',
     &         NF_DOUBLE,1, 0. _d 0 )

          err = NF_ENDDEF(fidadglobal(num_file))
          err = NF_CLOSE(fidadglobal(num_file))
          err = NF_OPEN(adfnameequincglo, NF_WRITE, 
     &         fidadglobal(num_file))

         ELSE
          err = NF_OPEN(adfnameequincglo, NF_WRITE, 
     &         fidadglobal(num_file))
         ENDIF
#endif
        ENDIF

!===========================================================
!       READ IHOP OBS NETCDF VARIABLES
!===========================================================
! NOTE: not applying reading in chunks
!       see pkg/profiles or pkg/obsfit for how to chunk.

! default values
        DO k=1,NOBSMAX_IHOP
         ihopObs_ind_glob(num_file,k) = 0
        ENDDO
        ObsNo_valid = 0

! read time
        err = NF_GET_VARA_DOUBLE(ncid,varid1a,1,
     &   ObsNo(num_file), tmpyymmdd)
        err = NF_GET_VARA_DOUBLE(ncid,varid1b,1,
     &   ObsNo(num_file), tmphhmmss)

        IF (err.NE.NF_NOERR) THEN
         WRITE(msgBuf,'(3A)')
     &   'IHOP_COST_INIT_FIXED: file ', ihopObsFile(1:IL),
     &   '.nc was not read properly: check yyyymmdd, hhmmss.'
         CALL PRINT_ERROR( msgBuf, myThid)
        ENDIF

! loop through number of observations, in obsfit this is a 'chunk'
        DO m=1,ObsNo(num_file)

!AV if obs starts before model run, will be ignored even if it ends after 
!AV model starts - ok for now      
         IF (( ( tmpyymmdd(m).GT.yymmddMin ).OR.(( tmpyymmdd(m).EQ.
     &     yymmddMin ).AND.( tmphhmmss(m).GE.hhmmssMin ))).AND.
     &   ( ( tmpyymmdd(m).LT.yymmddMax ).OR.(( tmpyymmdd(m).EQ.
     &     yymmddMax ).AND.( tmphhmmss(m).LE.hhmmssMax ))) ) THEN
          hh = int(tmphhmmss(m))/10000
          IF ( hh.LT.hoursPerDay ) THEN
           obsIsInRunTime=1
#ifdef ALLOW_CAL
           CALL cal_FullDate( int(tmpyymmdd(m)),int(tmphhmmss(m)),
     &      tmpdate,myThid )
           CALL cal_TimePassed( modelstartdate,tmpdate,tmpdiff,
     &      myThid )
           CALL cal_ToSeconds (tmpdiff,diffsecs,myThid)
#endif /*ALLOW_CAL*/
           diffsecs=diffsecs+nIter0*deltaTclock
          ELSE
! if tmphhmmss is out of range then disregard obs
           obsIsInRunTime=0
           diffsecs=-deltaTclock
           ObsNo_hh=ObsNo_hh+1
          ENDIF
         ELSE
          obsIsInRunTime=0
          diffsecs=-deltaTclock
         ENDIF

         IF (obsIsInRunTime.EQ.1) THEN
! if yes then store obs index position:
          ObsNo_valid = ObsNo_valid+1
          ihopObs_ind_glob(num_file,ObsNo_valid) = m

         ENDIF !if (obsIsInRunTime.EQ.1) then

! regardless of whether obs is in run time or not, 
! keep track of ihop obs index position
! and assign time to ihop obs
         tmp_ihopObs_time(m)=diffsecs

! check that maximum size was not reached:
         IF (ObsNo(num_file).GT.NOBSMAX_IHOP) THEN
          WRITE(msgBuf,'(3A)')
     &     'IHOP_COST_INIT_FIXED: file ', ihopObsFile(1:IL),
     &     '.nc was not read properly: increase NOBSMAX_IHOP.'
          CALL PRINT_ERROR( msgBuf, myThid)
         ENDIF

        ENDDO !m

! calibrate number of VALID obs datum
        ObsNo(num_file) = ObsNo_valid

! output ihop cost init info
        WRITE(msgbuf,'(A,I9)')
     &     '  # of ihop obs with erroneous HHMMSS values =',
     &     ObsNo_hh
        CALL print_message(
     &     msgbuf, standardmessageunit, SQUEEZE_RIGHT, myThid)

        WRITE(msgbuf,'(A,I9)')
     &     '  # of ihop obs within time period  =',
     &     ObsNo(num_file)
        CALL print_message(
     &     msgbuf, standardmessageunit, SQUEEZE_RIGHT, myThid)

! ==============================================================================
! Initiate tiled defaults for ihopObs before reading in from netcdf
        DO bj=1,nSy
         DO bi=1,nSx
          ObsNo_valid = 0
          ObsNo_tiled(num_file,bi,bj)=ObsNo(num_file)

          DO k=1,NOBSMAX_IHOP
           ihopObs_time(num_file,k,bi,bj) =    -999. _d 0
           ihopObs_lon(num_file,k,bi,bj)=      -999. _d 0
           ihopObs_lat(num_file,k,bi,bj)=      -999. _d 0
           ihopObs_depth(num_file,k,bi,bj)=    -999. _d 0
           ihopObs_uncert(num_file,k,bi,bj)=      1. _d 0
           ihopObs_ind_glob_tiled(num_file,k,bi,bj)=    0
          ENDDO !k

! Read ALL info using respective varid and save to tmp, NO CHUNKS!
          IF (weighIhopObs.EQ.1) THEN  
           err = NF_GET_VARA_DOUBLE(ncid,varid3,1,
     &      ObsNo_tiled(num_file,bi,bj), tmp_uncert2)
          ENDIF
          err = NF_GET_VARA_DOUBLE(ncid,varid4a,1,
     &     ObsNo_tiled(num_file,bi,bj), tmp_lon2)
          err = NF_GET_VARA_DOUBLE(ncid,varid4b,1,
     &     ObsNo_tiled(num_file,bi,bj), tmp_lat2)
          err = NF_GET_VARA_DOUBLE(ncid,varid4c,1,
     &     ObsNo_tiled(num_file,bi,bj), tmp_depth2)

          IF (err.NE.NF_NOERR) THEN
           WRITE(msgBuf,'(3A)')
     &      'IHOP_COST_INIT_FIXED: file ', ihopObsFile(1:IL),
     &      '.nc was not read properly: check x,y,z,uncert.'
           CALL PRINT_ERROR(msgBuf, myThid)
          ENDIF

! Find obs datum in MITgcm spatial grid
! loop through observations, NOTE: in obsfit this is a 'chunk'
          DO m=1,ObsNo_tiled(num_file,bi,bj)
           ihopObsInTile=1
      
! in a valid time?
           IF ((tmp_ihopObs_time(m).LT.0. _d 0) .OR.
     &        (tmp_ihopObs_time(m).GT. modelend-modelstart)) THEN
            ihopObsInTile=0
           ENDIF

! Determine whether obs datum is in current tile domain
           lon_cur=tmp_lon2(m)
           lat_cur=tmp_lat2(m)
           depth_cur=tmp_depth2(m)

           IF (ihopObsInTile.EQ.1) THEN
            IF (xC(sNx+1,1,bi,bj).LT.xC(1,1,bi,bj)) THEN
             tmp_lon=xC(sNx+1,1,bi,bj)+360. _d 0
            ELSE
             tmp_lon=xC(sNx+1,1,bi,bj)
            ENDIF

            IF ((xC(1,1,bi,bj).LE.lon_cur).AND.
     &       (tmp_lon.GT.lon_cur).AND.
     &       (yC(1,1,bi,bj).LE.lat_cur).AND.
     &       (yC(1,sNy+1,bi,bj).GT.lat_cur) ) THEN
!  do nothing 
            ELSEIF ((xC(sNx+1,1,bi,bj).LT.xC(1,1,bi,bj)).AND.
     &       (xC(1,1,bi,bj).LE.lon_cur+360. _d 0).AND.
     &       (tmp_lon.GT.lon_cur+360. _d 0).AND.
     &       (yC(1,1,bi,bj).LE.lat_cur).AND.
     &       (yC(1,sNy+1,bi,bj).GT.lat_cur)
     &       ) THEN
             lon_cur=lon_cur+360. _d 0         
            ELSE
!  not in this tile
             ihopObsInTile=0
            ENDIF
 
! Determine i,j,k to the south-west-down of obs datum
! set default values
            ihopObs_i=-10
            ihopObs_j=-10
            ihopObs_k1=-10
            ihopObs_k2=-10
            lon_1=-10
            lon_2=-10
            lat_1=-10
            lat_2=-10
            depth_1=-10
            depth_2=-10

            IF (ihopObsInTile.EQ.1) THEN
! Value of k, below the data point:
! case 1: above first depth level
             IF (-rC(1).GT.depth_cur) THEN
              ihopObs_k1=1
              ihopObs_k2=1
              depth_fac=1. _d 0
             ENDIF
! case 2: below last depth level
             IF (-rC(Nr).LE.depth_cur) THEN
              ihopObs_k1=Nr
              ihopObs_k2=Nr
              depth_fac=1. _d 0
             ENDIF
! case 3: between two depth levels
             DO k=1,Nr-1
              IF ((-rC(k).LE.depth_cur).AND.
     &          (-rC(k+1).GT.depth_cur)) THEN
               ihopObs_k1=k
               ihopObs_k2=k+1
               depth_1=-rC(k)
               depth_2=-rC(k+1)
               depth_fac=(depth_cur-depth_1)/(depth_2-depth_1)
              ENDIF
             ENDDO

             DO j=1,sNy+1
              DO i=1,sNx+1

               lon_tmp1=xC(i,j,bi,bj)
               lon_tmp2=xC(i+1,j,bi,bj)
               lat_tmp1=yC(i,j,bi,bj)
               lat_tmp2=yC(i,j+1,bi,bj)
          
! value of j, south of the obs datum:
               IF ((lat_tmp1.LE.lat_cur).AND.
     &          (lat_tmp2.GT.lat_cur)) THEN
                ihopObs_j=j
                lat_1=lat_tmp1
                lat_2=lat_tmp2
               ENDIF

! value of i, west of the obs datum:
               IF (xC(i+1,j,bi,bj).LT.xC(1,j,bi,bj)) THEN
                lon_tmp2=lon_tmp2+360
               ENDIF
               IF (xC(i,j,bi,bj).LT.xC(1,j,bi,bj)) THEN
                lon_tmp1=lon_tmp1+360
               ENDIF

               IF ((lon_tmp1.LE.lon_cur).AND.
     &          (lon_tmp2.GT.lon_cur)) THEN
                ihopObs_i=i
                lon_1=lon_tmp1
                lon_2=lon_tmp2
               ENDIF

              ENDDO !j
             ENDDO !i
            ENDIF !ihopObsInTile.EQ.1 2nd time

            IF ( ((ihopObs_i.EQ.-10).OR.(ihopObs_j.EQ.-10)) .OR.
     &           ((ihopObs_k1.EQ.-10).OR.(ihopObs_k2.EQ.-10)) ) THEN
             ihopObsInTile=0
            ENDIF

! 'Reshape' observed datum to tiled structures
            IF (ihopObsInTile.EQ.1) THEN
             ObsNo_valid = ObsNo_valid+1
             ihopObs_ind_glob_tiled(num_file,ObsNo_valid,bi,bj) = m

             ihopObs_time(num_file,ObsNo_valid,bi,bj) =
     &                                          tmp_ihopObs_time(m)
             IF (weighIhopObs.EQ.1) THEN  
              ihopObs_uncert(num_file,ObsNo_valid,bi,bj)=tmp_uncert2(m)
             ENDIF
             ihopObs_lon(num_file,ObsNo_valid,bi,bj)=lon_cur
             ihopObs_lat(num_file,ObsNo_valid,bi,bj)=lat_cur
             ihopObs_depth(num_file,ObsNo_valid,bi,bj)=depth_cur

             ihopObs_i_tiled(num_file,ObsNo_valid,bi,bj) = ihopObs_i
             ihopObs_j_tiled(num_file,ObsNo_valid,bi,bj) = ihopObs_j
             ihopObs_k_tiled(num_file,ObsNo_valid,bi,bj) = ihopObs_k1

            ENDIF !ihopObsInTile.EQ.1 3rd time

           ENDIF !ihopObsInTile.EQ.1 1st time

! ==============================================================================
! Check that maximum size was not reached:
           IF (m.GE.NOBSMAX_IHOP) then
            WRITE(msgBuf,'(3A)')
     &       'IHOP_COST_INIT_FIXED: file ', ihopObsFile(1:IL),
     &       '.nc was not read properly: increase NOBSMAX_IHOP.'
            CALL PRINT_ERROR( msgBuf, myThid)
           ENDIF

          ENDDO !m

! Number of valid ihop obs in the file
          ObsNo_tiled(num_file,bi,bj) = ObsNo_valid
       
          WRITE(msgbuf,'(A,I4,A,I4)')
     &       '  current tile is bi,bj                      =',
     &       bi,',',bj
          CALL print_message(
     &       msgbuf, standardmessageunit, SQUEEZE_RIGHT, myThid)

          WRITE(msgbuf,'(A,I9)')
     &       '  # of ihop obs within tile and time period  =',
     &       ObsNo_tiled(num_file,bi,bj)
          CALL print_message(
     &       msgbuf, standardmessageunit, SQUEEZE_RIGHT, myThid)

!===========================================================
! Create files for model counterparts to ihop obs
!===========================================================

          IF (ObsNo_tiled(num_file,bi,bj).GT.0) THEN
           iG=bi+(myXGlobalLo-1)/sNx ! Kludge until unstructered tiles
           jG=bj+(myYGlobalLo-1)/sNy ! Kludge until unstructered tiles

           JL  = ILNBLNK( ihopObsDir )

!IE:           IF (ihopObsDoNcOutput) THEN
! ALWAYS Write to nc file

           WRITE(fnameequinc,'(3A,I3.3,A,I3.3,A)')
     &      ihopObsDir(1:JL),ihopObsFile(1:IL),'.',iG,'.',jG,'.equi.nc'
           WRITE(adfnameequinc,'(4A,I3.3,A,I3.3,A)')
     &      ihopObsDir(1:JL),'ad',
     &      ihopObsFile(1:IL),'.',iG,'.',jG,'.equi.nc'

           inquire( file=fnameequinc, exist=exst )
           IF (.NOT.exst) THEN
      
            CALL ihop_cost_init_equifiles(num_file,
     &       fiddata_obs(num_file),fnameequinc,
     &       fidfwd_obs(num_file,bi,bj),ObsNo_tiled(num_file,bi,bj),
     &       bi,bj,myThid)
     
           ELSE
      
            err = NF_OPEN(fnameequinc,NF_WRITE,
     &        fidfwd_obs(num_file,bi,bj))
           ENDIF

#ifdef ALLOW_ADJOINT_RUN
           inquire( file=adfnameequinc, exist=exst )
           IF (.NOT.exst) THEN
            CALL ihop_cost_init_equifiles(num_file,
     &       fiddata_obs(num_file),
     &       adfnameequinc, fidadj_obs(num_file,bi,bj),
     &       ObsNo_tiled(num_file,bi,bj),
     &       bi,bj,myThid)
           ELSE
            err = NF_OPEN(adfnameequinc,NF_WRITE,
     &       fidadj_obs(num_file,bi,bj))
           ENDIF
#endif

          ENDIF

         ENDDO !bj
        ENDDO !bi
!===========================================================
       ELSE !if (IL.NE.0) then

        DO bj=1,nSy
         DO bi=1,nSx
          ObsNo(num_file)=0
          ObsNo_tiled(num_file,bi,bj)=0
          DO k=1,NOBSMAX_IHOP
           ihopObs_lon(num_file,k,bi,bj)=      -999. _d 0
           ihopObs_lat(num_file,k,bi,bj)=      -999. _d 0
           ihopObs_depth(num_file,k,bi,bj)=    -999. _d 0
           ihopObs_uncert(num_file,k,bi,bj)=      1. _d 0
           ihopObs_ind_glob_tiled(num_file,k,bi,bj)=    0
          ENDDO
         ENDDO !bj
        ENDDO !bi

       ENDIF !if (IL.NE.0) then
      ENDDO !      do num_file=1,NFILESMAX_IHOP
      
      _END_MASTER( myThid )
      _BARRIER

!===========================================================
      WRITE(msgbuf,'(A)') ' '
      CALL print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT, myThid)
      WRITE(msgbuf,'(A)')
     &'// ======================================================='
      CALL print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT, myThid)
      WRITE(msgbuf,'(A)')
     &'// ihop cost read obs data file >>> END <<<'
      CALL print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT, myThid)
      WRITE(msgbuf,'(A)')
     &'// ======================================================='
      CALL print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT, myThid)
      WRITE(msgbuf,'(A)') ' '
      CALL print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT, myThid)
      
#endif /* ALLOW_COST && ALLOW_IHOP*/

      RETURN
      END
