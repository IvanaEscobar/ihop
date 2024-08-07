#include "IHOP_OPTIONS.h"
!BOP
! !ROUTINE: BELLHOP
! !INTERFACE:
MODULE BELLHOP
  ! Written as a module to be used in ihop_*.F parts of the MITgcm package
  ! BELLHOP Beam tracing for ocean acoustics

  ! Copyright (C) 2009 Michael B. Porter

  ! This program is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.

  ! You should have received a copy of the GNU General Public License
  ! along with this program.  If not, see <http://www.gnu.org/licenses/>.

  ! First version (1983) originally developed with Homer Bucker, Naval Ocean 
  ! Systems Center

  
  USE ihop_mod,     only:   rad2deg, i, MaxN, Title, Beam, ray2D, istep,       &
                            NRz_per_range, afreq, SrcDeclAngle,                &
                            PRTFile, SHDFile, ARRFile, RAYFile, DELFile   
  USE initenvihop,  only:   initEnv, openOutputFiles, resetMemory
  USE angle_mod,    only:   Angles, ialpha
  USE srPos_mod,    only:   Pos
  USE ssp_mod,      only:   EvaluateSSP, HSInfo, Bdry, SSP, fT
  USE bdry_mod,     only:   initATI, initBTY, GetTopSeg, GetBotSeg, Bot, Top,  &
                            atiType, btyType, NatiPts, NbtyPts, iSmallStepCtr, &
                            IsegTop, IsegBot, rTopSeg, rBotSeg
  USE refCoef,      only:   readReflectionCoefficient,                         &
                            InterpolateReflectionCoefficient, ReflectionCoef,  &
                            RTop, RBot, NBotPts, NTopPts
  USE influence,    only:   InfluenceGeoHatRayCen,&! InfluenceSGB,             &
                            InfluenceGeoGaussianCart, InfluenceGeoHatCart,     &
                            ScalePressure
  USE beamPattern 
  USE writeRay,     only:   WriteRay2D, WriteDel2D
  USE arr_mod,      only:   WriteArrivalsASCII,WriteArrivalsBinary,MaxNArr,    &
                            Arr, NArr, U
  

!   !USES:
  IMPLICIT NONE
!   == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef ALLOW_CTRL
# include "CTRL_FIELDS.h"
#endif

!   == External Functions ==
    INTEGER  ILNBLNK
    EXTERNAL ILNBLNK

CONTAINS
  SUBROUTINE IHOP_INIT ( myTime, myIter, myThid )
  !     !INPUT/OUTPUT PARAMETERS:
  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    _RL, INTENT( IN )       :: myTime
    INTEGER, INTENT( IN )   :: myIter
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    LOGICAL, PARAMETER  :: Inline = .FALSE.
    INTEGER             :: iostat, iAllocStat, ierr
    INTEGER             :: jj 
    REAL                :: Tstart, Tstop
  ! added locally previously read in from unknown mod ... IEsco22
    CHARACTER ( LEN=2 ) :: AttenUnit
  ! For MPI writing: copying eeboot_minimal.F
    CHARACTER*(MAX_LEN_FNAM)    :: fNam
    CHARACTER*(6)               :: fmtStr
    INTEGER                     :: mpiRC, IL
  ! ===========================================================================
    INTEGER              :: iSeg
    INTEGER, PARAMETER   :: ArrivalsStorage = 2000, MinNArr = 10
  ! ===========================================================================
 
    ! Open the print file: template from eeboot_minimal.F
#ifdef IHOP_WRITE_OUT
    IF ( .NOT.usingMPI ) THEN
        WRITE(myProcessStr, '(I10.10)') myIter
        IL=ILNBLNK( myProcessStr )
        WRITE(fNam,'(4A)') TRIM(IHOP_fileroot),'.',myProcessStr(1:IL),'.prt'
        IF ( IHOP_dumpfreq.GE.0 ) &
         OPEN( PRTFile, FILE = fNam, STATUS = 'UNKNOWN', IOSTAT = iostat )
#ifdef ALLOW_USE_MPI
    ELSE ! using MPI
        CALL MPI_COMM_RANK( MPI_COMM_MODEL, mpiMyId, mpiRC )
        myProcId = mpiMyId
        IL = MAX(4,1+INT(LOG10(DFLOAT(nPx*nPy))))
        WRITE(fmtStr,'(2(A,I1),A)') '(I',IL,'.',IL,')'
        WRITE(myProcessStr,fmtStr) myProcId
        IL = ILNBLNK( myProcessStr )
        mpiPidIo = myProcId
        pidIO    = mpiPidIo

        IF( mpiPidIo.EQ.myProcId ) THEN
#  ifdef SINGLE_DISK_IO
         IF( myProcId.eq.0 ) THEN
#  endif
            IF (myIter.GE.0) THEN
                WRITE(fNam,'(4A,I10.10,A)') &
                    TRIM(IHOP_fileroot),'.',myProcessStr(1:IL),'.',myIter,'.prt'
            ELSE
                WRITE(fNam,'(4A)') &
                    TRIM(IHOP_fileroot),'.',myProcessStr(1:IL),'.prt'
            ENDIF

            IF ( IHOP_dumpfreq .GE. 0) &
             OPEN(PRTFile, FILE=fNam, STATUS='UNKNOWN', IOSTAT=iostat )
            IF ( iostat /= 0 ) THEN
                WRITE(*,*) 'ihop: IHOP_fileroot not recognized, ', &
                    TRIM(IHOP_fileroot)
                WRITE(msgBuf,'(A)') 'IHOP_INIT: Unable to recognize file'
                CALL PRINT_ERROR( msgBuf, myThid )
                STOP 'ABNORMAL END: S/R IHOP_INIT'
            END IF
#  ifdef SINGLE_DISK_IO
         END IF
#  endif
        END IF
# endif /* ALLOW_USE_MPI */
    END IF
#endif /* IHOP_WRITE_OUT */

    ! Reset memory
    CALL resetMemory() 
    ! save data.ihop, gcm SSP: REQUIRED
    CALL initEnv( myTime, myIter, myThid )
    ! AlTImetry: OPTIONAL, default is no ATIFile
    CALL initATI( Bdry%Top%HS%Opt( 5:5 ), Bdry%Top%HS%Depth, myThid )
    ! BaThYmetry: OPTIONAL, default is BTYFile
    CALL initBTY( Bdry%Bot%HS%Opt( 2:2 ), Bdry%Bot%HS%Depth, myThid )
    ! (top and bottom): OPTIONAL
    CALL readReflectionCoefficient( Bdry%Bot%HS%Opt( 1:1 ), &
                                    Bdry%Top%HS%Opt( 2:2 ), myThid ) 
    ! Source Beam Pattern: OPTIONAL, default is omni source pattern
    SBPFlag = Beam%RunType( 3:3 )
    CALL readPat( myThid )
    Pos%Ntheta = 1
    ALLOCATE( Pos%theta( Pos%Ntheta ), Stat = IAllocStat )
    IF ( IAllocStat/=0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'BELLHOP IHOP_INIT: failed allocation Pos%theta'
        CALL PRINT_ERROR( msgBuf, myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R  IHOP_INIT'
    ENDIF
    Pos%theta( 1 ) = 0.
  

! Allocate arrival and U variables on all MPI processes 
    SELECT CASE ( Beam%RunType( 5:5 ) )
    CASE ( 'I' )
       NRz_per_range = 1         ! irregular grid
    CASE DEFAULT
       NRz_per_range = Pos%NRz   ! rectilinear grid
    END SELECT
  
    IF ( ALLOCATED( U ) ) DEALLOCATE( U )
     SELECT CASE ( Beam%RunType( 1:1 ) )
     ! for a TL calculation, allocate space for the pressure matrix
     CASE ( 'C', 'S', 'I' )        ! TL calculation
          ALLOCATE ( U( NRz_per_range, Pos%NRr ), Stat = iAllocStat )
          IF ( iAllocStat/=0 ) THEN
#ifdef IHOP_WRITE_OUT
              WRITE(msgBuf,'(2A)') 'BELLHOP IHOP_INIT: ', & 
                             'Insufficient memory for TL matrix: reduce Nr*NRz'
              CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
              STOP 'ABNORMAL END: S/R IHOP_INIT'
          END IF
          U = 0.0
     CASE ( 'A', 'a', 'R', 'E', 'e' )   ! Arrivals calculation
          ALLOCATE ( U( 1,1 ), Stat = iAllocStat )   ! open a dummy variable
     CASE DEFAULT
          ALLOCATE ( U( 1,1 ), Stat = iAllocStat )   ! open a dummy variable
     END SELECT
  
     ! for an arrivals run, allocate space for arrivals matrices
     SELECT CASE ( Beam%RunType( 1:1 ) )
     CASE ( 'A', 'a', 'e' )
          ! allow space for at least MinNArr arrivals
          MaxNArr = MAX( ArrivalsStorage / ( NRz_per_range * Pos%NRr ), & 
                         MinNArr )
          ALLOCATE ( Arr( MaxNArr, Pos%NRr, NRz_per_range ), &
                     NArr( Pos%NRr, NRz_per_range ), Stat = iAllocStat )
          IF ( iAllocStat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
              WRITE(msgBuf,'(2A)') 'BELLHOP IHOP_INIT: ', & 
               'Not enough allocation for Arr; reduce ArrivalsStorage'
              CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
              STOP 'ABNORMAL END: S/R IHOP_INIT'
          END IF
     CASE DEFAULT
          MaxNArr = 1
          ALLOCATE ( Arr( 1, NRz_per_range, Pos%NRr ), &
                     NArr( Pos%NRr, NRz_per_range ), Stat = iAllocStat )
     END SELECT
  
     NArr( 1:Pos%NRr, 1:NRz_per_range ) = 0 ! IEsco22 unnecessary? NArr = 0 below
  
#ifdef IHOP_WRITE_OUT
     WRITE(msgBuf,'(A)') 
     ! In adjoint mode we do not write output besides on the first run
     IF (IHOP_dumpfreq.GE.0) &
       CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
#endif /* IHOP_WRITE_OUT */



! open all output files
    IF ( IHOP_dumpfreq .GE. 0 ) &
     CALL OpenOutputFiles( IHOP_fileroot, myTime, myIter, myThid )
  
    ! Run Bellhop solver on a single processor
    if (numberOfProcs.gt.1) then
! Use same single processID as IHOP COST package
!        if(myProcId.eq.(numberOfProcs-1)) then
        if(myProcId.eq.0) then
            CALL CPU_TIME( Tstart )
            CALL BellhopCore(myThid)
            CALL CPU_TIME( Tstop )
! Alternitavely, we can broadcast relevant info to all mpi processes Ask P. 
!#ifdef ALLOW_COST
!            ! Broadcast info to all MPI procs for COST function accumulation
!            print *, "escobar: broacasting from pid ", myProcId
!            CALL MPI_BCAST(i, 1, MPI_COMPLEX, myProcId, MPI_COMM_MODEL, ierr)
!
!#endif /* ALLOW_COST */
        endif
    else
        CALL CPU_TIME( Tstart )
        CALL BellhopCore(myThid)
        CALL CPU_TIME( Tstop )
    endif
  
#ifdef IHOP_WRITE_OUT
    IF ( IHOP_dumpfreq.GE.0 ) THEN
        ! print run time
        if (numberOfProcs.gt.1) then
            if(myProcId.ne.(numberOfProcs-1)) then
                WRITE(msgBuf,'(A,I4,A)') 'NOTE: Proc ',myProcId, &
                    " didn't run ihop"
                CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
            endif
        endif
        WRITE(msgBuf, '(A)' )
        CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
        WRITE(msgBuf, '(A,G15.3,A)' ) 'CPU Time = ', Tstop-Tstart, 's'
        CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
  
        ! close all files
        IF ( IHOP_dumpfreq .GE. 0) THEN
            SELECT CASE ( Beam%RunType( 1:1 ) )
            CASE ( 'C', 'S', 'I' )  ! TL calculation
               CLOSE( SHDFile )
            CASE ( 'A', 'a' )       ! arrivals calculation
               CLOSE( ARRFile )
            CASE ( 'R', 'E' )       ! ray and eigen ray trace
               CLOSE( RAYFile )
            CASE ( 'e' )
               CLOSE( RAYFile ) 
               CLOSE( ARRFile )
               IF ( writeDelay ) CLOSE( DELFile )
            END SELECT
  
            if (numberOfProcs.gt.1) then
                ! Erase prtfiles that aren't on procid = 0
                if(myProcId.ne.0) then
                    CLOSE(PRTFile, STATUS='DELETE')
                else
                    CLOSE(PRTFile)
                endif
            else
                CLOSE(PRTFile)
            endif
        ENDIF
    ENDIF
#endif /* IHOP_WRITE_OUT */
  
  RETURN
  END !SUBROUTINE IHOP_INIT
  
  ! **********************************************************************!
  SUBROUTINE BellhopCore( myThid )
  
  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    INTEGER              :: IBPvec( 1 ), ibp, is, iBeamWindow2, Irz1, Irec, &
                            NalphaOpt
    REAL    (KIND=_RL90) :: Amp0, DalphaOpt, xs( 2 ), RadMax, s, &
                            c, cimag, gradc( 2 ), crr, crz, czz, rho
  
    afreq = 2.0 * PI * IHOP_freq
  
    Angles%alpha  = Angles%alpha * deg2rad  ! convert to radians
    Angles%Dalpha = 0.0
    IF ( Angles%Nalpha > 1 ) THEN
         Angles%Dalpha = ( Angles%alpha( Angles%Nalpha ) - Angles%alpha( 1 ) ) &
                         / ( Angles%Nalpha - 1 )  ! angular spacing between beams
    ELSE
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'BELLHOP BellhopCore: ', & 
                      'Required: Nalpha>1, else add iSingle_alpha(see angleMod)'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R BellhopCore'
    END IF
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !         begin solve         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SourceDepth: DO is = 1, Pos%NSz
       xs = [ zeroRL, Pos%Sz( is ) ]   ! source coordinate, assuming source @ r=0
  
       SELECT CASE ( Beam%RunType( 1:1 ) )
       CASE ( 'C','S','I' ) ! TL calculation, zero out pressure matrix
          U = 0.0
       CASE ( 'A','a','e' )   ! Arrivals calculation, zero out arrival matrix
          NArr = 0
       END SELECT
  
       CALL EvaluateSSP(  xs, c, cimag, gradc, crr, crz, czz, rho, &
                          'TAB', myThid  )
  
       !!IESCO22: BEAM stuff !!
       RadMax = 5 * c / IHOP_freq  ! 5 wavelength max radius IEsco22: unused
       IF ( Beam%RunType( 1:1 ) == 'C' ) THEN ! for Coherent TL Run
       ! Are there enough rays?
          DalphaOpt = SQRT( c / ( 6.0 * IHOP_freq * Pos%Rr( Pos%NRr ) ) )
          NalphaOpt = 2 + INT( ( Angles%alpha( Angles%Nalpha ) &
                               - Angles%alpha( 1 ) ) / DalphaOpt )
#ifdef IHOP_WRITE_OUT
          IF ( Angles%Nalpha < NalphaOpt ) THEN
             WRITE( msgBuf, '(A,/,A,I10.4)' ) 'WARNING: Too few beams',&
                 'Nalpha should be at least = ', NalphaOpt
             ! In adjoint mode we do not write output besides on the first run
             IF (IHOP_dumpfreq.GE.0) &
                CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
          ENDIF
#endif /* IHOP_WRITE_OUT */
       ENDIF
       !!IESCO22: end BEAM stuff !!
  
       ! Trace successive beams
       DeclinationAngle: DO ialpha = 1, Angles%Nalpha
          ! take-off declination angle in degrees
          SrcDeclAngle = rad2deg * Angles%alpha( ialpha )
  
          ! Single ray run? then don't visit code below
          IF ( Angles%iSingle_alpha==0 .OR. ialpha==Angles%iSingle_alpha ) THEN
  
             !!IESCO22: BEAM stuff !!
             IBPvec = maxloc( SrcBmPat( :, 1 ), mask = SrcBmPat( :, 1 ) &
                      < SrcDeclAngle )  ! index of ray angle in beam pattern
             IBP    = IBPvec( 1 )
             IBP    = MAX( IBP, 1 )           ! don't go before beginning of table
             IBP    = MIN( IBP, NSBPPts - 1 ) ! don't go past end of table
             ! IEsco22: When a beam pattern isn't specified, IBP = 1
  
             ! linear interpolation to get amplitude
             s    = ( SrcDeclAngle  - SrcBmPat( IBP, 1 ) ) &
                    / ( SrcBmPat( IBP + 1, 1 ) - SrcBmPat( IBP, 1 ) )
             Amp0 = ( 1 - s ) * SrcBmPat( IBP, 2 ) + s * SrcBmPat( IBP + 1, 2 )
             ! IEsco22: When a beam pattern isn't specified, Amp0 = 0
  
             ! Lloyd mirror pattern for semi-coherent option
             IF ( Beam%RunType( 1:1 ) == 'S' ) &
                Amp0 = Amp0 * SQRT( 2.0 ) * ABS( SIN( afreq / c * xs( 2 ) &
                       * SIN( Angles%alpha( ialpha ) ) ) )
             !!IESCO22: end BEAM stuff !!
  
#ifdef IHOP_WRITE_OUT
             ! report progress in PRTFile (skipping some angles)
             IF ( MOD( ialpha - 1, max( Angles%Nalpha / 50, 1 ) ) == 0 ) THEN
                WRITE(msgBuf,'(A,I7,F10.2)') 'Tracing ray ', &
                       ialpha, SrcDeclAngle
                ! In adjoint mode we do not write output besides on the first run
                IF (IHOP_dumpfreq.GE.0) &
                    CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
                FLUSH( PRTFile )
             END IF
#endif /* IHOP_WRITE_OUT */
             
             ! Trace a ray, update ray2D structure
             CALL TraceRay2D( xs, Angles%alpha( ialpha ), Amp0, myThid )   
  
             ! Write the ray trajectory to RAYFile
             IF ( Beam%RunType(1:1) == 'R') THEN   
                CALL WriteRay2D( SrcDeclAngle, Beam%Nsteps )
                IF (writeDelay) CALL WriteDel2D( SrcDeclAngle, Beam%Nsteps )
             ELSE ! Compute the contribution to the field
                
                SELECT CASE ( Beam%Type( 1:1 ) )
                CASE ( 'g' )
                   CALL InfluenceGeoHatRayCen(    U, Angles%alpha( ialpha ), &
                                                  Angles%Dalpha, myThid )
                CASE ( 'B' )
                   CALL InfluenceGeoGaussianCart( U, Angles%alpha( ialpha ), &
                                                  Angles%Dalpha, myThid )
                CASE ( 'G','^' )
                   CALL InfluenceGeoHatCart(  U, Angles%alpha( ialpha ), &
                                              Angles%Dalpha, myThid )
                CASE DEFAULT !IEsco22: thesis is in default behavior
                   CALL InfluenceGeoHatCart(  U, Angles%alpha( ialpha ), &
                                              Angles%Dalpha, myThid )
                END SELECT
             END IF
  
          END IF
       END DO DeclinationAngle
  
       ! write results to disk
  
       SELECT CASE ( Beam%RunType( 1:1 ) )
       CASE ( 'C', 'S', 'I' )   ! TL calculation
          CALL ScalePressure( Angles%Dalpha, ray2D( 1 )%c, Pos%Rr, U, &
                              NRz_per_range, Pos%NRr, Beam%RunType, IHOP_freq )
          IRec = 10 + NRz_per_range * ( is - 1 )
          RcvrDepth: DO Irz1 = 1, NRz_per_range
             IRec = IRec + 1
             WRITE( SHDFile, REC = IRec ) U( Irz1, 1:Pos%NRr )
          END DO RcvrDepth
       CASE ( 'A', 'e' )             ! arrivals calculation, ascii
          CALL WriteArrivalsASCII(  Pos%Rr, NRz_per_range, Pos%NRr, &
                                    Beam%RunType( 4:4 ) )
       CASE ( 'a' )             ! arrivals calculation, binary
          CALL WriteArrivalsBinary( Pos%Rr, NRz_per_range, Pos%NRr, &
                                    Beam%RunType( 4:4 ) )
       END SELECT
  
    END DO SourceDepth
  
  RETURN
  END !SUBROUTINE BellhopCore
  
  ! **********************************************************************!
  
  SUBROUTINE TraceRay2D( xs, alpha, Amp0, myThid )
  
    ! Traces the beam corresponding to a particular take-off angle, alpha [rad]
  
    USE step,     only: Step2D
  
  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN ) :: xs( 2 )     ! coordinate of source
    REAL (KIND=_RL90), INTENT( IN ) :: alpha, Amp0 ! angle in rad, beam amp
    INTEGER           :: is, is1                   ! indices for ray step
    REAL (KIND=_RL90) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
    REAL (KIND=_RL90) :: dEndTop( 2 ), dEndBot( 2 ), TopnInt( 2 ), BotnInt( 2 ), &
                         ToptInt( 2 ), BottInt( 2 ), rayt(2), raytOld(2)
    ! Distances from ray beginning, end to top and bottom
    REAL (KIND=_RL90) :: DistBegTop, DistEndTop, DistBegBot, DistEndBot 
    REAL (KIND=_RL90) :: sss, declAlpha, declAlphaOld
    LOGICAL           :: RayTurn = .FALSE., continue_steps
  
    ! Initial conditions (IC)
    iSmallStepCtr = 0
    CALL EvaluateSSP( xs, c, cimag, gradc, crr, crz, czz, rho, &
                      'TAB', myThid )
    ray2D( 1 )%c         = c              ! sound speed at source [m/s]
    ray2D( 1 )%x         = xs             ! range and depth of source
    ray2D( 1 )%t         = [ COS( alpha ), SIN( alpha ) ] / c ! unit tangent / c
    ray2D( 1 )%p         = [ 1.0, 0.0 ]   ! IESCO22: slowness vector
    ! second component of qv is not supported in geometric beam tracing
    ! set I.C. to 0 in hopes of saving run time
    IF ( Beam%RunType( 2:2 ) == 'G' .or. Beam%RunType( 2:2 ) == 'B')  THEN 
        ray2D( 1 )%q = [ 0.0, 0.0 ]   ! IESCO22: geometric beam in Cartesian
    ELSE
        ray2D( 1 )%q = [ 0.0, 1.0 ]   ! IESCO22: ray centered coords 
    END IF
    ray2D( 1 )%tau       = 0.0
    ray2D( 1 )%Amp       = Amp0
    ray2D( 1 )%Phase     = 0.0
    ray2D( 1 )%NumTopBnc = 0
    ray2D( 1 )%NumBotBnc = 0
    ray2D( 1 )%NumTurnPt = 0
  
    ! IESCO22: update IsegTop, rTopSeg and IsegBot, rBotSeg in bdrymod.f90
    CALL GetTopSeg( xs(1), myThid )   ! find alimetry   segment above the source
    CALL GetBotSeg( xs(1), myThid )   ! find bathymetry segment below the source
  
    ! IESCO22: 'L' is long format. See BeadBTY s/r in bdrymod.f90. Default is to
    ! calculate cp, cs, and rho instead of reading them in
    IF ( atiType( 2 : 2 ) == 'L' ) THEN
       ! grab the geoacoustic info for the new segment
       Bdry%Top%HS%cp  = Top( IsegTop )%HS%cp   
       Bdry%Top%HS%cs  = Top( IsegTop )%HS%cs
       Bdry%Top%HS%rho = Top( IsegTop )%HS%rho
    END IF
    IF ( btyType( 2 : 2 ) == 'L' ) THEN
       Bdry%Bot%HS%cp  = Bot( IsegBot )%HS%cp
       Bdry%Bot%HS%cs  = Bot( IsegBot )%HS%cs
       Bdry%Bot%HS%rho = Bot( IsegBot )%HS%rho
    END IF
  
    CALL Distances2D( ray2D( 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, &
                                    dEndTop,          dEndBot, &
                                    Top( IsegTop )%n, Bot( IsegBot )%n, &
                                    DistBegTop,       DistBegBot )
  
    IF ( DistBegTop <= 0 .OR. DistBegBot <= 0 ) THEN
       Beam%Nsteps = 1
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') &
          'WARNING: TraceRay2D: The source is outside the domain boundaries'
       ! In adjoint mode we do not write output besides on the first run
       IF (IHOP_dumpfreq.GE.0) &
            CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
#endif /* IHOP_WRITE_OUT */
       RETURN       ! source must be within the domain
    END IF
  
    ! Trace the beam (Reflect2D increments the step index, is)
    is = 0
    continue_steps = .true.
    Stepping: DO istep = 1, MaxN - 1
       IF ( continue_steps ) THEN
         is  = is + 1 ! old step
         is1 = is + 1 ! new step forward
  
         CALL Step2D( ray2D( is ), ray2D( is1 ),  &
              Top( IsegTop )%x, Top( IsegTop )%n, &
              Bot( IsegBot )%x, Bot( IsegBot )%n, myThid )
  
         ! IESCO22: turning point check
         IF ( is > 1 ) THEN
            rayt    = ray2D(is1)%x - ray2D(is)%x
            raytOld = ray2D(is)%x  - ray2D(is-1)%x
            declAlpha    = ATAN2( rayt(2), rayt(1) ) 
            declAlphaOld = ATAN2( raytOld(2), raytOld(1) ) 
            RayTurn = ( declAlpha <= 0.0d0 .AND. declAlphaOld > 0.0d0 .OR. &
                        declAlpha >= 0.0d0 .AND. declAlphaOld < 0.0d0 )
            IF ( RayTurn) THEN
               ray2D( is1 )%NumTurnPt = ray2D( is )%NumTurnPt + 1
            END IF
         END IF
  
         ! New altimetry segment?
         IF ( ray2D( is1 )%x( 1 ) < rTopSeg( 1 ) .OR. &
              ray2D( is1 )%x( 1 ) > rTopSeg( 2 ) ) THEN
            CALL GetTopSeg( ray2D( is1 )%x( 1 ), myThid )
            IF ( atiType( 2 : 2 ) == 'L' ) THEN
               ! ATIFile geoacoustic info from new segment, cp
               Bdry%Top%HS%cp  = Top( IsegTop )%HS%cp   
               Bdry%Top%HS%cs  = Top( IsegTop )%HS%cs
               Bdry%Top%HS%rho = Top( IsegTop )%HS%rho
            END IF
         END IF
  
         ! New bathymetry segment?
         IF ( ray2D( is1 )%x( 1 ) < rBotSeg( 1 ) .OR. &
              ray2D( is1 )%x( 1 ) > rBotSeg( 2 ) ) THEN
            CALL GetBotSeg( ray2D( is1 )%x( 1 ), myThid )
            IF ( btyType( 2 : 2 ) == 'L' ) THEN
               ! BTYFile geoacoustic info from new segment, cp
               Bdry%Bot%HS%cp  = Bot( IsegBot )%HS%cp   
               Bdry%Bot%HS%cs  = Bot( IsegBot )%HS%cs
               Bdry%Bot%HS%rho = Bot( IsegBot )%HS%rho
            END IF
         END IF
  
         ! *** Reflections ***
         ! Tests ray at step is IS inside, and ray at step is+1 IS outside
         ! DistBeg is the distance at step is,   which is saved
         ! DistEnd is the distance at step is+1, which needs to be calculated
  
         CALL Distances2D( ray2D( is1 )%x,  &
             Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
             Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )
  
         ! IESCO22: Did new ray point cross top boundary? Then reflect
         IF ( DistBegTop > 0.0d0 .AND. DistEndTop <= 0.0d0 ) THEN 
  
            IF ( atiType == 'C' ) THEN ! curvilinear interpolation
               ! proportional distance along segment
               sss     = DOT_PRODUCT( dEndTop, Top( IsegTop )%t ) &
                         / Top( IsegTop )%Len
               ToptInt = ( 1-sss ) * Top( IsegTop   )%Nodet &
                         + sss     * Top( 1+IsegTop )%Nodet
               TopnInt = ( 1-sss ) * Top( IsegTop   )%Noden &
                         + sss     * Top( 1+IsegTop )%Noden
            ELSE
               TopnInt = Top( IsegTop )%n   ! normal is constant in a segment
               ToptInt = Top( IsegTop )%t
            END IF
  
            CALL Reflect2D( is, Bdry%Top%HS,    'TOP',  ToptInt,    TopnInt, &
                                Top( IsegTop )%kappa,   RTop,       NTopPTS, &
                                myThid ) 
  
            CALL Distances2D( ray2D( is+1 )%x, &
                Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
                Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )
  
         ! IESCO22: Did ray cross bottom boundary? Then reflect
         ELSE IF ( DistBegBot > 0.0d0 .AND. DistEndBot <= 0.0d0 ) THEN
  
            IF ( btyType == 'C' ) THEN ! curvilinear interpolation
               ! proportional distance along segment
               sss     = DOT_PRODUCT( dEndBot, Bot( IsegBot )%t ) &
                         / Bot( IsegBot )%Len
               BotnInt = ( 1-sss ) * Bot( IsegBot   )%Noden &
                         + sss     * Bot( 1+IsegBot )%Noden
               BottInt = ( 1-sss ) * Bot( IsegBot   )%Nodet &
                         + sss     * Bot( 1+IsegBot )%Nodet
            ELSE
               BotnInt = Bot( IsegBot )%n   ! normal is constant in a segment
               BottInt = Bot( IsegBot )%t
            END IF
  
            CALL Reflect2D( is, Bdry%Bot%HS,    'BOT',  BottInt,    BotnInt, &
                                Bot( IsegBot )%kappa,   RBot,       NBotPTS, &
                                myThid ) 
  
            CALL Distances2D( ray2D( is+1 )%x, &
                Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
                Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )
         END IF
  
         ! Has the ray left the box, lost its energy, escaped the boundaries, 
         ! or exceeded storage limit?
         ! IESCO22: Rewriting for debugging with gcov
         WRITE(msgBuf,'(A)') ' '
         IF ( ray2D( is+1 )%x( 1 ) > Beam%Box%r ) THEN
            WRITE(msgBuf,'(A)') 'TraceRay2D: ray left Box%r'
         ELSE IF ( ray2D( is+1 )%x( 1 ) < 0 ) THEN
            WRITE(msgBuf,'(A)') 'TraceRay2D: ray left Box r=0'
         ELSE IF ( ray2D( is+1 )%x( 2 ) > Beam%Box%z ) THEN 
            WRITE(msgBuf,'(A)') 'TraceRay2D: ray left Box%z'
         ELSE IF ( ABS( ray2D( is+1 )%Amp ) < 0.005 ) THEN
            WRITE(msgBuf,'(A)') 'TraceRay2D: ray lost energy'
         ELSE IF ( DistBegTop < 0.0 .AND. DistEndTop < 0.0 ) THEN 
            WRITE(msgBuf,'(A)') 'TraceRay2D: ray escaped top bound'
         ELSE IF ( DistBegBot < 0.0 .AND. DistEndBot < 0.0 ) THEN
            WRITE(msgBuf,'(A)') 'TraceRay2D: ray escaped bot bound'
         ELSE IF ( is >= MaxN - 3 ) THEN
            WRITE(msgBuf,'(2A)') 'WARNING: TraceRay2D: Check storage ',&
                                 'for ray trajectory'
         END IF
     
#ifdef IHOP_WRITE_OUT
         IF ( ( ray2D( is+1 )%x( 1 ) > Beam%Box%r ) .OR. & 
              ( ray2D( is+1 )%x( 1 ) < 0 ) .OR. &
              ( ray2D( is+1 )%x( 2 ) > Beam%Box%z ) .OR. &
              ( ABS( ray2D( is+1 )%Amp ) < 0.005 ) .OR. &
              ( DistBegTop < 0.0 .AND. DistEndTop < 0.0 ) .OR. &
              ( DistBegBot < 0.0 .AND. DistEndBot < 0.0 ) .OR. &
              ( is >= MaxN - 3 ) ) THEN
              IF ( IHOP_dumpfreq .GE. 0) &
                  CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
         ENDIF
#endif /* IHOP_WRITE_OUT */
         IF (INDEX(msgBuf, 'TraceRay2D').eq.1) THEN
             Beam%Nsteps = is+1
             continue_steps = .false.
         ELSE IF (INDEX(msgBuf, 'WARNING: TraceRay2D').eq.1) THEN
             Beam%Nsteps = is
             continue_steps = .false.
         END IF 

         DistBegTop = DistEndTop
         DistBegBot = DistEndBot
      END IF ! continue_steps 
    END DO Stepping
  
  RETURN
  END !SUBROUTINE TraceRay2D
  
  ! **********************************************************************!
  
  SUBROUTINE Distances2D( rayx, Topx, Botx, dTop, dBot, Topn, Botn, DistTop, &
                          DistBot )
  
    ! Calculates the distances to the boundaries
    ! Formula differs from JKPS because code applies outward pointing normals
  
    REAL (KIND=_RL90), INTENT( IN  ) :: rayx(2)          ! ray coordinate
    REAL (KIND=_RL90), INTENT( IN  ) :: Topx(2), Botx(2) ! top, bottom coordinate
    REAL (KIND=_RL90), INTENT( IN  ) :: Topn(2), Botn(2) ! top, bottom normal vector (outward)
    REAL (KIND=_RL90), INTENT( OUT ) :: dTop(2), dBot(2) ! vector pointing from top, bottom bdry to ray
    REAL (KIND=_RL90), INTENT( OUT ) :: DistTop, DistBot ! distance (normal to bdry) from the ray to top, bottom boundary
  
    dTop    = rayx - Topx  ! vector pointing from top    to ray
    dBot    = rayx - Botx  ! vector pointing from bottom to ray
    DistTop = -DOT_PRODUCT( Topn, dTop )
    DistBot = -DOT_PRODUCT( Botn, dBot )
  
  RETURN
  END !SUBROUTINE Distances2D
  
  ! **********************************************************************!
  
  SUBROUTINE Reflect2D( is, HS, BotTop, tBdry, nBdry, kappa, RefC, Npts, myThid )
  
  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    INTEGER,              INTENT( IN ) :: Npts ! unsued if there are no refcoef files
    REAL (KIND=_RL90),    INTENT( IN ) :: tBdry(2), nBdry(2)  ! Tangent and normal to the boundary
    REAL (KIND=_RL90),    INTENT( IN ) :: kappa ! Boundary curvature, for curvilinear grids
    CHARACTER (LEN=3),    INTENT( IN ) :: BotTop       ! bottom or top reflection
    TYPE( HSInfo ),       INTENT( IN ) :: HS           ! half-space properties
    TYPE(ReflectionCoef), INTENT( IN ) :: RefC( NPts ) ! reflection coefficient
    INTEGER,              INTENT( INOUT ) :: is
    INTEGER              :: is1
    REAL (KIND=_RL90)    :: c, cimag, gradc( 2 ), crr, crz, czz, &
                            rho ! derivatives of sound speed
    REAL (KIND=_RL90)    :: RM, RN, Tg, Th, rayt( 2 ), rayn( 2 ), &
                            rayt_tilde( 2 ), rayn_tilde( 2 ), cnjump, &
                            csjump  ! for curvature change
    REAL (KIND=_RL90)    :: ck, co, si, cco, ssi, pdelta, rddelta, sddelta, &
                            theta_bot ! for beam shift
    COMPLEX (KIND=_RL90) :: kx, kz, kzP, kzS, kzP2, kzS2, mu, f, g, y2, y4, &
                            Refl   ! for tabulated reflection coef.
    COMPLEX (KIND=_RL90) :: ch, a, b, d, sb, delta, ddelta ! for beam shift
    TYPE(ReflectionCoef) :: RInt
  
    is  = is + 1 ! old step
    is1 = is + 1 ! new step reflected (same x, updated basis vectors)
  
    Tg = DOT_PRODUCT( ray2D( is )%t, tBdry )  ! ray tan projected along boundary
    Th = DOT_PRODUCT( ray2D( is )%t, nBdry )  ! ray tan projected normal boundary
  
    ray2D( is1 )%NumTopBnc = ray2D( is )%NumTopBnc
    ray2D( is1 )%NumBotBnc = ray2D( is )%NumBotBnc
    ray2D( is1 )%x         = ray2D( is )%x
    ray2D( is1 )%t         = ray2D( is )%t - 2.0 * Th * nBdry ! change ray direction
  
    ! Calculate change in curvature, kappa
    ! Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).
  
    ! Get c
    CALL EvaluateSSP( ray2D( is )%x, c, cimag, gradc, crr, crz, czz, rho, & 
                      'TAB', myThid  )
  
    ! unmodified unit ray tangent and normal
    rayt = c * ray2D( is )%t                              ! unit tangent to ray
    rayn = [ -rayt( 2 ), rayt( 1 ) ]                      ! unit normal  to ray
  
    ! reflected unit ray tangent and normal
    rayt_tilde = c * ray2D( is1 )%t                       ! unit tangent to ray
    rayn_tilde = -[ -rayt_tilde( 2 ), rayt_tilde( 1 ) ]   ! unit normal  to ray
  
    ! get the jumps (this could be simplified, e.g. jump in rayt is 
    ! roughly 2 * Th * nbdry
    cnjump = -DOT_PRODUCT( gradc, rayn_tilde - rayn  )
    csjump = -DOT_PRODUCT( gradc, rayt_tilde - rayt )
    RN = 2 * kappa / c ** 2 / Th    ! boundary curvature correction
  
    IF ( BotTop == 'TOP' ) THEN
       ! cnjump changes sign because the (t,n) system of the top boundary has a 
       ! different sense to the bottom boundary
       cnjump = -cnjump
       RN     = -RN
    END IF
  
    RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
    RN = RN + RM * ( 2 * cnjump - RM * csjump ) / c ** 2
  
    SELECT CASE ( Beam%Type( 3 : 3 ) )
    CASE ( 'D' )
       RN = 2.0 * RN
    CASE ( 'Z' )
       RN = 0.0
    END SELECT
    
    ray2D( is1 )%c   = c
    ray2D( is1 )%tau = ray2D( is )%tau
    ray2D( is1 )%p   = ray2D( is )%p + ray2D( is )%q * RN
    ray2D( is1 )%q   = ray2D( is )%q
  
    ! account for phase change
  
    SELECT CASE ( HS%BC )
    CASE ( 'R' )                 ! rigid
       ray2D( is1 )%Amp   = ray2D( is )%Amp
       ray2D( is1 )%Phase = ray2D( is )%Phase
    CASE ( 'V' )                 ! vacuum
       ray2D( is1 )%Amp   = ray2D( is )%Amp
       ray2D( is1 )%Phase = ray2D( is )%Phase + PI
    CASE ( 'F' )                 ! file
       RInt%theta = rad2deg * ABS( ATAN2( Th, Tg ) )           ! angle of incidence (relative to normal to bathymetry)
       IF ( RInt%theta > 90 ) RInt%theta = 180. - RInt%theta  ! reflection coefficient is symmetric about 90 degrees
       CALL InterpolateReflectionCoefficient( RInt, RefC, Npts )
       ray2D( is1 )%Amp   = ray2D( is )%Amp * RInt%R
       ray2D( is1 )%Phase = ray2D( is )%Phase + RInt%phi
    CASE ( 'A', 'G' )     ! half-space
       kx = afreq * Tg    ! wavenumber in direction parallel      to bathymetry
       kz = afreq * Th    ! wavenumber in direction perpendicular to bathymetry
  
       ! notation below is a bit mis-leading
       ! kzS, kzP is really what I called gamma in other codes, and differs by a 
       ! factor of +/- i
       IF ( REAL( HS%cS ) > 0.0 ) THEN
          kzS2 = kx**2 - ( afreq / HS%cS )**2
          kzP2 = kx**2 - ( afreq / HS%cP )**2
          kzS  = SQRT( kzS2 )
          kzP  = SQRT( kzP2 )
          mu   = HS%rho * HS%cS**2
  
          y2 = ( ( kzS2 + kx**2 )**2 - 4.0D0 * kzS * kzP * kx**2 ) * mu
          y4 = kzP * ( kx**2 - kzS2 )
  
          f = afreq**2 * y4
          g = y2
       ELSE
          kzP = SQRT( kx**2 - ( afreq / HS%cP )**2 )
  
          ! Intel and GFortran compilers return different branches of the SQRT 
          ! for negative reals
          IF ( REAL( kzP ) == 0.0D0 .AND. AIMAG( kzP ) < 0.0D0 ) kzP = -kzP
          f   = kzP
          g   = HS%rho
       ENDIF
  
       ! complex reflection coef.
       Refl =  - ( rho*f - i * kz*g ) / ( rho*f + i*kz*g )   
  
       IF ( ABS( Refl ) < 1.0E-5 ) THEN   ! kill a ray that has lost its energy in reflection
          ray2D( is1 )%Amp   = 0.0
          ray2D( is1 )%Phase = ray2D( is )%Phase
       ELSE
          ray2D( is1 )%Amp   = ABS( Refl ) * ray2D(  is )%Amp
          ray2D( is1 )%Phase = ray2D( is )%Phase + &
                               ATAN2( AIMAG( Refl ), REAL( Refl ) )
  
          if ( Beam%Type( 4:4 ) == 'S' ) then   ! beam displacement & width change (Seongil's version)
             ch = ray2D( is )%c / conjg( HS%cP )
             co = ray2D( is )%t( 1 ) * ray2D( is )%c
             si = ray2D( is )%t( 2 ) * ray2D( is )%c
             ck = afreq / ray2D( is )%c
  
             a   = 2 * HS%rho * ( 1 - ch * ch )
             b   = co * co - ch * ch
             d   = HS%rho * HS%rho * si * si + b
             sb  = sqrt( b )
             cco = co * co
             ssi = si * si
  
             IF ( si /= 0.0 ) THEN
                delta = a * co / si / ( ck * sb * d )   ! Do we need an abs() on this???
             ELSE
                delta = 0.0
             END IF
  
             pdelta  = real( delta ) / ( ray2D( is )%c / co)
             ddelta  = -a / ( ck*sb*d ) - a*cco / ssi / (ck*sb*d) &
                       + a*cco / (ck*b*sb*d) &
                       -a*co / si / (ck*sb*d*d) &
                       * (2* HS%rho * HS%rho *si*co-2*co*si)
             rddelta = -real( ddelta )
             sddelta = rddelta / abs( rddelta )        
  
             ! next 3 lines have an update by Diana McCammon to allow a sloping 
             ! bottom . I think the formulas are good, but this won't be reliable 
             ! because it doesn't have the logic that tracks crossing into new 
             ! segments after the ray displacement.
             
             theta_bot = datan( tBdry( 2 ) / tBdry( 1 ))  ! bottom angle
             ray2D( is1 )%x( 1 ) = ray2D( is1 )%x( 1 ) + real( delta ) &
                                   * dcos( theta_bot )   ! range displacement
             ray2D( is1 )%x( 2 ) = ray2D( is1 )%x( 2 ) + real( delta ) &
                                   * dsin( theta_bot )   ! depth displacement
             ray2D( is1 )%tau    = ray2D( is1 )%tau + pdelta  ! phase change
             ray2D( is1 )%q      = ray2D( is1 )%q + sddelta * rddelta * si * c &
                                   * ray2D( is )%p   ! beam-width change
          endif
  
       ENDIF
  
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(2A)') 'HS%BC = ', HS%BC
       ! In adjoint mode we do not write output besides on the first run
       IF (IHOP_dumpfreq.GE.0) &
        CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
       WRITE(msgBuf,'(A)') 'BELLHOP Reflect2D: Unknown boundary condition type'
       CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
       STOP 'ABNORMAL END: S/R Reflect2D'
    END SELECT
  
    ! Update top/bottom bounce counter
    IF (BotTop == 'TOP') THEN
       ray2D( is+1 )%NumTopBnc = ray2D( is )%NumTopBnc + 1
    ELSE IF ( BotTop == 'BOT' ) THEN
       ray2D( is+1 )%NumBotBnc = ray2D( is )%NumBotBnc + 1
    ELSE
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(2A)') 'BELLHOP Reflect2D: ', & 
                            'no reflection bounce, but in relfect2d somehow'
       CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
       STOP 'ABNORMAL END: S/R Reflect2D'
    END IF
  
  RETURN
  END !SUBROUTINE Reflect2D
  
END MODULE BELLHOP
