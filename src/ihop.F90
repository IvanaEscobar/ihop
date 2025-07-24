#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: IHOP
MODULE IHOP
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   Written as a module to be used in ihop_*.F parts of the MITgcm package
!   IHOP Beam tracing for ocean acoustics
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Second version BELLHOP Copyright (C) 2009 Michael B. Porter
!
! First version (1983) originally developed with Homer Bucker, Naval Ocean
! Systems Center

! !USES:
  USE ihop_mod, only: oneCMPLX, PRTFile, SHDFile, ARRFile, RAYFile, DELFile
  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#ifdef ALLOW_USE_MPI
# include "EESUPPORT.h"
#endif
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef ALLOW_CTRL
# include "CTRL_FIELDS.h"
#endif

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC ihop_main
!=======================================================================

! == Module Variables == None

! == External Functions ==
  INTEGER  ILNBLNK
  EXTERNAL ILNBLNK
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R IHOP_MAIN
! S/R IHOPCore
! S/R TraceRay2D
! S/R Distances2D
! S/R Reflect2D
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: IHOP_MAIN
! !INTERFACE:
  SUBROUTINE IHOP_MAIN ( myTime, myIter, myThid )
! !DESCRIPTION:
!   Main routine for IHOP ray tracing

! !USES:
  USE ihop_init_diag, only: initPRTFile, openOutputFiles, resetMemory
  USE bdry_mod,       only: Bdry, writeBdry
  USE ssp_mod,        only: setSSP
  USE refCoef,        only: writeRefCoef 
  USE beampat,        only: writePat
  USE ihop_mod,       only: Beam

! !INPUT PARAMETERS:
! myTime   :: time in seconds
! myIter   :: iteration number
! myThid   :: my thread ID
  _RL, INTENT( IN )       :: myTime
  INTEGER, INTENT( IN )   :: myIter
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! Tstart, Tstop :: Timing variables
! mpiRC :: MPI return code
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  REAL :: Tstart, Tstop
#ifdef ALLOW_USE_MPI
  INTEGER :: mpiRC
  mpiRC  = 0
#endif

  Tstart = 0.
  Tstop  = 0.

  ! Reset memory and default values: REQUIRED
  CALL resetMemory()

  ! Only do IO and computation on a single process!
  IF ( .NOT.usingMPI .AND. IHOP_dumpfreq.GE.0 ) THEN
    myProcId=0

    ! open PRTFile
    CALL initPRTFile( myTime, myIter, myThid )

    ! write Top/Bot to PRTFile: REQUIRED
    CALL writeBdry( myThid )

    ! write refCoef: OPTIONAL
    CALL writeRefCoef( myThid ) 

    ! Source Beam Pattern: OPTIONAL, default is omni source pattern
    CALL writePat( myThid )

    ! Open output files
    CALL OpenOutputFiles( IHOP_fileroot, myTime, myIter, myThid )

#ifdef ALLOW_USE_MPI
  ELSE ! using MPI
    CALL MPI_COMM_RANK( MPI_COMM_MODEL, mpiMyId, mpiRC )
    myProcId = mpiMyId

    ! Hard coded write on single proc
    IF ( myProcId.EQ.0 .AND. IHOP_dumpfreq.GE.0 ) THEN
      ! open PRTFile
      CALL initPRTFile( myTime, myIter, myThid )

      ! write Top/Bot to PRTFile: REQUIRED
      CALL writeBdry( myThid )

      ! write refCoef: OPTIONAL
      CALL writeRefCoef( myThid ) 

      ! Source Beam Pattern: OPTIONAL, default is omni source pattern
      CALL writePat( myThid )

      ! Open output files:
      CALL OpenOutputFiles( IHOP_fileroot, myTime, myIter, myThid )

    ENDIF ! IF ( myProcId.EQ.0 ) 

#endif /* ALLOW_USE_MPI */
  ENDIF ! IF ( .NOT.usingMPI )

  ! set SSP%cmat from gcm SSP: REQUIRED
  CALL setSSP( myThid )

  ! Run IHOP solver on a single processor
  IF ( myProcId.EQ.0 ) THEN
    CALL CPU_TIME( Tstart )
    CALL IHOPCore( myThid )
    CALL CPU_TIME( Tstop  )
  ENDIF ! IF ( myProcId.EQ.0 )
#ifdef ALLOW_USE_MPI
  IF ( usingMPI ) CALL MPI_BARRIER(MPI_COMM_WORLD, mpiRC)
#endif

#ifdef IHOP_WRITE_OUT
  IF ( myProcId.EQ.0 .AND. IHOP_dumpfreq.GE.0 ) THEN
    WRITE(msgBuf, '(A)' )
    CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
    WRITE(msgBuf, '(A,G15.3,A)' ) 'CPU Time = ', Tstop-Tstart, 's'
    CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
    CLOSE(PRTFile)

    ! Close all files
    SELECT CASE ( Beam%RunType( 1:1 ) )
    CASE ( 'C', 'S', 'I' )  ! TL calculation
      CLOSE( SHDFile )
    CASE ( 'A', 'a' )       ! arrivals calculation
      CLOSE( ARRFile )
    CASE ( 'R', 'E' )       ! ray and eigen ray trace
      CLOSE( RAYFile )
      IF ( writeDelay ) CLOSE( DELFile )
    CASE ( 'e' )
      CLOSE( RAYFile )
      CLOSE( ARRFile )
      IF ( writeDelay ) CLOSE( DELFile )
    CASE DEFAULT
      STOP "ABNORMAL END: S/R IHOP_MAIN"
    END SELECT

  ENDIF ! IF ( myProcId.EQ.0 )
#endif /* IHOP_WRITE_OUT */

  RETURN
  END !SUBROUTINE IHOP_MAIN

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: IHOPCore
! !INTERFACE:
  SUBROUTINE IHOPCore( myThid )
! !DESCRIPTION:
!   IHOP core solver.

! !USES:
  USE ssp_mod,   only: evalSSP, iSegr  !RG
  USE angle_mod, only: Angles, iAlpha
  USE srPos_mod, only: Pos
  USE arr_mod,   only: WriteArrivalsASCII, WriteArrivalsBinary, nArr, U
  USE writeRay,  only: WriteRayOutput
  USE influence, only: InfluenceGeoHatRayCen, InfluenceGeoGaussianCart, &
                       InfluenceGeoHatCart, ScalePressure
  USE beampat,   only: NSBPPts, SrcBmPat
  USE ihop_mod,  only: Beam, ray2D, rad2deg, SrcDeclAngle, afreq, &
                       nRz_per_range, RAYFile, DELFile, maxN
! IESCO25: FOR TAF -RG
!  USE influence,  only: ratio1, rB
!  USE bdry_mod, only: bdry
!  USE arr_mod,  only: Arr

! !INPUT PARAMETERS:
! myThid :: my thread ID
  INTEGER, INTENT( IN )   :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! IBPvec :: Index of beam pattern
! ibp    :: Index of beam pattern
! is     :: Index of source depth
! iBeamWindow2 :: Index of beam window
! Irz1   :: Index of receiver depth
! iRec   :: Index of receiver
! nAlphaOpt :: Number of optimal angles
! nSteps :: Number of steps in ray tracing
! Amp0   :: Initial amplitude of the ray
! DalphaOpt :: Optimal angle step
! xs     :: Source coordinates
! RadMax :: Maximum radius for the beam
! s     :: Interpolation factor
! c     :: Sound speed
! cimag :: Imaginary part of sound speed
! gradc :: Gradient of sound speed
! crr   :: Radial component of sound speed
! crz   :: Vertical component of sound speed
! czz   :: Axial component of sound speed
! rho   :: Density
! tmpDelay :: Temporary delay array for ray tracing
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER :: IBPvec(1), ibp, is, iBeamWindow2, Irz1, iRec, &
             nAlphaOpt, nSteps
  REAL(KIND=_RL90) :: Amp0, DalphaOpt, xs(2), RadMax, s, &
                      c, cimag, gradc(2), crr, crz, czz, rho
  REAL(KIND=_RL90) :: tmpDelay(maxN)
!EOP

!$TAF init IHOPCore2 = static, Pos%nSZ*Angles%nAlpha

  afreq = 2.0 * PI * IHOP_freq

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !         begin solve         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SourceDepth: DO is = 1, Pos%nSZ
    xs = [ zeroRL, Pos%SZ( is ) ]  ! assuming source @ r=0

    ! Reset U and nArr for each new source depth
    SELECT CASE ( Beam%RunType( 1:1 ) )
    CASE ( 'C','S','I' ) ! TL calculation, zero out pressure matrix
      U = 0.0
      nArr = 0
    CASE ( 'A','a','e' )   ! Arrivals calculation, zero out arrival matrix
      U = 0.0
      nArr = 0
    CASE DEFAULT ! Ray tracing only
      U = 0.0
      nArr = 0
    END SELECT

    CALL evalSSP(  xs, c, cimag, gradc, crr, crz, czz, rho, myThid  )

    !!IESCO22: BEAM stuff !!
    RadMax = 5 * c / IHOP_freq  ! 5 wavelength max radius IEsco22: unused
    IF ( Beam%RunType( 1:1 ).EQ.'C' ) THEN ! for Coherent TL Run
      ! Are there enough rays?
      DalphaOpt = SQRT( c / ( 6.0 * IHOP_freq * Pos%RR( Pos%nRR ) ) )
      nAlphaOpt = 2 + INT( ( Angles%arad( Angles%nAlpha ) &
                            - Angles%arad( 1 ) ) / DalphaOpt )
#ifdef IHOP_WRITE_OUT
      IF ( Angles%nAlpha .LT. nAlphaOpt ) THEN
        WRITE( msgBuf, '(A,/,A,I10.4)' ) 'WARNING: Too few beams', &
          'nAlpha should be at least = ', nAlphaOpt
        ! In adjoint mode we do not write output besides on the first run
        IF (IHOP_dumpfreq.GE.0) &
          CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
      ENDIF
#endif /* IHOP_WRITE_OUT */
    ENDIF ! IF ( Beam%RunType( 1:1 ).EQ.'C' )
    !!IESCO22: end BEAM stuff !!

    ! Trace beams: IESCO25 MPI distribute this loop!
    DeclinationAngle: DO iAlpha = 1, Angles%nAlpha

!$TAF store ray2d,arr,nArr,u = IHOPCore2

      ! take-off declination angle in degrees
      SrcDeclAngle = Angles%adeg( iAlpha )

      ! Single ray run? then don't visit code below
      IF ( Angles%iSingle_alpha.EQ.0 .OR. iAlpha.EQ.Angles%iSingle_alpha ) THEN
        !!IESCO22: BEAM stuff !!
        IBPvec = maxloc( SrcBmPat( :, 1 ), mask=SrcBmPat( :, 1 ) &
                < SrcDeclAngle )  ! index of ray angle in beam pattern
        IBP    = IBPvec( 1 )
        IBP    = MAX( IBP, 1 )         ! don't go before beginning of table
        IBP    = MIN( IBP, NSBPPts-1 ) ! don't go past end of table
        ! IEsco22: When a beam pattern isn't specified, IBP = 1

        ! linear interpolation to get amplitude
        s    = ( SrcDeclAngle  - SrcBmPat( IBP, 1 ) ) &
              / ( SrcBmPat( IBP + 1, 1 ) - SrcBmPat( IBP, 1 ) )
        Amp0 = ( 1 - s ) * SrcBmPat( IBP, 2 ) + s * SrcBmPat( IBP+1, 2 )
        ! IEsco22: When a beam pattern isn't specified, Amp0 = 0

!$TAF store amp0,beam%runtype,beam%nsteps = IHOPCore2
! IESCO24: Store derived type by data type: Bdry from bdry_mod
! Scalar components:
!$TAF store bdry%top%hs%cp,bdry%top%hs%cs,bdry%top%hs%rho = IHOPCore2
! Fixed arrays:
! Allocatable arrays:

        ! Lloyd mirror pattern for semi-coherent option
        IF ( Beam%RunType( 1:1 ).EQ.'S' ) &
!$TAF store amp0,beam%runtype = IHOPCore2
          Amp0 = Amp0 * SQRT( 2.0 ) * ABS( SIN( afreq / c * xs( 2 ) &
                  * SIN( Angles%arad( iAlpha ) ) ) )
        !!IESCO22: end BEAM stuff !!

#ifdef IHOP_WRITE_OUT
          ! report progress in PRTFile (skipping some angles)
          IF ( MOD( iAlpha-1, MAX( Angles%nAlpha/50, 1 ) ).EQ.0 ) THEN
            WRITE(msgBuf,'(A,I7,F10.2)') 'Tracing ray ', &
              iAlpha, SrcDeclAngle
            ! In adjoint mode we do not write output besides on the first run
            IF (IHOP_dumpfreq.GE.0) &
              CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
            FLUSH( PRTFile )
          ENDIF
#endif /* IHOP_WRITE_OUT */

          ! Trace a ray, update ray2D structure
          CALL TraceRay2D( xs, Angles%arad( iAlpha ), Amp0, myThid )

          ! Write the ray trajectory to RAYFile
          IF ( Beam%RunType(1:1).EQ.'R') THEN
            nSteps = Beam%nSteps
            CALL WriteRayOutput( RAYFile, nSteps, &
              ray2D(1:nSteps)%x(1),    ray2D(1:nSteps)%x(2), &
              ray2D(nSteps)%NumTopBnc, ray2D(nSteps)%NumBotBnc )
            IF (writeDelay) THEN
              tmpDelay = REAL(ray2D(1:nSteps)%tau)
              CALL WriteRayOutput( DELFile, nSteps, &
                tmpDelay,                ray2D(1:nSteps)%x(2), &
                ray2D(nSteps)%NumTopBnc, ray2D(nSteps)%NumBotBnc )
            ENDIF

          ELSE ! Compute the contribution to the field
            SELECT CASE ( Beam%Type( 1:1 ) )
            CASE ( 'g' )
              CALL InfluenceGeoHatRayCen( U, Angles%Dalpha, myThid )
            CASE ( 'B' )
              CALL InfluenceGeoGaussianCart( U, Angles%Dalpha, &
                myThid )
            CASE ( 'G','^' )
              CALL InfluenceGeoHatCart( U, Angles%Dalpha, myThid )
            CASE DEFAULT !IEsco22: thesis is in default behavior
              CALL InfluenceGeoHatCart( U, Angles%Dalpha, myThid )
            END SELECT

          ENDIF ! IF ( Beam%RunType(1:1).EQ.'R')

      ENDIF ! IF ( Angles%iSingle_alpha.EQ.0 .OR. iAlpha.EQ.Angles%iSingle_alpha )

    ENDDO DeclinationAngle

    ! Write results to disk
    SELECT CASE ( Beam%RunType( 1:1 ) )
    CASE ( 'C', 'S', 'I' )   ! TL calculation
      CALL ScalePressure( Angles%Dalpha, ray2D( 1 )%c, Pos%RR, U, &
                          nRz_per_range, Pos%nRR, Beam%RunType, &
                          IHOP_freq )
      iRec = 10 + nRz_per_range * ( is-1 )
      RcvrDepth: DO Irz1 = 1, nRz_per_range
        iRec = iRec + 1
        WRITE( SHDFile, REC=iRec ) U( Irz1, 1:Pos%nRR )
      ENDDO RcvrDepth

    CASE ( 'A', 'e' ) ! arrivals calculation, ascii
      CALL WriteArrivalsASCII(  Pos%RR, nRz_per_range, Pos%nRR, &
                                Beam%RunType( 4:4 ) )
    CASE ( 'a' ) ! arrivals calculation, binary
      CALL WriteArrivalsBinary( Pos%RR, nRz_per_range, Pos%nRR, &
                                Beam%RunType( 4:4 ) )
    CASE DEFAULT ! ray trace only, NO ARRFile
      !CALL WriteArrivalsASCII( Pos%RR, nRz_per_range, Pos%nRR, &
      !                         Beam%RunType( 4:4 ) )
    END SELECT

  ENDDO SourceDepth

  RETURN
  END !SUBROUTINE IHOPCore

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: TraceRay2D
! !INTERFACE:
  SUBROUTINE TraceRay2D( xs, alpha, Amp0, myThid )
! !DESCRIPTION:
!   Traces a ray in 2D geometry, starting at source xs with take-off angle alpha [rad]

! !USES:
  USE ihop_mod, only: maxN, istep
  USE bdry_mod, only: GetTopSeg, GetBotSeg, atiType, btyType, Bdry, &
                      IsegTop, IsegBot, rTopSeg, rBotSeg, Top, Bot
  USE refCoef,  only: RTop, RBot, NBotPts, NTopPts
  USE step,     only: Step2D
  USE ihop_mod, only: Beam, ray2D, iSmallStepCtr
! IESCO25: FOR TAF -RG
  USE ssp_mod,  only: evalSSP, iSegr

! !INPUT PARAMETERS:
! xs     :: Source coordinates [m, m]
! alpha  :: Take-off angle [rad]
! Amp0   :: Initial beam amplitude
! myThid :: my thread ID
  REAL (KIND=_RL90), INTENT( IN ) :: xs(2)
  REAL (KIND=_RL90), INTENT( IN ) :: alpha, Amp0
  INTEGER,           INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! is, is1 :: indices for ray step
! c, cimag :: sound speed and its imaginary part
! gradc :: gradient of sound speed
! crr, crz, czz :: radial, vertical, and axial components of sound speed
! rho :: density
! dEndTop, dEndBot :: Distances from ray beginning, end to top and bottom
! TopnInt, BotnInt :: Interpolated normals at top and bottom
! ToptInt, BottInt :: Interpolated tangents at top and bottom
! rayt, raytOld :: ray tangents at current and previous step
! DistBegTop, DistEndTop :: Distances from ray beginning, end to top
! DistBegBot, DistEndBot :: Distances from ray beginning, end to bottom
! sss :: proportional distance along segment
! declAlpha, declAlphaOld :: Declination angles at current and previous step
! RayTurn :: Logical flag for ray turning point
! endRay :: Logical flag for ray end
! continue_steps :: Logical flag to continue ray steps
! reflect :: Logical flag for ray reflection
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  INTEGER           :: is, is1
  REAL (KIND=_RL90) :: c, cimag, gradc(2), crr, crz, czz, rho
  REAL (KIND=_RL90) :: dEndTop(2), dEndBot(2), TopnInt(2), BotnInt(2), &
                       ToptInt(2), BottInt(2), rayt(2), raytOld(2)
  REAL (KIND=_RL90) :: DistBegTop, DistEndTop, DistBegBot, DistEndBot
  REAL (KIND=_RL90) :: sss, declAlpha, declAlphaOld
  LOGICAL           :: RayTurn=.FALSE., endRay=.FALSE.
  LOGICAL           :: continue_steps, reflect
!EOP

!$TAF init TraceRay2D = static, maxN-1

  ! Initial conditions (IC)
  iSmallStepCtr = 0
  CALL evalSSP( xs, c, cimag, gradc, crr, crz, czz, rho, myThid )
  ray2D( 1 )%c         = c              ! sound speed at source [m/s]
  ray2D( 1 )%x         = xs             ! range and depth of source
  ray2D( 1 )%t         = [ COS( alpha ), SIN( alpha ) ] / c ! unit tangent / c
  ray2D( 1 )%p         = [ 1.0, 0.0 ]   ! IESCO22: slowness vector
  ! second component of qv is not supported in geometric beam tracing
  ! set I.C. to 0 in hopes of saving run time
!$TAF store beam%runtype = TraceRay2D
  IF ( Beam%RunType( 2:2 ).EQ.'G' .OR. Beam%RunType( 2:2 ).EQ.'B')  THEN
    ray2D( 1 )%q = [ 0.0, 0.0 ]   ! IESCO22: geometric beam in Cartesian
  ELSE
    ray2D( 1 )%q = [ 0.0, 1.0 ]   ! IESCO22: ray centered coords
  ENDIF

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
  IF ( atiType( 2:2 ).EQ.'L' ) THEN
    ! grab the geoacoustic info for the new segment
    Bdry%Top%HS%cp  = Top( IsegTop )%HS%cp
    Bdry%Top%HS%cs  = Top( IsegTop )%HS%cs
    Bdry%Top%HS%rho = Top( IsegTop )%HS%rho
  ENDIF

  IF ( btyType( 2:2 ).EQ.'L' ) THEN
    ! grab the geoacoustic info for the new segment
    Bdry%Bot%HS%cp  = Bot( IsegBot )%HS%cp
    Bdry%Bot%HS%cs  = Bot( IsegBot )%HS%cs
    Bdry%Bot%HS%rho = Bot( IsegBot )%HS%rho
  ENDIF

  CALL Distances2D( ray2D( 1 )%x, &
    Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
    Top( IsegTop )%n, Bot( IsegBot )%n, DistBegTop, DistBegBot )

! Source MUST be within the domain
  IF ( DistBegTop.LE.0 .OR. DistBegBot.LE.0 ) THEN
    Beam%Nsteps = 1
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') &
      'WARNING: TraceRay2D: The source is outside the domain boundaries'
    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.GE.0) &
      CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
#endif /* IHOP_WRITE_OUT */
  ELSE
    ! Trace the beam (Reflect2D increments the step index, is)
    is = 0
    continue_steps = .true.
    reflect=.false.

    Stepping: DO istep = 1, maxN-1
!$TAF store is,bdry,beam,continue_steps,distbegbot,distbegtop = TraceRay2D

      IF ( continue_steps ) THEN
        is  = is + 1 ! old step
        is1 = is + 1 ! new step forward

!$TAF store is,isegbot,isegtop,rbotseg,rtopseg = TraceRay2D
!$TAF store ray2d = TraceRay2D

        CALL Step2D( ray2D( is ), ray2D( is1 ), &
          Top( IsegTop )%x, Top( IsegTop )%n,   &
          Bot( IsegBot )%x, Bot( IsegBot )%n, myThid )

        ! IESCO22: turning point check
        IF ( is.GT.1 ) THEN
          rayt    = ray2D(is1)%x - ray2D(is)%x
          raytOld = ray2D(is )%x - ray2D(is-1)%x

          IF ( ALL(rayt.EQ.0.0) ) THEN
            declAlpha = 0.0
          ELSE
            declAlpha = ATAN2( rayt(2), rayt(1) )
          ENDIF

          IF ( ALL(raytOld.EQ.0.0) ) THEN
            declAlphaOld = 0.0
          ELSE
            declAlphaOld = ATAN2( raytOld(2), raytOld(1) )
          ENDIF

          RayTurn = ( declAlpha.LE.0.0d0 .AND. declAlphaOld.GT.0.0d0 .OR. &
                      declAlpha.GE.0.0d0 .AND. declAlphaOld.LT.0.0d0 )
          IF (RayTurn) ray2D( is1 )%NumTurnPt = ray2D( is )%NumTurnPt + 1

        ENDIF ! IF ( is.GT.1 )

        ! New altimetry segment?
        IF ( ray2D( is1 )%x( 1 ).LT.rTopSeg( 1 ) .OR. &
             ray2D( is1 )%x( 1 ).GT.rTopSeg( 2 ) ) THEN
          CALL GetTopSeg( ray2D( is1 )%x( 1 ), myThid )

          IF ( atiType( 2:2 ).EQ.'L' ) THEN
            ! ATIFile geoacoustic info from new segment, cp
            Bdry%Top%HS%cp  = Top( IsegTop )%HS%cp
            Bdry%Top%HS%cs  = Top( IsegTop )%HS%cs
            Bdry%Top%HS%rho = Top( IsegTop )%HS%rho
          ENDIF

        ENDIF

        ! New bathymetry segment?
        IF ( ray2D( is1 )%x( 1 ).LT.rBotSeg( 1 ) .OR. &
             ray2D( is1 )%x( 1 ).GT.rBotSeg( 2 ) ) THEN
          CALL GetBotSeg( ray2D( is1 )%x( 1 ), myThid )

          IF ( btyType( 2:2 ).EQ.'L' ) THEN
            ! BTYFile geoacoustic info from new segment, cp
            Bdry%Bot%HS%cp  = Bot( IsegBot )%HS%cp
            Bdry%Bot%HS%cs  = Bot( IsegBot )%HS%cs
            Bdry%Bot%HS%rho = Bot( IsegBot )%HS%rho
          ENDIF

        ENDIF

        ! *** Reflections ***
        ! Tests ray at step is IS inside, and ray at step is+1 IS outside
        ! DistBeg is the distance at step is,   which is saved
        ! DistEnd is the distance at step is+1, which needs to be calculated

        CALL Distances2D( ray2D( is1 )%x,  &
          Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
          Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

        ! IESCO22: Did new ray point cross top boundary? Then reflect
        IF ( DistBegTop.GT.0.0d0 .AND. DistEndTop.LE.0.0d0 ) THEN
          reflect=.true.

!$TAF store isegtop = TraceRay2D

          IF ( atiType.EQ.'C' ) THEN ! curvilinear interpolation
            ! proportional distance along segment
            sss     = DOT_PRODUCT( dEndTop, Top( IsegTop )%t ) &
                      / Top( IsegTop )%Len
            ToptInt = ( 1-sss ) * Top( IsegTop   )%Nodet &
                      + sss     * Top( 1+IsegTop )%Nodet
            TopnInt = ( 1-sss ) * Top( IsegTop   )%Noden &
                      + sss     * Top( 1+IsegTop )%Noden
          ELSE
            TopnInt = Top( IsegTop )%n   ! constant normal in a segment
            ToptInt = Top( IsegTop )%t
          ENDIF

!$TAF store is,isegtop = TraceRay2D
          CALL Reflect2D( is, Bdry%Top%HS, 'TOP', ToptInt, TopnInt, &
            Top( IsegTop )%kappa, RTop, NTopPTS, myThid )

!$TAF store is,isegbot = TraceRay2D
          CALL Distances2D( ray2D( is+1 )%x, &
            Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
            Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

        ! IESCO22: Did ray cross bottom boundary? Then reflect
        ELSEIF ( DistBegBot.GT.0.0d0 .AND. DistEndBot.LE.0.0d0 ) THEN
          reflect=.true.

!$TAF store isegbot = TraceRay2D
          IF ( btyType.EQ.'C' ) THEN ! curvilinear interpolation
            ! proportional distance along segment
            sss     = DOT_PRODUCT( dEndBot, Bot( IsegBot )%t ) &
                      / Bot( IsegBot )%Len
            BotnInt = ( 1-sss ) * Bot( IsegBot   )%Noden &
                      + sss     * Bot( 1+IsegBot )%Noden
            BottInt = ( 1-sss ) * Bot( IsegBot   )%Nodet &
                      + sss     * Bot( 1+IsegBot )%Nodet
          ELSE ! btyType not 'C'
            BotnInt = Bot( IsegBot )%n   ! normal is constant in a segment
            BottInt = Bot( IsegBot )%t
          ENDIF

!$TAF store is,isegbot = TraceRay2D
          CALL Reflect2D( is, Bdry%Bot%HS, 'BOT', BottInt, BotnInt, &
            Bot( IsegBot )%kappa, RBot, NBotPTS, myThid )

!$TAF store is,isegbot = TraceRay2D
          CALL Distances2D( ray2D( is+1 )%x, &
            Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
            Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

        ELSE
          ! Do not reflect
          reflect=.false.

        ENDIF ! IF ( DistBegTop .GT. 0.0d0 .AND. DistEndTop .LE. 0.0d0 )

        ! Has the ray left the box, lost its energy, escaped the boundaries,
        ! or exceeded storage limit?
        ! IESCO22: Rewriting for debugging with gcov
        WRITE(msgBuf,'(A)') ' '
        IF ( ray2D( is+1 )%x( 1 ).GT.Beam%Box%R ) THEN
          WRITE(msgBuf,'(A)') 'TraceRay2D: ray left Box%R'
          endRay=.TRUE.
        ELSEIF ( ray2D( is+1 )%x( 1 ).LT.0 ) THEN
          WRITE(msgBuf,'(A)') 'TraceRay2D: ray left Box r=0'
          endRay=.TRUE.
        ELSEIF ( ray2D( is+1 )%x( 2 ).GT.Beam%Box%Z ) THEN
          WRITE(msgBuf,'(A)') 'TraceRay2D: ray left Box%Z'
          endRay=.TRUE.
        ELSEIF ( ABS( ray2D( is+1 )%Amp ).LT.0.005 ) THEN
          WRITE(msgBuf,'(A)') 'TraceRay2D: ray lost energy'
          endRay=.TRUE.
        ELSEIF ( DistBegTop.LT.0.0 .AND. DistEndTop.LT.0.0 ) THEN
          WRITE(msgBuf,'(A)') 'TraceRay2D: ray escaped top bound'
          endRay=.TRUE.
        ELSEIF ( DistBegBot.LT.0.0 .AND. DistEndBot.LT.0.0 ) THEN
          WRITE(msgBuf,'(A)') 'TraceRay2D: ray escaped bot bound'
          endRay=.TRUE.
        ELSEIF ( is.GE.maxN-3 ) THEN
          WRITE(msgBuf,'(2A)') 'WARNING: TraceRay2D: Check storage ',&
                                'for ray trajectory'
          endRay=.TRUE.
        ELSE
          WRITE(msgBuf,'(A)')
          endRay=.FALSE.
        ENDIF

#ifdef IHOP_WRITE_OUT
        IF ( endRay ) THEN
          IF ( IHOP_dumpfreq.GE.0 ) &
            CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
        ENDIF
#endif /* IHOP_WRITE_OUT */

        IF (INDEX(msgBuf, 'TraceRay2D').EQ.1) THEN
          Beam%Nsteps = is+1
          continue_steps = .false.
        ELSEIF (INDEX(msgBuf, 'WARNING: TraceRay2D').EQ.1) THEN
          Beam%Nsteps = is
          continue_steps = .false.
        ELSE
          continue_steps = .true.
        ENDIF

        DistBegTop = DistEndTop
        DistBegBot = DistEndBot

      ENDIF ! continue_steps

    ENDDO Stepping

  ENDIF ! IF ( DistBegTop.LE.0 .OR. DistBegBot.LE.0 )

  RETURN
  END !SUBROUTINE TraceRay2D

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: Distances2D
! !INTERFACE:
  SUBROUTINE Distances2D( rayx, Topx, Botx, dTop, dBot, Topn, Botn, &
                          DistTop, DistBot )
! !DESCRIPTION:
!   Calculates the distances from a ray to the top and bottom boundaries

! !USES: None

! !INPUT PARAMETERS:
! rayx  :: Ray coordinate [m, m]
! Topx  :: Top boundary coordinate [m, m]
! Botx  :: Bottom boundary coordinate [m, m]
! Topn  :: Top boundary normal vector (outward) [m, m]
! Botn  :: Bottom boundary normal vector (outward) [m, m]
! DistTop, DistBot :: Distances from the ray to top and bottom boundaries
!                    (normal to boundary) [m]
  REAL (KIND=_RL90), INTENT( IN  ) :: rayx(2)
  REAL (KIND=_RL90), INTENT( IN  ) :: Topx(2), Botx(2)
  REAL (KIND=_RL90), INTENT( IN  ) :: Topn(2), Botn(2) 
  REAL (KIND=_RL90), INTENT( OUT ) :: dTop(2), dBot(2)
  REAL (KIND=_RL90), INTENT( OUT ) :: DistTop, DistBot
! !OUTPUT PARAMETERS: dTop, dBot, DistTop, DistBot

! !LOCAL VARIABLES: None
!EOP

  dTop    = rayx - Topx  ! vector pointing from top    to ray
  dBot    = rayx - Botx  ! vector pointing from bottom to ray
  DistTop = -DOT_PRODUCT( Topn, dTop )
  DistBot = -DOT_PRODUCT( Botn, dBot )

  RETURN
  END !SUBROUTINE Distances2D

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: Reflect2D
! !INTERFACE:
  SUBROUTINE Reflect2D( is, HS, BotTop, tBdry, nBdry, kappa, RefC, Npts, &
                        myThid )
! !DESCRIPTION:
!   Reflects a ray at the boundary, using the reflection coefficient
!   and the boundary properties. The ray is reflected in the same x, and
!   basis vectors are updated. The reflected ray is stored in ray2D(is+1).

! !USES:
  USE ssp_mod,  only: evalSSP
  USE bdry_mod, only: HSInfo
  USE refCoef,  only: ReflectionCoef, InterpolateReflectionCoefficient
  USE ihop_mod, only: Beam, ray2D, rad2deg, afreq

! !INPUT PARAMETERS:
! is     :: Step index of the ray to be reflected
! HS     :: Half-space properties at the boundary
! BotTop :: 'TOP' or 'BOT' for top or bottom boundary reflection
! tBdry  :: Tangent vector to the boundary at the reflection point [m, m]
! nBdry  :: Normal vector to the boundary at the reflection point [m, m]
! kappa  :: Boundary curvature, for curvilinear grids
! RefC   :: Reflection coefficient at the boundary
! Npts   :: Number of points in the reflection coefficient array
! myThid :: my thread ID
  INTEGER,              INTENT( INOUT ) :: is
  TYPE( HSInfo ),       INTENT( IN ) :: HS
  CHARACTER (LEN=3),    INTENT( IN ) :: BotTop
  REAL (KIND=_RL90),    INTENT( IN ) :: tBdry(2), nBdry(2), kappa
  TYPE(ReflectionCoef), INTENT( IN ) :: RefC( Npts )
  INTEGER,              INTENT( IN ) :: Npts
  INTEGER,              INTENT( IN ) :: myThid
! !OUTPUT PARAMETERS: is

! !LOCAL VARIABLES:
! msgBuf :: Informational/error message buffer
! is1    :: Step index of the reflected ray
! c, cimag :: Sound speed and its imaginary part
! gradc  :: Gradient of sound speed
! crr, crz, czz :: Radial, vertical, and axial components of sound speed
! rho    :: Density
! RM, RN :: Curvature change parameters
! Tg, Th :: Ray tangent projected along and normal to the boundary
! rayt, rayn :: Unit tangent and normal vectors to the ray
! rayt_tilde, rayn_tilde :: Reflected unit tangent and normal vectors to the ray
! cnjump, csjump :: Jumps in normal and tangential components of the
!                  sound speed gradient
! ck, co, si, cco, ssi :: Parameters for beam shift
! pdelta, rddelta, sddelta :: Parameters for beam shift
! theta_bot :: Beam shift angle for bottom boundary
! kz, kzP, kzS, kzP2, kzS2 :: Complex wave numbers for P and S waves
! mu, f, g, y2, y4 :: Complex coefficients for tabulated reflection coefficient
! Refl :: Reflection coefficient
! ch, a, b, d, sb, delta, ddelta :: Parameters for beam shift
! RInt :: Local derived type for reflection coefficient
  CHARACTER*(MAX_LEN_MBUF) :: msgBuf
  INTEGER              :: is1
  REAL (KIND=_RL90)    :: c, cimag, gradc( 2 ), crr, crz, czz, &
                          rho
  REAL (KIND=_RL90)    :: RM, RN, Tg, Th, rayt( 2 ), rayn( 2 ), &
                          rayt_tilde( 2 ), rayn_tilde( 2 ), cnjump, &
                          csjump
  REAL (KIND=_RL90)    :: ck, co, si, cco, ssi, pdelta, rddelta, sddelta, &
                          theta_bot
  COMPLEX (KIND=_RL90) :: kx, kz, kzP, kzS, kzP2, kzS2, mu, f, g, y2, y4, &
                          Refl
  COMPLEX (KIND=_RL90) :: ch, a, b, d, sb, delta, ddelta
  TYPE(ReflectionCoef) :: RInt

!$TAF init reflect2d1 = 'IHOPreflectray2d'

  ! Init default values for local derived type Rint
  Rint%R = 0.0
  Rint%phi = 0.0
  Rint%theta = -999.0

  ! increment stepping counters
  is  = is + 1 ! old step
  is1 = is + 1 ! new step reflected (same x, updated basis vectors)

!$TAF store ray2D(is)%t = reflect2d1
  Tg = DOT_PRODUCT( ray2D( is )%t, tBdry )  ! ray tan projected along boundary
  Th = DOT_PRODUCT( ray2D( is )%t, nBdry )  ! ray tan projected normal boundary

  ray2D( is1 )%NumTopBnc = ray2D( is )%NumTopBnc
  ray2D( is1 )%NumBotBnc = ray2D( is )%NumBotBnc
  ray2D( is1 )%x         = ray2D( is )%x
  ray2D( is1 )%t         = ray2D( is )%t - 2.0 * Th * nBdry ! change ray direction

  ! Calculate change in curvature, kappa
  ! Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).

  ! Get c
  CALL evalSSP( ray2D( is )%x, c, cimag, gradc, crr, crz, czz, rho, myThid )

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

  IF ( BotTop.EQ.'TOP' ) THEN
    ! cnjump changes sign because the (t,n) system of the top boundary has a
    ! different sense to the bottom boundary
    cnjump = -cnjump
    RN     = -RN
  END IF

  RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
  RN = RN + RM * ( 2 * cnjump - RM * csjump ) / c ** 2

  SELECT CASE ( Beam%Type( 3:3 ) )
  CASE ( 'D' )
    RN = 2.0 * RN
  CASE ( 'Z' )
    RN = 0.0
  CASE DEFAULT
    RN = RN
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
!$TAF store rint = reflect2d1
    ! angle of incidence (relative to normal to bathymetry)
    IF ( (Th.EQ.0.0) .AND. (Tg.EQ.0.0) ) THEN
      RInt%theta = 0.0
    ELSE
      RInt%theta = rad2deg * ABS( ATAN2( Th, Tg ) )
    ENDIF
    
    ! reflection coefficient is symmetric about 90 degrees
    IF ( RInt%theta.GT.90. ) RInt%theta = 180. - RInt%theta 

    CALL InterpolateReflectionCoefficient( RInt, RefC, Npts )
    ray2D( is1 )%Amp   = ray2D( is )%Amp * RInt%R
    ray2D( is1 )%Phase = ray2D( is )%Phase + RInt%phi

  CASE ( 'A', 'G' )     ! half-space
    kx = afreq * Tg    ! wavenumber in direction parallel      to bathymetry
    kz = afreq * Th    ! wavenumber in direction perpendicular to bathymetry

    ! notation below is a bit mis-leading
    ! kzS, kzP is really what I called gamma in other codes, and differs by a
    ! factor of +/- i
    IF ( REAL( HS%cS ).GT.0. ) THEN
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
      IF ( REAL( kzP ).EQ.0. .AND. AIMAG( kzP ).LT.0. ) kzP = -kzP
      f   = kzP
      g   = HS%rho

    ENDIF

    ! complex reflection coef.
    Refl =  - ( rho*f - oneCMPLX * kz*g ) / ( rho*f + oneCMPLX*kz*g )

    IF ( ABS( Refl ).LT.1.0E-5 ) THEN   ! kill a ray that has lost its energy in reflection
      ray2D( is1 )%Amp   = 0.0
      ray2D( is1 )%Phase = ray2D( is )%Phase

    ELSE
      ray2D( is1 )%Amp   = ABS( Refl ) * ray2D(  is )%Amp
      ray2D( is1 )%Phase = ray2D( is )%Phase + &
                            ATAN2( AIMAG( Refl ), REAL( Refl ) )

      IF ( Beam%Type( 4:4 ).EQ.'S' ) THEN   ! beam displacement & width change (Seongil's version)
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

        IF ( si.NE.0. ) THEN
            delta = a * co / si / ( ck * sb * d )   ! Do we need an abs() on this???
        ELSE
            delta = 0.0
        ENDIF

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
                          
      ENDIF ! IF ( Beam%Type( 4:4 ).EQ.'S' )

    ENDIF ! IF ( ABS( Refl ).LT.1.0E-5 )

  CASE DEFAULT
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'HS%BC = ', HS%BC
    ! In adjoint mode we do not write output besides on the first run
    IF ( IHOP_dumpfreq.GE.0 ) &
      CALL PRINT_MESSAGE(msgBuf, PRTFile, SQUEEZE_RIGHT, myThid)
    WRITE(msgBuf,'(A)') 'IHOP Reflect2D: Unknown boundary condition type'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R Reflect2D'
  END SELECT !CASE ( Beam%Type( 3:3 ) )

  ! Update top/bottom bounce counter
  IF (BotTop.EQ.'TOP') THEN
    ray2D( is+1 )%NumTopBnc = ray2D( is )%NumTopBnc + 1
  ELSEIF ( BotTop.EQ.'BOT' ) THEN
    ray2D( is+1 )%NumBotBnc = ray2D( is )%NumBotBnc + 1
  ELSE
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'IHOP Reflect2D: ', &
      'no reflection bounce, but in reflect2d somehow'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R Reflect2D'
    
  ENDIF

  RETURN
  END !SUBROUTINE Reflect2D

END MODULE IHOP
