!BOP
!     !ROUTINE: CTRL_SIZE.h
!     !INTERFACE:
!     #include "CTRL_SIZE.h"

!     !DESCRIPTION:
!     *================================================================*
!     | CTRL_SIZE.h
!     | o set maximum number of control variables
!     *================================================================*
!EOP

!     Generic control variable array dimension
!     ----------------------------------------
!
!     maxCtrlArr2D :: number of 2-D generic init. ctrl variables
!     maxCtrlArr3D :: number of 3-D generic init. ctrl variables
!     maxCtrlTim2D :: number of 2-D generic time-varying ctrl variables
!     maxCtrlProc  :: number of pre-processing options per ctrl variable

      INTEGER     maxCtrlArr2D
      PARAMETER ( maxCtrlArr2D = 1 )

      INTEGER     maxCtrlArr3D
      PARAMETER ( maxCtrlArr3D = 8 )

      INTEGER     maxCtrlTim2D
      PARAMETER ( maxCtrlTim2D = 1 )

      INTEGER     maxCtrlProc
      PARAMETER ( maxCtrlProc = 1 )

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
