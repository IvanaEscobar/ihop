!BOP
!     !ROUTINE: CTRL_FIELDS.h
!     !INTERFACE:
!     include "CTRL_FIELDS.h"
!     !DESCRIPTION:
!     \bv
!     *==============================================================*
!     | CTRL_FIELDS.h
!     | o Additional control fields (from main model or pkgs)
!     |   which are only used when specific CTRL options are defined.
!     *==============================================================*
!     | Note: Other model or pkg variables which can also be used
!     | independently from CTRL options remain in their respective
!     | header files.
!     *==============================================================*
!     \ev
!EOP

#ifdef ALLOW_BOTTOMDRAG_CONTROL
      COMMON /CTRL_FIELDS_BOTTOMDRAG/
     &                       bottomDragFld
      _RL  bottomDragFld (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
