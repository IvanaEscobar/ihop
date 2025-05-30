!
!BOP
!    !ROUTINE: RECIP_DYG_MACROS.h
!    !INTERFACE:
!    include RECIP_DYG_MACROS.h
!    !DESCRIPTION: \bv
!     *==========================================================*
!     | RECIP_DYG_MACROS.h
!     *==========================================================*
!     | These macros are used to reduce memory requirement and/or
!     | memory references when variables are fixed along a given
!     | axis or axes.
!     *==========================================================*
!     \ev
!EOP

#ifdef RECIP_DYG_CONST
#define  _recip_dyG(i,j,bi,bj) recip_dyG(1,1,1,1)
#endif

#ifdef RECIP_DYG_FX
#define  _recip_dyG(i,j,bi,bj) recip_dyG(i,1,bi,1)
#endif

#ifdef RECIP_DYG_FY
#define  _recip_dyG(i,j,bi,bj) recip_dyG(1,j,1,bj)
#endif

#ifndef _recip_dyG
#define  _recip_dyG(i,j,bi,bj) recip_dyG(i,j,bi,bj)
#endif
