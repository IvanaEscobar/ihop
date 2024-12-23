!
!BOP
!    !ROUTINE: HFACS_MACROS.h
!    !INTERFACE:
!    include HFACS_MACROS.h
!    !DESCRIPTION: \bv
!     *==========================================================*
!     | HFACS_MACROS.h
!     *==========================================================*
!     | These macros are used to reduce memory requirement and/or
!     | memory references when variables are fixed along a given
!     | axis or axes.
!     *==========================================================*
!     \ev
!EOP

#ifdef HFACS_CONST
#define  _hFacS(i,j,k,bi,bj) hFacS(1,1,1,1,1)
#endif

#ifdef HFACS_FX
#define  _hFacS(i,j,k,bi,bj) hFacS(i,1,1,bi,1)
#endif

#ifdef HFACS_FY
#define  _hFacS(i,j,k,bi,bj) hFacS(1,j,1,1,bj)
#endif

#ifdef HFACS_FXY
#define  _hFacS(i,j,k,bi,bj) hFacS(i,j,1,bi,bj)
#endif

#ifdef ALLOW_DEPTH_CONTROL
# define _hFacS(i,j,k,bi,bj) hFacS(i,j,k,bi,bj)*maskS(i,j,k,bi,bj)
#endif

#ifndef _hFacS
#define  _hFacS(i,j,k,bi,bj) hFacS(i,j,k,bi,bj)
#endif
