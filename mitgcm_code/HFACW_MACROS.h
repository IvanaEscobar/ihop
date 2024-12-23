!
!BOP
!    !ROUTINE: HFACW_MACROS.h
!    !INTERFACE:
!    include HFACW_MACROS.h
!    !DESCRIPTION: \bv
!     *==========================================================*
!     | HFACW_MACROS.h
!     *==========================================================*
!     | These macros are used to reduce memory requirement and/or
!     | memory references when variables are fixed along a given
!     | axis or axes.
!     *==========================================================*
!     \ev
!EOP

#ifdef HFACW_CONST
#define  _hFacW(i,j,k,bi,bj) hFacW(1,1,1,1,1)
#endif

#ifdef HFACW_FX
#define  _hFacW(i,j,k,bi,bj) hFacW(i,1,1,bi,1)
#endif

#ifdef HFACW_FY
#define  _hFacW(i,j,k,bi,bj) hFacW(1,j,1,1,bj)
#endif

#ifdef HFACW_FXY
#define  _hFacW(i,j,k,bi,bj) hFacW(i,j,1,bi,bj)
#endif

#ifdef ALLOW_DEPTH_CONTROL
# define _hFacW(i,j,k,bi,bj) hFacW(i,j,k,bi,bj)*maskW(i,j,k,bi,bj)
#endif

#ifndef _hFacW
#define  _hFacW(i,j,k,bi,bj) hFacW(i,j,k,bi,bj)
#endif
