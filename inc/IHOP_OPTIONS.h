#ifndef IHOP_OPTIONS_H
#define IHOP_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

!BOP
! !ROUTINE: IHOP_OPTINOS.h
! !INTERFACE:
! #include "IHOP_OPTIONS.h"

! !DESCRIPTION:
! CPP options file for IHOP package:
! Use this file for selecting options within package "ihop"
!EOP

#ifdef ALLOW_IHOP
! Place CPP define/undef flag here
#define IHOP_DEBUG

#define IHOP_WRITE_OUT

! Only consider 2D propagation
#undef IHOP_THREED

! to reduce memory storage, disable unused array with those CPP flags :
#define IHOP_3D_STATE
#define IHOP_2D_STATE
#define IHOP_TENDENCY

#undef IHOP_MULTIPLE_SOURCES

#undef IHOP_MULTIPLE_RECEIVER_DEPTHS
#undef IHOP_MULTIPLE_RECEIVER_RANGES

#endif /* ALLOW_IHOP */
#endif /* IHOP_OPTIONS_H */

!EH3 ;;; Local Variables: ***
!EH3 ;;; mode:fortran ***
!EH3 ;;; End: ***
