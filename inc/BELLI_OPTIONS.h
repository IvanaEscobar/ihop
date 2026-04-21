#ifndef BELLI_OPTIONS_H
#define BELLI_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

!BOP
! !ROUTINE: BELLI_OPTIONS.h
! !INTERFACE:
! #include "BELLI_OPTIONS.h"

! !DESCRIPTION:
! CPP options file for belli package:
! Use this file for selecting options within package "belli"
!EOP

#ifdef ALLOW_BELLI
! Place CPP define/undef flag here
#define BELLI_DEBUG

! Allow writing to PRTFile and RAY/DEL/ARR Files
#define BELLI_WRITE_OUT

! Three-dimensional sound propagation
#undef BELLI_THREED

! to reduce memory storage, disable unused array with those CPP flags :
#define BELLI_3D_STATE
#define BELLI_2D_STATE
#undef BELLI_TENDENCY

! a time-series of sound source propagation
#define BELLI_MULTIPLE_TIMES

! multiple sound sources
#undef BELLI_MULTIPLE_SOURCES

! multiple sound receivers; e.g. a VLA
#undef BELLI_MULTIPLE_RECEIVER_DEPTHS
#undef BELLI_MULTIPLE_RECEIVER_RANGES

!!! belli cost options !!!
! QoI as J = C^2 or J = cMat^2 at final timestep
#undef TEST_BELLI_COST
! QoI as J = cMat(t)^2, etc.
#undef TEST_BELLI_COST_INLOOP

#endif /* ALLOW_BELLI */
#endif /* BELLI_OPTIONS_H */

!EH3 ;;; Local Variables: ***
!EH3 ;;; mode:fortran ***
!EH3 ;;; End: ***
