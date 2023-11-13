!     ==================================================================
!     HEADER calendar
!     ==================================================================
!
!     o This header file contains variables that are used by the
!       calendar tool. The calendar tool can be used in the ECCO
!       SEALION release of the MITgcmUV.
!
!     started: Christian Eckert eckert@mit.edu  30-Jun-1999
!     changed: Christian Eckert eckert@mit.edu  17-Dec-1999
!              - restructured the original version in order to have a
!                better interface to the MITgcmUV.
!
!     ==================================================================
!     HEADER calendar
!     ==================================================================

!   - Parameters of the numerical model:
!
!     modelStart       :: start time of the numerical model.
!     modelStartDate   :: start date of the numerical model.
!     modelEnd         :: end   time of the numerical model.
!     modelEndDate     :: end   date of the numerical model.
!     modelStep        :: timestep of the numerical model.
!     modelIntSteps    :: number of timestep that are to be performed.
!     modelIter0       :: the numerical models initial timestep number.
!     modelIterEnd     :: the models last timestep number.
!     modelStepsperday :: number of model time steps per day (<- removed).

!   - Parameters used by the calendar:
!
!     refDate          :: first day of the Gregorian Calendar.
!     nMonthYear       :: number months in a year.
!     nDayMonth        :: days per month depending on the year being a leap
!                         year or not. If the Model calendar is used a 360
!                         days year with 30 days months is used instead.
!     nDaysNoLeap      :: number of days in a usual year.
!     nDaysLeap        :: number of days in a leap year.
!     nMaxDayMonth     :: maximum number of days in a years month.
!     hoursPerDay      :: number of hours   in a calendars day.
!     minutesPerDay    :: number of minutes in a calendars day.
!     minutesPerHour   :: number of minutes in a calendars hour.
!     secondsPerDay    :: number of seconds in a calendars day.
!     secondsPerHour   :: number of seconds in a calendars hour.
!     secondsPerMinute :: number of seconds in a calendars minute.
!     cal_setStatus    :: status of calendar parms setting (0=none, 3=fully set)

      INTEGER nMonthYear
      PARAMETER ( nMonthYear = 12 )

      COMMON /CALENDAR_RL/                                                                                                              &
     &                modelStart,                                                                                                       &
     &                modelEnd,                                                                                                         &
     &                modelStep
      _RL modelStart
      _RL modelEnd
      _RL modelStep

      COMMON /CALENDAR_I/                                                                                                               &
     &               refDate,                                                                                                           &
     &               nDayMonth,                                                                                                         &
     &               nDaysNoLeap,                                                                                                       &
     &               nDaysLeap,                                                                                                         &
     &               nMaxDayMonth,                                                                                                      &
     &               hoursPerDay,                                                                                                       &
     &               minutesPerDay,                                                                                                     &
     &               minutesPerHour,                                                                                                    &
     &               secondsPerDay,                                                                                                     &
     &               secondsPerHour,                                                                                                    &
     &               secondsPerMinute,                                                                                                  &
     &               modelStartDate,                                                                                                    &
     &               modelEndDate,                                                                                                      &
     &               modelIter0,                                                                                                        &
     &               modelIterEnd,                                                                                                      &
     &               modelIntSteps,                                                                                                     &
     &               cal_setStatus,                                                                                                     &
     &               startdate_1,                                                                                                       &
     &               startdate_2

      INTEGER refDate(4)
      INTEGER nDayMonth(nMonthYear,2)
      INTEGER nDaysNoLeap
      INTEGER nDaysLeap
      INTEGER nMaxDayMonth
      INTEGER hoursPerDay
      INTEGER minutesPerDay
      INTEGER minutesPerHour
      INTEGER secondsPerDay
      INTEGER secondsPerHour
      INTEGER secondsPerMinute

      INTEGER modelStartDate(4)
      INTEGER modelEndDate(4)
      INTEGER modelIter0
      INTEGER modelIterEnd
      INTEGER modelIntSteps

      INTEGER cal_setStatus
      INTEGER startdate_1
      INTEGER startdate_2

!   calendarDumps :: When set, approximate months (30-31 days) and years (360-372 days)
!                    for parameters chkPtFreq, pChkPtFreq, taveFreq, SEAICE_taveFreq,
!                    KPP_taveFreq, and freq in pkg/diagnostics are converted to exact
!                    calendar months and years.  Requires pkg/cal.
      COMMON /CALENDAR_L/                                                                                                               &
     &               calendarDumps,                                                                                                     &
     &               usingModelCalendar,                                                                                                &
     &               usingNoLeapYearCal,                                                                                                &
     &               usingJulianCalendar,                                                                                               &
     &               usingGregorianCalendar
      LOGICAL calendarDumps
      LOGICAL usingModelCalendar
      LOGICAL usingNoLeapYearCal
      LOGICAL usingJulianCalendar
      LOGICAL usingGregorianCalendar

!     theCalendar :: type of calendar to use; available:
!                    'model', 'gregorian' or 'noLeapYear'.
!     dayOfWeek   :: Week day number one is the week day of refDate.
!                    For the Gregorian calendar this is Friday, 15-Oct-1582.
!     monthOfYear :: Both available calendars are assumed to have twelve
!                    months.
      COMMON /CALENDAR_C/                                                                                                               &
     &                     theCalendar,                                                                                                 &
     &                     dayOfWeek,                                                                                                   &
     &                     monthOfYear
      CHARACTER*(20) theCalendar
      CHARACTER*(3) dayOfWeek(7)
      CHARACTER*(3) monthOfYear(nMonthYear)

