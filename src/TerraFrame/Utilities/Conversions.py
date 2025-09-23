# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import math

from .Time import Deltas
from .Time.JulianDate import JulianDate
from .Time.TimeScales import TimeScales


def any_to_tt(jd_a):
    """
    This function takes a Julian Date (JD) in UTC, TT, or TAI and
    converts it to TT. This convertion is leap second aware.

    Note, UT1 is not supported.

    :param jd_a: Julian Date in UTC, TT, or TAI
    :type jd_a: JulianDate
    :return: Julian Date in TAI
    :rtype: JulianDate
    """

    assert (isinstance(jd_a, JulianDate))

    match jd_a.time_scale:
        case TimeScales.TT:
            return jd_a
        case TimeScales.UTC:
            return utc_to_tt(jd_a)
        case TimeScales.TAI:
            return tai_to_tt(jd_a)
        case _:
            raise RuntimeError(f'Unsupported timescale in convertion to TT: '
                               f'{jd_a.time_scale}')


def utc_to_tt(jd_utc):
    """
    This function takes a Julian Date (JD) in UTC and converts it to TT. This
    convertion is leap second aware.

    :param jd_utc: Julian Date in UTC
    :type jd_utc: JulianDate
    :return: Julian Date in TT
    :rtype: JulianDate
    """

    assert (isinstance(jd_utc, JulianDate))
    assert jd_utc.time_scale == TimeScales.UTC

    jd_tai = utc_to_tai(jd_utc)

    jd_tt = tai_to_tt(jd_tai)
    jd_tt.time_scale = TimeScales.TT

    return jd_tt


def utc_to_tai(jd_utc):
    """
    This function takes a Julian Date (JD) in UTC and converts it to TAI. This
    convertion is leap second aware.

    :param jd_utc: Julian Date in UTC
    :type jd_utc: JulianDate
    :return: Julian Date in TAI
    :rtype: JulianDate
    """

    assert (isinstance(jd_utc, JulianDate))
    assert jd_utc.time_scale == TimeScales.UTC

    delta = Deltas.TaiUtcDelta().get_delta(jd_utc)

    jd_tai = jd_utc + seconds_to_days(delta)
    jd_tai.time_scale = TimeScales.TAI

    return jd_tai


def utc_to_ut1(jd_utc):
    """
    This function takes a Julian Date (JD) in UTC and converts it to UT1. This
    convertion is leap second aware.

    :param jd_utc: Julian Date in UTC
    :type jd_utc: JulianDate
    :return: Julian Date in UT1
    :rtype: JulianDate
    """

    assert (isinstance(jd_utc, JulianDate))
    assert jd_utc.time_scale == TimeScales.UTC

    # We convert to TAI as an intermediate to avoid UTC leap second ambiguity
    jd_tai = utc_to_tai(jd_utc)

    delta_ut1_utc = Deltas.Ut1UtcDelta().get_delta(jd_utc)
    delta_tai_utc = Deltas.TaiUtcDelta().get_delta(jd_utc)

    delta_ut1_tai = delta_ut1_utc - delta_tai_utc

    jd_ut1 = jd_tai + seconds_to_days(delta_ut1_tai)
    jd_ut1.time_scale = TimeScales.UT1

    return jd_ut1


def tai_to_ut1(jd_tai):
    """
    This function takes a Julian Date (JD) in TAI and converts it to UT1. This
    convertion is leap second aware.

    :param jd_tai: Julian Date in TAI
    :type jd_tai: JulianDate
    :return: Julian Date in UT1
    :rtype: JulianDate
    """

    assert (isinstance(jd_tai, JulianDate))
    assert jd_tai.time_scale == TimeScales.TAI

    jd_utc = tai_to_utc(jd_tai)

    delta_ut1_utc = Deltas.Ut1UtcDelta().get_delta(jd_utc)
    delta_tai_utc = Deltas.TaiUtcDelta().get_delta(jd_utc)

    delta_ut1_tai = delta_ut1_utc - delta_tai_utc

    jd_ut1 = jd_tai + seconds_to_days(delta_ut1_tai)
    jd_ut1.time_scale = TimeScales.UT1

    return jd_ut1


def tt_to_ut1(jd_tt):
    """
    This function takes a Julian Date (JD) in TT and converts it to UT1. This
    convertion is leap second aware.

    :param jd_tt: Julian Date in TT
    :type jd_tt: JulianDate
    :return: Julian Date in UT1
    :rtype: JulianDate
    """

    assert (isinstance(jd_tt, JulianDate))
    assert jd_tt.time_scale == TimeScales.TT

    jd_tai = tt_to_tai(jd_tt)

    jd_ut1 = tai_to_ut1(jd_tai)

    return jd_ut1


def tai_to_utc(jd_tai):
    """
    This function takes a Julian Date (JD) in TAI and converts it to UTC.
    This convertion is leap second aware. Though note that in general the UTC
    representation of a leap second is ambiguous.

    :param jd_tai: Julian Date in TAI
    :type jd_tai: JulianDate
    :return: Julian Date in UTC
    :rtype: JulianDate
    """

    assert (isinstance(jd_tai, JulianDate))
    assert jd_tai.time_scale == TimeScales.TAI

    lhs = Deltas.TaiUtcDeltaInverted()

    jd_utc = jd_tai - seconds_to_days(lhs.get_delta(jd_tai))
    jd_utc.time_scale = TimeScales.UTC

    return jd_utc


def tt_to_utc(jd_tt):
    """
    This function takes a Julian Date (JD) in TT and converts it to UTC.
    This convertion is leap second aware. Though note that in general the UTC
    representation of a leap second is ambiguous.

    :param jd_tt: Julian Date in TT
    :type jd_tt: JulianDate
    :return: Julian Date in UTC
    :rtype: JulianDate
    """

    assert (isinstance(jd_tt, JulianDate))
    assert jd_tt.time_scale == TimeScales.TT

    jd_tai = tt_to_tai(jd_tt)

    jd_utc = tai_to_utc(jd_tai)

    return jd_utc


def tt_to_tai(jd_tt):
    """
    This function takes a Julian Date (JD) in TT and converts it to TAI.

    :param jd_tt: Julian Date in TT
    :type jd_tt: JulianDate
    :return: Julian Date in TAI
    :rtype: JulianDate
    """
    assert (isinstance(jd_tt, JulianDate))
    assert jd_tt.time_scale == TimeScales.TT

    jd_tai = jd_tt - seconds_to_days(32.184)
    jd_tai.time_scale = TimeScales.TAI

    return jd_tai


def tai_to_tt(jd_tai):
    """
    This function takes a Julian Date (JD) in TAI and converts it to TT.

    :param jd_tai: Julian Date in TAI
    :type jd_tai: JulianDate
    :return: Julian Date in TT
    :rtype: JulianDate
    """
    assert (isinstance(jd_tai, JulianDate))
    assert jd_tai.time_scale == TimeScales.TAI

    jd_tt = jd_tai + seconds_to_days(32.184)
    jd_tt.time_scale = TimeScales.TT

    return jd_tt


def muas_to_rad(x):
    """
    This function converts the input from microarcsecond (muas) to radians

    :rtype: float
    :type x: float
    :param x: Input muas angle
    :return: Output angle in radians
    """

    val = x / 1e6  # microarcsecond to arcsecond

    return arcsec_to_rad(val)


def mas_to_rad(x):
    """
    This function converts the input from milliarcseconds (mas) to radians

    :rtype: float
    :type x: float
    :param x: Input mas angle
    :return: Output angle in radians
    """

    val = x / 1000.0  # milliarcsecond to arcsecond

    return arcsec_to_rad(val)


def arcsec_to_rad(x):
    """
    This function converts the input from milliarcseconds (mas) to radians

    :rtype: float
    :type x: float
    :param x: Input arcsecond angle
    :return: Output angle in radians
    """

    val = x / 3600.0  # arcsecond to degrees
    val *= math.pi / 180.0  # degrees to radians

    return val


def seconds_to_days(value):
    """
    This function takes a numeric input representing seconds and scales it to
    days. The convertion assumes there is exactly 86,400 seconds in a day.

    :param value: Input value in seconds
    :type value: int | float
    :return: Output value in days
    :rtype: int | float
    """

    return value / 86400


def days_to_centuries(value):
    """
    This function takes a numeric input representing days and scales it to
    centuries. The convertion assumes there is exactly 365 days in a year.

    :param value: Input value in days
    :type value: int | float
    :return: Output value in centuries
    :rtype: int | float
    """

    return value / (365 * 100)


def seconds_to_centuries(value):
    """
    This function takes a numeric input representing seconds and scales it to
    centuries. The convertion assumes there is exactly 365 days in a year.

    :param value: Input value in seconds
    :type value: int | float
    :return: Output value in centuries
    :rtype: int | float
    """

    value = days_to_centuries(seconds_to_days(value))

    return value