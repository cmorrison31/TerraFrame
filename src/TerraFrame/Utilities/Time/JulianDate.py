# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import math
from copy import copy, deepcopy
from TerraFrame.Utilities.Time.TimeScales import TimeScales
import datetime


class JulianBase:
    """
    This class is an implementation of a Julian date (JD) with an integer part
    and a floating point fractional part to maximize precision at the expense
    of computational speed and complexity for basic math operations.

    This base class does not assume a specific timescale or unit time. This
    will work for any time that can be represented by an integer whole and
    floating point fractional part.

    Conversion to float should be delayed as long as possible to retain the
    extra precision throughout calculations.
    """

    time_scale: TimeScales
    _integer_part: int
    _fraction_part: float

    def __init__(self, integer_part, fraction_part=0.0,
                 time_scale: TimeScales = TimeScales.UTC):
        if not isinstance(fraction_part, (float, int)):
            raise RuntimeError('The fraction part of a JulianDate must '
                               'be a floating point number.')

        self.time_scale = time_scale

        self._integer_part = 0
        self._fraction_part = 0.0

        self._add_number(integer_part)
        self._add_number(fraction_part)

    def integer_part(self):
        return self._integer_part

    def fraction_part(self):
        return self._fraction_part

    def _add_number(self, value):
        if isinstance(value, float):
            fractional_part, integer_part = math.modf(value)

            self._integer_part += int(integer_part)
            self._fraction_part = (math.copysign(fractional_part,
                                                 integer_part) + self._fraction_part)

            # Check if we need to roll over the fractional part
            if abs(self._fraction_part) >= 1.0:
                fractional_part, integer_part = math.modf(self._fraction_part)

                self._integer_part += int(integer_part)

                if fractional_part < 0:
                    self._integer_part -= 1
                    self._fraction_part = 1.0 - fractional_part
                else:
                    self._fraction_part = fractional_part

            if self._fraction_part < 0.0 and self._integer_part >= 0:
                self._integer_part -= 1
                self._fraction_part = 1.0 + self._fraction_part
            else:
                self._fraction_part = abs(self._fraction_part)

        elif isinstance(value, int):
            self._integer_part += value

            return

    def __copy__(self):
        return type(self)(self._integer_part, self._fraction_part,
                          time_scale=self.time_scale)

    def __deepcopy__(self, memo):
        id_self = id(self)
        _copy = memo.get(id_self)
        if _copy is None:
            _copy = type(self)(deepcopy(self._integer_part, memo),
                               deepcopy(self._fraction_part, memo),
                               deepcopy(self.time_scale, memo))
            memo[id_self] = _copy
        return _copy

    def __add__(self, other):
        if isinstance(other, float) or isinstance(other, int):
            value = copy(self)
            value._add_number(other)
            return value

        elif isinstance(other, JulianBase):
            value = copy(self)
            value._add_number(other._integer_part)
            value._add_number(other._fraction_part)
            return value

        else:
            return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, float):
            value = copy(self)
            value._add_number(other * -1.0)
            return value

        elif isinstance(other, int):
            value = copy(self)
            value._add_number(other * -1)
            return value

        elif isinstance(other, JulianBase):
            value = copy(self)
            value._add_number(other._integer_part * -1)
            value._add_number(other._fraction_part * -1.0)
            return value

        else:
            return NotImplemented

    def __rsub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return JulianDate(self._integer_part * other,
                              self._fraction_part * other, self.time_scale)

        elif isinstance(other, JulianBase):
            jd = JulianDate(0, 0.0, self.time_scale)

            frac1 = math.copysign(self._fraction_part, self._integer_part)
            frac2 = math.copysign(other._fraction_part, other._integer_part)

            jd._add_number(self._integer_part * other._integer_part)
            jd._add_number(frac1 * frac2)
            jd._add_number(self._integer_part * frac2)
            jd._add_number(frac1 * other._integer_part)

            return jd

        else:
            return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __float__(self):
        return float(self._integer_part) + self._fraction_part

    def __str__(self):
        return (f'{self._integer_part}' + '{:.6f}'.format(
            self._fraction_part).lstrip('-0')).rstrip('0')

    def __repr__(self):
        return ((f'{self._integer_part}' + '{:.17f}'.format(
            self._fraction_part).lstrip('-0')).rstrip('0') +
                f' {self.time_scale.name}')

    def __lt__(self, other):
        assert (isinstance(other, JulianBase))

        if self._integer_part < other._integer_part:
            return True
        elif self._integer_part > other._integer_part:
            return False
        else:
            if self._fraction_part < other._fraction_part:
                return True
            else:
                return False

    def __le__(self, other):
        assert (isinstance(other, JulianBase))

        if self._integer_part < other._integer_part:
            return True
        elif self._integer_part > other._integer_part:
            return False
        else:
            if self._fraction_part <= other._fraction_part:
                return True
            else:
                return False

    def __hash__(self):
        return hash((self._integer_part, self._fraction_part))

    def __eq__(self, other):
        assert (isinstance(other, JulianBase))

        return (self._integer_part == other._integer_part and
                self._fraction_part == other._fraction_part)


class JulianDate(JulianBase):
    """
    This class is a realization of a Julian Day Date (JD) using the JulianBase
    class.

    A Julian date is defined as the number of days since January 1, 4713, B.C.
    at 12 noon. It's very important to note that 0, and therefore whole day
    numbers, occur at noon and not midnight.
    """

    def __init__(self, integer_part, fraction_part=0.0,
                 time_scale=TimeScales.UTC):
        super().__init__(integer_part, fraction_part, time_scale)

    def day_fraction(self):
        return self._fraction_part

    def round_to_days(self):
        # Floor round to the nearest day keeping in mind the noon epoch
        new_integer = self._integer_part
        new_fraction = self._fraction_part

        if new_fraction < 0.5:
            new_integer -= 1

        new_fraction = 0.5

        return JulianDate(new_integer, new_fraction)


    def next_gregorian_day(self):
        new_integer = self._integer_part
        new_fraction = self._fraction_part

        if new_fraction< 0.5:
            new_fraction = 0.5
        else:
            new_integer += 1.0
            new_fraction = 0.5

        return JulianDate(new_integer, new_fraction)

    @staticmethod
    def j2000(time_scale=TimeScales.UTC):
        """
        This function returns the JD for the J2000 epoch

        The J2000 epoch is defined as JD 2451545.0 exactly.

        :param time_scale: Timescale to use for epoch value. Defaults to UTC.
        :type time_scale: TimeScales
        :return: J2000 epoch as a Julian Date
        :rtype: JulianDate
        """

        return JulianDate(2451545, 0.0,
                          time_scale=time_scale)


class JulianCentury(JulianBase):
    """
    This class is a realisation of a Julian Century Date (JC) using the
    JulianBase class.

    A Julian Century is defined as having exactly 36525 days. Just as with
    Julian Dates, the epoch is at January 1, 4713, B.C., 12 noon.
    """

    def __init__(self, integer_part, fraction_part=0.0):
        super().__init__(integer_part, fraction_part, TimeScales.TT)

    def century_fraction(self):
        return self._fraction_part

    def round_to_century(self):
        # Floor round to the nearest day keeping in mind the noon epoch
        if self._fraction_part < 0.5:
            self._fraction_part = 0.0
        else:
            self._fraction_part = 0.5

        return JulianDate(self._integer_part, self._fraction_part)

    @staticmethod
    def j2000():
        """
        This function returns a Julian century object at the J2000 epoch

        :return: Julian century date at J2000 epoch
        :rtype: JulianCentury
        """

        return JulianCentury(67, 0.11964407939767340849357)


def julian_date_now():
    dt_now = datetime.datetime.now(datetime.UTC)
    dt_ref = datetime.datetime(2000, 1, 1, 12, 0, 0, tzinfo=datetime.UTC)

    delta_days = (dt_now - dt_ref).total_seconds() / 86400.0
    jd = JulianDate.j2000() + JulianDate(delta_days)

    return jd


def julian_day_number_from_date(year, month, day):
    """
    This function converts a datetime object to a Julian Days Date (JD).

    The algorthm comes from Explanatory Supplement to the Astronomical Almanac,
    3rd. See section 15.11.3, page 618.

    :param year: Gregorian year
    :type year: int
    :param month: Gregorian month
    :type month: int
    :param day: Gregorian day
    :type day: int
    :return jdn: Julian day number in the same timescale as the input
    :rtype: int
    """

    m = 2
    n = 12
    r = 4
    p = 1461
    y = 4716
    q = 0
    j = 1401
    s = 153
    u = 5
    t = 2
    a = 184
    c = -38

    h = month - m
    g = year + y - (n - h) // n
    f = (h - 1 + n) % n
    e = (p * g + q) // r + day - 1 - j
    jdn = e + (s * f + t) // u
    jdn = jdn - (3 * ((g + a) // 100)) // 4 - c

    return jdn


def julian_date_from_datetime(year, month, day, hour=0, minute=0, second=0,
                              microsecond=0, time_scale=TimeScales.UTC):
    """
    This function converts a datetime object to a Julian Days Date (JD).

    The algorthm comes from Explanatory Supplement to the Astronomical Almanac,
    3rd. See section 15.11.3, page 618.

    :param year: Gregorian year
    :type year: int
    :param month: Gregorian month
    :type month: int
    :param day: Day of the month
    :type day: int
    :param hour: hour
    :type hour: int
    :param minute: minute
    :type minute: int
    :param second: second
    :type second: int
    :param microsecond: microseconds, range 0 -> 1000000
    :type microsecond: int
    :return jdn: Julian day number in the same timescale as the input
    :param time_scale: Timescale to use. Defaults to UTC.
    type time_scale: TimeScale
    :rtype: JulianDate
    """

    jdn = julian_day_number_from_date(year, month, day)

    # Add fractional part
    time_seconds = (hour * 3600 + minute * 60 + second + microsecond / 1e6)
    fractional_day = time_seconds / 86400

    # Shift epoch to noon while maximizing precision
    if fractional_day > 0.5:
        fractional_day -= 0.5
    elif fractional_day == 0.5:
        fractional_day = 0
    else:
        jdn -= 1.0
        fractional_day += 0.5

    return JulianDate(jdn, fractional_day, time_scale=time_scale)


def pydatetime_from_julian_date(jd):
    """

    :param jd:
    :type jd: JulianDate
    :return:
    """
    m = 2
    n = 12
    r = 4
    p = 1461
    y = 4716
    j = 1401
    s = 153
    v = 3
    u = 5
    b = 274277
    c = -38
    w = 2

    jd = jd + 0.5 # Account for noon epoch in JD
    frac = jd.fraction_part()
    jdn = jd.integer_part()

    f = jdn + j
    f = f + (((4 * jdn + b) // 146097) * 3) // 4 + c
    e = r * f + v
    g = e % p // r
    h = u * g + w
    days = h % s // u + 1
    months = (h // s + m) % n + 1
    years = e // p - y + (n + m - months) // n

    frac *= 24
    hours = math.floor(frac)
    frac -= hours
    frac *= 60
    minutes = math.floor(frac)
    frac -= minutes
    frac *= 60
    seconds = math.floor(frac)
    frac -= seconds
    frac *= 1e6
    microseconds = math.floor(frac)

    dt = datetime.datetime(year=years, month=months, day=days, hour=hours,
                           minute=minutes, second=seconds,
                           microsecond=microseconds)

    return dt


def julian_date_from_pydatetime(dt):
    """
    Converts a datetime object to a Julian Date (JD) object. See the called
    function for more information.

    :return: Julian Date

    :type dt: datetime.datetime
    :rtype: JulianDate
    """

    return julian_date_from_datetime(dt.year, dt.month, dt.day, dt.hour,
                                     dt.minute, dt.second, dt.microsecond)


def julian_day_datetime_to_century_datetime(jd):
    """
    This function takes in a Julian Day (JD) datetime and converts it to a
    Julian century datetime based on the J2000 epoch.

    :type jd: JulianDate
    :param jd: Datetime as a JD
    :return:
    """

    t = 1.0 / 36525.0 * (jd - JulianDate.j2000(time_scale=jd.time_scale))

    return t


def julian_date_to_modified_julian_date(jd):
    """
    This function takes in a Julian Date (JD) and converts it to a
    Modified Julian Date (MJD).

    :type jd: JulianDate
    :param jd: Time as a JD
    :return: MJD
    :rtype: JulianDate
    """

    mjd = jd - JulianDate(2400000.0, 0.5,
                          time_scale=jd.time_scale)

    return mjd


def modified_julian_date_to_julian_date(mjd):
    """
    This function takes in a Julian Date (JD) and converts it to a
    Modified Julian Date (MJD).

    :type mjd: JulianDate
    :param mjd: Time as a MJD
    :return: JD
    :rtype: JulianDate
    """

    jd = mjd + JulianDate(2400000, 0.5,
                          time_scale=mjd.time_scale)

    return jd
