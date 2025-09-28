# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import math
from abc import ABC, abstractmethod

from TerraFrame.Utilities import Time, Conversions
from TerraFrame.Utilities.Time import JulianDate


class EarthBase(ABC):
    def __init__(self):
        self._gm = 3.9860050e14  # m^3/s^2 per WGS84
        self._omega = 7.292115e-5  # rad/s per WGS84

    @property
    @abstractmethod
    def semi_major_axis(self):
        pass

    @property
    @abstractmethod
    def flattening(self):
        pass

    @property
    @abstractmethod
    def eccentricity(self):
        pass

    @property
    def gravitational_constant(self):
        return self._gm

    @property
    def mean_angular_velocity(self):
        return self._omega

    def angular_velocity(self, jd_ut1):
        """
        This function takes in a Julian Day (JD) datetime in UT1 and returns
        the angular velocity of the earth in a precessing reference frame per
        WGS84 specification.

        :param jd_ut1: JD in UT1
        :type jd_ut1: JulianDate.JulianDate
        :return:
        :rtype: float
        """

        jd_ut1_c = JulianDate.julian_day_datetime_to_century_datetime(jd_ut1)

        omega_prime = 7.2921151467e-5
        m = 7.086e-12 + 4.3e-15 * jd_ut1_c.hour

        omega_star = omega_prime * self._omega + m  # rad/s

        return omega_star

    def radius_geocentric(self, lat_geocentric):
        e2 = self.eccentricity() ** 2

        # (Radius / a)^2
        roa2 = (1.0 - e2) / (1.0 - e2 * math.cos(lat_geocentric)**2)

        r = self.semi_major_axis() * math.sqrt(roa2)

        return r

    def cartesian_from_lat_lon_alt(self, lat, lon, alt):
        """
        This function takes in a geodetic latitude, longitude, altitude and
        computes the cartesian x, y, and z coordinates. Altitude is height
        above ellipsoid.

        Note that geodetic latitude is not the same as geocentric latitude
        (longitude however, is the same). Geocentric latitude is the angle
        between the position vector from the origin and the equatorial plane.
        Geodetic latitude is the angle between the normal to the spheroid and
        the equatorial plane. If you don't know which you want, you want
        geodetic latitude. If no qualifier is given, 99% of the time latitude
        means geodetic latitude.

        :param lat: Geodetic latitude in radians
        :param lon: Geodetic longitude in radians
        :param alt: Height above ellipsoid in meters
        :return: x, y, z geodetic cartesian coordinates
        :rtype: float, float, float
        """

        n_phi = self.semi_major_axis() / math.sqrt(
            1.0 - self.eccentricity() ** 2 * math.sin(lat) ** 2)

        x = (n_phi + alt) * math.cos(lat) * math.cos(lon)
        y = (n_phi + alt) * math.cos(lat) * math.sin(lon)
        z = ((1.0 - self.eccentricity() ** 2) * n_phi + alt) * math.sin(lat)

        return x, y, z

    @staticmethod
    def cartesian_from_geocentric_lat_lon_radius(lat, lon, radius):
        """
        This function takes in a geocentric latitude, longitude, radius and
        computes the cartesian x, y, and z coordinates. Radius is distance
        from the ellipsoid center.

        Note that geodetic latitude is not the same as geocentric latitude
        (longitude however, is the same). Geocentric latitude is the angle
        between the position vector from the origin and the equatorial plane.
        Geodetic latitude is the angle between the normal to the spheroid and
        the equatorial plane. If you don't know which you want, you want
        geodetic latitude. If no qualifier is given, 99% of the time latitude
        means geodetic latitude.

        :param lat: Geodetic latitude in radians
        :param lon: Geodetic longitude in radians
        :param radius: Distance to the ellipsoid center
        :return: x, y, z geodetic cartesian coordinates
        :rtype: float, float, float
        """

        x = radius * math.cos(lat) * math.cos(lon)
        y = radius * math.cos(lat) * math.sin(lon)
        z = radius * math.sin(lat)

        return x, y, z

    def lat_lon_alt_from_cartesian(self, x, y, z):
        """
        This function takes in cartesian coordinates and calculates the
        geodetic latitude, longitude, and altitude. The algorithm used is
        from Fukushima 2006.

        Note that geodetic latitude is not the same as geocentric latitude
        (longitude however, is the same). Geocentric latitude is the angle
        between the position vector from the origin and the equatorial plane.
        Geodetic latitude is the angle between the normal to the spheroid and
        the equatorial plane. If you don't know which you want, you want
        geodetic latitude. If no qualifier is given, 99% of the time latitude
        means geodetic latitude.

        Reference:
        Fukushima, T. Transformation from Cartesian to Geodetic Coordinates
        Accelerated by Halley’s Method. J Geodesy 79, 689–693 (2006).
        https://doi.org/10.1007/s00190-006-0023-2

        :param x: Geocentric cartesian x coordinate
        :type x: float
        :param y: Geocentric cartesian y coordinate
        :type y: float
        :param z: Geocentric cartesian z coordinate
        :type z: float
        :return: latitude, longitude, altitude (height above ellipsoid)
        rtype: float, float, float
        """

        ec = math.sqrt(1.0 - self.eccentricity() ** 2)
        pl = math.sqrt(x ** 2 + y ** 2)

        a_inv = 1.0 / self.semi_major_axis()

        p = a_inv * pl
        zb = a_inv * ec * abs(z)
        e = self.eccentricity() ** 2

        sn1 = zb
        cn1 = ec * p

        # Just two iterations gets us to the micro-arcsecond accuracy level
        count = 0
        while count < 2:
            s = sn1
            c = cn1

            # It's symbol soup. No way around it. If you want this to make
            # sense, you'll have to look at the reference paper: Fukushima 2006.
            a = math.sqrt(s ** 2 + c ** 2)
            b = 1.5 * e * s * c ** 2 * ((p * s - zb * c) * a - e * s * c)
            f = p * a ** 3 - e * c ** 3
            d = zb * a ** 3 + e * s ** 3

            cn1 = f ** 2 - b * c
            sn1 = d * f - b * s

            count += 1

        s = sn1
        c = cn1

        cc = ec * c

        lon = math.atan2(y, x)
        lat = math.copysign(1.0, z) * math.atan2(s, cc)
        tmp = self.semi_major_axis() * math.sqrt(ec ** 2 * s ** 2 + cc ** 2)
        alt = ((pl * cc + abs(z) * s - tmp) / math.sqrt(cc ** 2 + s ** 2))

        return lat, lon, alt

    @staticmethod
    def geocentric_lat_lon_radius_from_cartesian(x, y, z):
        """
        This function takes in cartesian coordinates and calculates the
        geocentric latitude, longitude, and radius.

        Note that geodetic latitude is not the same as geocentric latitude
        (longitude however, is the same). Geocentric latitude is the angle
        between the position vector from the origin and the equatorial plane.
        Geodetic latitude is the angle between the normal to the spheroid and
        the equatorial plane. If you don't know which you want, you want
        geodetic latitude. If no qualifier is given, 99% of the time latitude
        means geodetic latitude.

        :param x: Geocentric cartesian x coordinate
        :type x: float
        :param y: Geocentric cartesian y coordinate
        :type y: float
        :param z: Geocentric cartesian z coordinate
        :type z: float
        :return: Geocentric latitude, longitude, radius (distance from
        geocentric origin)
        rtype: float, float, float
        """

        lat = math.atan2(z, math.sqrt(x ** 2 + y ** 2))
        lon = math.atan2(y, x)
        radius = math.sqrt(x ** 2 + y ** 2 + z ** 2)

        return lat, lon, radius

    def geodetic_to_geocentric(self, lat, lon, alt=0.0):
        """
        This function takes in geodetic latitude, longitude, and altitude and
        computes the corresponding geocentric latitude, longitude, and radius
        from ellipsoid center.

        :param lat: Geodetic latitude in radians
        :type lat: float
        :param lon: Geodetic longitude in radians
        :type lon: float
        :param alt: Height above ellipsoid in meters
        :type alt: float
        :return: Geocentric latitude, longitude, radius
        rtype: float, float, float
        """

        x, y, z = self.cartesian_from_lat_lon_alt(lat, lon, alt)

        lat_geo, lon_geo, r = self.geocentric_lat_lon_radius_from_cartesian(
            x, y, z)

        return lat_geo, lon_geo, r


class SphericalEarth(EarthBase):
    def __init__(self):
        super().__init__()
        # Mean radius of the three semi-axes per WGS84
        self._r = 6371008.7714  # meters.
        self._f = 0.0
        self._e = 0.0

    def semi_major_axis(self):
        return self._r

    def flattening(self):
        return self._f

    def eccentricity(self):
        return self._e


class WGS84Ellipsoid(EarthBase):
    def __init__(self):
        super().__init__()
        self._a = 6378137.0  # meters, per WGS84
        self._f = 1.0 / 298.257223563  # per WGS84
        self._e = math.sqrt(2 * self._f - self._f ** 2)

    def semi_major_axis(self):
        return self._a

    def flattening(self):
        return self._f

    def eccentricity(self):
        return self._e


def earth_rotation_angle(time):
    """
    This function computes the earth rotation angle at a given datetime in UT1.

    :param time: JulianDate.JulianDate in UT1
    :return: Earth rotation angle in radians
    :type time: JulianDate.JulianDate
    :rtype: float
    """

    assert (time.time_scale == Time.TimeScales.UT1)

    day_frac = time.day_fraction()
    tu = time - JulianDate.JulianDate.j2000()

    era = 2.0 * math.pi * (
            day_frac + 0.7790572732640 + float(0.00273781191135448 * tu))

    era = math.fmod(era, 2.0 * math.pi)

    return era


def earth_rotation_angle_derivative():
    """
    This function computes the time derivative of the earth rotation angle at a
    given datetime in UT1.

    :return: Earth rotation angle in radians/sec
    :rtype: float
    """

    dera_dt = 2.0 * math.pi * 1.00273781191135448
    dera_dt = Conversions.seconds_to_days(dera_dt)

    return dera_dt
