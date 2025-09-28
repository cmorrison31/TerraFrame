# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import math

from TerraFrame.Utilities import Conversions


def mean_anomaly_of_the_moon(time):
    """
    This function computes the mean anomaly of the Moon per IERS Conventions
    (2010)

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean anomaly of the Moon
    :rtype: float
    """

    # The polynomial is in arcseconds
    value = (485868.24903600005 + 1717915923.217800 * time +
             31.879200 * time ** 2 + 0.05163500 * time ** 3 -
             0.0002447000 * time ** 4)

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 360 * 60 * 60)

    value = Conversions.arcsec_to_rad(value)

    return value


def mean_anomaly_of_the_moon_derivative(time):
    """
    This function computes the mean anomaly of the Moon time derivative per
    IERS Conventions (2010)

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean anomaly of the Moon time derivative
    :rtype: float
    """

    # The polynomial is in arcseconds
    value = (1717915923.2178 + 63.7584 * time + 0.15490500000000001 *
             time ** 2 - 0.0009788 * time ** 3)

    value = Conversions.arcsec_to_rad(value)

    return value


def mean_anomaly_of_the_sun(time):
    """
    This function computes the mean anomaly of the Sun per IERS Conventions
    (2010)

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean anomaly of the Sun
    :rtype: float
    """

    # The polynomial is in arcseconds
    value = (1287104.793048 + 129596581.0481 * time - 0.5532 * time ** 2 +
             0.0001360 * time ** 3 - 0.00001149 * time ** 4)

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 360 * 60 * 60)

    value = Conversions.arcsec_to_rad(value)

    return value


def mean_anomaly_of_the_sun_derivative(time):
    """
    This function computes the mean anomaly of the Sun time derivative per
    IERS Conventions (2010)

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean anomaly of the Sun time derivative
    :rtype: float
    """

    # The polynomial is in arcseconds
    value = (129596581.0481 - 1.1064 * time + 0.000408 * time ** 2 -
             0.00004596 * time ** 3)

    value = Conversions.arcsec_to_rad(value)

    return value


def mean_longitude_moon_minus_ascending_node(time):
    """
    This function computes Mean longitude of the Moon minus the ascending node
    per IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of the Moon minus the ascending node
    :rtype: float
    """

    # The polynomial is in arcseconds
    value = (335779.526232 + 1739527262.8478 * time - 12.75120 * time ** 2 -
             0.001037 * time ** 3 + 0.00000417 * time ** 4)

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 360 * 60 * 60)

    value = Conversions.arcsec_to_rad(value)

    return value


def mean_longitude_moon_minus_ascending_node_derivative(time):
    """
    This function computes Mean longitude of the Moon minus the ascending node
    time derivative per IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of the Moon minus the ascending node time derivative
    :rtype: float
    """

    # The polynomial is in arcseconds
    value = (1739527262.8478-25.5024 * time - 0.003111 * time ** 2 +
             0.00001668 * time ** 3)

    value = Conversions.arcsec_to_rad(value)

    return value


def mean_elongation_of_the_moon_from_the_sun(time):
    """
    This function computes Mean Elongation of the Moon from the Sun per IERS
    Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean Elongation of the Moon from the Sun
    :rtype: float
    """

    # The polynomial is in arcseconds
    value = (1072260.7036920001 + 1602961601.2090 * time -
             6.3706 * time ** 2 + 0.006593 * time ** 3 - 0.00003169 * time ** 4)

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 360 * 60 * 60)

    value = Conversions.arcsec_to_rad(value)

    return value


def mean_elongation_of_the_moon_from_the_sun_derivative(time):
    """
    This function computes Mean Elongation of the Moon from the Sun time
    derivative per IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean Elongation of the Moon from the Sun time derivative
    :rtype: float
    """

    # The polynomial is in arcseconds
    value = (1602961601.209 - 12.7412 * time + 0.019779 * time ** 2 -
    0.00012676 * time ** 3)

    value = Conversions.arcsec_to_rad(value)

    return value


def mean_longitude_of_the_ascending_node_of_the_moon(time):
    """
    This function computes Mean Longitude of the Ascending Node of the Moon per
    IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of the Ascending Node of the Moon
    :rtype: float
    """

    # The polynomial is in arcseconds
    value = (450160.39803599997 - 6962890.5431 * time + 7.4722 * time ** 2 +
             0.0077020 * time ** 3 - 0.00005939 * time ** 4)

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 360 * 60 * 60)

    value = Conversions.arcsec_to_rad(value)

    return value


def mean_longitude_of_the_ascending_node_of_the_moon_derivative(time):
    """
    This function computes Mean Longitude of the Ascending Node of the Moon
    time derivative per IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean Longitude of the Ascending Node of the Moon time derivative
    :rtype: float
    """

    # The polynomial is in arcseconds
    value = (-6962890.5431 + 14.9444 * time + 0.023106 * time ** 2 -
             0.00023756 * time ** 3)

    value = Conversions.arcsec_to_rad(value)

    return value


def mean_longitude_of_mercury(time):
    """
    This function computes the mean longitude of Mercury per IERS Conventions
    (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of Mercury
    :rtype: float
    """

    # The polynomial is in radians
    value = 4.402608842 + 2608.7903141574 * time

    value = math.fmod(value, 2.0 * math.pi)

    return value


def mean_longitude_of_mercury_derivative():
    """
    This function computes the mean longitude of Mercury time derivative per
    IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :return: Mean longitude of Mercury time derivative
    :rtype: float
    """

    # The polynomial is in radians
    value = 2608.7903141574

    return value


def mean_longitude_of_venus(time):
    """
    This function computes the mean longitude of Venus per IERS Conventions
    (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of Venus
    :rtype: float
    """

    # The polynomial is in radians
    value = 3.176146697 + 1021.3285546211 * time

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 2.0 * math.pi)

    return value


def mean_longitude_of_venus_derivative():
    """
    This function computes the mean longitude of Venus time derivative per IERS
    Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :return: Mean longitude of Venus time derivative
    :rtype: float
    """

    # The polynomial is in radians
    value = 1021.3285546211

    return value


def mean_longitude_of_earth(time):
    """
    This function computes the mean longitude of Earth per IERS Conventions
    (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of Earth
    :rtype: float
    """

    # The polynomial is in radians
    value = 1.753470314 + 628.3075849991 * time

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 2.0 * math.pi)

    return value


def mean_longitude_of_earth_derivative():
    """
    This function computes the mean longitude of Earth time derivative per IERS
    Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :return: Mean longitude of Earth time derivative
    :rtype: float
    """

    # The polynomial is in radians
    value = 628.3075849991

    return value


def mean_longitude_of_mars(time):
    """
    This function computes the mean longitude of Mars per IERS Conventions
    (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of Mars
    :rtype: float
    """

    # The polynomial is in radians
    value = 6.203480913 + 334.0612426700 * time

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 2.0 * math.pi)

    return value


def mean_longitude_of_mars_derivative():
    """
    This function computes the mean longitude of Mars time derivative per IERS
    Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :return: Mean longitude of Mars time derivative
    :rtype: float
    """

    # The polynomial is in radians
    value = 334.0612426700

    return value


def mean_longitude_of_jupiter(time):
    """
    This function computes the mean longitude of Jupiter per IERS Conventions
    (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of Jupiter
    :rtype: float
    """

    # The polynomial is in radians
    value = 0.599546497 + 52.9690962641 * time

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 2.0 * math.pi)

    return value


def mean_longitude_of_jupiter_derivative():
    """
    This function computes the mean longitude of Jupiter time derivative per
    IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :return: Mean longitude of Jupiter time derivative
    :rtype: float
    """

    # The polynomial is in radians
    value = 52.9690962641

    return value


def mean_longitude_of_saturn(time):
    """
    This function computes the mean longitude of Saturn per IERS Conventions
    (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of Saturn
    :rtype: float
    """

    # The polynomial is in radians
    value = 0.874016757 + 21.3299104960 * time

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 2.0 * math.pi)

    return value


def mean_longitude_of_saturn_derivative():
    """
    This function computes the mean longitude of Saturn time derivative per
    IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :return: Mean longitude of Saturn time derivative
    :rtype: float
    """

    # The polynomial is in radians
    value = 21.3299104960

    return value


def mean_longitude_of_uranus(time):
    """
    This function computes the mean longitude of Uranus per IERS Conventions
    (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of Uranus
    :rtype: float
    """

    # The polynomial is in radians
    value = 5.481293872 + 7.4781598567 * time

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 2.0 * math.pi)

    return value


def mean_longitude_of_uranus_derivative():
    """
    This function computes the mean longitude of Uranus time derivative per
    IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :return: Mean longitude of Uranus time derivative
    :rtype: float
    """

    # The polynomial is in radians
    value = 7.4781598567

    return value


def mean_longitude_of_neptune(time):
    """
    This function computes the mean longitude of Neptune per IERS Conventions
    (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: Mean longitude of Neptune
    :rtype: float
    """

    # The polynomial is in radians
    value = 5.311886287 + 3.8133035638 * time

    # Take the modulus before convertion to maintain accuracy
    value = math.fmod(value, 2.0 * math.pi)

    return value


def mean_longitude_of_neptune_derivative():
    """
    This function computes the mean longitude of Neptune time derivative per
    IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :return: Mean longitude of Neptune time derivative
    :rtype: float
    """

    # The polynomial is in radians
    value = 3.8133035638

    return value


def general_precession_in_longitude(time):
    """
    This function computes the general precession in longitude per IERS
    Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: General precession in longitude
    :rtype: float
    """

    # The polynomial is in radians
    value = 0.02438175 * time + 0.00000538691 * time ** 2

    return value


def general_precession_in_longitude_derivative(time):
    """
    This function computes the general precession time derivative in longitude
    per IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: float
    :param time: Terrestrial time measured in Julian centuries.
    :return: General precession in longitude time derivative
    :rtype: float
    """

    # The polynomial is in radians
    value = 0.02438175 + 0.00001077382 * time

    return value
