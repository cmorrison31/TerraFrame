# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import sys

import numpy as np

import TerraFrame.Earth
from TerraFrame.Utilities import Time, Conversions
from TerraFrame.Utilities.Time.JulianDate import JulianDate


def r1(phi):
    """
    This function computes the R1 rotation matrix. As per Kaplan (2005), R1 is
    defined as:

    [A] rotation matrix to transform column 3-vectors from one cartesian
    coordinate system to another. Final system is formed by rotating original
    system about its own x-axis by angle phi (counterclockwise as viewed from
    the +x direction):

    Source:
    Kaplan, G. H., 2005, U.S. Naval Observatory Circular No. 179 (Washington:
    USNO), page xi

    :param phi: Rotation angle in radians
    :return: R1 matrix
    :type phi: float
    :rtype: np.ndarray
    """

    r = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(phi), np.sin(phi)],
                  [0.0, -np.sin(phi), np.cos(phi)]])

    return r


def dr1dt(phi, dphi_dt):
    """
    This function computes the time derivative of the R2 rotation matrix. As
    per Kaplan (2005), R2 is defined as:

    [A] rotation matrix to transform column 3-vectors from one cartesian
    coordinate system to another. Final system is formed by rotating original
    system about its own y-axis by angle φ (counterclockwise as viewed from
    the +y direction):

    Source:
    Kaplan, G. H., 2005, U.S. Naval Observatory Circular No. 179 (Washington:
    USNO), page xi

    :param phi: Input rotation angle in radians
    :param dphi_dt: The time derivative of phi where dt is in seconds
    :return: dR3dt matrix
    :type phi: float
    :type dphi_dt: float
    :rtype: np.ndarray
    """

    dr_dt = np.array(
        [[0.0, 0.0, 0.0],
         [0.0, -np.sin(phi) * dphi_dt, np.cos(phi) * dphi_dt],
         [0.0, -np.cos(phi) * dphi_dt, -np.sin(phi) * dphi_dt]]
    )

    return dr_dt


def r2(theta):
    """
    This function computes the R2 rotation matrix. As per Kaplan (2005), R2 is
    defined as:

    [A] rotation matrix to transform column 3-vectors from one cartesian
    coordinate system to another. Final system is formed by rotating original
    system about its own y-axis by angle φ (counterclockwise as viewed from
    the +y direction):

    Source:
    Kaplan, G. H., 2005, U.S. Naval Observatory Circular No. 179 (Washington:
    USNO), page xi

    :param theta: Input rotation angle in radians
    :return: R2 matrix
    :type theta: float
    :rtype: np.ndarray
    """

    r = np.array(
        [[np.cos(theta), 0.0, -np.sin(theta)],
         [0.0, 1.0, 0.0],
         [np.sin(theta), 0.0, np.cos(theta)]]
    )

    return r


def dr2dt(theta, dtheta_dt):
    """
    This function computes the time derivative of the R2 rotation matrix. As
    per Kaplan (2005), R2 is defined as:

    [A] rotation matrix to transform column 3-vectors from one cartesian
    coordinate system to another. Final system is formed by rotating original
    system about its own y-axis by angle φ (counterclockwise as viewed from
    the +y direction):

    Source:
    Kaplan, G. H., 2005, U.S. Naval Observatory Circular No. 179 (Washington:
    USNO), page xi

    :param theta: Input rotation angle in radians
    :param dtheta_dt: The time derivative of theta where dt is in seconds
    :return: dR3dt matrix
    :type theta: float
    :type dtheta_dt: float
    :rtype: np.ndarray
    """

    dr_dt = np.array(
        [[-np.sin(theta) * dtheta_dt, 0.0, -np.cos(theta) * dtheta_dt],
         [0.0, 0.0, 0.0],
         [np.cos(theta) * dtheta_dt, 0.0, -np.sin(theta) * dtheta_dt]]
    )

    return dr_dt


def r3(psi):
    """
    This function computes the R3 rotation matrix. As per Kaplan (2005), R3 is
    defined as:

    [A] rotation matrix to transform column 3-vectors from one cartesian
    coordinate system to another. Final system is formed by rotating original
    system about its own z-axis by angle φ (counterclockwise as viewed from
    the +z direction):

    Source:
    Kaplan, G. H., 2005, U.S. Naval Observatory Circular No. 179 (Washington:
    USNO), page xi

    :param psi: Input rotation angle in radians
    :return: R3 matrix
    :type psi: float
    :rtype: np.ndarray
    """

    r = np.array(
        [[np.cos(psi), np.sin(psi), 0.0],
         [-np.sin(psi), np.cos(psi), 0.0],
         [0.0, 0.0, 1.0]]
    )

    return r


def dr3dt(psi, dpsi_dt):
    """
    This function computes the time derivative of the R3 rotation matrix. As
    per Kaplan (2005), R3 is defined as:

    [A] rotation matrix to transform column 3-vectors from one cartesian
    coordinate system to another. Final system is formed by rotating original
    system about its own z-axis by angle φ (counterclockwise as viewed from
    the +z direction):

    Source:
    Kaplan, G. H., 2005, U.S. Naval Observatory Circular No. 179 (Washington:
    USNO), page xi

    :param psi: Input rotation angle in radians
    :param dpsi_dt: The time derivative of psi where dt is in seconds
    :return: dR3dt matrix
    :type psi: float
    :type dpsi_dt: float
    :rtype: np.ndarray
    """

    dr_dt = np.array(
        [[-np.sin(psi) * dpsi_dt, np.cos(psi) * dpsi_dt, 0.0],
         [-np.cos(psi) * dpsi_dt, -np.sin(psi) * dpsi_dt, 0.0],
         [0.0, 0.0, 0.0]]
    )

    return dr_dt


def euler_angles_from_transformation(t_m):
    """
    This function takes a transformation matrix and calculates the corresponding
    Tait–Bryan angles, following z-y′-x″ (intrinsic rotations).

    The angles are technically Tait–Bryan angles but are often called Euler
    angles. This function has been named to align with common usage.

    The Tait–Bryan angles and ordering in this function align with common
    usage in navigation and engineering.

    :param t_m: Transformation matrix
    :return: Array of Tait–Bryan angles, ordered z-y′-x″ (yaw-pitch-roll)
    :type t_m: np.ndarray
    :rtype: np.ndarray
    """

    pitch = np.asin(-t_m[0, 2])
    yaw = np.atan2(t_m[0, 1], t_m[0, 0])
    roll = np.atan2(t_m[1, 2], t_m[2, 2])

    return np.array([yaw, pitch, roll])


def transformation_from_euler(yaw, pitch, roll):
    """
    This function takes in Tait–Bryan angles, following z-y′-x″
    (intrinsic rotations), and creates the corresponding transformation matrix.

    The angles are technically Tait–Bryan angles but are often called Euler
    angles. This function has been named to align with common usage.

    The Tait–Bryan angles and ordering in this function align with common
    usage in navigation and engineering.

    :param yaw: Yaw rotation angle in radians
    :param pitch: Pitch rotation angle in radians
    :param roll: Roll rotation angle in radians
    :return: Transformation matrix
    :type yaw: float
    :type pitch: float
    :type roll: float
    :rtype: np.ndarray
    """

    t_m = r3(yaw) @ r2(pitch) @ r1(roll)

    return t_m


def angle_and_axis_from_transformation(t_m):
    """
    This function takes in a transformation matrix and computes the
    corresponding angle and axis of rotation.

    :param t_m: Transformation matrix
    :return: angle and axis of rotation
    :type t_m: np.ndarray
    :rtype: tuple(float, np.ndarray)
    """

    angle = np.arccos((np.trace(t_m) - 1.0) / 2.0)

    # Avoid dividing by zero if rotation is effectively none.
    if abs(angle) < sys.float_info.epsilon:
        angle = 0.0
        axis = np.array([1.0, 0, 0])

        return angle, axis

    axis = np.array(
        (t_m[2, 1] - t_m[1, 2], t_m[0, 2] - t_m[2, 0], t_m[1, 0] - t_m[0, 1]))

    axis /= (2.0 * np.sin(angle))

    return angle, axis


def calculate_s_prime(time):
    """
    This function computes the Terrestrial Intermediate Origin (TIO) locator
    called s' (or s prime) per IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :type time: JulianDate
    :param time: Terrestrial time measured in Julian centuries.
    :return: s prime
    :rtype: float
    """

    assert (time.time_scale == Time.TimeScales.TT)

    # This is an approximation good for the next century. See section 5.5.2 of
    # IERS Conventions (2010) for more context.
    s_prime = -47e-6 * float(time)

    s_prime = Conversions.arcsec_to_rad(s_prime)

    return s_prime


def calculate_s_prime_derivative():
    """
    This function computes the time derivative of the Terrestrial Intermediate
    Origin (TIO) locator called s' (or s prime) per IERS Conventions (2010).

    Note that technically, per IERS Conventions (2010), the input time should
    be Barycentric Dynamical Time (TDB) but the difference between TDB and TT is
    already small. Additionally, we do not consider effects outside the
    Geocentric Celestial Reference System (GCRS) and the primary driver of the
    TDB vs TT difference is earth's mean anomaly in its orbit. The error from
    this simplication is less than a microarcsecond in nutation.

    :return: s prime time derivative
    :rtype: float
    """

    # This is an approximation good for the next century. See section 5.5.2 of
    # IERS Conventions (2010) for more context.
    dt_s_prime_dt = Conversions.seconds_to_centuries(-47e-6)
    dt_s_prime_dt = Conversions.arcsec_to_rad(dt_s_prime_dt)

    return dt_s_prime_dt


def cirs_to_gcrs(x, y, s):
    """
    This function computes the transformation matrix from the Celestial
    Intermediate Reference System (CIRS) to the Geocentric Celestial
    Reference System (GCRS) per IERS Conventions (2010).

    x and y are coordinates of the Celestial Intermediate Pole (CIP) and s is
    the Celestial Intermediate Origin (CIO) locator parameter which provides
    the position of the CIO on the equator of the CIP.

    :type x: float
    :type y: float
    :type s: float
    :param x: X coordinate of the CIP
    :param y: Y coordinate of the CIP
    :param s: CIO location parameter
    :return: CGRS to CIRS transformation matrix
    :rtype: np.ndarray
    """

    # This should never be true in reality
    assert (1.0 - x ** 2 - y ** 2 > 0.0)

    # e and d formulas from Capitaine (2003)
    e = np.atan2(y, x)
    d = np.atan2(np.sqrt(x ** 2 + y ** 2), np.sqrt(1 - x ** 2 - y ** 2))

    t_gc = r3(-e) @ r2(-d) @ r3(e) @ r3(s)

    return t_gc


def cirs_to_gcrs_derivative(x, y, s, dx_dt, dy_dt, ds_dt):
    """
    This function computes the time derivative of the transformation matrix
    from the Celestial Intermediate Reference System (CIRS) to the Geocentric
    Celestial Reference System (GCRS) per IERS Conventions (2010).

    x and y are coordinates of the Celestial Intermediate Pole (CIP) and s is
    the Celestial Intermediate Origin (CIO) locator parameter which provides
    the position of the CIO on the equator of the CIP.

    :param x: X coordinate of the CIP
    :type x: float
    :param y: Y coordinate of the CIP
    :type y: float
    :param s: CIO location parameter
    :type s: float
    :param dx_dt: X coordinate of the CIP time derivative
    :type dx_dt: float
    :param dy_dt: Y coordinate of the CIP time derivative
    :type dy_dt: float
    :param ds_dt: CIO location parameter time derivative
    :type ds_dt: float
    :return: CGRS to CIRS transformation matrix time derivative
    :rtype: np.ndarray
    """

    # This should never be true in reality
    assert (1.0 - x ** 2 - y ** 2 > 0.0)

    # e and d formulas from Capitaine (2003)
    e = np.atan2(y, x)

    de_dt = (x * dy_dt - y * dx_dt)/(x**2 + y**2)

    d = np.atan2(np.sqrt(x ** 2 + y ** 2), np.sqrt(1 - x ** 2 - y ** 2))

    dd_dt = (np.sqrt(-1 - 1/(x**2 + y**2 - 1)) * (x * dx_dt + y * dy_dt) /
             (x**2 + y**2))

    dt_gc_dt = (r3(-e) @ r2(-d) @ r3(e) @ dr3dt(s, ds_dt) +
                r3(-e) @ r2(-d) @ dr3dt(e, de_dt) @ r3(s) +
                r3(-e) @ dr2dt(-d, -dd_dt) @ r3(e) @ r3(s) +
                dr3dt(-e, -de_dt) @ r2(-d) @ r3(e) @ r3(s))

    return dt_gc_dt


def earth_rotation_matrix(time):
    """
    This function computes the earth rotation matrix at a given datetime in UT1.

    :param time: JulianDate in UT1
    :return: TIRS to ITRS transformation matrix
    :type time: JulianDate
    :rtype: np.ndarray
    """

    assert (time.time_scale == Time.TimeScales.UT1)

    era = TerraFrame.Earth.earth_rotation_angle(time)

    r_era = r3(-era)

    return r_era


def earth_rotation_matrix_derivative(time):
    """
    This function computes the time derivative of the earth rotation matrix
    at a given datetime in UT1.

    :param time: JulianDate in UT1
    :return: TIRS to ITRS transformation matrix
    :type time: JulianDate
    :rtype: np.ndarray
    """

    assert (time.time_scale == Time.TimeScales.UT1)

    era = TerraFrame.Earth.earth_rotation_angle(time)
    dera_dt = TerraFrame.Earth.earth_rotation_angle_derivative()

    dr_era_dt = dr3dt(-era, -dera_dt)

    return dr_era_dt


def itrs_to_tirs(pm_x, pm_y, sp):
    """
    This function computes the transformation matrix from the International
    Terrestrial Reference System (ITRS) to the Terrestrial Intermediate
    Reference System (TIRS) per IERS Conventions (2010).

    pm_x and pm_y are coordinates of polar motion and sp (s') is the
    Terrestrial Intermediate Origin (TIO) locator parameter which provides
    the position of the TIO on the equator of the CIP.

    :type pm_x: float
    :type pm_y: float
    :type sp: float
    :param pm_x: Polar motion x coordinate
    :param pm_y: Polar motion y coordinate
    :param sp: TIO location parameter
    :return: TIRS to ITRS transformation matrix
    :rtype: np.ndarray
    """

    t_ti = (r3(-sp) @ r2(pm_x) @ r1(pm_y))

    return t_ti


def itrs_to_tirs_derivative(pm_x, pm_y, sp, dpm_x_dt, dpm_y_dt, dsp_dt):
    """
    This function computes the transformation matrix from the International
    Terrestrial Reference System (ITRS) to the Terrestrial Intermediate
    Reference System (TIRS) per IERS Conventions (2010).

    pm_x and pm_y are coordinates of polar motion and sp (s') is the
    Terrestrial Intermediate Origin (TIO) locator parameter which provides
    the position of the TIO on the equator of the CIP.

    :type pm_x: float
    :param pm_x: Polar motion x coordinate
    :type pm_y: float
    :param pm_y: Polar motion y coordinate
    :type sp: float
    :param sp: TIO location parameter
    :type dpm_x_dt: float
    :param pm_x: Polar motion x coordinate time derivative
    :type dpm_y_dt: float
    :param pm_y: Polar motion y coordinate time derivative
    :type dsp_dt: float
    :param sp: TIO location parameter time derivative
    :return: TIRS to ITRS transformation matrix time derivative
    :rtype: np.ndarray
    """

    dr1_dt = dr1dt(pm_y, dpm_y_dt)
    dr2_dt = dr2dt(pm_x, dpm_x_dt)
    dr3_dt = dr2dt(-sp, -dsp_dt)

    dt_ti_dt = (r3(-sp) @ r2(pm_x) @ dr1_dt + r3(-sp) @ dr2_dt @ r1(pm_y) +
                dr3_dt @ r2(pm_x) @ r1(pm_y))

    return dt_ti_dt
