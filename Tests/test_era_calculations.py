# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import random

import erfa
import numpy as np

import TerraFrame.Earth
from TerraFrame.Utilities import TransformationMatrices, Conversions
from TerraFrame.Utilities.Time import JulianDate


def test_era_angle_calculation():
    val = random.uniform(0, 100.0)
    jd_ut1 = (JulianDate.JulianDate.j2000(
        time_scale=JulianDate.TimeScales.UT1) + val)

    era = TerraFrame.Earth.earth_rotation_angle(jd_ut1)

    jd1, jd2 = jd_ut1.integer_part(), jd_ut1.fraction_part()
    era_a = erfa.era00(jd1, jd2)

    assert (np.abs(era - era_a) < 1e-10)


def test_era_transformation_calculation():
    val = random.uniform(0, 100.0)
    jd_ut1 = (JulianDate.JulianDate.j2000(
        time_scale=JulianDate.TimeScales.UT1) + val)

    t_e = TransformationMatrices.earth_rotation_matrix(jd_ut1)

    jd1, jd2 = jd_ut1.integer_part(), jd_ut1.fraction_part()
    era_a = erfa.era00(jd1, jd2)

    t_e_erfa = TransformationMatrices.r3(-era_a)

    assert (np.max(np.abs(t_e - t_e_erfa)) < 1e-10)


def test_era_transformation_derivative_calculation():
    val = random.uniform(0, 100.0)
    jd_ut1 = (JulianDate.JulianDate.j2000(
        time_scale=JulianDate.TimeScales.UT1) + val)

    dt = 1e-6

    t_e2 = TransformationMatrices.earth_rotation_matrix(jd_ut1 + dt)
    t_e1 = TransformationMatrices.earth_rotation_matrix(jd_ut1 - dt)

    d_te_dt_fd = (t_e2 - t_e1) / (2 * dt)
    d_te_dt_fd *= Conversions.seconds_to_days(1)

    d_te_dt = TransformationMatrices.earth_rotation_matrix_derivative(jd_ut1)

    error = np.abs(d_te_dt - d_te_dt_fd)

    assert (np.max(error) < 1e-13)
