# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import erfa
import numpy as np

from TerraFrame.PrecessionNutation import SeriesExpansion
from TerraFrame.Utilities.Time import JulianDate
from TerraFrame.Utilities import Conversions


def test_cip_calculation():
    se_cip_x = SeriesExpansion.cip_x()
    se_cip_y = SeriesExpansion.cip_y()
    se_cip_sxy2 = SeriesExpansion.cip_sxy2()

    n = 250
    frac = np.linspace(0, 365.25 * 0.25, n)

    cip_x = np.zeros((n,))
    cip_y = np.zeros((n,))
    cip_s = np.zeros((n,))
    cip_x_a = np.zeros((n,))
    cip_y_a = np.zeros((n,))
    cip_s_a = np.zeros((n,))

    for i, val in enumerate(frac):
        jd_tt = (JulianDate.JulianDate.j2000(
            time_scale=JulianDate.TimeScales.TT) + val)
        jdc_tt = JulianDate.julian_day_datetime_to_century_datetime(jd_tt)

        cip_x[i] = se_cip_x.compute(jdc_tt)
        cip_y[i] = se_cip_y.compute(jdc_tt)
        sxy2 = se_cip_sxy2.compute(jdc_tt)
        cip_s[i] = sxy2 - cip_x[i] * cip_y[i] / 2.0

        jd1, jd2 = jd_tt.integer_part(), jd_tt.fraction_part()

        # Get X, Y, s using IAU 2006/2000A model
        x, y, s = erfa.xys06a(jd1, jd2)

        cip_x_a[i] = x
        cip_y_a[i] = y
        cip_s_a[i] = s

    assert np.max(np.abs(cip_x - cip_x_a)) < 1e-10
    assert np.max(np.abs(cip_y - cip_y_a)) < 1e-10
    assert np.max(np.abs(cip_s - cip_s_a)) < 1e-10


def test_cip_derivatives():
    se_cip_x = SeriesExpansion.cip_x()
    se_cip_y = SeriesExpansion.cip_y()
    se_cip_sxy2 = SeriesExpansion.cip_sxy2()

    n = 10
    frac = np.linspace(0, 365.25 * 100, n)

    d_cip_x_dt = np.zeros((n,))
    d_cip_x_dt_fd = np.zeros((n,))
    d_cip_y_dt = np.zeros((n,))
    d_cip_y_dt_fd = np.zeros((n,))
    d_cip_s_dt = np.zeros((n,))
    d_cip_s_dt_fd = np.zeros((n,))

    dt = 1e-6

    for i, val in enumerate(frac):
        jd_tt = (JulianDate.JulianDate.j2000(
            time_scale=JulianDate.TimeScales.TT) + val)
        jdc_tt = JulianDate.julian_day_datetime_to_century_datetime(jd_tt)

        _, d_cip_x_dt[i] = se_cip_x.compute(jdc_tt, derivative=True)
        _, d_cip_y_dt[i] = se_cip_y.compute(jdc_tt, derivative=True)
        _, d_cip_s_dt[i] = se_cip_sxy2.compute(jdc_tt, derivative=True)

        f2 = se_cip_x.compute(jdc_tt + 1 * dt)
        f1 = se_cip_x.compute(jdc_tt - 1 * dt)

        d_cip_x_dt_fd[i] = (f2 - f1) / (2 * dt)
        # noinspection PyTypeChecker
        d_cip_x_dt_fd[i] = Conversions.seconds_to_centuries(d_cip_x_dt_fd[i])

        f2 = se_cip_y.compute(jdc_tt + 1 * dt)
        f1 = se_cip_y.compute(jdc_tt - 1 * dt)

        d_cip_y_dt_fd[i] = (f2 - f1) / (2 * dt)
        # noinspection PyTypeChecker
        d_cip_y_dt_fd[i] = Conversions.seconds_to_centuries(d_cip_y_dt_fd[i])

        f2 = se_cip_sxy2.compute(jdc_tt + 1 * dt)
        f1 = se_cip_sxy2.compute(jdc_tt - 1 * dt)

        d_cip_s_dt_fd[i] = (f2 - f1) / (2 * dt)
        # noinspection PyTypeChecker
        d_cip_s_dt_fd[i] = Conversions.seconds_to_centuries(d_cip_s_dt_fd[i])

    assert (np.max(np.abs(d_cip_x_dt - d_cip_x_dt_fd)) < 1e-13)
    assert (np.max(np.abs(d_cip_y_dt - d_cip_y_dt_fd)) < 1e-13)
    assert (np.max(np.abs(d_cip_s_dt - d_cip_s_dt_fd)) < 1e-13)
