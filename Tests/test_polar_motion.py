# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import random

import erfa
import numpy as np

from TerraFrame.Utilities import (TransformationMatrices, BulletinData,
                                  Conversions)
from TerraFrame.Utilities.Time import JulianDate


def test_s_prime():
    val = random.uniform(0, 100.0)
    jd_tt = (JulianDate.JulianDate.j2000(
        time_scale=JulianDate.TimeScales.TT) + val)
    jdc_tt = JulianDate.julian_day_datetime_to_century_datetime(jd_tt)

    s_p = TransformationMatrices.calculate_s_prime(jdc_tt)

    jd1, jd2 = jd_tt.integer_part(), jd_tt.fraction_part()
    s_p_a = erfa.sp00(jd1, jd2)

    assert (abs(s_p - s_p_a) < 1e-10)


def test_itrs_to_tirs():
    val = random.uniform(0, 100.0)
    pm_x = random.uniform(0, 1e-3)
    pm_y = random.uniform(0, 1e-3)

    jd_tt = (JulianDate.JulianDate.j2000(
        time_scale=JulianDate.TimeScales.TT) + val)
    jdc_tt = JulianDate.julian_day_datetime_to_century_datetime(jd_tt)

    s_p = TransformationMatrices.calculate_s_prime(jdc_tt)

    t_ti = TransformationMatrices.itrs_to_tirs(pm_x, pm_y, s_p)

    jd1, jd2 = jd_tt.integer_part(), jd_tt.fraction_part()
    s_p_a = erfa.sp00(jd1, jd2)

    # ERFA/SOFA computes the inverse transform, so we need to take the transpose
    t_it_erfa = erfa.pom00(pm_x, pm_y, s_p_a)
    t_ti_erfa = t_it_erfa.T

    assert (np.max(np.abs(t_ti - t_ti_erfa)) < 1e-10)


def test_itrs_to_tirs_derivative_calculation():
    val = random.uniform(0, 100.0)
    dt = 1e-6
    dtc = Conversions.seconds_to_centuries(dt)
    bd = BulletinData.BulletinData()

    jd_utc = (JulianDate.JulianDate.j2000(
        time_scale=JulianDate.TimeScales.UTC) + val)
    mjd_utc = JulianDate.julian_date_to_modified_julian_date(jd_utc)

    jd_tt = Conversions.any_to_tt(jd_utc)
    jdc_tt = JulianDate.julian_day_datetime_to_century_datetime(jd_tt)

    pm_x = bd.f_pm_x(float(mjd_utc))
    pm_y = bd.f_pm_y(float(mjd_utc))
    pm_x = Conversions.arcsec_to_rad(pm_x)
    pm_y = Conversions.arcsec_to_rad(pm_y)
    dpm_x_dt = bd.f_pm_x(float(mjd_utc), derivative=True)
    dpm_x_dt = Conversions.arcsec_to_rad(dpm_x_dt)
    dpm_x_dt = Conversions.seconds_to_days(dpm_x_dt)
    dpm_y_dt = bd.f_pm_y(float(mjd_utc), derivative=True)
    dpm_y_dt = Conversions.arcsec_to_rad(dpm_y_dt)
    dpm_y_dt = Conversions.seconds_to_days(dpm_y_dt)

    sp = TransformationMatrices.calculate_s_prime(jdc_tt)
    dsp_dt = TransformationMatrices.calculate_s_prime_derivative()

    dt_ti_dt = TransformationMatrices.itrs_to_tirs_derivative(pm_x, pm_y, sp,
        dpm_x_dt, dpm_y_dt, dsp_dt)

    pm_x2 = bd.f_pm_x(float(mjd_utc) + dt)
    pm_y2 = bd.f_pm_y(float(mjd_utc) + dt)
    pm_x2 = Conversions.arcsec_to_rad(pm_x2)
    pm_y2 = Conversions.arcsec_to_rad(pm_y2)
    sp2 = TransformationMatrices.calculate_s_prime(jdc_tt + dtc)

    pm_x1 = bd.f_pm_x(float(mjd_utc) - dt)
    pm_y1 = bd.f_pm_y(float(mjd_utc) - dt)
    pm_x1 = Conversions.arcsec_to_rad(pm_x1)
    pm_y1 = Conversions.arcsec_to_rad(pm_y1)
    sp1 = TransformationMatrices.calculate_s_prime(jdc_tt - dtc)

    f2 = TransformationMatrices.itrs_to_tirs(pm_x2, pm_y2, sp2)
    f1 = TransformationMatrices.itrs_to_tirs(pm_x1, pm_y1, sp1)

    dt_ti_dt_fd = (f2 - f1) / (2 * dt)
    dt_ti_dt_fd *= Conversions.seconds_to_days(1.0)

    error = np.abs(dt_ti_dt - dt_ti_dt_fd)

    assert (np.max(error) < 1e-9)
