# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import datetime

from scipy.differentiate import derivative

import TerraFrame.Utilities.Conversions
from TerraFrame.PrecessionNutation import SeriesExpansion
from TerraFrame.Utilities import (Conversions, Time, BulletinData,
                                  TransformationMatrices)
from TerraFrame.Utilities.Time.JulianDate import JulianDate


class CelestialTerrestrialTransformation:
    def __init__(self, user_polar_motion=True, user_nutation_corrections=True):
        self.se_cip_x = SeriesExpansion.cip_x()
        self.se_cip_y = SeriesExpansion.cip_y()
        self.se_cip_sxy2 = SeriesExpansion.cip_sxy2()

        self._user_polar_motion = user_polar_motion
        self._user_nutation_corrections = user_nutation_corrections

        if not self._user_polar_motion and not self._user_nutation_corrections:
            self.bd = None
        else:
            self.bd = BulletinData.BulletinData()

        # Cached results
        self.t_gi = None
        self.t_gc = None
        self.t_ct = None
        self.t_ti = None

    def itrs_to_gcrs(self, time):
        if isinstance(time, datetime.datetime):
            time = Time.JulianDate.julian_date_from_pydatetime(time)

        assert isinstance(time, JulianDate)

        jd_tt = Conversions.any_to_tt(time)

        if time.time_scale == Time.TimeScales.UTC:
            jd_utc = time
        else:
            jd_utc = Conversions.tt_to_utc(jd_tt)

        # We also need time in Modified Julian Date (MJD) for the Bulletin
        # corrections lookup table.
        mjd_utc = Time.JulianDate.julian_date_to_modified_julian_date(jd_utc)

        # We also need time in UT1 for the ERA
        jd_ut1 = Conversions.tt_to_ut1(jd_tt)

        # Time needs to be in Julian centuries
        jdc_tt = Time.JulianDate.julian_day_datetime_to_century_datetime(jd_tt)

        # For the given terrestrial time (TT), call the routines to obtain the
        # IAU 2006/2000A X and Y from series. Then calculate "s" which is the
        # CIO locator
        cip_x = self.se_cip_x.compute(jdc_tt)
        cip_y = self.se_cip_y.compute(jdc_tt)
        sxy2 = self.se_cip_sxy2.compute(jdc_tt)
        cip_s = sxy2 - cip_x * cip_y / 2.0

        # Any CIP corrections ∆X, ∆Y can now be applied, and the corrected
        # X, Y, and s can be used to construct the Celestial Intermediate
        # Reference System (CIRS) to Geocentric Celestial Reference System
        # (GCRS) matrix: CIRS -> GCRS.
        # Get corrections by interpolating in the IERS Bulletin A data
        if self._user_nutation_corrections:
            dx = self.bd.f_nc_dx(float(mjd_utc))
            dy = self.bd.f_nc_dy(float(mjd_utc))

            cip_x += Conversions.mas_to_rad(dx)
            cip_y += Conversions.mas_to_rad(dy)

        # Create the first transformation matrix
        t_gc = TransformationMatrices.cirs_to_gcrs(cip_x, cip_y, cip_s)

        # The Earth rotation matrix is the transformation from the Terrestrial
        # Intermediate Reference System (TIRS) to the Celestial Intermediate
        # Reference System (CIRS): TIRS -> CIRS.
        # This function uses normal JD time in UT1.
        t_ct = TransformationMatrices.earth_rotation_matrix(jd_ut1)

        # Given polar motion offsets pm_x and pm_y, along with the Terrestrial
        # Intermediate Origin (TIO) locator (s prime or sp), the International
        # Terrestrial Reference System (ITRS) to Terrestrial Intermediate
        # Reference System (TIRS) transformation matrix can be constructed:
        # ITRS -> TIRS.
        if self._user_polar_motion:
            pm_x = self.bd.f_pm_x(float(mjd_utc))
            pm_y = self.bd.f_pm_y(float(mjd_utc))
        else:
            pm_x = 0.0
            pm_y = 0.0

        sp = TransformationMatrices.calculate_s_prime(jdc_tt)

        pm_x = TerraFrame.Utilities.Conversions.arcsec_to_rad(pm_x)
        pm_y = TerraFrame.Utilities.Conversions.arcsec_to_rad(pm_y)

        t_ti = TransformationMatrices.itrs_to_tirs(pm_x, pm_y, sp)

        # Construct the final transformation matrix: ITRS -> GCRS
        t_gi = t_gc @ t_ct @ t_ti

        self.t_gi = t_gi
        self.t_gc = t_gc
        self.t_ct = t_ct
        self.t_ti = t_ti

        return t_gi

    def gcrs_to_itrs(self, time):
        t_gi = self.itrs_to_gcrs(time)

        return t_gi.T


    def itrs_to_gcrs_angular_vel(self, time):
        if isinstance(time, datetime.datetime):
            time = Time.JulianDate.julian_date_from_pydatetime(time)

        assert isinstance(time, JulianDate)

        jd_tt = Conversions.any_to_tt(time)

        if time.time_scale == Time.TimeScales.UTC:
            jd_utc = time
        else:
            jd_utc = Conversions.tt_to_utc(jd_tt)

        # We also need time in Modified Julian Date (MJD) for the Bulletin
        # corrections lookup table.
        mjd_utc = Time.JulianDate.julian_date_to_modified_julian_date(jd_utc)

        # We also need time in UT1 for the ERA
        jd_ut1 = Conversions.tt_to_ut1(jd_tt)

        # Time needs to be in Julian centuries
        jdc_tt = Time.JulianDate.julian_day_datetime_to_century_datetime(jd_tt)

        # For the given terrestrial time (TT), call the routines to obtain the
        # IAU 2006/2000A X and Y from series. Then calculate "s" which is the
        # CIO locator
        (cip_x, dcip_x_dt) = self.se_cip_x.compute(jdc_tt, derivative=True)
        (cip_y, dcip_y_dt) = self.se_cip_y.compute(jdc_tt, derivative=True)
        (sxy2, dsxy2_dt) = self.se_cip_sxy2.compute(jdc_tt, derivative=True)
        cip_s = sxy2 - cip_x * cip_y / 2.0

        dcip_s_dt = dsxy2_dt - 0.5 * cip_y * dcip_x_dt - 0.5 * cip_x * dcip_y_dt

        # Any CIP corrections ∆X, ∆Y can now be applied, and the corrected
        # X, Y, and s can be used to construct the Celestial Intermediate
        # Reference System (CIRS) to Geocentric Celestial Reference System
        # (GCRS) matrix: CIRS -> GCRS.
        # Get corrections by interpolating in the IERS Bulletin A data
        if self._user_nutation_corrections:
            dx = self.bd.f_nc_dx(float(mjd_utc))
            dy = self.bd.f_nc_dy(float(mjd_utc))

            cip_x += Conversions.mas_to_rad(dx)
            cip_y += Conversions.mas_to_rad(dy)

            ddx_dt = self.bd.f_nc_dx(float(mjd_utc), derivative=True)
            ddy_dt = self.bd.f_nc_dy(float(mjd_utc), derivative=True)

            dcip_x_dt += ddx_dt
            dcip_y_dt += ddy_dt

        # Create the first transformation matrix
        t_gc = TransformationMatrices.cirs_to_gcrs(cip_x, cip_y, cip_s)
        dt_gc_dt = (TransformationMatrices
                    .cirs_to_gcrs_derivative(cip_x, cip_y, cip_s,
                                             dcip_x_dt, dcip_y_dt, dcip_s_dt))

        # The Earth rotation matrix is the transformation from the Terrestrial
        # Intermediate Reference System (TIRS) to the Celestial Intermediate
        # Reference System (CIRS): TIRS -> CIRS.
        # This function uses normal JD time in UT1.
        t_ct = TransformationMatrices.earth_rotation_matrix(jd_ut1)

        dt_ct_dt = (TransformationMatrices
                    .earth_rotation_matrix_derivative(jd_ut1))

        # Given polar motion offsets pm_x and pm_y, along with the Terrestrial
        # Intermediate Origin (TIO) locator (s prime or sp), the International
        # Terrestrial Reference System (ITRS) to Terrestrial Intermediate
        # Reference System (TIRS) transformation matrix can be constructed:
        # ITRS -> TIRS.
        if self._user_polar_motion:
            pm_x = self.bd.f_pm_x(float(mjd_utc))
            pm_y = self.bd.f_pm_y(float(mjd_utc))
            dpm_x_dt = self.bd.f_pm_x(float(mjd_utc), derivative=True)
            dpm_y_dt = self.bd.f_pm_y(float(mjd_utc), derivative=True)
        else:
            pm_x = 0.0
            pm_y = 0.0
            dpm_x_dt = 0.0
            dpm_y_dt = 0.0

        sp = TransformationMatrices.calculate_s_prime(jdc_tt)
        dsp_dt = TransformationMatrices.calculate_s_prime_derivative()

        pm_x = TerraFrame.Utilities.Conversions.arcsec_to_rad(pm_x)
        pm_y = TerraFrame.Utilities.Conversions.arcsec_to_rad(pm_y)

        t_ti = TransformationMatrices.itrs_to_tirs(pm_x, pm_y, sp)
        dt_ti_dt = TransformationMatrices.itrs_to_tirs_derivative(
            pm_x, pm_y, sp,
            dpm_x_dt, dpm_y_dt, dsp_dt)

        # Construct the final transformation matrix: ITRS -> GCRS
        t_gi = t_gc @ t_ct @ t_ti

        self.t_gi = t_gi
        self.t_gc = t_gc
        self.t_ct = t_ct
        self.t_ti = t_ti

        return t_gi
