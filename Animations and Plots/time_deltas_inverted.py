# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import matplotlib.pyplot as plt
import numpy as np

from TerraFrame.Utilities import Conversions
from TerraFrame.Utilities import Time
from TerraFrame.Utilities.Time import JulianDate


def main():
    tai_utc = Time.Deltas.TaiUtcDelta()

    jd_base = (JulianDate.
               julian_date_from_datetime(2015, 7, 1, 0, 0, 33,
                                         time_scale=JulianDate.TimeScales.TAI))

    offsets = [0, 1, 1.1, 2, 2.2, 4]
    n = len(offsets)
    deltas_tai_utc_t = np.zeros((n,))
    deltas_tai_utc_u = np.zeros((n,))

    jds_tai = n * [jd_base]
    jds_tai_dt = n * [jd_base]
    jds_uct_dt = n * [jd_base]
    jds_utc = n * [jd_base]
    jds_ut1 = n * [jd_base]

    for i, jd in enumerate(jds_tai):
        jds_tai[i] = jd + Conversions.seconds_to_days(offsets[i])
        jds_utc[i] = Conversions.tai_to_utc(jds_tai[i])
        jds_ut1[i] = Conversions.tai_to_ut1(jds_tai[i])

        deltas_tai_utc_t[i] = float(jds_tai[i] - jds_utc[i]) * 86400
        deltas_tai_utc_u[i] = tai_utc.get_delta(jds_utc[i])
        jds_tai_dt[i] = JulianDate.pydatetime_from_julian_date(jds_tai[i])
        jds_uct_dt[i] = JulianDate.pydatetime_from_julian_date(jds_utc[i])

    f, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.plot(jds_uct_dt, deltas_tai_utc_t, marker='^')
    ax1.plot(jds_uct_dt, deltas_tai_utc_u, marker='s')

    ax1.set_xlabel('Date UTC')
    ax1.set_ylabel('Time Delta')

    ax1.legend(['TAI - UTC Manual', 'TAI - UTC Calculated'])

    ax2.plot(jds_ut1, jds_utc, marker='^', linestyle='--')
    ax2.plot(jds_ut1, jds_ut1, marker='s', linestyle='--')

    ax2.set_xlabel('Date UT1')
    ax2.set_ylabel('Time')

    ax2.legend(['UTC', 'UT1'])

    ax3.plot(jds_utc, jds_utc, marker='^', linestyle='--')
    ax3.plot(jds_utc, jds_tai, marker='d', linestyle='--')
    ax3.plot(jds_utc, jds_ut1, marker='s', linestyle='--')

    ax3.set_xlabel('Datetime UTC')
    ax3.set_ylabel('Datetime')
    ax3.legend(['UTC', 'TAI', 'UT1'])

    plt.show()


if __name__ == '__main__':
    main()
