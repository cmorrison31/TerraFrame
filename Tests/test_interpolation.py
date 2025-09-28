# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np

from TerraFrame.Utilities import Interpolation


def test_linear_interpolation():
    x = np.array([1, 2, 3])
    y = np.array([1, 4, 9])

    f = Interpolation.Interpolation1D(x, y)

    xv = [0.5, 1.5, 2.0, 2.5, 3.0, 3.5]
    yv = f(xv)

    y_answer = [1.0, 2.5, 4.0, 6.5, 9.0, 9.0]

    for i, v in enumerate(yv):
        assert abs(v - y_answer[i]) < 1e-10


def test_pchip_interpolation():
    x = np.array([-1, 1, 2, 3, 4])
    y = np.array([1, 1, 4, 9, 16])

    f = Interpolation.InterpolationPchip(x, y)

    xv = np.linspace(-1.5, 4.5, 10)
    yv = f(xv)
    yvp = f(xv, derivative=True)

    # Answers calculated using Scipy's Pchip routine.
    y_answer = [1.0, 1.0, 1.0, 1.0, 1.1354166666666663, 3.3437499999999982,
                6.239583333333333, 10.193672839506169, 15.616512345679013, 16.0]
    yp_answer = [0.0, 0.0, 0.0, 0.0, 1.5624999999999987, 4.0625,
                 5.104166666666666, 8.263888888888884, 4.3750000000000036, 0.0]

    for i, (v, vp) in enumerate(zip(yv, yvp)):
        assert (abs(v - y_answer[i]) < 1e-10)
        assert (abs(vp - yp_answer[i]) < 1e-10)
