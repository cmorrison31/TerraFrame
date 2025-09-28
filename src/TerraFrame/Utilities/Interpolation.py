# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np

from .Helpers import clamp, ensure_iterable


class Interpolation1D:
    """
    This is not a general 1D interpolation class. It is specialized to the
    specific usage requirements and patterns of TerraFrame.

    Since most queries will be near each other, there is an index cache. If
    the index is out of bounds, the first or last values are used.
    """

    def __init__(self, x, y):
        self._x = x
        self._y = y

        self._index_cache = None

    def _get_index(self, xv):
        index = None

        # Under nominal usage patterns, most queries will use the same index
        # with only the occasional change
        if self._index_cache is not None:
            x1 = self._x[clamp(self._index_cache - 1, 0, len(self._x) - 1)]
            x2 = self._x[clamp(self._index_cache, 0, len(self._x) - 1)]

            if x1 < xv <= x2:
                index = self._index_cache

        if index is None:
            index = np.searchsorted(self._x, xv)
            self._index_cache = index

        return index

    def __call__(self, xv, derivative=False):
        xv = ensure_iterable(xv)

        yv = len(xv) * [0.0]

        for i, v in enumerate(xv):
            index = self._get_index(v)

            # If we're out of bounds, return the first or last value
            if index == 0:
                if derivative:
                    yv[i] = 0.0
                else:
                    yv[i] = self._y[index]
                continue
            elif index >= len(self._x):
                if derivative:
                    yv[i] = 0.0
                else:
                    yv[i] = self._y[-1]
                continue

            y2 = self._y[index]
            y1 = self._y[index - 1]

            x2 = self._x[index]
            x1 = self._x[index - 1]

            m = (y2 - y1) / (x2 - x1)

            if not derivative:
                yv[i] = m * (xv[i] - x1) + y1
            else:
                yv[i] = m

        if len(yv) == 1:
            return yv[0]
        else:
            return yv


class InterpolationPchip:
    """
    This is not a general PCHIP interpolation class. It is specialized to the
    specific usage requirements and patterns of TerraFrame.

    This Pchip routine is based on the method laid out by Butland et al. in [1].

    Since most queries will be near each other, there is an index cache. If
    the index is out of bounds, the first or last values are used.

    References:
        [1] F., N., And, F., J., & Butland (1984). A Method for Constructing
            Local Monotone Piecewise Cubic Interpolants. Siam Journal on
            Scientific and Statistical Computing, 5, 300-304.
            DOI: https://doi.org/10.1137/0905021
    """

    def __init__(self, x, y):
        self._x = x
        self._y = y

        self._index_cache = None
        self._coefficient_cache = None

    def _get_index(self, xv):
        index = None

        # Under nominal usage patterns, most queries will use the same index
        # with only the occasional change
        if self._index_cache is not None:
            x1 = self._x[clamp(self._index_cache - 1, 0, len(self._x) - 1)]
            x2 = self._x[clamp(self._index_cache, 0, len(self._x) - 1)]

            if x1 < xv <= x2:
                index = self._index_cache

        if index is None:
            # If the index cache isn't correct, we need to clear the
            # coefficient cache.
            self._coefficient_cache = None
            index = np.searchsorted(self._x, xv)
            self._index_cache = index

        return index

    def _get_gi(self, index):
        if index == 0 or index >= len(self._x) - 1:
            return 0

        h2 = self._x[index + 1] - self._x[index + 0]
        h1 = self._x[index + 0] - self._x[index - 1]

        s2 = (self._y[index + 1] - self._y[index + 0]) / h2
        s1 = (self._y[index + 0] - self._y[index - 1]) / h1

        if s1 * s2 > 0:
            alpha = 1.0 / 3.0 * (1 + h2 / (h1 + h2))

            gi = s1 * s2 / (alpha * s2 + (1 - alpha) * s1)
        else:
            gi = 0

        return gi

    def __call__(self, xv, derivative=False):
        xv = ensure_iterable(xv)

        yv = len(xv) * [0.0]

        for i, v in enumerate(xv):
            index = self._get_index(v)

            # If we're out of bounds, return the first or last value
            if index == 0:
                if derivative:
                    yv[i] = 0.0
                else:
                    yv[i] = self._y[index]
                continue
            elif index >= len(self._x):
                if derivative:
                    yv[i] = 0.0
                else:
                    yv[i] = self._y[-1]
                continue

            if self._coefficient_cache is None:
                g1 = self._get_gi(index - 1)
                g2 = self._get_gi(index + 0)

                y2 = self._y[index]
                y1 = self._y[index - 1]

                x2 = self._x[index]
                x1 = self._x[index - 1]

                a = np.array((
                (1, x1, x1 ** 2, x1 ** 3),
                (1, x2, x2 ** 2, x2 ** 3),
                (0, 1, 2 * x1, 3 * x1 ** 2),
                (0, 1, 2 * x2, 3 * x2 ** 2)
                ))

                b = np.array((
                (y1, ),
                (y2, ),
                (g1, ),
                (g2, )
                ))

                coefficients = np.linalg.solve(a, b)

                self._coefficient_cache = coefficients.flatten()

            if not derivative:
                yv[i] = (self._coefficient_cache[0] +
                         self._coefficient_cache[1] * xv[i] +
                         self._coefficient_cache[2] * xv[i] ** 2 +
                         self._coefficient_cache[3] * xv[i] ** 3)
            else:
                yv[i] = (self._coefficient_cache[1] +
                         self._coefficient_cache[2] * 2 * xv[i] +
                         self._coefficient_cache[3] * 3 * xv[i] ** 2)

        if len(yv) == 1:
            return yv[0]
        else:
            return yv
