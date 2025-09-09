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

    def __call__(self, xv):
        xv = ensure_iterable(xv)

        yv = len(xv) * [0.0]

        for i, v in enumerate(xv):
            index = self._get_index(v)

            # If we're out of bounds, return the first or last value
            if index == 0:
                yv[i] = self._y[index]
                continue
            elif index >= len(self._x):
                yv[i] = self._y[-1]
                continue

            y2 = self._y[index]
            y1 = self._y[index - 1]

            x2 = self._x[index]
            x1 = self._x[index - 1]

            yv[i] = (y2 - y1) / (x2 - x1) * (xv[i] - x1) + y1

        if len(yv) == 1:
            return yv[0]
        else:
            return yv
