# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import math
from importlib import resources
from typing import Optional
import struct

import numpy as np
import numpy.typing as npt


class EGM2008:
    _data: Optional[npt.NDArray[np.float64]]

    # Data layout is Cbar(n, m), Sbar(n, m)
    # Coefficients are listed in triangular ordering (n increasing, m=0..n)
    _data = None

    def __init__(self, order=18):
        self.order = order

        # This is different from the WGS84 ellipsoid value and I think that's
        # intentional.
        self._a = 6378136.3  # meters.

        # Same comment about value difference.
        self._gm = 3.986004415e14  # m^3/s^2.

        self._load_data()

    def _load_data(self):
        """
        Load the EGM2008 coefficient data only if it hasn't been loaded
        already by another class instance.

        :return: None
        """

        if EGM2008._data is not None:
            return

        file_name = 'EGM2008_to200_TideFree.bin'
        file_path = resources.files("TerraFrame.Data").joinpath(file_name)

        row_format = struct.Struct("<HHdd")
        file_contents = []

        # noinspection PyTypeChecker
        with open(file_path, "rb") as f:
            while True:
                chunk = f.read(row_format.size)

                if not chunk:
                    break

                row = row_format.unpack(chunk)
                file_contents.append(row)

        EGM2008._data = np.array(file_contents)

        # If the requested order was zero, we read no data from the file
        if len(EGM2008._data) == 0:
            EGM2008._data = np.zeros((3, 4))

        # Check if the data starts and n = 2 and if so, add the implicit values
        elif EGM2008._data[0, 0] == 2:
            EGM2008._data = np.vstack((np.zeros((3, 4)), EGM2008._data))

        # Drop the n and m columns
        EGM2008._data = EGM2008._data[:, 2:]

        # Enforce C(0, 0) = 1 which isn't explicitly set by the EGM2008
        # coefficients but is expected by the Clenshaw summation
        EGM2008._data[0, 0] = 1

    @staticmethod
    def _get_index(n, m):
        """
        Given degree n and order m, this function returns the linear index
        into the coefficient array corresponding to the degree n and order m.

        :param n: degree
        :type n: int
        :param m: order
        :type m: int
        :return: index into EGM2008._data
        :rtype: int
        """

        index = n * (n + 1) // 2 + m

        return index

    def calculate(self, x, y, z):
        """
        This function calculates the gravitational potential and the three
        gradient terms in the ECEF frame. The input values are cartesian
        coordinates also in the ECEF frame.

        This function uses a double Clenshaw summation technique to calculate
        the spherical harmonic summations. This approach is fast, memory
        efficient, and numerically stable to the full degree and order of the
        EGM2008 model. All Legendre polynomial evaluations are avoided and
        just six trigonometric function calls are used.

        The approach used in this function was developed primarily from the
        equations and methodology outlined in [1]. Unless you're already very
        family with Clenshaw sums, this code will likely be unreadable without
        reading [1] first. The nature of the problem makes "variable name
        soup" unavoidable unfortunately.

        Sources:
        [1] Deakin, R.E., 1998, 'Derivatives of the earth’s potentials'.
            Geomatics Research Australasia, No.68, June, 1998, pp. 31-60.

        [2] Holmes, S., Featherstone, W. A unified approach to the Clenshaw
            summation and the recursive computation of very high degree and
            order normalised associated Legendre functions. Journal of
            Geodesy 76, 279–299 (2002).
            https://doi.org/10.1007/s00190-002-0216-2

        :param x: Cartesian x-coordinate in the ECEF frame
        :type x: float
        :param y: Cartesian y-coordinate in the ECEF frame
        :type y: float
        :param z: Cartesian z-coordinate in the ECEF frame
        :type z: float
        :return: Gravitational potential and the three gravitational
                 acceleration vector components in the ECEF frame.
        :rtype: float, float, float, float
        """

        # Note: all latitudes and longitudes in this function are geocentric,
        # not geodetic.

        longitude = math.atan2(y, x)
        r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
        q = self._a / r
        q2 = q ** 2
        t = z / r  # equivalent to sin(latitude)
        u = math.sqrt(x ** 2 + y ** 2) / r  # equivalent to cos(latitude)

        # Seed value for gravitational potential recurrence and it's partials
        v = 0
        dvdr = 0
        dvdl = 0  # l = longitude
        dvdt = 0  # t = sin(latitude)

        # Seed values for cosine and sine recurrence
        cl0 = math.cos(longitude)
        sl0 = math.sin(longitude)
        cl2 = math.cos((self.order + 2) * longitude)
        cl1 = math.cos((self.order + 1) * longitude)
        sl2 = math.sin((self.order + 2) * longitude)
        sl1 = math.sin((self.order + 1) * longitude)

        # Outer loop recurrence using Clenshaw's summation
        # We're tracking six summations:
        # 1. Gravitational potential, v
        # 2. Partial derivative of v w.r.t. spherical radius: dvdr
        # 3. Partial derivative of v w.r.t. longitude: dvdl
        # 4. Partial derivative of v w.r.t. t (which is sin(latitude): dvdt
        # 5. sin(m*longitude)
        # 6. cos(m*longitude)
        for m in range(self.order, -1, -1):
            # Inner loop recurrence seed values for gravitational potential
            # Note that dvdl re-uses the potential summation values
            sc2 = 0  # S, cosine m+2 term
            sc1 = 0  # S, cosine m+1 term
            sc = 0  # S, cosine m term
            ss2 = 0  # S, sine m+2 term
            ss1 = 0  # S, sine m+1 term
            ss = 0  # S, sine m term

            # Inner loop recurrence seed values for dvdr
            sc2_dr = 0  # S, cosine m+2 term
            sc1_dr = 0  # S, cosine m+1 term
            sc_dr = 0  # S, cosine m term
            ss2_dr = 0  # S, sine m+2 term
            ss1_dr = 0  # S, sine m+1 term
            ss_dr = 0  # S, sine m term

            # Inner loop recurrence seed values for dvdt
            sc_dt2 = 0  # S, cosine m+2 term
            sc_dt1 = 0  # S, cosine m+1 term
            sc_dt = 0  # S, cosine m term
            ss_dt2 = 0  # S, sine m+2 term
            ss_dt1 = 0  # S, sine m+1 term
            ss_dt = 0  # S, sine m term

            # Inner loop recurrence using Clenshaw's summation
            # We're tracking six summations:
            # 1. All v cosine related terms: c(n, m) * q^n * P(n, m)
            # 2. All v sine related terms: s(n, m) * q^n * P(n, m)
            # 3. The derivative w.r.t. t of all v cosine related terms
            # 4. The derivative w.r.t. t of all v sine related terms
            # 5. All dvdr cosine related terms
            # 6. All dvdr sine related terms
            #
            # Note that at no point is a legendre polynomial actually evaluated.
            # Isn't math crazy?
            for n in range(self.order, m - 1, -1):
                index = self._get_index(n, m)

                dadt = (math.sqrt((2 * n + 3) * (2 * n + 1) /
                                  ((n + m + 1) * (n - m + 1))) * q)
                a = dadt * t

                b = -q2 * math.sqrt((2 * n + 5) * (n + m + 1) * (n - m + 1) /
                                    ((2 * n + 1) * (n + m + 2) * (n - m + 2)))

                # dvdt update
                # Do this first since it uses sc1 and ss1
                sc_dt = a * sc_dt1 + b * sc_dt2 + dadt * sc1
                ss_dt = a * ss_dt1 + b * ss_dt2 + dadt * ss1

                sc_dt2 = sc_dt1
                sc_dt1 = sc_dt

                ss_dt2 = ss_dt1
                ss_dt1 = ss_dt

                # Gravitational potential update
                sc = a * sc1 + b * sc2 + EGM2008._data[index, 0]
                ss = a * ss1 + b * ss2 + EGM2008._data[index, 1]

                sc2 = sc1
                sc1 = sc

                ss2 = ss1
                ss1 = ss

                # dvdr update
                sc_dr = (a * sc1_dr + b * sc2_dr +
                         EGM2008._data[index, 0] * -1 * (n + 1))
                ss_dr = (a * ss1_dr + b * ss2_dr +
                         EGM2008._data[index, 1] * -1 * (n + 1))

                sc2_dr = sc1_dr
                sc1_dr = sc_dr

                ss2_dr = ss1_dr
                ss1_dr = ss_dr

            # Outer loop Clenshaw's summation for sine and cosine
            sl = 2 * cl0 * sl1 - sl2
            cl = 2 * cl0 * cl1 - cl2

            sl2 = sl1
            sl1 = sl

            cl2 = cl1
            cl1 = cl

            # We use a trick to get the fully normalized recurrence from [2] to
            # apply for all m.
            k = 2 if m > 0 else 1
            a = math.sqrt((2 * m + 3) / (k * (m + 1))) * u * q

            # Outer loop Clenshaw's summation for potential
            v = a * v + ss * sl + sc * cl

            # Outer loop Clenshaw's summation for dvdr
            dvdr = a * dvdr + ss_dr * sl + sc_dr * cl

            # Outer loop Clenshaw's summation for dvdl
            dvdl = a * dvdl + m * (ss * cl - sc * sl)

            # Outer loop Clenshaw's summation for dvdl
            dvdt = (a * dvdt + (sc_dt * cl + ss_dt * sl) -
                    (sc * cl + ss * sl) * t / u ** 2 * m)

        v = self._gm / r * v
        dvdr = self._gm / r ** 2 * dvdr
        dvdl = self._gm / r * dvdl
        dvdt = self._gm / r * dvdt

        # Convert from spherical coordinates to cartesian coordinates
        dvdx = u * cl0 * dvdr - u * t * cl0 / r * dvdt - sl0 / (u * r) * dvdl
        dvdy = u * sl0 * dvdr - u * t * sl0 / r * dvdt + cl0 / (u * r) * dvdl
        dvdz = t * dvdr + u ** 2 / r * dvdt

        return v, dvdx, dvdy, dvdz
