# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from collections.abc import Callable, Iterable
from typing import Optional

import numpy as np
import numpy.typing as npt

from .DataFiles import DataFileManager
from .Interpolation import InterpolationPchip


class BulletinData:
    """
    This class acts as a container and interpolation helper for IERS Bulletin A
    correction data. All data is loaded only once and shared between class
    instances.

    The data is held in a numpy array. Each row is a different time
    (monotonically increasing). The layout in a row is as follows:
        0. Modified Julian Date (MJD UTC)
        1. UT1-UTC (seconds)
        2. Polar motion x-coordinate (arcseconds)
        3. Polar motion y-coordinate (arcseconds)
        4. Nutation correction x-coordinate (arcseconds)
        5. Nutation correction y-coordinate (arcseconds)

    Since the UTC deltas are provided in UTC but at single day resolutions, we
    treat the UTC delta data as if it's a function of UT1.

    """
    data: Optional[npt.NDArray[np.float64]]

    f_pm_x: Optional[
        Callable[[float | Iterable[float]], float | Iterable[float]]]
    f_pm_y: Optional[
        Callable[[float | Iterable[float]], float | Iterable[float]]]
    f_nc_dx: Optional[
        Callable[[float | Iterable[float]], float | Iterable[float]]]
    f_nc_dy: Optional[
        Callable[[float | Iterable[float]], float | Iterable[float]]]

    data = None
    f_pm_x = None
    f_pm_y = None
    f_nc_dx = None
    f_nc_dy = None

    def __init__(self):
        self._file_name = r'finals.all.iau2000.txt'
        self._url = (r'https://datacenter.iers.org/data/latestVersion/'
                     r'finals.all.iau2000.txt')

        seconds_in_a_month = 30 * 24 * 60 * 60

        self._datafile_manager = DataFileManager(self._file_name, self._url,
                                                 seconds_in_a_month)

        self._parse_file()

        self._init_interpolants()

    def __len__(self):
        return BulletinData.data[:, 0].shape[0]

    def dat_file_path(self):
        return self._datafile_manager.file_path

    @staticmethod
    def modified_julian_dates(index):
        if BulletinData.data is not None:
            return BulletinData.data[index, 0]
        else:
            raise RuntimeError('BulletinData must be initialized first.')

    @staticmethod
    def ut1_utc_delta(index):
        if BulletinData.data is not None:
            return BulletinData.data[index, 1]
        else:
            raise RuntimeError('BulletinData must be initialized first.')

    def _init_interpolants(self):
        if BulletinData.f_pm_x is None:
            BulletinData.f_pm_x = InterpolationPchip(self.data[:, 0],
                                                  self.data[:, 2])

            BulletinData.f_pm_y = InterpolationPchip(self.data[:, 0],
                                                  self.data[:, 3])

            BulletinData.f_nc_dx = InterpolationPchip(self.data[:, 0],
                                                   self.data[:, 4])

            BulletinData.f_nc_dy = InterpolationPchip(self.data[:, 0],
                                                   self.data[:, 5])

    def _parse_file(self):
        # Don't reparse the file data
        if BulletinData.data is not None:
            return

        file_content = []

        self._datafile_manager.get_file()

        with open(self._datafile_manager.file_path, 'r', encoding='utf8') as f:
            file_content = f.readlines()

        data_tmp = []

        for line in file_content:
            line_data = np.zeros((6,))

            if len(line.strip()) >= 125:
                # The file format is fixed width with a strict
                # specification. See readme.finals2000A.txt on the IERS
                # website.

                # Modified Julian Date (MJD UTC)
                line_data[0] = float(line[7:15])

                # Use Bulletin A values
                try:
                    # UT1-UTC (sec. of time)
                    line_data[1] = float(line[154:165])

                    # Polar motion, x (arcseconds)
                    line_data[2] = float(line[134:144])

                    # Polar motion, y (arcseconds)
                    line_data[3] = float(line[144:154])

                    # Nutation correction, dx (milliarcseconds)
                    line_data[4] = float(line[165:175])

                    # Nutation correction, dy (milliarcseconds)
                    line_data[5] = float(line[175:185])
                except ValueError:
                    # UT1-UTC (sec. of time)
                    line_data[1] = float(line[58:68])

                    # Polar motion, x (arcseconds)
                    line_data[2] = float(line[18:27])

                    # Polar motion, y (arcseconds)
                    line_data[3] = float(line[37:46])

                    # Nutation correction, dx (milliarcseconds)
                    line_data[4] = float(line[97:106])

                    # Nutation correction, dy (milliarcseconds)
                    line_data[5] = float(line[116:125])

                data_tmp.append(line_data)

            else:
                # Skip invalid lines
                continue

        BulletinData.data = np.array(data_tmp)
