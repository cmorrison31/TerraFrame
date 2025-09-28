# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import re
from abc import ABC, abstractmethod

import numpy as np

from TerraFrame.PrecessionNutation import Arguments
from TerraFrame.Utilities import Conversions
from importlib import resources


class SeriesExpansion(ABC):
    def __init__(self, data_file_path):
        self.data_file_path = data_file_path
        self.data = []

        self._parse_file()

    def _parse_file(self):
        with open(self.data_file_path, 'r') as f:
            file_content = f.readlines()

            j = -1

            for line in file_content:
                # Search for a header line that set's the j value
                result = re.search(r'j\s*=\s*(\d+)\s+Number of terms',
                                   line.strip(), re.IGNORECASE)

                if result:
                    j = int(result.group(1))

                if j >= 0 and len(line.strip()) > 0:
                    try:
                        self.data.append([j, ] + [float(x)
                                                  for x in line.split()])
                    except ValueError:
                        pass

            self.data = np.array(self.data )

    @abstractmethod
    def compute(self, t):
        pass


class CipCoordinate(SeriesExpansion):
    def __init__(self, data_file_path, polynomial_coefficients):
        super().__init__(data_file_path)
        self._polynomial_coefficients = polynomial_coefficients

    def compute(self, t, derivative=False):
        """
        This function evaluates a celestial intermediate pole (CIP) series at
        the given Julian century date.

        :param t: Time as a Julian Century Date (jdc)
        :type t: JulianDate
        :param derivative: Compute the time derivative of the series
                           (default false)
        :type derivative: bool
        :return: The computed series value at time t and optionally the series
        time derivative
        :rtype: float | float, float
        """

        t = float(t)

        # units are micro-arcseconds
        poly_part = 0.0
        poly_part_dt = 0.0

        for j in range(len(self._polynomial_coefficients)):
            poly_part += self._polynomial_coefficients[j] * t ** j

        if derivative:
            for j in range(1, len(self._polynomial_coefficients)):
                poly_part_dt += self._polynomial_coefficients[j] * t ** (j - 1)

        # Initialize all the argument parameters. "argument" is the term that
        # IERS uses to refer to the input to the trigonometric functions. The
        # order is tightly coupled with the file format.
        arguments = np.zeros((14, ))
        d_arguments_dt = np.zeros((14, ))

        arguments[0] = Arguments.mean_anomaly_of_the_moon(t) # l
        arguments[1] = Arguments.mean_anomaly_of_the_sun(t) # l'
        arguments[2] = Arguments.mean_longitude_moon_minus_ascending_node(t) # F
        arguments[3] = Arguments.mean_elongation_of_the_moon_from_the_sun(t) # D
        arguments[4] = (
            Arguments.mean_longitude_of_the_ascending_node_of_the_moon(t)) # Ω
        arguments[5] = Arguments.mean_longitude_of_mercury(t) # L_Me
        arguments[6] = Arguments.mean_longitude_of_venus(t) # L_Ve
        arguments[7] = Arguments.mean_longitude_of_earth(t) # L_E
        arguments[8] = Arguments.mean_longitude_of_mars(t) # L_Ma
        arguments[9] = Arguments.mean_longitude_of_jupiter(t) # L_J
        arguments[10] = Arguments.mean_longitude_of_saturn(t) # L_Sa
        arguments[11] = Arguments.mean_longitude_of_uranus(t) # L_U
        arguments[12] = Arguments.mean_longitude_of_neptune(t) # L_Ne
        arguments[13] = Arguments.general_precession_in_longitude(t) # p_A

        if derivative:
            d_arguments_dt[0] = Arguments.mean_anomaly_of_the_moon_derivative(t)
            d_arguments_dt[1] = Arguments.mean_anomaly_of_the_sun_derivative(t)
            d_arguments_dt[2] = (
                Arguments
                .mean_longitude_moon_minus_ascending_node_derivative(t)
            )
            d_arguments_dt[3] = (
                Arguments
                .mean_elongation_of_the_moon_from_the_sun_derivative(t)
            )
            d_arguments_dt[4] = (
                Arguments
                .mean_longitude_of_the_ascending_node_of_the_moon_derivative(t)
            )
            d_arguments_dt[5] = Arguments.mean_longitude_of_mercury_derivative()
            d_arguments_dt[6] = Arguments.mean_longitude_of_venus_derivative()
            d_arguments_dt[7] = Arguments.mean_longitude_of_earth_derivative()
            d_arguments_dt[8] = Arguments.mean_longitude_of_mars_derivative()
            d_arguments_dt[9] = Arguments.mean_longitude_of_jupiter_derivative()
            d_arguments_dt[10] = Arguments.mean_longitude_of_saturn_derivative()
            d_arguments_dt[11] = Arguments.mean_longitude_of_uranus_derivative()
            d_arguments_dt[12] = (
                Arguments.mean_longitude_of_neptune_derivative())
            d_arguments_dt[13] = (
                Arguments
                .general_precession_in_longitude_derivative(t)
            )

        non_poly_part = 0.0
        d_non_poly_part_dt = 0.0

        for row in self.data:
            j = row[0]
            a_s = row[2]
            a_c = row[3]

            arg = np.linalg.vecdot(row[4:], arguments)

            arg_p = np.linalg.vecdot(row[4:], d_arguments_dt)

            tmp = a_s * np.sin(arg) + a_c * np.cos(arg)

            non_poly_part += tmp * t ** j

            if derivative:
                tmp2 = (a_s * np.cos(arg) - a_c * np.sin(arg)) * arg_p

                # Avoid a numerical issue with the generic formula when j = 0
                if j > 0:
                    d_non_poly_part_dt += j * t ** (j - 1) * tmp + t ** j * tmp2
                else:
                    d_non_poly_part_dt += tmp2


        total = poly_part + non_poly_part

        total_dt = poly_part_dt + d_non_poly_part_dt

        total = Conversions.muas_to_rad(total)
        total_dt = Conversions.muas_to_rad(total_dt)

        # Up till this point, we've been computing the derivative with
        # respect to centuries. Change it to be with respect to seconds.
        total_dt = Conversions.seconds_to_centuries(total_dt)

        if derivative:
            return total, total_dt
        else:
            return total


def cip_x(file_name=r'tab5.2a.txt'):
    file_path = resources.files("TerraFrame.Data").joinpath(file_name)
    return CipCoordinate(file_path, (-16617.0, 2004191898.0,
                                     -429782.9, -198618.34, 7.578, 5.9285))

def cip_y(file_name=r'tab5.2b.txt'):
    file_path = resources.files("TerraFrame.Data").joinpath(file_name)
    return CipCoordinate(file_path, (-6951.0, -25896.0,
                                     -22407274.7, 1900.59, 1112.526, 0.1358))

def cip_sxy2(file_name=r'tab5.2d.txt'):
    file_path = resources.files("TerraFrame.Data").joinpath(file_name)
    return CipCoordinate(file_path, (94.0, 3808.65,
                                     -122.68, -72574.11, 27.98, 15.62))
