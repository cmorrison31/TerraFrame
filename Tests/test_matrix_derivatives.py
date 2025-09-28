# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np

from TerraFrame.Utilities import TransformationMatrices


def test_r1_derivative():
    n = 25
    t = np.linspace(0, 24 * 60 * 60, n)
    psi = [x ** 1.2 for x in t]
    dpsidt = 1.2
    dt = 1e-4
    dpsi = dpsidt * dt

    error = np.zeros((n,))

    for i in reversed(range(n)):
        r1 = TransformationMatrices.r1(psi[i] - dpsi)
        r2 = TransformationMatrices.r1(psi[i] + dpsi)

        drdt_e = (r2 - r1) / (2 * dt)
        drdt = TransformationMatrices.dr1dt(psi[i], dpsidt)

        error[i] = np.max(np.abs(drdt_e - drdt))

    assert (max(error) < 1e-6)


def test_r2_derivative():
    n = 25
    t = np.linspace(0, 24 * 60 * 60, n)
    psi = [x ** 1.2 for x in t]
    dpsidt = 1.2
    dt = 1e-4
    dpsi = dpsidt * dt

    error = np.zeros((n,))

    for i in reversed(range(n)):
        r1 = TransformationMatrices.r2(psi[i] - dpsi)
        r2 = TransformationMatrices.r2(psi[i] + dpsi)

        drdt_e = (r2 - r1) / (2 * dt)
        drdt = TransformationMatrices.dr2dt(psi[i], dpsidt)

        error[i] = np.max(np.abs(drdt_e - drdt))

    assert (max(error) < 1e-6)


def test_r3_derivative():
    n = 25
    t = np.linspace(0, 24 * 60 * 60, n)
    psi = [x ** 1.2 for x in t]
    dpsidt = 1.2
    dt = 1e-4
    dpsi = dpsidt * dt

    error = np.zeros((n,))

    for i in reversed(range(n)):
        r1 = TransformationMatrices.r3(psi[i] - dpsi)
        r2 = TransformationMatrices.r3(psi[i] + dpsi)

        drdt_e = (r2 - r1) / (2 * dt)
        drdt = TransformationMatrices.dr3dt(psi[i], dpsidt)

        error[i] = np.max(np.abs(drdt_e - drdt))

    assert (max(error) < 1e-6)
