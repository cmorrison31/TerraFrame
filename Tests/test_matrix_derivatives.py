import numpy as np

from TerraFrame.Utilities import TransformationMatrices


def test_r3_derivative():
    n = 25
    t = np.linspace(0, 24 * 60 * 60, n)
    psi = [x ** 1.2 for x in t]
    dpsidt = 1.2
    dt = 1e-4
    dpsi = dpsidt * dt

    error = np.zeros((n,))

    for i in reversed(range(n)):
        r1 = TransformationMatrices.r3(psi[i])
        r2 = TransformationMatrices.r3(psi[i] + dpsi)

        drdt_e = (r2 - r1) / dt
        drdt = TransformationMatrices.dr3dt(psi[i], dpsidt)

        error[i] = np.max(np.abs(drdt_e - drdt))

    assert (max(error) < 1e-4)
