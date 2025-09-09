# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

import TerraFrame
from TerraFrame.Utilities import TransformationMatrices
from TerraFrame.Utilities.Time import JulianDate


def get_earth_image_and_coordinates():
    img = plt.imread(r'world.topo.bathy.200412.3x5400x2700.jpg')

    # define a grid matching the map size, subsample along with pixels
    theta = np.linspace(0, np.pi, img.shape[0])
    phi = np.linspace(0, 2 * np.pi, img.shape[1])

    count = 300  # keep 180 points along theta and phi
    theta_inds = np.linspace(0, img.shape[0] - 1, count).round().astype(int)
    phi_inds = np.linspace(0, img.shape[1] - 1, count).round().astype(int)
    theta = theta[theta_inds]
    phi = phi[phi_inds]
    img = img[np.ix_(theta_inds, phi_inds)]

    theta, phi = np.meshgrid(theta, phi)

    return img, theta, phi


def main():
    axis_len = 6371.0e3
    n = 500  # Number of time steps
    days = 7
    times = np.linspace(0, 1 * days, n)
    earth_radius = 6371.0e3

    gcrs_basis = np.eye(3)
    itrs_basis_itrs = np.eye(3)

    itrs_vecs: list[np.ndarray] = []
    datetime_list = []
    t_list = []

    base_time = JulianDate.JulianDate.j2000()

    ct = TerraFrame.CelestialTerrestrialTransformation()

    for t in times:
        time = base_time + t

        datetime_list.append(JulianDate.pydatetime_from_julian_date(time))

        t_gi = ct.itrs_to_gcrs(time)

        angles = TransformationMatrices.euler_angles_from_transformation(t_gi)

        scale = 8000
        angles[0] *= 0
        angles[1] *= scale
        angles[2] *= scale

        t_gi = TransformationMatrices.transformation_from_euler(*angles)

        t_list.append(t_gi)

        itrs_basis_cgrs = t_gi @ itrs_basis_itrs  # 3x3, each column is one axis
        itrs_vecs.append(itrs_basis_cgrs)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection=Axes3D.name)
    ax.set_xlim(-axis_len, axis_len)
    ax.set_ylim(-axis_len, axis_len)
    # noinspection PyUnresolvedReferences
    ax.set_zlim(-axis_len, axis_len)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    # noinspection PyUnresolvedReferences
    ax.set_zlabel('Z')
    ax.view_init(elev=35.26, azim=20)
    gcrs_origin = np.zeros((3,))

    gcrs_basis *= earth_radius * 1.7

    for i in range(len(itrs_vecs)):
        itrs_vecs[i] *= earth_radius * 1.5

    # Add the Earth
    img, theta, phi = get_earth_image_and_coordinates()
    x = earth_radius * np.sin(theta) * np.cos(phi)
    y = earth_radius * np.sin(theta) * np.sin(phi)
    z = earth_radius * np.cos(theta)

    # The earth image needs to be rotated 180 degrees to align with (0, 0)
    # longitude and latitude
    tmp_a = -180 * np.pi / 180
    t_image = np.array(
        ((np.cos(tmp_a), np.sin(tmp_a), 0), (-np.sin(tmp_a), np.cos(tmp_a), 0),
         (0, 0, 1)), dtype='float')

    # GCRS: fixed basis vectors
    ax.quiver(gcrs_origin[0], gcrs_origin[1], gcrs_origin[2], gcrs_basis[0, 0],
              gcrs_basis[1, 0], gcrs_basis[2, 0], color='r', label='GCRS X'),
    ax.quiver(gcrs_origin[0], gcrs_origin[1], gcrs_origin[2], gcrs_basis[0, 1],
              gcrs_basis[1, 1], gcrs_basis[2, 1], color='g', label='GCRS Y'),
    ax.quiver(gcrs_origin[0], gcrs_origin[1], gcrs_origin[2], gcrs_basis[0, 2],
              gcrs_basis[1, 2], gcrs_basis[2, 2], color='b', label='GCRS Z')

    # ITRS: animated basis vectors
    itrs_quivers = [ax.quiver(gcrs_origin[0], gcrs_origin[1], gcrs_origin[2],
                              itrs_vecs[0][0, 0], itrs_vecs[0][1, 0],
                              itrs_vecs[0][2, 0], color='r', label='ITRS X',
                              linestyle='dashed'),
                    ax.quiver(gcrs_origin[0], gcrs_origin[1], gcrs_origin[2],
                              itrs_vecs[0][0, 1], itrs_vecs[0][1, 1],
                              itrs_vecs[0][2, 1], color='g', label='ITRS Y',
                              linestyle='dashed'),
                    ax.quiver(gcrs_origin[0], gcrs_origin[1], gcrs_origin[2],
                              itrs_vecs[0][0, 2], itrs_vecs[0][1, 2],
                              itrs_vecs[0][2, 2], color='b', label='ITRS Z',
                              linestyle='dashed')]

    ax.legend(loc='upper right')

    # noinspection PyUnresolvedReferences
    earths = [ax.plot_surface(x.T, y.T, z.T, facecolors=img / 255, cstride=1,
                              rstride=1, shade=False, edgecolor='none',
                              zorder=10), ]

    x_n = np.copy(x)
    y_n = np.copy(y)
    z_n = np.copy(z)

    # Animation Function
    def update(frame):
        earths[0].remove()

        # Rotate the earth mesh
        for xi in range(0, len(x)):
            for yj in range(0, len(x)):
                vec = np.array((x[xi, yj], y[xi, yj], z[xi, yj]))
                vec = np.matmul(np.transpose(np.copy(t_image)), vec)
                vec = np.matmul(t_list[frame], vec)

                x_n[xi, yj] = vec[0]
                y_n[xi, yj] = vec[1]
                z_n[xi, yj] = vec[2]

        for q in itrs_quivers:
            q.remove()
        itrs_quivers[0] = ax.quiver(gcrs_origin[0], gcrs_origin[1],
                                    gcrs_origin[2], itrs_vecs[frame][0, 0],
                                    itrs_vecs[frame][1, 0],
                                    itrs_vecs[frame][2, 0], color='r',
                                    label='ITRS X', linestyle='dashed')
        itrs_quivers[1] = ax.quiver(gcrs_origin[0], gcrs_origin[1],
                                    gcrs_origin[2], itrs_vecs[frame][0, 1],
                                    itrs_vecs[frame][1, 1],
                                    itrs_vecs[frame][2, 1], color='g',
                                    label='ITRS Y', linestyle='dashed')
        itrs_quivers[2] = ax.quiver(gcrs_origin[0], gcrs_origin[1],
                                    gcrs_origin[2], itrs_vecs[frame][0, 2],
                                    itrs_vecs[frame][1, 2],
                                    itrs_vecs[frame][2, 2], color='b',
                                    label='ITRS Z', linestyle='dashed')

        dt = datetime_list[frame].replace(microsecond=0)
        ax.set_title(f"Precession, Nutation, & Polar Motion over {days} "
                     f"Days\nScale {scale}x, Earth Rotation Suppressed"
                     f"\nUTC:"
                     f" {dt.isoformat()}", fontfamily="monospace", pad=0)

        # noinspection PyUnresolvedReferences
        earths[0] = ax.plot_surface(x_n.T, y_n.T, z_n.T, facecolors=img / 255,
                                    cstride=1, rstride=1, shade=False,
                                    edgecolor='none', zorder=10)

        return itrs_quivers

    ax_lim = (-1.5 * axis_len, 1.5 * axis_len)

    ax.set_xlim(ax_lim)
    ax.set_ylim(ax_lim)
    # noinspection PyUnresolvedReferences
    ax.set_zlim(ax_lim)
    ax.set_aspect('equal', 'box')

    # Create animation. The "ani" variable must exist even if it's not used.
    # noinspection PyUnusedLocal
    ani = FuncAnimation(fig, update, frames=len(times), interval=1/30*1e3)

    ani.save('animation.webm', fps=30)

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
