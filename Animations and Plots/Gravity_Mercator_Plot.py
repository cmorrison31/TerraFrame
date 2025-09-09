from TerraFrame.GravitationModel.GravitationalPotential import EGM2008
from TerraFrame import Earth
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def main():
    order = 5
    alt = 0
    egm = EGM2008(order)
    wgs84 = Earth.WGS84Ellipsoid()

    num_points = 200
    lats_deg = np.linspace(-90, 90, num_points)
    lons_deg = np.linspace(-180, 180, num_points)

    lats_rad = np.radians(lats_deg)
    lons_rad = np.radians(lons_deg)

    calcs = np.zeros((num_points, num_points))
    calcs_g = np.zeros((num_points, num_points))

    for i, lat_rad in enumerate(lats_rad):
        print(f'{i}/{num_points}')

        for j, lon_rad in enumerate(lons_rad):
            lat_geoc, lon_geoc, r = wgs84.geodetic_to_geocentric(lat_rad,
                                                                 lon_rad, alt)
            v, dvdx, dvdy, dvdz = egm.calculate(lat_geoc, lon_geoc, r)
            g = np.sqrt(dvdx ** 2 + dvdy ** 2 + dvdz ** 2)
            calcs_g[i, j] = g
            calcs[i, j] = v

    fig = plt.figure(figsize=(14, 6))

    # First subplot: Mercator projection with coastlines
    ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.Mercator())
    im = ax1.imshow(
        calcs_g,
        cmap="hot",
        interpolation="nearest",
        extent=[lons_deg[0], lons_deg[-1], lats_deg[0], lats_deg[-1]],
        transform=ccrs.PlateCarree(),
        origin="lower",
        aspect="auto"
    )
    ax1.coastlines(resolution="110m", linewidth=1)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax1.set_title("EGM2008 Gravitational Acceleration (Mercator projection)")
    fig.colorbar(im, ax=ax1, orientation="vertical", shrink=0.6,
                 label="Gravitational Acceleration")

    # Second subplot: line plot at longitude = 0
    ax2 = fig.add_subplot(1, 2, 2)
    lon0_idx = np.argmin(np.abs(lons_deg - 0))
    ax2.plot(lats_deg, calcs[:, lon0_idx])
    ax2.set_xlabel("Latitude (deg)")
    ax2.set_ylabel("Potential")
    ax2.set_title("Longitude = 0")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
