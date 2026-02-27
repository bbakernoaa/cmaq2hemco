import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import hvplot.xarray
from cmaq2hemco.utils import regrid_dataset, _add_bounds_to_cmaq


def run_viz_example():
    """Example of AERO Protocol Visualization."""
    # 1. Generate dummy data
    nx, ny = 40, 40
    lon = np.linspace(-125, -70, nx)
    lat = np.linspace(25, 50, ny)
    LON, LAT = np.meshgrid(lon, lat)
    data = np.exp(-((LON + 95) ** 2 + (LAT - 38) ** 2) / 20)

    ds = xr.Dataset(
        {"emis": (("ROW", "COL"), data.astype("f"))},
        coords={"ROW": np.arange(ny), "COL": np.arange(nx)},
    )
    ds["lon"] = (("ROW", "COL"), LON, {"units": "degrees_east"})
    ds["lat"] = (("ROW", "COL"), LAT, {"units": "degrees_north"})
    # Instead of explicit 'crs' attribute, use IOAPI metadata
    ds.attrs["GDTYP"] = 2
    ds.attrs["P_ALP"] = 33.0
    ds.attrs["P_BET"] = 45.0
    ds.attrs["P_GAM"] = -97.0
    ds.attrs["XCENT"] = -97.0
    ds.attrs["YCENT"] = 40.0
    ds.attrs["XORIG"] = -1000000.0
    ds.attrs["YORIG"] = -1000000.0
    ds.attrs["XCELL"] = 50000.0
    ds.attrs["YCELL"] = 50000.0
    ds = _add_bounds_to_cmaq(ds)

    # 2. Regrid
    elat = np.linspace(25, 50, 21)
    elon = np.linspace(-125, -70, 21)
    clat = (elat[1:] + elat[:-1]) / 2
    clon = (elon[1:] + elon[:-1]) / 2
    target = xr.Dataset(
        {
            "lat": (["lat"], clat, {"units": "degrees_north"}),
            "lon": (["lon"], clon, {"units": "degrees_east"}),
        }
    )
    target["lat_b"] = (
        ["lat_v", "lon_v"],
        (elat[:, None] * np.ones((1, elon.size))).astype("f"),
    )
    target["lon_b"] = (
        ["lat_v", "lon_v"],
        (np.ones((elat.size, 1)) * elon[None, :]).astype("f"),
    )

    ds_regridded = regrid_dataset(ds, target, method="conservative")

    # TRACK A: Publication (Matplotlib + Cartopy)
    plt.figure(figsize=(10, 5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ds_regridded.emis.plot(
        ax=ax, transform=ccrs.PlateCarree(), x="lon", y="lat", cmap="viridis"
    )
    ax.coastlines()
    ax.set_title("Regridded Emissions (Conservative) - Track A")
    plt.savefig("examples/viz_track_a.png")
    print("Saved examples/viz_track_a.png")

    # TRACK B: Exploration (HvPlot)
    # Note: This is meant for interactive environments like Jupyter
    hv_plot = ds_regridded.emis.hvplot.quadmesh(
        x="lon",
        y="lat",
        geo=True,
        rasterize=True,
        cmap="viridis",
        title="Regridded Emissions (Conservative) - Track B",
    )
    hvplot.save(hv_plot, "examples/viz_track_b.html")
    print("Saved examples/viz_track_b.html")


if __name__ == "__main__":
    run_viz_example()
