import numpy as np
import xarray as xr
import pytest
import os
from cmaq2hemco.utils import regrid_dataset


def create_mock_ds(lazy=False):
    """Create a mock source dataset in LC projection."""
    nx, ny = 10, 10
    lon = np.linspace(-100, -90, nx)
    lat = np.linspace(30, 40, ny)
    LON, LAT = np.meshgrid(lon, lat)

    data = np.random.rand(1, ny, nx).astype("f")
    ds = xr.Dataset(
        {"emis": (("time", "ROW", "COL"), data)},
        coords={
            "time": [np.datetime64("2022-01-01")],
            "ROW": np.arange(ny),
            "COL": np.arange(nx),
        },
    )
    ds["lon"] = (("ROW", "COL"), LON, {"units": "degrees_east"})
    ds["lat"] = (("ROW", "COL"), LAT, {"units": "degrees_north"})
    ds.attrs["crs"] = (
        "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
    )
    ds.attrs["XCELL"] = 10000.0
    ds.attrs["YCELL"] = 10000.0

    if lazy:
        ds = ds.chunk({"ROW": 5, "COL": 5})
    return ds


def create_target_grid():
    """Create a target regular grid."""
    elat = np.linspace(31, 39, 5)
    elon = np.linspace(-99, -91, 5)
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
    return target


@pytest.mark.skipif(os.environ.get("SKIP_ESMPY") == "1", reason="esmpy not available")
def test_regrid_backend_agnostic():
    """Verify conservative re-gridding logic works identically for NumPy and Dask arrays."""
    from cmaq2hemco.utils import _add_bounds_to_cmaq

    ds_eager = _add_bounds_to_cmaq(create_mock_ds(lazy=False))
    ds_lazy = _add_bounds_to_cmaq(create_mock_ds(lazy=True))
    target = create_target_grid()

    # Run regridding
    out_eager = regrid_dataset(
        ds_eager, target, method="conservative", reuse_weights=False
    )
    out_lazy = regrid_dataset(
        ds_lazy, target, method="conservative", reuse_weights=False
    )

    # Assertions
    xr.testing.assert_allclose(out_eager, out_lazy.compute(), atol=1e-6)
    assert "history" in out_eager.attrs
    assert "Regridded using xregrid" in out_eager.attrs["history"]

    # Verify laziness
    assert out_lazy.emis.chunks is not None


@pytest.mark.skipif(os.environ.get("SKIP_ESMPY") == "1", reason="esmpy not available")
def test_weight_reuse():
    """Verify weight reuse functionality."""
    from cmaq2hemco.utils import _add_bounds_to_cmaq

    ds = _add_bounds_to_cmaq(create_mock_ds(lazy=False))
    target = create_target_grid()
    weights_file = "test_weights.nc"
    if os.path.exists(weights_file):
        os.remove(weights_file)

    try:
        # First pass: generate weights
        out1 = regrid_dataset(
            ds,
            target,
            method="conservative",
            reuse_weights=True,
            weights_file=weights_file,
        )
        assert os.path.exists(weights_file)

        # Second pass: reuse weights
        out2 = regrid_dataset(
            ds,
            target,
            method="conservative",
            reuse_weights=True,
            weights_file=weights_file,
        )

        xr.testing.assert_allclose(out1, out2)
    finally:
        if os.path.exists(weights_file):
            os.remove(weights_file)
