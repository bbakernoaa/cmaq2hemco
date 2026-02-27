"""
Utilities for converting CMAQ-ready emissions to HEMCO supported files.
Follows the AERO Protocol 🍃⚡.
"""

from typing import Any, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr

# Optional dependencies
try:
    import xregrid
    from xregrid import Regridder
except ImportError:
    xregrid = None

xr.set_options(keep_attrs=True)

__all__ = [
    "plumerise_briggs",
    "open_date",
    "gd2hemco",
    "pt2hemco",
    "pt2gd",
    "merge",
    "to_ioapi",
    "getmw",
    "se_file",
    "gd2matrix",
    "gd2hemco_fast",
    "unitconvert",
    "hemco_area",
    "symlinks",
    "gd_file",
    "regrid_dataset",
]

# https://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids
# default vertical edge levels in meters above ground level.
_defez = np.array(
    [
        -6.00,
        123.0,
        254.0,
        387.0,
        521.0,
        657.0,
        795.0,
        934.0,
        1075.0,
        1218.0,
        1363.0,
        1510.0,
        1659.0,
        1860.0,
        2118.0,
        2382.0,
        2654.0,
        2932.0,
        3219.0,
        3665.0,
        4132.0,
        4623.0,
        5142.0,
        5692.0,
        6277.0,
        6905.0,
        7582.0,
        8320.0,
        9409.0,
        10504,
        11578,
        12633,
        13674,
        14706,
        15731,
        16753,
        17773,
        18807,
        19855,
        20920,
        22004,
        23108,
        24240,
        25402,
        26596,
        27824,
        29085,
        30382,
        31716,
        33101,
        34539,
        36030,
        37574,
        39173,
        40825,
        42529,
        44286,
        46092,
        47946,
        49844,
        51788,
        53773,
        55794,
        57846,
        59924,
        62021,
        64129,
        66245,
        68392,
        70657,
        73180,
        76357,
        80581,
    ],
    dtype="f",
)
# default vertical mid points in s = (p - ptop) / (p - p0) eta levels
_deflevs = np.array(
    [
        0.99250002413,
        0.97749990013,
        0.962499776,
        0.947499955,
        0.93250006,
        0.91749991,
        0.90249991,
        0.88749996,
        0.87249996,
        0.85750006,
        0.842500125,
        0.82750016,
        0.8100002,
        0.78750002,
        0.762499965,
        0.737500105,
        0.7125001,
        0.6875001,
        0.65625015,
        0.6187502,
        0.58125015,
        0.5437501,
        0.5062501,
        0.4687501,
        0.4312501,
        0.3937501,
        0.3562501,
        0.31279158,
        0.26647905,
        0.2265135325,
        0.192541016587707,
        0.163661504087706,
        0.139115,
        0.11825,
        0.10051436,
        0.085439015,
        0.07255786,
        0.06149566,
        0.05201591,
        0.04390966,
        0.03699271,
        0.03108891,
        0.02604911,
        0.021761005,
        0.01812435,
        0.01505025,
        0.01246015,
        0.010284921,
        0.008456392,
        0.0069183215,
        0.005631801,
        0.004561686,
        0.003676501,
        0.002948321,
        0.0023525905,
        0.00186788,
        0.00147565,
        0.001159975,
        0.00090728705,
        0.0007059566,
        0.0005462926,
        0.0004204236,
        0.0003217836,
        0.00024493755,
        0.000185422,
        0.000139599,
        0.00010452401,
        7.7672515e-05,
        5.679251e-05,
        4.0142505e-05,
        2.635e-05,
        1.5e-05,
    ],
    dtype="d",
)


_levattrs = dict(
    long_name="hybrid level at midpoints ((A/P0)+B)",
    units="level",
    axis="Z",
    positive="up",
    standard_name="atmosphere_hybrid_sigma_pressure_coordinate",
    formula_terms="a: hyam b: hybm p0: P0 ps: PS",
)
_latattrs = dict(
    long_name="Latitude",
    units="degrees_north",
    axis="Y",
    bounds="lat_bnds",
)
_lonattrs = dict(
    long_name="Longitude",
    units="degrees_east",
    axis="X",
    bounds="lon_bnds",
)
_reftime = "1970-01-01 00:00:00"
_timeattrs = dict(
    long_name="Time",
    units=f"hours since {_reftime}",
    calendar="gregorian",
    axis="T",
)


def getmw(key: str, gc: str = "cb6r5_ae7_aq", nr: str = "cb6r5hap_ae7_aq") -> float:
    """
    Get the molecular weight (kg/mol) for a chemical mechanism species.

    The species may be an explicit or lumped species, so weights are mechanism
    specific.

    Parameters
    ----------
    key : str
        Species key in mechanism.
    gc : str, default 'cb6r5_ae7_aq'
        Name of gas-phase chemical mechanism.
    nr : str, default 'cb6r5hap_ae7_aq'
        Name of non-reactive gas-phase mechanism (typically haps).

    Returns
    -------
    mw : float
        Molecular weight of species in kg/mol.

    Raises
    ------
    KeyError
        If species key is not found in the mechanism files.
    """
    import io
    import os
    import re

    import pandas as pd
    import requests

    mwpath = f"cmaq_{gc}_molwt.csv"
    fillin = {
        "CH4": 16.0,  # from ECH4
        "ETHYLBENZ": 106.2,  # from XYLMN
        "BENZ": 78.1,  # from BENZENE
        "NH3": 17.0,  # from NR_{mech}.nml
    }
    if not os.path.exists(mwpath):
        mwdfs = []
        for prfx, mech in [("GC", gc), ("NR", nr)]:
            url = (
                "https://raw.githubusercontent.com/USEPA/CMAQ/refs/heads/main/"
                f"CCTM/src/MECHS/{mech}/{prfx}_{mech}.nml"
            )
            r = requests.get(url, timeout=10)
            txtlines = r.text.split("\n")
            datlines = [
                _l for _l in txtlines if _l.startswith("'") or _l.startswith("!")
            ]
            datlines[0] = datlines[0].replace("!", "") + ","
            dat = "\n".join(datlines).replace("'", "")
            dat = re.sub("[ ]+,", ",", dat)
            rmwdf = pd.read_csv(io.StringIO(dat), index_col=False)
            rmwdf.columns = [k for k in rmwdf.columns]
            mwdfs.append(rmwdf.set_index("SPECIES"))
        mwdf = pd.concat(mwdfs)
        for newk, mw in fillin.items():
            if newk not in mwdf.index:
                mwdf.loc[newk, "MOLWT"] = mw
        mwdf[["MOLWT"]].to_csv(mwpath, index=True)
    mwdf = pd.read_csv(mwpath, index_col=0)
    try:
        mw = mwdf.loc[key, "MOLWT"] / 1e3
    except KeyError:
        raise KeyError(f"{key} not found in {mwpath}")
    return float(mw)


def plumerise_briggs(
    stkdm: Union[float, xr.DataArray],
    stkvel: Union[float, xr.DataArray],
    stktk: Union[float, xr.DataArray],
    pres_a: Union[float, xr.DataArray] = 101325.0,
    temp_a: Union[float, xr.DataArray] = 288.15,
    u: Union[float, xr.DataArray] = 2.5,
    x: Union[float, xr.DataArray] = 6000.0,
    theta_lapse: Optional[Union[float, xr.DataArray]] = None,
    F: Optional[Union[float, xr.DataArray]] = None,
) -> Union[float, xr.DataArray]:
    """
    Briggs (1969, 1971, 1974) equations of Plume Rise.

    Approximates calculations used by CMAQ and SMOKE to calculate plume rise.
    On pg 868, Table 18.4 gives 7 equations -- 3 for stable conditions and
    4 for neutral/unstable conditions.

    Parameters
    ----------
    stkdm : float or xr.DataArray
        Diameter of stack opening (m).
    stkvel : float or xr.DataArray
        Velocity of stack gas at opening (m/s).
    stktk : float or xr.DataArray
        Temperature of gas at opening (K).
    pres_a : float or xr.DataArray, default 101325.
        Pressure (Pa) of ambient environment.
    temp_a : float or xr.DataArray, default 288.15
        Temperature (K) of ambient environment.
    u : float or xr.DataArray, default 2.5
        Wind speed (m/s).
    x : float or xr.DataArray, default 6000.
        Distance from stack at which plume-rise is calculated (m).
    theta_lapse : float or xr.DataArray, optional
        Potential Temperature Gradient (dtheta / dz) in K/m.
    F : float or xr.DataArray, optional
        Buoyancy parameter (m4/s3).

    Returns
    -------
    dz : float or xr.DataArray
        Plume rise height of centerline.
    """
    g = 9.80665
    # Buoyancy flux parameter
    if F is None:
        dT = np.maximum(0, stktk - temp_a)
        F = g * stkdm**2 / 4 * stkvel * dT / stktk

    if theta_lapse is None or (
        isinstance(theta_lapse, (float, int)) and theta_lapse <= 0.005
    ):
        dz_neutral_loF_loX = (1.6 * F ** (1 / 3.0)) * x ** (2 / 3.0) / u
        dz_neutral_loF_hiX = (21.4 * F ** (3 / 4.0)) / u
        dz_neutral_hiF_hiX = (38.7 * F ** (3 / 5.0)) / u

        dz_loF = xr.where(
            x < (49 * F ** (5 / 8.0)), dz_neutral_loF_loX, dz_neutral_loF_hiX
        )
        dz_hiF = xr.where(
            x < (119 * F ** (2 / 5.0)), dz_neutral_loF_loX, dz_neutral_hiF_hiX
        )
        dz = xr.where(F < 55, dz_loF, dz_hiF)
    else:
        S2 = theta_lapse * g / temp_a
        dz_stable_1 = (2.4 * (F / S2)) ** (1 / 3.0) / u ** (1 / 3.0)
        dz_stable_2 = 5.0 * F ** (1 / 4.0) * S2 ** (-3 / 8.0)
        dz_stable_3 = 1.6 * F ** (1 / 3.0) * x ** (2 / 3.0) / u
        dz = np.minimum(dz_stable_1, dz_stable_2)
        dz = np.minimum(dz, dz_stable_3)

    return dz


def open_date(date: Any, tmpl: str, bucket: str, cache: bool = True) -> xr.Dataset:
    """
    Open all files for specific date from AWS S3.

    Parameters
    ----------
    date : str or datetime-like
        Date parsable by pandas.to_datetime.
    tmpl : str
        strftime template for date file (e.g., MCIP/12US1/GRIDCRO2D.12US1.35L.%y%m%d).
    bucket : str
        S3 Bucket name.
    cache : bool, default True
        Store file to prevent re-downloading.

    Returns
    -------
    ds : xarray.Dataset
        Opened dataset.

    Raises
    ------
    IOError
        If the file does not exist in the specified S3 bucket.
    """
    import gzip
    import io
    import os

    import boto3
    import cmaqsatproc as csp
    from botocore import UNSIGNED
    from botocore.client import Config

    client = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    date_dt = pd.to_datetime(date)
    path = date_dt.strftime(tmpl)
    dest = os.path.join(bucket, path)

    if cache:
        if not os.path.exists(dest):
            res = client.list_objects(Bucket=bucket, Prefix=path)
            if len(res.get("Contents", [])) != 1:
                raise IOError(f"{path} does not exist in {bucket}")
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            client.download_file(bucket, path, dest)

        if dest.endswith(".gz"):
            with gzip.open(dest, "rb") as gz:
                bdy = io.BytesIO(gz.read())
            f = csp.open_ioapi(bdy, engine="scipy")
        else:
            f = csp.open_ioapi(dest)
    else:
        res = client.list_objects(Bucket=bucket, Prefix=path)
        if len(res.get("Contents", [])) != 1:
            raise IOError(f"{path} does not exist in {bucket}")
        obj = client.get_object(Bucket=bucket, Key=path)
        bdy_data = obj["Body"].read()
        if path.endswith(".gz"):
            bdy_data = gzip.decompress(bdy_data)
        bdy = io.BytesIO(bdy_data)
        f = csp.open_ioapi(bdy, engine="scipy")

    return f


def pt2hemco(
    path: str,
    pf: xr.Dataset,
    elat: np.ndarray,
    elon: np.ndarray,
    ez: Optional[np.ndarray] = None,
    nk: int = 11,
    temp_a: float = 288.15,
    pres_a: float = 101325,
    u: float = 2.5,
    verbose: int = 0,
) -> "hemcofile":
    """
    Convert a point source file to a hemco-ready file.

    Allocates emissions to level based on stack height (STKHT) and Briggs
    plume rise.

    Parameters
    ----------
    path : str
        Path to save as HEMCO file.
    pf : xarray.Dataset
        Point source dataset from se_file.
    elat : np.ndarray
        Edge latitudes for regular grid.
    elon : np.ndarray
        Edge longitudes for regular grid.
    ez : np.ndarray, optional
        Edge altitudes in meters for vertical structure.
    nk : int, default 11
        Number of vertical levels to use if ez not specified.
    temp_a : float, default 288.15
        Ambient temperature (K) for plume rise.
    pres_a : float, default 101325
        Ambient pressure (Pa) for plume rise.
    u : float, default 2.5
        Wind speed (m/s) for plume rise.
    verbose : int, default 0
        Verbosity level.

    Returns
    -------
    outf : hemcofile
        HEMCO file object.
    """
    import warnings

    assert pf["lat"].min() >= elat.min()
    assert pf["lat"].max() <= elat.max()
    assert pf["lon"].min() >= elon.min()
    assert pf["lon"].max() <= elon.max()

    if ez is None:
        ez = _defez[: nk + 1]
    nk = ez.size - 1
    nr = elat.size - 1
    nc = elon.size - 1
    clat = (elat[1:] + elat[:-1]) / 2
    clon = (elon[1:] + elon[:-1]) / 2
    ilat = np.arange(nr)
    ilon = np.arange(nc)
    dx = pf.attrs.get("XCELL", 12000.0) / 2

    ris = pd.cut(pf["lat"], bins=elat, labels=ilat).astype("i")
    cis = pd.cut(pf["lon"], bins=elon, labels=ilon).astype("i")

    if "STKHT" not in pf or pf["STKHT"][:].isnull().all():
        warnings.warn("STKHT missing or all null; 2d output")
        nk = 1
        kis = np.zeros(pf.sizes["stack"], dtype="i")
    else:
        iz = np.arange(nk)
        dz = plumerise_briggs(
            pf["STKDM"],
            pf["STKVE"],
            pf["STKTK"],
            temp_a=temp_a,
            pres_a=pres_a,
            u=u,
            x=dx,
        )
        z = pf["STKHT"] + dz
        z = np.minimum(np.maximum(z, ez[0]), ez[-1])
        kis = pd.cut(z, bins=ez, labels=iz, include_lowest=True).astype("i")

    clev = _deflevs[:nk]
    tis = (
        ((pf.time - pf.time.min()).dt.total_seconds() / 3600).round(0).astype("i").data
    )
    pf["ti"] = ("time",), tis
    pf["ki"] = ("stack",), kis
    pf["ri"] = ("stack",), ris
    pf["ci"] = ("stack",), cis

    nt = pf.sizes["time"]
    datakeys = [
        k
        for k, v in pf.data_vars.items()
        if (
            k not in ("TFLAG", "lon", "lat", "ti", "ki", "ri", "ci") and len(v.dims) > 1
        )
    ]
    outf = hemcofile(
        path, pf.time, clat, clon, lev=clev, varkeys=datakeys, attrs=pf.attrs
    )
    area = hemco_area(elat, elon)
    outf.addvar("AREA", area, units="m2", dims=("lat", "lon"))

    for dk in datakeys:
        if verbose > 0:
            print(f"Processing {dk}")
        df = pf[["ti", "ki", "ri", "ci", dk]].to_dataframe()
        df = df.loc[df[dk] != 0]
        if df.empty:
            continue
        vals = df.groupby(["ti", "ki", "ri", "ci"], as_index=False).sum()
        tmp = np.zeros((nt, nk, nr, nc), dtype="f")
        tmp[vals.ti, vals.ki, vals.ri, vals.ci] = vals[dk].values
        attrs = {k: v for k, v in pf[dk].attrs.items()}
        unit = attrs.get("units", "unknown").strip()
        tmp, unit = unitconvert(dk, tmp, unit, area=area)
        attrs["units"] = unit
        outf.addvar(dk, tmp, **attrs)

    return outf


def regrid_dataset(
    ds: xr.Dataset,
    target_grid: xr.Dataset,
    method: str = "conservative",
    reuse_weights: bool = True,
    weights_file: str = "weights.nc",
    **kwargs: Any,
) -> xr.Dataset:
    """
    Regrid a Dataset using xregrid following AERO protocol.

    Parameters
    ----------
    ds : xarray.Dataset
        Source dataset.
    target_grid : xarray.Dataset
        Target grid dataset containing 'lat' and 'lon' (and boundaries for
        conservative).
    method : str, default 'conservative'
        Regridding method (bilinear, conservative, nearest_s2d, etc.).
    reuse_weights : bool, default True
        Whether to save and reuse weights from disk.
    weights_file : str, default 'weights.nc'
        Path to weights file.
    **kwargs : Any
        Additional arguments for xregrid.Regridder.

    Returns
    -------
    out_ds : xarray.Dataset
        Regridded dataset with updated scientific provenance in history.

    Raises
    ------
    ImportError
        If xregrid is not installed.
    """
    if xregrid is None:
        raise ImportError("xregrid is required for this operation.")

    import os

    # Handle weight reuse
    if reuse_weights and os.path.exists(weights_file):
        regridder = Regridder.from_weights(
            weights_file, ds, target_grid, method=method, **kwargs
        )
    else:
        regridder = Regridder(
            ds,
            target_grid,
            method=method,
            reuse_weights=reuse_weights,
            filename=weights_file,
            **kwargs,
        )

    out_ds = regridder(ds)

    # Update history for provenance
    history = out_ds.attrs.get("history", "")
    ts = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
    history += f"\n{ts}: Regridded using xregrid method={method}"
    out_ds.attrs["history"] = history

    return out_ds


def gd2hemco(
    path: str,
    gf: xr.Dataset,
    elat: np.ndarray,
    elon: np.ndarray,
    method: str = "conservative",
    reuse_weights: bool = True,
    weights_file: str = "weights.nc",
    verbose: int = 0,
) -> "hemcofile":
    """
    Regrid a CMAQ gridded file to HEMCO using xregrid (AERO Protocol).

    Parameters
    ----------
    path : str
        Output file path.
    gf : xarray.Dataset
        Source gridded emissions dataset.
    elat : np.ndarray
        Edge latitudes for target grid.
    elon : np.ndarray
        Edge longitudes for target grid.
    method : str, default 'conservative'
        Regridding method.
    reuse_weights : bool, default True
        Whether to reuse weights.
    weights_file : str, default 'weights.nc'
        File to save/load weights.
    verbose : int, default 0
        Verbosity level.

    Returns
    -------
    outf : hemcofile
        HEMCO file object.
    """
    if "lat" not in gf or "lon" not in gf:
        gf = gd_file(gf)

    # Create target grid
    clat = (elat[1:] + elat[:-1]) / 2
    clon = (elon[1:] + elon[:-1]) / 2
    target_grid = xr.Dataset(
        {
            "lat": (["lat"], clat, _latattrs),
            "lon": (["lon"], clon, _lonattrs),
        }
    )
    if method == "conservative":
        target_grid["lat_b"] = (
            ["lat_v", "lon_v"],
            (elat[:, None] * np.ones((1, elon.size))).astype("f"),
        )
        target_grid["lon_b"] = (
            ["lat_v", "lon_v"],
            (np.ones((elat.size, 1)) * elon[None, :]).astype("f"),
        )
        target_grid.lat.attrs["units"] = "degrees_north"
        target_grid.lon.attrs["units"] = "degrees_east"

        # Ensure source has bounds
        if "lat_b" not in gf:
            gf = _add_bounds_to_cmaq(gf)

    rgf = regrid_dataset(
        gf,
        target_grid,
        method=method,
        reuse_weights=reuse_weights,
        weights_file=weights_file,
    )

    clev = _deflevs[:1]
    datakeys = [
        k for k in rgf.data_vars if k not in ("TFLAG", "lon", "lat", "lat_b", "lon_b")
    ]
    outf = hemcofile(
        path, rgf.time, clat, clon, lev=clev, varkeys=datakeys, attrs=rgf.attrs
    )
    area = hemco_area(elat, elon)
    outf.addvar("AREA", area, units="m2", dims=("lat", "lon"))

    for dk in datakeys:
        if verbose > 0:
            print(f"Adding {dk}")
        attrs = {k: v for k, v in rgf[dk].attrs.items()}
        unit = attrs.get("units", "unknown").strip()
        # Unit conversion might be needed if source was not flux
        tmp, unit = unitconvert(dk, rgf[dk].expand_dims("lev", axis=1), unit, area=area)
        attrs["units"] = unit
        outf.addvar(dk, tmp, **attrs)

    return outf


def _add_bounds_to_cmaq(ds: xr.Dataset) -> xr.Dataset:
    """
    Add approximate bounds to a CMAQ dataset for xregrid.

    Parameters
    ----------
    ds : xarray.Dataset
        Input CMAQ dataset in projected coordinates.

    Returns
    -------
    ds : xarray.Dataset
        Dataset with 'lat_b' and 'lon_b' corner coordinates added.
    """
    import pyproj

    proj = pyproj.Proj(ds.crs)
    # Use coordinates directly if possible, else assume from attributes
    x = ds.COL
    y = ds.ROW
    dx = ds.attrs.get("XCELL", 12000.0)
    dy = ds.attrs.get("YCELL", 12000.0)

    # Edge coordinates in projected space
    xe = np.linspace(x.data[0] - dx / 2, x.data[-1] + dx / 2, len(x) + 1)
    ye = np.linspace(y.data[0] - dy / 2, y.data[-1] + dy / 2, len(y) + 1)

    XE, YE = np.meshgrid(xe, ye)
    LONE, LATE = proj(XE, YE, inverse=True)

    ds["lat_b"] = (("lat_v", "lon_v"), LATE.astype("f"))
    ds["lon_b"] = (("lat_v", "lon_v"), LONE.astype("f"))
    ds.lat.attrs["units"] = "degrees_north"
    ds.lon.attrs["units"] = "degrees_east"
    return ds


def gd2hemco_fast(
    path: str,
    gf: xr.Dataset,
    elat: np.ndarray,
    elon: np.ndarray,
    verbose: int = 0,
) -> "hemcofile":
    """
    Fast regridding using xregrid conservative method (AERO Protocol).

    Replaces old bilinear implementation for improved mass conservation.

    Parameters
    ----------
    path : str
        Output file path.
    gf : xarray.Dataset
        Source gridded emissions dataset.
    elat : np.ndarray
        Edge latitudes for target grid.
    elon : np.ndarray
        Edge longitudes for target grid.
    verbose : int, default 0
        Verbosity level.

    Returns
    -------
    outf : hemcofile
        HEMCO file object.
    """
    return gd2hemco(path, gf, elat, elon, method="conservative", verbose=verbose)


def gd2matrix(gf: xr.Dataset, elat: np.ndarray, elon: np.ndarray) -> pd.DataFrame:
    """
    Create a pixel to pixel fractional mapping matrix (Legacy).

    Used for fractional area overlap interpolation before xregrid integration.

    Parameters
    ----------
    gf : xarray.Dataset
        Input gridded dataset presenting the csp.geodf interface.
    elat : np.ndarray
        Edges of grid latitudes in degrees_north.
    elon : np.ndarray
        Edges of grid longitudes in degrees_east.

    Returns
    -------
    gdf : pd.DataFrame
        Mapping from ROW/COL to lon/lat cells with fraction of mass per
        ROW/COL assigned to lon/lat.
    """
    import warnings

    import geopandas as gpd
    from shapely.geometry import box

    clat = (elat[:-1] + elat[1:]) / 2
    clon = (elon[:-1] + elon[1:]) / 2
    hdx = np.diff(elon).mean() / 2
    hdy = np.diff(elat).mean() / 2
    qgeodf = gf.csp.geodf
    qgeodf["original_area"] = qgeodf.geometry.area
    if hdx == hdy:
        latj = np.arange(clat.size)
        loni = np.arange(clon.size)
        LONI, LATJ = np.meshgrid(loni, latj)
        LON, LAT = np.meshgrid(clon, clat)
        hgeodf = gpd.GeoDataFrame(
            dict(
                lat=LAT.ravel(), lon=LON.ravel(), lati=LATJ.ravel(), loni=LONI.ravel()
            ),
            geometry=gpd.points_from_xy(LON.ravel(), LAT.ravel()),
            crs=4326,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hgeodf["geometry"] = hgeodf["geometry"].buffer(hdx, cap_style="square")
    else:
        geoms = []
        lats, latis, lons, lonis = [], [], [], []
        for lati, lat in enumerate(clat):
            for loni, lon in enumerate(clon):
                lats.append(lat)
                latis.append(lati)
                lons.append(lon)
                lonis.append(loni)
                geoms.append(box(lon - hdx, lat - hdy, lon + hdx, lat + hdy))
        hgeodf = gpd.GeoDataFrame(
            dict(lat=lats, lon=lons, lati=latis, loni=lonis), geometry=geoms, crs=4326
        )
    ol = gpd.overlay(qgeodf.reset_index(), hgeodf.to_crs(qgeodf.crs))
    ol["intx_area"] = ol.geometry.area
    ol["fraction"] = ol["intx_area"] / ol["original_area"]
    ol["ri"] = ol["ROW"].astype("i")
    ol["ci"] = ol["COL"].astype("i")
    return ol.set_index(["ROW", "COL", "lati", "loni"])


def unitconvert(
    key: str,
    val: Union[np.ndarray, xr.DataArray],
    unit: str,
    area: Optional[Union[np.ndarray, xr.DataArray]] = None,
    inplace: bool = True,
) -> Tuple[Union[np.ndarray, xr.DataArray], str]:
    """
    Perform unit conversion to kg/m2/s.

    Parameters
    ----------
    key : str
        Species name to get molecular weight.
    val : np.ndarray or xr.DataArray
        Values in input units to be converted.
    unit : str
        Input unit string (e.g., 'moles/s', 'g/s/m2').
    area : np.ndarray or xr.DataArray, optional
        Area for each pixel (m2) if input is not per unit area.
    inplace : bool, default True
        If True, do the conversion in-place.

    Returns
    -------
    outval : np.ndarray or xr.DataArray
        Converted values in kg/m2/s.
    outunit : str
        Final unit string ('kg/m2/s').

    Raises
    ------
    AssertionError
        If '/s' is not in the input unit.
    """
    unit = unit.strip()
    outunit_parts = []
    factor = 1.0
    assert "/s" in unit
    inunit = unit.replace("/s", "")
    gps = ("g/s", "g/s/m**2", "g/m**2/s", "g/m2/s", "g/s/m2")
    nps = ("moles/s", "moles/s/m**2", "moles/m**2/s", "moles/s/m2", "moles/m2/s")

    if unit in gps:
        factor /= 1000.0
        outunit_parts.append("kg")
    elif unit in nps:
        try:
            mw = getmw(key)
            factor *= mw
            outunit_parts.append("kg")
        except KeyError as e:
            print(f"**WARNING: {key} in {unit} not converted to kg: {e}")
            outunit_parts.append(inunit)
    else:
        print(f"**WARNING: {key} [{unit}] not converted to kg: {unit} unknown")
        outunit_parts.append(inunit)

    if "m**2" not in unit and area is not None:
        factor /= area
        outunit_parts.append("/m**2")
    elif "/m**2" in unit:
        outunit_parts.append("/m**2")
    outunit_parts.append("/s")
    outunit = "".join(outunit_parts).replace("/s/m", "/m2/s")

    if inplace:
        val *= factor
        outval = val
    else:
        outval = val * factor

    return outval, outunit


def merge(fs: List[xr.Dataset], bf: Optional[xr.Dataset] = None) -> xr.Dataset:
    """
    Combine many datasets into a single dataset.

    Uses the mass from all datasets and the vertical structure of the
    tallest one.

    Parameters
    ----------
    fs : list of xarray.Dataset
        List of dataset objects to merge.
    bf : xarray.Dataset, optional
        Basis dataset for coordinates. If None, uses the tallest.

    Returns
    -------
    bf : xarray.Dataset
        Merged dataset with updated history.
    """
    import copy

    if bf is None:
        bf = sorted([(f.sizes.get("lev", 1), f) for f in fs])[-1][1][
            ["time", "lev", "lat", "lon"]
        ]
    fattrs = copy.deepcopy(bf.attrs)
    vattrs = {k: copy.deepcopy(v.attrs) for k, v in bf.data_vars.items()}
    for f in fs:
        for k in f.data_vars:
            if k not in bf:
                bf[k] = f[k]
                vattrs[k] = copy.deepcopy(f[k].attrs)
            else:
                bf[k] = bf[k] + f[k].data
            bf[k].attrs.update(vattrs[k])
    bf.attrs.update(fattrs)
    # Update history
    ts = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
    bf.attrs["history"] = bf.attrs.get("history", "") + f"\n{ts}: Merged datasets."
    return bf


def pt2gd(
    pf: xr.Dataset,
    nr: int,
    nc: int,
    ez: Optional[np.ndarray] = None,
    vgtyp: int = -9999,
    vgtop: float = 5000.0,
    vglvls: Optional[np.ndarray] = None,
    byvar: bool = True,
) -> xr.Dataset:
    """
    Convert point file to a CMAQ gridded emission file.

    Assigns vertical layers using Briggs plume rise calculations.

    Parameters
    ----------
    pf : xarray.Dataset
        Point source dataset with 'time' and 'stack' dimensions.
    nr : int
        Number of rows for the output grid.
    nc : int
        Number of columns for the output grid.
    ez : np.ndarray, optional
        Edges of the vertical array being created (m).
    vgtyp : int, default -9999
        Vertical grid type using IOAPI parameters.
    vgtop : float, default 5000.0
        Top of the vertical grid in Pascals.
    vglvls : np.ndarray, optional
        Edges of the vertical grid in terrain-following system.
    byvar : bool, default True
        Perform calculations by variable to reduce memory load.

    Returns
    -------
    outf : xarray.Dataset
        Gridded dataset representing point source mass.
    """
    ns = pf.sizes["stack"]
    nt = pf.sizes["time"]
    if ez is None:
        nz = 1
        kis = np.zeros((ns,), dtype="i")
    else:
        nz = ez.size - 1
        dz = plumerise_briggs(pf["STKDM"], pf["STKVE"], pf["STKTK"])
        dz = xr.where(np.isnan(dz), 0, dz)
        zs = pf["STKHT"] + dz
        zs = np.minimum(ez.max(), np.maximum(ez.min(), zs))
        kis = pd.cut(zs, bins=ez, labels=np.arange(nz, dtype="i")).astype("i")

    ex = np.arange(nc + 1) * pf.XCELL + pf.XORIG
    ey = np.arange(nr + 1) * pf.YCELL + pf.YORIG
    ris = pd.cut(pf["YLOCA"], bins=ey, labels=np.arange(nr, dtype="i"))
    cis = pd.cut(pf["XLOCA"], bins=ex, labels=np.arange(nc, dtype="i"))

    outf = xr.Dataset()
    outf.attrs.update(pf.attrs)
    outf.attrs["NCOLS"] = nc
    outf.attrs["NROWS"] = nr
    outf.attrs["NLAYS"] = nz
    if vglvls is None:
        vglvls = np.arange(nz + 1)
    outf.attrs["VGLVLS"] = vglvls
    outf.attrs["VGTYP"] = vgtyp
    outf.attrs["VGTOP"] = float(vgtop)

    datakeys = [k for k in pf.data_vars if k not in ("TFLAG",)]

    pf["ti"] = ("time",), pf.time.dt.hour.data
    pf["ki"] = ("stack",), kis
    pf["ri"] = ("stack",), ris
    pf["ci"] = ("stack",), cis

    for dk in datakeys:
        if len(pf[dk].dims) == 1:
            continue
        tmp = np.zeros((nt, nz, nr, nc), dtype="f")
        df = pf[["ti", "ki", "ri", "ci", dk]].to_dataframe()
        df = df.loc[df[dk] != 0].dropna()
        if not df.empty:
            vals = df.groupby(["ti", "ki", "ri", "ci"], as_index=False).sum()
            tmp[vals.ti, vals.ki, vals.ri, vals.ci] = vals[dk].values

        outf[dk] = (("TSTEP", "LAY", "ROW", "COL"), tmp, pf[dk].attrs)

    return outf


def se_file(sf: xr.Dataset, ef: xr.Dataset) -> xr.Dataset:
    """
    Combine stack group file and emission line file.

    Parameters
    ----------
    sf : xarray.Dataset
        CMAQ stack group dataset (dimensions include ROW).
    ef : xarray.Dataset
        CMAQ emission line dataset (dimensions include ROW, TSTEP).

    Returns
    -------
    ef : xarray.Dataset
        Combined dataset with dimensions ('time', 'stack').
    """
    ef = ef.rename(ROW="stack", TSTEP="time")
    sf = sf.isel(TSTEP=0, LAY=0, COL=0, drop=True).rename(ROW="stack")
    if "TFLAG" in sf:
        del sf["TFLAG"]
    if "TFLAG" in ef:
        del ef["TFLAG"]
    for k in sf.data_vars:
        if k not in ef:
            ef[k] = sf[k]
    ef["lat"] = sf["LATITUDE"]
    ef["lon"] = sf["LONGITUDE"]
    del ef["LATITUDE"]
    del ef["LONGITUDE"]
    return ef


def hemco_area(elat: np.ndarray, elon: np.ndarray, R: float = 6371007.2) -> np.ndarray:
    """
    Calculate area of grid cells on a sphere.

    Parameters
    ----------
    elat : np.ndarray
        Edges that define a regular grid in degrees north.
    elon : np.ndarray
        Edges that define a regular grid in degrees east.
    R : float, default 6371007.2
        Radius of the Earth (meters).

    Returns
    -------
    A : np.ndarray
        Area in square meters with dimensions (nlat, nlon).
    """
    dx = np.diff(elon).mean()
    a = dx * np.pi / 180 * R**2 * np.diff(np.sin(np.radians(elat)))
    A = a.astype("f")[:, None].repeat(elon.size - 1, 1)
    return A


class hemcofile:
    """HEMCO NetCDF file writer (Eager)."""

    def __init__(
        self,
        path: str,
        time: Any,
        lat: np.ndarray,
        lon: np.ndarray,
        lev: Optional[np.ndarray] = None,
        varkeys: Optional[List[str]] = None,
        attrs: Optional[dict] = None,
    ):
        """
        Initialize HEMCO NetCDF file.

        Parameters
        ----------
        path : str
            Path of the NetCDF file to be created.
        time : Any
            Times of the output file in UTC.
        lat : np.ndarray
            Latitudes for midpoints of grid (degrees_north).
        lon : np.ndarray
            Longitudes for midpoints of grid (degrees_east).
        lev : np.ndarray, optional
            Vertical layers for the destination file.
        varkeys : list of str, optional
            List of variables to pre-define.
        attrs : dict, optional
            Attributes for the output file.
        """
        import netCDF4

        self.nc = netCDF4.Dataset(path, mode="ws", format="NETCDF4_CLASSIC")
        if attrs:
            for pk, pv in attrs.items():
                self.nc.setncattr(pk, pv)
        self.nc.createDimension("time", None)
        if lev is not None:
            self.nc.createDimension("lev", lev.size)
        self.nc.createDimension("lat", lat.size)
        self.nc.createDimension("lon", lon.size)

        tv = self.nc.createVariable("time", "d", ("time",))
        for k, v in _timeattrs.items():
            tv.setncattr(k, v)
        timec = (
            (pd.to_datetime(time) - pd.to_datetime(_reftime)).total_seconds() / 3600.0
        ).astype("d")
        tv[: timec.size] = timec

        if lev is not None:
            levv = self.nc.createVariable("lev", "d", ("lev",))
            for k, v in _levattrs.items():
                levv.setncattr(k, v)
            levv[:] = lev

        latv = self.nc.createVariable("lat", "d", ("lat",))
        for k, v in _latattrs.items():
            latv.setncattr(k, v)
        latv[:] = lat

        lonv = self.nc.createVariable("lon", "d", ("lon",))
        for k, v in _lonattrs.items():
            lonv.setncattr(k, v)
        lonv[:] = lon

        if varkeys:
            for vk in varkeys:
                self.defvar(vk)

    def defvar(self, vk: str, dims: Optional[Tuple[str, ...]] = None, **attrs: Any):
        """
        Define a variable using HEMCO expectations.

        Parameters
        ----------
        vk : str
            Name of the output variable.
        dims : tuple of str, optional
            Named dimensions of the output variable.
        **attrs : Any
            Variable attributes.
        """
        nr = len(self.nc.dimensions["lat"])
        nc = len(self.nc.dimensions["lon"])
        chunkdefsizes = dict(time=1, lev=1, lat=nr, lon=nc)
        if dims is None:
            dims = (
                ("time", "lev", "lat", "lon")
                if "lev" in self.nc.dimensions
                else ("time", "lat", "lon")
            )
        chunks = [chunkdefsizes[dk] for dk in dims]
        vv = self.nc.createVariable(
            vk, "f", dims, chunksizes=chunks, zlib=True, complevel=1
        )
        vv.setncattr("standard_name", vk)
        vv.setncattr("units", "unknown")
        for pk, pv in attrs.items():
            vv.setncattr(pk, pv)

    def addvar(
        self,
        key: str,
        vals: np.ndarray,
        dims: Optional[Tuple[str, ...]] = None,
        **attrs: Any,
    ):
        """
        Add data to a variable, defining it if necessary.

        Parameters
        ----------
        key : str
            Name of the variable.
        vals : np.ndarray
            Values to add.
        dims : tuple of str, optional
            Named dimensions.
        **attrs : Any
            Variable attributes.
        """
        if key not in self.nc.variables:
            self.defvar(key, dims=dims, **attrs)
        vv = self.nc.variables[key]
        for pk, pv in attrs.items():
            vv.setncattr(pk, pv)
        nt = vals.shape[0]
        vv[:nt] = vals

    def close(self):
        """Close the NetCDF file."""
        if hasattr(self, "nc") and self.nc:
            self.nc.close()

    def __del__(self):
        self.close()


def to_ioapi(ef: xr.Dataset, path: str, **wopts: Any):
    """
    Write Dataset to a NetCDF file following IOAPI conventions.

    Parameters
    ----------
    ef : xarray.Dataset
        Emission dataset with time dimension to allow constructing TFLAG.
    path : str
        Output file path.
    **wopts : Any
        Write options for xarray.to_netcdf.
    """
    wopts.setdefault("mode", "ws")
    wopts.setdefault("format", "NETCDF4_CLASSIC")
    wopts.setdefault("unlimited_dims", ("TSTEP",))
    if "TFLAG" not in ef:
        nt = ef.sizes.get("TSTEP", 1)
        # Approximate time if missing
        time_vals = (
            ef.time
            if "time" in ef
            else pd.to_datetime([ef.attrs.get("SDATE", 2022001)] * nt)
        )
        date = time_vals.strftime("%Y%j").astype("i")
        time = time_vals.strftime("%H%M%S").astype("i")
        tflag = (
            xr.DataArray(
                np.array([date, time]).T,
                dims=("TSTEP", "DATE-TIME"),
            )
            .expand_dims(VAR=np.arange(len(ef.data_vars)))
            .transpose("TSTEP", "VAR", "DATE-TIME")
        )
        ef["TFLAG"] = tflag
    ef.to_netcdf(path, **wopts)


def symlinks(
    tmpl: str,
    dates: Union[pd.Series, str],
    datetype: Optional[str] = None,
    verbose: int = 0,
):
    """
    Create symlinks for date-based files.

    Parameters
    ----------
    tmpl : str
        strftime template for paths.
    dates : pd.Series or str
        Source and destination dates mapping.
    datetype : str, optional
        Key for selecting the date column from a CSV file.
    verbose : int, default 0
        Verbosity level.

    Returns
    -------
    links : list of str
        List of created symlink paths.
    """
    import os

    if isinstance(dates, str):
        df = pd.read_csv(dates)
        df.columns = [k.strip() for k in df.columns]
        df[datetype] = pd.to_datetime(df[datetype])
        df[df.columns[0]] = pd.to_datetime(df[df.columns[0]])
        dates_ser = df[[datetype, df.columns[0]]].set_index(df.columns[0])[datetype]
    else:
        dates_ser = dates

    links = []
    for outdate, indate in dates_ser.items():
        src = indate.strftime(tmpl)
        dst = outdate.strftime(tmpl)
        if not os.path.exists(dst):
            if not os.path.exists(src):
                continue
            os.symlink(src, dst)
            links.append(dst)
    return links


def gd_file(ef: xr.Dataset) -> xr.Dataset:
    """
    Add lon/lat coordinates and rename dimensions for gridded CMAQ datasets.

    Parameters
    ----------
    ef : xarray.Dataset
        Input dataset with 'crs' attribute and ROW/COL dimensions.

    Returns
    -------
    ef : xarray.Dataset
        Modified dataset with 'lon', 'lat', and 'time' meta-variables.
    """
    import pyproj

    if "LAY" in ef.dims:
        ef = ef.isel(LAY=0, drop=True)
    if "TSTEP" in ef.dims:
        ef = ef.rename(TSTEP="time")

    proj = pyproj.Proj(ef.crs)
    Y, X = xr.broadcast(ef.ROW, ef.COL)
    LON, LAT = proj(X.data, Y.data, inverse=True)
    ef["lon"] = (("ROW", "COL"), LON, dict(units="degrees_east", long_name="longitude"))
    ef["lat"] = (("ROW", "COL"), LAT, dict(units="degrees_north", long_name="latitude"))
    if "TFLAG" in ef:
        del ef["TFLAG"]
    return ef
