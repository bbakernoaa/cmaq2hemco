__all__ = ["open_gdemis", "open_ptemis"]

import xarray as xr

xr.set_options(keep_attrs=True)


def open_gdemis(date, sector, cache=True):
    """
    Open gridded (gd) emission file for a date/sector combination

    Arguments
    ---------
    date : datetime
        Datetime object that supports strftime interface
    sector : str
        Named premerged sector
    cache : bool
        If True, store the input file locally, if false read in memory.

    Returns
    -------
    ef : xarray.Dataset
        xarray Dataset with lon, lat, and time meta variables.
    """
    from .utils import find_s3_file, gd_file, open_date

    bucket = "epa-2022-modeling-platform"
    root = "emis/2022v2/"

    if sector == "merged_nobeis_norwc":
        search_pattern = f"{sector}/emis_mole_all"
    else:
        search_pattern = f"premerged_area/{sector}/emis_mole_{sector}"

    # Use fuzzy=True for sectors that may only have representative dates
    actual_key = find_s3_file(date, bucket, root, search_pattern, fuzzy=True)

    if not actual_key:
        # Fallback to older 2022v1
        root_v1 = "emis/2022v1/"
        actual_key = find_s3_file(date, bucket, root_v1, search_pattern, fuzzy=True)

    if not actual_key:
        raise IOError(f"Could not find S3 file for {sector} on {date}")

    ef = open_date(date, actual_key, bucket, cache=cache)
    ef = gd_file(ef)

    return ef


def open_ptemis(date, sector, cache=True):
    """
    Open point (pt) emission file for a date/sector combination.
    This interface adds metadata from the stack groups file to
    create a single self-describing data file for subsequent processing.


    Arguments
    ---------
    date : datetime
        Datetime object that supports strftime interface
    sector : str
        Named point sector
    cache : bool
        If True, store the input file locally, if false read in memory.

    Returns
    -------
    ef : xarray.Dataset
        xarray Dataset with lon, lat, and time meta variables.
    """
    from .utils import find_s3_file, open_date, se_file

    bucket = "epa-2022-modeling-platform"
    root = "emis/2022v2/"

    # Find stack groups file
    sg_key = find_s3_file(date, bucket, root, f"{sector}/stack_groups_{sector}", fuzzy=True)
    if not sg_key:
        # Fallback to v1
        sg_key = find_s3_file(
            date, bucket, "emis/2022v1/", f"{sector}/stack_groups_{sector}", fuzzy=True
        )
    # Search for emissions file (inline emissions are more likely to have daily files)
    epath = find_s3_file(date, bucket, root, f"{sector}/inln_mole_{sector}", fuzzy=True)
    if not epath:
        epath = find_s3_file(
            date, bucket, "emis/2022v1/", f"{sector}/inln_mole_{sector}", fuzzy=True
        )

    if not epath:
        raise IOError(f"Could not find emissions file for {sector} on {date}")

    ef = open_date(date, epath, bucket, cache=cache).isel(LAY=0, COL=0, drop=True)
    ef = se_file(sf, ef)
    return ef
