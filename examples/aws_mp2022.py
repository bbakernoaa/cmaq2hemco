import argparse
import pandas as pd
import numpy as np
import os
from joblib import Parallel, delayed
from cmaq2hemco.mp2022 import open_gdemis, open_ptemis
from cmaq2hemco.utils import pt2hemco, gd2hemco
from cmaq2hemco.mechrc.cb6r5_ae7_aq import writeconfig

debug = False
# Define grid by edges
elat = np.linspace(15, 65, 501)
elon = np.linspace(-135, -50, 851)
# excluding burn and biogenics
# openburn
# beis4
gkeys = """rwc
rail
onroad_gas
onroad_diesel
onroad_ca_adj_gas
onroad_ca_adj_diesel
np_solvents
np_oilgas
nonroad_gas
nonroad_diesel
nonpt
mexico_onroad
livestock
canada_ptdust_adj
canada_onroad
canada_og2D
canada_afdust_adj
airports
afdust_adj""".split()

# ignoring fires to avoid duplication
# ptagfire
# ptfire-wild
# ptfire-rx
# ptfire_othna
pkeys = """ptegu
canmex_point
cmv_c1c2_12
cmv_c3_12
pt_oilgas
ptnonipm""".split()


def process_gkey(date, gkey, elat, elon):
    outpath = f"epa2022v1/{gkey}/{gkey}_{date:%Y-%m-%d}_epa2022v1_hc_22m.nc"
    if os.path.exists(outpath):
        return
    print(date, gkey)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    try:
        # open_gdemis could fail if date is not available (eg., only weekday)
        gf = open_gdemis(date, gkey)
    except Exception as e:
        print(f"**WARNING:: Skipping {date} {gkey}: {e}")
        return
    # using xregrid conservative regridding
    # All gkeys share the same grid, so weights.nc can be reused.
    rgf = gd2hemco(
        outpath, gf, elat, elon, method="conservative", weights_file="weights.nc"
    )
    del rgf, gf


def process_pkey(date, pkey, elat, elon):
    outpath = f"epa2022v1/{pkey}/{pkey}_{date:%Y-%m-%d}_epa2022v1_hc_22m.nc"
    if os.path.exists(outpath):
        return
    print(date, pkey)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    try:
        # open_ptemis could fail if date is not available (eg., only weekday)
        pf = open_ptemis(date, pkey)
    except Exception as e:
        print(f"**WARNING:: Skipping {date} {pkey}: {e}")
        return
    rpf = pt2hemco(outpath, pf, elat, elon)  # apply plume rise
    del rpf, pf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process EPA 2022 Emissions to HEMCO format.")
    parser.add_argument("--start", type=str, help="Start date (YYYY-MM-DD)", default=None)
    parser.add_argument("--end", type=str, help="End date (YYYY-MM-DD)", default=None)
    parser.add_argument("--date", type=str, help="Single date to process (YYYY-MM-DD)", default=None)
    
    args = parser.parse_args()

    if args.date:
        dates = pd.date_range(args.date, args.date, freq="D")
    elif args.start and args.end:
        dates = pd.date_range(args.start, args.end, freq="D")
    elif args.start:
        dates = pd.date_range(args.start, args.start, freq="D")
    else:
        # Default to full year if no arguments provided
        dates = pd.date_range("2022-01-01", "2022-12-31", freq="D")

    # Pre-generate weights.nc to avoid race conditions in parallel
    if gkeys:
        first_date = dates[0]
        first_gkey = gkeys[0]
        print(f"Pre-generating weights using {first_gkey} on {first_date}")
        process_gkey(first_date, first_gkey, elat, elon)

    # Use Parallel to speed up processing (tune n_jobs as needed)
    for date in dates: 
        for gkey in gkeys:
            process_gkey(date,gkey, elat, elon)

    # Parallel(n_jobs=3)(
    #     delayed(process_pkey)(date, pkey, elat, elon)
    #     for date in dates
    #     for pkey in pkeys
    # )

    for sector in gkeys + pkeys:
        hcpath = f"epa2022v1/{sector}/HEMCO_{sector}.rc"
        secttmpl = f"epa2022v1/{sector}/{sector}_%Y-%m-%d_epa2022v1_hc_22m.nc"
        # Only write config if the directory exists (implying some files were processed)
        if os.path.exists(os.path.dirname(hcpath)):
            writeconfig(hcpath, 2022, sector, secttmpl)
