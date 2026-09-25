#!/usr/bin/env python3
import os
import sys
import numpy as np

try:
    import xarray as xr
except ImportError:
    print("ERROR: Python package 'xarray' is not installed.", file=sys.stderr)
    sys.exit(1)

try:
    import netCDF4
except ImportError:
    print("ERROR: Python package 'netCDF4' is not installed.", file=sys.stderr)
    sys.exit(1)

EXP_METADATA = {
    "CTRL": {
        "basic_state": "p1",
        "region_1": "p1",
        "region_2": "p1",
        "region_3": "p1",
        "forcing_tag": "p1",
        "description": "Reference experiment: P1 basic state and P1 SST forcing over all regions"
    },
    "EXP_B": {
        "basic_state": "p2",
        "region_1": "p1",
        "region_2": "p1",
        "region_3": "p1",
        "forcing_tag": "p1",
        "description": "Basic-state effect: P2 basic state while retaining P1 SST-forcing"
    },
    "EXP_S-123": {
        "basic_state": "p1",
        "region_1": "p2",
        "region_2": "p2",
        "region_3": "p2",
        "forcing_tag": "p2",
        "description": "Total SST-forcing effect: P2 forcing over all three regions with P1 basic state"
    },
    "EXP_B_S-123": {
        "basic_state": "p2",
        "region_1": "p2",
        "region_2": "p2",
        "region_3": "p2",
        "forcing_tag": "p2",
        "description": "Full P2 configuration: P2 basic state and P2 SST forcing over all three regions"
    },
    "EXP_B_S-1": {
        "basic_state": "p2",
        "region_1": "p2",
        "region_2": "p1",
        "region_3": "p1",
        "forcing_tag": "p2_s1",
        "description": "P2 SST-forcing change in Region 1 under P2 atmospheric basic state"
    },
    "EXP_B_S-12": {
        "basic_state": "p2",
        "region_1": "p2",
        "region_2": "p2",
        "region_3": "p1",
        "forcing_tag": "p2_s12",
        "description": "Cumulative P2 SST-forcing changes in Regions 1 & 2 under P2 atmospheric basic state"
    }
}

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    tmp_nc_dir = os.path.join(base_dir, "_tmp_nc")
    output_nc_path = os.path.join(base_dir, "combined_ncep_r2_moist_lbm.nc")

    experiments = list(EXP_METADATA.keys())
    datasets = []
    first_ds_dims = None

    print(f"Opening and validating temporary NetCDF files for {len(experiments)} experiments...")
    for exp in experiments:
        nc_file = os.path.join(tmp_nc_dir, f"{exp}.nc")
        if not os.path.exists(nc_file):
            print(f"ERROR: NetCDF file does not exist: {nc_file}", file=sys.stderr)
            sys.exit(1)

        try:
            ds = xr.open_dataset(nc_file, decode_times=False)
        except Exception as e:
            print(f"ERROR: Cannot open NetCDF file {nc_file}: {e}", file=sys.stderr)
            sys.exit(1)

        current_dims = {k: v for k, v in ds.sizes.items()}
        if first_ds_dims is None:
            first_ds_dims = current_dims
        else:
            if current_dims != first_ds_dims:
                print(f"ERROR: Dimension mismatch in {exp}!", file=sys.stderr)
                print(f"Expected: {first_ds_dims}, got: {current_dims}", file=sys.stderr)
                sys.exit(1)

        meta = EXP_METADATA[exp]
        ds = ds.assign_coords(
            experiment=exp,
            basic_state=meta["basic_state"],
            region_1=meta["region_1"],
            region_2=meta["region_2"],
            region_3=meta["region_3"],
            forcing_tag=meta["forcing_tag"]
        )
        ds.attrs["description"] = meta["description"]
        ds = ds.expand_dims("experiment")
        datasets.append(ds)

    print(f"Concatenating {len(datasets)} datasets along 'experiment' dimension...")
    try:
        combined_ds = xr.concat(datasets, dim="experiment")
    except Exception as e:
        print(f"ERROR: Failed to concatenate datasets: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Saving combined dataset to {output_nc_path}...")
    try:
        combined_ds.to_netcdf(output_nc_path, format="NETCDF4")
    except Exception as e:
        print(f"ERROR: Failed to save final combined NetCDF: {e}", file=sys.stderr)
        sys.exit(1)

    print("SUCCESS: Combined NetCDF saved successfully.")

    # Validation
    try:
        ds_val = xr.open_dataset(output_nc_path, decode_times=False)
        print("\n--- Summary of Combined NetCDF ---")
        print(ds_val)
        print("\n--- Dimensions ---")
        for dim, size in ds_val.sizes.items():
            print(f"  {dim}: {size}")
        print("\n--- Experiments Present ---")
        print(ds_val["experiment"].values)
    except Exception as e:
        print(f"ERROR during validation: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
