"""
process_wind.py
---------------
Functions for loading ERA5 NetCDF wind data, selecting a target grid point,
and deriving wind speed and direction from u/v components.

All logic here was validated in notebooks/01_era5_exploration.ipynb.
"""

import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

def load_era5_year(year: int, data_dir: str, lat: float, lon: float) -> pd.DataFrame:
    filepath = Path(data_dir) / f"era5_sydney_{year}.nc"

    if not filepath.exists():
        raise FileNotFoundError(
            f"ERA5 file not found: {filepath}\n"
            f"Run fetch_era5.py to download missing data."
        )
    
    ds = xr.open_dataset(filepath)
    ds_point = ds.sel(latitude=lat, longitude=lon, method="nearest")

    sel_lat = float(ds_point.latitude)
    sel_lon = float(ds_point.longitude)

    # Extract arrays — convert to numpy immediately, release xarray
    u    = ds_point.u100.values.astype(np.float64)
    v    = ds_point.v100.values.astype(np.float64)
    time = ds_point.valid_time.values
    ds.close()

    wind_speed = np.sqrt(u**2 + v**2)
    wind_direction = (270 - np.degrees(np.arctan2(v, u))) % 360

    df = pd.DataFrame({
        "u100":           u,
        "v100":           v,
        "wind_speed":     wind_speed,
        "wind_direction": wind_direction
    }, index=pd.to_datetime(time))
    df.index.name = "time"

    return df

def load_era5_multiyear(years: list, data_dir: str,
                         lat: float, lon: float) -> pd.DataFrame:
    frames = []
    for year in years:
        df_year = load_era5_year(year, data_dir, lat, lon)
        frames.append(df_year)
        print(f"  Loaded {year}: {len(df_year)} timesteps | "
              f"Mean speed: {df_year['wind_speed'].mean():.2f} m/s")

    df_all = pd.concat(frames).sort_index()

    print(f"\n  Combined: {len(df_all)} total timesteps")
    print(f"  Period:   {df_all.index[0]} to {df_all.index[-1]}")
    print(f"  Overall mean speed: {df_all['wind_speed'].mean():.2f} m/s")

    return df_all

def add_temporal_features(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()  # avoid mutating the input

    df["month"] = df.index.month
    df["season"] = df.index.month.map({
        12: "Summer", 1: "Summer",  2: "Summer",
         3: "Autumn", 4: "Autumn",  5: "Autumn",
         6: "Winter", 7: "Winter",  8: "Winter",
         9: "Spring", 10: "Spring", 11: "Spring"
    })

    return df

def basic_diagnostics(df: pd.DataFrame) -> dict:
    speeds = df["wind_speed"]

    diagnostics = {
        "n_hours":          len(speeds),
        "mean_ms":          speeds.mean(),
        "median_ms":        speeds.median(),
        "std_ms":           speeds.std(),
        "max_ms":           speeds.max(),
        "min_ms":           speeds.min(),
        "pct_below_cutin":  (speeds < 3.0).mean()  * 100,
        "pct_above_rated":  (speeds > 12.0).mean() * 100,
        "pct_above_cutout": (speeds > 25.0).mean() * 100,
        "missing_values":   speeds.isnull().sum()
    }

    print("=== WIND SPEED DIAGNOSTICS ===")
    print(f"  Hours:              {diagnostics['n_hours']:,}")
    print(f"  Mean:               {diagnostics['mean_ms']:.2f} m/s")
    print(f"  Median:             {diagnostics['median_ms']:.2f} m/s")
    print(f"  Std deviation:      {diagnostics['std_ms']:.2f} m/s")
    print(f"  Max:                {diagnostics['max_ms']:.2f} m/s")
    print(f"  % below cut-in:     {diagnostics['pct_below_cutin']:.1f}%")
    print(f"  % above rated:      {diagnostics['pct_above_rated']:.1f}%")
    print(f"  % above cut-out:    {diagnostics['pct_above_cutout']:.1f}%")
    print(f"  Missing values:     {diagnostics['missing_values']}")

    return diagnostics

