"""
wind_rose.py
------------
Functions for directional wind analysis including sector frequency
statistics, power rose calculation, and wind rose visualisation.

All logic validated in notebooks/01_era5_exploration.ipynb.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from windrose import WindroseAxes
from pathlib import Path

# Standard 16-sector compass labels
SECTOR_LABELS = [
    "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
    "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"
]

SEASON_COLOURS = {
    "Summer": "#FF6B6B",
    "Autumn": "#FFA07A",
    "Winter": "#4169E1",
    "Spring": "#32CD32"
}

def assign_sectors(wind_direction: np.ndarray,
                   n_sectors: int = 16) -> np.ndarray:
    sector_width = 360.0 / n_sectors

    # Shift by half sector width so sector boundaries fall between labels
    # e.g. N sector covers 348.75°–11.25° not 0°–22.5°
    shifted = (wind_direction + sector_width / 2) % 360
    indices = (shifted / sector_width).astype(int) % n_sectors

    return indices

def directional_statistics(df: pd.DataFrame,
                            n_sectors: int = 16) -> pd.DataFrame:
    sector_width   = 360.0 / n_sectors
    sector_centres = np.arange(0, 360, sector_width)
    sector_indices = assign_sectors(df["wind_direction"].values, n_sectors)

    rows = []
    for i, (label, centre) in enumerate(zip(SECTOR_LABELS[:n_sectors],
                                             sector_centres)):
        mask  = sector_indices == i
        count = mask.sum()
        freq  = count / len(df) * 100

        if count > 0:
            speeds       = df["wind_speed"][mask]
            mean_spd     = speeds.mean()
            max_spd      = speeds.max()
            pct_rated    = (speeds > 12.0).mean() * 100
            mean_u_cubed = (speeds**3).mean()
            energy_index = (freq / 100) * mean_u_cubed
        else:
            mean_spd = max_spd = pct_rated = mean_u_cubed = energy_index = 0.0

        rows.append({
            "label":            label,
            "centre_deg":       centre,
            "frequency_pct":    freq,
            "mean_speed":       mean_spd,
            "max_speed":        max_spd,
            "pct_above_rated":  pct_rated,
            "energy_index":     energy_index
        })

    df_sectors = pd.DataFrame(rows)

    # Normalise energy index to percentage of total
    total_energy = df_sectors["energy_index"].sum()
    df_sectors["energy_pct"] = (
        df_sectors["energy_index"] / total_energy * 100
        if total_energy > 0 else 0.0
    )

    return df_sectors

def print_directional_table(df_sectors: pd.DataFrame) -> None:
    print("=== DIRECTIONAL ANALYSIS ===")
    print(f"{'Sector':<6} {'Centre':>7} {'Freq':>8} {'Mean U':>8} "
          f"{'Max U':>7} {'> Rated':>9} {'Energy':>8}")
    print("-" * 62)

    for _, row in df_sectors.iterrows():
        bar = "█" * int(row["frequency_pct"] / 1.5)
        print(f"{row['label']:<6} {row['centre_deg']:>6.1f}°  "
              f"{row['frequency_pct']:>7.1f}%  "
              f"{row['mean_speed']:>7.2f}   "
              f"{row['max_speed']:>6.2f}  "
              f"{row['pct_above_rated']:>8.1f}%  "
              f"{row['energy_pct']:>6.1f}%  {bar}")
        
    # Key findings
    dominant_freq   = df_sectors.loc[df_sectors["frequency_pct"].idxmax()]
    dominant_energy = df_sectors.loc[df_sectors["energy_pct"].idxmax()]

    print(f"\n  Most frequent:      {dominant_freq['label']} "
          f"({dominant_freq['centre_deg']:.1f}°) - "
          f"{dominant_freq['frequency_pct']:.1f}% of hours")
    print(f"  Most energetic:     {dominant_energy['label']} "
          f"({dominant_energy['centre_deg']:.1f}°) - "
          f"{dominant_energy['energy_pct']:.1f}% of energy flux")

    if dominant_freq["label"] != dominant_energy["label"]:
        print(f"\n  Note: Most frequent and most energetic directions differ.")
        print(f"  The cubic U³ weighting shifts energy dominance toward "
              f"higher-speed sectors.")
        
def plot_wind_rose(df: pd.DataFrame,
                   title: str = "",
                   save_path: str = None) -> None:
    
    fig = plt.figure(figsize=(9, 9))
    ax  = WindroseAxes.from_ax(fig=fig)

    ax.bar(df["wind_direction"],
           df["wind_speed"],
           normed=True,
           opening=0.9,
           edgecolor="white",
           linewidth=0.5,
           bins=[0, 3, 7, 12, 18, 25],
           nsector=16)

    ax.set_legend(title="Wind Speed (m/s)",
                  loc="lower right",
                  fontsize=9)
    
    full_title = title if title else "Wind Rose"
    ax.set_title(full_title, fontsize=12, pad=20)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved: {save_path}")

    plt.show()

def plot_seasonal_roses(df: pd.DataFrame,
                        title: str = "",
                        save_path: str = None) -> None:
    if "season" not in df.columns:
        raise ValueError(
            "DataFrame missing 'season' column. "
            "Run add_temporal_features() first."
        )

    seasons   = ["Summer", "Autumn", "Winter", "Spring"]
    positions = [(0, 0), (0, 1), (1, 0), (1, 1)]

    fig, axes = plt.subplots(2, 2, figsize=(14, 14),
                              subplot_kw=dict(projection="windrose"))

    for season, (row, col) in zip(seasons, positions):
        mask        = df["season"] == season
        season_data = df[mask]
        n_hours     = mask.sum()
        mean_spd    = season_data["wind_speed"].mean()

        ax = WindroseAxes.from_ax(ax=axes[row, col])
        ax.bar(season_data["wind_direction"],
               season_data["wind_speed"],
               normed=True,
               opening=0.9,
               edgecolor="white",
               linewidth=0.5,
               bins=[0, 3, 7, 12, 18, 25],
               nsector=16)

        ax.set_title(f"{season}\n"
                     f"n={n_hours:,} hrs | "
                     f"Mean: {mean_spd:.2f} m/s",
                     fontsize=11, pad=15)
        
    full_title = title if title else "Seasonal Wind Roses"
    plt.suptitle(full_title, fontsize=13, y=1.01)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved: {save_path}")

    plt.show()