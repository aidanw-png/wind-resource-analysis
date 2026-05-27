"""
aep.py
------
Functions for defining turbine power curves, applying hub height wind
shear corrections, and computing Annual Energy Production (AEP) from
Weibull distribution parameters.

All logic validated in notebooks/01_era5_exploration.ipynb.

Reference turbine: NREL 5MW (Jonkman et al., 2009, NREL/TP-500-38060)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, integrate
from scipy.interpolate import interp1d
from pathlib import Path
import math

# ── Power curve data ────────────────────────────────────────────────────────
# NREL 5MW Reference Turbine
# Source: Jonkman et al. (2009) NREL/TP-500-38060, Table 3-2
# Wind speed (m/s) : Power output (kW)
NREL_5MW_CURVE = {
    "wind_speed": [
         0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,
         9.0, 10.0, 11.0, 11.4, 12.0, 13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 25.01
    ],
    "power_kw": [
           0,    0,    0,   40,  178,  403,  737, 1116, 1582,
        2100, 2665, 3248, 3500, 3500, 3500, 3500, 3500, 3500,
        3500, 3500, 3500, 3500, 3500, 3500, 3500, 3500, 3500,    0
    ],
    "hub_height_m":   90.0,
    "rotor_diameter_m": 126.0,
    "rated_power_kw": 3500.0,
    "cut_in_ms":       3.0,
    "rated_ms":       11.4,
    "cut_out_ms":     25.0
}

def build_power_curve(curve_data: dict = None) -> interp1d:
    if curve_data is None:
        curve_data = NREL_5MW_CURVE

    u = np.array(curve_data["wind_speed"])
    p = np.array(curve_data["power_kw"])

    power_fn = interp1d(
        u, p,
        kind="linear",
        bounds_error=False,
        fill_value=0.0
    )

    return power_fn

def calculate_aep(k_hub: float,
                  lam_hub: float,
                  power_fn: interp1d,
                  availability: float = 0.95,
                  rated_power_mw: float = 5.0,
                  u_min: float = 0.01,
                  u_max: float = 26.0,
                  n_points: int = 10000) -> dict:
    hours_per_year = 8760

    # Integration grid
    u_range = np.linspace(u_min, u_max, n_points)

    # Weibull PDF at hub height
    pdf = stats.weibull_min.pdf(u_range, c=k_hub, loc=0, scale=lam_hub)

    # Power curve evaluated at each speed
    power = power_fn(u_range)

    # Integrand: expected power contribution at each speed
    integrand = power * pdf

    # Numerical integration using Simpson's rule
    mean_power_kw = integrate.simpson(integrand, x=u_range)

    # Scale to annual energy
    gross_aep_mwh = mean_power_kw * hours_per_year / 1000
    gross_aep_gwh = gross_aep_mwh / 1000
    net_aep_mwh   = gross_aep_mwh * availability
    net_aep_gwh   = net_aep_mwh / 1000

    # Capacity factor - actual output / maximum possible output
    max_annual_mwh = rated_power_mw * hours_per_year
    cf_gross = gross_aep_mwh / max_annual_mwh
    cf_net   = net_aep_mwh   / max_annual_mwh

    result = {
        "mean_power_kw":          mean_power_kw,
        "gross_aep_mwh":          gross_aep_mwh,
        "gross_aep_gwh":          gross_aep_gwh,
        "net_aep_mwh":            net_aep_mwh,
        "net_aep_gwh":            net_aep_gwh,
        "capacity_factor_gross":  cf_gross,
        "capacity_factor_net":    cf_net,
        "availability":           availability,
        "k_hub":                  k_hub,
        "lam_hub":                lam_hub
    }

    print("=== AEP RESULTS ===")
    print(f"  Mean power:             {mean_power_kw:.1f} kW")
    print(f"  Gross AEP:              {gross_aep_gwh:.3f} GWh/year")
    print(f"  Net AEP ({availability:.0%} avail):   "
          f"{net_aep_gwh:.3f} GWh/year")
    print(f"  Gross capacity factor:  {cf_gross:.1%}")
    print(f"  Net capacity factor:    {cf_net:.1%}")

    return result

def plot_aep_analysis(k_hub: float,
                      lam_hub: float,
                      power_fn: interp1d,
                      aep_result: dict,
                      title: str = "",
                      save_path: str = None) -> None:
    
    u_range   = np.linspace(0.01, 26.0, 10000)
    pdf       = stats.weibull_min.pdf(u_range, c=k_hub, loc=0, scale=lam_hub)
    power     = power_fn(u_range)
    integrand = power * pdf

    mean_hub  = lam_hub * math.gamma(1 + 1/k_hub)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # --- Left: AEP integrand ---
    ax = axes[0]
    ax.plot(u_range, integrand, color="steelblue",
            linewidth=2.0, label="P(U) x f(U)")
    ax.fill_between(u_range, integrand, alpha=0.3, color="steelblue")
    ax.fill_between(u_range, integrand,
                    where=(u_range < 3.0),
                    alpha=0.5, color="lightgrey",
                    label="Below cut-in")
    ax.axvline(3.0,    color="green",  linewidth=1.2,
               linestyle=":", label="Cut-in (3 m/s)")
    ax.axvline(11.4,   color="purple", linewidth=1.2,
               linestyle=":", label="Rated (11.4 m/s)")
    ax.axvline(25.0,   color="red",    linewidth=1.2,
               linestyle=":", label="Cut-out (25 m/s)")
    ax.axvline(mean_hub, color="black", linewidth=1.2,
               linestyle="--",
               label=f"Mean ({mean_hub:.1f} m/s)")
    ax.set_xlabel("Wind Speed at Hub Height (m/s)")
    ax.set_ylabel("Power x Probability Density (kW per m/s)")
    ax.set_title("AEP Integrand - Energy Contribution by Speed")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 27)

    # --- Right: Cumulative AEP ---
    ax2 = axes[1]
    cum = np.cumsum(integrand) * 26.0 / 10000
    cum_frac = cum / cum[-1]

    ax2.plot(u_range, cum_frac * 100,
             color="steelblue", linewidth=2.0)

    for target, colour in [(50, "orange"), (80, "red")]:
        idx = np.searchsorted(cum_frac, target / 100)
        if idx < len(u_range):
            spd = u_range[idx]
            ax2.axhline(target, color=colour, linewidth=1.0,
                        linestyle="--", alpha=0.7)
            ax2.axvline(spd, color=colour, linewidth=1.0,
                        linestyle="--", alpha=0.7,
                        label=f"{target}% AEP below {spd:.1f} m/s")
            
    ax2.axvline(3.0,  color="green",  linewidth=1.2,
                linestyle=":", label="Cut-in (3 m/s)")
    ax2.axvline(11.4, color="purple", linewidth=1.2,
                linestyle=":", label="Rated (11.4 m/s)")
    ax2.set_xlabel("Wind Speed at Hub Height (m/s)")
    ax2.set_ylabel("Cumulative AEP (%)")
    ax2.set_title("Cumulative AEP Distribution by Wind Speed")
    ax2.legend(fontsize=8)
    ax2.set_xlim(0, 27)
    ax2.set_ylim(0, 100)

    cf_net = aep_result["capacity_factor_net"]
    gross  = aep_result["gross_aep_gwh"]
    net    = aep_result["net_aep_gwh"]

    full_title = (f"AEP Analysis - Gross: {gross:.3f} GWh/year | "
                  f"Net: {net:.3f} GWh/year | CF: {cf_net:.1%}")
    if title:
        full_title = f"{full_title}\n{title}"

    plt.suptitle(full_title, fontsize=11, y=1.02)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved: {save_path}")

    plt.show()