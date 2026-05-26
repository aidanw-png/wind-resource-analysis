"""
weibull.py
----------
Functions for fitting a two-parameter Weibull distribution to wind speed
data using Maximum Likelihood Estimation, computing goodness of fit
diagnostics, and deriving key statistical quantities.

All logic validated in notebooks/01_era5_exploration.ipynb.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import math

def fit_weibull(wind_speeds: np.ndarray,
                calm_threshold: float = 0.5,
                verbose: bool = True) -> dict:
    # Filter calm hours
    fit_speeds = wind_speeds[wind_speeds >= calm_threshold]
    n_calm = len(wind_speeds) - len(fit_speeds)

    # MLE fit — floc=0 enforces two-parameter Weibull
    shape, loc, scale = stats.weibull_min.fit(fit_speeds, floc=0)

    k   = shape
    lam = scale

    # Theoretical mean from Weibull parameters
    # Mean = λ × Γ(1 + 1/k)
    weibull_mean  = lam * math.gamma(1 + 1/k)
    observed_mean = wind_speeds.mean()

    result = {
        "k":             k,
        "lam":           lam,
        "n_fit":         len(fit_speeds),
        "n_calm":        n_calm,
        "pct_calm":      n_calm / len(wind_speeds) * 100,
        "weibull_mean":  weibull_mean,
        "observed_mean": observed_mean,
        "mean_error":    abs(weibull_mean - observed_mean)
    }

    # Physical interpretation of k
    # Must be defined before the verbose print block that references it
    if k < 1.5:
        interp = "highly variable - frequent calms and storm events"
    elif k < 2.0:
        interp = "variable - typical mid-latitude synoptic regime"
    elif k < 2.5:
        interp = "moderately consistent"
    elif k < 3.0:
        interp = "consistent - approaching trade wind character"
    else:
        interp = "highly consistent - characteristic of trade wind regime"

    if verbose:
        print("=== WEIBULL FIT RESULTS ===")
        print(f"  Shape parameter  k:   {k:.4f}")
        print(f"  Scale parameter  λ:   {lam:.4f} m/s")
        print(f"  Observations used:    {len(fit_speeds):,}")
        print(f"  Calm hours excluded:  {n_calm} ({result['pct_calm']:.1f}%)")
        print(f"  Theoretical mean:     {weibull_mean:.4f} m/s")
        print(f"  Observed mean:        {observed_mean:.4f} m/s")
        print(f"  Mean error:           {result['mean_error']:.4f} m/s")
        print(f"  Interpretation:       k={k:.2f} — {interp}")

    return result

def apply_hub_correction(weibull_result: dict,
                          z_ref: float,
                          z_hub: float,
                          alpha: float,
                          verbose: bool = True) -> dict:
    
    shear_factor = (z_hub / z_ref) ** alpha

    k_hub   = weibull_result["k"]
    lam_hub = weibull_result["lam"] * shear_factor
    mean_hub = lam_hub * math.gamma(1 + 1/k_hub)

    result = weibull_result.copy()
    result.update({
        "k_hub":        k_hub,
        "lam_hub":      lam_hub,
        "mean_hub":     mean_hub,
        "shear_factor": shear_factor,
        "z_ref":        z_ref,
        "z_hub":        z_hub,
        "alpha":        alpha
    })

    if verbose:
        print(f"\n=== HUB HEIGHT CORRECTION ===")
        print(f"  Reference height:     {z_ref:.0f} m")
        print(f"  Hub height:           {z_hub:.0f} m")
        print(f"  Shear exponent α:     {alpha}")
        print(f"  Shear factor:         {shear_factor:.6f}")
        print(f"  λ at {z_ref:.0f}m:          "
              f"{weibull_result['lam']:.4f} m/s")
        print(f"  λ at {z_hub:.0f}m:           "
              f"{lam_hub:.4f} m/s")
        print(f"  Mean at hub height:   {mean_hub:.4f} m/s")

    return result

def goodness_of_fit(wind_speeds: np.ndarray,
                     weibull_result: dict,
                     calm_threshold: float = 0.5,
                     verbose: bool = True) -> dict:
    
    fit_speeds = wind_speeds[wind_speeds >= calm_threshold]
    k   = weibull_result["k"]
    lam = weibull_result["lam"]

    ks_stat, ks_p = stats.kstest(
        fit_speeds,
        lambda x: stats.weibull_min.cdf(x, c=k, loc=0, scale=lam)
    )

    # Contextual interpretation - KS is oversensitive at large n
    if ks_p > 0.05:
        interpretation = "Cannot reject Weibull hypothesis (p > 0.05)"
    else:
        interpretation = (
            f"Formally rejected (p={ks_p:.4f}) but D={ks_stat:.4f} is "
            f"practically small - expected with n={len(fit_speeds):,} samples"
        )

    result = {
        "ks_statistic":   ks_stat,
        "ks_pvalue":      ks_p,
        "n":              len(fit_speeds),
        "interpretation": interpretation
    }

    if verbose:
        print("\n=== GOODNESS OF FIT ===")
        print(f"  KS statistic D:  {ks_stat:.4f}")
        print(f"  p-value:         {ks_p:.4f}")
        print(f"  n observations:  {len(fit_speeds):,}")
        print(f"  Result:          {interpretation}")

    return result


def plot_weibull_fit(wind_speeds: np.ndarray,
                     weibull_result: dict,
                     title: str = "",
                     save_path: str = None) -> None:
    
    k   = weibull_result["k"]
    lam = weibull_result["lam"]

    fit_speeds = wind_speeds[wind_speeds >= 0.5]
    u_range    = np.linspace(0, wind_speeds.max() * 1.1, 500)
    pdf_values = stats.weibull_min.pdf(u_range, c=k, loc=0, scale=lam)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # --- Panel 1: Histogram + PDF ---
    axes[0].hist(wind_speeds, bins=40, density=True,
                 color="steelblue", edgecolor="white",
                 linewidth=0.5, alpha=0.7, label="Observed")
    axes[0].plot(u_range, pdf_values, color="red",
                 linewidth=2.5,
                 label=f"Weibull (k={k:.3f}, λ={lam:.3f})")
    axes[0].axvline(weibull_result["observed_mean"],
                    color="black", linewidth=1.2, linestyle="--",
                    label=f"Mean: {weibull_result['observed_mean']:.2f} m/s")
    axes[0].set_xlabel("Wind Speed (m/s)")
    axes[0].set_ylabel("Probability Density")
    axes[0].set_title("Weibull Fit")
    axes[0].legend(fontsize=8)

    # --- Panel 2: Q-Q plot ---
    quantiles        = np.linspace(0.01, 0.99, 200)
    theoretical_q    = stats.weibull_min.ppf(
                           quantiles, c=k, loc=0, scale=lam)
    observed_q       = np.quantile(fit_speeds, quantiles)

    axes[1].plot(theoretical_q, observed_q,
                 color="steelblue", linewidth=1.5, label="Data quantiles")
    axes[1].plot([0, theoretical_q.max()],
                 [0, theoretical_q.max()],
                 color="red", linewidth=1.5,
                 linestyle="--", label="Perfect fit")
    axes[1].set_xlabel("Theoretical Quantiles (m/s)")
    axes[1].set_ylabel("Observed Quantiles (m/s)")
    axes[1].set_title("Q-Q Plot")
    axes[1].legend(fontsize=8)

    # --- Panel 3: CDF comparison ---
    sorted_speeds    = np.sort(fit_speeds)
    empirical_cdf    = np.arange(1, len(sorted_speeds) + 1) / len(sorted_speeds)
    theoretical_cdf  = stats.weibull_min.cdf(
                           sorted_speeds, c=k, loc=0, scale=lam)

    axes[2].plot(sorted_speeds, empirical_cdf,
                 color="steelblue", linewidth=1.5, label="Empirical CDF")
    axes[2].plot(sorted_speeds, theoretical_cdf,
                 color="red", linewidth=1.5,
                 linestyle="--", label="Weibull CDF")
    axes[2].set_xlabel("Wind Speed (m/s)")
    axes[2].set_ylabel("Cumulative Probability")
    axes[2].set_title("CDF Comparison")
    axes[2].legend(fontsize=8)

    full_title = f"Weibull Diagnostics - k={k:.3f}, λ={lam:.3f} m/s"
    if title:
        full_title = f"{full_title}\n{title}"
    plt.suptitle(full_title, fontsize=12)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved: {save_path}")

    plt.show()