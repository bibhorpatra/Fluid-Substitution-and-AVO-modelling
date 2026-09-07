"""
synthetic_well.py
-------------------
Builds a physically-consistent synthetic type-well log dataset for a
shaly-sand Tertiary clastic sequence, representative of a shallow,
under-compacted North Sea-type reservoir setting (e.g. Utsira/Hordaland
Group-style high-porosity Tertiary sands).

No proprietary or measured well data is used anywhere in this project --
every curve here is generated from standard empirical/compaction relations
(Athy's law porosity-depth trend, a quartz-clay VRH mineral mixture, the
Krief dry-frame relation, and Gassmann fluid substitution), calibrated to
literature-typical parameter ranges for this depositional setting. This is
explicitly a synthetic "type well", not a claimed real-field dataset.

Reference:
    Athy, L. F., 1930, Density, porosity, and compaction of sedimentary
    rocks: AAPG Bulletin, 14, 1-24.
"""

import numpy as np
import pandas as pd

from . import fluid_properties as fl
from . import rock_physics_model as rpm

# ---------------------------------------------------------------------------
# Well geometry
# ---------------------------------------------------------------------------
DEPTH_TOP_M = 2000.0
DEPTH_BASE_M = 2200.0
DEPTH_STEP_M = 0.1524  # 0.5 ft, a typical wireline log sample rate

RESERVOIR_TOP_M = 2140.0
RESERVOIR_BASE_M = 2156.0

# Secondary, non-reservoir silty streaks included only for log realism.
SILT_STREAKS = [(2033.0, 2039.0), (2172.0, 2178.0)]

SALINITY_PPM = 35_000.0
OIL_API = 32.0
GAS_GRAVITY = 0.65

RANDOM_SEED = 42


def _logistic_step(z, center, width, low, high, invert=False):
    """Smooth sigmoidal transition between `low` and `high` around `center`."""
    s = 1.0 / (1.0 + np.exp(-(z - center) / (width / 4.0)))
    if invert:
        s = 1.0 - s
    return low + (high - low) * s


def _build_vshale(depth):
    rng = np.random.default_rng(RANDOM_SEED)

    background = 0.82 + 0.05 * np.sin(2 * np.pi * depth / 150.0)

    # Main reservoir sand body: a smooth, log-resolution-limited transition
    # down into a clean sand at the top and back up to shale at the base.
    top_step = _logistic_step(depth, RESERVOIR_TOP_M, 1.0, 0.82, 0.08)
    base_step = _logistic_step(depth, RESERVOIR_BASE_M, 1.0, 0.08, 0.82)
    midpoint = (RESERVOIR_TOP_M + RESERVOIR_BASE_M) / 2.0
    reservoir_shape = np.where(depth < midpoint, top_step, base_step)
    in_reservoir = (depth >= RESERVOIR_TOP_M - 3) & (depth <= RESERVOIR_BASE_M + 3)
    vsh = np.where(in_reservoir, reservoir_shape, background)

    # Thin non-pay silty streaks elsewhere in the section, for log realism.
    for lo, hi in SILT_STREAKS:
        center, width = (lo + hi) / 2.0, (hi - lo)
        streak = 0.82 - 0.47 * np.exp(-0.5 * ((depth - center) / (width / 2.5)) ** 2)
        near_streak = (depth >= lo - 2) & (depth <= hi + 2)
        vsh = np.where(near_streak, streak, vsh)

    # Correlated (band-limited) noise for a realistic log texture.
    from scipy.ndimage import gaussian_filter1d
    noise = gaussian_filter1d(rng.normal(0, 1, size=depth.size), sigma=3.0)
    noise = 0.035 * noise / np.std(noise)
    vsh = np.clip(vsh + noise, 0.02, 0.98)
    return vsh


def _athy_porosity_trend(depth, phi0, decay_per_m):
    return phi0 * np.exp(-decay_per_m * depth)


def _build_porosity(depth, vshale):
    rng = np.random.default_rng(RANDOM_SEED + 1)
    phi_sand_trend = _athy_porosity_trend(depth, phi0=0.45, decay_per_m=0.00012)
    phi_shale_trend = _athy_porosity_trend(depth, phi0=0.55, decay_per_m=0.00035)
    phi = (1.0 - vshale) * phi_sand_trend + vshale * phi_shale_trend

    from scipy.ndimage import gaussian_filter1d
    noise = gaussian_filter1d(rng.normal(0, 1, size=depth.size), sigma=3.0)
    noise = 0.010 * noise / np.std(noise)
    phi = np.clip(phi + noise, 0.02, 0.42)
    return phi


def build_synthetic_well():
    """
    Returns a pandas DataFrame with one row per depth sample, containing the
    lithology/porosity model, in-situ P/T, mineral and dry-frame moduli, the
    fully brine-saturated ("in-situ / pre-drill") elastic log, and -- within
    the reservoir interval only -- the Gassmann-substituted oil-case and
    gas-case elastic logs.
    """
    depth = np.arange(DEPTH_TOP_M, DEPTH_BASE_M + DEPTH_STEP_M / 2, DEPTH_STEP_M)

    vshale = _build_vshale(depth)
    porosity = _build_porosity(depth, vshale)

    k_min, mu_min, rho_min = rpm.mineral_moduli(vshale)
    k_dry, mu_dry = rpm.krief_dry_frame(k_min, mu_min, porosity)

    temp_c, pres_pa = fl.in_situ_conditions(depth)
    rho_brine, k_brine = fl.brine_properties(temp_c, pres_pa, SALINITY_PPM)

    k_sat_brine = rpm.gassmann_saturated_bulk_modulus(k_dry, k_min, k_brine, porosity)
    rho_sat_brine = rpm.saturated_density(porosity, rho_min, rho_brine)
    vp_brine = rpm.vp_from_moduli(k_sat_brine, mu_dry, rho_sat_brine)
    vs_brine = rpm.vs_from_moduli(mu_dry, rho_sat_brine)

    df = pd.DataFrame({
        "DEPTH_M": depth,
        "VSHALE": vshale,
        "PHIE": porosity,
        "TEMP_C": temp_c,
        "PRESSURE_PA": pres_pa,
        "K_MINERAL_PA": k_min,
        "MU_MINERAL_PA": mu_min,
        "RHO_MINERAL_KGM3": rho_min,
        "K_DRY_PA": k_dry,
        "MU_DRY_PA": mu_dry,
        "RHO_BRINE_FLUID_KGM3": rho_brine,
        "K_BRINE_FLUID_PA": k_brine,
        "VP_BRINE_MPS": vp_brine,
        "VS_BRINE_MPS": vs_brine,
        "RHO_BRINE_KGM3": rho_sat_brine,
    })
    df["AI_BRINE"] = rpm.acoustic_impedance(df["VP_BRINE_MPS"], df["RHO_BRINE_KGM3"])
    df["VPVS_BRINE"] = df["VP_BRINE_MPS"] / df["VS_BRINE_MPS"]

    # --- Fluid substitution within the reservoir interval only -------------
    reservoir_mask = ((df["DEPTH_M"] >= RESERVOIR_TOP_M) &
                       (df["DEPTH_M"] <= RESERVOIR_BASE_M) &
                       (df["VSHALE"] < 0.2)).to_numpy()

    for label, api, gravity in (("OIL", OIL_API, None), ("GAS", None, GAS_GRAVITY)):
        vp_col, vs_col, rho_col = f"VP_{label}_MPS", f"VS_{label}_MPS", f"RHO_{label}_KGM3"
        df[vp_col] = np.nan
        df[vs_col] = np.nan
        df[rho_col] = np.nan

        idx = np.where(reservoir_mask)[0]
        if idx.size == 0:
            continue

        t_res = df["TEMP_C"].to_numpy()[idx]
        p_res = df["PRESSURE_PA"].to_numpy()[idx]
        if label == "OIL":
            rho_fl_out, k_fl_out = fl.dead_oil_properties(t_res, p_res, api_gravity=api)
        else:
            rho_fl_out, k_fl_out = fl.gas_properties(t_res, p_res, gas_gravity=gravity)

        k_dry_r = df["K_DRY_PA"].to_numpy()[idx]
        k_min_r = df["K_MINERAL_PA"].to_numpy()[idx]
        mu_r = df["MU_DRY_PA"].to_numpy()[idx]
        phi_r = df["PHIE"].to_numpy()[idx]
        rho_brine_r = df["RHO_BRINE_FLUID_KGM3"].to_numpy()[idx]
        rho_in_r = df["RHO_BRINE_KGM3"].to_numpy()[idx]

        k_sat_out = rpm.gassmann_saturated_bulk_modulus(k_dry_r, k_min_r, k_fl_out, phi_r)
        rho_out = rho_in_r - phi_r * (rho_brine_r - rho_fl_out)
        vp_out = rpm.vp_from_moduli(k_sat_out, mu_r, rho_out)
        vs_out = rpm.vs_from_moduli(mu_r, rho_out)

        df.loc[idx, vp_col] = vp_out
        df.loc[idx, vs_col] = vs_out
        df.loc[idx, rho_col] = rho_out

    df["AI_OIL"] = rpm.acoustic_impedance(df["VP_OIL_MPS"], df["RHO_OIL_KGM3"])
    df["VPVS_OIL"] = df["VP_OIL_MPS"] / df["VS_OIL_MPS"]
    df["AI_GAS"] = rpm.acoustic_impedance(df["VP_GAS_MPS"], df["RHO_GAS_KGM3"])
    df["VPVS_GAS"] = df["VP_GAS_MPS"] / df["VS_GAS_MPS"]

    df.attrs["reservoir_top_m"] = RESERVOIR_TOP_M
    df.attrs["reservoir_base_m"] = RESERVOIR_BASE_M
    return df


def reservoir_average_properties(df, fluid, zone_top, zone_base):
    """Average Vp, Vs, rho over a clean sub-interval (avoids boundary/tuning
    transition samples), for use as a single AVO-modelling layer property."""
    mask = (df["DEPTH_M"] >= zone_top) & (df["DEPTH_M"] <= zone_base)
    vp = df.loc[mask, f"VP_{fluid}_MPS"].mean()
    vs = df.loc[mask, f"VS_{fluid}_MPS"].mean()
    rho = df.loc[mask, f"RHO_{fluid}_KGM3"].mean()
    return vp, vs, rho
