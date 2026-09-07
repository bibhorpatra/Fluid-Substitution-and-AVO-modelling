"""
rock_physics_template.py
--------------------------
Constructs a Vp/Vs-vs-Acoustic-Impedance rock physics template (RPT) at
fixed reservoir pressure/temperature conditions, parametrised by porosity,
for the brine / oil / gas saturation cases and for a shale reference trend.
Also synthesizes sparse "core plug" calibration samples with measurement
scatter, in the style of a log-and-core calibrated RPT.
"""

import numpy as np
import pandas as pd

from . import rock_physics_model as rpm


def _forward_model(vshale, porosity, k_fluid, rho_fluid):
    vshale = np.full_like(porosity, vshale, dtype=float) if np.isscalar(vshale) else vshale
    k_min, mu_min, rho_min = rpm.mineral_moduli(vshale)
    k_dry, mu_dry = rpm.krief_dry_frame(k_min, mu_min, porosity)
    k_sat = rpm.gassmann_saturated_bulk_modulus(k_dry, k_min, k_fluid, porosity)
    rho_sat = rpm.saturated_density(porosity, rho_min, rho_fluid)
    vp = rpm.vp_from_moduli(k_sat, mu_dry, rho_sat)
    vs = rpm.vs_from_moduli(mu_dry, rho_sat)
    ai = rpm.acoustic_impedance(vp, rho_sat)
    return vp, vs, rho_sat, ai


def build_template_curves(porosity_range, sand_vshale, shale_vshale,
                           brine_props, oil_props, gas_props):
    """
    porosity_range : ndarray of porosities to sweep.
    *_props : (rho_fluid, k_fluid) tuples at the representative reservoir P/T.
    Returns a dict of DataFrames: 'brine_sand', 'oil_sand', 'gas_sand', 'shale'.
    """
    rho_b, k_b = brine_props
    rho_o, k_o = oil_props
    rho_g, k_g = gas_props

    curves = {}
    for label, vshale, (rho_fl, k_fl) in (
        ("brine_sand", sand_vshale, (rho_b, k_b)),
        ("oil_sand", sand_vshale, (rho_o, k_o)),
        ("gas_sand", sand_vshale, (rho_g, k_g)),
        ("shale", shale_vshale, (rho_b, k_b)),
    ):
        vp, vs, rho, ai = _forward_model(vshale, porosity_range, k_fl, rho_fl)
        curves[label] = pd.DataFrame({
            "PHIE": porosity_range, "VP_MPS": vp, "VS_MPS": vs,
            "RHO_KGM3": rho, "AI": ai, "VPVS": vp / vs,
        })
    return curves


def synthesize_core_samples(n_samples, porosity_range, sand_vshale, fluid_props,
                             porosity_noise_std=0.012, measurement_noise_frac=0.015,
                             seed=7):
    """
    Sparse synthetic "core plug" samples along the sand trend, with random
    porosity draws and small independent measurement noise on Vp, Vs, rho --
    representative of the scatter seen when calibrating a log-derived RPT
    against discrete core measurements.
    """
    rng = np.random.default_rng(seed)
    rho_fl, k_fl = fluid_props

    phi_samples = rng.uniform(porosity_range.min(), porosity_range.max(), n_samples)
    phi_samples = np.clip(phi_samples + rng.normal(0, porosity_noise_std, n_samples),
                           porosity_range.min(), porosity_range.max())

    vp, vs, rho, ai = _forward_model(sand_vshale, phi_samples, k_fl, rho_fl)
    vp *= 1 + rng.normal(0, measurement_noise_frac, n_samples)
    vs *= 1 + rng.normal(0, measurement_noise_frac * 1.3, n_samples)
    rho *= 1 + rng.normal(0, measurement_noise_frac * 0.4, n_samples)
    ai = vp * rho

    return pd.DataFrame({
        "PHIE_core": phi_samples, "VP_core": vp, "VS_core": vs,
        "RHO_core": rho, "AI_core": ai, "VPVS_core": vp / vs,
    })
