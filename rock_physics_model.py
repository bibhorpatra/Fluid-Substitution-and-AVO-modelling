"""
rock_physics_model.py
----------------------
Mineral-mixing, dry-frame, and Gassmann fluid-substitution equations for a
shaly-sand clastic system.

Mineral moduli are mixed with the Voigt-Reuss-Hill (VRH) average (Hill,
1952). The dry-rock frame is modelled with the Krief et al. (1990)
critical-porosity-style relation, which links dry-frame moduli to the
mineral moduli through porosity alone -- appropriate for a first-pass,
log-based rock-physics workflow where no core-calibrated contact-cement
model is available. Fluid substitution follows Gassmann (1951).

References:
    Hill, R., 1952, The elastic behaviour of a crystalline aggregate:
        Proc. Phys. Soc. London, A65, 349-354.
    Krief, M., Garat, J., Stellingwerff, J., and Ventre, J., 1990,
        A petrophysical interpretation using the velocities of P and S
        waves (full waveform sonic): The Log Analyst, 31, 355-369.
    Gassmann, F., 1951, Uber die elastizitat poroser medien:
        Vierteljahrsschrift der Naturforschenden Gesellschaft in Zurich,
        96, 1-23.
"""

import numpy as np

# Mineral end-member elastic moduli and densities, typical values reported
# in the rock-physics literature (e.g. Mavko, Mukerji & Dvorkin, 2009, The
# Rock Physics Handbook).
QUARTZ = dict(K=36.6e9, MU=45.0e9, RHO=2650.0)   # Pa, Pa, kg/m3
CLAY = dict(K=20.9e9, MU=6.85e9, RHO=2580.0)      # Pa, Pa, kg/m3

KRIEF_EXPONENT_A = 3.0  # standard default exponent in the Krief relation


def voigt_reuss_hill(f1, m1, f2, m2):
    """VRH average of a two-mineral mixture with fractions f1, f2 (f1+f2=1)
    and moduli m1, m2."""
    voigt = f1 * m1 + f2 * m2
    reuss = 1.0 / (f1 / m1 + f2 / m2)
    return 0.5 * (voigt + reuss)


def mineral_moduli(vshale):
    """
    Effective mineral bulk modulus, shear modulus, and density for a
    quartz-clay mixture, given the shale volume fraction (Vshale).
    """
    vshale = np.asarray(vshale, dtype=float)
    vsand = 1.0 - vshale
    k_min = voigt_reuss_hill(vsand, QUARTZ["K"], vshale, CLAY["K"])
    mu_min = voigt_reuss_hill(vsand, QUARTZ["MU"], vshale, CLAY["MU"])
    rho_min = vsand * QUARTZ["RHO"] + vshale * CLAY["RHO"]
    return k_min, mu_min, rho_min


def krief_dry_frame(k_mineral, mu_mineral, porosity, a=KRIEF_EXPONENT_A):
    """
    Krief et al. (1990) dry-frame moduli:
        K_dry / K_mineral = mu_dry / mu_mineral = (1 - phi) ** (a / (1 - phi))
    """
    porosity = np.clip(np.asarray(porosity, dtype=float), 1e-6, 0.999)
    exponent = a / (1.0 - porosity)
    factor = (1.0 - porosity) ** exponent
    return k_mineral * factor, mu_mineral * factor


def gassmann_saturated_bulk_modulus(k_dry, k_mineral, k_fluid, porosity):
    """
    Gassmann (1951) fluid-substitution equation for the saturated-rock bulk
    modulus:

        K_sat = K_dry + (1 - K_dry/K_min)^2 / (phi/K_fl + (1-phi)/K_min - K_dry/K_min^2)

    Shear modulus is fluid-independent under Gassmann's assumptions
    (mu_sat = mu_dry) and is not recomputed here.
    """
    porosity = np.asarray(porosity, dtype=float)
    numerator = (1.0 - k_dry / k_mineral) ** 2
    denominator = (porosity / k_fluid) + ((1.0 - porosity) / k_mineral) - (k_dry / k_mineral**2)
    return k_dry + numerator / denominator


def saturated_density(porosity, rho_mineral, rho_fluid):
    """Bulk density of the fluid-saturated rock."""
    porosity = np.asarray(porosity, dtype=float)
    return (1.0 - porosity) * rho_mineral + porosity * rho_fluid


def vp_from_moduli(k_sat, mu, rho):
    return np.sqrt((k_sat + 4.0 / 3.0 * mu) / rho)


def vs_from_moduli(mu, rho):
    return np.sqrt(mu / rho)


def acoustic_impedance(vp, rho):
    return vp * rho


def lambda_rho(vp, vs, rho):
    """Lambda-Rho (incompressibility) attribute: rho * (Vp^2 - 2 Vs^2)."""
    return rho * (vp**2 - 2.0 * vs**2)


def mu_rho(vs, rho):
    """Mu-Rho (rigidity) attribute: rho * Vs^2."""
    return rho * vs**2


def gassmann_fluid_substitution(vp_in, vs_in, rho_in, porosity, vshale,
                                 k_fluid_in, rho_fluid_in,
                                 k_fluid_out, rho_fluid_out, a=KRIEF_EXPONENT_A):
    """
    Full forward + inverse Gassmann fluid-substitution workflow starting
    from measured/modelled logs saturated with an initial fluid, ending
    with the rock re-saturated with a target fluid.

    Steps
    -----
    1. Estimate mineral moduli (VRH, quartz-clay mix) from Vshale.
    2. Back out K_sat(in) and mu from the input Vp, Vs, rho logs.
    3. Remove the initial fluid effect (inverse Gassmann) to get K_dry.
       (mu is assumed fluid-independent and is reused directly.)
    4. Re-saturate with the target fluid (forward Gassmann) to get K_sat(out).
    5. Recompute rho, Vp, Vs for the new fluid.

    Using measured logs to back out K_dry (rather than only the Krief model)
    means the substitution honours the actual measured/modelled elastic
    state of the rock; the Krief relation is used only as a consistency
    check and for extrapolating the rock-physics template beyond the
    logged interval (see rock_physics_template.py-equivalent routine below).

    Returns a dict of the output-fluid Vp, Vs, rho, K_sat, mu, K_dry.
    """
    k_min, mu_min, rho_min = mineral_moduli(vshale)

    mu = rho_in * vs_in**2
    k_sat_in = rho_in * vp_in**2 - 4.0 / 3.0 * mu

    # Inverse Gassmann: recover K_dry from K_sat with the initial fluid.
    numerator = k_sat_in * (porosity * k_min / k_fluid_in + 1.0 - porosity) - k_min
    denominator = porosity * k_min / k_fluid_in + k_sat_in / k_min - 1.0 - porosity
    k_dry = numerator / denominator

    # Forward Gassmann with the target fluid. The rock frame (mineral volume)
    # is unchanged by fluid substitution, so the new bulk density is simply
    # the input density with the pore-filling fluid swapped out:
    #   rho_in  = (1-phi) rho_mineral + phi rho_fluid_in
    #   rho_out = (1-phi) rho_mineral + phi rho_fluid_out
    #           = rho_in - phi (rho_fluid_in - rho_fluid_out)
    k_sat_out = gassmann_saturated_bulk_modulus(k_dry, k_min, k_fluid_out, porosity)
    rho_out = rho_in - porosity * (rho_fluid_in - rho_fluid_out)

    vp_out = vp_from_moduli(k_sat_out, mu, rho_out)
    vs_out = vs_from_moduli(mu, rho_out)

    return dict(vp=vp_out, vs=vs_out, rho=rho_out, k_sat=k_sat_out, mu=mu,
                k_dry=k_dry, k_mineral=k_min, mu_mineral=mu_min, rho_mineral=rho_min)
