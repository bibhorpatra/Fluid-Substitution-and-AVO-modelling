"""
fluid_properties.py
--------------------
Pore-fluid density and bulk-modulus estimation for brine, oil, and gas at
in-situ reservoir pressure and temperature.

Brine and pure-water properties use the Batzle & Wang (1992) empirical
relations exactly as implemented in the open-source `bruges` rock-physics
library (bruges.rockphysics.fluids), which reproduces equations 27-29 of
that paper. Oil and gas properties are estimated with simplified,
clearly-flagged representative PVT relations, since no measured PVT report
(GOR, bubble point, gas composition) exists for this synthetic case. The
simplifications are documented at the point of use and are only applied to
the hydrocarbon phase being substituted in, never to the brine baseline.

Reference:
    Batzle, M., and Wang, Z., 1992, Seismic properties of pore fluids:
    Geophysics, 57(11), 1396-1408.
"""

import numpy as np
from bruges.rockphysics import fluids as _bw

GAS_CONSTANT = 8.3144598  # J / (mol K)


def brine_properties(temperature_c, pressure_pa, salinity_ppm=35_000.0):
    """
    Brine density (kg/m3) and bulk modulus (Pa) via Batzle & Wang (1992),
    eqs. 27b and 29, as implemented in bruges.

    Parameters
    ----------
    temperature_c : float or ndarray
        In-situ temperature, degrees C.
    pressure_pa : float or ndarray
        In-situ pore pressure, Pa.
    salinity_ppm : float
        NaCl salinity in parts per million (default 35,000 ppm, ~seawater).

    Returns
    -------
    rho : ndarray, kg/m3
    k   : ndarray, Pa
    """
    salinity_frac = salinity_ppm * 1e-6
    # NOTE: bruges.rockphysics.fluids.rho_brine returns density in g/cc
    # (its salinity-correction coefficients from Batzle & Wang, 1992, eq.
    # 27b are calibrated in that unit), despite a docstring that says
    # kg/m3 -- confirmed by comparing it against rho_water(), which does
    # return kg/m3 as documented. Convert to kg/m3 here so every fluid and
    # mineral property in this project is consistently SI (kg/m3, Pa).
    rho_g_cc = _bw.rho_brine(temperature_c, pressure_pa, salinity_frac)
    rho = rho_g_cc * 1000.0
    v = _bw.v_brine(temperature_c, pressure_pa, salinity_frac)
    k = rho * v**2
    return rho, k


def dead_oil_properties(temperature_c, pressure_pa, api_gravity=32.0,
                         reference_vp=1370.0, thermal_expansion=7.0e-4,
                         compressibility=5.0e-4):
    """
    Simplified dead (gas-free) oil density and bulk modulus at reservoir
    conditions.

    Density is computed exactly from the API-gravity definition at stock-tank
    conditions, then adjusted to reservoir temperature and pressure with a
    simple linear thermal-expansion / compressibility correction (typical
    coefficients for a medium-gravity dead oil). Acoustic velocity is taken
    as a representative literature value for a dead oil of this gravity at
    reservoir conditions (consistent with the range reported by Batzle &
    Wang, 1992, and standard rock-physics reference tables), rather than
    the full live-oil Batzle-Wang correlation, since no measured GOR or gas
    gravity is available for this synthetic case.

    Parameters
    ----------
    temperature_c : float
        In-situ temperature, degrees C.
    pressure_pa : float
        In-situ pore pressure, Pa.
    api_gravity : float
        Stock-tank oil API gravity (default 32 deg API, a medium-gravity oil).
    reference_vp : float
        Representative dead-oil acoustic velocity at reservoir conditions, m/s.
    thermal_expansion : float
        Linear thermal expansion coefficient, per degree C.
    compressibility : float
        Linear compressibility coefficient, per MPa.

    Returns
    -------
    rho : float, kg/m3
    k   : float, Pa
    """
    rho_stock = (141.5 / (131.5 + api_gravity)) * 1000.0  # kg/m3, exact API definition
    pressure_mpa = pressure_pa * 1e-6
    rho = rho_stock * (1.0 - thermal_expansion * (temperature_c - 15.6)
                        + compressibility * pressure_mpa)
    k = rho * reference_vp**2
    return rho, k


def gas_properties(temperature_c, pressure_pa, gas_gravity=0.65, heat_capacity_ratio=1.25):
    """
    Simplified reservoir-gas density and bulk modulus under an ideal-gas
    approximation, consistent with the ideal-gas density relation used by
    bruges.rockphysics.fluids.rho_gas.

    Density: rho = M P / (R T)                      (ideal gas law)
    Modulus: K = gamma * P                            (ideal-gas adiabatic bulk modulus)

    Parameters
    ----------
    temperature_c : float
    pressure_pa : float
    gas_gravity : float
        Gas specific gravity relative to air (default 0.65, typical for a
        light hydrocarbon reservoir gas dominated by methane).
    heat_capacity_ratio : float
        Effective ratio of specific heats, gamma (default 1.25).

    Returns
    -------
    rho : float, kg/m3
    k   : float, Pa
    """
    molecular_weight_g_mol = gas_gravity * 28.97  # g/mol (28.97 g/mol = air)
    # NOTE: bruges.rockphysics.fluids.rho_gas takes molecular weight in
    # g/mol but, because it internally rescales pressure from Pa to MPa
    # without a compensating factor, returns a value 1000x too small to be
    # kg/m3 (verified against the ideal-gas law for air at STP, which
    # should give ~1.225 kg/m3). Correct for that here.
    rho = _bw.rho_gas(temperature_c, pressure_pa, molecular_weight_g_mol) * 1000.0
    k = heat_capacity_ratio * pressure_pa
    return rho, k


def in_situ_conditions(depth_m, surface_temp_c=8.0, geothermal_gradient_c_per_km=30.0,
                        pressure_gradient_mpa_per_km=10.0):
    """
    Simple hydrostatic pressure / linear geothermal gradient model for a
    normally-pressured North Sea Tertiary clastic sequence.

    Returns (temperature_c, pressure_pa).
    """
    temperature_c = surface_temp_c + geothermal_gradient_c_per_km * (depth_m / 1000.0)
    pressure_pa = pressure_gradient_mpa_per_km * (depth_m / 1000.0) * 1e6
    return temperature_c, pressure_pa
