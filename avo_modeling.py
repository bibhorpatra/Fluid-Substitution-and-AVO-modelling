"""
avo_modeling.py
-----------------
Pre-stack synthetic seismic (angle-gather) modelling and AVO
intercept/gradient classification.

Reflectivity uses the Shuey (1985) three-term approximation to the
linearized Aki & Richards (1980) P-P reflection coefficient, which is
accurate for angles of incidence up to roughly 30-35 degrees -- adequate
for the near/mid/far stack ranges modelled here. AVO classification uses
the two-term intercept/gradient convention of Rutherford & Williams (1989),
extended with Castagna & Swan's (1997) Class IV designation.

References:
    Aki, K., and Richards, P. G., 1980, Quantitative Seismology: W. H. Freeman.
    Shuey, R. T., 1985, A simplification of the Zoeppritz equations:
        Geophysics, 50(4), 609-614.
    Rutherford, S. R., and Williams, R. H., 1989, Amplitude-versus-offset
        variations in gas sands: Geophysics, 54(6), 680-688.
    Castagna, J. P., and Swan, H. W., 1997, Principles of AVO crossplotting:
        The Leading Edge, 16(4), 337-342.
"""

import numpy as np


def ricker_wavelet(frequency_hz, dt_s=0.001, length_s=0.128):
    """Zero-phase Ricker wavelet, standard closed form."""
    t = np.arange(-length_s, length_s + dt_s / 2, dt_s)
    arg = (np.pi * frequency_hz * t) ** 2
    w = (1.0 - 2.0 * arg) * np.exp(-arg)
    return t, w


def shuey_reflectivity(vp1, vs1, rho1, vp2, vs2, rho2, angles_deg):
    """
    Three-term Shuey (1985) reflectivity vs. angle for a single interface
    (medium 1 over medium 2). Returns (R(theta), intercept A, gradient B,
    curvature C) evaluated at the supplied angles (degrees).
    """
    theta = np.deg2rad(np.asarray(angles_deg, dtype=float))

    vp_avg, vs_avg, rho_avg = (vp1 + vp2) / 2.0, (vs1 + vs2) / 2.0, (rho1 + rho2) / 2.0
    d_vp, d_vs, d_rho = vp2 - vp1, vs2 - vs1, rho2 - rho1

    A = 0.5 * (d_vp / vp_avg + d_rho / rho_avg)
    B = (0.5 * (d_vp / vp_avg)
         - 2.0 * (vs_avg / vp_avg) ** 2 * (d_rho / rho_avg + 2.0 * d_vs / vs_avg))
    C = 0.5 * (d_vp / vp_avg)

    R = A + B * np.sin(theta) ** 2 + C * (np.sin(theta) ** 2) * (np.tan(theta) ** 2)
    return R, A, B, C


def classify_avo(intercept, gradient):
    """
    Rutherford & Williams (1989) / Castagna & Swan (1997) AVO class from
    the two-term intercept (A) and gradient (B).
    """
    if intercept > 0.02 and gradient < 0:
        return "Class I (high-impedance sand, amplitude decreases with offset)"
    if abs(intercept) <= 0.02:
        return "Class II / IIp (near-zero impedance contrast, possible polarity reversal with offset)"
    if intercept < -0.02 and gradient < 0:
        return "Class III (low-impedance sand, amplitude increases with offset -- classic bright-spot)"
    if intercept < -0.02 and gradient > 0:
        return "Class IV (low-impedance sand, amplitude decreases in magnitude with offset)"
    return "Unclassified"


def build_angle_gather(interfaces, angles_deg, time_axis_s, wavelet):
    """
    Build a synthetic angle gather from one or more reflecting interfaces.

    Parameters
    ----------
    interfaces : list of dict
        Each dict has keys vp1, vs1, rho1 (layer above), vp2, vs2, rho2
        (layer below), and twt_s (two-way time of that interface).
    angles_deg : array-like
        Incidence angles to model (degrees).
    time_axis_s : ndarray
        Regularly-sampled time axis for the output trace (seconds).
    wavelet : ndarray
        Zero-phase wavelet (e.g. from ricker_wavelet), same dt as time_axis_s.

    Returns
    -------
    gather : ndarray, shape (len(time_axis_s), len(angles_deg))
    top_interface_ABC : tuple (A, B, C) for the first interface in the list,
        for AVO intercept/gradient classification.
    """
    dt = time_axis_s[1] - time_axis_s[0]
    n_t = time_axis_s.size
    n_ang = len(angles_deg)

    reflectivity = np.zeros((n_t, n_ang))
    top_abc = None
    for k, iface in enumerate(interfaces):
        R, A, B, C = shuey_reflectivity(iface["vp1"], iface["vs1"], iface["rho1"],
                                         iface["vp2"], iface["vs2"], iface["rho2"],
                                         angles_deg)
        if k == 0:
            top_abc = (A, B, C)
        spike_idx = int(round((iface["twt_s"] - time_axis_s[0]) / dt))
        if 0 <= spike_idx < n_t:
            reflectivity[spike_idx, :] += R

    gather = np.zeros((n_t, n_ang))
    for j in range(n_ang):
        gather[:, j] = np.convolve(reflectivity[:, j], wavelet, mode="same")
    return gather, top_abc


def stack_gather(gather, angles_deg, angle_range):
    lo, hi = angle_range
    sel = (np.asarray(angles_deg) >= lo) & (np.asarray(angles_deg) <= hi)
    return gather[:, sel].mean(axis=1)
