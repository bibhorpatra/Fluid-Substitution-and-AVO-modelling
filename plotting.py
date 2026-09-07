"""
plotting.py
-------------
All figure-generation routines for the fluid-substitution / AVO workflow.
Kept deliberately separate from the physics modules so the analysis script
and notebook can call a small, documented plotting API.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

FLUID_COLORS = {"BRINE": "#1f6fb2", "OIL": "#2e8b57", "GAS": "#c0392b", "SHALE": "#8d7860"}
FLUID_LABELS = {"BRINE": "Brine sand (in-situ)", "OIL": "Oil sand (substituted)",
                "GAS": "Gas sand (substituted)", "SHALE": "Shale"}

plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.dpi": 140,
    "font.size": 9.5,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "axes.axisbelow": True,
})


def _shade_reservoir(ax, top, base, orientation="y"):
    if orientation == "y":
        ax.axhspan(top, base, color="#f2c14e", alpha=0.15, zorder=0)
    else:
        ax.axvspan(top, base, color="#f2c14e", alpha=0.15, zorder=0)


def plot_log_tracks(df, reservoir_top, reservoir_base, out_path):
    fig, axes = plt.subplots(1, 4, figsize=(9.5, 8.5), sharey=True)
    depth = df["DEPTH_M"]

    axes[0].plot(df["VSHALE"], depth, color="#6b4f3a", lw=0.8)
    axes[0].set_xlabel("Vshale (v/v)")
    axes[0].set_xlim(0, 1)

    axes[1].plot(df["PHIE"], depth, color="#1f6fb2", lw=0.8)
    axes[1].set_xlabel("Porosity (v/v)")
    axes[1].set_xlim(0, 0.45)

    axes[2].plot(df["VP_BRINE_MPS"], depth, color="#222222", lw=0.8, label="Vp")
    axes[2].plot(df["VS_BRINE_MPS"], depth, color="#c0392b", lw=0.8, label="Vs")
    axes[2].set_xlabel("Velocity (m/s)")
    axes[2].legend(loc="lower right", fontsize=7, framealpha=0.85)

    axes[3].plot(df["RHO_BRINE_KGM3"] / 1000.0, depth, color="#2e8b57", lw=0.8)
    axes[3].set_xlabel("Bulk density (g/cc)")

    for ax in axes:
        _shade_reservoir(ax, reservoir_top, reservoir_base)
        ax.invert_yaxis()
    axes[0].set_ylabel("Depth (m)")
    fig.suptitle("Synthetic type-well: lithology, porosity, and in-situ (brine) elastic logs", y=0.995)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_depth_trends(df, reservoir_top, reservoir_base, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 8.0), sharey=True)
    depth = df["DEPTH_M"]

    axes[0].plot(df["AI_BRINE"] / 1000.0, depth, color=FLUID_COLORS["BRINE"], lw=0.9,
                 label=FLUID_LABELS["BRINE"])
    axes[0].plot(df["AI_OIL"] / 1000.0, depth, color=FLUID_COLORS["OIL"], lw=1.4,
                 label=FLUID_LABELS["OIL"])
    axes[0].plot(df["AI_GAS"] / 1000.0, depth, color=FLUID_COLORS["GAS"], lw=1.4,
                 label=FLUID_LABELS["GAS"])
    axes[0].set_xlabel(r"Acoustic impedance, AI (10$^3$ m/s $\cdot$ g/cc)")
    axes[0].legend(loc="upper left", fontsize=7, framealpha=0.85)

    axes[1].plot(df["VPVS_BRINE"], depth, color=FLUID_COLORS["BRINE"], lw=0.9)
    axes[1].plot(df["VPVS_OIL"], depth, color=FLUID_COLORS["OIL"], lw=1.4)
    axes[1].plot(df["VPVS_GAS"], depth, color=FLUID_COLORS["GAS"], lw=1.4)
    axes[1].set_xlabel("Vp/Vs ratio")

    for ax in axes:
        _shade_reservoir(ax, reservoir_top, reservoir_base)
        ax.invert_yaxis()
    axes[0].set_ylabel("Depth (m)")
    fig.suptitle("Depth trends: in-situ brine case vs. Gassmann-substituted oil/gas case", y=0.995)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_rock_physics_template(df, template_curves, core_samples, reservoir_top,
                                reservoir_base, out_path):
    fig, ax = plt.subplots(figsize=(7.2, 6.2))

    for key, style in (("shale", dict(color=FLUID_COLORS["SHALE"], ls="--")),
                        ("brine_sand", dict(color=FLUID_COLORS["BRINE"])),
                        ("oil_sand", dict(color=FLUID_COLORS["OIL"])),
                        ("gas_sand", dict(color=FLUID_COLORS["GAS"]))):
        curve = template_curves[key]
        ax.plot(curve["AI"] / 1000.0, curve["VPVS"], lw=2.0,
                label=f"{key.replace('_', ' ').title()} trend (porosity 0.05-0.40)", **style)

    reservoir = df[(df["DEPTH_M"] >= reservoir_top) & (df["DEPTH_M"] <= reservoir_base)]
    background = df[(df["DEPTH_M"] < reservoir_top) | (df["DEPTH_M"] > reservoir_base)]
    ax.scatter(background["AI_BRINE"] / 1000.0, background["VPVS_BRINE"], s=4, c="#999999",
               alpha=0.35, label="Background shale/silt log samples", zorder=2)
    ax.scatter(reservoir["AI_BRINE"] / 1000.0, reservoir["VPVS_BRINE"], s=8,
               c=FLUID_COLORS["BRINE"], alpha=0.6, label="Reservoir sand log samples (in-situ)",
               zorder=3)

    ax.scatter(core_samples["AI_core"] / 1000.0, core_samples["VPVS_core"], marker="D", s=32,
               facecolor="white", edgecolor="black", linewidth=0.8, zorder=4,
               label="Core plug calibration samples")

    ax.set_xlabel(r"Acoustic impedance, AI (10$^3$ m/s $\cdot$ g/cc)")
    ax.set_ylabel("Vp/Vs ratio")
    ax.set_title("Rock physics template: lithology and fluid-type discrimination")
    ax.legend(fontsize=7.5, framealpha=0.9, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_lambda_mu_rho(points, out_path):
    """points: dict[label] = (lambda_rho, mu_rho) in Pa*kg/m3 units."""
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    for label, (lr, mr) in points.items():
        ax.scatter(mr / 1e9, lr / 1e9, s=90, color=FLUID_COLORS.get(label, "black"),
                   edgecolor="black", linewidth=0.6, zorder=3, label=FLUID_LABELS.get(label, label))
    ax.set_xlabel(r"$\mu\rho$  (GPa $\cdot$ g/cc)")
    ax.set_ylabel(r"$\lambda\rho$  (GPa $\cdot$ g/cc)")
    ax.set_title(r"$\lambda\rho$ vs. $\mu\rho$: reservoir sand fluid discrimination")
    ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_vp_porosity(df, template_curves, reservoir_top, reservoir_base, out_path):
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    reservoir = df[(df["DEPTH_M"] >= reservoir_top) & (df["DEPTH_M"] <= reservoir_base)]
    background = df[(df["DEPTH_M"] < reservoir_top) | (df["DEPTH_M"] > reservoir_base)]

    ax.scatter(background["PHIE"], background["VP_BRINE_MPS"], s=4, c="#999999", alpha=0.35,
               label="Background shale/silt (brine)")
    ax.scatter(reservoir["PHIE"], reservoir["VP_BRINE_MPS"], s=10, c=FLUID_COLORS["BRINE"],
               alpha=0.7, label="Reservoir sand (brine, in-situ)")
    ax.scatter(reservoir["PHIE"], reservoir["VP_GAS_MPS"], s=10, c=FLUID_COLORS["GAS"],
               alpha=0.7, label="Reservoir sand (gas, substituted)")

    for key, style in (("brine_sand", dict(color=FLUID_COLORS["BRINE"])),
                        ("gas_sand", dict(color=FLUID_COLORS["GAS"]))):
        curve = template_curves[key]
        ax.plot(curve["PHIE"], curve["VP_MPS"], lw=1.6, ls="--", **style)

    ax.set_xlabel("Porosity (v/v)")
    ax.set_ylabel("Vp (m/s)")
    ax.set_title("Vp-porosity compaction and fluid-substitution trend")
    ax.legend(fontsize=7.5, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_angle_gathers(gathers, time_axis_s, angles_deg, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 6.2), sharey=True)
    t_ms = time_axis_s * 1000.0
    vmax = max(np.abs(g).max() for g in gathers.values())

    for ax, (label, gather) in zip(axes, gathers.items()):
        extent = [angles_deg[0], angles_deg[-1], t_ms[-1], t_ms[0]]
        im = ax.imshow(gather, aspect="auto", extent=extent, cmap="seismic",
                        vmin=-vmax, vmax=vmax)
        ax.set_xlabel("Incidence angle (deg)")
        ax.set_title(FLUID_LABELS[label])
    axes[0].set_ylabel("Two-way time (ms)")
    cbar = fig.colorbar(im, ax=axes, shrink=0.7, pad=0.02)
    cbar.set_label("Reflection amplitude")
    fig.suptitle("Synthetic pre-stack angle gathers, top- and base-reservoir reflections", y=0.99)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_near_mid_far(stacks, time_axis_s, window_ms, out_path):
    """stacks: dict[fluid_label][stack_label] -> 1D trace array."""
    t_ms = time_axis_s * 1000.0
    mask = (t_ms >= window_ms[0]) & (t_ms <= window_ms[1])

    stack_labels = list(next(iter(stacks.values())).keys())
    fig, axes = plt.subplots(1, len(stack_labels), figsize=(8.0, 6.0), sharey=True)
    for ax, stack_label in zip(axes, stack_labels):
        for fluid_label, stack_dict in stacks.items():
            ax.plot(stack_dict[stack_label][mask], t_ms[mask],
                    color=FLUID_COLORS[fluid_label], lw=1.4, label=FLUID_LABELS[fluid_label])
        ax.axvline(0, color="black", lw=0.5)
        ax.set_xlabel(stack_label)
        ax.invert_yaxis()
    axes[0].set_ylabel("Two-way time (ms)")
    axes[0].legend(fontsize=7, framealpha=0.9, loc="upper left")
    fig.suptitle("Near / mid / far angle stacks at the reservoir interval", y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_avo_crossplot(points, out_path):
    """points: dict[label] -> (intercept A, gradient B)."""
    fig, ax = plt.subplots(figsize=(6.2, 6.0))

    a_line = np.linspace(-0.4, 0.4, 50)
    ax.plot(a_line, -a_line, color="#999999", lw=1.0, ls=":", label="Background trend (B = -A)")
    ax.axhline(0, color="#cccccc", lw=0.8)
    ax.axvline(0, color="#cccccc", lw=0.8)

    ax.axvspan(-0.4, -0.02, color=FLUID_COLORS["GAS"], alpha=0.05)
    ax.axvspan(0.02, 0.4, color=FLUID_COLORS["BRINE"], alpha=0.05)

    for label, (a, b) in points.items():
        ax.scatter(a, b, s=110, color=FLUID_COLORS[label], edgecolor="black", linewidth=0.7,
                   zorder=3, label=f"{FLUID_LABELS[label]} (top of reservoir)")

    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-0.4, 0.4)
    ax.set_xlabel("Intercept, A")
    ax.set_ylabel("Gradient, B")
    ax.set_title("AVO intercept-gradient crossplot, top-of-reservoir reflection")
    ax.legend(fontsize=8, framealpha=0.9, loc="lower left")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
