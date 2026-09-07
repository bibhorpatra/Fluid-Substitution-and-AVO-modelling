# Fluid Substitution & Pre-Stack Seismic Modeling

A Gassmann fluid-substitution and AVO rock-physics workflow for a shaly-sand clastic reservoir, implemented from first principles in Python.

This project builds a synthetic well through a high-porosity Tertiary clastic reservoir, performs Gassmann fluid substitution to predict how the rock's elastic properties change when the brine-filled sand is instead charged with oil or gas, calibrates a Vp/Vs-vs-acoustic-impedance rock physics template against log and core data, and generates pre-stack synthetic angle gathers to evaluate the resulting AVO response.

**No proprietary or measured field dataset is used anywhere in this repository.** Every log curve is generated from standard, cited empirical and petrophysical relations — Athy's law compaction, a quartz-clay Voigt-Reuss-Hill mineral mixture, the Krief et al. (1990) dry-frame model, Batzle & Wang (1992) fluid properties, and Gassmann (1951) fluid substitution — calibrated to literature-typical parameter ranges for a shallow, unconsolidated North Sea-type Tertiary clastic setting. It is a synthetic "type well," not a claimed real-field dataset.

## Workflow

1. **Synthetic type-well construction** — a 200 m (2000-2200 m TVD) log section with a 16 m clean reservoir sand (2140-2156 m), built from an Athy's-law porosity-depth compaction trend and a Vshale-controlled quartz-clay mineral mixture.
2. **Gassmann fluid substitution** — the in-situ, fully brine-saturated reservoir sand is substituted with oil (32° API dead oil) and gas (0.65 specific gravity) at in-situ pressure/temperature, holding the rock frame fixed.
3. **Rock physics template** — a porosity-swept forward model traces brine/oil/gas-sand and shale trends in Vp/Vs-vs-acoustic-impedance space, calibrated against the well log samples and a set of synthetic core-plug measurements.
4. **Pre-stack AVO modelling** — Shuey (1985) three-term reflectivity, a 30 Hz Ricker wavelet, and synthetic angle gathers (0-35°) for each fluid case, with near/mid/far stacks and an intercept-gradient AVO classification.

## Results

Reservoir sand elastic properties by fluid case (2142-2154 m clean-sand average), against a shale cap of Vp = 2584 m/s, Vs = 1144 m/s, ρ = 2.167 g/cc:

| Fluid case | Vp (m/s) | Vs (m/s) | ρ (g/cc) | AI (m/s·g/cc) | Vp/Vs | AVO class (top of reservoir) |
|---|---|---|---|---|---|---|
| Brine (in-situ) | 2929 | 1637 | 2.093 | 6131 | 1.790 | **Class I** — A = +0.045, B = −0.280 |
| Oil (substituted) | 2804 | 1660 | 2.036 | 5709 | 1.690 | **Class II/IIp** — A = +0.010, B = −0.324 |
| Gas (substituted) | 2681 | 1765 | 1.800 | 4826 | 1.519 | **Class III** — A = −0.074, B = −0.390 |

Gas substitution lowers Vp by ~8%, density by ~14%, and Vp/Vs by ~15% relative to the in-situ brine case — the classic gas signature, since pore-fluid softening affects the bulk modulus far more than it affects the rock frame, while Vs is comparatively fluid-insensitive.

The oil case is the interesting one: its impedance sits almost exactly on top of the shale cap's, driving the AVO intercept to within a few hundredths of zero — a genuinely subtle anomaly that a single near-offset stack would likely miss, and that only becomes diagnostic once the full-angle gradient is examined. The gas case, by contrast, gives an unambiguous Class III bright-spot response: strong negative intercept reinforced by a strong negative gradient, i.e. amplitude growing markedly with offset.

### Log tracks and lithology model

![Log tracks](figures/01_log_tracks.png)

### Depth trends: in-situ brine vs. Gassmann-substituted oil/gas

![Depth trends](figures/02_depth_trends.png)

### Rock physics template (Vp/Vs vs. acoustic impedance)

![Rock physics template](figures/03_rock_physics_template.png)

### Lambda-Rho / Mu-Rho fluid discrimination

Under Gassmann's assumptions the shear modulus — and therefore μρ, since `μρ = ρVs² = μ` identically — is fluid-independent, so the brine/oil/gas points fall at essentially the same μρ while spreading widely in λρ (Goodway et al., 1997). Shale, with different mineralogy and frame properties, separates cleanly on the μρ axis.

![Lambda-Mu-Rho](figures/04_lambda_mu_rho.png)

### Vp-porosity compaction and fluid trend

![Vp-porosity](figures/05_vp_porosity.png)

### Synthetic pre-stack angle gathers

Top- and base-of-reservoir reflections, 0-35° incidence, 30 Hz Ricker wavelet.

![Angle gathers](figures/06_angle_gathers.png)

### Near / mid / far angle stacks

![Near mid far stacks](figures/07_near_mid_far_stacks.png)

### AVO intercept-gradient crossplot

![AVO crossplot](figures/08_avo_crossplot.png)

## Repository structure

```
├── run_analysis.py              # End-to-end driver: run this to reproduce everything
├── src/
│   ├── fluid_properties.py      # Batzle & Wang (1992) brine/oil/gas PVT properties
│   ├── rock_physics_model.py    # VRH mineral mixing, Krief dry-frame, Gassmann substitution
│   ├── synthetic_well.py        # Synthetic type-well log construction
│   ├── rock_physics_template.py # Porosity-swept RPT curves + synthetic core calibration
│   ├── avo_modeling.py          # Shuey reflectivity, Ricker wavelet, angle gathers, AVO classification
│   └── plotting.py              # All figure-generation routines
├── notebooks/
│   └── 01_fluid_substitution_avo_analysis.ipynb   # Narrated, executed end-to-end notebook
├── figures/                     # All output figures (PNG), regenerated by run_analysis.py
├── data/
│   ├── synthetic_well_logs.csv  # Full synthetic well log dataset
│   ├── synthetic_well.las       # Same data as a LAS 2.0 file
│   └── core_calibration_samples.csv
├── report/
│   └── Fluid_Substitution_AVO_Modeling_Report.pdf # Full written report
├── results_summary.json         # Machine-readable numeric results
└── requirements.txt
```

## Methodology notes and assumptions

- **Mineralogy**: a quartz (K = 36.6 GPa, μ = 45.0 GPa, ρ = 2650 kg/m³) / clay (K = 20.9 GPa, μ = 6.85 GPa, ρ = 2580 kg/m³) Voigt-Reuss-Hill mixture, weighted by Vshale.
- **Dry-frame model**: Krief et al. (1990), `K_dry/K_min = μ_dry/μ_min = (1-φ)^(A/(1-φ))`, with A = 3.
- **In-situ conditions**: hydrostatic pressure gradient (10.0 MPa/km) and a 30 °C/km geothermal gradient from an 8 °C surface temperature — typical of a normally-pressured North Sea Tertiary section.
- **Brine properties**: Batzle & Wang (1992), via the open-source `bruges` library, at 35,000 ppm NaCl salinity.
- **Oil properties**: a simplified dead-oil (gas-free) model — exact API-gravity density definition, a linear temperature/pressure density correction, and a representative reservoir-condition acoustic velocity — used because no measured PVT report (GOR, bubble point) exists for this synthetic case. See the docstring in `src/fluid_properties.py` for the full derivation and the reasoning behind this simplification.
- **Gas properties**: an ideal-gas approximation (density via the ideal gas law, adiabatic bulk modulus `K = γP` with γ = 1.25), consistent with the same approximation used internally by the brine/gas density functions this project depends on.
- **AVO modelling**: the Shuey (1985) three-term approximation is accurate to roughly 30-35°, which bounds the angle range modelled here (0-35°). A constant overburden velocity (2400 m/s) is used only for the depth-to-time conversion that positions the two reflectors on the synthetic trace — it does not affect the reflectivity modelling itself.
- **Limitations**: this is a 1D, weak-contrast, uniform (Reuss-average) saturation model. It does not capture patchy saturation effects, anisotropy, attenuation/dispersion, or a full nonlinear Zoeppritz solution beyond ~35°. It is intended as a rock-physics and AVO-modelling demonstration, not a field-ready interpretation.

## Reproducing this project

```bash
pip install -r requirements.txt
python run_analysis.py
```

This regenerates every figure in `figures/`, the data products in `data/`, and `results_summary.json` from scratch (fixed random seeds, so results are deterministic). The notebook in `notebooks/` mirrors the same workflow with narrative commentary; re-execute it with:

```bash
jupyter nbconvert --to notebook --execute --inplace notebooks/01_fluid_substitution_avo_analysis.ipynb
```

## References

- Aki, K., and Richards, P. G., 1980, *Quantitative Seismology*. W. H. Freeman.
- Athy, L. F., 1930, Density, porosity, and compaction of sedimentary rocks. *AAPG Bulletin*, 14, 1-24.
- Batzle, M., and Wang, Z., 1992, Seismic properties of pore fluids. *Geophysics*, 57(11), 1396-1408.
- Castagna, J. P., and Swan, H. W., 1997, Principles of AVO crossplotting. *The Leading Edge*, 16(4), 337-342.
- Gassmann, F., 1951, Über die Elastizität poröser Medien. *Vierteljahrsschrift der Naturforschenden Gesellschaft in Zürich*, 96, 1-23.
- Goodway, B., Chen, T., and Downton, J., 1997, Improved AVO fluid detection and lithology discrimination using Lambda-Mu-Rho (LMR). *67th SEG Annual Meeting, Expanded Abstracts*.
- Hill, R., 1952, The elastic behaviour of a crystalline aggregate. *Proc. Phys. Soc. London*, A65, 349-354.
- Krief, M., Garat, J., Stellingwerff, J., and Ventre, J., 1990, A petrophysical interpretation using the velocities of P and S waves. *The Log Analyst*, 31, 355-369.
- Rutherford, S. R., and Williams, R. H., 1989, Amplitude-versus-offset variations in gas sands. *Geophysics*, 54(6), 680-688.
- Shuey, R. T., 1985, A simplification of the Zoeppritz equations. *Geophysics*, 50(4), 609-614.

## License

MIT — see [LICENSE](LICENSE).
