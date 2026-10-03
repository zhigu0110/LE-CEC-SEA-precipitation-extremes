# Selected LE-CEC analysis functions

This limited English code subset contains calculations that match the manuscript's physical decomposition, mechanism relationships, and historical observational evaluation. It is not the full LE-CEC constraint or a complete figure reproduction pipeline.

| Module                    | Scope                                                                               |
| ------------------------- | ----------------------------------------------------------------------------------- |
| physical_decomposition.py | Historical-baseline thermodynamic, dynamic, and residual terms in Equation 2        |
| mechanism_statistics.py   | Paired inter-SMILE Pearson correlations and regional area means for Figures 3 and 4 |
| observation_metrics.py    | Weighted spatial COR, RMSE, RSD, and PBIAS for historical evaluation                |

The functions were adapted from scripts 04, 05, 06, and 09. File loading and plotting were separated from the scientific calculations. Input shape and finite-support checks were added; correlation significance uses the actual paired sample count. No headline result values or synthetic-data fallback are included.

Supply verified ensemble-mean fields for the manuscript's 12 SMILEs, historical 1985-2014 and future 2070-2099 under SSP5-8.5/RCP8.5. Use identical grids, units, model pairing, and explicitly defined spatial support. Saturation humidity must be prepared and validated independently; the questionable sea-level-pressure preprocessing is excluded. Unit conversions for figure display belong to the caller.

The residual is retained for closure checks even though the manuscript's four mechanism panels display only thermo_term, dynamic_term, beta_bar, and sens_q.

Dependencies are listed in requirements.txt. Array-level synthetic checks passed, including decomposition closure, finite-pair correlations, weighted means, and observation metrics. These checks do not verify published numerical results. Published numerical results and final figures have not been independently reproduced by this package. This repository is distributed under the MIT license; see LICENSE.

## Paths and execution

All data and output paths are configured in config.py:

```python
data_path = "./data/"
output_path = "./results/"
```

Run from the package directory after preparing the inputs described in data/README.md:

```bash
python -m pip install -r requirements.txt
python run_analysis.py
```

The runner exports decomposition fields, four mechanism correlation fields with valid sample counts, a regional model summary, and optional observation metrics. Input/output archives are ignored by Git. No personal absolute paths or font paths are included. These portable array inputs require verified preprocessing; they do not replace the original NetCDF archive or establish full manuscript reproduction.
