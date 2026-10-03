"""Run the selected manuscript calculations using explicit array inputs."""
from pathlib import Path
import csv
import json
import numpy as np
from config import data_path, output_path, input_file, observation_file
from physical_decomposition import decompose
from mechanism_statistics import paired_correlation, area_mean
from observation_metrics import spatial_metrics


def main():
    source = Path(data_path) / input_file
    target = Path(output_path)
    required = ("rx_his", "rx_future", "q_his", "q_future", "tas_his", "tas_future",
                "latitude", "longitude", "mask", "model_names")
    with np.load(source, allow_pickle=False) as data:
        missing = sorted(set(required) - set(data.files))
        if missing:
            raise ValueError(f"Missing input arrays: {missing}")
        arrays = {key: data[key] for key in required}
    names = arrays["model_names"]
    if names.ndim != 1 or len(names) != 12 or len(set(names.tolist())) != 12:
        raise ValueError("Exactly 12 distinct SMILE names are required.")
    expected = (12, len(arrays["latitude"]), len(arrays["longitude"]))
    for key in required[:6]:
        if arrays[key].shape != expected:
            raise ValueError(f"{key} must have shape {expected}.")
    mask = arrays["mask"]
    if mask.dtype != bool or mask.shape != expected[1:]:
        raise ValueError("The common regional mask must be boolean and match the grid.")
    fields = decompose(*(arrays[key] for key in required[:6]))
    target.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(target / "physical_decomposition.npz", **fields,
                        latitude=arrays["latitude"], longitude=arrays["longitude"],
                        mask=mask, model_names=names)
    keys = ("thermo_term", "dynamic_term", "beta_bar", "sens_q")
    correlations = {}
    for key in keys:
        relation = paired_correlation(fields["rx1day_his"], fields[key], axis=0)
        correlations.update({key + "_" + metric: np.where(mask, values,
                             0 if metric == "n" else np.nan)
                             for metric, values in relation.items()})
    np.savez_compressed(target / "mechanism_correlations.npz", **correlations,
                        latitude=arrays["latitude"], longitude=arrays["longitude"], mask=mask)
    # Use common finite support across all models and plotted variables.
    regional_mask = mask.copy()
    for field in fields.values():
        regional_mask &= np.all(np.isfinite(field), axis=0)
    if not regional_mask.any():
        raise ValueError("No common finite regional support across the 12 SMILEs.")
    summary = {key: area_mean(value, arrays["latitude"], regional_mask)
               for key, value in fields.items()}
    with (target / "mechanism_model_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["model"] + list(summary))
        writer.writerows([[name] + [float(v[i]) for v in summary.values()]
                         for i, name in enumerate(names)])
    obs_source = Path(data_path) / observation_file
    if obs_source.exists():
        with np.load(obs_source, allow_pickle=False) as obs:
            climatology = obs["rx1day_climatology"]
            if not (np.array_equal(obs["latitude"], arrays["latitude"])
                    and np.array_equal(obs["longitude"], arrays["longitude"])):
                raise ValueError("Observation coordinates must exactly match the model grid.")
        weights = np.broadcast_to(np.cos(np.deg2rad(arrays["latitude"]))[:, None], mask.shape)
        eval_mask = regional_mask & np.isfinite(climatology)
        metrics = {str(name): spatial_metrics(arrays["rx_his"][i], climatology,
                                              weights, eval_mask)
                   for i, name in enumerate(names)}
        (target / "observation_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(f"Selected calculations saved to {target}. No combined constraint was computed.")


if __name__ == "__main__":
    main()
