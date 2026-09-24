"""Deterministic constructed-array study; no field data, FDTD, or training."""

import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path

import numpy as np

from research_evaluation_contract import waveform_metrics, energy_metrics, interval, configuration_labels
from research_operator_contract import catalogue, apply_configuration, ConfigUnavailable


def metrics(y, s, mask=None):
    if mask is None:
        mask = np.ones(s.shape, dtype=bool)
    return waveform_metrics(y, s, mask, reference_kind="complete_clean", state="mixed", scope="complete")["metrics"]


def pulse(samples, centers, width):
    u = (np.arange(samples)[:, None] - np.asarray(centers)[None, :]) / width
    return np.where(np.abs(u) <= 4, (1 - 2 * u**2) * np.exp(-u**2), 0.)


def shift_zero(x, count):
    y = np.zeros_like(x)
    if count == 0:
        return x.copy()
    y[count:] = x[:-count]
    return y


def local_errors(y, s, width):
    # All windows fixed by reference geometry, before any candidate outputs.
    squared_error = np.sum((y - s)**2, axis=0)
    squared_reference = np.sum(s**2, axis=0)
    kernel = np.ones(width)
    numerator = np.convolve(squared_error, kernel, mode="valid")
    denominator = np.convolve(squared_reference, kernel, mode="valid")
    valid = denominator > 0
    return np.sqrt(numerator[valid] / denominator[valid])


def nrmse_bounds(yhat, shat, ey, es):
    if not np.isfinite([ey, es]).all() or min(ey, es) < 0:
        raise ValueError("nonnegative_finite_error_budgets_required")
    r, b = np.linalg.norm(yhat - shat), np.linalg.norm(shat)
    if b <= es:
        return {"lower": None, "upper": None, "reason": "reference_norm_not_separated_from_zero"}
    return {"lower": float(max(0., r - ey - es) / (b + es)),
            "upper": float((r + ey + es) / (b - es)), "reason": None}


def run():
    checks, damage, weak_rows, candidate_rows, labels = [], [], [], [], []

    def check(name, ok, **details):
        if not bool(ok):
            raise AssertionError(name)
        checks.append({"id": name, "passed": True, **details})

    for width in (3, 6, 12):
        s = pulse(256, np.full(256, 128.), width)
        variations = [("amplitude", a, a*s) for a in (1., .95, .9, .75, .5, 0., -1.)]
        variations += [("shift", k, shift_zero(s, k)) for k in (0, 1, 2, 4, 8, 16)]
        for count in (0, 1, 2, 4, 8, 16, 64, 256):
            y = s.copy()
            start = min(120, 256-count)
            y[:, start:start+count] = 0
            variations.append(("delete_traces", count, y))
        for kind, value, y in variations:
            m = metrics(y, s)
            row = {"width_samples": width, "damage": kind, "amount": value, **m}
            for window in (1, 8, 16, 32):
                row[f"max_local_D_w{window}"] = float(local_errors(y, s, window).max())
            damage.append(row)
            if kind == "amplitude":
                check(f"amplitude_w{width}_{value}", np.isclose(m["nrmse"], abs(value-1))
                      and np.isclose(m["shape_residual"], 0, atol=1e-12))
            if kind == "delete_traces":
                check(f"deletion_w{width}_{value}", np.isclose(m["nrmse"], np.sqrt(value/256))
                      and np.isclose(m["amplitude_error"], value/256)
                      and np.isclose(row["max_local_D_w1"], 1. if value else 0.))
    residual = max(abs(r["nrmse"]**2-r["amplitude_error"]**2-r["shape_residual"]**2) for r in damage)
    check("projection_identity", residual < 1e-12, maximum_residual=residual)
    edge = np.zeros((8, 2)); edge[-1] = 1
    check("shift_has_no_wraparound", not shift_zero(edge, 1).any())

    strong = pulse(256, np.full(256, 64.), 3)
    weak = pulse(256, np.full(256, 160.), 3)
    weak_mask = np.zeros(strong.shape, bool); weak_mask[144:177] = True
    for ratio in (.01, .1, .5):
        s = strong + ratio*weak
        for count in (1, 4, 16, 256):
            y = s.copy(); start = min(120, 256-count)
            y[:, start:start+count] -= ratio*weak[:, start:start+count]
            whole, event = metrics(y, s), metrics(y, ratio*weak, weak_mask)
            weak_rows.append({"weak_strong_ratio": ratio, "deleted_traces": count,
                              "whole_D": whole["nrmse"], "weak_event_D": event["nrmse"],
                              "weak_max_trace_D": float(local_errors(y*weak_mask, ratio*weak, 1).max())})
            expected = ratio*np.sqrt(count/256)/np.sqrt(1+ratio**2)
            check(f"weak_energy_formula_{ratio}_{count}", np.isclose(whole["nrmse"], expected))

    # Independent geometric perturbation check, not a probabilistic confidence interval.
    rng = np.random.default_rng(20260924)
    bounds_trials = 1000
    for _ in range(bounds_trials):
        shat = rng.normal(size=12); yhat = shat + rng.normal(size=12)
        es = float(rng.uniform(0, .8)*np.linalg.norm(shat))
        ey = float(rng.uniform(0, .8)*np.linalg.norm(shat))
        ds, dy = rng.normal(size=12), rng.normal(size=12)
        ds *= es*rng.uniform()/np.linalg.norm(ds)
        dy *= ey*rng.uniform()/np.linalg.norm(dy)
        bounds = nrmse_bounds(yhat, shat, ey, es)
        actual = np.linalg.norm(yhat+dy-shat-ds)/np.linalg.norm(shat+ds)
        if not bounds["lower"]-1e-12 <= actual <= bounds["upper"]+1e-12:
            raise AssertionError("uncertainty_enclosure")
    check("uncertainty_enclosure", True, trials=bounds_trials, seed=20260924)
    bounds = nrmse_bounds(np.array([1.3]), np.array([1.]), .02, .05)
    check("one_dimensional_bounds_attained", np.isclose(bounds["lower"], (1.28-1.05)/1.05)
          and np.isclose(bounds["upper"], (1.32-.95)/.95))
    check("uncertain_zero_reference_abstains", nrmse_bounds(np.ones(2), np.ones(2), 0, 2)["upper"] is None)
    uncertainty_examples = []
    for es in (0., .01, .03, .09, .5, 1.):
        b = nrmse_bounds(np.array([.92]), np.array([1.]), 0., es)
        uncertainty_examples.append({"nominal_D": .08, "reference_error_bound": es, **b,
                                     "fixture_limit": .1,
                                     "decision": "undetermined" if b["upper"] is None else
                                     "pass" if b["upper"] <= .1 else
                                     "fail" if b["lower"] > .1 else "undetermined"})

    # Gain diagnostic at the same sample, on known separate signal/noise components.
    gain_rows = []
    for gain in (1., 2., 4.):
        signal, noise = 2., .1
        ratio_before = signal**2/noise**2
        ratio_after = (gain*signal)**2/(gain*noise)**2
        gain_rows.append({"gain": gain, "noise_energy_ratio": gain**2,
                          "component_power_ratio_before": ratio_before, "component_power_ratio_after": ratio_after})
    check("pointwise_gain_keeps_component_snr", all(np.isclose(r["component_power_ratio_before"], r["component_power_ratio_after"]) for r in gain_rows))

    traces = np.arange(64)
    flat = pulse(192, np.full(64, 100.), 3)
    curved = pulse(192, 100+12*np.cos(2*np.pi*traces/64), 3)
    lateral = np.exp(-((traces-31.5)/6)**2); lateral -= lateral.mean()
    local = flat*lateral[None, :]
    early = 2*pulse(192, np.full(64, 30.), 3)
    weak_layer = .1*pulse(192, np.full(64, 145.), 3)
    scenes = [("flat_early", [flat], early), ("curved_early", [curved], early),
              ("local_early", [local], early), ("strong_weak_early", [flat, weak_layer], early),
              ("flat_clean", [flat], np.zeros_like(flat)), ("flat_overlap", [flat], 2*flat)]
    saved_arrays = {}
    for name, components, clutter in scenes:
        s = sum(components); x = s+clutter
        masks = []
        for component in components:
            occupied = np.flatnonzero(np.any(component != 0, axis=1))
            mask = np.zeros(x.shape, bool)
            mask[max(0, occupied[0]-8):min(192, occupied[-1]+9)] = True
            masks.append(mask)
        clutter_mask = np.zeros(x.shape, bool); clutter_mask[10:51] = True
        has_pure_clutter = name.endswith("early")
        records = []
        saved_arrays[f"{name}_reference"], saved_arrays[f"{name}_input"] = s, x
        saved_arrays[f"{name}_clutter"] = clutter
        for i, (component, mask) in enumerate(zip(components, masks)):
            saved_arrays[f"{name}_event{i}"] = component
            saved_arrays[f"{name}_mask{i}"] = mask
        for config in [c for c in catalogue() if c["end_gain"] == 1]:
            cid = config["id"]
            try:
                y = apply_configuration(x, cid)["output"]
            except ConfigUnavailable as exc:
                candidate_rows.append({"scene": name, "configuration": cid, "available": False, "reason": str(exc)})
                records.append({"id": cid, "available": False, "metrics": {}, "failure_reason": str(exc)})
                continue
            event_metrics = [metrics(y, s, mask) for mask in masks]
            d = max(m["nrmse"] for m in event_metrics)
            a = max(m["amplitude_error"] for m in event_metrics)
            r = energy_metrics(y, x, clutter_mask)["energy_ratio"] if has_pure_clutter else None
            candidate_rows.append({"scene": name, "configuration": cid, "available": True,
                                   "reason": None, "D_max_event": d, "A_max_event": a,
                                   "clutter_energy_ratio": r, "whole_D": metrics(y, s)["nrmse"]})
            for event_id, m in enumerate(event_metrics):
                candidate_rows[-1][f"event{event_id}_D"] = m["nrmse"]
                candidate_rows[-1][f"event{event_id}_A"] = m["amplitude_error"]
            saved_arrays[f"{name}_{cid}"] = y
            # Snapping only machine-size residues in exact constructed fixtures.
            losses = {"D": interval(0. if d < 1e-12 else d), "A": interval(0. if a < 1e-12 else a)}
            if r is not None:
                losses["R"] = interval(0. if r < 1e-12 else r)
            records.append({"id": cid, "available": True, "metrics": losses})
        objectives = ["D", "A", "R"] if has_pure_clutter else ["D", "A"]
        for limit in (None, .05, .1, .25, .5):
            label = configuration_labels(records, limits={"D": limit, "A": limit}, objectives=objectives)
            labels.append({"scene": name, "fixture_only_limit": limit, **label})
            if limit is None:
                check(f"null_limits_{name}", label["status"] == "partial_only" and label["unique_label"] is None)
    check("candidate_rows_complete", len(candidate_rows) == 42)
    check("unmodified_flat_mean_loss", np.isclose(next(r["D_max_event"] for r in candidate_rows if r["scene"] == "flat_clean" and r["configuration"] == "B1_G1_BG"), .25))
    # Same observed array admits mutually incompatible target/clutter assignments.
    observed = 3*flat
    check("overlap_nonidentifiability", np.array_equal(flat+2*flat, 3*flat+np.zeros_like(flat)))
    identifiability = {"observation": "X=3*flat", "world_A_target": "flat", "world_A_clutter": "2*flat",
                      "world_B_target": "3*flat", "world_B_clutter": "zero",
                      "same_input": True, "identity_D_world_A": metrics(observed, flat)["nrmse"],
                      "identity_D_world_B": metrics(observed, observed)["nrmse"]}
    result = {"schema": "constructed-damage-pilot/1", "evidence_level": "constructed_arrays_only",
              "field_data_used": False, "fdtd_executed": False, "training_executed": False,
              "physical_thresholds_frozen": False, "python": platform.python_version(), "numpy": np.__version__,
              "checks": checks, "damage": damage, "weak_events": weak_rows, "candidates": candidate_rows,
              "labels": labels, "uncertainty_examples": uncertainty_examples, "gain": gain_rows,
              "nonidentifiability": identifiability,
              "limitations": ["no physical time/depth units", "six designed scenes, no population inference",
                              "no random noise in candidate study", "fixture thresholds are not acceptance criteria",
                              "uncertainty bounds require independently justified input error budgets",
                              "candidate differences on constructed mixtures are not FDTD evidence"]}
    return result, saved_arrays


def plot(result, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with plt.rc_context({"font.size": 10, "svg.fonttype": "none"}):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
        rows = [r for r in result["damage"] if r["width_samples"] == 6 and r["damage"] == "delete_traces" and r["amount"] <= 16]
        for key, label, marker in [("nrmse", "Whole event", "o"), ("max_local_D_w16", "Max fixed 16-trace window", "s"), ("max_local_D_w1", "Max individual trace", "^")]:
            axes[0, 0].plot([r["amount"] for r in rows], [r[key] for r in rows], marker=marker, label=label)
        axes[0, 0].set(title="A. Local deletion can be diluted", xlabel="Deleted traces (of 256)", ylabel="NRMSE", ylim=(0, 1.05))
        axes[0, 0].legend(fontsize=8)
        for width, marker in zip((3, 6, 12), ("o", "s", "^")):
            rows = [r for r in result["damage"] if r["width_samples"] == width and r["damage"] == "shift"]
            axes[0, 1].plot([r["amount"] for r in rows], [r["nrmse"] for r in rows], marker=marker, label=f"Width parameter {width}")
        axes[0, 1].set(title="B. Shift sensitivity depends on pulse width", xlabel="Shift (samples; no depth calibration)", ylabel="NRMSE", ylim=(0, 1.9))
        axes[0, 1].legend(fontsize=8)
        for name, marker in (("flat_early", "o"), ("curved_early", "s"), ("local_early", "^")):
            rows = [r for r in result["candidates"] if r["scene"] == name and r["available"]]
            axes[1, 0].scatter([r["D_max_event"] for r in rows], [r["clutter_energy_ratio"] for r in rows], label=name, marker=marker, s=50)
        axes[1, 0].set(title="C. All available q=1 candidates; no ranking", xlabel="Max protected-event NRMSE", ylabel="Pure-clutter energy ratio (lower is less)", xlim=(-.03, 1.1), ylim=(-.03, 1.1))
        axes[1, 0].legend(fontsize=8)
        rows = result["uncertainty_examples"][:4]
        xx = np.arange(len(rows))
        # Draw interval endpoints directly: subtracting rounded 0.08 can create
        # a negative machine-size error-bar length for the exact point interval.
        axes[1, 1].vlines(xx, [r["lower"] for r in rows], [r["upper"] for r in rows], color="#0072B2")
        axes[1, 1].scatter(xx, [r["nominal_D"] for r in rows], color="#0072B2", marker="o")
        axes[1, 1].scatter(xx, [r["lower"] for r in rows], color="#0072B2", marker="_")
        axes[1, 1].scatter(xx, [r["upper"] for r in rows], color="#0072B2", marker="_")
        axes[1, 1].axhline(.1, linestyle="--", color="black", label="Illustrative limit 0.1")
        axes[1, 1].set_xticks(xx, [str(r["reference_error_bound"]) for r in rows])
        axes[1, 1].set(title="D. Deterministic bounds, not confidence intervals", xlabel="Reference norm error bound (reference norm = 1)", ylabel="True NRMSE enclosure", ylim=(0, .22))
        axes[1, 1].legend(fontsize=8)
        fig.suptitle("Constructed-array pilot: diagnostics only; no field data or FDTD", fontsize=14)
        fig.savefig(output / "damage_pilot.png", dpi=180, facecolor="white")
        fig.savefig(output / "damage_pilot.svg", facecolor="white")
        plt.close(fig)
    return matplotlib.__version__


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true", help="Optional figures; requires a working Matplotlib installation")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    result, arrays = run()
    for key in ("damage", "weak_events", "candidates", "uncertainty_examples", "gain"):
        rows = result[key]
        fields = list(dict.fromkeys(k for row in rows for k in row))
        with (args.output_dir / f"{key}.csv").open("x", newline="", encoding="utf-8-sig") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)
    np.savez_compressed(args.output_dir / "candidate_arrays.npz", **arrays)
    result["matplotlib"] = plot(result, args.output_dir) if args.plot else None
    result["figure_status"] = "generated" if args.plot else "not_requested_use_complete_csv_tables"
    source_files = [Path(__file__), Path(__file__).with_name("research_evaluation_contract.py"), Path(__file__).with_name("research_operator_contract.py"), Path(__file__).resolve().parents[1]/"docs/research/2026-09-24_damage_pilot_plan.md"]
    result["source_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
    result["artifact_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in args.output_dir.iterdir() if p.is_file()}
    (args.output_dir / "results.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"checks_passed": len(result["checks"]), "damage_rows": len(result["damage"]), "candidate_rows": len(result["candidates"]), "output": str(args.output_dir)}))


if __name__ == "__main__":
    main()
