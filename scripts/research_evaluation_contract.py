"""Research metrics and partial configuration labels; no physical reference inference.

Callers supply independently defined reference provenance, masks, and constraints.
No filtering, field reading, learned score, or automatic interpretation is performed.
"""

import numpy as np


def _window(estimate, reference, mask):
    y, s, m = np.asarray(estimate), np.asarray(reference), np.asarray(mask)
    if y.shape != s.shape or m.shape != s.shape or s.ndim != 2 or m.dtype.kind != "b":
        raise ValueError("matching_2d_arrays_and_boolean_mask_required")
    if y.dtype.kind != "f" or s.dtype.kind != "f" or not np.isfinite(y).all() or not np.isfinite(s).all():
        raise ValueError("finite_real_floating_arrays_required")
    if not m.any():
        raise ValueError("empty_evaluation_window")
    return y[m].astype(float), s[m].astype(float)


def waveform_metrics(estimate, reference, mask, *, reference_kind, state, scope="event"):
    """No per-output scaling/alignment. Missing evidence is null, not zero loss.

    complete_clean: constructed full retained response; mixed allowed only for whole response.
    isolated_event: independently established absolute event reference.
    paired_contrast: evaluates a response difference only; never grants absolute preservation.
    """
    kinds = {"complete_clean", "isolated_event", "paired_contrast"}
    states = {"isolated", "mixed", "ambiguous", "absent", "numerically_unresolved"}
    if reference_kind not in kinds or state not in states or scope not in {"event", "complete", "contrast"}:
        raise ValueError("unknown_reference_contract")
    y, s = _window(estimate, reference, mask)
    eligible = state == "isolated" or (reference_kind == "complete_clean" and state == "mixed" and scope == "complete")
    contrast = reference_kind == "paired_contrast" and scope == "contrast" and state in {"isolated", "mixed"}
    if reference_kind == "paired_contrast":
        eligible = False
    if not (eligible or contrast):
        return {"available": False, "absolute_preservation_eligible": False,
                "reason": f"reference_state_or_scope:{state}:{reference_kind}:{scope}", "metrics": None}
    scale = float(np.max(np.abs(s)))
    if scale == 0:
        return {"available": False, "absolute_preservation_eligible": False,
                "reason": "zero_reference_use_negative_control", "metrics": None}
    # Numerical scaling shared by both arrays, not fitted to the candidate output.
    yn, sn = y / scale, s / scale
    ss = float(sn @ sn)
    alpha = float((yn @ sn) / ss)
    yy = float(yn @ yn)
    values = {"nrmse": float(np.linalg.norm(yn - sn) / np.sqrt(ss)),
              "amplitude_factor": alpha, "amplitude_error": abs(alpha - 1),
              "shape_residual": float(np.linalg.norm(yn - alpha * sn) / np.sqrt(ss)),
              "signed_cosine": None if yy == 0 else float((yn @ sn) / np.sqrt(yy * ss))}
    if not all(v is None or np.isfinite(v) for v in values.values()):
        raise ValueError("nonfinite_metric")
    return {"available": True, "absolute_preservation_eligible": eligible,
            "reason": "contrast_diagnostic_only" if contrast else None, "metrics": values}


def energy_metrics(output, baseline, mask):
    """Caller must establish pure-interference or absent-target status separately."""
    y, x = _window(output, baseline, mask)
    before, after = float(np.mean(x * x)), float(np.mean(y * y))
    if not np.isfinite(before) or not np.isfinite(after):
        raise ValueError("nonfinite_energy")
    return {"baseline_mean_square": before, "output_mean_square": after,
            "energy_ratio": None if before == 0 else after / before,
            "ratio_reason": "zero_baseline_report_absolute_energy" if before == 0 else None}


def interval(value, lower=None, upper=None):
    """Point values are permitted only when uncertainty is explicitly declared exact."""
    lo, hi = value if lower is None else lower, value if upper is None else upper
    if not all(np.isfinite(v) for v in (value, lo, hi)) or not 0 <= lo <= value <= hi:
        raise ValueError("finite_nonnegative_ordered_loss_interval_required")
    return {"value": float(value), "lower": float(lo), "upper": float(hi)}


def configuration_labels(records, *, limits, objectives):
    """Feasibility + robust Pareto labels; no arbitrary scalar weights.

    records: id, available, metrics={name: interval or None}, optional failure reason.
    limits: required metric -> maximum allowed loss; None means not calibrated.
    objectives: fixed, common, lower-is-better metric names for this label context.
    Failure dominates unknown in constraint evaluation. Incomparable is not equivalent.
    """
    if not records or not limits or not objectives or len(set(objectives)) != len(objectives):
        raise ValueError("nonempty_records_limits_and_unique_objectives_required")
    ids = [r["id"] for r in records]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate_configuration_id")
    if any(v is not None and (not np.isfinite(v) or v < 0) for v in limits.values()):
        raise ValueError("invalid_constraint_limit")
    assessed = []
    for record in records:
        if not isinstance(record["available"], bool):
            raise ValueError("boolean_availability_required")
        for item in record["metrics"].values():
            if item is not None:
                interval(item["value"], item["lower"], item["upper"])
        gates = {}
        for key, maximum in limits.items():
            value = record["metrics"].get(key)
            gates[key] = ("undetermined" if maximum is None or value is None else
                          "pass" if value["upper"] <= maximum else
                          "fail" if value["lower"] > maximum else "undetermined")
        status = ("unavailable" if not record["available"] else
                  "infeasible" if "fail" in gates.values() else
                  "undetermined" if "undetermined" in gates.values() else "feasible")
        assessed.append({"id": record["id"], "feasibility": status, "constraints": gates,
                         "metrics": record["metrics"], "failure_reason": record.get("failure_reason")})
    complete = [r for r in assessed if r["feasibility"] == "feasible"
                and all(r["metrics"].get(k) is not None for k in objectives)]
    edges = []
    for a in complete:
        for b in complete:
            if a["id"] == b["id"]:
                continue
            no_worse = all(a["metrics"][k]["upper"] <= b["metrics"][k]["lower"] for k in objectives)
            better = any(a["metrics"][k]["upper"] < b["metrics"][k]["lower"] for k in objectives)
            if no_worse and better:
                edges.append({"better": a["id"], "worse": b["id"]})
    dominated = {e["worse"] for e in edges}
    front = [r["id"] for r in complete if r["id"] not in dominated]
    partial = any(r["feasibility"] == "undetermined" or
                  (r["feasibility"] == "feasible" and r not in complete) for r in assessed)
    if partial:
        status = "partial_only"
    elif not complete:
        status = "no_admissible_candidate"
    else:
        status = "unique_admissible" if len(front) == 1 else "multiple_admissible"
    return {"status": status, "candidate_records": assessed, "nondominated_supported_set": front,
            "robust_preference_pairs": edges, "unique_label": front[0] if status == "unique_admissible" else None,
            "all_candidates_accounted_for": not partial,
            "objectives": list(objectives), "constraints": dict(limits),
            "missing_metrics_are_not_zero": True, "incomparable_is_not_equivalent": True}
