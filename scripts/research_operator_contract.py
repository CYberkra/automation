"""Dense-array research prototype, not a validated GPR processing library.

Only fixed shared gain, partial row-mean subtraction, and leading-SVD removal,
plus two out-of-catalogue mechanism references (local mean, robust PCA).
No field reader, FDTD execution, learning, ROI selection, or physical inversion.
"""

import numpy as np


VERSION = "research-operators/0.1"
VERSION_02 = "research-operators/0.2"
GAP_RTOL = 1e-8  # Proposed numerical guard, not a physical acceptance threshold.

# v0.2 adds the RPCA background kind as a lambda axis (0.5/1.0/2.0 x standard
# PCP lambda0 = 1/sqrt(max(shape))), end_gain 1, BG order only (frozen as
# operator_catalogue_v0.2.json, 2026-09-29; v0.1 remains the default so all
# v0.1-era checks and frozen hashes stay valid).
RPCA_LAM_FACTORS = (0.5, 1.0, 2.0)


class ConfigUnavailable(ValueError):
    """A valid catalogue entry cannot be applied to this particular input."""


def catalogue(version="0.1"):
    backgrounds = [("identity", None), ("mean", 0.25), ("mean", 0.5),
                   ("mean", 1.0), ("svd", 1), ("svd", 2), ("svd", 3)]
    if version == "0.2":
        backgrounds = backgrounds + [("rpca", f) for f in RPCA_LAM_FACTORS]
    elif version != "0.1":
        raise ValueError("unknown_catalogue_version")
    rows = []
    for index, (kind, parameter) in enumerate(backgrounds):
        for end_gain in (1, 2, 4):
            if kind == "rpca" and end_gain > 1:
                continue  # RPCA joins with end_gain 1 / BG order only in v0.2
            for order in (["BG", "GB"] if kind == "svd" and end_gain > 1 else ["BG"]):
                rows.append({"id": f"B{index}_G{end_gain}_{order}", "background": kind,
                             "parameter": parameter, "end_gain": end_gain, "order": order})
    return rows


def _background(x, kind, parameter, gap_rtol):
    if kind == "identity":
        return x.copy(), {"kind": kind}
    if kind == "mean":
        # Scale before averaging to avoid overflow from summing finite values.
        scale = float(np.max(np.abs(x)))
        mean = np.zeros((x.shape[0], 1)) if scale == 0 else (x / scale).mean(axis=1, keepdims=True) * scale
        return x - parameter * mean, {"kind": kind, "lambda": parameter}
    if kind == "rpca":
        # parameter is the lambda factor relative to the standard PCP choice
        # lambda0 = 1/sqrt(max(shape)); see rpca_control for the solver.
        lam = float(parameter) / (max(x.shape) ** 0.5)
        try:
            y, diag = rpca_control(x, lam)
        except (FloatingPointError, np.linalg.LinAlgError) as exc:
            raise ConfigUnavailable("numeric_failure") from exc
        return y, {"kind": kind, "lam_factor": float(parameter),
                   "lam_over_lambda0": diag["lam_over_standard"],
                   "status": diag["status"], "rank_L": diag["rank_L"],
                   "iterations": diag["iterations"],
                   "rel_residual": diag["rel_residual"]}
    k = int(parameter)
    if k >= min(x.shape):
        raise ConfigUnavailable("svd_rank_must_be_less_than_minimum_dimension")
    scale = float(np.max(np.abs(x)))
    if scale == 0:
        return x.copy(), {"kind": kind, "k": k, "status": "zero_input", "effective_k": 0}
    u, singular, vt = np.linalg.svd(x / scale, full_matrices=False)
    relative = singular / singular[0]
    rank_floor = max(x.shape) * np.finfo(np.float64).eps
    numerical_rank = int(np.count_nonzero(relative > rank_floor))
    effective_k = min(k, numerical_rank)
    gap = float(relative[k - 1] - relative[k])
    # Ties entirely below the numerical rank floor carry no resolved component.
    if k < numerical_rank and gap <= gap_rtol:
        raise ConfigUnavailable("svd_cutoff_gap_unresolved")
    removed = (u[:, :effective_k] * singular[:effective_k]) @ vt[:effective_k, :]
    return (x / scale - removed) * scale, {
        "kind": kind, "k": k, "effective_k": effective_k,
        "relative_singular_values": relative.tolist(), "cutoff_gap_relative_to_s1": gap,
        "gap_rtol": gap_rtol, "numerical_rank_floor": rank_floor,
        "status": "rank_redundant" if k > numerical_rank else "ok",
    }


def apply_configuration(x, config_id, *, mask=None, domain="time_real", gap_rtol=GAP_RTOL,
                        catalogue_version="0.1"):
    """Return Y and actual step arrays. BG pre_gain differs from GB audit_before_gain.

    Axes are [sample, trace]. The caller defines the entire input window in advance.
    Invalid masked cells are preserved only for the all-identity configuration.
    Other configurations require a dense, fully valid window; no imputation occurs.
    """
    a = np.asarray(x)
    if domain != "time_real" or a.ndim != 2 or min(a.shape) < 2:
        raise ValueError("require_time_real_2d_sample_trace_with_each_dimension_at_least_2")
    if a.dtype.kind != "f" or not np.isfinite(a).all():
        raise ValueError("require_finite_real_floating_input")
    if not np.isfinite(gap_rtol) or not 0 < gap_rtol < 1:
        raise ValueError("invalid_gap_guard")
    matches = [c for c in catalogue(catalogue_version) if c["id"] == config_id]
    if not matches:
        raise ValueError("unknown_configuration")
    config = matches[0]
    if mask is not None:
        mask = np.asarray(mask)
        if mask.dtype.kind != "b" or mask.shape != a.shape:
            raise ValueError("mask_must_be_boolean_and_match_shape")
        if not mask.all() and (config["background"] != "identity" or config["end_gain"] != 1):
            raise ConfigUnavailable("nonidentity_requires_dense_valid_window")
    a = a.astype(np.float64, copy=True)
    curve = np.power(float(config["end_gain"]), np.arange(a.shape[0]) / (a.shape[0] - 1))[:, None]
    steps = []
    current = a.copy()
    pre_gain = None
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        try:
            for op in config["order"]:
                before = current.copy()
                if op == "B":
                    current, diagnostics = _background(current, config["background"], config["parameter"], gap_rtol)
                else:
                    current = curve * current
                    diagnostics = {"kind": "fixed_shared_gain", "end_gain": config["end_gain"], "axis": 0}
                if not np.isfinite(current).all():
                    raise ConfigUnavailable("nonfinite_output")
                steps.append({"operator": op, "input": before, "output": current.copy(),
                              "diagnostics": diagnostics})
                if op == "B" and config["order"] == "BG":
                    pre_gain = current.copy()
        except (FloatingPointError, np.linalg.LinAlgError) as exc:
            raise ConfigUnavailable("numeric_failure") from exc
    return {"version": VERSION_02 if catalogue_version == "0.2" else VERSION,
            "config": config, "output": current,
            "pre_gain": pre_gain, "audit_before_gain": current / curve,
            "gain_curve": curve[:, 0].copy(), "steps": steps,
            "mask": None if mask is None else mask.copy(), "input_dtype": str(np.asarray(x).dtype),
            "compute_dtype": "float64", "domain": domain, "axes": ["sample", "trace"]}


def rpca_control(x, lam, max_iter=500, tol=1e-7):
    """Robust PCA / principal component pursuit, out-of-catalogue mechanism reference.

    Decomposes the [sample, trace] array as x = L + S with L low-rank (trace-
    invariant background) and S sparse (events/anomalies) via the inexact
    augmented-Lagrange method of Lin, Chen & Ma (2009); deterministic, float64.
    Returns (x - L, diagnostics): the sparse output is the background-suppressed
    radargram. lam is the sparsity weight (standard choice 1/sqrt(max(shape))).
    """
    a = np.array(x, dtype=np.float64, copy=True)
    if a.ndim != 2 or min(a.shape) < 2 or not np.isfinite(a).all():
        raise ValueError("invalid_control_array")
    if not np.isfinite(lam) or lam <= 0:
        raise ValueError("invalid_lambda")
    scale = float(np.max(np.abs(a)))
    if scale == 0:
        return a.copy(), {"kind": "rpca", "lam": lam, "status": "zero_input",
                          "rank_L": 0, "iterations": 0, "rel_residual": 0.0}
    m = a / scale
    fro = float(np.linalg.norm(m))
    norm2 = float(np.linalg.svd(m, compute_uv=False)[0])
    dual = max(norm2, float(np.linalg.norm(m, np.inf)) / lam)
    y = m / dual
    L = np.zeros_like(m)
    S = np.zeros_like(m)
    mu = 1.25 / norm2
    mu_bar = mu * 1e7
    rho = 1.5
    rel = None
    it = 0
    for it in range(1, max_iter + 1):
        u, sv, vt = np.linalg.svd(m - S + y / mu, full_matrices=False)
        sv_thr = np.maximum(sv - 1.0 / mu, 0.0)
        L = (u * sv_thr) @ vt
        t2 = m - L + y / mu
        S = np.sign(t2) * np.maximum(np.abs(t2) - lam / mu, 0.0)
        r = m - L - S
        y = y + mu * r
        mu = min(mu * rho, mu_bar)
        rel = float(np.linalg.norm(r) / fro)
        if rel < tol:
            break
    out = m - L
    return out * scale, {"kind": "rpca", "lam": lam,
                         "lam_over_standard": round(lam * (max(a.shape) ** 0.5), 6),
                         "status": "converged" if rel is not None and rel < tol else "max_iter",
                         "rank_L": int(np.count_nonzero(np.linalg.svd(L, compute_uv=False) > 1e-10 * max(L.shape))),
                         "iterations": it, "rel_residual": rel,
                         "sparse_fraction": float(np.count_nonzero(S) / S.size)}


def local_mean_control(x, width, strength=1.0):
    """Independent spatially adaptive reference, outside the 27-entry catalogue.

    Truncated and count-normalized edges; no implicit global mean subtraction.
    Mechanism helper only. Public API/mask/provenance support is not implemented.
    """
    a = np.array(x, dtype=np.float64, copy=True)
    if a.ndim != 2 or not np.isfinite(a).all():
        raise ValueError("invalid_control_array")
    if isinstance(width, bool) or not isinstance(width, int) or width < 3 or width % 2 == 0 or width > a.shape[1]:
        raise ValueError("require_odd_integer_window_between_3_and_trace_count")
    if not np.isfinite(strength) or not 0 <= strength <= 1:
        raise ValueError("invalid_strength")
    result = a.copy()
    half = width // 2
    for j in range(a.shape[1]):
        result[:, j] -= strength * a[:, max(0, j - half):min(a.shape[1], j + half + 1)].mean(axis=1)
    return result
