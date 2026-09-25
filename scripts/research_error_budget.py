"""Output error budget certificates for constructed-array operator chains.

Every budget number is computed from the actual arrays: correctly-rounded
math.fsum recomputation plus counted IEEE half-ulp slacks for arithmetic,
and Wedin/Davis-Kahan sin-Theta perturbation bounds for leading-SVD removal.
No epsilon-times-a-constant guesses, no field data, no FDTD, no training.
A budget certifies numerical identifiability of a direction only; never
physical accuracy. FDTD-sourced inputs are refused a budget.
"""

import math

import numpy as np

from research_operator_contract import catalogue, apply_configuration, ConfigUnavailable


U = np.finfo(np.float64).eps / 2.0  # IEEE 754 binary64 unit roundoff.


def _half_ulp(a):
    return np.abs(np.spacing(np.asarray(a, dtype=np.float64))) / 2.0


def _exact_construction(provenance, input_error_bound):
    if provenance in ("fdtd_unresolved", "fdtd_template"):
        return None, "fdtd_numerically_unresolved"
    if provenance != "constructed_exact":
        return None, "unknown_input_provenance"
    if not np.isfinite(input_error_bound) or input_error_bound < 0:
        raise ValueError("finite_nonnegative_input_error_bound_required")
    if input_error_bound > 0:
        return float(input_error_bound), "caller_declared_input_budget"
    return 0.0, "exact_construction"


def _mean_certificate(x, lam, y):
    """Bound ||y - (x - lam * exact_row_mean(x))||_F elementwise.

    Mirrors the operator path xs = x/scale, m = mean(xs, axis=1) * scale.
    fsum is correctly rounded; pairwise-vs-fsum difference is measured.
    """
    scale = float(np.max(np.abs(x)))
    n = x.shape[1]
    eb = np.zeros_like(y)
    if scale == 0:
        return eb
    xs = x / scale
    m_op = xs.mean(axis=1)
    m_err = np.empty(x.shape[0])
    for r in range(x.shape[0]):
        s_ref = math.fsum(xs[r].tolist())
        m_ref = s_ref / n
        delta_pair = abs(m_op[r] - m_ref)
        inbound = (delta_pair + float(_half_ulp(m_ref)) + float(_half_ulp(s_ref)) / n
                   + float(np.sum(_half_ulp(xs[r]))) / n)
        m_err[r] = scale * inbound + float(_half_ulp(m_op[r] * scale))
    lam_mul = 0.0 if math.frexp(lam)[0] == 0.5 else 1.0  # power-of-two lambda is exact
    for r in range(x.shape[0]):
        lm = lam * m_op[r] * scale
        eb[r, :] = lam * m_err[r] + lam_mul * float(_half_ulp(lm))
    eb += _half_ulp(y)
    return eb


def _gain_certificate(curve, y):
    if np.all(curve == 1.0):
        return np.zeros_like(y), True
    return _half_ulp(y).copy(), False


def _accurate_dot(a, b):
    """fsum-based dot with a counted rounding-error bound."""
    products = (np.asarray(a, dtype=np.float64) * np.asarray(b, dtype=np.float64)).tolist()
    value = math.fsum(products)
    error = float(np.sum(_half_ulp(np.asarray(products)))) + float(_half_ulp(value))
    return value, error


def _residual_bound(x, u, s, vt, i):
    """Upper bound on max(||X v_i - s_i u_i||, ||X^T u_i - s_i v_i||)."""
    n, m = x.shape
    worst = 0.0
    for mat, vec, unit in ((x, vt[i], u[:, i]), (x.T, u[:, i], vt[i])):
        count = mat.shape[0]
        w, werr = np.empty(count), np.empty(count)
        for j in range(count):
            row = mat[j]
            dot, derr = _accurate_dot(row, vec)
            p = s[i] * unit[j]
            w[j] = dot - p
            werr[j] = derr + float(_half_ulp(p)) + float(_half_ulp(w[j]))
        sq = (np.abs(w) + werr) ** 2
        total = math.fsum(sq.tolist())
        norm = math.sqrt(total + float(np.sum(_half_ulp(sq))) + float(_half_ulp(total)))
        worst = max(worst, norm)
    return worst

def svd_perturbation_bound(residuals, singular_values, k):
    """Per kept component: b_i = r_i + 2*sqrt(2)*sigma_i*sinTheta_i.

    sinTheta_i <= sqrt(2) r_i / gap_i (Wedin/Davis-Kahan sin-Theta theorem,
    standard matrix perturbation theory; not a GPR physical threshold).
    Missing when a gap is nonpositive or a bound is vacuous (sinTheta >= 1).
    """
    s = np.asarray(singular_values, dtype=np.float64)
    terms, details = [], []
    for i in range(k):
        r = float(residuals[i])
        left = s[i - 1] - s[i] if i > 0 else math.inf
        right = s[i] - s[i + 1]
        gap = float(min(left, right))
        if not gap > 0 or not np.isfinite(gap):
            return None, "svd_subspace_bound_vacuous"
        sin_theta = math.sqrt(2.0) * r / gap
        if not sin_theta < 1.0:
            return None, "svd_subspace_bound_vacuous"
        b = r + 2.0 * math.sqrt(2.0) * s[i] * sin_theta
        terms.append(b)
        details.append({"component": i, "residual": r, "gap": gap,
                        "sin_theta_bound": sin_theta, "output_perturbation_bound": b})
    return {"total": float(math.fsum(terms)), "components": details}, None


def _svd_certificate(x, k, y):
    """Arithmetic certificate plus subspace perturbation budget for SVD removal.

    Mathematical semantics: y_math = x - scale * R_math(xs_stored), where R_math
    is the exact rank-effk truncated-SVD component of the stored float64 xs.
    Terms: xs vs x/scale rounding, matmul-vs-fsum reconstruction, subtraction
    and rescaling roundings, and the sin-Theta perturbation of R_math.
    """
    scale = float(np.max(np.abs(x)))
    if scale == 0:
        return np.zeros_like(y), 0.0, [], None
    xs = x / scale
    u, s, vt = np.linalg.svd(xs, full_matrices=False)
    rank_floor = max(x.shape) * np.finfo(np.float64).eps
    numerical_rank = int(np.count_nonzero(s / s[0] > rank_floor))
    eff = min(k, numerical_rank)
    us = u[:, :eff] * s[:eff]
    removed = us @ vt[:eff, :]
    recon = np.zeros_like(y)
    if eff:
        product_arrays = []
        for q in range(eff):
            product_arrays.append(us[:, q][:, None] * vt[q][None, :])
        for index in np.ndindex(x.shape):
            terms = [float(pa[index]) for pa in product_arrays]
            ref = math.fsum(terms)
            slack = sum(float(_half_ulp(t)) for t in terms) + float(_half_ulp(ref))
            slack += sum(float(_half_ulp(us[index[0], q])) * abs(vt[q, index[1]])
                         for q in range(eff))
            recon[index] = abs(removed[index] - ref) + slack
        residuals = [_residual_bound(xs, u, s, vt, i) for i in range(eff)]
        perturb, missing = svd_perturbation_bound(residuals, s, eff)
        if missing:
            return None, None, [], missing
        perturb_total = scale * perturb["total"]
        perturb_details = perturb["components"]
    else:
        perturb_total, perturb_details = 0.0, []
    diff = xs - removed
    eb = scale * (recon + _half_ulp(xs) + _half_ulp(diff)) + _half_ulp(y)
    return eb, perturb_total, perturb_details, None


def configuration_output_budget(x, config_id, *, mask=None,
                                input_provenance="constructed_exact", input_error_bound=0.0):
    """Computed L2 output error budget certificate for one catalogue configuration.

    Returns available=False with a machine-readable reason for FDTD-sourced or
    unknown provenance, unavailable operators, and vacuous SVD bounds.
    The window budget is conservative: per-element arithmetic of the last step
    is restricted to the window; earlier-step and perturbation parts are
    full-array bounds, valid for any window since ||E[mask]|| <= ||E||.
    """
    e_in, certificate = _exact_construction(input_provenance, input_error_bound)
    if e_in is None:
        return {"available": False, "reason": certificate, "budget_full": None,
                "budget_window": None, "certificates": []}
    matches = [c for c in catalogue() if c["id"] == config_id]
    if not matches:
        raise ValueError("unknown_configuration")
    config = matches[0]
    try:
        result = apply_configuration(x, config_id)
    except ConfigUnavailable as exc:
        return {"available": False, "reason": f"operator_unavailable:{exc}",
                "budget_full": None, "budget_window": None, "certificates": []}
    certificates = ["input:" + certificate]
    carry = e_in
    eb_current, perturb = np.zeros_like(result["output"]), 0.0
    for step in result["steps"]:
        diagnostics = step["diagnostics"]
        kind = diagnostics["kind"]
        if carry > 0 and kind == "svd":
            return {"available": False, "reason": "svd_input_error_not_certified",
                    "budget_full": None, "budget_window": None, "certificates": certificates}
        if kind == "identity":
            eb_step, pert_step, sigma_max, note = np.zeros_like(x), 0.0, 1.0, "identity_exact_copy"
        elif kind == "mean":
            eb_step = _mean_certificate(step["input"], diagnostics["lambda"], step["output"])
            pert_step, sigma_max = 0.0, 1.0
            note = "mean_exact_recomputation_fsum"
        elif kind == "svd":
            eb_step, pert_step, details, missing = _svd_certificate(
                step["input"], diagnostics["k"], step["output"])
            if missing:
                return {"available": False, "reason": missing, "budget_full": None,
                        "budget_window": None, "certificates": certificates}
            sigma_max, note = 1.0, "svd_arith_fsum+wedin_davis_kahan_sin_theta"
            if details:
                certificates.append({"svd_components": details})
        elif kind == "fixed_shared_gain":
            curve = np.power(float(diagnostics["end_gain"]),
                             np.arange(x.shape[0]) / (x.shape[0] - 1))[:, None]
            eb_step, exact_one = _gain_certificate(curve, step["output"])
            pert_step = 0.0
            sigma_max = float(diagnostics["end_gain"])
            note = "gain_by_exact_one" if exact_one else "gain_counted_half_ulp"
        else:
            raise ValueError(f"unknown_step_kind:{kind}")
        carry = sigma_max * (carry + float(np.sqrt(np.sum(eb_current**2))) + perturb)
        eb_current, perturb = eb_step, pert_step
        certificates.append(note)
    window_last = float(np.sqrt(np.sum(eb_current[mask]**2))) if mask is not None \
        else float(np.sqrt(np.sum(eb_current**2)))
    full = carry + float(np.sqrt(np.sum(eb_current**2))) + perturb
    budget = {"available": True, "reason": None, "budget_full": full,
              "budget_window": carry + perturb + window_last,
              "output": result["output"], "certificates": certificates,
              "config": config}
    return budget


def difference_budget(budget_a, budget_b, norm_a, norm_b, norm_difference):
    """Triangle-inequality budget e_diff = e_a + e_b plus amplification diagnostic."""
    values = (budget_a, budget_b, norm_a, norm_b, norm_difference)
    if not all(np.isfinite(v) and v >= 0 for v in values):
        raise ValueError("finite_nonnegative_budgets_and_norms_required")
    return {"budget": budget_a + budget_b, "certificate": "triangle_inequality",
            "cancellation_amplification": (None if norm_difference == 0 else
                                           (norm_a + norm_b) / norm_difference)}
