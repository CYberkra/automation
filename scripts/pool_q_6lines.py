# -*- coding: utf-8 -*-
"""Pool per-trace q across the 6 field lines -> merged percentiles."""
import glob
import os

import numpy as np

qs = []
for f in sorted(glob.glob(r"E:\automation_djh\artifacts_check\ringing_field\*_q.npy")):
    q = np.load(f)
    qs.append(q)
    print("%-12s n=%5d  p25=%.4f p50=%.4f p75=%.4f" % (
        os.path.basename(f), len(q),
        np.percentile(q, 25), np.percentile(q, 50), np.percentile(q, 75)))

p = np.concatenate(qs)
print("pooled      n=%5d  p10=%.4f p25=%.4f p50=%.4f p75=%.4f p90=%.4f" % (
    (len(p),) + tuple(np.percentile(p, x) for x in (10, 25, 50, 75, 90))))
