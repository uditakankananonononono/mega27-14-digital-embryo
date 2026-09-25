"""statsmodels OLS cross-check of the dosage-coupled length-scale effect (DCLS).
Complements the Spearman analysis in results/bcd_windows.json with a parametric fit
on the same [0, 0.5] EL window."""
import json, multiprocessing as mp, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from digitalembryo.bcd import fit_single, load_live_gradients
from bcd_sdd_test import DOSAGE

def _work(rec):
    x, y = rec["x"], rec["intensity"]
    m = (x >= 0.0) & (x < 0.5)
    if m.sum() < 15:
        return None
    s = fit_single(x[m], y[m])
    return {"line": rec["line"], "lam": s["lam"], "dose": DOSAGE.get(rec["line"])}

def main():
    import statsmodels.api as sm
    recs = load_live_gradients()
    with mp.Pool(2) as pool:
        rows = [r for r in pool.map(_work, recs) if r and r["dose"]]
    dose = np.array([r["dose"] for r in rows], float)
    lam = np.array([r["lam"] for r in rows], float)
    X = sm.add_constant(dose)
    m = sm.OLS(lam, X).fit(cov_type='HC3')
    out = {
        "model": "OLS lambda ~ bcd_dose, HC3 robust SE, fit window [0,0.5]EL",
        "n_gradients": int(m.nobs),
        "slope_per_dose": float(m.params[1]),
        "slope_ci95": [float(v) for v in m.conf_int()[1]],
        "p_value_slope": float(m.pvalues[1]),
        "r2": float(m.rsquared),
        "crosscheck": "Spearman rho=0.356, p=7.9e-16 on same window (results/bcd_windows.json)",
    }
    json.dump(out, open("results/dcls_ols.json", "w"), indent=1)
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
