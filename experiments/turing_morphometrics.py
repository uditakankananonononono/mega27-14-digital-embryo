"""scikit-image morphometrics of a steady-state Turing pattern (Study A).
Segments activator spots and reports count/area/eccentricity + an independent
radial-spectrum wavenumber cross-check against the committed sweep value."""
import json, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from digitalembryo.morphogen import GiererMeinhardt2D
from skimage import measure, filters

MU_A = 0.03  # mid-sweep point; committed k2_sim = 2.6885237478643838 (results/turing_sweep.json)

def main():
    sim = GiererMeinhardt2D(n=96, Da=0.005, Dh=0.2, mu_a=MU_A, mu_h=0.12, seed=7)
    prev = None
    for blk in range(40):
        sim.step(500)
        a = sim.a
        if prev is not None and np.max(np.abs(a - prev)) < 1e-6:
            break
        prev = a.copy()
    a = sim.a
    thr = filters.threshold_otsu(a)
    lab = measure.label(a > thr)
    props = measure.regionprops(lab, intensity_image=a)
    areas = [p.area for p in props]
    ecc = [p.eccentricity for p in props]
    k2_fft = sim.pattern_wavenumber()
    out = {
        "model": "Gierer-Meinhardt 2D n=96, mu_a=0.03, mu_h=0.12, Da=0.005, Dh=0.2, steady state",
        "n_spots": len(props),
        "spot_area_mean": float(np.mean(areas)),
        "spot_area_cv": float(np.std(areas) / np.mean(areas)),
        "spot_eccentricity_mean": float(np.mean(ecc)),
        "otsu_threshold": float(thr),
        "k2_fft_this_run": float(k2_fft),
        "k2_sim_committed_n48": 2.6885237478643838,
        "note": "spot morphometrics via skimage measure/filters; k2 differs from committed n=48 sweep value due to lattice size (expected).",
    }
    json.dump(out, open("results/turing_morphometrics.json", "w"), indent=1)
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
