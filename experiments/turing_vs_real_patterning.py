"""Study A vs real embryonic patterning (verdict #3): can lateral-inhibition
(Turing) dynamics with MEASURED Drosophila blastoderm parameters produce
pair-rule stripe-scale periodic patterns?

Measured anchors (Gregor et al. 2007, Cell 130:141-152; full text fetched and
quoted in the paper): Bcd-GFP cortical diffusion D = 0.30 +/- 0.09 um^2/s
(direct pattern photobleaching, 21 curves, 4 embryos), transport-model estimate
0.37 +/- 0.05, unfertilized eggs 0.35 +/- 0.12; gradient forms in ~1 hr;
nuclear cycles 9-12 min; embryo ~500 um. eve stripe biology: 7 stripes driven
by 5 discrete enhancers (Bothma et al. 2014 / eLife reviews).

Turing wavelength (activator-inhibitor, geometric-mean form, O(1) kinetic
constant noted):  lambda* = 2 pi (D_a D_h / (mu_a mu_h))^{1/4}
Stripe-scale target: ~50 um (7 stripes over the ~65% EL trunk of a ~500 um egg).
"""
import json
import numpy as np

D_A_MEAS = 0.30   # um^2/s, Gregor 2007 direct measurement
def lambda_star(d_a, d_h, mu):
    return 2 * np.pi * (d_a * d_h / mu**2) ** 0.25

grid = []
for lifetime_min in (10, 30, 60, 180, 540):   # 540 min = SDD-consistent lifetime from l=100um, D=0.3
    mu = 1.0 / (lifetime_min * 60.0)
    for ratio in (2, 5, 10, 20):              # D_h / D_a
        lam = lambda_star(D_A_MEAS, D_A_MEAS * ratio, mu)
        grid.append({"lifetime_min": lifetime_min, "Dh_over_Da": ratio,
                     "lambda_star_um": round(float(lam), 1)})

# required activator diffusion for stripe-scale pattern at 1-hr lifetime
mu_1h = 1.0 / 3600.0
target = 50.0 / (2 * np.pi)
d_eff_needed = target**4 * mu_1h**2     # (D_a D_h)
d_a_needed = np.sqrt(d_eff_needed / 10) # with D_h/D_a = 10

out = {
 "measured_anchors": {"D_bcd_cortical_um2_s": 0.30, "D_sd": 0.09,
    "gradient_formation_hr": 1.0, "embryo_length_um": 500,
    "source": "Gregor et al. 2007 Cell 130:141-152"},
 "stripe_scale_target_um": 50.0,
 "wavelength_grid": grid,
 "min_lambda_star_um": min(g["lambda_star_um"] for g in grid),
 "max_lambda_star_um": max(g["lambda_star_um"] for g in grid),
 "d_activator_needed_for_50um_at_1h_lifetime": round(float(d_a_needed), 5),
 "fold_below_measured": round(D_A_MEAS / float(d_a_needed), 1),
 "reading": "measured morphogen mobility puts any Turing wavelength at embryo scale, not stripe scale",
}
json.dump(out, open("results/turing_vs_real_patterning.json", "w"), indent=1)
for g in grid:
    print(g)
print("D_a needed for 50 um stripes (1h lifetime):", out["d_activator_needed_for_50um_at_1h_lifetime"],
      "um^2/s =", out["fold_below_measured"], "x below measured")
