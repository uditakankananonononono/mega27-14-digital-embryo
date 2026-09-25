"""Symbolic verification of paper derivations + agreement with committed numerics."""
import json, os, sys
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
import verify_math as vm

def test_gm_steady_state():
    sol, f, g = vm.gm_steady_state()
    a, h, rho, mu_a, mu_h = sp.symbols('a h rho mu_a mu_h')
    assert sp.simplify(f.subs(sol)) == 0 and sp.simplify(g.subs(sol)) == 0
    assert sp.simplify(sol[a] - mu_h/mu_a) == 0
    assert sp.simplify(sol[h] - mu_h*rho/mu_a**2) == 0

def test_gm_dispersion_matches_committed_numerics():
    d = vm.gm_dispersion()
    mu_a, mu_h, Da, Dh = sp.symbols('mu_a mu_h Da Dh')
    # parameter set used in experiments/morphogen demo (results/results.json -> A_turing.analytic)
    sub = {mu_a: sp.Rational('0.06'), mu_h: sp.Rational('0.12'), Da: sp.Rational('0.005'), Dh: sp.Rational('0.2')}
    tr = float(d['trace'].subs(sub)); det = float(d['det'].subs(sub)); B = float(d['B'].subs(sub))
    here = os.path.dirname(__file__)
    committed = json.load(open(os.path.join(here, '..', 'results', 'results.json')))['A_turing']['analytic']
    assert abs(tr - committed['trace']) < 1e-9
    assert abs(det - committed['det']) < 1e-9
    assert abs(B - committed['B']) < 1e-9
    band = sorted(float(r.subs(sub)) for r in d['k2_roots'])
    assert abs(band[0] - committed['k2_min']) < 1e-6
    assert abs(band[1] - committed['k2_max']) < 1e-6
    assert committed['k2_min'] < committed['sim_dominant_k2'] < committed['k2_max']

def test_sdd_solution():
    s = vm.sdd_steady_state()
    assert s['ode_ok'] and s['flux_ok']

def test_threshold_readout():
    t = vm.threshold_readout_fisher()
    lam, sig, cT = sp.symbols('lam sig cT')
    assert sp.simplify(t['sigma_x'] - lam*sig/cT) == 0
