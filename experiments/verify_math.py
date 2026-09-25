"""Symbolically verify the derivations used in the paper (sympy).
Run standalone or via tests/test_math_verify.py."""
import sympy as sp

def gm_steady_state():
    # Gierer-Meinhardt (rescaled): da/dt = rho*a^2/h - mu_a a + Da nabla^2 a
    #                               dh/dt = rho*a^2   - mu_h h + Dh nabla^2 h
    a, h, rho, mu_a, mu_h = sp.symbols('a h rho mu_a mu_h', positive=True)
    f = rho*a**2/h - mu_a*a
    g = rho*a**2 - mu_h*h
    # homogeneous steady state: f=0 => h = rho*a/mu_a ; g=0 => h = rho*a**2/mu_h
    sol = sp.solve([f, g], [a, h], dict=True)
    sol = [s for s in sol if s[a] != 0][0]
    return sol, f, g

def gm_dispersion():
    sol, f, g = gm_steady_state()
    a, h, rho, mu_a, mu_h, Da, Dh, k2 = sp.symbols('a h rho mu_a mu_h Da Dh k2', positive=True)
    J = sp.Matrix([[sp.diff(f, a), sp.diff(f, h)],[sp.diff(g, a), sp.diff(g, h)]]).subs(sol)
    J = J.applyfunc(sp.simplify)
    M = J - sp.diag(Da*k2, Dh*k2).as_explicit() if False else sp.Matrix([[Da*k2,0],[0,Dh*k2]])
    char = sp.expand(M.charpoly().as_expr())  # sigma^2 - tr sigma + det
    tr = sp.simplify(J.trace()); det = sp.simplify(J.det())
    # instability band: det(M) < 0  =>  Da*Dh*k2^2 - (Dh*J11 + Da*J22)k2 + det(J) < 0
    B = sp.simplify(Dh*J[0,0] + Da*J[1,1])
    k2_roots = sp.solve(sp.Eq(Da*Dh*k2**2 - B*k2 + det, 0), k2)
    k2max = sp.simplify(B/(2*Da*Dh))
    return dict(steady=sol, J=J, trace=tr, det=det, B=B, k2_roots=k2_roots, k2max=k2max, charpoly=char)

def sdd_steady_state():
    # SDD: production at x=0 with flux J0, diffusion D, degradation gamma:
    #   -D c'' + gamma c = 0 on x>0,  -D c'(0) = J0  =>  c = (J0/sqrt(D gamma)) exp(-x sqrt(gamma/D))
    x, D, gam, J0 = sp.symbols('x D gam J0', positive=True)
    lam = sp.sqrt(D/gam)
    A = J0/(D*sp.sqrt(gam/D))
    c = A*sp.exp(-x/lam)
    ode_ok = sp.simplify(-D*sp.diff(c, x, 2) + gam*c) == 0
    flux_ok = sp.simplify(-D*sp.diff(c, x).subs(x, 0) - J0) == 0
    return dict(lam=lam, c=c, ode_ok=ode_ok, flux_ok=flux_ok)

def threshold_readout_fisher():
    # c(x) = A exp(-x/lam); estimator: read c at x with Gaussian noise sigma_c (abs),
    # position estimate from threshold crossing of c = cT:  x_hat = lam ln(A/cT)
    # error propagation: dx/dc = lam/c  =>  sigma_x = lam * sigma_c / cT
    A, lam, cT, sig = sp.symbols('A lam cT sig', positive=True)
    x_hat = lam*sp.log(A/cT)
    dxdA = sp.simplify(sp.diff(x_hat, A))          # lam/A -> rel. error of A maps to abs. error of x
    sigma_x = sp.simplify(lam * sig / cT)
    return dict(x_hat=x_hat, dxdA=dxdA, sigma_x=sigma_x)

if __name__ == '__main__':
    d = gm_dispersion()
    print('steady state:', d['steady'])
    print('trace:', d['trace'], ' det:', d['det'])
    print('B =', sp.simplify(d['B']))
    print('k2 band roots:', [sp.nsimplify(r) for r in d['k2_roots']])
    print('k2max (peak of dispersion):', d['k2max'])
    s = sdd_steady_state(); print('SDD ok:', s['ode_ok'], s['flux_ok'], 'lambda =', s['lam'])
    t = threshold_readout_fisher(); print('threshold readout:', t)
