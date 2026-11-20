"""Independent numerical checks. Run from any directory with NumPy and pandas."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SAVED = ROOT / "outputs/saved"
EMP = SAVED / "empirical"
METHODS = ["ss", "kernel", "quadratic", "changepoint"]
OBJECTS = ["x", "v", "a", "average_a"]
RECORDS = ["NASA", "NOAA", "HadCRU", "Berk", "ERA5", "Mean"]


def check_close(actual, expected, atol=1e-8):
    np.testing.assert_allclose(actual, expected, rtol=0, atol=atol)


def refit_benchmarks(t, y, z):
    mask = t >= 1970
    t, y, z = t[mask], y[mask], z[mask]
    u = t - t.mean()
    base = np.column_stack([np.ones(len(t)), u, z])
    quadratic = np.column_stack([base, u**2])
    b = np.linalg.lstsq(quadratic, y, rcond=None)[0]
    q = [b[0]+b[1]*u[-1]+b[-1]*u[-1]**2, 10*(b[1]+2*b[-1]*u[-1]), 200*b[-1], 200*b[-1]]
    # Residualize each candidate hinge against the shared regressors.
    knots = t[60:-60]
    hinges = np.maximum(t[:, None] - knots, 0)
    basis = np.linalg.qr(base, mode="reduced")[0]
    residual = hinges - basis @ (basis.T @ hinges)
    k = np.argmax((residual.T @ y)**2 / np.sum(residual**2, axis=0))
    b = np.linalg.lstsq(np.column_stack([base, hinges[:, k]]), y, rcond=None)[0]
    c = [b[0]+b[1]*u[-1]+b[-1]*hinges[-1, k], 10*(b[1]+b[-1]), 0, 10*b[-1]*(knots[k]>t[-121])]
    return q, c, knots[k]


def local_quadratic(t, values):
    fitted = np.empty((len(t), 3, values.shape[1]))
    for i, ti in enumerate(t):
        u = (t-ti)/25
        keep = np.abs(u)<1
        u = u[keep]
        weights = np.sqrt((1-np.abs(u)**3)**3)
        x = np.column_stack([np.ones(len(u)), u, u*u])*weights[:,None]
        fitted[i] = np.linalg.lstsq(x, values[keep]*weights[:,None], rcond=None)[0]
    return fitted


def check_empirical():
    raw = pd.read_csv(ROOT / "data/Monthly_HadISST.csv")
    fit = pd.read_csv(EMP / "Empirical_Results.csv").query("Specification == 'monthly'")
    assert len(fit)==96 and not fit.duplicated(["Dataset","Method","Object"]).any()
    base = raw.t.between(1991,2021,inclusive="left")
    values = raw[RECORDS].sub(raw.loc[base,RECORDS].mean()) + .88
    mean = values[RECORDS[:4]].mean(axis=1, skipna=False)
    later = raw.t>=1940
    mean.loc[later] = values.loc[later,RECORDS[:5]].mean(axis=1,skipna=False)
    check_close(mean, values.Mean)
    cols = []
    for name,degree in [("NINO34",2),("volc2",1),("ssn",1)]:
        v = (raw[name]-raw[name].mean())/raw[name].std()
        lags = np.column_stack([v.shift(j) for j in range(13)])
        for power in range(degree+1):
            cols.append(lags @ ((np.arange(13)/12)**power))
    z = np.column_stack(cols)
    for record in RECORDS:
        keep = np.isfinite(values[record]) & np.isfinite(z).all(axis=1)
        t,y,reg = raw.t.to_numpy()[keep],values[record].to_numpy()[keep],z[keep]
        def target(method):
            return fit[(fit.Dataset==record)&(fit.Method==method)].set_index("Object").loc[OBJECTS]
        q, c, knot = refit_benchmarks(t,y,reg)
        check_close(q,target("quadratic").Estimate)
        check_close(c,target("changepoint").Estimate)
        check_close(knot,target("changepoint").BreakDate.iloc[0])
        stacked = np.column_stack([y,reg])
        smooth = local_quadratic(t,stacked)
        partial = stacked-smooth[:,0,:]
        beta = np.linalg.lstsq(partial[:,1:],partial[:,0],rcond=None)[0]
        coeff = smooth[:,:,0]-smooth[:,:,1:]@beta
        level,rate,accel = coeff[:,0],coeff[:,1]*10/25,coeff[:,2]*200/25**2
        check_close([level[-1],rate[-1],accel[-1],rate[-1]-rate[-121]],target("kernel").Estimate)
        for m in METHODS:
            sub = "validated_monthly" if m=="ss" else "reporting"
            path = pd.read_csv(EMP/sub/f"{record}_{m}_paths.csv")
            check_close(path.t.iloc[-1],t[-1])
            check_close([*path.iloc[-1][["x","v","a"]],path.v.iloc[-1]-path.v.iloc[-121]],target(m).Estimate)
        rate_table = fit[(fit.Dataset==record)&(fit.Object=="v")].set_index("Method").Estimate
        assert rate_table[["ss","kernel"]].min() > rate_table[["quadratic","changepoint"]].max()
    velocity = pd.read_csv(EMP/"reporting/VelocityModel_Comparison.csv").set_index("Dataset")
    assert len(velocity)==5
    for m in ["quadratic","changepoint"]:
        rate = fit[(fit.Method==m)&(fit.Object=="v")].set_index("Dataset").Estimate
        assert (velocity.AlternativeRate > rate.loc[velocity.index]).all()
    check_close(velocity.AlternativeRate.mean(),.389,atol=.0005)
    meanfit = fit[fit.Dataset=="Mean"].pivot(index="Method",columns="Object",values="Estimate")
    endpoint = fit[fit.Dataset=="Mean"].EndTime.unique().item()
    for threshold, expected in [(1.5,[2030.9,2028.7,2028.5,2028.4,2027.5,2027.5]),(2,[2046.7,2042.9,2039.7,2037.5,2037.5,2035.7])]:
        for mode, m, date in zip(["fit","fit","flat","fit","flat","fit"],["quadratic","changepoint","ss","ss","kernel","kernel"],expected):
            row=meanfit.loc[m]; a=0 if mode=="flat" else row.a; gap=threshold-row.x
            year=endpoint+20*gap/(row.v+np.sqrt(row.v**2+2*a*gap))
            check_close(year,date,atol=.05)
    print("PASS: 18 independent monthly fits, all 24 path endpoints and paired decade contrasts, lower-order smoother, 12 crossing dates.")


def check_simulations():
    sim=SAVED/"simulations"
    accuracy=pd.read_csv(sim/"accuracy/AccuracySummary.csv")
    assert len(accuracy)==80 and accuracy.Valid.eq(1000).all()
    wins=accuracy.loc[accuracy.groupby(["DGP","Estimand"]).RMSE.idxmin()]
    assert wins.Method.eq("ss").sum()==16 and wins.query("Estimand == 'v_T'").Method.eq("ss").all()
    for folder, methods in [("smooth_coverage",METHODS),("reference_coverage",["quadratic","changepoint"])]:
        directory=sim/folder
        if folder=="smooth_coverage":
            rows=pd.read_csv(directory/"CoverageReplications.csv")
            summary=pd.read_csv(directory/"CoverageSummary.csv")
            tails=["LowerMiss","UpperMiss"]
        else:
            rows=pd.concat([pd.read_csv(directory/f"{m}_truth.csv") for m in methods])
            summary=pd.read_csv(directory/"Summary.csv")
            tails=["LowerAboveTruth","UpperBelowTruth"]
        assert len(rows)==400*4*len(methods) and len(summary)==4*len(methods)
        for _,s in summary.iterrows():
            g=rows[(rows.Method==s.Method)&(rows.Estimand==s.Estimand)]
            assert set(g.Replication)==set(range(1,401)) and g.Success.all()
            if s.Method=="changepoint" and s.Estimand=="a_T":
                assert g.Low.isna().all() and g.High.isna().all() and np.isnan(s.Coverage)
                continue
            assert np.isfinite(g[["Truth","Low","High"]]).all().all() and (g.Low<=g.High).all()
            left, right=g.Truth<g.Low,g.Truth>g.High
            p=(~(left|right)).mean()
            check_close([p,np.sqrt(p*(1-p)/400),left.mean(),right.mean()],s[["Coverage","CoverageMCSE",*tails]].to_numpy(dtype=float))
    print("PASS: accuracy ranking and coverage, sampling errors and both tails reconstructed from 9,600 estimator/estimand/replication rows.")


if __name__ == "__main__":
    check_empirical()
    check_simulations()
