"""All experiments, checks and the report."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import engine, plotting
from .engine import ROOT
from .report import md_table

MODES = {
    "linear_zero": "linear zero rates",
    "log_linear": "log-linear discount factors",
    "monotone_cubic": "monotone cubic",
}


def _check(checks, area, name, ok, detail, warn=False):
    checks.append(
        {"area": area, "check": name, "status": "PASS" if ok else ("WARN" if warn else "FAIL"), "detail": detail}
    )


def run(cfg_path=None, out=None, fast: bool = False) -> dict:
    cfg = yaml.safe_load(open(cfg_path or ROOT / "configs" / "default.yaml"))
    out = Path(out) if out else ROOT / "reports"
    (out / "figures").mkdir(parents=True, exist_ok=True)
    if not engine.available():
        raise SystemExit("C++ library not built: run `make build` first")
    checks: list[dict] = []
    hw = cfg["hull_white"]

    # 1. bootstrap under three interpolations
    curves, rows = {}, {}
    for mode, label in MODES.items():
        c = engine.call("curve", mode=mode, step=0.05)
        f = np.array(c["rows"]["forward"])
        curves[label] = f
        rows[label] = {
            "sweeps": int(c["passes"]),
            "largest repricing error": c["max_error"],
            "largest forward jump, bp": 1e4 * np.abs(np.diff(f[:, 1])).max(),
            "lowest forward, %": 100 * f[:, 1].min(),
            "highest forward, %": 100 * f[:, 1].max(),
        }
        pillars = pd.DataFrame(c["rows"]["pillar"], columns=["maturity, years", "zero rate", "discount factor"])
    boot = pd.DataFrame(rows).T
    boot["sweeps"] = boot["sweeps"].astype(int)
    _check(
        checks,
        "curve",
        "every quote repriced under every interpolation",
        boot["largest repricing error"].max() < 1e-12,
        f"largest error {boot['largest repricing error'].max():.1e}",
    )
    _check(
        checks,
        "curve",
        "the cubic removes the forward jumps",
        boot.loc[MODES["monotone_cubic"], "largest forward jump, bp"]
        < 0.25 * boot.loc[MODES["log_linear"], "largest forward jump, bp"],
        f"{boot.loc[MODES['monotone_cubic'], 'largest forward jump, bp']:.1f} bp between points 0.05 years apart, against {boot.loc[MODES['log_linear'], 'largest forward jump, bp']:.1f} (log-linear) and {boot.loc[MODES['linear_zero'], 'largest forward jump, bp']:.1f} (linear zero)",
    )
    _check(
        checks,
        "curve",
        "forward rates stay positive",
        boot["lowest forward, %"].min() > 0,
        f"lowest {boot['lowest forward, %'].min():.2f}%",
    )

    # 2. bond analytics
    b = engine.call("bond", mode="linear_zero", shifts=cfg["shifts_bp"], **cfg["bond"])
    krd = pd.DataFrame(b["rows"]["krd"], columns=["pillar, years", "key-rate duration"])
    sh = pd.DataFrame(
        b["rows"]["shift"], columns=["yield shift, bp", "exact price", "duration only", "duration and convexity"]
    ).set_index("yield shift, bp")
    sh["error, duration only"] = sh["duration only"] - sh["exact price"]
    sh["error, with convexity"] = sh["duration and convexity"] - sh["exact price"]
    krd_gap = abs(krd["key-rate duration"].sum() - b["effective_duration"])
    _check(
        checks,
        "bond",
        "key-rate durations add up to the effective duration",
        krd_gap < 1e-5,
        f"sum {krd['key-rate duration'].sum():.6f} against {b['effective_duration']:.6f}",
    )
    big = sh.loc[[s for s in sh.index if abs(s) == 100]]
    _check(
        checks,
        "bond",
        "convexity removes most of the error at 100 bp",
        (big["error, with convexity"].abs() < 0.05 * big["error, duration only"].abs()).all(),
        f"{big['error, duration only'].abs().max():.3f} down to {big['error, with convexity'].abs().max():.4f} per 100 of face",
    )
    _check(
        checks,
        "bond",
        "duration alone always understates the price",
        bool((sh["error, duration only"] < 0).all()),
        "negative error for every shift",
    )

    # 3. Nelson-Siegel
    ns = engine.call("nelson")
    nsf = np.array(ns["rows"]["fit"])
    _check(
        checks,
        "curve",
        "Nelson-Siegel fits the zero rates within 10 bp",
        ns["max_error"] < 1e-3,
        f"largest error {1e4 * ns['max_error']:.1f} bp, root mean square {1e4 * ns['rmse']:.1f} bp",
    )

    # 4. Hull-White: closed forms against the tree
    steps = cfg["tree_steps"][: 5 if fast else None]
    h = engine.call("hw", steps=steps, **hw, **cfg["swaption"])
    tr = pd.DataFrame(
        h["rows"]["tree"],
        columns=[
            "steps per year",
            "half width",
            "curve fit error",
            "bond option",
            "caplet",
            "payer",
            "receiver",
            "Bermudan payer",
        ],
    )
    tr["steps per year"] = tr["steps per year"].astype(int)
    conv = pd.DataFrame(
        {
            "steps per year": tr["steps per year"],
            "nodes across": 2 * tr["half width"].astype(int) + 1,
            "bond option": tr["bond option"] / h["bond_option"] - 1,
            "caplet": tr["caplet"] / h["caplet"] - 1,
            "swaption": tr["payer"] / h["payer"] - 1,
            "tree payer - receiver": tr["payer"] - tr["receiver"],
            "Bermudan / European": tr["Bermudan payer"] / tr["payer"],
        }
    ).set_index("steps per year")
    _check(
        checks,
        "Hull-White",
        "tree reproduces the curve",
        tr["curve fit error"].max() < 1e-12,
        f"largest relative error of a discount factor {tr['curve fit error'].max():.1e}",
    )
    _check(
        checks,
        "Hull-White",
        "closed-form payer and receiver satisfy parity",
        abs(h["parity_gap"]) < 1e-13,
        f"gap {h['parity_gap']:.1e}",
    )
    _check(
        checks,
        "Hull-White",
        "at-the-money parity holds on the tree",
        conv["tree payer - receiver"].abs().max() < 1e-12,
        f"largest gap {conv['tree payer - receiver'].abs().max():.1e}",
    )
    final = conv.iloc[-1][["bond option", "caplet", "swaption"]].abs().max()
    _check(
        checks,
        "Hull-White",
        "finest tree within 0.1% of the closed forms",
        final < (5e-3 if fast else 1e-3),
        f"largest relative error {final:.1e} at {conv.index[-1]} steps per year",
    )
    _check(
        checks,
        "Hull-White",
        "Bermudan worth more than European",
        conv["Bermudan / European"].min() > 1.05,
        f"ratio {conv['Bermudan / European'].iloc[-1]:.3f}",
    )

    # 5. implied volatility surface and calibration
    v = pd.DataFrame(
        engine.call("vols", **hw, **cfg["surface"])["rows"]["vol"], columns=["expiry", "tenor", "vol", "rate"]
    )
    surface = v.pivot(index="expiry", columns="tenor", values="vol")
    surface.index.name, surface.columns = "expiry \\ tenor", [f"{t:g}" for t in surface.columns]
    exact = engine.call("calibrate", **hw, **cfg["surface"])
    noisy = engine.call("calibrate", noise=cfg["calibration_noise"], **hw, **cfg["surface"])
    cal = pd.DataFrame(
        {
            "exact quotes": {
                "mean reversion": exact["fitted_a"],
                "sigma": exact["fitted_sigma"],
                "volatility error, bp": 1e4 * exact["rmse_vol"],
            },
            f"quotes perturbed by {cfg['calibration_noise']:.0%}": {
                "mean reversion": noisy["fitted_a"],
                "sigma": noisy["fitted_sigma"],
                "volatility error, bp": 1e4 * noisy["rmse_vol"],
            },
        }
    )
    _check(
        checks,
        "calibration",
        "parameters recovered from their own volatilities",
        abs(exact["fitted_a"] / hw["a"] - 1) < 1e-3 and abs(exact["fitted_sigma"] / hw["sigma"] - 1) < 1e-4,
        f"a {exact['fitted_a']:.6f}, sigma {exact['fitted_sigma']:.7f}",
    )
    _check(
        checks,
        "calibration",
        "stable under perturbed quotes",
        abs(noisy["fitted_sigma"] / hw["sigma"] - 1) < 0.05 and abs(noisy["fitted_a"] - hw["a"]) < 0.03,
        f"a {noisy['fitted_a']:.4f}, sigma {noisy['fitted_sigma']:.5f}",
    )
    _check(
        checks,
        "calibration",
        "volatility falls with expiry and tenor",
        bool((np.diff(surface.to_numpy(), axis=0) < 0).all() and (np.diff(surface.to_numpy(), axis=1) < 0).all()),
        f"{surface.iloc[0, 0]:.1%} to {surface.iloc[-1, -1]:.1%}",
    )

    # 6. short-rate models against Monte Carlo
    sr = cfg["short_rate"]
    s = engine.call("shortrate", **{k: (max(sr[k] // 10, 2000) if fast and k == "paths" else sr[k]) for k in sr})
    mc = pd.DataFrame(
        [["Vasicek", *r] for r in s["rows"]["vasicek"]] + [["Cox-Ingersoll-Ross", *r] for r in s["rows"]["cir"]],
        columns=["model", "maturity", "closed form", "Monte Carlo", "standard error"],
    )
    mc["difference, bp of price"] = 1e4 * (mc["Monte Carlo"] - mc["closed form"])
    mc["z"] = (mc["Monte Carlo"] - mc["closed form"]) / mc["standard error"]
    vz = mc.loc[mc["model"] == "Vasicek", "z"].abs().max()
    _check(checks, "short rate", "Vasicek formula within Monte Carlo error", vz < 4, f"largest |z| {vz:.2f}")
    cir_gap = mc.loc[mc["model"] != "Vasicek", "difference, bp of price"].abs().max()
    _check(
        checks,
        "short rate",
        "Cox-Ingersoll-Ross formula within 2 bp of an Euler simulation",
        cir_gap < 2,
        f"largest gap {cir_gap:.2f} bp of price",
    )

    plotting.overview(curves, nsf, conv, surface, out / "figures" / "overview.png")

    val = pd.DataFrame(checks)
    counts = val["status"].value_counts().to_dict()
    L = [
        "# CurveKit report",
        "",
        "Times in years. Quotes in `configs/quotes.txt`: two deposits, one forward rate agreement, nine swaps out to 30 years.",
        "",
        "## 1. Bootstrap",
        "",
        md_table(boot, 2),
        "",
        "Forward jump: largest change of the instantaneous forward rate between two points 0.05 years apart.",
        "",
        "Pillars (monotone cubic):",
        "",
        md_table(pillars, 6, index=False),
        "",
        f"## 2. Bond: {cfg['bond']['maturity']:g} years, {cfg['bond']['coupon']:.1%} coupon",
        "",
        f"Price {b['price']:.4f}, yield {b['yield']:.4%}. Macaulay duration {b['macaulay']:.4f}, modified {b['modified']:.4f}, convexity {b['convexity']:.3f}, value of a basis point {b['dv01']:.5f}. "
        f"Effective duration on the curve {b['effective_duration']:.4f}.",
        "",
        md_table(krd[krd["key-rate duration"] != 0], 4, index=False),
        "",
        "Repricing after a yield shift:",
        "",
        md_table(sh, 4),
        "",
        "## 3. Nelson-Siegel fit",
        "",
        f"Level {ns['level']:.4%}, slope {ns['slope']:.4%}, curvature {ns['curvature']:.4%}, decay time {ns['tau']:.2f} years. Root mean square error {1e4 * ns['rmse']:.2f} bp, largest {1e4 * ns['max_error']:.2f} bp.",
        "",
        f"## 4. Hull-White, a = {hw['a']}, sigma = {hw['sigma']}",
        "",
        f"Option expiring in {cfg['swaption']['expiry']:g} years on a {cfg['swaption']['tenor']:g}-year swap at the forward rate {h['atm']:.4%}. Closed forms: payer swaption {h['payer']:.6f}, "
        f"zero-coupon bond option {h['bond_option']:.6f}, caplet {h['caplet']:.6f} (Black volatility {h['black_caplet_vol']:.2%}).",
        "",
        "Relative error of the tree against the closed forms:",
        "",
        md_table(conv, 5),
        "",
        "## 5. Swaption volatilities implied by the model, and calibration back",
        "",
        md_table(surface, 4),
        "",
        md_table(cal, 6),
        "",
        "## 6. Short-rate models against Monte Carlo",
        "",
        md_table(mc, 6, index=False),
        "",
        "## Validation",
        "",
        f"{counts.get('PASS', 0)} PASS, {counts.get('WARN', 0)} WARN, {counts.get('FAIL', 0)} FAIL. See `validation.md`.",
        "",
    ]
    (out / "report.md").write_text("\n".join(L))
    (out / "validation.md").write_text("# Validation checks\n\n" + md_table(val, index=False) + "\n")
    metrics = {
        "bootstrap": boot.to_dict("index"),
        "bond": {k: v for k, v in b.items() if k != "rows"},
        "nelson_siegel": {k: v for k, v in ns.items() if k != "rows"},
        "hull_white": {k: v for k, v in h.items() if k != "rows"},
        "tree": {str(k): v for k, v in conv.to_dict("index").items()},
        "surface": surface.to_dict(),
        "calibration": cal.to_dict(),
        "short_rate": mc.to_dict("records"),
        "checks": checks,
    }
    (out / "metrics.json").write_text(json.dumps(metrics, indent=1, default=float))
    return metrics
