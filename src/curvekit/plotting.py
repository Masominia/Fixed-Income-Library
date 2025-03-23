"""Static figures for the report (matplotlib, PNG for the README)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Validated categorical order (fixed, never cycled) + recessive ink.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
NEUTRAL = "#8a8985"
CRITICAL = "#c4271a"

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "axes.labelcolor": INK2,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "legend.frameon": False,
        "lines.linewidth": 1.6,
        "figure.dpi": 150,
        "axes.axisbelow": True,
    }
)


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def overview(curves, ns_fit, conv, surface, path):
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.6))
    a = ax[0, 0]
    for i, (label, f) in enumerate(curves.items()):
        a.plot(f[:, 0], 100 * f[:, 1], color=SERIES[i], lw=1.4, label=label)
    a.set(title="Instantaneous forward rate by interpolation", xlabel="years", ylabel="%")
    a.legend()

    a = ax[0, 1]
    f = next(reversed(curves.values()))
    a.plot(f[1:, 0], 100 * f[1:, 2], color=SERIES[0], label="bootstrapped zero rate")
    a.plot(ns_fit[:, 0], 100 * ns_fit[:, 2], "o", color=SERIES[1], ms=4, label="Nelson-Siegel at the pillars")
    a.set(title="Zero curve and its three-factor fit", xlabel="years", ylabel="%")
    a.legend()

    a = ax[1, 0]
    for i, col in enumerate(["bond option", "caplet", "swaption"]):
        a.loglog(conv.index, conv[col].abs(), "o-", color=SERIES[i], ms=4, label=col)
    a.set(title="Hull-White tree against closed forms", xlabel="time steps per year", ylabel="relative error")
    a.legend()

    a = ax[1, 1]
    for i, tenor in enumerate(surface.columns):
        a.plot(surface.index, 100 * surface[tenor], "o-", color=SERIES[i], ms=4, label=f"{tenor}-year swap")
    a.set(title="Black volatility implied by Hull-White swaption prices", xlabel="option expiry, years", ylabel="%")
    a.legend()
    _save(fig, path)
