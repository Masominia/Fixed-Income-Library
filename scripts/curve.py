"""Bootstrap a quote file and print the curve.

python scripts/curve.py [quotes.txt] [linear_zero|log_linear|monotone_cubic]
"""

import sys

from curvekit import engine

if __name__ == "__main__":
    quotes = sys.argv[1] if len(sys.argv) > 1 else None
    mode = sys.argv[2] if len(sys.argv) > 2 else "monotone_cubic"
    c = engine.call("curve", quotes, mode=mode, step=1.0)
    print(f"{mode}: {int(c['passes'])} sweeps, largest repricing error {c['max_error']:.1e}")
    print(f"{'maturity':>9} {'zero':>8} {'discount':>10}")
    for t, z, df in c["rows"]["pillar"]:
        print(f"{t:>9.2f} {z:>8.4%} {df:>10.6f}")
