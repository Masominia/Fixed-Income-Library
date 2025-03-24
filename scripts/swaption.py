"""Price a swaption under Hull-White: closed form, tree, and the Bermudan version.

python scripts/swaption.py expiry tenor [a] [sigma] [steps per year]
"""

import sys

from curvekit import engine

if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    expiry, tenor = float(sys.argv[1]), float(sys.argv[2])
    a, sigma = (float(x) for x in (sys.argv[3:5] if len(sys.argv) >= 5 else (0.05, 0.01)))
    steps = int(sys.argv[5]) if len(sys.argv) > 5 else 64
    h = engine.call("hw", expiry=expiry, tenor=tenor, a=a, sigma=sigma, steps=[steps])
    tree = h["rows"]["tree"][0]
    print(f"forward swap rate   {h['atm']:.4%}")
    print(f"payer, closed form  {h['payer']:.6f}")
    print(f"payer, tree         {tree[5]:.6f}")
    print(f"Bermudan payer      {tree[7]:.6f}")
