"""Python side of the C++ library: runs the binary and parses its output."""

from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "cpp" / "build" / "curvekit"
QUOTES = ROOT / "configs" / "quotes.txt"


def available() -> bool:
    return BINARY.exists()


def _value(text: str):
    try:
        return float(text)
    except ValueError:
        return text


def call(command: str, quotes: Path | None = None, **params) -> dict:
    """`call("hw", a=0.05, sigma=0.01)` -> {"payer": ..., "rows": {"tree": [[...], ...]}}"""
    args = [f"{k}={','.join(map(str, v)) if isinstance(v, (list, tuple)) else v}" for k, v in params.items()]
    out = subprocess.run(
        [str(BINARY), command, str(quotes or QUOTES), *args], capture_output=True, text=True, check=True
    ).stdout
    res: dict = {"rows": {}}
    for line in out.splitlines():
        f = line.split()
        if f[0] == "row":
            res["rows"].setdefault(f[1], []).append([_value(x) for x in f[2:]])
        else:
            res[f[0]] = _value(f[1])
    return res
