import pytest

from curvekit import engine
from curvekit.pipeline import run


@pytest.mark.skipif(not engine.available(), reason="C++ library not built")
def test_fast_pipeline(tmp_path):
    m = run(out=tmp_path, fast=True)
    assert not [c for c in m["checks"] if c["status"] == "FAIL"]
    assert (tmp_path / "report.md").exists() and (tmp_path / "figures" / "overview.png").exists()
    assert m["bootstrap"]["monotone cubic"]["largest repricing error"] < 1e-12
