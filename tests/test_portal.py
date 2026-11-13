import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(path):
    with open(path) as fh:
        return json.load(fh)


def test_portal_builds_from_the_result_tables(tmp_path):
    out = tmp_path / "site"
    subprocess.run([sys.executable, str(ROOT / "portal" / "build_portal.py"), "--out", str(out),
                    "--commit", "test"], check=True)
    index = load(out / "data" / "index.json")
    assert len(index) > 15000
    mpz = next(r for r in index if r["s"] == "MPZ")
    gene = load(out / "data" / "genes" / f"{mpz['k']}.json")[mpz["e"]]
    states = load(out / "data" / "states.json")
    assert len(states["states"]) == len(gene["spec"]) == 31
    assert len(states["deconv_states"]) == len(gene["deconv"]) == 16
    top = max(zip(gene["spec"], [s["state"] for s in states["states"]]))[1]
    assert top == "myelinating Schwann"
    assert gene["truth"]["category"] == "CMT" and gene["truth"]["tier"] == 1
    meta = load(out / "data" / "meta.json")
    assert meta["commit"] == "test" and meta["genes"] == len(index)
    assert (out / "index.html").exists() and (out / "tables" / "backtest_summary.tsv").exists()
