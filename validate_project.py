"""Integrity and synthetic helper tests only; never reports PBMC biological results."""
import ast
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parent
results = []
def check(name, action):
    action()
    results.append({"check": name, "passed": True})

source = (ROOT / "analysis.py").read_text(encoding="utf-8")
tree = ast.parse(source)
helper_names = {"qc_decisions", "sparse_fingerprint", "partition_comparison", "propose_identity"}
helper_nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in helper_names]
namespace = dict(np=np, pd=pd, sp=sp, hashlib=hashlib, adjusted_rand_score=adjusted_rand_score)
exec(compile(ast.Module(body=helper_nodes, type_ignores=[]), "extracted_helpers", "exec"), namespace)
cfg = dict(min_genes=200, max_genes_exclusive=2500, max_mt_pct_exclusive=5.0,
           stable_fraction=.90, split_fraction=.10)

def test_python_and_notebook():
    for path in ROOT.glob("*.py"):
        compile(path.read_text(encoding="utf-8"), str(path), "exec")
    nb = json.loads((ROOT / "Single_Cell_Immune_Profiling.ipynb").read_text(encoding="utf-8"))
    assert nb["nbformat"] == 4 and nb["nbformat_minor"] == 5
    assert nb["metadata"]["kernelspec"]["name"] == "python3"
    assert len({cell["id"] for cell in nb["cells"]}) == len(nb["cells"])
    code = []
    for i, cell in enumerate(nb["cells"]):
        assert cell["cell_type"] in {"markdown", "code"}
        assert isinstance(cell["metadata"], dict) and isinstance(cell["source"], list)
        assert all(isinstance(line, str) for line in cell["source"])
        if cell["cell_type"] == "code":
            text = "".join(cell["source"])
            compile(text, f"cell_{i}", "exec")
            assert cell["execution_count"] is None and cell["outputs"] == []
            code.append(text)
    assert len(code) == 11
    assert ast.dump(ast.parse("\n".join(code))) == ast.dump(ast.parse(source))
    requirements_node = next(node for node in tree.body if isinstance(node, ast.Assign)
                             and any(isinstance(t, ast.Name) and t.id == "REQUIREMENTS" for t in node.targets))
    pinned = ast.literal_eval(requirements_node.value)
    file_pins = [line for line in (ROOT / "requirements.txt").read_text().splitlines()
                 if line and not line.startswith("#")]
    assert pinned == file_pins

def test_qc_boundaries():
    obs = pd.DataFrame({"n_genes_by_counts": [199, 200, 2499, 2500, 900, 900, 2500],
                        "pct_counts_mt": [1, 1, 4.999, 1, 5, 4.999, 5]})
    got = namespace["qc_decisions"](obs, cfg)
    assert got.qc_pass.tolist() == [False, True, True, False, False, True, False]
    assert got.iloc[-1].fail_high_genes and got.iloc[-1].fail_mt
    assert int((~got.qc_pass).sum()) == 4
    assert int(got.drop(columns="qc_pass").sum().sum()) == 5  # overlapping failures

def test_fingerprint_preservation():
    counts = sp.csr_matrix(np.array([[0, 4, 2], [3, 0, 1]], dtype=np.float32))
    stored = counts.copy()
    before = namespace["sparse_fingerprint"](stored)
    transformed = counts.copy()
    transformed.data = np.log1p(transformed.data)
    assert namespace["sparse_fingerprint"](stored) == before
    assert namespace["sparse_fingerprint"](transformed) != before
    assert np.allclose(stored.data, np.rint(stored.data))

def test_partition_permutation():
    low = pd.Series(["a", "a", "b", "b", "c", "c"])
    high = pd.Series(["7", "7", "2", "2", "0", "0"])
    counts, fractions, stability, ari = namespace["partition_comparison"](low, high, cfg)
    assert ari == 1.0
    assert (stability.change == "stable one-to-one").all()
    assert counts.to_numpy().sum() == len(low)
    assert np.allclose(fractions.sum(axis=1), 1)

def test_partition_split_merge_and_alignment():
    low = pd.Series(["a"]*4 + ["b"]*4)
    split = pd.Series(["x"]*2 + ["y"]*2 + ["z"]*4)
    _, _, stability, ari = namespace["partition_comparison"](low, split, cfg)
    assert 0 < ari < 1
    assert stability.set_index("low_cluster").loc["a", "change"].startswith("split")
    assert stability.set_index("low_cluster").loc["b", "change"] == "stable one-to-one"
    merged = pd.Series(["x"]*8)
    _, _, stability, ari = namespace["partition_comparison"](low, merged, cfg)
    assert ari == 0
    assert (stability.change == "boundary shift / merge").all()
    assert (stability.dominant_fraction == 1).all()
    assert (stability.dominant_child_purity == .5).all()
    try:
        namespace["partition_comparison"](low, split.set_axis(range(8, 16)), cfg)
    except AssertionError:
        pass
    else:
        raise AssertionError("Mismatched cell indices were not rejected")

def evidence(active):
    genes = ["CD3D", "CD3E", "IL7R", "CD8A", "MS4A1", "CD79A", "LYZ", "CD14", "FCGR3A", "MS4A7", "NKG7", "GNLY", "FCER1A", "CST3", "PPBP", "PF4"]
    return pd.DataFrame({"gene": genes, "fraction_detected": [.8 if g in active else .02 for g in genes],
                         "mean_difference": [1.0 if g in active else -.1 for g in genes]})

def test_annotation_conservative():
    propose = namespace["propose_identity"]
    label, status, _ = propose(evidence(["MS4A1", "CD79A"]))
    assert label == "B-like" and "provisional" in status
    assert propose(evidence(["MS4A1"]))[1] == "uncertain"
    assert propose(evidence(["MS4A1", "CD79A", "CD3D", "CD3E"]))[1] == "mixed/uncertain"
    assert propose(evidence(["CD3D", "CD3E", "NKG7", "GNLY"]))[0] == "T-like (cytotoxic program)"
    assert propose(evidence(["NKG7", "GNLY"]))[0] == "NK-like"
    assert propose(evidence(["LYZ", "CD14", "FCER1A", "CST3"]))[1] == "mixed/uncertain"
    assert propose(evidence(["PPBP", "PF4"]))[0] == "Platelet-like"
    missing_t = evidence(["NKG7", "GNLY"]).query("gene not in ['CD3D', 'CD3E']")
    assert "uncertain" in propose(missing_t)[0]  # missing counterevidence is not negative evidence

for name, action in [
    ("All Python and 11 notebook code cells compile; notebook structure, empty outputs, source equivalence, and pins agree", test_python_and_notebook),
    ("Synthetic QC exact boundaries and overlapping failure counts", test_qc_boundaries),
    ("Synthetic sparse count copy stays intact after expression transformation", test_fingerprint_preservation),
    ("Synthetic ARI is invariant to permuted cluster names", test_partition_permutation),
    ("Synthetic split/merge, purity, and cell-index checks", test_partition_split_merge_and_alignment),
    ("Synthetic multi-marker, mixed, cytotoxic, missing-marker, and insufficient-evidence cases", test_annotation_conservative),
]:
    check(name, action)

report = {"status": "passed", "scope": "Code integrity and synthetic helper smoke tests only",
    "python": platform.python_version(), "checks": results,
    "pbmc_analysis_executed": False, "scanpy_integration_tested": False,
    "blocked_reason": "Scanpy/AnnData/graph dependencies absent; pip cannot resolve pypi.org in this environment",
    "notebook_code_cells": 11,
    "remaining": "Run notebook cells 1-10 in fresh Colab Python 3.11/3.12; cell 11 downloads results",
    "synthetic_data_are_not_biological_results": True}
(ROOT / "validation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))
