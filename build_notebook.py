"""Build a self-contained nbformat 4.5 notebook from the readable percent script."""
from pathlib import Path
import json
import re
import hashlib

ROOT = Path(__file__).resolve().parent
source = (ROOT / "analysis.py").read_text(encoding="utf-8")
parts = re.split(r"^# %%([^\n]*)\n", source, flags=re.MULTILINE)
cells = []
for i in range(1, len(parts), 2):
    kind, text = parts[i].strip(), parts[i+1].rstrip() + "\n"
    if kind == "[markdown]":
        text = "\n".join(line[2:] if line.startswith("# ") else line[1:] if line.startswith("#") else line
                         for line in text.splitlines()) + "\n"
        cell = {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(keepends=True)}
    else:
        compile(text, f"notebook-cell-{len(cells)+1}", "exec")
        cell = {"cell_type": "code", "execution_count": None, "metadata": {},
                "outputs": [], "source": text.splitlines(keepends=True)}
    cell["id"] = hashlib.sha256(f"{len(cells)}{text}".encode()).hexdigest()[:12]
    cells.append(cell)
notebook = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.12"},
    "colab": {"name": "Single_Cell_Immune_Profiling.ipynb", "provenance": []}}, "nbformat": 4, "nbformat_minor": 5}
destination = ROOT / "Single_Cell_Immune_Profiling.ipynb"
destination.write_text(json.dumps(notebook, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"Built {destination.name}: {sum(c['cell_type']=='code' for c in cells)} code cells; all compile.")
