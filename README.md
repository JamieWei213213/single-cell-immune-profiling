# Single-Cell Immune Profiling with Scanpy

A reproducible, explainable analysis of public PBMC3k sequencing counts, with a fixed-preprocessing Leiden-resolution experiment.

## Current status — read first

**Implemented, but the real PBMC analysis has not been executed in this delivery environment.** Python 3.12.10 is available; Scanpy, AnnData, igraph, and leidenalg were not installed. A project-local dependency installation failed because the execution environment could not resolve `pypi.org`. No quality-control numbers, clusters, biological findings, ARI score, marker tables, or PBMC plots have been fabricated. The source dimensions below are published metadata, not measurements from our run.

The notebook contains all code and explanations, with no external project-file dependency. Local code checks and synthetic helper smoke tests are recorded in `validation_report.json`; these do **not** validate a full Scanpy run. Dependency integration, data download, graph construction, figure generation, and `.h5ad` exports still require the Colab run below.

## Question and scope

Which broad immune-cell expression patterns are present, and how much do cluster memberships change between Leiden resolutions 0.4 and 0.8 when preprocessing and the graph are held fixed? This demonstrates computational bioinformatics, exploratory analysis, and sensitivity assessment. It is not a trained clinical predictor or wet-lab project.

## Run in Google Colab

1. Open [Google Colab](https://colab.research.google.com/), choose **File → Upload notebook**, and upload `Single_Cell_Immune_Profiling.ipynb`.
2. Use a fresh Python runtime. **CPU is sufficient**; a selected GPU is not used by these Scanpy calls. The pinned stack targets Python 3.11/3.12; cell 1 explicitly checks this. If Colab changes its default Python version, select a compatible runtime version where available or use a local Python 3.12 environment.
3. Run **cell 1**, the package installation cell. If the runtime reports a binary/ABI conflict after replacing preinstalled packages, restart the session and rerun from cell 1. Do not continue from a half-imported environment.
4. Run **cells 2–10 in order**, or use **Runtime → Run all**. Cell 4 downloads raw counts automatically. The first Numba compilation may take several minutes. CPU runtime and memory use are not benchmarked here; the dense PCA workspace is limited to the variable genes.
5. Cell 10 must report **COMPLETE**, and `run_status.json` must contain `"pbmc_analysis_complete": true`. Read generated `OBSERVATIONS.md`, the annotated stability table, both marker dot plots, and the uncertain labels.
6. **Cell 11** downloads the generated results ZIP in Colab. Also use **File → Download → .ipynb** to retain the executed notebook with its visible outputs. A Colab runtime is temporary; download results before disconnecting.

All biological analysis cells still need to run. There is no API key, private data upload, Drive mount, paid service, or manual cluster-number mapping requirement.

## Run locally

From this project directory in Windows PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe analysis.py
```

On Linux/macOS, use `python3.12 -m venv .venv`, `.venv/bin/python -m pip install -r requirements.txt`, then `.venv/bin/python analysis.py`. The script and notebook contain the same analysis. `python build_notebook.py` regenerates the notebook from the script. `python validate_project.py` performs local integrity and synthetic helper checks using NumPy, pandas, SciPy, and scikit-learn. Notebook structure is checked directly; this is not a full notebook-kernel or official nbformat-schema execution test.

## Dataset, provenance, and terms

[10x Genomics, 3k PBMCs from a Healthy Donor](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0), published May 26, 2016, contains **one sample from one healthy human donor** and **2,700 detected cells**. Cell Ranger 1.1.0 produced the counts. Scanpy's [official raw loader](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html) supplies 2,700 cells × 32,738 genes via `sc.datasets.pbmc3k()`; the compressed download is about 5.9 MB. These are raw molecule counts for upstream-called cells, not FASTQ reads or all empty droplets. We do not use `pbmc3k_processed()`.

The dataset page explicitly states **CC BY 4.0**. Credit 10x Genomics, link the [source](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0) and [license](https://creativecommons.org/licenses/by/4.0/), and identify derived changes when sharing. The pipeline writes this attribution into each results folder. Scanpy software has a separate BSD-3-Clause license; third-party libraries retain their own licenses. See `REFERENCES.md` for official evidence, workflow credit, and package release links checked October 20, 2026.

## Methods and decisions

| Stage | Implemented behavior | Reason / caveat |
|---|---|---|
| Input | Download official raw counts; verify shape, finite nonnegative integer values; save input `.h5ad` | Avoid accidental use of processed expression |
| Cell QC | Keep ≥200 and <2,500 detected genes; mitochondrial fraction <5% | Dataset-specific reference-style heuristics; not universal clinical cutoffs |
| Count diagnostics | Report total counts, detected genes, mitochondrial percentage; histogram and scatter before/after | No separate total-count threshold; failure categories overlap |
| Gene QC | Keep genes detected in ≥3 retained cells | Reduce extremely sparse features; may discard rare biology |
| Count preservation | Copy filtered counts to `layers['counts']`; retain whole input separately | Hash checks detect mutation; `.raw` is deliberately unused |
| Expression | Normalize each cell to 10,000 counts across retained genes; natural `log1p` | Reduce depth differences and dynamic range; does not solve all bias |
| Features | 2,000 Seurat-flavor dispersion HVGs on log data | Limit noise/computation; all genes retained for marker ranking |
| PCA | Scale a separate HVG copy, clip at 10; 40 PCs | Balance gene contributions; diagnostic variance plot; not optimized after viewing results |
| Graph | 10 nearest neighbors in PCA space | Define local similarity before clustering; graph computed once |
| Leiden | Explicit `leidenalg` backend, undirected weighted graph, seed 42, convergence iterations | Compare resolutions 0.4 and 0.8 with all else fixed |
| UMAP | One shared seeded embedding for both clusterings | Visualization only; does not establish identity or lineage |
| Markers | Wilcoxon with tie correction and BH adjustment, cluster versus rest, all retained log-expression genes | Exploratory rankings; same data used for clustering and testing |
| Annotation | Multi-marker rules; export supporting expression/fractions and competing lineages | Provisional proposals; uncertain and mixed labels preserved |
| Sensitivity | ARI, count cross-tab, row fractions, retention plus child purity | Do not compare numeric labels directly; groups need not be nested |

Baseline 0.4 is a declared broad-summary view, **not a claim of optimal resolution**. The experiment asks what changes at 0.8. Stable is defined as ≥90% parent retention and ≥90% child purity. A substantial split has at least two children each receiving ≥10% of the parent's cells. These are operational reporting rules, not statistical confidence thresholds.

## Reproduced workflow versus additions

Credit the official [Scanpy PBMC3k reference workflow](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html), which current documentation calls a legacy workflow, and the [Seurat PBMC3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) for marker interpretation. This project follows the broad reference analysis sequence. It does **not** claim exact numerical replication: gene-filter ordering, fixed HVG count, explicit Leiden backend/settings, no regression of counts/mitochondrial fraction, and modern package versions can alter results.

Project additions are the count-integrity checks, auditable QC waterfall, conservative reproducible marker proposals, explicit uncertainty, same-graph two-resolution comparison, group-level stability descriptions, one-cell trace, run manifest, and export round-trip checks.

## Actual observations

**Pending execution.** It is not yet known from this project's run how many cells pass QC, how many clusters appear, what markers rank highest, or which populations split. Do not substitute numbers from Scanpy's tutorial. The notebook generates `OBSERVATIONS.md` directly from actual runtime tables, including proposed identities linked to stable/split groups. The output documents measured membership changes; those changes alone do not verify biological populations.

## Files and outputs

| Delivered file | Purpose |
|---|---|
| `Single_Cell_Immune_Profiling.ipynb` | Standalone Colab notebook, 11 code cells |
| `analysis.py` | Matching readable executable pipeline |
| `requirements.txt` | Computational package pins |
| `build_notebook.py` | Regenerate notebook from the percent script |
| `validate_project.py`, `validation_report.json` | Integrity and clearly synthetic smoke checks |
| `INTERVIEW_GUIDE.md` | Architecture, data flow, tradeoffs, explanations, 8 Q&As |
| `RESUME_BULLETS.md` | Two truthful current-status bullets |
| `REFERENCES.md` | Source links, terms, and attribution |

After execution, `pbmc3k_runs/run_<timestamp>/` contains:

- `pbmc3k_input_counts.h5ad`, `pbmc3k_processed.h5ad`.
- `tables/cell_annotations.csv`, QC summaries and waterfall, gene metadata, PCA variance.
- Full compressed marker tables and top-20 tables at both resolutions, marker-panel evidence, cluster annotation proposals.
- Membership count/fraction cross-tabs and annotated cluster-stability results.
- PNG and PDF QC, PCA, UMAP, dot-plot, and overlap figures.
- `OBSERVATIONS.md`, `one_cell_trace.json`, `sensitivity.json`, configuration, environment, a resolved package freeze, attribution, and completion status.

## Limitations and next steps

One donor is one biological sample; thousands of cells are not thousands of independent donors. There are no cancer patients, treatment arms, or outcomes here. This dataset cannot support disease effects, treatment response, clinical accuracy, or population-level generalization. Clusters and UMAP positions do not independently establish identities or relationships. Marker rankings are circular exploratory summaries after clustering, not confirmatory differential expression. There is no protein assay, flow-cytometry label, pathology annotation, or wet-lab validation.

Thresholds can remove legitimate cells. A high gene count can suggest but does not prove a doublet. Formal doublet/ambient-RNA correction is not implemented. RNA dropout and technical variation can produce missing or misleading markers. The heuristic marker rules are intentionally conservative and may leave valid groups unassigned or label overlapping myeloid programs as mixed. Seed, QC, HVG, PCA, and neighbor-count sensitivity are not tested by this two-resolution experiment. Cross-platform numerical differences can persist despite fixed seeds and pinned main packages. A real replicated study would need additional donors, independent annotation evidence, appropriate sample-level statistics, and broader sensitivity checks.
