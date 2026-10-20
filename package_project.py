"""Package only intentional deliverables; omit environments, caches, and unrelated projects."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
import json

ROOT = Path(__file__).resolve().parent
names = ["Single_Cell_Immune_Profiling.ipynb", "README.md", "requirements.txt", "analysis.py",
         "INTERVIEW_GUIDE.md", "RESUME_BULLETS.md", "REFERENCES.md", "build_notebook.py",
         "validate_project.py", "validation_report.json", "package_project.py"]
manifest = {name: {"bytes": (ROOT/name).stat().st_size,
    "sha256": hashlib.sha256((ROOT/name).read_bytes()).hexdigest()} for name in names}
(ROOT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
archive = ROOT.parent / "outputs" / "Single_Cell_Immune_Profiling.zip"
archive.parent.mkdir(parents=True, exist_ok=True)
with ZipFile(archive, "w", ZIP_DEFLATED) as out:
    for name in [*names, "MANIFEST.json"]:
        out.write(ROOT / name, f"Single_Cell_Immune_Profiling/{name}")
with ZipFile(archive) as check:
    assert check.testzip() is None
    assert len(check.namelist()) == len(names) + 1
    for name, record in manifest.items():
        assert hashlib.sha256(check.read(f"Single_Cell_Immune_Profiling/{name}")).hexdigest() == record["sha256"]
print(f"Verified archive: {archive} | {archive.stat().st_size:,} bytes | {len(names)+1} files")
