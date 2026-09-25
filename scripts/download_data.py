#!/usr/bin/env python3
"""Download real external datasets for mega27-14 (digital embryo).

Provenance is recorded in data/DATA_MANIFEST.md. Every accession-level
record (one embryo gradient, one BioModels model, ...) counts as one
dataset per the program's uniform counting rule.
"""
import hashlib
import json
import pathlib
import urllib.request

ZENODO_LIU2013 = "https://zenodo.org/api/records/4942019/files/{key}/content"
FILES = {
    "liu2013/LiveImaging.mat": (
        ZENODO_LIU2013.format(key="Liu%20et%20al.%20-%20LiveImaging_data_structure.mat"),
        None,  # md5 filled after first download, then enforced
    ),
    "liu2013/Immunofluorescence.mat": (
        ZENODO_LIU2013.format(key="Liu%20et%20al.%20-%20Immunofluorescence_data_structure.mat"),
        None,
    ),
    "liu2013/README_IF.txt": (
        ZENODO_LIU2013.format(key="README_for_Liu%20et%20al.%20-%20Immunofluorescence_data_structure.txt"),
        None,
    ),
    "liu2013/README_Live.txt": (
        ZENODO_LIU2013.format(key="README_for_Liu%20et%20al.%20-%20LiveImaging_data_structure.txt"),
        None,
    ),
}

def fetch(url: str, dest: pathlib.Path) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=120) as r, open(dest, "wb") as f:
        f.write(r.read())
    return hashlib.md5(dest.read_bytes()).hexdigest()

def main() -> None:
    root = pathlib.Path(__file__).resolve().parent.parent / "data" / "raw"
    report = {}
    for rel, (url, md5) in FILES.items():
        dest = root / rel
        if dest.exists():
            digest = hashlib.md5(dest.read_bytes()).hexdigest()
            status = "cached"
        else:
            digest = fetch(url, dest)
            status = "downloaded"
        if md5 and digest != md5:
            raise SystemExit(f"md5 mismatch for {rel}: {digest} != {md5}")
        report[rel] = {"url": url, "md5": digest, "status": status}
        print(f"{status:>10} {rel} md5={digest}")
    (root / "download_report.json").write_text(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
