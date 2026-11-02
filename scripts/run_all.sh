#!/usr/bin/env bash
# Full reproduction: acquisition -> manifest -> models -> robustness -> single cell -> figures.
# Requires: Python 3.12/3.13, fonts-urw-base35 (Nimbus Sans), ~1 GB free disk, NCBI HTTPS access.
set -euo pipefail
cd "$(dirname "$0")/.."
[ -x .venv/bin/python ] || { python3 -m venv .venv && .venv/bin/pip install -q -r requirements.txt fonttools==4.60.1; }
PY=.venv/bin/python
$PY scripts/acquire_geo.py
$PY scripts/parse_soft.py GSE225224 GSE198895 GSE281091 GSE320302 GSE233163 GSE309905 > /dev/null
$PY scripts/download_supplementary.py
for a in GSE225224 GSE198895 GSE309905; do mkdir -p data/raw/$a/x && tar -xf data/raw/$a/${a}_RAW.tar -C data/raw/$a/x; done
$PY scripts/build_manifest.py
$PY scripts/run_bulk_de.py
$PY scripts/sensitivity.py
$PY scripts/sc_pseudobulk.py
$PY scripts/sc_descriptive.py
$PY scripts/context_gse309905.py
$PY scripts/make_figures.py
