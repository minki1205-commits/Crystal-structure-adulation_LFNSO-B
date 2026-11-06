#!/usr/bin/env bash
set -euo pipefail
cd /workspace/Crystal-structure-adulation_LFNSO-B
export MPLCONFIGDIR=/workspace/lfnso-cache/matplotlib
mkdir -p "$MPLCONFIGDIR"
python3 -m venv /workspace/lfnso-venv
/workspace/lfnso-venv/bin/python -m pip install -r requirements.lock.txt
/workspace/lfnso-venv/bin/python -m pip check
/workspace/lfnso-venv/bin/python -m pytest -q
