#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
MPLBACKEND=Agg python main.py
python -m pytest -q
