#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
MPLBACKEND=Agg python main.py
MPLBACKEND=Agg python reproduce.py
python -m pytest -q
