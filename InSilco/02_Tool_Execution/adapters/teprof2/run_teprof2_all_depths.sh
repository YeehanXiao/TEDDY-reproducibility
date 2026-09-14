#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${1:-$PWD}"
cd "${REPO_ROOT}"
exec bash InSilco/02_Tool_Execution/06_run_teprof2.sh
