#!/usr/bin/env bash
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

run_stage() {
  local script="$1"
  echo "===== ${script} ====="
  Rscript "${HERE}/${script}"
}

run_stage 01_compute_and_queue.R || exit $?
run_stage 08_export_plot_data.R
run_stage 02_render_core_gsea_plots.R
run_stage 03_render_functional_distance_plots.R
run_stage 04_render_expanded_evidence_plots.R
run_stage 05_render_summary_lfcpv_plots.R
run_stage 07_finalize_plot_status.R
