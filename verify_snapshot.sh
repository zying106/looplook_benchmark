#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Checking R syntax"
while IFS= read -r script; do
  Rscript -e 'parse(file=commandArgs(TRUE)[1])' "$script" >/dev/null
  echo "OK  ${script#"${root_dir}/"}"
done < <(find "${root_dir}/analysis" "${root_dir}/package_snapshot/R" -type f -name '*.R' | sort)

echo "Checking shell syntax"
while IFS= read -r script; do
  bash -n "$script"
  echo "OK  ${script#"${root_dir}/"}"
done < <(find "${root_dir}" -type f -name '*.sh' | sort)

echo "Checking for bundled caches and generated analysis files"
if find "${root_dir}" -type f \( -name '*.RData' -o -name '*.rds' -o -name '*.pdf' -o -name '*.xlsx' \) | grep -q .; then
  echo "Unexpected generated files found" >&2
  exit 1
fi

echo "Checking for local metadata and personal analysis paths"
if find "${root_dir}" -type f -name '.DS_Store' | grep -q .; then
  echo "Unexpected .DS_Store files found" >&2
  exit 1
fi
if rg -n '/home/rstudio|~/workspace' \
  "${root_dir}/analysis/figure_scripts" \
  "${root_dir}/analysis/pipeline_v13" \
  "${root_dir}/analysis/looplook_benchmark_v12_plot_data_master.R"; then
  echo "Personal analysis paths found" >&2
  exit 1
fi

echo "Checking paper-default benchmark settings"
rg -q '^n_iterations <- 300L$' "${root_dir}/analysis/looplook_benchmark_v12_plot_data_master.R"
rg -q '^sample_sizes <- 200L$' "${root_dir}/analysis/looplook_benchmark_v12_plot_data_master.R"
rg -q '^primary_sample_size <- 200L$' "${root_dir}/analysis/looplook_benchmark_v12_plot_data_master.R"
rg -q '^primary_hop <- 0L$' "${root_dir}/analysis/looplook_benchmark_v12_plot_data_master.R"

echo "Snapshot structural checks passed"
