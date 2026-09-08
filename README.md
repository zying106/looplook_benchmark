# Looplook manuscript analysis

This repository contains the analysis code used to evaluate Looplook target-gene assignment in LPS141 cells. It integrates processed chromatin-loop, ChIP-seq, RNA-seq, ATAC-seq, and histone-mark data, compares Looplook assignments with ChIPseeker, and generates the numerical results and figures reported in the manuscript.

## Contents

- `analysis/target_generation/`: generates Looplook target-gene assignments across four TSS windows.
- `analysis/looplook_benchmark_v12_plot_data_master.R`: runs the primary benchmark and downstream statistical analyses.
- `analysis/figure_scripts/`: recreates manuscript figures from benchmark outputs.
- `analysis/pipeline_v13/`: runs the benchmark and queued figure rendering.
- `package_snapshot/`: Looplook v0.99.19 source used by the workflow.
- `manifests/`: figure-to-script mapping and reported analysis parameters.

## Primary benchmark settings

The reported benchmark uses 16 hop0 assignment modes, 300 resampling iterations, and 200 genes per iteration. The modes combine four annotation pipelines with two target scopes (`all` and `promoter`) and two assignment policies (`strict` and `filled`).

Analyses were performed using Looplook, and the reproducible implementation is publicly available as version 0.99.19.

## Running the analysis

The workflow starts from processed loop, peak, expression, and chromatin-track files. Configure the dataset-specific filenames in the target-generation and benchmark scripts, then set the required directory variables:

```bash
export LOOPLOOK_TARGET_DATA_DIR=/path/to/processed_inputs
export LOOPLOOK_TARGET_OUT_DIR=/path/to/target_assignments
Rscript analysis/target_generation/generate_multi_tss_data_NO_bedpe.R

export LOOPLOOK_DATA_BASE=/path/to/benchmark_inputs
export LOOPLOOK_RDATA_BASE=/path/to/target_assignments
export LOOPLOOK_OUT_BASE=/path/to/benchmark_results
export LOOPLOOK_CACHE_POLICY=resume
export LOOPLOOK_PLOT_POLICY=queue
Rscript analysis/looplook_benchmark_v12_plot_data_master.R
bash analysis/pipeline_v13/run_pipeline.sh --render-only
```

