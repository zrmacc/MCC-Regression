# MCC regression replication

Zachary R. McCaw <br>
Updated: 2026-11-07

Code and display artifacts for *Pseudo-value-Based Mean Cumulative Count Regression*. Figure and table numbers match the current manuscript and supplement. Figure 5 (ORATORIO) is supplied by the clinical collaborators; participant data and clinical analysis code are not distributed here.

All reported regression inference uses HC3-adjusted standard errors with the second-order pseudo-value correction and normal Wald intervals and tests. The corrected subject contribution is `h11_hc3 + h12`; its uncentered outer product uses denominator n, with no additional residual-degrees-of-freedom multiplier.

## Software and validation

Use an MCC source build containing the tested HC3-plus-second-order implementation, commit `5dbe0a02f5d5a520497b4501439f58fed61c22c1`. A version number alone does not establish this capability. The simulation drivers run a small fit to verify `PseudoReg()`, its method string, and the manual corrected covariance before starting.

For example, install the pinned source in a separate library and run the capability tests:

```sh
mkdir -p .Rlib
export R_LIBS_USER="$PWD/.Rlib"
R CMD INSTALL --library="$R_LIBS_USER" /path/to/pinned/MCC
Rscript scripts/validation/validate_corrected_inference.R
```

Additional R packages are arrow, dplyr, ggplot2, glue, optparse, tidyr, gridExtra, and jsonlite. The command above assumes that the local MCC source is at the pinned commit. The package also reports its own dependencies during installation. `type="AUC"` denotes AUMCF; `type="MCF"` denotes MCF.

Scripts locate this checkout from their own path and accept absolute input/output paths. New raw results are written under the Git-ignored `simulations/` directory. Sample sizes are per arm. The retained historical pseudo-value concordance study is single-arm and is not included in the current supplement.

## Display map

| Display | Simulation scripts in `scripts/sims/` | Display script | Supplied output |
|---|---|---|---|
| Figure 1: bias | `sim_estimation.R`, `sim_estimation_truth.R` | `plot_figure_1.R` | `fig1a.pdf`, `fig1b.pdf` |
| Figure 2: coverage | same corrected run and truth as Figure 1 | `plot_figure_2.R` | `fig2a.pdf`, `fig2b.pdf` |
| Figure 3: null QQ | `sim_type1.R` | `plot_figure_3.R` | `fig3.pdf`, `fig3.png` |
| Figure 4: power | `sim_power.R` | `plot_figure_4.R` | `fig4.pdf` |
| Web Figure 1: frailty bias and coverage | `sim_frailty.R`, `sim_frailty_truth.R` | `plot_figure_S2.R` | `figS2a.pdf`, `figS2b.pdf` |
| Web Figure 2: terminal-event MCF | `sim_terminal_event.R`, `sim_terminal_event_truth.R` | `plot_figure_S3.R` | `figS3a.pdf`, `figS3b.pdf` |
| Web Figure 3: terminal-event AUMCF | same corrected run as Web Figure 2 | `plot_figure_S4.R` | `figS4a.pdf`, `figS4b.pdf` |
| Table 1: type I error | `sim_type1.R` | `make_table_1.R` | `tab1.tex` |
| Web Table 1 / Table S1: SE calibration | same corrected run and truth as Figures 1–2 | `make_table_S1.R` | `tabS1.tex` |

Web Table 1 reports the empirical Monte Carlo SD of the AUMCF treatment estimates, RMS corrected SE, and corrected 95% coverage in 27 cells with 2,000 replicates each. It uses the same estimates, SEs, and normal intervals as Figure 2. Empirical SD is a simulation benchmark. The historical concordance display (`plot_figure_S1.R`, `figS1a.png`, and `figS1b.png`) is retained but is not an active manuscript figure. Historical `S2`–`S4` script/output names remain unchanged; they generate the current Web Figures 1–3.

Figures 1 and 2 use 10 × 8-inch panels, included at full manuscript width on continued figure pages so their increased height does not reduce text size. Supplemental bias and coverage panels use the same 10 × 8-inch canvas and continued-page layout. Terminal-event figures export separate `a` (bias) and `b` (coverage) PDFs for the manuscript; combined `figS3.pdf`/`figS4.pdf` exports remain available at 10 × 16 inches. Figure 3 uses an 8 × 8-inch canvas, and Figure 4 uses a taller 9 × 7-inch canvas. `MainFigureBaseSize()` scales base type with source width, giving the same printed text size across the generated figures: 11.7 points at a 6.5-inch manuscript width.

All confidence-interval coverage axes display 80%–100%, with ticks every 5 percentage points and a dashed 95% reference. Coordinate zoom preserves every plotted estimate and Monte Carlo interval rather than dropping observations through scale limits. Bias limits remain symmetric around zero. Caption, subcaption, and alt-text spacing is single within the continued bias/coverage figures. Coefficients use β for both estimands. Figure numbers and source-to-manuscript filenames are defined by `DisplayAssetMap()`; validation groups continued panels before checking figure numbers, labels, and byte-identical asset copies. Figure 5 remains collaborator-supplied.

## Running the analyses

From this checkout, use one R process and one BLAS thread per independent job:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# A small capability/simulation test, written outside the production directories.
Rscript scripts/sims/sim_estimation.R --testing --n 50 --tau 2 --type AUC --beta_name beta_A --beta_value 0 --reps 2 --workers 1 --outdir simulations/pilot --batch pilot

# Full designs are long-running simulations.
Rscript scripts/sims/sim_estimation_truth.R --n_ground 10000 --reps 1000
Rscript scripts/sims/sim_estimation.R --reps 2000 --seed 2026090401 --workers 1 --batch replication --resume
Rscript scripts/sims/sim_type1.R --reps 100000 --workers 1 --batch replication --resume
Rscript scripts/sims/sim_power.R --reps 500 --seed 2026090404 --workers 1 --batch replication --resume
Rscript scripts/sims/sim_pseudovalue.R
Rscript scripts/sims/sim_frailty_truth.R --reps 100 --seed 20260904 --workers 1
Rscript scripts/sims/sim_frailty.R --reps 2000 --seed 2026090402 --workers 1 --batch replication --resume
Rscript scripts/sims/sim_terminal_event_truth.R --n_ground 10000 --reps 1000
Rscript scripts/sims/sim_terminal_event.R --n_blocks 1 --block_id 1 --reps 2000 --batch replication --resume
```

For independent baseline, frailty, and power jobs, filter the published design using CLI options and use separate output paths or batch names. Canonical configuration indices and per-replicate seeds remain fixed under filtering. Resume validates the method, cell metadata, and complete replicate IDs before skipping a saved cell. Type I error additionally supports `--rep_start` for independent replicate batches; terminal-event runs support cell blocks and saved shards. Keep new pilot outputs separate from production.

After the corrected simulations and their truth inputs are available, regenerate the displays:

```sh
Rscript scripts/reproduce_displays.R
```

The driver invokes all inferential figure/table scripts, regenerates the corresponding bias panels, and writes a generation manifest and command logs under `provenance/`. The historical concordance outputs are not regenerated or copied into the manuscript. Each inferential loader rejects missing or incorrect covariance metadata, duplicate keys, incomplete designs, and incomplete replicate counts. Where coefficient and SE fields exist, the loaders also verify the normal-Wald interval endpoints and p-values.

For existing saved production data, supply `--data_root`, `--event_truth`, `--terminal_truth`, and `--frailty_truth` explicitly. Use `--out_dir`, `--table_dir`, and `--manifest` to stage output before replacing supplied artifacts. An optional `--manuscript_dir` synchronizes the generated files into the existing manuscript asset map, including the caption-free longtable body of Table S1. Before generation or copying, the driver checks the active figure numbers, labels, and filenames. After synchronization, it also checks byte equality of the generated and manuscript assets. Commented-out figures are excluded. The generation manifest records the active map and hashes of the input data, output files, and plotting sources. The supplement source owns that table's caption, label, and surrounding environment. Individual scripts still accept `--help`; table scripts require `--overwrite` when replacing an existing table.

To check the supplied figures against an existing manuscript checkout without regenerating them:

```sh
Rscript scripts/validation/validate_display_assets.R --manuscript_dir /path/to/manuscript --figure_dir figures
```

The check includes all ten generated PDF assets and verifies the reference to the separately supplied clinical figure. The manuscript checkout must contain `main.tex`, `supplement.tex`, and `figures/`.

## Reproducibility

Raw production results and frozen truth are not bundled. The supplied displays were generated from validated saved results; their original generation record is `provenance/current_displays.json`. The readability refresh of Figures 1–4 is recorded in `provenance/figure_readability_20261106.json`; it reuses the same validated inputs and changes only plotting typography, dimensions, and manuscript layout. The matching supplemental refresh is recorded in `provenance/supplement_readability_20261107.json`; the saved input files and plotted numerical summaries were verified unchanged. The full display-driver audit, active figure map, source hashes, and reproduction checks are recorded in `provenance/display_alignment_20261107.json`. New Monte Carlo results need not be byte-identical, and PDF timestamps can differ even when the numerical summaries agree. Keep each truth-generation run separate and use an explicit truth file or directory; do not pool unrelated runs.

The historical concordance script retains a fixed grid traversal, including unused tau=1 and n=400 cells, to reproduce the plotted cells' random sequence at seed 3141. This display concerns pseudo-value agreement rather than regression inference.

The ORATORIO clinical collaborators maintain the application analysis. The replication scripts neither generate its forest plot nor provide a local clinical-data rerun.
