# TRACE FTICR Data

Curated data package and reproducible analysis scripts for TRACE porewater FT-ICR molecular composition data.

## Directory Layout

```text
data/
  raw/
    merged_output/
      Test_Processed-Unprocessed_Mol.csv
      Test_Processed-Unprocssed_Data.csv.gz
    sample_key/
      60398_Sihi_Porewater_July82024_DS.csv
  cleaned/
    fticr_all.csv
    trace_fticr_all.csv.gz
  derived/
    fticr_chemodiversity_metrics.csv
    monthly_covariates/
      fticr_monthly_covariates.csv
      fticr_monthly_covariates_metadata.csv
scripts/
  fticr_integration/
    *.R
    *.Rmd
output/                         # generated locally; ignored by Git
  fticr_integration/
```

## Data Contents

`data/raw/merged_output/Test_Processed-Unprocessed_Mol.csv` contains molecular formula annotations and calculated chemical properties.

`data/raw/merged_output/Test_Processed-Unprocssed_Data.csv.gz` contains the large FTICR formula-by-sample data matrix as a gzip-compressed CSV. The filename preserves the upstream spelling.

`data/raw/sample_key/60398_Sihi_Porewater_July82024_DS.csv` contains the sample key used to map FTICR sample IDs to TRACE plot, treatment, depth, and date metadata.

`data/cleaned/fticr_all.csv` and `data/cleaned/trace_fticr_all.csv.gz` are cleaned TRACE FTICR exports used by downstream analyses.

`data/derived/fticr_chemodiversity_metrics.csv` contains sample-level chemodiversity metrics.

## Monthly TRACE Covariates

`data/derived/monthly_covariates/fticr_monthly_covariates.csv` is a complete plot-by-calendar-month grid for May 2017 through April 2024. It has 504 rows (six plots by 84 months) and 84 columns and can be joined to FT-ICR samples using `plot` and `year_month`. Existing FT-ICR-derived columns are excluded to prevent outcome leakage.

The environmental profile contains monthly mean temperature and VWC at 0-10, 20-30, and 40-50 cm from the central CS655 sensors, plus monthly mean plot-level soil CO2 flux. It also retains selected porewater chemistry, nutrient-core and root-soil-core measurements, nutrient-core microbial biomass, and root stock, growth, and mortality at 0-10, 10-20, and 20-30 cm.

Sparse biogeochemical measurements are not interpolated or imputed. Multiple observations within a plot-month are represented by their arithmetic mean; blank cells indicate that no source observation was available.

`data/derived/monthly_covariates/fticr_monthly_covariates_metadata.csv` defines every output column, including its source variable and file, units, depth, aggregation, data type, and coverage.

Rebuild both files with:

```bash
Rscript scripts/fticr_integration/15_build_monthly_covariates.R
```

The script reads the hourly and daily `joined_data_flux_final_2` derivatives from `/Users/jrweverk/Documents/TRACE_MCMC` by default. Set `TRACE_MCMC_DIR` to use a different checkout.

## Analysis Outputs

The scripts in `scripts/fticr_integration/` generate interpolation tables, temporal-distance and semivariance tests, continuous-time autocorrelation summaries, temporal-driver models, calendar-time trends, cross-plot convergence tests, and figures.

These regenerated products are written to or retained under `output/fticr_integration/`, which is intentionally ignored by Git. This keeps the shared repository focused on source data, final collaborator datasets, and the code needed to reproduce analyses.

Run the calendar-time trend and cross-plot convergence analysis with:

```bash
Rscript scripts/fticr_integration/14_calendar_time_trends_and_cross_plot_convergence.R
```

Its figures are written to `output/fticr_integration/figures/` without embedded titles or subtitles.

## Compressed CSVs

Two CSVs are stored as `.csv.gz` files so the repository can be pushed to GitHub without Git LFS:

- `data/raw/merged_output/Test_Processed-Unprocssed_Data.csv.gz`
- `data/cleaned/trace_fticr_all.csv.gz`

They can be decompressed with:

```bash
gunzip -k data/raw/merged_output/Test_Processed-Unprocssed_Data.csv.gz
gunzip -k data/cleaned/trace_fticr_all.csv.gz
```

The calendar-time analysis can read the compressed formula-by-sample matrix through the system `gzip` command when `R.utils` is unavailable. It prefers the decompressed CSV when both versions are present.
