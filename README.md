# TRACE FTICR Data

Curated FTICR data package for sharing TRACE porewater FT-ICR molecular composition data and derived temporal-analysis outputs with collaborators.


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
    fticr_integration_output/
      *.csv
scripts/
  fticr_integration/
    *.R
    *.Rmd
```

## Data Contents

`data/raw/merged_output/Test_Processed-Unprocessed_Mol.csv` contains molecular formula annotations and calculated chemical properties.

`data/raw/merged_output/Test_Processed-Unprocssed_Data.csv.gz` contains the large FTICR formula-by-sample data matrix as a gzip-compressed CSV. The filename preserves the upstream spelling.

`data/raw/sample_key/60398_Sihi_Porewater_July82024_DS.csv` contains the sample key used to map FTICR sample IDs to TRACE plot, treatment, depth, and date metadata.

`data/cleaned/fticr_all.csv` and `data/cleaned/trace_fticr_all.csv.gz` are cleaned TRACE FTICR exports used by downstream analyses.

`data/derived/fticr_chemodiversity_metrics.csv` contains sample-level chemodiversity metrics.

`data/derived/fticr_integration_output/` contains derived FTICR integration outputs, including interpolation tables, formula-count temporal metrics, pairwise temporal distance data, semivariance slope tests, CAR(1) autocorrelation summaries, model coefficients, driver-analysis outputs, and plot-ready CSVs.

`scripts/fticr_integration/` contains the R scripts and interactive R Markdown workbook used to generate the current temporal variability and interpolation analyses.

## Compressed CSVs

Two CSVs are stored as `.csv.gz` files so the repository can be pushed to GitHub without Git LFS:

- `data/raw/merged_output/Test_Processed-Unprocssed_Data.csv.gz`
- `data/cleaned/trace_fticr_all.csv.gz`

They can be decompressed with:

```bash
gunzip -k data/raw/merged_output/Test_Processed-Unprocssed_Data.csv.gz
gunzip -k data/cleaned/trace_fticr_all.csv.gz
```
