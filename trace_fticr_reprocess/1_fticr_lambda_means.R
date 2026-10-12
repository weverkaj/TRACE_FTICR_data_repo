
library(tidyverse)

source("getLambda.R")

# FTICR-MS processed data ----


data <- read_csv("output/Trace_ICR_data_SN12_sigma8_final.csv")

mol <- read_csv("output/Trace_ICR_mol_SN12_sigma8_final.csv")


## derived dir

out_dir = file.path("output/derived_clean")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Remove sample 194 that did not calibrated well

## Formulas detected only in the poorly calibrated sample

bad_sample <- "Sihi_60398_194"

data_long <- data |>
  pivot_longer(-Sample_ID, names_to = "Molecular_Formula", values_to = "intensity") |>
  filter(intensity > 0) # zeros are absent formulas

formulas_only_bad <- data_long |>
  filter(n() == 1, Sample_ID == bad_sample, .by = Molecular_Formula) |>
  pull(Molecular_Formula)

mol_ft <- mol |>
  filter(!Molecular_Formula %in% formulas_only_bad)

data_long <- data_long |>
  filter(!Molecular_Formula %in% formulas_only_bad) %>% 
  filter(Sample_ID != bad_sample)


## Lambda

# get_compositions() expects a `MolForm` column
lambda_comp <- mol_ft |>
  rename(MolForm = Molecular_Formula) |>
  get_compositions()

# positions follow the first 12 values returned by getThermoStoich()
thermo_cols <- c(delGcox0 = 1, delGd0 = 2, lambda0 = 6, delGcox = 7, delGd = 8, lambda = 12)

index_cols <- c("OC", "HC", "NOSC", "GFE", "DBE", "AI_mod")


lambda_mol <- get_lambda(lambda_comp$chemical_compositions)[, thermo_cols] |>
  as_tibble(.name_repair = ~ names(thermo_cols)) |>
  mutate(Molecular_Formula = lambda_comp$formulas, .before = 1) |>
  drop_na(lambda)


## Joining to mol

mol_ready <- mol %>% 
  left_join(lambda_mol)

write_csv(
  mol_ready,
  file.path(out_dir, "trace_mol_ready.csv")
)

## Mean, medians, quartiles

data_long <- data_long |>
  inner_join(mol_ready, by = "Molecular_Formula")

ready_summary <- data_long |>
  summarise(
    n_formulas = n(),
    across(
      c(names(thermo_cols), all_of(index_cols)),
      list(
        mean = mean,
        wmean = \(x) weighted.mean(x, intensity),
        median = median,
        q25 = \(x) quantile(x, 0.25),
        q75 = \(x) quantile(x, 0.75)
      ),
      .names = "{.col}_{.fn}"
    ),
    .by = Sample_ID
  )


## Save summary


write_csv(
  ready_summary,
  file.path(out_dir, "trace_mol_summary_sample.csv")
)






