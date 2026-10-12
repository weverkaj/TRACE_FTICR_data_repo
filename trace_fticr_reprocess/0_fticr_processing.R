library(furrr)
library(readxl)
library(tidyverse)


# FTICR-MS data ----

fticr_files_all <- list.files('Sihi_TRACE_output/', pattern = '.csv',
                              full.names = TRUE) 

blank_files <- fticr_files_all[str_detect(fticr_files_all, 'QC_Blank')]

fticr_files_samples <- fticr_files_all[!str_detect(fticr_files_all, 'QC|SRFA|Blank')]

plan(multisession, workers = 12)

blank_dataframes <- future_map(blank_files, function(path){
  
  print(path)
  if(str_detect(path, '.csv')){
    data <- read_csv(path) 
  } else if(str_detect(path, '.xlsx')){
    data <- read_excel(path, sheet = 'Sheet1') 
  }
  
  data_filt <- data %>% 
    filter(`Ion Type` == 'de-protonated') %>% 
    drop_na(`Molecular Formula`) %>% 
    rename(Molecular_Formula = `Molecular Formula`) %>% 
    mutate(filename = tools::file_path_sans_ext(basename(path)),
           qc_date = str_remove(filename, '.*_Fir_'),
           qc_date = str_remove(qc_date, '_300SA.*')) %>% 
    select(Molecular_Formula, qc_date, filename) %>% 
    distinct()
  
  return(data_filt)
  
}) %>% reduce(rbind)

blank_formulas <- blank_dataframes %>%
  group_by(Molecular_Formula) %>% 
  count() %>% 
  filter(n >= 4) %>% 
  select(Molecular_Formula) %>% 
  distinct()

data_frames <- future_map(fticr_files_samples, function(path){
  
  # Rounding data to be able to merge peaks that are highly similar to one another,
  # For masses with molecular formula, keeping only those with the best score
  
  
  data <- read_csv(path,
                   col_types = paste0(c(rep('d', 17), 'cccddcdddddddd'),
                                      collapse = ''))
  
  if('P' %notin% colnames(data)) data <- mutate(data, P = NA)
  if('S' %notin% colnames(data)) data <- mutate(data, S = NA)
  
  data_filt_blanks <- data %>% 
    drop_na(`Molecular Formula`) %>% 
    rename(Molecular_Formula = `Molecular Formula`, 
           mz = `Calibrated m/z`) %>% 
    mutate(filename = tools::file_path_sans_ext(basename(path))) %>% 
    filter(!Molecular_Formula %in% blank_formulas$Molecular_Formula,
           `Is Isotopologue` == 0,
           `Ion Type` == 'de-protonated') %>% 
    select(Sample_ID = filename, mz, `Peak Height`, 
           `Confidence Score`, `m/z Error (ppm)`, Molecular_Formula, 
           C, H, N, O, P, S) %>% 
    group_by(mz) %>% 
    slice_max(n = 1, order_by = `Confidence Score`) %>% 
    ungroup()
  
  data_ready <- data_filt_blanks %>% 
    ## To avoid floting point issues with R
    mutate(
      mz = mz * 1000000,
      mz = ceiling(mz),
      mz = mz / 1000000,
      mz = format(mz, 6))
  
  peak_abundances <- data_ready %>% 
    select(Molecular_Formula, Sample_ID, Intensity = `Peak Height`) %>% 
    group_by(Molecular_Formula, Sample_ID) %>% 
    summarise(Intensity = sum(Intensity)) %>%
    ungroup()
  
  peak_annot <- data_ready %>% 
    select(-c(`Peak Height`, `Confidence Score`, `m/z Error (ppm)`, -mz)) %>% 
    distinct() %>% 
    select(Molecular_Formula, everything()) %>% 
    mutate(across(!Molecular_Formula, ~ifelse(is.na(.x), 0, .x)))
  
  res <- list(abundances = peak_abundances,
              annotation = peak_annot)
})

all_peak_abundances <- map(data_frames, ~.x$abundances) %>% reduce(rbind) 

peak_abundances_derep <- all_peak_abundances %>% 
  mutate(Sample_ID = str_remove(Sample_ID, '_r[1-4].*')) %>% 
  group_by(Sample_ID, Molecular_Formula) %>% 
  mutate(rep_count = n()) %>% 
  group_by(Sample_ID) %>% 
  filter(rep_count == max(rep_count)) %>%
  group_by(Sample_ID, Molecular_Formula) %>%
  summarise(Intensity = ceiling(mean(Intensity))) %>% 
  pivot_wider(names_from = Molecular_Formula, values_from = Intensity, values_fill = 0)

write_csv(peak_abundances_derep, 'output/Trace_ICR_data_SN12_sigma8_final.csv')

all_annotations <- map(data_frames, ~.x$annotation) %>% reduce(plyr::rbind.fill) %>% 
  filter(Molecular_Formula %in% colnames(peak_abundances_derep)) %>% 
  select(-Sample_ID, -mz) %>% 
  distinct()

all_annotations_ready <- all_annotations %>%  
  mutate(OC = O / C,
         HC = H / C,
         NOSC = 4 - ((4*C + H - 3*N - 2*O + 5*P - 2*S)/
                       C),
         GFE = 60.3 - 28.5*NOSC,
         DBE = 1 + 0.5*(2*C - H + N + P),
         AI = (1 + C - O - S - 0.5*(H + N + P))/
           (C - O - S - N - P),
         AI_mod = (1 + C - 0.5*O - S -0.5*(N + P + H))/
           (C - 0.5*O - N - S - P),
         AI_mod = ifelse(AI_mod < 0, 0, AI_mod),
         AI_mod = ifelse(is.infinite(AI_mod), 0, AI_mod),
         AI = ifelse(is.infinite(AI), 0, AI),
         AI = ifelse(AI < 0, 0, AI),
         heteroatom = str_remove(Molecular_Formula, 'C\\d+ H\\d+ ')) %>% 
  mutate(Class = case_when(between(OC, 0, 0.3) & between(HC, 1.5, 2.5) ~ 'Lipid',
                           between(OC, 0, 0.125) & between(HC, 0.8, 1.5) ~ 'Unsat. HC',
                           between(OC, 0, 0.95) & between(HC, 0.2, 0.8) ~ 'Cond. HC',
                           between(OC, 0.3, 0.55) & between(HC, 1.5, 2.3) ~ 'Protein',
                           between(OC, 0.55, 0.7) & between(HC, 1.5, 2.2) ~ 'Amino sugar',
                           between(OC, 0.7, 1.5) & between(HC, 1.5, 2.5) ~ 'Carbohydrate',
                           between(OC, 0.125, 0.65) & between(HC, 0.8, 1.5) ~ 'Lignin',
                           between(OC, 0.65, 1.1) & between(HC, 0.8, 1.5) ~ 'Tannin'))

write_csv(all_annotations_ready, 'output/Trace_ICR_mol_SN12_sigma8_final.csv')
