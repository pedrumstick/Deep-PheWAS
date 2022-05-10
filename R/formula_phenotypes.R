#! /usr/bin/env Rscript
'Creates formula phenotypes using data from previously generated phenotype scripts. Uses the PheWAS_manifest file as a guide to creating phenotypes and will only create phenotypes where the required data is availble, there is no need to edit the PheWAS_manifest if the user does not have the full dataset required. 

The default save location is located within the Deep_PheWAS files the user can choose to use different folder/file name. If specified folders will be created if not already present. The PheWAS_manifest which guides the phenotype creation can be chosen manually if it has been edited to add/change phenotypes, it is recommended to read the user guide if wanting to edit the PheWAS_manifest. 

Usage:
    formula_phenotypes.R (--min_data=<FILE>) [--data_field_phenotypes=<FILE> --sex_info=<FILE> --PheWAS_manifest_overide=<FILE> --phenotype_save_file=<FILE>]
    
Options:
    -h --help                                 Show this screen.
    -v --version                              Show version.
    
    Mandated inputs
    --min_data=<FILE>                         Full file path of the min_data file for UKB, see minimum_data script for format.
  
    Options
    --data_field_phenotypes=<FILE>            Full file path of the data_field.RDS file produced by data_field_phenotypes.R, only used for eGFR phenotype.
    --sex_info=<FILE>                         Full file path of the combined_sex file a file containing sex information (0=female, 
                                              1=male). Used only for eGFR field_ID phenotype.
    --PheWAS_manifest_overide=<FILE>          Full file path of the alternative PheWAS_manifest file.
    --phenotype_save_file=<FILE>              Full file path for the save file for the generated RDS used for phenotype creation. 
                                              [default: data/phenotypes/formula_phenotypes.RDS]

' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.1 formula_phenotypes.R')

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(parallel))
suppressMessages(library(lubridate))
library(here)

# Load in defining variables -----------------------------------------------------
# data_field_phenotypes
if(!is.null(arguments$data_field_phenotypes)){
data_field_phenotypes <- readRDS(arguments$data_field_phenotypes)
current_IDs <- names(data_field_phenotypes)
} else {
  current_IDs <- c(NA)
}
# minimum data
min_data <- fread(arguments$min_data)
current_field_id <- gsub("-.*", "\\1", colnames(min_data))
# PheWAS manifest
if(is.null(arguments$PheWAS_manifest_overide)){
  PheWAS_manifest <- fread(here("data","PheWAS_manifest.csv"))
} else {
  PheWAS_manifest <- fread(arguments$PheWAS_manifest_overide)
}
# save location
if(arguments$phenotype_save_file=="data/phenotypes/formula_phenotypes.RDS"){
  phenotype_save_location <- here("data","phenotypes","formula_phenotypes.RDS")
  if(!dir.exists(here("data","phenotypes"))){
    dir.create(here("data","phenotypes"))
  }
} else {
  phenotype_save_location <- arguments$phenotype_save_file
  new_folder <- str_remove(arguments$phenotype_save_file,basename(arguments$phenotype_save_file))
  if(!dir.exists(new_folder)){
    dir.create(new_folder)}
}
# empty list
formula_phenotypes <- list()
# eGFR --------------------------------------------------------------------
# eGFR uses cystatin field ID but a slightly more involved formula so is sectioned off here, will only run if the appropriate 
# PheWAS_ID is available from the data
if ("Q0800" %in% PheWAS_manifest$PheWAS_ID & "Q0513" %in% current_IDs & !is.null(arguments$sex_info)) {
  # Subset biomarker data to baseline cystatin-C measure to calculate eGFR taken from here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4398023/
  # sex info
  # 0= female, 1=male
  sex <- fread(arguments$sex_info)
  eGFR <- data_field_phenotypes[c("Q0513")] %>% 
    reduce(full_join,by="eid") %>% 
    rename(cystatin=Q0513,age=Q0513_age) %>% 
    inner_join(sex) %>% 
    drop_na(sex,cystatin) %>% 
    mutate(eGFR=ifelse(cystatin>0.8 & sex==0,
                       (133*(cystatin/0.8)^-1.328*(0.996^age)*0.932),ifelse(
                         cystatin>0.8 & sex==1,
                         (133*(cystatin/0.8)^-1.328*(0.996^age)),ifelse(
                           cystatin<=0.8 & sex==0,
                           (133*(cystatin/0.8)^-0.499*(0.996^age)*0.932),(133*(cystatin/0.8)^-0.499*(0.996^age)*0.932))))) %>% 
    select(eid,Q2020=eGFR,earliest_date,Q2020_age=age)
  eGFR_list <- list(Q2020=eGFR)
  formula_phenotypes <- append(formula_phenotypes,eGFR_list)
} 
# QTc ---------------------------------------------------------------------
if("21003" %in% current_field_id & "22333" %in% current_field_id & "53" %in% current_field_id & "22331" %in% current_field_id) {
  ## QTc
  age <- min_data %>% 
    select(eid,matches("^21003-")) %>% 
    select(1,age=2) %>% 
    drop_na()
  date <- min_data %>% 
    select(eid,matches("^53-")) %>% 
    select(1,earliest_date=2) %>% 
    drop_na() %>% 
    mutate(earliest_date=ymd(earliest_date))
  QT <- min_data %>% 
    select(eid,matches("^22331-")) %>% 
    select(1,QT=2) %>% 
    drop_na()
  RR <- min_data %>% 
    select(eid,matches("^22333-")) %>% 
    select(1,RR=2) %>% 
    drop_na()
  QTc <- QT %>% 
    left_join(RR, by="eid") %>% 
    left_join(age, by="eid") %>% 
    left_join(date,by = "eid") %>% 
    mutate(qtc_interval_bazett = QT/sqrt(RR/1000)) %>% 
    select(eid,Q1500=qtc_interval_bazett,earliest_date,Q1500_age=age)
  QTc_list <- list(Q1500=QTc)
  formula_phenotypes <- append(formula_phenotypes,QTc_list)  
}
# WHR ---------------------------------------------------------------------
if("21003" %in% current_field_id & "48" %in% current_field_id & "53" %in% current_field_id & "49" %in% current_field_id){
  ## WHR
  age <- min_data %>% 
    select(eid,matches("^21003-")) %>% 
    select(1,age=2) %>% 
    drop_na()
  date <- min_data %>% 
    select(eid,matches("^53-")) %>% 
    select(1,earliest_date=2) %>% 
    drop_na() %>% 
    mutate(earliest_date=ymd(earliest_date))
  waist <- min_data %>% 
    select(eid,matches("^48-")) %>% 
    select(1,waist=2) %>% 
    drop_na()
  hip <- min_data %>% 
    select(eid,matches("^49-")) %>% 
    select(1,hip=2) %>% 
    drop_na()
  WHR <- waist %>% 
    left_join(hip,by="eid") %>% 
    left_join(age, by="eid") %>% 
    left_join(date,by = "eid") %>% 
    mutate(WHR=waist / hip) %>% 
    select(eid,Q0007=WHR,earliest_date,Q0007_age=age)
  WHR_list <- list(Q0007=WHR)
  formula_phenotypes <- append(formula_phenotypes,WHR_list)
}

# LF ratio ---------------------------------------------------------------
if("53" %in% current_field_id & "21003" %in% current_field_id & "3062" %in% current_field_id & "20151" %in% current_field_id & "3063" %in% current_field_id & "20150" %in% current_field_id){
# need date of assessment 
date_assesment <- min_data %>% 
  select(eid,matches("^53-"))
date_cols <-  colnames(date_assesment)[-1]
# and age at assessment
age_assesment <-  min_data %>% 
  select(eid,matches("^21003-")) 
age_cols <- colnames(age_assesment)[-1]
# FVC all measure and single best measure
FVC_all <- min_data %>% 
  select(eid,matches("^3062-"))
FVC_best <- min_data %>% 
  select(eid,matches("^20151-")) %>% 
  drop_na() %>% 
  rename(FVC=2)
# col names used to combine values
data_field_FVC <- colnames(FVC_all)[-1]
# find date position for best FVC value
FVC_mod <- FVC_all %>% 
  mutate(results_list=transpose(across(all_of(data_field_FVC)))) %>% 
  left_join(FVC_best) %>% 
  drop_na(FVC) %>% 
  mutate(position_best=unlist(map2(results_list,FVC,~match(.y,unlist(.x)))),
         matching_date_position_first=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_FVC[position_best])))+1,) %>% 
  select(eid,FVC,FVC_position=matching_date_position_first)
# FEV1 all measures and best measure
FEV1_all <- min_data %>% 
  select(eid,matches("^3063-"))
FEV1_best <- min_data %>% 
  select(eid,matches("^20150-")) %>% 
  drop_na() %>% 
  rename(FEV1=2)
# col names used for combining
data_field_FEV1 <- colnames(FEV1_all)[-1]
# finding date/age position
FEV1_mod <- FEV1_all %>% 
  mutate(results_list=transpose(across(all_of(data_field_FEV1)))) %>% 
  left_join(FEV1_best) %>% 
  drop_na(FEV1) %>% 
  mutate(position_best=unlist(map2(results_list,FEV1,~match(.y,unlist(.x)))),
         matching_date_position_first=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_FEV1[position_best])))+1,) %>% 
  select(eid,FEV1,FEV1_position=matching_date_position_first)
# making ratio and finding date and age based on shared position
ratio <- FEV1_mod %>% 
  left_join(FVC_mod, by="eid") %>% 
  mutate(divergence=ifelse(FEV1_position==FVC_position,1,0)) %>% 
  filter(divergence==1) %>% 
  mutate(ratio=FEV1/FVC) %>% 
  select(eid,ratio,position=FEV1_position) %>% 
  left_join(date_assesment, by="eid") %>% 
  mutate(across(all_of(date_cols), ~gsub(" .*", "", .x))) %>% 
  mutate(across(all_of(date_cols),as.character)) %>% 
  mutate(date_list = transpose(across(all_of(date_cols))),
         earliest_date=map2(date_list,position, function(x,y) x[[y]])) %>% 
  left_join(age_assesment, by="eid") %>% 
  mutate(age_list=transpose(across(all_of(age_cols))),
         pheno_age = unlist(map2(age_list,position, function(x,y) x[[y]]))) %>% 
  select(eid,ratio,earliest_date,pheno_age) %>% 
  mutate(earliest_date=ymd(earliest_date)) %>% 
  set_names(c("eid","Q0123","earliest_date","Q0123_age"))
# make list and append
ratio_list <- list(Q0123=ratio)
formula_phenotypes <- append(formula_phenotypes,ratio_list)
}
saveRDS(formula_phenotypes,phenotype_save_location)