#! /usr/bin/env Rscript
'Creates primary care quantitative phenotypes. Uses the PheWAS_manifest file as a guide to creating phenotypes. If using non UK-Biobank data the GP clinical data would need to match the format of UK Biobank data and use one or both of Read V2 or Read V3 codes. As for other scripts the user has the option for parallel computing by inputting a value for N_codes. 

The default save location is located within the Deep_PheWAS files the user can choose to use different folder/file name. If specified folders will be created if not already present. The PheWAS_manifest which guides the phenotype creation can be chosen manually if it has been edited to add/change phenotypes, it is recommended to read the user guide if wanting to edit the PheWAS_manifest. 

Usage:
    primary_care_quantitiative_phenotypes.R (--GPC=<FILE> --DOB=<FILE>) [--N_cores=<number> --PheWAS_manifest_overide=<FILE> --phenotype_save_file=<FILE>]
    
Options:
    -h --help                                 Show this screen.
    -v --version                              Show version.
    
    Mandated inputs
    --GPC=<FILE>                              Full file path of the GP clinical data file from UKB.
    --DOB=<FILE>                              Full file path of date of birth file. File containing eid and DOB columns. DOB does not have to be exact day
                                              Default DOB file created for UK Biobank data from month of birth (field-ID 52), and year of birth (field-ID 34), 
                                              to create a pseudo date of birth.
    Options
    --phenotype_save_file=<FILE>              Full file path for the save filefor the generated RDS used for phenotype creation. 
                                              [default: data/phenotypes/PQP.RDS]
    --N_cores=<number>                        Number of cores requested if wanting to use parallel computing.
    --PheWAS_manifest_overide=<FILE>          Full file path of the alternative PheWAS_manifest file.
  
' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.1 primary_care_quantitiative_phenotypes.R')

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(parallel))
suppressMessages(library(lubridate))
library(here)

# functions ---------------------------------------------------------------
quant_primary_care <- function(a,b,c,d,e) {
  # message to track progress
  message(a)
  # create an age column for that phenotype
  PheWAS_ID_age <- paste0(a,"_age")
  # read in codes from list
  codes <-  fread(here("data","PQP_codes",b))
  if(c==0) {
    # splitting by Read code Version
    read_V2 <- codes %>% 
      filter(type=="read_V2") %>% 
      pull(read_code)
    read_V3 <- codes %>% 
      filter(type=="read_V3") %>% 
      pull(read_code)
    # extracting values
    prim_care_values <- prim_care %>%
      filter(read_2 %in% read_V2 | read_3 %in% read_V3) %>%
      group_by(eid) %>%
      filter(value1 != "",
             !grepl("^OPR|\\^",value1),
             !grepl("^OPR|\\^",value2),
             !grepl("^OPR|\\^",value3)) %>% 
      mutate(date=event_dt,
             value1=as.numeric(value1)) %>%
      arrange(date) %>%
      slice(1) %>%
      ungroup() %>% 
      left_join(DOB) %>% 
      mutate(age=round(time_length(x = difftime(date,DOB),unit = "years"),digits = 0)) %>% 
      select(eid,value1,age) %>% 
      set_names(c("eid",a,PheWAS_ID_age))
    return(prim_care_values)
  } else if (c==1) {
    # splitting by Read code Version
    read_V2 <- codes %>% 
      filter(type=="read_V2") %>% 
      pull(read_code)
    read_V3 <- codes %>% 
      filter(type=="read_V3") %>% 
      pull(read_code)
    # extracting values
    prim_care_values <- prim_care %>%
      filter(read_2 %in% read_V2 | read_3 %in% read_V3) %>%
      group_by(eid) %>%
      filter(value1 != "",
             !grepl("^OPR|\\^",value1),
             !grepl("^OPR|\\^",value2),
             !grepl("^OPR|\\^",value3)) %>% 
      mutate(date=event_dt,
             value1=as.numeric(value1)) %>% 
      filter(value1>d,value1<e) %>% 
      arrange(date) %>%
      slice(1) %>%
      ungroup() %>% 
      left_join(DOB) %>% 
      mutate(age=round(time_length(x = difftime(date,DOB),unit = "years"),digits = 0)) %>% 
      select(eid,value1,age) %>% 
      set_names(c("eid",a,PheWAS_ID_age))
    return(prim_care_values)
  } 
}
# Load in defining variables -----------------------------------------------------
# primary care data
prim_care <- fread(arguments$GPC)
# DOB
DOB <- fread(arguments$DOB)
# PheWAS manifest
if(is.null(arguments$PheWAS_manifest_overide)){
  PheWAS_manifest <- fread(here("data","PheWAS_manifest.csv"))
} else {
  PheWAS_manifest <- fread(arguments$PheWAS_manifest_overide)
}
# N_Cores
if(!is.null(arguments$N_cores)) {
  N_cores <- as.numeric(arguments$N_cores)
} else {
  N_cores <- NA
}
# save location
if(arguments$phenotype_save_file=="data/phenotypes/PQP.RDS"){
  phenotype_save_location <- here("data","phenotypes","PQP.RDS")
  if(!dir.exists(here("data","phenotypes"))){
    dir.create(here("data","phenotypes"))
  }
} else {
  phenotype_save_location <- arguments$phenotype_save_file
  new_folder <- str_remove(arguments$phenotype_save_file,basename(arguments$phenotype_save_file))
  if(!dir.exists(new_folder)){
    dir.create(new_folder)}
}
# selects phenotypes based on the phenotype manifest.
primary_care_phewas_ID_info <- PheWAS_manifest %>% 
  filter(category=="primary_care_quantitative_phenotype")
# Run functions -----------------------------------------------------------
if (is.numeric(N_cores)) {
  primary_care_quantitiative_data <- mcmapply(quant_primary_care,
                                              primary_care_phewas_ID_info$PheWAS_ID,
                                              primary_care_phewas_ID_info$primary_care_code_list,
                                              primary_care_phewas_ID_info$limits,
                                              primary_care_phewas_ID_info$lower_limit,
                                              primary_care_phewas_ID_info$upper_limit,
                                              SIMPLIFY = F,mc.cores = N_cores,USE.NAMES = T)
} else {
  primary_care_quantitiative_data <- mapply(quant_primary_care,
                                            primary_care_phewas_ID_info$PheWAS_ID,
                                            primary_care_phewas_ID_info$primary_care_code_list,
                                            primary_care_phewas_ID_info$limits,
                                            primary_care_phewas_ID_info$lower_limit,
                                            primary_care_phewas_ID_info$upper_limit,
                                            SIMPLIFY = F,USE.NAMES = T)
  
}
saveRDS(primary_care_quantitiative_data,phenotype_save_location)