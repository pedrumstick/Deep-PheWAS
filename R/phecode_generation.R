#!/usr/bin/env Rscript
'Uses pre-loaded mapping files to first map health care data from ICD10 and ICD9 code lists to phecodes and separately for use in other phenotypes the range_IDs (exclusions from being controls). Original mapping files can be found https://phewascatalog.org/phecodes. 

User inputs location of health_data and sex_info files (see user guide for format) to run with default settings. 

If using user generated files/data then the user can control the labels assigned for ICD10 and ICD9 data in health_data using the ICD10 and ICD9 inputs. If wanting only phecodes or range_IDs select no_range_ID and no_phecode respectively, the save locations of the range_ID and phecode list object can also be inputted. As for other scripts the user has the option for parallel computing by inputting a value for N_codes, see user guide for more detail on required file formatting.

Usage:
    phecode_generation.R (--health_data=<FILE> --sex_info=<FILE>) [--control_exclusions=<FILE> --N_cores=<number> --no_phecodes --no_range_ID --ICD10=<text_comma> --ICD9=<text_comma> --phecode_save_file=<FILE> --range_ID_save_file=<FILE>]
    
Options:
    -h --help                                 Show this screen.
    -v --version                              Show version.
    
    Mandated inputs
    --health_data=<FILE>                      Full file path of the tabdata file for UKB, see user guide for file format.
    
    --sex_info=<FILE>                         Full file path of the combined_sex file a file containing participant ID and sex information 
                                              (0=female, 1=male).
                                              
    Options
    --control_exclusions=<FILE>               Full file path of the optional control exclusions file.
    --N_cores=<number>                        Number of cores requested if wanting to use parallel computing.
    --no_range_ID                             Select if not wanting to run the range_ID function, will save all IDs for each range_ID 
                                              described by the phecodes ~233 range IDs.
                                              
    --no_phecodes                             Select if not wanting to extract the phecodes as phenotypes for analysis. 
    
    --ICD10=<text_comma>                      Comma separated text representing the labels used in health_data for ICD10 values if 
                                              no ICD10 values use NA or any text not used as a source in health data. If not used default 
                                              for ICD10 (ICD10,cancer,MD) will be used.
                                              
    --ICD9=<text_comma>                       Comma separated text representing the labels used in health_data for ICD9 values if no
                                              ICD9 value use NA or any text not used as a source in health data. If not used default 
                                              for ICD9 (ICD9) will be used.
                                              
    --phecode_save_file=<FILE>                Full file path for the save file for the phecode phenotypes list object. 
                                              [default: data/phenotypes/phecodes.RDS]
    --range_ID_save_file=<FILE>               Full file path to the folder used to store the range_ID list object. 
                                              [default: data/phenotypes/range_ID.RDS]
                                              
    
' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.1 pheocode_generation.R')

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(parallel))
suppressMessages(library(lubridate))
library(here)

# Functions ---------------------------------------------------------------
phecoding <- function(a,b,c,d,e,f,g,h,i) {
  # j is sex, a small number of phecodes are sex specific and so the function needs to account for this.  
  if (is.na(i)) {
    # searching for data
    phecode <- phecodes_mapped %>%
      filter(phecode==b)
    # making cases
    cases_epi_date <- phecode %>% 
      group_by(eid) %>% 
      summarise(earliest_date=min(date))
    cases_epi <- phecode %>%
      count(eid,phecode) %>% 
      left_join(cases_epi_date) %>% 
      mutate(!!a := n) %>% 
      select(eid,{{a}},earliest_date)
    case_id_epi <- cases_epi %>% 
      select(1) %>% 
      pull()
    #excluded are those in the exclusion ranges so these cannot be controls.
    excluded_epi_df <- phecodes_mapped %>%
      filter(!eid %in% case_id_epi) %>% 
      filter(between(phecode,c,d) | between(phecode,e,f) | between(phecode,g,h)) %>% 
      mutate(!!a := NA,
             earliest_date=NA,
             earliest_date=ymd(earliest_date)) %>%  
      distinct(eid, .keep_all = T) %>% 
      select(eid,{{a}},earliest_date)
    # vector to use later
    excluded_epi <- excluded_epi_df %>% 
      select(eid) %>% 
      pull()
    #controls are then selected from the remaining population sex label has virtually everyone in UK Biobank and so is the population 
    #sample for the epidemiological definition
    controls_epi <- combined_sex %>% 
      filter(!eid %in% case_id_epi) %>% 
      filter(!eid %in% excluded_epi) %>% 
      filter(!eid %in% control_exclusions) %>%
      mutate(!!a := 0,
             earliest_date =NA,
             earliest_date=ymd(earliest_date)) %>% 
      select(eid,{{a}},earliest_date)
    # final output
    phecode_status <- cases_epi %>% 
      bind_rows(controls_epi, excluded_epi_df)
    return(phecode_status)
  } else {
    #exactly as above but filtering for sex
    phecode <- phecodes_mapped %>%
      filter(sex==i) %>%
      filter(phecode==b)
    # making cases
    cases_epi_date <- phecode %>% 
      group_by(eid) %>% 
      summarise(earliest_date=min(date))
    cases_epi <- phecode %>%
      count(eid,phecode) %>% 
      left_join(cases_epi_date) %>% 
      mutate(!!a := n) %>% 
      select(eid,{{a}},earliest_date)
    # vector used later
    case_id_epi <- cases_epi %>% 
      select(1) %>% 
      pull()
    # now exclusions
    excluded_epi_df <- phecodes_mapped %>%
      filter(sex==i) %>% 
      filter(between(phecode,c,d) | between(phecode,e,f) | between(phecode,g,h)) %>% 
      mutate(!!a := NA,
             earliest_date=NA,
             earliest_date=ymd(earliest_date)) %>%  
      distinct(eid, .keep_all = T) %>% 
      select(eid,{{a}},earliest_date)
    # vector used later
    excluded_epi <- excluded_epi_df %>% 
      select(eid) %>% 
      pull()
    #controls are then selected from the remaining population sex label has virtually everyone in UK Biobank and so is the population 
    #sample for the epidemiological definition
    controls_epi <- combined_sex %>%  
      filter(sex==i) %>% 
      filter(!eid %in% case_id_epi) %>% 
      filter(!eid %in% excluded_epi) %>% 
      filter(!eid %in% control_exclusions) %>%
      mutate(!!a := 0,
             earliest_date =NA,
             earliest_date=ymd(earliest_date)) %>% 
      select(eid,{{a}},earliest_date)
    # final output
    phecode_status <- cases_epi %>% 
      bind_rows(controls_epi, excluded_epi_df)
    return(phecode_status)
  }
}  
exclusion_generator <- function(a,b,c,d,e,f,g) {
  #use the largest data file which includes the mortality and cancer data.
  exclusions <- phecodes_mapped %>% 
    filter(between(phecode,b,c) | between(phecode,d,e) | between(phecode,f,g)) %>% 
    mutate(!!a := 1) %>% 
    distinct(eid)
  return(exclusions)
}
# Load in and define variables -----------------------------------------------------------------
# exclusions
if(!is.null(arguments$exclusions)) {
  exclusions <- fread(arguments$exclusions, header = F) %>% 
    pull()
} else {
  exclusions <- data.frame(V1=NA) %>% 
    pull()
}
if(!is.null(arguments$N_cores)) {
  N_cores <- as.numeric(arguments$N_cores)
} else {
  N_cores <- NA
}
if(arguments$phecode_save_file=="data/phenotypes/phecodes.RDS") {
  phecode_save_location <- here("data","phenotypes","phecodes.RDS")
  if(!dir.exists(here("data","phenotypes"))){
    dir.create(here("data","phenotypes"))
  }
} else {
  phecode_save_location <- arguments$phecode_save_file
  new_folder <- str_remove(arguments$phecode_save_file,basename(arguments$phecode_save_file))
  if(!dir.exists(new_folder)){
    dir.create(new_folder)
  } 
}
if(arguments$range_ID_save_file=="data/phenotypes/range_ID.RDS"){
  range_ID_save_location <- here("data","phenotypes","range_ID.RDS")
  if(!dir.exists(here("data","phenotypes"))){
    dir.create(here("data","phenotypes"))
  } 
} else {
  range_ID_save_location <- arguments$range_ID_save_file
  new_folder <- str_remove(arguments$range_ID_save_file,basename(arguments$range_ID_save_file))
  if(!dir.exists(new_folder)){
    dir.create(new_folder)}
}
# read in the phecode mapping files for ICD10 and ICD9
phecode_map_ICD10 <- fread(here("data","phecode_map_rollup_ICD10"))
phecode_map_ICD9 <- fread(here("data","phecode_map_rollup_ICD9"))

# phecode definitions for mapping to phecodes
phecode_definitions <- fread(here("data","phecode_definitions"))

# optional file for control exclusions due to lack of suitable follow-up
if(is.null(arguments$control_exclusions)) {
  control_exclusions <- data.frame(list("eid"="")) %>% 
    pull(eid)
} else {
  control_exclusions <- fread(arguments$control_exclusions) %>% 
    select(1) %>% 
    pull()
}
# cases and controls can only be defined where recorded sex is available
combined_sex <- fread(arguments$sex_info)
# health data
health_data <- fread(arguments$health_data)
# Map to hospital data to phecodes -----------------------------------------------------
# source names for ICD10 records and ICD9 records respectively
if(is.null(arguments$ICD10)) {
  ICD10 <- c("ICD10_1","ICD10_2","ICD10_3","MD","cancer")
} else {
  ICD10 <- c(unlist(strsplit(arguments$ICD10,",")))
}
if(is.null(arguments$ICD9)) {
  ICD9 <- c("ICD9_1","ICD9_2","ICD9_3")
} else {
  ICD9 <- c(unlist(strsplit(arguments$ICD9,",")))
}
# creating an IC9 and ICD10 specific records for searching
all_ICD10 <- health_data %>% 
  filter(source %in% ICD10)
all_ICD9 <- health_data %>% 
  filter(source %in% ICD9)
## map to the phecodes separately for the different data sources ICD10 and 9 (phecodes have different maps)
ICD10_mapped <- left_join(all_ICD10,phecode_map_ICD10, by=c("code"="ICD10")) %>% 
  select(eid,phecode=PHECODE,date) %>% 
  drop_na()
ICD9_mapped <- left_join(all_ICD9,phecode_map_ICD9, by=c("code"="icd9")) %>% 
  select(eid,phecode,date) %>% 
  drop_na()
## then join the ICD10 and ICD9 data together and add in sex information via combined sex file
phecodes_mapped <- bind_rows(ICD10_mapped,ICD9_mapped) %>% 
  left_join(combined_sex)
# Phecode phenotype function -----------------------------------------------------
# these are the variables that are used in the mapply function
phewas_ID <- phecode_definitions[["phecode_name"]]
code_of_interest <- phecode_definitions[["phecode"]]
l_1 <- phecode_definitions[["l1"]]
u_1 <- phecode_definitions[["u1"]]
l_2 <- phecode_definitions[["l2"]] 
u_2 <- phecode_definitions[["u2"]] 
l_3 <- phecode_definitions[["l3"]] 
u_3 <- phecode_definitions[["u3"]] 
sex <- phecode_definitions[["sex_num"]]

if (isFALSE(arguments$no_phecodes)) {
  if (is.numeric(N_cores)) {
    epi_phenotypes <- mcmapply(phecoding,
                               phewas_ID,code_of_interest,l_1,u_1,l_2,u_2,l_3,u_3,sex, 
                               SIMPLIFY = FALSE,USE.NAMES = T,mc.cores = N_cores) 
  } else {
    epi_phenotypes <- mapply(phecoding,
                             phewas_ID,code_of_interest,l_1,u_1,l_2,u_2,l_3,u_3,sex, 
                             SIMPLIFY = FALSE,USE.NAMES = T)  
  }
  # save
  saveRDS(epi_phenotypes,phecode_save_location)
}
# Exclusion range function ------------------------------------------------
# Now making a file that contains just the 233 groups of exclusions. This is then used in the primary care derived curated phenotypes where appropriate
if (isFALSE(arguments$no_range_ID)) {
  # this uses the previously defined phecode definitions as input
  exclude_definitions <- phecode_definitions %>% 
    group_by(range_ID) %>% 
    summarise(l1=max(l1),u1=max(u1),l2=max(l2),u2=max(u2),l3=max(l3),u3=max(u3))
  l1_e <- exclude_definitions[["l1"]]
  u1_e <- exclude_definitions[["u1"]]
  l2_e <- exclude_definitions[["l2"]]
  u2_e <- exclude_definitions[["u2"]]
  l3_e <- exclude_definitions[["l3"]]
  u3_e <- exclude_definitions[["u3"]]
  names_e <-exclude_definitions[["range_ID"]]
  # run
  range_exclusions <- mapply(exclusion_generator,names_e,l1_e,u1_e,l2_e,u2_e,l3_e,u3_e, SIMPLIFY = F, USE.NAMES = T)
  # create an empty exclusion range
  eid <- ""
  R_0 <- data.frame(eid)
  # append for final output saved as list object
  range_exclusions_edit <- append(range_exclusions,list(R_0=R_0))
  saveRDS(range_exclusions_edit,range_ID_save_location)
}
