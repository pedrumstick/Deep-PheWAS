#! /usr/bin/env Rscript
'Combines tab data sourced from UK-Biobank into a single minimum_tab_data.gz file that is used throughout the Deep-PheWAS pipeline. UK-Biobank field-ids can be added as required by editing the fields-minium.txt file found by default in data/fields-minimum.txt. Can process tab data in the following format option used by UK-Biobank r,csv or txt default is either txt or csv, to select r version use optional flag R_format. Further details of the final file produced, and inputs can be found in the user guide. As for other scripts the user has the option for parallel computing by inputting a value for N_codes.

Usage:
    minimum_tab_data_edit.R (--data_folder=<FOLDER> | --data_files=<FILES>) [--r_format --data_field_ID=<FILE> --data_name_pattern=<text> --N_cores=<number> --save_loc=<FILE> --exclusions=<FILE>]
    
Options:
    -h --help                     Show this screen.
    -v --version                  Show version.
    
    Mandated inputs
    --data_folder=<FOLDER>        Location of all the tab data files that will be searched.
    Or
    --data_files=<FILES>          Comma seperated input of full file path to data files.
    
    Options
    --r_format                    Input if data from UK-Biobank has been downloaded using the r option in ukbconv. If using the txt or csv option no 
                                  input is required (effective default).
    
    --data_field_ID=<FILE>        Full path to the file containing the field_IDs required for Deep-PheWAS if not using default file. Is a plain text
                                  file with no header one field-ID per row. Field-ID is a numeric value, example field-ID 54 is UK biobank assesment
                                  centre.
                                  [default: fields-minimum.txt]
    
    --data_name_pattern=<text>    Common pattern for isolating tab_data files if in directory with mixed files.
    --N_cores=<number>            Number of cores requested if wanting to use parallel computing.
    --save_loc=<FILE>             Full path to save file location for minimum_tab_data, defaults to minimum_tab_data.gz located in 
                                  the data/ folder. [default: minimum_tab_data.gz]
                                  
    --exclusions=<FILE>           Full path to an exclusions file (single column of IDs) that documents any IDs to be excluded from the analysis.
    
' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.1 minimum_data.R')

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(parallel))
library(here)

# Functions ---------------------------------------------------------------
# simple function that searches for the tab names in the data and then extracts them as a list
tab_data_clean <- function (x) {
  # finding column names in data set 
  tabdata_colname <- colnames(fread(x, nrows=0))
  colname_search <- str_subset(tabdata_colname,paste0(headings,collapse = "|"))
  colname_search <- colname_search[!is.na(colname_search)]
  # and extracting only required columns
  tabdata <- fread(x,select = colname_search) %>% 
    rename(eid=1) %>% 
    filter(!eid %in% exclusions$V1)
  return(tabdata)
}
# Load in and define variables -----------------------------------------------------------------
# N_cores
if(!is.null(arguments$N_cores)) {
  N_cores <- as.numeric(arguments$N_cores)
} else {
  N_cores <- NA
}
# required data_fields
if(arguments$data_field_ID!="fields-minimum.txt") {
  field_ID_headings <- fread(arguments$data_field_ID, header = F) } else {
    field_ID_headings <- fread(here("data","fields-minimum.txt"), header = F)
  }
adding_eid <- data.frame(conversion=c("eid","f.eid"))
# selecting data format for UK-Biobank only
if(arguments$r_format){
  headings <- field_ID_headings %>% 
    mutate(conversion=paste0("f.",V1,".")) %>% 
    bind_rows(adding_eid) %>% 
    pull(conversion)
} else {
  headings <- field_ID_headings %>% 
    mutate(conversion=paste0("^",V1,"-")) %>%
    bind_rows(adding_eid) %>%
    pull(conversion)
}
# exclusions
if(!is.null(arguments$exclusions)) {
  exclusions <- fread(arguments$exclusions, header = F)
} else {
  exclusions <- data.frame(V1=NA)
}
# selecting files to search
if(!is.null(arguments$data_folder)) {
if(is.null(arguments$data_name_pattern)) {
  tab_files <- list.files(path = arguments$data_folder, full.names = T)
} else {
  tab_files <- list.files(path = arguments$data_folder, pattern = arguments$data_name_pattern, full.names = T)
}
} else if(!is.null(arguments$data_files)){
  tab_files <- unlist(str_split(arguments$data_files,pattern = ","))
}
# Run functions -------------------------------------------------------
if(is.numeric(N_cores)){
  min_tab_data <- mclapply(tab_files, tab_data_clean,mc.cores = N_cores)
} else {
  min_tab_data <- mapply(tab_data_clean,tab_files, SIMPLIFY = F)
}
# Use purr to map and join Worst case scenario numbers don't match but it jst introduces NA's.
  minimum_tab <- min_tab_data %>% 
    reduce(full_join, by="eid")
  if(arguments$r_format) {
  colnames(minimum_tab) <- str_replace(colnames(minimum_tab),"f\\.","")
  colnames(minimum_tab) <- str_replace(colnames(minimum_tab),"\\.","-")
  }
  if(arguments$save_loc=="minimum_tab_data.gz") {
    fwrite(minimum_tab, here("data","minimum_tab_data.gz"), sep = "\t", na = NA)
  } else {
    fwrite(minimum_tab, arguments$save_loc, sep = "\t", na = NA)
  }
