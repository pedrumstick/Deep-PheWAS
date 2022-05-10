#! /usr/bin/env Rscript
'Creates Data Field phenotypes from UK-Biobank data_fields. Uses the PheWAS_manifest file as a guide to creating phenotypes and will only create phenotypes where the required data_field data is available, there is no need to edit the PheWAS_manifest if the user does not have the full dataset required. If using non UK-Biobank data data_fields would need to be mapped if similar measures were available. 

Only required data is min_data. As for other scripts the user has the option for parallel computing by inputting a value for N_codes. 

The default save location is located within the Deep_PheWAS files the user can choose to use different folder/file name. If specified folders will be created if not already present. The PheWAS_manifest which guides the phenotype creation can be chosen manually if it has been edited to add/change phenotypes, it is recommended to read the user guide if wanting to edit the PheWAS_manifest. 

Usage:
    data_field_phenotypes.R (--min_data=<FILE>) [--N_cores=<number> --PheWAS_manifest_overide=<FILE> --phenotype_save_file=<FILE> --append_file=<FILE>]
    
Options:
    -h --help                                 Show this screen.
    -v --version                              Show version.
    
    Mandated inputs
    --min_data=<FILE>                         Full file path of the tabdata file for UKB, see minimum_data script for format.
  
    Options
    --N_cores=<number>                        Number of cores requested if wanting to use parallel computing.
    --PheWAS_manifest_overide=<FILE>          Full file path of the alternative PheWAS_manifest file.
    --phenotype_save_file=<FILE>              Full file path for the save filefor the generated RDS used for phenotype creation. 
                                              [default: data/phenotypes/data_field_phenotypes.RDS]
    --append_file=<FILE>                      Full file path of RDS object containing complete field_ID phenotypes, will read object in and compare complete                                                   field_ID phenotypes saved in the object with those planned for devlopment. It will then only attempt to create
                                              phenotypes not present in the RDS object. This allows phenotypes to be added into the larger file without re-doing                                               the process for all field_ID phenotypes.
' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.1 data_field_phenotypes.R')

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
data_field_extraction <- function (a,b,c,d,e,f,g,h,i,j,k,l,m,n) {
  # a slightly different method is needed when considering traits that have a separate QC flag, we need to be able to 
  # remove the measurement that has been flagged and replace with an NA before we evaluate the measurements.
  message(a)
  if(l=="first"){
    position_col <- c("first_result","matching_date_position_first")
  } else if(l=="last"){
    position_col <- c("last_result","matching_date_position_last")
  } else if(l=="min"){
    position_col <- c("min_result","matching_date_position_min")
  } else if(l=="max"){
    position_col <- c("max_result","matching_date_position_max")
  }
  
  age_name <- paste0(a,"_age")
  
  if(d=="binary") {
    # extract columns for phenotype
    per_data_field <- min_data %>% 
      select(eid,matches(c))
    # vector used later
    data_field_colnames <- colnames(per_data_field)[-1]
    # extract date column
    date_col <- min_data %>% 
      select(eid,matches(j))
    # vector used later
    date_col_names <- colnames(date_col)[-1]
    # need to extract cases and controls separately, to achieve this need to identify unique values and set all non-cases or non-controls as NA
    all_values <- unique(as.vector(as.matrix(per_data_field[,-1])))
    case_values <- as.numeric(unlist(strsplit(e,",")))
    control_values <- as.numeric(unlist(strsplit(f,",")))
    # leave only values that represent cases
    na_values_controls <- all_values[- which(all_values %in% case_values)]
    # and controls
    na_values_cases <- all_values[- which(all_values %in% control_values)]
    # turn non-cases/controls into NA  
    just_cases <- per_data_field %>% 
      mutate(across(everything(),~replace(., . %in% na_values_controls,NA)))
    just_controls <- per_data_field %>% 
      mutate(across(everything(),~replace(., . %in% na_values_cases,NA)))
    # create cases
    cases <- just_cases %>% 
      unite(col=results,2:length(colnames(.)),sep = ",",na.rm = T,remove = F) %>% 
      mutate(results=na_if(results,"")) %>% 
      drop_na(results) %>% 
      mutate(results_list=transpose(across(all_of(data_field_colnames))))
    # too computationally inefficent to create all variations so crete only the one that is required
    if(l=="first"){
      cases_extracted <- cases %>%
        mutate(first_result= unlist(map(results_list, ~ first(na.omit(.x)))),
               position_first=unlist(map2(results_list,first_result, ~match(.y,unlist(.x)))),
               matching_date_position_first=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_first])))+1,)
    } else if(l=="last"){
      cases_extracted <- cases %>%
        mutate(last_result= unlist(map(results_list, ~ last(na.omit(.x)))),
               position_last=unlist(map2(results_list,last_result, ~match(.y,unlist(.x)))),
               matching_date_position_last=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_last])))+1,)
    } else if(l=="min"){
      cases_extracted <- cases %>%
        mutate(min_result= unlist(map(results_list, ~ min(unlist(.x),na.rm = T))),
               position_min=unlist(map2(results_list,min_result, ~match(.y,unlist(.x)))),
               matching_date_position_min=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_min])))+1,)
    } else if(l=="max"){
      cases_extracted <- cases %>%
        mutate(max_result= unlist(map(results_list, ~ max(unlist(.x),na.rm = T))),
               position_max=unlist(map2(results_list,max_result, ~match(.y,unlist(.x)))),
               matching_date_position_max=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_max])))+1,)
    }
    cases_extracted_complete <- cases_extracted %>% 
      select(eid,all_of(position_col)) %>% 
      rename(eid=1,pheno=2,position_match=3) %>% 
      mutate(pheno=1) %>% 
      drop_na() %>% 
      left_join(date_col, by="eid") %>% 
      mutate(across(all_of(date_col_names), ~gsub(" .*", "", .x))) %>% 
      mutate(across(all_of(date_col_names),as.character)) %>% 
      mutate(date_list = transpose(across(all_of(date_col_names))),
             earliest_date=map2(date_list,position_match, function(x,y) x[[y]])) %>%  
      select(eid,pheno,earliest_date) %>% 
      mutate(earliest_date=ymd(earliest_date),
             pheno=1) %>% 
      set_names(c("eid",a,"earliest_date"))
    # create controls
    controls <- just_controls %>% 
      unite(col=results,2:length(colnames(.)),sep = ",",na.rm = T,remove = F) %>% 
      mutate(results=na_if(results,"")) %>% 
      drop_na(results) %>% 
      mutate(results_list=transpose(across(all_of(data_field_colnames))))
    # too computationally inefficent to create all variations so crete only the one that is required
    if(l=="first"){
      controls_extracted <- controls %>%
        mutate(first_result= unlist(map(results_list, ~ first(na.omit(.x)))),
               position_first=unlist(map2(results_list,first_result, ~match(.y,unlist(.x)))),
               matching_date_position_first=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_first])))+1,)
    } else if(l=="last"){
      controls_extracted <- controls %>%
        mutate(last_result= unlist(map(results_list, ~ last(na.omit(.x)))),
               position_last=unlist(map2(results_list,last_result, ~match(.y,unlist(.x)))),
               matching_date_position_last=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_last])))+1,)
    } else if(l=="min"){
      controls_extracted <- controls %>%
        mutate(min_result= unlist(map(results_list, ~ min(unlist(.x),na.rm = T))),
               position_min=unlist(map2(results_list,min_result, ~match(.y,unlist(.x)))),
               matching_date_position_min=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_min])))+1,)
    } else if(l=="max"){
      controls_extracted <- controls %>%
        mutate(max_result= unlist(map(results_list, ~ max(unlist(.x),na.rm = T))),
               position_max=unlist(map2(results_list,max_result, ~match(.y,unlist(.x)))),
               matching_date_position_max=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_max])))+1,)
    }
    controls_extracted_complete <- controls_extracted %>% 
      select(eid,all_of(position_col)) %>% 
      rename(eid=1,pheno=2,position_match=3) %>% 
      mutate(pheno=1) %>% 
      drop_na() %>% 
      left_join(date_col, by="eid") %>% 
      mutate(across(all_of(date_col_names), ~gsub(" .*", "", .x))) %>% 
      mutate(across(all_of(date_col_names),as.character)) %>% 
      mutate(date_list = transpose(across(all_of(date_col_names))),
             earliest_date=map2(date_list,position_match, function(x,y) x[[y]])) %>%  
      select(eid,pheno,earliest_date) %>% 
      mutate(earliest_date=ymd(earliest_date),
             pheno=0) %>% 
      set_names(c("eid",a,"earliest_date")) 
    # union of cases and controls to ensure no repeat ID's
    union_cases_controls <- cases_extracted_complete %>% 
      inner_join(controls_extracted_complete, by="eid") %>% 
      pull(eid)
    if(m==1) {
      # exclude crossover cases/controls from being either cases or controls
      phenotype_extraction_final <- cases_extracted_complete %>% 
        bind_rows(controls_extracted_complete) %>% 
        filter(!eid %in% union_cases_controls)
    } else if(m==0){
      #exclude the crossover only from controls to ensure that there are no duplicate ids
      controls_modified <- controls_extracted_complete %>% 
        filter(!eid %in% union_cases_controls)
      phenotype_extraction_final <- cases_extracted_complete %>% 
        bind_rows(controls_modified)
    }
    return(phenotype_extraction_final)
  } else {
    if(is.na(b) || ncol(min_data %>% select(eid,matches(na.omit(b)))) <2){
      # extract columns for phenotype
      per_data_field <- min_data %>% 
        select(eid,matches(c))
      # vector used later
      data_field_colnames <- colnames(per_data_field)[-1]
      # extract date column
      date_col <- min_data %>% 
        select(eid,matches(j))
      # vector used later
      date_col_names <- colnames(date_col)[-1]
      # selecting single value for quantitative measure where multiple options are available,
      # joins with age at assessment to get accurate age when measure taken.
      phenotype_extraction <- per_data_field %>% 
        unite(col=results,2:length(colnames(.)),sep = ",",na.rm = T,remove = F) %>% 
        mutate(results=na_if(results,"")) %>% 
        drop_na(results) %>% 
        mutate(results_list = transpose(across(all_of(data_field_colnames))))
      
      # too computationally inefficent to create all variations so crete only the one that is required
      if(l=="first"){
        phenotype_extracted <- phenotype_extraction %>%
          mutate(first_result= unlist(map(results_list, ~ first(na.omit(.x)))),
                 position_first=unlist(map2(results_list,first_result, ~match(.y,unlist(.x)))),
                 matching_date_position_first=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_first])))+1,)
      } else if(l=="last"){
        phenotype_extracted <- phenotype_extraction %>%
          mutate(last_result= unlist(map(results_list, ~ last(na.omit(.x)))),
                 position_last=unlist(map2(results_list,last_result, ~match(.y,unlist(.x)))),
                 matching_date_position_last=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_last])))+1,)
      } else if(l=="min"){
        phenotype_extracted <- phenotype_extraction %>%
          mutate(min_result= unlist(map(results_list, ~ min(unlist(.x),na.rm = T))),
                 position_min=unlist(map2(results_list,min_result, ~match(.y,unlist(.x)))),
                 matching_date_position_min=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_min])))+1,)
      } else if(l=="max"){
        phenotype_extracted <- phenotype_extraction %>%
          mutate(max_result= unlist(map(results_list, ~ max(unlist(.x),na.rm = T))),
                 position_max=unlist(map2(results_list,max_result, ~match(.y,unlist(.x)))),
                 matching_date_position_max=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_max])))+1,)
      }
      phenotype_extracted_complete <- phenotype_extracted %>% 
        select(eid,all_of(position_col)) %>% 
        rename(eid=1,pheno=2,position_match=3) %>% 
        drop_na() 
      # create either date and age column or just date column depending on input
      if(is.na(k)){
        phenotype_extraction_final <- phenotype_extracted_complete %>% 
          left_join(date_col, by="eid") %>% 
          mutate(across(all_of(date_col_names), ~gsub(" .*", "", .x))) %>% 
          mutate(across(all_of(date_col_names),as.character)) %>% 
          mutate(date_list = transpose(across(all_of(date_col_names))),
                 earliest_date=map2(date_list,position_match, function(x,y) x[[y]])) %>% 
          select(eid,pheno,earliest_date) %>% 
          mutate(earliest_date=ymd(earliest_date)) %>% 
          set_names(c("eid",a,"earliest_date"))
        return(phenotype_extraction_final)
      } else {
        # extract age column
        age_col <- min_data %>% 
          select(eid,matches(k))
        # vector used later
        age_col_names <- colnames(age_col)[-1]
        # add in age specific column
        phenotype_extraction_final <- phenotype_extracted_complete %>% 
          left_join(date_col, by="eid") %>% 
          mutate(across(all_of(date_col_names), ~gsub(" .*", "", .x))) %>% 
          mutate(across(all_of(date_col_names),as.character)) %>% 
          left_join(age_col, by="eid") %>% 
          mutate(date_list=transpose(across(all_of(date_col_names))),
                 earliest_date = map2(date_list,position_match, function(x,y) x[[y]])) %>% 
          mutate(age_list=transpose(across(all_of(age_col_names))),
                 pheno_age = unlist(map2(age_list,position_match, function(x,y) x[[y]]))) %>% 
          select(eid,pheno,earliest_date,pheno_age) %>% 
          mutate(earliest_date=ymd(earliest_date)) %>% 
          set_names(c("eid",a,"earliest_date",age_name))
      }
      return(phenotype_extraction_final)
      
    } else {
      # creating a variable used for naming a column
      PheWAS_ID_age <- paste0(a,"_age")
      # extract date column
      date_col <- min_data %>% 
        select(eid,matches(j))
      # vector used later
      date_col_names <- colnames(date_col)[-1]
      # searching for all cols with the QC flag field ID
      flag_ID <- min_data %>% 
        select(eid,matches(b))
      # can select which QC flags to use as filters if not discriminating then use 'all' in PheWAS_manifest
      if(n=="all") {
        just_QC <- flag_ID
      } else {
        # need to extract cases and controls separately, to achieve this need to identify unique values and set all non-cases or non-controls as NA
        all_values <- unique(as.vector(as.matrix(flag_ID[,-1])))
        flag_values <- as.numeric(unlist(strsplit(n,",")))
        # leave only values that represent non-accepted QC values
        na_values <- all_values[- which(all_values %in% flag_values)]
        # turn non-accepted qc flags into NAs
        just_QC <- per_data_field %>% 
          mutate(across(everything(),~replace(., . %in% na_values,NA)))
      }
      # converts any non-empty or NA string to the value of 1 as it is not important for this analysis why a flag is present.
      just_QC[,2:ncol(just_QC)][just_QC[,2:ncol(just_QC)] != ""] <- 1
      # converts to numeric
      just_QC <- just_QC %>% 
        mutate(across(c(2:ncol(.)),as.numeric))
      
      # in some cases the QC flag has more measures than there are measures for the variable it is 
      # flagging so there are more QC flags than measures, to account for this the flag measures 
      # are combined into a single value per equivalent quantitative measure. They are then re-coded as the specific 
      # cause of the QC flag is not relevant for our analysis.
      flag_cols_vector <- colnames(flag_ID)[-1]
      # use the vector of flag ID col names and identify how many visits this represents
      unique_visit <- paste0(b,unique(as.numeric(gsub(".*-(.+)\\..*", "\\1", flag_cols_vector))))
      # then select each group of cols that represent each visit and convert into a single value if any value present
      split <- lapply(unique_visit, function (x) test <- just_QC %>% 
                        select(1,matches(x)))
      divide <- lapply(split, function (x) edit <- x %>% 
                         mutate(results = rowSums(across(-eid),na.rm = T)) %>%  
                         select(1,results)) 
      divided <- divide %>% 
        reduce(left_join,by=c("eid")) %>% 
        mutate(filter = rowSums(across(-eid),na.rm = T)) %>% 
        filter(filter>0) %>% 
        select(-filter)
      # convert into a list that should be in equal length to the data_field the flag is used for
      flag_edit <- divided %>% 
        mutate(across(all_of(colnames(.)[-1]), ~na_if(.,0))) %>% 
        mutate(results_list=transpose(across(all_of(colnames(.)[-1]))),
               flags=map(results_list, function(x) which(!is.na(unlist(x))))) %>% 
        select(eid,flags)
      # searching for all cols with the data_field
      per_data_field <- min_data %>% 
        select(eid,matches(c))
      # vector needed later
      data_field_colnames <- colnames(per_data_field)[-1]
      # edits the original values based on the presence of a flag ID and then selects the first valid value
      ID_edit <- per_data_field %>% 
        unite(col=results,2:length(colnames(.)),sep = ",",na.rm = T,remove = F) %>% 
        mutate(results=na_if(results,"")) %>% 
        drop_na(results) %>% 
        left_join(flag_edit) %>% 
        mutate(results_list=transpose(across(all_of(data_field_colnames))),
               new_results=map2(results_list,flags, function(x,y) replace(unlist(x),list = unlist(y),values = NA)))
      # too computationally inefficient to create all variations so create only the one that is required
      if(l=="first"){
        phenotype_extracted <- ID_edit %>%
          mutate(first_result= unlist(map(new_results, ~ first(na.omit(.x)))),
                 position_first=unlist(map2(new_results,first_result, ~match(.y,unlist(.x)))),
                 matching_date_position_first=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_first])))+1,)
      } else if(l=="last"){
        phenotype_extracted <- ID_edit %>%
          mutate(last_result= unlist(map(new_results, ~ last(na.omit(.x)))),
                 position_last=unlist(map2(new_results,last_result, ~match(.y,unlist(.x)))),
                 matching_date_position_last=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_last])))+1,)
      } else if(l=="min"){
        phenotype_extracted <- ID_edit %>%
          mutate(min_result= unlist(map(new_results, ~ min(unlist(.x),na.rm = T))),
                 position_min=unlist(map2(new_results,min_result, ~match(.y,unlist(.x)))),
                 matching_date_position_min=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_min])))+1,)
      } else if(l=="max"){
        phenotype_extracted <- ID_edit %>%
          mutate(max_result= unlist(map(new_results, ~ max(unlist(.x),na.rm = T))),
                 position_max=unlist(map2(new_results,max_result, ~match(.y,unlist(.x)))),
                 matching_date_position_max=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[position_max])))+1,)
      }
      phenotype_extracted_complete <- phenotype_extracted %>% 
        select(eid,all_of(position_col)) %>% 
        rename(eid=1,pheno=2,position_match=3) %>% 
        drop_na() 
      # used to add named_age or not as inputted
      if(is.na(k)){
        phenotype_extraction_final <- phenotype_extracted_complete %>% 
          left_join(date_col, by="eid") %>% 
          mutate(across(all_of(date_col_names), ~gsub(" .*", "", .x))) %>% 
          mutate(across(all_of(date_col_names),as.character)) %>% 
          mutate(date_list = transpose(across(all_of(date_col_names))),
                 earliest_date=map2(date_list,position_match, function(x,y) x[[y]])) %>% 
          select(eid,pheno,earliest_date) %>% 
          mutate(earliest_date=ymd(earliest_date)) %>% 
          set_names(c("eid",a,"earliest_date"))
        return(phenotype_extraction_final)
      } else {
        # extract age column
        age_col <- min_data %>% 
          select(eid,matches(k))
        # vector used later
        age_col_names <- colnames(age_col)[-1]
        # add in age specific column
        phenotype_extraction_final <- phenotype_extracted_complete %>% 
          left_join(date_col, by="eid") %>% 
          mutate(across(all_of(date_col_names), ~gsub(" .*", "", .x))) %>% 
          mutate(across(all_of(date_col_names),as.character)) %>% 
          left_join(age_col, by="eid") %>% 
          mutate(date_list=transpose(across(all_of(date_col_names))),
                 earliest_date = map2(date_list,position_match, function(x,y) x[[y]])) %>% 
          mutate(age_list=transpose(across(all_of(age_col_names))),
                 pheno_age = unlist(map2(age_list,position_match, function(x,y) x[[y]]))) %>% 
          select(eid,pheno,earliest_date,pheno_age) %>% 
          mutate(earliest_date=ymd(earliest_date)) %>% 
          set_names(c("eid",a,"earliest_date",age_name))
      }
      return(phenotype_extraction_final)
    }
  }
}
combining_data_field <- function (a,b,c,d,e,f) {
  message(a)
  # selecting min or max values
  if(b=="min"){
    position_col <- c("min_result","position_min")
  } else if(b=="max"){
    position_col <- c("max_result","position_max")
  }
  # defining age column name is required
  age_name <- paste0(a,"_age")
  # value to search for
  phewas_ID_edit <- paste0(a,".")
  # combining values
  combinees <- PheWAS_manifest %>% 
    filter(str_detect(PheWAS_ID,phewas_ID_edit)) %>% 
    pull(PheWAS_ID)
  combining_variables <- data_field_variables[combinees] %>% 
    reduce(full_join,by="eid") 
  #getting combinations of column names to split columns by
  ID_colnames <- colnames(combining_variables)[-1]
  ID_columns_no_date <- ID_colnames[!grepl("earliest_date", ID_colnames)]
  ID_colnames_just_date <- ID_colnames[!grepl(a, ID_colnames)]
  ID_columns_age <- ID_colnames[grepl("age", ID_colnames)]
  non_data_columns <- c(ID_colnames_just_date,ID_columns_age)
  # just the data
  ID_colnames_just_data <- ID_colnames [! ID_colnames %in% non_data_columns]
  # need to recode other columns that are not dates
  combined_variables <-combining_variables %>% 
    mutate(across(all_of(ID_columns_no_date), as.double))
  # selecting the component parts to recombine later
  just_data <- combined_variables %>% 
    select(eid,all_of(ID_colnames_just_data))
  just_dates <- combined_variables %>% 
    select(eid,all_of(ID_colnames_just_date)) %>% 
    mutate(across(all_of(ID_colnames_just_date), as.character))
  just_age <- combined_variables %>% 
    select(eid,all_of(ID_columns_age))
  # creating columns to select max and min and record position
  selecting_data <- just_data %>% 
    mutate(results_list=transpose(across(all_of(ID_colnames_just_data))),
           max_result=unlist(map(results_list, ~ max(unlist(.x),na.rm = T))),
           min_result=unlist(map(results_list, ~ min(unlist(.x),na.rm = T))),
           position_min=unlist(map2(results_list,min_result, ~match(.y,unlist(.x)))),
           position_max=unlist(map2(results_list,max_result, ~match(.y,unlist(.x))))) %>% 
    select(eid,all_of(position_col)) %>% 
    rename(eid=1,pheno=2,position_match=3)
  # option to add in age column as required, otherwise combine with just_dates and select correct value         
  if(f=="age"){
    phenotype_extraction_final <- selecting_data %>% 
      left_join(just_dates, by="eid") %>% 
      mutate(date_list = transpose(across(all_of(ID_colnames_just_date))),
             earliest_date=map2(date_list,position_match, function(x,y) x[[y]])) %>% 
      select(eid,pheno,earliest_date) %>% 
      mutate(earliest_date=ymd(earliest_date)) %>% 
      set_names(c("eid",a,"earliest_date"))
    return(phenotype_extraction_final)
  } else if(f=="named"){
    phenotype_extraction_final <- selecting_data %>% 
      left_join(just_dates, by="eid") %>% 
      mutate(date_list = transpose(across(all_of(ID_colnames_just_date))),
             earliest_date=map2(date_list,position_match, function(x,y) x[[y]])) %>% 
      left_join(just_age, by="eid") %>% 
      mutate(age_list = transpose(across(all_of(ID_columns_age))),
             age=unlist(map2(age_list,position_match, function(x,y) x[[y]]))) %>% 
      select(eid,pheno,earliest_date,age) %>% 
      mutate(earliest_date=ymd(earliest_date)) %>% 
      set_names(c("eid",a,"earliest_date",age_name))
    return(phenotype_extraction_final) 
  }
}
# Load in defining variables -----------------------------------------------------
min_data <- fread(arguments$min_data)
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
# Select if wanting to append existing data_field phenotype file
if(!is.null(arguments$append_file)){
  data_field_variables <- readRDS(arguments$append_file)
} else {
  data_field_variables <- list()
}
# save location
if(arguments$phenotype_save_file=="data/phenotypes/data_field_phenotypes.RDS"){
  phenotype_save_location <- here("data","phenotypes","data_field_phenotypes.RDS")
  if(!dir.exists(here("data","phenotypes"))){
    dir.create(here("data","phenotypes"))
  }
} else {
  phenotype_save_location <- arguments$phenotype_save_file
  new_folder <- str_remove(arguments$phenotype_save_file,basename(arguments$phenotype_save_file))
  if(!dir.exists(new_folder)){
    dir.create(new_folder)}
}
# Data_field_phenotype function -------------------------------------------------------
# check to see if all the data is available and editing phenotype generation based on missing data and check to see if any existing phenotypes in data_field_variables
current_data_field <- gsub("-.*", "\\1", colnames(min_data))
current_IDs <- names(data_field_variables)
data_field <- PheWAS_manifest %>%  
  drop_na(field_code) %>% 
  filter(field_code %in% current_data_field,
         !PheWAS_ID %in% current_IDs) %>% 
  mutate(search_id=paste0("^",field_code,"-"),
         QC_flag_ID=ifelse(is.na(QC_flag_ID),QC_flag_ID,paste0("^",QC_flag_ID,"-")),
         date_code=ifelse(is.na(date_code),date_code,paste0("^",date_code,"-")),
         age_code=ifelse(is.na(age_code),age_code,paste0("^",age_code,"-"))) 
# function to extract data_field phenotypes
if (is.numeric(N_cores)) {
  data_field_variables_created <- mcmapply(data_field_extraction,
                                           data_field$PheWAS_ID,
                                           data_field$QC_flag_ID,
                                           data_field$search_id,
                                           data_field$analysis,
                                           data_field$case_code,
                                           data_field$control_code,
                                           data_field$limits,
                                           data_field$lower_limit,
                                           data_field$upper_limit,
                                           data_field$date_code,
                                           data_field$age_code,
                                           data_field$first_last_max_min_value,
                                           data_field$exclude_case_control_crossover,
                                           data_field$QC_filter_values,
                                           SIMPLIFY = F,
                                           mc.cores = N_cores,
                                           USE.NAMES = T)
} else {
  data_field_variables_created <- mapply(data_field_extraction,
                                         data_field$PheWAS_ID,
                                         data_field$QC_flag_ID,
                                         data_field$search_id,
                                         data_field$analysis,
                                         data_field$case_code,
                                         data_field$control_code,
                                         data_field$limits,
                                         data_field$lower_limit,
                                         data_field$upper_limit,
                                         data_field$date_code,
                                         data_field$age_code,
                                         data_field$first_last_max_min_value,
                                         data_field$exclude_case_control_crossover,
                                         data_field$QC_filter_values,
                                         SIMPLIFY = F,
                                         USE.NAMES = T)
}
# append to existing list either empty or with previously created phenotypes
data_field_variables <- append(data_field_variables,data_field_variables_created)
# Combined_data_field_phenotype function ---------------------------------
data_field_comb <- PheWAS_manifest %>% 
  filter(quant_combination==1) 
# function
if (length(data_field_comb >0)) {
  filed_ID_combined <- mapply(combining_data_field,
                              data_field_comb$PheWAS_ID,
                              data_field_comb$first_last_max_min_value,
                              data_field_comb$limits,
                              data_field_comb$lower_limit,
                              data_field_comb$upper_limit,
                              data_field_comb$age_col,
                              SIMPLIFY = F,
                              USE.NAMES = T)
  data_field_variables <- append(data_field_variables,filed_ID_combined)
}
# save file
saveRDS(data_field_variables,phenotype_save_location)