#! /usr/bin/env Rscript
'Using code lists provided to define concepts. Concepts can then be combined to form cases and controls for curated-phenotypes. A concept is a single homogenous group of codes that all describe a single disease, symptom, or related group of medicines. 

There are two outputs, one saves the summarised concept saved by default as concepts.RDS. THis contains a single record per individual that summarises the total number of occurances any of the codes making up that concept have occurrred in their records. This is then used in combination with other phenotypes/concepts to make composite phenotypes. The second file is all_dates.RDS (deafult), this contains records for every code within each participants records, each participant will have N records for N occurances of any of the codes of the concpt within their record. THis file is used to filter per event rather than the aggregate events which are summarised in the concept.RDS file.

UK Biobank prescription data is unusual in that it does not contain consistent codes (namely BNF codes) for most of the data. As such most of searching is done through using key search terms. A small subsection requires using Read_V2 drug codes. Because of this major difference with the other types of data, prescription concepts have their own function and would need to be used with caution in a non-UK Biobank setting (see user guide).

Users can choose to define only concepts that use clinical date or prescription data or both at least one must be inputted. The default save location is located within the Deep_PheWAS files the user can choose to use different folders. If specified folders will be created if not already present. The PheWAS_manifest which guides the phenotype creation can be chosen manually if it has been edited to add/change phenotypes, it is recommended to read the user guide if wanting to edit the PheWAS_manifest.

Usage:
    concept_creation.R --GPP=<FILE> [--health_data=<FILE> --concept_save_file=<FILE> --all_dates_save_file=<FILE>  --PheWAS_manifest_overide=<FILE> --code_list_folder=<FOLDER>]
    concept_creation.R --health_data=<FILE> [--GPP=<FILE> --concept_save_file=<FILE> --all_dates_save_file=<FILE>  --PheWAS_manifest_overide=<FILE> --code_list_folder=<FOLDER>]
    
Options:
    -h --help                                 Show this screen.
    -v --version                              Show version.
    
    Options
    --health_data=<FILE>                      Full file path of the tabdata file for UK-Biobank. Selecting will generate 
                                              clinical concepts.
    --GPP=<FILE>                              Full file path of the GP prescription data from UK-Biobank. Selecting will 
                                              generate prescription concepts.
                                              
    --concept_save_file=<FILE>                Full file path for the save file for the generated concepts RDS used for concept creation. 
                                              [default: data/phenotypes/concepts.RDS]
    --all_dates_save_file=<FILE>              Full file path for the save file for the generated all_dates.RDS used for per-event combinations of concepts. 
                                              [default: data/phenotypes/all_dates.RDS]
    --PheWAS_manifest_overide=<FILE>          Full file path of the alternative PheWAS_manifest file.
    --code_list_folder=<FOLDER>               Full file for the folder containing code lists.
                                              [default: data/concept_codes]
    
' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.1 concept_creation.R')

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(parallel))
suppressMessages(library(lubridate))
suppressMessages(library(stringi))
library(here)

# Functions ---------------------------------------------------------------
clinical_code_lookup <- function (x,z) {
  message(x)
  # this is the rated code list for each concept
  code_list_alpha <- fread(z)
  concept_name <- x
  sources_alpha <- unique(code_list_alpha$Source)
  # variable if statement to isolate presence or absence of ICD 10/9 or both
  if(any(str_detect("ICD10$",sources_alpha),na.rm = T) & any(str_detect("ICD9$",sources_alpha),na.rm = T)){
    # filter to remove ICD10 or 9 or both
    no_ICD <- code_list_alpha %>% 
      filter(Source!="ICD10", Source!="ICD9")
    # create a seperate data frame for each position of ICD code in 10/9 or both
    ICD10_1 <- code_list_alpha %>% 
      filter(Source=="ICD10") %>% 
      mutate(Source="ICD10_1")
    ICD10_2 <- code_list_alpha %>% 
      filter(Source=="ICD10") %>% 
      mutate(Source="ICD10_2")
    ICD10_3 <- code_list_alpha %>% 
      filter(Source=="ICD10") %>% 
      mutate(Source="ICD10_3")
    ICD9_1 <- code_list_alpha %>% 
      filter(Source=="ICD9") %>% 
      mutate(Source="ICD9_1")
    ICD9_2 <- code_list_alpha %>% 
      filter(Source=="ICD9") %>% 
      mutate(Source="ICD9_2")
    ICD9_3 <- code_list_alpha %>% 
      filter(Source=="ICD9") %>% 
      mutate(Source="ICD9_3")
    # combine into a new code_list variable that is used throughout function
    code_list <- no_ICD %>% 
      bind_rows(ICD10_1,ICD10_2,ICD10_3,ICD9_1,ICD9_2,ICD9_3)
    sources <- unique(code_list$Source)
  } else if(any(str_detect("ICD10$",sources_alpha),na.rm = T)) {
    no_ICD <- code_list_alpha %>% 
      filter(Source!="ICD10")
    ICD10_1 <- code_list_alpha %>% 
      filter(Source=="ICD10") %>% 
      mutate(Source="ICD10_1")
    ICD10_2 <- code_list_alpha %>% 
      filter(Source=="ICD10") %>% 
      mutate(Source="ICD10_2")
    ICD10_3 <- code_list_alpha %>% 
      filter(Source=="ICD10") %>% 
      mutate(Source="ICD10_3")
    code_list <- no_ICD %>% 
      bind_rows(ICD10_1,ICD10_2,ICD10_3)
    sources <- unique(code_list$Source)
  } else if(any(str_detect("ICD9$",sources_alpha),na.rm = T)){
    no_ICD <- code_list_alpha %>% 
      filter(Source!="ICD9")
    ICD9_1 <- code_list_alpha %>% 
      filter(Source=="ICD9") %>% 
      mutate(Source="ICD9_1")
    ICD9_2 <- code_list_alpha %>% 
      filter(Source=="ICD9") %>% 
      mutate(Source="ICD9_2")
    ICD9_3 <- code_list_alpha %>% 
      filter(Source=="ICD9") %>% 
      mutate(Source="ICD9_3")
    code_list <- no_ICD %>% 
      bind_rows(ICD9_1,ICD9_2,ICD9_3)
    sources <- unique(code_list$Source)
  } else {
    code_list <- code_list_alpha
    sources <-  sources_alpha
  }
  
  concept_source_data_extraction <- function(a) {
    
    all_codes <- code_list %>% 
      filter(Source==a) %>% 
      mutate(Code=as.character(Code))
    
    included_codes  <- all_codes %>% 
      filter(Source==a,Decision==1)
    
    heath_data_source <- health_data %>% 
      filter(source==a)
    
    # and looks those up in the database 
    all_codes_extracted <- heath_data_source %>% 
      filter(code %in% all_codes$Code)  
    
    included_codes_extracted <- heath_data_source %>% 
      filter(code %in% included_codes$Code) 
    
    ### wide date
    min_date_label <- paste0("first_code_",a)
    
    wide_date <- included_codes_extracted %>% 
      group_by(eid) %>% 
      summarise(!!min_date_label:=min(date))
    
    #### wide file
    N_source_type <- included_codes_extracted %>% 
      count(eid) %>% 
      rename(!!a:=n)
    
    #### all dates
    all_dates <- included_codes_extracted %>% 
      select(eid,date)
    
    #### TN
    included_ID <- unique(included_codes_extracted$eid)
    included_ID_N <-length(included_ID)
    type_of_code <- a
    half_TN <- data.frame(included_ID_N,type_of_code) %>% 
      select(Source=type_of_code,Total_participants=included_ID_N)
    
    #### CTID
    N_codes <- all_codes_extracted %>% 
      group_by(code) %>% 
      count() %>% 
      mutate(Source={{a}},) %>% 
      rename(N_Codes=n)
    
    N_codes_ID <- all_codes_extracted %>% 
      group_by(code) %>% 
      summarise(N_Participants= n_distinct(eid))
    
    # and identifies how many participants have only one code out of the code list in their records  
    IDs_only_one_code <- all_codes_extracted[all_codes_extracted$eid %in% names(which(table(all_codes_extracted$eid) < 2)), ]
    only_one_code <- IDs_only_one_code %>% 
      group_by(code) %>% 
      summarise(N_Participants_with_1_code= n_distinct(eid))
    
    CTID_primer <- all_codes %>% 
      left_join(N_codes,by=c("Code"="code","Source")) %>% 
      left_join(N_codes_ID,by=c("Code"="code")) %>% 
      left_join(only_one_code,by=c("Code"="code")) %>% 
      select(Code,Description,Source,Decision,N_Codes,N_Participants,N_Participants_with_1_code)
    
    output <- list(wide_date=wide_date,wide_file=N_source_type,all_dates=all_dates,half_TN=half_TN,included_IDs=included_ID,CTID_primer=CTID_primer)
    
    return(output)
    
  }
  
  per_source_data <- lapply(sources,concept_source_data_extraction)
  
  ## Wide dates all file primary file used for phenotype generation
  wide_dates <- lapply(per_source_data,'[[',"wide_date") %>% 
    reduce(full_join) %>% 
    mutate(earliest_date = pmin(!!! rlang::syms(names(.)[2:ncol(.)]), na.rm = T)) %>% 
    select(eid,earliest_date)
  
  wide_count <- lapply(per_source_data,'[[',"wide_file") %>% 
    reduce(full_join) %>% 
    rowwise() %>% 
    mutate(any_code=sum(!!! rlang::syms(names(.)[2:ncol(.)]), na.rm = T))
  
  WA <- wide_dates %>% 
    left_join(wide_count)
  
  WA[is.na(WA)] <- 0
  
  WDA <- WA %>% 
    select(eid,!!x:=any_code,earliest_date)
  
  ##  all dates
  AD <- lapply(per_source_data,'[[',"all_dates") %>% 
    reduce(bind_rows)
  
  ## AID
  AID <- data.frame(length(unique(lapply(per_source_data,'[[',"included_IDs") %>% 
                                    reduce(c)))) %>% 
    rename(N_participants=1)
  
  ##  TN
  half_TN <- lapply(per_source_data,'[[',"half_TN") %>% 
    reduce(bind_rows)
  
  included_ID <- lapply(per_source_data,'[[',"included_IDs")
  
  calc_minus_source_ID <- function(x) {
    current_ID <- unlist(included_ID[x])
    alternate_ID <- unlist(included_ID[-x])
    exclusive_ID <- sum((!current_ID %in% alternate_ID),na.rm = T)
    
    return(exclusive_ID)
  }
  range_of_x <- c(1:length(included_ID))
  only_source_ID <- lapply(range_of_x,calc_minus_source_ID) %>% 
    reduce(c) %>% 
    as.data.frame() %>% 
    rename(N_Participants_only_this_code_source=1)
  
  TN <- half_TN %>% 
    bind_cols(only_source_ID)
  
  ## CTID
  CTID_primer <- lapply(per_source_data,'[[',"CTID_primer") %>% 
    reduce(bind_rows) 
  CTID_primer[is.na(CTID_primer)] <- 0
  
  CTID <- CTID_primer %>% 
    arrange(desc(Decision),desc(N_Participants),Source)
  
  #  fwrite(WA,paste0(concept_save_location,"/",paste0(concept_name,"_WA.txt")), sep = "\t")
  #  fwrite(WDA,paste0(concept_save_location,"/",paste0(concept_name,"_WDA.txt")),quote = T, na = NA, sep = "\t")
  #  fwrite(AID,paste0(concept_save_location,"/",paste0(concept_name,"_AID.txt")), sep = "\t")
  #  fwrite(AD,paste0(concept_save_location,"/",paste0(concept_name,"_AD.txt")), sep = "\t")
  #  fwrite(CTID,paste0(concept_save_location,"/",paste0(concept_name,"_CTID.txt")),quote = T, na = NA, sep = "\t")
  
  new_return <- list(WDA=WDA,AD=AD)  
  return(new_return)
}
drug_code_lookup <- function(x,y,z) {
  message(x)
  ## creating a name for the files
  both_names <- x
  
  drug_list_BNF_DMD <- fread(z,sep = ",") %>% 
    pull()
  
  drug_list_V2 <- fread(y) %>% 
    filter(Decision==1) %>% 
    select(1) %>% 
    pull()
  
  ## variable used for regex search function used below
  string_regex_options=stri_opts_regex(case_insensitive = T)
  
  ## use these names in a string search, this not only searches for the keywords but also reports which keywords are then found in each row. I use this
  ## to add a keyword for describing the number of participants found for each broad keyword.
  drug_text_search <- GP_P %>% 
    filter(str_detect(drug_name, regex((paste(drug_list_BNF_DMD, collapse="|")),ignore_case = TRUE))) %>% 
    mutate(keyword_or_V2code=(stri_extract_first_regex(drug_name,paste(drug_list_BNF_DMD, collapse = '|'),simplify = T,opts_regex = string_regex_options)),
           keyword_or_V2code=as.character(keyword_or_V2code),
           keyword_or_V2code=str_to_title(keyword_or_V2code),
           Source="BNF/DMD")
  
  ## now use read_v2 codes to find the small percentage of codes with only read_v2 identification this only searches in those with read_v2 codes and not those eids
  ## already found within the keyword search above
  drug_V2_search <- GP_P %>% 
    filter(read_2!="NA", !row_number %in% drug_text_search$row_number, read_2 %in% drug_list_V2) %>% 
    mutate(keyword_or_V2code=read_2,
           Source="Read_V2") %>% 
    left_join(V2_drugs, by=c("read_2"="read_code")) %>% 
    select(-status_flag)
  
  ## combine for a list of results
  all_drug_results <- bind_rows(drug_text_search,drug_V2_search) %>% 
    distinct(row_number, .keep_all = TRUE)
  
  ## create all_id
  all_id_drugs <- list(length(unique(all_drug_results$eid)))
  
  ## This creates the codes total ID's or CTID files for prescription data, it does this for the keyword/Read_v2 code
  ## first it counts all incidents of the codes
  CTID_ish <- all_drug_results %>% 
    group_by(keyword_or_V2code) %>% 
    count() %>% 
    rename(N_Codes=n)
  
  ## then counts the N of unique eids i.e. how many particpants in total this represents
  CTID2_ish <- all_drug_results %>% 
    group_by(keyword_or_V2code) %>%
    summarise(N_Participants=n_distinct(eid))
  
  ## identifies those with only one code in the record
  prescription_only_one_code <- all_drug_results %>% 
    group_by(eid) %>% 
    summarise(N_eid=n()) %>% 
    filter(N_eid <2)
  
  ## and uses the info above to show how many for each keyword/Read_v2 code have only a single code
  prescription_ids_ooc <- all_drug_results %>% 
    filter(eid %in% prescription_only_one_code$eid) %>% 
    group_by(keyword_or_V2code) %>% 
    summarise(N_Particicpants_with_1_code= n_distinct(eid))
  
  ## calculates how many prescriptions for each code there is on average
  CTID3_ish <- all_drug_results %>% 
    group_by(keyword_or_V2code) %>%
    summarise(Mean_prescriptions_per_participant=round((n()/n_distinct(eid)),digits = 1))
  
  ## this does the same as the above but removing those with only one prescription in their records
  CTID4_ish <- all_drug_results %>% 
    filter(!eid %in% prescription_only_one_code$eid) %>% 
    group_by(keyword_or_V2code) %>%
    summarise(`Mean_prescriptions_per_participant*`=round((n()/n_distinct(eid)),digits=1))
  
  ## the two below add in data on description and Source for the final table
  read_descriptions <- drug_V2_search %>% 
    group_by(keyword_or_V2code) %>% 
    summarise(terms_x=unique(term_description),Source=unique(Source))
  
  drug_list_description <- fread(z) %>% 
    mutate(Source="BNF/DMD",
           Description=drug_unique) %>% 
    rename(keyword_or_V2code=1)
  
  ## making the final CTID table
  CTID_drugs <- CTID_ish %>% 
    left_join(CTID2_ish)  %>% 
    left_join(prescription_ids_ooc) %>% 
    left_join(CTID3_ish) %>% 
    left_join(CTID4_ish)  %>%
    left_join(read_descriptions) %>% 
    left_join(drug_list_description, by="keyword_or_V2code") %>% 
    mutate(Source=ifelse(is.na(Source.x),Source.y,Source.x),
           Description=ifelse(is.na(Description),terms_x,Description),
           Decision=1) %>% 
    select(1,Description,Source,Decision,N_Codes,N_Participants,N_Particicpants_with_1_code,Mean_prescriptions_per_participant,`Mean_prescriptions_per_participant*`) %>% 
    mutate_all(~replace(., is.na(.), 0))
  
  ## now formatting dates
  all_drugs_with_date <- all_drug_results %>% 
    filter(issue_date!="NA") %>% 
    rename(code_date =issue_date) %>% 
    mutate(code_date=ymd(code_date))
  all_drugs_with_date <- all_drugs_with_date %>% 
    drop_na(code_date) 
  
  ## create earliest and latest date of prescription
  all_drugs_time <- all_drugs_with_date %>% 
    group_by(eid) %>% 
    summarize(earliest_date=min(code_date), any_code=max(code_date), any_code=n()) 
  
  ## all dates file
  all_dates_drugs <- all_drugs_with_date %>% 
    select(eid,code_date)
  
  ## to make equivalent of all_wide_dates and all_wide file need to create a separation in the codes into BNF/DMD and Read_V2
  ## do for BNF
  BNF_DMD <- all_drug_results %>% 
    filter(Source=="BNF/DMD")
  BNF_DMD_wide <- BNF_DMD %>% 
    count(eid,keyword_or_V2code) %>% 
    group_by(keyword_or_V2code) %>%
    mutate(row = row_number()) %>% 
    pivot_wider(names_from = keyword_or_V2code, values_from = n) %>% 
    select(-row) %>% 
    group_by(eid) 
  
  BNF_DMD_wide[is.na(BNF_DMD_wide)] <- 0
  
  BNF_DMD_wide <- BNF_DMD_wide %>%
    summarise_all(list(max)) %>% 
    mutate(BNF_DMD = rowSums(.[-1],na.rm = TRUE))
  
  BNF_DMD_wide_date <- BNF_DMD %>% 
    group_by(eid) %>% 
    summarise(BNF_DMD_first_code=min(issue_date))
  
  ## Now for Read_V2
  Read_V2 <- all_drug_results %>% 
    filter(Source=="Read_V2") 
  
  Read_V2_wide <- Read_V2 %>% 
    count(eid,keyword_or_V2code) %>% 
    group_by(keyword_or_V2code) %>%
    mutate(row = row_number()) %>% 
    pivot_wider(names_from = keyword_or_V2code, values_from = n) %>% 
    select(-row) %>% 
    group_by(eid) 
  
  Read_V2_wide[is.na(Read_V2_wide)] <- 0
  
  Read_V2_wide <- Read_V2_wide %>%
    summarise_all(list(max)) %>% 
    mutate(Read_V2 = rowSums(.[-1],na.rm = TRUE))
  
  Read_V2_wide_date <- Read_V2 %>% 
    group_by(eid) %>% 
    summarise(Read_V2_first_code=min(issue_date))
  
  ## Now combine them
  wide_drugs <- BNF_DMD_wide %>% 
    left_join(Read_V2_wide) %>% 
    select(eid,BNF_DMD,Read_V2) 
  
  wide_drugs[is.na(wide_drugs)] <- 0
  
  wide_drugs <- wide_drugs %>% 
    mutate(any_code = BNF_DMD + Read_V2)
  
  ##now with dates added in
  wide_drugs_date <- wide_drugs %>% 
    left_join(BNF_DMD_wide_date) %>% 
    left_join(Read_V2_wide_date) %>% 
    mutate(earliest_date=ifelse(BNF_DMD_first_code==0, Read_V2_first_code, BNF_DMD_first_code)) %>% 
    select(eid,earliest_date,BNF_DMD,Read_V2, any_code) %>% 
    filter(earliest_date!="NA")
  
  wide_dates_all <- wide_drugs_date %>% 
    select(eid,!!x:=any_code,earliest_date)
  
  #Total Numbers
  BNF_DMD_ID <- length(unique(BNF_DMD$eid))
  Read_V2_ID <- length(unique(Read_V2$eid))
  
  ID_N <- c(BNF_DMD_ID,Read_V2_ID)
  just_ID_N <- c(BNF_DMD_ID,Read_V2_ID)
  numbers_row <- c("BNF/DMD","Read_V2")
  
  total_numbers_drugs <- setNames(data.frame(matrix(ncol = 1, nrow = 2)), "flib") %>% 
    add_column(ID_N,just_ID_N,numbers_row) %>% 
    select(Source=numbers_row,Total_Participants=ID_N,N_Participants_only_this_code_source=just_ID_N)
  
  #  fwrite(wide_drugs_date,paste0(concept_save_location,"/",paste0(both_names,"_WA.txt")), quote = TRUE, na=NA)
  #  fwrite(all_dates_drugs,paste0(concept_save_location,"/",paste0(both_names,"_AD.txt")), quote = TRUE, na=NA)
  #  fwrite(CTID_drugs,paste0(concept_save_location,"/",paste0(both_names,"_CTID.txt")), quote = TRUE, na=NA)
  #  fwrite(all_id_drugs,paste0(concept_save_location,"/",paste0(both_names,"_AID.txt")), quote = TRUE, na=NA)
  #  fwrite(total_numbers_drugs,paste0(concept_save_location,"/",paste0(both_names,"_TN.txt")), quote = TRUE, na=NA)
  #  fwrite(wide_dates_all,paste0(concept_save_location,"/",paste0(both_names,"_WDA.txt")), quote = TRUE, na=NA)
  
  new_return <- list(WDA=wide_dates_all,AD=all_dates_drugs)
  
  return(new_return)
}
# Load in defining variables -----------------------------------------------------
# creating empty lists to fill with function results
concepts <- list()
all_dates <- list()
# defining save locations for
# concept.RDS
if(arguments$concept_save_file=="data/phenotypes/concepts.RDS"){
  phenotype_save_location <- here("data","phenotypes","concepts.RDS")
  if(!dir.exists(here("data","phenotypes"))){
    dir.create(here("data","phenotypes"))
  }
} else {
  phenotype_save_location <- arguments$concept_save_file
  new_folder <- str_remove(arguments$concept_save_file,basename(arguments$concept_save_file))
  if(!dir.exists(new_folder)){
    dir.create(new_folder)}
}
# all dates
if(arguments$all_dates_save_file=="data/phenotypes/all_dates.RDS"){
  all_dates_save_location <- here("data","phenotypes","all_dates.RDS")
  if(!dir.exists(here("data","phenotypes"))){
    dir.create(here("data","phenotypes"))
  }
} else {
  all_dates_save_location <- arguments$all_dates_save_file
  new_folder <- str_remove(arguments$all_dates_save_file,basename(arguments$all_dates_save_file))
  if(!dir.exists(new_folder)){
    dir.create(new_folder)}
}
# PheWAS_manifest
if(is.null(arguments$PheWAS_manifest_overide)){
  PheWAS_manifest <- fread(here("data","PheWAS_manifest.csv"))
} else {
  PheWAS_manifest <- fread(arguments$PheWAS_manifest_overide)
}
# code list folder location
if(arguments$code_list_folder=="data/concept_codes"){
  code_list_folder <- here("data","concept_codes")
} else {
  code_list_folder <- arguments$code_list_folder
}
# Clinical_code_concepts --------------------------------------------------
if(!is.null(arguments$health_data)) {
  # health_data
  health_data <- fread(arguments$health_data)
  # Clinical search codes works by searching for all clinical code list files in specified folder
  concept_clinical_codes <- data.frame(names=(str_remove(list.files(path = code_list_folder, 
                                                                    pattern = "codes_rated.csv$"),"_codes_rated.csv")))%>% 
    left_join(PheWAS_manifest, by=c("names"="concept_name")) %>% 
    drop_na(PheWAS_ID) %>% 
    select(names) %>% 
    pull()
  # gets PheWAS_ID which is used to save and then access the files later from joining with the PheWAS_manifest
  concept_clinical_codes_PheWAS_ID <- data.frame(names=(str_remove(list.files(path = code_list_folder,
                                                                              pattern = "codes_rated.csv$"),"_codes_rated.csv"))) %>% 
    left_join(PheWAS_manifest, by=c("names"="concept_name")) %>%  
    drop_na(PheWAS_ID) %>% 
    select(PheWAS_ID) %>% 
    pull()
  message(concept_clinical_codes)
  # running the clinical concept function
  concept_clinical_codes_list <- lapply(concept_clinical_codes, function(x) paste0(code_list_folder,"/",x,"_codes_rated.csv"))
  clinical_concepts <- mapply(clinical_code_lookup,
                              concept_clinical_codes_PheWAS_ID,
                              concept_clinical_codes_list,
                              SIMPLIFY = F,
                              USE.NAMES = T)
  # is a list with two outputs, retrive WDA which saves as concepts.RDS and AD which saves as all_dates.RDS
  WDA <- lapply(clinical_concepts,'[[',"WDA")
  AD <- lapply(clinical_concepts,'[[',"AD")
  new_names <- paste0(base::names(AD),"_AD")
  base::names(AD) <- new_names
  # appending lists
  concepts <- append(concepts,WDA)
  all_dates <-append(all_dates,AD)
  # for running individual code list lookups
  #l <- concept_clinical_codes[5]
  #m <- concept_clinical_codes_list[5]
  #mapply(clinical_code_lookup,l,m)
}
# Prescription concepts ---------------------------------------------------
if(!is.null(arguments$GPP)) {
  # prescription data
  GP_P <- fread(arguments$GPP, na.strings = "") %>% 
    select(eid,read_2,issue_date,drug_name) %>% 
    mutate(row_number = 1:n())
  V2_drugs <- fread(here("data","V2_drugs.csv"))
  # Prescription search codes works by searching for all clinical code list files in specified folder
  prescription_search_terms <- data.frame(names=(str_remove(list.files(path = code_list_folder, pattern = "DMD.csv$"),"_BNF_DMD.csv"))) %>% 
    left_join(PheWAS_manifest, by=c("names"="concept_name")) %>% 
    drop_na(PheWAS_ID) %>% 
    select(names) %>% 
    pull()
  # and extracting the PheWAS ID by joining with the PheWAS_mainfest
  prescription_search_terms_PheWAS_ID <- data.frame(names=(str_remove(list.files(path = code_list_folder, pattern = "DMD.csv$"),"_BNF_DMD.csv"))) %>% 
    left_join(PheWAS_manifest, by=c("names"="concept_name")) %>% 
    drop_na(PheWAS_ID) %>% 
    select(PheWAS_ID) %>% 
    pull()
  # list of file locations
  prescription_search_terms_list <- lapply(prescription_search_terms, function(x) paste0(code_list_folder,"/",x,"_BNF_DMD.csv"))
  # V2 files here unique quirk of the UK Biobank data
  V2_search_terms <- as.list(list.files(path = code_list_folder,pattern = "V2_rated.csv$"))
  V2_search_terms_list <- lapply(V2_search_terms, function(x) paste0(code_list_folder,"/",x))
  # run the function
  drug_concepts <- mapply(drug_code_lookup,
                          prescription_search_terms_PheWAS_ID,
                          V2_search_terms_list,
                          prescription_search_terms_list,
                          SIMPLIFY = F,
                          USE.NAMES = T)
  # has two outputs as list need to extract and then append to lists
  WDA <- lapply(drug_concepts,'[[',"WDA")
  AD <- lapply(drug_concepts,'[[',"AD")
  new_names <- paste0(base::names(AD),"_AD")
  base::names(AD) <- new_names
  # append to lists
  concepts <- append(concepts,WDA)
  all_dates <-append(all_dates,AD)
}
# Save RDS ----------------------------------------------------------------
saveRDS(concepts,phenotype_save_location)
saveRDS(all_dates,all_dates_save_location) 