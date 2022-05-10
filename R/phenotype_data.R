#! /usr/bin/env Rscript
'Combines health data from UK-Biobank into a single long file with eid, code, date, source. Also creates seven additonal files used in various stages in the pipeline. 

Current sources of data for the health records are primary care clinical codes Read V2 and V3 (V2,V3), ICD9 and 10 codes (ICD9,ICD10), UK biobank self-report non-cancer illness codes (SR), UK biobank self-report operation codes (SROP), OPCS-4 codes (OPCS). Primary care prescription data from UKB does not lend itself to be combined in this way and so is edited and saved as a separate file alondside an edited version of the GP clinical data, and list of IDs of those with GP data.

Alongside self-reported non-cancer codes and the cancer registry data the --min_data flag will also generate combined_sex and call_rate_kinship, both files are used further in the pipeline and are explained in more detail in the user guide. Finally the control_populations file is created this is used to define two population for 

All inputs are optional, the script will combine whichever of the allowable files are inputted, to generate the full list of phenotypes availble in the phenotype_manifest.csv all inputs must be provided, the script will produce no output if none of the options are selected. In cases where a data field-ID is availble but the corresponding date field-ID is unavailable neither will be extracted. All ouputs are saved in data within the Deep_PheWAS file structure.

Usage:
    phenotype_data.R [--min_data=<FILE> --GPC=<FILE> --GPP=<FILE> --hesin_diag=<FILE> --HESIN=<FILE> --hesin_oper=<FILE> --death_cause=<FILE> --death=<FILE> --exclusions=<FILE> --king_coef=<FILE> --save_location=<FOLDER>]
    
Options:
    -h --help                 Show this screen.
    -v --version              Show version.
    
    Options
    --min_data=<FILE>         Full file path of the minimum_tab_data file for UK-Biobank, see minimum_tab_data_edit.R for 
                              instructions and user guide for details of the file format.
    --GPC=<FILE>              Full file path of the GP clinical data file from UK-Biobank.
    --GPP=<FILE>              Full file path of the GP prescription data from UK-Biobank.
    --hesin_diag=<FILE>       Full file path of hesin_diag file from UK-Biobank.
    --HESIN=<FILE>            Full file path of HESIN file from UK-Biobank.
    --hesin_oper=<FILE>       Full file path of the hesin_oper file from UK-Biobank.
    --death_cause=<FILE>      Full file path of the death_cause file from UK-Biobank.
    --death=<FILE>            Full file path of the death file from UK-Biobank contains date of death info.
    --exclusions=<FILE>       Full file path of the exclusions file that documents the latest exclusion list for 
                              participants from UK-Biobank.
    --king_coef=<FILE>        Full file path of related data file containing King coefficient scores for related ID pairs 
                              from UK-Biobank.
    --save_location=<FOLDER>  Full file path for the common folder to save created files deafults to data/. 
                              [default: data/]
' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.1 phenotype_data.R')

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(lubridate))
library(here)

# Functions ---------------------------------------------------------------
SR_data_out <-  function(a,b,c,d,e) {
  if (b=="NC") {
    field_id_date <- paste0("^20008-",a)
    field_id_code <- paste0("^20002-",a)
  } else if (b=="OP") {
    field_id_date <- paste0("^20010-",a)
    field_id_code <- paste0("^20004-",a)
  }
  
  field_id_assesment <- paste0("^53-",a) 
  field_id_assesment_selection <- paste0("53-",a,".0")
  
  SR_dates <- get(c) %>% 
    select(1,matches(field_id_date)) %>% 
    pivot_longer(c(2:length(colnames(.))),names_to = "type",values_to = "year_dx") %>% 
    rowid_to_column() %>% 
    drop_na() %>% 
    select(rowid,year_dx)
  
  SR_visit <- get(d) %>% 
    select(1,matches(field_id_code)) %>% 
    pivot_longer(c(2:length(colnames(.))),names_to = "type",values_to = "code") %>% 
    rowid_to_column() %>% 
    drop_na() %>% 
    left_join(select(assesment_centre_date,1,matches(field_id_assesment))) %>% 
    left_join(SR_dates,by=c("rowid")) %>% 
    select(eid,code,date_of_dx=year_dx,date_of_visit={{field_id_assesment_selection}}) %>% 
    mutate(date_of_dx = ifelse(date_of_dx == -1, NA, ifelse(date_of_dx == -3, NA, date_of_dx)),
           date_of_dx = ymd(round_date(date_decimal(date_of_dx), unit = "day")),
           date_of_visit = ymd(date_of_visit),
           date = coalesce(date_of_dx,date_of_visit),
           source=e,
           code=as.character(code)) %>% 
    select(eid,code,date,source)
  
  return(SR_visit)
  
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
# min_data
if(!is.null(arguments$min_data)){
  tab_data <- fread(arguments$min_data) %>% 
    filter(!eid %in% exclusions)
  available_tab_data <- colnames(tab_data)
}
# empty df
all_phenotype_data <- data.frame(remove=NA)
# the assessment centre date
if(length(str_which(available_tab_data,"^53-"))>0) {
  assesment_centre_date <- tab_data %>% 
    select(1,matches("^53-"))
}
# save location
if(arguments$save_location == "data/") {
  save_location <- paste0(here("data"),"/")
} else {
  save_location <- arguments$save_location
}
# UK Biobank self-reported non-cancer diagnosis  --------------------------------------
# self-report non-cancer codes
if(!is.null(arguments$min_data)) {
  if(length(str_which(available_tab_data,"^20002-"))>0 & length(str_which(available_tab_data,"^20008-"))>0) {
    SR_NC <- tab_data %>% 
      select(1,matches("^20002-")) 
    # the year of diagnosis
    SR_NC_year_dx <- tab_data %>% 
      select(1,matches("^20008-")) 
    # calculating N of visits
    SR_NC_visits <- as.numeric(unique(gsub(".*-(.+)\\..*","\\1",colnames(SR_NC_year_dx)[-1])))
    # function to process the data
    SR_NC_data <- mapply(SR_data_out,SR_NC_visits,MoreArgs=list(b="NC",c="SR_NC_year_dx",d="SR_NC",e="SR"), SIMPLIFY = F) %>% 
      reduce(bind_rows)
    # combine
    all_phenotype_data <- all_phenotype_data  %>% 
      bind_rows(SR_NC_data)
  }
  # UK Biobank self-report operations  -----------------------------------------------------------
  if(length(str_which(available_tab_data,"^20004-"))>0 & length(str_which(available_tab_data,"^20010-"))>0) {
    SR_OP <- tab_data %>% 
      select(1,matches("^20004-")) 
    # the date of operation
    SR_OP_year_dx <- tab_data %>% 
      select(1,matches("^20010-")) 
    # calculating N of visits
    SR_OP_visits <- as.numeric(unique(gsub(".*-(.+)\\..*","\\1",colnames(SR_OP_year_dx)[-1])))
    # function to process the data
    SR_OP_data <- mapply(SR_data_out,SR_OP_visits,MoreArgs=list(b="OP",c="SR_OP_year_dx",d="SR_OP",e="SROP"), SIMPLIFY = F) %>% 
      reduce(bind_rows)
    # combine
    all_phenotype_data <- all_phenotype_data  %>% 
      bind_rows(SR_OP_data)
  }
  # Cancer_registry ---------------------------------------------------------
  # split into date and diagnostic codes
  if(length(str_which(available_tab_data,"^40005-"))>0 & length(str_which(available_tab_data,"^40006-"))>0) {
    # date
    cancer_reg_date <- tab_data %>% 
      select(1,matches("^40005-")) %>% 
      pivot_longer(c(2:length(colnames(.))),names_to = "type",values_to = "year_dx") %>% 
      rowid_to_column() %>% 
      drop_na() %>% 
      select(rowid,year_dx)
    # data
    cancer_reg_data <- tab_data %>% 
      select(1,matches("^40006-")) %>% 
      mutate_all(na_if,"") %>% 
      pivot_longer(c(2:length(colnames(.))),names_to = "type",values_to = "code") %>% 
      rowid_to_column() %>% 
      drop_na() %>% 
      left_join(cancer_reg_date, by=c("rowid")) %>% 
      drop_na() %>% 
      mutate(code=ifelse(nchar(code) > 4, (str_extract(code, "^.{4}")), code),
             date=ymd(year_dx),
             source="cancer") %>% 
      select(eid,code,date,source)
    # combine
    all_phenotype_data <- all_phenotype_data  %>% 
      bind_rows(cancer_reg_data)
  }
  # combined_sex ---------------------------------------------------------------------
  ## sex 1 = male, 0 = female
  if(length(str_which(available_tab_data,"^22001-"))>0 & length(str_which(available_tab_data,"^31-"))>0) {
    # genetic sex
    gen_sex <- tab_data %>% 
      select(1,matches("^22001-")) %>% 
      rename(eid=1,gen_sex=2)
    # self-report sex
    reported_sex <- tab_data %>% 
      select(1,matches(("^31-"))) %>% 
      rename(eid=1,reported_sex=2)
    # select genetic sex unless no genetic sex reported
    combined_sex <- gen_sex %>% 
      full_join(reported_sex) %>% 
      mutate(sex=ifelse(is.na(gen_sex),reported_sex,gen_sex)) %>% 
      drop_na(sex) %>% 
      filter(!eid %in% exclusions) %>% 
      select(eid,sex)
    # write
    fwrite(combined_sex,paste0(save_location,"/combined_sex"),sep = "\t")
  }
# Loss to followup --------------------------------------------------------
  if(length(str_which(available_tab_data,"^190-0.0"))>0) {
    loss_to_follow_up <- tab_data %>% 
      select(1,`190-0.0`) %>% 
      filter(`190-0.0`!=1) %>% 
      select(eid)
    #write
    fwrite(loss_to_follow_up, paste0(save_location,"/control_exclusions"), na = NA, sep = "\t")
  }
  # Kinship_call_rate -------------------------------------------------------
  if(length(str_which(available_tab_data,"^22005-"))>0 & !is.null(arguments$king_coef)) {
    # KINGS coefficent
    related <- fread(arguments$king_coef) %>% 
      select(ID1,ID2,Kinship)
    # call rate
    call_r <- tab_data %>% 
      select(1,matches("^22005-")) %>% 
      rename(missingness=`22005-0.0`)
    # combine to give missingness information to each of ID, this is done as when removing related individuals you remove the one with the highest missingness
    call_rate_kinship <- related %>% 
      left_join(call_r,by=c("ID1"="eid")) %>%
      rename(missingness_ID1=missingness) %>% 
      left_join(call_r,by=c("ID2"="eid")) %>% 
      rename(missingness_ID2=missingness) %>% 
      filter(Kinship>=0.0884) %>% 
      mutate(lower_missing = if_else(missingness_ID1 < missingness_ID2, 1, 2)) %>% 
      filter(!ID1 %in% exclusions,!ID2 %in% exclusions)
    # write
    fwrite(call_rate_kinship,here("data","related_callrate"), sep = "\t")
  }
}
# GP data -----------------------------------------------------------------
# simple edit to remove exclusions and format date
if(!is.null(arguments$GPC)) {
  # GP clinical data
  GP_C <- fread(arguments$GPC, na.strings = "") %>% 
    filter(!eid %in% exclusions) %>% 
    mutate(event_dt=dmy(event_dt))
  # need to retain an edited copy as further phenotypes require this information
  fwrite(GP_C,paste0(save_location,"/GP_C_edit.txt.gz"), na = NA, quote = TRUE)
  # GP_ID file
  GP_ID <- GP_C %>% 
    distinct(eid)
  fwrite(GP_ID,paste0(save_location,"/GP_C_ID.txt.gz"), na = NA, quote = TRUE)
  # split V2 and V3 as source
  GP_read_V2 <- GP_C %>%
    drop_na(read_2) %>% 
    mutate(source="V2") %>% 
    select(eid,code=read_2,date=event_dt,source) %>% 
    drop_na()
  GP_read_V3 <- GP_C %>%
    drop_na(read_3) %>% 
    mutate(source="V3") %>% 
    select(eid,code=read_3,date=event_dt,source) %>%
    drop_na()
  # combine
  all_phenotype_data <- all_phenotype_data  %>% 
    bind_rows(GP_read_V3,GP_read_V2)
}
## Need to edit GP_P data to make DMD codes a character for searching prescription data is used separately in the later functions and so is saved as an edited file
if(!is.null(arguments$GPP)) {
  # GP prescription data
  GP_P <- fread(arguments$GPP, na.strings = "") %>% 
    filter(!eid %in% exclusions) %>% 
    mutate(issue_date=dmy(issue_date)) %>% 
    mutate(issue_date=format(issue_date, "%Y-%m-%d"),
           dmd_code=as.character(dmd_code))
  # write
  fwrite(GP_P,paste0(save_location,"/GP_P_edit.txt.gz"), na = NA, quote = TRUE)
}
# HES ---------------------------------------------------------------------
# HES data needs to be combined from 2 tables, HES_diag and HESIN, HESIN has dates, diag has codes
if(!is.null(arguments$hesin_diag) & !is.null(arguments$HESIN)) {
  HES <- fread(arguments$hesin_diag)
  HES_dates <- fread(arguments$HESIN) %>% 
    select(eid, ins_index, epistart,admidate) %>% 
    na_if("") %>% 
    mutate(dated=ifelse((is.na(epistart)==T & is.na(admidate)==T), NA,ifelse(is.na(epistart)==T, admidate, epistart))) %>% 
    mutate(dates=dmy(dated)) %>% 
    select(eid, ins_index, dates) %>% 
    drop_na()
  # combine
  HES_all <- HES %>% 
    left_join(HES_dates) %>% 
    drop_na(dates) %>% 
    mutate_all(na_if,"") %>% 
    select(-diag_icd9_nb,-diag_icd10_nb,-arr_index) %>% 
    filter(!eid %in% exclusions) %>% 
    mutate(diag_icd10= ifelse(nchar(diag_icd10) > 4, (str_extract(diag_icd10, "^.{4}")), diag_icd10),
           diag_icd9= ifelse(nchar(diag_icd9) > 5, (str_extract(diag_icd9, "^.{5}")), diag_icd9)) 
  # split by ICD10/9 and by poistion primary secondary and external (1,2,3)
  # Primary
  ICD10_HES_first <-  HES_all %>%
    filter(level==1) %>%  
    mutate(source="ICD10_1") %>% 
    select(eid,code=diag_icd10,date=dates,source,grouping_variable=ins_index) %>% 
    drop_na(code)
  # Secondary
  ICD10_HES_second <-  HES_all %>%
    filter(level==2) %>% 
    mutate(source="ICD10_2") %>% 
    select(eid,code=diag_icd10,date=dates,source,grouping_variable=ins_index) %>% 
    drop_na(code)
  # External
  ICD10_HES_third <-  HES_all %>%
    filter(level==3) %>%  
    mutate(source="ICD10_3") %>% 
    select(eid,code=diag_icd10,date=dates,source,grouping_variable=ins_index) %>% 
    drop_na(code)
  # Primary
  ICD9_HES_first <-  HES_all %>%
    filter(level==1) %>% 
    mutate(source="ICD9_1") %>% 
    select(eid,code=diag_icd9,date=dates,source,grouping_variable=ins_index) %>% 
    drop_na(code) %>% 
    mutate(code=as.character(code))
  # Secondary
  ICD9_HES_second <-  HES_all %>%
    filter(level==2) %>% 
    mutate(source="ICD9_2") %>% 
    select(eid,code=diag_icd9,date=dates,source,grouping_variable=ins_index) %>% 
    drop_na(code) %>% 
    mutate(code=as.character(code))
  # External
  ICD9_HES_third <-  HES_all %>%
    filter(level==3) %>% 
    mutate(source="ICD9_3") %>% 
    select(eid,code=diag_icd9,date=dates,source,grouping_variable=ins_index) %>% 
    drop_na(code) %>% 
    mutate(code=as.character(code))
  # combine
  all_phenotype_data <- all_phenotype_data  %>% 
    bind_rows(ICD10_HES_first,ICD10_HES_second,ICD10_HES_third,ICD9_HES_first,ICD9_HES_second,ICD9_HES_third)
}
# HES_op ------------------------------------------------------------------
if(!is.null(arguments$hesin_oper)) {
  HES_op <- fread(arguments$hesin_oper)
  HES_op_data <- HES_op %>% 
    filter(!eid %in% exclusions) %>% 
    left_join(HES_dates) %>% 
    drop_na(dates) %>% 
    mutate(source="OPCS") %>% 
    select(eid,code=oper4,date=dates,source) %>% 
    drop_na(code)
  # combine
  all_phenotype_data <- all_phenotype_data  %>% 
    bind_rows(HES_op_data)
}
# Mortality data ----------------------------------------------------------
if(!is.null(arguments$death_cause) & !is.null(arguments$death)) {
  cause_of_death <- fread(arguments$death_cause)
  date_of_death <- fread(arguments$death)
  MD <- cause_of_death %>% 
    left_join(date_of_death) %>% 
    mutate(date_of_death=dmy(date_of_death)) %>% 
    drop_na(date_of_death)
  ## check to see if there is anyone with multiple dates of death
  MD_date_check <- MD %>% 
    group_by(eid) %>% 
    summarise(differnce=max(date_of_death)-min(date_of_death)) %>% 
    filter(differnce != 0) %>% 
    pull(eid)
  # edit ICD10 codes the same as HES to allow maximum mapping and exclude the ids with multiple dates of deaths 
  MD_data <- MD %>% 
    filter(!eid %in% MD_date_check) %>% 
    filter(!eid %in% exclusions)  %>% 
    mutate(code= ifelse(nchar(cause_icd10) > 4, (str_extract(cause_icd10, "^.{4}")), cause_icd10),
           source="MD") %>% 
    select(eid,code,date=date_of_death,source) %>% 
    drop_na(code)
  # combine
  all_phenotype_data <- all_phenotype_data  %>% 
    bind_rows(MD_data)
}
## combine and save
all_phenotype_data <- all_phenotype_data %>% 
  select(-remove)
fwrite(all_phenotype_data,paste0(save_location,"/health_records.txt.gz"),sep = "\t",na = NA)
# Control populations -----------------------------------------------------
if(all(c("combined_sex","GP_C_ID.txt.gz") %in% list.files(here("data")))){
  # making populations from combined sex (as should eb everybody)
  combined_sex_edit <- fread(here("data","combined_sex")) %>% 
    mutate(all_pop=1) %>% 
    select(eid,all_pop)
  # and GP data (required to make effect controls and minimise misclassification of controls)
  GP_ID_edit <- fread(here("data","GP_C_ID.txt.gz")) %>% 
    mutate(primary_care_pop=1) %>% 
    select(eid,primary_care_pop)
  # combine
  control_populations <- GP_ID_edit %>% 
    full_join(combined_sex_edit)
  # write
  fwrite(control_populations,paste0(save_location,"/control_populations"), sep = "\t")
}
# DOB ---------------------------------------------------------------------
if(length(str_which(available_tab_data,"^52-"))>0 & length(str_which(available_tab_data,"^34-"))>0){
  # making approximate DOB using year and month of birth set day at 15th
  # year
  YOB <- tab_data %>% 
    select(eid, matches("^34-")) %>% 
    rename(YOB=2)
  # month
  MOB <- tab_data %>% 
    select(eid, matches("^52-")) %>% 
    rename(MOB=2)
  # DOB
  DOB <- YOB %>% 
    left_join(MOB) %>% 
    mutate(day=15,
           DOB=make_date(day=day,month=MOB,year=YOB)) %>% 
    select(eid,DOB)
# write
  fwrite(DOB,paste0(save_location,"/DOB"))
  }