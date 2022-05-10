#! /usr/bin/env Rscript
'
Takes the created phenotypes and filters for case number and optionally for relatedness. Can filter the whole sample or perform filtering by group. These groups are classically ancestries produces a list for each grouping. Can also create sex specific phenotypes (1 each for male and female), these can be selected in subsequent scripts if not wanting all sex_split phenotypes. At each filtering stage will save an R object list containing filtered phenotypes as a list of lists.

All inputs are optional, by default the script will filter for default case numbers in the phenotypes marked for inclusion by the PheWAS_manifest and save the output, edit these default values with quantitative_Case_N and binary_Case_N inputs and use no_case_N_save if not wanting to save the output of the initial case number filtering. If wanting to filter by group use groupings input, if wanting to filter for relatedness then use relate_remove and kinship inputs, to save to none default name use relate_remove_save_name. T add in additional sex specific phenotypes choose the sex_split input alongside the sex_info file and what value male and females are coded in using male and female inputs.

Usage:
    phenotype_preparation.R [--phenotype_folder=<FOLDER> --groupings=<file> --quantitative_Case_N=<number> --binary_Case_N=<number> --no_case_N_save --case_N_filtered_save_name=<name> --relate_remove --relate_remove_save_name=<name> --kinship_file=<file> --sex_split --sex_info=<file> --male=<number> --female=<number> --PheWAS_manifest_overide=<FILE>]
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    
    
    
    Filtering by group
    --groupings=<FILE>                      Full file path of file containing group information used for seperating out the population,
                                            normally representing ancestry. If selected all filtering will be performed on a per group level,
                                            If not used will use all participants and label output lists as all_population.
    
    Case number filtering options
    --phenotype_folder=<FOLDER>             Full file path of the folder containing the phenotypes RDS files.
                                            [default: data/phenotypes/]
    --quantitative_Case_N=<number>          Number that represents the minimum number of cases for quantitative phenotype inclusion. 
                                            [default: 100]
    --binary_Case_N=<number>                Number that represents the minimum number of cases for binary phenotype inclusion. 
                                            [default: 50]
    --no_case_N_save                        Input if not wanting to save the output of the filtering by number of cases.
    --case_N_filtered_save_name=<FILE>      Full file path for save location of the case_N filtered R list. If not inputted will save to default. 
                                            [default: analysis/phenotype/preperation/all_pheno_epi_case_N_filtered]
    Realte remove filtering                                        
    --relate_remove                         Input option to perform relate remove functions to produce list of phenotypes filtered for 
                                            relateness and number of cases.
    --relate_remove_save_name=<FILE>        Full file path for save location of the relate_remove R list created when selected relate_remove option.
                                            If not inputted will save to default name in analysis/phenotype_preparation folder. 
                                            [default: analysis/phenotype/preperation/relate_remove_pheno_N_filtered]
    --kinship_file=<FILE>                   Full file path to the kinship file required if wanting to filter by relateness. 
                                            See user guide for required format. [default: data/related_callrate]
    
    Sex split phenotype options                                        
    --sex_split                             Input option if wanting to create sex specific phenotypes, will derive and add a male and female version
                                            of each phenotype. Further scripts allow user to select specific phenotypes for inclusion.
    --sex_info=<FILE>                       Full file path for file containing sex data, see user guide for format. [default: data/combined_sex]
    --male=<number>                         Value of male sex coded in table must be a number. [default: 1]
    --female=<number>                       Value of female sex coded in table must be a number. [default: 0]
    --PheWAS_manifest_overide=<FILE>        Full file path of the alternative PheWAS_manifest file.
   
' -> doc

suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.1 phenotype_preparation')

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
library(here)

# Functions ---------------------------------------------------------------
sum_stats <- function(x,y,z) {
  name <- str_remove(z,"_male|_female")
  
  type <- PheWAS_manifest %>% 
    filter(PheWAS_ID=={{name}}) %>% 
    select(analysis) %>% 
    pull()
  
  #message(paste(type))
  
  if (type=="quant") {
    
    cases <- x %>% 
      drop_na(2)
    
    n_cases <- nrow(cases)
    n_control <- NA
    
    phecode_name <- z
    stats_name <- paste0("N_cases_",y)
    control_name <- paste0("N_control_",y)
    
    
    summary_sample_stats <- data.frame(phecode_name,n_cases,n_control) %>%
      rename(PheWAS_ID=1,!!stats_name:=2,!!control_name:=3)
    
  } else {
    
    ## epi cases and controls
    cases <- x
    n_cases <- length(which(cases[,2]>0))
    n_control <- length(which(cases[,2]==0))
    
    phecode_name <- z
    stats_name <- paste0("N_cases_",y)
    control_name <- paste0("N_control_",y)
    
    summary_sample_stats <- data.frame(phecode_name,n_cases,n_control) %>%
      rename(PheWAS_ID=1,!!stats_name:=2,!!control_name:=3)
    
  }
  
  return(summary_sample_stats)
}
group_id_creater <- function(x) {
  ID <- groupings %>% 
    filter(group==x) %>% 
    pull(eid)
  return(ID)
}
relate_remove <- function(x,y) {
  
  name <- str_remove(y,"_male|_female")
  
  type <- PheWAS_manifest %>% 
    filter(PheWAS_ID=={{name}}) %>% 
    select(analysis) %>% 
    pull()
  
  pheno <- x %>% 
    rename(any_code=2)
  
  if (type=="quant") {
    
    rel<-call_rate_kinships
    
    pheno<-x %>% 
      drop_na()
    
    # Remove all rows from the relatedness file where (1) one or both IDs are not in the phenotype 
    #(1)
    rel <- call_rate_kinships %>% 
      filter(ID1 %in% pheno$eid & ID2 %in% pheno$eid)
    
    # Create identifier for whether pair contains a duplicated ID or not and partition into two data frames (independent pairs and pairs containing duplicated IDs)
    
    AllID<-c(rel$ID1,rel$ID2)
    dupID<-AllID[which(duplicated(AllID))]
    rel$dup<-"NO"
    rel$dup[which(rel$ID1%in%dupID|rel$ID2%in%dupID)]<-"YES"
    
    indpairs<-rel[which(rel$dup=="NO"),]
    duppairs<-rel[which(rel$dup=="YES"),]
    
    # How many times does each duplicated ID arise?
    
    repeat_counter <- function(x) {
      
      duplicated_count_ID1 <- rel %>% 
        filter(dup=="YES") %>% 
        select(ID=ID1)
      
      duplicated_count_ID2 <- rel %>% 
        filter(dup=="YES") %>% 
        select(ID=ID2)
      
      duplicated_count_all <- duplicated_count_ID1 %>% 
        bind_rows(duplicated_count_ID2) %>% 
        group_by(ID) %>% 
        count()
      
      duplicated_n <- duplicated_count_all %>% 
        filter(n>x) %>% 
        pull(ID)
      
      return(duplicated_n)
      
    }
    
    over_6 <- repeat_counter(5)
    rel <- rel %>% 
      filter(!ID1 %in% over_6 & !ID2 %in% over_6)
    
    over_5 <- repeat_counter(4)
    rel <- rel %>% 
      filter(!ID1 %in% over_5 & !ID2 %in% over_5)
    
    over_4 <- repeat_counter(3)
    rel <- rel %>% 
      filter(!ID1 %in% over_4 & !ID2 %in% over_4)
    
    over_3 <- repeat_counter(2)
    rel <- rel %>% 
      filter(!ID1 %in% over_3 & !ID2 %in% over_3)
    
    over_2 <- repeat_counter(1)
    rel <- rel %>% 
      filter(!ID1 %in% over_2 & !ID2 %in% over_2)
    
    exlcusion_over <- c(over_6,over_5,over_4,over_3,over_2)
    
    remaining_ID <- rel %>% 
      mutate(excluded_ID = case_when(lower_missing==1 ~ ID2,
                                     is.na(missingness_ID2) ~ ID2,
                                     lower_missing==2 ~ ID1,
                                     is.na(missingness_ID1) ~ ID1)) %>% 
      pull(excluded_ID)
    all_exclusions <- c(remaining_ID,exlcusion_over)  
    
    pheno_new <- pheno %>% 
      filter(!eid %in% all_exclusions)
    
    return(pheno_new)
    
  } else {
    
    ## now related removed
    cases <- pheno %>% 
      filter(any_code>0) %>% 
      pull(eid)
    
    controls <- pheno %>% 
      filter(any_code==0) %>% 
      pull(eid)
    
    ## need to remove related people for each phenotype
    ## edit to remove exclusion from kinship file
    call_rate_kinship <- call_rate_kinships %>% 
      filter(ID1 %in% pheno$eid & ID2 %in% pheno$eid)
    
    ## remove the case with the lowest call rate
    case_case <- call_rate_kinship %>% 
      filter(ID1 %in% cases, ID2 %in% cases) %>% 
      mutate(exclusion=if_else(lower_missing==1,ID2,ID1)) %>% 
      select(exclusion) %>% 
      pull()
    
    ## edit the callrate_kinship file to remove any excluded cases
    call_rate_kinship <- call_rate_kinship %>% 
      filter(!ID1 %in% case_case) %>% 
      filter(!ID2 %in% case_case)
    
    ## now remove controls, here controls are defined simply as not being exclusions or cases
    case_control <- call_rate_kinship %>% 
      filter(ID1 %in% cases & ID2 %in% controls) %>% 
      mutate(exclusion=ID2) %>% 
      select(exclusion) %>% 
      pull()
    
    ## edit the callrate_kinship file to remove any excluded controls
    call_rate_kinship <- call_rate_kinship %>% 
      filter(!ID1 %in% case_control) %>% 
      filter(!ID2 %in% case_control)
    
    control_case <- call_rate_kinship %>% 
      filter(ID1 %in% controls & ID2 %in% cases) %>% 
      mutate(exclusion=ID1) %>% 
      select(exclusion) %>% 
      pull()
    
    ## edit the callrate_kinship file to remove any excluded controls 
    call_rate_kinship <- call_rate_kinship %>% 
      filter(!ID1 %in% control_case) %>% 
      filter(!ID2 %in% control_case)
    
    ## 
    control_control <- call_rate_kinship %>% 
      filter(ID1 %in% controls & ID2 %in% controls) %>% 
      mutate(exclusion=if_else(lower_missing==1,ID2,ID1)) %>% 
      select(exclusion) %>% 
      pull()
    
    ## combining the vectors
    excluded_cases <- unique(c(case_case))
    excluded_controls <- unique(c(case_control,control_case,control_control))
    removed_ID <- c(excluded_cases,excluded_controls)
    
    related_removed_final <- pheno %>% 
      filter(!eid %in% removed_ID) 
    
    return(related_removed_final)
  } 
}


# Load in defining variables ----------------------------------------------
if(is.null(arguments$PheWAS_manifest_overide)){
  PheWAS_manifest <- fread(here("data","PheWAS_manifest.csv"))
} else {
  PheWAS_manifest <- fread(arguments$PheWAS_manifest_overide)
}
if(arguments$phenotype_folder=="data/phenotypes/"){
  phenotype_folder <- here("data","phenotypes")
} else {
  phenotype_folder <- arguments$phenotype_folder
}

quant_min_cases <- as.numeric(arguments$quantitative_Case_N)
binary_min_cases <- as.numeric(arguments$binary_Case_N)


# load optional arguments
if(!is.null(arguments$groupings)) {
  groupings <- fread(arguments$groupings)
}

if(arguments$kinship_file=="data/related_callrate"){
  call_rate_kinships <- fread(here("data","related_callrate"))  
  } else {
    call_rate_kinships <- fread(arguments$kinship_file)
    }

# create directory

if(arguments$case_N_filtered_save_name=="analysis/phenotype/preperation/all_pheno_epi_case_N_filtered") {
  if(!dir.exists(here("analysis","phenotypes","phenotype_preparation"))){
    dir.create(here("analysis","phenotypes","phenotype_preparation"),recursive = T)
  }
  case_N_save_name <- here("analysis","phenotype","preperation","all_pheno_epi_case_N_filtered")
} else {
  case_N_save_name <- arguments$case_N_filtered_save_name
  new_folder <- str_remove(arguments$case_N_filtered_save_name,basename(arguments$case_N_filtered_save_name))
  if(!dir.exists(new_folder)){
    dir.create(new_folder,recursive = T)
  }
  
  }


# read in all the files in the phenotypes folder and combine them.
files <- paste0(phenotype_folder,"/",list.files(phenotype_folder))

phenotypes <- files %>% 
  map(readRDS) %>% 
  reduce(c)

# first edit the phenotypes to those included in the phenotype manifest for inclusion 
included_phenotypes <- PheWAS_manifest %>% 
  filter(included_in_analysis == 1) 

all_phenotypes <- phenotypes[names(phenotypes) %in% included_phenotypes$PheWAS_ID == TRUE]

# optional addition of sex_split phenotypes
if(arguments$sex_split){
  
  male_N <- as.numeric(arguments$male)
  female_N <- as.numeric(arguments$female)
  
  if(arguments$sex_info=="data/combined_sex"){
    sex_info <- fread(here("data","combined_sex")) %>% 
      rename(eid=1) }
    else{
      sex_info <- fread(arguments$sex_info)
    }
    
    male <- sex_info %>% 
      filter(sex==male_N) %>% 
      pull(eid)
    female <- sex_info %>% 
      filter(sex==female_N) %>% 
      pull(eid)
  
  sex_specific_female <- function(x,y) {
    female_pheno <- y %>%
      ungroup() %>% 
      filter(eid %in% female) %>% 
      rename(!!x:=2)
    return(female_pheno)
  }
  
  sex_specific_male <- function(x,y) {
    male_pheno <- y %>% 
      ungroup() %>% 
      filter(eid %in% male) %>% 
      rename(!!x:=2)
    return(male_pheno)
  }
  
  pre_sex_split <- PheWAS_manifest %>% 
    filter(sex=="Male"|sex=="Female")
  
  all_phenotypes_sex_edit<-all_phenotypes[names(all_phenotypes) %in% pre_sex_split$PheWAS_ID == F]
  
  male_names <- paste0(names(all_phenotypes_sex_edit),"_male")
  female_names <- paste0(names(all_phenotypes_sex_edit),"_female")
  
  female_phenotypes <- mapply(sex_specific_female,female_names,all_phenotypes_sex_edit, SIMPLIFY = F, USE.NAMES = T)
  male_phenotypes <- mapply(sex_specific_male,male_names,all_phenotypes_sex_edit, SIMPLIFY = F, USE.NAMES = T)
  
  all_phenotypes <- c(all_phenotypes,female_phenotypes,male_phenotypes)
  
}

# now split the phenotypes into ancestry or any other division
if(!is.null(arguments$groupings)){
  
  groups <- discard(unique(groupings$group),is.na)
  group_ID <- mapply(group_id_creater, groups,USE.NAMES = T,SIMPLIFY = F)
  
  spliting_by_population <- function(x) {
    
    per_pheno_ancestry <- function (x,y) {
      
      split_pheno <- x %>%
        ungroup() %>% 
        filter(eid %in% y)
      
      return(split_pheno)
    }
    
    per_population <- mapply(per_pheno_ancestry,all_phenotypes,SIMPLIFY = F, USE.NAMES = T,MoreArgs = list(y=x))
    
    return(per_population)
  }
  
  per_population_phenotypes <- mapply(spliting_by_population,group_ID,SIMPLIFY = F,USE.NAMES = T)
  
} else {
  groups <- c("all_population")
  per_population_ungroup <- all_phenotypes %>% 
    ungroup()
  per_population_phenotypes <- list(all_pop=per_population_ungroup)
}

# Run functions -----------------------------------------------------------
# now run summary stats on those phenotypes using whatever groupings selected here we use ancestry
summary_all_pops <- function(x,y,z) {
  
  population_name <- y
  pheno_names <- names(x)
  
  summary_N <- mapply(sum_stats,x,z=pheno_names,SIMPLIFY = F,USE.NAMES = T,MoreArgs = list(y=population_name)) %>% 
    reduce(bind_rows)
  
  return(summary_N)
  
} 

per_population_summary <- mapply(summary_all_pops,per_population_phenotypes,groups,SIMPLIFY = F,USE.NAMES = T)

# now exclude those with <N cases <N cases quantitative. 
selecting_phenotype_cases_N <- function(x) {
  
  PheWAS_ID_type <- PheWAS_manifest %>% 
    select(PheWAS_ID,analysis)
  
  filtering_process <- x %>% 
    rename(cases=2) %>% 
    mutate(ID_2=str_remove(PheWAS_ID,"_male|_female")) %>% 
    left_join(PheWAS_ID_type, by=c("ID_2"="PheWAS_ID")) %>%
    mutate(sufficent_cases=ifelse(analysis=="quant" & 
                                    cases >= quant_min_cases,1,ifelse(analysis=="binary" &
                                                                      cases >= binary_min_cases, 1, 0))) %>% 
    filter(sufficent_cases>0)
  
  N_cases_filtered <- x %>% 
    filter(PheWAS_ID %in% filtering_process$PheWAS_ID) %>% 
    pull(PheWAS_ID)
  
  return(N_cases_filtered)
}

trimmed_PheWAS_ID <- mapply(selecting_phenotype_cases_N,per_population_summary,SIMPLIFY = F,USE.NAMES = T)

reducing_phenotypes <- function(x,y) {
  
  x <- x[names(x) %in% y == TRUE]
  return(x)
}

epi_case_trimmed <- mapply(reducing_phenotypes,per_population_phenotypes,trimmed_PheWAS_ID,SIMPLIFY = F)

if(isFALSE(arguments$no_case_N_save)) {
      saveRDS(epi_case_trimmed,case_N_save_name)
  }
# now to remove related individuals within each population.
if(arguments$relate_remove) {
  relate_remove_all_pops <- function(x) {
    
    PheWAS_ID <- names(x)
    
    relate_remove_step <- mapply(relate_remove,x,PheWAS_ID,SIMPLIFY = F,USE.NAMES = T)
    
    return(relate_remove_step)
    
  } 
  
  relate_remove_phenotypes <- mapply(relate_remove_all_pops,epi_case_trimmed,SIMPLIFY = F,USE.NAMES = T)
  
  RR_population_summary <- mapply(summary_all_pops,relate_remove_phenotypes,groups,SIMPLIFY = F,USE.NAMES = T)
  RR_trimmed_PheWAS_ID <- mapply(selecting_phenotype_cases_N,RR_population_summary,SIMPLIFY = F,USE.NAMES = T)
  RR_case_trimmed <- mapply(reducing_phenotypes,relate_remove_phenotypes,RR_trimmed_PheWAS_ID,SIMPLIFY = F)
  
  if(arguments$relate_remove_save_name=="analysis/phenotype/preperation/relate_remove_pheno_N_filtered") {
    if(!dir.exists(here("analysis","phenotype","preperation"))) {
      dir.create(here("analysis","phenotype","preperation"))
    }
    relate_removed_save <- here("analysis","phenotype","relate_remove_pheno_N_filtered")
    
    }else{
      relate_removed_save <- arguments$relate_remove_save_name
      new_folder <- str_remove(arguments$relate_remove_save_name,basename(arguments$relate_remove_save_name))
      if(!dir.exists(new_folder)){
        dir.create(new_folder,recursive = T)
      }
    }
  saveRDS(RR_case_trimmed,relate_removed_save)
}

