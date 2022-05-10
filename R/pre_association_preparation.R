#! /usr/bin/env Rscript
'
Final preparation step for the phenotypes prior to association testing, this will create tables for the phenotypes which are then used in the association tests. The tables can be split by group if groupings input was used in phenotype_preparation.R and the quantitative traits transformed inverse normal rank transformation. 

The user must input a location of a phenotype RDS file which will be an output of the phenotype_preparation.R script, a folder name to store the phenotype tables and a root name for the phenotype table saved files. If selecting a grouping file this must be the same as the one used in phenotype_preparartion.R. IVNT will perform inverse normal rank transformation on the quantitative phenotypes. If sex split phenotypes were created in the phenotype_preparation.R script then the user can choose to keep all the phenotypes created (all_sex) or provide a list for inclusion (sex_list). 

Usage:
    pre_association_preparation.R (--phenotypes=<FILE> --folder_name=<text> --file_name=<text>) [--IVNT --groupings=<file> PheWAS_manifest_overide=<FILE>] [--all_sex | --sex_list=<file>]
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    
    Mandatory inputs
    --phenotypes=<FILE>               Full file path of RDS file containing phenotypes processed with phenotype_preparation.R 
    --folder_name=<FOLDER>            Full file path of the folder where phenotype tables will be saved. 
    --file_name=<text>                Common name of the saved files, saved as grouping_file_name. Example if grouping was EUR and file name 
                                      IVNT_relate_removed then final name would be EUR_IVNT_relate_removed.
    
    Options
    --IVNT                            Input if wanting to perform inverse normal transformation on quantitative phenotypes.
                                      
    --groupings=<FILE>                Full file path of file containing group information used for separating out the population,
                                      normally representing ancestry. If selected all functions will be performed on a per group level,
                                      If not selected will use all participants and label output lists as all_pop.
                                      
    --all_sex                         Input if wanting to include all sex specific phenotypes.
    --sex_list=<FILE>                 Input full file path of a file containing sex specific phenotypes for inclusion. Format is a single column
                                      with any header, listing the sex specific phenotypes. Requires sex-specific names, example P350_male not 
                                      just P350.
    --PheWAS_manifest_overide=<FILE>  Full file path of the alternative PheWAS_manifest file.
    
' -> doc


suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.1 pre_association_preparation.R')

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(RNOmni))
suppressMessages(library(parallel))
library(here)


# Functions ---------------------------------------------------------------
join_cols_pheno <- function (x) {
  per_list <- x
  per_list <- per_list[names(per_list) %in% included_phenotypes$PheWAS_ID == TRUE]
  
  PheWAS_ID <<- names(per_list)
  final_phenotypes <- mapply(converted_for_association,per_list,PheWAS_ID,SIMPLIFY = F,USE.NAMES = T) 
  
  if (arguments$IVNT) {
    
    PheWAS_IDs <- names(final_phenotypes)
    transformed_phenoptypes_final <- mapply(IVNT_transformation,final_phenotypes,PheWAS_IDs,SIMPLIFY = F,USE.NAMES = T) %>% 
      reduce(full_join)
    
    return(transformed_phenoptypes_final)
    
  } else {
    
    final_phenotypes_tabled <- final_phenotypes %>% 
      reduce(full_join)
    
    return(final_phenotypes_tabled)
    
  }
  
}
converted_for_association <- function (x,y) {
  name <- str_remove(y,"_male|_female")
  analysis_type <- PheWAS_manifest %>% 
    filter(PheWAS_ID=={{name}}) %>% 
    select(analysis) %>% 
    pull()
  
  age_col <- PheWAS_manifest %>% 
    filter(PheWAS_ID=={{name}}) %>% 
    select(age_col) %>% 
    pull()
  
  if (analysis_type == "binary") {
    phenotype <- x %>% 
      rename(any_code=2) %>% 
      mutate(!!y:=ifelse(is.na(any_code),NA,ifelse(any_code>0,1,0))) %>% 
      select(eid,{{y}}) %>% 
      drop_na()
    return(phenotype)
    
  } else {
    if(str_detect(y,"_male|_female") & str_remove(y,"_male|_female") %in% PheWAS_ID) {
      desired_col_names <- c("eid",y)
    } else {
      if(age_col=="named") {
        desired_col_names <- c("eid",y,paste0(name,"_age"))
      } else {
        desired_col_names <- c("eid",y)
      }
    }
    phenotype <- x %>% 
      rename(!!y:=2) %>% 
      select(any_of(desired_col_names)) %>% 
      drop_na()
    
    return(phenotype)
  }
}
IVNT_transformation <- function(x,y) {
  name <- str_remove(y,"_male|_female")
  
  age_col <- PheWAS_manifest %>% 
    filter(PheWAS_ID=={{name}}) %>% 
    select(age_col) %>% 
    pull()
  
  type <- PheWAS_manifest %>% 
    filter(PheWAS_ID=={{name}}) %>% 
    select(transformation) %>% 
    pull()
  
  if (is.na(type)) {
    return(x)
  } else if (type == "IVNT") {
    
    pheno <- x %>% 
      rename(any_code=2)
    pheno_col <- pheno %>% 
      pull(any_code)
    
    message(y)
    
    IVNT <- RankNorm(pheno_col)
    
    IVNT_col <- replace(pheno[,2],values = IVNT) %>% 
      rename(!!y:=1)
    
    if(age_col=="named"){
      desired_col_names <- c("eid",y,paste0(name,"_age"))
    } else {
      desired_col_names <- c("eid",y)
    }
    pheno_new <- pheno %>%
      bind_cols(IVNT_col) %>% 
      select(any_of(desired_col_names))
    
    return(pheno_new)
  }
}
saving_files <- function(a,b) {
  
  phenos_table <- a 
  name <- paste0(b,"_",arguments$file_name)
  fwrite(phenos_table,paste0(arguments$folder_name,"/",name),sep = "\t",na = NA)
}


# Load in defining variables ----------------------------------------------
if(is.null(arguments$PheWAS_manifest_overide)){
  PheWAS_manifest <- fread(here("data","PheWAS_manifest.csv"))
} else {
  PheWAS_manifest <- fread(arguments$PheWAS_manifest_overide)
}
phenotype_list <- readRDS(arguments$phenotypes)
if(!is.null(arguments$groupings)){
  population_grouping <- fread(arguments$groupings) 
  file_prefix <- discard(unique(population_grouping$group),is.na)
} else {
  file_prefix <- c("all_pop")
}

if(isFALSE(arguments$all_sex) & is.null(arguments$sex_list)) {
  included_phenotypes <- PheWAS_manifest %>% 
    filter(included_in_analysis == 1) 
  
} else if(arguments$all_sex) {
  included_phenotypes_manifest <- PheWAS_manifest %>% 
    filter(included_in_analysis == 1) %>% 
    select(PheWAS_ID) %>% 
    mutate(male=paste0(PheWAS_ID,"_male"),
           female=paste0(PheWAS_ID,"_female"))
  included_phenotypes <- (data.frame(PheWAS_ID=c(included_phenotypes_manifest$PheWAS_ID,
                                                 included_phenotypes_manifest$male,
                                                 included_phenotypes_manifest$female)))
  
} else if(!is.null(arguments$sex_list)) {
  
  included_phenotypes_manifest <- PheWAS_manifest %>% 
    filter(included_in_analysis == 1) %>% 
    select(PheWAS_ID)
  
  sex_phenos <- fread(arguments$sex_list) %>% 
    rename(PheWAS_ID=1) 
  
  included_phenotypes <- sex_phenos %>% 
    bind_rows(included_phenotypes_manifest)
}

if(!dir.exists(arguments$folder_name)){
  dir.create(arguments$folder_name,recursive = T)
}

# Run functions -----------------------------------------------------------
final_phenotypes_columns <- mapply(join_cols_pheno,phenotype_list,SIMPLIFY = F,USE.NAMES = T)
mapply(saving_files,final_phenotypes_columns,file_prefix, SIMPLIFY = F)
