#! /usr/bin/env Rscript
'Runs association analysis in using plink2. Takes the variants extracted using extract_snps.R and performs regression analysis on the availble phenotypes per inputted group. 

Requires location of a folder where analysis will be hosted, a comma separated list of phenotype files derived from pre_association_preparation.R, a covariate file edited for use in plink (see user guide), and the location of the variants for association prepared with extract_SNPs.R. By default, it will run association tests on all included phenotypes in the phenotype files used and use those phenotype files to assign group names, all analysis is then performed per group. Phenotypes can be specified using one of three flags phenotype_overide, exclude_sex_phenotypes or sex_include_file, group names can be specified using the group_name_overide flag. Results are saved in a R list with option to save tables per-group of the combined raw plink results. With a large number of variants the regression analysis can take a long time, to make for more efficient analysis the phenotypes for any one or more of the analysed groups can be split using the split_group input. When a group is split the phenotype files are split into smaller chunks, the analysis can then be performed in a cluster enviroment.To perfomr the analysis this same script can be used again but with the split_analysis flag used to amend how results are saved. Results are then combined with  split_combine.R. If splitting analysis tables and figures should not be compiled until split_combine.R is complete.   

Usage:
PheWAS_association_PLINK.R  (--analysis_folder=<FOLDER> --covariate=<file> --phenotype_files=<FILE> --variants_for_association=<FILE> --analysis_name=<name>) [--group_name_overide=<text>] [--split_analysis --N_cores=<number> --PheWAS_manifest_overide=<FILE> --plink_exe=<command> --save_plink_tables --split_group=<name> --N_quant_split=<number> --N_binary_split=<number> --model=<text> --check_existing_results] [--phenotype_overide=<FILE> | --exclude_sex_phenotypes | --sex_include_file=<file>]
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    
    Mandatory inputs
    --analysis_folder=<FOLDER>          Full path to the folder location of the folder hosting this analysis.
    --phenotype_files=<FILE>            Comma separated full file paths of the phenotype data. For example 
                                        /home/phenotypes/EUR_pheno,/home/phenotypes/AFR_pheno 
    --covariate=<file>                  Full file path to the covariate file edited for use in plink, see user guide.
    --variants_for_association=<FILE>   Full file path to the pgen file for the variants used for association, generated with extract_SNP.R                         --analysis_name=<name>              Name for the analysis, is used later in saving tables, so should distignuish between other analyses.
    
    Options
    --group_name_overide=<text>         Comma separated text input for group names. By default group names are created by extracting all text 
                                        before the first underscore (_) from the file names (note not the full path but the files names) of 
                                        the inputs from phenotype_files. Example /home/phenotypes/EUR_pheno would extract EUR as the group name. 
                                        The order for names must match the order of the input for phenotype_files.
    
    --N_cores=<number>                  Number of cores requested if wanting to use parallel computing.
    --PheWAS_manifest_overide=<FILE>    Full file path of the alternative PheWAS_manifest file.
    --plink_exe=<command>               Command to execute plink2 program [default: plink2]
    --save_plink_tables                 Select if wanting to save the raw plink results per group as a table always saved in 
                                        analysis_folder/assoiation_results/group/group_plink_results_raw
                                        
    --phenotype_overide=<FILE>          Full file path to a txt file containing single column of phenotypes used to filter for inclusion, 
                                        default uses PheWAS_manifest alongside any additonal sex specific phenotypes created in the previous stages.
    --exclude_sex_phenotypes            Select if wanting to exclude the sex specific phenotypes from the analysis.
    --sex_include_file=<file>           Full file path to list of desired included sex specific phenotypes.
  
    --split_group=<name>                Comma separated groups that require splitting for more efficient analysis. Does not run the analysis per group 
                                        but splits and saves the phenotypes into smaller files which can then be analysed using cluster based computing. 
                                        The number of phenotypes in each split file is dependent on the type (binary or quantitative) and is set using
                                        N_quant_split and N_binary_split. Split files are saved in analysis_folder/group_split with group being the group
                                        name. A file is created in analysis_folder/group_split saved as group_split_guide with group once again being the 
                                        group name inputted. This file lists the file names of the splits, which can be used to guide distributed 
                                        computing analysis.
    --N_quant_split=<number>            Number of quantitative phenotypes per split file. [default: 200]
    --N_binary_split=<number>           Number of binary phenotypes per split file. [default: 80]
    --split_analysis                    Select if the input being analysed is from the split_group command, only changes the saving and attempted combining of                                          files.
    --model=<text>                      Genetic model to use for analysis, can be one of genotypic, hethom, dominant, recessive, hetonly. If not one of 
                                        these options will be converted to blank.
    --check_existing_results            Input if wanting to re-run part of the analysis whilst checking for existing results in the plink_results folder. Used                                          primarily where a run has aborted part way through.
    
' -> doc

suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.1 PheWAS_association_PLINK.R')

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
library(here)
library(lubridate)


# Functions ---------------------------------------------------------------
association_tests <- function(y,x) {
  if(any(str_detect(split_groups,y))) {
    #create directory for the split files
    dir.create(paste0(arguments$analysis_folder,"/","split_",y))
    # edit phenotype file id much like covariates  
    phenotype_edit <- fread(x) %>% 
      rename(IID=1) %>% 
      filter(IID %in% psam$IID)
    # define split values for different types of phenotypes    
    split_number_binary <- as.numeric(arguments$N_binary_split)
    split_number_quant <- as.numeric(arguments$N_quant_split)
    # separate phenotypes into quant and binary    
    binary_phenotypes <- phenotype_edit %>% 
      select(1,any_of(binary_ID))
    quant_phenotypes <- phenotype_edit %>% 
      select(1,any_of(quant_ID))
    # now split PheWAS_ID by the values selected    
    binary_column_splits <- split(colnames(binary_phenotypes)[-1], ceiling(seq_along(colnames(binary_phenotypes)[-1])/split_number_binary))
    quant_column_splits <- split(colnames(quant_phenotypes)[-1], ceiling(seq_along(colnames(quant_phenotypes)[-1])/split_number_quant))
    # combine    
    column_splits <- c(binary_column_splits,quant_column_splits)
    # use PheWAS_ID splits to select off columns    
    split_phenotypes_df <- lapply(column_splits, function(x) split_pheno <- phenotype_edit %>% 
                                    select(1,any_of(x)))
    # create a file for each of the splits 
    N_splits <- length(split_phenotypes_df)
    file_names <- paste0(y,"_",seq(1:N_splits))
    # simple save function
    saving_files <- function(x,y,z) {
      file_to_save <- x
      fwrite(file_to_save,paste0(arguments$analysis_folder,"/","split_",z,"/",y),sep = "\t", na = NA)
    }
    mapply(saving_files,split_phenotypes_df,file_names,y)
    # guide to help with distributed computing
    guide_file <- data.frame(list.files(paste0(arguments$analysis_folder,"/","split_",y),full.names = T))
    fwrite(guide_file,paste0(arguments$analysis_folder,"/",y,"_split_guide"),col.names = F)
    
    return()
  }
  #save location for the plink results 
  save_location <- paste0(arguments$analysis_folder,"/","plink_results","/",y)
  group_name <- y
  #check existing files for phenotypes
  all_pheno_results <- list.files(paste0(arguments$analysis_folder,"/","plink_results"))
  convert_PheWAS_ID <- str_remove(str_remove(all_pheno_results[str_detect(all_pheno_results,y)],paste0(y,".")),".glm.*")
  # ability to just run pheno-group test that havn't already been completed
  if(arguments$check_existing_results){
    remaining_binary <- setdiff(binary_ID,convert_PheWAS_ID)
    remaining_quatitative <- setdiff(quant_ID,convert_PheWAS_ID)
    all_remaining <- c(remaining_binary,remaining_quatitative)
  } else {
    all_remaining <- c(binary_ID,quant_ID)
  }
  # edit phenotype file id much like covariates  
  phenotype_edit <- fread(x) %>% 
    rename(IID=1) %>% 
    filter(IID %in% psam$IID) %>%
    select(IID,any_of(all_remaining))
  # exit if no phenotypes to test
  if(ncol(phenotype_edit)==1){
    if(isFALSE(arguments$split_analysis)){
      # combine results
      results_files <-  list.files(path = paste0(arguments$analysis_folder,"/","plink_results"), pattern = paste0(y,"(.*?)hybrid|",y,"(.*?)linear"))
      remove_string <- c(".glm.logistic.hybrid", paste0(y,"."), ".glm.linear")
      results_file_location <- paste0(paste0(arguments$analysis_folder,"/","plink_results","/",results_files))
      plink_file_edit_join <- function(x,y) {
        name <- data.frame(str_remove_all(y, paste(remove_string, collapse = "|"))) %>% 
          rename(PheWAS_ID=1)
        per_result <- fread(x) %>% 
          bind_cols(name) %>% 
          mutate(across(everything(), as.character),
                 table_save=arguments$analysis_name)
        return(per_result)  
      }
      # join the plink results into a single table
      plink_results <- mapply(plink_file_edit_join,results_file_location,results_files,SIMPLIFY = F) %>% 
        rbindlist(.,fill = T)
      # option to save raw results per group
      if(arguments$save_plink_tables){
        save_results_table_name <- paste0(y,"_",arguments$analysis_name,"_plink_results_raw.csv")
        fwrite(plink_results,paste0(arguments$analysis_folder,"/","association_results","/",y,"/",save_results_table_name), na = NA)
      }
      return(plink_results)
    } else {
      return()
    }
  }
  # seperate into binary and quant this is for neater results
  binary_pheno <- phenotype_edit %>% 
    mutate(`#FID`=IID) %>% 
    select(`#FID`,IID,any_of(binary_ID))
  # quant
  quant_pheno <- phenotype_edit %>% 
    mutate(`#FID`=IID) %>%
    select(`#FID`,IID,any_of(quant_ID))
  # save locations for temp phenotypes
  binary_pheno_file <- paste0(arguments$analysis_folder,"/","temp_pheno","/",basename(x),"_temp_binary")
  quant_pheno_file <- paste0(arguments$analysis_folder,"/","temp_pheno","/",basename(x),"_temp_quant")
  # selecting whether to do binary or quant or both types of analysis
  if(ncol(quant_pheno)==2) {
    # write temp binary file
    fwrite(binary_pheno,binary_pheno_file, sep = "\t", na = NA, quote = F)
    # run the regression in plink2
    system(paste0(plink_exe," --pfile ",association_variants," --glm sex cols=+totallele,+a1count,+a1countcc,+a1countcc,+a1freq,+machr2,+a1freqcc hide-covar ",model," --pheno ",binary_pheno_file," --covar ",covariates," --ci 0.95  --1 --vif 90 --out ",save_location," --covar-variance-standardize"))
    # remove the binary file  
    file.remove(binary_pheno_file)
  } else if(ncol(binary_pheno)==2) {
    fwrite(quant_pheno,quant_pheno_file, sep="\t", na = NA,quote = F)
    # run the regression in plink2
    system(paste0(plink_exe," --pfile ",association_variants," --glm sex cols=+totallele,+a1count,+a1countcc,+a1countcc,+a1freq,+machr2,+a1freqcc hide-covar ",model," --pheno ",quant_pheno_file," --covar ",covariates," --ci 0.95 --1 --vif 90 --out ",save_location," --covar-variance-standardize"))
    # remove the quant file  
    file.remove(quant_pheno_file)
  } else {
    # write temp phenotype files
    fwrite(binary_pheno,binary_pheno_file, sep = "\t", na = NA, quote = F)
    fwrite(quant_pheno,quant_pheno_file, sep="\t", na = NA, quote = F)
    # run the regression in plink2
    system(paste0(plink_exe," --pfile ",association_variants," --glm sex cols=+totallele,+a1count,+a1countcc,+a1countcc,+a1freq,+machr2,+a1freqcc hide-covar ",model," --pheno ",binary_pheno_file," --covar ",covariates," --ci 0.95 --1 --vif 90 --out ",save_location," --covar-variance-standardize"))
    
    system(paste0(plink_exe," --pfile ",association_variants," --glm sex cols=+totallele,+a1count,+a1countcc,+a1countcc,+a1freq,+machr2,+a1freqcc hide-covar ",model," --pheno ",quant_pheno_file," --covar ",covariates," --ci 0.95 --1 --vif 90 --out ",save_location," --covar-variance-standardize"))
    # remove the temp phenotype files 
    file.remove(binary_pheno_file)
    file.remove(quant_pheno_file)
  }
  if(isFALSE(arguments$split_analysis)){
    # combine results
    results_files <-  list.files(path = paste0(arguments$analysis_folder,"/","plink_results"), pattern = paste0(y,"(.*?)hybrid|",y,"(.*?)linear"))
    remove_string <- c(".glm.logistic.hybrid", paste0(y,"."), ".glm.linear")
    results_file_location <- paste0(paste0(arguments$analysis_folder,"/","plink_results","/",results_files))
    plink_file_edit_join <- function(x,y) {
      name <- data.frame(str_remove_all(y, paste(remove_string, collapse = "|"))) %>% 
        rename(PheWAS_ID=1)
      per_result <- fread(x) %>% 
        bind_cols(name) %>% 
        mutate(across(everything(), as.character),
               table_save=arguments$analysis_name)
      return(per_result)  
    }
    # join the plink results into a single table
    plink_results <- mapply(plink_file_edit_join,results_file_location,results_files,SIMPLIFY = F) %>% 
      rbindlist(.,fill = T)
    # option to save raw results per group
    if(arguments$save_plink_tables){
      save_results_table_name <- paste0(y,"_",arguments$analysis_name,"_plink_results_raw.csv")
      fwrite(plink_results,paste0(arguments$analysis_folder,"/","association_results","/",y,"/",save_results_table_name), na = NA)
    }
    return(plink_results)
  } else {
    return()
  }
}

# Load in defining variables ----------------------------------------------
if(is.null(arguments$PheWAS_manifest_overide)){
  PheWAS_manifest <- fread(here("data","PheWAS_manifest.csv"))
} else {
  PheWAS_manifest <- fread(arguments$PheWAS_manifest_overide)
}
# creating analysis folder
if(!dir.exists(arguments$analysis_folder)) {
  dir.create(arguments$analysis_folder,recursive = T)
}
#create folder for temp phenotypes
if(!dir.exists(paste0(arguments$analysis_folder,"/","temp_pheno"))){
  dir.create(paste0(arguments$analysis_folder,"/","temp_pheno"),recursive = T)
}
# assigning groups and phenotype_files dependent on input
phenotype_files <- unlist(str_split(arguments$phenotype_files, ","))
if(!is.null(arguments$group_name_overide)) {
  groups <- unlist(str_split(arguments$group_name_overide, ","))
  lapply(groups, function(x) dir.create(paste0(arguments$analysis_folder,"/","association_results","/",x),recursive = T))
} else {
  groups <-gsub("_.*", "", basename(unlist(strsplit(arguments$phenotype_files,","))))
  lapply(groups, function(x) dir.create(paste0(arguments$analysis_folder,"/","association_results","/",x),recursive = T))
}

if(!is.null(arguments$split_group)) {
  split_groups <- unlist(str_split(arguments$split_group,pattern = ","))
} else {
  split_groups <- ""
}

# variants for association edit to allow input into plink2
association_variants <- str_remove(arguments$variants_for_association,".pgen")
# defining all_phenos
if(!is.null(arguments$phenotype_overide)) {
  phenotype_overide <- fread(arguments$phenotype_overide) %>% 
    pull(1)
  # incorporate age as 
  age_phenos <- PheWAS_manifest %>% 
    filter(age_col!="age") %>% 
    filter(PheWAS_ID %in% phenotype_overide) %>% 
    mutate(age_col_complete = ifelse(age_col=="named",paste0(PheWAS_ID,"_age"),"age")) %>% 
    select(PheWAS_ID=age_col_complete)
  # make sex phenotypes names  
  sex_phenos <- PheWAS_manifest %>% 
    filter(included_in_analysis==1) %>%
    mutate(male=paste0(PheWAS_ID,"_male"),
           female=paste0(PheWAS_ID,"_female"),
           male_analysis=analysis,female_analysis=analysis,
           male_age_col=age_col,female_age_col=age_col) %>% 
    select(PheWAS_ID,male,female,analysis,male_analysis,female_analysis,age_col,male_age_col,female_age_col)
  # make into data.frame  
  all_phenos_complete <- data.frame(PheWAS_ID=c(sex_phenos$PheWAS_ID,sex_phenos$male,sex_phenos$female), analysis=c(sex_phenos$analysis,sex_phenos$male_analysis,sex_phenos$female_analysis), age_col=c(sex_phenos$age_col,sex_phenos$male_age_col,sex_phenos$female_age_col))
  # combine filtering for override file  
  all_phenos <- all_phenos_complete %>%
    filter(PheWAS_ID %in% phenotype_overide) %>% 
    mutate(include=1) %>% 
    bind_rows(age_phenos)
} else if(arguments$exclude_sex_phenotypes) {
  # no sex phenotypes just use PheWAS_manifest
  age_phenos <- PheWAS_manifest %>% 
    filter(age_col!="age") %>% 
    filter(included_in_analysis==1) %>% 
    mutate(age_col_complete = ifelse(age_col=="named",paste0(PheWAS_ID,"_age"),"age")) %>% 
    select(PheWAS_ID=age_col_complete)
  
  all_phenos <- PheWAS_manifest %>% 
    filter(included_in_analysis==1) %>% 
    mutate(include=1) %>% 
    select(PheWAS_ID,analysis,age_col,include) %>% 
    bind_rows(age_phenos)
} else if(is.null(arguments$sex_include_file)) {
  age_phenos <- PheWAS_manifest %>% 
    filter(age_col!="age") %>% 
    filter(included_in_analysis==1) %>% 
    mutate(age_col_complete = ifelse(age_col=="named",paste0(PheWAS_ID,"_age"),"age")) %>% 
    select(PheWAS_ID=age_col_complete)
  
  sex_phenos <- PheWAS_manifest %>% 
    filter(included_in_analysis==1) %>%
    mutate(male=paste0(PheWAS_ID,"_male"),
           female=paste0(PheWAS_ID,"_female"),
           male_analysis=analysis,female_analysis=analysis,
           male_age_col=age_col,female_age_col=age_col) %>% 
    select(PheWAS_ID,male,female,analysis,male_analysis,female_analysis,age_col,male_age_col,female_age_col)
  
  all_phenos_complete <- data.frame(PheWAS_ID=c(sex_phenos$PheWAS_ID,sex_phenos$male,sex_phenos$female), analysis=c(sex_phenos$analysis,sex_phenos$male_analysis,sex_phenos$female_analysis), age_col=c(sex_phenos$age_col,sex_phenos$male_age_col,sex_phenos$female_age_col))
  
  all_phenos <- all_phenos_complete %>% 
    mutate(include=1) %>% 
    bind_rows(age_phenos)
} else {
  age_phenos <- PheWAS_manifest %>% 
    filter(age_col!="age") %>% 
    filter(included_in_analysis==1) %>% 
    mutate(age_col_complete = ifelse(age_col=="named",paste0(PheWAS_ID,"_age"),"age")) %>% 
    select(PheWAS_ID=age_col_complete)
  
  sex_phenos_list <- fread(arguments$sex_include_file) %>% 
    left_join(PheWAS_manifest) %>% 
    mutate(male=paste0(PheWAS_ID,"_male"),
           female=paste0(PheWAS_ID,"_female"),
           male_analysis=analysis,female_analysis=analysis,
           male_age_col=age_col,female_age_col=age_col) %>% 
    select(PheWAS_ID,male,female,analysis,male_analysis,female_analysis,age_col,male_age_col,female_age_col)
  
  sex_phenos <- data.frame(PheWAS_ID=c(sex_phenos_list$male,sex_phenos_list$female), analysis=c(sex_phenos_list$male_analysis,sex_phenos_list$female_analysis), age_col=c(sex_phenos_list$male_age_col,sex_phenos_list$female_age_col))
  
  all_phenos <- PheWAS_manifest %>% 
    filter(included_in_analysis==1) %>% 
    bind_rows(sex_phenos) %>% 
    mutate(include=1) %>% 
    select(PheWAS_ID,analysis,age_col,include) %>% 
    bind_rows(age_phenos)
  
}
# create required folder
if(!dir.exists(paste0(arguments$analysis_folder,"/","plink_results"))) {
  dir.create(paste0(arguments$analysis_folder,"plink_results"))
}
# split phenotypes
binary_ID <- all_phenos %>% 
  filter(analysis=="binary") %>% 
  pull(PheWAS_ID)
quant_ID <- all_phenos %>% 
  filter(analysis=="quant") %>% 
  pull(PheWAS_ID)
# psam required to match sample id for testing
psam <- fread(paste0(association_variants,".psam"))
# edited to match psam and required col names
edited_covartiate <- fread(arguments$covariate) %>% 
  filter(IID %in% psam$IID)
fwrite(edited_covartiate,paste0(arguments$analysis_folder,"/","edited_covariate"),na = NA,quote = F,sep = "\t")
# location for plink2 commands
covariates <- paste0(arguments$analysis_folder,"/","edited_covariate")
plink_exe <- arguments$plink_exe
if(!is.null(arguments$model)){
  model <- arguments$model
  acceptable_input <- c("recessive","dominant","hetonly","genotypic","hethom")
  if(!model %in% acceptable_input){
    model <- ""
  }
} else {
  model <- ""
}

# Run functions -----------------------------------------------------------
all_results <- mapply(association_tests,groups,phenotype_files,SIMPLIFY = F, USE.NAMES = T)
#save results
if(arguments$split_analysis) {
} else {
  saveRDS(all_results,paste0(arguments$analysis_folder,"/","association_results","/",arguments$analysis_name,"_association_results_list.gz"))
}